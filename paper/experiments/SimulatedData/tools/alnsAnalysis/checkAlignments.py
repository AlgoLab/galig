#!/usr/bin/python3

import sys, os, argparse
import glob

from Bio import SeqIO
import gffutils
import pysam

BaseFold = ""
RefsFold = ""
GTFsFold = ""
SamplesFold = ""
ResultsFold = ""

class Gene:
    def __init__(self, id, strand, start, end):
        self.id = id
        self.strand = strand
        self.start = start
        self.end = end
        self.transcripts = {}

    def addTranscript(self, transcript):
        self.transcripts[transcript.id] = transcript

class Transcript:
    def __init__(self, id, strand, start, end, exons):
        self.id = id
        self.strand = strand
        self.start = start
        self.end = end
        self.exons = exons
        self.length = sum([e-s+1 for (s,e) in self.exons])

#############
# GTF utils #
#############
def openGTF(gtfPath):
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath),
                                 keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(gtfPath,
                                 dbfn="{}.db".format(gtfPath),
                                 force=True, keep_order=True,
                                 disable_infer_genes=True,
                                 disable_infer_transcripts=True,
                                 merge_strategy='merge',
                                 sort_attribute_values=True)
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath), keep_order=True)
    return gtf

def extractGene(gtf, genes):
    for g in gtf.features_of_type('gene'):
        gene = Gene(g.id, g.strand, g.start, g.end)
        for tr in gtf.children(g, featuretype='transcript', order_by='start'):
            exons = list(gtf.children(tr, featuretype='exon', order_by='start'))
            transcript = Transcript(tr.id, tr.strand, tr.start, tr.end, [(ex.start, ex.end) for ex in exons])
            gene.addTranscript(transcript)
        genes[g.id] = gene

def buildTr2GeneMap(genes):
    tr2gene = {}
    for geneID,gene in genes.items():
        for transcript in gene.transcripts:
            tr2gene[transcript] = geneID
    return tr2gene

'''
 Return True if the position pos is inside the gene considered in the analysis
 (only the gene with AS events are considered)
'''
def isInsideCoveredGenes(pos, genes, coveredGenes):
    for gene in [genes[geneID] for geneID in coveredGenes]:
        if pos in range(gene.start, gene.end):
            return True
    return False

def getCoveredPositions(transcript, fragmentStart, fragmentEnd, readStrand, readLen):
    strand = transcript.strand
    startOffset = 0
    endOffset = 0
    coveredPositions = []
    pos = 0    

    if fragmentStart<=0:
        startOffset = abs(fragmentStart)
    if fragmentEnd>=transcript.length:
        endOffset = fragmentEnd - transcript.length

    for (s,e) in transcript.exons if strand is '+' else transcript.exons[::-1]:
        if pos > fragmentEnd:
            break
        for p in range(s,e+1) if strand is '+' else list(range(s,e+1))[::-1]:
            if fragmentStart <= pos <= fragmentEnd:
                if strand is '+':
                    coveredPositions.append(p-1)
                else:
                    coveredPositions.insert(0, p-1)
            pos+=1

    if strand is '+' and readStrand is '+':
        if startOffset > 0:
            externalPositions = list(range(transcript.start-startOffset-1, transcript.start-1))
            # print("## + + 1",len(externalPositions), len(coveredPositions), externalPositions[-3:], coveredPositions[:3])
            coveredPositions = externalPositions + coveredPositions
        if endOffset > 0:
            externalPositions = list(range(transcript.end, transcript.end + 100 - len(coveredPositions)))
            # print("## + + 2", len(coveredPositions), len(externalPositions), coveredPositions[-3:], externalPositions[:3])
            coveredPositions = coveredPositions + externalPositions
        return coveredPositions[:100]
    elif strand is '+' and readStrand is '-':
        if startOffset > 0:
            externalPositions = list(range(transcript.start - 100 + len(coveredPositions), transcript.start-1))
            # print("## + - 1", len(externalPositions), len(coveredPositions), externalPositions[-3:], coveredPositions[:3])
            coveredPositions = externalPositions + coveredPositions
        if endOffset > 0:
            externalPositions = list(range(transcript.end, transcript.end + endOffset))
            # print("## + - 2", len(coveredPositions), len(externalPositions), coveredPositions[-3:], externalPositions[:3])
            coveredPositions = coveredPositions + externalPositions
        return coveredPositions[-100:]

    elif strand is '-' and readStrand is '+':
        if startOffset > 0:
            externalPositions = list(range(transcript.end, transcript.end + startOffset))
            # print("## - + 1", len(coveredPositions), len(externalPositions), coveredPositions[-3:], externalPositions[:3])
            coveredPositions = coveredPositions + externalPositions
        if endOffset > 0:
            externalPositions = list(range(transcript.start - 100 + len(coveredPositions), transcript.start - 1))
            # print("## - + 2", len(externalPositions), len(coveredPositions), externalPositions[-3:], coveredPositions[:3])
            coveredPositions = externalPositions + coveredPositions
        return coveredPositions[-100:]
    else:
        if startOffset > 0:
            externalPositions = list(range(transcript.end, transcript.end + 100 - len(coveredPositions)))
            # print("##", len(coveredPositions), len(externalPositions), coveredPositions[-3:], externalPositions[:3])
            coveredPositions = coveredPositions + externalPositions
        if endOffset > 0:
            externalPositions = list(range(transcript.start - endOffset - 1 - 1, transcript.start - 1))
            # print("##", len(externalPositions), len(coveredPositions), externalPositions[-3:], coveredPositions[:3])
            coveredPositions = externalPositions + coveredPositions
        return coveredPositions[:100]

def checkAln(aln, readLen, gene, alnsNumDict, alnsMismatchesDict, alnsClipsDict, tool):
    seqName = aln.query_name
    trID = seqName.split(':')[2]
    fragmentStart = int(seqName.split(':')[5])
    fragmentEnd = int(seqName.split(':')[6])
    readStrand = '+' if seqName.split(':')[7][0] is 'S' else '-'

    # Number of alignments per read
    alnsNumDict[seqName] = alnsNumDict[seqName] + 1 if seqName in alnsNumDict else 1

    # Mismatches per alignment
    nm = aln.get_cigar_stats()[0][10]
    alnsMismatchesDict[nm] = alnsMismatchesDict[nm] + 1 if nm in alnsMismatchesDict else 1

    # Clips per alignment
    clips = aln.get_cigar_stats()[0][4] + aln.get_cigar_stats()[0][5]
    alnsClipsDict[clips] = alnsClipsDict[clips] + 1 if clips in alnsClipsDict else 1

    # Is a good alignment?
    coveredPositions = []
    for blockList in [list(range(block[0], block[1])) for block in aln.get_blocks()]:
        coveredPositions.extend(blockList)
    coveredPositions = set(coveredPositions)
    realCoveredPositions = set(getCoveredPositions(gene.transcripts[trID], fragmentStart, fragmentEnd, readStrand, readLen))
    
    clips = aln.get_cigar_stats()[0][4]
    insertions = aln.get_cigar_stats()[0][1]
    deletions = aln.get_cigar_stats()[0][2]
    exactCovered = len(coveredPositions & realCoveredPositions)
    wrongCovered = max(0, readLen - exactCovered - clips - insertions - deletions)

    if wrongCovered is 0:
        return 1
    if exactCovered > 0:
        return 2
    return 0

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', dest='base_folder', help='Data folder', required=True)
    parser.add_argument('-c', '--chr', dest='chrom', help='Chromosome', required=True)
    parser.add_argument('-s', '--sample', dest='sample', help='Sample', required=True)

    args = parser.parse_args()

    BaseFold = "/home/prj_gralign/data/BMC-exp-sim"
    RefsFold = os.path.join(args.base_folder, "Refs")
    GTFsFold = os.path.join(args.base_folder, "RealAnnotations")
    SamplesFold = os.path.join(args.base_folder, "Samples")
    ResultsFold = os.path.join(args.base_folder, "Results")
    outFold = os.path.join(ResultsFold, args.sample, args.chrom, "alignments_analysis")

    fastaPath = os.path.join(RefsFold, "{}.fasta".format(args.chrom))
    RefSeq = list(SeqIO.parse(fastaPath, "fasta"))[0]

    gtfPaths = glob.glob(os.path.join(GTFsFold, args.chrom, '*.gtf'))
    genes = {}
    for gtfPath in gtfPaths:
        gtf = openGTF(gtfPath)
        extractGene(gtf, genes)
    tr2gene = buildTr2GeneMap(genes)
    
    samplePath = os.path.join(SamplesFold, args.sample, "{}.fasta".format(args.chrom))
    readsPerTranscript = {}
    readsLen = {}
    for read in SeqIO.parse(samplePath, "fasta"):
        trID = read.id.split(':')[2]
        readsPerTranscript[trID] = readsPerTranscript[trID] + 1 if trID in readsPerTranscript else 1
        readsLen[read.id.split('/')[0]] = len(read.seq)

    coveredGenes = set()
    asgalTotal = 0
    asgalGood = 0
    asgalPartially = 0
    asgalAlnsNum = {}
    asgalAlnsMismatches = {}
    asgalAlnsClips = {}
    asgalSamPaths = glob.glob(os.path.join(ResultsFold, args.sample, args.chrom, 'asgal_annot', '*/*.sam'))
    for asgalSamPath in asgalSamPaths:
        geneID = os.path.splitext(os.path.basename(asgalSamPath))[0]
        coveredGenes.add(geneID)
        asgalSam = pysam.AlignmentFile(asgalSamPath, "r")
        for aln in asgalSam.fetch():
            if not aln.is_secondary:
                readLen = readsLen[aln.query_name.split('/')[0]]
                seqName = aln.query_name
                trID = seqName.split(':')[2]
                geneID = tr2gene[trID]
                isGood = checkAln(aln, readLen, genes[geneID], asgalAlnsNum, asgalAlnsMismatches, asgalAlnsClips, 'ASGAL')
                if isGood == 1:
                    asgalGood += 1
                elif isGood == 2:
                    asgalPartially += 1
                asgalTotal += 1
        asgalSam.close()

    totalConsideredReads = sum([readsPerTranscript[trID] for gene in coveredGenes for trID in genes[gene].transcripts])

    starTotal = 0
    starGood = 0
    starPartially = 0
    starAlnsNum = {}
    starAlnsMismatches = {}
    starAlnsClips = {}
    starBamPath = os.path.join(ResultsFold, args.sample, args.chrom, 'star-aln', '{}_Aligned.sortedByCoord.out.bam'.format(args.chrom))
    starBam = pysam.AlignmentFile(starBamPath, "rb")
    for aln in starBam.fetch():
        pos = aln.reference_start + 1
        if not aln.is_secondary and isInsideCoveredGenes(pos, genes, coveredGenes):
            readLen = readsLen[aln.query_name] #STAR remove trailing /int
            seqName = aln.query_name
            trID = seqName.split(':')[2]
            geneID = tr2gene[trID]
            isGood = checkAln(aln, readLen, genes[geneID], starAlnsNum, starAlnsMismatches, starAlnsClips, 'STAR')
            if isGood == 1:
                starGood += isGood
            elif isGood == 2:
                starPartially += 1
            starTotal += 1
    starBam.close()

    if not os.path.isdir(outFold):
        os.mkdir(outFold)

    with open(os.path.join(outFold, "alignmentsAccuracy.csv".format(args.chrom)), 'w') as f:
        f.write("Chrom,Tool,TotalReads,GoodAlns,PartAlns,BadAlns,TotalAlns\n")
        f.write("{},{},{},{},{},{},{},{}\n".format(args.sample, args.chrom, 'asgal', totalConsideredReads, asgalGood, asgalPartially, asgalTotal-asgalGood-asgalPartially, asgalTotal))
        f.write("{},{},{},{},{},{},{},{}\n".format(args.sample, args.chrom, 'star', totalConsideredReads, starGood, starPartially, starTotal-starGood-starPartially, starTotal))

    with open(os.path.join(outFold, "alignmentsStatistics.csv".format(args.chrom)), 'w') as f:
        f.write("Chrom,Tool,Type,TotalReads,NumberOfType,NumberAlignments\n")
        # Number of unique alignments
        hist = {}
        for seqName,w in asgalAlnsNum.items():
            hist[w] = hist[w]+1 if w in hist else 1
        for k in sorted(hist):
            f.write("{},{},{},{},{},{}\n".format(args.sample, args.chrom, 'ASGAL', 'uniqueAl', k, hist[k]))
        hist = {}
        for seqName,w in starAlnsNum.items():
            hist[w] = hist[w]+1 if w in hist else 1
        for k in sorted(hist):
            f.write("{},{},{},{},{},{}\n".format(args.sample, args.chrom, 'STAR', 'uniqueAl', k, hist[k]))

        # Mismatches
        for k in sorted(asgalAlnsMismatches):
            f.write("{},{},{},{},{},{}\n".format(args.sample, args.chrom, 'ASGAL', 'mismatches',  k, asgalAlnsMismatches[k]))
        for k in sorted(starAlnsMismatches):
            f.write("{},{},{},{},{},{}\n".format(args.sample, args.chrom, 'STAR', 'mismatches', k, starAlnsMismatches[k]))

        #Clips
        for k in sorted(asgalAlnsClips):
            f.write("{},{},{},{},{},{}\n".format(args.sample, args.chrom, 'ASGAL', 'clips', k, asgalAlnsClips[k]))
        for k in sorted(starAlnsClips):
            f.write("{},{},{},{},{},{}\n".format(args.sample, args.chrom, 'STAR', 'clips', k, starAlnsClips[k]))

if __name__ == '__main__':
    main()
