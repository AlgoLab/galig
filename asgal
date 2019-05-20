#!/usr/bin/python3

import sys, os
import argparse
import time
import subprocess
import glob

import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffutils
import pysam

# CONST #############################################################################################################
refsFold = "refs"
annosFold = "annos"
samplesFold = "samples"
asgalFold = "ASGAL"
salmonIndexFold = "salmon/salmon_index"
salmonOutFold = "salmon/salmon_out"
logsFold = "logs"
bar_length = 50

# UTILITIES #########################################################################################################
'''
Functions to:
 - get current time
 - print to stderr with current time
 - print error and exit
 - print progress bar
'''
def getTime():
    return time.strftime('[ %b %d, %Y - %l:%M:%S%p ]')

def eprint(*args, **kwargs):
    print(getTime(), *args, file=sys.stderr, **kwargs)

def eprint_error(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    sys.exit(1)

def print_bar(BAR, i, n):
    sys.stderr.write("[{}] {}/{}\r".format(''.join(BAR), min(i,n), n))
    sys.stderr.flush()

# REFERENCE #########################################################################################################
'''
Function to split input reference in more references,
one for each chromosome.
'''
def splitReference(args):
    eprint("Splitting input reference...")
    outFold = os.path.join(args.outputPath, refsFold)
    if not os.path.isdir(outFold):
        os.makedirs(outFold)
    for record in SeqIO.parse(open(args.refPath), "fasta"):
        outFasta = open(os.path.join(outFold, "{}.fa".format(record.id)), "w")
        SeqIO.write(record, outFasta, "fasta")
        outFasta.close()
    eprint("Done.")

# GTF ###############################################################################################################
'''
Function to open a GTF file
'''
def openGTF(gtfPath, verbose=True):
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath),
                                 keep_order=True)
    except ValueError:
        if verbose:
            eprint("Indexing...")
        gtf = gffutils.create_db(gtfPath,
                                 dbfn="{}.db".format(gtfPath),
                                 force=True, keep_order=True,
                                 disable_infer_genes=True,
                                 disable_infer_transcripts=True,
                                 merge_strategy='merge',
                                 sort_attribute_values=True)
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath), keep_order=True)
    return gtf
'''
Function to:
 1. split input GTF in smaller GTFs, one for each gene (only if multi mode)
 2. build dictionary chr->genes and transcript->gene
'''
def splitAnnotation(args):
    eprint("Opening input annotation...")
    gtf = openGTF(args.annoPath)

    if args.multiMode:
        eprint("Splitting input annotation...")
        outFold = os.path.join(args.outputPath, annosFold)
        if not os.path.isdir(outFold):
            os.makedirs(outFold)
    else:
        eprint("Reading input annotation...")

    genes = list()
    chr_genes_dict = {}
    tr_gene_dict = {}

    nGenes = len(list(gtf.features_of_type('gene')))
    i = 0
    BAR = [' ' for i in range(0,bar_length)]
    print_bar(BAR, 0, nGenes)
    for gene in gtf.features_of_type('gene'):
        chrom = gene.seqid
        geneID = gene.id.split('.')[0]
        genes.append(geneID)
        chr_genes_dict[chrom] = chr_genes_dict[chrom] + 1 if chrom in chr_genes_dict else 1
        if args.multiMode:
            outPath = os.path.join(outFold, "{}.gtf".format(geneID))
            with open(outPath, 'w') as out:
                out.write(str(gene) + "\n")
                for transcript in gtf.children(gene, featuretype='transcript', order_by='start'):
                    transcriptID = transcript.id.split('.')[0]
                    tr_gene_dict[transcriptID] = geneID
                    out.write(str(transcript) + "\n")
                    for exon in gtf.children(transcript, featuretype='exon', order_by='start'):
                        out.write(str(exon) + "\n")
        i+=1
        if i >= nGenes-1:
            BAR = ['#' for i in range(0,bar_length)]
            i = nGenes - 1
        else:
            index = int(i*bar_length/nGenes)
            BAR = ['#' for i in range(0, index+1)] + [' ' for i in range(index+1,bar_length)]
        print_bar(BAR, i+1, nGenes)
    print("", file=sys.stderr)
    eprint("Done.")
    return genes, chr_genes_dict, tr_gene_dict

# SALMON ############################################################################################################
'''
Function to run Salmon and SAMtools.
'''
def runSalmon(args):
    if not os.path.isdir(os.path.join(args.outputPath, logsFold)):
        os.makedirs(os.path.join(args.outputPath, logsFold))
    salmonIndexLog = os.path.join(args.outputPath, logsFold, "salmon_index.log")
    salmonQuantLog = os.path.join(args.outputPath, logsFold, "salmon_quant.log")
    samtoolsLog = os.path.join(args.outputPath, logsFold, "samtools.log")
    salmonSam = os.path.join(args.outputPath, "salmon", "salmon.sam")
    salmonTmpBam = os.path.join(args.outputPath, "salmon", "salmon.sort.bam")
    salmonBam = os.path.join(args.outputPath, "salmon", "salmon.bam")

    salmon = "salmon"
    try:
        subprocess.call(salmon,
                        stdout=open("/dev/null", 'w'),
                        stderr=open("/dev/null", 'w'))
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            salmon = os.path.join(WP, "salmon", "bin", "salmon")

    salmonIndex = os.path.join(args.outputPath, salmonIndexFold)
    if not os.path.isdir(salmonIndex):
        os.makedirs(salmonIndex)
    salmonOut = os.path.join(args.outputPath, salmonOutFold)
    if not os.path.isdir(salmonOut):
        os.makedirs(salmonOut)
    eprint("Running Salmon indexing...")
    subprocess.run([salmon, "index",
                    "-t", args.transPath,
                    "-i", salmonIndex],
                   stdout=open(salmonIndexLog, 'w'), stderr=open(salmonIndexLog, 'w'))
    eprint("Done.")
    if args.sample2Path == '-':
        eprint("Running Salmon quasi-mapping on single-end sample...")
        subprocess.run([salmon, "quant",
                        "-i", salmonIndex,
                        "-l", "A",
                        "-r", args.sample1Path,
                        "-o", salmonOut,
                        ("--writeMappings="+salmonSam)],
                       stdout=open(salmonQuantLog, 'w'),
                       stderr=open(salmonQuantLog, 'w'))
        eprint("Done.")
    else:
        eprint("Running Salmon quasi-mapping on paired-end sample...")
        subprocess.run([salmon, "quant",
                        "-i", salmonIndex,
                        "-l", "A",
                        "-1", args.sample1Path,
                        "-2", args.sample2Path,
                        "-o", salmonOut,
			"--validateMappings",
                        ("--writeMappings="+salmonSam),
                        "--writeUnmappedNames"],
                       stdout=open(salmonQuantLog, 'w'),
                       stderr=open(salmonQuantLog, 'w'))
        eprint("Done.")
    # SAMtools post-processing
    eprint("Running SAMtools postprocessing...")
    subprocess.run(["samtools", "view",
                    "-Sb", salmonSam],
                   stdout=open(salmonTmpBam, 'w'),
                   stderr=open(samtoolsLog, 'w'))
    subprocess.run(["samtools", "sort",
                    salmonTmpBam],
                   stdout=open(salmonBam, 'w'),
                   stderr=open(samtoolsLog, 'a'))
    subprocess.run(["samtools", "index",
                    salmonBam],
                   stdout=open(samtoolsLog, 'a'),
                   stderr=open(samtoolsLog, 'a'))
    os.remove(salmonSam)
    os.remove(salmonTmpBam)
    eprint("Done.")

# Salmon Post-Processing ############################################################################################
'''
Function to parse the unmapped_names.txt file produced by Salmon
when run with paired-end sample.
'''
def parse_unmapped_file(unmapped_filename, sample_filenames):
    unmapped_tags = set()
    with open(unmapped_filename) as unmapped_file:
        for line in unmapped_file:
            tag, unmap_type = line.strip('\n').split(' ')
            if unmap_type in ['m1', 'm2']:
                file_suff = "2" if unmap_type == 'm1' else "1"
                unmapped_tags.add(tag[:-1] + file_suff)

    unmapped_reads = {}
    for sample_fn in sample_filenames:
        eprint("Parsing {} to retrieve reads...".format(sample_fn))
        isZipped = False
        isFastq = False
        fnames = sample_fn.split('.')
        if fnames[-1] == 'gz':
            isZipped = True
            if fnames[-2] in ['fastq', 'fq']:
                isFastq = True
            elif fnames[-2] in ['fasta', 'fa']:
                isFastq = False
            else:
                pass
        else:
            if fnames[-1] in ['fastq', 'fq']:
                isFastq = True
            elif fnames[-1] in ['fasta', 'fa']:
                isFastq = False
            else:
                pass

        file_handle = gzip.open(sample_fn, "rt") if isZipped else open(sample_fn, "r")
        fType = "fastq" if isFastq else "fasta"
        for record in SeqIO.parse(file_handle, fType):
            if record.id in unmapped_tags:
                unmapped_reads[record.id[:-2]] = record.seq
    eprint("Done.")
    return unmapped_reads

'''
Given a dict {gene_idx : [read_idxs]} and a dict {read_idx : fasta_record},
save a fasta file for each gene_idx containing the corresponding reads
'''
def save_to_fasta(gene_reads_dict, reads, out_dir):
    for gene_idx, read_idxs in gene_reads_dict.items():
        file_handle = open(os.path.join(out_dir, "{}.fa".format(gene_idx)), 'a')
        records = []
        for read_idx in read_idxs:
            records.append(reads[read_idx])
        SeqIO.write(records, file_handle, "fasta")
        file_handle.close()

'''
Function to split Salmon BAM into smaller fastas,
one for each gene covered by at least one read
'''
def split_bam(aln_filename, genes, tr_gene_dict, out_dir, unmapped_reads):
    chunk_size = 250000
    gene_reads_dict = {}
    reads = {}
    all_reads = {}
    with pysam.AlignmentFile(aln_filename, 'rb') as aln_file:
        for aln in aln_file:
            if aln.is_secondary or aln.is_unmapped:
                continue
            tr_name = aln_file.get_reference_name(aln.reference_id).split('.')[0]
            if tr_name not in tr_gene_dict:
                continue
            if tr_gene_dict[tr_name] in genes:
                gene_idx = tr_gene_dict[tr_name]
                possibly_not_uniq_read_idx = aln.query_name
                if possibly_not_uniq_read_idx not in all_reads:
                    all_reads[possibly_not_uniq_read_idx] = 1
                else:
                    all_reads[possibly_not_uniq_read_idx] += 1
                read_idx = possibly_not_uniq_read_idx + "." + str(all_reads[possibly_not_uniq_read_idx])
                seq = SeqRecord(Seq(aln.query_sequence), id=read_idx, name=read_idx, description="From Salmon alignment file")
                if read_idx not in reads:
                    reads[read_idx] = seq
                gene_reads_dict[gene_idx] = gene_reads_dict[gene_idx] | set([read_idx]) if gene_idx in gene_reads_dict else set([read_idx])
                if possibly_not_uniq_read_idx in unmapped_reads:
                    # Mate is unmapped
                    all_reads[possibly_not_uniq_read_idx] += 1
                    read_idx = possibly_not_uniq_read_idx + "." + str(all_reads[possibly_not_uniq_read_idx])
                    seq = SeqRecord(unmapped_reads[possibly_not_uniq_read_idx], id=read_idx, name=read_idx, description="From Salmon alignment file")
                    if read_idx not in reads:
                        reads[read_idx] = seq
                    gene_reads_dict[gene_idx] = gene_reads_dict[gene_idx] | set([read_idx])
            if len(reads) == chunk_size:
                save_to_fasta(gene_reads_dict, reads, out_dir)
                reads = {}
                gene_reads_dict = {}
    if len(reads) != 0:
        save_to_fasta(gene_reads_dict, reads, out_dir)

def splitSalmon(args, genes, tr_gene_dict):
    salmonBam = os.path.join(args.outputPath, "salmon", "salmon.bam")
    salmonOut = os.path.join(args.outputPath, salmonOutFold)

    outFold = os.path.join(args.outputPath, samplesFold)
    if not os.path.isdir(outFold):
        os.makedirs(outFold)

    unmapped_reads = {}
    if args.sample2Path != '-':
        # Paired-end
        eprint("Retrieving unmapped reads...")
        unmapped_reads = parse_unmapped_file(os.path.join(salmonOut, "aux_info", "unmapped_names.txt"),
                                             list([args.sample1Path, args.sample2Path]))
        eprint("Unmapped reads that will be remapped: {}".format(len(unmapped_reads)))
        eprint("Done.")
    eprint("Splitting Salmon BAM...")
    split_bam(salmonBam, genes, tr_gene_dict, outFold, unmapped_reads)
    eprint("Done.")

# ASGAL #############################################################################################################
'''
Function to run ASGAL on the considered genes.
'''
def runASGAL(args, genes, chr_genes_dict):
    samples_paired = []
    if args.multiMode:
        refs = os.path.join(args.outputPath, refsFold)
        annos = os.path.join(args.outputPath, annosFold)
        samples = glob.glob(os.path.join(os.path.join(args.outputPath, samplesFold), '*.fa'))
        outFold = os.path.join(args.outputPath, asgalFold)
        if not os.path.isdir(outFold):
            os.makedirs(outFold)
        logFold = os.path.join(args.outputPath, logsFold, "ASGAL")
        if not os.path.isdir(logFold):
            os.makedirs(logFold)
    else:
        gene = genes[0]
        ref = args.refPath
        anno = args.annoPath
        samples = [args.sample1Path]
        outFold = os.path.join(args.outputPath)
        if args.pairedEnd:
            # NOTE: [:-2] removes -1 and -2 from read names, these will be added later where necessary
            out = os.path.join(outFold, "{}-{}".format((os.path.basename(samples[0]).split('.')[0])[:-2], gene))
        else:
            out = os.path.join(outFold, "{}-{}".format(os.path.basename(samples[0]).split('.')[0], gene))
        logFold = os.path.join(args.outputPath)
        if args.pairedEnd:
            # NOTE: [:-2] removes -1 and -2 from read names, these will be added later where necessary
            log = os.path.join(outFold, "{}-{}.log".format((os.path.basename(samples[0]).split('.')[0])[:-2], gene))
        else:
            log = os.path.join(outFold, "{}-{}.log".format(os.path.basename(samples[0]).split('.')[0], gene))


    if args.pairedEnd:
        samples_paired = [args.sample2Path]

    if not args.pairedEnd:    # Single-end sample

        eprint("Running ASGAL on {} {}...".format(len(samples), 'gene' if len(samples) == 1 else 'genes'))
        BAR = [' ' for i in range(0,bar_length)]
        print_bar(BAR, 0, len(samples))
        i = 0

        for sample in samples:
            if args.multiMode:
                gene = os.path.basename(sample)[:-3]
                anno = os.path.join(annos, "{}.gtf".format(gene))
                gtf = openGTF(anno, verbose=False)
                chrom = list(gtf.features_of_type('gene'))[0].seqid
                ref = os.path.join(refs, "{}.fa".format(chrom))
                out = os.path.join(outFold, "{}".format(gene))
                log = os.path.join(logFold, "{}".format(gene))

            allevents_flag = "--allevents" if args.allevents else ""


            asgal_CMDs = {'align': ["{}/bin/SpliceAwareAligner".format(WP),
                                    "-g", ref,
                                    "-a", anno,
                                    "-s", sample,
                                    "-l", args.l,
                                    "-e", args.e,
                                    "-o", "{}.mem".format(out)],
                          'sam' : ["python3",
                                   "{}/scripts/formatSAM.py".format(WP),
                                   "-m", "{}.mem".format(out),
                                   "-g", ref,
                                   "-a", anno,
                                   "-e", args.e,
                                   "-o", "{}.sam".format(out)],
                          'events' : ["python3",
                                      "{}/scripts/detectEvents.py".format(WP),
                                      "-g", ref,
                                      "-a", anno,
                                      "-m", "{}.mem".format(out),
                                      "-o", "{}.events.csv".format(out),
                                      "-e", args.e,
                                      "-w", args.w]}
            if args.allevents:
                asgal_CMDs['events'].append("--allevents")

            subprocess.run(asgal_CMDs['align'],
                           stdout=open(log, 'w'),
                           stderr=open(log, 'w'))
            subprocess.run(asgal_CMDs['sam'],
                           stdout=open(log, 'w'),
                           stderr=open(log, 'w'))
            subprocess.run(asgal_CMDs['events'],
                           stdout=open(log, 'w'),
                           stderr=open(log, 'w'))
            i+=1
            index = int(i*bar_length/len(samples))
            if i == len(genes)-1:
                BAR = ['#' for i in range(0,bar_length)]
            else:
                BAR = ['#' for i in range(0, index+1)] + [' ' for i in range(index+1,bar_length)]
                print_bar(BAR, i+1, len(samples))

    else:   # Paired-end sample

        flt = args.flt

        log = os.path.join(logFold, "{}-{}.log".format((os.path.basename(samples[0]).split('.')[0])[:-2], gene))
        eprint("Running ASGAL on {} {}...".format(len(samples), 'read-pair' if len(samples) == 1 else 'read-pairs'))
        BAR = [' ' for i in range(0,bar_length)]
        print_bar(BAR, 0, len(samples))
        i = 0

        for sample1,sample2 in zip(samples,samples_paired):
            if args.multiMode:
                gene = os.path.basename(sample)[:-3]
                anno = os.path.join(annos, "{}.gtf".format(gene))
                gtf = openGTF(anno, verbose=False)
                chrom = list(gtf.features_of_type('gene'))[0].seqid
                ref = os.path.join(refs, "{}.fa".format(chrom))
                out = os.path.join(outFold, "{}".format(gene))
                log = os.path.join(logFold, "{}".format(gene))

            allevents_flag = "--allevents" if args.allevents else ""

            asgal_CMDs_paired = {'align': ["{}/bin/SpliceAwareAlignerPaired".format(WP),
                                    "-g", ref,
                                    "-a", anno,
                                    "-1", sample1,
                                    "-2", sample2,
                                    "-f", flt,
                                    "-l", args.l,
                                    "-e", args.e,
                                    "-o", "{}".format(out)],
                                'sam' : ["python3",
                                   "{}/scripts/formatSAMPaired.py".format(WP),
                                   "-1", "{}-1.mem".format(out),
                                   "-2", "{}-2.mem".format(out),
                                   "-g", ref,
                                   "-a", anno,
                                   "-e", args.e,
                                   "-o", "{}".format(out)],
                                'events' : ["python3",
                                   "{}/scripts/detectEventsPaired.py".format(WP),
                                   "-g", ref,
                                   "-a", anno,
                                   "-1", "{}-1.mem".format(out),
                                   "-2", "{}-2.mem".format(out),
                                   "-o", "{}.events.csv".format(out),
                                   "-e", args.e,
                                   "-w", args.w]}
            if args.allevents:
                asgal_CMDs_paired['events'].append("--allevents")

            subprocess.run(asgal_CMDs_paired['align'],
                           stdout=open(log, 'w'),
                           stderr=open(log, 'w'))
            subprocess.run(asgal_CMDs_paired['sam'],
                           stdout=open(log, 'w'),
                           stderr=open(log, 'w'))
            subprocess.run(asgal_CMDs_paired['events'],
                           stdout=open(log, 'w'),
                           stderr=open(log, 'w'))
            i+=1
            index = int(i*bar_length/len(samples))
            if i == len(genes)-1:
                BAR = ['#' for i in range(0,bar_length)]
            else:
                BAR = ['#' for i in range(0, index+1)] + [' ' for i in range(index+1,bar_length)]
            print_bar(BAR, i+1, len(samples))

    print("", file=sys.stderr)
    eprint("Done.")

# INPUTS CHECK ######################################################################################################
'''
Function to check all inputs (if files exist, file extensions...).
'''
def checkInputs(args):
    # Reference
    if os.path.isfile(args.refPath):
        if args.refPath.split('.')[-1] not in ['fa', 'fasta']:
            eprint_error("\nUnknown extension for reference genome, is it a .fa or .fasta file? Halting...\n")
    else:
        eprint_error("\nReference genome {} not found. Halting...\n".format(args.refPath))

    # Annotation
    if os.path.isfile(args.annoPath):
        if args.annoPath[-3:] != 'gtf':
            eprint_error("\nUnknown extension for annotation, is it a .gtf file? Halting...\n")
    else:
        eprint_error("\nAnnotation {} not found. Halting...\n".format(args.annoPath))

    # Transcripts
    if args.transPath != '-':
        if os.path.isfile(args.transPath):
            if args.transPath.split('.')[-1] == 'gz':
                if args.transPath.split('.')[-2] not in ['fa', 'fasta']:
                    eprint_error("\nUnknown extension for transcripts file, is it a .fa.gz or .fasta.gz file? Halting...\n")
            else:
                if args.transPath.split('.')[-1] not in ['fa', 'fasta']:
                    eprint_error("\nUnknown extension for transcripts file, is it a .fa or .fasta file? Halting...\n")
        else:
            eprint_error("\nTranscripts file {} not found. Halting...\n".format(args.annoPath))

    # Samples
    sample = args.sample1Path
    if os.path.isfile(sample):
        if sample.split('.')[-1] == 'gz':
            if sample.split('.')[-2] not in ['fa', 'fasta', 'fq', 'fastq']:
                eprint_error("\nUnknown extension for sample 1, is it a .fa.gz/.fasta.gz/.fq.gz/.fastq.gz file? Halting...\n")
        else:
            if sample.split('.')[-1] not in ['fa', 'fasta', 'fq', 'fastq']:
                eprint_error("\nUnknown extension for sample 1, is it a .fa/.fasta/.fq/.fastq file? Halting...\n")
    else:
        eprint_error("\nSample 1 {} not found. Halting...\n".format(sample))

    sample = args.sample2Path
    if sample != '-':
        if os.path.isfile(sample):
            if sample.split('.')[-1] == 'gz':
                if sample.split('.')[-2] not in ['fa', 'fasta', 'fq', 'fastq']:
                    eprint_error("\nUnknown extension for sample 2, is it a .fa.gz/.fasta.gz/.fq.gz/.fastq.gz file? Halting...\n")
            else:
                if sample.split('.')[-1] not in ['fa', 'fasta', 'fq', 'fastq']:
                    eprint_error("\nUnknown extension for sample 2, is it a .fa/.fasta/.fq/.fastq file? Halting...\n")
        else:
            eprint_error("\nSample 2 {} not found. Halting...\n".format(sample))

    # L, e, w
    try:
        args.l = str(int(args.l))
    except ValueError:
        eprint_error("\n l value must be an integer. Halting...\n")
    try:
        args.e = str(int(args.e))
    except ValueError:
        eprint_error("\n e value must be an integer. Halting...\n")
    try:
        args.w = str(int(args.w))
    except ValueError:
        eprint_error("\n w value must be an integer. Halting...\n")

    # Other checks
    if args.sample2Path != '-' and args.transPath == '-' and not args.pairedEnd:
        eprint_error("\nIf you pass two samples, you have to pass also the transcripts of the genes. Halting...\n")
    if args.multiMode and args.transPath == '-':
        eprint_error("\nIn multi mode, you have to pass the transcripts of the genes. Halting...\n")

# main ##############################################################################################################
def main():
    global WP
    WP = '/'.join(sys.argv[0].split('/')[:-1])
    if len(WP) == 0:
        WP = '.'

    parser = argparse.ArgumentParser(description='ASGAL - Alternative Splicing Graph ALigner',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-g', '--genome', dest='refPath', help='Path to genome', required=True)
    parser.add_argument('-a', '--annotation', dest='annoPath', help='Path to annotation', required=True)
    parser.add_argument('-s', '--sample', dest='sample1Path', help='Path to sample (1)', required=True)
    parser.add_argument('-o', '--output', dest='outputPath', help='Path to output folder', required=True)

    parser.add_argument('-s2', '--sample2', dest='sample2Path', help='Path to sample (2)', required=False, default='-')
    parser.add_argument('-t', '--transcripts', dest='transPath', help='Path to transcripts', required=False, default='-')

    parser.add_argument('-l', '--L', dest='l', help='MEMs length', required=False, default="15")
    parser.add_argument('-e', '--erate', dest='e', help='Error rate', required=False, default="3")
    parser.add_argument('-w', '--support', dest='w', help='Minimum intron coverage', required=False, default="3")
    parser.add_argument('--allevents', dest='allevents', help='Use this if you want to detect all events, also annotated ones', required=False, action='store_true')
    parser.add_argument('--multi', dest='multiMode', help='Use this to run ASGAL in genome-wide mode', required=False, action='store_true')

    parser.add_argument('--paired', dest='pairedEnd', help='Use this to run ASGAL with paired-end reads', required=False, action='store_true')
    parser.add_argument('-f', '--flt', dest='flt', help='Fragment Library Type', required=False, default="ISF")

    args = parser.parse_args()

    checkInputs(args)

    if not os.path.isdir(args.outputPath):
        os.makedirs(args.outputPath)

    genes, chr_genes_dict, tr_gene_dict = splitAnnotation(args)
    if args.multiMode:
        splitReference(args)

    if args.multiMode:
        runSalmon(args)
        splitSalmon(args, genes, tr_gene_dict)

    runASGAL(args, genes, chr_genes_dict)

if __name__ == '__main__':
    main()
