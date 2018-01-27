import sys

import gffutils

#Read gtf
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

#Extract transcripts from gtf
def extractFromGTF(gtf):
    transcripts = {}
    strand = ""

    for g in gtf.features_of_type('gene'):
        strand = g.strand
        for tr in gtf.children(g, featuretype='transcript', order_by='start'):
            transcript = []
            for ex in gtf.children(tr, featuretype='exon', order_by='start'):
                transcript += [ex.start]
                transcript += [ex.end]
            transcripts[tr.id.split(".")[0]] = transcript
    return strand, transcripts

def printEv(chrom, t, p1, p2):
    print("{} {} {} {}".format(chrom, t, p1, p2))

def main():
    truthPath = sys.argv[1]
    gtfPath = sys.argv[2]
    bedPath = sys.argv[3]

    strand, transcripts = extractFromGTF(openGTF(gtfPath))

    novelTranscripts = []
    for line in open(bedPath).readlines():
        tID = line.strip("\n").split("\t")[-1].split("--")[0]
        novelTranscripts.append(transcripts[tID])

    for line in open(truthPath).readlines():
        line = line.strip("\n").split(" ")
        if line[1] == "ES":
            chrom, evT, p0, p1 = line
            p0,p1 = int(p0)-1, int(p1)+1
            isValidEvent = False
            for T in novelTranscripts:
                if p0 in T:
                    i = T.index(p0)
                    if T[i+1] == p1:
                        isValidEvent = True
                        break
            if isValidEvent:
                printEv(chrom, evT, p0+1, p1-1)
        elif line[1] == "IR":
            chrom, evT, p0, p1 = line
            p0,p1 = int(p0)-1, int(p1)+1
            isValidEvent = False
            for T in novelTranscripts:
                if p0 in T:
                    i = T.index(p0)
                    if T[i+1] == p1:
                        isValidEvent = True
                        break
            if isValidEvent:
                printEv(chrom, evT, p0+1, p1-1)
        elif line[1] == "A3":
            chrom, evT, p0, p1, p2 = line
            p0,p1,p2 = int(p0), int(p1), int(p2)
            if strand == '+':
                P0,P1 = 0,0
                isValidEvent = False
                for T in novelTranscripts:
                    if p2 in T and p0 in T:
                        isValidEvent = True
                        evT = evT+"e"
                        P0,P1 = p0,p2
                        break
                    elif p2 in T and p1 in T:
                        isValidEvent = True
                        evT = evT+"i"
                        P0,P1 = p1,p2
                        break
                if isValidEvent:
                    printEv(chrom, evT, P0-1, P1+1)
            else:
                P0,P1 = 0,0
                isValidEvent = False
                for T in novelTranscripts:
                    if p2 in T and p0 in T:
                        isValidEvent = True
                        evT = evT+"e"
                        P0,P1 = p0,p2
                        break
                    elif p2 in T and p1 in T:
                        isValidEvent = True
                        evT = evT+"i"
                        P0,P1 = p1,p2
                        break
                if isValidEvent:
                    printEv(chrom, evT, P0+1, P1-1)
        elif line[1] == "A5":
            chrom, evT, p0, p1, p2 = line
            p0,p1,p2 = int(p0), int(p1), int(p2)
            if strand == '+':
                P0,P1 = 0,0
                isValidEvent = False
                for T in novelTranscripts:
                    if p0 in T and p2 in T:
                        isValidEvent = True
                        evT = evT+"i"
                        P0,P1 = p0,p2
                        break
                    elif p1 in T and p2 in T:
                        isValidEvent = True
                        evT = evT+"e"
                        P0,P1 = p1,p2
                        break
                if isValidEvent:
                    printEv(chrom, evT, P0+1, P1-1)
            else:
                P0,P1 = 0,0
                isValidEvent = False
                for T in novelTranscripts:
                    if p0 in T and p2 in T:
                        isValidEvent = True
                        evT = evT+"i"
                        P0,P1 = p0,p2
                        break
                    elif p1 in T and p2 in T:
                        isValidEvent = True
                        evT = evT+"e"
                        P0,P1 = p1,p2
                        break
                if isValidEvent:
                    printEv(chrom, evT, P0-1, P1+1)

if __name__ == '__main__':
    main()
