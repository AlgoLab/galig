import sys

import gffutils

#Extract transcripts from gtf
def extractFromGTF(gtfPath):
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

    strand = ""
    for g in gtf.features_of_type('gene'):
        strand = g.strand
    return strand

def printEv(chrom, t, p1, p2):
    print("{} {} {} {}".format(chrom, t, p1, p2))

def main():
    truthPath = sys.argv[1]
    gtfPath = sys.argv[2]
    strand = extractFromGTF(gtfPath)

    for line in open(truthPath).readlines():
        line = line.strip("\n").split(" ")
        if line[1] == "ES":
            chrom, evT, p0, p1 = line
            p0,p1 = int(p0), int(p1)
            printEv(chrom, evT, p0, p1)
        elif line[1] == "IR":
            chrom, evT, p0, p1 = line
            p0,p1 = int(p0), int(p1)
            printEv(chrom, evT, p0, p1)
        elif line[1] == "A3":
            chrom, evT, p0, p1, p2 = line
            p0,p1,p2 = int(p0), int(p1), int(p2)
            if strand == '+':
                printEv(chrom, evT, p0-1, p2+1)
                printEv(chrom, evT, p1-1, p2+1)
            else:
                printEv(chrom, evT, p0+1, p2-1)
                printEv(chrom, evT, p1+1, p2-1)
        elif line[1] == "A5":
            chrom, evT, p0, p1, p2 = line
            p0,p1,p2 = int(p0), int(p1), int(p2)
            if strand == '+':
                printEv(chrom, evT, p0+1, p2-1)
                printEv(chrom, evT, p1+1, p2-1)
            else:
                printEv(chrom, evT, p0-1, p2+1)
                printEv(chrom, evT, p1-1, p2+1)

if __name__ == '__main__':
    main()
