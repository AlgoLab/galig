import os, sys

from Bio import SeqIO

chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
        "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

def main():
    refPath = sys.argv[1]
    dirPath = sys.argv[2]

    print("* Reading FASTA")
    FullRef = SeqIO.index(refPath, "fasta")

    if not(os.path.isdir(dirPath)):
        os.mkdir(dirPath)

    print("* Splitting FASTA")
    for c in chrs:
        out = "{}/{}.fasta".format(dirPath, c)
        ref = FullRef[c]
        SeqIO.write(ref, out, "fasta")

if __name__ == '__main__':
    main()
