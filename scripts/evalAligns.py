import sys, os, re

from Bio import SeqIO
import gffutils

def main():
    sam_path = sys.argv[1]
    fa_path = sys.argv[2]

    gene = os.path.basename(sam_path)[:-4]

    TP = 0
    FP = 0
    FN = 0

    found_ids = []
    with open(sam_path, 'r') as sam:
        for line in sam:
            if line[0] != '@':
                spl_line = line.split('\t')
                identifier, strand, start, cigar = spl_line[0], int(spl_line[1])==0, int(spl_line[3]), spl_line[5]
                found_ids.append(identifier)
                regex = re.search(".*_(.+)$", identifier)
                real_start = int(regex.group(1))

                if start == real_start+1:
                    TP+=1
                else:
                    FP+=1

    real = 0
    for record in SeqIO.parse(fa_path, "fasta"):
        real+=1

    FN = real-TP
    if FN<0:
        FN=0

    P = TP/(TP+FP)
    R = TP/(TP+FN)

    print(gene, real, len(found_ids), TP, FP, FN, P, R, sep=",")

if __name__ == '__main__':
    main()
