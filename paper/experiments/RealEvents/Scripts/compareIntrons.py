import sys

import pysam

def main():
    eventsPath = sys.argv[1]
    bamPath = sys.argv[2]
    chrom = sys.argv[3]
    gene = sys.argv[4]

    FoundIntrons = set()
    for line in open(eventsPath).readlines():
        if line[0:4] == "Type":
            continue
        line = line.strip("\n").split(",")
        p1,p2 = int(line[1]),int(line[2])
        FoundIntrons.add((p1,p2))

    RealIntrons = set()
    bam = pysam.AlignmentFile(bamPath, "rb")
    for record in bam.fetch():
        if "N" in record.cigarstring:
            Is = []
            i = 0
            for (f,_) in record.cigartuples:
                if f==3:
                    Is.append(i)
                if f == 0:
                    i+=1
            blocks = record.get_blocks()
            for i in Is:
                RealIntrons.add((blocks[i-1][1]+1, blocks[i][0]))
    print(chrom, gene, len(FoundIntrons & RealIntrons), len(FoundIntrons), len(RealIntrons), sep=",")

if __name__ == '__main__':
    main()
