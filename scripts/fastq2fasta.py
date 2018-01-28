import sys, os

from Bio import SeqIO

def main():
    fastqPath = sys.argv[1]
    fastaPath = sys.argv[2]

    if not os.path.exists(fastqPath):
        print("No FASTQ file at", fastqPath)
        sys.exit(1)
    
    with open(fastqPath, "rU") as fastq:
        with open(fastaPath, "w") as fasta:
            sequences = SeqIO.parse(fastq, "fastq")
            SeqIO.write(sequences, fasta, "fasta")

if __name__ == '__main__':
    main()
