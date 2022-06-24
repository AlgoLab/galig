import sys
from Bio import SeqIO


def main():
    fa_path = sys.argv[1]
    idxs = set()
    for record in SeqIO.parse(fa_path, "fasta"):
        if record.id in idxs:
            continue
        idxs.add(record.id)
        SeqIO.write(record, sys.stdout, "fasta")


if __name__ == "__main__":
    main()
