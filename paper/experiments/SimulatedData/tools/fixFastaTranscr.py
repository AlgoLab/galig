from Bio import SeqIO
import sys

file_in = sys.argv[1]

file_out = sys.argv[2]

fasta_sequences = SeqIO.parse(open(file_in), 'fasta')

with open(file_out, "w") as out:
    for fasta in fasta_sequences:
        name, seq = fasta.id, fasta.seq.tostring()
        newname = name.split("|")[1]
        out.write(">" + newname + "\n" + seq + "\n")
