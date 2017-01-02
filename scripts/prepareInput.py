import sys, os

from BCBio import GFF
from Bio import SeqIO

def reverseAndComplement(text):
  text = text[::-1]
  nucl_bases_dict = {"A":"T", "T":"A", "C":"G", "G":"C", "a":"t", "t":"a", "c":"g", "g":"c"}
  new_text = ""
  for elem in text:
    try:
      new_text += nucl_bases_dict[elem]
    except KeyError:
      new_text += elem
  return new_text

def extractInfo(genomic, gene_annotation):
  genomic_dict = SeqIO.index(genomic, "fasta")
  in_handle = open(gene_annotation)

  T = "|"
  e_lens = []
  added_exs = []
  exs_pos = []
  gene_name = ""
  
  for rec in GFF.parse(in_handle, limit_info=dict(gff_type = ["gene", "transcript", "exon"])):
    chromosome = rec.id
    ref_genomic = genomic_dict[chromosome]
    for gene in rec.features:
      gene_name = gene.id
      for transcript in gene.sub_features:
        for exon in transcript.sub_features:
          if(not(exon.id in added_exs)):
            start = exon.location.start
            end = exon.location.end
            strand = exon.location.strand
            exon_seq = ref_genomic[start:end].seq
            if strand < 0:
              exon_seq = reverseAndComplement(exon_seq)
            T += str(exon_seq) + "|"
            e_lens.append(str(len(str(exon_seq))))
            added_exs.append(exon.id)
            exs_pos.append("{},{}".format(int(start),int(end)))

  in_handle.close()

  f = open("./tmp/gene_name", "w")
  f.write(gene_name)
  f.close()
  
  f = open("./tmp/e_pos", "w")
  f.write("{}".format("\n".join(exs_pos)))
  f.close()
  
  f = open("./tmp/T.fa", "w")
  f.write(">T\n{}\n".format(T))
  f.close()
  
  f = open("./tmp/e_lens", "w")
  f.write("{}".format("\n".join(e_lens)))
  f.close()
  
  edg_string = ""
  for i in range(0,len(added_exs)):
    for j in range(0,len(added_exs)):
      if i != j:
        if int(exs_pos[i].split(",")[1]) <= int(exs_pos[j].split(",")[0]):
          edg_string += "{},{}\n".format(i+1,j+1)
    
  f = open("./tmp/edges", "w")
  f.write("{}".format(edg_string))
  f.close()

def main(genomic, gene_annotation):
  print("Extracting splicing graph (T.fa, e_lens, edges)...")
  extractInfo(genomic, gene_annotation)

if __name__ == '__main__':
  #Genomic, annotation, rna-seq reads
  main(sys.argv[1], sys.argv[2])
