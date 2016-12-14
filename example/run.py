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

  exons_number = 0
  T = "|"
  e_lens = []
  added_exs = []
  
  for rec in GFF.parse(in_handle, limit_info=dict(gff_type = ["gene", "transcript", "exon"])):
    chromosome = rec.id
    ref_genomic = genomic_dict[chromosome]
    for gene in rec.features:
      for transcript in gene.sub_features:
        for exon in transcript.sub_features:
          if(not(exon.id in added_exs)):
            exons_number += 1
            start = exon.location.start
            end = exon.location.end
            strand = exon.location.strand
            exon_seq = ref_genomic[start:end].seq
            if strand < 0:
              exon_seq = reverseAndComplement(exon_seq)
            T += str(exon_seq) + "|"
            e_lens.append(str(len(str(exon_seq))))
            added_exs.append(exon.id)
  in_handle.close()
  
  f = open("./tmp/T.fa", "w")
  f.write(">T\n{}\n".format(T))
  f.close()
  
  f = open("./tmp/e_lens", "w")
  f.write("{}".format("\n".join(e_lens)))
  f.close()

def savePatternFile(pattern):
  f = open("./tmp/P.fa", "w")
  f.write("> P\n{}\n".format(pattern))
  f.close()


def main(genomic, gene_annotation, rna_seqs, L, k):
  os.system("rm T.fa*")
  os.system("rm P.fa")
  os.system("rm e_lens")
  os.system("rm mems")
  
  extractInfo(genomic, gene_annotation)
  
  patterns_fa = list(SeqIO.parse("{}".format(rna_seqs), "fasta"))
  
  for pattern in patterns_fa:
    p = pattern.seq
    
    f = open("patterns.fa")
    f.write("{}\n".foramt(str(p)))
    f.close()
    
    savePatternFile(str(p))
    p_len = len(p)
    os.system("./bMEM/backwardMEM -l={} ./tmp/T.fa ./tmp/P.fa > ./tmp/mems".format(L))
    os.system("./bin/main ./tmp/mems ./tmp/e_lens {} {} {}".format(L, k, p_len))
  
  os.system("rm stopwatch.txt")
  os.system("rm lock.txt")
    

if __name__ == '__main__':
  main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
