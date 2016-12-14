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
  
  for rec in GFF.parse(in_handle, limit_info=dict(gff_type = ["gene", "transcript", "exon"])):
    chromosome = rec.id
    ref_genomic = genomic_dict[chromosome]
    for gene in rec.features:
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
            exs_pos.append([int(start),int(end)])

  in_handle.close()
  
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
        if exs_pos[i][1] <= exs_pos[j][0]:
          edg_string += "{},{}\n".format(i+1,j+1)
    
  f = open("./tmp/edges", "w")
  f.write("{}".format(edg_string))
  f.close()

def savePatternFile(pattern):
  f = open("./tmp/P.fa", "w")
  f.write("> P\n{}\n".format(pattern))
  f.close()
  
  f = open("patterns", "a")
  f.write("{}\n".format(pattern))
  f.close()


def main(genomic, gene_annotation, rna_seqs, L, k):
  os.system("rm T.fa*")
  os.system("rm P.fa")
  os.system("rm e_lens")
  os.system("rm mems")
  os.system("rm patterns")
  os.system("touch patterns")
  
  extractInfo(genomic, gene_annotation)
  
  patterns_fa = list(SeqIO.parse("{}".format(rna_seqs), "fasta"))
  c = 0
  log = open("avanzamento.log", "w")
  for pattern in patterns_fa:
    c += 1
    p = pattern.seq
    savePatternFile(str(p))
    p_len = len(p)
    os.system("./bMEM/backwardMEM -l={} ./tmp/T.fa ./tmp/P.fa > ./tmp/mems".format(L))
    print("../bin/main ./tmp/mems ./tmp/e_lens ./tmp/edges {} {} {}".format(L, k, p_len))
    os.system("../bin/main ./tmp/mems ./tmp/e_lens ./tmp/edges {} {} {}".format(L, k, p_len))
    if c % 100 == 0:
       log.write("{0}\n".format(c))
       log.flush()
  os.system("rm stopwatch.txt")
  os.system("rm lock.txt")
    

if __name__ == '__main__':
  main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
