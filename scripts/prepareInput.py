import sys, os

from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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

def extractInfo(genomic, gene_annotation, rna, out):
  genomic_dict = SeqIO.index(genomic, "fasta")
  in_handle = open(gene_annotation)

  T = "|"
  e_lens = []
  added_exs = []
  exs_pos = []
  chromosome = ""
  chromo_len = 0
  gene_name = ""
  real_edges = []
  real_edg_string = ""
  
  for rec in GFF.parse(in_handle, limit_info=dict(gff_type = ["gene", "transcript", "exon"])):
    chromosome = rec.id
    ref_genomic = genomic_dict[chromosome]
    #chromo_len = len(ref_genomic)
    for gene in rec.features:
      gene_name = gene.id
      for transcript in gene.sub_features:
        last_id = -1
        for exon in transcript.sub_features:
          if(not(exon.id in added_exs)):
            start = exon.location.start
            end = exon.location.end
            strand = exon.location.strand
            exon_seq = ref_genomic[start:end].seq
            #if strand < 0:
            #  exon_seq = reverseAndComplement(exon_seq)
            T += str(exon_seq) + "|"
            e_lens.append(str(len(str(exon_seq))))
            added_exs.append(exon.id)
            exs_pos.append("{},{}".format(int(start),int(end)))
          exs_index = added_exs.index(exon.id)+1
          if last_id != -1:
            if not([last_id, exs_index] in real_edges):
              real_edges.append([last_id, exs_index])
              real_edg_string += "{},{}\n".format(last_id, exs_index)
          last_id = exs_index
  
  in_handle.close()

  f = open("{}/tmp/Ps.fa".format(out), "w")
  for record in SeqIO.parse(rna, "fasta"):
    id1 = record.id
    seq1 = record.seq
    record = SeqRecord(seq1, id=id1, description="")
    SeqIO.write(record, f, "fasta")
    id2 = id1 + "'"
    seq2 = reverseAndComplement(seq1)
    record = SeqRecord(Seq(seq2), id=id2, description="")
    SeqIO.write(record, f, "fasta")
  f.close()

  f = open("{}/tmp/gene_info".format(out), "w")
  f.write("{}\n{}\n{}".format(chromosome, chromo_len, gene_name))
  f.close()
  
  f = open("{}/tmp/e_pos".format(out), "w")
  f.write("{}".format("\n".join(exs_pos)))
  f.close()
  
  f = open("{}/tmp/T.fa".format(out), "w")
  f.write(">T\n{}\n".format(T))
  f.close()
  
  f = open("{}/tmp/e_lens".format(out), "w")
  f.write("{}".format("\n".join(e_lens)))
  f.close()

  f = open("{}/tmp/real_edges".format(out), "w")
  f.write("{}".format(real_edg_string))
  f.close()

  exs_pos_list = []
  for elem in exs_pos:
    elem_list = elem.split(",")
    exs_pos_list.append([int(elem_list[0]), int(elem_list[1])])
  edg_string = ""
  for i in range(0,len(added_exs)):
    for j in range(0,len(added_exs)):
      if i != j and not([i+1,j+1] in real_edges):
        if exs_pos_list[i][1] <= exs_pos_list[j][0]:
          edg_string += "{},{}\n".format(i+1,j+1)
    
  f = open("./tmp/added_edges", "w")
  f.write("{}".format(edg_string))
  f.close()

def main(genomic, gene_annotation, rna, out):
  extractInfo(genomic, gene_annotation, rna, out)

if __name__ == '__main__':
  #Genomic, annotation, outs
  main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
