import os, sys

from Bio import SeqIO
import gffutils

def openGTF(gtfPath):
    try:
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath),
                                 keep_order=True)
    except ValueError:
        gtf = gffutils.create_db(gtfPath,
                                 dbfn="{}.db".format(gtfPath),
                                 force=True, keep_order=True,
                                 disable_infer_genes=True,
                                 disable_infer_transcripts=True,
                                 merge_strategy='merge',
                                 sort_attribute_values=True)
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath),
                                 keep_order=True)
    return gtf

def main():
    refPath = sys.argv[1]
    GTFsFold = sys.argv[2]
    OutFold = sys.argv[3]

    ref = list(SeqIO.parse(refPath, "fasta"))[0]

    for f in os.listdir(GTFsFold):
        if f[-4:] == ".gtf":
            gene = os.path.basename(f)[:-4]
            gtfPath = os.path.join(GTFsFold, f)
            gtf = openGTF(gtfPath)
            
            out = os.path.join(OutFold, "{}.fasta".format(gene))
            
            for gene in gtf.features_of_type('gene'):
                name = gene.attributes['gene_id'][0]
                s = gene.start
                e = gene.end
                seq = ref[s-1:e]
                SeqIO.write(seq, out, "fasta")

if __name__ == '__main__':
    main()
