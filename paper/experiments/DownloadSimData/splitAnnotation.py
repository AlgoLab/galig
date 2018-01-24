import sys, os

import gffutils

def main():
    gtfPath = sys.argv[1]
    gtfsDir = sys.argv[2]

    print("* Reading GTF")
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
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath), keep_order=True)

    print("* Iterating GTF")
    if not os.path.exists(gtfsDir):
        os.makedirs(gtfsDir)
    for gene in gtf.features_of_type('gene'):
        chrom = gene.seqid
        geneID = gene.id
        geneDir = "{}/{}".format(gtfsDir, chrom)
        if not os.path.exists(geneDir):
            os.makedirs(geneDir)
        out = open("{}/{}.gtf".format(geneDir,geneID), "w")
        out.write(str(gene) + "\n")
        for tr in gtf.children(gene, featuretype='transcript', order_by='start'):
            out.write(str(tr) + "\n")
            for ex in gtf.children(tr, featuretype='exon', order_by='start'):
                out.write(str(ex) + "\n")
        out.close()

if __name__ == '__main__':
    main()
