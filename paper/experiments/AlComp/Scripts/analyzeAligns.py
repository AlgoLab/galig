import sys, os, time

from Bio import SeqIO
import gffutils

#Read gtf
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
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath), keep_order=True)
    return gtf

#Extract transcript IDs from gtf
def getTrIDs(gtf):
    trIDs = set()

    for g in gtf.features_of_type('gene'):
        for tr in gtf.children(g, featuretype='transcript', order_by='start'):
            trIDs.add(tr.attributes['transcript_id'][0])
    return trIDs

def extractTrID(h):
    return h.split(":")[2]

def main():
    samplePath = sys.argv[1]
    GTFsFold = sys.argv[2]
    SAMsFold = sys.argv[3]

    sample = [extractTrID(x.id) for x in list(SeqIO.parse(samplePath, "fasta"))]

    for f in os.listdir(GTFsFold):
        gtfPath = ""
        gtf = ""
        if f[-4:] == ".gtf":
            gene = os.path.basename(f)[:-4]
            gtfPath = os.path.join(GTFsFold, f)
            gtf = openGTF(gtfPath)
            trIDs = getTrIDs(gtf)
            Truth = 0
            for s in sample:
                if s in trIDs:
                    Truth+=1

            TP = 0
            FP = 0

            #for record in open(os.path.normpath(SAMsFold + "/" + gene.split(".")[0] + ".sam")):
            for record in open(os.path.normpath(SAMsFold + "/" + gene + ".sam")):
                if record[0] != "@":
                    rID, flag, *rest = record.strip("\n").split("\t")
                    rID = extractTrID(rID)
                    if flag in ["0", "16"]:
                        if rID in trIDs:
                            TP += 1
                        else:
                            FP += 1
            FN = Truth - TP
            print(gene, TP, FP, FN, sep=",")

if __name__ == '__main__':
    main()
