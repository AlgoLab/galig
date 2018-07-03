import sys, os, argparse

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
        gtf = gffutils.FeatureDB("{}.db".format(gtfPath), keep_order=True)
    return gtf

def contains(gtf, transcript, p1, p2):
    Ss = sorted([exon.start-1 for exon in list(gtf.children(transcript, featuretype='exon', order_by='start'))[1:]])
    Es = sorted([exon.end+1 for exon in list(gtf.children(transcript, featuretype='exon', order_by='start'))[:-1]])
    introns = list(zip(Es, Ss))
    if (p1,p2) in introns:
        return True
    return False

def extractTruth():
    eventsPath = sys.argv[1]
    for line in open(eventsPath).readlines():
        info = line.strip("\n").split("\t")[0]
        if info.startswith("ENSG"):
            geneEv,_,p1,p2,_ = info.split(':')
            gene = geneEv.split(";")[0]
            p1 = int(p1.split("-")[0])+1
            p2 = int(p2.split("-")[1])-1
        print(gene,p1,p2, sep=',')

def main():
    truthPath = sys.argv[1]
    GTFsPath = sys.argv[2]
    outFold = sys.argv[3]

    for line in open(truthPath).readlines()[1:-1]:
        info = line.strip("\n").split("\t")[7]
        if info.startswith("ENSG"):
            geneEv,_,p1,p2,_ = info.split(':')
            geneName = geneEv.split(";")[0]
            p1 = int(p1.split("-")[0])+1
            p2 = int(p2.split("-")[1])-1
            gtfPath = GTFsPath + "/" + geneName + ".gtf"
            if not os.path.isfile(gtfPath):
                continue
            p1,p2 = int(p1),int(p2)
            gtf = openGTF(gtfPath)
            gene = list(gtf.features_of_type('gene'))[0]
            transcripts = list(gtf.children(gene, featuretype='transcript', order_by='start'))

            toAdd = []
            for transcript in transcripts:
                if contains(gtf, transcript, p1, p2):
                    toAdd.append(0)
                else:
                    toAdd.append(1)

            out = None
            if sum(toAdd) > 0:
                outPath = os.path.join(outFold, "{}.gtf".format(geneName))
                if not os.path.exists(os.path.dirname(outPath)):
                    os.makedirs(os.path.dirname(outPath))
                out = open(outPath, 'w')
                out.write("{}\n".format(str(gene)))
            i = 0
            for transcript in transcripts:
                if toAdd[i] == 1:
                    out.write("{}\n".format(str(transcript)))
                    for exon in list(gtf.children(transcript, featuretype='exon', order_by='start')):
                        out.write("{}\n".format(str(exon)))
                i+=1
            if sum(toAdd) > 0:
                out.close()

if __name__ == '__main__':
    main()
