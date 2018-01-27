import sys

import gffutils

# Opens gtf (gffutils wrapper)
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

# Extracts transcripts and introns from gtf
def extractFromGTF(gtf):
    strand = "+"
    exons = set()
    for g in gtf.features_of_type('gene'):
        strand = g.strand
        for tr in gtf.children(g, featuretype='transcript', order_by='start'):
            E = list(gtf.children(tr, featuretype='exon', order_by='start'))

            exons_ = set([(ex.start, ex.end) for ex in E])
            exons = exons | exons_
            
    return strand, exons

def getUnannotated(exons, p1, p2, flag):
    p1,p2 = int(p1),int(p2)
    L = []
    if flag: #starts
        L = [e[0] for e in exons]
    else: #ends
        L = [e[1] for e in exons]
    if p1 in L:
        return True,(p2)
    if p2 in L:
        return True,(p1)
    return False,(p1,p2)

def main():
    splOut = sys.argv[1]
    gtfPath = sys.argv[2]

    gtf = openGTF(gtfPath)
    strand, exons = extractFromGTF(gtf)

    for line in open(splOut):
        ev = line.strip("\n").split("\t")
        if ev[0] != "contig":
            if ev[2][:4] == "exon" or ev[2][:4] == "mult":
                p1,p2 = int(ev[5]), int(ev[8])
                print("ES", p1+1, p2-1, "-", "-", sep=",")
            elif ev[2][:5] == "alt_5":
                if ev[1] == '+':
                    f,p = getUnannotated(exons, ev[7], ev[9], False)
                    if f:
                        p1,p2 = p+1,int(ev[4])-1
                        print("A5", p1, p2, "-", "-", sep=",")
                    else:
                        p1,p2 = p[0]+1,int(ev[4])-1
                        print("A5", p1, p2, "-", "-", sep=",")
                        p1,p2 = p[0]+1,int(ev[4])-1
                        print("A5", p1, p2, "-", "-", sep=",")
                else:
                    f,p = getUnannotated(exons, ev[6], ev[8], True)
                    if f:
                        p1,p2 = int(ev[5])+1,p-1
                        print("A5", p1, p2, "-", "-", sep=",")
                    else:
                        p1,p2 = int(ev[5])+1,p[0]-1
                        print("A5", p1, p2, "-", "-", sep=",")
                        p1,p2 = int(ev[5])+1,p[1]-1
                        print("A5", p1, p2, "-", "-", sep=",")
            elif ev[2][:5] == "alt_3":
                if ev[1] == '+':
                    f,p = getUnannotated(exons, ev[6], ev[8], True)
                    if f:
                        p1,p2 = int(ev[5])+1,p-1
                        print("A3", p1, p2, "-", "-", sep=",")
                    else:
                        p1,p2 = int(ev[5])+1,p[0]-1
                        print("A3", p1, p2, "-", "-", sep=",")
                        p1,p2 = int(ev[5])+1,p[1]-1
                        print("A3", p1, p2, "-", "-", sep=",")
                else:
                    f,p = getUnannotated(exons, ev[7], ev[9], False)
                    if f:
                        p1,p2 = p+1,int(ev[4])-1
                        print("A3", p1, p2, "-", "-", sep=",")
                    else:
                        p1,p2 = p[0]+1,int(ev[4])-1
                        print("A3", p1, p2, "-", "-", sep=",")
                        p1,p2 = p[1]+1,int(ev[4])-1
                        print("A3", p1, p2, "-", "-", sep=",")
            elif ev[2][:4] == "intr":
                p1,p2 = int(ev[6]), int(ev[7])
                print("IR", p1+1, p2-1, "-", "-", sep=",")

if __name__ == '__main__':
    main()
