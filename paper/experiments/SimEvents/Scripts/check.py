import sys, os

import gffutils

def getUnannotated(exons, strand, evT, pos1, pos2):
    AnnS = [ex[0] for ex in exons]
    AnnE = [ex[1] for ex in exons]
    if strand == '+':
        if evT == 'A3':
            if pos1 in AnnS and pos2 not in AnnS:
                return pos2
            elif pos1 not in AnnS and pos2 in AnnS:
                return pos1
            else:
                return False
        else:
            if pos1 in AnnE and pos2 not in AnnE:
                return pos2
            elif pos1 not in AnnE and pos2 in AnnE:
                return pos1
            else:
                return False
    else:
        if evT == 'A3':
            if pos1 in AnnE and pos2 not in AnnE:
                return pos2
            elif pos1 not in AnnE and pos2 in AnnE:
                return pos1
            else:
                return False
        else:
            if pos1 in AnnS and pos2 not in AnnS:
                return pos2
            elif pos1 not in AnnS and pos2 in AnnS:
                return pos1
            else:
                return False

#Build truth combining infoPath and truthPath
def buildTruth(eventID, infoPath, truthPath, exons, strand):
    IDs = []
    for line in open(infoPath).readlines():
        uID, *eIDs = line.strip("\n").split(" ")
        if uID == eventID:
            realeIDs = []
            for e in eIDs:
                realeIDs.append(e)
            IDs = realeIDs
            break

    truthTypes = []
    truth = {'ES':[], 'A3':[], 'A5':[], 'IR':[]}
    totalTruth = {'ES':[], 'A3':[], 'A5':[], 'IR':[]}
    for line in open(truthPath).readlines():
        if line[0] == "#":
            continue
        info = line.strip("\n").split(" ")
        evT = info[1]
        if evT[0] == 'A':
            pos1, pos2, pos3, uID  = info[2], info[3], info[4], info[5]
            pos1, pos2, pos3 = int(pos1), int(pos2), int(pos3)
            if uID+"e" in IDs:
                if abs(pos3-pos1)<abs(pos3-pos2):
                    truth[evT[0:2]].append((min(pos1, pos3), max(pos1, pos3)))
                    truthTypes.append(evT+"e")
                    totalTruth[evT[0:2]].append((min(pos2, pos3), max(pos2, pos3)))
                else:
                    truth[evT[0:2]].append((min(pos2, pos3), max(pos2, pos3)))
                    truthTypes.append(evT+"e")
                    totalTruth[evT[0:2]].append((min(pos1, pos3), max(pos1, pos3)))
            elif uID+"i" in IDs:
                if abs(pos3-pos1)<abs(pos3-pos2):
                    truth[evT[0:2]].append((min(pos2, pos3), max(pos2, pos3)))
                    truthTypes.append(evT+"i")
                    totalTruth[evT[0:2]].append((min(pos1, pos3), max(pos1, pos3)))
                else:
                    truth[evT[0:2]].append((min(pos1, pos3), max(pos1, pos3)))
                    truthTypes.append(evT+"i")
                    totalTruth[evT[0:2]].append((min(pos2, pos3), max(pos2, pos3)))
            else:
                totalTruth[evT[0:2]].append((min(pos1, pos3), max(pos1, pos3)))
                totalTruth[evT[0:2]].append((min(pos2, pos3), max(pos2, pos3)))
        else:
            pos1, pos2, uID  = info[2], info[3], info[4]
            pos1,pos2 = int(pos1),int(pos2)
            if uID in IDs:
                truth[evT].append((min(pos1, pos2), max(pos1, pos2)))
                truthTypes.append(evT)
            else:
                totalTruth[evT].append((min(pos1, pos2), max(pos1, pos2)))

    return truth, truthTypes, totalTruth

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

#Extract strand and exons from gtf
def extractFromGTF(gtf):
    chrom = ""
    strand = ""
    exons = set()
    introns = set()

    for g in gtf.features_of_type('gene'):
        chrom = g.seqid
        strand = g.strand
        for tr in gtf.children(g, featuretype='transcript', order_by='start'):
            exs = list(gtf.children(tr, featuretype='exon', order_by='start'))
            exons = exons | set([(ex.start, ex.end) for ex in exs])

            ExSs = [ex.start-1 for ex in exs]
            ExEs = [ex.end+1 for ex in exs]
            I = set(list(zip(ExEs[:-1], ExSs[1:])))
            introns = introns | I

    return chrom, strand, exons, introns

def isExternal(t, pos1, pos2, strand, exons):
    if t == '3':
        if strand == '+':
            if pos1 in [e[0] for e in exons]:
                if pos1 < pos2:
                    return False
                else:
                    return True
            elif pos2 in [e[0] for e in exons]:
                if pos1 < pos2:
                    return True
                else:
                    return False
        elif strand == '-':
            if pos1 in [e[1] for e in exons]:
                if pos1 < pos2:
                    return True
                else:
                    return False
            elif pos2 in [e[1] for e in exons]:
                if pos1 < pos2:
                    return False
                else:
                    return True
    elif t == '5':
        if strand == '+':
            if pos1 in [e[1] for e in exons]:
                if pos1 < pos2:
                    return True
                else:
                    return False
            elif pos2 in [e[1] for e in exons]:
                if pos1 < pos2:
                    return False
                else:
                    return True
        elif strand == '-':
            if pos1 in [e[0] for e in exons]:
                if pos1 < pos2:
                    return False
                else:
                    return True
            elif pos2 in [e[0] for e in exons]:
                if pos1 < pos2:
                    return True
                else:
                    return False

def main():
    eventsPath = sys.argv[1]
    infoPath = sys.argv[2]
    truthPath = sys.argv[3]
    gtfPath = sys.argv[4]

    gtf = openGTF(gtfPath)
    chrom, strand, exons, introns = extractFromGTF(gtf)

    eventID = os.path.basename(eventsPath)[:-7]
    gene = os.path.basename(os.path.dirname(eventsPath))

    truth, truthTypes, totalTruth = buildTruth(eventID, infoPath, truthPath, exons, strand)

    TP = 0
    FP = 0
    FPevs = {'ES':[], 'A3':[], 'A5':[], 'IR':[]}
    checked = {'ES':[], 'A3':[], 'A5':[], 'IR':[]}

    for line in open(eventsPath).readlines():
        if line[0:4] == "Type":
            continue
        evT, p1, p2, _, _ = line.strip("\n").split(",")
        p1,p2 = int(p1)-1,int(p2)+1
        pos1,pos2 = min(p1,p2), max(p1,p2)
        if (pos1,pos2) in truth[evT]:
            if (pos1,pos2) not in checked[evT]:
                TP+=1
                checked[evT].append((pos1,pos2))
        else:
            if (pos1+1, pos2-1) not in introns:
                if (pos1,pos2) not in totalTruth[evT]:
                    FPevs[evT].append((pos1,pos2))
                    FP+=1

    FN = len(truthTypes)-TP
    t = "-"
    if len(truthTypes) == 1:
        t = truthTypes[0]

    fpES = "."
    if len(FPevs['ES']) != 0:
        fpES = "/".join([":".join([str(x1+1),str(x2-1)]) for (x1,x2) in FPevs['ES']])
    fpA3 = "."
    if len(FPevs['A3']) != 0:
        fpA3 = "/".join([":".join([str(x1+1),str(x2-1)]) for (x1,x2) in FPevs['A3']])
    fpA5 = "."
    if len(FPevs['A5']) != 0:
        fpA5 = "/".join([":".join([str(x1+1),str(x2-1)]) for (x1,x2) in FPevs['A5']])
    fpIR = "."
    if len(FPevs['IR']) != 0:
        fpIR = "/".join([":".join([str(x1+1),str(x2-1)]) for (x1,x2) in FPevs['IR']])
    print(chrom,gene,eventID,len(truthTypes),t,TP,FN,FP,fpES,fpA3,fpA5,fpIR,sep=',')

if __name__ == '__main__':
    main()
