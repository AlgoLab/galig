import sys, os, itertools, operator

from Bio import SeqIO
import gffutils

from BitVector import BitVector

# Computes edit distance between two strings
def editDistance(s1, s2):
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

# L -> (l0,l1), (l1,l2), (l2, l3), ...
def pairwise(L):
    L0, L1 = itertools.tee(L)
    next(L1, None)
    return zip(L0,L1)

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

# Extracts strand, transcripts and introns from gtf
def extractFromGTF(gtf):
    strand = "+"
    introns = set()
    transcripts = {}
    for g in gtf.features_of_type('gene'):
        strand = g.strand
        for tr in gtf.children(g, featuretype='transcript', order_by='start'):
            exons = list(gtf.children(tr, featuretype='exon', order_by='start'))

            transcript = [(ex.start, ex.end) for ex in exons]
            transcripts.update({tr.id:transcript})

            introns_ = set(zip([ex.end+1 for ex in exons[:-1]], [ex.start-1 for ex in exons[1:]]))
            introns = introns | introns_
    return strand, transcripts, introns

# Extracts text and exon positions from "index" file
def extractFromInfoFile(infoPath):
    lines = open(infoPath).readlines()
    text = lines[1].strip("\n")
    exons = [(int(p[0]), int(p[1])) for p in [pos.split(",") for pos in lines[4].strip("\n").split()]]
    return text, exons

# Extracts elements from one line of the custom alignments file
def readLine(line):
    # 0: strand
    # 1: ID
    # 2: errors
    # 3 to -1: mems
    # -1: read
    line = line.strip("\n").strip(" ").split(" ")
    strand = line[0]
    readID = line[1]
    err = int(line[2])
    mems = line[3:-1]
    read = line[-1]
    return readID, err, mems, read

# Removes annotated introns
def filterAnnotated(newIntrons, annIntrons, tresh):
    introns = {}
    allIntrons = set()
    for (p1,p2),w in newIntrons.items():
        if (p1,p2) not in annIntrons:
            introns.update({(p1,p2):w})
            allIntrons.add((p1,p2))
        else:
            allIntrons.add((p1,p2))
            if w >= tresh:
                print("# A {} {}".format(p1,p2))
            else:
                print("# B {} {}".format(p1,p2))
    return introns, allIntrons

# Filters new introns that are not sufficiently covered
def filterLowCovered(newIntrons, tresh):
    introns = {}
    for (p1,p2),w in newIntrons.items():
        if w >= tresh:
            introns.update({(p1,p2):w})
        else:
            print("# W {} {}".format(p1,p2))
    return introns

# Reconciliate a given intron (start/end) position with respect to the input pattern
def reconciliate(pos, ref, patt, isStart):
    off = 1
    MaxOff = 3
    while off<=MaxOff:
        newPos = pos-off
        if isStart:
            newPatt = str(ref[newPos-1:newPos+1].seq)
        else:
            newPatt = str(ref[newPos-2:newPos].seq)
        if newPatt == patt:
            return True,newPos

        newPos = pos+off
        if isStart:
            newPatt = str(ref[newPos-1:newPos+1].seq)
        else:
            newPatt = str(ref[newPos-2:newPos].seq)
        if newPatt == patt:
            return True,newPos
        off+=1
    return False,-1

# Cleans new introns basing on canonical patterns (reconciliation)
def reconciliateIntrons(newIntrons, ref, strand):
    introns = {}
    for (p1,p2),w in newIntrons.items():
        ipS = str(ref[p1-1:p1+1].seq)
        ipE = str(ref[p2-2:p2].seq)
        if strand == '+':
            canIP = {"GT":["AG"], "GC":["AG"]}
            canIPrev = {"AG":["GT", "GC"]}
        else:
            canIP = {"CT":["AC", "GC"]}
            canIPrev = {"AC":["CT"], "GC":["CT"]}
        if ipS in canIP.keys() and ipE in canIP[ipS]:
            key = (p1,p2)
            introns[key] = introns[key]+w if key in introns else w
        elif ipS in canIP.keys():
            newPos = float('inf')
            for rightIP in canIP[ipS]:
                found,pos = reconciliate(p2, ref, rightIP, False)
                if found:
                    if abs(pos-p2) < abs(newPos-p2):
                        newPos = pos
            if newPos == float('inf'):
                newPos = p2
            key = (p1,newPos)
            introns[key] = introns[key]+w if key in introns else w
        elif ipE in canIPrev:
            newPos = float('inf')
            for rightIP in canIPrev[ipE]:
                found,pos = reconciliate(p1, ref, rightIP, True)
                if found:
                    if abs(pos-p1) < abs(newPos-p1):
                        newPos = pos
            if newPos == float('inf'):
                newPos = p1
            key = (newPos,p2)
            introns[key] = introns[key]+w if key in introns else w
        else:
            off = 1
            MaxOff = 3
            while off <= MaxOff:
                newP1, newP2 = p1-1, p2-1
                ipS = str(ref[newP1-1:newP1+1].seq)
                ipE = str(ref[newP2-2:newP2].seq)
                if ipS in canIP.keys() and ipE in canIP[ipS]:
                    key = (newP1,newP2)
                    introns[key] = introns[key]+w if key in introns else w
                    break

                newP1, newP2 = p1+1, p2+1
                ipS = str(ref[newP1-1:newP1+1].seq)
                ipE = str(ref[newP2-2:newP2].seq)
                if ipS in canIP.keys() and ipE in canIP[ipS]:
                    key = (newP1,newP2)
                    introns[key] = introns[key]+w if key in introns else w
                    break
                off+=1
    return introns

def isClose(Ps,p):
    for p_ in Ps:
        if abs(p-p_) <= 100-30: ###readLen - 2*L
            return True
    return False

def isInsideExon(Es,p):
    for (s,e) in Es:
        if s <= p <= e:
            return True
    return False

#######################
def getExonsCloseTo(Es,p):
    exsS = set()
    exsE = set()
    for (s,e) in Es:
        if 0 <= s-p <= 100-30: ###readLen - 2*L
            exsS.add((s,e))
        elif 0 <= p-e <= 100-30: ###readLen - 2*L
            exsE.add((s,e))
    return (list(exsS),list(exsE))

def getExonsContaining(Es,p):
    exs = set()
    for (s,e) in Es:
        if s <= p <= e:
            exs.add((s,e))
    return list(exs)
#######################

def existsIntronEndingAt(introns, p):
    for (s,e) in introns:
        if e == p:
            return True
    return False

def existsIntronStartingAt(introns, p):
    for (s,e) in introns:
        if s == p:
            return True
    return False

def getSuccIntrons(Introns, I):
    Introns.sort()
    i = Introns.index(I)

    Succ = set()
    if i<len(Introns)-1:
        found = False
        while not found and i<len(Introns)-1:
            i+=1
            (s,e) = Introns[i]
            if I[1]<s:
                found = True
                Succ.add((s, e))

        flag = True
        while flag and i<len(Introns)-1:
            i+=1
            (s_,e_) = Introns[i]
            if s<=s_<=e:
                Succ.add((s_, e_))
            else:
                flag = False
    return list(Succ)

def getPrecIntrons(Introns, I):
    Introns.sort(key=operator.itemgetter(1))
    i = Introns.index(I)

    Prec = set()
    if i>0:
        found = False
        while not found and i>0:
            i-=1
            (s,e) = Introns[i]
            if e<I[0]:
                found = True
                Prec.add((s, e))

        flag = True
        while flag and i>0:
            i-=1
            (s_,e_) = Introns[i]
            if s<=e_<=e:
                Prec.add((s_, e_))
            else:
                flag = False
    return list(Prec)

# Extract events from novel introns
def checkNewIntrons(newIntrons, allIntrons, strand, transcripts):
    allIntrons = list(allIntrons)
    allIntrons.sort()
    events = {'ES': {}, 'A3': {}, 'A5': {}, 'IR': {}}
    for (p1,p2),w in newIntrons.items():
        precIntrons = getPrecIntrons(allIntrons, (p1,p2))
        succIntrons = getSuccIntrons(allIntrons, (p1,p2))
        ES = False
        IR = False
        A3 = False
        A5 = False
        for trID,exons in transcripts.items():
            tranSt,tranEnd = exons[0][0], exons[-1][1]
            starts = [ex[0]-1 for ex in exons]
            ends = [ex[1]+1 for ex in exons]

            #Checking exon skippings...
            if p1 in ends and p2 in starts:
                i1 = ends.index(p1)
                i2 = starts.index(p2)
                if i1 != i2-1:
                    key = (p1,p2,w)
                    if key not in events['ES']:
                        events['ES'][key] = []
                    events['ES'][key].append(trID)
                    ES = True

            #Checking intron retentions...
            for (s,e) in exons:
                if s < p1 < p2 < e:
                    if (len(precIntrons) == 0 or s-1 in [i[1] for i in precIntrons]) and (len(succIntrons) == 0 or e+1 in [i[0] for i in succIntrons]):
                        key = (p1,p2,w)
                        if key not in events['IR']:
                            events['IR'][key] = []
                        events['IR'][key].append(trID)
                        IR = True

            #Given (p1,p2), checking alt splice sites:
            if p1 in ends and p1 != tranEnd and p2 not in starts and p2 < tranEnd:
                exonsContaining = getExonsContaining(exons,p2)
                exonsCloseS,_ = getExonsCloseTo(exons,p2)
                #if isInsideExon(exons,p2) or isClose(starts,p2):
                if len(exonsContaining)>0 or len(exonsCloseS)>0:
                    if len(succIntrons) == 0 or any([i[0] in ends for i in succIntrons]) or tranEnd in [e[1] for e in exonsCloseS + exonsContaining]:
                        if strand == '+':
                            t = 'A3'
                            A3 = True
                        else:
                            t = 'A5'
                            A5 = True
                        key = (p1,p2,w)
                        if key not in events[t]:
                            events[t][key] = []
                        events[t][key].append(trID)

            if p1 not in ends and p2 in starts and p2 != tranSt and p1 > tranSt:
                exonsContaining = getExonsContaining(exons,p1)
                _,exonsCloseE = getExonsCloseTo(exons,p1)
                #if isInsideExon(exons,p1) or isClose(ends,p1):
                if len(exonsContaining)>0 or len(exonsCloseE)>0:
                    if len(precIntrons) == 0 or any([i[1] in starts for i in precIntrons]) or tranSt in [e[0] for e in exonsCloseE + exonsContaining]:
                        if strand == '+':
                            t = 'A5'
                            A5 = True
                        else:
                            t = 'A3'
                            A3 = True
                        key = (p1,p2,w)
                        if key not in events[t]:
                            events[t][key] = []
                        events[t][key].append(trID)
        print("# {}-{}: {}, {}, {}, {}".format(p1, p2, ES, IR, A3, A5))
    for t,evs in events.items():
        for (p1,p2,w),trs in evs.items():
            print(t,p1,p2,w,"/".join(trs))

def main():
    #Cmd line arguments
    refPath = sys.argv[1]      #Reference
    gtfPath = sys.argv[2]      #GTF
    memsPath = sys.argv[3]     #MEMs file
    errRate = int(sys.argv[4]) #Error rate
    tresh = int(sys.argv[5])   #New introns treshold
    if len(sys.argv)<7:
        reconciliate = True
    else:
        reconciliate = bool(int(sys.argv[6]))

    Ref = list(SeqIO.parse(refPath, "fasta"))[0]

    gtf = openGTF(gtfPath)
    strand, transcripts, annIntrons = extractFromGTF(gtf)

    text, exons = extractFromInfoFile(gtfPath + ".sg")
    BitV = BitVector(text)

    newIntrons = {}
    ass = {}
    for line in open(memsPath, 'r').readlines():
        readID, err, mems, read = readLine(line)

        if len(mems) > 1:
            for mem1,mem2 in pairwise(mems):
                # Remove ( and ) from mem and cast to int
                mem1 = [int(x) for x in mem1[1:-1].split(",")]
                id1 = BitV.rank(mem1[0] - 1)

                mem2 = [int(x) for x in mem2[1:-1].split(",")]
                id2 = BitV.rank(mem2[0] - 1)

                Poverlap = mem2[1]-mem1[1]-mem1[2]
                if id1 == id2: #MEMs inside the same exon
                    Toverlap = mem2[0]-mem1[0]-mem1[2]
                    # Intron Retention
                    if Poverlap <= 0 and Toverlap > 0:
                        gap = Toverlap+abs(Poverlap)-1
                        if gap > 0:
                            pos1 = exons[id1-1][0] + mem1[0] + mem1[2] - BitV.select(id1) + Poverlap
                            pos2 = pos1 + gap
                            key = (pos1, pos2)
                            newIntrons[key] = newIntrons[key]+1 if key in newIntrons else 1
                else: #MEMs on different exons
                    offset1 = BitV.select(id1 + 1) - (mem1[0] + mem1[2])
                    offset2 = mem2[0] - BitV.select(id2)-1
                    if Poverlap <= 0:
                        Poverlap = abs(Poverlap)
                        # No gap on P: possible Int.Alt.S.S.
                        if offset1 == 0:
                            offset2 += Poverlap
                        else: #anyway, not only if offset2 == 0 !!! maybe this is wrong
                            offset1 += Poverlap
                        pos1 = exons[id1-1][1] - offset1 + 1
                        pos2 = exons[id2-1][0] + offset2 - 1
                        key = (pos1, pos2)
                        newIntrons[key] = newIntrons[key]+1 if key in newIntrons else 1
                    else:
                        #Gap on P
                        if offset1 == 0 and offset2 == 0:
                            #No gap on T -> possible Ext.Alt.S.S.
                            intronStart, intronEnd  = exons[id1-1][1] + 1, exons[id2-1][0] - 1
                            intronString = Ref[intronStart-1:intronEnd-1] #-1 since strings are indexed starting from 0, gtf from 1
                            readGapString = read[mem1[1]+mem1[2]-1:mem2[1]-1]
                            maxErr = round(len(read)*errRate/100)
                            err1 = editDistance(readGapString, intronString[:len(readGapString)])
                            err2 = editDistance(readGapString, intronString[len(intronString)-len(readGapString):])

                            pos1, pos2 = -1, -1
                            if err1 <= err2:
                                # We check the start of the intron
                                if err1 + err <= maxErr:
                                    pos1 = intronStart + Poverlap
                                    pos2 = intronEnd
                            else:
                                #We check the end of the intron
                                if err2 + err <= maxErr:
                                    pos1 = intronStart
                                    pos2 = intronEnd - Poverlap

                            if pos1 != -1 and pos2 != -1:
                                key = (pos1, pos2)
                                newIntrons[key] = newIntrons[key]+1 if key in newIntrons else 1

    # for (p1,p2),w in newIntrons.items():
    #     print("{}-{} : {}".format(p1,p2,w))
    # print("")

    newIntrons, allIntrons = filterAnnotated(newIntrons, annIntrons, tresh)
    if reconciliate:
        newIntrons = reconciliateIntrons(newIntrons, Ref, strand)
    newIntrons, allIntrons_ = filterAnnotated(newIntrons, annIntrons, tresh)
    allIntrons = allIntrons | allIntrons_
    newIntrons = filterLowCovered(newIntrons, tresh)

    # for (p1,p2),w in newIntrons.items():
    #     print("{}-{} : {}".format(p1,p2,w))
    # print("")
    checkNewIntrons(newIntrons, allIntrons, strand, transcripts)

if __name__ == '__main__':
    main()
