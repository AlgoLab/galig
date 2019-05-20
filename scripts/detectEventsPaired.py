#!/usr/bin/python3

import sys, os, itertools, operator
import argparse
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

# TODO: extract directly from FASTA + GTF
# Extracts text and exon positions from "index" file
def extractFromInfoFile(infoPath):
    lines = open(infoPath).readlines()
    text = lines[1].strip("\n")
    exons = [(int(p[0]), int(p[1])) for p in [pos.split(",") for pos in lines[4].strip("\n").split()]]
    return text, exons

# Extracts elements from one line of the spliced graph-alignments file
def readLine(line):
    mapped = True
    strand = ""
    readID = ""
    err = 0
    mems = []
    read = ""

    line = line.strip("\n").strip(" ").split(" ")
    # 0: MAPPED/UNMAPPED
    mapped = True if line[0] == "MAPPED" else False
    placeholder = True if line[0] == "PLACEHOLDER" else False

    if mapped:
        # 1: strand
        # 2: ID
        # 3: errors
        # 4 to -1: mems
        # -1: read
        strand = line[1]
        readID = line[2]
        err = int(line[3])
        mems = line[4:-1]
        read = line[-1]
    elif not mapped and not placeholder:
        # 1: ID
        # 2: read
        readID = line[1]
        read = line[-1]
    else: # not mapped and placeholder
        # 1: ID
        readID = line[1]
    return placeholder, mapped, strand, readID, err, mems, read

# Removes annotated introns
def filterAnnotated(newIntrons, annIntrons):
    introns = {}
    annFoundIntrons = {}
    for (p1,p2),w in newIntrons.items():
        if (p1,p2) not in annIntrons:
            introns.update({(p1,p2):w})
        else:
            annFoundIntrons.update({(p1,p2):w})
    return introns, annFoundIntrons

# Filters new introns that are not sufficiently covered
def filterLowCovered(introns, tresh):
    filtIntrons = {}
    for (p1,p2),w in introns.items():
        if w >= tresh:
            filtIntrons.update({(p1,p2):w})
        #else:
        #    print("# W {} {}".format(p1,p2))
    return filtIntrons

# Reconciliate a given intron (start/end) position with respect to the input pattern
def reconciliate(pos, ref, patt, isStart):
    off = 1
    MaxOff = 3 # TODO: this one could be a parameter
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

# Cleans introns basing on canonical patterns (reconciliation)
def reconciliateIntrons(introns, ref, strand):
    recIntrons = {}
    for (p1,p2),w in introns.items():
        intronPref = str(ref[p1-1:p1+1].seq)
        intronSuff = str(ref[p2-2:p2].seq)
        if strand == '+':
            canIP = {"GT":["AG"], "GC":["AG"]}
            canIPrev = {"AG":["GT", "GC"]}
        else:
            canIP = {"CT":["AC", "GC"]}
            canIPrev = {"AC":["CT"], "GC":["CT"]}
        if intronPref in canIP.keys() and intronSuff in canIP[intronPref]:
            key = (p1,p2)
            recIntrons[key] = recIntrons[key]+w if key in recIntrons else w
        elif intronPref in canIP.keys():
            newPos = float('inf')
            for acceptedIntronicPattern in canIP[intronPref]:
                found,pos = reconciliate(p2, ref, acceptedIntronicPattern, False)
                if found:
                    if abs(pos-p2) < abs(newPos-p2):
                        newPos = pos
            if newPos == float('inf'):
                newPos = p2
            key = (p1,newPos)
            recIntrons[key] = recIntrons[key]+w if key in recIntrons else w
        elif intronSuff in canIPrev:
            newPos = float('inf')
            for acceptedIntronicPattern in canIPrev[intronSuff]:
                found,pos = reconciliate(p1, ref, acceptedIntronicPattern, True)
                if found:
                    if abs(pos-p1) < abs(newPos-p1):
                        newPos = pos
            if newPos == float('inf'):
                newPos = p1
            key = (newPos,p2)
            recIntrons[key] = recIntrons[key]+w if key in recIntrons else w
        else:
            off = 1
            MaxOff = 3
            while off <= MaxOff:
                newP1, newP2 = p1-1, p2-1
                intronPref = str(ref[newP1-1:newP1+1].seq)
                intronSuff = str(ref[newP2-2:newP2].seq)
                if intronPref in canIP.keys() and intronSuff in canIP[intronPref]:
                    key = (newP1,newP2)
                    recIntrons[key] = recIntrons[key]+w if key in recIntrons else w
                    break

                newP1, newP2 = p1+1, p2+1
                intronPref = str(ref[newP1-1:newP1+1].seq)
                intronSuff = str(ref[newP2-2:newP2].seq)
                if intronPref in canIP.keys() and intronSuff in canIP[intronPref]:
                    key = (newP1,newP2)
                    recIntrons[key] = recIntrons[key]+w if key in recIntrons else w
                    break
                off+=1
    return recIntrons

# Returns all the exons starting and ending close to a given position
def getExonsCloseTo(Es,p):
    exsS = set()
    exsE = set()
    for (s,e) in Es:
        if 0 <= s-p <= 100-30: # TODO: set 100-30 to readLen - 2*L
            exsS.add((s,e))
        elif 0 <= p-e <= 100-30: # TODO: set 100-30 to readLen - 2*L
            exsE.add((s,e))
    return (list(exsS),list(exsE))

# Returns all the exons containing a position
def getExonsContaining(Es,p):
    exs = set()
    for (s,e) in Es:
        if s <= p <= e:
            exs.add((s,e))
    return list(exs)

# Returns True if there exists an intron ending at the given position, False otherwise
def existsIntronEndingAt(introns, p):
    for (s,e) in introns:
        if e == p:
            return True
    return False

# Returns True if there exists an intron starting at the given position, False otherwise
def existsIntronStartingAt(introns, p):
    for (s,e) in introns:
        if s == p:
            return True
    return False

# Returns all the (possibly overlapping) introns that follow the given intron
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

# Returns all the (possibly overlapping) introns that precede the given intron
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

# Extracts events from introns
def checkNewIntrons(newIntrons, allIntrons, strand, transcripts):
    # TODO: use idmp to make more accurate guesses

    allIntrons = list(allIntrons)
    allIntrons.sort()
    events = {'ES': {}, 'A3': {}, 'A5': {}, 'IR': {}}
    for (p1,p2),w in newIntrons.items():
        # Getting introns preceding and following the considered intron
        precIntrons = getPrecIntrons(allIntrons, (p1,p2))
        succIntrons = getSuccIntrons(allIntrons, (p1,p2))

        # For each transcript...
        for trID,exons in transcripts.items():
            tranSt,tranEnd = exons[0][0], exons[-1][1]
            intronsEnds = [ex[0]-1 for ex in exons]
            intronsStarts = [ex[1]+1 for ex in exons]

            # Checking exon skippings
            if p1 in intronsStarts and p2 in intronsEnds:
                i1 = intronsStarts.index(p1)
                i2 = intronsEnds.index(p2)
                if i1 != i2-1:
                    key = (p1,p2,w)
                    if key not in events['ES']:
                        events['ES'][key] = []
                    events['ES'][key].append(trID)

            # Checking intron retentions
            for (s,e) in exons:
                if s < p1 < p2 < e:
                    if (len(precIntrons) == 0 or s-1 in [i[1] for i in precIntrons]) and (len(succIntrons) == 0 or e+1 in [i[0] for i in succIntrons]):
                        key = (p1,p2,w)
                        if key not in events['IR']:
                            events['IR'][key] = []
                        events['IR'][key].append(trID)

            # Checking alternative splice sites
            if p1 in intronsStarts and p1 != tranEnd and p2 not in intronsEnds and p2 < tranEnd:
                exonsContaining = getExonsContaining(exons,p2)
                exonsCloseS,_ = getExonsCloseTo(exons,p2)
                if len(exonsContaining)>0 or len(exonsCloseS)>0:
                    if len(succIntrons) == 0 or any([i[0] in intronsStarts for i in succIntrons]) or tranEnd in [e[1] for e in exonsCloseS + exonsContaining]:
                        if strand == '+':
                            t = 'A3'
                        else:
                            t = 'A5'
                        key = (p1,p2,w)
                        if key not in events[t]:
                            events[t][key] = []
                        events[t][key].append(trID)

            if p1 not in intronsStarts and p2 in intronsEnds and p2 != tranSt and p1 > tranSt:
                exonsContaining = getExonsContaining(exons,p1)
                _,exonsCloseE = getExonsCloseTo(exons,p1)
                if len(exonsContaining)>0 or len(exonsCloseE)>0:
                    if len(precIntrons) == 0 or any([i[1] in intronsEnds for i in precIntrons]) or tranSt in [e[0] for e in exonsCloseE + exonsContaining]:
                        if strand == '+':
                            t = 'A5'
                        else:
                            t = 'A3'
                        key = (p1,p2,w)
                        if key not in events[t]:
                            events[t][key] = []
                        events[t][key].append(trID)

                #check new splicing event

    return events

# Printing events (TODO: they can be printed when found, but maybe the dict could be useful for some analysis)
def printEvents(events, outPath):
    out = open(outPath, 'w')
    out.write("Type,Start,End,Support,Transcripts\n")
    for t,evs in events.items():
        for (p1,p2,w),trs in evs.items():
            out.write("{},{},{},{},{}\n".format(t,p1,p2,w,"/".join(trs)))

def extractIntrons(memsPath, Ref, exons, BitV, errRate, onlyPrimary):
    introns = {}
    lastID = ""
    for line in open(memsPath, 'r').readlines():
        placeholder, mapped, alStrand, readID, err, mems, read = readLine(line)
        if onlyPrimary:
            if readID == lastID:
                continue
            lastID = readID

        if mapped:
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
                                introns[key] = introns[key]+1 if key in introns else 1
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
                            introns[key] = introns[key]+1 if key in introns else 1
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
                                    introns[key] = introns[key]+1 if key in introns else 1
    return introns

# Merge two introns-dict into one
# NOTE: an introns-dict has the following structure: {(intronStart, intronEnd): number_of_reads}
def mergeIntrons(introns1, introns2):
    introns = {}
    for (p1,p2),w in introns1.items():
            introns[(p1,p2)] = w
    for (p1,p2),w in introns2.items():
            if (p1,p2) not in introns.keys():
                introns[(p1,p2)] = w
            else:
                introns[(p1,p2)] += w
    return introns

def main(memsPath1, memsPath2, refPath, gtfPath, errRate, tresh, outPath, allevents):
    #TODO: add this as cmd line parameter
    onlyPrimary = False

    # Reading reference genome
    Ref = list(SeqIO.parse(refPath, "fasta"))[0]

    # Reading annotation
    gtf = openGTF(gtfPath)
    strand, transcripts, annIntrons = extractFromGTF(gtf)

    # Extracting text and exons from "index" file
    text, exons = extractFromInfoFile(gtfPath + ".sg")
    BitV = BitVector(text)

    # Extracting introns from spliced graph-alignments
    introns1 = extractIntrons(memsPath1, Ref, exons, BitV, errRate, onlyPrimary)
    introns2 = extractIntrons(memsPath2, Ref, exons, BitV, errRate, onlyPrimary)
    introns = mergeIntrons(introns1, introns2)

    # Cleaning introns
    if not(allevents):
        newIntrons, annFoundIntrons = filterAnnotated(introns, annIntrons)
        newIntrons = reconciliateIntrons(newIntrons, Ref, strand)
        newIntrons, annFoundIntrons_ = filterAnnotated(newIntrons, annIntrons)
        newIntrons = filterLowCovered(newIntrons, tresh)
        annFoundIntrons = filterLowCovered(annFoundIntrons, tresh)
        annFoundIntrons_ = filterLowCovered(annFoundIntrons_, tresh)
        allIntronsKey = set(newIntrons.keys()) | set(annFoundIntrons.keys())
    else:
        newIntrons = reconciliateIntrons(introns, Ref, strand)
        newIntrons = filterLowCovered(newIntrons, tresh)
        allIntronsKey = newIntrons.keys()

    # Extracting events from introns
    events = checkNewIntrons(newIntrons, allIntronsKey, strand, transcripts)
    printEvents(events, outPath)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Detects alternative splicing events from splice-aware alignments to a splicing graph")
    parser.add_argument('-g', '--genome', required=True, help='FASTA input file containing the reference')
    parser.add_argument('-a', '--annotation', required=True, help='GTF input file containing the gene annotation')
    parser.add_argument('-1', '--mems1', required=True, help='input file containing the alignments to the splicing graph')
    parser.add_argument('-2', '--mems2', required=True, help='input file containing the alignments to the splicing graph')
    parser.add_argument('-o', '--output', required=True, help='SAM output file')
    #parser.add_argument('-i', '--idmp', required=True, help='inner distance between mate-pairs')
    parser.add_argument('-e', '--erate', required=False, default=3, type=int, help='error rate of alignments (from 0 to 100, default: 3)')
    parser.add_argument('-w', '--support', required=False, default=3, type=int, help='minimum number of reads needed to confirm an event (default: 3)')
    parser.add_argument('--allevents', required=False, action='store_true', help='output all events, not only the novel ones')

    args = parser.parse_args()

    memsPath1 = args.mems1
    memsPath2 = args.mems2
    refPath = args.genome
    gtfPath = args.annotation
    errRate = args.erate
    tresh = args.support
    outPath = args.output
    allevents = args.allevents
    #idmp = args.idmp

    main(memsPath1, memsPath2, refPath, gtfPath, errRate, tresh, outPath, allevents)
