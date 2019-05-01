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

# Extracts text and exon positions from "index" file
def extractFromInfoFile(infoPath):
    lines = open(infoPath).readlines()

    line1 = lines[0].split(" ")
    ref, refLen = line1[0], int(line1[1])

    text = lines[1].strip("\n")

    exons = [(int(p[0]), int(p[1])) for p in [pos.split(",") for pos in lines[4].strip("\n").split()]]

    return ref, refLen, text, exons

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


def extractMEMs(mems):
    MEMs = []
    for mem in mems:
        # Remove ( and ) from mem and cast to int
        mem = [int(x) for x in mem[1:-1].split(",")]
        MEMs.append(mem)
    return MEMs

def getStart(mem, bv, exPos):
    i = bv.rank(mem[0])
    return exPos[i-1][0] + (mem[0] - bv.select(i)) - 1

def getFlagPaired(mapped1, strand1, readID1, mapped2, strand2, readID2, read1=True):
    f = 1 #since it's a paired-end read

    # SEE: http://seqanswers.com/forums/showthread.php?t=17314

    if mapped1 and mapped2:  # Both reads mapped
        if read1: # setting flag for read coming from file 1
            if strand1 == "-":
                if strand2 == "-": # - - read1
                    f = 113
                else:  # - + read1
                    f = 81
            else:
                if strand2 == "-": # + - read1
                    f = 97
                else: # + + read1
                    f = 65

        else:   # setting flag for read coming from file 2
            if strand2 == "-":
                if strand1 == "-": # - - read2
                    f = 177
                else:  # + - read2
                    f = 145
            else:
                if strand1 == "-": # - + read2
                    f = 161
                else: # + + read2
                    f = 129

    elif mapped1 and not mapped2:  # only read1 mapped
        if read1:
            f = 73
        else:
            f = 133
    elif not mapped1 and mapped2:  # only read2 mapped
        if read1:
            f = 69
        else:
            f = 137
    else:   # not mapped1 and not mapped2
        if read1:
            f = 77
        else:   # read2
            f = 141

    return f

# Reconciliate a given intron (s:e)
def reconciliate(s, e, RefSeq):
    off = 0
    MaxOff = 3 # TODO: this one could be a parameter
    while off<=MaxOff:
        Iseq = RefSeq[s-off-1:e-off].seq
        if Iseq[:2] in ["GT", "GC"] and Iseq[-2:] == "AG":
            return -off
        Iseq = RefSeq[s+off-1:e+off].seq
        if Iseq[:2] in ["GT", "GC"] and Iseq[-2:] == "AG":
            return off
        off+=1
    return ""


def getCIGAR(mems, RefSeq, bv, exPos, read, errRate, err):
    m = len(read)
    cigarList = []
    CIGAR = ""
    i = 0
    while i<len(mems): # if mems == [] (unmapped read), len(mems) == 0, therefore the loop will be skipped
        N = 0
        M = 0
        if i == 0:
            if mems[i][1] != 1:
                initialClips = mems[i][1]-1
                cigarList.append([initialClips, 'S'])
            cigarList.append([mems[i][2], 'M'])
        else:
            ##################################################################
            id1 = bv.rank(mems[i-1][0] - 1)
            id2 = bv.rank(mems[i][0] - 1)
            if id1 == id2:
                errors_P = mems[i][1] - mems[i-1][1] - mems[i-1][2]
                errors_T = mems[i][0] - mems[i-1][0] - mems[i-1][2]
                #------------------------------------------------
                if errors_P == 0:
                    if errors_T < 0:
                        #0 on P, - on T
                        cigarList.append([abs(errors_T), 'I'])
                        cigarList.append([mems[i][2]-abs(errors_T), 'M'])
                    elif errors_T > 0:
                        #0 on P, + on T
                        cigarList.append([errors_T, 'D'])
                        cigarList.append([mems[i][2], 'M'])
                #------------------------------------------------
                elif errors_P < 0:
                    if errors_T == 0:
                        #- on P, 0 on T
                        cigarList.append([abs(errors_P), 'D'])
                        cigarList.append([mems[i][2]-abs(errors_P), 'M'])
                    elif errors_T < 0:
                        #- on P, - on T
                        errs = abs(errors_P) - abs(errors_T)
                        if errs > 0:
                            #More on P
                            cigarList.append([errs, 'D'])
                        elif errs < 0:
                            #More on T
                            cigarList.append([abs(errs), 'I'])
                        cigarList.append([mems[i][2]-max(abs(errors_P), abs(errors_T)), 'M'])
                    elif errors_T > 0:
                        #- on P, + on T
                        I = (exPos[id1-1][0]+mems[i-1][0]+mems[i-1][2]-bv.select(id1)-1 - abs(errors_P), exPos[id1-1][0]+mems[i-1][0]+mems[i-1][2]-bv.select(id1)-2 + errors_T)
                        Iseq = RefSeq[I[0]-1:I[1]].seq
                        if Iseq[:2] in ["GT", "GC"] and Iseq[-2:] == "AG":
                            cigarList[-1][0] = cigarList[-1][0] - abs(errors_P)
                            N = abs(errors_P) + errors_T
                            if N>0:
                                cigarList.append([N, 'N'])
                            else:
                                if cigarList[-1][1] is 'M':
                                    M = cigarList[-1][0]
                                    cigarList.pop(-1)
                            cigarList.append([M + mems[i][2], 'M'])
                        else:
                            N = abs(errors_P) + errors_T
                            if N>0:
                                cigarList.append([N, 'N'])
                            else:
                                if cigarList[-1][1] is 'M':
                                    M = cigarList[-1][0]
                                    cigarList.pop(-1)
                            cigarList.append([M + mems[i][2]-abs(errors_P), 'M'])
                #------------------------------------------------
                elif errors_P > 0:
                    if errors_T == 0:
                        #+ on P, 0 on T
                        cigarList.append([errors_P, 'I'])
                        cigarList.append([mems[i][2], 'M'])
                    elif errors_T < 0:
                        #+ on P, - on T
                        cigarList.append([abs(errors_P) + abs(errors_T), 'I'])
                        cigarList.append([mems[i][2]-abs(errors_T), 'M'])
                    elif errors_T > 0:
                        #+ on P, + on T
                        if errors_P == errors_T:
                            #Same: all SUBs
                            cigarList[-1][0] = cigarList[-1][0] + errors_P + mems[i][2]
                        elif errors_P > errors_T:
                            #More on P: some INSs, some SUBs
                            cigarList.append([errors_P - errors_T, 'I'])
                            cigarList.append([errors_T + mems[i][2], 'M'])
                        else:
                            #More on T: some DELs, some SUBs
                            cigarList.append([errors_T - errors_P, 'D'])
                            cigarList.append([errors_P + mems[i][2], 'M'])
            ##################################################################
            else:
                errors_P = mems[i][1] - mems[i-1][1] - mems[i-1][2]
                errors_T1 = bv.select(bv.rank(mems[i-1][0]) + 1) - mems[i-1][0] - mems[i-1][2]
                errors_T2 = mems[i][0] - bv.select(bv.rank(mems[i][0])) - 1
                intron = exPos[id2-1][0] - exPos[id1-1][1] - 1
                #------------------------------------------------
                if errors_P == 0:
                    N = intron + errors_T1 + errors_T2
                    if N>0:
                        cigarList.append([N, 'N'])
                    else:
                        if cigarList[-1][1] is 'M':
                            M = cigarList[-1][0]
                            cigarList.pop(-1)
                    cigarList.append([M + mems[i][2], 'M'])
                #------------------------------------------------
                elif errors_P < 0:
                    if errors_T1 == 0 and errors_T2 != 0:
                        N = intron + errors_T1 + errors_T2 + abs(errors_P)
                        if N>0:
                            cigarList.append([N, 'N'])
                        else:
                            if cigarList[-1][1] is 'M':
                                M = cigarList[-1][0]
                                cigarList.pop(-1)
                        cigarList.append([M + mems[i][2]-abs(errors_P), 'M'])
                    elif errors_T1 != 0 and errors_T2 == 0:
                        N = 0
                        cigarList[-1][0] = cigarList[-1][0] - abs(errors_P)
                        if cigarList[-1][0]<=0:
                            M = cigarList[-1][0]
                            cigarList.pop(-1)
                            if cigarList[-1][1] is 'N':
                                cigarList[-2][0] = cigarList[-2][0] + M
                                M = 0
                                N = cigarList[-1][0]
                                cigarList.pop(-1)
                        N += intron + errors_T1 + errors_T2 + abs(errors_P)
                        if N>0:
                            cigarList.append([N, 'N'])
                        else:
                            if cigarList[-1][1] is 'M':
                                M = cigarList[-1][0]
                                cigarList.pop(-1)
                        cigarList.append([M + mems[i][2], 'M'])
                    else:
                        I = (exPos[id1-1][0]+mems[i-1][0]+mems[i-1][2]-bv.select(id1)-1 - abs(errors_P), exPos[id1-1][0]+mems[i-1][0]+mems[i-1][2]-bv.select(id1) + intron+errors_T1+errors_T2-2)
                        Iseq = RefSeq[I[0]-1:I[1]].seq
                        if Iseq[:2] in ["GT", "GC"] and Iseq[-2:] == "AG":
                            cigarList[-1][0] = cigarList[-1][0] - abs(errors_P)
                            N = intron + errors_T1 + errors_T2 + abs(errors_P)
                            if N>0:
                                cigarList.append([N, 'N'])
                            else:
                                if cigarList[-1][1] is 'M':
                                    M = cigarList[-1][0]
                                    cigarList.pop(-1)
                            cigarList.append([M + mems[i][2], 'M'])
                        else:
                            N = intron + errors_T1 + errors_T2 + abs(errors_P)
                            if N>0:
                                cigarList.append([N, 'N'])
                            else:
                                if cigarList[-1][1] is 'M':
                                    M = cigarList[-1][0]
                                    cigarList.pop(-1)
                            cigarList.append([M + mems[i][2]-abs(errors_P), 'M'])
                #------------------------------------------------
                else:
                    I = (exPos[id1-1][0]+mems[i-1][0]+mems[i-1][2]-bv.select(id1)-1, exPos[id1-1][0]+mems[i-1][0]+mems[i-1][2]-bv.select(id1)-2 + intron+errors_T1+errors_T2)
                    intronString = RefSeq[I[0]-1:I[1]]
                    readGapString = read[mems[i-1][1]+mems[i-1][2]-1:mems[i][1]-1]
                    maxErr = round(len(read)*errRate/100)
                    err1 = editDistance(readGapString, intronString[:len(readGapString)])
                    err2 = editDistance(readGapString, intronString[len(intronString)-len(readGapString):])
                    if err1 <= err2 and err1 + err <= maxErr:
                        cigarList[-1][0] = cigarList[-1][0] + len(readGapString)
                        N = intron + errors_T1 + errors_T2 - len(readGapString)
                        if N <= 0:
                            cigarList[-1][0] = cigarList[-1][0] + mems[i][2]
                        else:
                            if N>0:
                                cigarList.append([N, 'N'])
                            else:
                                if cigarList[-1][1] is 'M':
                                    M = cigarList[-1][0]
                                    cigarList.pop(-1)
                            cigarList.append([M + mems[i][2], 'M'])
                    elif err2 <= err1 and err2 + err <= maxErr:
                        N = intron + errors_T1 + errors_T2 - len(readGapString)
                        if N <= 0:
                            cigarList[-1][0] = cigarList[-1][0] + len(readGapString) + mems[i][2]
                        else:
                            if N>0:
                                cigarList.append([N, 'N'])
                            else:
                                if cigarList[-1][1] is 'M':
                                    M = cigarList[-1][0]
                                    cigarList.pop(-1)
                            cigarList.append([M + len(readGapString)+mems[i][2], 'M'])
                    else:
                        cigarList.append([abs(errors_P), 'I'])
                        N = intron + errors_T1 + errors_T2
                        if N>0:
                            cigarList.append([N, 'N'])
                        else:
                            if cigarList[-1][1] is 'M':
                                M = cigarList[-1][0]
                                cigarList.pop(-1)
                        cigarList.append([M + mems[i][2], 'M'])
        i+=1
    if mems != []: # same as before
        finalDels = m - mems[-1][1] - mems[-1][2] + 1
        if finalDels != 0:
            cigarList.append([finalDels, 'S'])

        CIGAR = ""
        for (n,l) in cigarList:
            CIGAR += str(n) + l
    else:   # mems == [] (Unmapped slignment)
        CIGAR = str(len(read)) + "X" # Mismatches only
    return CIGAR

def getEnd(mem, bv, exPos):
    exonN = bv.rank(mem[0])

    # get starting position (on reference)
    exonStartingPos = exPos[exonN-1][0]

    # get offset from bitvector
    exonStartingPos_onT = bv.select(exonN)
    offset = mem[0] - exonStartingPos_onT + 1;

    # find the end
    # NOTE: offset is the same on reference and bitvector
    end = exonStartingPos + offset + mem[2]

    return end

def getTlen(start1, mems2, bv, exPos):
    lastMem2 = mems2[-1]

    # get ending position (in the reference)
    end2 = getEnd(lastMem2, bv, exPos)

    return end2 - start1

def getIdmp(start2, mems1, bv, exPos):
    lastMem1 = mems1[-1]

    # get ending position (in the reference)
    end1 = getEnd(lastMem1, bv, exPos)

    return start2 - end1

# Get transcript-based idmp
def getTranscriptIdmp(transcripts, mems1, mems2, bv, exPos):
    tIdmp = 0

    # get last mem from read1 and first mem from read2
    m1 = mems1[-1]
    m2 = mems2[0]

    # find exon for m1 and m2
    id1 = bv.rank(m1[0]) - 1
    id2 = bv.rank(m2[0]) - 1

    # find positions on bitvector
    start1 = m1[0]
    start2 = m2[0]
    len1 = m1[2]
    len2 = m2[2]

    if id1 == id2: # m1 and m2 are on same exon
        distance = start2 - (start1 + len1)
        tIdmp += distance
    else:
        # for now only consecutive exons
        # TODO: add non-consecutive exons
        consecutiveExons = False

        # check in all transcripts if the exons are consecutive
        # NOTE: two exons are consecutive if they appear next to each other
        #  in a transcript (an exon is represented as (startPosReference, endPosReference) )
        for _,exons in transcripts.items():
            start1Reference = getStart(m1,bv,exPos)
            start2Reference = getStart(m2,bv,exPos)
            end1Reference = getEnd(m1,bv,exPos)
            end2Reference = getEnd(m2,bv,exPos)
            for (s1, e1), (s2, e2) in pairwise(exons):
                if s1 == start1Reference and e1 == end1Reference and s2 == start2Reference and e2 == end2Reference:
                    consecutiveExons = True
                    exon1EndingPosBv = bv.select(id1+1)
                    exon2startingPosBv = bv.select(id2)
                    distance = exon1EndingPosBv - (start1 + len1) + (start2 - exon2startingPosBv) - 1
                    tIdmp += distance
                    break

        if not consecutiveExons:
            pass

    return tIdmp

# Extracts transcripts from gtf
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
    return transcripts

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

# L -> (l0,l1), (l1,l2), (l2, l3), ...
def pairwise(L):
    L0, L1 = itertools.tee(L)
    next(L1, None)
    return zip(L0,L1)


def main(memsPath1, memsPath2, refPath, gtfPath, errRate, outPath):
    RefSeq = list(SeqIO.parse(refPath, "fasta"))[0]

    ref, refLen, text, exPos = extractFromInfoFile(gtfPath + ".sg")
    bv = BitVector(text)

    # Reading annotation
    gtf = openGTF(gtfPath)
    transcripts = extractFromGTF(gtf)

    # Open output sam file
    out = open(outPath, "w")
    out.write("@HD\tVN:1.4\n")
    out.write("@SQ\tSN:{}\tLN:{}\n".format(ref, refLen))

    # Open output stats file
    out_stats = open(outPath+".alignsinfo.txt", "w")

    lastMapped1 = False
    lastID1 = ""
    lastStart1 = -1
    lastCigar1 = ""
    lastStrand1 = ""

    lastMapped2 = False
    lastID2 = ""
    lastStart2 = -1
    lastCigar2 = ""
    lastStrand2 = ""
    file1 = open(memsPath1)
    file2 = open(memsPath2)

    line1 = file1.readline()
    line2 = file2.readline()

    count_reads1 = 0
    count_reads2 = 0

    idmp = 0
    tIdmp = 0
    count_mapped_pairs = 0
    count_mapped1 = 0
    count_mapped2 = 0
    count_primary_allignments = 0
    count_secondary_allignments = 0
    count_unmapped_reads1 = 0
    count_unmapped_reads2 = 0

    pos_tlen = []

    # 'merge' the 2 files containing mems together
    # NOTE-1: reads in .mem files are sorted!!!
    # NOTE-2: the .mem files have the same length
    while line1!='' and line2!='':
        placeholder1, mapped1, strand1, readID1, err1, mems1, read1 = readLine(line1)
        placeholder2, mapped2, strand2, readID2, err2, mems2, read2 = readLine(line2)

        rnext1 = "*"
        pnext1 = 0
        tlen1 = 0
        flag1 = 0
        if not placeholder1:
            mems1 = extractMEMs(mems1)
            start1 = getStart(mems1[0], bv, exPos) if mapped1 else 0
            flag1 = getFlagPaired(mapped1, strand1, readID1, mapped2, strand2, readID2, read1=True)
            cigar1 = getCIGAR(mems1, RefSeq, bv, exPos, read1, errRate, err1)

        rnext2 = "*"
        pnext2 = 0
        tlen2 = 0
        flag2 = 0
        if not placeholder2:
            mems2 = extractMEMs(mems2)
            start2 = getStart(mems2[0], bv, exPos) if mapped2 else 0
            flag2 = getFlagPaired(mapped1, strand1, readID1, mapped2, strand2, readID2, read1=False)
            cigar2 = getCIGAR(mems2, RefSeq, bv, exPos, read2, errRate, err2)


        '''
            SEE: https://samtools.github.io/hts-specs/SAMv1.pdf
            For a unmapped paired-end or mate-pair read whose mate is mapped, the unmapped read should
            have RNAME and ***POS*** identical to its mate.
        '''

        if flag1 == 69:             # not mapped1 and mapped2
            start1 = start2
            #pnext2 = start2
            pnext2 = "*"
        if flag2 == 133:            # mapped1 and not mapped2
            start2 = start1
            #pnext1 = start1
            pnext1 = "*"

        # calculate rnext, pnext and tlen
        if mapped1:
            pnext1 = start2
            if mapped2:
                rnext1 = "="
                if len(mems2) == 1:
                    tlen1 = abs(start2 + len(read2) - start1)
                else:
                    tlen1 = getTlen(start1, mems2, bv, exPos)
            else:
                tlen1 = 0
        if mapped2:
            pnext2 = start1
            if mapped1:
                rnext2 = "="
                tlen2 = -tlen1
            else:
                tlen2 = 0

        # Calculate stats
        if not placeholder1:
            count_reads1 += 1
        if not placeholder2:
            count_reads2 += 1

        if mapped1 and mapped2:
            count_mapped1 += 1
            count_mapped2 += 1
            if len(mems1) == 1:
                idmp += abs(start2 - start1 + len(read1))
            else:
                idmp += getIdmp(start2, mems1, bv, exPos)
            count_mapped_pairs += 1
            if readID1 != lastID1 and readID2 != lastID2:
                count_primary_allignments += 1
                pos_tlen.append(tlen1)
            else:
                count_secondary_allignments += 1
            tIdmp += getTranscriptIdmp(transcripts, mems1, mems2, bv, exPos)
        elif mapped1: # not mapped2
            count_mapped1 += 1
            if not placeholder2:
                count_unmapped_reads2 += 1
        elif mapped2: # not mapped1
            count_mapped2 += 1
            if not placeholder1:
                count_unmapped_reads1 += 1
        else: # not mapped1 and not mapped2
            if not placeholder1:
                count_unmapped_reads1 += 1
            if not placeholder2:
                count_unmapped_reads2 += 1

        #Same alignment is not output twice
        if (readID1 != lastID1 or start1 != lastStart1 or cigar1 != lastCigar1) and not placeholder1:
            lastMapped1 = mapped1
            lastID1 = readID1
            lastStart1 = start1
            lastCigar1 = cigar1
            lastStrand1 = strand1
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\n".format(readID1, flag1, ref, start1, 255, cigar1, rnext1, pnext1, tlen1, read1, "*", err1))
        #Same alignment is not output twice
        if (readID2 != lastID2 or start2 != lastStart2 or cigar2 != lastCigar2) and not placeholder2:
            lastMapped2 = mapped2
            lastID2 = readID2
            lastStart2 = start2
            lastCigar2 = cigar2
            lastStrand2 = strand2
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\n".format(readID2, flag2, ref, start2, 255, cigar2, rnext2, pnext2, tlen2, read2, "*", err2))

        line1 = file1.readline()
        line2 = file2.readline()

    # all mems have been read

    # find the average
    idmp /= count_mapped_pairs
    tIdmp /= count_mapped_pairs

    out_stats.write("Count mapped1: " + str(count_mapped1) + "/" + str(count_reads1) + "\n")
    out_stats.write("Count mapped2: " + str(count_mapped2) + "/" + str(count_reads2) + "\n")
    out_stats.write("Count unmapped reads1: " + str(count_unmapped_reads1) + "\n")
    out_stats.write("Count unmapped reads2: " + str(count_unmapped_reads2) + "\n")
    out_stats.write("Count mapped pairs: " + str(count_mapped_pairs) + "\n")
    out_stats.write("Count primary allignments: " + str(count_primary_allignments) + "\n")
    out_stats.write("Count secondary allignments: " + str(count_secondary_allignments) + "\n")
    out_stats.write("idmp: " + str(idmp) + "\n")
    out_stats.write("tidmp: " + str(tIdmp) + "\n")
    out_stats.write("average tlen: " + str(sum(pos_tlen)/count_primary_allignments) + "\n")

    out_stats.close()
    out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Converts alignments to a splicing graph to alignments to a reference genome (in SAM format)")
    parser.add_argument('-g', '--genome', required=True, help='FASTA input file containing the reference')
    parser.add_argument('-a', '--annotation', required=True, help='GTF input file containing the gene annotation')
    parser.add_argument('-1', '--mems1', required=True, help='input file containing the alignments to the splicing graph')
    parser.add_argument('-2', '--mems2', required=True, help='input file containing the alignments to the splicing graph')
    parser.add_argument('-o', '--output', required=True, help='SAM output file')
    parser.add_argument('-e', '--erate', required=False, default=3, type=int, help='error rate of alignments (from 0 to 100, default: 3)')
    args = parser.parse_args()
    memsPath1 = args.mems1
    memsPath2 = args.mems2
    refPath = args.genome
    gtfPath = args.annotation
    errRate = args.erate
    outPath = args.output
    main(memsPath1, memsPath2, refPath, gtfPath, errRate, outPath)
