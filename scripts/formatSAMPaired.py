#!/usr/bin/python3

import sys, os
import argparse
from Bio import SeqIO

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
    return strand, readID, err, mems, read

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

def getFlagPaired(strand1, readID1, strand2, readID2, read1=True):
    f = 1 #since it's a paired-end read

    # SEE: http://seqanswers.com/forums/showthread.php?t=17314

    if readID1 == readID2:  # Both reads mapped
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

        else:   #setting flag for read coming from file 2
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

    else:  #Other read not mapped
        if read1:
            f = 73
        else:
            f = 137
    # NOTE: There is no "both reads umapped" case since we are starting from mems

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
    while i<len(mems):
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
    finalDels = m - mems[-1][1] - mems[-1][2] + 1
    if  finalDels != 0:
        cigarList.append([finalDels, 'S'])

    CIGAR = ""
    for (n,l) in cigarList:
        CIGAR += str(n) + l
    return CIGAR

def main(memsPath1, memsPath2, refPath, gtfPath, errRate, outPath):
    RefSeq = list(SeqIO.parse(refPath, "fasta"))[0]

    ref, refLen, text, exPos = extractFromInfoFile(gtfPath + ".sg")
    bv = BitVector(text)

    out = open(outPath, "w")
    out.write("@HD\tVN:1.4\n")
    out.write("@SQ\tSN:{}\tLN:{}\n".format(ref, refLen))

    lastID1 = ""
    lastStart1 = -1
    lastCigar1 = ""
    lastStrand1 = ""

    lastID2 = ""
    lastStart2 = -1
    lastCigar2 = ""
    lastStrand2 = ""

    file1 = open(memsPath1)
    file2 = open(memsPath2)

    line1 = file1.readline()
    line2 = file2.readline()

    # 'merge' the 2 files containing mems together
    # NOTE-1: reads in .mem files are sorted!!!
    # NOTE-2: .mems file may NOT have the same length!!!
    while line1!='' and line2!='':
        strand1, readID1, err1, mems1, read1 = readLine(line1)
        strand2, readID2, err2, mems2, read2 = readLine(line2)
        mems1 = extractMEMs(mems1)
        mems2 = extractMEMs(mems2)

        start1 = getStart(mems1[0], bv, exPos)
        start2 = getStart(mems2[0], bv, exPos)
        flag1 = getFlagPaired(strand1, readID1, strand2, readID2, read1=True)
        flag2 = getFlagPaired(strand1, readID1, strand2, readID2, read1=False)
        cigar1 = getCIGAR(mems1, RefSeq, bv, exPos, read1, errRate, err1)
        cigar2 = getCIGAR(mems2, RefSeq, bv, exPos, read2, errRate, err2)

        #Same alignment is not output twice
        if readID1 != lastID1 or start1 != lastStart1 or cigar1 != lastCigar1:
            lastID1 = readID1
            lastStart1 = start1
            lastCigar1 = cigar1
            lastStrand1 = strand1
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\n".format(readID1, flag1, ref, start1, 255, cigar1, "*", 0, 0, read1, "*", err1))
        #Same alignment is not output twice
        if readID2 != lastID2 or start2 != lastStart2 or cigar2 != lastCigar2:
            lastID2 = readID2
            lastStart2 = start2
            lastCigar2 = cigar2
            lastStrand2 = strand2
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\n".format(readID2, flag2, ref, start2, 255, cigar2, "*", 0, 0, read2, "*", err2))

        # keep looping over the first .mem file until ReadID1 changes
        line1 = file1.readline();
        if line1 != '':
            strand1, readID1, err1, mems1, read1 = readLine(line1)
        while line1 != '' and readID1 == lastID1:
            mems1 = extractMEMs(mems1)
            start1 = getStart(mems1[0], bv, exPos)
            flag1 = getFlagPaired(strand1, readID1, strand2, readID2, read1=True)
            cigar1 = getCIGAR(mems1, RefSeq, bv, exPos, read1, errRate, err1)

            if readID1 != lastID1 or start1 != lastStart1 or cigar1 != lastCigar1:
                lastID1 = readID1
                lastStart1 = start1
                lastCigar1 = cigar1
                lastStrand1 = strand1
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\n".format(readID1, flag1, ref, start1, 255, cigar1, "*", 0, 0, read1, "*", err1))

            line1 = file1.readline()
            strand1, readID1, err1, mems1, read1 = readLine(line1)
            # in order to end the while loop, readID1 (and therefore all the other -1 variables) have to change
            # their previous values will be in the last- variables

        # keep looping over the second .mem file until ReadID2 changes
        line2 = file2.readline();
        if line2 != '':
            strand2, readID2, err2, mems2, read2 = readLine(line2)
        while line2 != '' and readID2 == lastID2:
            mems2 = extractMEMs(mems2)
            start2 = getStart(mems2[0], bv, exPos)
            flag2 = getFlagPaired(lastStrand1, lastID1, strand2, readID2, read1=False) #lastStrand1 and lastID1 are needed since strand1 could have changed
            cigar2 = getCIGAR(mems2, RefSeq, bv, exPos, read2, errRate, err2)

            if readID2 != lastID2 or start2 != lastStart2 or cigar2 != lastCigar2:
                lastID2 = readID2
                lastStart2 = start2
                lastCigar2 = cigar2
                lastStrand2 = strand2
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\n".format(readID2, flag2, ref, start2, 255, cigar2, "*", 0, 0, read2, "*", err2))

            line2 = file2.readline()
            strand2, readID2, err2, mems2, read2 = readLine(line2)
            # in order to end the while loop, readID2 (and therefore all the other -2 variables) have to change
            # their previous values will be in the last- variables

        # by this point, all reads having the same readID in both files have been written in the SAM file


    # by this point, either line1 or line2 will be = '' (or both if both .mem files have the same length)

    # -1.mem file not yet finished
    while line1!='':
        strand1, readID1, err1, mems1, read1 = readLine(line1)
        mems1 = extractMEMs(mems1)

        start1 = getStart(mems1[0], bv, exPos)
        flag1 = getFlagPaired(strand1, readID1, lastStrand2, lastID2, read1=True) #lastStrand2 and lastID2 are needed since strand1 could have changed

        cigar = getCIGAR(mems1, RefSeq, bv, exPos, read1, errRate, err1)
        if readID1 != lastID1 or start1 != lastStart1 or cigar1 != lastCigar1:
            lastID1 = readID1
            lastStart1 = start1
            lastCigar1 = cigar1
            lastStrand1 = strand1
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\n".format(readID1, flag1, ref, start1, 255, cigar1, "*", 0, 0, read1, "*", err1))
        line1 = file1.readline()

    # -2.mem file not yet finished
    while line2!='':
        strand2, readID2, err2, mems2, read2 = readLine(line2)
        mems2 = extractMEMs(mems2)

        start2 = getStart(mems2[0], bv, exPos)
        flag2 = getFlagPaired(lastStrand1, lastID1, strand2, readID2, read1=False)

        cigar = getCIGAR(mems2, RefSeq, bv, exPos, read2, errRate, err2)
        if readID2 != lastID2 or start2 != lastStart2 or cigar2 != lastCigar2:
            lastID2 = readID2
            lastStart2 = start2
            lastCigar2 = cigar2
            lastStrand2 = strand2
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\n".format(readID2, flag2, ref, start2, 255, cigar2, "*", 0, 0, read2, "*", err2))
        line2 = file2.readLine()

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
