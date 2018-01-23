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

def getFlag(strand, readID, lastID):
    f = 0
    if strand == "-":
        if lastID == readID:
            f = 272
        else:
            f = 16
    else:
        if lastID == readID:
            f=256
        else:
            f = 0
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
                            cigarList.append([abs(errors_P) + errors_T, 'N'])
                            cigarList.append([mems[i][2], 'M'])
                        else:
                            cigarList.append([abs(errors_P) + errors_T, 'N'])
                            cigarList.append([mems[i][2]-abs(errors_P), 'M'])
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
                    cigarList.append([intron + errors_T1 + errors_T2, 'N'])
                    cigarList.append([mems[i][2], 'M'])
                #------------------------------------------------
                elif errors_P < 0:
                    if errors_T1 == 0 and errors_T2 != 0:
                        cigarList.append([intron + errors_T1 + errors_T2 + abs(errors_P), 'N'])
                        cigarList.append([mems[i][2]-abs(errors_P), 'M'])
                    elif errors_T1 != 0 and errors_T2 == 0:
                        cigarList[-1][0] = cigarList[-1][0] - abs(errors_P)
                        cigarList.append([intron + errors_T1 + errors_T2 + abs(errors_P), 'N'])
                        cigarList.append([mems[i][2], 'M'])
                    else:
                        I = (exPos[id1-1][0]+mems[i-1][0]+mems[i-1][2]-bv.select(id1)-1 - abs(errors_P), exPos[id1-1][0]+mems[i-1][0]+mems[i-1][2]-bv.select(id1) + intron+errors_T1+errors_T2-2)
                        Iseq = RefSeq[I[0]-1:I[1]].seq
                        if Iseq[:2] in ["GT", "GC"] and Iseq[-2:] == "AG":
                            cigarList[-1][0] = cigarList[-1][0] - abs(errors_P)
                            cigarList.append([intron + errors_T1 + errors_T2 + abs(errors_P), 'N'])
                            cigarList.append([mems[i][2], 'M'])
                        else:
                            cigarList.append([intron + errors_T1 + errors_T2 + abs(errors_P), 'N'])
                            cigarList.append([mems[i][2]-abs(errors_P), 'M'])
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
                            cigarList.append([N, 'N'])
                            cigarList.append([mems[i][2], 'M'])
                    elif err2 <= err1 and err2 + err <= maxErr:
                        N = intron + errors_T1 + errors_T2 - len(readGapString)
                        if N <= 0:
                            cigarList[-1][0] = cigarList[-1][0] + len(readGapString) + mems[i][2]
                        else:
                            cigarList.append([N, 'N'])
                            cigarList.append([len(readGapString)+mems[i][2], 'M'])
                    else:
                        cigarList.append([abs(errors_P), 'I'])
                        cigarList.append([intron + errors_T1 + errors_T2, 'N'])
                        cigarList.append([mems[i][2], 'M'])
        i+=1
    finalDels = m - mems[-1][1] - mems[-1][2] + 1
    if  finalDels != 0:
        cigarList.append([finalDels, 'S'])

    CIGAR = ""
    for (n,l) in cigarList:
        CIGAR += str(n) + l
    return CIGAR

def main(memsPath, refPath, gtfPath, errRate, outPath):
    RefSeq = list(SeqIO.parse(refPath, "fasta"))[0]

    ref, refLen, text, exPos = extractFromInfoFile(gtfPath + ".sg")
    bv = BitVector(text)

    out = open(outPath, "w")
    out.write("@HD\tVN:1.4\n")
    out.write("@SQ\tSN:{}\tLN:{}\n".format(ref, refLen))

    lastID = ""
    lastStart = -1
    lastCigar = ""
    for line in open(memsPath).readlines():
        strand, readID, err, mems, read = readLine(line)
        mems = extractMEMs(mems)

        start = getStart(mems[0], bv, exPos)
        flag = getFlag(strand, readID, lastID)

        cigar = getCIGAR(mems, RefSeq, bv, exPos, read, errRate, err)
        #Same alignment is not output twice
        if readID != lastID or start != lastStart or cigar != lastCigar:
            lastID = readID
            lastStart = start
            lastCigar = cigar
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\n".format(readID, flag, ref, start, 255, cigar, "*", 0, 0, read, "*", err))

    out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Converts alignments to a splicing graph to alignments to a reference genome (in SAM format)")
    parser.add_argument('-g', '--genome', required=True, help='FASTA input file containing the reference')
    parser.add_argument('-a', '--annotation', required=True, help='GTF input file containing the gene annotation')
    parser.add_argument('-m', '--mems', required=True, help='input file containing the alignments to the splicing graph')
    parser.add_argument('-o', '--output', required=True, help='SAM output file')
    parser.add_argument('-e', '--erate', required=False, default=3, type=int, help='error rate of alignments (from 0 to 100, default: 3)')
    args = parser.parse_args()
    memsPath = args.mems
    refPath = args.genome
    gtfPath = args.annotation
    errRate = args.erate
    outPath = args.output
    main(memsPath, refPath, gtfPath, errRate, outPath)
