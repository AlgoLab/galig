import sys
#from Bio import SeqIO

from BitVector import BitVector

class SAMFormatter:
    def __init__(self, outFile, indexFile):#, rna_seqs_file):
        self.outFile = outFile

        #Output reader
        self.outs = []
        for line in open(outFile).readlines():
            line = line.strip('\n')
            if line != "":
                l = line.split(" ")
                self.outs.append([l[0], l[1], l[2], l[3:-1], l[-1]])

        #Index reader
        with open(indexFile) as o:
            line = o.readline()
            i = 0
            while line:
                if i==0:
                    (self.reference, self.ref_length, self.strand) = line[:-1].split(" ")
                if i==1:
                    self.text = line[:-1]
                if i==4:
                    pos_S = line[:-2].split(" ")
                    self.pos = []
                    for p in pos_S:
                        x = p.split(",")
                        self.pos.append([int(x[0]), int(x[1])])
                i+=1
                line = o.readline()

        #Rna-Seqs extraction
        #self.rna_seqs = SeqIO.index(rna_seqs_file, "fasta")

        #Bit Vector Setup
        self.bv = BitVector(self.text)

    def format(self):
        out = open(self.outFile + ".sam", "w")
        out.write("@HD\tVN:1.4\n")
        out.write("@SQ\tSN:{}\tLN:{}\n".format(self.reference, self.ref_length))
        lastID = ""
        lastStart = ""
        lastCigar = ""
        for (strand, rID, err, mems, read) in self.outs:
            memsList = self.extractMEMs(mems)
            if strand == "-":
                if lastID == rID:
                    f = 272
                else:
                    f = 16
                #read = self.reverse_and_complement(read)
            else:
                if lastID == rID:
                    f=256
                else:
                    f = 0
            start = self.getStart(memsList[0])
            cigar = self.getCIGAR(memsList, len(read))
            #Same alignment is not output twice
            if rID != lastID or start != lastStart or cigar != lastCigar:
                lastID = rID
                lastStart = start
                lastCigar = cigar
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNM:i:{}\n".format(rID, f, self.reference, start, 255, cigar, "*", 0, 0, read, "*", err))
        out.close()

    #Utils
    def extractMEMs(self, mems_S):
        mems = []
        for s in mems_S:
            string_mem = s[1:-1].split(",")
            mems.append([int(string_mem[0]), int(string_mem[1]), int(string_mem[2])])
        return mems

    def getStart(self, mem):
        id = self.bv.rank(mem[0])
        return self.pos[id-1][0] + (mem[0] - self.bv.select(id)) - 1

    def getEnd(self, mem):
        id = self.bv.rank(mem[0])
        return self.pos[id-1][0] + (mem[0] + mem[2] - 1 - self.bv.select(id)) - 1

    def reverse_and_complement(self, text):
        text = text[::-1]
        nucl_bases_dict = {"A":"T", "T":"A", "C":"G", "G":"C", "a":"t", "t":"a", "c":"g", "g":"c", "N":"N", "n":"n"}
        new_text = ""
        for elem in text:
            new_text += nucl_bases_dict[elem]
        return new_text

    def getCIGAR(self, mems, m):
        CIGAR = ""
        i = 0
        while i<len(mems):
            if i == 0:
                if mems[i][1] != 1:
                    initial_clips = mems[i][1]-1
                    CIGAR += "{}S".format(initial_clips)
                    poss_text = self.text[mems[i][0]-initial_clips-1:mems[i][0]-1]
                CIGAR += "{}M".format(mems[i][2])
            else:
                ##################################################################
                id1 = self.bv.rank(mems[i-1][0])
                id2 = self.bv.rank(mems[i][0])
                if id1 == id2:
                    errors_P = mems[i][1] - mems[i-1][1] - mems[i-1][2]
                    errors_T = mems[i][0] - mems[i-1][0] - mems[i-1][2]
                    #------------------------------------------------
                    if errors_P == 0:
                        if errors_T < 0:
                            #0 on P, - on T
                            CIGAR += "{}I".format(abs(errors_T))
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_T))
                        elif errors_T > 0:
                            #0 on P, + on T
                            CIGAR += "{}D".format(abs(errors_T))
                            CIGAR += "{}M".format(mems[i][2])
                    #------------------------------------------------
                    elif errors_P < 0:
                        if errors_T == 0:
                            #- on P, 0 on T
                            CIGAR += "{}D".format(abs(errors_P))
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                        elif errors_T < 0:
                            #- on P, - on T
                            errs = abs(errors_P) - abs(errors_T)
                            if errs > 0:
                                #More on P
                                CIGAR += "{}D".format(errs)
                            elif errs < 0:
                                #More on T
                                CIGAR += "{}I".format(abs(errs))
                            CIGAR += "{}M".format(mems[i][2]-max(abs(errors_P), abs(errors_T)))
                        elif errors_T > 0:
                            #- on P, + on T
                            CIGAR += "{}D".format(abs(errors_P) + errors_T)
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                    #------------------------------------------------
                    elif errors_P > 0:
                        if errors_T == 0:
                            #+ on P, 0 on T
                            CIGAR += "{}I".format(errors_P)
                            CIGAR += "{}M".format(mems[i][2])
                        elif errors_T < 0:
                            #+ on P, - on T
                            CIGAR += "{}I".format(abs(errors_P + abs(errors_T)))
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_T))
                        elif errors_T > 0:
                            #+ on P, + on T
                            if errors_P == errors_T:
                                #Same: all SUBs
                                prev_matches = 0
                                f1 = True
                                j1 = 0
                                while f1:
                                    j1-=1
                                    if j1<=-len(CIGAR) or CIGAR[:-1][j1] in ["I", "D", "N", "S"]:
                                        f1 = False
                                        prev_matches += int(CIGAR[:-1][j1+1:])
                                CIGAR = CIGAR[:j1] + "{}M".format(prev_matches + errors_P + mems[i][2])
                            elif errors_P > errors_T:
                                #More on P: some INSs, some SUBs
                                CIGAR += "{}I".format(errors_P - errors_T)
                                CIGAR += "{}M".format(errors_T + mems[i][2])
                            else:
                                #More on T: some DELs, some SUBs
                                CIGAR += "{}D".format(errors_T - errors_P)
                                CIGAR += "{}M".format(errors_P + mems[i][2])
                ##################################################################
                else:
                    errors_P = mems[i][1] - mems[i-1][1] - mems[i-1][2]
                    errors_T1 = self.bv.select(self.bv.rank(mems[i-1][0]) + 1) - mems[i-1][0] - mems[i-1][2]
                    errors_T2 = mems[i][0] - self.bv.select(self.bv.rank(mems[i][0])) - 1
                    intron = self.pos[id2-1][0] - self.pos[id1-1][1] - 1
                    #------------------------------------------------
                    if errors_P == 0:
                        CIGAR += "{}N".format(intron + errors_T1 + errors_T2)
                        CIGAR += "{}M".format(mems[i][2])
                    #------------------------------------------------
                    elif errors_P < 0:
                        CIGAR += "{}N".format(intron + errors_T1 + errors_T2 + abs(errors_P))
                        CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                    #------------------------------------------------
                    else:
                        CIGAR += "{}I".format(abs(errors_P))
                        CIGAR += "{}N".format(intron + errors_T1 + errors_T2)
                        CIGAR += "{}M".format(mems[i][2])
            i+=1
        final_dels = m - mems[-1][1] - mems[-1][2] + 1
        if  final_dels != 0:
            CIGAR += "{}S".format(final_dels)
        return CIGAR


if __name__ == '__main__':
    outFile = sys.argv[1]
    indexFile = sys.argv[2]
    #reads = sys.argv[3]
    sf = SAMFormatter(outFile, indexFile)
    sf.format()
