import sys
from Bio import SeqIO

from BitVector import BV

class SAMFormatter:
    def __init__(self, out_file, index_file, rna_seqs_file):
        self.out_file = out_file

        #Output reader
        with open(out_file) as o:
            self.outs = []
            for line in o.read().split("\n"):
                if line != "":
                    l = line.split(" ")
                    l.remove("")
                    self.outs.append([l[0], l[1], l[2], l[3:]])

        #Index reader
        with open(index_file) as o:
            line = o.readline()
            i = 0
            while line:
                if i==0:
                    (self.reference, self.ref_length) = line[:-1].split(" ")
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
        self.rna_seqs = SeqIO.index(rna_seqs_file, "fasta")

        #Bit Vector Setup
        self.bv = BV(self.text)

    def format(self):
        out_mems = open(self.out_file + "_in_SAM", "w")
        out = open(self.out_file + ".sam", "w")
        out.write("@HD\tVN:1.4\n")
        out.write("@SQ\tSN:{}\tLN:{}\n".format(self.reference, self.ref_length))
        last_start = ""
        last_cigar = ""
        for (strand, p_id, err, mems) in self.outs:
            mems_list = self.extractMEMs(mems)
            rna_seq = self.rna_seqs[p_id].seq
            if strand == "-":
                f = 16
                rna_seq = self.reverse_and_complement(rna_seq)
            else:
                f = 0
            start = self.getStart(mems_list[0])
            end = self.getEnd(mems_list[-1])
            cigar = self.getCIGAR(mems_list, len(rna_seq))
            #Same alignments is not output
            if start != last_start or cigar != last_cigar:
                last_start = start
                last_cigar = cigar
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tER:A:{}\tFP:A:{}\n".format(p_id, f, self.reference, start, 255, cigar, "*", 0, 0, rna_seq, "*", err, end))
                out_mems.write("{} {} {} {}\n".format(strand, p_id, err, ' '.join(mems)))
        out.close()
        out_mems.close()

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
                            #Case 1
                            #print("1")
                            CIGAR += "{}I".format(abs(errors_T))
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_T))
                        elif errors_T > 0:
                            #Case 2
                            #print("2")
                            CIGAR += "{}D".format(abs(errors_T))
                            CIGAR += "{}M".format(mems[i][2])
                    #------------------------------------------------
                    elif errors_P < 0:
                        if errors_T == 0:
                            #Case 3
                            #print("3")
                            CIGAR += "{}D".format(abs(errors_P))
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                        elif errors_T < 0:
                            #Case 4
                            errs = abs(errors_P) - abs(errors_T)
                            if errs > 0:
                                #Case 4a
                                #print("4a")
                                CIGAR += "{}D".format(errs)
                            elif errs < 0:
                                #Case 4b
                                #print("4b")
                                CIGAR += "{}I".format(abs(errs))
                            CIGAR += "{}M".format(mems[i][2]-max(abs(errors_P), abs(errors_T)))
                        elif errors_T > 0:
                            #Case 5
                            #print("5")
                            CIGAR += "{}D".format(abs(errors_P) + errors_T)
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                    #------------------------------------------------
                    elif errors_P > 0:
                        if errors_T == 0:
                            #Case 6
                            #print("6")
                            CIGAR += "{}I".format(errors_P)
                            CIGAR += "{}M".format(mems[i][2])
                        elif errors_T < 0:
                            #Case 7
                            #print("7")
                            CIGAR += "{}I".format(abs(errors_P + abs(errors_T)))
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_T))
                        elif errors_T > 0:
                            #Case 8
                            if errors_P == errors_T:
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
                                CIGAR += "{}I".format(errors_P - errors_T)
                                CIGAR += "{}M".format(errors_T + mems[i][2])
                            else:
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
                        '''
                        if errors_T1 == 0 and errors_T2 == 0:
                            #print("9")
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                                CIGAR += "{}M".format(mems[i][2])
                            else:
                                prev_matches = 0
                                f1 = True
                                j1 = 0
                                while f1:
                                    j1-=1
                                    if j1<=-len(CIGAR) or CIGAR[:-1][j1] in ["I", "D", "N", "S"]:
                                        f1 = False
                                        prev_matches += int(CIGAR[:-1][j1+1:])
                                CIGAR = CIGAR[:j1] + "{}M".format(prev_matches + mems[i][2])
                        elif errors_T1 > 0 and errors_T2 == 0:
                            #print("10")
                            CIGAR += "{}D".format(errors_T1)
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}M".format(mems[i][2])
                        elif errors_T1 == 0 and errors_T2 > 0:
                            #print("11")
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}D".format(errors_T2)
                            CIGAR += "{}M".format(mems[i][2])
                        else:
                            #print("12")
                            if intron != 0:
                                CIGAR += "{}D".format(errors_T1)
                                CIGAR += "{}N".format(intron)
                                CIGAR += "{}D".format(errors_T2)
                            else:
                                CIGAR += "{}D".format(errors_T1 + errors_T2)
                            CIGAR += "{}M".format(mems[i][2])
                        '''
                    #------------------------------------------------
                    elif errors_P < 0:
                        CIGAR += "{}D".format(abs(errors_P))
                        CIGAR += "{}N".format(intron + errors_T1 + errors_T2 - 1)
                        CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                        '''
                        if errors_T1 == 0 and errors_T2 == 0:
                            #print("13")
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}D".format(abs(errors_P))
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                        elif errors_T1 > 0 and errors_T2 == 0:
                            #print("14")
                            if intron != 0:
                                CIGAR += "{}D".format(errors_T1)
                                CIGAR += "{}N".format(intron)
                                CIGAR += "{}D".format(abs(errors_P))
                            else:
                                CIGAR += "{}D".format(errors_T1 + abs(errors_P))
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                        elif errors_T1 == 0 and errors_T2 > 0:
                            #print("15")
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}D".format(abs(errors_P) + errors_T2)
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                        else:
                            #print("16")
                            if intron != 0:
                                CIGAR += "{}D".format(errors_T1)
                                CIGAR += "{}N".format(intron)
                                CIGAR += "{}D".format(abs(errors_P) + errors_T2)
                            else:
                                CIGAR += "{}D".format(errors_T1 + abs(errors_P) + errors_T2)
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                        '''
                    #------------------------------------------------
                    else:
                        #CIGAR += "{}I".format(abs(errors_P))
                        #CIGAR += "{}N".format(intron + errors_T1 + errors_T2 - 1)
                        #CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                        '''
                        if errors_T1 + errors_T2 == errors_P:
                            prev_matches = 0
                            f1 = True
                            j1 = 0
                            while f1:
                                j1-=1
                                if j1<=-len(CIGAR) or CIGAR[:-1][j1] in ["I", "D", "N", "S"]:
                                    f1 = False
                                    prev_matches += int(CIGAR[:-1][j1+1:])
                            CIGAR = CIGAR[:j1] + "{}M".format(prev_matches + errors_T1)
                            CIGAR += "{}N".format(intron)
                            CIGAR += "{}M".format(errors_T2 + mems[i][2])
                        elif errors_T1 + errors_T2 < errors_P:
                        '''
                        if errors_T1 == 0 and errors_T2 == 0:
                            CIGAR += "{}I".format(errors_P)
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}M".format(mems[i][2])
                        elif errors_T1 > 0 and errors_T2 == 0:
                            if errors_P == errors_T1:
                                prev_matches = 0
                                f1 = True
                                j1 = 0
                                while f1:
                                    j1-=1
                                    if j1<=-len(CIGAR) or CIGAR[:-1][j1] in ["I", "D", "N", "S"]:
                                        f1 = False
                                        prev_matches += int(CIGAR[:-1][j1+1:])
                                if intron != 0:
                                    CIGAR = CIGAR[:j1] + "{}M".format(prev_matches + errors_P)
                                    CIGAR += "{}N".format(intron)
                                    CIGAR += "{}M".format(mems[i][2])
                                else:
                                    CIGAR = CIGAR[:j1] + "{}M".format(prev_matches + errors_P + mems[i][2])
                            elif errors_P > errors_T1:
                                CIGAR += "{}I".format(errors_P - errors_T1)
                                if intron != 0:
                                    CIGAR += "{}M".format(errors_T1)
                                    CIGAR += "{}N".format(intron)
                                    CIGAR += "{}M".format(mems[i][2])
                                else:
                                    CIGAR += "{}M".format(errors_T1 + mems[i][2])
                            else:
                                prev_matches = 0
                                f1 = True
                                j1 = 0
                                while f1:
                                    j1-=1
                                    if j1<=-len(CIGAR) or CIGAR[:-1][j1] in ["I", "D", "N", "S"]:
                                        f1 = False
                                        prev_matches += int(CIGAR[:-1][j1+1:])
                                CIGAR =  CIGAR[:j1] + "{}M".format(prev_matches + errors_P)
                                CIGAR += "{}D".format(errors_T1 - errors_P)
                                if intron != 0:
                                    CIGAR += "{}N".format(intron)
                                CIGAR += "{}M".format(mems[i][2])
                        elif errors_T1 == 0 and errors_T2 > 0:
                            if errors_P == errors_T2:
                                if intron != 0:
                                    CIGAR += "{}N".format(intron)
                                    CIGAR += "{}M".format(errors_P + mems[i][2])
                                else:
                                    prev_matches = 0
                                    f1 = True
                                    j1 = 0
                                    while f1:
                                        j1-=1
                                        if j1<=-len(CIGAR) or CIGAR[:-1][j1] in ["I", "D", "N", "S"]:
                                            f1 = False
                                            prev_matches += int(CIGAR[:-1][j1+1:])
                                    CIGAR = CIGAR[:j1] + "{}M".format(prev_matches + errors_P + mems[i][2])
                            elif errors_P > errors_T2:
                                if intron != 0:
                                    CIGAR += "{}N".format(intron)
                                CIGAR += "{}I".format(errors_P - errors_T2)
                                CIGAR += "{}M".format(errors_T2 + mems[i][2])
                            else:
                                if intron != 0:
                                    CIGAR += "{}N".format(intron)
                                CIGAR += "{}D".format(errors_T2 - errors_P)
                                CIGAR += "{}M".format(errors_P + mems[i][2])
                        else:
                            if errors_P == errors_T1 + errors_T2:
                                prev_matches = 0
                                f1 = True
                                j1 = 0
                                while f1:
                                    j1-=1
                                    if j1<=-len(CIGAR) or CIGAR[:-1][j1] in ["I", "D", "N", "S"]:
                                        f1 = False
                                        prev_matches += int(CIGAR[:-1][j1+1:])
                                if intron != 0:
                                    CIGAR = CIGAR[:j1] + "{}M".format(prev_matches + errors_T1)
                                    CIGAR += "{}N".format(intron)
                                    CIGAR += "{}M".format(errors_T2 + mems[i][2])
                                else:
                                    CIGAR = CIGAR[:j1] + "{}M".format(prev_matches + errors_P + mems[i][2])
                            elif errors_P > errors_T1 + errors_T2:
                                CIGAR += "{}I".format(errors_P - errors_T1 - errors_T2)
                                if intron != 0:
                                    CIGAR += "{}M".format(errors_T1)
                                    CIGAR += "{}N".format(intron)
                                    CIGAR += "{}M".format(errors_T2 + mems[i][2])
                                else:
                                    CIGAR = CIGAR[:j1] + "{}M".format(errors_T1 + errors_T2 + mems[i][2])
                            else:
                                tot_ins = errors_T1 + errors_T2 - errors_P
                                if tot_ins <= errors_T1:
                                    CIGAR += "{}D".format(tot_ins)
                                    if intron != 0:
                                        if errors_T1 - tot_ins > 0:
                                            CIGAR += "{}M".format(errors_T1 - tot_ins)
                                        CIGAR += "{}N".format(intron)
                                        CIGAR += "{}M".format(errors_T2 + mems[i][2])
                                    else:
                                        CIGAR += "{}M".format(errors_T1 - tot_ins + errors_T2 + mems[i][2])
                                else:
                                    if intron != 0:
                                        CIGAR += "{}D".format(errors_T1)
                                        tot_ins -= errors_T1
                                        CIGAR += "{}N".format(intron)
                                        CIGAR += "{}D".format(tot_ins)
                                        CIGAR += "{}M".format(errors_T2 - tot_ins + mems[i][2])
                                    else:
                                        CIGAR += "{}D".format(tot_ins)
                                        CIGAR += "{}M".format(errors_T2 - tot_ins + mems[i][2])
            i+=1
        final_dels = m - mems[-1][1] - mems[-1][2] + 1
        if  final_dels != 0:
            CIGAR += "{}S".format(final_dels)
        return CIGAR


if __name__ == '__main__':
    sf = SAMFormatter(sys.argv[1], sys.argv[2], sys.argv[3])
    sf.format()
