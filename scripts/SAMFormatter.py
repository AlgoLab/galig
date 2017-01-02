import sys
from Bio import SeqIO

from BitVector import BV

class SAMFormatter:
    def __init__(self, out_file, rna_seqs_file):
        self.out_file = out_file
        #Output reader
        with open(out_file) as o:
            self.outs = []
            for line in o.read().split("\n"):
                if line != "":
                    l = line.split(" ")
                    self.outs.append((l[0],[self.extractMEM(m) for m in l[1:-2]]))

        #Gene_name extraction
        with open("./tmp/gene_name") as g:
            self.gene_name = g.read()

        #Rna-Seqs extraction
        ### !!!BAD!!!
        rna_seqs = list(SeqIO.parse(rna_seqs_file, "fasta"))
        self.rna_seqs = {}
        for elem in rna_seqs:
            self.rna_seqs.update({elem.id:elem.seq})
        print(self.rna_seqs)
        
        #Bit Vector Setup
        with open("./tmp/T.fa") as t:
            text = t.read().split("\n")[1]
            self.bv = BV(text)

        #Exons Position Setup
        with open("./tmp/e_pos") as f:
            self.exs_pos = []
            for line in f.read().split("\n"):
                if line != "":
                    self.exs_pos.append([int(n) for n in line.split(",")])

    def format(self):
        out = open(self.out_file + ".sam", "w")
        out.write("@HD\n")
        out.write("@SQ\n")
        for (p_id, mems) in self.outs:
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(p_id, 0, self.gene_name, self.getStart(mems[0]), 255, self.getCIGAR(mems), "*", 0, 0, self.rna_seqs[p_id], "*"))

    #Utils
    def extractMEM(self, string):
        string_mem = string[1:-1].split(",")
        return [int(string_mem[0]), int(string_mem[1]), int(string_mem[2])]

    def getStart(self, mem):
        id = self.bv.rank(mem[0])
        return self.exs_pos[id-1][0] + (mem[0] - self.bv.select(id))

    def getCIGAR(self, mems):
        CIGAR = ""
        i = 0
        while i<len(mems):
            if i == 0:
                if mems[i][1] != 1:
                    CIGAR += "{}D".format(mems[i][1]-1)
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
                            CIGAR += "{}D".format(abs(errors_T))
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_T))
                        elif errors_T > 0:
                            #Case 2
                            #print("2")
                            CIGAR += "{}I".format(abs(errors_T))
                            CIGAR += "{}M".format(mems[i][2])
                    #------------------------------------------------
                    elif errors_P < 0:
                        if errors_T == 0:
                            #Case 3
                            #print("3")
                            CIGAR += "{}I".format(abs(errors_P))
                            CIGAR += "{}M".format(mems[i][2-abs(errors_P)])
                        elif errors_T < 0:
                            #Case 4
                            errs = abs(errors_P) - abs(errors_T)
                            if errs > 0:
                                #Case 4a
                                #print("4a")
                                CIGAR += "{}I".format(errs)
                            elif errs < 0:
                                #Case 4b
                                #print("4b")
                                CIGAR += "{}D".format(abs(errs))
                                CIGAR += "{}M".format(mems[i][2]-max(abs(errors_P), abs(errors_T)))
                        elif errors_T > 0:
                            #Case 5
                            #print("5")
                            CIGAR += "{}I".format(abs(errors_P) + errors_T)
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                    #------------------------------------------------
                    elif errors_P > 0:
                        if errors_T == 0:
                            #Case 6
                            #print("6")
                            CIGAR += "{}D".format(errors_P)
                            CIGAR += "{}M".format(mems[i][2])
                        elif errors_T < 0:
                            #Case 7
                            #print("7")
                            CIGAR += "{}D".format(abs(errors_P + abs(errors_T)))
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_T))
                        elif errors_T > 0:
                            #Case 8
                            CIGAR += "{}D".format(abs(errors_P))
                            CIGAR += "{}I".format(abs(errors_T))
                            CIGAR += "{}M".format(mems[i][2])
                ##################################################################
                else:
                    errors_P = mems[i][1] - mems[i-1][1] - mems[i-1][2]
                    errors_T1 = self.bv.select(self.bv.rank(mems[i-1][0]) + 1) - mems[i-1][0] - mems[i-1][2]
                    errors_T2 = mems[i][0] - self.bv.select(self.bv.rank(mems[i][0])) - 1
                    intron = self.exs_pos[id2-1][0] - self.exs_pos[id1-1][1]
                    #------------------------------------------------
                    if errors_P == 0:
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
                                    if CIGAR[:-1][j1] in ["I", "D", "N"]:
                                        f1 = False
                                        prev_matches += int(CIGAR[:-1][j1+1:])
                                CIGAR = CIGAR[:j1] + "{}M".format(prev_matches + mems[i][2])
                        elif errors_T1 > 0 and errors_T2 == 0:
                            #print("10")
                            CIGAR += "{}I".format(errors_T1)
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}M".format(mems[i][2])
                        elif errors_T1 == 0 and errors_T2 > 0:
                            #print("11")
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}I".format(errors_T2)
                            CIGAR += "{}M".format(mems[i][2])
                        else:
                            #print("12")
                            CIGAR += "{}I".format(errors_T1)
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}I".format(errors_T2)
                            CIGAR += "{}M".format(mems[i][2])
                    #------------------------------------------------
                    elif errors_P < 0:
                        if errors_T1 == 0 and errors_T2 == 0:
                            #print("13")
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}I".format(abs(errors_P))
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                        elif errors_T1 > 0 and errors_T2 == 0:
                            #print("14")
                            CIGAR += "{}I".format(abs(errors_T1))
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}I".format(abs(errors_P))
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                        elif errors_T1 == 0 and errors_T2 > 0:
                            #print("15")
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}I".format(abs(errors_P) + errors_T2)
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                        else:
                            #print("16")
                            CIGAR += "{}I".format(abs(errors_T1))
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}I".format(abs(errors_P) + errors_T2)
                            CIGAR += "{}M".format(mems[i][2]-abs(errors_P))
                    #------------------------------------------------
                    else:
                        if errors_T1 == 0 and errors_T2 == 0:
                            #print("17")
                            CIGAR += "{}D".format(errors_P)
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}M".format(mems[i][2])
                        elif errors_T1 > 0 and errors_T2 == 0:
                            #print("18")
                            CIGAR += "{}D".format(errors_P)
                            CIGAR += "{}I".format(errors_T1)
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}M".format(mems[i][2])
                        elif errors_T1 == 0 and errors_T2 > 0:
                            #print("19")
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}D".format(errors_P)
                            CIGAR += "{}I".format(errors_T2)
                            CIGAR += "{}M".format(mems[i][2])
                        else:
                            #print("20")
                            CIGAR += "{}D".format(errors_P)
                            CIGAR += "{}I".format(errors_T1)
                            if intron != 0:
                                CIGAR += "{}N".format(intron)
                            CIGAR += "{}I".format(errors_T2)
                            CIGAR += "{}M".format(mems[i][2])
            i+=1
        return CIGAR

if __name__ == '__main__':
    #Outfile, rna_seqs
    f = SAMFormatter(sys.argv[1], sys.argv[2])
    f.format()
