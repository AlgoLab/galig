import sys, os

from Bio import SeqIO
from graphviz import Digraph

from BitVector import BitVector
from utils import *

def editDistance(s1, s2):
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

class SplicingGraph:
    def __init__(self, ref, infoPath, errRate): #reads,
        infofile = open(infoPath)
        lines = infofile.readlines()
        chromo = lines[0].strip("\n").split(" ")[0]
        self.ref = list(SeqIO.parse(ref, "fasta"))[0].seq
        #self.reads = SeqIO.index(reads, "fasta")
        if lines[0].strip("\n").split(" ")[2] == '+':
            self.GeneStrand = True
        else:
            self.GeneStrand = False
        self.ExonNames = lines[5].strip("\n").split()
        self.Text = lines[1]
        self.BV = BitVector(self.Text)
        pos = [(int(p[0]), int(p[1])) for p in [pos.split(",") for pos in lines[4].strip("\n").split()]]
        self.ExonPos = dict(zip(self.ExonNames, pos))
        self.labels = []
        self.nodes = {}
        self.edges = {}
        self.newEdges = {}
        adjMatrix = [[int(elem) for elem in row.split()] for row in lines[3].strip("\n")[:-2].split(";")]
        for name in self.ExonNames:
            self.addNode(name)
        r = 0
        for row in adjMatrix:
            c = 0
            for elem in row:
                if elem == 1:
                    self.addEdge(r,c,'e')
                if elem == 2:
                    self.addEdge(r,c,'n')
                c+=1
            r+=1
        self.errRate = errRate

    #Augmenting
    def augments(self, memsPath, onlyPrimFlag):
        '''
        - comps:   {newPos: weight}
        '''
        self.IntA3 = {}
        self.IntA5 = {}
        self.ExtA3 = {}
        self.ExtA5 = {}
        self.introns = {}
        lastID = ""
        for line in open(memsPath, 'r').readlines():
            # 0: strand
            # 1: ID
            # 2: errors
            # 3 to -1: mems
            # -1: read
            info = line.strip("\n").strip(" ").split(" ")
            strand = info[0]
            readID = info[1]
            err = int(info[2])
            read = info[-1]
            usedExons = set([self.BV.rank(int(mem[1:-1].split(',')[0]) - 1) for mem in info[3:-1]])

            #Only the first alignment is considered (primary one)
            if onlyPrimFlag and lastID == readID:
                continue
            lastID = readID

            #The weight of all exons used by alignments are incremented by 1
            for e in usedExons:
                self.incrementNode(e)

            #If number of MEMs > 1
            if len(info[3:-1]) > 1:
                for currMEM, nextMEM in pairwise(info[3:-1]):
                    # Remove ( and ) from mem and cast to int
                    currMEM = [int(elem) for elem in currMEM[1:-1].split(",")]
                    currIndex = self.BV.rank(currMEM[0] - 1)
                    currLabel = self.getNodeLabel(currIndex)

                    nextMEM = [int(elem) for elem in nextMEM[1:-1].split(",")]
                    nextIndex = self.BV.rank(nextMEM[0] - 1)
                    nextLabel = self.getNodeLabel(nextIndex)

                    if currIndex != nextIndex: #MEMs are aligned to different exons
                        #The edge is incremented by 1
                        self.incrementEdge(currIndex, nextIndex)

                        #Overlap on P
                        overlap = nextMEM[1]-currMEM[1]-currMEM[2]
                        #Gap on first exon end
                        currOffset = self.BV.select(currIndex + 1) - (currMEM[0] + currMEM[2])
                        #Gap on second exon start
                        nextOffset = nextMEM[0] - self.BV.select(nextIndex)-1

                        #Now we check for internal alt splice sites
                        if overlap <= 0:
                            #  --- TODO -----------
                            # - Canonical Patterns -
                            #  --------------------
                            overlap = abs(overlap)
                            #Values needed for IR
                            pos1 = 0
                            pos2 = 0
                            # First we check for alt splice site on "current_index" exon
                            if currOffset > 0:
                                currOffset += overlap
                                annPos = self.ExonPos[currLabel][1] #end of exon
                                newPos = annPos-currOffset
                                if not self.isExonBound(newPos):
                                    if self.GeneStrand: #Strand + -> AS 5'
                                        self.IntA5[newPos] = self.IntA5[newPos]+1 if newPos in self.IntA5 else 1
                                    else: #Strand + -> AS 3'
                                        self.IntA3[newPos] = self.IntA3[newPos]+1 if newPos in self.IntA3 else 1
                                pos1 = newPos
                            else:
                                pos1 = self.ExonPos[currLabel][1]
                            # Now we check for alt splice site on "next_index" exon
                            if nextOffset > 0:
                                nextOffset += overlap
                                annPos = self.ExonPos[nextLabel][0] #start of exon
                                newPos = annPos+nextOffset
                                if not self.isExonBound(newPos):
                                    if self.GeneStrand: #Strand + -> AS 3'
                                        self.IntA3[newPos] = self.IntA3[newPos]+1 if newPos in self.IntA3 else 1
                                    else:
                                        self.IntA5[newPos] = self.IntA5[newPos]+1 if newPos in self.IntA5 else 1
                                pos2 = newPos
                            else:
                                pos2 = self.ExonPos[nextLabel][0]
                            # Now we check for intron retention
                            if self.isInsideExon(pos1, pos2):
                                key = (pos1, pos2)
                                self.introns[key] = self.introns[key]+1 if key in self.introns else 1

                        #Insertions -> External alt splice sites
                        else:
                            #print(currOffset, nextOffset)
                            if currOffset == 0 and nextOffset == 0:
                                #There is an insertions
                                IntronStart = self.ExonPos[self.ExonNames[currIndex-1]][1] + 1
                                IntronEnd = self.ExonPos[self.ExonNames[nextIndex-1]][0]
                                intron = self.ref[IntronStart-1:IntronEnd-1] #-1 since strings are indexed starting from 0, gtf from 1
                                readGap = read[currMEM[1]+currMEM[2]-1:nextMEM[1]-1]
                                maxErr = round(len(read)*self.errRate/100)
                                err1 = editDistance(readGap, intron[:len(readGap)])
                                err2 = editDistance(readGap, intron[len(intron)-len(readGap):])
                                if err1<=err2:
                                    #We check the start of the intron
                                    if err1 + err <= maxErr:
                                        annPos = self.ExonPos[currLabel][1] #end of curr exon (start of intron - 1)
                                        newPos = annPos + len(readGap)
                                        if not self.isInsideExon(annPos+1, newPos) and not self.isExonBound(newPos):
                                            if self.GeneStrand: #Strand + -> AS 5'
                                                self.ExtA5[(annPos, newPos)] = self.ExtA5[(annPos, newPos)]+1 if (annPos, newPos) in self.ExtA5 else 1
                                            else: #Strand + -> AS 3'
                                                self.ExtA3[(annPos, newPos)] = self.ExtA3[(annPos, newPos)]+1 if (annPos, newPos) in self.ExtA3 else 1
                                else:
                                    #We check the end of the intron
                                    if err2 + err <= maxErr:
                                        annPos = self.ExonPos[nextLabel][0] #start of next exon (end of intron + 1)
                                        newPos = annPos - len(readGap)
                                        if not self.isInsideExon(newPos, annPos-1) and not self.isExonBound(newPos):
                                            if self.GeneStrand: #Strand + -> AS 3'
                                                self.ExtA3[(annPos, newPos)] = self.ExtA3[(annPos, newPos)]+1 if (annPos, newPos) in self.ExtA3 else 1
                                            else: #Strand + -> AS 5'
                                                self.ExtA5[(annPos, newPos)] = self.ExtA5[(annPos, newPos)]+1 if (annPos, newPos) in self.ExtA5 else 1
                    else: #mems are aligned to same exon
                        Poverlap = nextMEM[1]-currMEM[1]-currMEM[2]
                        Toverlap = nextMEM[0]-currMEM[0]-currMEM[2]
                        #Intron Retention
                        if Poverlap <= 0 and Toverlap > 0:
                            gap = Toverlap+abs(Poverlap)-1
                            if gap > 0:
                                label = self.getNodeLabel(currIndex)
                                pos1 = self.ExonPos[label][0]+currMEM[0]+currMEM[2]-self.BV.select(currIndex)+Poverlap
                                pos2 = pos1 + gap
                                key = (pos1, pos2)
                                self.introns[key] = self.introns[key]+1 if key in self.introns else 1

    def printEv(self, evType, label1, label2, pos1, pos2, info, w):
        print("{} {} {} {} {} {} {}".format(evType, label1, label2, pos1, pos2, info, w))

    #Return True if pos is an exon boundary, False otherwise
    def isExonBound(self, pos):
        for (s,e) in self.ExonPos.values():
            if pos == s or pos == e:
                return True
        return False

    #Return True if [pos1..pos2] is inside at least one exon, False otherwise
    def isInsideExon(self, pos1, pos2):
        for (s,e) in self.ExonPos.values():
            if s <= pos1 < pos2 <= e:
                return True
        return False

    def isRealES(self, pos1, pos2):
        for (s,e) in self.ExonPos.values():
            if pos1 <= s < e <= pos2:
                return True
        return False

    # #Return an integer (0..100) representing the percentage of real intron in [pos1..pos2]
    # def intronRelativePerc(self,pos1,pos2):
    #     print("---")
    #     x = [1]*(pos2-pos1)
    #     print(x)
    #     for (s,e) in self.ExonPos.values():
    #         if s < pos1 < pos2 < e:
    #             print(0)
    #             #return 0
    #         elif pos1 <= s < e <= pos2:
    #             print(1)
    #             x[s-pos1:pos2-e+1] = [0]*(e-s)
    #         elif s < pos1 < e <= pos2:
    #             print(2)
    #             print(s, pos1, e, pos2)
    #             print(0, pos2-pos1, e-pos1)
    #             x[:e-pos1] = [0]*(e-pos1)
    #         elif pos1 <= s < pos1 < e:
    #             print(3)
    #             x[s-pos1:] = [0]*(pos2-s)
    #         print(x)
    #     print(sum(x)/len(x)*100)

    ## EVENTS
    ##########
    def extractEvents(self):
        self.extractES()
        self.extractIntAS()
        self.extractExtAS()
        self.extractIR()

    def extractES(self):
        combinedES = {}
        for (n1,n2),w in self.newEdges.items():
            label1 = self.getNodeLabel(n1)
            label2 = self.getNodeLabel(n2)
            pos1 = self.ExonPos[label1][1]
            pos2 = self.ExonPos[label2][0]
            if (pos1,pos2) not in combinedES:
                combinedES.update({(pos1,pos2):0})
            combinedES[(pos1,pos2)] = combinedES[(pos1,pos2)] + w
        for (pos1,pos2),w in combinedES.items():
            if w >= conf:
                if self.isRealES(pos1, pos2):
                    self.printEv("ES", "-", "-", pos1, pos2, "-", w)

    def extractExtAS(self):
        for (annPos,newPos),w in self.ExtA5.items():
            if w >= conf:
                self.printEv("A5", "-", "-", annPos, newPos, "-", w)
        for (annPos,newPos),w in self.ExtA3.items():
            if w >= conf:
                self.printEv("A3", "-", "-", annPos, newPos, "-", w)

    def extractIntAS(self):
        combinedA5 = set()
        for newPos,w in self.IntA5.items():
            if w >= conf:
                for label,(p1,p2) in self.ExonPos.items():
                    if p1<newPos<p2:
                        annPos = self.ExonPos[label][1] if self.GeneStrand else self.ExonPos[label][0]
                        if (annPos, newPos) not in combinedA5:
                            combinedA5.add((annPos, newPos))
                            self.printEv("A5", "-", "-", annPos, newPos, "-", w)

        combinedA3 = set()
        for newPos,w in self.IntA3.items():
            if w >= conf:
                for label,(p1,p2) in self.ExonPos.items():
                    if p1<newPos<p2:
                        annPos = self.ExonPos[label][0] if self.GeneStrand else self.ExonPos[label][1]
                        if (annPos, newPos) not in combinedA3:
                            combinedA3.add((annPos, newPos))
                            self.printEv("A3", "-", "-", annPos, newPos, "-", w)

    def extractIR(self):
        #{(exID, posT, len) : weight}
        for (pos1,pos2),w in self.introns.items():
            if w >= conf:
                self.printEv("IR", "-", "-", pos1, pos2, "-", w)

    #Nodes Methods
    def addNode(self, label, w=0):
        if label != "" and label not in self.labels:
            self.labels.append(label)
            self.nodes.update({len(self.labels):[label, w]})

    def incrementNode(self, n1, w=1):
        if n1 in self.nodes:
            self.nodes[n1][1] += w

    def getNodeLabel(self, n):
        return self.nodes[n][0]

    def getNodeW(self, n):
        return self.nodes[n][1]

    #Edges Methods
    def addEdge(self, n1, n2, t, w=0):
        if t == 'e':
            if (n1,n2) not in self.edges:
                self.edges.update({(n1,n2):w})
        elif t == 'n':
            if (n1,n2) not in self.newEdges:
                self.newEdges.update({(n1,n2):w})

    def incrementEdge(self, n1, n2, w=1):
        if (n1,n2) in self.edges:
            self.edges[(n1,n2)] += w
        elif (n1,n2) in self.newEdges:
            self.newEdges[(n1,n2)] += w

    def isEdge(self, n1, n2):
        if (n1,n2) in self.edges:
            return True
        elif (n1,n2) in self.newEdges:
            return True
        else:
            return False

    def isNewEdge(self, n1, n2):
        if (n1,n2) in self.newEdges:
            return True
        else:
            return False

    ## Utils
    #########
    def clean(self, nTresh = 0, eTresh = 0):
        edgesToRemove = []
        for label in self.labels:
            if label == "":
                continue
            index = self.labels.index(label)+1
            if self.getNodeW(index) <= nTresh:
                self.labels[index-1] = ""
                self.nodes.pop(index)
                for (n1,n2),cov in self.edges.items():
                    if index == n1 or index == n2:
                        edgesToRemove.append((n1,n2))
                for (n1,n2),cov in self.newEdges.items():
                    if index == n1 or index == n2:
                        edgesToRemove.append((n1,n2))

        for (n1,n2),w in self.edges.items():
            if w <= eTresh:
                edgesToRemove.append((n1,n2))
        for (n1,n2),w in self.newEdges.items():
            if w <= eTresh:
                edgesToRemove.append((n1,n2))
        for edge in edgesToRemove:
            try:
                self.edges.pop(edge)
                continue
            except KeyError:
                pass
            try:
                self.newEdges.pop(edge)
                continue
            except KeyError:
                pass

    def getAdjMatrix(self):
        A = [[0 for x in range(0, len(self.labels))] for y in range(0, len(self.labels))]
        for (n1,n2) in self.edges:
            if n1 != n2:
                A[n1-1][n2-1] = 1
        for (n1,n2) in self.newEdges:
            A[n1-1][n2-1] = 1
        return A

    def print(self):
        print("### NODES ###")
        for label in self.labels:
            if label != "":
                index = self.labels.index(label)+1
                print("{} {}: {}".format(index, label, self.getNodeW(index)))
        print("### EDGES ###")
        for (n1,n2),c in self.edges.items():
            print("{}->{}: {}".format(n1, n2, c))
        print("### NEW EDGES ###")
        for (n1,n2),c in self.newEdges.items():
           print("{}->{}: {}".format(n1, n2, c))

    def save(self, name):
        print("Saving...")
        g = Digraph('G', filename="{}.gv".format(name))#os.path.join(OUT, "graph.gv"))
        g.attr('node', shape='circle')
        for label in self.labels:
            if label != "":
                index = self.labels.index(label)
                n_label = "{} ({})".format(label, self.getNodeW(index+1))
                g.node(n_label)
        for (n1,n2),cov in self.edges.items():
            n1_label = "{} ({})".format(self.labels[n1-1], self.getNodeW(n1))
            n2_label = "{} ({})".format(self.labels[n2-1], self.getNodeW(n2))
            g.edge(n1_label, n2_label, label=str(cov), color = "black")
        for (n1,n2),cov in self.newEdges.items():
            n1_label = "{} ({})".format(self.labels[n1-1], self.getNodeW(n1))
            n2_label = "{} ({})".format(self.labels[n2-1], self.getNodeW(n2))
            g.edge(n1_label, n2_label, label=str(cov), color = "red")
        g.render()

def main():
    global conf
    ref = sys.argv[1]      #reference
    infoPath = sys.argv[2] #index
    memsPath = sys.argv[3] #MEMs
    errRate = int(sys.argv[4])
    conf = int(sys.argv[5])
    onlyPrimFlag = bool(sys.argv[6])

    G = SplicingGraph(ref, infoPath, errRate)
    G.augments(memsPath, onlyPrimFlag)
    G.extractEvents()

if __name__ == '__main__':
    main()
