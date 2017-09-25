import sys, os

from Bio import SeqIO
from graphviz import Digraph

from BitVector import BitVector
from utils import *

#Confidence values for each event
conf = 1

class SplicingGraph:
    def __init__(self, ref, reads, infoPath):
        infofile = open(infoPath)
        lines = infofile.readlines()
        chromo = lines[0].strip("\n").split(" ")[0]
        self.ref = list(SeqIO.parse(ref, "fasta"))[0].seq
        self.reads = SeqIO.index(reads, "fasta")
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

    #Augmenting
    def augments(self, memsPath):
        '''
        - comps:   {exID: {offset : weight}}
        - ins:     {(exID1,exID2,size) : weight}
        - SNVs:    {(exID, posT, posP) : weight} (posT: distance from label beginning)
        - introns: {(exID, posT, len) : weight} (posT: distance from label beginning)
        '''
        self.comp3 = {}
        self.comp5 = {}
        self.insertions = {}
        self.SNVs = {}
        self.introns = {}
        with open(memsPath, 'r') as out:
            for line in out:
                # 0: strand
                # 1: ID
                # 2: errors
                # 3+: mems
                align = line.strip("\n").strip(" ").split(" ")
                readID = align[1]
                usedExons = set([self.BV.rank(int(mem[1:-1].split(',')[0]) - 1) for mem in align[3:]])
                for e in usedExons:
                    self.incrementNode(e)
                if len(align[3:]) <= 1:
                    continue
                for currMEM, nextMEM in pairwise(align[3:]):
                    # Remove ( and ) from mem and cast to int
                    currMEM = [int(elem) for elem in currMEM[1:-1].split(",")]
                    nextMEM = [int(elem) for elem in nextMEM[1:-1].split(",")]
                    currIndex = self.BV.rank(currMEM[0] - 1)
                    nextIndex = self.BV.rank(nextMEM[0] - 1)
                    if currIndex != nextIndex: #mems are aligned to different exons
                        #edge is incremented
                        self.incrementEdge(currIndex, nextIndex)

                        overlap = nextMEM[1]-currMEM[1]-currMEM[2]
                        currOffset = self.BV.select(currIndex + 1) - (currMEM[0] + currMEM[2])
                        nextOffset = nextMEM[0] - self.BV.select(nextIndex)-1

                        #Competings
                        if overlap <= 0:
                            #  --- TODO -----------
                            # - Canonical Patterns -
                            #  --------------------
                            overlap = abs(overlap)
                            # First we check for 5' competings in exon "current_index"
                            if currOffset > 0:
                                currOffset += overlap
                                if currIndex not in self.comp5:
                                    self.comp5[currIndex] = {}
                                self.comp5[currIndex][currOffset] = self.comp5[currIndex][currOffset]+1 if currOffset in self.comp5[currIndex] else 1
                            # Now we check for 3' competings in exon "next_index"
                            if nextOffset > 0:
                                nextOffset += overlap
                                if nextIndex not in self.comp3:
                                    self.comp3[nextIndex] = {}
                                self.comp3[nextIndex][nextOffset] = self.comp3[nextIndex][nextOffset]+1 if nextOffset in self.comp3[nextIndex] else 1
                        #Insertions and SNVs
                        else:
                            if currOffset == 0 and nextOffset == 0:
                                #Insertions
                                IntronStart = self.ExonPos[self.ExonNames[currIndex-1]][1] + 1
                                IntronEnd = self.ExonPos[self.ExonNames[nextIndex-1]][0]
                                intron = self.ref[IntronStart:IntronEnd]
                                read = self.reads[readID].seq
                                readGap = read[currMEM[1]+currMEM[2]:nextMEM[1]]
                                if readGap == intron[:len(readGap)]:
                                    if currIndex not in self.comp5:
                                        self.comp5[currIndex] = {}
                                    self.comp5[currIndex][-len(readGap)] = self.comp5[currIndex][-len(readGap)]+1 if -len(readGap) in self.comp5[currIndex] else 1
                                elif readGap == intron[len(intron)-len(readGap):]:
                                    if nextIndex not in self.comp3:
                                        self.comp3[nextIndex] = {}
                                    self.comp3[nextIndex][-len(readGap)] = self.comp3[nextIndex][-len(readGap)]+1 if -len(readGap) in self.comp3[nextIndex] else 1
                                else:
                                    key = (currIndex, nextIndex, overlap)
                                    self.insertions[key] = self.insertions[key]+1 if key in self.insertions else 1
                            else:
                                #SNV
                                if overlap == currOffset + nextOffset:
                                    key = (currIndex, currMEM[0]+currMEM[2]-self.BV.select(currIndex), currMEM[1]+currMEM[2])
                                    self.SNVs[key] = self.SNVs[key]+1 if key in self.SNVs else 1
                                    key = (nextIndex, nextMEM[0]-1-self.BV.select(nextIndex), nextMEM[1]-1)
                                    self.SNVs[key] = self.SNVs[key]+1 if key in self.SNVs else 1
                    else: #mems are aligned to same exon
                        Poverlap = nextMEM[1]-currMEM[1]-currMEM[2]
                        Toverlap = nextMEM[0]-currMEM[0]-currMEM[2]
                        #NewIntron
                        if Poverlap == 0 and Toverlap > 0:
                            key = (currIndex, currMEM[1]+currMEM[2]-self.BV.select(currIndex), Toverlap)
                            self.introns[key] = self.introns[key]+1 if key in self.introns else 1
                        #SNV
                        if Poverlap == Toverlap:
                            key1 = (currIndex, currMEM[0]+currMEM[2]-self.BV.select(currIndex), currMEM[1]+currMEM[2])
                            key2 = (currIndex, nextMEM[0]-1-self.BV.select(currIndex), nextMEM[1]-1)
                            self.SNVs[key1] = self.SNVs[key1]+1 if key1 in self.SNVs else 1
                            if key1 != key2:
                                self.SNVs[key2] = self.SNVs[key2]+1 if key2 in self.SNVs else 1

    def printEv(self, evType, label1, label2, pos1, pos2, info, w):
        print("{} {} {} {} {} {} {}".format(evType, label1, label2, pos1, pos2, info, w))

    ## EVENTS
    ##########
    def extractEvents(self):
        self.extractES()
        self.extractC(self.comp5, True)
        self.extractC(self.comp3, False)
        self.extractMEE()

    def extractES(self):
        for (n1,n2),w in self.newEdges.items():
            if w >= conf:
                label1 = self.getNodeLabel(n1)
                label2 = self.getNodeLabel(n2)
                self.printEv("ES", label1, label2, self.ExonPos[label1][1], self.ExonPos[label2][0], "-", w)

    def extractC(self,CompDict,tFlag):
        for exid,comp in CompDict.items():
            label = self.getNodeLabel(exid)
            for offset,w in comp.items():
                if w >= conf:
                    if self.GeneStrand:
                        t = "A5" if tFlag else "A3"
                        annPos = self.ExonPos[label][1] if tFlag else self.ExonPos[label][0]
                        newPos = annPos-offset if tFlag else annPos+offset
                    else:
                        t = "A3" if tFlag else "A5"
                        annPos = self.ExonPos[label][1] if tFlag else self.ExonPos[label][0]
                        newPos = annPos-offset if tFlag else annPos+offset
                    self.printEv(t, label, "-", annPos, newPos, "-", w)

    def extractMEE(self):
        self.clean(conf,conf)
        A = self.getAdjMatrix()
        alreadyChecked = []
        for node1 in [node for node in self.nodes]:
            for node2 in [node for node in self.nodes if node != node1]:
                if (node1,node2) in alreadyChecked or (node2,node1) in alreadyChecked:
                    continue
                alreadyChecked.append((node1,node2))
                if not(self.isEdge(node1,node2)):
                    if sum(getCommonAncestors(A, node1-1, node2-1))>0 and sum(getCommonDescendants(A, node1-1, node2-1))>0:
                        w = 2
                        label1 = self.getNodeLabel(node1)
                        label2 = self.getNodeLabel(node2)
                        if w >= conf:
                            print("MEE {} {} {} {} {} {} {}".format(label1, label2, self.ExonPos[label1][0], self.ExonPos[label1][1], self.ExonPos[label2][0], self.ExonPos[label2][1], w))

    ## INFO
    ########
    def extractInfo(self):
        self.extractIns()
        self.extractSNV()
        self.extractIntrons()

    def extractIns(self):
        #{(exID1,exID2,size) : weight}
        for (n1,n2,l),w in self.insertions.items():
            if w >= conf:
                label1 = self.getNodeLabel(n1)
                label2 = self.getNodeLabel(n2)
                self.printEv("NE", label1, label2, self.ExonPos[label1][1], self.ExonPos[label2][0], l, w)

    def extractSNV(self):
        #{(exID, posT, posP) : weight}
        for (n,posT,posP),w in self.SNVs.items():
            if w >= conf:
                label = self.getNodeLabel(n)
                self.printEv("SNV", label, "-", self.ExonPos[label][0]+posT, "-", posP, w)

    def extractIntrons(self):
        #{(exID, posT, len) : weight}
        for (n,posT,l),w in self.introns.items():
            label = self.getNodeLabel(n)
            self.printEv("IR", label, "-", self.ExonPos[label][0]+posT, self.ExonPos[label][0]+posT+l, "-", w)
    
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
    ref = sys.argv[1]      #reference
    reads = sys.argv[2]    #reads
    infoPath = sys.argv[3] #index
    memsPath = sys.argv[4] #MEMs
    G = SplicingGraph(ref, reads, infoPath)
    G.augments(memsPath)
    G.extractEvents()
    G.extractInfo()

if __name__ == '__main__':
    main()
