import sys, os

from BitVector import BitVector
from SplicingGraph import SplicingGraph
from utils import *

#Confidence values for each event
ES_conf = 30
C_conf = 10
MEE_conf = 50

def create_adj_matrix_from_line(line):
    adj_matrix = [[int(elem) for elem in row.split()] for row in line[:-2].split(";")]
    return adj_matrix

def readInfoPath(info_path):
    infofile = open(info_path)
    lines = infofile.readlines()
    adj_matrix = create_adj_matrix_from_line(lines[3].strip("\n"))
    names = lines[5].strip("\n").split()
    text = lines[1]
    BV = BitVector(text)
    return names, BV, adj_matrix

def initializeSG(names, adj_matrix):
    for name in names:
        G.addNode(name)
    r = 0
    for row in adj_matrix:
        c = 0
        for elem in row:
            if elem == 1:
                G.addEdge(r,c,'e')
            if elem == 2:
                G.addEdge(r,c,'n')
            c+=1
        r+=1

def extractInfoFromOut(out_path):
    competings = {}
    with open(out_path, 'r') as out:
        for line in out:
            # 0: strand
            # 1: ID
            # 2: errors
            # 3+: mems
            align = line.strip("\n").strip(" ").split(" ")
            used_exons = set([BV.rank(int(mem[1:-1].split(',')[0]) - 1) for mem in align[3:]])
            for e in used_exons:
                G.incrementNode(e)
            if len(align[3:]) <= 1:
                continue
            for current_mem, next_mem in pairwise(align[3:]):
                # Remove ( and ) from mem
                current_mem = [int(elem) for elem in current_mem[1:-1].split(",")]
                next_mem = [int(elem) for elem in next_mem[1:-1].split(",")]
                current_index = BV.rank(current_mem[0] - 1)
                next_index = BV.rank(next_mem[0] - 1)
                if current_index is not next_index:
                    # mems are aligned to different exons
                    G.incrementEdge(current_index, next_index)
                    # Check if there is some competing
                    # First we check for 5' competings in exon "current_index"
                    current_offset = BV.select(current_index + 1) - (current_mem[0] + current_mem[2])
                    if current_offset > 0:
                        if current_index not in competings:
                            competings[current_index] = {}
                        current_key = (True, next_index, current_offset)
                        competings[current_index][current_key] = competings[current_index][current_key] + 1 if current_key in competings[current_index] else 1
                    # Now we check for 3' competings in exon "next_index"
                    next_offset = next_mem[0]-1 - BV.select(next_index)
                    if next_offset > 0:
                        if next_index not in competings:
                            competings[next_index] = {}
                        next_key = (False, current_index, next_offset)
                        competings[next_index][next_key] = competings[next_index][next_key] + 1 if next_key in competings[next_index] else 1
    return competings

def checkCompetings(competings):
    confirmed_competings = []
    competings_for_position = {}
    for first_exon in competings:
        for competing_prop in competings[first_exon]:
            k = (first_exon, competing_prop[0], competing_prop[2])
            competings_for_position[k] = competings_for_position[k] + competings[first_exon][competing_prop] if k in competings_for_position else competings[first_exon][competing_prop]

    for first_exon in competings:
        for competing_prop in competings[first_exon]:
            side, second_exon, offset = competing_prop
            if competings_for_position[(first_exon, side, offset)] > C_conf:
                new_node_id = G.splitNode(first_exon, side, offset, competings[first_exon][competing_prop])
                if side is True:
                    source_exon_tmp, target_exon_tmp = first_exon, second_exon
                    G.addEdge(new_node_id, target_exon_tmp, G.getEdgeType(source_exon_tmp, target_exon_tmp), competings[first_exon][competing_prop])
                    ###
                    for parent in G.getParents(first_exon):
                        G.addEdge(parent, new_node_id, G.getEdgeType(parent, first_exon), G.getEdgeWeight(parent, first_exon))
                    ###
                else:
                    source_exon_tmp, target_exon_tmp = second_exon, first_exon
                    G.addEdge(source_exon_tmp, new_node_id, G.getEdgeType(source_exon_tmp, target_exon_tmp), competings[first_exon][competing_prop])
                    ###
                    if first_exon in G.outLists:
                        for son in G.outLists[first_exon]:
                            G.addEdge(new_node_id, son, G.getEdgeType(first_exon, son), G.getEdgeWeight(first_exon, son))
                    ###
                G.decrementEdge(source_exon_tmp, target_exon_tmp, competings[first_exon][competing_prop])

    # TODO: questo e' chiaramente un filtro            
    for c in competings_for_position:
        if competings_for_position[c] > C_conf:
            confirmed_competings.append([c[1], c[0], c[2], competings_for_position[c]])
    return confirmed_competings

def checkESEvents():
    for (n1,n2),w in G.new_edges.items():
        if w > ES_conf:
            print("ES", G.getLabel(n1), G.getLabel(n2), w, sep=",")

def checkCEvents(confirmed_competings):
    for [side, exid1, offset, cov] in confirmed_competings:    
        t = "5'" if side else "3'"
        if G.isSource(exid1):
            t+='S'
        elif G.isSink(exid1):
            t+='E'
        print(t, G.getLabel(exid1), cov, offset, sep=",")

def checkMEEEvents():
    A = G.getAdjMatrix()
    found_MEEs = []
    checked = []
    '''
    for node1 in G.nodes:
        for node2 in [node for node in G.nodes if node != node1]:
    '''
    for node1 in [node for node in G.nodes if node not in G.newNodes]:
        for node2 in [node for node in G.nodes if node not in G.newNodes and node != node1]:
            if (node1, node2) in checked or (node2, node1) in checked:
                continue
            checked.append((node1, node2))
            if not(G.isEdge(node1, node2)):
                if (node1, node2) in found_MEEs:
                    pass
                elif sum(getCommonAncestors(A, node1-1, node2-1))>0 and sum(getCommonDescendants(A, node1-1, node2-1))>0:
                    w=sum([c for _,c in G.getParents(node1).items()])
                    w+=sum([c for _,c in G.getParents(node2).items()])
                    if node1 in G.outLists:
                        w+=sum([c for _,c in G.outLists[node1].items()])
                    if node2 in G.outLists:
                        w+=sum([c for _,c in G.outLists[node2].items()])
                    print("MEE", G.getLabel(node1), G.getLabel(node2), w, sep=",")
                    found_MEEs.append((node1, node2))

def main():
    global G, BV, text

    info_path = sys.argv[1]
    out_path = sys.argv[2]

    names, text, BV, adj_matrix = readInfoPath(info_path)

    # G = SplicingGraph(info_path)
    # G = SplicingGraph(names, adj_matrix)
    G = SplicingGraph()
    initializeSG(names, adj_matrix)

    G.buildOutLists()
    competings = extractInfoFromOut(out_path)
    confirmed_competings = checkCompetings(competings)

    G.clean(0,0) # G.filter() def filter(self, x = 0, y = 0 ):
    G.buildOutLists()
    G.save("graph")
    checkESEvents()                     # G.checkESEvents()
    checkCEvents(confirmed_competings)  # G.check...
    checkMEEEvents()


if __name__ == '__main__':
    main()
