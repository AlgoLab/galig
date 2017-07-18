import sys, os

from BitVector import BitVector
from SplicingGraph import SplicingGraph
from utils import *

#Confidence values for each event
ES_conf = 1
C_conf = 5
MEE_conf = 1

def create_adj_matrix_from_line(line):
    adj_matrix = [[int(elem) for elem in row.split()] for row in line[:-2].split(";")]
    return adj_matrix

def readInfoPath(info_path):
    infofile = open(info_path)
    lines = infofile.readlines()
    if lines[0].strip("\n").split(" ")[2] == '+':
        strand = True
    else:
        strand = False
    adj_matrix = create_adj_matrix_from_line(lines[3].strip("\n"))
    names = lines[5].strip("\n").split()
    text = lines[1]
    BV = BitVector(text)
    pos = [(int(p[0]), int(p[1])) for p in [pos.split(",") for pos in lines[4].strip("\n").split()]]
    EPos = dict(zip(names, pos))
    return strand, names, text, BV, adj_matrix, EPos

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

                    overlap = abs(next_mem[1]-current_mem[1]-current_mem[2])

                    #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
                    #!# ADD CANONICAL PATTERNS
                    #!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#!#
                    # First we check for 5' competings in exon "current_index"
                    current_offset = BV.select(current_index + 1) - (current_mem[0] + current_mem[2])
                    if current_offset > 0:
                        current_offset += overlap
                        if current_index not in competings:
                            competings[current_index] = {}
                        current_key = (True, next_index, current_offset)
                        competings[current_index][current_key] = competings[current_index][current_key] + 1 if current_key in competings[current_index] else 1
                    # Now we check for 3' competings in exon "next_index"
                    next_offset = next_mem[0] - BV.select(next_index)-1
                    if next_offset > 0:
                        current_offset -= overlap
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

def checkESEvents(out):
    ESs = []
    for (n1,n2),w in G.new_edges.items():
        if w > ES_conf:
            n1_label = G.getNodeLabel(n1).split("-")[0]
            n2_label = G.getNodeLabel(n2).split("-")[0]
            if (n1_label,n2_label,w) not in ESs:
                ESs.append((n1_label, n2_label, w))

    for (n1_label,n2_label,w) in ESs:
        out.write("ES {} {} {} {} {}\n".format(n1_label, n2_label, EPos[n1_label][1], EPos[n2_label][0], w))

def checkCEvents(confirmed_competings, out, strand):
    for [side, exid1, offset, cov] in confirmed_competings:
        n_label = G.getNodeLabel(exid1)
        if strand:
            t = "5'" if side else "3'"
            ann_pos = EPos[n_label][1] if side else EPos[n_label][0]
            new_pos = ann_pos-offset if side else ann_pos+offset
        else:
            t = "3'" if side else "5'"
            ann_pos = EPos[n_label][1] if side else EPos[n_label][0]
            new_pos = ann_pos-offset if side else ann_pos+offset
        if G.isSource(exid1):
            t+='S'
        elif G.isSink(exid1):
            t+='E'

        out.write("{} {} {} {} {}\n".format(t, n_label, ann_pos, new_pos, cov))

def checkMEEEvents(out):
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
                    n1_label = G.getNodeLabel(node1)
                    n2_label = G.getNodeLabel(node2)
                    out.write("MEE {} {} {} {} {} {}\n".format(n1_label, n2_label, EPos[n1_label][0], EPos[n1_label][1], EPos[n2_label][0], EPos[n2_label][1], w))
                    found_MEEs.append((node1, node2))

def main():
    global G, BV, text, EPos

    info_path = sys.argv[1]
    out_path = sys.argv[2]
    out_file = sys.argv[3]

    strand, names, text, BV, adj_matrix, EPos = readInfoPath(info_path)

    # G = SplicingGraph(info_path)
    # G = SplicingGraph(names, adj_matrix)
    G = SplicingGraph()
    initializeSG(names, adj_matrix)

    G.buildOutLists()
    competings = extractInfoFromOut(out_path)
    confirmed_competings = checkCompetings(competings)

    G.clean(0,0) # G.filter() def filter(self, x = 0, y = 0 ):
    G.buildOutLists()
    G.save(out_file)
    
    out = open(out_file, 'w')
    checkESEvents(out)                     # G.checkESEvents()
    checkCEvents(confirmed_competings, out, strand)  # G.check...
    checkMEEEvents(out)
    out.close()

if __name__ == '__main__':
    main()
