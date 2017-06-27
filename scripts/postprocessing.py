import sys, os

from BitVector import BitVector
from SplicingGraph import SplicingGraph

#Confidence values for each event
ES_conf = 30
C_conf = 10
MEE_conf = 50

def main():
    info_path = sys.argv[1]
    out_path = sys.argv[2]

    with open(info_path) as info:
        i = 0
        for line in info:
            if i == 1:
                text = line
            if i == 3:
                adj_matrix = []
                for row in line[:-2].split(";"):
                    int_row = []
                    for elem in row.split(" "):
                        if elem != "":
                            int_row.append(int(elem))
                    if int_row != []:
                        adj_matrix.append(int_row)
            if i == 5:
                names = line[:-2].split(" ")
            i+=1
    BV = BitVector(text)

    G = SplicingGraph()
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
    G.save("1")

    competings = {}
    with open(out_path, 'r') as out:
        for line in out:
            '''
             0: strand
             1: ID
             2: errors
             3+: mems
            '''
            align = line[:-2].split(" ")
            used_exons = []
            last_exid = -1
            last_mem = (-1,-1,-1)
            for mem in align[3:]:
                m = mem[1:-1].split(',')
                (t, p, l) = (int(m[0]), int(m[1]), int(m[2]))
                curr_exid = BV.rank(t-1)
                if curr_exid not in used_exons:
                    G.incrementNode(curr_exid)
                    used_exons.append(curr_exid)
                if last_exid != -1:
                    if last_exid != curr_exid:
                        G.incrementEdge(last_exid, curr_exid)
                        (last_t, last_p, last_l) = (last_mem[0], last_mem[1], last_mem[2])
                        last_comp = BV.select(BV.rank(last_t-1)+1) - (last_t+last_l)
                        curr_comp = t-1 - BV.select(BV.rank(t-1))
                        if last_comp > 0:
                            if last_exid not in competings:
                                competings.update({last_exid:{}})
                            if (True, curr_exid, last_comp) not in competings[last_exid]:
                                competings[last_exid].update({(True, curr_exid, last_comp):0})
                            competings[last_exid][(True, curr_exid, last_comp)] += 1
                        if curr_comp > 0:
                            if curr_exid not in competings:
                                competings.update({curr_exid:{}})
                            if (False, last_exid, curr_comp) not in competings[curr_exid]:
                                competings[curr_exid].update({(False, last_exid, curr_comp):0})
                            competings[curr_exid][(False, last_exid, curr_comp)] += 1
                last_exid = curr_exid
                last_mem = (t,p,l)
    G.save("2")
    confirmed_competings = []
    for exid1, comps in competings.items():
        for (type,exid2,length),cov in comps.items():
            if cov > C_conf:
                if type:
                    #5'
                    G.decrementNode(exid1, w=cov)
                    G.decrementNode(exid2, w=cov)
                    G.decrementEdge(exid1, exid2, w=cov)
                    label = "{}_{}_{}".format(G.getLabel(exid1), 0, length)
                    G.addNode(label, w=cov, comp3=0, comp5=length)
                    new_exid = G.getIndex(label)
                    G.addEdge(exid1, new_exid, 'f')
                    t = G.getType(exid1, exid2)
                    G.addEdge(new_exid, exid2, t, w=cov)
                    confirmed_competings.append([type, exid1, length, cov])
                else:
                    #3'
                    G.decrementNode(exid1, w=cov)
                    G.decrementNode(exid2, w=cov)
                    G.decrementEdge(exid2, exid1, w=cov)
                    label = "{}_{}_{}".format(G.getLabel(exid1), length, 0)
                    G.addNode(label, w=cov, comp3=length, comp5=0)
                    new_exid = G.getIndex(label)
                    G.addEdge(exid1, new_exid, 'f')
                    t = G.getType(exid2, exid1)
                    G.addEdge(exid2, new_exid, t, w=cov)
                    confirmed_competings.append([type, exid1, length, cov])
    G.clean(0, 0)
    G.save("3")

    #Exons Skipping
    for (n1,n2),w in G.new_edges.items():
        if w > ES_conf:
            print("ES", G.getLabel(n1), G.getLabel(n2), w, sep=",")

    #Competing
    for [type, exid1, length, cov] in confirmed_competings:    
        if type:
            t = "5'"
            if G.isSource(exid1):
                t+='S'
            if G.isSink(exid1):
                t+='E'
            print(t, G.getLabel(exid1), cov, length, sep=",")
        else:
            t = "3'"
            if G.isSource(exid1):
                t+='S'
            if G.isSink(exid1):
                t+='E'
            print(t, G.getLabel(exid1), cov, length, sep=",")

    #Mutually Exclusive Exons
    G.buildOutLists()
    G.print()
    found_MEEs = []
    for node in G.nodes:
        if node in G.outLists and node not in G.newNodes:
            fake_sons = G.getFakeSons(node)
            sons_of_fake_sons = {}
            for fake_son in fake_sons:
                sons_of_fake_sons.update(G.outLists[fake_son])
            all_sons = {}
            all_sons.update(G.outLists[node])
            all_sons.update(sons_of_fake_sons)
            for out1,w1 in all_sons.items():
                for out2,w2 in G.outLists[node].items():
                    if out1 != out2:
                        if out1 in G.outLists and out2 in G.outLists and out1 not in G.outLists[out2] and out2 not in G.outLists[out1] and ((out1, out2) not in found_MEEs and (out2, out1) not in found_MEEs):
                            print(node, out1, out2)
                            try:
                                fake_sons_of_sons_1 = {}
                                for son in G.outLists[out1]:
                                    fake_sons_of_sons_1.update(G.getFakeSons(son))
                                fake_sons_of_sons_2 = {}
                                for son in G.outLists[out2]:
                                    fake_sons_of_sons_2.update(G.getFakeSons(son))
                                print(set(G.outLists[out1]), set(fake_sons_of_sons_1), set(G.outLists[out2]), set(fake_sons_of_sons_2))
                                common_sons = set.intersection(set.union(set(G.outLists[out1]), set(fake_sons_of_sons_1)), set.union(set(G.outLists[out2]), set(fake_sons_of_sons_2)))
                                #common_sons = set.intersection(set(G.outLists[out1]), set(G.outLists[out2]))
                            except KeyError:
                                continue
                            print(common_sons)
                            if len(common_sons) > 0:
                                w = 0
                                for son in common_sons:
                                    print(node, out1, out2, son)
                                    w += w1+w2+G.getEdgeWeight(out1, son)+G.getEdgeWeight(out2, son)
                                print("MEE", G.getLabel(out1), G.getLabel(out2), w, sep=",")
                                found_MEEs.append((out1, out2))
    
if __name__ == '__main__':
    main()
