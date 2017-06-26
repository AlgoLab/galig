import sys, os

from BitVector import BitVector
from SplicingGraph import SplicingGraph

#Confidence values for each event
ES_conf = 50
C_conf = 2000
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
    print(competings)
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
                    G.addEdge(new_exid, exid2, 'e', w=cov)
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
                    G.addEdge(exid2, new_exid, 'e', w=cov)
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
            print("5'", G.getLabel(exid1), cov, length, sep=",")
        else:
            print("3'", G.getLabel(exid1), cov, length, sep=",")
    
if __name__ == '__main__':
    main()

'''

import numpy as np

import matplotlib
matplotlib.use('Agg')

def main():
    info_path = sys.argv[1]
    

    ES_COV = 1
    
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
                exs_name = line[:-2].split(" ")
            i+=1
    BV = BitVector(text)
    
    OUT = os.path.dirname(out_path) + "/sgal/"
    try:
        os.mkdir(OUT)
    except FileExistsError:
        pass

    nodes = {}
    edges = {}
    out_edges = {}
    comps_3 = {}
    comps_5 = {}
    with open(out_path, 'r') as out:
        for line in out:
            align = line[:-2].split(" ")
             #0: strand
             #1: ID
             #2: errors
             #3+: mems
            last_ex = -1
            last_mem = (-1,-1,-1)
            for mem in align[3:]:
                m = mem[1:-1].split(',')
                (t, p, l) = (int(m[0]), int(m[1]), int(m[2]))
                curr_ex = BV.rank(t-1)
                if curr_ex not in nodes:
                    nodes.update({curr_ex:0})
                    comps_list_5 = []
                    comps_list_3 = []
                    ex_len = BV.select(curr_ex+1) - BV.select(curr_ex) + 1 - 2
                    for i in range(0,ex_len):
                        comps_list_3.append(0)
                        comps_list_5.append(0)
                    comps_3.update({curr_ex:comps_list_3})
                    comps_5.update({curr_ex:comps_list_5})
                nodes[curr_ex] = nodes[curr_ex]+1

                if last_ex != -1:
                    if last_ex != curr_ex:
                        edge = (last_ex, curr_ex)
                        if edge not in edges:
                            edges.update({edge:0})
                        edges[edge] = edges[edge]+1

                        if last_ex not in out_edges:
                            out_edges.update({last_ex:[]})
                        if curr_ex not in out_edges[last_ex]: 
                            out_edges[last_ex].append(curr_ex)

                        (last_t, last_p, last_l) = (last_mem[0], last_mem[1], last_mem[2])
                        last_comp = last_t+last_l-1 - BV.select(BV.rank(last_t-1)) - 1
                        curr_comp = (t-1) - BV.select(BV.rank(t-1))
                        comps_5[last_ex][last_comp] += 1
                        comps_3[curr_ex][curr_comp] += 1
                last_ex = curr_ex
                last_mem = (t,p,l)

    #3' alternative
    for ex_id,coverage in comps_3.items():
        i = 1
        for pos in coverage[1:-1]:
            if pos != 0:
                print("3'", exs_name[ex_id-1], i, pos, sep=",")
            i+=1

    #5' alternative
    for ex_id,coverage in comps_5.items():
        i = 1
        for pos in coverage[1:-1]:
            if pos != 0:
                print("5'", exs_name[ex_id-1], i, pos, sep=",")
            i+=1

    #Exon Skipping
    for ids,coverage in edges.items():
        if adj_matrix[ids[0]][ids[1]] == 2:
            if coverage > ES_COV:
                print("ES", exs_name[ids[0]-1], exs_name[ids[1]-1], coverage, sep=",")

    #Mutually Exclusive Exons
    found_MEEs = []
    for node in nodes:
        if node in out_edges:
            for out1 in out_edges[node]: 
                for out2 in out_edges[node]:
                    if out1 != out2:
                        if out1 not in out_edges[out2] and out2 not in out_edges[out1] and ((out1, out2) not in found_MEEs and (out2, out1) not in found_MEEs):
                            try:
                                common_sons = set.intersection(set(out_edges[out1]), set(out_edges[out2]))
                            except KeyError:
                                continue
                            if len(common_sons) > 0:
                                w = 0
                                for son in common_sons:
                                    w += edges[(node, out1)]+edges[(node, out2)]+edges[(out1, son)]+edges[(out2, son)]
                                print("MEE", exs_name[out1-1], exs_name[out2-1], w, sep=",")
                                found_MEEs.append((out1, out2))
        
    with open(os.path.join(OUT, "graph.log"), 'w') as out:
        nodes_line = ""
        for k,v in nodes.items():
            nodes_line += "{}:{}:{} ".format(k,exs_name[k-1],v)
        out.write(nodes_line[:-1])
        out.write("\n")
        edges_line = ""
        for k,v in edges.items():
            edges_line += "{}-{}:{} ".format(k[0],k[1],v)
        out.write(edges_line[:-1])
        out.write("\n")

    # with open(os.path.join(OUT, "comps.log"), 'w') as out:
    #     for k,v in comps.items():
    #         name = exs_name[k-1]
    #         out.write("{} {} ".format(k, name))
    #         hist_line = ""
    #         for i in v:
    #             hist_line += str(i) + ","
    #         out.write(hist_line[:-1])
    #         out.write("\n")

    # for k,comp in compsnn.items():
    #     if sum(comp) != 0:
    #         name = A.getExonName(k)
    #         x = range(len(comp))
    #         plt.clf()
    #         plt.bar(x, comp, 1, color="blue")
    #         plt.title(name)
    #         plt.xlabel("Position")
    #         plt.ylabel("No Alignments")
    #         plt.xlim(0,len(comp))
    #         plt.ylim(0,max(comp))
    #         plt.savefig(os.path.join(OUT, "{}.png".format(name)))

    g = Digraph('G', filename=os.path.join(OUT, "graph.gv"))
    g.attr('node', shape='circle')
    for node,count in nodes.items():
        g.node(exs_name[node-1] + ',' + str(count))
    for edge in edges:
        g.edge(str(exs_name[edge[0]-1]) + ',' + str(nodes[edge[0]]), str(exs_name[edge[1]-1]) + ',' + str(nodes[edge[1]]), label=str(edges[edge]))
    g.render()

if __name__ == "__main__":
    main()
'''
