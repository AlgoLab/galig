import sys, os

from utils import Annotation

from graphviz import Digraph
import numpy as np

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

def main():
    genomic_path = sys.argv[1]
    annotation_path = sys.argv[2]
    out_path = sys.argv[3]
    OUT = os.path.dirname(out_path)

    A = Annotation(genomic_path, annotation_path)

    nodes = {}
    edges = {}
    comps = {}
    with open(out_path, 'r') as out:
        for line in out:
            align = line[:-2].split(" ")
            '''
             0: strand
             1: ID
             2: errors
             3+: mems
            '''
            last_ex = -1
            last_mem = (-1,-1,-1)
            for mem in align[3:]:
                m = mem[1:-1].split(',')
                (t, p, l) = (int(m[0]), int(m[1]), int(m[2]))
                curr_ex = A.rank(t-1)
                if curr_ex not in nodes:
                    nodes.update({curr_ex:0})
                    comps_list = []
                    for i in range(0,A.getExonLength(curr_ex)):
                        comps_list.append(0)
                    comps.update({curr_ex:comps_list})
                nodes[curr_ex] = nodes[curr_ex]+1

                if last_ex != -1:
                    if last_ex != curr_ex:
                        edge = (last_ex, curr_ex)
                        if edge not in edges:
                            edges.update({edge:0})
                        edges[edge] = edges[edge]+1

                        (last_t, last_p, last_l) = (last_mem[0], last_mem[1], last_mem[2])
                        last_comp = last_t+last_l-1 - A.select(A.rank(last_t-1)) - 1
                        curr_comp = (t-1) - A.select(A.rank(t-1))
                        comps[last_ex][last_comp] += 1
                        comps[curr_ex][curr_comp] += 1
                last_ex = curr_ex
                last_mem = (t,p,l)

    with open(os.path.join(OUT, "graph.log"), 'w') as out:
        nodes_line = ""
        for k,v in nodes.items():
            nodes_line += "{}:{} ".format(k,v)
        out.write(nodes_line[:-1])
        out.write("\n")
        edges_line = ""
        for k,v in edges.items():
            edges_line += "{}-{}:{} ".format(k[0],k[1],v)
        out.write(edges_line[:-1])
        out.write("\n")

    with open(os.path.join(OUT, "comps.log"), 'w') as out:
        for k,v in comps.items():
            name = A.getExonName(k)
            out.write("{} {} ".format(k, name))
            hist_line = ""
            for i in v:
                hist_line += str(i) + ","
            out.write(hist_line[:-1])
            out.write("\n")

    for k,comp in comps.items():
        if sum(comp) != 0:
            name = A.getExonName(k)
            x = range(len(comp))
            plt.clf()
            plt.bar(x, comp, 1, color="blue")
            plt.title(name)
            plt.xlabel("Position")
            plt.ylabel("No Alignments")
            plt.xlim(0,len(comp))
            plt.ylim(0,max(comp))
            plt.savefig(os.path.join(OUT, "{}.png".format(name)))

    g = Digraph('G', filename=os.path.join(OUT, "graph.gv"))
    g.attr('node', shape='circle')
    for node,count in nodes.items():
        g.node(A.getExonName(node) + ',' + str(count))
    for edge in edges:
        g.edge(str(A.getExonName(edge[0])) + ',' + str(nodes[edge[0]]), str(A.getExonName(edge[1])) + ',' + str(nodes[edge[1]]), label=str(edges[edge]))
    g.render()

if __name__ == "__main__":
    main()
