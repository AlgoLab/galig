import os, sys
from BitVector import BV

def main(out_path):
    genomic = "/data/Toxoplasma/genomic.fasta"
    annotation = "/data/Toxoplasma/Genes_200/"
    with open("./tmp/T.fa") as t:
        text = t.read().split("\n")[1]
        bv = BV(text)
    with open("./tmp/added_edges") as a:
        edges = {}
        es = a.read().split("\n")
        for e in es:
            if e != "":
                edges.update({e:0})

    p_ids = {}
    p_edges = {}
    p_edges_use = {}
    p_edges_count = {}
    p_errs = {}
    if len(edges) > 0:
        with open(out_path) as o:
            for line in o:
                if line != "":
                    strings = line.split(" ")
                    p_id, mems_string, err = strings[0], strings[1:-1], int(strings[-1])
                    if p_id not in p_ids:
                        p_ids.update({p_id:1})
                        p_edges.update({p_id:[]})
                        p_edges_use.update({p_id:[]})
                        p_edges_count.update({p_id:0})
                        p_errs.update({p_id:err})
                    else:
                        p_ids[p_id] = p_ids[p_id] + 1
                    mems = []
                    for mem in mems_string:
                        if mem != "" and mem != "\n":
                            m = mem[1:-1].split(",")
                            mems.append({"t":int(m[0]), "p":int(m[1]), "l":int(m[2])})
                    if len(mems) > 1:
                        i = 0
                        use = 0
                        while i<len(mems)-1:
                            edge = str(bv.rank(mems[i]["t"] - 1)) + "," + str(bv.rank(mems[i+1]["t"] - 1))
                            if edge in edges:
                                edges[edge] = edges[edge] + 1
                                p_edges[p_id].append(edge)
                                p_edges_count[p_id] = p_edges_count[p_id] + 1
                                use+=1
                            i+=1
                        if use != 0:
                            p_edges_use[p_id].append(str(use))

    f = open("{}_NewEdges".format(out_path), "w")
    for p_id in p_ids:
        if p_edges[p_id] != []:
            f.write("{}\t{}\t{}\n".format(p_id, p_errs[p_id], " ".join(p_edges[p_id])))
    f.close()

if __name__ == '__main__':
    main(sys.argv[1])
