import os, sys
from BitVector import BV

def main(out_path):
    genomic = "/data/Toxoplasma/genomic.fasta"
    annotation = "/data/Toxoplasma/Genes_200/"
    with open("./tmp/T.fa") as t:
        text = t.read().split("\n")[1]
        bv = BV(text)
    edges = {}
    with open("./tmp/real_edges") as a:
        es = a.read().split("\n")
        for e in es:
            if e != "":
                edges.update({e:0})
    with open("./tmp/added_edges") as a:
        es = a.read().split("\n")
        for e in es:
            if e != "":
                edges.update({e:0})

    f = open("{}_Edges".format(out_path), "w")
    #if len(edges) > 0:
    with open(out_path) as o:
        for line in o:
            curr_edges = []
            if line != "":
                strings = line.split(" ")
                p_id, mems_string, err = strings[0], strings[1:-1], int(strings[-1])
                mems = []
                for mem in mems_string:
                    if mem != "" and mem != "\n":
                        m = mem[1:-1].split(",")
                        mems.append({"t":int(m[0]), "p":int(m[1]), "l":int(m[2])})
                if len(mems) > 1:
                    i = 0
                    while i<len(mems)-1:
                        e1 = str(bv.rank(mems[i]["t"] - 1))
                        e2 = str(bv.rank(mems[i+1]["t"] - 1))
                        if e1 != e2:
                            curr_edges.append("{},{}".format(e1, e2))
                        i+=1
                    if curr_edges == []:
                        curr_edges.append(str(bv.rank(mems[1]["t"] - 1)))
                    f.write("{}\t{}\t{}\n".format(p_id, err, " ".join(curr_edges)))
    f.close()

if __name__ == '__main__':
    main(sys.argv[1])
