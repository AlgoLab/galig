import sys

from BitVector import BV

#MEMs file reader
def extract_output(in_file):
    outs = []
    with open(in_file) as o:
        for line in o.read().split("\n"):
            if line != "":
                l = line.split(" ")
                l.remove("")
                row = []
                if l[0] == '+':
                    row.append(True)
                else:
                    row.append(False)
                row.append(l[1])
                row.append(int(l[2]))
                mems = []
                for mem in l[3:]:
                    mems.append([int(x) for x in mem[1:-1].split(",")])
                row.append(mems)
                outs.append(row)
    return outs

#Index reader
def read_index(index_file):
    with open(index_file) as o:
        line = o.readline()
        i = 0
        while line:
            if i==0:
                (reference, ref_len) = line[:-1].split(" ")
            if i==1:
                text = line[:-1]
            if i==4:
                pos_S = line[:-2].split(" ")
                ExPos = []
                for p in pos_S:
                    x = p.split(",")
                    ExPos.append([int(x[0]), int(x[1])])
            i+=1
            line = o.readline()
    return reference, ref_len, text, ExPos

def main():
    in_file = sys.argv[1]
    index_file = sys.argv[2]

    outs = extract_output(in_file)
    reference, ref_len, text, ExPos = read_index(index_file)
    bv = BV(text)

    ExsN = len(ExPos)
    starts = [float('infinity') for i in range(ExsN)]
    ends = [-1 for i in range(ExsN)]
    ExsCov = {}
    edges = {}

    for out in outs:
        err = out[2]

        UsedExs = []
        
        (t,p,l) = out[3][0]
        index = bv.rank(t)
        UsedExs.append(index)
        if index in ExsCov:
            ExsCov.update({index:ExsCov[index]+1})
        else:
            ExsCov.update({index:1})

        start = ExPos[index-1][0]+(t-bv.select(index))
        if start < starts[index-1]:
            starts[index-1] = ExPos[index-1][0]+(t-bv.select(index)-1)
        (old_t,old_p,old_l) = (t,p,l)

        for (t,p,l) in out[3][1:]:
            old_index = bv.rank(old_t)
            index = bv.rank(t)
            if index not in UsedExs:
                UsedExs.append(index)
                if index in ExsCov:
                    ExsCov.update({index:ExsCov[index]+1})
                else:
                    ExsCov.update({index:1})
            if old_index != index:
                end = (old_t+old_l-1 - bv.select(old_index)-1) + ExPos[old_index-1][0]
                if end > ends[old_index-1]:
                    ends[old_index-1] = end
                start = ExPos[index-1][0]+(t-bv.select(index)-1)
                if start < starts[index-1]:
                    starts[index-1] = start

                key = str(old_index) + "_" + str(index)
                old_w = 0
                if key in edges:
                    old_w = edges[key]
                if err == 0:
                    w = 0
                else:
                    w = 1/err
                edges.update({key:old_w+w})
            (old_t,old_p,old_l) = (t,p,l)
        old_index = bv.rank(old_t)
        end = (old_t+old_l-1 - bv.select(old_index)-1) + ExPos[old_index-1][0]
        if end > ends[old_index-1]:
            ends[old_index-1] = end
    print(ExsCov)
    print(starts)
    print(ends)
    print(edges)

    i = 0
    while i<len(starts):
        s = starts[i]
        e = ends[i]
        (sp, ep) = ExPos[i]
        if s != float('infinity') and s-sp != 0:
            print("Exon {}: {} S".format(i+1, s-sp))
        if e != -1 and ep-e != 0:
            print("Exon {}: {} E".format(i+1, ep-e))
        i+=1            

if __name__ == '__main__':
    main()
