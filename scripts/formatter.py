import os, sys
from BitVector import BV

def main(out_path):
    with open("./tmp/T.fa") as t:
        text = t.read().split("\n")[1]
        bv = BV(text)

    f = open("{}_withRef".format(out_path), "w")
    with open(out_path) as o:
        for line in o:
            if line != "":
                strings = line.split(" ")
                p_id, mems_string, err = strings[0], strings[1:-1], int(strings[-1])
                mems = []
                for mem in mems_string:
                    if mem != "" and mem != "\n":
                        m = mem[1:-1].split(",")
                        mems.append({"t":int(m[0]), "p":int(m[1]), "l":int(m[2])})
                T = text[mems[0]["t"]-1:mems[0]["t"]+mems[0]["l"] - 1]
                t = mems[0]["t"]+mems[0]["l"]
                i = 1
                while i<len(mems):
                    if bv.rank(mems[i-1]["t"] - 1) != bv.rank(mems[i]["t"] - 1):
                        T += text[t:bv.select(bv.rank(mems[i-1]["t"] - 1) + 1)]
                        T += text[bv.select(bv.rank(mems[i]["t"] - 1)):mems[i]["t"]+mems[i]["l"] - 1]
                    else:
                        T += text[t:mems[i]["t"]+mems[i]["l"] - 1]
                    t = mems[i]["t"]+mems[i]["l"]
                    i+=1

                f.write("{}\t{}\t{}\t{}\n".format(p_id, " ".join(mems_string[:-1]), err, T))
    f.close()

if __name__ == '__main__':
    main(sys.argv[1])
