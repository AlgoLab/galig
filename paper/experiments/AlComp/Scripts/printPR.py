import sys

def main():
    fpath = sys.argv[1]

    TP = 0
    FP = 0
    FN = 0

    for line in open(fpath).readlines():
        gene,tp,fp,fn = line.strip("\n").split(",")
        TP+=int(tp)
        FP+=int(fp)
        FN+=int(fn)
    P = round(TP/(TP+FP), 3)
    R = round(TP/(TP+FN), 3)
    F = round(2*(P*R)/(P+R),3)

    print(P, R, F, sep="\t")

if __name__ == '__main__':
    main()
