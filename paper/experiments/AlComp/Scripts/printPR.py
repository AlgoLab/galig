import sys

def printPRF(ev, dats):
    if dats['TP']+dats['FP'] == 0:
        P = 0
    else:
        P = round(dats['TP']/(dats['TP']+dats['FP']),3)
    if dats['TP']+dats['FN'] == 0:
        R = 0
    else:
        R = round(dats['TP']/(dats['TP']+dats['FN']),3)
    if P+R != 0:
        F = round(2*(P*R)/(P+R),3)
    else:
        F = 0
    print(ev,P,R,F,sep=',')

def main():
    fpath = sys.argv[1]

    TP = 0
    FP = 0
    FN = 0

    for line in open(fpath).readlines():
        chrom,gene,tp,fp,fn = line.strip("\n").split(",")
        TP+=int(tp)
        FP+=int(fp)
        FN+=int(fn)
    P = round(TP/(TP+FP), 3)
    R = round(TP/(TP+FN), 3)
    F = round(2*(P*R)/(P+R),3)
    print("P", "R", "F", sep="\t")
    print(P, R, F, sep="\t")

        

if __name__ == '__main__':
    main()
