import sys, os

def main():
    eventsPath = sys.argv[1]
    truthPath = sys.argv[2]
    allEventsPath = sys.argv[3]
    chrom = sys.argv[4]
    gene = sys.argv[5]
    
    total = 0
    Truth = {'ES':[], 'A3':[], 'A5':[], 'IR':[]}
    for line in open(truthPath).readlines():
        _,evT,p1,p2 = line.strip("\n").split(" ")
        p1,p2 = int(p1),int(p2)
        p1,p2 = min(p1,p2), max(p1,p2)
        Truth[evT[0:2]].append((p1,p2))
        total+=1

    allEvents = {'ES':[], 'A3':[], 'A5':[], 'IR':[]}
    for line in open(allEventsPath).readlines():
        _,evT,p1,p2 = line.strip("\n").split(" ")
        p1,p2 = int(p1),int(p2)
        p1,p2 = min(p1,p2), max(p1,p2)
        allEvents[evT[0:2]].append((p1,p2))

    TP = 0
    FP = 0
    FN = 0
    FPs = {'ES':set(), 'A3':set(), 'A5':set(), 'IR':set()}
    for line in open(eventsPath).readlines():
        if line[0:4] == "Type":
            continue
        line = line.strip("\n").split(",")
        p1,p2 = int(line[1]), int(line[2])
        p1,p2 = min(p1,p2), max(p1,p2)
        if (p1,p2) in Truth[line[0]]:
            TP+=1
        else:
            if (p1,p2) not in allEvents[line[0]]:
                FP+=1
                FPs[line[0]].add("{}-{}".format(str(p1),str(p2)))
    FN = total-TP
    print(chrom, gene, TP, FP, FN, " ".join(FPs['ES']), " ".join(FPs['IR']), " ".join(FPs['A3']), " ".join(FPs['A5']), sep=",")

if __name__ == '__main__':
    main()
