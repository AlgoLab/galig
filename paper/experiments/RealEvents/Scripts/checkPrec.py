import sys, os

def main():
    eventsPath = sys.argv[1]
    truthPath = sys.argv[2]
    allEventsPath = sys.argv[3]
    chrom = sys.argv[4]
    gene = sys.argv[5]

    Truth = {'ES':[], 'A3':[], 'A5':[], 'IR':[]}
    for line in open(truthPath).readlines():
        _,evT,p1,p2 = line.strip("\n").split(" ")
        p1,p2 = int(p1),int(p2)
        p1,p2 = min(p1,p2), max(p1,p2)
        Truth[evT[0:2]].append((p1,p2))

    allEvents = {'ES':[], 'A3':[], 'A5':[], 'IR':[]}
    for line in open(allEventsPath).readlines():
        _,evT,p1,p2 = line.strip("\n").split(" ")
        p1,p2 = int(p1),int(p2)
        p1,p2 = min(p1,p2), max(p1,p2)
        allEvents[evT[0:2]].append((p1,p2))

    for line in open(eventsPath).readlines():
        if line[0:4] == "Type":
            continue
        line = line.strip("\n").split(",")
        p1,p2 = int(line[1]), int(line[2])
        p1,p2 = min(p1,p2), max(p1,p2)
        if (p1,p2) in Truth[line[0]]:
            print(chrom,gene,line[0],p1,p2,1)
        else:
            if (p1,p2) not in allEvents[line[0]]:
                print(chrom,gene,line[0],p1,p2,0)

if __name__ == '__main__':
    main()
