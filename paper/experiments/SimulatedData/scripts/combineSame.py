import sys

def main():
    f = sys.argv[1]
    Sames = []
    Equals = set()
    Total = set()
    for line in open(f).readlines():
        x,y,d = line.strip("\n").split(" ")
        Total.add(x)
        Total.add(y)
        if d == "0":
            Equals.add(x)
            Equals.add(y)
            found = False
            for S in Sames:
                if x in S or y in S:
                    found = True
                    S.add(x)
                    S.add(y)
            if not found:
                Sames.append(set([x,y]))
    i = 0
    for S in Sames:
        i+=1
        print("E{}".format(i), " ".join(list(S)))
    for t in Total:
        if t not in Equals:
            i+=1
            print("E{}".format(i), t)
            

if __name__ == '__main__':
    main()
