import sys

def splitTime(time):
    t = time.split(":")
    if len(t)<3:
        m = int(t[0])
        s,d = int(t[1].split(".")[0]), int(t[1].split(".")[1])
        if d > 50:
            s+=1
        s+=60*m
    else:
        h = int(t[0])
        m = int(t[1])
        s = int(t[2])
        s+=3600*h+60*m
    if s == 0:
        s = 1
    return s

def main():
    fpath = sys.argv[1]
    Ts = []
    for line in open(fpath).readlines():
        dat = line.strip("\n")
        Ts.append(splitTime(dat))

    print(sum(Ts), min(Ts), max(Ts))#, mean(Ts))

if __name__ == '__main__':
    main()
