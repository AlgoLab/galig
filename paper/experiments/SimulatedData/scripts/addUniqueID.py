import sys, os

def main():
    f = sys.argv[1]
    events = {"ES":0, "A3":0, "A5":0, "IR":0}
    for line in open(f).readlines():
        ev = line.strip("\n").split(" ")
        if ev[1] in events:
            events[ev[1]]+=1
            ev.append("{}{}".format(ev[1], events[ev[1]]))
            print(" ".join(ev))

if __name__ == '__main__':
    main()
