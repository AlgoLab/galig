import sys, os

def main(files):
    ES = {}
    AS3 = {}
    AS5 = {}
    MEE = {}
    samples = []
    for f in files:
        sample = os.path.basename(f)
        samples.append(sample)
        events = [x.strip("\n").split(" ") for x in open(f).readlines()]
        for event in events:
            ev = "{}-{}".format(event[1], event[2])
            if event[0] == "ES":
                if ev not in ES:
                    ES[ev] = []
                ES[ev].append(sample)
            elif event[0] == "3'":
                if ev not in AS3:
                    AS3[ev] = []
                AS3[ev].append(sample)
            elif event[0] == "5'":
                if ev not in AS5:
                    AS5[ev] = []
                AS5[ev].append(sample)
            elif event[0] == "MEE":
                if ev not in MEE:
                    MEE[ev] = []
                MEE[ev].append(sample)

    print("Event,Info,{}".format(",".join(samples)))
    for e,ss in ES.items():
        line = "ES," + e + ","
        for s in samples:
            if s in ss:
                line += "1,"
            else:
                line += "0,"
        print(line[:-1])

    for e,ss in AS3.items():
        line = "AS3," + e + ","
        for s in samples:
            if s in ss:
                line += "1,"
            else:
                line += "0,"
        print(line[:-1])

    for e,ss in AS5.items():
        line = "AS5," + e + ","
        for s in samples:
            if s in ss:
                line += "1,"
            else:
                line += "0,"
        print(line[:-1])

    for e,ss in MEE.items():
        line = "MEE," + e + ","
        for s in samples:
            if s in ss:
                line += "1,"
            else:
                line += "0,"
        print(line[:-1])

if __name__ == '__main__':
    main(sys.argv[1:])
