import sys, re

class Feature:
    def __init__(self, line):
        line_tmp = line.split("\t")
        self.ref = line_tmp[0]
        start = int(line_tmp[3])
        end = int(line_tmp[4])
        self.strand = line_tmp[6]

        attributes = line_tmp[8].strip("\n")

        structure_group = re.findall("structure \"([0-9\^\-]*,[0-9\^\-]*)\";", attributes)
        self.structure = structure_group[0]

        all_splices = re.findall("splice_chain \"([0-9\^\-,]+)\";", attributes)[0].replace(",", "")
        splices_group = re.findall("([0-9]+[-\^])", all_splices)
        self.splices = []
        for splice in splices_group:
            self.splices.append(int(splice[:-1]))

        if self.strand == '+':
            self.flanks = (start, end)
            self.splices.sort()
        else:
            self.flanks = (end, start)
            self.splices.sort(reverse = True)

    def getTranscripts(self):
        return self.transcripts

    def getFlank1(self):
        return self.flanks[0]

    def getFlank2(self):
        return self.flanks[1]

    def getSplice(self, i):
        return int(self.splices[i-1])

    def getSplices(self):
        splices = []
        for splice in self.splices:
            splices.append(splice)
        return splices

    def printRecord(self):
        f = "{} {}".format(self.getFlank1(), self.getFlank2())
        s = self.getSplices()
        print("------------------")
        print(self.structure)
        print(f)
        print(s)
        print("------------------")

def printES(f,p1,p2):
    print("{0} {1} {2} {3}".format(f.ref, "ES", min(p1,p2)+1, max(p1,p2)-1))

def printIR(f,i,j):
    i = f.getSplice(i)
    j = f.getSplice(j)
    print("{0} {1} {2} {3}".format(f.ref, "IR", min(i,j)+1, max(i,j)-1))

def printA3(f, i, j):
    i = f.getSplice(i)
    j = f.getSplice(j)
    print("{0} {1} {2} {3} {4}".format(f.ref, "A3", min(i,j), max(i,j), f.getFlank1()))

def printA5(f, i, j):
    i = f.getSplice(i)
    j = f.getSplice(j)
    print("{0} {1} {2} {3} {4}".format(f.ref, "A5", min(i,j), max(i,j), f.getFlank2()))

def main():
    events_path = sys.argv[1]
    for line in open(events_path, 'r'):
        f = Feature(line)
        if f.structure[0:4] == "0,1-":
            printES(f, f.getFlank1(), f.getFlank2())
        elif f.structure == "1^,2^":
            printA5(f, 1, 2)
        elif f.structure == "1-,2-":
            printA3(f, 1, 2)
        elif f.structure == "0,1^2-":
            printIR(f,1,2)
        elif f.structure == "0,1^2-3^4-":
            printIR(f,1,2)
            printIR(f,3,4)
        elif f.structure == "0,1^2-3^4-5^6-":
            printIR(f,1,2)
            printIR(f,3,4)
            printIR(f,5,6)
        elif f.structure == "0,1^2-3^4-5^6-7^8-":
            printIR(f,1,2)
            printIR(f,3,4)
            printIR(f,5,6)
            printIR(f,7,8)
        elif f.structure == "1^2-,3^4-":
            printIR(f,1,2)
            printIR(f,3,4)
        else:
            pass

if __name__ == '__main__':
    main()
