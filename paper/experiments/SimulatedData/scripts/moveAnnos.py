import sys, os

def main():
    infoPath = sys.argv[1]
    folder = sys.argv[2]

    for line in open(infoPath).readlines():
        newName, oldName, *toDel = line.strip("\n").split(" ")
        os.rename(folder + "/" + oldName + ".gtf", folder + "/" + newName + ".gtf")
        for n in toDel:
            os.remove(folder + "/" + n + ".gtf")

if __name__ == '__main__':
    main()
