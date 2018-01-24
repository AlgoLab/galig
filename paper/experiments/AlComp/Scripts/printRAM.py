import sys

def main():
    fpath = sys.argv[1]
    Rs = []
    for line in open(fpath).readlines():
        dat = int(line.strip("\n"))
        Rs.append(dat)

    print(max(Rs), min(Rs))

if __name__ == '__main__':
    main()
