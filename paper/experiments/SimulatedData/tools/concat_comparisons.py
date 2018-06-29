import sys

for f in sys.argv[1:]:
    with open(f, "r") as res:
        res.readline()
        for l in res.readlines():
            print(l.rstrip())
