import sys

from BitVector import BV

def extractMEMs(string):
    mems = []
    string_mems = string.split(" ")
    for s_mem in string_mems:
        mems.append([int(n) for n in s_mem[1:-1].split(",")])
    return mems

def main(path):
    with open("../tmp/T.fa") as t:
        text = t.read().split("\n")[1]
        bv = BV(text) #bv.rank, bv.select
    print(bv)
    
    with open("../tmp/e_pos") as f:
        exs_pos = []
        for line in f.read().split("\n"):
            if line != "":
                exs_pos.append([int(n) for n in line.split(",")])

    print(exs_pos)

    mems = extractMEMs(path)
    print(mems)
    i = 0
    CIGAR = ""
    matches = 0
    while i<len(mems):
        if i == 0:
            if mems[i][1] != 1:
                print("{}S".format(mems[i][1]))
                matches = 0
            matches += mems[i][2]
        else:
            if bv.rank(mems[i-1][0]) ==  bv.rank(mems[i][0]):
                errors_P = mems[i][1] - mems[i-1][1] - mems[i-1][2]
                errors_T = mems[i][0] - mems[i-1][0] - mems[i-1][2]
                #------------------------------------------------
                if errors_P == 0:
                    if errors_T < 0:
                        #Case 1
                        print("1")
                        if matches != 0:
                            print("{}M".format(matches))
                            matches = 0
                        print("{}D".format(abs(errors_T)))
                        matches = mems[i][2]-abs(errors_T)
                    elif errors_T > 0:
                        #Case 2
                        print("2")
                        if matches != 0:
                            print("{}M".format(matches))
                            matches = 0
                        print("{}I".format(errors_T))
                        matches = mems[i][2]
                #------------------------------------------------
                elif errors_P < 0:
                    if errors_T == 0:
                        #Case 3
                        print("3")
                        if matches != 0:
                            print("{}M".format(matches))
                            matches = 0
                        print("{}I".format(abs(errors_P)))
                        matches = mems[i][2]-abs(errors_P)
                    elif errors_T < 0:
                        #Case 4
                        if matches != 0:
                            print("{}M".format(matches))
                            matches = 0
                        errs = abs(errors_P) - abs(errors_T)
                        if errs > 0:
                            #Case 4a
                            print("4a")
                            print("{}I".format(errs))
                        elif errs < 0:
                            #Case 4b
                            print("4b")
                            print("{}D".format(abs(errs)))
                        matches = mems[i][2]-max(abs(errors_P), abs(errors_T))
                    elif errors_T > 0:
                        #Case 5
                        print("5")
                        if matches != 0:
                            print("{}M".format(matches))
                            matches = 0
                        errs = abs(errors_P) + errors_T
                        print("{}I".format(abs(errs)))
                        matches += mems[i][2]-abs(errors_P)
                #------------------------------------------------
                elif errors_P > 0:
                    if errors_T == 0:
                        #Case 6
                        print("6")
                        if matches != 0:
                            print("{}M".format(matches))
                            matches = 0
                        print("{}D".format(errors_P))
                        matches = mems[i][2]
                    elif errors_T < 0:
                        #Case 7
                        print("7")
                        if matches != 0:
                            print("{}M".format(matches))
                            matches = 0
                        errs = errors_P + abs(errors_T)
                        print("{}D".format(abs(errs)))
                        matches += mems[i][2]-abs(errors_T)
                    elif errors_T > 0:
                        #Case 8
                        if matches != 0:
                            print("{}M".format(matches))
                            matches = 0
                        errs = abs(errors_P) - abs(errors_T)
                        if errs > 0:
                            #Case 8a
                            print("8a")
                            print("{}I".format(errs))
                        elif errs < 0:
                            #Case 8b
                            print("8b")
                            print("{}D".format(abs(errs)))
                        matches = mems[i][2]-max(abs(errors_P), abs(errors_T))
                        

            else:
                pass
        i+=1
    if matches != 0:
        print("{}M".format(matches))

if __name__ == '__main__':
    main("(2,1,3) (5,6,3)")
