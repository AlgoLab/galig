import itertools

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

#################################
### ANCESTORS AND DESCENDANTS ###
#################################

#Utils
def bitOr(lists):
    res = lists[0]
    for l in lists[1:]:
        res = [i | j for i,j in zip(res,l)]
    return res

def bitAnd(list1, list2):
    return [i & j for i,j in zip(list1, list2)]

def getColumn(A, c):
    return [row[c] for row in A]

#Ancestors
def getAncestors(A, x):
    ancestors = []
    used = []
    toCheck = [getColumn(A,x)]
    used.append(x)
    while len(toCheck)>0:
        curr_c = toCheck.pop(0)
        ancestors.append(curr_c)
        r = 0
        for elem in curr_c:
            if elem == 1 and r not in used:
                used.append(r)
                next_c = getColumn(A,r)
                if sum(next_c) > 0:
                    toCheck.append(next_c)
            r+=1
    return bitOr(ancestors)

def getVector(x, n):
    return [0 for i in range(0,x)] + [1] + [0 for i in range(x+1, n)]

def getCommonAncestors(A, x, y):
    anc_x = getAncestors(A, x)
    anc_y = getAncestors(A, y)
           
    if sum(bitAnd(anc_x, getVector(y, len(A)))) == 0 and sum(bitAnd(anc_y, getVector(x, len(A)))) == 0:
        return bitAnd(anc_x, anc_y)
    return [0 for i in range(0,len(A))]

#Descendants
def getDescendants(A, x):
    descendants = []
    used = []
    toCheck = [A[x]]
    used.append(x)
    while len(toCheck)>0:
        curr_r = toCheck.pop(0)
        descendants.append(curr_r)
        c = 0
        for elem in curr_r:
            if elem == 1 and c not in used:
                used.append(c)
                next_r = A[c]
                if sum(next_r) > 0:
                    toCheck.append(next_r)
            c+=1
    return bitOr(descendants)

def getCommonDescendants(A, x, y):
    des_x = getDescendants(A, x)
    des_y = getDescendants(A, y)

    if sum(bitAnd(des_x, getVector(y, len(A)))) > 0 and sum(bitAnd(des_y, getVector(x, len(A)))) > 0:
        return bitAnd(des_x, des_y)
    return bitAnd(des_x, des_y)
