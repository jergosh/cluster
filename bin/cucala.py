import copy
import random
from numpy import mean, sqrt, zeros
from scipy import stats

def statscanMW(vec):
    v = stats.rankdata(vec)
    n = len(v)
    
    tab = zeros([n*(n-1)/2, 3])
    ind = 0
    ind_max = None
    P_max = 0.0

    for i in range(0, n-1):
        sumrank = v[i]

        for j in range(i+1, n):
            sumrank += v[j]
            tab[ind, 0] = i
            tab[ind, 1] = j
            m = j - i + 1

            if m < n:
                tab[ind, 2] = (sumrank-m*(n+1)/2)/sqrt(m*(n-m)*(n+1)/12)
            else:
                tab[ind, 2] = 0

            if tab[ind, 2] > P_max:
                P_max = tab[ind, 2]
                ind_max = ind

            ind += 1

    return tab[ind_max, ]

def signscanMW(vec, niter):
    v = stats.rankdata(vec)
    n = len(v)
    pval = 0

    res = statscanMW(vec)
    print res
    stat = res[2]

    for i in range(niter):
        vec_sim = random.sample(vec, n)
        res_sim = statscanMW(vec_sim)
        print res_sim
        if stat > res_sim[2]:
            pval += 1

    pval = float(pval)/(niter+1)

    return [res, 1-pval]

dist = lambda c1, c2: reduce(lambda x1, x2: x1+x2, [ (c[0]-c[1])**2 for c in zip(c1, c2) ])
# dist_naive = lambda c1, c2: (c1[0]-c2[0])**2+(c1[1]-c2[1])**2
    

# The data is locations + marks -- how to best represent it, 
# bearing in mind that we shall need to permute the marks
def order_dists(coords):
    """
    Preprocessing: calculating the distances to be able to iterate more 
    efficiently
    """

    n = len(coords)
    dists = [None]*n

    for i, c1 in enumerate(coords):
        ds = [ (j, dist(c1, c2)) for j, c2 in enumerate(coords) ]
        dists[i] = sorted(ds, key=lambda k: k[1])

    return dists 

def MWnaive(coords, marks):
    assert len(coords) == len(marks)

    n = len(coords)
    maxI = 0.0
    maxCoords = None
    maxR = 0.0
    maxV = []
    E = mean(marks)

    for i, c1 in enumerate(coords):
        for j, c2 in enumerate(coords):
            if marks[j] == 0:
                continue

            R = dist(c1, c2)
            m = 0
            SV = 0
            v = []

            for k, c3 in enumerate(coords):
                if dist(c1, c3) <= R:
                    m += 1
                    SV += marks[k]
                    v.append(k)

            I = 0
            if m < n:
                I = (SV-m*E)/sqrt(m*(n-m))
            if I > maxI:
                maxI = I
                maxCoords = coords[i]
                maxR = R
                maxV = v

    return maxI, maxCoords, maxR, maxV
            

def MWcont(coords, marks, dists):
    assert len(coords) == len(marks)

    n = len(coords)
    maxI = 0.0
    maxCoords = None
    maxR = 0.0
    maxV = []
    E = mean(marks)

    for i, c1 in enumerate(coords):
        m = 0
        SV = 0
        v = []
        for j, R in dists[i]:
            SV += marks[j]
            m += 1

            v.append(j)

            if marks[j] == 0:
                continue

            if m < n:
                I = (SV-m*E)/sqrt(m*(n-m))

                if I > maxI:
                    maxI = I
                    maxCoords = coords[i]
                    maxR = R
                    maxV = copy.copy(v)

    return maxI, maxCoords, maxR, maxV
      
def signMWcont(coords, marks, dists, niter):
    maxI, maxCoords, maxR, maxV = MWcont(coords, marks, dists)

    pval = 1
    for i in range(niter):
        print i,
        marks_p = random.sample(marks, len(marks))
        I, c, R, v =  MWcont(coords, marks_p, dists)
        print I
        if I >= maxI:
            pval += 1

    pval = float(pval) / (niter+1)
    return maxI, maxCoords, maxR, maxV, pval
    

def main():
    coords = []
    marks = []
    datfile = open("data/chorley.dat")
    
    datfile.readline() # Skip header

    for l in datfile:
        f = [ float(f) for f in l.rstrip().split() ]
        coords.append(f[0:2])
        marks.append(f[2])

    dists = order_dists(coords)
    # print MWnaive(coords, marks)
    maxI, maxCoords, maxR, maxV, pval = signMWcont(coords, marks, dists)
    print maxI, maxCoords, maxR, maxV, pval

if __name__ == "__main__":
    main()
