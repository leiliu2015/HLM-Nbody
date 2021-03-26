#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import h5py
from numpy import *
from scipy.special import erf
from scipy.integrate import romb
from itertools import permutations
from itertools import combinations

# n-body contact probability given a Gaussian form factor with specific VPs
def K2pnmGaus_VPs(K, n, rc2, vps):
    n = int(n)
    d = sum(K, axis=0)
    D = diag(d)
    km= D - K
    vps = list(vps)
    #
    km= km[1:,1:]
    Dk= linalg.det(km)
    N = len(K)
    pnm = zeros((N, N))*nan
    fDTA = 3./rc2
    for i in range(0, N):
        for j in range(i+1, N):
            if (not i in vps) and (not j in vps):
                A = [i,j] + vps
                DTA = zeros((N, N))
                for u in A:
                    for v in A:
                        if u == v:
                            DTA[u,u] = n-1.
                        else:
                            DTA[u,v] = -1.
                mDk = linalg.det(km + fDTA*DTA[1:,1:])
                pn = (Dk/mDk)**(1.5)
                pnm[i,j] = pn
                pnm[j,i] = pn
    return pnm

# Save a matrix of size n*n
def saveMx(fn, xy, ct):
    fw = open(fn, 'w')
    fw.write(ct)
    n = len(xy)
    for i in range(n):
        lt = ''
        for v in xy[i]:
            lt += "%11s "%('NaN') if isnan(v) else "%11.5e "%(v)
        fw.write(lt+'\n')
    fw.close()

# Calculate p_{n} at specific VPs assuming a Gaussian form factor
if not len(sys.argv) >= 3:
    print('usage:: python mc4c_nbodyVP.py xxx.K_fit VP1 [VP2...]')
    sys.exit()
fk = str(sys.argv[1])
vps= sort(list(map(int, sys.argv[2:])))
nvp= len(vps)
ndx= nvp+2
xi2= 1.0**2 # r^2_c

# Main
if True:
    # Read K-matrix
    if not os.path.isfile(fk):
        print('Cannot find '+fk)
        sys.exit()
    else:
        K_fit = []
        with open(fk) as fr:
            for line in fr:
                if not line[0] == '#':
                    lt = line.strip()
                    lt = lt.split()
                    K_fit.append( list(map(float, lt)) )
        K_fit = array(K_fit)
        N = len(K_fit)

    # Compute p_{n}
    P = K2pnmGaus_VPs(K_fit, ndx, xi2, vps)
    fo= fk+".q%d.vp"%(ndx)
    for vp in vps:
        fo += "%d-"%(vp)
    saveMx(fo[:-1]+'.txt', P, "#P_%d N %d min: %11.5e max: %11.5e\n"%(ndx, N, nanmin(P), nanmax(P)))



