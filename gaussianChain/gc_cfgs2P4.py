#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import h5py
from numpy import *
from numba import jit
from itertools import combinations
from itertools import permutations
import warnings
warnings.filterwarnings('ignore')

# Compute p_{i,j,k,l} based on 3D coordinates assuming a Heaviside-step form factor,
# whose analytical solution is not available yet. @Apr.18.2020 by LL
if not len(sys.argv) == 4:
    print('usage:: python gc_cfgs2P4.py N freq h5q[0/1]')
    sys.exit()
N = int(sys.argv[1])
fq = int(sys.argv[2])
h5q = int(sys.argv[3])

ndx= 4
fdx= 1 # Heaviside step form factor
xi2= 1.0
fo = "g%d.xyz.p4" % (N)

# Compute binary contact tensor based on 3D coordinates
@jit(nopython=True)
def xyz2P4(xyz, N, xi2, pair, nqua):
    B = zeros((N, N))
    for i in range(0, N):
        for j in range(i+1, N):
            v = xyz[j] - xyz[i]
            r2= dot(v,v)
            if r2 < xi2:
                B[i,j] = 1
                B[j,i] = 1
    C = zeros(nqua) + 1
    for q in range(0, nqua):
        for p in range(0, 6):
            i = pair[q,p,0]
            j = pair[q,p,1]
            if B[i,j] == 0:
                C[q] = 0
                break
    return C

# Save contact probability tensor in a TXT file
def saveTsTx(fn, ts, ct):
    d = len(shape(ts))
    N = shape(ts)[0]
    dv= d - 2 # dimension of viewpoints
    #
    fw= open(fn, 'w')
    fw.write(ct)
    bx= 0
    for index in combinations(list(range(N)), dv):
        mx = ts[index]
        ct = '#vp: '
        for i in index:
            ct += "%d " % (i)
        ct+= "min: %11.5e " % (nanmin(mx))
        ct+= "max: %11.5e " % (nanmax(mx))
        ct+= "block-index: %d " % (bx)
        fw.write(ct+'\n')
        bx += 1
        for i in range(0, N):
            lt = ''
            for j in range(0, N):
                lt += "%11s "%('NaN') if isnan(mx[i,j]) else "%11.5e "%(mx[i,j])
            fw.write(lt+'\n')
        fw.write('\n\n')
    fw.close()

# Save contact probability tensor in a H5 file
def saveTsH5(fn, ts):
    fw = h5py.File(fn, 'w')
    fw.create_dataset('cm', data=ts, dtype='f')
    fw.close()

# Main
if True:
    fx = "g%d.xyz.h5" % (N)
    if not os.path.isfile(fx):
        print('Cannot find '+fx)
        sys.exit()
    else:
        pair = []
        for ijkl in list(combinations(list(range(N)), 4)):
            pair.append( list(combinations(ijkl, 2)) )
        pair = array(pair, dtype=int)
        nqua = shape(pair)[0]
        cqua = zeros(nqua)

        fr= h5py.File(fx, 'r')
        fmx = fr.get('xyz').shape[0]
        nsample = 0
        for f in range(0, fmx, fq):
            R = fr.get('xyz')[f]
            cqua += xyz2P4(R, N, xi2, pair, nqua)
            nsample += 1
            #
            if nsample%100 == 0: print("3D2P4: %5.1f"%(100.*f/float(fmx)) + ' %')
        fr.close()

        P = zeros((N,N,N,N))*nan
        q = 0
        for ijkl in list(combinations(list(range(N)), 4)):
            p4 = cqua[q]/float(nsample)
            q += 1
            for index in permutations(ijkl):
                P[index] = p4

        if h5q:
            saveTsH5(fo+'.h5',  P)
        else:
            saveTsTx(fo+'.txt', P, "#N: %d fdx: %d freq: %d\n"%(N, fdx, fq))
