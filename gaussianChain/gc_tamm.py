#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import h5py
from numpy import *
from itertools import combinations
set_printoptions(precision=3, linewidth=200)
import warnings
warnings.filterwarnings('ignore')

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

# Calculate the boost to triplet contacts correlator, namely f_{2}(i,j,k) = p_{i,j,k}/p_{i,j}/p_{j,k} at i<j<k,
# which has benn analyzed by M.V.Tamm in PRE 99,032501 (2019).
if not len(sys.argv) == 4:
    print('usage:: python gc_tamm.py N fdx[0/1] k0')
    sys.exit()
N  = int(sys.argv[1])
fdx= int(sys.argv[2])
k0 = float(sys.argv[3])

xi2= 1.0**2 # r^2_c
pq = 'p' if fdx else 'q'
aq = True # re-compute f_{2}(i,j,k) by taking advantage of Gaussian chain
lb = "g%d"%(N)
fy = lb + ".tamm-F%d"%(fdx)

if (fdx==0) and (aq):
    # Compute f_{2}(i,v,j) based on the property of Gaussian chain: sigma_{i,j} = min(i,j)/k_{0}
    F2 = zeros((N,N,N))*nan
    for v in range(0, N):
        for i in range(0, N):
            for j in range(0, N):
                if (i<v) and (v<j):
                    s1 = v-i
                    s2 = j-v
                    f2 = ((k0*k0 + 6.*k0*(s1+s2) + 27.*s1*s2) / (k0*k0 + 3.*k0*(s1+s2) + 9.*s1*s2))**(-1.5)
                    F2[v,i,j] = f2
                    F2[v,j,i] = f2
else:
    # Read p_{i,j}
    fx = lb + ".%s2.h5" % (pq)
    if not os.path.isfile(fx):
        print('Cannot find '+fx)
        sys.exit()
    else:
        fr = h5py.File(fx, 'r')
        P2 = fr.get('cm')[()]
        fr.close()
    # Read p_{i,j,k}
    fx = lb + ".%s3.h5" % (pq)
    if not os.path.isfile(fx):
        print('Cannot find '+fx)
        sys.exit()
    else:
        fr = h5py.File(fx, 'r')
        P3 = fr.get('cm')[()]
        fr.close()
    # Compute f_{2}(i,v,j)
    F2 = zeros((N,N,N))*nan
    for v in range(0, N):
        for i in range(0, N):
            for j in range(0, N):
                if (i<v) and (v<j):
                    f2 = P3[v,i,j]/P2[v,i]/P2[v,j]
                    F2[v,i,j] = f2
                    F2[v,j,i] = f2
saveTsTx(fy+'.txt', F2, "#N: %d aq: %d\n"%(N, aq))

