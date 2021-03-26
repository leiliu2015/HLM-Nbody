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

# Save an array of size n*m
def saveLg(fn, xy, ct):
    fw = open(fn, 'w')
    fw.write(ct)
    n = len(xy)
    m = shape(xy)[1] if xy.ndim==2 else 0
    for i in range(n):
        if m == 0:
            lt = "%11s "%('NaN') if isnan(xy[i]) else "%11.5e"%(xy[i])
        else:
            lt = ''
            for v in xy[i]:
                lt += "%11s "%('NaN') if isnan(v) else "%11.5e "%(v)
        fw.write(lt+'\n')
    fw.close()

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

# Repeat Tanay's analysis @Apr07.2020 by LL.
if not len(sys.argv) == 5:
    print('usage:: python gc_tanay.py N fdx[0/1] vp_index[0,N) TDX[1/2]')
    sys.exit()
N  = int(sys.argv[1])
fdx= int(sys.argv[2])
vp = int(sys.argv[3])
tx = int(sys.argv[4]) # Different way to sample {vp,i} in the randomized triplets set
normalizeq = False

pq = 'p' if fdx else 'q'
lb = "g%d"%(N)
vp_lb = "vp%02d"%(vp)
fy = lb + ".tanay-F%d"%(fdx)

if True:
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

    # Compute marginal probability p(vp, i)
    if tx == 1:
        mp = zeros(N)
        for i in range(0, N):
            mp[i] = nansum(P3[vp,i,:])
        mp[vp] = nan
    else:
        mp = P2[vp,:]
    #
    smp = nansum(mp)
    vsp = zeros((N, 3))
    vsp[:,0] = arange(N)
    vsp[:,1] = mp*1.0
    vsp[:,2] = mp/smp
    saveLg(fy+".%s.T%d.vpsoiProfile"%(vp_lb, tx), vsp, '#i '+('p_{vp,i}' if tx==2 else 'sum_j(p_{vp,i,j})')+' frac\n')

    # Compute hub-score
    rad = zeros((N, N))*nan
    for i in range(0, N):
        for j in range(0, N):
            rad[i,j] = (mp[i]*mp[j] + mp[i]*P2[i,j] + mp[j]*P2[i,j])/3.0
    if normalizeq:
        srad = nansum(rad)
        rad = rad/srad*nansum(P3[vp])
    hs = P3[vp]/rad
    saveMx(fy+".%s.T%d.hs.txt"%(vp_lb, tx), hs, "#N: %d min: %11.5e max: %11.5e\n"%(N, nanmin(hs), nanmax(hs)))

