#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import h5py
from numpy import *

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

# Calculate p_{i,j,vp} at specific promoters
if not len(sys.argv) == 3:
    print('usage:: python sprite_3bodyPromoter-chromatixTADs.py TAD-id[0-38] xxx.p3.h5')
    sys.exit()
tdx = int(sys.argv[1])
fx = str(sys.argv[2])

# Read TAD info
fs = 'sprite-data/CHROMATIX-suppTable3b.txt'
if not os.path.isfile(fs):
    print('Cannot find '+fs)
    sys.exit()
else:
    TS = "TAD%d:"%(tdx)
    sp = 0
    with open(fs) as f:
        for line in f:
            if line == '\n':
                sp = 0
            else:
                if TS in line:
                    sp = 1
                    lt = line.strip().split()
                    cx, se = lt[1].split(':')
                    roi_gs, roi_ge = list(map(int, se.split('-')))
                    N = int(lt[-1])
                    bts = [] # bead type
                else:
                    if (sp==1) and (not line[0]=='#'):
                        lt = line.strip().split()
                        bts.append( list(map(int, lt[1:4])) )
    bts = array(bts).astype(int)
    #
    pdx = where(bts[:,1]==1)[0] # promoter indices
    PN = len(pdx) # number of promoters
    if PN == 0:
        print("No promoter is found in TAD%d" % (tdx))
        sys.exit()
    else:
        pbx = zeros(PN).astype(int) # promoter-block indices
        for i in range(1, PN):
            if pdx[i] == pdx[i-1]+1:
                pbx[i] = pbx[i-1]
            else:
                pbx[i] = pbx[i-1]+1
        PBN = max(pbx)+1 # number of promoter-blocks

# Read p_{ijk}
if not os.path.isfile(fx):
    print('Cannot find '+fx)
    sys.exit()
else:
    p3_vp = zeros((N, N, PN))
    fr = h5py.File(fx, 'r')
    for p in range(0, PN):
        p3_vp[:,:,p] = fr.get('cm')[pdx[p]]
    fr.close()

# Read p^{expectation}_{ijk}
fz = fx[:-3]+'.exp.h5'
if not os.path.isfile(fz):
    print('Cannot find '+fz)
    sys.exit()
else:
    fr = h5py.File(fz, 'r')
    p3_exp = fr.get('p3_exp')[()]
    fr.close()
    #
    z3_vp = zeros((N, N, PN))*nan
    for p in range(0, PN):
        k = pdx[p]
        for i in range(0, N):
            for j in range(i+1, N):
                if (not i==k) and (not j==k):
                    s = abs(array([i-j, i-k, j-k]))
                    mis = min(s)
                    mas = max(s)
                    z3_vp[i,j,p] = p3_exp[mas, mis, 0]
                    z3_vp[j,i,p] = p3_exp[mas, mis, 0]

# Save TXT (merge consecutive promoters into promoter-blocks)
for b in range(0, PBN):
    ps = where(pbx==b)[0]
    bs = min(pdx[ps])
    be = max(pdx[ps])
    fo = fx[:-3]+".vp%03d-%03d" % (bs, be)
    #
    P = zeros((N, N))
    for p in ps:
        P += p3_vp[:,:,p]
    P /= float(len(ps))
    saveMx(fo+'.pm', P, "#N: %d min: %11.5e max: %11.5e\n"%(N, nanmin(P), nanmax(P)))
    #
    E = zeros((N, N))
    for p in ps:
        E += z3_vp[:,:,p]
    E /= float(len(ps))
    saveMx(fo+'.em', E, "#N: %d min: %11.5e max: %11.5e\n"%(N, nanmin(E), nanmax(E)))
    #
    Z = P/E
    saveMx(fo+'.sm', Z, "#N: %d min: %11.5e max: %11.5e\n"%(N, nanmin(Z), nanmax(Z)))






