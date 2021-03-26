#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import h5py
from numpy import *
from itertools import combinations

# Compute z_{ijk} of all triplets, and count their statistics by bioinformatic features
if not len(sys.argv) == 3:
    print('usage:: python sprite_3body_ENPS.py TAD-id[0-38] xxx.p3.h5')
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
    bts = array(bts).astype(int) #E/P/SE

# Compute z_{ijk}
fz = fx[:-3]+'.exp.h5'
if not os.path.isfile(fz):
    print('Cannot find '+fz)
    sys.exit()
else:
    fr = h5py.File(fz, 'r')
    p3_exp = fr.get('p3_exp')[()]
    fr.close()
#
abb = ['E','P','S','#'] # Enhancer/Promoter/Super-enhancer/Nothing
cbs = ['EPS','#PS','#ES','#EP', \
       '##E','##P','##S', \
       '#EE','EEP','EES', \
       '#PP','EPP','PPS', \
       '#SS','ESS','PSS', \
       '###','EEE','PPP','SSS','ALL'] # Combinations with characters in order

M = len(cbs)
zs= [[] for i in range(M)]
#
if not os.path.isfile(fx):
    print('Cannot find '+fx)
    sys.exit()
else:
    fr = h5py.File(fx, 'r')
    for i,j,k in combinations(range(N), 3):
        s = abs(array([i-j, i-k, j-k]))
        mis = min(s)
        mas = max(s)
        z = fr.get('cm')[i,j,k]/p3_exp[mas, mis, 0]
        zs[-1].append(z)
        if not isnan(z):
            #
            us = [abb[q] for q in range(0,3) if bts[i,q]==1]
            if us == []:
                us = [abb[3]]
            vs = [abb[q] for q in range(0,3) if bts[j,q]==1]
            if vs == []:
                vs = [abb[3]]
            ws = [abb[q] for q in range(0,3) if bts[k,q]==1]
            if ws == []:
                ws = [abb[3]]
            #
            for u in us:
                for v in vs:
                    for w in ws:
                        uvw = ''.join(sort([u,v,w]))
                        if not uvw in cbs:
                            print("%s is not included"%(uvw))
                            sys.exit()
                        else:
                            zs[cbs.index(uvw)].append( z )
    fr.close()
    # Statistics
    fw = open(fx[:-3]+'.ENPS', 'w')
    for m in range(M):
        lt = "%s " % (cbs[m])
        if zs[m] == []:
            lt += "N: %7d z: %12s %12s" % (0, 'NaN', 'NaN')
        else:
            zm = array(zs[m])
            zm_ln = log(zm[zm>0])
            lt += "N: %7d z: %11.5e %11.5e log(z): %+12.5e %11.5e" % (len(zm), nanmean(zm), nanstd(zm), mean(zm_ln), std(zm_ln))
        fw.write(lt+'\n')
    fw.close()

