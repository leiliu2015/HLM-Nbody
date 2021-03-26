#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
from numpy import *
set_printoptions(precision=3, linewidth=200)

# Subtract the association z-score among SOIs.
if not len(sys.argv) == 3:
    print('usage:: python mc4c_h5association_soi2zm.py soi.txt xxx.zm')
    sys.exit()
fs = str(sys.argv[1])
fz = str(sys.argv[2])
cq = True # Line:48

if not os.path.isfile(fs):
    print('Cannot find '+fs)
    sys.exit()
else:
    soi_lb = []
    soi_bs = []
    soi_be = []
    with open(fs) as f:
        for line in f:
            if not line[0] == '#':
                lt = line.strip().split()
                soi_lb.append( str(lt[0]) )
                soi_bs.append( str(lt[1]) )
                soi_be.append( str(lt[2]) )
    soi_bs = array(soi_bs, dtype=int)
    soi_be = array(soi_be, dtype=int)
    soi_N  = len(soi_lb)

if not os.path.isfile(fz):
    print('Cannot find '+fz)
    sys.exit()
else:
    zm = []
    with open(fz) as f:
        for line in f:
            if not line[0] == '#':
                lt = line.strip().split()
                zm.append( list(map(float, lt)) )
    zm = array(zm)
    N  = len(zm)

soi_zm = zeros((soi_N, soi_N))
for i in range(soi_N):
    for j in range(soi_N):
        if cq:
            soi_zm[i,j] = zm[soi_bs[i], soi_bs[j]]
        else:
            zs = []
            for u in range(soi_bs[i], soi_be[i]+1):
                for v in range(soi_bs[j], soi_be[j]+1):
                    zs.append( zm[u,v] )
            zs = array(zs)
            if (sum(isfinite(zs))==0) or (i==j):
                soi_zm[i,j] = nan
            else:
                soi_zm[i,j] = nanmean(zs)

fw = open(fz+'-soi2', 'w')
ct = "#fs: %s\n" % (fs)
ct+= "#fz: %s\n" % (fz)
ct+= "#shape: %d\n" % (soi_N)
ct+= 'xxx'
for lb in soi_lb:
    ct += "\t%s" % (lb)
fw.write(ct+'\n')
for i in range(soi_N):
    lt = "%s" % (soi_lb[i])
    for j in range(soi_N):
        lt += "\t%12s "%('NaN') if isnan(soi_zm[i,j]) else "\t%+12.5e "%(soi_zm[i,j])
    fw.write(lt+'\n')
fw.close()

