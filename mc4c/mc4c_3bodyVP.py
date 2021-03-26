#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import h5py
from numpy import *

def fragSpan(gs, ge, gz, rs, nb):
    ws = zeros(nb)
    bs = (gs-gz)/float(rs)
    be = (ge-gz)/float(rs)
    ns = max(int(floor(bs)), 0)
    ne = min(int(ceil(be)), nb)
    iw = [[], []] # index, weight
    if ns < ne:
        for n in range(ns, ne):
            if n < bs:
                if be < (n+1):
                    ws[n] = be-bs
                else:
                    ws[n] = (n+1)-bs
            else:
                if be < (n+1):
                    ws[n] = be-n
                else:
                    ws[n] = 1.0
        if sum(ws) > 0:
            ws = ws/sum(ws)
        for n in range(0, nb):
            if ws[n] > 0:
                iw[0].append(n)
                iw[1].append(ws[n])
    else:
        if (ns==ne) and (be>bs) and (ns<nb):
            iw[0].append(ns)
            iw[1].append(1.)
    return iw

# Calculate p_{i,j,vp} at a specific viewpoint
if not len(sys.argv) == 6:
    print('usage:: python mc4c_3bodyVP.py xxx.p3.h5 g_s(ROI) res[bp] g_s(VP) g_e(VP)')
    sys.exit()
fx = str(sys.argv[1])
gz = int(sys.argv[2])
rs = int(sys.argv[3])
gs = int(sys.argv[4])
ge = int(sys.argv[5])

if not os.path.isfile(fx):
    print('Cannot find '+fx)
    sys.exit()
else:
    fr = h5py.File(fx, 'r')
    N  = fr.get('cm').shape[0]
    vp = fragSpan(gs, ge, gz, rs, N)[0]
    p3 = zeros((N, N))
    for p in vp:
        p3 += fr.get('cm')[p]
    fr.close()
    #
    fw = open(fx[:-3]+".vp%d-%d.txt" % (gs, ge), 'w')
    ct = "#N: %d VP: " % (N)
    for p in vp:
        ct += "%d " % (p)
    ct+= "min: %11.5e max: %11.5e" % (nanmin(p3), nanmax(p3))
    fw.write(ct+'\n')
    for i in range(0, N):
        lt = ''
        for j in range(0, N):
            lt += "%11s "%('NaN') if isnan(p3[i,j]) else "%11.5e "%(p3[i,j])
        fw.write(lt+'\n')
    fw.close()
