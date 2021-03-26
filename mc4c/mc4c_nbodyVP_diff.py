#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
from numpy import *

if not len(sys.argv) == 3:
    print('usage:: python mc4c_nbodyVP_diff.py cellLine1.txt cellLine2.txt')
    sys.exit()
fx = [str(sys.argv[1]), str(sys.argv[2])]

pns= []
for f in fx:
    if not os.path.isfile(f):
        print('Cannot find '+f)
        sys.exit()
    else:
        pns.append([])
        with open(f) as fr:
            for line in fr:
                if not line[0]=='#':
                    lt = line.strip().split()
                    pns[-1].append( list(map(float, lt)) )

pns = array(pns)
dpn = pns[0] - pns[1]
dpn2= dpn/pns[1]
N = len(dpn)

# Absolute difference
fw = open(fx[0]+'-diff', 'w')
ct = "#fa: %s\n" % (fx[0])
ct+= "#fb: %s\n" % (fx[1])
ct+= "#shape: %d fa-fb min: %+12.5e max: %+12.5e\n" % (N, nanmin(dpn), nanmax(dpn))
fw.write(ct)
for i in range(N):
    lt = ''
    for j in range(N):
        lt += " %12s "%('NaN') if isnan(dpn[i,j]) else " %+12.5e "%(dpn[i,j])
    fw.write(lt+'\n')
fw.close()

