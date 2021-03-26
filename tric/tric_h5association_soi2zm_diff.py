#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
from numpy import *
set_printoptions(precision=3, linewidth=200)

# Compute z(cellLine1) - z(cellLine2) among sites of interests
if not len(sys.argv) == 3:
    print('usage:: python diff.py cellLine1.zm-soi2 cellLine2.zm-soi2')
    sys.exit()
fx = [str(sys.argv[1]), str(sys.argv[2])]

soi_zm = []
soi_lb = []
for f in fx:
    if not os.path.isfile(f):
        print('Cannot find '+f)
        sys.exit()
    else:
        soi_zm.append([])
        with open(f) as fr:
            for line in fr:
                if not line[0]=='#':
                    lt = line.strip().split()
                    if not lt[0] == 'xxx':
                        soi_zm[-1].append( list(map(float, lt[1:])) )
                    else:
                        if soi_lb == []:
                            soi_lb = lt[1:]

soi_zm = array(soi_zm)
soi_N  = shape(soi_zm)[1]
soi_dz = soi_zm[0] - soi_zm[1]

fw = open(fx[0]+'-diff', 'w')
ct = "#fa: %s\n" % (fx[0])
ct+= "#fb: %s\n" % (fx[1])
ct+= "#shape: %d fa-fb\n" % (soi_N)
ct+= 'xxx'
for lb in soi_lb:
    ct += "\t%s" % (lb)
fw.write(ct+'\n')
for i in range(soi_N):
    lt = "%s" % (soi_lb[i])
    for j in range(soi_N):
        lt += "\t%12s "%('NaN') if isnan(soi_dz[i,j]) else "\t%+12.5e "%(soi_dz[i,j])
    fw.write(lt+'\n')
fw.close()
