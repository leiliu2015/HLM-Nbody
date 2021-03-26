#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
from numpy import *

# Average xxx.log over different replicas
if not len(sys.argv) == 4:
    print("usage:: python averageLogs.py dataDir Nxx-Fx-Lx-Mx numberOfReplicas")
    sys.exit()
ddir=str(sys.argv[1])
fg = str(sys.argv[2])
nr = int(sys.argv[3])

c_trajs = [] # {cost, pearsonr, time}
for rdx in range(nr):
    fx = "%s/%s-r%d.log"%(ddir, fg, rdx)
    if not os.path.isfile(fx):
        print('Cannot find '+fx)
    else:
        c_trajs.append([])
        with open(fx) as fr:
            for line in fr:
                if not line[0] == '#':
                    lt = line.strip()
                    lt = lt.split()
                    c_trajs[-1].append( list(map(float, lt)) )
c_trajs = array(c_trajs)
iterations = shape(c_trajs)[1]

fw = open("%s/%s-rx.log"%(ddir, fg), 'w')
fw.write('#<cost> min_c max_c <pearsonr> min_r max_r <time>\n')
for i in range(iterations):
    ys = c_trajs[:,i,0]
    lt = "%11.5e %11.5e %11.5e   " % (mean(ys), min(ys), max(ys))
    ys = c_trajs[:,i,1]
    lt+= "%11.5e %11.5e %11.5e   " % (mean(ys), min(ys), max(ys))
    ys = c_trajs[:,i,2]
    lt+= "%11.5e" % (mean(ys))
    fw.write(lt+'\n')
fw.close()

