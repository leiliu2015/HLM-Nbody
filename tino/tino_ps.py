#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
from numpy import *

# Calculate P(s) from p_{ij}
def calPs(cm, N):
    ps = zeros(N)
    ns = zeros(N)
    for i in range(0, N):
        for j in range(i, N):
            pij = cm[i,j]
            if pij > 0:
                ps[j-i] += pij
                ns[j-i] += 1.
    for i in range(0, N):
        if ns[i] > 0:
            ps[i] = ps[i]/ns[i]
        else:
            ps[i] = nan
    return ps

# Print P(s)_{obs} & P(s)_{fit}
def savePs(fn, pss, N, ct):
    fw = open(fn, 'w')
    fw.write(ct)
    for i in range(N):
        lt = "%4d " % (i)
        for m in range(2):
            lt += "%11s "%('NaN') if isnan(pss[m][i]) else "%11.4e"%(pss[m][i])
        fw.write(lt+'\n')
    fw.close()
    return

# Print pij_{obs} & pij_{fit}
def saveMx(fn, mxs, N, ct):
    fw = open(fn, 'w')
    fw.write(ct)
    for i in range(N):
        lt = ''
        for j in range(N):
            if i==j:
                lt += "%11s "%('NaN')
            elif i>j:
                lt += "%11s "%('NaN') if isnan(mxs[0][i,j]) else "%11.4e"%(mxs[0][i,j])
            else:
                lt += "%11s "%('NaN') if isnan(mxs[1][i,j]) else "%11.4e"%(mxs[1][i,j])
        fw.write(lt+'\n')
    fw.close()
    return

# P_{obs}(s) versus P_{fit}(s)
if not len(sys.argv) == 4:
    print('usage:: python tino_ps.py P_{obs} fitDataDir P_{fit}')
    sys.exit()
fhic = str(sys.argv[1])
ddir = str(sys.argv[2])
ffit = str(sys.argv[3])

# P_{obs}(s)
if not os.path.isfile(fhic):
    print('Cannot find '+fhic)
    sys.exit()
else:
    P_obs = []
    with open(fhic) as fr:
        for line in fr:
            if not line[0] == '#':
                lt = line.strip()
                lt = lt.split()
                P_obs.append( list(map(float, lt)) )
    P_obs = array(P_obs)
    N = len(P_obs)
    #
    ps_obs = calPs(P_obs, N)

# P_{fit}(s)
fx = "%s/%s"%(ddir, ffit)
if not os.path.isfile(fx):
    print('Cannot find '+fx)
else:
    P_fit = []
    with open(fx) as fr:
        for line in fr:
            if not line[0] == '#':
                lt = line.strip()
                lt = lt.split()
                P_fit.append( list(map(float, lt)) )
    P_fit = array(P_fit)
    #
    ps_fit = calPs(P_fit, N)
    # Save
    ct = "#N: %d\n#fa: %s\n#fb: %s\n" % (N, fhic, ffit)
    savePs(fx+'.vp', [ps_obs,ps_fit], N, ct)
    saveMx(fx+'.vm', [P_obs,P_fit], N, ct)

