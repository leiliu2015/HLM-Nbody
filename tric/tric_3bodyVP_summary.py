#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import h5py
from numpy import *
from numba import jit
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# Asscending order
def rankArray(x):
    ai = argsort(x)
    rk = empty_like(ai)
    rk[ai] = arange(len(x))
    return rk

# N_k, r_{2k}, and rho_{k}, which are required to compute rho_s by using Eq.(7)
def calSCC(cm, N, smx, lnq, r2q):
    #
    nk = zeros(smx)
    rok= zeros(smx)
    r2k= zeros(smx)
    cs = zeros(smx)*nan
    #
    for s in range(1, smx):
        x = []
        y = []
        for i in range(0, N-s):
            j = i+s
            cas = cm[0,i,j]
            cbs = cm[1,i,j]
            if (cas>0) and (cbs>0): # excluding NaN or either c_{ij}==0, different from HiCRep
                if lnq:
                    x.append(log10(cas))
                    y.append(log10(cbs))
                else:
                    x.append(cas)
                    y.append(cbs)
        x = array(x)
        y = array(y)
#        #
#        if var(x) == 0:
#            x[0] *= (1.0+1.0e-9)
#        if var(y) == 0:
#            y[0] *= (1.0+1.0e-9)
        #
        n = len(x)
        if n>1:
            if (var(x)==0) or (var(y)==0):
                pass
            else:
                cs[s] = pearsonr(x,y)[0]
        if n>0:
            nk[s] = n*1.0
            #
            if (var(x)>0) and (var(y)>0):
                rok[s] = corrcoef(x,y)[0,1]
            else:
                rok[s] = 0.0
            #
            if r2q:
                xi = rankArray(x)/float(n)
                yi = rankArray(y)/float(n)
                r2k[s] = sqrt(var(xi)*var(yi))
            else:
                r2k[s] = sqrt(var(x)*var(y))
    #
    A = 0.0 # sum_{k} N_{k}*r_{2k}
    B = 0.0 # sum_{k} N_{k}*r_{2k}*rho_{k}
    C = 0.0 # sum_{k} (N_{k}*r_{2k})^2 * (1-rho_{k}^2)^2 / N_{k}-3
    #
    for s in range(1, smx-3):
        if nk[s] > 3:
            nr = nk[s]*r2k[s]
            A += nr
            B += nr*rok[s]
            C += (nr**2)*((1.0 - rok[s]**2)**2)/(nk[s] - 3.0)
    #
    rhos = B/A # Eq.(7)
    rhov = C/(A**2) # Eq.(9)
    return [rhos, rhov, cs]

# pearson & spearman
def calCorr(cm, N, smx, lnq):
    x = []
    y = []
    for i in range(0, N):
        for j in range(i+1, N):
            s = j-i
            if s < smx:
                cas = cm[0,i,j]
                cbs = cm[1,i,j]
                if (cas>0) and (cbs>0):
                    if lnq:
                        x.append(log10(cas))
                        y.append(log10(cbs))
                    else:
                        x.append(cas)
                        y.append(cbs)
    x = array(x)
    y = array(y)
    #
    pc, pp = pearsonr(x,y)
    sc, sp = spearmanr(x,y)
    sx = sum(y*x)/sum(x*x) # best fit: y=sx*x
    #
    if pp < 1.0e-32: pp = 0.0
    if sp < 1.0e-32: sp = 0.0
    return [pc,pp,sc,sp,sx]

# compare Two contact probability matrices
# by usint the Stratum-adjusted Correlation Coefficient (SCC, HiCRep, Genome Res. 2017). @Feb.01.2020 by LL
#
if not len(sys.argv) >= 3:
    print('usage:: python tric_3bodyVP_summary.py q3_vp.tric q3_vp.hlm [lnq=0] [dq=0:deleteLowTrustRegionAtHba]')
    sys.exit()
#
fx = [str(sys.argv[1]), str(sys.argv[2])]
fo = fx[-1][:-4]+'.scc'
#
fq = [0, 0] # 0/1: txt:h5
lnq= 0 if len(sys.argv)==3 else int(sys.argv[3]) # 0/1: c_{ij}/log10(c_{ij})
dq = 0 if len(sys.argv)==4 else int(sys.argv[4]) # Hba1 and Hba2 are indistinguishable bioinformatically, which results in artifitial contacts around Hba 
r2q= 1 # 0/1: Eq.(5)/Eq.(10)

#Hba-a1	91	91	P
#Hbq1b	93	93	P
#Hba-a2	97	98	P
#Hbq1a	99	100	P
bs_del = 91; be_del=100

# readin
for f in fx:
    if not os.path.isfile(f):
        print('cannot find '+f)
        sys.exit()
#
cm = []
if sum(fq) == 0:
    for i in range(0, 2):
        cm.append([])
        fr = open(fx[i], 'r')
        for line in fr.readlines():
            if not line[0] == '#':
                lt = line.strip()
                lt = lt.split()
                cm[-1].append( list(map(float, lt)) )
        fr.close()
    #
    if not len(cm[0]) == len(cm[1]):
        print('check matrix dimensions')
        sys.exit()
    else:
        cm = array(cm)
        cm_bk = cm*1.
        N  = len(cm[0])
        smx= N # 5Mb/resolution, recommended by HiCRep.
        if dq:
            for k in range(2):
                cm[k, bs_del:(be_del+1), bs_del:(be_del+1)] = nan
        # print([nanmax(cm[0]), nanmax(cm[1])])
else:
    sys.exit()

# correlations
scc, scv, scs = calSCC(cm, N, smx, lnq, r2q)
pc,pp,sc,sp,sx = calCorr(cm, N, smx, lnq)

# save
ct = "#fa: %s\n#fb: %s\n" % (fx[0], fx[1])
ct+= "#smx: %4d lnq: %d r2q: %d   " % (smx, lnq, r2q)
ct+= "scc: %+12.5e %11.5e   " % (scc,scv)
ct+= "pearson: %+12.5e %11.5e   " % (pc, pp)
ct+= "spearman: %+12.5e %11.5e   " % (sc, sp)
#
fw = open(fo, 'w')
fw.write(ct+'\n')
for s in range(1, smx):
    if isnan(scs[s]):
        fw.write("%4d %12s\n" % (s, 'NaN'))
    else:
        fw.write("%4d %+12.5e\n" % (s, scs[s]))
fw.close()

ct = "#fa: %s\n#fb: %s\n" % (fx[0], fx[1])
ct+= "#fa/fb: %+12.5e deleteQ[%d-%d]: %d" % (sx, bs_del, be_del, dq)
fw = open(fo[:-4]+'.vm', 'w')
fw.write(ct+'\n')
for i in range(0, N):
    lt = ''
    for j in range(0, N):
        if i==j:
            lt += "%11s "%('NaN')
        elif i>j:
            lt += "%11s "%('NaN') if isnan(cm_bk[0][i,j]) else "%11.4e"%(cm_bk[0][i,j]*sx)
        else:
            lt += "%11s "%('NaN') if isnan(cm_bk[1][i,j]) else "%11.4e"%(cm_bk[1][i,j])
    fw.write(lt+'\n')
fw.close()

