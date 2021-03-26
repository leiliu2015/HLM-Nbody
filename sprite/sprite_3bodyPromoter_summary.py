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

# compare Two contact probability matrices via {scc, pearson, spearman} @Jul.05.2020 by LL
if not len(sys.argv) >= 3:
    print('usage:: python sprite_3bodyPromoter_summary.py TAD-id[0-38] xxx.p3.h5 [lnq=0]')
    sys.exit()
tdx = int(sys.argv[1])
fh = str(sys.argv[2])
fo = fh[:-3]+'.summary'
#
lnq= 0 if len(sys.argv)==3 else int(sys.argv[3]) # 0/1: c_{ij}/log10(c_{ij})
r2q= 1 # 0/1: Eq.(5)/Eq.(10)

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
                    smx = N # 5Mb/resolution, recommended by HiCRep.
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

fw = open(fo, 'w')
apx = ['pm', 'em']
for b in range(0, PBN):
    ps = where(pbx==b)[0]
    bs = min(pdx[ps])
    be = max(pdx[ps])
    ct = "TAD%02d vp: %03d-%03d smx: %4d lnq: %d r2q: %d   " % (tdx, bs, be, smx, lnq, r2q)
    #
    for u in range(0, 2):
        ct+= "   [%s] " % (apx[u])
        fx = ["sprite-data/TAD%02d.p3.vp%03d-%03d.txt"%(tdx, bs, be), \
              fh[:-3]+".vp%03d-%03d.%s"%(bs, be, apx[u])]
        # Read p_{ijk}
        for f in fx:
            if not os.path.isfile(f):
                print('cannot find '+f)
                sys.exit()
        #
        cm = []
        if True:
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
        #
        scc, scv, scs = calSCC(cm, N, smx, lnq, r2q)
        pc,pp,sc,sp,sx = calCorr(cm, N, smx, lnq)
        #
        ct+= "scc: %+12.5e %11.5e " % (scc,scv)
        ct+= "pearson: %+12.5e %11.5e " % (pc, pp)
        ct+= "spearman: %+12.5e %11.5e " % (sc, sp)
        #
        if u == 0:
            p3_spt = []
            p3_hlm = []
            for i in range(0, N):
                for j in range(i+1, N):
                    if not isnan(cm[0,i,j]+cm[1,i,j]):
                        p3_spt.append( cm[0,i,j] )
                        p3_hlm.append( cm[1,i,j] )
            p3_spt = array(p3_spt)
            p3_hlm = array(p3_hlm)
            #
            Q = 4
            M = int(len(p3_spt)/float(Q))
            idx = argsort(p3_hlm)
            #
            fc = open(fh[:-3]+".vp%03d-%03d.%s-qt"%(bs, be, apx[u]), 'w')
            for q in range(0, Q):
                p3 = p3_hlm[ idx[q*M:min(q*M+M,len(p3_spt))] ]
                p3_m = mean(p3)
                p3_s =  std(p3)
                p3_q = [quantile(p3, i*0.25) for i in range(0,5)]
                lt = "q: %d   hlm: %11.5e %11.5e qt: " % (q, p3_m, p3_s)
                for i in range(5):
                    lt += "%11.5e " % (p3_q[i])
                #
                p3 = p3_spt[ idx[q*M:min(q*M+M,len(p3_spt))] ]
                p3_m = mean(p3)
                p3_s =  std(p3)
                p3_q = [quantile(p3, i*0.25) for i in range(0,5)]
                lt += "  sprite: %11.5e %11.5e qt: " % (p3_m, p3_s)
                for i in range(5):
                    lt += "%11.5e " % (p3_q[i])
                fc.write(lt+'\n')
            fc.close()
    fw.write(ct+'\n')
fw.close()

