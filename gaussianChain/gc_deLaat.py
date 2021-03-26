#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import h5py
import networkx as nx
from numpy import *
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

def findMaxCliques(xyz, N, fdx):
    G = nx.Graph()
    G.add_nodes_from(range(N))
    for i in range(0, N):
        for j in range(i+1, N):
            v = xyz[j] - xyz[i]
            r2= dot(v,v)
            if fdx == 0:
                if random.random() < exp(-1.5*r2):
                    G.add_edge(i,j)
            else:
                if r2 < 1:
                    G.add_edge(i,j)
    return list(nx.find_cliques(G))

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

# Try to repeat de Laat's analysis based on 3D chromatin coordinates @Mar25.2020 by LL
if not len(sys.argv) == 5:
    print('usage:: python mc4c_h5association.py xxx.xyz.h5 fdx[0/1] vp_index[0,N) freq')
    sys.exit()
fx = str(sys.argv[1])
fdx= int(sys.argv[2])
vp = [int(sys.argv[3])]
fq = int(sys.argv[4])

rds = 1274; random.seed(rds)
npos_min = 100
nperm = 200
vp_lb = "vp%02d"%(vp[0])
fy = fx[:-3]+".deLaat-F%d"%(fdx)

opq = True # Overall profile
vsq = True # Three-body viewpoint association

if not os.path.isfile(fx):
    print('Cannot find '+fx)
    sys.exit()
else:
    cpx2= [] # Valid concatemers/complexes/reads in the units of bin with k>=2
    cpx3= [] # Valid concatemers/complexes/reads in the units of bin with k>=3
    fr = h5py.File(fx, 'r')
    NF, nb = fr.get('xyz').shape[:2]
    for f in range(0, NF, fq):
        R = fr.get('xyz')[f]
        maxCliques = findMaxCliques(R, nb, fdx)
        if False:
            for read in maxCliques:
                if sum(isin(read, vp)):
                    if len(read) >= 2:
                        cpx2.append(read)
                    if len(read) >= 3:
                        cpx3.append(read)
        else:
            for maxClique in maxCliques:
                if sum(isin(maxClique, vp)):
                    read = []
                    for i in maxClique:
                        if not i in vp:
                            read.append(i)
                    if len(read) >= 1:
                        cpx2.append(maxClique)
                    if len(read) >= 2:
                        cpx3.append(maxClique)
    fr.close()

    if opq:
        # Overall profile
        op2 = zeros(nb)
        for read in cpx2:
            op2[read] += 1.0
        op2[vp] = nan
        ncpx2 = len(cpx2)*1.0
        #
        op3 = zeros(nb)
        for read in cpx3:
            op3[read] += 1.0
        op3[vp] = nan
        ncpx3 = len(cpx3)*1.0
        #
        fw = open(fy+".%s.overallProfile"%(vp_lb), 'w')
        ct = "#MQ>=NaN totalReads(k>=2): %d totalReads(k>=3): %d\n" % (ncpx2, ncpx3)
        ct+= '#i gs[bin] ge[bin] reads(k>=2) fraction(k>=2) read(k>=3) fraction(k>=3)\n'
        fw.write(ct)
        for i in range(nb):
            lt = "%d %d %d " % (i, i, i)
            lt+= 'NaN ' if isnan(op2[i]) else "%d "%(op2[i])
            lt+= 'NaN ' if isnan(op2[i]) else "%11.5e "%(op2[i]/ncpx2)
            lt+= 'NaN ' if isnan(op3[i]) else "%d "%(op3[i])
            lt+= 'NaN ' if isnan(op3[i]) else "%11.5e "%(op3[i]/ncpx3)
            fw.write(lt+'\n')
        fw.close()
        #
        siz= [len(read) for read in cpx2] # Size of reads in units of bin
        sic= bincount(siz)
        fw = open(fy+".%s.readSize"%(vp_lb), 'w')
        for i in range(2, len(sic)):
            fw.write("%d %d\n" % (i, sic[i]))
        fw.close()

    if vsq:
        # Association analysis
        ncpx3 = len(cpx3)*1.0
        frq_pos = zeros((nb, nb)) # Set SOI along the second axis
        frq_neg = zeros((nb, nb, 2)) # average/std
        N_POS  = zeros(nb) # Number of reads in positive set for each SOI
        for soi in range(nb):
            if not soi in vp:
                qs = array([1 if soi in read else 0 for read in cpx3])
                q_pos = where(qs==1)[0]
                q_neg = where(qs==0)[0]
                n_pos = len(q_pos)
                # positive set
                for q in q_pos:
                    frq_pos[cpx3[q], soi] += 1.0
                # negative set
                frq_neg_perm = zeros((nb, nperm))
                for p in range(nperm):
                    q_neg_perm = random.permutation(q_neg)
                    for q in q_neg_perm[:n_pos]: # down-sampling
                        frq_neg_perm[ random.permutation(cpx3[q])[1:], p] += 1.0 # fragment number compensation
                for i in range(nb):
                    frq_neg[i, soi, 0] = mean(frq_neg_perm[i])
                    frq_neg[i, soi, 1] =  std(frq_neg_perm[i])
                #
                for i in vp+[soi]:
                    frq_pos[i, soi] = nan
                    frq_neg[i, soi, :] = nan
                N_POS[soi] = n_pos
            else:
                frq_pos[:, soi] = nan
                frq_neg[:, soi, :] = nan
                N_POS[soi] = nan
            print("z: %3d / %3d" % (soi, nb))
        #
        zm = zeros((nb, nb))*nan
        for soi in range(nb):
            if (not soi in vp) and (N_POS[soi] >= npos_min):
                for i in range(0, nb):
                    z = (frq_pos[i,soi] - frq_neg[i,soi,0])/frq_neg[i,soi,1]
                    if not isnan(z):
                        zm[i,soi] = z
        saveMx(fy+".%s.zm"%(vp_lb), zm, "#k>=3 MQ>=NaN n_pos>=%d\n#i-soi-z\n"%(npos_min))
        #
        fw = open(fy+".%s.vpsoiProfile"%(vp_lb), 'w')
        fw.write("#k>=3 MQ>=NaN totalReads: %d\n"%(ncpx3))
        for soi in range(nb):
            n_pos = N_POS[soi]
            ct = "#soi: %d %d-%d " % (soi, soi, soi)
            ct+= 'n_pos: NaN' if isnan(n_pos) else "n_pos: %d"%(n_pos)
            ct+= '\n#i pos neg_ave neg_std pos/n_pos*100.0 neg_ave/n_pos*100.0 neg_std/n_pos*100.0 z-score\n'
            fw.write(ct)
            for i in range(nb):
                fp = frq_pos[i,soi]
                fnm= frq_neg[i,soi,0]
                fns= frq_neg[i,soi,1]
                z  = zm[i,soi]
                #
                lt = "%d " % (i)
                lt+= 'NaN ' if isnan(fp) else "%d "%(fp)
                lt+= "%11s "%('NaN') if isnan(fnm) else "%11.5e "%(fnm)
                lt+= "%11s "%('NaN') if isnan(fns) else "%11.5e "%(fns)
                #
                lt+= "%11s "%('NaN') if isnan(fp) else "%11.5e "%(fp*100.0/n_pos)
                lt+= "%11s "%('NaN') if isnan(fnm) else "%11.5e "%(fnm*100.0/n_pos)
                lt+= "%11s "%('NaN') if isnan(fns) else "%11.5e "%(fns*100.0/n_pos)
                #
                lt+= "%11s "%('NaN') if isnan(z) else "%11.5e "%(z)
                fw.write(lt+'\n')
            fw.write('\n\n')
        fw.close()
