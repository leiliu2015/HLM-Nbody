#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import h5py
from numpy import *
from scipy.special import erf
from scipy.integrate import romb
set_printoptions(precision=3, linewidth=200)

# Calculate the 'cooperativity', namely p_{jk|ij}/p_{jk} and p_{jk|!ij}/p_{jk} at i<j<k,
# which has benn measured by X. Zhuang in Science 362, 419 (2018).
if not len(sys.argv) == 4:
    print('usage:: python zhuang2018_K2CP.py xxx.K_fit fdx[0/1] CTCF-profile')
    sys.exit()
fk = str(sys.argv[1])
fdx= int(sys.argv[2])
fb = str(sys.argv[3])
xi2= 1.0**2 # r^2_c
NL = 2**5 + 1 # Number of integration points between [0, 2*rc=2] for each r_{ij} at fdx==1

# Onami's method. Note that here GAM = 1/(2*\gamma) = <r^2_{ij}>,
# which is different from \gamma used in Hyeon's paper in 2019.
def K2GAM(K):
    # K to Laplacian matrix
    d = sum(K, axis=0)
    D = diag(d)
    L = D - K
    # Eigenvalues and eigenvectors
    lam, Q = linalg.eigh(L)
    inv_lam = 1 / lam   # inverse of the eigenvalues
    inv_lam[0] = 0
    inv_Lam = diag(inv_lam)
    # L to M
    M = dot(Q, dot(inv_Lam, Q.T))
    # M to Σ^2
    M_diag = diag(diag(M))
    Ones = ones(shape(K))
    A = dot(M_diag, Ones)
    GAM = A + A.T - 2 * M
    return GAM

# 2-body contact probability
def GAM2P(GAM, fdx):
    if fdx == 0:
        P = (1.+3.*GAM/xi2)**(-1.5)
    else:
        imx = eye(N)
        g = (0.5*xi2)* (1./(GAM + imx) - imx) # To avoid division by 0 at the diagonal
        P = erf(sqrt(g)) - (2.0/sqrt(pi))*sqrt(g)*exp(-g) + imx
    P += diag(zeros(N)*nan)
    return P

# p_{ij,jk}, with a order of indices i<j<k, given a gaussian form factor
def K2p2mGaus(K, rc2):
    d = sum(K, axis=0)
    D = diag(d)
    km= D - K
    N = len(K)
    S = zeros((N, N))
    S[1:,1:] = linalg.inv(km[1:,1:])
    p2m = zeros((N, N, N))*nan
    for i in range(0, N):
        for j in range(i+1, N):
            for k in range(j+1, N):
                a = S[i,i] + S[j,j] - 2.*S[i,j]
                b = S[i,j] + S[j,k] - S[i,k] - S[j,j]
                c = S[j,j] + S[k,k] - 2.*S[j,k]
                z = a*c - b**2
                p2m[i,j,k] = (1. + 3.*(a+c)/rc2 + 9.*z/(rc2**2))**(-1.5)
    return p2m

# p_{ij,jk}, with a order of indices i<j<k, given a step form factor
def K2p2mStep(K, NL, rc):
    xxx= 1.0e-12
    d = sum(K, axis=0)
    D = diag(d)
    km= D - K
    N = len(km)
    G = 0.5/K2GAM(K)
    #
    dl = rc/(NL-1)
    rijs=linspace(0, rc, NL)
    p2m = zeros((N, N, N))*nan
    #
    for i in range(0, N):
        for j in range(i+1, N):
            nij = list(range(0, N))
            nij.remove(i)
            nij.remove(j)
            nij = array(nij).astype(int)
            # partial kirhoff matrix [delete i,j]
            pkm = km[nij,:][:,nij]
            pDk = linalg.det(pkm) # integration over x_p, p \neq {i,j}
            # x+sft will give a quadratic compact form x^{T}*A*x of the potential energy at the conditions of x_i=1, x_j=0
            sft = dot( km[i,nij], linalg.inv(pkm) )
            #
            gij= G[i,j]
            Gij = 4.0*(gij**1.5)/sqrt(pi)
            #print km
            #
            for k in range(j+1, N): # i<j<k
                if (not i==k) and (not j==k):
                    sfk = sft[ where(nij==k)[0][0] ]
                    #
                    nijk= list(nij)
                    nijk.remove(k)
                    nijk= array(nijk).astype(int)
                    pdk = linalg.det(km[nijk,:][:,nijk]) # integration over x_p, p \neq {i,j,k}
                    Ck  = pDk/pdk/2.
                    Gk  = (pi/Ck)**1.5
                    #print [pDk, pdk, Ck, Gk, sft]
                    #
                    sCk= sqrt(Ck)
                    qL = zeros(NL) # to compute p(r_ij<=rc & r_jk<=rc)
                    for l in range(0, NL):
                        rij = rijs[l] # given r_{ij}
                        sfz = -sfk*rij # exp(-C*(z-sfz)^2)
                        #
                        # A: \int_{-rc}^{rc} exp(-C*(z-sfz)^2)
                        za = (-rc-sfz)*sCk
                        zb = ( rc-sfz)*sCk
                        if za*zb < 0:
                            A = erf(abs(za)) + erf(abs(zb))
                        else:
                            A = abs( erf(abs(za)) - erf(abs(zb)) )
                        A = A*sqrt(pi)/2.0/sCk
                        #
                        # B: exp[-C(sfz^2 + rc^2)] * sinh(2*C*sfz*rc) / (C*sfz)
                        if abs(sfz) < xxx:
                            B = exp(-Ck*(sfz**2 + rc**2)) * 2.0*rc
                        else:
                            B = exp(-Ck*(sfz**2 + rc**2)) * sinh(2.0*Ck*sfz*rc) / (Ck*sfz)
                        # p(r_jk<=rc) at the conditions of {x_i=y_i=x_j=y_j=0, z_i=rij, z_j=0}
                        qL[l] = ( A-B )*pi/Ck/Gk
                        # p(r_ij<=rc & r_jk<=rc) = \int_{0}^{rc} dr_{ij} [ p(r_ij) * (r_jk<=rc | r_ij) ]
                        prij = Gij * (rij**2) * exp(-gij*rij*rij)
                        #
                        qL[l]=prij*qL[l]
                    #
                    pijk = romb(qL, dx=dl)
                    p2m[i,j,k] = pijk
    return p2m

# Main
if True:
    # Read K-matrix
    if not os.path.isfile(fk):
        print('Cannot find '+fk)
        sys.exit()
    else:
        K_fit = []
        with open(fk) as fr:
            for line in fr:
                if not line[0] == '#':
                    lt = line.strip()
                    lt = lt.split()
                    K_fit.append( list(map(float, lt)) )
        K_fit = array(K_fit)
        N = len(K_fit)

    # Read CTCF binding profile
    if not os.path.isfile(fb):
        print('Cannot find '+fb)
        sys.exit()
    else:
        bf = []
        with open(fb) as fr:
            for line in fr:
                if not line[0] == '#':
                    lt = line.strip()
                    lt = lt.split()
                    bf.append( 1 if float(lt[3])>0 else 0 )
        bf = array(bf)
        if not len(bf) == N:
            print('Check N')
            sys.exit()

    # Compute p_{ij}
    P = GAM2P(K2GAM(K_fit), fdx)
    # Compute p_{ij,jk}
    if fdx:
        CP = K2p2mStep(K_fit, NL, sqrt(xi2))
    else:
        CP = K2p2mGaus(K_fit, xi2)

    # Save results in TXT
    res = []
    for i in range(0, N):
        for j in range(i+1, N):
            for k in range(j+1, N):
                pijjk = CP[i,j,k]
                pij = P[i,j]
                pjk = P[j,k]
                cp1 = pijjk/pij
                cp0 = (pjk - pijjk)/(1.-pij)
                res.append([i,j,k,pijjk,pij,pjk,cp1,cp0])
    res = array(res)
    idx = argsort(res[:,5]) # sort by p_{jk}

    fy = fk+'.zhuang2018'
    fw = open(fy+'.txt', 'w')
    ct = "#N: %d fdx: %d NL: %d\n"%(N, fdx, NL)
    ct+= '#i j k p_{ij,jk} p_{ij} p_{jk} p_{jk|ij} p_{jk|!ij} ctcf-triplet-index'
    fw.write(ct+'\n')
    cti= 0
    for t in idx:
        i = int(res[t,0])
        j = int(res[t,1])
        k = int(res[t,2])
        lt= "%3d %3d %3d   " % (i,j,k)
        for p in res[t,3:]:
            lt += "%11.5e "%(p)
        #
        if (bf[i]+bf[j]+bf[k])==3:
            lt += "%d" % (cti)
            cti+= 1
        fw.write(lt+'\n')
    fw.close()
