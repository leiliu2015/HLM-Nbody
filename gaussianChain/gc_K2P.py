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
from scipy.linalg import toeplitz
from itertools import permutations
from itertools import combinations
set_printoptions(precision=3, linewidth=200)
import warnings
warnings.filterwarnings('ignore')

# Calculate p_{n} with n>=2
if not len(sys.argv) == 6:
    print('usage:: python gc_K2P.py N fdx[0/1] n[>=2] h5q[0/1] k0')
    sys.exit()
N  = int(sys.argv[1])
fdx= int(sys.argv[2])
ndx= int(sys.argv[3])
h5q= int(sys.argv[4])
k0=float(sys.argv[5])
xi2= 1.0**2 # r^2_c
NL = 2**5 + 1 # Number of integration points between [0, 2*rc=2] for each r_{ij} at fdx==1

# Check avaible choices of {fdx,n,gq}
if fdx:
    if ndx > 3:
        print('Analytical solution of p_n with n>3 and fdx=1 is not available.')
        sys.exit()
else:
    gq = True # Compute the general solution of q_n or not
    if (gq == False) and (ndx > 4):
        print('Analytical solution of q_n with n>4 is only available at gq=True.')
        sys.exit()

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

# 3-body contact probability given a heaviside-step form factor
def K2p3mStep(K, NL, rc):
    xxx= 1.0e-12
    d = sum(K, axis=0)
    D = diag(d)
    km= D - K
    N = len(km)
    G = 0.5/K2GAM(K)
    #
    dl = rc/(NL-1)
    rijs= linspace(0, rc, NL)
    p3m = zeros((N, N, N))*nan
    for i in range(0, N):
        # print("p3m: %3d / %3d" % (i, N))
        for j in range(i+1, N):
            nij = list(range(0, N))
            nij.remove(i)
            nij.remove(j)
            nij = array(nij).astype(int)
            # partial kirhoff matrix [delete i,j]
            pkm = km[nij,:][:,nij]
            pDk = linalg.det(pkm) # integration over x_p, p \neq {i,j}
            # x+sft will give a quadratic compact form x^{T}*A*x of the potential energy at the conditions of x_i=1, x_j=-1
            sft = dot( km[i,nij]-km[j,nij], linalg.inv(pkm) )
            #
            gij= G[i,j]
            Gij = 4.0*(gij**1.5)/sqrt(pi)
            for k in range(j+1, N):
                sfk = sft[ where(nij==k)[0][0] ]
                #
                nijk= list(nij)
                nijk.remove(k)
                nijk= array(nijk).astype(int)
                pdk = linalg.det(km[nijk,:][:,nijk]) # integration over x_p, p \neq {i,j,k}
                Ck  = pDk/pdk/2.
                Gk  = (pi/Ck)**1.5
                # print([pDk, pdk, Ck, Gk, sfk])
                #
                sCk= sqrt(Ck)
                qL = zeros(NL) # to compute p(r_ij<=rc & r_ik<=rc & r_jk<=rc)
                for l in range(0, NL):
                    rij = rijs[l] # given r_{ij}
                    hrij= rij/2.0
                    D   = rc - hrij # z_{k} in [-D,D]
                    sfz = -sfk*hrij # exp(-C*(z-sfz)^2)
                    u   = rij + 2.0*sfz
                    v   = rij - 2.0*sfz
                    ##print [l, hrij, sfz, -D-sfz, D-sfz]
                    #
                    # A: \int_{-D}^{D} exp(-C*(z-sfz)^2)
                    za = (-D-sfz)*sCk
                    zb = ( D-sfz)*sCk
                    if za*zb < 0:
                        A = erf(abs(za)) + erf(abs(zb))
                    else:
                        A = abs( erf(abs(za)) - erf(abs(zb)) )
                    A = A*sqrt(pi)/2.0/sCk
                    #
                    # B: exp[C*(rij^2 - 4sfz^2)/4] / C*(rij+2sfz) * [1 - exp(C*D*(rij+2sfz))]
                    if abs(u) < xxx:
                        B =-D
                    else:
                        B = exp(Ck*u*v/4.0) / Ck / u * (1.0 - exp(Ck*D*u))
                    # C: exp[C*(rij^2 - 4sfz^2)/4] / C*(rij-2sfz) * [1 - exp(C*D*(rij-2sfz))]
                    if abs(v) < xxx:
                        C =-D
                    else:
                        C = exp(Ck*u*v/4.0) / Ck / v * (1.0 - exp(Ck*D*v))
                    # p(r_ik<=rc & r_jk<=rc) at the conditions of {x_i=y_i=x_j=y_j=0, z_i=-z_j=rij/2}
                    qL[l] = ( A+(B+C)*exp(-Ck*rc*rc) )*pi/Ck/Gk
                    # p(r_ij<=rc & r_ik<=rc & r_jk<=rc) = \int_{0}^{rc} dr_{ij} [ p(r_ij) * (r_ik<=rc & r_jk<=rc | r_ij) ]
                    prij = Gij * (rij**2) * exp(-gij*rij*rij)
                    #
                    qL[l]=prij*qL[l]
                #
                pijk = romb(qL, dx=dl)
                for index in permutations([i,j,k]):
                    p3m[index] = pijk
    return p3m

# 3-body contact probability given a Gaussian form factor
def K2p3mGaus(K, rc2):
    d = sum(K, axis=0)
    D = diag(d)
    km= D - K
    N = len(K)
    S = zeros((N, N))
    S[1:,1:] = linalg.inv(km[1:,1:])
    p3m = zeros((N, N, N))*nan
    for i,j,k in combinations(list(range(N)), 3):
        a = S[i,i] + S[j,j] - 2.*S[i,j]
        b = S[i,j] + S[j,k] - S[i,k] - S[j,j]
        c = S[j,j] + S[k,k] - 2.*S[j,k]
        z = a*c - b**2
        p3= (1. + 6.*(a+b+c)/rc2 + 27.*z/(rc2**2))**(-1.5)
        for index in permutations([i,j,k]):
            p3m[index] = p3
    return p3m

# 4-body contact probability given a Gaussian form factor
def K2p4mGaus(K, rc2):
    d = sum(K, axis=0)
    D = diag(d)
    km= D - K
    N = len(K)
    S = zeros((N, N))
    S[1:,1:] = linalg.inv(km[1:,1:])
    p4m = zeros((N, N, N, N))*nan
    LA4 = array([[3,2,1],[2,4,2],[1,2,3]])*3./rc2
    for i,j,k,l in combinations(list(range(N)), 4):
        a = S[i,i] + S[j,j] - 2.*S[i,j]
        b = S[i,j] + S[j,k] - S[i,k] - S[j,j]
        c = S[j,j] + S[k,k] - 2.*S[j,k]
        d = S[j,k] + S[k,l] - S[j,l] - S[k,k]
        e = S[k,k] + S[l,l] - 2.*S[k,l]
        f = S[i,k] + S[j,l] - S[i,l] - S[j,k]
        W4= array([[a,b,f],[b,c,d],[f,d,e]])
        p4= ( linalg.det(W4) * linalg.det(linalg.inv(W4) + LA4) )**(-1.5)
        for index in permutations([i,j,k,l]):
            p4m[index] = p4
    return p4m

# n-body contact probability given a Gaussian form factor
def K2pnmGaus(K, n, rc2):
    n = int(n)
    d = sum(K, axis=0)
    D = diag(d)
    km= D - K
    #
    km= km[1:,1:]
    Dk= linalg.det(km)
    N = len(K)
    pnm = zeros([N for i in range(n)])*nan
    fDTA = 3./rc2
    for A in combinations(list(range(N)), n):
        DTA = zeros((N, N))
        for u in A:
            for v in A:
                if u == v:
                    DTA[u,u] = n-1.
                else:
                    DTA[u,v] = -1.
        mDk = linalg.det(km + fDTA*DTA[1:,1:])
        pn = (Dk/mDk)**(1.5)
        for index in permutations(A):
            pnm[index] = pn
    return pnm

# Save contact probability tensor in a TXT file
def saveTsTx(fn, ts, ct):
    d = len(shape(ts))
    N = shape(ts)[0]
    dv= d - 2 # dimension of viewpoints
    #
    fw= open(fn, 'w')
    fw.write(ct)
    bx= 0
    for index in combinations(list(range(N)), dv):
        mx = ts[index]
        ct = '#vp: '
        for i in index:
            ct += "%d " % (i)
        ct+= "min: %11.5e " % (nanmin(mx))
        ct+= "max: %11.5e " % (nanmax(mx))
        ct+= "block-index: %d " % (bx)
        fw.write(ct+'\n')
        bx += 1
        for i in range(0, N):
            lt = ''
            for j in range(0, N):
                lt += "%11s "%('NaN') if isnan(mx[i,j]) else "%11.5e "%(mx[i,j])
            fw.write(lt+'\n')
        fw.write('\n\n')
    fw.close()

# Save contact probability tensor in a H5 file
def saveTsH5(fn, ts):
    fw = h5py.File(fn, 'w')
    fw.create_dataset('cm', data=ts, dtype='f')
    fw.close()

# Main
if True:
    # Define K-matrix
    ks = zeros(N)
    ks[1] = k0
    K_fit = toeplitz(ks)

    # Compute p_{n}
    if fdx:
        if ndx == 2:
            P = GAM2P(K2GAM(K_fit), fdx)
        else:
            P = K2p3mStep(K_fit, NL, sqrt(xi2))
    else:
        if gq:
            P = K2pnmGaus(K_fit, ndx, xi2)
        else:
            if ndx == 2:
                P = GAM2P(K2GAM(K_fit), fdx)
            elif ndx==3:
                P = K2p3mGaus(K_fit, xi2)
            else:
                P = K2p4mGaus(K_fit, xi2)

    # Save results
    fo = "g%d.%s%d" % (N, 'p' if fdx else 'q', ndx)
    if h5q:
        saveTsH5(fo+'.h5',  P)
    else:
        saveTsTx(fo+'.txt', P, "#N: %d fdx: %d "%(N, fdx)+("NL: %d"%(NL) if fdx else "gq: %d"%(gq))+'\n')

