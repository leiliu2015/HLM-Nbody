#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import time
import warnings
from numpy import *
from scipy.special import erf
from scipy.stats import pearsonr
set_printoptions(precision=3, linewidth=200)
warnings.filterwarnings('ignore')

# Compare different minimization methods
if not len(sys.argv) == 6:
    print("usage:: python tino_nan.py fdx[0/1] ldx[0/1/2/3] mdx[0/1/2/3/4] rdx[>=0] normalized-HiC-Contact-Matrix")
    sys.exit()
# Options
# Form factor 0/1: F(x) = exp(-3*x^2 / 2*r^2_c) or F(x) = StepHeaviside(r_c - x)
fdx = int(sys.argv[1])
# Cost function:
# 0: {\sum_{ij} (p_{ij} - p^{HiC}_{ij})^2}
# 1: {\sum_{ij} [log10(p_{ij}) - log10(p^{HiC}_{ij})]^2}
# 2: 1 - pearsonr(log10(p_{ij}), log10(p^{HiC}_{ij}))
# 3: \sum_{ij} [<r^2_{ij}>/<r^{2,HiC}_{ij}> - 1]^2
# where i<j, and p_{ij} is the pairwise contact probability between the i-th and j-th monomers.
ldx = int(sys.argv[2])
# Minimization method:
# 0: steepest descent single variable update
# 1: stochastic single variable update (Onami's choice)
# 2: gradient descent (Orland's choice)
# 3: RMS-Prop
# 4: ADAM
mdx = int(sys.argv[3])
# Random seed index
rdx = int(sys.argv[4])
# HiC observation
fhic= str(sys.argv[5])

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

# Form factor
def GAM2P(GAM, fdx):
    if fdx == 0:
        P = (1.+3.*GAM/xi2)**(-1.5)
    else:
        imx = eye(N)
        g = (0.5*xi2)* (1./(GAM + imx) - imx) # To avoid division by 0 at the diagonal
        P = erf(sqrt(g)) - (2.0/sqrt(pi))*sqrt(g)*exp(-g) + imx
    return P

# Pearsonr(GAM)
def GAM2pearson(GAM, fdx):
    lnP = log10( GAM2P(GAM, fdx) )
    return pearsonr(lnP[trilMask], lnP_obs[trilMask])[0]

# Cost function L(GAM)
def GAM2cost(GAM, fdx, ldx):
    if ldx == 0:
        P = GAM2P(GAM, fdx)
        lv= nansum((P - P_obs)**2)/2.
    elif ldx==1:
        P = GAM2P(GAM, fdx)
        lv= nansum((log10(P) - lnP_obs)**2)/2.
    elif ldx==2:
        lnP = log10(GAM2P(GAM, fdx))
        if True:
            sF1 = sum(nanMask* lnP)/2.
            sF2 = sum(nanMask* lnP**2)/2.
            sFO = nansum(lnP*lnP_obs)/2.
            lv= 1. - (npc*sFO - sF1*sO1)/sqrt(npc*sF2 - sF1**2)/vOn
        else:
            lv= 1. - pearsonr(lnP[trilMask], lnP_obs[trilMask])[0]
    elif ldx==3:
        lv = nansum(WS* (GAM - GAM_obs)**2 )/2.
    return lv

# Calculate dL(GAM)/dGAM, called only by the function K2dJ()
def dcost_dGAM(GAM, fdx, ldx):
    if ldx == 0:
        P = GAM2P(GAM, fdx)
        if fdx == 0:
            dF = -4.5/xi2 * (1.+3.*GAM/xi2)**(-2.5)
        else:
            invG = 1./(GAM+eye(N)) - eye(N)
            dF = - xi2**(1.5)/sqrt(2.*pi) * exp(-0.5*xi2*invG) * (invG**(2.5))
        Z = (P - P_obs)*dF
    elif ldx==1:
        P = GAM2P(GAM, fdx)
        lnP = log10(P)
        if fdx == 0:
            dF = -4.5/log(10)/(3.*GAM+xi2)
        else:
            invG = 1./(GAM+eye(N)) - eye(N)
            dF = - xi2**(1.5)/sqrt(2.*pi) * exp(-0.5*xi2*invG) * (invG**(2.5))
            dF = dF/log(10)/P
        Z = (lnP - lnP_obs)*dF
    elif ldx==2:
        P = GAM2P(GAM, fdx)
        lnP = log10(P)
        sF1 = sum(nanMask* lnP)/2.
        sF2 = sum(nanMask* lnP**2)/2.
        sFO = nansum(lnP*lnP_obs)/2.
        vFn = sqrt(npc*sF2 - sF1**2)
        lva = (npc*sFO - sF1*sO1)/(vFn**2)
        lvb = 1./(2.*vFn*vOn)
        lvc = lvb*(sO1 - lva*sF1)
        if fdx == 0:
            dF = -4.5/log(10)/(3.*GAM+xi2)
        else:
            invG = 1./(GAM+eye(N)) - eye(N)
            dF = - xi2**(1.5)/sqrt(2.*pi) * exp(-0.5*xi2*invG) * (invG**(2.5))
            dF = dF/log(10)/P
        Z = -(lvb*npc*(lnP_obs - lva*lnP) - lvc)*dF
    elif ldx==3:
        Z = WS*(GAM - GAM_obs)
    return Z

# Calculate dL/dK with matrix calculus, based on Orland's derivation.
def K2dJ(K, fdx, ldx):
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
    #
    Z = dcost_dGAM(GAM, fdx, ldx)
    Z[nanMaskN] = 0
    #
    Y = 2.*(diag(sum(Z, axis=0)) - Z)
    BET = dot(dot(M, Y), M)
    #
    BET_diag = diag(diag(BET))
    H = dot(BET_diag, Ones)
    dJ= -0.5*(H + H.T - 2*BET)
    return dJ

# A subroutine only called by the function zbrent()
def cpij(gij, rc):
    cp = -2.0*rc*sqrt(gij/pi)*exp(-gij*rc*rc)
    cp+= erf(sqrt(gij)*rc)
    return cp

# Calculate gamma_{ij} based on p_{ij} and r_c
# in the case of a Heaviside Step form factor, only called by the function P2GAM()
def zbrent(rc, cp, x1, x2):
    #default
    ITMAX = 500
    EPS = 1.0e-16 # machine float-point precision
    tol = 1.0e-16 # convergent criteria
    #
    a = x1
    b = x2
    c = x2
    #
    fa= cpij(a, rc) - cp
    fb= cpij(b, rc) - cp
    #
    if (fa>0 and fb>0) or (fa<0 and fb<0):
        print('Root must be bracketed in zbrent')
        sys.exit()
    fc= fb
    #
    for iter in range(0, ITMAX):
        if (fb>0 and fc>0) or (fb<0 and fc<0):
            c = a
            fc= fa
            d = b-a
            e = b-a
        if (abs(fc) < abs(fb)):
            a = b
            b = c
            c = a
            fa= fb
            fb= fc
            fc= fa
        tol1 = 2.0*EPS*abs(b) + 0.5*tol
        xm= (c-b)/2.0
        if (abs(xm)<=tol1 or fb==0):
            return 0.5/b # GAM = 0.5/\gamma
        #
        if (abs(e)>=tol1 and abs(fa)>abs(fb)):
            s = fb/fa
            if (a==c):
                p = 2.0*xm*s
                q = 1.0-s
            else:
                q = fa/fc
                r = fb/fc
                p = s*(2.0*xm*q*(q-r) - (b-a)*(r-1.0))
                q = (q-1.0)*(r-1.0)*(s-1.0)
            #
            if (p>0):
                q = -q
            p = abs(p)
            min1 = 3.0*xm*q - abs(tol1*q)
            min2 = abs(e*q)
            #
            if (2.0*p < min(min1, min2)):
                e = d
                d = p/q
            else:
                d = xm
                e = d
        else:
            d = xm
            e = d
        #
        a = b
        fa= fb
        if (abs(d) > tol1):
            b += d
        else:
            # nrutil.h - SIGN(a,b)
            if xm >= 0:
                b += abs(tol1)
            else:
                b -= abs(tol1)
        #
        fb = cpij(b, rc) - cp
    #
    print('Maximum double of iterations exceeded in zbrent\n')
    sys.exit()

# Calculate GAM from p_ij
def P2GAM(P, fdx):
    if fdx == 0:
        GAM = xi2/3.0*(P**(-2.0/3)-1.)
    else:
        GAM = zeros((N, N))
        for i in range(0, N):
            for j in range(i+1, N):
                if isnan(P[i,j]):
                    GAM[i,j] = GAM[j,i] = nan
                else:
                    GAM[i,j] = GAM[j,i] = zbrent(sqrt(xi2), P[i,j], 0, 100)
            if isnan(P[i,i]):
                GAM[i,i] = nan
    return GAM

# Initialize k_ij
def Init_K(N, INIT_K0):
    K = zeros((N, N))
    for i in range(1, N):
        j = i-1
        K[i,j] = K[j,i] = INIT_K0
    return K

# Choose {i,i+1} randomly
def Random_Int_Pair_Adjacent(N):
    i = random.randint(N-1)
    j = i + 1
    return i, j

# Choose {i,j} with j>i+1 randomly
def Random_Int_Pair_Nonadj(N, nps):
    k = random.randint(nps)
    v = sqrt((2*N-3)**2 - 8.*k)/2.
    i = int(N-1.5-v)
    j = int(k-(2*N-i-3)*i/2+i+2)
    return i, j

# Choose m pairs {i,i+1} randomly
def Random_Int_Pairs_Adjacent(N, m):
    k_mask = zeros((N, N))
    for i in random.choice(N-1, size=m, replace=False):
        k_mask[i,i+1] = k_mask[i+1,i] = 1
    return k_mask

# Choose m pairs {i,j} with j>i+1 randomly
def Random_Int_Pairs_Nonadj(N, nps, m):
    k_mask = zeros((N, N))
    for k in random.choice(nps, size=m, replace=False):
        v = sqrt((2*N-3)**2 - 8.*k)/2.
        i = int(N-1.5-v)
        j = int(k-(2*N-i-3)*i/2+i+2)
        k_mask[i,j] = k_mask[j,i] = 1
    return k_mask

# Save an array of size n*m
def saveLg(fn, xy, ct):
    fw = open(fn, 'w')
    fw.write(ct)
    n = len(xy)
    m = shape(xy)[1] if xy.ndim==2 else 0
    for i in range(n):
        if m == 0:
            lt = "%11s "%('NaN') if isnan(xy[i]) else "%11.5e"%(xy[i])
        else:
            lt = ''
            for v in xy[i]:
                lt += "%11s "%('NaN') if isnan(v) else "%11.5e "%(v)
        fw.write(lt+'\n')
    fw.close()

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

# Steepest descent single variable update
def sd_svu(KZ, ITERATIONS=100, STEP=25, EPS=1.0e-7):
    beta = 0.9
    paras_fit = "%4d\t%4d\t%f\t%f" % (ITERATIONS, STEP, EPS, beta)
    # Preparation
    K = KZ
    G = K2GAM(K)
    C = GAM2cost(G, fdx, ldx)
    c_traj = zeros((ITERATIONS+1, 3)); c_traj[0] = C, GAM2pearson(G, fdx), time.time(); # {cost, 1-pearsonr, time}
    # Minimization
    for iteration in range(ITERATIONS):
        s = 0
        for step in range(STEP):
            dJ = K2dJ(K, fdx, ldx)
            i, j = unravel_index(argmax(abs(dJ)), (N,N))
            tmp_kij = K[i,j]
            K[i,j] = K[i,j] - dJ[i,j]*EPS
            K[j,i] = K[j,i] - dJ[i,j]*EPS
            tmp_C = GAM2cost(K2GAM(K), fdx, ldx)
            if tmp_C < C:
                C = tmp_C
                s+= 1
                if spq==2: print("1\t%d\tC\t%04d\t%f\tK[%4d,%4d]\t%f" % (iteration, step, C, i, j, abs(dJ[i,j])))
            else:
                K[i,j] = tmp_kij
                K[j,i] = tmp_kij
                EPS *= beta
                if spq==2: print("0\t%d\tC\t%04d\t%f\tK[%4d,%4d]\t%f" % (iteration, step, C, i, j, abs(dJ[i,j])))
        if spq==1: print("%d\t%f\t%.3f" % (iteration, C, s/STEP))
        #
        G = K2GAM(K)
        c_traj[iteration+1] = C, GAM2pearson(G, fdx), time.time()
    return [K, GAM2P(G,fdx), c_traj, paras_fit]

# Stochastic single variable update
def rs_svu(KZ, ITERATIONS=100, STEP1=2000, STEP2=5000, ALPHA1=0.002, ALPHA2=0.0001):
    paras_fit = "%4d\t%4d\t%4d\t%f\t%f" % (ITERATIONS, STEP1, STEP2, ALPHA1, ALPHA2)
    # Preparation
    K = KZ
    G = K2GAM(K)
    C = GAM2cost(G, fdx, ldx)
    c_traj = zeros((ITERATIONS+1, 3)); c_traj[0] = C, GAM2pearson(G, fdx), time.time(); # {cost, 1-pearsonr, time}
    # Minimization
    for iteration in range(ITERATIONS):
        s = 0
        for step1 in range(STEP1):
            i, j = Random_Int_Pair_Adjacent(N)
            tmp_kij = K[i,j]
            dkij = random.uniform(0, ALPHA1)
            if G[i,j] > GAM_obs[i,j]:
                K[i,j] += dkij
                K[j,i] += dkij
            else:
                K[i,j] -= dkij
                K[j,i] -= dkij
            tmp_G = K2GAM(K)
            tmp_C = GAM2cost(tmp_G, fdx, ldx)
            if tmp_C < C:
                G = tmp_G.copy()
                C = tmp_C
                s+= 1
                if spq==2: print("1\t%d\tA\t%04d\t%f\tK[%4d,%4d]\t%f" % (iteration, step1, C, i, j, abs(dkij)))
            else:
                K[i,j] = tmp_kij
                K[j,i] = tmp_kij
                if spq==2: print("0\t%d\tA\t%04d\t%f\tK[%4d,%4d]\t%f" % (iteration, step1, C, i, j, abs(dkij)))
        if spq==1: print("%d\t%f\t%.3f" % (iteration, C, s/STEP1))
        #
        s = 0
        for step2 in range(STEP2):
            i, j = Random_Int_Pair_Nonadj(N, nps)
            tmp_kij = K[i,j]
            dkij = random.uniform(0, ALPHA2)
            if G[i,j] > GAM_obs[i,j]:
                K[i,j] += dkij
                K[j,i] += dkij
            else:
                K[i,j] -= dkij
                K[j,i] -= dkij
            tmp_G = K2GAM(K)
            tmp_C = GAM2cost(tmp_G, fdx, ldx)
            if tmp_C < C:
                G = tmp_G.copy()
                C = tmp_C
                s+= 1
                if spq==2: print("1\t%d\tB\t%04d\t%f\tK[%4d,%4d]\t%f" % (iteration, step2, C, i, j, abs(dkij)))
            else:
                K[i,j] = tmp_kij
                K[j,i] = tmp_kij
                if spq==2: print("0\t%d\tB\t%04d\t%f\tK[%4d,%4d]\t%f" % (iteration, step2, C, i, j, abs(dkij)))
        if spq==1: print("%d\t%f\t%.3f" % (iteration, C, s/STEP2))
        #
        c_traj[iteration+1] = C, GAM2pearson(G, fdx), time.time()
    return [K, GAM2P(G,fdx), c_traj, paras_fit]

# Gradient descent
def gd(KZ, ITERATIONS=100, STEP1=25, STEP2=25, NK1=80, NK2=200, EPS1=1.0e-2, EPS2=1.0e-7):
    paras_fit = "%4d\t%4d\t%4d\t%4d\t%4d\t%f\t%f" % (ITERATIONS, STEP1, STEP2, NK1, NK2, EPS1, EPS2)
    # Preparation
    K = KZ
    G = K2GAM(K)
    C = GAM2cost(G, fdx, ldx)
    c_traj = zeros((ITERATIONS+1, 3)); c_traj[0] = C, GAM2pearson(G, fdx), time.time(); # {cost, 1-pearsonr, time}
    # Minimization
    for iteration in range(ITERATIONS):
        for step1 in range(STEP1):
            k_mask = Random_Int_Pairs_Adjacent(N, NK1)
            v = K2dJ(K, fdx, ldx)*k_mask*EPS1
            K = K - v
            G = K2GAM(K)
            C = GAM2cost(G, fdx, ldx)
            if spq==2: print("1\t%d\tA\t%04d\t%f\t%f" % (iteration, step1, C, abs(v).max()))
        if isnan(C): sys.exit()
        if spq==1: print("%d\t%f\t%.3f" % (iteration, C, 1.))
        #
        for step2 in range(STEP2):
            k_mask = Random_Int_Pairs_Nonadj(N, nps, NK2)
            v = K2dJ(K, fdx, ldx)*k_mask*EPS2
            K = K - v
            G = K2GAM(K)
            C = GAM2cost(G, fdx, ldx)
            if spq==2: print("1\t%d\tB\t%04d\t%f\t%f" % (iteration, step2, C, abs(v).max()))
        if isnan(C): sys.exit()
        if spq==1: print("%d\t%f\t%.3f" % (iteration, C, 1.))
        #
        c_traj[iteration+1] = C, GAM2pearson(G, fdx), time.time()
    return [K, GAM2P(G,fdx), c_traj, paras_fit]

def rms_prop(KZ, ITERATIONS=100, STEP1=5, STEP2=5, NK1=80, NK2=200, EPS1=0.05, EPS2=0.00005):
    beta = 0.9; epsilon = 1.0e-12;
    paras_fit = "%4d\t%4d\t%4d\t%4d\t%4d\t%f\t%f\t%f\t%e" % (ITERATIONS, STEP1, STEP2, NK1, NK2, EPS1, EPS2, beta, epsilon)
    # Preparation
    K = KZ
    G = K2GAM(K)
    C = GAM2cost(G, fdx, ldx)
    c_traj = zeros((ITERATIONS+1, 3)); c_traj[0] = C, GAM2pearson(G, fdx), time.time(); # {cost, 1-pearsonr, time}
    # Minimization
    grad_sq1 = zeros((N, N)); grad_sq2 = zeros((N, N))
    for iteration in range(ITERATIONS):
        C_bk = C
        for step1 in range(STEP1):
            k_mask = Random_Int_Pairs_Adjacent(N, NK1)
            #
            g1 = K2dJ(K, fdx, ldx)*k_mask
            grad_sq1 = beta*grad_sq1 + (1-beta)*g1*g1
            v = g1/sqrt(grad_sq1 + epsilon)
            #
            K = K - v*EPS1
            C = GAM2cost(K2GAM(K), fdx, ldx)
            if spq==2: print("1\t%d\tA\t%04d\t%f\t%f\t%e" % (iteration, step1, C, abs(v).max(), EPS1))
        if isnan(C): sys.exit()
        if spq==1: print("%d\t%f\t%.3f" % (iteration, C, 1.))
        if C >= C_bk: EPS1 *= 0.2
        #
        C_bk = C
        for step2 in range(STEP2):
            k_mask = Random_Int_Pairs_Nonadj(N, nps, NK2)
            #
            g2 = K2dJ(K, fdx, ldx)*k_mask
            grad_sq2 = beta*grad_sq2 + (1-beta)*g2*g2
            v = g2/sqrt(grad_sq2 + epsilon)
            #
            K = K - v*EPS2
            C = GAM2cost(K2GAM(K), fdx, ldx)
            if spq==2: print("1\t%d\tB\t%04d\t%f\t%f\t%e" % (iteration, step2, C, abs(v).max(), EPS2))
        if isnan(C): sys.exit()
        if spq==1: print("%d\t%f\t%.3f" % (iteration, C, 1.))
        if C >= C_bk: EPS2 *= 0.2
        #
        G = K2GAM(K)
        c_traj[iteration+1] = C, GAM2pearson(G, fdx), time.time()
    return [K, GAM2P(G,fdx), c_traj, paras_fit]

def adam(KZ, ITERATIONS=100, STEP1=5, STEP2=5, NK1=80, NK2=200, EPS1=0.05, EPS2=0.00005):
    gamma = 0.9; beta = 0.99; epsilon = 1.0e-12;
    paras_fit = "%4d\t%4d\t%4d\t%4d\t%4d\t%f\t%f\t%f\t%f\t%e" % (ITERATIONS, STEP1, STEP2, NK1, NK2, EPS1, EPS2, gamma, beta, epsilon)
    # Preparation
    K = KZ
    G = K2GAM(K)
    C = GAM2cost(G, fdx, ldx)
    c_traj = zeros((ITERATIONS+1, 3)); c_traj[0] = C, GAM2pearson(G, fdx), time.time(); # {cost, 1-pearsonr, time}
    # Minimization
    v1 = zeros((N, N)); grad_sq1 = zeros((N, N))
    v2 = zeros((N, N)); grad_sq2 = zeros((N, N))
    EPS1Z = EPS1; EPS2Z = EPS2;
    DEPS1 = EPS1*0.01; DEPS2 = EPS2*0.01
    for iteration in range(ITERATIONS):
        C_bk = C
        for step1 in range(STEP1):
            k_mask = Random_Int_Pairs_Adjacent(N, NK1)
            #
            g1 = K2dJ(K, fdx, ldx)*k_mask
            v1 = gamma*v1 + (1-gamma)*g1
            grad_sq1 = beta*grad_sq1 + (1-beta)*g1*g1
            v_hat1 = v1/(1-gamma)
            grad_sq_hat1 = grad_sq1/(1-beta)
            v = v_hat1/sqrt(grad_sq_hat1 + epsilon)
            #
            K = K - v*EPS1
            C = GAM2cost(K2GAM(K), fdx, ldx)
            if spq==2: print("1\t%d\tA\t%04d\t%f\t%f\t%e" % (iteration, step1, C, abs(v).max(), EPS1))
        if isnan(C): sys.exit()
        if spq==1: print("%d\t%f\t%.3f" % (iteration, C, 1.))
        if C >= C_bk:
            EPS1 -= DEPS1
            EPS1 = max(EPS1, DEPS1)
        else:
            EPS1 += DEPS1
            EPS1 = min(EPS1, EPS1Z)
        #
        C_bk = C
        for step2 in range(STEP2):
            k_mask = Random_Int_Pairs_Nonadj(N, nps, NK2)
            #
            g2 = K2dJ(K, fdx, ldx)*k_mask
            v2 = gamma*v2 + (1-gamma)*g2
            grad_sq2 = beta*grad_sq2 + (1-beta)*g2*g2
            v_hat2 = v2/(1-gamma)
            grad_sq_hat2 = grad_sq2/(1-beta)
            v = v_hat2/sqrt(grad_sq_hat2 + epsilon)
            #
            K = K - v*EPS2
            C = GAM2cost(K2GAM(K), fdx, ldx)
            if spq==2: print("1\t%d\tB\t%04d\t%f\t%f\t%e" % (iteration, step2, C, abs(v).max(), EPS2))
        if isnan(C): sys.exit()
        if spq==1: print("%d\t%f\t%.3f" % (iteration, C, 1.))
        if C >= C_bk:
            EPS2 -= DEPS2
            EPS2 = max(EPS2, DEPS2)
        else:
            EPS2 += DEPS2
            EPS2 = min(EPS2, EPS2Z)
        #
        G = K2GAM(K)
        c_traj[iteration+1] = C, GAM2pearson(G, fdx), time.time()
    return [K, GAM2P(G,fdx), c_traj, paras_fit]

if True:
    xi2 = 1.0**2 # r^2_c
    # Read random seed
    fx = 'mdRandomSeeds'
    if not os.path.isfile(fx):
        print('Cannot find '+fx)
        sys.exit()
    else:
        lineno = 0
        with open(fx) as fr:
            for line in fr:
                if lineno == rdx:
                    rds = int(line.strip())
                    break
                lineno += 1
        random.seed(rds)

    # Read Hi-C assuring that p_ii = 1, and there is no ZERO
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
        # Check the validity of p_{ij}
        if sum(P_obs==0):
            print("%s contains Zero" %(fhic))
            sys.exit()
        #
        N = len(P_obs)
        # P_obs[20,:] = nan; P_obs[:,20] = nan # Set NaN iterms manually !
        nanMask = (~isnan(P_obs)).astype(float)
        nanMaskN= isnan(P_obs)
        # This should not happen, but...
        if fdx == 1:
            for i in range(0, N):
                for j in range(0, N):
                    if P_obs[i,j] > 1:
                        P_obs[i,j] = 1.0

    # Observation relevant quantities
    npc = N*(N-1)/2.- isnan(triu(P_obs,1)).sum()
    nps = int((N-1)*(N-2)/2)
    lnP_obs = log10(P_obs)
    GAM_obs = P2GAM(P_obs, fdx)
    if ldx==2:
        sO1 = nansum(lnP_obs)/2.
        sO2 = nansum(lnP_obs**2)/2.
        vOn = sqrt(npc*sO2 - sO1**2)
    elif ldx==3:
        WS = (1./(GAM_obs**2 + eye(N)) - eye(N)) / npc
    trilMask = where( (tril(ones((N,N)))+nanMaskN)==0 ) # Non-NaN matrix indices of j>i

    # Minimization
    K_fit = Init_K(N, INIT_K0=0.3)
    spq = 0 # 0/1/2: noScreenPrint/ screenPrint/ screenPrintAtEveryUpdate
    if mdx == 0:
        K_fit, P_fit, c_traj, paras_fit = sd_svu(K_fit, ITERATIONS=10, STEP=700, EPS=1.0e-3) # STEP=npc
    elif mdx==1:
        K_fit, P_fit, c_traj, paras_fit = rs_svu(K_fit, ITERATIONS=10, STEP1=2000, STEP2=5000, ALPHA1=0.002, ALPHA2=0.0001)
    elif mdx==2:
        K_fit, P_fit, c_traj, paras_fit = gd(K_fit, ITERATIONS=1000, STEP1=5, STEP2=5, NK1=80, NK2=200, EPS1=1.0e-2, EPS2=1.0e-7)
    elif mdx==3:
        K_fit, P_fit, c_traj, paras_fit = rms_prop(K_fit, ITERATIONS=1000, STEP1=5, STEP2=5, NK1=80, NK2=200, EPS1=1.0e-2 if ldx==2 else 5.0e-2, EPS2=1.0e-5 if ldx==2 else 5.0e-5)
    elif mdx==4:
        adam_its = max(1000, round(N/10.)*100)
        adam_sp1 = round(N*0.8)
        adam_sp2 = round(N*2.0)
        K_fit, P_fit, c_traj, paras_fit = adam(K_fit, ITERATIONS=adam_its, STEP1=5, STEP2=5, NK1=adam_sp1, NK2=adam_sp2, EPS1=1.0e-2 if ldx==2 else 5.0e-2, EPS2=1.0e-5 if ldx==2 else 5.0e-5)

    # Save results
    dataDir = fhic[:fhic.rfind('.')]
    os.makedirs(dataDir, exist_ok=True)
    fo = "%s/N%d-F%d-L%d-M%d-r%d"%(dataDir, N, fdx, ldx, mdx, rdx)

    c_traj[:,2] = c_traj[:,2] - c_traj[0,2]
    saveLg(fo+'.log', c_traj, "#%s\n#cost pearsonr(log10(p)) systemTime\n"%(paras_fit))

    saveMx(fo+'.K_fit', K_fit, "#K_fit N %d min: %11.5e max: %11.5e\n"%(N, K_fit.min(), K_fit.max()))

    P_fit[nanMaskN] = nan; fmx = nan_to_num(P_fit, copy=True)
    saveMx(fo+'.P_fit', P_fit, "#P_fit N %d min: %11.5e max: %11.5e\n"%(N, min(fmx[fmx>0]), max(fmx[fmx<1])))

    # fmx = nan_to_num(P_obs, copy=True) # NaN to 0
    # saveMx(fo+'.P_obs', P_obs, "#P_obs N %d min: %11.5e max: %11.5e\n"%(N, min(fmx[fmx>0]), max(fmx[fmx<1])))
