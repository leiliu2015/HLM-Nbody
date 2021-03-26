#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import h5py
from numpy import *
from scipy.linalg import toeplitz
import warnings
warnings.filterwarnings('ignore')

# Produce 3D chromatin conformations from K-matrix, based on Onami's work.
if not len(sys.argv) == 4:
    print('usage:: python gc_K2cfg.py N NSamples k0')
    sys.exit()
N  = int(sys.argv[1])
nc = int(sys.argv[2])
k0=float(sys.argv[3])
rds= 1274; random.seed(rds)

if True:
    # Define K-matrix
    ks = zeros(N)
    ks[1] = k0
    K_fit = toeplitz(ks)
    # K to Laplacian matrix
    d = sum(K_fit, axis=0)
    D = diag(d)
    L = D - K_fit
    # Eigenvalues and eigenvectors
    lam, Q = linalg.eigh(L)

    fo = h5py.File("g%d"%(N) + '.xyz.h5', 'w')
    xyz= fo.create_dataset('xyz', (nc, N, 3), dtype='f')
    for c in range(0, nc):
        # Collective coordinates X
        X = zeros((N, 3))
        for k in range(3):
            X[1:,k] = sqrt(1.0/lam[1:])*random.randn(N-1) # PHi-C used sqrt(1.0/3/lam[1:])
        # X -> 3D coordinates R
        R = zeros((N, 3))
        for k in range(3):
            R[:,k] = dot(Q, X[:,k])
        xyz[c,:,:] = R
    fo.close()
