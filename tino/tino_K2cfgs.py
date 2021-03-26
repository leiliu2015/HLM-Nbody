#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import h5py
from numpy import *

# Produce 3D chromatin conformations from K-matrix, based on Onami's work.
if not len(sys.argv) == 3:
    print('usage:: python tino_K2cfg.py xxx.K_fit NSamples')
    sys.exit()
fk = str(sys.argv[1])
nc = int(sys.argv[2])
rds= 1274; random.seed(rds)

# Read K-matrix
if True:
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
    # K to Laplacian matrix
    d = sum(K_fit, axis=0)
    D = diag(d)
    L = D - K_fit
    # Eigenvalues and eigenvectors
    lam, Q = linalg.eigh(L)

    fo = h5py.File(fk[:-6] + '.xyz.h5', 'w')
    xyz= fo.create_dataset('xyz', (nc, N, 3), dtype='f')
    for c in range(0, nc):
        # Collective coordinates X
        X = zeros((N, 3))
        for k in range(3):
            # X[1:,k] = sqrt(1.0/3/lam[1:])*random.randn(N-1)
            X[1:,k] = sqrt(1.0/lam[1:])*random.randn(N-1)
        # X -> 3D coordinates R
        R = zeros((N, 3))
        for k in range(3):
            R[:,k] = dot(Q, X[:,k])
        xyz[c,:,:] = R
    fo.close()
