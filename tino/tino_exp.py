#
# Copyright (C) 2020 Lei Liu & Changbong Hyeon
# This file is part of HLM-tino.
#
import os
import sys
import h5py
from numpy import *

# Save a list of matrices, with a shape of N*N*M, in a TXT file
def saveLmTx(fn, mxs, cts):
    N, M = shape(mxs)[1:3]
    fw = open(fn, 'w')
    fw.write(cts[0])
    for m in range(M):
        fw.write(cts[m+1])
        fw.write("#min: %+12.5e max: %+12.5e\n"%(nanmin(mxs[:,:,m]), nanmax(mxs[:,:,m])))
        for i in range(0, N):
            lt = ''
            for j in range(0, N):
                lt += "%12s "%('NaN') if isnan(mxs[i,j,m]) else "%+11.5e "%(mxs[i,j,m])
            fw.write(lt+'\n')
        fw.write('\n\n')
    fw.close()

# Compute p_{3}(m,t) where m=min(j-i,k-j) and t=k-i for i<j<k
if not len(sys.argv) == 3:
    print('usage:: python tino_exp.py xxx.[p/q]3.h5 saveTXT/H5[0/1]')
    sys.exit()
fx = str(sys.argv[1])
h5q= int(sys.argv[2])

if True:
    if not os.path.isfile(fx):
        print('Cannot find '+fx)
        sys.exit()
    else:
        fr = h5py.File(fx, 'r')
        p3 = fr.get('cm')[()]
        fr.close()
        N = len(p3)
        # Expected 3-body contact probability as a function of major and minor spans, {mean, std, #ofTriplets}
        p3_exp = zeros((N, N, 3))*nan
        for mas in range(2, N):
            for mis in range(1, int(mas/2.)+1):
                ps = []
                for i in range(0, N-mas):
                    k = i+mas
                    ji= i+mis
                    jk= k-mis
                    if ji == jk:
                        ps.append(p3[i,ji,k])
                    else:
                        ps.append(p3[i,ji,k])
                        ps.append(p3[i,jk,k])
                ps = array(ps)
                p3_exp[mas, mis, 0] = nanmean(ps)
                p3_exp[mas, mis, 1] = nanstd(ps)
                p3_exp[mas, mis, 2] = sum(isfinite(ps))
        #
        if h5q:
            fw = h5py.File(fx[:-3]+'.exp.h5', 'w')
            fw.create_dataset('p3_exp', data=p3_exp, dtype='f')
            fw.close()
        else:
            saveLmTx(fx[:-3]+'.exp.txt', p3_exp, ["#N: %d\n"%(N), '#mean\n', '#std\n', '#ofTriplets\n'])
