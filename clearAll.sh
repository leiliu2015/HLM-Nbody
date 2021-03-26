#!/bin/bash

fds=(numeric_mESC numeric_GM12878 chromatinTracing tric mc4c sprite)
sfs=(Bonev_ES_chr8-42100-44500kb_res25kb GSE63525_GM12878_chr5-90000-100000kb_res50kb GSE63525_IMR90_chr21-29370000-30600000_res30kb E*_alpha *_CTCF TAD1[09])
for ((fdx=0; "$fdx" < 6; fdx++))
do
    fd=${fds[$fdx]}
    cd ${fd}
    if [ -f mdRandomSeeds ]
    then
        echo clean ./${fd}
        rm -r mdRandomSeeds ${sfs[$fdx]}
    fi
    cd ../
done

cd ./gaussianChain
if [ -f g20.p2.h5 ]
then
   echo clean ./gaussianChain
   rm *h5 *zm *txt *Profile *readSize
fi
cd ../
