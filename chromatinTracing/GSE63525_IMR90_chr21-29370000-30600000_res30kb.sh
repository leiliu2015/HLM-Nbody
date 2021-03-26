#!/bin/bash

fhic=GSE63525_IMR90_chr21-29370000-30600000_res30kb.txt
ddir=GSE63525_IMR90_chr21-29370000-30600000_res30kb
ctcf=wgEncodeAwgTfbsSydhImr90CtcfbIggrabUniPk.narrowPeak.chr21-29370000-30600000-res30000
numberOfReplicas=1 # 5
cp ../tino/mdRandomSeeds ./

mdx=4
for ((fdx=0; "$fdx" <= "1"; fdx++))
do
    if (("$fdx" == "0"))
    then
        ldx=1
    else
        ldx=3
    fi
    #
    for ((rdx=0; "$rdx" < "$numberOfReplicas"; rdx++))
    do
        lb=F${fdx}-L${ldx}-M${mdx}-r${rdx}
        echo ${lb}
        python ../tino/tino_nan.py ${fdx} ${ldx} ${mdx} ${rdx} ${fhic}
    done
done

fits=(0 0)
mdx=4
for ((fdx=0; "$fdx" <= "1"; fdx++))
do
    if (("$fdx" == "0"))
    then
        ldx=1
    else
        ldx=3
    fi
    #
    lb=F${fdx}-L${ldx}-M${mdx}-r${fits[$fdx]}
    echo $lb
    python ../tino/tino_ps.py ${fhic} ${ddir} N41-${lb}.P_fit
    python zhuang2018_K2CP.py ${ddir}/N41-${lb}.K_fit 0 ${ctcf}
done
