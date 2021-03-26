#!/bin/bash

codesDir=../tino
experDir=/home/kapok/Downloads/triplets/CHROMATIX
rois=(19 10)
Ns=(89 104)
cp ../tino/mdRandomSeeds ./

mdx=4
rmx=1
for ((idx=0; "$idx" < "2"; idx++))
do
    tdx=${rois[$idx]}
    N=${Ns[$idx]}
    printf -v ddir "TAD%02d" ${tdx}
    fhic=${ddir}.txt
    for ((fdx=0; "$fdx" < "1"; fdx++))
    do
        if (("$fdx" == "0"))
        then
            ldx=1
            pq=q
        else
            ldx=3
            pq=p
        fi
        #
        for ((rdx=0; "$rdx" < "$rmx"; rdx++))
        do
            lb=N${N}-F${fdx}-L${ldx}-M${mdx}-r${rdx}
            echo ${ddir} ${lb}
            python ${codesDir}/tino_nan.py ${fdx} ${ldx} ${mdx} ${rdx} ${fhic}
        done
        #
        lb=N${N}-F${fdx}-L${ldx}-M${mdx}-r0
        python ${codesDir}/tino_ps.py ${fhic} ${ddir} ${lb}.P_fit
        python ${codesDir}/tino_K2P.py ${ddir}/${lb}.K_fit ${fdx} 3
        python ${codesDir}/tino_exp.py ${ddir}/${lb}.K_fit.${pq}3.h5 1
        python sprite_3bodyPromoter.py ${tdx} ${ddir}/N${N}-F${fdx}-L${ldx}-M${mdx}-r0.K_fit.${pq}3.h5
        python sprite_3bodyPromoter_summary.py ${tdx} ${ddir}/N${N}-F${fdx}-L${ldx}-M${mdx}-r0.K_fit.${pq}3.h5
        python sprite_3body_ENPS.py ${tdx} ${ddir}/N${N}-F${fdx}-L${ldx}-M${mdx}-r0.K_fit.${pq}3.h5
    done
done
