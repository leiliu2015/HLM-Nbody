#!/bin/bash

codesDir=../tino
experDir=./tric-data
rois=(ESC_alpha ERY_alpha)
cells=(ES C57)
Ns=(150 150)
cp ../tino/mdRandomSeeds ./

mdx=4
rmx=1
for ((idx=0; "$idx" < "2"; idx++))
do
    fhic=${rois[$idx]}.txt
    ddir=${rois[$idx]}
    N=${Ns[$idx]}
    for ((fdx=0; "$fdx" < "1"; fdx++))
    do
        if (("$fdx" == "0"))
        then
            ldx=1
        else
            ldx=3
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
    done
done

vps=(75 68 91 99 110 112)
vpNames=(R2 HS-39 Hba-a1 Hbq1a HS+44 HS+48)
gss=(32151060 32137176 32182969 32199491 32220918 32224323)
ges=(32151883 32137426 32183821 32200139 32221720 32227298)
gz=32000000
N=150
for ((fdx=0; "$fdx" < "1"; fdx++))
do
    if (("$fdx" == "0"))
    then
        lb=N${N}-F0-L1-M4-r0
        pq=q
    else
        lb=N${N}-F1-L3-M4-r0
        pq=p
    fi
    for ((vdx=0; "$vdx" < "6"; vdx++))
    do
        vp=${vps[$vdx]}
        gs=${gss[$vdx]}
        ge=${ges[$vdx]}
        printf -v vs "%03d" ${vp}
        for ((idx=0; "$idx" < "2"; idx++))
        do
            ddir=${rois[$idx]}
            #
            echo ${ddir} ${lb} ${vs} ${gs} ${ge}
            python tric_3bodyVP.py ${ddir}/${lb}.K_fit.${pq}3.h5 ${gz} 2000 ${gs} ${ge}
            if (("$vdx" < 2))
            then
                python tric_3bodyVP_summary.py ${experDir}/GSE107755_${cells[$idx]}_${vpNames[$vdx]}_chr11-32000000-32300000_res2000.txt ${ddir}/${lb}.K_fit.${pq}3.vp${gs}-${ge}.txt 0 1 
            fi
            python tric_h5association_soi2zm.py ${experDir}/GSE107755_soi.chr11-32000000-32300000_res2000 ${ddir}/${lb}.K_fit.${pq}3.vp${gs}-${ge}.txt
        done
        python tric_h5association_soi2zm_diff.py ERY_alpha/${lb}.K_fit.${pq}3.vp${gs}-${ge}.txt-soi2 ESC_alpha/${lb}.K_fit.${pq}3.vp${gs}-${ge}.txt-soi2
    done
done
