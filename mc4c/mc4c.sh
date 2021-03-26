#!/bin/bash

codesDir=../tino
experDir=./mc4c-data
rois=(Hap1_CTCF WaplKO_CTCF)
Ns=(128 128)
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

for ((fdx=0; "$fdx" < "1"; fdx++))
do
    #WPL_WTD
    python ./mc4c_3bodyVP.py Hap1_CTCF/N128-F0-L1-M4-r0.K_fit.q3.h5 120800000 10000 121126666 121145791
    python ./mc4c_3bodyVP_summary.py ${experDir}/MC4C_WPL-WTD.res10000.cm3 Hap1_CTCF/N128-F0-L1-M4-r0.K_fit.q3.vp121126666-121145791.txt
    #WPL_WTC
    python ./mc4c_3bodyVP.py Hap1_CTCF/N128-F0-L1-M4-r0.K_fit.q3.h5 120800000 10000 121941808 121960933
    python ./mc4c_3bodyVP_summary.py ${experDir}/MC4C_WPL-WTC.res10000.cm3 Hap1_CTCF/N128-F0-L1-M4-r0.K_fit.q3.vp121941808-121960933.txt
    #WPL_KOD
    python ./mc4c_3bodyVP.py WaplKO_CTCF/N128-F0-L1-M4-r0.K_fit.q3.h5 120800000 10000 121126666 121145791
    python ./mc4c_3bodyVP_summary.py ${experDir}/MC4C_WPL-KOD.res10000.cm3 WaplKO_CTCF/N128-F0-L1-M4-r0.K_fit.q3.vp121126666-121145791.txt
    #WPL_KOC
    python ./mc4c_3bodyVP.py WaplKO_CTCF/N128-F0-L1-M4-r0.K_fit.q3.h5 120800000 10000 121941808 121960933
    python ./mc4c_3bodyVP_summary.py ${experDir}/MC4C_WPL-KOC.res10000.cm3 WaplKO_CTCF/N128-F0-L1-M4-r0.K_fit.q3.vp121941808-121960933.txt
done

vps=(7 19 29 30 33 55 63 66 81 102 115 116 119)
gz=120800000
N=128
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
    for ((vdx=0; "$vdx" < "13"; vdx++))
    do
        vp=${vps[$vdx]}
        gs=$(($vp*10000))
        gs=$(($gs+$gz))
        ge=$(($gs+10000))
        printf -v vs "%03d" ${vp}
        for ((idx=0; "$idx" < "2"; idx++))
        do
            ddir=${rois[$idx]}
            echo ${ddir} ${lb} ${vdx}
            python ./mc4c_3bodyVP.py ${ddir}/${lb}.K_fit.${pq}3.h5 ${gz} 10000 ${gs} ${ge}
            python ./mc4c_h5association_soi2zm.py ${experDir}/GSM249387x_CTCF.chr8-120800000-122080000-res10000 ${ddir}/${lb}.K_fit.${pq}3.vp${gs}-${ge}.txt
        done
        python ./mc4c_h5association_soi2zm_diff.py WaplKO_CTCF/${lb}.K_fit.${pq}3.vp${gs}-${ge}.txt-soi2 Hap1_CTCF/${lb}.K_fit.${pq}3.vp${gs}-${ge}.txt-soi2
    done
done

# q_{4} at VPs {E, K}
for ((idx=0; "$idx" < "2"; idx++))
do
    ddir=${rois[$idx]}
    lb=N128-F0-L1-M4-r0
    echo $ddir
    python ./mc4c_nbodyVP.py ${ddir}/${lb}.K_fit 33 115
    python ./mc4c_h5association_soi2zm.py ${experDir}/GSM249387x_CTCF.chr8-120800000-122080000-res10000 ${ddir}/${lb}.K_fit.q4.vp33-115.txt
done
python ./mc4c_nbodyVP_diff.py WaplKO_CTCF/${lb}.K_fit.q4.vp33-115.txt Hap1_CTCF/${lb}.K_fit.q4.vp33-115.txt
python ./mc4c_h5association_soi2zm_diff.py WaplKO_CTCF/${lb}.K_fit.q4.vp33-115.txt-soi2 Hap1_CTCF/${lb}.K_fit.q4.vp33-115.txt-soi2 1

# q_{5} at VPs {E, H, K}
for ((idx=0; "$idx" < "2"; idx++))
do
    ddir=${rois[$idx]}
    lb=N128-F0-L1-M4-r0
    echo $ddir
    python ./mc4c_nbodyVP.py ${ddir}/${lb}.K_fit 33 66 115
    python ./mc4c_h5association_soi2zm.py ${experDir}/GSM249387x_CTCF.chr8-120800000-122080000-res10000 ${ddir}/${lb}.K_fit.q5.vp33-66-115.txt
done
python ./mc4c_nbodyVP_diff.py WaplKO_CTCF/${lb}.K_fit.q5.vp33-66-115.txt Hap1_CTCF/${lb}.K_fit.q5.vp33-66-115.txt
python ./mc4c_h5association_soi2zm_diff.py WaplKO_CTCF/${lb}.K_fit.q5.vp33-66-115.txt-soi2 Hap1_CTCF/${lb}.K_fit.q5.vp33-66-115.txt-soi2 1

