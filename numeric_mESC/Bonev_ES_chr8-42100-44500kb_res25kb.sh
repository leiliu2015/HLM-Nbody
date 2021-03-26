#!/bin/bash

fhic=Bonev_ES_chr8-42100-44500kb_res25kb.txt
ddir=Bonev_ES_chr8-42100-44500kb_res25kb
numberOfReplicas=1 # 5
cp ../tino/mdRandomSeeds ./

for ((fdx=0; "$fdx" <= "1"; fdx++))
do
    for ((ldx=0; "$ldx" <= "3"; ldx++))
    do
        for ((mdx=1; "$mdx" <= "4"; mdx++))
        do
            for ((rdx=0; "$rdx" < "$numberOfReplicas"; rdx++))
            do
                lb=F${fdx}-L${ldx}-M${mdx}-r${rdx}
                echo ${lb}
                python ../tino/tino_nan.py ${fdx} ${ldx} ${mdx} ${rdx} ${fhic}
            done
            python ../tino/averageLogs.py ${ddir} N97-F${fdx}-L${ldx}-M${mdx} ${numberOfReplicas}
        done
    done
done

fits=(F1-L1-M4-r0 F1-L2-M4-r0 F1-L3-M4-r0)
for ((i=0; "$i" < "3"; i++))
do
    fit=${fits[$i]}
    python ../tino/tino_ps.py ${fhic} ${ddir} N97-${fit}.P_fit
done
