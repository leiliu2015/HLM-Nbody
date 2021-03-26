#!/bin/bash

N=20
k=3.0
vp=7

# 3D coordinates will be used by gc_cfgs2P4.py & gc_deLaat.py
if [ ! -f ./g${N}.xyz.h5 ]
then
    python ./gc_K2cfgs.py ${N} 100000 ${k}
fi

for ((fdx=0; "$fdx" < "2"; fdx++))
do
    for ((h5q=0; "$h5q" < "2"; h5q++))
    do
        python ./gc_K2P.py ${N} ${fdx} 2 ${h5q} ${k}
        python ./gc_K2P.py ${N} ${fdx} 3 ${h5q} ${k}
        if (("$h5q" == "0"))
        then
            if (("$fdx" == "0"))
            then
                python ./gc_K2P.py ${N} ${fdx} 4 ${h5q} ${k}
            else
                python ./gc_cfgs2P4.py ${N} 1 ${h5q}
            fi
        fi
    done
    #
    python ./gc_tamm.py ${N} ${fdx} ${k}
    python ./gc_tanay.py ${N} ${fdx} ${vp} 2
    python ./gc_deLaat.py g${N}.xyz.h5 ${fdx} ${vp} 10
done

