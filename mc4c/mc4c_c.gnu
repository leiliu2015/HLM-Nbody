N=13
set xrange [-0.5:N-0.5]
set yrange [-0.5:N-0.5]
set size square
set tics front out scale 0
unset key

# p_{ijk}(WalpKO) - p_{ijk}(WT)
set cbrange [-0.0007:0.0007]
set palette defined(0'blue',0.5'white',1'red')
unset colorbox

vp(i) = i==0 ? 7 : (i==1 ? 19 : (i==2 ? 29 : (i==3 ? 30 : (i==4 ? 33 : (i==5 ? 55 : (i==6 ? 63 : (i==7 ? 66 : (i==8 ? 81 : (i==9 ? 102 : (i==10 ? 115 : (i==11 ? 116 : 119)))))))))))
gs(i) = vp(i)*10000 + 120800000
ge(i) = gs(i) + 10000
fx(i) = sprintf("WaplKO_CTCF/N128-F0-L1-M4-r0.K_fit.q3.vp%d-%d.txt-soi2-diff",gs(i),ge(i))

set object 1 rect from graph 0,0 to graph 1,1 fc rgb 'black'
set term wxt 0 size 800,600
set multiplot layout 3,4 spacing 0,0
do for [i=0:11] {
  set title 'P_{dWAPL} - P_{WT}' offset 0,-0.8
  plot fx(i) matrix rowheader columnheader u 1:2:3 w image
}
unset multiplot

