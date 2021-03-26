load '../tino/piyg.pal'
fa = 'Hap1_CTCF/N128-F0-L1-M4-r0.K_fit.q3.vp121126666-121145791.vm'
fb = 'Hap1_CTCF/N128-F0-L1-M4-r0.K_fit.q3.vp121941808-121960933.vm'
fc = 'WaplKO_CTCF/N128-F0-L1-M4-r0.K_fit.q3.vp121126666-121145791.vm'
fd = 'WaplKO_CTCF/N128-F0-L1-M4-r0.K_fit.q3.vp121941808-121960933.vm'

N = 128
set xrange [-0.5:N-0.5]
set yrange [-0.5:N-0.5]
set xtics 0,50,150; set mxtics 5
set ytics 0,50,150; set mytics 5
do for [i=0:2] {
  set xtics add(''i*50)
  set ytics add(''i*50)
} 
set size square
set tics front out scale 0.3
set cbtics in
unset key

set cbrange [-4.5:-2.5]; set cbtics -6,1,0; set cbtics in
set term wxt 0 size 600,600
set multiplot layout 2,2
set object 1 rectangle from graph 0,0 to graph 1,1 fc rgb 'black' lw 0

set title 'HAP1 WT (VP: E)'
plot fa matrix u 1:2:($3>0 ? log10($3) : NaN) w image
set title 'HAP1 dWAPL (VP: E)'
plot fc matrix u 1:2:($3>0 ? log10($3) : NaN) w image
set title 'HAP1 WT (VP: K)'
plot fb matrix u 1:2:($3>0 ? log10($3) : NaN) w image
set title 'HAP1 dWAPL (VP: K)'
plot fd matrix u 1:2:($3>0 ? log10($3) : NaN) w image

unset multiplot
