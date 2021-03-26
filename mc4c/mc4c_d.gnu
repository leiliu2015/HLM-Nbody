N  = 13
set xrange [-0.5:N-0.5]
set yrange [-0.5:N-0.5]
set size square
set tics front out scale 0.75 nomirror
set palette defined(0'blue', 0.5'white', 1'red')
unset key

set object 1 rect from graph 0,0 to graph 1,1 fc rgb 'black'
set multiplot layout 1,2 spacing 3,0

fa = 'WaplKO_CTCF/N128-F0-L1-M4-r0.K_fit.q4.vp33-115.txt-soi2-diff'
set cbrange [-5e-6:+5e-6]; set cbtics in
set title 'P_{dWAPL} - P_{WT} (VP: E,K)'
plot fa matrix rowheader columnheader u 1:2:($1==$2 ? NaN : $3) w image pixels

fb = 'WaplKO_CTCF/N128-F0-L1-M4-r0.K_fit.q5.vp33-66-115.txt-soi2-diff'
set cbrange [-2e-8:+2e-8]; set cbtics in
set title 'P_{dWAPL} - P_{WT} (VP: E,H,K)'
plot fb matrix rowheader columnheader u 1:2:($1==$2 ? NaN : $3) w image pixels
unset multiplot


