load '../tino/piyg.pal'
gs(vp) = vp==75 ? 32151060 : (vp==68 ? 32137176 : NaN)
ge(vp) = vp==75 ? 32151883 : (vp==68 ? 32137426 : NaN)
f(c,vp) = sprintf("%s_alpha/N150-F0-L1-M4-r0.K_fit.q3.vp%d-%d.vm", c==0 ? 'ESC' : 'ERY', gs(vp), ge(vp))

N = 150
#set xrange [-0.5:N-0.5]
#set yrange [-0.5:N-0.5]
set xrange [40:120]
set yrange [40:120]
set size square
set tics front out scale 0.3
set cbtics in
unset key
set xtics 0,30,150 nomirror; set mxtics 3
set ytics 0,30,150 nomirror; set mytics 3

set cbrange [-5.5:-1.5]; set cbtics -6,1,0; set cbtics in
set multiplot layout 2,2
set object 1 rectangle from graph 0,0 to graph 1,1 fc rgb 'black' lw 0

set title 'ESC R2'
plot f(0,75) matrix u 1:2:($3>0 ? log10($3) : NaN) w image
set title 'ERY R2'
plot f(1,75) matrix u 1:2:($3>0 ? log10($3) : NaN) w image
set title 'ESC HS-39'
plot f(0,68) matrix u 1:2:($3>0 ? log10($3) : NaN) w image
set title 'ERY HS-39'
plot f(1,68) matrix u 1:2:($3>0 ? log10($3) : NaN) w image

unset multiplot
