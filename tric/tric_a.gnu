load '../tino/piyg.pal'
fdx=0; ldx=1; mdx=4
fx(cdx) = sprintf("%s_alpha/N150-F%d-L%d-M%d-r0.P_fit.vm",cdx==0 ? 'ESC' : 'ERY', fdx,ldx,mdx)
N  = 150
rs = 0.002
gz = 32.0

set size square
set tics front out scale 0.75
unset key
set xrange [0:N]
set yrange [0:N]

i2g(i) = gz+i*rs
set xtics 1,50,N nomirror offset 0,0.4; set mxtics 5;
set ytics 1,50,N nomirror offset 0.8,0; set mytics 5;
do for [i=1:N:50] {
  set xtics add(sprintf("%.1f",i2g(i)) i)
  set ytics add(sprintf("%.1f",i2g(i)) i)
}

set cbrange [-3:-1]; set cbtics -5,1,0 in scale 0.5; set mcbtics 5
do for [i=-3:0:1] {
  set cbtics add(sprintf("10^{%d}",i) i)
}

set multiplot layout 1,2
set title 'ESC'
plot fx(0) matrix u ($1+0.5):($2+0.5):($3>0 ? log10($3) : NaN) w image notitle
set title 'ERY'
plot fx(1) matrix u ($1+0.5):($2+0.5):($3>0 ? log10($3) : NaN) w image notitle
unset multiplot

