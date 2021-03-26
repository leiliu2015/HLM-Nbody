ddir = 'GSE63525_GM12878_chr5-90000-100000kb_res50kb'
set size square
set multiplot layout 1,2

fa(f,l,m,r) = sprintf("%s/N200-F%d-L%d-M%d-r%d.P_fit.vp",ddir,f,l,m,r)
set logscale x 10
set logscale y 10
set xlabel 's [Mb]'
set ylabel 'P(s)' offset 3,0
mlw = 2
N = 200
rs = 0.05
set xrange [rs:200*rs]
set yrange [5e-5:0.4]
set ytics add('0.1'0.1, '10^{-2}'0.01, '10^{-3}'0.001, '10^{-4}'0.0001, '10^{-5}'0.00001)
lab(i) = i==0 ? 'Hi-C' : sprintf("L_{%d}",i)
set key samplen 0.3
plot fa(1,1,4,0) u ($1*rs):2 w p ls 4 lw mlw pt 6 t lab(0), \
     fa(1,1,4,0) u ($1*rs):3 w l ls 1 lw mlw t lab(1), \
     fa(1,2,4,0) u ($1*rs):3 w l ls 8 lw mlw dt (1,.75) t lab(2), \
     fa(1,3,4,0) u ($1*rs):3 w l ls 3 lw mlw dt (2,1) t lab(3)
unset xrange
unset yrange
unset logscale
unset xlabel
unset ylabel

load '../tino/jet.pal'
fb(f,l,m,r) = sprintf("%s/N200-F%d-L%d-M%d-r%d.P_fit.vm",ddir,f,l,m,r)
N = 200
rs= 0.05
set xrange [-0.5:N]
set yrange [-0.5:N]
set size square
set tics front out scale 0.3
i2g(i) = 90.0+i*rs
set xtics 0,40,200 nomirror; set mxtics 2;
set ytics 0,40,200 nomirror; set mytics 2;
do for [i=0:200:40] {
  set xtics add(sprintf("%d",i2g(i)) i)
  set ytics add(sprintf("%d",i2g(i)) i)
}
unset key
set cbrange [-3:-0.5]; set cbtics -5,1,0 in scale 0.5; set mcbtics 5
do for [i=-3:0:1] {
  set cbtics add(sprintf("10^{%d}",i) i)
}
set cbtics offset -0.5,0
set label 1 'Hi-C' at first 5,90 textcolor rgb 'white' front 
set label 2 'L_{1}' at first 82,10 textcolor rgb 'white' front 
plot fb(1,1,4,0) matrix u 1:2:($3>0? log10($3) : NaN) w image

unset multiplot
