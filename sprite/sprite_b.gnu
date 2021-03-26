tdx = 10
vps = 22
vpe = 24
N = 104
gz= 11.69
mlw = 1
mps = 1
set term wxt 0 size 800, 600
set multiplot layout 2, 3

# p_{ij}
unset key
set size square
load '../tino/piyg.pal'
fdx=0; ldx=1; mdx=4
fb = sprintf("TAD%02d/N%d-F%d-L%d-M%d-r0.P_fit.vm",tdx,N,fdx,ldx,mdx)
rs = 0.005
set xrange [0:N]
set yrange [0:N]
i2g(i) = gz+i*rs
set xtics 1,40,N nomirror offset 0,0.4; set mxtics 2;
set ytics 1,40,N nomirror offset 0.8,0; set mytics 2;
do for [i=1:N:40] {
  set xtics add(sprintf("%.1f",i2g(i)) i)
  set ytics add(sprintf("%.1f",i2g(i)) i)
}
set cbrange [-2.:-1.0]; set cbtics -5,1,0 in scale 0.5; set mcbtics 5
set title 'p_{ij}'
plot fb matrix u ($1+0.5):($2+0.5):($1==$2 ? NaN : log10($3)) w image pixels

# p_{ijk}
fa = sprintf("sprite-data/TAD%02d.p3.vp%03d-%03d.txt", tdx, vps, vpe)
fb = sprintf("TAD%02d/N%d-F%d-L%d-M%d-r0.K_fit.q3.vp%03d-%03d.pm", tdx, N, fdx, ldx, mdx, vps, vpe)
set cbrange [-4.3:-2.7]; set cbtics -6,1,0; set cbtics in scale 0.5; set mcbtics 2;
set title 'p_{ijk}'
plot fa matrix u ($1+0.5):($2+0.5):($1<$2 ? $3*1.6/40-4.3 : NaN) w image pixels, \
     fb matrix u ($1+0.5):($2+0.5):($1>$2 ? log10($3) : NaN) w image pixels

# \overline{p}_{ijk}
fb = sprintf("TAD%02d/N%d-F%d-L%d-M%d-r0.K_fit.q3.vp%03d-%03d.em", tdx, N, fdx, ldx, mdx, vps, vpe)
set cbrange [-4.3:-2.7]; set cbtics -6,1,0; set cbtics in scale 0.5; set mcbtics 2;
set title 'p@^{expected}_{ijk}'
plot fb matrix u ($1+0.5):($2+0.5):($1==$2 ? NaN : log10($3)) w image pixels

# Z_{ijk}
fb = sprintf("TAD%02d/N%d-F%d-L%d-M%d-r0.K_fit.q3.vp%03d-%03d.sm", tdx, N, fdx, ldx, mdx, vps, vpe)
set cbrange [0.5:2.5]; set cbtics -6,2,0; set cbtics in scale 0.5; set mcbtics 2;
set palette defined(0'#3b4cc0',0.2'white',1'#b40426')
set title 'Z_{ijk}'
plot fb matrix u ($1+0.5):($2+0.5):($1==$2 ? NaN : $3) w image pixels

# quantilized p_{ijk}
fx = sprintf("TAD%02d/N%d-F%d-L%d-M%d-r0.K_fit.q3.vp%03d-%03d.pm-qt", tdx, N, fdx, ldx, mdx, vps, vpe)
set xrange [0:0.0017]
set yrange [10:28]
set xtics 0,0.001,0.003;
set ytics 0,10,30
set xlabel 'p_{ijk}'
set ylabel 'SPRITE Freq'
set border 3
set size ratio 0.7
unset title
plot fx u 9:18:8:10:17:19 w xyerror lw mlw pt 7 ps mps notitle

# ENPS
set datafile commentschars "%"
set tics nomirror
fx= sprintf("TAD%02d/N%d-F0-L1-M4-r0.K_fit.q3.ENPS", tdx, N)
cr= 0.380468
set size ratio 0.3
set xrange [-0.5:19.5]
set yrange [0.1:2.3]
set xtics rotate 90 scale 0.2
set ytics 0,1,4; set mytics 2
set ylabel 'Z' offset 2, 0
unset xlabel
set border 2
set arrow 1 from graph 0, first 1 to graph 1, first 1 nohead lc rgb 'black' lw 1 back
set obj 2 rectangle from graph 0, first 1-cr to graph 1, first 1+cr
set obj 2 fc rgb '#E6E6E6' fs solid 1.0 noborder back
plot fx u 0:($5> 1+cr? $5 : NaN):6:xticlabels(1) w error lw mlw lc rgb 'red'  pt 7 ps mps notitle, \
     fx u 0:($5<=1+cr? $5 : NaN):6:xticlabels(1) w error lw mlw lc rgb 'blue' pt 7 ps mps notitle

unset multiplot
