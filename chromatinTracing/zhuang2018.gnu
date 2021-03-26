set term wxt 0 size 500,600
set multiplot layout 3,2
load '../tino/moreland.pal'

# p_{ij}
dz = 'GSE63525_IMR90_chr21-29370000-30600000_res30kb'
fx = sprintf("%s/N41-F0-L1-M4-r0.P_fit.vm",dz)
gx = sprintf("%s/N41-F1-L3-M4-r0.P_fit.vm",dz)
N = 41
set xrange [-0.5:N-0.5]
set yrange [-0.5:N-0.5]
set size square
set tics front out scale 0.3
unset key
set palette defined(0'blue', 0.5'white', 1'red')
set cbrange [-2.5:-0.5]; set cbtics in;
plot fx matrix u 1:2:($3>0 ? log10($3) : NaN) w image
plot gx matrix u 1:2:($3>0 ? log10($3) : NaN) w image

# p_{jk|ij} of ctcf triplets
dz = 'GSE63525_IMR90_chr21-29370000-30600000_res30kb'
fx = sprintf("%s/N41-F0-L1-M4-r0.K_fit.zhuang2018.txt",dz)
gx = sprintf("%s/N41-F1-L3-M4-r0.K_fit.zhuang2018.txt",dz)
set xrange [0:370]
set yrange [-0.02:0.28]
set xlabel 'CTCF triplet index'
set ylabel 'Contact probability'
set size ratio 0.7
set tics front in
set xtics 0,100,500
set ytics 0,0.2,1
set key left top samplen 0.3
set border 3; set xtics nomirror; set ytics nomirror
mlw = 2
mps = 0.5
plot fx u 9:7 w p ls 1 lc rgb 'red' pt 7 ps mps t '+', \
     fx u 9:8 w p ls 1 lc rgb 'blue' pt 7 ps mps t '-', \
     fx u 9:6 w l ls 8 lw mlw lc rgb 'cyan' t '0'
plot gx u 9:7 w p ls 1 lc rgb 'red' pt 7 ps mps t '+', \
     gx u 9:8 w p ls 1 lc rgb 'blue' pt 7 ps mps t '-', \
     gx u 9:6 w l ls 8 lw mlw lc rgb 'cyan' t '0'

# p_{jk|ij} of triplets
dz = 'GSE63525_IMR90_chr21-29370000-30600000_res30kb'
fx = sprintf("%s/N41-F0-L1-M4-r0.K_fit.zhuang2018.txt",dz)
gx = sprintf("%s/N41-F1-L3-M4-r0.K_fit.zhuang2018.txt",dz)
set xrange [0:11000]
set yrange [0.001:1]
set xlabel 'triplet index'
set ylabel 'Contact probability'
set size ratio 0.7
set tics front in
set xtics 0,4000,12000
set key left top samplen 0.3
set logscale y 10
set ytics ('10^{-1}'0.1,'10^{-2}'0.01,'10^{-3}'0.001)
set border 3; set xtics nomirror; set ytics nomirror
mlw = 2
mps = 0.5
fq = 5
plot fx u ($0*fq):7 every fq w p ls 1 lc rgb 'red' pt 7 ps mps t '+', \
     fx u ($0*fq):8 every fq w p ls 1 lc rgb 'blue' pt 7 ps mps t '-', \
     fx u 0:6 w l ls 8 lw mlw lc rgb 'cyan' t '0'
plot gx u ($0*fq):7 every fq w p ls 1 lc rgb 'red' pt 7 ps mps t '+', \
     gx u ($0*fq):8 every fq w p ls 1 lc rgb 'blue' pt 7 ps mps t '-', \
     gx u 0:6 w l ls 8 lw mlw lc rgb 'cyan' t '0'

unset multiplot
