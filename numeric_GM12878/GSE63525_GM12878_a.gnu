ddir = 'GSE63525_GM12878_chr5-90000-100000kb_res50kb'
lg(f,l,m) = sprintf("%s/N200-F%d-L%d-M%d-rx.log",ddir,f,l,m)

mlw = 2
xsf = 60.0
set logscale x 10
set logscale y 10
set xrange []
set ytics add('10^{2}'100,'10^{3}'1000,'10^{4}'10000,'10^{-2}'0.01,'10^{-3}'0.001)
ya(l) = l==0 ? 0.2 : (l==1 ? 300 : (l==2 ? 0.01 : 0.02))
yb(l) = l==0 ? 50 : (l==1 ? 10000 : (l==2 ? 0.1 : 40))
set xlabel 'Time [min]'
unset key
set multiplot layout 2,2
do for [ldx=0:3] {
  set yrange [ya(ldx):yb(ldx)]
  set ylabel sprintf("L_{%d}",ldx) offset 2,0
  plot for [m=1:4] lg(0,ldx,m) u ($7/xsf):1 every ::2 w l ls m lw mlw t sprintf("F%d-M%d",0,m), \
       for [m=1:4] lg(1,ldx,m) u ($7/xsf):1 every ::2 w l ls m lw mlw dt 3 t sprintf("F%d-M%d",1,m)
}
unset multiplot
