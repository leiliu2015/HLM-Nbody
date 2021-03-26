ddir = 'Bonev_ES_chr8-42100-44500kb_res25kb'
lg(f,l,m) = sprintf("%s/N97-F%d-L%d-M%d-rx.log",ddir,f,l,m)

mlw = 2
xsf = 60.0
set logscale x 10
set logscale y 10
set xrange []
set ytics add('10^{2}'100,'10^{3}'1000,'10^{4}'10000,'10^{-2}'0.01)
ya(l) = l==0 ? 0.01 : (l==1 ? 10 : (l==2 ? 0.005 : 0.007))
yb(l) = l==0 ? 15 : (l==1 ? 3000 : (l==2 ? 0.1 : 20))
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
