N = 18
set xrange [-0.5:N-0.5]
set yrange [-0.5:N-0.5]
set size square
set tics out scale 0.1
unset key
set cbtics in
set xtics nomirror
set ytics nomirror
set xtics rotate 90
set term wxt 0 size 1100,900
set object 1 rectangle from graph 0,0 to graph 1,1 fc rgb 'black' lw 0 # black background
set multiplot layout 2,3

# compare p_{3}
vpName(vp) = vp==75 ? 'R2' : (vp==68 ? 'HS-39' : (vp==99 ? 'Hbq1a' : (vp==110 ? 'HS+44' : (vp==112 ? 'HS+48' : (vp==91 ? 'Hba-a1' : NaN)))))
gs(vp) = vp==75 ? 32151060 : (vp==68 ? 32137176 : (vp==99 ? 32199491 : (vp==110 ? 32220918 : (vp==112 ? 32224323 : (vp==91 ? 32182969 : NaN)))))
ge(vp) = vp==75 ? 32151883 : (vp==68 ? 32137426 : (vp==99 ? 32200139 : (vp==110 ? 32221720 : (vp==112 ? 32227298 : (vp==91 ? 32183821 : NaN)))))
f(vp) = sprintf("ERY_alpha/N150-F0-L1-M4-r0.K_fit.q3.vp%d-%d.txt-soi2-diff", gs(vp), ge(vp))
g(vp) = sprintf("ERY_alpha/N150-F1-L3-M4-r0.K_fit.q3.vp%d-%d.txt-soi2-diff", gs(vp), ge(vp))
set cbrange [-0.0005:+0.0075]; set palette defined(0'blue', 0.0625'white', 1'red')
set title sprintf("p_{Erythroid} - p_{ES} [%s]", vpName(75))
plot f(75) matrix rowheader columnheader u 1:2:3 w image
set title sprintf("p_{Erythroid} - p_{ES} [%s]", vpName(68))
plot f(68) matrix rowheader columnheader u 1:2:3 w image
set title sprintf("p_{Erythroid} - p_{ES} [%s]", vpName(91))
plot f(91) matrix rowheader columnheader u 1:2:3 w image
set title sprintf("p_{Erythroid} - p_{ES} [%s]", vpName(99))
plot f(99) matrix rowheader columnheader u 1:2:3 w image
set title sprintf("p_{Erythroid} - p_{ES} [%s]", vpName(110))
plot f(110) matrix rowheader columnheader u 1:2:3 w image
set title sprintf("p_{Erythroid} - p_{ES} [%s]", vpName(112))
plot f(112) matrix rowheader columnheader u 1:2:3 w image
unset multiplot
