N = 20
vp= 7
fdx = 0
set xrange [-0.5:N-0.5]
set yrange [-0.5:N-0.5]
set size square
set tics front out scale 0.3
set cbtics in
unset key

set multiplot layout 2,3

load '../tino/piyg.pal'
set title 'p(i,j)'
fx = sprintf("g%d.%s2.txt",N,fdx==0 ? 'q' : 'p')
plot fx matrix u 1:2:($3>0 ? log10($3) : NaN) w image
set title 'p(i,j,7)'
fx = sprintf("g%d.%s3.txt",N,fdx==0 ? 'q' : 'p')
plot fx index vp matrix u 1:2:($3>0 ? log10($3) : NaN) w image
set title 'p(i,j,6,15)'
if (fdx==0) {
  fx = sprintf("g%d.%s4.txt",N,fdx==0 ? 'q' : 'p')
  plot fx index 107 matrix u 1:2:($3>0 ? log10($3) : NaN) w image
} else {
  fx = sprintf("g%d.xyz.%s4.txt",N,fdx==0 ? 'q' : 'p')
  plot fx index 107 matrix u 1:2:($3>0 ? log10($3) : NaN) w image
}

load '../tino/moreland.pal'
set object 1 rect from first -0.5,vp-0.5 to first N-0.5,vp+0.5 fc rgb 'black' fs noborder front
set object 2 rect from first vp-0.5,-0.5 to first vp+0.5,N-0.5 fc rgb 'black' fs noborder front
set title 'p(i,7,k)/p(i,7)p(7,k)'
fx = sprintf("g%d.tamm-F%d.txt",N,fdx)
plot fx index vp matrix u 1:2:3 w image

set title 'H(i,j,7)'
fx = sprintf("g%d.tanay-F%d.vp%02d.T2.hs.txt",N,fdx,vp)
plot fx matrix u 1:2:3 w image

set title 'A(i,j,7)'
set cbrange [-6:6]
fx = sprintf("g%d.xyz.deLaat-F%d.vp%02d.zm",N,fdx,vp)
plot fx matrix u 1:2:3 w image

unset multiplot

