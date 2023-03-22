#Common parameters

set xlabel "X"
set ylabel "Y"
set zlabel "H"
set term eps size 5.5,4.5
set autoscale fix
set key off


# 3D plots

set size 1.1,1.1
set pm3d
set hidden3d
unset colorbox
set origin -0.01, -0.05

set output "test1-3d.eps"
splot "test1.txt" using 1:2:5

set output "test2-3d.eps"
splot "test2.txt" using 1:2:5

set output "test3-3d.eps"
splot "test3.txt" using 1:2:5

set output "test4-3d.eps"
splot "test4.txt" using 1:2:5


#2D plots

unset size
unset origin
unset hidden3d
set pm3d map
set colorbox
set size square
set lmargin at screen 0.05;
set rmargin at screen 0.9;
set bmargin at screen 0.1;
set tmargin at screen 0.95;

set output "test1-2d.eps"
splot "test1.txt" using 1:2:5

set output "test2-2d.eps"
splot "test2.txt" using 1:2:5

set output "test3-2d.eps"
splot "test3.txt" using 1:2:5

set output "test4-2d.eps"
splot "test4.txt" using 1:2:5