#!/usr/local/bin/gnuplot -persist
# Gnuplot script file for comparing iPBS obtained results
# to Gouy-Chapman (gouy19a, chapman13a) and Debye Hueckel theory (see e.g. debye23a)

set  autoscale                        # scale axes automatically
set encoding utf8
unset log
set log
unset label                            # remove any previous labels
set xtic auto                         # set xtics automatically
set ytic auto                          # set ytics automatically

lambda1 = 0.4
lambda2 = .1
lambda3 = .2
lambda4 = .125
lambda5 = 1.
kappa1 = 1./lambda1
kappa2 = 1./lambda2
kappa3 = 1./lambda3
kappa4 = 1./lambda4
kappa5 = 1./lambda5
bjerrum = 1.
radius1=2.2
radius2=50.
radius3=20.
radius4=100.
radius5=1.
f1(x)= 4*pi*bjerrum*x*(radius1/(1+kappa1*radius1))
f2(x)= 4*pi*bjerrum*x*(radius2/(1+kappa2*radius2))
f3(x)= 4*pi*bjerrum*x*(radius3/(1+kappa3*radius3))
f4(x)= 4*pi*bjerrum*x*(radius4/(1+kappa4*radius4))
f5(x)= 4*pi*bjerrum*x*(radius5/(1+kappa5*radius5))
g1(x)=2*asinh(2.*pi*bjerrum*lambda1*x)
g2(x)=2*asinh(2.*pi*bjerrum*lambda2*x)
g3(x)=2*asinh(2.*pi*bjerrum*lambda3*x)
g4(x)=2*asinh(2.*pi*bjerrum*lambda4*x)
g5(x)=2*asinh(2.*pi*bjerrum*lambda5*x)

set terminal postscript eps color enhanced "Helvetica" 16
set out 'surface_potential_gouy.eps'

set title "iPBS solution for {/Symbol k}{/Times-Italic a} >> 1 \n compared to Gouy-Chapman theory (solid) \n \
          and the linear Debye-Hückel surface potential (dashed)" # font "Helvetica,16"
set xlabel 'surface charge density {/Symbol s}' # font "Helvetica,14"
set key bottom right
set ylabel 'reduced surface potential {/Symbol f}_{/Times-Italic S} ({/Times-Italic k}_B {/Times-Italic T} / {/Times-Italic e})' # font "Helvetica,14"
plot [.005:60000][.007:40] 'gouy_potential_k01_a01.dat' u 1:(-$2) title '{/Symbol k} = 1 {/Times-Italic a} = 1' pt 7 lc 1, f5(x) lw 1 lt 2 lc 1 title '', g5(x) lt 1 lc 1 title '',\
     'gouy_potential_k10_a05.dat' u 1:(-$2) title '{/Symbol k} = 10 {/Times-Italic a} = 5' pt 7 lc 2, f2(x) lw 1 lt 2 lc 2 title '', g2(x) lt 1 lc 2 title '',\
     'gouy_potential_k02.2_a_02.2.dat' u 1:(-$2) title '{/Symbol k} = 2.2, {/Times-Italic a} = 2.2' pt 7 lc 3, f1(x) lw 1 lt 2 lc 3 title '', g1(x) lt 1 lc 3 title '' 
     #'gouy_potential_k05_a20.dat' u 1:(-$2) title '{/Symbol k} = 5 {/Times-Italic a} = 20' pt 7 lc 3, f3(x) lw 1 lt 9 lc 3 title '', g3(x) lt 1 lc 3 title '',\
     #'gouy_potential_k08_a100.dat' u 1:(-$2) title '{/Symbol k} = 8 {/Times-Italic a} = 100' pt 7 lc 4, f4(x) lw 1 lt 9 lc 4 title '', g4(x) lt 1 lc 4 title '',\
     #'gouy_potential_k01_a01.dat' u 1:(-$2) title '{/Symbol k} = 1 {/Times-Italic a} = 1' pt 7 lc 5, f5(x) lw 1 lt 9 lc 5 title '', g5(x) lt 1 lc 5 title ''
