#!/usr/local/bin/gnuplot -persist
# Gnuplot script file for comparing iPBS obtained results
# to Debye Hueckel theory (see e.g. debye23a)

set  autoscale                        # scale axes automatically
set encoding utf8
unset log
set log
unset label                            # remove any previous labels
set xtic (1,2,5,10)                          # set xtics automatically
set ytic auto                          # set ytics automatically

lambda = 5
kappa = 1./lambda
charge = 1
bjerrum = 1
radius = 5
A = charge * bjerrum * exp(kappa*radius) / (1 + kappa*radius)
f(x)= A * exp ( - kappa * x ) / x

set terminal postscript eps color enhanced "Helvetica" 16

set out 'debye_hueckel.eps'
set multiplot;
set origin 0.0,0.0;
set size 1,1;  
set title "iPBS solution and Debye-Hückel potential" # font "Helvetica,16"
set xlabel 'distance {/Times-Italic x} / {/Symbol l}' # font "Helvetica,14"
set ylabel 'reduced potential {/Symbol f} ({/Times-Italic k}_B {/Times-Italic T} / {/Times-Italic e})' # font "Helvetica,14"
plot [1:10][5e-8:0.2] "../ipbs_solution.dat" u ((sqrt($1**2+$2**2))/lambda):(-$3) every 5 ls 1 pt 2 ps 1 title 'iPBS solution' , f(x*lambda) ls -1 lw 3 title 'Debye-Hückel potential' 
set origin 0.15,0.1;
set size 0.5,0.5;  
unset log;
set title 'reproduction of the surface potential'
set xtic auto
set xlabel '' 
set ylabel ''
plot [0.9999:1.01][.097:0.1005] "../ipbs_solution.dat" u ((sqrt($1**2+$2**2))/lambda):(-$3) every 10 ls 1 pt 2 ps 1 title '' , f(x*lambda) ls -1 lw 3 title '' 
unset multiplot
