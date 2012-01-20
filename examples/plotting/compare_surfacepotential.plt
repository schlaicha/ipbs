#!/usr/local/bin/gnuplot -persist
# Gnuplot script file for comparing iPBS obtained results
# to Debye Hueckel theory (see e.g. debye23a)

set  autoscale                        # scale axes automatically
set encoding utf8
# unset log
set log
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically

charge_1 = 1
charge_2 = .1
charge_3 = 10
f(x)= charge_1 / (1+x)
g(x)= charge_2 / (1+x)
h(x)= charge_3 / (1+x)

set terminal postscript eps color enhanced "Helvetica" 16
set out 'surface_potential.eps'

set multiplot;
set origin 0.0,0.0;
set size 1,1;  

set title "Surface potentials of iPBS and Debye-Hückel theory" # font "Helvetica,16"
set xlabel '{/Symbol k} {/Italic a}' # font "Helvetica,14"
set ylabel 'surface potential {/Symbol f}_{/Italic S} ({/Times-Italic k}_B {/Times-Italic T} / {/Times-Italic e})' # font "Helvetica,14"
plot [][8e-4:12] "../data/surface_potential_Z0.1.dat" u 1:(-$2) pt 7 title '' , g(x) ls 2 lw 2 lt 2 title '{/Times-Italic Z} = 0.1', \
           "../data/surface_potential_Z1.dat" u 1:(-$2) pt 7 title '' , f(x) ls 3 lw 2 lt 3 title '{/Times-Italic Z} = 1', \
           "../data/surface_potential_Z10.dat" u 1:(-$2) ls -1 pt 2 ps 1 title '' , h(x) ls 4 lw 2 lt 4 title '{/Times-Italic Z} = 10'

unset log
set log x
set ytic (0.094, 0.099)
set origin 0.12, 0.12;
set xlabel ''
set ylabel ''
set title ''
set size 0.4, 0.35;
plot [1e-2:1e-1]  "../data/surface_potential_Z0.1.dat" u 1:(-$2) pt 7 title '' , g(x) ls 2 lw 2 lt 2 title '' 

unset multiplot
