#!/usr/local/bin/gnuplot -persist
# Gnuplot script file for plotting data in file "plane_analytic.dat"
# against the analytical solution of two infinite charged planes
# both with constant surface charge and iterated chare

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Electrostatic potential between two charged planes" #font "Helvetica,12"
set xlabel "{/Times-Italic x}/{/Symbol l}, where the origin is located in the center between both planes" #font "Helvetica,10"
set ylabel 'potential {/Symbol f} ({/Times-Italic k}_B {/Times-Italic T} / {/Times-Italic e})' # font "Helvetica,14"
#set key left bottom

charge_density = 0.00795
kappa = 1
dist = 1
dist2 = 1
f(x) = 4*pi*charge_density*exp(-(x-dist-+dist))
g(x) = 1./tanh(kappa*dist)*4*pi*charge_density*cosh(x)/cosh(kappa*dist) 
h(x) = 4*pi*charge_density/kappa*((2+kappa*dist2)/(1+(1+kappa*dist2)*tanh(kappa*dist)))*cosh(kappa*x)/cosh(kappa*dist)
i(x) = 4*pi*charge_density/kappa*(2-(2+kappa*dist2)/(1+kappa*dist2+1/tanh(kappa*dist)))*exp(-kappa*(x-dist-dist2))
cut(r,y) = (r>1?-10:y)

set terminal postscript eps color enhanced "Helvetica" 16
set out 'plane_analytic.eps'
plot [:5][0:.15] './plane_analytic.dat' u  ($1>0?$1:(-$1)):(cut($2,$3)) lt 1 pt 7 title "",'./plane_analytic.dat' u  ($1>0?$1:(-$1)):(cut($2,-$4)) lt 0 pt 7 title "",\
    f(x) ls -1 title "", g(x) ls -1 title "constant surface charge", h(x) lt 1 title "" , i(x) lt 1 title "effective charge"
