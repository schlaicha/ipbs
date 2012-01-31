#!/usr/local/bin/gnuplot -persist

reset
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Integrated charge within a spherical shell of radius r around a colloid" #font "Helvetica,12"
set xlabel "radius r / {/Symbol l}" #font "Helvetica,10"
set ylabel 'Q (r) / Ze' # font "Helvetica,14"
set key right bottom

f(x)=-exp(1)/(1+.1*10)*4*pi/10*(1+0.1*x*10)*exp(-.1*x*10)+4*pi/10

set terminal postscript eps color enhanced "Helvetica" 16
set out 'chargecomp.eps'
#set log x
plot f(x)/(4*pi/10) title "Integrated DH charge density" lt -1,\
'charge2.dat' u ($1/10):($2/(4*pi/10)) title "a = 10 nm, {/Symbol s} = 0.001 e/nm^2, {/Symbol l} = 10 nm" lc 1 pt 5,\
 'charge3.dat' u ($1/40):($2/255) w lp title 'a = 40 nm, Z = 255 e, {/Symbol l} = 40 nm' pt 7 lc 3
