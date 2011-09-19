#!/usr/local/bin/gnuplot -persist
# Gnuplot script file for plotting data in file "dlvo_sigma*_force.dat"
# against the results of the DLVO theory

reset
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Force between two identical spherical particles" #font "Helvetica,12"
set xlabel "{/Times-Italic x}/{/Symbol l} (center-center)" #font "Helvetica,10"
set ylabel 'F {/Symbol l} / ({/Times-Italic k}_B {/Times-Italic T})' # font "Helvetica,14"
#set key left bottom

s = 0.00922
a = 1
l_B = 1
k = 1
f(x)=(4*pi*s*a**2)**2*(exp(k*a)/(1+k*a))**2*(exp(-k*x)/x**2+k*exp(-k*x)/x)

set terminal postscript eps color enhanced "Helvetica" 16
set out 'dlvo.eps'
set log y
plot 'dlvo_sigma0.01_force.dat' title "", f(x) title "DLVO force fit" lt -1
