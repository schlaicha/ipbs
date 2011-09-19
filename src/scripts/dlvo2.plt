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

s1 = 0.01#*.893
s2 = 0.1#*.86
a = 1.5
l_B = 1
k = 1
f(x)=(4*pi*s1*a**2)**2*(exp(k*a)/(1+k*a))**2*(exp(-k*x)/x**2+k*exp(-k*x)/x)
g(x)=(4*pi*s2*a**2)**2*(exp(k*a)/(1+k*a))**2*(exp(-k*x)/x**2+k*exp(-k*x)/x)

set terminal postscript eps color enhanced "Helvetica" 16
set out 'dlvo2.eps'
set log y
plot 'dlvo_sigma0.01_without.dat' lt 1 pt 1 title "without osmotic pressure", 'dlvo_sigma0.01_with.dat' lt 1 pt 2 title "with osmotic pressure", f(x) title "DLVO force {/Symbol s} = 0.01" lt 1,\
'dlvo_sigma0.1_without.dat' lt 3 pt 1 title "", 'dlvo_sigma0.1_with.dat' lt 3 pt 2 title "", g(x) title "DLVO force {/Symbol s} = 0.1 " lt 3
