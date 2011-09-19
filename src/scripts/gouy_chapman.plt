#!/usr/local/bin/gnuplot -persist
# Gnuplot script file for plotting data in file "ipbs_solution.dat"
# against the analytical solution of an infinite charged plane
# (Gouy-Chapman model)

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Electrostatic potential in the Gouy-Chapman model" #font "Helvetica,12"
set xlabel "distance {/Times-Italic d}/{/Symbol l} from the boarder of the plane" #font "Helvetica,10"
set ylabel 'potential {/Symbol f} ({/Times-Italic k}_B {/Times-Italic T} / {/Times-Italic e})' # font "Helvetica,14"
#set key left bottom

set log x
charge_density = .88
z = 1.
lambda = 1.
phi_0 = 2. / z * asinh (charge_density * lambda * z *2.*pi)
f(x) = 2. * log ( (1+tanh(phi_0/4.)*exp(-x / lambda)) / (1-tanh(phi_0/4.)*exp(- x / lambda)))
cut(r,y) = (r>4?0:y)
cut2(x,z) = (abs(x)>3?0:z)

set terminal postscript eps color enhanced "Helvetica" 16
set out 'wall.eps'
plot [9e-3:12][2e-2:] './wall2d.dat' u ($1-.99):(-cut($2,$3)) ev 5 title "cylindrical symmetry (disc)", \
'./wall2dcart.dat' u ($1-.99):(-cut($2,$3)) ev 5 title "cartesian symmetry (stick)",'./wallflat.dat' u 2:(-cut2($1,$3)) ev 20 title "cartesian symmtry, flat wall",\
f(x) lt -1 title "Gouy-Chapman potential"

#plot "../ipbs_solution.dat" u 1:(-$3) title 'solution using iPBS algorithm' pt 1, "../ipbs_comparison.dat" u 1:(-$4) title 'constant surface charge' pt 1, f(x) title 'Gouy-Chapman model' ls 3 lw 3
