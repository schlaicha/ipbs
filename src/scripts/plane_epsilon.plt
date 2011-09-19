#!/usr/local/bin/gnuplot -persist
# Gnuplot script file for plotting data in file "plane_analytic.dat"
# against the analytical solution of two infinite charged planes
# both with constant surface charge and iterated chare

reset
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "Electrostatic potential between two charged planes with dielectric missmatch" #font "Helvetica,12"
set xlabel "{/Times-Italic x}/{/Symbol l}, where the origin is located in the center between both planes" #font "Helvetica,10"
set ylabel 'potential {/Symbol f} ({/Times-Italic k}_B {/Times-Italic T} / {/Times-Italic e})' # font "Helvetica,14"
#set key left bottom

charge_density = 0.00795
kappa = 1
dist = 0.8
dist2 = 5

eps_a = 1
eps_b = 1

eps_bad = 1.5

f(x) = 4*pi*charge_density/kappa/eps_a*((cosh(kappa*dist)*eps_a*(-dist2+dist1)*kappa-2*eps_bad))/(eps_a*(-dist2+dist))*cosh(kappa*x)/cosh(kappa*dist)
g(x) = 4*pi*charge_density/kappa*(2-(2+kappa*dist2)/(1+kappa*dist2+1/tanh(kappa*dist)))*exp(-kappa*(x-dist-dist2))

h(x) = 4*pi*charge_density/kappa/eps_a*((2*eps_b+kappa*dist2*eps_a)/(1+(eps_b+eps_a*kappa*dist2)*tanh(kappa*dist)))*cosh(kappa*x)/cosh(kappa*dist)
i(x) = 4*pi*charge_density/kappa*(2-(2+kappa*dist2)/(1+kappa*dist2+1/tanh(kappa*dist)))*exp(-kappa*(x-dist-dist2))
cut(r,y) = (r>1?-10:y)

set terminal postscript eps color enhanced "Helvetica" 16
set out 'plane_epsilon.eps'
#plot [:dist2+5][0:.2] './without/ipbs_solution.dat' u  ($1>0?$1:(-$1)):(cut($2,$3)) lt 0 pt 7 title "",'./with_bad/ipbs_solution.dat' u  ($1>0?$1:(-$1)):(cut($2,$3)) lt 1 pt 7 title "",\
#    f(x) ls 1 title "", g(x) lt 1 title "bad dielectric", h(x) lt -1 title "" , i(x) ls -1 title "without contrast"
 plot [:8][0.001:.25] './with_good/ipbs_solution.dat' u ($1>0?$1:-$1):(cut($2,$3)) ls 1 title "", 0.1274*cosh(x) ls 1 title "", .1336*exp(-x+3.5) ls 1 title "{/Symbol e}_1 / {/Symbol e}_2 = 1/10 ",\
 './without/ipbs_solution.dat' u ($1>0?$1:-$1):(cut($2,$3)) ls 2 title "", .1556*cosh(x) ls -1 lt 2 title "", .1188*exp(-x+3.5) ls 2 title "{/Symbol e}_1 / {/Symbol e}_2 = 1 ",\
 "./with_bad/ipbs_solution.dat" u ($1>0?$1:-$1):(cut($2,$3)) lt 3 title "", .1/tanh(0.5)*cosh(x)/cosh(0.5) lt 3 title "", 0.1*exp(-x+3.5) lt 3 title "{/Symbol e}_1 / {/Symbol e}_2 = 80"
