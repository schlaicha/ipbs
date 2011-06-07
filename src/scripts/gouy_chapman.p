# Gnuplot script file for plotting data in file "ipbs_solution.dat"
# against the analytical solution of an infinite charged plane
# (Gouy-Chapman model)

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "electrostatic potential in the Gouy-Chapman model" font "Helvetica,12"
set xlabel "distance from plane" font "Helvetica,10"
set ylabel "potential (reduced units)" font "Helvetica,10"

set log x
charge_density = 0.3
z = 1
lambda = 1
phi_0 = 2 / z * asinh (charge_density * lambda * z / sqrt(8))
f(x) = 2 * log ( (1+tanh(phi_0/4)*exp(- sqrt(1) / lambda * x)) / (1-tanh(phi_0/4)*exp(- sqrt(1) / lambda * x)))

plot "ipbs_solution.dat" u 2:3 title 'solution using iPBS algorithm' pt 1, "reference_step_0.dat" u 2:3 title 'solution using standard FE with constant flux' pt 1, f(x) title 'Gouy-Chapman model' ls 3 lw 3
