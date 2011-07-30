# Gnuplot script file for comparing iPBS obtained results
# to Debye Hueckel theory (see e.g. debye23a)

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "iPBS solution fit to Debye-Hueckel potential" font "Helvetica,16"
set xlabel "radial coordinate" font "Helvetica,12"
set ylabel "potential (reduced units)" font "Helvetica,12"

lambda = 1.0
kappa = 1 / lambda
#charge = 31.416
charge = 188.5
bjerrum = 0.7
radius = 5
A = charge * bjerrum * exp(kappa*radius) / (1 + kappa*radius)
f(x)= A * exp ( - kappa * x ) / x
set pointsize 1

plot [5:50] "../ipbs_step_7.dat" u (sqrt($1*$1+$2*$2)):(-$3) title 'iPBS solution' pt 5 ps 2, "../reference_step_0.dat" u (sqrt($1*$1+$2*$2)):(-$3) title 'finite element solution' pt 3 ps 1, f(x) title 'debye hueckel potential' ls 3 lw 3
