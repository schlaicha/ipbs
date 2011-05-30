# Gnuplot script file for plotting data in file "reference.dat"
# and comparing it to the analytical solution

set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set title "iPBS solution fit to Debye-Hueckel potential" font "Helvetica,16"
set xlabel "radial coordinate" font "Helvetica,12"
set ylabel "potential (reduced units)" font "Helvetica,12"

#set key 0.01,100
#set label "Yield Point" at 0.003,260
#set arrow from 0.0028,250 to 0.003,280
set xr [5:15]
set yr [1e-5:]
set log y
f(x)=a*exp(-x)/x
fit f(x) 'reference_step_7.dat' u (sqrt($1*$1+$2*$2)):(-$3) via a
set pointsize 1

plot "reference_step_7.dat" u (sqrt($1*$1+$2*$2)):(-$3) title 'finite element solution' pt 7, f(x) title 'debye hueckel fit' ls 3 lw 2
