#!/usr/bin/python

import numpy
import scipy
from scipy import interpolate
import pylab
from getData import *
import glob
import os

def dens(phi):
  kappa=0.025
  bjerrum=.96
  dens=numpy.sinh(-phi)
  rho=kappa**2/(4*numpy.pi*bjerrum)
  rho = rho * dens
  return(rho)

savefile=open("charge2.dat", "w")

radius=40

#distances=[11., 13., 17., 21.2, 22.46, 23.82, 25.24, 26.76, 28.36, 30.06, 31.86, 33.78, 35.8, 37.96, 40.24, 42.64, 45.2, 47.92, 50.8, 53.84, 57.08, 60.5, 64.14, 67.98, 72.06, 76.38, 80.96, 85.82, 98.]

path=""
filename="ipbs_solution.vtu"

xmin=-400.
xmax=400.
ymin=0.
ymax=400.
resolution=2

x,y,sol=getData(filename, xmin, xmax, ymin, ymax, resolution, resolution)
phi = scipy.interpolate.RectBivariateSpline(x,y,sol, kx=3, ky=3)
numpy.savetxt("sol.txt", sol.transpose())

#distances=[ 180. ]

#for dist in distances:
for dist in range(45, 395, 10):
  dist=float(dist)
  #now calculate the radial distribution
  alpha=0
  alphasteps= 10
  intsteps = 1800
  alpha_step = numpy.pi/alphasteps
  intlen = dist - radius
  #print intlen
  alpha = 0
  int_step=intlen/intsteps
  data=numpy.zeros(alphasteps)
  alphas=numpy.zeros(alphasteps)
    
  ofile=open("tester.dat", "w")
  for i in range(alphasteps):
    intpos = numpy.array([radius*numpy.cos(alpha),radius*numpy.sin(alpha)])
    dir =  numpy.array([numpy.cos(alpha),numpy.sin(alpha)]) * int_step
    #print 'alpha = ', alpha, ' intpos = ', intpos, ' dir = ', dir
    n = 0
    for j in range(intsteps):
      r = numpy.sqrt((intpos[0])**2+intpos[1]**2)
      phi1 = phi(intpos[0],intpos[1])[0][0]
      d = dens(phi1)
      if j == 0:
        d/=2
      n = n+d*r**2*int_step*2*numpy.pi
      ofile.write(  str(r) + " " + str(phi1) + " " + str(d) +" " +  str(n) + "\n")
      #ofile.write( str(intpos[0]) + " " + str(intpos[1]) + " " + str(phi(intpos[0],intpos[1])[0][0]) + "\n")
      #print intpos, n
      intpos += dir
    alphas[i]=alpha
    data[i]=n
    alpha += alpha_step
    n=0
  ofile.close()
  pylab.plot(alphas, -data, "o-")
  
  sum=0
  for i in alphas:
    sum += data[i] * numpy.sin(i) * alpha_step
    #savefile.write(str(i) + " " + str(data[i]) + " " + str(sum) + "\n")
  print sum
  
  savefile.write( str(dist) + " " + str(sum) + "\n")
savefile.close()


pylab.show()
