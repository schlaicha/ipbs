#!/usr/bin/python

import numpy
import scipy
from scipy import interpolate
import pylab
from getData import *

def dens(phi):
  kappa=0.1
  bjerrum=1.0
  dens=numpy.sinh(-phi)
  rho=kappa**2/(4*numpy.pi*bjerrum)
  rho = rho * dens
  return(rho)

radius=10
dist='40.240'
path="two_colloid_dist_" + str(dist)
filename="two_colloid_dist_" + str(dist) + "/ipbs_solution.vtu"
dist=float(dist)

xmin=-dist
xmax=0
ymin=0
ymax=dist/2.
resolution=.125

x,y,sol=getData(filename, xmin, xmax, ymin, ymax, resolution, resolution)

print x.shape, y.shape, sol.shape

phi = scipy.interpolate.RectBivariateSpline(x,y,sol, kx=3, ky=3)
print phi(2,3)
numpy.savetxt("sol.txt", sol.transpose())
#rho=scipy.interpolate.interp2d(x,y,dens(sol), kind='cubic')
#rho=scipy.interpolate.RectBivariateSpline(x,y,sol)
#x = numpy.linspace(xmin,xmax,100)# x coordinates
#y = numpy.linspace(ymin,ymax,100)# y coordinates
#z = solutiondata(x,y)
#print z
print phi(x,y).shape


CSF = pylab.contourf(x.flatten(),y.flatten(),dens(phi(x,y).transpose()))
CS = pylab.contour(x,y,phi(x,y).transpose())
cbar = pylab.colorbar(CSF, format='%g')
#cbar.add_lines(CS)
pylab.title('effective ion density')
pylab.xlabel('z',fontsize=12)
pylab.ylabel('r',fontsize=12)
pylab.figure()


#now calculate the radial distribution
alpha=0
alphasteps=100
intsteps = 100

alpha_step = numpy.pi/alphasteps
intlen = dist/2 - radius
int_step=intlen/intsteps
data=numpy.zeros(alphasteps)
alphas=numpy.zeros(alphasteps)

ofile=open("tester.dat", "w")
for i in range(alphasteps):
  intpos = numpy.array([radius*numpy.cos(alpha)-dist/2,radius*numpy.sin(alpha)])
  dir =  numpy.array([numpy.cos(alpha),numpy.sin(alpha)]) * int_step
  #print 'alpha = ', alpha, ' intpos = ', intpos, ' dir = ', dir
  n = 0
  for j in range(intsteps):
    n = n+dens(phi(intpos[0],intpos[1]))*((intpos[0]+dist/2)**2+intpos[1]**2)*int_step*2*numpy.pi
    ofile.write( str(intpos[0]) + " " + str(intpos[1]) + " " + str(dens(phi(intpos[0],intpos[1]))[0][0]*((intpos[0]+dist/2)**2+intpos[1]**2)*int_step*2*numpy.pi) + "\n")
    #print intpos, n
    intpos += dir
  alphas[i]=alpha
  data[i]=n
  alpha += alpha_step
  n=0
ofile.close()
#pylab.ylim(0,0.5)
pylab.plot(alphas, -data, "o-")

sum=0
for i in alphas:
  sum += data[i] * numpy.sin(i) * alpha_step
print sum

pylab.show()

ofile = open(path+"/radial_density.dat", "w")
for i in range(0, len(alphas)):
  ofile.write(str(alphas[i]) + " " + str(data[i]) + "\n")
ofile.close()

#data=numpy.zeros((alphasteps, intsteps ))
#rs=numpy.zeros(intsteps)
#alpha=0
#for i in range(alphasteps):
#  sum=0
#  r=10
#  alphas[i]=alpha
#  for j in range(intsteps):
#    r = r + int_step
#    intpos = array([radius*cos(alpha) - dist/2,radius*sin(alpha)])
#    sum+= -rho(intpos[0],intpos[1])*r**2 * int_step
#    data[i,j]=sum
#    rs[j] = r
#
#  alpha+=alpha_step
#
#pylab.plot(alphas, data.transpose()[intsteps-1])
#pylab.show()
