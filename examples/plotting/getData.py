#!/usr/bin/python

import numpy
import scipy
from scipy import interpolate
import vtktools
#import pylab

def getData(file, xmin=float('nan'), xmax=float('nan'), ymin=float('nan'), ymax=float('nan'), step_x=1,step_y=1):
  u=vtktools.vtu(file)

  print numpy.isnan(xmin), numpy.isnan(xmax), numpy.isnan(ymin), numpy.isnan(ymax)

  if numpy.isnan(xmin): 
    xmin=u.ugrid.GetBounds()[0]
    print 'xmin = ', xmin
  if numpy.isnan(ymin):
    ymin=u.ugrid.GetBounds()[2]
    print 'ymin = ', ymin
  if numpy.isnan(xmax):
    xmax=u.ugrid.GetBounds()[1]
    print 'xmax = ', xmax
  if numpy.isnan(ymax):
    ymax=u.ugrid.GetBounds()[3]
    print 'ymax = ', ymax
  
  Xlist = numpy.arange(xmin,xmax,step_x)# x coordinates
  Ylist = numpy.arange(ymin,ymax,step_y)# y coordinates
  [X0,Y0] = scipy.meshgrid(Xlist,Ylist)
  X0=X0.transpose()
  Y0=Y0.transpose()
  print Xlist.shape, Ylist.shape, X0.shape, Y0.shape
  Z0 = 0.0*Y0 # This is 2d so z is an array of zeros.
  X = numpy.reshape(X0,(numpy.size(X0),))
  Y = numpy.reshape(Y0,(numpy.size(Y0),))
  Z = numpy.reshape(Z0,(numpy.size(Z0),))
  pts = zip(X,Y,Z)
  pts = vtktools.arr(pts)
  # create arrays of velocity and temperature values at the desired points
  sol = u.ProbeData(pts, 'solution')
  print sol.shape, Xlist.shape, Ylist.shape
  #temperature_structured = u.ProbeData(pts, 'Temperature')
  numpy.savetxt("pts.dat", zip(X,Y,sol))

  sol = sol.reshape((numpy.size(Xlist),numpy.size(Ylist)))
  x=[]
  y=[]
  z=[]
  for i in range(len(Xlist)):
    for j in range(len(Ylist)):
      #print i,j
      x.append(X0[i,j])
      y.append(Y0[i,j])
      z.append(sol[i,j])
  numpy.savetxt("pts2.dat", numpy.array(zip(x, y, z)))

  #data = scipy.interpolate.RectBivariateSpline(Xlist,Ylist,sol)
  
  # return Xlist, Ylist, sol
  return Xlist,Ylist,sol
