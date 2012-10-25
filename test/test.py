#!/usr/bin/python
import subprocess
import numpy
import matplotlib.pyplot as p
import pylab

def runtest(test):
  p = subprocess.Popen(test, shell=True)
  retval = p.wait()
  return retval

def checkmsh(meshname):
  p = subprocess.Popen('[ -f grids/' + str(meshname) + '.msh ]', shell=True)
  if p.wait() == 0:
    p = subprocess.Popen('[ grids/' + str(meshname) + '.geo -ot grids/' + str(meshname) + '.msh ]', shell=True)
    return p.wait()
  else: return 1

def domesh(meshname):
  p = subprocess.Popen('gmsh -2 $(pwd)/grids/' + str(meshname) + '.geo ]', shell=True)
  retval = p.wait()
  return retval

def getData2d(datafile):
  data = numpy.loadtxt(datafile)
  xdata = data[:,0]
  ydata = data[:,1]
  sol = data[:,2]
  return xdata, ydata, sol

def dh(r):
  s = 0.001
  R = 3
  k = 1
  lb = 1
  return (4*numpy.pi*R**2*s*numpy.exp(k*R)/(1+k*R)*numpy.exp(-k*r)/r)

def test_2d_dh():
  print "Testing 2D solution for Debye-Hueckel approx"
  print "===================================================================\n"
  if (checkmsh("sphere") == 0):
    print "using existing meshfile"
  else:
    print "calling GMSH"
    g = domesh("sphere")
    if g == 0:
      print "call was successful..."
    else: 
      print "GMSH failed"
      return 1
    print "now calling ipbs"
  cmd = '$(pwd)/../src/ipbs test2d_sphere.cfg'
  success = runtest(cmd)
  if success == 0:
    print 'IPBS has succesfully finished'
  else:
    print 'IPBS failed'
    return 1
  x, y, sol = getData2d("dh2D_solution.dat")
  radial = numpy.sqrt(x**2 + y**2)
  # plot the IPBS results
  p.plot(radial, sol,'.')
  # analytical approximation
  xvals = pylab.arange(numpy.min(radial),numpy.max(radial),0.02)
  p.plot(xvals, dh(xvals), "r-", linewidth=2)
  p.xlabel("radial coordinate")
  p.ylabel("potential")
  p.title("IPBS solution vs. Debye-Hueckel approx. (2D)")
  dev = numpy.subtract( sol, dh(radial) )
  reldev = dev / dh(radial)
  p.text( 10, .002,"standard deviation: " + str(numpy.std(dev)) )
  p.axes([.45,.45,.35,.35])
  p.plot(radial, reldev, '.')
  p.title("relative deviation")

test_2d_dh()
p.show()
