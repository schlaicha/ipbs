#!/usr/bin/python
# -*- encoding: utf8 -*-

# show the absolute difference between the iPBS solution of a single sphere in reduced symmetry
# and the constant surface flux (Neumann) reference solution

import sys
import matplotlib
matplotlib.use('pdf')
import numpy
import pylab
from pylab import savefig
from scipy import *
from matplotlib.delaunay import *
from matplotlib import rc

def readfile(file):
    x, y ,ref, sol, dif, rel = loadtxt(file, unpack=True) 
    return x, y, ref, sol, dif, rel

def defgrid(x,y):
    paso = 1000
    intervalx = (x.max()-x.min())/paso
    intervaly = (y.max()-y.min())/paso
    xvalues = arange(x.min(),x.max(),intervalx)
    yvalues = arange(y.min(), y.max(),intervaly)
    xi, yi = meshgrid(xvalues,yvalues)
    return xi,yi

def interpolacion(x,y,z,xi,yi):
    # triangulate data
    tri = Triangulation(x,y)
    # interpolate data
    interp = tri.nn_interpolator(z)
    zi = interp(xi,yi)
    return zi

#do the acutal plotting

def plotear(xi,yi,zi):
    # mask inner circle
    interior = sqrt((xi**2) + (yi**2)) < 5.0 
    zi[interior] = ma.empty
    pylab.figure(figsize=(16,10))
    levels = [9e-4, 9.6e-4, 9.8e-4, 1.0e-3, 1.1e-3, 1.2e-3, 1.4e-3, 1.15e-3, 1.05e-3, 1.005e-3, 9.9e-4, 1.3e-3, 1.03e-3, 1.05e-5, 1.04e-3, 1.09e-3]
    labels = [1.2e-3, 1.4e-3, 1.15e-3, 1.3e-3, 1.03e-3, 1.04e-3, 1.09e-3, 1.1e-3]
    #levels = numpy.linspace(9.e-4, 1.1e-3, num=18, endpoint=True)
    #levels = np.logspace(-8,-1, num=20, endpoint=True, base=10)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    rc('text', usetex=True)
    #rc('font', family='serif')
    CS = pylab.contour(xi,yi,zi, levels, colors='black', linewidths=1.0, lynestiles='solid')
    pylab.xlim((-6.0,6.0))
    pylab.ylim((0, 6.0))
    pylab.clabel(CS,labels, fmt='%1.2e', rightside_up='true', colors='black', inline_spacing=5, fontsize=20, inline=1)
    #pylab.colorbar(CS, spacing='proportional', format='%1.2e')
    pylab.title('iPBS relative error', fontsize=25)
    pylab.xlabel(r'$\displaystyle z/\lambda_D$',fontsize=20)
    pylab.ylabel(r'$\displaystyle r/\lambda_D$',fontsize=20)
    #pylab.show()
    savefig('ipbs_rel')
    return xi, yi, zi

def main():
    filename = sys.argv[1]
    x, y, ref, sol, dif, rel = readfile(filename)
    z = abs(rel)
    xi, yi = defgrid(x,y)
    zi = interpolacion(x,y,z,xi,yi)
    return plotear(xi,yi,zi)


if __name__ == "__main__":
    if len (sys.argv) == 2:
        xi, yi, zi=main ()
    else:
        print "How to use: %s data_file" % sys.argv[0]

