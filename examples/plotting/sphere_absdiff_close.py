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
    paso = 100
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
    interior = sqrt((xi**2) + (yi**2)) < 1.0 
    zi[interior] = ma.zeros
    pylab.figure(figsize=(16,10))
    #levels = [log(1.625e-3) ,log(3.75e-3), log(7.5e-3), log(1.5e-4), log(3e-4), log(6e-4), log(1.2e-3), log(2.4e-3), log(4.8e-3), log(9.6e-3)]
    levels = numpy.linspace(1e-9, 1e-3, 30, endpoint=True)
    #levels = numpy.logspace(1e-10, 1e-4, 20, endpoint=True)
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    rc('text', usetex=True)
    #rc('font', family='serif')
    CS = pylab.contour(xi,yi,zi, levels, colors='black', linewidths=1.0, lynestiles='solid')
    pylab.xlim((-10,10))
    pylab.ylim((0, 10))
    pylab.clabel(CS,levels[0::2], fmt='%1.2e', rightside_up='true', colors='black', inline_spacing=5, fontsize=20, inline=1)
    pylab.title('iPBS absolute error', fontsize=25)
    pylab.xlabel(r'$\displaystyle z/\lambda_D$',fontsize=20)
    pylab.ylabel(r'$\displaystyle r/\lambda_D$',fontsize=20)
    #pylab.show()
    savefig('ipbs_abs_close')
    return xi, yi, zi

def main():
    filename = sys.argv[1]
    x, y, ref, sol, dif, rel = readfile(filename)
    #z = log(abs(dif))
    z = abs(dif)
    x_l = x/5.
    y_l = y/5.
    xi, yi = defgrid(x_l,y_l)
    zi = interpolacion(x_l,y_l,z,xi,yi)
    plotear(xi,yi,zi)

if __name__ == "__main__":
    if len (sys.argv) == 2:
        main ()
    else:
        print "How to use: %s data_file" % sys.argv[0]
