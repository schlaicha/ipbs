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
    interior = sqrt((xi**2) + (yi**2)) < 1.0 
    zi[interior] = ma.zeros
    pylab.figure(figsize=(16,10))
    levels = [log(1e-9), log(5e-9), log(1e-8), log(5e-8), log(1e-7), log(5e-7), log(1e-6), log(5e-6), log(1e-5), log(5e-5), log(1e-4), log(2e-4), log(4e-4), log(8e-4), log(1.2e-3), log(2.4e-3), log(4.8e-3), log(9.6e-3)]
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    rc('text', usetex=True)
    scale = list()
    counter = 0
    labels = dict()
    for i in levels:
      scale.insert(counter,str(pow(e,i)))
      labels[levels[counter]] = scale[counter]
      counter = counter + 1
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    CS = pylab.contour(xi,yi,zi, levels, colors='black', lynestiles='solid')
    print str(labels)
    pylab.clabel(CS,levels[0::2], fmt=labels, fontsize=20, inline=1)
    #p.clabel(CS,levels[1::2], fmt={log(1e-9):'1e-9', log(5e-9):'5e-9', log(1e-8):'1e-8', log(5e-8):'5e-8', log(1e-7):'1e-7', log(5e-7):'5e-7', log(1e-6):'1e-6', log(5e-6):'5e-8', log(1e-5):'1e-5', log(5e-5):'5e-5', log(1e-4):'1e-4', log(2e-4):'2e-4', log(4e-4):'4e-4', log(8e-4):'8e-4', log(1.2e-3):'1.2e-3', log(2.4e-3):'2.4e-3', log(4.8e-3):'4.8e-3'}, fontsize=9, inline=1)
    pylab.xlim((-7.0,7.0))
    pylab.ylim((0, 7.0))
    pylab.title('iPBS absolute error', fontsize=25)
    pylab.xlabel(r'$\displaystyle z/\lambda_D$',fontsize=20)
    pylab.ylabel(r'$\displaystyle r/\lambda_D$',fontsize=20)
    #pylab.show()
    savefig('ipbs_abs_far')
    return xi, yi, zi

def main():
    filename = sys.argv[1]
    x, y, ref, sol, dif, rel = readfile(filename)
    z = log(abs(dif))
    #z = abs(rel)
    xi, yi = defgrid(x/5.,y/5.)
    zi = interpolacion(x,y,z,xi,yi)
    plotear(xi,yi,zi)

if __name__ == "__main__":
    if len (sys.argv) == 2:
        main ()
    else:
        print "How to use: %s data_file" % sys.argv[0]
