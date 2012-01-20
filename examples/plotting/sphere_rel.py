#!/usr/bin/python
# -*- encoding: utf8 -*-

# show the absolute difference between the iPBS solution of a single sphere in reduced symmetry
# and the constant surface flux (Neumann) reference solution

import sys
from scipy import *
from matplotlib.delaunay import *
import matplotlib as mp
import pylab as p
import matplotlib.cm as cm
from matplotlib import colors
import numpy as np

def readfile(file):
    x, y ,ref, sol, dif, rel = loadtxt(file, unpack=True) 
    print x[:10]
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
    #exterior = sqrt((xi**2) + (yi**2)) > 15.0
    #zi[exterior] = ma.empty
    p.figure(figsize=(16,10))
    levels = [log(1.1e-3), log(1.12e-3), log(1.14e-3), log(1.18e-3), log(1.26e-3), log(1.42e-3), log(1.74e-3), log(1.29e-3)]
    #levels = [log(3.25e-3),log(3.5e-3), log(4.0e-3), log(4.1e-3),log(4.2e-3),log(4.4e-3), log(4.6e-3)]
    #levels = [8e-4, 1e-3, 1.15e-3, 1.20e-3, 1.25e-3, 1.3e-3, 1.4e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3, 3.5e-3, 4e-3, 4.5e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 1e-2] 
    #levels = [7e-4, 8e-4, 9e-4, 1e-3, 1.1e-3, 1.12e-3, 1.125e-3, 1.13e-3 ,1.14e-3, 1.16e-3, 1.2e-3, 1.22e-3, 1.24e-3, 1.26e-3, 1.28e-3, 1.3e-3, 1.4e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3, 3.5e-3, 4e-3, 4.5e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3, 1e-2] 
    #levels = [log(2.8125e-7),log(5.625e-7), log(1.125e-6) , log(3.25e-6), log(.75e-5), log(1.25e-4), log(2.5e-4), log(0.005), log(0.01), log(0.02), log(0.04), log(0.08), log(0.16), log(0.32), log(0.64)]
    #levels = [log(.64), log(.32), log(.16), log(.08), log(.04), log(.02), log(.01), log(0.005), log(0.0025)]
    #levels = np.linspace(7.e-4, 1e-2, num=15, endpoint=True)
    #levels = np.logspace(-8,-1, num=20, endpoint=True, base=10)
    scale = list()
    counter = 0
    labels = dict()
    for i in levels:
      scale.insert(counter,str('%1.2e' %pow(e,i)))
      #scale.insert(counter,str(i))
      labels[levels[counter]] = scale[counter]
      counter = counter + 1
    #CSF = p.contourf(xi,yi,zi, levels, cmap=cm.jet)
    mp.rcParams['contour.negative_linestyle'] = 'solid'
    CS = p.contour(xi,yi,zi, levels, colors='black', lynestiles='solid')
    print str(labels)
    p.clabel(CS,levels, fmt=labels, fontsize=12, inline=1, color='black')
    #p.clabel(CS)
    p.xlim((-6,6))
    p.ylim((0, 6))
    #p.clabel(CS,levels[1::2], fmt='%1.1e', fontsize=14, inline=1, color='black')
    #p.clabel(CS,levels[1::2], fmt={log(1e-9):'1e-9', log(5e-9):'5e-9', log(1e-8):'1e-8', log(5e-8):'5e-8', log(1e-7):'1e-7', log(5e-7):'5e-7', log(1e-6):'1e-6', log(5e-6):'5e-8', log(1e-5):'1e-5', log(5e-5):'5e-5', log(1e-4):'1e-4', log(2e-4):'2e-4', log(4e-4):'4e-4', log(8e-4):'8e-4', log(1.2e-3):'1.2e-3', log(2.4e-3):'2.4e-3', log(4.8e-3):'4.8e-3'}, fontsize=9, inline=1)
    p.title('iPBS relative error')
    p.xlabel('z-coordinate',fontsize=14)
    p.ylabel('radial coordinate',fontsize=14)
    # add a vertical bar with the color values
    #cbar = p.colorbar(CSF, ticks=levels, format='%e')
    #cbar.ax.set_yticklabels(scale)
    #cbar.ax.set_ylabel('iPBS relative error',fontsize=12)
    #cbar.add_lines(CS)
    return xi, yi, zi

def main():
    filename = sys.argv[1]
    x, y, ref, sol, dif, rel = readfile(filename)
    #z = log(abs(rel))
    z = abs(rel)
    xi, yi = defgrid(x,y)
    zi = interpolacion(x,y,z,xi,yi)
    return plotear(xi,yi,zi)


if __name__ == "__main__":
    if len (sys.argv) == 2:
        xi, yi, zi=main ()
    else:
        print "How to use: %s data_file" % sys.argv[0]
