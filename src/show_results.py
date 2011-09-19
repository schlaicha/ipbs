#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
#import numpy as np
from scipy import *
from matplotlib.delaunay import *
import pylab as p
import matplotlib.cm as cm
from matplotlib import colors

def readfile(file):
    x, y ,ref, sol, dif, rel = loadtxt(file, unpack=True) 
    return x, y, ref, sol, dif, rel

def defgrid(x,y):
    paso = 1200.
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
    zi[interior] = ma.zeros
    #exterior = sqrt((xi**2) + (yi**2)) > 25.0
    #zi[exterior] = ma.empty
    p.figure(figsize=(16,10))
    levels = [log(.75e-5), log(1.25e-4), log(2.5e-4), log(0.005), log(0.01), log(0.02), log(0.04), log(0.08), log(0.16), log(0.32), log(0.64)]
    scale = list()
    counter = 0
    for i in levels:
      scale.insert(counter, pow(e,i))
      counter = counter + 1
      #scale[i]=10^i
    #  scale[i] = pow(10,levels[i])
    #levels = [0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56]
    #levels = [log(0.005), log(0.01), log(0.02),log(0.04), log(0.08), log(0.16), log(0.32), log(0.64), log(1.28), log(2.56)]
    #levels = [zi.min() , zi.max() , (zi.max()-zi.min())/10]
    #v = np.linspace(0., 1., 25, endpoint=True)
    #v = np.logspace(1e-8,1e-1, num=10, endpoint=True, base=2)
    #levels = (0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.)
    CSF = p.contourf(xi,yi,zi, levels, cmap=cm.jet)
    #CSF = p.contourf(xi,yi,zi)
    CS = p.contour(xi,yi,zi, levels)
    p.clabel(CS,levels, fmt={log(.75e-5):'hallo', log(1.25e-4):'3', log(2.5e-4):'4', log(0.005):'5', log(0.01):'6', log(0.02):'7', log(0.04):'8', log(0.08):'9', log(0.16):'10', log(0.32):'11', log(0.64):'12'})
    p.title('iPBS relative error')
    p.xlabel('x-coordinate',fontsize=12)
    p.ylabel('y-coordinate',fontsize=12)
    # add a vertical bar with the color values
    cbar = p.colorbar(CSF, ticks=levels, format='%g')
    cbar.ax.set_yticklabels(scale)
    cbar.ax.set_ylabel('Relative error in %',fontsize=12)
    cbar.add_lines(CS)
    p.show()

def main():
    filename = sys.argv[1]
    x, y, ref, sol, dif, rel = readfile(filename)
    z = log(abs(dif))
    #z = abs(rel)
    xi, yi = defgrid(x,y)
    zi = interpolacion(x,y,z,xi,yi)
    plotear(xi,yi,zi)

if __name__ == "__main__":
    if len (sys.argv) == 2:
        main ()
    else:
        print "How to use: %s data_file" % sys.argv[0]
