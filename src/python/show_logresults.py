#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
from scipy import *
from matplotlib.delaunay import *
from matplotlib import colors, ticker
import pylab as p
import numpy as np

def readfile(file):
    x, y ,ref, sol, dif, rel = loadtxt(file, unpack=True) 
    return x, y, ref, sol, dif, rel

def defgrid(x,y):
    paso = 600.
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
    zi[interior] = ma.ones
    p.figure(figsize=(16,10))
    #levels = [zi.min() , zi.max() , (zi.max()-zi.min())/10]
    levels = [.75e-5, 1.25e-4, 2.5e-4, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64]
    #v = np.logspace(1e-8,1e-4, num=50, endpoint=True, base=2.0)
    #levels = [1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,1E-2,1E-1,1E0,1E1]
    CSF = p.contourf(xi,yi,zi,levels,norm=colors.LogNorm())
    #CSF = p.contourf(xi,yi,zi,levels)
    CS = p.contour(xi,yi,zi,levels)
    p.clabel(CS)
    p.title('IPBS vs Neumann B.C. (absolute difference)')
    p.ylabel('radial coordinate r',fontsize=12)
    p.xlabel('z-coordinate',fontsize=12)
    # add a vertical bar with the color values
    cbar = p.colorbar(CSF)
    cbar.ax.set_ylabel('electrostatic potential (reduced units)',fontsize=12)
    #cbar.add_lines(CS)
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
