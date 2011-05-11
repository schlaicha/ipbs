#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
from scipy import *
from matplotlib.delaunay import *
import pylab as p

def readfile(file):
    x, y ,ref, sol, dif, rel = loadtxt(file, unpack=True) 
    return x, y, ref, sol, dif, rel

def defgrid(x,y):
    paso = 1600.
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
    p.figure(figsize=(16,10))
    #levels = [zi.min() , zi.max() , (zi.max()-zi.min())/10]
    #levels = [1E-8,1E-7,1E-6,1E-5,1E-4,1E-3,5E-2,1E-2]
    #CSF = p.contourf(xi,yi,zi,levels)
    CSF = p.contourf(xi,yi,zi)
    CS = p.contour(xi,yi,zi)
    p.clabel(CS)
    p.title('iPBS absolute difference')
    p.xlabel('x-coordinate',fontsize=12)
    p.ylabel('y-coordinate',fontsize=12)
    # add a vertical bar with the color values
    cbar = p.colorbar(CSF)
    cbar.ax.set_ylabel('Solution difference',fontsize=12)
    cbar.add_lines(CS)
    p.show()

def main():
    filename = sys.argv[1]
    x, y, ref, sol, dif, rel = readfile(filename)
    #z = abs((ref - sol))
    z = dif
    xi, yi = defgrid(x,y)
    zi = interpolacion(x,y,z,xi,yi)
    plotear(xi,yi,zi)

if __name__ == "__main__":
    if len (sys.argv) == 2:
        main ()
    else:
        print "How to use: %s data_file" % sys.argv[0]
