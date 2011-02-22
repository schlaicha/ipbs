#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
from scipy import *
from matplotlib.delaunay import *
import pylab as p

def readfile(file):
    x, y ,z = loadtxt(file, unpack=True) 
    return x, y, z

def defgrid(x,y):
    paso = 800.
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
    #interior = sqrt((xi**2) + (yi**2)) < 1.0 - 1E-4
    #zi[interior] = ma.masked
    p.figure(figsize=(16,10))
    #levels = [zi.min() , zi.max() , (zi.max()-zi.min())/10]
    #levels = [0,1E-6, 1E-5,1E-4,1E-3]
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
    x, y, z = readfile(filename)
    xi, yi = defgrid(x,y)
    zi = interpolacion(x,y,z,xi,yi)
    plotear(xi,yi,zi)

if __name__ == "__main__":
    if len (sys.argv) == 2:
        main ()
    else:
        print "How to use: %s data_file1" % sys.argv[0]
