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
    interior = sqrt((xi**2) + (yi**2)) < 10.0 
    zi[interior] = ma.zeros
    p.figure(figsize=(16,10))
    #levels = [zi.min() , zi.max() , (zi.max()-zi.min())/10]
    #levels = [0 ,0.005, 0.01,0.015 , 0.02, 0.03, 0.04, 0.05]
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
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    x1, y1, z1 = readfile(filename1)
    x2, y2, z2 = readfile(filename2)
    z = abs(z2 - z1)
    xi, yi = defgrid(x1,y1)
    zi = interpolacion(x1,y1,z,xi,yi)
    plotear(xi,yi,zi)

if __name__ == "__main__":
    if len (sys.argv) == 3:
        main ()
    else:
        print "How to use: %s data_file1 %data_file2 " % sys.argv[0]
