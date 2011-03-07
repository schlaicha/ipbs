#!/usr/bin/python
# -*- encoding: utf8 -*-

import sys
from scipy import *
from matplotlib.delaunay import *
import pylab as p
from matplotlib import colors,ticker,pyplot

def readfile(file):
    x, y ,z = loadtxt(file, unpack=True) 
    return x, y, z

def defgrid(x,y):
    paso = 1000.
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
    interior1 = sqrt(((xi+1.5)**2) + (yi**2)) < 1.0 
    interior2 = sqrt(((xi-1.5)**2) + (yi**2)) < 1.0
    zi[interior1] = ma.masked
    zi[interior2] = ma.masked
    p.figure(figsize=(16,10))
    pyplot.jet()
    max=2.8
    min=0.4
    steps = 50
    levels=list()
    labels=list()
    for i in range(0,steps):
	levels.append(int((max-min)/steps*100*i)*0.01+min)
    for i in range(0,steps/2):
	labels.append(levels[2*i])
    CSF = p.contourf(xi,yi,zi,levels,norm=colors.LogNorm())
    CS = p.contour(xi,yi,zi,levels, format='%.3f', labelsize='18')
    p.clabel(CS,labels,inline=1,fontsize=9)
    p.title('electrostatic potential of two spherical colloids, R=lambda/3',fontsize=24)
    p.xlabel('z-coordinate (3*lambda)',fontsize=18)
    p.ylabel('radial coordinate r (3*lambda)',fontsize=18)
    # add a vertical bar with the color values
    cbar = p.colorbar(CSF,ticks=labels,format='%.3f')
    cbar.ax.set_ylabel('potential (reduced units)',fontsize=18)
    cbar.add_lines(CS)
    p.show()

def main():
    filename = sys.argv[1]
    x, y, z = readfile(filename)
    xi, yi = defgrid(x,y)
    zi = interpolacion(x,y,-z,xi,yi)
    plotear(xi,yi,zi)

if __name__ == "__main__":
    if len (sys.argv) == 2:
        main ()
    else:
        print "How to use: %s data_file1" % sys.argv[0]
