#!/usr/bin/env python

from pygame import *
from math import *
from random import *
import sys
from time import sleep
import collections
# import pdb;
import os
from numpy import fromfile
from numpy import sin,pi
from os import system
import numpy as np

execfile('globals.py')                    # global simulation parameters.
####################################################################
# Display stuff
# screen init
# Set in globals.py:
# Nsites=128
# Cellsize = 4
# Width = Cellsize*Nsites
# Height = Cellsize*Nsites

ncount = 0                                # number of calls to trace
palette = []                              # color map array of RGB triples
cell = None                               # will be set to a Surface in gr_init
screen = None                             # will be set to a Surface in gr_init

# set window position
def gr_init():
    global palette, cell, screen
    os.environ['SDL_VIDEO_WINDOW_POS'] = "%d,%d" % (20,300)
    init()                                    # pygame
    display.set_caption("Heredity")
    screen = display.set_mode([Width, Height])
    draw.rect(screen, [10, 10, 10],(0, 0 , Width, Height + 1), 0)
    display.update()
    cell = Surface((Cellsize,Cellsize))

    # rainbow colormap:
    norm = lambda x: min(max(int((x+1)*128),0),255)
    s = lambda t: sin(2*pi*t)
    spec = lambda t: (norm(s(t*0.9+0.2)), norm(s(t*0.9+0.9)), norm(s(t*0.9+0.5)))
    palette = tuple(spec(x/256.) for x in range(256))

# display lattice as a line (scrolling)
def gr_disp(lattice):            # eachcolval = [key,activityvalue]
    global ncount, cell,screen,latmn,latmx
    # do the scroll:
    # convert to ints (color map entries)
    pvals = [int(255*((x-latmn)/(latmx-latmn))) for x in lattice]   # color idx into palette
    pvals = [min([x,255]) for x in pvals]
    pvals = [max([x,0]) for x in pvals]
    # compute y
    if ncount<Nsites:      # first, don't scroll
        y = Cellsize*ncount
    else:                                 #  then scroll when full
        y = Height-Cellsize           # bottom
        screen.scroll(0,-1*Cellsize)    # shift up
    for i in range(len(pvals)):         # go through lattice sites, changing x, same y, for blit
        col = pvals[i]
        cell.fill(palette[col])
        x = Cellsize*i+1
        screen.blit(cell,(x,y))
    ncount += 1
    display.update()

# End display stuff
####################################################################



def main():
    gr_init()
    lattice = [x/float(Nsites) for x in range(Nsites)]
    for i in range(500):
        gr_disp(lattice)
#        lattice = [lattice[(i+1)%Nsites] for i in range(len(lattice))]

if __name__=='__main__':
    main()

