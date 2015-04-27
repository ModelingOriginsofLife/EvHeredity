#!/opt/local/bin/python


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


argv = sys.argv
if len(argv)>1:
    datfile = argv[1]
    datout = open(datfile,"w")
else:
    datfile = None
    
# set window position
os.environ['SDL_VIDEO_WINDOW_POS'] = "%d,%d" % (20,300)



Nsites=128
Cellsize = 4
Width = Cellsize*Nsites
Height = Cellsize*Nsites
ncount = 0                                #number of calls to trace

lattice = np.zeros(Nsites)

# screen init
init()
screen = display.set_mode([Width, Height])
display.set_caption("Activity")
draw.rect(screen, [10, 10, 10],(0, 0 , Width, Height + 1), 0)
display.update()
cell = Surface((Cellsize,Cellsize))

# rainbow colormap:
norm = lambda x: min(max(int((x+1)*128),0),255)
s = lambda t: sin(2*pi*t)
spec = lambda t: (norm(s(t*0.9+0.2)), norm(s(t*0.9+0.9)), norm(s(t*0.9+0.5)))
palette = tuple(spec(x/256.) for x in range(256))

for i in range(20):
    draw.rect(screen, [10, 10, 10],(0, 0 , Width, Height + 1), 0)
    display.update()
    draw.rect(screen, [100, 100, 100],(0, 0 , Width, Height + 1), 0)
    display.update()
    
