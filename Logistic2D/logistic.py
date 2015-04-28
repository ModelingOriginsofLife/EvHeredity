#!/usr/bin/env python

import numpy as np
import sys
from graph import *

execfile('globals.py')

np.random.seed(10)

# bit locations for fixed pts and targets:
bits0 = range(0,Nbits) # 4 bits for fix0
bits1 = range(Nbits,2*Nbits) # 4 bits for fix1
bitsT = range(2*Nbits,3*Nbits)
rseq = [1*(i%3==0) for i in range(Nvar)]
rseq[2] = 1
rseq[7] = 1
#print rseq

def get0(bits):
    rtn=0
    for i in range(Nbits):
        rtn = rtn<<1 | (bits[bits0[i]]^rseq[bits0[i]]) # NG hack to mix things up a bit
    rtn = rtn / pow(2,Nbits)
    return rtn
def get1(bits):
    rtn=0
    for i in range(Nbits):
        rtn = rtn<<1 | (bits[bits1[i]]^rseq[bits1[i]])
    rtn = rtn / pow(2,Nbits)
    return rtn
def getT(bits):
    rtn=0
    targbits = [Target>>i & 1 for i in range(Nbits)]
    trybits = [bits[i] for i in bitsT]
    rtn = sum([x==y for x,y in zip(targbits,trybits)]) /float(len(targbits)) # 1-Hamming dist
    return rtn

class Site:
    logis = []                         # state
    fix = []                           # 2 fixed pts
    W = 0.0                            # weight of my fixed pt to nbrs
    bits = [] # will be 0,1 according to whether state is nearest fix[0] or fix[1]
    targ = []
	
    def __init__(self):
        self.logis = np.random.rand(Nvar)
        self.fix = [0.25,0.75]
        self.bits = [0]*Nvar 

    def Bits(self):
        for i in range(Nvar):
            x = self.logis[i]
            d0 = abs(x-self.fix[0])
            d1 = abs(x-self.fix[1])
            if(d0>d1):
                self.bits[i] = 1
            else:
                self.bits[i] = 0
                
    def Fix(self):
        self.Bits()
        self.fix[0] = get0(self.bits)
        self.fix[1] = get1(self.bits)
        
    def Setw(self,targ):
#        assert targ<16                  #4 bits
		
        targbits = [targ>>i & 1 for i in range(Nbits)]
        trybits = [self.bits[i] for i in bitsT]
        self.W = max(0,Wbias*(sum([x==y for x,y in zip(targbits,trybits)]) /float(len(targbits)) - 0.5))# 1-Hamming dist

    def Iterate(self):
        for i in range(Nvar):
            for _ in range(Niter):
                self.logis[i] = R*self.logis[i]*(1-self.logis[i])

    def Contract(self):
        for i in range(Nvar):
            for _ in range(Ncontract):
                #self.logis[i] = self.logis[i] - Crate*(self.logis[i] - self.fix[self.bits[i]])
                self.logis[i] = self.logis[i] - Crate*(self.logis[i] - self.targ[i])
            if self.logis[i] < 0:
                self.logis[i] = 0.0001
            if self.logis[i] >1 :
                self.logis[i] = 0.9999

        
# nearest nbr
# Nbrs =[[i-1,i+1] for i in range(Nsites)]
# Nbrs[0][0] = Nsites-1                     # circular boundary...
# Nbrs[Nsites-1][1] = 0

# next to nearest nbr

def mk2dnbrs():
    global Nbrs
    Nbrs =[[i-Nside,i+1,i+Nside,i-1] for i in range(Nsites)]
    # top row upper nbr:
    for i in range(Nside):
        Nbrs[i][0] = Nsites-Nside+i
    # bottom row lower nbr:
    for i in range(Nside):
        ii = Nsites-Nside+i
        Nbrs[ii][2] = i
    # L col
    for i in range(Nside):
        ii = Nside*i
        Nbrs[ii][3] = ii+Nside-1
    # R col
    for i in range(Nside):
        ii = (i+1)*Nside - 1
        Nbrs[ii][1] = ii-Nside+1    
    # upper L corner L nbr
    Nbrs[0][3] = Nside-1
    # lower R corner R nbr
    Nbrs[Nsites-1][1] = Nsites-Nside
    
mk2dnbrs()    

# For testing the 2d neighbors by taking a close look...
#~ for i in range(len(Nbrs)):
     #~ print i,i/Nside,i%Nside,Nbrs[i]


class Lattice:
    sites = []

    def __init__(self):
        self.sites = [Site() for _ in range(Nsites)]
    # establish fixed points from bit vector
    def Fix(self):
        for i in range(Nsites):
            self.sites[i].Fix()

    # update fixed points of each site with weighted sum of nbr fixed pts
    def Fixw(self):
        # first compute into temporary storate
        ttemp = [ [0 for j in range(Nvar)] for i in range(Nsites)]
        wtmp = [ [0,0] for i in range(Nsites)] # self.sites[i].fix
        for i in range(Nsites):
            for j in range(2):            # j over fix
                Wnorm = 1                 # for current site
                for nbr in Nbrs[i]:
                    Wnorm += self.sites[nbr].W
                fixcur = self.sites[i].fix[j] / Wnorm
                tcur = [self.sites[i].fix[k]/Wnorm for k in self.sites[i].bits]                
                for nbr in Nbrs[i]:
                    fixcur += self.sites[nbr].W * self.sites[nbr].fix[j] / Wnorm
                    for idx in range(len(tcur)):
						tcur[idx] += self.sites[nbr].fix[ self.sites[nbr].bits[idx] ] * self.sites[nbr].W / Wnorm
                    assert fixcur>=0
                    assert fixcur<=1
                wtmp[i][j] = fixcur
                ttemp[i] = tcur
        # now copy into the real thing
        for i in range(Nsites):
			self.sites[i].targ = ttemp[i]			
            #~ for j in range(2):            # j over fix
                #~ self.sites[i].fix[j] = wtmp[i][j]

    def Update(self):
        for site in self.sites:           # freewheel
            site.Iterate()                
        self.Fix()                        # establish fixed pts
        for site in self.sites:           # compute weights according to Target match
            site.Setw(Target)                   
        self.Fixw()                       # new fixd pts from coupling
        for site in self.sites:           # contract toward new fixed pts
            site.Contract()

    def Colors(self):
        # fitness
        if Coltype==0:
            return [getT(self.sites[i].bits) for i in range(Nsites)]  
        # ave of fixed pts
        if Coltype==1:
            return [(self.sites[i].fix[0]+self.sites[i].fix[1])/2 for i in range(Nsites)]
#        return [self.sites[i].logis[0] for i in  range(Nsites)]
            
def main():
    gr_init()
    foo = Lattice()
    for i in range(500):
        gr_disp(foo.Colors())
        foo.Update()
        print foo.sites[0].fix

if __name__=='__main__':
    main()
