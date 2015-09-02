#!/usr/bin/env python

import numpy as np
import sys
from graph import *
from numba import jit

execfile('globals.py')

#np.random.seed(10)

# bit locations for fixed pts and target:
bits0x = range(0,Nbits) # Nbits bits for fix0x
bits0y = [Nbits + x for x in range(0,Nbits)] # Nbits bits for fix0y
bits1x = [2*Nbits + x for x in range(0,Nbits)] # Nbits bits for fix1x
bits1y = [3*Nbits + x for x in range(0,Nbits)]# Nbits bits for fix1y
bitsT = [4*Nbits + x for x in range(0,Nbits)] # Nbits bits for target
rseq = [1*(i%3==0) for i in range(Nvar)]   # for NG hack to mix things up a bit
rseq[2] = 1
rseq[7] = 1
#print rseq

@jit
def getbits(bits,rng):
    rtn=0
    for i in rng:
#        rtn = rtn <<1 | bits[i]
        rtn = rtn<<1 | (bits[i]^rseq[i]) # NG hack to mix things up a bit
    rtn = 2*((rtn / pow(2,Nbits)) - 0.5)
    return rtn

@jit
def getT(bits):
    rtn=0
    for i in bitsT:
        rtn = rtn<<1 | (bits[i]^rseq[i]) # NG hack to mix things up a bit
    rtn = rtn / pow(2,Nbits)
    return rtn
    

class Site:
    statex = []                         # state
    statey = []                         # state
    fix = []                           # 2 fixed pts
    W = 0.0                            # weight of my fixed pt to nbrs
    bits = [] # will be 0,1 according to whether state is nearest fix[0] or fix[1]
    targx = []
    targy = []
    def __init__(self):
        self.statex = [x-0.5 for x in np.random.rand(Nvar)]
        self.statey = [x-0.5 for x in np.random.rand(Nvar)]
        self.fix = [[0.2,0.2],[0.4,0.4]]
        self.bits = [0]*Nvar 

    def Bits(self):                       # set self.bits according to which fixed pt is nearest
        for i in range(Nvar):
            x = self.statex[i]
            y = self.statey[i]
            xfix = self.fix[0][0]
            yfix = self.fix[0][1]
            d0 = sqrt((x-xfix)*(x-xfix))
            xfix = self.fix[1][0]
            yfix = self.fix[1][1]
            d1 = sqrt((x-xfix)*(x-xfix))
            if(d0>d1):
                self.bits[i] = 1
            else:
                self.bits[i] = 0
                
    def Fix(self):
        self.Bits()
        self.fix[0][0] = getbits(self.bits,bits0x)
        self.fix[0][1] = getbits(self.bits,bits0y)
        self.fix[1][0] = getbits(self.bits,bits1x)
        self.fix[1][1] = getbits(self.bits,bits1y)
        
    def Setw(self,targ):
#        assert targ<16                  #4 bits
		
        targbits = [targ>>i & 1 for i in range(Nbits)]
        trybits = [self.bits[i] for i in bitsT]
        self.W = max(0,Wbias*(sum([x==y for x,y in zip(targbits,trybits)]) /float(len(targbits)) - 0.5))# 1-Hamming dist

    def Iterate(self):
        for i in range(Nvar):
            for _ in range(Niter):
                x = self.statex[i]
                y = self.statey[i]
                self.statex[i] = 1 - alpha * x * x + y
                self.statey[i] = beta * x

    def Contract(self):
        for i in range(Nvar):
            for _ in range(Ncontract):
                #self.state[i] = self.state[i] - Crate*(self.state[i] - self.fix[self.bits[i]])
                self.statex[i] = self.statex[i] - Crate*(self.statex[i] - self.targx[i])
                self.statey[i] = self.statey[i] - Crate*(self.statey[i] - self.targy[i])
            if (abs(self.statex[i]) > 1.5) | (abs(self.statey[i]) > 1.5): # throw it back if escaping
                self.statex[i] = 0.63     #near the fixed pt.
                self.statey[i] = 0.189

        
# nearest nbr
# Nbrs =[[i-1,i+1] for i in range(Nsites)]
# Nbrs[0][0] = Nsites-1                     # circular boundary...
# Nbrs[Nsites-1][1] = 0

# nearest nbr
Nbrs =[[i-1,i+1] for i in range(Nsites)]
Nbrs[0][0] = Nsites-1                     # circular boundary...
Nbrs[Nsites-1][1] = 0

# next to nearest nbr
# Nbrs =[[i-2,i-1,i+1,i+2] for i in range(Nsites)]
# Nbrs[0][0] = Nsites-2                     # circular boundary...
# Nbrs[0][1] = Nsites-1                     # circular boundary...
# Nbrs[1][0] = Nsites-1
# Nbrs[Nsites-2][3] = 0
# Nbrs[Nsites-1][2] = 0
# Nbrs[Nsites-1][3] = 1

###################
# for 2d:
# def mk2dnbrs():
#     global Nbrs
#     Nbrs =[[i-Nside,i+1,i+Nside,i-1] for i in range(Nsites)]
#     # top row upper nbr:
#     for i in range(Nside):
#         Nbrs[i][0] = Nsites-Nside+i
#     # bottom row lower nbr:
#     for i in range(Nside):
#         ii = Nsites-Nside+i
#         Nbrs[ii][2] = i
#     # L col
#     for i in range(Nside):
#         ii = Nside*i
#         Nbrs[ii][3] = ii+Nside-1
#     # R col
#     for i in range(Nside):
#         ii = (i+1)*Nside - 1
#         Nbrs[ii][1] = ii-Nside+1    
#     # upper L corner L nbr
#     Nbrs[0][3] = Nside-1
#     # lower R corner R nbr
#     Nbrs[Nsites-1][1] = Nsites-Nside
    
# mk2dnbrs()    

# # For testing the 2d neighbors by taking a close look...
# #~ for i in range(len(Nbrs)):
#      #~ print i,i/Nside,i%Nside,Nbrs[i]


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
        ##################################################
        # for x coord of both fixed pts:
        # first compute into temporary storage
        ttemp = [ [0 for j in range(Nvar)] for i in range(Nsites)]
        wtmp = [ [0,0] for i in range(Nsites)] # self.sites[i].fix
        for i in range(Nsites):
            for j in range(2):            # j over fix
                Wnorm = 1                 # for current site
                for nbr in Nbrs[i]:
                    Wnorm += self.sites[nbr].W
                fixcur = self.sites[i].fix[j][0] / Wnorm   #[0] => x coord
                tcur = [self.sites[i].fix[k][0]/Wnorm for k in self.sites[i].bits]                
                for nbr in Nbrs[i]:
                    fixcur += self.sites[nbr].W * self.sites[nbr].fix[j][0] / Wnorm
                    for idx in range(len(tcur)):
						tcur[idx] += self.sites[nbr].fix[ self.sites[nbr].bits[idx] ][0] * self.sites[nbr].W / Wnorm
                    assert fixcur>=-1.5
                    assert fixcur<=1.5
                wtmp[i][j] = fixcur
                ttemp[i] = tcur
        # now copy into the real thing
        for i in range(Nsites):
			self.sites[i].targx = ttemp[i]			
        ##################################################
        # for y coord of both fixed pts:
        # first compute into temporary storage
        ttemp = [ [0 for j in range(Nvar)] for i in range(Nsites)]
        wtmp = [ [0,0] for i in range(Nsites)] # self.sites[i].fix
        for i in range(Nsites):
            for j in range(2):            # j over fix
                Wnorm = 1                 # for current site
                for nbr in Nbrs[i]:
                    Wnorm += self.sites[nbr].W
                fixcur = self.sites[i].fix[j][1] / Wnorm   #[1] => y coord
                tcur = [self.sites[i].fix[k][1]/Wnorm for k in self.sites[i].bits]                
                for nbr in Nbrs[i]:
                    fixcur += self.sites[nbr].W * self.sites[nbr].fix[j][1] / Wnorm
                    for idx in range(len(tcur)):
						tcur[idx] += self.sites[nbr].fix[ self.sites[nbr].bits[idx] ][1] * self.sites[nbr].W / Wnorm
                    assert fixcur>=-1.5
                    assert fixcur<=1.5
                wtmp[i][j] = fixcur
                ttemp[i] = tcur
        # now copy into the real thing
        for i in range(Nsites):
			self.sites[i].targy = ttemp[i]			

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
#        return [self.sites[i].state[0] for i in  range(Nsites)]
            
def main():
    gr_init()
    foo = Lattice()
    for i in range(500):
        gr_disp(foo.Colors())
        foo.Update()

if __name__=='__main__':
    main()
        
