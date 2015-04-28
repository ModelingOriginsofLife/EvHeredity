

##############################################################
# Dynamics
Target = 21               # random # between 0 and 2^Nbits-1
Nbits = 12
Niter = 6                # number of iterations for freewheel updates
Ncontract = 1             # number of iterations for contraction updates
R = 3.7                   # R value for logistic update
Crate = 1               # contraction rate
Wbias = 1.0                 # bias toward fitness of nbr.

##############################################################
# Lattice
Nside = 32                     # side of 2d lattice
Nsites=Nside*Nside              # total lattice size
Nvar = 3*Nbits                       # number of internal dynamic variables


##############################################################
# Graphics
Cellsize = 4
Width = Cellsize*Nside
Height = Cellsize*Nside
latmn = 0                                 # lattice values latmin and less will be color 0
latmx = 1.0                               # lattice values latmax and more will be color 255
Coltype = 0                               # 1 => mean of fixed pts, 0 => fitness

