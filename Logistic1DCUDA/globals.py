

##############################################################
# Dynamics
Target = 21               # random # between 0 and 2^Nbits-1
Nbits = 12
Niter = 6                # number of iterations for freewheel updates
Ncontract = 1             # number of iterations for contraction updates
R = 3.7                   # R value for logistic update
Crate = 0.99              # contraction rate
Wbias = 1.0                 # bias toward fitness of nbr.

##############################################################
# Lattice
Nsites=64                      # total lattice size
Nvar = 3*Nbits                 # number of internal dynamic variables


##############################################################
# Graphics
Cellsize = 8
Width = Cellsize*Nsites
Height = Cellsize*Nsites
latmn = 0                                 # lattice values latmin and less will be color 0
latmx = 1.0                               # lattice values latmax and more will be color 255
Coltype = 0                               # 1 => mean of fixed pts, 0 => fitness

