

##############################################################
# Dynamics
Target = 21               # random # between 0 and 2^Nbits-1
Nbits = 12
Niter = 6                # number of iterations for freewheel updates
Ncontract = 1             # number of iterations for contraction updates
alpha = 1.3               # Henon
beta = 0.3                # params
Crate = 0.78              # contraction rate
Wbias = 2.0                 # bias toward fitness of nbr.

##############################################################
# Lattice
Nsites=128                      # total lattice size
Nvar = 5*Nbits                 # number of internal dynamic variables (must be at least 5*Nbits for 2 fixed pts+targ


##############################################################
# Graphics
Cellsize = 4
Width = Cellsize*Nsites
Height = Cellsize*Nsites
latmn = 0                                 # lattice values latmin and less will be color 0
latmx = 1.0                               # lattice values latmax and more will be color 255
Coltype = 0                               # 1 => mean of fixed pts, 0 => fitness

