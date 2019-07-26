import numpy as np 
import matplotlib.pyplot as plt
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB
import math

# Configure one microring
from components.couplers.DC_from_data import DC_from_data
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import flattop4
from components.MRF import MRF

num_rings 	= 4 		# Number of cascaded rings [-]
radius 		= 2.5e-6 	# Radius of the ring resonators [m]
neff 		= 2.4449 	# Effective index of the ring waveguide []
ng 			= 4.18 		# Group index of the ring waveguide []
alpha_wg	= 3. 		# Losses of the ring waveguide [dB/cm]
k_vec 		= flattop4(0.1)
loss_coupler = [0.99,0.99,0.99,0.99,0.99]
couplers = [virtual_DC(k_vec[0], loss_coupler[0]),  
			virtual_DC(k_vec[1], loss_coupler[1]), 
			virtual_DC(k_vec[2], loss_coupler[2]), 
			virtual_DC(k_vec[3], loss_coupler[3]), 
			virtual_DC(k_vec[4], loss_coupler[4])]


# Build the Microring filter and measure transmission across the wavelength range
mrf = MRF('', num_rings, radius, neff, ng, alpha_wg, couplers, [1, 0, 0])

# Input Signal
wavelength  = np.linspace(1540e-9, 1560e-9, 1000)
inputPower  = np.zeros(1000)

print(np.abs(mrf.TMM(1555e-9, E_in=0, E_add=0, E_drop=1, E_thru=0)))

# 5 RINGS
# (I A D T)(1 0 0 0) (X T X X)(0 ~1 ~0 0) 
# (I A D T)(0 1 0 0) (X X D X)(0 0 ~1 0)
# (I A D T)(0 0 1 0) (A X X X)(~1 0 0 ~0)
# (I A D T)(0 0 0 1) (X X X I)(0 0 0 ~1) 
# (IADT) (IDTA) (ATDI)

# 4 RINGS
# (I A D T)(1 0 0 0) (X T X X)(0 ~1 0 ~0) 
# (I A D T)(0 1 0 0) (X X D X)(0 0 ~1 0)
# (I A D T)(0 0 1 0) (A X X X)(~1 0 ~0 0)
# (I A D T)(0 0 0 1) (X X X I)(0 0 0 ~1) 
# (IADT) (IDTA) (ATDI)