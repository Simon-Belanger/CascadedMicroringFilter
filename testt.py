"""
Test used to make sure the ports attribution in the TMM method is working right.

Author      : Simon Bélanger-de Villers (simon.belanger-de-villers.1@ulaval.ca)
Created     : July 24th 2019
Last edited : July 24th 2019
"""

import numpy as np 
import matplotlib.pyplot as plt
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB
import math

# Configure one microring
from components.couplers.virtual_DC import virtual_DC
from components.MRF import MRF

ei = 1
ed = 0
et = 0
ea = 0

# One ring
num_rings 	= 1			# Number of cascaded rings [-]
radius 		= 2.5e-6 	# Radius of the ring resonators [m]
neff 		= 2.4449 	# Effective index of the ring waveguide []
ng 			= 4.18 		# Group index of the ring waveguide []
alpha_wg	= 3. 		# Losses of the ring waveguide [dB/cm]
k_vec 		= [0.1, 0.1]
loss_coupler = [0.99,0.99]
couplers = [virtual_DC(k_vec[0], loss_coupler[0]),  
			virtual_DC(k_vec[1], loss_coupler[1])]


# Build the Microring filter and measure transmission across the wavelength range
mrf = MRF('', num_rings, radius, neff, ng, alpha_wg, couplers, [1])

# Input Signal
wavelength  = np.linspace(1540e-9, 1560e-9, 1000)
inputPower  = np.zeros(1000)

print(mrf.TMM(1555e-9, E_in=ei, E_add=ea, E_drop=ed, E_thru=et))


# 2 Rings
num_rings 	= 2			# Number of cascaded rings [-]
radius 		= 2.5e-6 	# Radius of the ring resonators [m]
neff 		= 2.4449 	# Effective index of the ring waveguide []
ng 			= 4.18 		# Group index of the ring waveguide []
alpha_wg	= 3. 		# Losses of the ring waveguide [dB/cm]
k_vec 		= [0.1, 0.05, 0.1]
loss_coupler = [0.99,0.99, 0.99]
couplers = [virtual_DC(k_vec[0], loss_coupler[0]),  
			virtual_DC(k_vec[1], loss_coupler[1]),
			virtual_DC(k_vec[2], loss_coupler[2])]


# Build the Microring filter and measure transmission across the wavelength range
mrf = MRF('', num_rings, radius, neff, ng, alpha_wg, couplers, [1])

# Input Signal
wavelength  = np.linspace(1540e-9, 1560e-9, 1000)
inputPower  = np.zeros(1000)

#print(mrf.TMM(1555e-9, E_in=ei, E_add=ea, E_drop=ed, E_thru=et))

# 3 Rings
num_rings 	= 3			# Number of cascaded rings [-]
radius 		= 2.5e-6 	# Radius of the ring resonators [m]
neff 		= 2.4449 	# Effective index of the ring waveguide []
ng 			= 4.18 		# Group index of the ring waveguide []
alpha_wg	= 3. 		# Losses of the ring waveguide [dB/cm]
k_vec 		= [0.1, 0.05, 0.05, 0.1]
loss_coupler = [0.99,0.99, 0.99, 0.99]
couplers = [virtual_DC(k_vec[0], loss_coupler[0]),  
			virtual_DC(k_vec[1], loss_coupler[1]),
			virtual_DC(k_vec[2], loss_coupler[2]),
			virtual_DC(k_vec[3], loss_coupler[3])]


# Build the Microring filter and measure transmission across the wavelength range
mrf = MRF('', num_rings, radius, neff, ng, alpha_wg, couplers, [1])

# Input Signal
wavelength  = np.linspace(1540e-9, 1560e-9, 1000)
inputPower  = np.zeros(1000)

print(mrf.TMM(1555e-9, E_in=ei, E_add=ea, E_drop=ed, E_thru=et))

# 4 Rings
num_rings 	= 4			# Number of cascaded rings [-]
radius 		= 2.5e-6 	# Radius of the ring resonators [m]
neff 		= 2.4449 	# Effective index of the ring waveguide []
ng 			= 4.18 		# Group index of the ring waveguide []
alpha_wg	= 3. 		# Losses of the ring waveguide [dB/cm]
k_vec 		= [0.1, 0.05, 0.001, 0.05, 0.1]
loss_coupler = [0.99,0.99, 0.99, 0.99, 0.99]
couplers = [virtual_DC(k_vec[0], loss_coupler[0]),  
			virtual_DC(k_vec[1], loss_coupler[1]),
			virtual_DC(k_vec[2], loss_coupler[2]),
			virtual_DC(k_vec[3], loss_coupler[3]),
			virtual_DC(k_vec[4], loss_coupler[4])]


# Build the Microring filter and measure transmission across the wavelength range
mrf = MRF('', num_rings, radius, neff, ng, alpha_wg, couplers, [1])

# Input Signal
wavelength  = np.linspace(1540e-9, 1560e-9, 1000)
inputPower  = np.zeros(1000)

#print(mrf.TMM(1555e-9, E_in=ei, E_add=ea, E_drop=ed, E_thru=et))
