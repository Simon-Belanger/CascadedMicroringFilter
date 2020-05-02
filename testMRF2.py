"""
Script used to test a single microring filter that consists of multiple rings.

Author      : Simon Belanger-de Villers (simon.belanger-de-villers.1@ulaval.ca)
Created     : July 24th 2019
Last edited : July 24th 2019
"""

from components.couplers.DC_from_data import DC_from_data
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import flattop5
from components.MRF import MRF
import numpy as np

# Parameters for the rings
num_rings 	= 2 		# Number of cascaded rings [-]
radius 		= 2.5e-6 	# Radius of the ring resonators [m]
neff 		= 2.4449 	# Effective index of the ring waveguide []
ng 			= 4.18 		# Group index of the ring waveguide []
alpha_wg	= 3. 		# Losses of the ring waveguide [dB/cm]

# Parameters for the couplers (num_rings+1)
couplers = [virtual_DC(0.1,0.99),  
			virtual_DC(0.05,0.99), 
			virtual_DC(0.1,0.99)]


# Build the Microring filter and measure transmission across the wavelength range
mrf = MRF('', num_rings, radius, neff, ng, alpha_wg, couplers, [1, 0, 0])
mrf.sweep(np.linspace(1500,1600,1000)*1e-9)