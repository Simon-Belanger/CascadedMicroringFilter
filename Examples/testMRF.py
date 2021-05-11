"""
Script used to test a single microring filter that consists of multiple rings.

Author      : Simon Bélanger-de Villers (simon.belanger-de-villers.1@ulaval.ca)
Created     : July 24th 2019
Last edited : July 24th 2019
"""
import sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from components.couplers.DC_from_data import DC_from_data
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import flatTopCoupling
from components.MRF import MRF
import numpy as np

# Parameters for the rings
num_rings 	= 5 		# Number of cascaded rings [-]
radius 		= 2.5e-6 	# Radius of the ring resonators [m]
neff 		= 2.4449 	# Effective index of the ring waveguide []
ng 			= 4.18 		# Group index of the ring waveguide []
alpha_wg	= 3. 		# Losses of the ring waveguide [dB/cm]

# Parameters for the couplers (num_rings+1)
k_vec = flatTopCoupling(num_rings, 0.5)
loss_coupler = [0.98,0.99,0.99,0.99,0.99,0.99]
couplers = [virtual_DC(k_vec[0],loss_coupler[0]),  
			virtual_DC(k_vec[1],loss_coupler[1]), 
			virtual_DC(k_vec[2],loss_coupler[2]), 
			virtual_DC(k_vec[3],loss_coupler[3]), 
			virtual_DC(k_vec[4],loss_coupler[4]), 
			virtual_DC(k_vec[5],loss_coupler[5])]


# Build the Microring filter and measure transmission across the wavelength range
wavelength=np.linspace(1530e-9, 1570e-9, 1000)
mrf = MRF(num_rings, radius, neff, ng, alpha_wg, couplers, [1, 0, 0])
mrf.clipping=False
outFields = mrf.TMM(wavelength, np.ones(len(wavelength)))
mrf.plotTransmission(wavelength, outFields)