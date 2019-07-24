from components.couplers.DC_from_data import DC_from_data
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import flattop5
from components.MRF import MRF
import numpy as np

# Parameters for the rings
num_rings 	= 5 		# Number of cascaded rings [-]
radius 		= 2.5e-6 	# Radius of the ring resonators [m]
neff 		= 2.4449 	# Effective index of the ring waveguide []
ng 			= 4.18 		# Group index of the ring waveguide []
alpha_wg	= 3. 		# Losses of the ring waveguide [dB/cm]

# Parameters for the couplers (num_rings+1)
k_vec = flattop5(0.1)
loss_coupler = [0.98,0.99,0.99,0.99,0.99,0.99]
couplers = [virtual_DC(k_vec[0],loss_coupler[0]),  
			virtual_DC(k_vec[1],loss_coupler[1]), 
			virtual_DC(k_vec[2],loss_coupler[2]), 
			virtual_DC(k_vec[3],loss_coupler[3]), 
			virtual_DC(k_vec[4],loss_coupler[4]), 
			virtual_DC(k_vec[5],loss_coupler[5])]


# Build the Microring filter and measure transmission across the wavelength range
mrf = MRF('', num_rings, radius, neff, ng, alpha_wg, couplers, [1, 0, 0])
mrf.sweep(np.linspace(1500,1600,1000)*1e-9)