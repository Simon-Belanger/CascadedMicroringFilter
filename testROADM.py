"""
This File is used to simulate the response of a 4 channels ROADM made with 
cascaded microring filters.

Author      : Simon Bélanger-de Villers
Created     : July 24th 2019
Last edited : July 24th 2019
"""
import numpy as np 
import matplotlib 
matplotlib.use("Qt5agg") 
import matplotlib.pyplot as plt

# Configure one microring
from components.couplers.DC_from_data import DC_from_data
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import flattop5
from components.MRF import MRF

num_rings 	= 5 		# Number of cascaded rings [-]
radius 		= 2.5e-6 	# Radius of the ring resonators [m]
neff 		= 2.4449 	# Effective index of the ring waveguide []
ng 			= 4.18 		# Group index of the ring waveguide []
alpha_wg	= 3. 		# Losses of the ring waveguide [dB/cm]
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

# Input Signal
wavelength  = np.linspace(1500e-9, 1600e-9, 1000)
inputPower  = np.ones(1000)
plt.plot(wavelength*1e9, inputPower)
plt.xlabel('Wavelength [nm]')
plt.ylabel('Power [a.u]')
plt.show()

# Filterd by filter #1

# Filterd by filter #2

# Filterd by filter #3

# Filterd by filter #4

# Through component
