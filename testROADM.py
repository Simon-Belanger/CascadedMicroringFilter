"""
This File is used to simulate the response of a 4 channels ROADM made with 
cascaded microring filters.

Author      : Simon Bélanger-de Villers (simon.belanger-de-villers.1@ulaval.ca)
Created     : July 24th 2019
Last edited : July 24th 2019
"""
import numpy as np 
import matplotlib.pyplot as plt
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB
import math

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
loss_coupler = [0.99,0.99,0.99,0.99,0.99,0.99]
couplers = [virtual_DC(k_vec[0], loss_coupler[0]),  
			virtual_DC(k_vec[1], loss_coupler[1]), 
			virtual_DC(k_vec[2], loss_coupler[2]), 
			virtual_DC(k_vec[3], loss_coupler[3]), 
			virtual_DC(k_vec[4], loss_coupler[4]), 
			virtual_DC(k_vec[5], loss_coupler[5])]


# Build the Microring filter and measure transmission across the wavelength range
mrf = MRF('', num_rings, radius, neff, ng, alpha_wg, couplers, [1, 0, 0])

# Input Signal
wavelength  = np.linspace(1540e-9, 1560e-9, 1000)
inputPower  = np.zeros(1000)

#TODO Edge coupler loss and wavelength dependance
#TODO reconfigure MRF so it accepts input for the thru/drop ports instead of just 0
#TODO PSR to split into input and Add
#TODO PSR ro combine add and drop
#TODO Circulator to have thru/drop distinct from Input/Add
#TODO Model Power taps 
#TODO Model VOAs

# Filterd by filter #1
[E_drop1, E_thru1] = mrf.TMM_v2(wavelength, powerdB2Field(inputPower), np.zeros(1000))
#mrf.plot_transmission(wavelength, E_drop1, E_thru1)

# Filterd by filter #2
mrf.apply_phase_tuning([math.pi/4, math.pi/4, math.pi/4, math.pi/4, math.pi/4])
[E_drop2, E_thru2] = mrf.TMM_v2(wavelength, E_thru1, np.zeros(1000))
#mrf.plot_transmission(wavelength, E_drop2, E_thru2)

# Filterd by filter #3
mrf.apply_phase_tuning([2*math.pi/4, 2*math.pi/4, 2*math.pi/4, 2*math.pi/4, 2*math.pi/4])
[E_drop3, E_thru3] = mrf.TMM_v2(wavelength, E_thru2, np.zeros(1000))
#mrf.plot_transmission(wavelength, E_drop3, E_thru3)

# Filterd by filter #4
mrf.apply_phase_tuning([3*math.pi/4, 3*math.pi/4, 3*math.pi/4, 3*math.pi/4, 3*math.pi/4])
[E_drop4, outputField] = mrf.TMM_v2(wavelength, E_thru3, np.zeros(1000))
#mrf.plot_transmission(wavelength, E_drop4, outputField)

# Device output transmission spectrum
plt.plot(wavelength*1e9, inputPower, label="Input")
plt.plot(wavelength*1e9, field2PowerdB(E_drop1), label="Channel 1")
plt.plot(wavelength*1e9, field2PowerdB(E_drop2), label="Channel 2")
plt.plot(wavelength*1e9, field2PowerdB(E_drop3), label="Channel 3")
plt.plot(wavelength*1e9, field2PowerdB(E_drop4), label="Channel 4")
plt.plot(wavelength*1e9, field2PowerdB(outputField), label="Through")
plt.xlabel('Wavelength [nm]')
plt.ylabel('Power [a.u]')
plt.legend()
plt.show()