"""
This File is used to simulate the response of a 4 channels ROADM made with 
cascaded microring filters.

Author      : Simon Bélanger-de Villers (simon.belanger-de-villers.1@ulaval.ca)
Created     : July 24th 2019
Last edited : October 17th 2019
"""
# Things to do when I have time
#TODO Edge coupler loss and wavelength dependance
#TODO reconfigure MRF so it accepts input for the thru/drop ports instead of just 0
#TODO PSR to split into input and Add
#TODO PSR to combine add and drop
#TODO Circulator to have thru/drop distinct from Input/Add
#TODO Model Power taps 
#TODO Model VOAs
#TODO Algorithm to tune it 
#TODO make instances of the filter instead of using the same exact one
#TODO Easier way to tune the filters 
#TODO Phase difference with waveguide imbalance + tuning
#TODO Use EMPy integrated mode solver to perform waveguide characterization
#TODO Implement the couplers from JLT Bahadori with EMPy mode solver

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
k_vec 		= flattop5(0.1)
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
nonepower = np.zeros(1000, dtype=complex)

# Filterd by filter #1
outFields1 = mrf.TMM(wavelength, powerdB2Field(inputPower))

# Filterd by filter #2
mrf.apply_phase_tuning([math.pi/4, math.pi/4, math.pi/4, math.pi/4, math.pi/4])
outFields2 = mrf.TMM(wavelength, outFields1.through)

# Filterd by filter #3
mrf.apply_phase_tuning([2*math.pi/4, 2*math.pi/4, 2*math.pi/4, 2*math.pi/4, 2*math.pi/4])
outFields3 = mrf.TMM(wavelength, outFields2.through)

# Filterd by filter #4
mrf.apply_phase_tuning([3*math.pi/4, 3*math.pi/4, 3*math.pi/4, 3*math.pi/4, 3*math.pi/4])
outFields4 = mrf.TMM(wavelength, outFields3.through)

# Device output transmission spectrum
plt.plot(wavelength*1e9, inputPower, label="Input")
plt.plot(wavelength*1e9, field2PowerdB(outFields1.drop), label="Channel 1")
plt.plot(wavelength*1e9, field2PowerdB(outFields2.drop), label="Channel 2")
plt.plot(wavelength*1e9, field2PowerdB(outFields3.drop), label="Channel 3")
plt.plot(wavelength*1e9, field2PowerdB(outFields4.drop), label="Channel 4")
plt.plot(wavelength*1e9, field2PowerdB(outFields4.through), label="Through")
plt.xlabel('Wavelength [nm]')
plt.ylabel('Power [dB]')
plt.legend()
plt.show()