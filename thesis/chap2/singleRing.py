"""
Show the transmission and phase of a single ring resonator filter.
"""

import numpy as np 
import matplotlib.pyplot as plt
import math, sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB, field2PowerLinear
from misc.Filter import *
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import *
from components.MRF import parentMRF
from components.ring import Ring
import matplotlib.font_manager as fm

# General parameters
pdfFilename = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter2/singleRing.pdf'
wavelength=np.linspace(1500e-9, 1600e-9, 10000)

"""
Case 1 : 2 rings notVernier filter
"""
num_rings 	= 1 		# Number of cascaded rings [-]
radius 		= 2.5e-6 	# Radius of the ring resonators [m]
neff 		= 2.4449 	# Effective index of the ring waveguide []
ng 			= 4.18 		# Group index of the ring waveguide []
alpha_wg	= 3. 		# Losses of the ring waveguide [dB/cm]

rings    = [Ring(radius, neff, ng, alpha_wg) for i in range(num_rings)]
k_vec = [0.1, 0.1]
loss_coupler = [1, 1]
couplers = [virtual_DC(k_vec[0],loss_coupler[0]),
			virtual_DC(k_vec[1],loss_coupler[1])]
mrf = parentMRF(rings, couplers)
mrf.clipping = False
outFieldsNV = mrf.TMM(wavelength, np.ones(len(wavelength)))

mrf.plotTransmission(wavelength, outFieldsNV)