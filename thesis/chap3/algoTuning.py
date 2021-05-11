"""
Experiment with central wavelength tuning of the MRF
"""

import numpy as np 
import matplotlib.pyplot as plt
import math, sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.Filter import *
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import *
from components.MRF import idealMRF
from components.ring import Ring
from active_tuning.Algorithms import *
import matplotlib.font_manager as fm
from math import pi
import random
import pickle

# Build the dual ring filter/
filterOrder 	= 2 		# Number of cascaded rings [-]
radius 		    = 2.5e-6 	# Radius of the ring resonators [m]
neff 		    = 2.4449 	# Effective index of the ring waveguide []
ng 			    = 4.18 		# Group index of the ring waveguide []
alpha_wg	    = 3. 		# Losses of the ring waveguide [dB/cm]
powerCoupling   = 0.1       # Power coupling to the bus []
crosstalkMat    = [1., 0.5]  # Crosstalk matrix indices [normalized]
mrf = idealMRF(filterOrder, radius, neff, ng, alpha_wg, powerCoupling, crosstalkMat)
mrf.clipping = False

# Random center but well defined delta 
maximumInterval = pi/2  
delta           = pi/4
center = random.uniform(-maximumInterval/2, maximumInterval/2)
mrf.manufacturingValue([center-delta, center+delta])


# Plot the transmission spectrum before wavelength tuning
wavelength  = np.linspace(1500e-9, 1600e-9, 1000)
outFields   = mrf.TMM(wavelength, np.ones(len(wavelength)))
mrf.plotTransmission(wavelength, outFields)


# Tuning Coordinates method
cdr = coordinatesDescent(mrf, numIter=5, plot_maps=False)
#pickle.dump( cdr, open( "cdrTest.p", "wb" ) )
#cdr = pickle.load( open( "cdrTest.p", "rb" ) )


# Plot the convergence for each ring / iteration
#plotCoordinatesDescentIteration(cdr)


print(np.matrix([[1, 0.15, 0],[0.15, 1, 0.15],[0, 0.15, 1]])*np.matrix([[-1, 1, 1],[0, sqrt(2), -sqrt(2)],[1, 1, 1]]))
print(np.matrix([[-1, 1, 1],[0, sqrt(2), -sqrt(2)],[1, 1, 1]])*np.matrix([[1, 0.15, 0],[0.15, 1, 0.15],[0, 0.15, 1]]))

#plotConvergence(cdr['phase'], cdr['power'])