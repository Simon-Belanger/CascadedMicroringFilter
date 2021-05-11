"""
Plot the response of a 5 rings filter for increasing delta phi random
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
import matplotlib.font_manager as fm
from math import pi

# General parameters
pdfFilename = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter2/singleRing.pdf'


# Build the dual ring filter
filterOrder 	= 5 		# Number of cascaded rings [-]
radius 		    = 2.5e-6 	# Radius of the ring resonators [m]
neff 		    = 2.4449 	# Effective index of the ring waveguide []
ng 			    = 4.18 		# Group index of the ring waveguide []
alpha_wg	    = 100. 		# Losses of the ring waveguide [dB/cm]
powerCoupling   = 0.4       # Power coupling to the bus []
mrf = idealMRF(filterOrder, radius, neff, ng, alpha_wg, powerCoupling)
mrf.clipping = True

wavelength=np.linspace(1530e-9, 1550e-9, 1000)
for phiR in [0, pi/16, pi/8, pi/4, pi/2]:
    # Randomize it's response
    mrf.manufacturing(phiR)
    outFields = mrf.TMM(wavelength, np.ones(len(wavelength)))
    #mrf.plotTransmission(wavelength, outFieldsNV)


    ra = ringAnalyisGeneral(wavelength, transmissionDrop=20*np.log10(outFields.drop))      # Normalise the data using ringAnalysis (findpeaks)
    cw = ra.measureLongitudinalModes('wavelength')
    plt.plot(wavelength, 20*np.log10(outFields.drop))
plt.show()



if False:
    fig, (leftAx, rightAx) = plt.subplots(ncols=2, figsize=(10, 4), sharey=True, gridspec_kw={'width_ratios': [2, 1]})
    leftAx.set_ylabel('Insertion loss [dB]', font=overpassFont)
    # Left Ax : Tuning spectrum
    leftAx.plot(lm5U.read_param('wvl')*1e9, lm5U.applySmoothingFilter(lm5U.read_param('il')), color='grey', label='As-fabricated')
    leftAx.plot(lm5T.read_param('wvl')*1e9, lm5T.applySmoothingFilter(lm5T.read_param('il')), color='crimson', label='After tuning')
    leftAx.set_xlim([1550, 1592])
    leftAx.set_xlabel('Wavelength [nm]', font=overpassFont)
    leftAx.grid(color='grey', linestyle='--', linewidth=0.25)
    leftAx.legend(prop=overpassFont, frameon=False, loc='upper right', bbox_to_anchor=(0.92, 0.9))
    # Right Ax : Relative wavelength
    raU = ringAnalyisGeneral(lm5U.read_param('wvl'), lm5U.applySmoothingFilter(lm5U.read_param('il')))      # Normalise the data using ringAnalysis (findpeaks)
    raT = ringAnalyisGeneral(lm5T.read_param('wvl'), lm5T.applySmoothingFilter(lm5T.read_param('il')))      # Normalise the data using ringAnalysis (findpeaks)
    cwu, cwt = raU.measureLongitudinalModes('wavelength'), raT.measureLongitudinalModes('wavelength')
    print(raT.measureLongitudinalModesOrder(2.4449, 5e-6))
    rightAx.plot(lm5U.read_param('wvl')*1e9-cwu[2]*1e9, lm5U.applySmoothingFilter(lm5U.read_param('il')), color='grey')
    rightAx.plot(lm5T.read_param('wvl')*1e9-cwt[1]*1e9+ajustment, lm5T.applySmoothingFilter(lm5T.read_param('il')), color='crimson')    
    rightAx.set_xlim([-3, 3])
    rightAx.set_xlabel('Relative wavelength, ( $\lambda-\lambda_m$ ) [nm]', font=overpassFont)
    rightAx.grid(color='grey', linestyle='--', linewidth=0.25)
    # Second x axis
    rightAxF = rightAx.secondary_xaxis('top', functions=(deltaNmToGHz, deltaGHzToNm))
    rightAxF.set_xlabel('Relative frequency, ( $\\nu-\\nu_m$ ) [GHz]', font=overpassFont)
    # Rectangle for inset
    leftAx.add_patch(Rectangle((cwt[1]*1e9-4/2, leftAx.get_ylim()[0]), 4, leftAx.get_ylim()[1]-leftAx.get_ylim()[0], fill=None, linestyle='-.'))
    plt.subplots_adjust(wspace=.1)
    plt.savefig(pdfDirectory + 'roadm5um_tuning.pdf')
    plt.show()