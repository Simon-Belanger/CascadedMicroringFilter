import numpy as np 
import matplotlib.pyplot as plt
import math, sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from components.MRF import idealMRF

"""
Test used to make sure the ports attribution in the TMM method is working right no matter which number of rings.

Author      : Simon Belanger-de Villers (simon.belanger-de-villers.1@ulaval.ca)
Created     : July 24th 2019
Last edited : August 13th 2020
"""

# Input Signal for all tests
wavelength  = np.linspace(1540e-9, 1545e-9, 1000)
inputPower  = np.ones(1000)

# One ring
mrf = idealMRF(1, 2.5e-6, 2.4449, 4.18, 3., 0.1, [1])
outFields = mrf.TMM(wavelength, inputPower)
print('1 ring')
mrf.plotTransmission(wavelength, outFields)

# 2 Rings
mrf = idealMRF(2, 2.5e-6, 2.4449, 4.18, 3., 0.1, [1])
outFields = mrf.TMM(wavelength, inputPower)
print('2 rings')
mrf.plotTransmission(wavelength, outFields)

# 3 Rings
mrf = idealMRF(3, 2.5e-6, 2.4449, 4.18, 3., 0.1, [1])
outFields = mrf.TMM(wavelength, inputPower)
print('3 rings')
mrf.plotTransmission(wavelength, outFields)

# 4 Rings
mrf = idealMRF(4, 2.5e-6, 2.4449, 4.18, 3., 0.1, [1])
outFields = mrf.TMM(wavelength, inputPower)
print('4 rings')
mrf.plotTransmission(wavelength, outFields)

# 5 rings
mrf = idealMRF(5, 2.5e-6, 2.4449, 4.18, 3., 0.1, [1])
mrf.manufacturing(math.pi/50)
outFields = mrf.TMM(wavelength, inputPower)
print('5 rings')
mrf.plotTransmission(wavelength, outFields)
mrf.plotPhaseResponse(wavelength, mrf.measurePhaseResponse(wavelength, outFields)[0], mrf.measurePhaseResponse(wavelength, outFields)[1])
mrf.plotGroupDelay(wavelength, mrf.measureGroupDelay(wavelength, outFields)[0], mrf.measureGroupDelay(wavelength, outFields)[1])