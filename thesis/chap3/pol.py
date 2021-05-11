"""
Analysis of polarization of the PSR sagnac loop
"""
import sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB, attenuationCoefficientFromWaveguideLosses
from misc.constants import c
from math import pi
import cmath as cm
import numpy as np
import matplotlib.pyplot as plt

# Methods
def getPropagationConstant(wavelength, effectiveIndex, groupIndex, dispersion):
    " Getter for the propagation constant, dispersive. "
    refWavelength = 1550e-9

    dOmega = ((2*pi*c)/float(wavelength) - (2*pi*c)/refWavelength)
    firstOrder  = (2*pi*c)/refWavelength * (float(effectiveIndex)/c)
    secondOrder = float(groupIndex)/c
    thirdOrder  = 0.5 * (- float(dispersion) * 1e-6 * refWavelength**2 / 2 * pi * c )
    return firstOrder + secondOrder * dOmega + thirdOrder * dOmega**2

# Properties of the waveguide for both TE and TM modes
modeTE = {'effectiveIndex': 2.44, 'groupIndex': 4.18, 'dispersion': 0, 'wgLosses': 3} #530
modeTM = {'effectiveIndex': 1.78, 'groupIndex': 3.805, 'dispersion': 0, 'wgLosses': 2} # -2.1e4


# PSR transfer coefficients
TEBranch = {'TE-TE': -1, 'TE-TM': -1000, 'TM-TE': -1000,'TM-TM': -1000}
TMBranch = {'TE-TE': -1000, 'TE-TM': -1000, 'TM-TE': -1, 'TM-TM':-1000}

inputField      = 1         # Electric field amplitude [V/m]
inputTEratio    = 0.5       # TE ratio of the input field (0 if purely TM, 1 if purely TE and in between for some of both)
wavelengthArray = np.linspace(1500e-9, 1600e-9, 1000)
initLength      = 200e-6
length          = 2500e-6    # Length of the sagnac loop [m]

output = {'TE': [], 'TM': [], 'Total': []}
beta = {'TE':[], 'TM':[]}
test = []
for wavelength in wavelengthArray:
    # Mode propagation properties
    betaTE = getPropagationConstant(wavelength, modeTE['effectiveIndex'], modeTE['groupIndex'], modeTE['dispersion'])
    betaTM = getPropagationConstant(wavelength, modeTM['effectiveIndex'], modeTM['groupIndex'], modeTM['dispersion'])
    beta['TE'].append(betaTE); beta['TM'].append(betaTM)
    alphaTE = attenuationCoefficientFromWaveguideLosses(modeTE['wgLosses'], 'field')
    alphaTM = attenuationCoefficientFromWaveguideLosses(modeTM['wgLosses'], 'field')


    # Input fields
    inputTE = inputField * cm.exp(-1j*0) * inputTEratio
    inputTM = inputField * cm.exp(-1j*0) * (1-inputTEratio)

    # Initial propagation
    inputTE = inputTE * cm.exp(-1j * (betaTE - 1j * alphaTE) * initLength)
    inputTM = inputTM * cm.exp(-1j * (betaTM - 1j * alphaTM) * initLength)

    # After the polarization splitting
    topTE = inputTE * powerdB2Field(TEBranch['TE-TE']) + inputTM * powerdB2Field(TEBranch['TM-TE'])
    topTM = inputTE * powerdB2Field(TEBranch['TE-TM']) + inputTM * powerdB2Field(TEBranch['TM-TM'])
    botTE = inputTE * powerdB2Field(TMBranch['TE-TE']) + inputTM * powerdB2Field(TMBranch['TM-TE'])
    botTM = inputTE * powerdB2Field(TMBranch['TE-TM']) + inputTM * powerdB2Field(TMBranch['TM-TM'])

    # Propagation
    betaTE = getPropagationConstant(wavelength, modeTE['effectiveIndex'], modeTE['groupIndex'], modeTE['dispersion'])
    betaTM = getPropagationConstant(wavelength, modeTE['effectiveIndex'], modeTE['groupIndex'], modeTE['dispersion'])
    alphaTE = attenuationCoefficientFromWaveguideLosses(modeTE['wgLosses'], 'field')
    alphaTM = attenuationCoefficientFromWaveguideLosses(modeTE['wgLosses'], 'field')

    topTE = topTE * cm.exp(-1j * (betaTE - 1j * alphaTE) * length)
    topTM = topTM * cm.exp(-1j * (betaTE - 1j * alphaTE) * length)
    botTE = botTE * cm.exp(-1j * (betaTE - 1j * alphaTE) * length)
    botTM = botTM * cm.exp(-1j * (betaTE - 1j * alphaTE) * length)

    # After the polarization recombining
    outputTE = topTE * powerdB2Field(TMBranch['TE-TE']) + topTM * powerdB2Field(TMBranch['TE-TM']) + botTE * powerdB2Field(TEBranch['TE-TE']) + botTM * powerdB2Field(TEBranch['TE-TM'])
    outputTM = topTE * powerdB2Field(TMBranch['TM-TE']) + topTM * powerdB2Field(TMBranch['TM-TM']) + botTE * powerdB2Field(TEBranch['TE-TE']) + botTM * powerdB2Field(TEBranch['TE-TM'])
    outputTE = outputTE * cm.exp(-1j * (betaTE - 1j * alphaTE) * initLength)
    outputTM = outputTM * cm.exp(-1j * (betaTM - 1j * alphaTM) * initLength)

    # Output fields
    output['TE'].append(field2PowerdB(outputTE))
    output['TM'].append(field2PowerdB(outputTM))
    output['Total'].append(field2PowerdB(outputTE + outputTM))

#plt.plot(wavelengthArray*1e9, np.asarray(test),label='test')
plt.plot(wavelengthArray*1e9, np.asarray(output['TE']),label='output TE')
plt.plot(wavelengthArray*1e9, np.asarray(output['TM']),label='output TM')
plt.plot(wavelengthArray*1e9, np.asarray(output['Total']),label='Interference')
plt.plot(wavelengthArray*1e9, np.asarray(output['TE'])+np.asarray(output['TM']),label='TE+TM')
plt.axis([1500, 1600, -100, 0])
plt.legend()
plt.show()

# Beta plot
#x = wavelengthArray*1e9
##x = c/wavelengthArray*1e-12
#plt.plot(x, beta['TE'], label='TE')
#plt.plot(x, beta['TM'], label='TM')
#plt.grid()
#plt.show()