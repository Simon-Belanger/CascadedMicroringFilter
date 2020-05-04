""" All-Pass Microring Modulator analytical model from W. Bogaerts paper 'Silicon microring resonators' ."""

import sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from analytical.MRR_AP import *
from misc.dataIO import saveData, loadData
import os, math
import numpy as np
import matplotlib.pyplot as plt

## Things to do
# TODO: Group index from mode solver/ effective index method
# TODO: Evaluate the losses from bends + propag + pn junction
# TODO: Corner Analysis 
# TODO: Photoconcuctive heater in the model
# TODO: Interpolate the data from the triangular mesh on a rectilinear grid and apply the effective index method (https://matplotlib.org/gallery/images_contours_and_fields/triinterp_demo.html#sphx-glr-gallery-images-contours-and-fields-triinterp-demo-py)
# TODO: Critical coupling is for a given bias (Operation point)? Or is it measured for 0V Bias?
# TODO: Doped waveguide class that has a waveguide and a phase shifter cross section

class MRM_Static(MRR_AP):
    """ Static model of a microring modulator. 

    Example : see mainMRM.py
    """
    wavelengthRange = {'start':1552e-9, 'stop':1553e-9, 'pts':1000}     # Wavelength range for sweeps
    biasSweep       = np.linspace(0, -4.0, 5)                                # Bias array for sweep
    offsetGap       = 10e-9           # Variation in gap in order to achieve undercoupling/overcoupling [m]
    wavelength      = 1550e-9

    def __init__(self, radius, loss_per_cm, n_eff, n_g, gap, pnCoverage, coupler, phaseshifter):

        if radius != 10e-6:
            print('Error: the coupler data is only for radius = 10 microns.');exit()

        # Properties
        self.radius = radius
        
        # Components
        self.coupler    = coupler
        self.phaseShifter = phaseshifter

        # Waveguides
        self.pnCoverage     = pnCoverage
        self.alpha          = attenuationCoefficientFromWaveguideLosses(loss_per_cm, type='power')
        self.waveguide      = waveguideMRM(effectiveIndex=n_eff, groupIndex=n_g, attenuationCoefficient=self.alpha, length=self.roundtripLength-self.modLength)
        self.dopedWaveguide = waveguideMRM(effectiveIndex=self.phaseShifter.intrinsicEffectiveIndex, groupIndex=n_g, attenuationCoefficient=self.phaseShifter.intrinsicAttenuationCoefficient, length=self.modLength)

        self.gap = gap # Needs to be after the phase shifter definition in case CC is asked it needs to know the losses

    # Properties / Attributes
    @property
    def modLength(self):
        " Getter for the PN Junction modulator length [m]. "
        return self.roundtripLength * self.pnCoverage / 100

    @property
    def gap(self):
        " Getter for the directional coupler gap [m]. "
        return self._gap
    @gap.setter
    def gap(self, gap):
        " Setter for the directional coupler gap [m]. "
        self.setBias(0) # TODO:  To have the losses of the modulator included (intrinsic vs perturbative)
        criticalGap = self.findCriticalCouplingCondition()
        if gap == 'CC': # Critically Coupled
            self._gap = criticalGap
        elif gap == 'UC': # Under Coupled
            self._gap = criticalGap + self.offsetGap
        elif gap == 'OC': # Over Coupled
            self._gap = criticalGap - self.offsetGap
        else:
            self._gap = gap
        self.setCrossCoupling(self._gap)

    @property
    def roundtripTransmission(self):
        """ Obtain the roundtrip losses from the loss coefficient and the roundtrip length of the ring resonator. """
        return self.waveguide.getAmplitudeTransmission() * self.dopedWaveguide.getAmplitudeTransmission()

    @property
    def roundtripPhase(self):
        """ Obtain the roundtrip phase at a particular wavelength from the effective index and the roundtrip length of
        the ring modulator. """
        return self.waveguide.getPhaseShift(self.wavelength) + self.dopedWaveguide.getPhaseShift(self.wavelength)

    @property
    def freeSpectralRange(self):
        return self.wavelength**2/((self.waveguide.groupIndex * self.waveguide.length)+(self.dopedWaveguide.groupIndex * self.dopedWaveguide.length))

    # Methods
    def setCrossCoupling(self, gap):
        """ Set the cross coupling coefficient in power. """
        self.C = self.coupler.measureCouplingCoefficient(gap, self.wavelength)

    def setBias(self, bias):
        """ Set the reverse bias applied to the PN junction. """

        # Interpolate the complex effective index to the selected reverse bias
        neffInterp = np.interp(-bias, -self.phaseShifter.phaseShifter['V'], self.phaseShifter.phaseShifter['neff'])

        # Set the new effective index and losses
        self.dopedWaveguide.effectiveIndex = neffInterp.real
        loss_per_cm = -0.2 * np.log10(np.exp(1))*(-2 * math.pi * neffInterp.imag/self.wavelength) # [dB/cm]
        self.dopedWaveguide.attenuationCoefficient = attenuationCoefficientFromWaveguideLosses(loss_per_cm, type='power')

    def findCriticalCouplingCondition(self):
        """ Find the gap for which the cross-coupling coefficient satisfies the critical coupling condition. """
        return self.coupler.reverseInterpolateCouplingCoefficient((1 - self.roundtripTransmission**2))

    def plotTransmissionVsBias(self):
        " Measure the MRM transmission vs DV Reverse Bias "
        plt.figure(1);plt.title('Transmission');plt.xlabel('Wavelength [nm]');plt.ylabel('Transmission [dB]')
        plt.figure(2);plt.title('Phase');plt.xlabel('Wavelength [nm]');plt.ylabel('Phase [rad]')
        for Vi in self.biasSweep:
            self.setBias(Vi)
            wvl, power = self.sweep_power_transmission(self.wavelengthRange['start'], self.wavelengthRange['stop'], self.wavelengthRange['pts'])
            plt.figure(1);plt.plot(wvl*1e9, power, label='{:0.1f} V'.format(Vi))
            wvl, phase = self.sweep_phase(self.wavelengthRange['start'], self.wavelengthRange['stop'], self.wavelengthRange['pts'])
            plt.figure(2);plt.plot(wvl*1e9, phase, label='{:0.1f} V'.format(Vi))
        plt.figure(1);plt.legend();plt.grid();plt.xlim(self.wavelengthRange['start']*1e9, self.wavelengthRange['stop']*1e9);plt.savefig('transmission.pdf')
        plt.figure(2);plt.legend();plt.grid();plt.xlim(self.wavelengthRange['start']*1e9, self.wavelengthRange['stop']*1e9);plt.savefig('phase.pdf')
        plt.show()

    def measureOpticalModulationAmplitude(self):
        pass

    def getQfactor(self, lambdaRes=1550e-9):
        """ Return the Q factor for the cavity for a given resonance wavelength. """
        return (math.pi * self.waveguide.groupIndex * self.roundtripLength * math.sqrt(self.r * self.roundtripTransmission))/(lambdaRes * (1 - self.r * self.roundtripTransmission))

    def measureOpticalBandwidthVsBias(self):
        """ Measure the optical bandwidth of the modulator [GHz]. """
        BW = []
        for Vi in self.biasSweep:
            self.setBias(Vi)
            BW.append(self.getOpticalBandwidth(self.wavelength))
        return BW # [GHz]

    def plotOpticalBandwidthVsBias(self):
        """ Plot the optical bandwidth vs reverse bias. """
        plt.plot(self.biasSweep, self.measureOpticalBandwidthVsBias())
        plt.xlabel('Reverse Bias [V]');plt.ylabel('Optical Bandwidth [GHz]')
        plt.grid();plt.savefig('opticalBandwidth.pdf');plt.show()