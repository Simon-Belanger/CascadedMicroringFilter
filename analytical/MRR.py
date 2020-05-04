""" General Microring Resonator analytical model from W. Bogaerts paper 'Silicon microring resonators' ."""

from math import cos,sin,atan,pi,exp,sqrt
import numpy as np
import cmath, math
import matplotlib.pyplot as plt
import sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.utils import attenuationCoefficientFromWaveguideLosses
from misc.constants import c
from components.waveguides.waveguideMRM import waveguideMRM

class MRR(object):
    """ General Microring Resonator."""

    refWavelength   = 1550e-9 # Reference wavelength for the computation of the propagation constant (taylor expansion) [m]

    def __init__(self, radius, loss_per_cm, n_eff, n_g):
        """ Constructor for the ring resonator class. """

        # Main properties
        self.wavelength = 1550e-9       # Wavelength @ which the transmission will be measured [m]
        self.radius     = radius        # Ring radius [m]
        self.loss       = loss_per_cm   # Waveguide losses [dB/cm]
        self.n_eff      = n_eff         # Effective index of the waveguide []
        self.n_g        = n_g           # Group index of the waveguide []
        # Measured properties
        self.alpha = attenuationCoefficientFromWaveguideLosses(loss_per_cm, type='power')
        self.waveguide = waveguideMRM(effectiveIndex=n_eff, groupIndex=n_g, attenuationCoefficient=self.alpha, length=self.roundtripLength)

    @property
    def wavelength(self):
        """ Wavelength at which the transmission/phase will be measured. """
        return self._wavelength
    @wavelength.setter
    def wavelength(self, wavelength):
        self._wavelength = wavelength

    @property
    def radius(self):
        """ Radius of the microring resonator. """
        return self._radius
    @radius.setter
    def radius(self, radius):
        self._radius = radius

    @property
    def roundtripLength(self):
        """ Waveguide length in the resonator. """
        return 2 * pi * self.radius

    @property
    def roundtripLoss(self):
        """ Obtain the field roundtrip losses from the loss coefficient and the roundtrip length of the ring resonator. 
            Rename to roundtripTransmission. """
        return self.waveguide.getAmplitudeTransmission()

    @property
    def roundtripPhase(self):
        """ Obtain the roundtrip phase at a particular wavelength from the effective index and the roundtrip length of
        the ring resonator. """
        return self.waveguide.getPhaseShift(self.wavelength)

    def get_field_transmission(self, wavelength, E_in):
        """ This method will be defined in the inherited classes 'MRR_AP' and 'MRR_AD'. """
        pass

    def sweep_power_transmission(self, lambda_min, lambda_max, lambda_points):
        """ Obtain the power transmission spectrum of the ring resonator. Implemented in subclasses. """
        pass

    def sweep_field_transmission(self, lambda_min, lambda_max, lambda_points):
        """ Obtain the field transmission spectrum of the ring resonator."""
        lambda_0 = np.linspace(lambda_min, lambda_max, lambda_points)
        E_t = []
        for wvl in lambda_0:
            E_t.append(self.get_field_transmission(wvl, 1))
        return lambda_0, np.asarray(E_t)

    def plot_field_transmission(self, lambda_0, E_t):
        """ Plot the drop/thru port transmission spectrum. """
        plt.plot(lambda_0 * 1e9, 10 * np.log10(abs(E_t) ** 2))
        plt.xlabel('Wavelength (nm)', fontsize=14)
        plt.ylabel('Power transmission (dB)', fontsize=14)
        plt.show()

    def plot_phase(self, lambda_0, E_t):
        """ Plot the drop/thru port transmission spectrum. """
        plt.plot(lambda_0 * 1e9, np.unwrap(np.angle(E_t))/ np.pi)
        plt.xlabel('Wavelength (nm)', Fontsize=14)
        plt.ylabel('Phase (/pi)', Fontsize=14)
        plt.show()

    def sweep_phase(self, lambda_min, lambda_max, lambda_points):
        """ Obtain the field transmission spectrum of the ring resonator."""
        lambda_0 = np.linspace(lambda_min, lambda_max, lambda_points)
        Phi = []
        for wvl in lambda_0:
            Phi.append(self.get_effective_phase_shift(wvl))
        return lambda_0, np.asarray(Phi)

    def measure_power_thru(self, E_t):
        return 10 * np.log10(abs(E_t) ** 2)

    def measure_phase_thru(self, E_t):
        return np.angle(E_t)/ np.pi

    @staticmethod
    def get_r(C):
        """ Obtain the field thru-coupling coefficient 'r' from the power cross-coupling coefficient 'C'. """
        return sqrt(1 - C)

    def getQfactor(self):
        """ This method will be defined in the inherited classes 'MRR_AP' and 'MRR_AD'. """
        pass

    def getOpticalBandwidth(self, wavelength=1550e-9):
        """ This method returns the optical bandwidth of a microring resonator that is dependant on the Q factor. [GHz]"""
        omega_0 = 2 * pi * c / wavelength
        return omega_0 / (2 * math.pi * self.getQfactor()) * 1e-9 # [GHz]