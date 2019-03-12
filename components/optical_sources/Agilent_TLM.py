from scipy.signal import gaussian
import matplotlib.pyplot as plt
from math import sqrt, pi, exp, log
import numpy as np

class Agilent_TLM(object):
    """
    Class to generate the spectral response of the Agilent 81600B Tunable Laser Module (TLM).
    With associated noise and perturbative effects.
    """

    output_power        = 0         # Laser output power [dBm]

    central_wavelength  = 1550e-9   # Selected central wavelength for the laser [m]
    min_wavelength      = 1440e-9   # Maximum wavelength. [m]
    max_wavelength      = 1640e-9   # Minimum wavelength. [m]
    res_wavelength     = 0.1e-12  # Wavelength resolution: smallest selectable wavelength increment/decrement. [m]

    linewidth           = 100e3     # Static 3dB linewidth (FWHM) [Hz] 100 kHz
    linewidth_eff       = 50e6      # Dynamic 3dB linewidth [Hz] 50 MHz

    SSR                 = 60        # Sidemode Suppression Ratio: Ratio of optical power in the main mode t optical power
                                    # of the highest sidemode [dB]

    SSSER               = 70        # Signal to Source spontanoeus emission ratio : Ratio of signal po [dB/nm]
                                    #
    c                   = 299792458 # Velocity of light in vacuum [m/s]

    def __init__(self, plot_spectrum=False):
        self.plot_spectrum = plot_spectrum

    def set_output_power(self, output_power):
        """Set the output power for the laser."""
        self.output_power = output_power

    def set_central_wavelength(self, wavelength):
        """Set the central wavelength for the laser."""
        self.central_wavelength = wavelength

    def set_wavelength_range(self, wavelength_start, wavelength_stop, wavelength_step):
        self.wavelength_start = wavelength_start
        self.wavelength_stop = wavelength_stop
        self.wavelength_step = wavelength_step

    def gaussian_linewidth(self):
        """Model the spectral linewidth as a gaussian function with parameters"""

        # Laser line parameters
        mu = self.central_wavelength    # central wavelength (mean)
        sigma = self.FWHM_wavelength()/(2*sqrt(2*log(2))) # waist (standard deviation)

        # Spectral response of the laser
        wavelength = np.arange(self.wavelength_start, self.wavelength_stop, self.wavelength_step)
        P = []
        for i in wavelength:
            P.append(1e-3*10**(self.output_power/10)*exp(-((i-mu)**2)/(2*sigma**2)) + 1e-3*10**((self.output_power-self.SSSER)/10))

        if self.plot_spectrum==True:
            plt.plot(wavelength*1e9,10*np.log10(np.asarray(P)/1e-3))
            plt.xlabel("Wavelength [nm]")
            plt.ylabel("Power [dB]")
            plt.title("Input signal from the laser")
            plt.show()

        return wavelength, P

    def FWHM_wavelength(self):
        """Convert the FWHM from frequency to wavelength."""
        return self.central_wavelength**2/self.c * self.linewidth


    def generate_SSE(self):
        """Add Source Spontaneous Emission (SSE) to the spectrum. """

        pass

if __name__ == "__main__":
    laser = Agilent_TLM(plot_spectrum=True)
    laser.set_central_wavelength(1540e-9)
    laser.set_wavelength_range(1535e-9, 1545e-9, 0.01e-9)
    laser.gaussian_linewidth()