""" Microring Resonator analytical model from W. Bogaerts paper 'Silicon microring resonators' ."""

from math import cos,sin,atan,pi,exp,sqrt
import numpy as np
import cmath
import matplotlib.pyplot as plt

class MRR(object):
    """ General Microring Resonator."""

    ref_wavelength = 1550e-9    # Reference wavelength for the computation of the propagation constant [m]
                                # For the taylor expansion
    c = 299792458               # Velocity of light in vaccum [m/s]

    def __init__(self, radius, loss_per_cm, n_eff, n_g):
        """ Constructor for the ring resonator class. """

        self.R = radius         # Ring radius [m]
        self.loss = loss_per_cm # Waveguide losses [dB/cm]
        self.n_eff = n_eff      # Effective index of the waveguide []
        self.n_g = n_g          # Group index of the waveguide []

        self.get_roundtrip_length()
        self.get_alpha()
        self.get_roundtrip_loss()

    def get_roundtrip_length(self):
        """ Obtain the roundtrip length from the radius of the ring resonator. """
        self.L_rt = 2 * pi * self.R  # Round-trip length [m]

    def get_alpha(self):
        """ Obtain the loss coefficient from the waveguide losses in dB/cm. """
        self.alpha = 23 * self.loss  # Linear loss coefficient [m-1]

    def get_roundtrip_loss(self):
        """ Obtain the roundtrip losses from the loss coefficient and the roundtrip length of the ring resonator. """
        self.a = sqrt(exp(- self.alpha * self.L_rt))

    def get_phase(self, lambda_0):
        """ Obtain the roundtrip phase at a particular wavelength from the effective index and the roundtrip length of
        the ring resonator. """
        beta = (2 * pi * self.n_eff) / self.ref_wavelength + (2 * pi * self.n_g) * (1/lambda_0 - 1/self.ref_wavelength) # Propagation constant[m-1]
        phi = beta * self.L_rt                      # Optical phase [rad]
        return phi

    def get_field_transmission(self, lambda_0, E_in):
        """ This method will be defined in the inherited classes 'MRR_AP' and 'MRR_AD'. """
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
        plt.xlabel('Wavelength (nm)', Fontsize=14)
        plt.ylabel('Power transmission (dB)', Fontsize=14)
        plt.show()

    def plot_phase(self, lambda_0, E_t):
        """ Plot the drop/thru port transmission spectrum. """
        plt.plot(lambda_0 * 1e9, np.unwrap(np.angle(E_t))/ np.pi)
        plt.xlabel('Wavelength (nm)', Fontsize=14)
        plt.ylabel('Phase (/pi)', Fontsize=14)
        plt.show()

    def measure_power_thru(self, E_t):
        return 10 * np.log10(abs(E_t) ** 2)

    def measure_phase_thru(self, E_t):
        return np.angle(E_t)/ np.pi

    @staticmethod
    def get_r(C):
        """ Obtain the field thru-coupling coefficient 'r' from the power cross-coupling coefficient 'C'. """
        return sqrt(1 - C)



class MRR_AP(MRR):
    """ All-Pass Microring Resonator"""
    def __init__(self, radius, loss_per_cm, n_eff, n_g, C):
        """ Constructor for the all-pass ring resonator class. """
        super().__init__(radius, loss_per_cm, n_eff, n_g)

        self.C = C                  # Power cross-coupling coefficient for the coupler [W/W]
        self.r = self.get_r(self.C) # Field thru-coupling coefficient for the coupler [W/W]

    def get_field_transmission(self, lambda_0, E_in):
        phi = self.get_phase(lambda_0)
        return E_in * cmath.exp(1j * (cmath.pi + phi)) * (self.a - self.r * cmath.exp(-1j*phi))/(1 - self.a * self.r * cmath.exp(1j*phi))

    def get_effective_phase_shift(self, lambda_0):
        phi = self.get_phase(lambda_0)
        return pi + phi + atan((self.r * sin(phi))/(self.a - self.r * cos(phi))) + atan((self.r * self.a * sin(phi))/(1 - self.r * self.a * cos(phi)))

    def get_power_transmission(self, lambda_0, I_in):
        """ Obtain the power transmission to the thru port. """
        phi = self.get_phase(lambda_0)
        return (self.a ** 2 - 2 * self.a * self.r * cos(phi) + self.r ** 2) / (1 - 2 * self.a * self.r * cos(phi) + (self.r * self.a) ** 2)

    def critical_coupling_condition(self):
        """ Find the power cross-coupling coefficient which satisfies the critical coupling condition. """

        print("The coupling coefficient required to achieve critical coupling condition is {}".format(1 - self.a**2))


class MRR_AD(MRR):
    """Add-Drop Microring Resonator"""

    def __init__(self, radius, loss_per_cm, n_eff, n_g, C1, C2):
        """ Constructor for the add-drop ring resonator class. """
        super().__init__(radius, loss_per_cm, n_eff, n_g)

        self.C1 = C1    # Power cross-coupling coefficient for the input coupler [W/W]
        self.C2 = C2    # Power cross-coupling coefficient for the output coupler [W/W]

        self.r1 = self.get_r(self.C1) # Field thru-coupling coefficient for the input coupler [W/W]
        self.r2 = self.get_r(self.C2) # Field thru-coupling coefficient for the output coupler [W/W]

    def get_power_transmission(self, wavelength):
        """ Obtain the power transmission to the drop and thru port. """

        phi = self.get_phase(wavelength)

        denom = 1 - 2 * self.r1 * self.r2 * self.a * cos(phi) + (self.r1 * self.r2 * self.a)**2             # Common denominator
        T_p = (self.r2**2 * self.a**2 - 2 * self.r1 * self.r2 * self.a * cos(phi) + self.r1**2) / denom     # Thru
        T_d = ((1 - self.r1**2) * (1 - self.r2**2) * self.a) / denom                                        # Drop

        return T_p, T_d

    def plot_transmission_spectrum(self, wavelength_list):
        """ Plot the power transmission spectrum for the drop and thru port. """

        # Wavelength sweep
        T_d, T_t = [], []
        for wavelength in wavelength_list:
            T, D = self.get_power_transmission(wavelength)
            T_t.append(T); T_d.append(D)

        # Plot the figure
        plt.plot(wavelength_list * 1e9, 10 * np.log10(T_d), label='Drop')
        plt.plot(wavelength_list * 1e9, 10 * np.log10(T_t), label='Thru')
        plt.xlabel('Wavelength (nm)', Fontsize=14)
        plt.ylabel('Power transmission (dB)', Fontsize=14)
        plt.legend()
        plt.show()

    def critical_coupling_condition(self):
        """ Find the power cross-coupling coefficient for the input coupler which satisfies the critical coupling condition. """

        print("The coupling coefficient required to achieve critical coupling condition is {}".format(1 - (self.r2 * self.a)**2))

if __name__ == "__main__":

    #AP1 = MRR_AP(20e-6, 2, 2.5, 4.5, 0.5)
    #wvl, field = AP1.sweep_field_transmission(1540e-9, 1560e-9, 1000)
    #AP1.plot_field_transmission(wvl, field)
    #AP1.plot_phase(wvl, field)

    #AP1.critical_coupling_condition()

    #phase = []
    #for wvl in np.linspace(1540e-9, 1560e-9, 1000):
    #    phase.append(AP1.get_effective_phase_shift(wvl))

    #plt.plot(np.linspace(1540e-9, 1560e-9, 1000),phase)
    #plt.show()

    AD1 = MRR_AD(radius=2.5e-6, loss_per_cm=2, n_eff=2.5, n_g=4.5, C1=0.20057784425772285, C2=0.2)

    AD1.plot_transmission_spectrum(np.linspace(1500e-9,1600e-9,1000))

    AD1.critical_coupling_condition()


