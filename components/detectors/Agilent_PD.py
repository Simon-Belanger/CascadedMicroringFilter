import numpy as np
import matplotlib.pyplot as plt
import random

class Agilent_PD(object):

    name = 'Agilent 81635A Dual Power Sensor Module'

    # wavelength_range  = 800-1650 nm
    # power_range       = +10 to –80 dBm
    # Noise (peak-to-peak) = < 20 pW
    # Averaging time (minimal) = 100µs

    # FEATURES
    # Thermal noise and Shot Noise to understand and implement
    # Averaging time more tricky to implement


    # Parameters
    saturation              = 10        # Max power measurable by the PD [dBm]
    noise_floor_mean        = -80       # Mean value of the noise floor [dBm]
    noise_floor_amplitude   = 100e-10    # Amplitude of the noise peak 2 peak

    def __init__(self, name_channel1='Detector 1', name_channel2='Detector 2'):
        self.name_channel1 = name_channel1
        self.name_channel2 = name_channel2

    def measure_power(self, Input_power):
        """ Measure the power and return the digital value. """

        Output_power = []
        for ii in Input_power:
            if ii<self.noise_floor_mean:    # Noise Floor
                Output_power.append(self.gen_noise_floor())
            elif ii>self.saturation:          # Saturation
                Output_power.append(self.saturation)
            else:                           # Assuming perfect responsivity (no offset nor attenuation)
                Output_power.append(ii)
        return Output_power

    def gen_noise_floor(self):
        """
        Generate a random noise floor value in dBm with probabilistic parameters such as:

        mean        : noise floor mean value given by noise_floor_mean [dB]
        Noise_pp    : noise floor peak 2 peak value given by noise_floor_amplitude [dB]

        """
        return self.mW_2_dB(random.uniform(self.dB_2_mW(self.noise_floor_mean) - self.noise_floor_amplitude / 2,
                                    self.dB_2_mW(self.noise_floor_mean) + self.noise_floor_amplitude / 2))

    def plot_transmission(self, wavelength, optical_in1, optical_in2):
        """
        Plot the transmission spectrum for the 2 channels.

        Inputs:
            wavelength      : Wavelength grid for the input/add signals [nx1 list]
            optical_in_1    : Electric field (complex) at optical input #1 vs wavelength [nx1 list]
            optical_in_2    : Electric field (complex) at optical input #2 vs wavelength [nx1 list]
        """

        # Convert complex field to power for both inputs
        P_1 = self.field_2_power(optical_in1)
        P_2 = self.field_2_power(optical_in2)

        # PD filtering effect (saturation, noise, etc)
        P_1_detected = self.measure_power(P_1)
        P_2_detected = self.measure_power(P_2)

        # Plot
        plt.plot(wavelength * 1e9, P_1_detected, label=self.name_channel1)
        plt.plot(wavelength * 1e9, P_2_detected, label=self.name_channel2)
        plt.xlabel('wavelength (nm)', Fontsize=14)
        plt.ylabel('Power transmission (dB)', Fontsize=14)
        plt.legend(loc='upper right')
        plt.axis([np.min(wavelength)*1e9, np.max(wavelength)*1e9, -90., 0.])
        plt.show()

    @staticmethod
    def field_2_power(E):
        """ Convert the electric field amplitude E to optical power. """
        return 10 * np.log10(abs(E) ** 2)

    @staticmethod
    def mW_2_dB(power_mW):
        """ Convert the power from mW to dB. """
        return 10 * np.log10(power_mW)

    @staticmethod
    def dB_2_mW(power_dB):
        """ Convert the power from mW to dB. """
        return 10 ** (power_dB/10)

    def measure_total_power_v2(self, wavelength, optical_in1, optical_in2):
        """ Get the total power measured by the PD with the spectral density data."""

        # Convert complex field to power for both inputs
        P_1 = self.field_2_power(optical_in1)
        P_2 = self.field_2_power(optical_in2)

        # PD filtering effect (saturation, noise, etc)
        P_1_detected = self.measure_power(P_1)
        P_2_detected = self.measure_power(P_2)

        #Linear power density
        P_1_detected_lin = np.power(10, np.asarray(P_1_detected)/10)
        P_2_detected_lin = np.power(10, np.asarray(P_2_detected)/10)

        #Total power
        Total_P1 = 0
        Total_P2 = 0

        delta_wavelength = wavelength[1]-wavelength[0]
        for i in range(len(wavelength)):
            Total_P1 += P_1_detected_lin[i] * delta_wavelength
            Total_P2 += P_2_detected_lin[i] * delta_wavelength

        return self.mW_2_dB(2*Total_P1), self.mW_2_dB(2*Total_P2)

if __name__ == '__main__':
    pass