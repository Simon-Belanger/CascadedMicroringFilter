"""
Model for a Cascaded Microring Filter using the Transfer Matrix Method (TMM)
    - Bidirectionnal (4x4 transfer matrices)

Author      : Simon Bélanger-de Villers (simon.belanger-de-villers.1@ulaval.ca)
Created     : 2018
Last edited : 26-02-2019
"""

# Future modifications:
# TODO : Be able to import different types of waveguides (simple analytical, MODE simulation, etc).
# TODO : Be able to use different types of couplers (simple analytical, FDTD simulations).
# TODO : Each ring should have its own ring object with basic properties such as length round trip, losses, phase ...
# TODO : Clean the code.
# TODO : Publish on github and make documentation for it. Do version history.
# TODO : Compare with EMPy for code structure and functionnalities


import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
from misc.utils import t2s, thermal_crosstalk
from math import pi, sqrt, log10, exp
import cmath as cm
import random
from components.phase_shifter import heater_basic
from components.ring import Ring


# Microring Filter Class, generates a model for a high-order microring filter
class MRF(object):
    """"""
    c = 299792458  # Velocity of light in vaccum [m/s]

    def __init__(self, name='', num_rings=5, radius=2.5e-6, neff=2.4449, alpha_wg=3., couplers=[None, None, None, None, None, None], crosstalk_coeff=[1, 0., 0.]):
        """ Constructor for the cascaded microring filter object. """

        self.name = name

        # Rings
        self.Rings = [Ring(radius, neff, alpha_wg)]*num_rings

        # Couplers
        self.C_Mat = []
        for c in couplers:
            self.C_Mat.append(c.T)

        # Additional phase (zero by default)
        self.bias = np.zeros((num_rings,))
        self.phase_deviation = np.zeros((num_rings,))
        self.desired_tuning_phase = np.zeros((num_rings,))  # Desired phase tuning if there was no thermal crosstalk
        self.actual_tuning_phase = np.zeros((num_rings,))   # Actual phase tuning after thermal crosstalk

        # Phase shifting
        self.phaseshifters = [heater_basic(phase_efficiency=20e-3, resistance=600)] * num_rings
        self.phase_coupling_matrix = thermal_crosstalk(num_rings, crosstalk_coeff[0], crosstalk_coeff[1], crosstalk_coeff[2])

        # For algorithms
        self.NM_phase, self.NM_power = [], []

    def get_total_phase(self, wavelength):
        """
        Get the modulus for the total phase for the cascaded ring.
        """

        inner_product = 0
        for ring_id in range(len(self.Rings)):
            inner_product += (self.Rings[ring_id].get_total_phase(wavelength))**2
        # If the remaining phase is 0 (perfectly tuned) prevents errors
        try:
            return log10(sqrt(inner_product)%(2*pi))
        except ValueError:
            return -100

    def manufacturing(self, var):
        """ Add a random phase perturbation that can be attributed to manufacturing variability. """
        for ring_id in range(len(self.Rings)):
            self.Rings[ring_id].set_phase_deviation(random.uniform(-var, var))


    def TMM(self, wavelength, E_in, E_add):
        """ Get the total transfer matrix for the cascaded microring filter. """

        # Transfer matrix multiplication
        P_list = []
        for ring_id in range(len(self.Rings)):
            P_list.append(self.Rings[ring_id].get_propagation_matrix(wavelength))
        C_list = self.C_Mat
        listmat = [None] * (len(P_list) + len(C_list))
        listmat[0:len(listmat):2] = C_list
        listmat[1:len(listmat) - 1:2] = P_list

        H = listmat[0]
        for el in listmat[1::]:
            H = H * el

        # Convert T matrix to S matrix
        S = t2s(H)
        E_out = S * np.matrix([[E_in], [0], [0], [E_add]])

        return E_out

    def TMM_v2(self, wavelength, E_in, E_add):
        """
        Calculate the response of the MRF system to an input signal over a broad wavelength range.


        Inputs:
            wavelength  : Wavelength grid for the input/add signals [nx1 list]
            E_in        : Electric field at the input port vs wavelength [nx1 list]
            E_add       : Electric field at the add port vs wavelength [nx1 list]

        Outputs:
            E_through   : Electric field at the through port vs wavelength [nx1 list]
            E_drop      : Electric field at the drop port vs wavelength [nx1 list]

        """

        E_thru = np.zeros(len(wavelength), dtype=complex)  # through-port response
        E_drop = np.zeros(len(wavelength), dtype=complex)  # drop-port response

        for ii in range(len(wavelength)):
            E_out = self.TMM(wavelength[ii], E_in[ii], E_add[ii])

            E_thru[ii] = E_out[1]  # through-port field response
            E_drop[ii] = E_out[2]  # drop-port field response

        return E_drop, E_thru

    def apply_bias(self, ring_id, voltage):
        """ Apply voltage to the corresponding phase shifter and update the tuning phase variable accordingly.
            Phase crosstalk is considered through the phase coupling matrix.

            Inputs:
                ring_id : Heater for which voltage is applied
                voltage : Voltage to apply to the i-th heater [V]
        """

        # Measure desired phase shift
        self.desired_tuning_phase[ring_id] = self.phaseshifters[ring_id].apply_voltage(voltage)

        # Obtain actual phase shift by including phase crosstalk
        self.actual_tuning_phase = np.asarray(np.squeeze(self.phase_coupling_matrix * np.transpose(np.asmatrix(self.desired_tuning_phase))))[0]

    def apply_tuning_phase(self, ring_id, phase):
        """ Apply voltage to the corresponding phase shifter and update the tuning phase variable accordingly. """
        self.actual_tuning_phase[ring_id] = phase

    def measure_power(self, lambda_0):
        """Measure the power coming out of the drop port at wavelength lambda_0."""

        # Get the field at the four ports
        E = self.TMM(lambda_0, 1, 0)
        P_drop = float(10 * np.log10(abs(E[2]) ** 2))
        P_thru = float(10 * np.log10(abs(E[1]) ** 2))

        # Effect of photodetector Noise floor
        noise_floor_mean        = -80
        noise_floor_amplitude   = 2
        if  P_drop < noise_floor_mean:  # Noise Floor
            P_drop = noise_floor_mean + random.uniform(-noise_floor_amplitude/2,noise_floor_amplitude/2)
        if  P_thru < noise_floor_mean:  # Noise Floor
            P_thru = noise_floor_mean + random.uniform(-noise_floor_amplitude / 2, noise_floor_amplitude / 2)

        return P_drop, P_thru

    def test_MRF(self, wavelength, bias_array):
        """"""
        # Apply bias
        for ring_id in range(len(self.Rings)):
            self.apply_bias(ring_id, bias_array[ring_id])
        # Store total phase in an array
        self.NM_phase.append(self.get_total_phase(wavelength))
        self.NM_power.append(self.measure_power(wavelength)[0])

        # Return power
        return self.measure_power(wavelength)[0]

    def minimize_MRF(self, bias_array):
        """"""
        return -1 * self.test_MRF(self.target_wavelength, bias_array)

    def sweep(self, wavelength=np.linspace(1530,1560,1000)*1e-9, E_in=1, E_add=0, plot_results=True):
        """"""
        E_thru = np.zeros(len(wavelength), dtype=complex)  # through-port response
        E_drop = np.zeros(len(wavelength), dtype=complex)  # drop-port response

        for ii in range(len(wavelength)):
            E_out = self.TMM(wavelength[ii], E_in, E_add)

            E_thru[ii] = E_out[1]  # through-port field response
            E_drop[ii] = E_out[2]  # drop-port field response

        if plot_results == True:
            self.plot_transmission(wavelength, E_drop, E_thru)

        return E_drop, E_thru

    def phase_response(self, wavelength, E_in, E_add, plot_results=False):
        """"""
        E_thru, E_drop = self.TMM(wavelength, E_in, E_add);
        phi_thru = np.angle(E_thru)
        phi_drop = np.angle(E_drop)

        if plot_results == True:
            self.plot_phase_response(wavelength, phi_drop, phi_thru)

        return phi_drop, phi_thru

    def group_delay(self, wavelength, E_in, E_add, plot_results=False):
        """"""
        omega = self.c / wavelength * 2*pi
        phi_thru, phi_drop = self.phase_response(wavelength, E_in, E_add)
        tau_thru = - np.diff(phi_thru) / np.diff(omega)
        tau_drop = - np.diff(phi_drop) / np.diff(omega)

        if plot_results == True:
            self.plot_group_delay(wavelength, tau_drop, tau_thru)
        return tau_thru, tau_drop

    @staticmethod
    def plot_transmission(wavelength, E_drop, E_thru):
        """ Plot the drop/thru port transmission spectrum. """
        plt.plot(wavelength * 1e9, 10 * np.log10(abs(E_drop) ** 2), label='Drop')
        plt.plot(wavelength * 1e9, 10 * np.log10(abs(E_thru) ** 2), label='Thru')
        plt.xlabel('wavelength (nm)', Fontsize=14)
        plt.ylabel('Power transmission (dB)', Fontsize=14)
        plt.legend(loc='upper right')
        plt.show()

    @staticmethod
    def plot_phase_response(wavelength, phi_drop, phi_thru):
        """ Plot the drop/thru port phase response. """
        plt.plot(wavelength * 1e9, phi_drop / np.pi, label='Drop')
        plt.plot(wavelength * 1e9, phi_thru / np.pi, label='Thru')
        plt.xlabel('wavelength (nm)', Fontsize=14)
        plt.ylabel('Phase(\pi)', Fontsize=14)
        plt.legend()
        plt.show()

    @staticmethod
    def plot_group_delay(wavelength, tau_drop, tau_thru):
        """ Plot the drop/thru port group delay. """
        plt.plot(lambda_0[0:len(wavelength) - 1] * 1e9, tau_drop * 1e12, label='Drop')
        plt.plot(lambda_0[0:len(wavelength) - 1] * 1e9, tau_thru * 1e12, label='Thru')
        plt.xlabel('wavelength (nm)', Fontsize=14)
        plt.ylabel('Group delay (ps)', Fontsize=14)
        plt.legend()
        plt.show()

if __name__ == '__main__':

    # Build the couplers
    from components.couplers.virtual_DC import virtual_DC
    from components.couplers.apodization import flattop5

    k_vec = flattop5(0.1)
    loss_coupler = [0.98,0.99,0.99,0.99,0.99,0.99]
    couplers = [virtual_DC(k_vec[0],loss_coupler[0]), virtual_DC(k_vec[1],loss_coupler[1]), virtual_DC(k_vec[2],loss_coupler[2]), virtual_DC(k_vec[3],loss_coupler[3]), virtual_DC(k_vec[4],loss_coupler[4]), virtual_DC(k_vec[5],loss_coupler[5])]

    # Build the Microring filter
    mrf = MRF( name='', num_rings=5, R=2.5e-6, couplers=couplers, alpha_wg=3., crosstalk_coeff=[1, 0.75, 0.1])

    # Sweep the device
    #E_drop, E_thru = mrf.sweep(lambda_0=np.linspace(1530,1560,1000)*1e-9, E_in=1, E_add=0, plot_results=True)
    #P_drop = 10 * np.log10(abs(E_drop))
    #P_thru = 10 * np.log10(abs(E_thru))
    # For data analysis (to remove soon)
    #np.savetxt("/Users/simonbelanger/Documents/UL/Silicon_Photonics/01_Projets/PaperLukas_dec2018/02_experimental_data/02_analyzed/simulation_data.txt",(np.linspace(1530,1560,1000)*1e-9,P_thru,P_drop))

    # Sweep the device with real laser
    from components.optical_sources.Agilent_TLM import Agilent_TLM
    laser = Agilent_TLM()
    laser.set_central_wavelength(1550e-9)
    laser.set_wavelength_range(1500e-9,1600e-9,0.1e-9)
    wavelength, P   = laser.gaussian_linewidth()
    E_drop, E_thru = mrf.TMM_v2(wavelength, np.sqrt(P), np.zeros(len(wavelength)))

    # Measure the power using the PD
    from components.detectors.Agilent_PD import Agilent_PD
    PD = Agilent_PD()
    PD.plot_transmission(wavelength, E_drop, E_thru)



