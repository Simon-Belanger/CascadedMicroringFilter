import numpy as np
import matplotlib.pyplot as plt


class Microring_AD(object):
    """ Microring resonator model
        by Wei Shi -- wei.j.shi@gmail.com
        created in 2012 and modified in May 2013
        Class defined in May 2017
        conversion to python by Simon Bélanger-de Villers in August 2018
    """

    c = 299792458  # speed of light m/s

    def __init__(self, R=10e-6, k1_amp=0.1, k2_amp=0.1, Tc1=1., Tc2=1., alpha_wg=3.):

        self.R = R                          # Ring radius
        self.k1_amp = k1_amp                # coupler 1 (bus WG): amplitude of the cross-coupling coefficient
        self.k2_amp = k2_amp                # coupler 2 (add-drop WG): amplitude of the cross-coupling coefficient
        self.Tc1 = Tc1                      # coupling loss: t*t'-k*k'=1 in the lossless case
        self.Tc2 = Tc2                      # coupling loss: t*t'-k*k'=1 in the lossless case
        self.alpha_wg = alpha_wg            # Propagation losses dB/cm

        self.k1 = self.get_k(k1_amp)        # coupler 1 (bus WG): cross-coupling coefficient
        self.k2 = self.get_k(k2_amp)        # coupler 2 (add-drop WG): cross-coupling coefficient
        self.t1 = self.get_t(k1_amp, Tc1)   # coupler 1: straight-through coupling coefficient
        self.t2 = self.get_t(k2_amp, Tc2)   # coupler 2: straight-through coupling coefficient
        self.get_L_rt()                     # round-trip length
        self.get_a()                        # round-trip field loss

    @staticmethod
    def get_k(k_amp):
        return k_amp * 1j

    @staticmethod
    def get_t(k_amp, Tc):
        return np.sqrt(Tc - k_amp ** 2)

    def get_L_rt(self):
        self.L_rt = 2 * np.pi * self.R;

    def get_a(self):
        self.a = np.sqrt(10 ** (-self.alpha_wg * self.L_rt * 100 / 10))

    def get_beta(self, lambda_0):

        # MODE data
        neff_data = np.array([2.4992, 2.4725, 2.4449, 2.4163, 2.3868])  # calculated using MODE
        lambda_data = np.array([1500, 1525, 1550, 1575, 1600]) * 1e-9  # wavelength points in the MODE simulation

        neff = np.interp(lambda_0, lambda_data, neff_data)
        beta = neff * 2 * np.pi / lambda_0
        return beta

    def analytic(self, lambda_0):
        A = self.a ** 2
        beta = self.get_beta(lambda_0)
        Thru_analytic = np.zeros(len(lambda_0))
        for i in range(len(lambda_0)):
            Thru_analytic[i] = (self.t1 ** 2 + A * self.t2 ** 2 - 2 * self.t1 * self.t2 * self.a * np.cos(
                beta[i] * self.L_rt)) / (
                                       1 + A * self.t1 ** 2 * self.t2 ** 2 - 2 * self.a * self.t1 * self.t2 * np.cos(
                                   beta[i] * self.L_rt))  # add-drop filter
        return Thru_analytic

    def TMM(self, lambda_0, E_in, E_add, plot_results=False):
        """"""
        C1 = 1 / self.k1 * np.matrix([[-self.t1, 1], [-self.Tc1, self.t1]])  # coupler matrix
        C2 = 1 / self.k2 * np.matrix([[-self.t2, 1], [-self.Tc2, self.t2]])

        beta = self.get_beta(lambda_0)  # propagation constant

        E_thru = np.zeros(len(lambda_0), dtype=complex)  # through-port response
        E_drop = np.zeros(len(lambda_0), dtype=complex)  # drop-port response

        for ii in range(len(lambda_0)):
            P = np.matrix([[0, np.sqrt(self.a) * np.exp(-1j * beta[ii] * self.L_rt / 2)],
                           [1 / (np.sqrt(self.a) * np.exp(-1j * beta[ii] * self.L_rt / 2)), 0]])  # ring cavity matrix
            H = C1 * P * C2  # total transfer matrix
            E_drop[ii] = (E_in - H[0, 0] * E_add) / H[0, 1]  # through-port field response
            E_thru[ii] = H[1, 0] * E_add + H[1, 1] * E_drop[ii]  # drop-port field response

        if plot_results == True:
            self.plot_TMM(lambda_0, E_drop, E_thru)

        return E_drop, E_thru

    def plot_TMM(self, lambda_0, E_drop, E_thru):

        plt.plot(lambda_0 * 1e9, 10 * np.log10(abs(E_drop) ** 2), label='Drop')
        plt.plot(lambda_0 * 1e9, 10 * np.log10(abs(E_thru) ** 2), label='Thru')
        plt.xlabel('wavelength (nm)', Fontsize=14)
        plt.ylabel('Power transmission (dB)', Fontsize=14)
        plt.legend(loc='upper right')
        plt.show()

    def phase_response(self, lambda_0, E_in, E_add, plot_results=False):
        """"""
        E_thru, E_drop = self.TMM(lambda_0, E_in, E_add);
        phi_thru = np.angle(E_thru);
        phi_drop = np.angle(E_drop);

        if plot_results == True:
            self.plot_phase_response(lambda_0, phi_drop, phi_thru)

        return phi_drop, phi_thru

    def plot_phase_response(self, lambda_0, phi_drop, phi_thru):

        plt.plot(lambda_0 * 1e9, phi_drop / np.pi, label='Drop')
        plt.plot(lambda_0 * 1e9, phi_thru / np.pi, label='Thru')
        plt.xlabel('wavelength (nm)', Fontsize=14)
        plt.ylabel('Phase(\pi)', Fontsize=14)
        plt.legend()
        plt.show()

    def group_delay(self, lambda_0, E_in, E_add, plot_results=False):
        omega = self.c / lambda_0 * 2 * np.pi;
        phi_thru, phi_drop = self.phase_response(lambda_0, E_in, E_add)
        tau_thru = - np.diff(phi_thru) / np.diff(omega);
        tau_drop = - np.diff(phi_drop) / np.diff(omega);

        if plot_results == True:
            self.plot_group_delay(lambda_0, tau_drop, tau_thru)
        return tau_thru, tau_drop

    def plot_group_delay(self, lambda_0, tau_drop, tau_thru):

        plt.plot(lambda_0[0:len(lambda_0) - 1] * 1e9, tau_drop * 1e12, label='Drop')
        plt.plot(lambda_0[0:len(lambda_0) - 1] * 1e9, tau_thru * 1e12, label='Thru')
        plt.xlabel('wavelength (nm)', Fontsize=14)
        plt.ylabel('Group delay (ps)', Fontsize=14)
        plt.legend()
        plt.show()

if __name__ == "__main__":

    # Wavelength range
    lambda_0 = np.linspace(1555, 1565, 1000) * 1e-9

    # Instantiate Microring object
    mrr = Microring_AD()

    # Get the power transmission for the drop port and for the through port over the wavelength range.
    E_drop, E_thru = mrr.TMM(lambda_0, 1, 0, plot_results=True)

    # Get the phase response for the drop port and for the through port over the wavelength range.
    phi_drop, phi_thru = mrr.phase_response(lambda_0, 1, 0, plot_results=True)

    # Get the group delay for the drop port and for the through port over the wavelength range.
    tau_thru, tau_drop = mrr.group_delay(lambda_0, 1, 0, plot_results=True)