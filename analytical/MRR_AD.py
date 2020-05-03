""" Add-Drop Microring Resonator analytical model from W. Bogaerts paper 'Silicon microring resonators' ."""

from analytical.MRR import MRR

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