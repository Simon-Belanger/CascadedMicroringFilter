""" All-Pass Microring Resonator analytical model from W. Bogaerts paper 'Silicon microring resonators' ."""

from analytical.MRR import *

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
        """ Get the effective phase shift induced by the resonator ."""
        phi = self.get_phase(lambda_0)
        return pi + phi + atan((self.r * sin(phi))/(self.a - self.r * cos(phi))) + atan((self.r * self.a * sin(phi))/(1 - self.r * self.a * cos(phi)))

    def get_power_transmission(self, lambda_0, I_in):
        """ Obtain the power transmission to the thru port. """
        phi = self.get_phase(lambda_0)
        return (self.a ** 2 - 2 * self.a * self.r * cos(phi) + self.r ** 2) / (1 - 2 * self.a * self.r * cos(phi) + (self.r * self.a) ** 2)

    def sweep_power_transmission(self, lambda_min, lambda_max, lambda_points):
        """ Obtain the power transmission spectrum of the ring resonator. Implemented in subclasses. """
        lambda_0 = np.linspace(lambda_min, lambda_max, lambda_points)
        E_t = []
        for wvl in lambda_0:
            E_t.append(self.measure_power_thru(self.get_power_transmission(wvl, 1)))
        return lambda_0, np.asarray(E_t)

    def critical_coupling_condition(self):
        """ Find the power cross-coupling coefficient which satisfies the critical coupling condition. """

        print("The coupling coefficient required to achieve critical coupling condition is {}".format(1 - self.a**2))

    def getQfactor(self, lambdaRes=1550e-9):
        """ Return the Q factor for the cavity for a given resonance wavelength. """
        return (math.pi * self.n_g * self.L_rt * math.sqrt(self.r * self.a))/(lambdaRes * (1 - self.r * self.a))