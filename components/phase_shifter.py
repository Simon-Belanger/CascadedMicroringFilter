"""
Thermo-Optic Phase Shifters (TOPS) models.

Author: Simon Belanger-de Villers
Created : 2019
Last Edited : March 14th 2021
"""
from math import pi, sqrt

# TODO : Have the resistance change with applied power (nonlinear)
# TODO : Make a resistance array with crosstalk in order to predict the ohmic crosstalk
# TODO : Model the phase_efficiency with physical parameters (coverage, temperature spread, etc)
# TODO : Specify different types of phase efficiencies (rad/mW, for 2pi phaseshift etc)
# TODO : Make use of MaxPower 

class heaterBasic(object):
    """ Basic model for a thermal phase shifter with constant resistance and constant phase efficiency. 

        Arguments:
            phaseEfficiency : Input power which results in a pi phase shift (P_pi)      [W]
            resistance      : Electrical resistance of the heater.                      [Ohms]
            maxPower        : Input power which causes breakage of the device.          [W]
    """

    def __init__(self, phaseEfficiency, resistance):
        """ Constructor for the heater object. """
        self.phaseEfficiency    = phaseEfficiency
        self.resistance         = resistance

    def getPhaseShift(self, inputPower):
        """ Get the phase-shift [rad] with the applied heater power [W]. """
        return pi * inputPower/self.phaseEfficiency

    def getHeaterPower(self, inputVoltage):
        """ Get the heater power [W] from the applied voltage [V] assuming ohmic power. """
        return inputVoltage**2/self.resistance

    def applyVoltage(self, inputVoltage):
        """ Apply voltage to the heater [V] and measure the resulting phase shift [rad]. """
        return self.getPhaseShift(self.getHeaterPower(inputVoltage))

    def getFsrVoltage(self):
        """ Get the voltage which allows to sweep one full free spectral range (FSR). """
        return sqrt(2*self.phaseEfficiency * self.resistance)

if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    # Sweep the voltage and measure the heater power
    PS = heaterBasic(phaseEfficiency=20e-3, resistance=100)
    V, P = np.linspace(0, 10, 10), []
    for v in V:
        P.append(PS.getHeaterPower(v))
    plt.plot(V, np.asarray(P)*1e3)
    plt.xlabel('Input voltage, $V$ [V]')
    plt.ylabel('Input power, $P$ [mW]')
    plt.grid()
    plt.show()

    # Sweep the voltage and measure the phase shift
    V, phase = np.linspace(0,5,10), []
    for v in V:
        phase.append(PS.applyVoltage(v))
    plt.plot(V, np.asarray(phase)/pi)
    plt.xlabel('Input voltage, $V$ [V]')
    plt.ylabel('Phase shift, $\Delta\phi$ [$\pi$ rad]')
    plt.grid()
    plt.show()

    # Return the maximal voltage that can be applied
    print(PS.getFsrVoltage())

