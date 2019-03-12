from math import pi, sqrt

class heater_basic(object):
    """ Basic model for a thermal phase shifter with constant resistance and constant phase efficiency. """

    # TODO : Have the resistance change with applied power
    # TODO : Make a resistance array with crosstalk in order to predict the ohmic crosstalk
    # TODO : Model the phase_efficiency with physical parameters (coverage, temperature spread, etc)

    def __init__(self, phase_efficiency, resistance):
        """ Constructor for the heater object. """
        self.phase_efficiency = phase_efficiency    # Thermo_optic phase efficiency : heater power required to get a pi phase-shift. [W]
        self.resistance = resistance                # Resistance : resistance of the heater [Ohm].

    def get_phase_shift(self, heater_power):
        """ Get the phase-shift [rad] with the applied heater power [W]. """
        return pi/self.phase_efficiency * heater_power

    def get_heater_power(self, voltage):
        """ Get the heater power [W] from the applied voltage [V] assuming ohmic power. """
        return voltage**2/self.resistance

    def apply_voltage(self, voltage):
        """ Apply voltage to the heater [V] and measure the resulting phase shift [rad]. """
        return self.get_phase_shift(self.get_heater_power(voltage))

    def get_FSR_voltage(self):
        """ Get the voltage which allows to sweep one full free spectral range (FSR). """
        return sqrt(2*self.phase_efficiency * self.resistance)


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt

    # Sweep the voltage and measure the heater power
    PS = heater_basic(20e-3, 600)
    V, P = np.linspace(0,10,10), []
    for v in V:
        P.append(PS.get_heater_power(v))
    plt.plot(V, np.asarray(P)*1e3)
    plt.xlabel('Voltage [V]')
    plt.ylabel('Heater Power [mW]')
    plt.grid()
    plt.show()

    # Sweep the voltage and measure the phase shift
    V, phase = np.linspace(0,5,10), []
    for v in V:
        phase.append(PS.apply_voltage(v))
    plt.plot(V, np.asarray(phase)/pi)
    plt.xlabel('Voltage [V]')
    plt.ylabel('Phase [$\pi$ rad]')
    plt.grid()
    plt.show()

    # Return the maximal voltage that can be applied
    print(PS.get_FSR_voltage())

