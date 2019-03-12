import cmath as cm
from math import pi,exp
import numpy as np

class Ring(object):
    """ Ring Class, generates a model for a single ring and it's parameters. """

    def __init__(self, ring_radius, effective_index, alpha_wg):
        " Constructor for the ring object. "

        self.ring_radius        = ring_radius       # Radius of the ring [m]
        self.effective_index    = effective_index   # Effective index of the waveguide [-]
        self.alpha_wg           = alpha_wg          # Propagation losses [dB/cm]

        # Phase
        self.phase_deviation    = 0                 # Phase offset due to manufacturing variability [rad]
        self.tuning_phase       = 0                 # Phase added through a tuning mechanism [rad]

    def get_roundtrip_length(self):
        " Getter for the roundtrip length attribute of the resonator [m]. "
        return 2*pi * self.ring_radius

    def get_loss_coefficient(self):
        " Getter for the loss coefficient attribute of the resonator [m-1]. "
        return 0.23 / 2 * 100 * self.alpha_wg

    def get_complex_phase(self, wavelength):
        " Getter for the complex phase attribute. "
        return (2*pi * self.effective_index / wavelength - 1j * self.get_loss_coefficient()) * self.get_roundtrip_length()

    def get_roundtrip_phase(self, wavelength):
        """ Setter for the roundtrip phase attribute.
        An input wavelength that is on resonance will yield a 0 (mod 2pi) roundtrip phase.
        """
        return (self.get_complex_phase(wavelength).real) % (2*pi)

    def set_phase_deviation(self, phase_deviation):
        " Setter for the phase deviation attribute. "
        self.phase_deviation = phase_deviation % (2*pi)

    def get_phase_deviation(self):
        " Getter for the phase deviation attribute. "
        return self.phase_deviation % (2*pi)

    def set_tuning_phase(self, tuning_phase):
        " Setter for the tuning phase attribute. "
        self.tuning_phase = tuning_phase % (2*pi)

    def get_tuning_phase(self):
        " Getter for the tuning phase attribute. "
        return self.tuning_phase % (2*pi)

    def get_total_phase(self, wavelength):
        " Getter for the total phase attribute. "
        return (self.get_roundtrip_phase(wavelength) + self.get_phase_deviation() + self.get_tuning_phase()) % (2*pi)

    def get_propagation_matrix(self, wavelength):
        " Getter for the propagation matrix attribute. "

        # Total phase = intrinsic phase + phase deviation + tuning phase
        phi = self.get_complex_phase(wavelength) + self.get_phase_deviation() + self.get_tuning_phase()

        return np.matrix([[cm.exp(-1j * phi/2),    0,                  0,                          0                   ],
                          [0,                   cm.exp(-1j * phi/2),   0,                          0                   ],
                          [0,                   0,                  1/cm.exp(-1j * phi/2),         0                   ],
                          [0,                   0,                  0,                          1/cm.exp(-1j * phi/2)]])