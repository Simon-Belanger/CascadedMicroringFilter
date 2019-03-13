import cmath as cm
from math import pi,exp
import numpy as np
from misc.utils import s2t

class Ring(object):
    """ Ring Class, generates a model for a single ring and it's parameters. """

    def __init__(self, ring_radius, effective_index, alpha_wg):
        " Constructor for the ring object. "

        self.ring_radius        = ring_radius       # Radius of the ring [m]
        self.effective_index    = effective_index   # Effective index of the waveguide [-]
        self.group_index        = 4.18              # Group index of the waveguide [-]
        self.ref_wavelength     = 1550e-9           # Reference wavelength for taylor expansion [m]
        self.alpha_wg           = alpha_wg          # Propagation losses [dB/cm]

        # Phase
        self.phase_deviation    = 0.                 # Phase offset due to manufacturing variability [rad]
        self.tuning_phase       = 0.                 # Phase added through a tuning mechanism [rad]

    def get_roundtrip_length(self):
        " Getter for the roundtrip length attribute of the resonator [m]. "
        return 2*pi * self.ring_radius

    def get_loss_coefficient(self):
        " Getter for the loss coefficient attribute of the resonator [m-1]. "
        return 0.23 / 2 * 100 * self.alpha_wg

    def get_propagation_constant_old(self, wavelength):
        " Getter for the complex phase attribute. "
        return 2*pi * self.effective_index / wavelength

    def get_propagation_constant(self, wavelength):
        " Getter for the propagation constant, dispersive. "
        return 2*pi*(self.effective_index/self.ref_wavelength + self.group_index*(1/wavelength - 1/self.ref_wavelength))

    def get_roundtrip_phase(self, wavelength):
        """ Getter for the complex phase attribute.
                An input wavelength that is on resonance will yield a 0 (mod 2pi) roundtrip phase."""
        return (self.get_propagation_constant(wavelength) * self.get_roundtrip_length()) % (2*pi)

    def get_imaginary_phase(self):
        " Getter for the complex phase attribute. "
        return self.get_loss_coefficient() * self.get_roundtrip_length()

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

    def get_transfer_matrix(self, wavelength):
        " Getter for the propagation matrix attribute. "

        # TODO : Have a single waveguide (2x2) and create a map ports function in order to build a 4x4

        # Total phase = intrinsic phase + phase deviation + tuning phase
        phi = self.get_total_phase(wavelength) - 1j * self.get_imaginary_phase()

        return s2t(np.matrix([[0,0,cm.exp(-1j * phi/2),0],
                              [0,0,0,cm.exp(-1j * phi/2)],
                              [cm.exp(-1j * phi/2),0,0,0],
                              [0,cm.exp(-1j * phi/2),0,0]]))