import cmath as cm
from math import pi,exp
import numpy as np
import sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc import constants
from misc.utils import s2t

class Ring(object):
    """ Ring Class, generates a model for a single ring and it's parameters. """

    def __init__(self, ringRadius, effectiveIndex, groupIndex, alpha_wg):
        " Constructor for the ring object. "

        self.ringRadius         = ringRadius        # Radius of the ring [m]
        self.effectiveIndex     = effectiveIndex    # Effective index of the waveguide [-]
        self.groupIndex         = groupIndex        # Group index of the waveguide [-]
        self.dispersion         = 0                 # Dispersion coefficient [ps nm-1 km-1]
        self.refWavelength      = 1550e-9           # Reference wavelength for taylor expansion [m]
        self.alpha_wg           = alpha_wg          # Propagation losses [dB/cm]

        # Phase
        self._phaseDeviation    = 0.                # Phase offset due to manufacturing variability [rad]
        self._tuningPhase       = 0.                # Phase added through a tuning mechanism [rad]

    # Properties / Attributes
    @property
    def roundtripLength(self):
        " Getter for the roundtrip length attribute of the resonator [m]. "
        return 2 * pi * self.ringRadius

    @property
    def lossCoefficient(self):
        " Getter for the loss coefficient attribute of the resonator [m-1]. "
        return 0.23 / 2 * 100 * self.alpha_wg

    @property
    def powerLossPerRoundtrip(self):
        " Get the fraction of power lost per roundtrip because of the propagation losses. a^2"
        return 1 - exp(-self.lossCoefficient * self.roundtripLength)**2

    @property
    def imaginaryPhase(self):
        " Getter for the complex phase attribute. "
        return self.lossCoefficient * self.roundtripLength

    @property
    def phaseDeviation(self):
        " Getter for the phase deviation attribute which is the phase attributed to randomness. "
        return self._phaseDeviation % (2*pi)
    @phaseDeviation.setter
    def phaseDeviation(self, phaseDeviation):
        " Setter for the phase deviation attribute which is the phase attributed to randomness. "
        self._phaseDeviation = phaseDeviation

    @property
    def tuningPhase(self):
        " Getter for the tuning phase attribute. "
        return self._tuningPhase % (2*pi)
    @tuningPhase.setter
    def tuningPhase(self, tuningPhase):
        " Setter for the tuning phase attribute. "
        self._tuningPhase = tuningPhase

    @property
    def freeSpectralRange(self):
        " Getter for the theoretical free specral range of the ring object. "
        return self.refWavelength**2 / ( self.groupIndex * self.roundtripLength)

    # Methods
    def getPropagationConstant(self, wavelength):
        " Getter for the propagation constant, dispersive. "
        dOmega = ((2*pi*constants.c)/float(wavelength) - (2*pi*constants.c)/self.refWavelength)
        firstOrder  = (2*pi*constants.c)/self.refWavelength * (float(self.effectiveIndex)/constants.c)
        secondOrder = float(self.groupIndex)/constants.c
        thirdOrder  = 0.5 * (- float(self.dispersion) * 1e-6 * self.refWavelength**2 / 2 * pi * constants.c )
        return firstOrder + secondOrder * dOmega + thirdOrder * dOmega**2

    def getRoundtripPhase(self, wavelength):
        """ Getter for the complex phase attribute.
                An input wavelength that is on resonance will yield a 0 (mod 2pi) roundtrip phase."""
        return (self.getPropagationConstant(wavelength) * self.roundtripLength) % (2*pi)

    def getTotalPhase(self, wavelength):
        " Getter for the total phase attribute. "
        return (self.getRoundtripPhase(wavelength) + self.phaseDeviation + self.tuningPhase) % (2*pi)

    def get_transfer_matrix(self, wavelength):
        " Getter for the propagation matrix attribute. "

        # Total phase = intrinsic phase + phase deviation + tuning phase
        phi = self.getTotalPhase(wavelength) - 1j * self.imaginaryPhase

        return s2t(np.matrix([[0, 0, cm.exp(-1j * phi/2), 0],
                              [0, 0, 0, cm.exp(-1j * phi/2)],
                              [cm.exp(-1j * phi/2), 0, 0, 0],
                              [0, cm.exp(-1j * phi/2), 0, 0]]))

    def resonanceWavelengths(self):
        " Find the resonant wavelengths of the ring. "
        resonantWavelengths = []
        for m in range(0, 10):
            resonantWavelengths.append(2*pi*self.ringRadius*self.effectiveIndex/m)
        return resonantWavelengths


if __name__ == '__main__':

    r = Ring(2.5e-6, 2.5, 4, 3.)
    print(r.getTotalPhase(1550e-9))