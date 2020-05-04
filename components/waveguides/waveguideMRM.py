""" Basic class for a waveguide to be used with the Microring Modulator."""

from math import pi, exp, sqrt

# TODO : Include waveguide geometric parameters and make it possible to simulate from MODE
# TODO : Make a children of this waveguide with more parameters

class waveguideMRM(object):
    refWavelength = 1550e-9 # Reference wavelength for the computation of the propagation constant (taylor expansion) [m]
    def __init__(self, effectiveIndex=2.5, groupIndex=4.5, attenuationCoefficient=46, length=1e-6):
        self.effectiveIndex         = effectiveIndex
        self.groupIndex             = groupIndex
        self.attenuationCoefficient = attenuationCoefficient
        self.length                 = length

    # Methods
    def getPhaseShift(self, wavelength):
        """ Simple non dispersive propagation in the waveguide. """
        beta = (2 * pi * self.effectiveIndex) / self.refWavelength + (2 * pi * self.groupIndex) * (1/wavelength - 1/self.refWavelength) # Propagation constant[m-1]
        return beta * self.length

    def getAmplitudeTransmission(self, wavelength):
        """ Simple Amplitude transmission (field) in the waveguide. """
        return exp(-1/2 * self.attenuationCoefficient * self.length)

    def solveEigenModes(self):
        """ Solve for the Eigenmodes of this waveguide. Obtain effective index, group index and losses. """
        pass