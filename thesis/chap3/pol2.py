"""
Work in progress. Implementation of Jones matrix for an optical system.
"""
from numpy import matrix, cos, sin, exp
import sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.utils import attenuationCoefficientFromWaveguideLosses
from math import pi
import scipy.constants as cte
from numpy.linalg import norm

def getPropagationConstant(wavelength, effectiveIndex, groupIndex, dispersion):
    " Getter for the propagation constant, dispersive. "
    refWavelength = 1550e-9

    dOmega = ((2*pi*cte.c)/float(wavelength) - (2*pi*cte.c)/refWavelength)
    firstOrder  = (2*pi*cte.c)/refWavelength * (float(effectiveIndex)/cte.c)
    secondOrder = float(groupIndex)/cte.c
    thirdOrder  = 0.5 * (- float(dispersion) * 1e-6 * refWavelength**2 / 2 * pi * cte.c )
    return firstOrder + secondOrder * dOmega + thirdOrder * dOmega**2

"""
Basic classes : Jones vector for field  and Jones matrix for system
"""
class jonesVector(object):
    " Jones vector describing the electric field as polarized light. "
    def __init__(self, Ex, Ey):
        self.Ex = Ex
        self.Ey = Ey

    @property
    def vector(self):
        return matrix([[self.Ex], [self.Ey]])

    @property
    def magnitude(self):
        return norm(self.vector)

    def propagate(self, jonesMatrix):
        " Propagate the given field into the given optical element given by the jones matrix. "
        outputVector = jonesMatrix.mat * self.vector
        return jonesVector(outputVector[0, 0], outputVector[1, 0])

class jonesVectorAngle(jonesVector):
    """ Jones vector describing the electric field as polarized light, where theta and phi are parameters of the poincare shpere. """
    def __init__(self, theta, phi):
        self.theta  = theta
        self.phi    = phi
    
    @property
    def Ex(self):
        return cos(self.theta) * exp(1j * self.phi)

    @property
    def Ey(self):
        return sin(self.theta) * exp(-1j * self.phi)

class jonesMatrix(object):
    " Jones matrix describing a general optical system. "
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    @property
    def mat(self):
        return matrix([[self.a , self.b], [self.c, self.d]])

"""
Actual components
"""
class buffer(jonesMatrix):
    " Jones matrix describing a buffer ."
    def __init__(self, length, betaX, alphaX, betaY, alphaY):

        a = 1
        b = 0
        c = 0
        d = 1
        super().__init__(a, b, c, d)

class waveguide(jonesMatrix):
    " Jones matrix describing a simple birefringent waveguide ."
    def __init__(self, length, betaX, alphaX, betaY, alphaY):

        a = exp((-1j * betaX - alphaX) * length)
        b = 0
        c = 0
        d = exp((-1j * betaY - alphaY) * length)
        super().__init__(a, b, c, d)

class polSplitterRotator(object):
    """ Polarization splitter rotator transfer matrix. Very simple model. 
    
    Fields reference:
            E_o      =   PSR        *   E_i

    # Forward direction

       [[   Ecw,x  ],   [[a, b],
        [   Ecw,y  ],    [c, d],    [[   Ei,x  ],
        [   Eccw,x ],    [e, f],     [   Ei,y ]]
        [   Eccw,y ]]    [g, h]]

        Clockwise branch            : Ei_x --> Ecw,x    Splitting
        Counterclockwise branch     : Ei_y --> Eccw,x   Splitting + Rotation

    # Backward direction

       [[   Ecw,x  ],   [[a, b],
        [   Ecw,y  ],    [c, d],    [[   Ei,x  ],
        [   Eccw,x ],    [e, f],     [   Ei,y ]]
        [   Eccw,y ]]    [g, h]]

        Clockwise branch            : Ecw_x     --> Ei_y    Splitting + rotation
        Counterclockwise branch     : Eccw,x    --> Ei_x    Splitting
    """
    def __init__(self, transmission=1, crosstalk=0):

        self.t = transmission
        self.c = crosstalk

    @property
    def forwardMat(self):
        return matrix([[self.t, self.c],[self.c, self.c],[self.c, self.t],[self.c, self.c]])

    @property
    def backwardMat(self):
        return matrix([[self.c, self.c, self.t, self.c],[self.t, self.c, self.c, self.c]])
    
    def propagateForward(self, E_i):
        outputVector = self.forwardMat * E_i.vector
        E_cw    = jonesVector(outputVector[0, 0], outputVector[1, 0])
        E_ccw   = jonesVector(outputVector[2, 0], outputVector[3, 0])
        return E_cw, E_ccw

    def propagateBackward(self, E_cw, E_ccw):
        outputVector = self.backwardMat * matrix(([E_cw.Ex], [E_cw.Ey], [E_ccw.Ex], [E_ccw.Ey]))
        E_r    = jonesVector(outputVector[0, 0], outputVector[1, 0])
        return E_r

"""
Actual code
"""
# Properties of the waveguide for both TE and TM modesvelength =
wavelength  = 1550e-9
length      = 100e-6
modeTE = {'effectiveIndex': 2.44, 'groupIndex': 4.18, 'dispersion': 0, 'wgLosses': 300} #530
modeTM = {'effectiveIndex': 1.78, 'groupIndex': 3.805, 'dispersion': 0, 'wgLosses': 200} # -2.1e4

betaTE  = getPropagationConstant(wavelength, modeTE['effectiveIndex'], modeTE['groupIndex'], modeTE['dispersion'])
betaTM  = getPropagationConstant(wavelength, modeTM['effectiveIndex'], modeTM['groupIndex'], modeTM['dispersion'])
alphaTE = attenuationCoefficientFromWaveguideLosses(modeTE['wgLosses'], 'field')
alphaTM = attenuationCoefficientFromWaveguideLosses(modeTM['wgLosses'], 'field')


# Full System model
if True:

    # Declare components
    PSR = polSplitterRotator(transmission=1, crosstalk=0)   # Lossless PSR perfect
    L   = waveguide(length, betaTE, alphaTE, betaTM, alphaTM)

    # Input Field
    E_i = jonesVectorAngle(theta=0, phi=0)

    # Initial -> Loop (PSR = TMat of the PSR)
    E_cw, E_ccw = PSR.propagateForward(E_i)

    # Porpagation inside the loop (L = TMat of the loop)
    E_cw2    = E_cw.propagate(L)
    E_ccw2   = E_ccw.propagate(L)

    # Loop --> Output (PSR = TMat of the PSR)
    E_r = PSR.propagateBackward(E_cw2, E_ccw2)
    print(E_r.magnitude)
