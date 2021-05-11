"""
Model for a Cascaded Microring Filter using the Transfer Matrix Method (TMM)
    - Bidirectionnal (4x4 transfer matrices)

Author      : Simon Belanger-de Villers (simon.belanger-de-villers.1@ulaval.ca)
Created     : 2018
Last edited : August 13th 2020
"""

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt
from misc.utils import t2s, listmat_multiply
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB
from misc.Crosstalk import build_crosstalk_matrix
from math import pi, sqrt, log10, exp
import cmath as cm
import random
from components.phase_shifter import heaterBasic
from components.ring import Ring
from components.couplers.apodization import *
from components.couplers.virtual_DC import virtual_DC
from misc import constants

# Future modifications:
# TODO : Be able to import different types of waveguides (simple analytical, MODE simulation, etc).
# TODO : Consider dispersion (at least first order)
# TODO : make a waveguide object and associate it to the ring object.
# TODO : Be able to use different types of couplers (simple analytical, FDTD simulations).
# TODO : Make a OADM class that would be parent of the MRF class with filter results and plotting methods.
# TODO : Make documentation for the code on github. readme
# TODO : Compare with EMPy for code structure and functionnalities
# TODO : Pass parameters as **kwargs instead (set defaults in the function)


# Matrix used to flip ports for parity
flipMat = lambda x : np.matrix([[x, 1-x, 0, 0], [1-x, x, 0, 0], [0, 0, x, 1-x], [0, 0, 1-x, x]])

# Default MRF class, rings and couplers can be different from one another
class parentMRF(object):
    """ Microring Filter Class, generates a model for a high-order microring filter. """

    def __init__(self, rings=[None]*5, couplers=[None]*6, crosstalk_coeff=[1, 0., 0., 0., 0.]):
        """ Constructor for the cascaded microring filter object. """
        
        try:
            self.Rings      = rings     # Ring resonator list
            self.couplers   = couplers  # Directionnal couplers list
            if not(len(self.couplers) == len(self.Rings) + 1):
                raise filterError('There should be one coupler more than the number of rings.')
        except filterError as err:
            print('filterError: {0}'.format(err));exit()


        # Additional phase (zero by default)
        self.bias                   = np.zeros(len(self.Rings))  # Bias applied to each phase shifter [V]
        self.desired_tuning_phase   = np.zeros(len(self.Rings))  # Desired phase tuning if there was no thermal crosstalk [rad]
        self.actual_tuning_phase    = np.zeros(len(self.Rings))  # Actual phase tuning after thermal crosstalk [rad]

        # Phase shifting
        self.phaseshifters          = [heaterBasic(phaseEfficiency=20e-3, resistance=600)] * len(self.Rings)
        if len(self.Rings) > 1:
            phaseCouplingMat  = build_crosstalk_matrix(crosstalk_coeff, len(self.Rings))

            # Normalization of the matrix
            row_sums = phaseCouplingMat.sum(axis=1)
            self.phase_coupling_matrix = phaseCouplingMat / row_sums[:, np.newaxis]
        else:
            self.phase_coupling_matrix  = 1

        # For algorithms
        self.NM_phase, self.NM_power = [], []
        self.target_wavelength = 1550e-9 # Default value

        # For measurements
        self.clipping = True

    def get_total_phase(self, wavelength):
        """ Get the modulus for the total phase for the cascaded ring. """
        inner_product = 0
        for ring in self.Rings:
            inner_product += (ring.getTotalPhase(wavelength))**2
        try:
            return log10(sqrt(inner_product)%(2*pi))
        except ValueError:
            return -100

    def manufacturing(self, maximumInterval):
        """ Add a random phase perturbation for each ring that can be attributed to manufacturing variability. """
        for ring in self.Rings:
            ring.phaseDeviation = random.uniform(-maximumInterval/2, maximumInterval/2)

    def manufacturingValue(self, phaseArray):
        """ Add a random phase perturbation for each ring that can be attributed to manufacturing variability. """
        
        if len(phaseArray)!=len(self.Rings):
            print('error');exit()
        
        for ring, phase in zip(self.Rings, phaseArray):
            ring.phaseDeviation = phase

    def TMM(self, wavelength, E_in, E_add=None, E_drop=None, E_thru=None):
        """
        Calculate the response of the microring filter to an input signal. The signal can be single value 
        or a list/ndarray.

        Inputs:
            wavelength  : Wavelength grid for the input/add signals [nx1 list]
            E_in        : Electric field at the input port vs wavelength [nx1 list]
            E_add       : Electric field at the add port vs wavelength [nx1 list]
            E_through   : Electric field at the through port vs wavelength [nx1 list]
            E_drop      : Electric field at the drop port vs wavelength [nx1 list]

        """

        if type(wavelength) in [np.ndarray, list]: # Wavelength array

            # Initialise drop, thru and add ports inputs if it has not been already set. 
            if E_thru is None:
                E_thru = np.zeros(len(wavelength), dtype=complex)  # through-port input
            if E_drop is None:
                E_drop = np.zeros(len(wavelength), dtype=complex)  # drop-port input
            if E_add is None:
                E_add = np.zeros(len(wavelength), dtype=complex)  # add-port input

            # Run the function for each wavelength
            outFields = fields(E_in=None, E_drop=np.zeros(len(wavelength), dtype=complex), E_thru=np.zeros(len(wavelength), dtype=complex), E_add=None)
            for ii in range(len(wavelength)):
                outputFields = self.TMM(wavelength[ii], E_in[ii], E_add[ii], E_drop[ii], E_thru[ii])
                outFields.through[ii], outFields.drop[ii] = outputFields.through, outputFields.drop
            return outFields

        elif type(wavelength) in [int, float, np.float64]: # Single wavelength

            # Initialise drop, thru and add ports inputs if it has not been already set. 
            if E_thru is None:
                E_thru = 0.  # through-port input
            if E_drop is None:
                E_drop = 0.  # drop-port input
            if E_add is None:
                E_add = 0.  # add-port input

            listMat     = self.buildTotalTransferMatrix(wavelength) # Order the matrices
            outputMat   = t2s(listmat_multiply(listMat) * flipMat(len(self.Rings)%2)) * np.matrix([[E_thru], [E_in], [E_add], [E_drop]]) # Matrix multiplication
            outputFields = fields(E_in=outputMat[1,0], E_drop=outputMat[3,0], E_thru=outputMat[0,0], E_add=outputMat[2,0])
            if self.clipping == True:
                outputFields = fields(E_in=self.clipMeasuredPower(outputMat[1,0]), E_drop=self.clipMeasuredPower(outputMat[3,0]), E_thru=self.clipMeasuredPower(outputMat[0,0]), E_add=self.clipMeasuredPower(outputMat[2,0]))
            return outputFields

        else:
            print('The input data format is not supported. Must be either numpy array, list, int or float.')

    def buildTotalTransferMatrix(self, wavelength):
        """ Order the matrices for T matrix multiplication by interleaving them. """
        listmat = [None] * (len(self.Rings) + len(self.couplers))
        listmat[0:len(listmat):2]   = [coupler.get_transfer_matrix(wavelength)  for coupler in self.couplers]
        listmat[1:len(listmat)-1:2] = [ring.get_transfer_matrix(wavelength)     for ring    in self.Rings]
        return listmat

    def apply_phase_tuning(self, phaseList):
        """ Apply the tuning phase to all rings. """
        self.actual_tuning_phase = phaseList
        for ring, phase in zip(self.Rings, phaseList):
            ring.tuningPhase = phase

    def applyPhaseTuningCrosstalk(self, phaseList):
        """ Apply the phase tuning but consider the crosstalk. """
        self.apply_phase_tuning(np.asarray(np.squeeze(self.phase_coupling_matrix * np.transpose(np.asmatrix(phaseList))))[0])

    def apply_bias(self, ringID, voltage):
        """ Apply voltage to the corresponding phase shifter and update the tuning phase variable accordingly.
            Phase crosstalk is considered through the phase coupling matrix.

            Inputs:
                ring_id : Heater for which voltage is applied
                voltage : Voltage to apply to the i-th heater [V]
        """

        # Measure desired phase shift
        self.bias[ringID]                  = voltage
        self.desired_tuning_phase[ringID]  = self.phaseshifters[ringID].applyVoltage(voltage)

        # Obtain actual phase shift by including phase crosstalk
        self.applyPhaseTuningCrosstalk(self.desired_tuning_phase)

    def applyPower(self, ringID, inputPower):
        """ Apply DC power to the corresponding phase shifter and update the tuning phase variable accordingly.
            Phase crosstalk is considered through the phase coupling matrix.

            Inputs:
                ring_id     : Heater for which DC power is applied
                inputPower  : DC Power to apply to the i-th heater [W]
        """

        # Measure desired phase shift
        self.bias[ringID]                  = np.sqrt(inputPower*self.phaseshifters[ringID].resistance)
        self.desired_tuning_phase[ringID]  = self.phaseshifters[ringID].getPhaseShift(inputPower)

        # Obtain actual phase shift by including phase crosstalk
        self.applyPhaseTuningCrosstalk(self.desired_tuning_phase)

    @staticmethod
    def clipMeasuredPower(measuredField):
        """ Clip the measured power if it is lower than the threshold of the noise floor. """
        noiseFloorMean, noiseFloorAmplitude        = -80, 2     #[dBm]

        if  field2PowerdB(measuredField) < noiseFloorMean:
            measuredField = powerdB2Field(noiseFloorMean + random.uniform(-noiseFloorAmplitude/2, noiseFloorAmplitude/2))
        return measuredField

    def objectiveFunction(self, wavelength, voltageArray):
        """
        function to optimize for the MRF.
        """
        # Apply bias
        for ringID in range(len(self.Rings)):
            self.apply_bias(ringID, voltageArray[ringID])
        # Store total phase in an array
        self.NM_phase.append(self.get_total_phase(wavelength))
        self.NM_power.append(10 * np.log10( np.abs(self.TMM(wavelength, E_in=1.).drop)**2 ))

        # Return power
        return self.TMM(wavelength, E_in=1.).drop

    def minimize_MRF(self, voltageArray):
        """"""
        return -1 * self.objectiveFunction(self.target_wavelength, voltageArray)

    @staticmethod
    def measurePhaseResponse(wavelength, fields):
        """ Measure the phase of the Drop/Through fields. """
        phiThru, phiDrop = np.angle(fields.through), np.angle(fields.drop)
        return phiThru, phiDrop

    @staticmethod
    def measureGroupDelay(wavelength, fields):
        """"""
        omega = constants.c / wavelength * 2*pi
        tauThru = - np.diff(np.angle(fields.through)) / np.diff(omega)
        tauDrop = - np.diff(np.angle(fields.drop)) / np.diff(omega)
        return tauThru, tauDrop

    @staticmethod
    def plotTransmission(wavelength, fields):
        """ Plot the drop/thru port transmission spectrum. """
        plt.plot(wavelength * 1e9, 10 * np.log10(abs(fields.drop) ** 2), label='Drop', linewidth=2)
        plt.plot(wavelength * 1e9, 10 * np.log10(abs(fields.through) ** 2), label='Through', linewidth=2)
        plt.xlabel('wavelength (nm)', fontsize=14);plt.ylabel('Power transmission (dB)', fontsize=14)
        plt.grid(); plt.legend(loc='upper right'); plt.show()

    @staticmethod
    def plotPhaseResponse(wavelength, phiThru, phiDrop, unit='rad'):
        """ Plot the drop/thru port phase response. """
        if unit=='rad':
            plt.plot(wavelength * 1e9, -np.unwrap(2*phiDrop) / 2, label='Drop', linewidth=2)
            plt.plot(wavelength * 1e9, -np.unwrap(2*phiThru) / 2, label='Through', linewidth=2)
            plt.xlabel('wavelength (nm)', fontsize=14); plt.ylabel('Phase (rad)', fontsize=14)
        elif unit=='pi':
            plt.plot(wavelength * 1e9, -np.unwrap(2*phiDrop) / 2*np.pi, label='Drop', linewidth=2)
            plt.plot(wavelength * 1e9, -np.unwrap(2*phiThru) / 2*np.pi, label='Through', linewidth=2)
            plt.xlabel('wavelength (nm)', fontsize=14); plt.ylabel('Phase (\\pi)', fontsize=14)
        plt.grid(); plt.legend(); plt.show()

    @staticmethod
    def plotGroupDelay(wavelength, tauThru, tauDrop):
        """ Plot the drop/thru port group delay. """
        plt.plot(wavelength[0:len(wavelength) - 1] * 1e9, tauDrop * 1e12, label='Drop', linewidth=2)
        plt.plot(wavelength[0:len(wavelength) - 1] * 1e9, tauThru * 1e12, label='Through', linewidth=2)
        plt.xlabel('wavelength (nm)', fontsize=14); plt.ylabel('Group delay (ps)', fontsize=14)
        plt.grid(); plt.legend(); plt.show()

# Cascaded ring filter where the rings are all the same
class MRF(parentMRF):
    " MRF object for which the ring radius "
    def __init__(self, num_rings=5, radius=2.5e-6, neff=2.4449, ng=4.18, alpha_wg=3., couplers=[None, None, None, None, None, None], crosstalk_coeff=[1, 0., 0., 0., 0.]):
        
        # Make all identical ring resonator objects from input parameters
        rings      = [Ring(radius, neff, ng, alpha_wg) for i in range(num_rings)]

        # Create the MRF object with the given attributes
        super().__init__(rings, couplers, crosstalk_coeff)

# MRF that exhibits a perfect Butterworth response
class idealMRF(MRF):
    " MRF that exhibits a perfect Butterworth response using perfect couplers. "
    def __init__(self, num_rings=5, radius=2.5e-6, neff=2.4449, ng=4.18, alpha_wg=3., powerCoupling=0.1, crosstalk_coeff=[1, 0., 0., 0., 0.]):
        
        # Make the couplers using perfect directional couplers (lossless, no wavelength dependance)
        couplers        = []
        for kappa, loss in zip( flatTopCoupling(order=num_rings, powerCoupling=powerCoupling), [1.]*(num_rings+1) ):
            couplers.append( virtual_DC(kappa, loss) )

        # Create the MRF object with the given attributes
        super().__init__(num_rings, radius, neff, ng, alpha_wg, couplers, crosstalk_coeff)

class fields(object):
    " This datatype contains 4 fields the 4 output fields of a microring filter."
    def __init__(self, E_in, E_drop, E_thru, E_add):
        self.input 		= E_in
        self.through 	= E_thru
        self.add 		= E_add
        self.drop 		= E_drop

    def asMatrix(self):
        # Make a 4x1 matrix out of all the fields
        return np.matrix([[self.through], [self.input], [self.add], [self.drop]])

    def __repr__(self):
        # Dummy printing function for debugging
        return "Input Field \t= {0:.4e}\nDrop Field \t= {1:.4e}\nThrough Field \t= {2:.4e}\nAdd Field \t= {3:.4e}\n".format(np.abs(self.input), np.abs(self.drop), np.abs(self.through), np.abs(self.add))

# Exceptions
class filterError(Exception): 
    pass