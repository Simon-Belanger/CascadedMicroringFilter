""" All-Pass Microring Modulator analytical model from W. Bogaerts paper 'Silicon microring resonators' ."""

import sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from analytical.MRR_AP import *
from misc.dataIO import saveData, loadData
import os, math
import numpy as np
import matplotlib.pyplot as plt

## Things to do
# TODO: Group index from mode solver/ effective index method
# TODO: Evaluate the losses from bends + propag + pn junction
# TODO: Corner Analysis 
# TODO: Photoconcuctive heater in the model
# TODO: Interpolate the data from the triangular mesh on a rectilinear grid and apply the effective index method (https://matplotlib.org/gallery/images_contours_and_fields/triinterp_demo.html#sphx-glr-gallery-images-contours-and-fields-triinterp-demo-py)
# TODO: Critical coupling is for a given bias (Operation point)? Or is it measured for 0V Bias?

class MRM_Static(MRR_AP):
    """ Static model of a microring modulator. 

    Example : see mainMRM.py
    """
    _modLength = 1 # Dummy value
    wavelengthRange = {'start':1552e-9, 'stop':1553e-9, 'pts':1000}     # Wavelength range for sweeps
    offsetGap = 10e-9           # Variation in gap in order to achieve undercoupling/overcoupling [m]

    def __init__(self, radius, loss_per_cm, n_eff, n_g, gap, pnCoverage, coupler, phaseshifter):
        C               = 200e-9 # Dummy value
        self.wavelength      = 1550e-9
        self.alphaMod   = 0
        self.coupler    = coupler
        self.phaseShifter = phaseshifter
        self.intrinsicLoss = None           # Propagation loss of the doped waveguide [dB/cm]
        super().__init__(radius, loss_per_cm, n_eff, n_g, C)
        self.modLength  = pnCoverage

        if radius != 10e-6:
            print('Error: the coupler data is only for radius = 10 microns.');exit()

        self.gap = gap

    # Properties / Attributes
    @property
    def modLength(self):
        " Getter for the PN Junction modulator length [m]. "
        return self._modLength
    @modLength.setter
    def modLength(self, pnCoverage):
        " Setter for the PN Junction modulator length [m]. "
        self._modLength = self.L_rt * pnCoverage / 100

    @property
    def gap(self):
        " Getter for the directional coupler gap [m]. "
        return self._gap
    @gap.setter
    def gap(self, gap):
        " Setter for the directional coupler gap [m]. "
        self.setBias(0) # TODO:  To have the losses of the modulator included (intrinsic vs perturbative)
        criticalGap = self.findCriticalCouplingCondition()
        if gap == 'CC': # Critically Coupled
            self._gap = criticalGap
        elif gap == 'UC': # Under Coupled
            self._gap = criticalGap + self.offsetGap
        elif gap == 'OC': # Over Coupled
            self._gap = criticalGap - self.offsetGap
        else:
            self._gap = gap
        self.setCrossCoupling(self._gap)

    # Methods
    def get_phase(self, lambda_0):
        """ Obtain the roundtrip phase at a particular wavelength from the effective index and the roundtrip length of
        the ring resonator. """
        beta0   = (2 * math.pi * self.n_eff) / self.ref_wavelength + (2 * math.pi * self.n_g) * (1/lambda_0 - 1/self.ref_wavelength) # Propagation constant for undoped wg [m-1]
        betaMod = (2 * math.pi * self.n_effMod) / self.ref_wavelength + (2 * math.pi * self.n_g) * (1/lambda_0 - 1/self.ref_wavelength) # Propagation constant for doped wg [m-1]
        phi     = (betaMod * self.modLength) + (beta0 * (self.L_rt - self.modLength)) # Optical phase [rad]
        return phi

    def get_roundtrip_loss(self):
        """ Obtain the roundtrip losses from the loss coefficient and the roundtrip length of the ring resonator. """
        baseLoss    = math.exp(- self.alpha * (self.L_rt - self.modLength))   # Intrinsic loss 
        modLoss     = math.exp(- 23 * self.alphaMod * self.modLength)         # PN junction related loss
        self.a      = math.sqrt(baseLoss * modLoss)

    def getCrossCoupling(self):
        """ Get the value of the cross coupling coefficient """
        return self.C

    def setCrossCoupling(self, gap):
        """ Set the cross coupling coefficient in power. """
        self.C = self.coupler.measureCouplingCoefficient(gap, self.wavelength)
        self.r = self.get_r(self.C)

    def setBias(self, bias):
        """ Set the reverse bias applied to the PN junction. """

        # Interpolate the complex effective index to the selected reverse bias
        neffInterp = np.interp(-bias, -self.phaseShifter.phaseShifter['V'], self.phaseShifter.phaseShifter['neff'])

        # Set the new effective index and losses
        self.n_effMod       = neffInterp.real
        self.alphaMod       = -0.2 * np.log10(np.exp(1))*(-2 * math.pi * neffInterp.imag/self.wavelength) # [dB/cm]
        self.get_roundtrip_loss()

    def findCriticalCouplingCondition(self):
        """ Find the gap for which the cross-coupling coefficient satisfies the critical coupling condition. """
        targetGap = self.coupler.reverseInterpolateCouplingCoefficient((1 - self.a**2))
        #print("The gap required to achieve critical coupling condition is {} nm".format(targetGap*1e9))
        return targetGap

    def measureTransmissionMRM(self):
        " Measure the MRM transmission vs DV Reverse Bias "

        bias = np.linspace(0, -4.0, 10)

        Q, BW = [], []
        fig1 = plt.figure(1);plt.title('Transmission');plt.xlabel('Wavelength [nm]');plt.ylabel('Transmission [dB]')
        fig2 = plt.figure(2);plt.title('Phase');plt.xlabel('Wavelength [nm]');plt.ylabel('Phase [rad]')
        for Vi in bias:
            self.setBias(Vi)
            wvl, power = self.sweep_power_transmission(self.wavelengthRange['start'], self.wavelengthRange['stop'], self.wavelengthRange['pts'])
            plt.figure(1);plt.plot(wvl*1e9, power, label='{:0.1f} V'.format(Vi))
            wvl, phase = self.sweep_phase(self.wavelengthRange['start'], self.wavelengthRange['stop'], self.wavelengthRange['pts'])
            plt.figure(2);plt.plot(wvl*1e9, phase, label='{:0.1f} V'.format(Vi))
            # Measure the Q factor 
            Q.append(self.getQfactor(1550e-9))
            BW.append(self.getOpticalBandwidth(wavelength=1550e-9))
        plt.figure(1);plt.legend();plt.grid();plt.figure(2);plt.legend();plt.grid()
        plt.show()

        return bias, Q, BW
"""
    def measureTotalBandwidth(VRCBW, RCBW, VoptBW, optBW):
        " Measure the total bandwidth from RC bandwidth and optical bandwidth. "
        # Interpolate the data to a given bias 
        outBias = np.linspace(0, 4, 100)
        outRCBW = np.interp(outBias, VRCBW, RCBW)
        outOptBW = np.interp(outBias, VoptBW, optBW)

        tbw = []
        for rcbw,obw in zip(outRCBW, outOptBW):
            tbw.append(math.sqrt(1/(1/rcbw**2 + 1/obw**2)))
        return outBias, tbw
"""