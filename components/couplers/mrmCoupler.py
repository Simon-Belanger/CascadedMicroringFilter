""" Directional coupler that is used for MRM circuit simulation.
"""
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.dataIO import saveData, loadData
from scipy import interpolate

# TODO : Obtain measurements from Simulation/ Do everything in FDTD
# TODO : Convergence testing
# TODO : Perform FDTD simulation and return the results
# TODO : Make it so it is not hardcoded 
# TODO : Add input parameters for FDTD simulations instead of being only a datafile load

class mrmDC(object):
    def __init__(self, couplerFilename):
        """ Constructor for the directionnal coupler object."""
        self.coupler = loadData(couplerFilename)

    def simulateCouplerUsing3DFDTD(self):
        """ Perform a 3D FDTD Simulation to measure the cross-coupling coefficient for a given waveguide geometry. """
        pass

    def plotCouplingCoefficient(self):
        """ Plot the coupling coefficient vs gap for all the wavelengths present in the datafile. [HARDCODED]"""
        gap     = np.asarray(self.coupler['1500']['gap'])
        kappa   = np.asarray([np.abs(self.coupler['1500']['k'])**2, np.abs(self.coupler['1548']['k'])**2, np.abs(self.coupler['1600']['k'])**2])
        plt.plot(gap*1e9, kappa[0,:], 'ok',label='1500 nm (3DFDTD)')
        plt.plot(gap*1e9, kappa[1,:], 'or',label='1548 nm (3DFDTD)')
        plt.plot(gap*1e9, kappa[2,:], 'ob',label='1600 nm (3DFDTD)')
        newgap = np.linspace(gap[0].real, gap[-1].real, 100)
        plt.plot(newgap*1e9, self.interpolateGapAndWavelength(newgap, 1500e-9),'-k', label='1500 nm (Interp.)')
        plt.plot(newgap*1e9, self.interpolateGapAndWavelength(newgap, 1548e-9),'-r', label='1548 nm (Interp.)')
        plt.plot(newgap*1e9, self.interpolateGapAndWavelength(newgap, 1600e-9),'-b', label='1600 nm (Interp.)')
        plt.xlabel('Gap [nm]');plt.ylabel('Power cross coupling coefficient [-]')
        plt.legend();plt.grid()
        plt.savefig('couplingCoefficient.pdf');plt.show()

    def interpolateGapAndWavelength(self, targetGap, targetWavelength):
        """ Perform 2D interpolation of the power cross-coupling coefficient (k^2) vs the gap [m] and the wavelength [m] [HARDCODED]"""
        wvl     = np.asarray([1500,1548,1600]) * 1e-9
        gap     = np.asarray(self.coupler['1500']['gap'])
        kappa   = np.asarray([np.abs(self.coupler['1500']['k'])**2, np.abs(self.coupler['1548']['k'])**2, np.abs(self.coupler['1600']['k'])**2])
        f       = interpolate.interp2d(gap, wvl, kappa)
        return np.squeeze(f(targetGap, targetWavelength))

    def reverseInterpolateCouplingCoefficient(self, powerCrossCouplingCoefficient):
        """ Find the gap corresponding """
        # TODO : Interpolate with the wavelength as well instead of choosing 1548 nm by default
        # TODO : makes sure the code always work (boundary conditions, etc)

        gap     = np.asarray(self.coupler['1548']['gap'])
        C       = np.abs(self.coupler['1548']['k'])**2
        return np.interp(-powerCrossCouplingCoefficient, -C, gap).real

    def measureCouplingCoefficient(self,gap, wavelength):
        """ Dummy function for now. """
        return self.interpolateGapAndWavelength(gap, wavelength)

if __name__ == '__main__':
    # Example usage

    import os
    filename = os.getcwd()+'/../../MRM/data/directionalCoupler/coupler1'
    dirCoup = mrmDC(filename)

    # Plot the coupling coefficient vs gap - For figures, reports, etc
    dirCoup.plotCouplingCoefficient()

    # Get the power cross-coupling coefficient for a given gap/wavelength combination
    print(dirCoup.interpolateGapAndWavelength(200e-9, 1550e-9))

