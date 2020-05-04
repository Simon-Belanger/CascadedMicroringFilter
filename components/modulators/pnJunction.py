""" Class that implements a PN junction that is simulated in Lumerical DEVICE. """

import sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.lumerical import *
from misc.dataIO import saveData, loadData
import os, math
import matplotlib.pyplot as plt
import numpy as np
import meshio
import matplotlib.tri as mtri

# Lumerical DEVICE Parameters
neffSweepDefaults = {'edge_per_wvl': 3, 'poly_order': 3, 'refinement': 20}
BWsweepDefaults = {'up_abs_tol': 0.0001, 'up_rel_tol': 1e-6, 'res_abs_tol': 0.0001, 'min_step': 0.02, 'max_step': 4.0, 'abs_lim': 0.001,
                    'rel_lim': 0.001, 'min_edge': 0.025, 'max_edge': 0.05, 't_ramp': 1.0, 't_hold': 5.0}

class pnJunction(object):
    """ 
        Example code to plot the phase shifter efficiency and introduced losses
            phaseshifter    = pnJunction('filename')
            phaseshifter.plotPhaseShifterEfficiency()
            phaseshifter.plotPropagationLosses()
    """

    wavelength = 1550e-9    # Target wavelength for simulations [m]

    def __init__(self, phaseShifterFilename=None):

        if phaseShifterFilename==None:
            print('Phase Shifter must be simulated. Not supported yet.')
        else:
            self.phaseShifter = loadData(phaseShifterFilename)
        # Results that are either simulated or loaded from file
        self.resultBias                     = None
        self.resultEffectiveIndex           = None
        self.resultPhaseShifterEfficiency   = None
        self.resultCapacitance              = None
        self.resultResistance               = None
        self.resultBandwidth                = None
        self.resultLosses                   = None
    
    def extractResultsFromFile(self):
        """ Extract the Reverse Bias and effective index data from a data file. """
        pass

    def measurePhaseShifterEfficiency(self):
        " Measure Vpi * Lpi, the phase shifter efficiency [V * cm]. "
        VpiLpi = []
        V = (-np.asarray(self.phaseShifter['V'])).tolist(); neff = self.phaseShifter['neff'].real
        for V_i in V:
            try:
                #print(V_i)
                VpiLpi.append((1e2 * self.wavelength * V_i)/(2*(neff[V.index(V_i)] - neff[0]))) # [V * cm]
            except ZeroDivisionError:
                VpiLpi.append(None)
        self.resultBias                     = V
        self.resultPhaseShifterEfficiency   = np.asarray(VpiLpi).real

    def plotPhaseShifterEfficiency(self):
        " Plot Vpi * Lpi, the phase shifter efficiency vs the reverse bias applied. "
        self.measurePhaseShifterEfficiency()
        plt.plot(self.resultBias, self.resultPhaseShifterEfficiency, '-ok')
        plt.xlabel('Reverse Bias [V]');plt.ylabel('Phase Shifter efficiency VpiLpi [V*cm]')
        plt.grid();plt.savefig('phaseShifterEfficiency.pdf');plt.show()

    def measurePropagationLosses(self):
        """ Measure the losses of the modulator [dB/cm]. """
        pass

    def plotPropagationLosses(self):
        """ Plot the propagation losses vs the reverse bias. """
        loss = -0.2 * np.log10(np.exp(1))*(-2 * math.pi * np.asarray(self.phaseShifter['neff'].imag)/self.wavelength)
        plt.plot(-1*self.phaseShifter['V'], loss, '-sk') # dB/cm
        plt.xlabel('Reverse bias [V]');plt.ylabel('Propagation loss [dB/cm]')
        plt.grid();plt.savefig('propagationLosses.pdf');plt.show()

    def simulateEffectiveIndexVsBiasVoltage(self, wavelength, bias):
        """ Simulate a PN Junction in Lumerical DEVICE using CHARGE + FEEM.
        Returns the modulators effective index vs reverse bias voltage."""

        # TODO : Interpolate the results on the input bias

        with lumapi.DEVICE(hide=False) as device:
            # Pass sweep parameters to DEVICE instance 
            device.putv("V", self.formatBiasForDEVICE(bias))
            device.putv("lambda_0", wavelength)
            device.putv("edge_per_wvl", neffSweepDefaults['edge_per_wvl'])
            device.putv("poly_order", neffSweepDefaults['poly_order'])
            device.putv("refinement", neffSweepDefaults['refinement'])

            # Run the script in DEVICE
            device.feval(os.getcwd() + "/deviceSimul/lumericalScriptFiles/" + "getneff.lsf")

            # Retrieve the results from DEVICE
            V       = np.squeeze(device.getv('V'))      # Voltage on anode (Reverse Bias)                 [V]
            neff    = np.squeeze(device.getv('neff'))   # Effective index vs voltage on anode             [-] 

        return V, neff

    def measureBandwidth(self, device, bias):
        """  """
        # Put deault options
        putOptions(BWsweepDefaults, device)
        device.putv('squarepulse', False)
        device.putv('conv_criteria', 'update')
        device.putv('Datafilename', '')

        BW = []
        for Vi in bias:
            device.putv('V_step', Vi)
            device.feval(os.getcwd() + "/deviceSimul/lumericalScriptFiles/" + "getBW.lsf")
            BW.append(device.getv('BW'))
        return BW

    def visualizeJunction(self, filename, outputGraphics=False, showMesh = False):
        " Visualize the Carrier concentration data initially -> The mesh. "

        # TODO : Do the same to plot the junction carrier density vs voltage for Electric simulations

        mesh = meshio.read(filename)
        x = mesh.points
        triangles = mesh.cells[0][1]
        N = mesh.point_data['N']
        triang = mtri.Triangulation(x[:,0]*1e6, x[:,1]*1e6, triangles)

        if outputGraphics==True:
            plt.figure(figsize=(15, 15))
            if showMesh:    # Show the triangular Mesh on top of the data
                plt.triplot(triang,'-k', alpha=0.3)
            plt.tricontourf(triang, N, 16)
            plt.gca().set_aspect('equal')
            plt.xlabel('x [µm]');plt.ylabel('y [µm]')
            plt.colorbar()
            plt.show()

        return x[:,0], x[:,1], triangles, N

    @staticmethod
    def formatBiasForDEVICE(bias):
        " Format the bias so it gets accepted in Lumerical DEVICE."
        return bias.astype(float)

if __name__ == '__main__':
    pass