"""
Simulation files that controls and handles all the simulations in lumerical DEVICE and Interconnect

Author      : Simon Belanger-de Villers (simon.belanger-de-villers.1@ulaval.ca)
Created     : February 20 2020
Last edited : Febroary 20 2020
"""

import matplotlib.pyplot as plt
import numpy as np
import os 
import math

def plotIVCurve(filename, saveToPDF=True):
    """ Get the IV data for a PN junction from Lumerical DEVICE and plot it. """

    V, I = [], []
    f = open(filename,'r')
    for line in f.read().split('\n')[0:-1]:
        V.append(float(line.split('\t')[0]))
        I.append(float(line.split('\t')[1]))
    f.close()
    V,I = np.asarray(V),np.asarray(I)

    # Plot the IV curve
    plt.plot(V, I/1e-3, '-k', linewidth=1.5)
    plt.xlabel('Bias [V]')
    plt.ylabel('Current [mA]')
    plt.grid()
    plt.axis([-1,1,-0.1,1.1])
    if saveToPDF:
        plt.savefig('figures/IVCurve.pdf')
    plt.show()

def plotCapacitanceCurve(filename, saveToPDF=True):
    """ Get the capacitance data for a PN junction from Lumerical DEVICE and plot it. """

    V, C = [], []
    f = open(filename,'r')
    for line in f.read().split('\n')[0:-1]:
        V.append(float(line.split('\t')[0]))
        C.append(float(line.split('\t')[1]))
    f.close()
    V,C = np.asarray(V),np.asarray(C)

    plt.plot(V, C, '-k', linewidth=1.5)
    plt.xlabel('Bias [V]')
    plt.ylabel('Capacitance [pF cm-1]')
    plt.grid()
    #plt.axis([-1,1,-0.1,1.1])
    if saveToPDF:
        plt.savefig('figures/capacitance.pdf')
    plt.show()

def plotBandwidthFromRC(capFilename, resFilename, saveToPDF=True):
    """ Plot the 3dB RC bandwidth from Rslab measurement and capacitance measurement. 
        measuredRes in [ohm cm]
    """
    # Extract capacitance vs bias voltage
    V, C = [], []
    f = open(capFilename,'r')
    for line in f.read().split('\n')[0:-1]:
        V.append(float(line.split('\t')[0]))
        C.append(float(line.split('\t')[1]))
    f.close()
    V,C = np.asarray(V),np.asarray(C)

    # Extract slab resistance
    f = open(resFilename,'r')
    for line in f.read().split('\n')[0:-1]:
        R = float(line.split('\t')[0])
    f.close()

    R = measuredRes # [ohm cm]
    plt.plot(V, 1e-9/(2*math.pi*R*C*1e-12), '-k', linewidth=1.5)
    plt.xlabel('Bias [V]')
    plt.ylabel('3dB RC Bandwidth [GHz]')
    plt.grid()
    #plt.axis([-1,1,-0.1,1.1])
    plt.savefig('figures/bandwidth.pdf')
    plt.show()

def plotOpticalProperties(filename, saveToPDF=True):
    """ Extract the effective index + propagation losses [dB/m] from FEEM simulations. """

    wavelength = 1550e-9 # [m]

    V, neff = [], []
    f = open(filename,'r')
    for line in f.read().split('\n')[0:-1]:
        V.append(float(line.split('\t')[0]))
        neff.append(float(line.split('\t')[1]) + 1j*float(line.split('\t')[2]))
    f.close()
    V, neff = np.asarray(V), np.asarray(neff)

    # Plot effective index vs reverse bias
    plt.plot(-1*V, np.real(neff),'-sk')
    plt.xlabel('Reverse bias [V]')
    plt.ylabel('Effective index')
    plt.grid()
    plt.show()

    # Plot propagation losses vs reverse bias
    loss = -0.2 * np.log10(np.exp(1))*(-2 * math.pi * np.imag(neff)/wavelength)
    plt.plot(-1*V, loss, '-sk') # dB/cm
    plt.xlabel('Reverse bias [V]')
    plt.ylabel('Propagation loss [dB/cm]')
    plt.grid()
    plt.show()

    # Measurement of Vpi (without considering the PN junction length/PN filling factor)
    R = 10e-6   # Ring radius [m]
    L = 2 * math.pi * R
    L = 10000e-6
    k0 = 2 * math.pi /wavelength
    
    dneff_pi = math.pi/(k0 * L)
    print(dneff_pi)

    # Plot the dneff
    dneff = np.real(neff-neff[0])
    plt.plot(-1*V, dneff,'-sk')
    plt.axhline(y=dneff_pi)
    plt.xlabel('Reverse bias [V]')
    plt.ylabel('Effective index change [x1E-4]')
    plt.grid()
    plt.show()

    # Plot the efficiency of the pn junction
    Lpi     = wavelength/dneff/0.02
    VpiLpi  = -V*Lpi
    plt.plot(-V, VpiLpi,'-sk')
    plt.plot()
    plt.xlabel('Reverse bias [V]')
    plt.ylabel('VpiLpi [V*m]')
    plt.grid()
    plt.show()


if __name__ == '__main__':

    # Inputs
    # IVCurve.txt -> I vs V for the PN junction
    # capacitance.txt -> C vs V for the PN junction
    # Resistance -> Rslab does not change with resistance

    plotIVCurve('IVCurve.txt', saveToPDF=True)
    plotCapacitanceCurve('capacitance.txt', saveToPDF=True)
    plotBandwidthFromRC('capacitance.txt', 'resistance.txt', saveToPDF=True)
    plotOpticalProperties('opticalProps.txt', saveToPDF=True)