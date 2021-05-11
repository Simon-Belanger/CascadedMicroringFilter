"""
mainDetuning.py
Script to generate the figures for the objective function mapping. With and without thermal crosstalk.

Author: Simon Belanger-de Villers
Created : March 14th 2021
"""

# TODO : Write a function that draws the detuning map for N rings
# TODO : Draw the detuning map for a coupled and uncoupled case

import numpy as np 
import matplotlib.pyplot as plt
import math, sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB, field2PowerLinear
from misc.Filter import *
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import *
from components.MRF import idealMRF
from components.ring import Ring
import matplotlib.font_manager as fm
from math import pi
import random
from detuningMapsFuncs import *

overpassFont = fm.FontProperties(family = 'Overpass', 
                                fname = '/Users/simonbelanger/Library/Fonts/overpass-semibold.otf', 
                                size = 12)

mrDef = {'fo':2, 'r': 2.5e-6, 'n':2.4449, 'ng':4.18, 'alf': 3., 'k': 0.1} # Ring default parameters
"""
1a - Detuning V1 : Map the objective function
"""
if False:
    # General parameters
    pdfFilename = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter4/detuning.pdf'
    wavelength=np.linspace(1540e-9, 1580e-9, 10000)

    # Build the dual ring filter from defaults
    crosstalkMat    = [1., 0.]  # Crosstalk matrix indices [normalized]
    mrf = idealMRF(mrDef['fo'], mrDef['r'], mrDef['n'], mrDef['ng'], mrDef['alf'], mrDef['k'], crosstalkMat)
    mrf.clipping = False

    # Randomize it's response
    #mrf.manufacturing(pi/4)
    #outFieldsNV = mrf.TMM(wavelength, np.ones(len(wavelength)))
    #mrf.plotTransmission(wavelength, outFieldsNV)

    # Do some phase tuning
    #mrf.manufacturingValue([-pi/8, pi/8])
    #outFieldsNV = mrf.TMM(wavelength, np.ones(len(wavelength)))
    #mrf.plotTransmission(wavelength, outFieldsNV)

    # Quad 1 (max)
    #mrf.manufacturingValue([pi/8, pi/8])
    # Quad 2 (not max)
    #mrf.manufacturingValue([-pi/8, pi/8])
    # Quad 3 (max)
    #mrf.manufacturingValue([-pi/8, -pi/4])
    # Quad 4 (not max)
    #mrf.manufacturingValue([pi/8, -pi/8])


    # Random center but well defined delta 
    #maximumInterval = pi/2
    #delta = pi/2
    #center = random.uniform(-maximumInterval/2, maximumInterval/2)
    #mrf.manufacturingValue([center-delta, center+delta])

    # Phase Domain 
    phiMax = 2*pi
    phi1, phi2 = np.linspace(0, phiMax, 100), np.linspace(0, phiMax, 100)
    objFunc = mapObjectiveFunction(mrf, phi1, phi2, targetWavelength=1541.9e-9)
    plotObjectiveFunction(phi1, phi2, objFunc, varType='phase')

    # Power Domain
    Pmax = 2*pi
    P1, P2 = np.linspace(0, Pmax, 100), np.linspace(0, Pmax, 100)
    objFunc = mapObjectiveFunction(mrf, P1, P2, targetWavelength=1541.9e-9)
    plotObjectiveFunction(P1, P2, objFunc, varType='power')



"""
1b - Detuning V2 : Map the objective function 
"""
if False:
    # General parameters
    pdfFilename = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter4/detuning.pdf'
    wavelength=np.linspace(1540e-9, 1580e-9, 10000)

    # Build the dual ring filter from defaults
    crosstalkMat    = [1., 0.5]  # Crosstalk matrix indices [normalized]
    mrf = idealMRF(mrDef['fo'], mrDef['r'], mrDef['n'], mrDef['ng'], mrDef['alf'], mrDef['k'], crosstalkMat)
    mrf.clipping = False

    # Random center but well defined delta 
    maximumInterval, delta = pi/2, pi/2
    center = random.uniform(-maximumInterval/2, maximumInterval/2)
    mrf.manufacturingValue([center-delta, center+delta])

    # Results
    mapObjectiveFunction(mrf, targetWavelength=1541.9e-9)


"""
2 - Feasible region vs initial conditions for the random phase shift
"""
if False:
    ran = pi/4
    fig, ax = plt.subplots()
    p1 = np.linspace(0, 2*pi, 100)
    p2 = p1
    ax.plot(p1, p2)
    p3 = p1 + ran
    ax.plot(p1, p2)
    p4 = p1 - ran
    ax.plot(p1, p2)
    ax.fill_between(p1, p3, p4, color='lightgrey')
    plotFeasibleRegion(ax, p1Min=0, p1Max=2*pi, p2Min=0, p2Max=2*pi, coupling=0.5)
    t = np.linspace(0,pi/2, 100)
    radius = pi
    ax.plot(radius*np.cos(t), radius*np.sin(t))
    radius = 2*pi
    ax.plot(radius*np.cos(t), radius*np.sin(t))
    ax.set_xlim([0, 2*pi])
    ax.set_ylim([0, 2*pi])
    ax.set_aspect('equal', 'box')
    plt.show()

"""
3 - Detuning map as a function of the coupling for various values of coupling. Multiple plots as pdfs
"""
if False:
    # General parameters
    pdfDir = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter3/objectiveFunctions/'

    # Build the dual ring filter from defaults
    for k in [0., 0.2, 0.4, 0.6, 0.8, 1.0]:
        crosstalkMat    = [1., k]  # Crosstalk matrix indices [normalized]
        mrf = idealMRF(mrDef['fo'], mrDef['r'], mrDef['n'], mrDef['ng'], mrDef['alf'], mrDef['k'], crosstalkMat)
        mrf.clipping = False

        # Plot the detuning maps for
        Pmax = 2*pi
        P1, P2 = np.linspace(0, Pmax, 100), np.linspace(0, Pmax, 100)
        objFunc = mapObjectiveFunction(mrf, P1, P2, targetWavelength=1541.9e-9)
        plotObjectiveFunction(P1, P2, objFunc, pdfDir + '{0:.0f}.pdf'.format(k*100), varType='power')

"""
3 - Detuning map as a function of the coupling for various values of coupling. Single plot
"""
if False:
    pdfDir = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter3/objectiveFunctions/'

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10, 6))
    for k, ax, cap in zip([0., 0.2, 0.4, 0.6, 0.8, 1.0], axes.flat, ['a) $\mu=0.0$', 'b) $\mu=0.2$', 'c) $\mu=0.4$', 'd) $\mu=0.6$', 'e) $\mu=0.8$', 'f) $\mu=1.0$']):

        crosstalkMat    = [1., k]  # Crosstalk matrix indices [normalized]
        mrf = idealMRF(mrDef['fo'], mrDef['r'], mrDef['n'], mrDef['ng'], mrDef['alf'], mrDef['k'], crosstalkMat)
        mrf.clipping = False

        # Get the detuning map
        Pmax = 2*pi
        P1, P2 = np.linspace(0, Pmax, 100), np.linspace(0, Pmax, 100)
        objFunc = mapObjectiveFunction(mrf, P1, P2, targetWavelength=1541.9e-9)

        # Plot the detuning maps for
        X, Y = np.meshgrid(P1/pi, P2/pi, copy=False, indexing='ij')
        cs      = ax.contourf(X, Y, objFunc, vmin=-60., vmax=0.)
        ax.contour(cs, colors='k')
        ax.set_aspect("equal") # Fix the aspect ratio. Use 'equal' for scaled and 'auto' for x=y
        #cbar    = fig.colorbar(cs, ax=ax)            
        #cbar.ax.set_ylabel('Insertion loss [dB]', font=overpassFont)
        ax.set_xlabel(r'$\frac{P_1}{P_{\pi, 1}}$', font=overpassFont)
        ax.set_ylabel(r'$\frac{P_2}{P_{\pi, 2}}$', font=overpassFont)
        # Plot the constraints equations
        #plotFeasibleRegion(ax, p1Min=0, p1Max=2, p2Min=0, p2Max=2, coupling=mrf.phase_coupling_matrix[0,1])
        #plt.savefig(filename)
        #plt.show()
        fig.tight_layout()
        fig.subplots_adjust(right=0.9)
        cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(cs, cax=cbar_ax)
        cbar.ax.set_ylabel('Insertion loss [dB]', font=overpassFont)

        ax.set_title('{}'.format(cap))
    plt.savefig(pdfDir + 'multi.pdf')
    plt.show()

"""
3 - Detuning map with feasible region on top
"""
if False:
    # Création du frame de la figure
    fig = plt.figure(figsize=(10, 6))
    gs = fig.add_gridspec(3, 4, wspace=0.5, hspace=1.)
    ax1 = fig.add_subplot(gs[:2, 0:2])
    ax1.set_title(r'a) $T_D(P_1, P_2)$')
    ax2 = fig.add_subplot(gs[:2, 2:4])
    ax2.set_title(r'b) $T_D(P_1, P_2)$')
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.set_title(r'c) $\phi_1(P_1, P_2)$')
    ax4 = fig.add_subplot(gs[2, 1])
    ax4.set_title(r'd) $\phi_2(P_1, P_2)$')
    ax5 = fig.add_subplot(gs[2, 2])
    ax5.set_title(r'e) $\phi_1(P_1, P_2)$')
    ax6 = fig.add_subplot(gs[2, 3])
    ax6.set_title(r'f) $\phi_2(P_1, P_2)$')

    # Création du cas uncoupled 
    k=0
    crosstalkMat    = [1., k]  # Crosstalk matrix indices [normalized]
    mrf = idealMRF(mrDef['fo'], mrDef['r'], mrDef['n'], mrDef['ng'], mrDef['alf'], mrDef['k'], crosstalkMat)
    mrf.clipping = False

    # Get the detuning map
    Pmax = 2*pi
    P1, P2 = np.linspace(0, Pmax, 100), np.linspace(0, Pmax, 100)
    objFunc, PS1, PS2 = mapObjectiveFunction(mrf, P1, P2, targetWavelength=1541.9e-9)

    # a) Objective Function
    X, Y = np.meshgrid(P1/pi, P2/pi, copy=False, indexing='ij')
    cs      = ax1.contourf(X, Y, objFunc, vmin=-60., vmax=0.)
    ax1.contour(cs, colors='k')
    ax1.set_aspect("equal") # Fix the aspect ratio. Use 'equal' for scaled and 'auto' for x=y
    #cbar    = fig.colorbar(cs, ax=ax)            
    #cbar.ax.set_ylabel('Insertion loss [dB]', font=overpassFont)
    ax1.set_xlabel(r'$\frac{P_1}{P_{\pi, 1}}$', font=overpassFont)
    ax1.set_ylabel(r'$\frac{P_2}{P_{\pi, 2}}$', font=overpassFont)

    # b) Phase shift 1 vs Pow
    X, Y = np.meshgrid(P1/pi, P2/pi, copy=False, indexing='ij')
    cs      = ax3.contourf(X, Y, PS1)
    ax3.contour(cs, colors='k')
    ax3.set_aspect("equal") # Fix the aspect ratio. Use 'equal' for scaled and 'auto' for x=y
    #cbar    = fig.colorbar(cs, ax=ax)            
    #cbar.ax.set_ylabel('Insertion loss [dB]', font=overpassFont)
    ax3.set_xlabel(r'$\frac{P_1}{P_{\pi, 1}}$', font=overpassFont)
    ax3.set_ylabel(r'$\frac{P_2}{P_{\pi, 2}}$', font=overpassFont)

    # c) Phase shift 2 vs Pow
    X, Y = np.meshgrid(P1/pi, P2/pi, copy=False, indexing='ij')
    cs      = ax4.contourf(X, Y, PS2)
    ax4.contour(cs, colors='k')
    ax4.set_aspect("equal") # Fix the aspect ratio. Use 'equal' for scaled and 'auto' for x=y
    #cbar    = fig.colorbar(cs, ax=ax)            
    #cbar.ax.set_ylabel('Insertion loss [dB]', font=overpassFont)
    ax4.set_xlabel(r'$\frac{P_1}{P_{\pi, 1}}$', font=overpassFont)
    ax4.set_ylabel(r'$\frac{P_2}{P_{\pi, 2}}$', font=overpassFont)


    # Création du cas coupled
    k=0.8
    crosstalkMat    = [1., k]  # Crosstalk matrix indices [normalized]
    mrf = idealMRF(mrDef['fo'], mrDef['r'], mrDef['n'], mrDef['ng'], mrDef['alf'], mrDef['k'], crosstalkMat)
    mrf.clipping = False

    # Get the detuning map
    Pmax = 2*pi
    P1, P2 = np.linspace(0, Pmax, 100), np.linspace(0, Pmax, 100)
    objFunc, PS1, PS2 = mapObjectiveFunction(mrf, P1, P2, targetWavelength=1541.9e-9)

    # a) Objective Function
    X, Y = np.meshgrid(P1/pi, P2/pi, copy=False, indexing='ij')
    cs      = ax2.contourf(X, Y, objFunc, vmin=-60., vmax=0.)
    ax2.contour(cs, colors='k')
    ax2.set_aspect("equal") # Fix the aspect ratio. Use 'equal' for scaled and 'auto' for x=y
    #cbar    = fig.colorbar(cs, ax=ax)            
    #cbar.ax.set_ylabel('Insertion loss [dB]', font=overpassFont)
    ax2.set_xlabel(r'$\frac{P_1}{P_{\pi, 1}}$', font=overpassFont)
    ax2.set_ylabel(r'$\frac{P_2}{P_{\pi, 2}}$', font=overpassFont)

    # b) Phase shift 1 vs Pow
    X, Y = np.meshgrid(P1/pi, P2/pi, copy=False, indexing='ij')
    cs      = ax5.contourf(X, Y, PS1)
    ax5.contour(cs, colors='k')
    ax5.set_aspect("equal") # Fix the aspect ratio. Use 'equal' for scaled and 'auto' for x=y
    #cbar    = fig.colorbar(cs, ax=ax)            
    #cbar.ax.set_ylabel('Insertion loss [dB]', font=overpassFont)
    ax5.set_xlabel(r'$\frac{P_1}{P_{\pi, 1}}$', font=overpassFont)
    ax5.set_ylabel(r'$\frac{P_2}{P_{\pi, 2}}$', font=overpassFont)

    # c) Phase shift 2 vs Pow
    X, Y = np.meshgrid(P1/pi, P2/pi, copy=False, indexing='ij')
    cs      = ax6.contourf(X, Y, PS2)
    ax6.contour(cs, colors='k')
    ax6.set_aspect("equal") # Fix the aspect ratio. Use 'equal' for scaled and 'auto' for x=y
    #cbar    = fig.colorbar(cs, ax=ax)            
    #cbar.ax.set_ylabel('Insertion loss [dB]', font=overpassFont)
    ax6.set_xlabel(r'$\frac{P_1}{P_{\pi, 1}}$', font=overpassFont)
    ax6.set_ylabel(r'$\frac{P_2}{P_{\pi, 2}}$', font=overpassFont)

    fig.tight_layout()
    plt.show()

"""
3 - Feasible region in the phase domain and in the power domain
"""
if False:
    # General parameters
    pdfFilename = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter4/detuning.pdf'
    wavelength=np.linspace(1540e-9, 1580e-9, 10000)

    # Build the dual ring filter from defaults
    k = 0.5
    crosstalkMat    = [1., k]  # Crosstalk matrix indices [normalized]
    mrf = idealMRF(mrDef['fo'], mrDef['r'], mrDef['n'], mrDef['ng'], mrDef['alf'], mrDef['k'], crosstalkMat)
    mrf.clipping = False

    phiMax = 2*pi
    phi1, phi2 = np.linspace(0, phiMax, 100), np.linspace(0, phiMax, 100)
    objFunc, PS1, PS2 = mapObjectiveFunction(mrf, phi1, phi2, targetWavelength=1541.9e-9)

    # Plots
    fig, ax = plt.subplots(figsize=(10, 6))

    X, Y = np.meshgrid(phi1/pi, phi2/pi, copy=False, indexing='ij')
    cs      = ax.contourf(X, Y, objFunc, vmin=-60., vmax=0.)
    ax.contour(cs, colors='k')
    ax.set_aspect("equal") # Fix the aspect ratio. Use 'equal' for scaled and 'auto' for x=y
    #cbar    = fig.colorbar(cs, ax=ax)
    #cbar.ax.set_ylabel('Insertion loss [dB]', font=overpassFont)
    ax.set_xlabel(r'$\frac{P_1}{P_{\pi, 1}}$', labelpad=0, font=overpassFont)
    ax.set_ylabel(r'$\frac{P_2}{P_{\pi, 2}}$', labelpad=0, font=overpassFont)
    # Plot the constraints equations
    plotFeasibleRegion(ax, p1Min=0, p1Max=2, p2Min=0, p2Max=2, coupling=mrf.phase_coupling_matrix[0, 1])
    plt.savefig(pdfFilename)
    plt.show()

