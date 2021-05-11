"""
Functions used to plot and analyze detuning maps aka objective function of the optimization problem

Author: Simon Bélanger-de Villers
Date : March 14th 2021
"""
import numpy as np
import matplotlib.pyplot as plt
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB, field2PowerLinear
from math import pi
import matplotlib.font_manager as fm

overpassFont = fm.FontProperties(family = 'Overpass', 
                                fname = '/Users/simonbelanger/Library/Fonts/overpass-semibold.otf', 
                                size = 12)

def mapObjectiveFunction(mrf, x1, x2, targetWavelength=1550-9):
    """ Measure and map the objective function of a 2 rings mrf for an input set of variables x1 and x2.

    """

    # TODO :Measure the objective function (in terms of power)
    # TODO :Measure the objective function (in terms of voltage)

    # Initialize the arrays
    objFunc         = np.zeros((len(x1), len(x2))) 
    phaseShift1     = np.zeros((len(x1), len(x2)))
    phaseShift2     = np.zeros((len(x1), len(x2)))

    # Measure the objective function (in terms of phase)
    for i in range(len(x1)):
        for j in range(len(x2)):
            mrf.applyPhaseTuningCrosstalk([x1[i], x2[j]])
            outFields = mrf.TMM(targetWavelength, 1.0)
            objFunc[i,j]        = field2PowerdB(outFields.drop)
            phaseShift1[i,j]    = mrf.actual_tuning_phase[0]
            phaseShift2[i,j]    = mrf.actual_tuning_phase[1]
    return objFunc, phaseShift1, phaseShift2
    # Plot type #2 : IL vs P1 & P2, normalized to Ppi


    #plt.contourf(x, y, phaseDetuningMapP1)
    #plt.show()

    #plt.contourf(x, y, phaseDetuningMapP2)
    #plt.show()

    #plt.contourf(phaseDetuningMapP1/pi, phaseDetuningMapP2/pi, phaseDetuningMap, vmin=-60., vmax=0.)
    #plt.show()

    # Plot in phase coordinates
    #fig, ax =plt.subplots()
    #ax.contourf(phaseDetuningMapP1/pi, phaseDetuningMapP2/pi, phaseDetuningMap, vmin=-60., vmax=0.)
    #plotFeasibleRegion(ax, p1Min=0, p1Max=2, p2Min=0, p2Max=2, coupling=mrf.phase_coupling_matrix[0,1])
    #plt.show()

    # Misc
    #plt.contourf(x, y, phaseDetuningMapP1)
    #plt.show()
    #plt.contourf(x, y, phaseDetuningMapP2)
    #plt.show()

def plotObjectiveFunction(x1, x2, objFunc, filename, varType):
    """[summary]

    Args:
        X ([type]): [description]
        Y ([type]): [description]
        objFunc ([type]): [description]
    """

    # Phase domain
    if varType=='phase':
        X, Y = np.meshgrid(x1/pi, x2/pi, copy=False, indexing='ij')
        fig, ax = plt.subplots()
        cs      = ax.contourf(X, Y, objFunc, vmin=-60., vmax=0.)
        cbar    = fig.colorbar(cs, ax=ax)
        ax.contour(cs, colors='k')
        ax.set_aspect("equal") # Fix the aspect ratio. Use 'equal' for scaled and 'auto' for x=y
        cbar.ax.set_ylabel('Insertion loss [dB]', font=overpassFont)
        ax.set_xlabel(r'$\Delta\phi_1$ [$\pi$]', labelpad=0, font=overpassFont)
        ax.set_ylabel(r'$\Delta\phi_2$ [$\pi$]', labelpad=0, font=overpassFont)
        # Plot the constraints equations
        #plotFeasibleRegion(ax, p1Min=0, p1Max=2, p2Min=0, p2Max=2, coupling=mrf.phase_coupling_matrix[0,1])
        plt.savefig(filename)
        plt.show()

    # Phase domain
    if varType=='power':
        X, Y = np.meshgrid(x1/pi, x2/pi, copy=False, indexing='ij')
        fig, ax = plt.subplots()
        cs      = ax.contourf(X, Y, objFunc, vmin=-60., vmax=0.)
        ax.contour(cs, colors='k')
        ax.set_aspect("equal") # Fix the aspect ratio. Use 'equal' for scaled and 'auto' for x=y
        #cbar    = fig.colorbar(cs, ax=ax)
        #cbar.ax.set_ylabel('Insertion loss [dB]', font=overpassFont)
        ax.set_xlabel(r'$\frac{P_1}{P_{\pi, 1}}$', labelpad=0, font=overpassFont)
        ax.set_ylabel(r'$\frac{P_2}{P_{\pi, 2}}$', labelpad=0, font=overpassFont)
        # Plot the constraints equations
        #plotFeasibleRegion(ax, p1Min=0, p1Max=2, p2Min=0, p2Max=2, coupling=mrf.phase_coupling_matrix[0,1])
        plt.savefig(filename)
        plt.show()

def plotFeasibleRegion(ax, p1Min=0, p1Max=2*pi, p2Min=0, p2Max=2*pi, coupling=0.5):
    """
        The feasible region is defined from the constraints imposed on the outputs of the DC sources.
    """
    k = coupling
    try:
        x = np.linspace(-p2Max, 2*p2Max)
        # Constraint functions
        y = k*x
        ax.plot(x, y, linestyle='--', color='black')
        y = 1/k * x
        ax.plot(x, y, linestyle='--', color='black')
        y = k*x + p2Max * (1-k**2)
        ax.plot(x, y, linestyle='--', color='black')
        y = 1/k*x + p1Max * (k**2-1)/k
        ax.plot(x, y, linestyle='--', color='black')

        # Constraint intersections
        p1 = [0, p2Max*k, p1Max, p1Max + p2Max*k]
        p2 = [0, p1Max, p2Max*k, p2Max + p1Max*k]
        ax.plot(p1, p2,'o', color='black')

        # Region
        p1 = [0, p2Max*k, p1Max + p2Max*k, p1Max, 0]
        p2 = [0, p1Max, p2Max + p1Max*k, p2Max*k, 0]
        ax.fill(p1, p2,'black', alpha=0.5)

    except ZeroDivisionError:
        # Constraint functions
        ax.hlines(p2Min, -100, 100, linestyle='--', color='black')
        ax.hlines(p2Max, -100, 100, linestyle='--', color='black')
        ax.vlines(p1Min, -100, 100, linestyle='--', color='black')
        ax.vlines(p1Max, -100, 100, linestyle='--', color='black')

        # Constraint intersections
        ax.plot([p1Min, p1Max, p1Max, p1Min],[p2Min, p2Min, p2Max, p2Max],'o', color='black')

        # Region
        ax.fill([p1Min, p1Max, p1Max, p1Min, p1Min],[p2Min, p2Min, p2Max, p2Max, p2Min],'black', alpha=0.5)

    ax.set_xlim([-1, 1.1*(p1Max + p2Max*k)])
    ax.set_ylim([-1, 1.1*(p1Max + p2Max*k)])