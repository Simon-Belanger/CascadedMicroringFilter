"""
Script used to generate figures for the thesis.

Author      : Simon Belanger-de Villers (simon.belanger-de-villers.1@ulaval.ca)
Created     : August 13th 2020
Last edited : October 23rd 2020
"""
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm 
import math, sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB, field2PowerLinear
from misc.Filter import ringAnalyisGeneral
from misc import constants
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import *
from components.MRF import MRF, idealMRF
import matplotlib.font_manager as fm

overpassFont = fm.FontProperties(family = 'Overpass', 
                                fname = '/Users/simonbelanger/Library/Fonts/overpass-semibold.otf', 
                                size = 12)

# TODO : Compare TMM with analytical model from J. Poon and J. C.C. Mak
# TODO : plot SBSR vs coupling
# TODO : Plot ripples amplitude vs coupling (or max loss in band)
# TODO : Confirm roll off does not depend on coupling, only on number of rings
# TODO : Make the model vectorial
# TODO : Implement Ripples_T: Ripples amplitude [dB]
# TODO : Implement Ripples_D : Ripples amplitude [dB]
# TODO : Implement Roll-off (dB/decade) to implement later : wavelength to frequency, plot and fit 

def plotSpectralResponseVsCouplingCoefficient(wavelength, couplingValues, pdfFilename, plotType='nm'):
	"""
	Response vs coupler coupling coefficient for 5 rings networks
	"""

	cmap = cm.get_cmap('jet')
	for C in couplingValues:
		mrf = idealMRF(5, 2.5e-6, 2.4449, 4.18, 3., C, [1])
		outFields = mrf.TMM(wavelength, np.ones(len(wavelength)))
		mrf.clipping = False
		outFields = mrf.TMM(wavelength, np.ones(len(wavelength)))
		lineColor = cmap(C/max(couplingValues))
		if plotType == 'nm':
			xAxis = wavelength * 1e9
		elif plotType == 'GHz':
			xAxis = constants.c / wavelength * 1e-9
		plt.plot(xAxis, 10 * np.log10(abs(outFields.drop) ** 2), linewidth=2, label='{}'.format(C), color=lineColor)
		plt.plot(xAxis, 10 * np.log10(abs(outFields.through) ** 2), linewidth=2, linestyle='dashdot', color=lineColor)

	if plotType == 'nm':
		plt.xlabel('Wavelength (nm)', fontsize=14)
	elif plotType == 'GHz':
		plt.xlabel('Frequency (GHz)', fontsize=14)
	plt.ylabel('Power transmission (dB)', fontsize=14)
	plt.grid()
	#plt.ylim([-80, 0])
	plt.legend(loc='upper right')
	plt.savefig(pdfFilename)
	plt.show()

def measureBandwidthVsCouplingCoefficient(radius, wavelength, order, couplingValues):
	' Response vs coupler coupling coefficient for 5 rings networks. ' 

	bandwidth = []
	for C in couplingValues:
		outFields = idealMRF(order, radius, 2.4449, 4.18, 3., C, [1]).TMM(wavelength, np.ones(len(wavelength)))
		anal = ringAnalyisGeneral(wavelength, field2PowerLinear(outFields.drop), field2PowerLinear(outFields.through))      # Normalise the data using ringAnalysis (findpeaks)
		try:
			bandwidth.append(anal.measureBandwidth()[0]*1e9)
		except IndexError:
			bandwidth.append(0)

	return couplingValues, bandwidth


"""
1) Response vs coupler coupling coefficient for 5 rings networks (fsr)
"""
if False:
    #Script for respvskappa
    wavelength  	= np.linspace(1500e-9, 1600e-9, 1000)
    couplingValues 	= np.linspace(0.1, 1, 2)
    pdfFilename 	= '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter2/respvskappa.pdf'
    plotType 		= 'nm' # nm or GHz
    plotSpectralResponseVsCouplingCoefficient(wavelength, couplingValues, pdfFilename, plotType)

"""
2) Response vs coupler coupling coefficient for 5 rings networks (peak)
"""
if False:
    #Script for respvskappa_zoom
    wavelength  	= np.linspace(1535e-9, 1545e-9, 1000)
    couplingValues 	= np.linspace(0.1, 1, 2)
    pdfFilename 	= '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter2/respvskappa_zoom.pdf'
    plotType 		= 'GHz' # nm or GHz
    plotSpectralResponseVsCouplingCoefficient(wavelength, couplingValues, pdfFilename, plotType)

"""
3) Plot filter 3dB bandwidth vs coupling coefficient for different number of rings, R=2.5 microns
"""

if False:
    pdfFilename 	= '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter2/bandwidthVkappa.pdf'

    wavelength  	= np.linspace(1500e-9, 1600e-9, 1000)
    couplingValues 	= np.linspace(0.001, 0.2, 10)

    # Get the results
    C1, B1 = measureBandwidthVsCouplingCoefficient(2.5e-6, wavelength, 1, couplingValues)
    C2, B2 = measureBandwidthVsCouplingCoefficient(2.5e-6, wavelength, 2, couplingValues)
    C3, B3 = measureBandwidthVsCouplingCoefficient(2.5e-6, wavelength, 3, couplingValues)
    C4, B4 = measureBandwidthVsCouplingCoefficient(2.5e-6, wavelength, 4, couplingValues)
    C5, B5 = measureBandwidthVsCouplingCoefficient(2.5e-6,wavelength, 5, couplingValues)

    # Plot the figure
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(C1, np.asarray(B1)*100/0.8, marker='s', label='Order 1')
    ax.plot(C2, np.asarray(B2)*100/0.8, marker='s', label='Order 2')
    ax.plot(C3, np.asarray(B3)*100/0.8, marker='s', label='Order 3')
    ax.plot(C4, np.asarray(B4)*100/0.8, marker='s', label='Order 4')
    ax.plot(C5, np.asarray(B5)*100/0.8, marker='s', label='Order 5')
    ax.set_xlabel(r'Power coupling coefficient, $\kappa^2$ []', font=overpassFont)
    ax.set_ylabel(r'3dB bandwidth [GHz]', font=overpassFont)
    ax.grid(color='grey', linestyle='--', linewidth=0.25)
    ax.legend(title='Width [nm]', prop=overpassFont, frameon=True, loc='upper left')
    plt.savefig(pdfFilename)
    plt.show()

"""
3) Plot filter 3dB bandwidth vs coupling coefficient for different number of rings, R=5 microns
"""

if True:
    pdfFilename 	= '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter2/bandwidthVkappa5.pdf'

    wavelength  	= np.linspace(1500e-9, 1600e-9, 1000)
    couplingValues 	= np.linspace(0.001, 0.2, 10)

    # Get the results
    C1, B1 = measureBandwidthVsCouplingCoefficient(5e-6, wavelength, 1, couplingValues)
    C2, B2 = measureBandwidthVsCouplingCoefficient(5e-6, wavelength, 2, couplingValues)
    C3, B3 = measureBandwidthVsCouplingCoefficient(5e-6, wavelength, 3, couplingValues)
    C4, B4 = measureBandwidthVsCouplingCoefficient(5e-6, wavelength, 4, couplingValues)
    C5, B5 = measureBandwidthVsCouplingCoefficient(5e-6,wavelength, 5, couplingValues)

    # Plot the figure
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(C1, np.asarray(B1)*100/0.8, marker='s', label='Order 1')
    ax.plot(C2, np.asarray(B2)*100/0.8, marker='s', label='Order 2')
    ax.plot(C3, np.asarray(B3)*100/0.8, marker='s', label='Order 3')
    ax.plot(C4, np.asarray(B4)*100/0.8, marker='s', label='Order 4')
    ax.plot(C5, np.asarray(B5)*100/0.8, marker='s', label='Order 5')
    ax.set_xlabel(r'Power coupling coefficient, $\kappa^2$ []', font=overpassFont)
    ax.set_ylabel(r'3dB bandwidth [GHz]', font=overpassFont)
    ax.grid(color='grey', linestyle='--', linewidth=0.25)
    ax.legend(title='Width [nm]', prop=overpassFont, frameon=True, loc='upper left')
    plt.savefig(pdfFilename)
    plt.show()
