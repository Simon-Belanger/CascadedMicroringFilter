"""
    fdtdDirectionalCouplersAnalysis.py
This script will be used to plot all the figures regarding the directional couplers

Plan
 - 5 µm bent coupler : kappa vs gap
 - 2.5 µm bent coupler : kappa vs gap
 - 5 µm interring coupler : kappa vs gap
 - 2.5 µm interring coupler : kappa vs gap

Simon Belanger-de Villers
 December 3rd 2020
"""

import numpy as np
import pickle
import os, sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from components.couplers import dataHandling
import matplotlib.pyplot as plt
from scipy import interpolate

# Plot the transmission of the coupler vs the coupler gap
def plotTransmissionVsGap(data, targetWavelength=1550e-9, options={'log':False}):
    'Plot the coupling vs gap for the given target wavelength'
    kappa = []; tr = []; total = []
    for index in range(0, len(data['gap [m]'])):
        f = interpolate.interp1d(data['scattering'][index]['wavelength'], np.abs(data['scattering'][index]['S_drop'])**2)
        kappa.append(f(targetWavelength))
        f2 = interpolate.interp1d(data['scattering'][index]['wavelength'], np.abs(data['scattering'][index]['S_through'])**2)
        tr.append(f2(targetWavelength))
        total.append(f(targetWavelength) + f2(targetWavelength))
    if options['log'] ==True:
        kappa, tr  = 10*np.log(kappa), 10*np.log(tr)
    plt.plot(np.asarray(data['gap [m]'])*1e9, kappa, marker='o', color='lightcoral', markersize=10, fillstyle='none', markeredgewidth = 3, linewidth=3, label='Cross-coupling')
    plt.plot(np.asarray(data['gap [m]'])*1e9, tr,marker='s', color='seagreen', markersize=10, fillstyle='none', markeredgewidth = 3, linewidth=3, label='Thru-Coupling')
    #plt.plot(np.asarray(data['gap [m]'])*1e9, total, color='blue', fillstyle='none', markeredgewidth = 3, linewidth=3, label='Total')
    plt.xlabel('Gap [nm]', fontsize=16); plt.ylabel('Power coupling coefficient', fontsize=16);plt.grid();#plt.legend()
    plt.xticks(fontsize=16);plt.yticks(fontsize=16)
    #plt.axis([min(data['gap [m]'])*1e9, max(data['gap [m]'])*1e9, 0, 1])
    #plt.show()

def plotTransmissionVsGapInter(data, targetWavelength=1550e-9, options={'log':False}):
    'Plot the coupling vs gap for the given target wavelength'
    kappa = []; tr = []; total = []
    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
    for index in range(0, len(data['gap [m]'])):
        f = interpolate.interp1d(data['scattering'][index]['wavelength'], np.abs(data['scattering'][index]['S_drop'])**2)
        kappa.append(f(targetWavelength))
        f2 = interpolate.interp1d(data['scattering'][index]['wavelength'], np.abs(data['scattering'][index]['S_through'])**2)
        tr.append(f2(targetWavelength))
        total.append(f(targetWavelength) + f2(targetWavelength))
    if options['log'] ==True:
        kappa, tr  = 10*np.log(kappa), 10*np.log(tr)
    ax2.plot(np.asarray(data['gap [m]'])*1e9, kappa, marker='o', color='lightcoral', markersize=10, fillstyle='none', markeredgewidth = 3, linewidth=3, label='Cross-coupling')
    ax1.plot(np.asarray(data['gap [m]'])*1e9, tr,marker='s', color='seagreen', markersize=10, fillstyle='none', markeredgewidth = 3, linewidth=3, label='Thru-Coupling')
    #plt.plot(np.asarray(data['gap [m]'])*1e9, total, color='blue', fillstyle='none', markeredgewidth = 3, linewidth=3, label='Total')
    plt.xlabel('Gap [nm]', fontsize=16);ax1.grid();ax2.grid();#plt.legend()
    fig.text(0.02, 0.5, 'Power coupling coefficient', fontsize=16, va='center', rotation='vertical')
    xticks = np.linspace(180, 360, 5)
    #plt.xticks(xticks)
    #ax2.set_xticklabels(xticks, fontsize=16)
    ax2.set_yticks(np.linspace(0.00, 0.02, 5))
    ax1.set_yticks(np.linspace(0.95, 1, 5))
    plt.setp(ax2.get_xticklabels(), Fontsize=16)
    plt.setp(ax1.get_yticklabels(), Fontsize=16)
    plt.setp(ax2.get_yticklabels(), Fontsize=16)
    #plt.xticks(xticks)
    #ax2.set_xticklabels(xticks, fontsize=16)
    #plt.axis([min(data['gap [m]'])*1e9, max(data['gap [m]'])*1e9, 0, 1])
    #plt.show()


    
# Location to save pdfs
dataDirectory   = os.getcwd() + "/thesis/chap2/data.nosync/"
pdfDirectory 	= '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter2/couplerData/'

"""
Graph #1 : transmission for 2.5 um bent coupler
"""
if False:
    dat = pickle.load(open(dataDirectory + "DC_Bent_R2_5.pickle", "rb") )
    plotTransmissionVsGap(dat, 1550e-9, {'log': False})
    #plt.savefig(pdfDirectory + 'DC_Bent_R2_5.pdf')
    plt.show()

"""
Graph #2 : transmission for 5 um bent coupler
"""
if False:
    dat = pickle.load(open(dataDirectory + "DC_Bent_R5.pickle", "rb" ))
    plotTransmissionVsGap(dat, 1550e-9, {'log': False})
    #plt.savefig(pdfDirectory + 'DC_Bent_R5.pdf')
    plt.show()

"""
Graph #3 : transmission for 2.5 um interring coupler
"""
if False:
    dat = pickle.load(open(dataDirectory + "DC_seriesrings_R2_5.pickle", "rb") )
    plotTransmissionVsGapInter(dat, 1550e-9, {'log': False})
    #plt.savefig(pdfDirectory + 'DC_seriesrings_R2_5.pdf')
    plt.show()

"""
Graph #4 : transmission for 5 um interring coupler (hardcoded parce que données erronées)
"""
if True:
    from scipy.io import loadmat

    data = loadmat(dataDirectory + "DC_seriesring_R5_V2.mat")
    g = np.squeeze(data['gap_interRing_1550'])
    kappa = np.squeeze(data['kappa_interRing_1550'])
    tr = np.squeeze(data['t_interRing_1550'])

    fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
    ax2.plot(np.asarray(g), kappa, marker='o', color='lightcoral', markersize=8, fillstyle='none', markeredgewidth = 2, linewidth=2, label='Cross-coupling')
    ax1.plot(np.asarray(g), tr, marker='s', color='seagreen', markersize=8, fillstyle='none', markeredgewidth = 2, linewidth=2, label='Thru-Coupling')
    plt.xlabel('Gap [nm]', fontsize=12);ax1.grid();ax2.grid();#plt.legend()
    fig.text(0.02, 0.5, 'Power coupling coefficient', fontsize=14, va='center', rotation='vertical')
    xticks = np.linspace(180, 360, 5)
    ax2.set_xlim([180, 300])
    ax1.set_ylim([0.97, 1])
    ax2.set_ylim([0, 0.03])
    #plt.xticks(xticks)
    #ax2.set_xticklabels(xticks, fontsize=16)
    ax2.set_yticks(np.linspace(0.00, 0.03, 5))
    ax1.set_yticks(np.linspace(0.97, 1, 5))
    plt.setp(ax2.get_xticklabels(), Fontsize=12)
    #plt.setp(ax1.get_yticklabels(), Fontsize=16)
    #plt.setp(ax2.get_yticklabels(), Fontsize=16)
    #plt.xticks(xticks)
    #ax2.set_xticklabels(xticks, fontsize=16)
    #plt.axis([min(data['gap [m]'])*1e9, max(data['gap [m]'])*1e9, 0, 1])
    plt.savefig(pdfDirectory + 'DC_seriesrings_R5_V2.pdf')
    plt.show()

# Old and obsolete
if False:
    dat = pickle.load(open(dataDirectory + "DC_seriesrings_R5.pickle", "rb" ))
    plotTransmissionVsGapInter(dat, 1550e-9, {'log': False})
    #plt.savefig(pdfDirectory + 'DC_seriesrings_R5.pdf')
    plt.show()

"""
Curve fitting of the coupling coefficient vs gap for all 4 cases
"""
if False:
    dat = pickle.load(open(dataDirectory + "DC_seriesrings_R5.pickle", "rb" ))
    kappa = []; tr = []; total = []
    for index in range(0, len(dat['gap [m]'])):
        f = interpolate.interp1d(dat['scattering'][index]['wavelength'], np.abs(dat['scattering'][index]['S_drop']))
        kappa.append(f(1550e-9))
        f2 = interpolate.interp1d(dat['scattering'][index]['wavelength'], np.abs(dat['scattering'][index]['S_through']))
        tr.append(f2(1550e-9))
    plt.plot(np.asarray(dat['gap [m]'])*1e9, np.asarray(kappa)**2)
    plt.show()