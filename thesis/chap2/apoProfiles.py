"""
Apodization profiles for the coupling of CROWs

Simon Bélanger-de Villers
April 2nd 2021

"""

import sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from components.couplers.apodization import *
import numpy as np

pdfDirectory = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter4/'

# Fig 1
#fig, ax = plt.subplots()
#displayApodizationCurves_couplingFixed(0.1, ax)
#plt.savefig(pdfDirectory + 'couplingFixed.pdf')
#plt.show()

"""
Default
"""
if False:
    fig, ax = plt.subplots()
    displayApodizationCurves_orderFixed(5, np.linspace(0.01, 1, 5), ax)
    plt.savefig(pdfDirectory + 'orderFixed.pdf')
    plt.show()

"""
1 - Profile d'apodization optimal selon B.E. Little vs ce qui a été fait pour R=5um
"""
if True:

    kExp    = [0.552, 0.123, 0.106, 0.106, 0.123, 0.552]        # Kappa expérimental
    kLilPow = flatTopCoupling(5, kExp[0]**4)         # Kappa expérimental
    kLil    = flatTopCoupling(5, kExp[0]**2)                    # Kappa selon B.E. Little
    kExpA, kLilA = np.asarray(kExp), np.asarray(kLil)

    # Plot
    fig, axes = plt.subplots(nrows=3, sharex=True)

    ## 2 courbes
    axes[0].plot(range(1,7),kExpA, color='k', linestyle='solid', label="Exp.")
    axes[0].plot(range(1,7),kLilA, color='k', linestyle='dashed', label="B.E. Little")
    axes[0].plot(range(1,7),kLilPow, color='k', linestyle='dashed', label="B.E. Little pow")
    ## différence, erreur absolue
    axes[1].plot(range(1,7), kExpA-kLilA, color='k', linestyle='solid', label="B.E. Little")
    ## différence, erreur relative
    axes[2].plot(range(1,7), (kExpA/kLilA-1)*100, color='k', linestyle='solid', label="B.E. Little")


    #ax.set_aspect(3)
    axes[2].set_xlabel(r'Coupler position, $i$')
    axes[0].set_ylabel(r'Field coupling coefficient, $\kappa_i$')
    axes[1].set_ylabel(r'Absolute error, ')
    axes[2].set_ylabel(r'Relative error, [%]')
    axes[0].set_xlim([1,6])
    axes[0].legend()
    plt.savefig(pdfDirectory+'apodization5.pdf')
    plt.show()

