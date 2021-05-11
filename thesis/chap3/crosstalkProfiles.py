import numpy as np 
import matplotlib.pyplot as plt
import math, sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.Crosstalk import *
import matplotlib.font_manager as fm

pdfFilename = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter3/crosstalkExp.pdf'

overpassFont = fm.FontProperties(family = 'Overpass', 
                                fname = '/Users/simonbelanger/Library/Fonts/overpass-semibold.otf', 
                                size = 12)

"""
Graph 1 : Plot the values of eta for different values of the coupling strength
"""
if True:

    # Plot the profile
    num_rings           = 5
    coupling_strength   = 1
    x = np.linspace(0, num_rings - 1, num_rings)
    fig, ax = plt.subplots()
    for xi in [0,0.3, 0.5, 1, 2, 5, 10, 100]:
        coupling_strength   = xi
        coeff = exponential_crosstalk(num_rings, coupling_strength)
        ax.plot(x, coeff, marker='s',label='{:.1f}'.format(coupling_strength))
        ax.set_xticks(x)
    ax.legend()
    ax.set_xlabel('Index, $i$', labelpad=0, font=overpassFont)
    ax.set_ylabel('Thermal crosstalk coefficient, $\eta_i$', labelpad=0, font=overpassFont)
    plt.savefig(pdfFilename)
    plt.show()