import numpy as np 
import matplotlib.pyplot as plt
import math, sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import *
from components.MRF import idealMRF
import matplotlib.font_manager as fm
"""
Script used to compare the response of rings vs the order of the filter
Thesis
"""

pdfFilename = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter2/ringOrder.pdf'

wavelength  = np.linspace(1538e-9, 1546e-9, 1000)

overpassFont = fm.FontProperties(family = 'Overpass', 
                                fname = '/Users/simonbelanger/Library/Fonts/overpass-semibold.otf', 
                                size = 12)

overpassFontLeg = fm.FontProperties(family = 'Overpass', 
                                fname = '/Users/simonbelanger/Library/Fonts/overpass-semibold.otf', 
                                size = 10)

# Filter parameters
radius 		    = 2.5e-6 	# Radius of the ring resonators [m]
neff 		    = 2.4449 	# Effective index of the ring waveguide []
ng 			    = 4.18 		# Group index of the ring waveguide []
alpha_wg	    = 3. 		# Losses of the ring waveguide [dB/cm]
powerCoupling   = 0.3       # PowerCoupling of the first ring resonator to the bus 

fig, ax = plt.subplots()
for filterOrder, color, annotationPosition,labelExp in zip([1,2,3,4,5], ['b', 'g','r','c','m'], [-11, -23, -42, -65, -86.5], ['st', 'nd', 'rd', 'th', 'th']):
    mrf = idealMRF(filterOrder, radius, neff, ng, alpha_wg, powerCoupling, [1])
    mrf.clipping = False
    outFields = mrf.TMM(wavelength, np.ones(len(wavelength)))
    ax.plot(wavelength * 1e9, 10 * np.log10(abs(outFields.drop) ** 2), color=color,linewidth=2)
    ax.plot(wavelength * 1e9, 10 * np.log10(abs(outFields.through) ** 2), color=color, linewidth=2, linestyle='dashed')
    ax.annotate(str(filterOrder)+'$^{'+labelExp+'}$ Order', [1544.5, annotationPosition], font=overpassFontLeg,color=color) #backgroundcolor='w'
ax.set_xlabel('Wavelength (nm)', font=overpassFont)
ax.set_ylabel('Power transmission (dB)', font=overpassFont)
ax.grid()
fig.savefig(pdfFilename)
plt.show()

