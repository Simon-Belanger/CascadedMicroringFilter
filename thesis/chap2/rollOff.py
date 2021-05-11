import numpy as np 
import matplotlib.pyplot as plt
import math, sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB, field2PowerLinear
from misc.Filter import *
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import *
from components.MRF import idealMRF
import matplotlib.font_manager as fm

"""
Script used to compare the response of rings vs the order of the filter
Thesis
"""

pdfFilename = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter2/rollOff.pdf'

wavelength  = np.linspace(1530e-9, 1575e-9, 10000)

#font = fm.FontProperties(family = 'Overpass', fname = '/Users/simonbelanger/Library/Fonts/overpass-semibold.otf')

# Filter parameters
radius 		    = 2.5e-6 	# Radius of the ring resonators [m]
neff 		    = 2.4449 	# Effective index of the ring waveguide []
ng 			    = 4.18 		# Group index of the ring waveguide []
alpha_wg	    = 3. 		# Losses of the ring waveguide [dB/cm]
powerCoupling   = 0.1       # PowerCoupling of the first ring resonator to the bus 

fig, ax = plt.subplots()
for filterOrder, color, annotationPosition, theo in zip([1,2,3,4,5], ['b', 'g','r','c','m'], [-10,-30,-50,-70,-80], [20,40,60,80,100]):
    mrf = idealMRF(filterOrder, radius, neff, ng, alpha_wg, powerCoupling, [1])
    outFields = mrf.TMM(wavelength, np.ones(len(wavelength)))
    anal = ringAnalyis(wavelength, outFields)
    anal.displayProps()
    nf, p = anal.measureRollOff()
    #ax.plot(np.log10(nf), p, color=color,linewidth=2.5)
    ax.plot(nf, p, color=color,linewidth=2.5)
    #ax.annotate('{} dB/decade'.format(theo), [1e11, annotationPosition], fontsize=20, color=color)
ax.set_xscale('log')
ax.set_xlabel('f - f_0 (GHz)', fontsize=14)
ax.set_ylabel('Power transmission (dB)', fontsize=14)
ax.set_xlim(1e9, 1e12)
ax.set_ylim(-60, 0)
ax.grid()
fig.savefig(pdfFilename)
plt.show()

# 1 check if roll-off is accurate for 1 rings
# 2 plot the roll-off for 5 rings
# Verify if the roll off is the same for different values of kappa, R, neff etc