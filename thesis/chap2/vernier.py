"""
Compare a Vernier filter to a non-Vernier filter. Show the extinction ratio that is bad.
"""

import numpy as np 
import matplotlib.pyplot as plt
import math, sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB, field2PowerLinear
from misc.Filter import *
from components.couplers.virtual_DC import virtual_DC
from components.couplers.apodization import *
from components.MRF import parentMRF
from components.ring import Ring
import matplotlib.font_manager as fm

# General parameters
pdfFilename = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter2/vernier/vernier.pdf'
wavelength=np.linspace(1540e-9, 1580e-9, 10000)

overpassFont = fm.FontProperties(family = 'Overpass', 
                                fname = '/Users/simonbelanger/Library/Fonts/overpass-semibold.otf', 
                                size = 12)

"""
Case 1 : 2 rings notVernier filter
"""
num_rings 	= 2 		# Number of cascaded rings [-]
radius 		= 2.5e-6 	# Radius of the ring resonators [m]
neff 		= 2.4449 	# Effective index of the ring waveguide []
ng 			= 4.18 		# Group index of the ring waveguide []
alpha_wg	= 3. 		# Losses of the ring waveguide [dB/cm]

rings    = [Ring(radius, neff, ng, alpha_wg) for i in range(num_rings)]
print(rings[0].freeSpectralRange*1e9)
print(rings[1].freeSpectralRange*1e9)
k_vec = flatTopCoupling(2, 0.01)
loss_coupler = [1, 1, 1]
couplers = [virtual_DC(k_vec[0],loss_coupler[0]),
			virtual_DC(k_vec[1],loss_coupler[1]),
			virtual_DC(k_vec[2],loss_coupler[2])]
mrf = parentMRF(rings, couplers, [1, 0, 0])
mrf.clipping = False
outFieldsNV = mrf.TMM(wavelength, np.ones(len(wavelength)))

"""
Case 2: 2 rings Vernier filter
"""
index = [4, 5]
radius 		= [2.5e-6*index[0], 2.5e-6*index[1]] 	# Radius of the ring resonators [m]
neff 		= 2.4449 	# Effective index of the ring waveguide []
ng 			= 4.18 		# Group index of the ring waveguide []
alpha_wg	= 3. 		# Losses of the ring waveguide [dB/cm]
rings    = [Ring(radius[0], neff, ng, alpha_wg),
            Ring(radius[1], neff, ng, alpha_wg)]
print(rings[0].freeSpectralRange*1e9)
print(rings[1].freeSpectralRange*1e9)
k_vec = flatTopCoupling(2, 0.01)
loss_coupler = [1, 1, 1]
couplers = [virtual_DC(k_vec[0],loss_coupler[0]),
			virtual_DC(k_vec[1],loss_coupler[1]),
			virtual_DC(k_vec[2],loss_coupler[2])]
mrf = parentMRF(rings, couplers, [1, 0, 0])
mrf.clipping = False
outFieldsV = mrf.TMM(wavelength, np.ones(len(wavelength)))

""" 
Plot everything
"""
colors = ['tab:blue', 'tab:red']
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(wavelength * 1e9, 10 * np.log10(abs(outFieldsNV.drop) ** 2), color=colors[0],linewidth=2)
ax.plot(wavelength * 1e9, 10 * np.log10(abs(outFieldsNV.through) ** 2), color=colors[0], linewidth=2, linestyle='dashed')
ax.plot(wavelength * 1e9, 10 * np.log10(abs(outFieldsV.drop) ** 2), color=colors[1],linewidth=2)
ax.plot(wavelength * 1e9, 10 * np.log10(abs(outFieldsV.through) ** 2), color=colors[1], linewidth=2, linestyle='dashed')
ax.set_xlabel('Wavelength, $\lambda$ [nm]', font=overpassFont)
ax.set_ylabel('Power transmission [dB]', font=overpassFont)
ax.grid(color='grey', linestyle='--', linewidth=0.25)
ax.set_xlim([1540, 1580])
plt.savefig(pdfFilename)
plt.show()
