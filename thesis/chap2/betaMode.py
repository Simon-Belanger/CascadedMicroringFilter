" Plot the propagation constant beta over a wavelength span. "

import numpy as np 
import matplotlib.pyplot as plt
import math, sys
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.utils import field2PowerdB, powerdB2Field, powerdB2PowerLinear, powerLinear2powerdB, field2PowerLinear
from components.ring import Ring

ringRadius = 2.5e-6
effectiveIndex = 2.4449
groupIndex = 4.18
alpha_wg = 0
ring = Ring(ringRadius, effectiveIndex, groupIndex, alpha_wg)

beta = []
wavelengthArray = np.linspace(1500e-9, 1600e-9, 100)
for wavelength in wavelengthArray:
    beta.append(ring.getPropagationConstant(wavelength))

plt.plot(wavelengthArray, beta)
plt.show()