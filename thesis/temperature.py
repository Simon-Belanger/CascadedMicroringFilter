"""
Script used to produce graphs for the temperature dependance of Silicon around 1530 nm (Extrapolate??)

Simon Belanger-de Villers
28 Septembre 2020
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

overpassFont = fm.FontProperties(family = 'Overpass', 
                                fname = '/Users/simonbelanger/Library/Fonts/overpass-semibold.otf', 
                                size = 12)

def betaFromTemperature(temperature):
    " Data from J. A. McCaulley et al., PHYSICAL REVIEW B, 1994"

    # Coefficients
    coeffs = {'1150': [-0.6176e-5, 1.58632e-7, 2.56342e-9, -1.77411e-11, 4.72221e-14 -5.70088e-17, 25.8917e-21],
                    '1310': [-3.6137e-5, 8.65085e-7, -3.83712e-9, 1.00556e-11, -1.49840e-14, 1.18078e-17, -3.80552e-21],
                    '1530': [-3.7239e-5, 8.61435e-7 ,-3.72908e-9, 0.92278e-11, -1.27065e-14, 0.91077e-17, -2.64153e-21],
                    '2390': [-4.9739e-5, 8.91462e-7, -3.47969e-9, 0.71330e-11, -0.72036e-14, 0.30047e-17, -0.21432e-21]}

    label = '1150'
    T = temperature
    beta = coeffs[label][0] + coeffs[label][1]*T + coeffs[label][2]*T**2 + coeffs[label][3]*T**3 + coeffs[label][4]*T**4 + coeffs[label][5]*T**5 + coeffs[label][6]*T**6


pdfFilename 	= '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter3/temperatureSilicon.pdf'

# silicon absolute refractive index vs temperature @ 1500 nm
temperature     = [30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 295]
absoluteIndex   = [3.45309, 3.45319, 3.45340, 3.45373, 3.45417, 3.45471, 3.45535, 3.45609, 3.46100, 3.46760, 3.47540, 3.48314]

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(temperature, absoluteIndex, linewidth=2, marker='s')
ax.set_xlabel('Temperature [K]', font=overpassFont)
ax.set_ylabel('Material refractive index', font=overpassFont)
ax.set_xlim([min(temperature), max(temperature)])
ax.set_ylim([min(absoluteIndex), max(absoluteIndex)])
ax.grid(color='grey', linestyle='--', linewidth=0.25)
plt.savefig(pdfFilename)
plt.show()

# silicon dn/dT vs temperature @ 1500 nm
if False:
    temperature     = [30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 250, 295]
    indexDiff = [5.86E-06, 1.61E-05, 2.64E-05, 3.67E-05, 4.70E-05, 5.72E-05, 6.75E-05, 7.78E-05, 1.19E-04, 1.42E-04, 1.66E-04, 1.87E-04]

    plt.plot(temperature, indexDiff)
    plt.plot(temperature[1:], np.diff(absoluteIndex)/np.diff(temperature))
    plt.show()