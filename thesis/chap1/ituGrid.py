"""
International Telecommunication association grid for DWDM and CWDM

DWDM Grids: 12.5, 25, 50, 100 GHz

Not very useful for the thesis.
Date: September 23rd 2020
"""

import matplotlib.pyplot as plt
import numpy as np

c = 299792458

centralWavelengthsCWDM = np.asarray([1271, 1291, 1311, 1331, 1351, 1371, 1391, 1411, 1431, 1451, 1471, 1491, 1511, 1531, 1511, 1531, 1551, 1571, 1591, 1611])*1e-9

def listGridCentralFrequenciesDWDM(channelSpacing):
    """ Return a list of the central frequencies of WDM channels. 
    Standard channel spacing for DWDM is 12.5, 25, 50, 100, 200, 300, 400, etc..."""
    frequencyList = []
    for n in range(-10, 10+1):
        frequencyList.append(193.1e12 + n * channelSpacing)
    return np.asarray(frequencyList)

def wavelength(frequency):
    return 2.99792458e8/frequency

def frequency(wavelength):
    return 2.99792458e8/wavelength

freqs = listGridCentralFrequenciesDWDM(25e9)
plt.stem(wavelength(freqs)*1e9, [0.9]*len(freqs), color='r')
plt.stem(centralWavelengthsCWDM*1e9, [0.9]*len(centralWavelengthsCWDM), color='r')
axes = plt.gca()
axes.set_ylim([0,1])
plt.show()

print((frequency(1261e-9)-frequency(1621e-9)))
print(frequency(1261e-9))
print(frequency(1621e-9))
