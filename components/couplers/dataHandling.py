"""
Methods used to convert Legacy coupler simulation data in .mat format to python. This format was used in 2017 but since then it has been replaced by files
that are 100% handled in lumerical FDTD Solutions.

Author: Simon Belanger-de Villers
Date:   October 20th 2020
"""

import numpy as np
import pickle, scipy, os, sys, h5py
import matplotlib.pyplot as plt
from scipy import interpolate

# Extract a single field
def convertMatlabComplexDoublToPython(h5Variable):
    'Convert a complex matlab h5 variable to a python complex variable. '
    complexList = []
    for points in np.squeeze(h5Variable[:]):
        complexList.append(points[0] + points[1]*1j)
    return np.asarray(complexList, dtype=complex)

# Extract all fields for a given file
def scatteringParameterVariableFromMatlab(filename):
    """Extract the data from a Sparameter matlab datafile and save to dict.
    
    Example : 
    data    = scatteringParameterVariableFromMatlab('<fileLocation/.../Sparam.mat')
    
    """
    data = {}
    with h5py.File(filename, 'r') as f:
        data['S_add']       = convertMatlabComplexDoublToPython(f['S_add']['S'])
        data['S_drop']      = convertMatlabComplexDoublToPython(f['S_drop']['S'])
        data['S_input']     = convertMatlabComplexDoublToPython(f['S_input']['S'])
        data['S_through']   = convertMatlabComplexDoublToPython(f['S_through']['S'])
        data['wavelength']  = np.squeeze(np.asarray(f['S_add']['lambda']))
    return data

# For seriesRings
def scatteringParameterVariableFromMatlabSeries(filename):
    'Extract the data from a Sparameter matlab datafile and save to dict.'
    data = {}
    with h5py.File(filename, 'r') as f:

        data['S_add']       = convertMatlabComplexDoublToPython(f['S_41']['S'])
        data['S_drop']      = convertMatlabComplexDoublToPython(f['S_31']['S'])
        data['S_input']     = convertMatlabComplexDoublToPython(f['S_11']['S'])
        data['S_through']   = convertMatlabComplexDoublToPython(f['S_21']['S'])
        data['wavelength']  = np.squeeze(np.asarray(f['S_11']['lambda']))
    return data

# Extract all files from a directory
def pickleFromMultipleFiles(directory):
    """ 
    Example:
    dat = pickleFromMultipleFiles(<fileLocation/.../>)
    pickle.dump(dat, open( "file.pickle", "wb" ) )
    """

    # Save gap dependant data
    gapList = [180, 200, 220, 240, 260, 280, 300, 320, 340]
    dataDict = {'gap [m]': (np.asarray(gapList)*1e-9).tolist(), 'scattering': []}
    for gap in gapList:
        filename = os.getcwd() + '/DC_seriesrings_R5/SIMUL/DC_SERIESRINGS_' + str(gap) + '/Sparam.mat'
        print(filename)
        data = scatteringParameterVariableFromMatlabSeries(filename)
        dataDict['scattering'].append(data)

    # Save other parameters 
    dataDict['type']                    = 'seriesrings DC'
    dataDict['ring radius [m]']         = 5e-6
    dataDict['ring width [m]']          = 450e-9
    #dataDict['bus width [m]']           = 340e-9
    dataDict['waveguide height [m]']    = 220e-9
    #dataDict['coupling arc [deg]']      = 110
    return dataDict