import sys, os
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
import scipy.io as sio
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches, rcParams
import matplotlib.font_manager as fm
from scipy import constants


overpassFont = fm.FontProperties(family = 'Overpass', 
                                fname = '/Users/simonbelanger/Library/Fonts/overpass-semibold.otf', 
                                size = 12)

# Global parameters
dataDir         = os.getcwd() + '/thesis/chap3/data/v1/'
pdfDirectory    = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter3/'

"""
1 - temperature profile across the rings. Comparison between solutions from the 2D and from the 3D heat diffusion 
    equation 

    dans les deux cas P = 0.01 W = 10 mW

    
"""
if True:    
    # 2D data
    f = h5py.File(dataDir + 'heatDiffusion2D.mat', 'r')
    yPos2D = np.squeeze(np.array(f.get('y')))
    temp2D = np.squeeze(np.array(f.get('T')))
    # 3D data
    f = h5py.File(dataDir + 'heatDiffusion3D.mat', 'r')
    yPos3D = np.squeeze(np.array(f.get('y')))
    temp3D = np.squeeze(np.array(f.get('T')))

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(yPos2D, temp2D)
    ax.plot(yPos3D, temp3D)
    ax.set_xlabel(r'y Position [$\mu$m]', font=overpassFont)
    ax.set_ylabel(r'Temperature, $T$ [K]',font=overpassFont)
    #ax.set_xlim([wvl.min()*1e9,wvl.max()*1e9])
    #ax.set_ylim([-0.6,0.6])
    plt.savefig(pdfDirectory + 'heatProfile2D3D.pdf')
    plt.show()