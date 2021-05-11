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
dataDir         = os.getcwd() + '/thesis/chap2/data.nosync/'
pdfDirectory    = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter2/waveguide/'
teFilename      = dataDir + 'TEmode_W450nm_T220nm_R2500nm.mat'
tmFilename      = dataDir + 'TMmode_W450nm_T220nm_R2500nm.mat'
wgWidth, wgThick = 0.45, 0.22

def formatField(field, x, y):
    """ Conversion of a matrix of tuples to a matrix of complex numbers """
    fieldComplex = np.zeros((len(y), len(x)), dtype=complex)
    for i in range(len(y)):
        for j in range(len(x)):
            fieldComplex[i,j] = field[i, j][0] +1j*field[i, j][1]
    return fieldComplex

def formatArray(array):
    """ Conversion of a matrix of tuples to a matrix of complex numbers """
    arrayComplex = np.zeros(len(array), dtype=complex)
    for i in range(len(array)):
            arrayComplex[i] = array[i][0] +1j*array[i][1]
    return arrayComplex

def loadFieldProfiles(filename):
    """ Load the electric and magnetic fields from the matlab file"""
    f = h5py.File(filename, 'r')

    x = np.squeeze(np.array(f.get('x')))
    y = np.squeeze(np.array(f.get('y')))
    Ex = formatField(np.array(f.get('Ex')), x, y)
    Ey = formatField(np.array(f.get('Ey')), x, y)
    Ez = formatField(np.array(f.get('Ez')), x, y)
    Hx = formatField(np.array(f.get('Hx')), x, y)
    Hy = formatField(np.array(f.get('Hy')), x, y)
    Hz = formatField(np.array(f.get('Hz')), x, y)

    return x,y,Ex,Ey,Ez,Hx,Hy,Hz

def loadDispersionData(filename):
    """ Load the dispersion data from the matlab file"""
    g = h5py.File(filename, 'r')

    f   = np.squeeze(np.array(g.get('f')))
    f_D     = np.squeeze(np.array(g.get('f_D')))
    D       = np.squeeze(np.array(g.get('D')))
    f_vg    = np.squeeze(np.array(g.get('f_vg')))
    vg      = np.squeeze(np.array(g.get('vg')))
    loss    = np.squeeze(np.array(g.get('loss')))
    beta    = formatArray(np.squeeze(np.array(g.get('beta'))))
    neff    = formatArray(np.squeeze(np.array(g.get('neff'))))

    return f, f_D, D, f_vg, vg, loss, beta, neff

"""
1 - Test bidon
"""
if False:
    # Load the data from a mat file containing the fields
    dir = os.getcwd() + '/thesis/chap2/data.nosync/'
    x,y,Ex,Ey,Ez,Hx,Hy,Hz = loadFieldProfiles(dir + 'testModeProfiles.mat')

    wgWidth, wgThick = 0.45, 0.22




    Y, X = np.meshgrid(y*1e6, x*1e6, copy=False, indexing='ij')

    fig, ax = plt.subplots()

    field = abs(Ex)**2+abs(Ey)**2+abs(Ez)**2

    cs      = ax.contourf(X, Y, field/field.max(), levels = 10, cmap=plt.cm.hot, vmin=0, vmax=1)
    ax.contour(cs, colors='k',linewidths=0.5)
    wgRect  = patches.Rectangle((-wgWidth/2, -wgThick/2), wgWidth, wgThick, linewidth=1, linestyle='dashed', edgecolor='white', facecolor='none')
    ax.add_patch(wgRect)
    cbar    = fig.colorbar(cs, ax=ax)
    ax.set_aspect("equal")
    cbar.ax.set_ylabel('Intensity (normalized)', font=overpassFont)
    ax.set_xlabel(r'x position [$\mu$m]', font=overpassFont)
    ax.set_ylabel(r'y position [$\mu$m]',font=overpassFont)
    ax.set_xlim([-0.8,0.8])
    ax.set_ylim([-0.8,0.8])
    #plt.savefig(filename)
    plt.show()

"""
2 - E field Intensity TE
"""
if False:    
    x,y,Ex,Ey,Ez,Hx,Hy,Hz = loadFieldProfiles(teFilename)
    
    # field to plot 
    field = abs(Ex)**2+abs(Ey)**2+abs(Ez)**2

    
    fig, ax = plt.subplots()
    Y, X = np.meshgrid(y*1e6, x*1e6, copy=False, indexing='ij')
    cs      = ax.contourf(X, Y, field, levels = 10, cmap=plt.cm.Reds, vmin=0, vmax=1)
    ax.contour(cs, colors='k',linewidths=0.5)
    wgRect  = patches.Rectangle((-wgWidth/2, -wgThick/2), wgWidth, wgThick, linewidth=1, linestyle='dashed', edgecolor='white', facecolor='none')
    ax.add_patch(wgRect)
    cbar    = fig.colorbar(cs, ax=ax)
    ax.set_aspect("equal")
    cbar.ax.set_ylabel('Intensity ', font=overpassFont)
    ax.set_xlabel(r'x position [$\mu$m]', font=overpassFont)
    ax.set_ylabel(r'y position [$\mu$m]',font=overpassFont)
    ax.set_xlim([-0.6,0.6])
    ax.set_ylim([-0.6,0.6])
    plt.savefig(pdfDirectory + 'TE_EfieldIntensity.pdf')
    plt.show()

"""
2 - H field Intensity TE
"""
if False:    
    x,y,Ex,Ey,Ez,Hx,Hy,Hz = loadFieldProfiles(teFilename)
    
    # field to plot 
    field = abs(Hx)**2+abs(Hy)**2+abs(Hz)**2

    
    fig, ax = plt.subplots()
    Y, X = np.meshgrid(y*1e6, x*1e6, copy=False, indexing='ij')
    cs      = ax.contourf(X, Y, field, levels = 10, cmap=plt.cm.Blues)
    ax.contour(cs, colors='k',linewidths=0.5)
    wgRect  = patches.Rectangle((-wgWidth/2, -wgThick/2), wgWidth, wgThick, linewidth=1, linestyle='dashed', edgecolor='white', facecolor='none')
    ax.add_patch(wgRect)
    cbar    = fig.colorbar(cs, ax=ax)
    ax.set_aspect("equal")
    cbar.ax.set_ylabel('Intensity (normalized)', font=overpassFont)
    ax.set_xlabel(r'x position [$\mu$m]', font=overpassFont)
    ax.set_ylabel(r'y position [$\mu$m]',font=overpassFont)
    ax.set_xlim([-0.6,0.6])
    ax.set_ylim([-0.6,0.6])
    plt.savefig(pdfDirectory + 'TE_HfieldIntensity.pdf')
    plt.show()

    """
3 - E field Intensity TM
"""
if False:    
    x,y,Ex,Ey,Ez,Hx,Hy,Hz = loadFieldProfiles(tmFilename)
    
    # field to plot 
    field = abs(Ex)**2+abs(Ey)**2+abs(Ez)**2

    
    fig, ax = plt.subplots()
    Y, X = np.meshgrid(y*1e6, x*1e6, copy=False, indexing='ij')
    cs      = ax.contourf(X, Y, field/field.max(), levels = 10, cmap=plt.cm.Reds, vmin=0, vmax=1)
    ax.contour(cs, colors='k',linewidths=0.5)
    wgRect  = patches.Rectangle((-wgWidth/2, -wgThick/2), wgWidth, wgThick, linewidth=1, linestyle='dashed', edgecolor='white', facecolor='none')
    ax.add_patch(wgRect)
    cbar    = fig.colorbar(cs, ax=ax)
    ax.set_aspect("equal")
    cbar.ax.set_ylabel('Intensity (normalized)', font=overpassFont)
    ax.set_xlabel(r'x position [$\mu$m]', font=overpassFont)
    ax.set_ylabel(r'y position [$\mu$m]',font=overpassFont)
    ax.set_xlim([-0.6,0.6])
    ax.set_ylim([-0.6,0.6])
    plt.savefig(pdfDirectory + 'TM_EfieldIntensity.pdf')
    plt.show()

"""
4 - H field Intensity TM
"""
if False:    
    x,y,Ex,Ey,Ez,Hx,Hy,Hz = loadFieldProfiles(tmFilename)
    
    # field to plot 
    field = abs(Hx)**2+abs(Hy)**2+abs(Hz)**2

    
    fig, ax = plt.subplots()
    Y, X = np.meshgrid(y*1e6, x*1e6, copy=False, indexing='ij')
    cs      = ax.contourf(X, Y, field, levels = 10, cmap=plt.cm.Blues)
    ax.contour(cs, colors='k',linewidths=0.5)
    wgRect  = patches.Rectangle((-wgWidth/2, -wgThick/2), wgWidth, wgThick, linewidth=1, linestyle='dashed', edgecolor='white', facecolor='none')
    ax.add_patch(wgRect)
    cbar    = fig.colorbar(cs, ax=ax)
    ax.set_aspect("equal")
    cbar.ax.set_ylabel('Intensity (normalized)', font=overpassFont)
    ax.set_xlabel(r'x position [$\mu$m]', font=overpassFont)
    ax.set_ylabel(r'y position [$\mu$m]',font=overpassFont)
    ax.set_xlim([-0.6,0.6])
    ax.set_ylim([-0.6,0.6])
    plt.savefig(pdfDirectory + 'TM_HfieldIntensity.pdf')
    plt.show()

"""
5 - Dispersion TE
"""
if  True:    

    teDispFilename      = dataDir + 'TEmodeDispersion_W450nm_T220nm_R2500nm.mat'
    f = h5py.File(teDispFilename, 'r')

    freq = np.squeeze(np.array(f.get('f')))
    wvl = constants.c/freq
    neff = formatArray(np.squeeze(np.array(f.get('neff'))))

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(wvl*1e9, np.real(neff))
    ax.set_xlabel(r'Wavelength, $\lambda$ [nm]', font=overpassFont)
    ax.set_ylabel(r'Effective index, $n_{eff}$',font=overpassFont)
    ax.set_xlim([wvl.min()*1e9,wvl.max()*1e9])
    #ax.set_ylim([-0.6,0.6])
    plt.savefig(pdfDirectory + 'TE_Dispersion.pdf')
    plt.show()

"""
5 - Dispersion TM
"""
if True:    

    teDispFilename      = dataDir + 'TMmodeDispersion_W450nm_T220nm_R2500nm.mat'
    f = h5py.File(teDispFilename, 'r')

    freq = np.squeeze(np.array(f.get('f')))
    wvl = constants.c/freq
    neff = formatArray(np.squeeze(np.array(f.get('neff'))))

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(wvl*1e9, np.real(neff))
    ax.set_xlabel(r'Wavelength, $\lambda$ [nm]', font=overpassFont)
    ax.set_ylabel(r'Effective index, $n_{eff}$',font=overpassFont)
    ax.set_xlim([wvl.min()*1e9,wvl.max()*1e9])
    #ax.set_ylim([-0.6,0.6])
    plt.savefig(pdfDirectory + 'TM_Dispersion.pdf')
    plt.show()

"""
Autre
"""
if False:

    filename      = dataDir + 'TEmodeDispersion_W450nm_T220nm_R2500nm.mat'
    f, f_D, D, f_vg, vg, loss, beta, neff = loadDispersionData(filename)

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(constants.c/f_vg, constants.c/vg)
    plt.show()

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(constants.c/f_D, D)
    plt.show()

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(constants.c/f, loss)
    plt.show()

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(constants.c/f, beta)
    plt.show()
    #ax.set_xlabel(r'Wavelength, $\lambda$ [nm]', font=overpassFont)
    #ax.set_ylabel(r'Effective index, $n_{eff}$',font=overpassFont)
    #ax.set_xlim([wvl.min()*1e9,wvl.max()*1e9])
    ##ax.set_ylim([-0.6,0.6])
    #plt.savefig(pdfDirectory + 'TM_Dispersion.pdf')