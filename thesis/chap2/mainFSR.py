"""
Script used to load the FSR relationship, from matlab
"""
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

pdfDirectory    = '/Users/simonbelanger/Documents/UL/Silicon_Photonics/07_Thesis/ulthese.nosync/fig/chapter2/fsr/'


overpassFont = fm.FontProperties(family = 'Overpass', 
                                fname = '/Users/simonbelanger/Library/Fonts/overpass-semibold.otf', 
                                size = 12)


radius = np.asarray([2.5e-6, 3.5e-6, 5.0e-6, 6.0e-6, 7.0e-6, 8.0e-6, 9.0e-06, 10.0e-06])
groupIndex = np.asarray([4.19544842336028, 4.18958856124690, 4.18605521051581, 4.18500058797412, 4.18436105059425, 4.18394452709435, 4.18365831195662, 4.18345326263329])
groupIndex2020 = np.asarray([3.351167164065220,4.192870347310389,4.190540076493132,4.189846797319420,4.189426851461823,4.189153531085577,4.188965802002462,4.188831351329833])

wavelength = 1500e-9

fsr = lambda groupIndex, radius, wavelength : wavelength**2 / (2 * np.pi * radius * groupIndex)

### FSR vs radius for different waveguide widths 
# group Index does not depend on ring radius
# !!! For narrow waveguides, does ng depend on radius
# V1
if False:
    # DATA
    widths      = np.linspace(200e-9, 800e-9, 6)
    groupIndex  = np.asarray([2.0911 ,4.2999, 4.2953, 4.0923, 3.9676, 3.8934])
    radius      = np.linspace(2e-6, 10e-6, 10)

    # Plot
    for w, ng in zip(widths[1:-1], groupIndex[1:-1]):
        plt.plot(radius*1e6, fsr(radius, ng, wavelength)*1e9, marker='s', mfc='none', mew=2, ms=8, linewidth=2, label='{0:.2f}'.format(w*1e9))
    plt.xlabel('Ring radius, R [$\mu$m]')
    plt.ylabel('Free spectral range, FSR [nm]')
    plt.grid()
    plt.legend()
    plt.savefig(pdfDirectory + 'fsrVSradiusVSwidth.pdf')
    plt.show()

# V2
if True:
    # DATA
    widths      = np.linspace(200e-9, 800e-9, 6)
    groupIndex  = np.asarray([2.0911 ,4.2999, 4.2953, 4.0923, 3.9676, 3.8934])
    radius      = np.linspace(2e-6, 10e-6, 10)

    # Plot
    colorCycle = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red']
    fig, Ax = plt.subplots(figsize=(10, 4))
    Ax.axhline(35, color ='k', linestyle='dashed',lw = 1) 
    Ax.axvline(2.5, color ='k', linestyle='dashed',lw = 1) 
    for w, ng, col in zip(widths[1:-1], groupIndex[1:-1], colorCycle):
        Ax.plot(radius*1e6, fsr(radius, ng, wavelength)*1e9, color=col, marker='s', mfc='none', mew=2, ms=8, linewidth=2, label='{0:.2f}'.format(w*1e9))
    Ax.legend(title='Width [nm]', prop=overpassFont, frameon=True, loc='upper right')
    Ax.grid(color='grey', linestyle='--', linewidth=0.25)
    Ax.set_xlabel('Ring radius, R [$\mu$m]', font=overpassFont)
    Ax.set_ylabel('Free spectral range, FSR [nm]', font=overpassFont)
    Ax.set_xlim([2, 10])
    plt.savefig(pdfDirectory + 'fsrVSradiusVSwidth.pdf')
    plt.show()



#### OBSOLETE

# Plot group index vs ring radius
if False:
    plt.plot(radius*1e6, groupIndex, marker='s', mfc='none', mec='red', color='red', mew=2, ms=8, linewidth=2)
    plt.xlabel('Ring radius, R [$\mu$m]')
    plt.ylabel('Free spectral range, FSR [nm]')
    plt.grid()
    plt.savefig('groupIndexVSradius.pdf')
    plt.show()

# Plot FSR vs ring radius
if False:
    plt.plot(radius*1e6, fsr(radius, groupIndex, wavelength)*1e9, marker='s', mfc='none', mec='red', color='red', mew=2, ms=8, linewidth=2)
    plt.plot(radius*1e6, fsr(radius, groupIndex2020, wavelength)*1e9, marker='s', mfc='none', mec='b', color='red', mew=2, ms=8, linewidth=2)
    plt.plot(radius*1e6, fsr(radius, groupIndex[-1], wavelength)*1e9, color='black', mew=2, ms=8, linewidth=2)
    plt.xlabel('Ring radius, R [$\mu$m]')
    plt.ylabel('Free spectral range, FSR [nm]')
    plt.grid()
    plt.savefig('fsrVSradius.pdf')
    plt.show()

### FSR vs wavelength for 500 nm wg
if False:
    radius = np.linspace(2e-6, 10e-6, 10)
    ng=4.18
    plt.plot(radius*1e6, fsr(radius, ng, 1530e-9)*1e9, marker='s', mfc='none', mew=2, ms=8, linewidth=2, label='{}'.format(1530))
    plt.plot(radius*1e6, fsr(radius, ng, 1550e-9)*1e9, marker='s', mfc='none', mew=2, ms=8, linewidth=2, label='{}'.format(1550))
    plt.plot(radius*1e6, fsr(radius, ng, 1570e-9)*1e9, marker='s', mfc='none', mew=2, ms=8, linewidth=2, label='{}'.format(1570))
    plt.xlabel('Ring radius, R [$\mu$m]')
    plt.ylabel('Free spectral range, FSR [nm]')
    plt.grid()
    plt.legend()
    plt.savefig('fsrVSradiusVSwavelength.pdf')
    plt.show()

### FSR strip vs rib
if False:
    radius = np.linspace(2e-6, 10e-6, 10)
    plt.plot(radius*1e6, fsr(radius, 4.18, 1550e-9)*1e9, marker='s', mfc='none', mew=2, ms=8, linewidth=2, label='strip')
    plt.plot(radius*1e6, fsr(radius, 3.85, 1550e-9)*1e9, marker='s', mfc='none', mew=2, ms=8, linewidth=2, label='rib')
    plt.xlabel('Ring radius, R [$\mu$m]')
    plt.ylabel('Free spectral range, FSR [nm]')
    plt.grid()
    plt.legend()
    plt.savefig('fsrVSradiusVSwg.pdf')
    plt.show()