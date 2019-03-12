"Possibly a useless script since Bogaerts script has been written, should check 21/02/2019"

from math import cos
import math
import numpy as np
import matplotlib.pyplot as plt

R = 2.5e-6 # Ring radius [m]
loss = 3 # Waveguide losses [dB/cm]
n_eff = 2.5 # Effective index of the waveguide []
n_g = 4.5 # Group index of the waveguide []
C = 0.9 # Cross-coupling coefficient [W/W]

"""
To get the ring response.
"""
def get_t(C):
    t = math.sqrt(1-C)
    return t

def get_roundtrip_length(R):
    L_rt = 2 * math.pi * R # Round-trip length [m]
    return L_rt

def get_alpha(loss):
    alpha = 23 * loss # Linear loss coefficient [m-1]
    return alpha

def get_roundtrip_loss(loss):
    alpha = get_alpha(loss) # loss coefficient [m-1]
    L_rt = get_roundtrip_length(R) # Round-trip length [m]
    a = math.sqrt(math.exp(- alpha * L_rt))
    return a

def get_phase(lambda_0):
    beta = (2 * math.pi * n_eff)/lambda_0 # Propagation constant[m-1]
    L_rt = get_roundtrip_length(R) # Round-trip length [m]
    phi = beta * L_rt # Optical phase [rad]
    return phi

def get_transmission(lambda_0):
    a = get_roundtrip_loss(loss)
    phi = get_phase(lambda_0)
    t = get_t(C)
    return (a**2 - 2 * a * t * cos(phi) + t**2)/(1 - 2 * a * t * cos(phi) + (t * a)**2)

def sweep(lambda_min, lambda_max, lambda_points):
        """"""
        lambda_0 = np.linspace(lambda_min, lambda_max, lambda_points)
        T_t = np.zeros(len(lambda_0))

        for ii in range(len(lambda_0)):
            T_t[ii] = get_transmission(lambda_0[ii])

        plot_transmission(lambda_0, T_t)

def plot_transmission(lambda_0, T_t):
    """ Plot the drop/thru port transmission spectrum. """
    plt.plot(lambda_0 * 1e9, 10*np.log10(abs(T_t)**2))
    plt.xlabel('Wavelength (nm)', Fontsize=14)
    plt.ylabel('Power transmission (dB)', Fontsize=14)
    plt.show()

"""
To get ring parameters.
"""
def get_lambda_res(lambda_min, lambda_max):
    L_rt = get_roundtrip_length(R) # Round-trip length [m]

    # wavelengths taht are smaller than lambda_max but larger than lambda_min
    m = np.array(range(1000)) # Integer for number each order of resonance
    lambda_res = []
    for i in m[1:]:
        lambda_i = (n_eff * L_rt) / i
        if (lambda_i>=lambda_min) and (lambda_i<=lambda_max):
            lambda_res.append(lambda_i)
    return lambda_res

def get_fsr(lambda_c):
    L_rt = get_roundtrip_length(R) # Round-trip length [m]
    fsr = lambda_c**2/(n_g * L_rt) # Free spectral range [m]
    return fsr

def get_ER(C):
    a = get_roundtrip_loss(loss)
    t = get_t(C)
    T_t = ((t + a)/(1 + t * a))**2 # Maximum transmission (power floor  ) []
    R_min = ((t - a)/(1 - t * a))**2 # Minimum transmission []
    return 10 * math.log(T_t/R_min)

def get_FWHM(C):
    a = get_roundtrip_loss(loss)
    t = get_t(C)
    L_rt = get_roundtrip_length(R)
    if (math.pi * n_g * L_rt * math.sqrt(t * a)) == 0.:
        FWHM = 0
    else:
        FWHM = ((1 - t * a) * (get_lambda_res(1500e-9, 1600e-9)[0])**2)/(math.pi * n_g * L_rt * math.sqrt(t * a))
    return FWHM

def plot_ER(C):
    ER = []
    for c in C:
        ER.append(get_ER(c))

    # Plot
    plt.plot(ER, C)
    plt.ylabel('Coupling ratio', Fontsize=14)
    plt.xlabel('Extinction Ratio (dB)', Fontsize=14)
    plt.show()

def plot_FWHM(C):
    FWHM = []
    for c in C:
        FWHM.append(get_FWHM(c))

    # Plot
    plt.plot(np.asarray(FWHM)*1e9, C)
    plt.ylabel('Coupling ratio', Fontsize=14)
    plt.xlabel('FWHM (nm)', Fontsize=14)
    plt.show()

#plot_ER(np.linspace(0,1,100)[1:-1])
#plot_FWHM(np.linspace(0,1,100)[0:-2])

#sweep(1550e-9,1560e-9,1000)