import math
import cmath
from misc.utils import *
import matplotlib.pyplot as plt

# Important constants [global]
c = 299792458 # Velocity of light in vacuum [m/s]

# Waveguide attrivutes
L       = 1e-6      # Waveguide length [m]
n_eff   = 1         # Waveguide effective index
n_g     = 1         # Waveguide group index
D       = 2         # Waveguide dispersion coefficient [s m-1 m-1]
loss    = 0         # Waveguide propagation losses [dB/cm]
lambda_c= 1550e-9   # Central wavelength for taylor expansion [m]

def get_alpha(loss):
    loss_per_meter = loss * 100
    alpha = (0.23/2) * loss_per_meter # Linear loss coefficient (Field) [m-1]
    return alpha

def ang_freq(lambda_i):
    return (2 * math.pi * c) / lambda_i

def propagation_constant(lambda_0):
    omega = ang_freq(lambda_0)
    omega_c = ang_freq(lambda_c)
    d_omega = omega - omega_c
    beta_0 = (omega_c * n_eff) / c
    beta_1 = n_g / c
    beta_2 = -(D * lambda_c**2)/(2 * math.pi * c)
    beta = beta_0 + beta_1 * d_omega + 1/2 * beta_2 * d_omega**2
    return beta

def complex_phase(lambda_0):
    """"""
    beta = propagation_constant(lambda_0)
    alpha = get_alpha(loss)
    return (beta - alpha * 1j) * L

# Waveguide Scattering matrix
def wg_Smat(lambda_0):
    """"""
    phi = complex_phase(lambda_0)
    return np.matrix([[0, cmath.exp(-1j * phi)],
                      [cmath.exp(-1j * phi), 0]])

def TMM(lambda_0, E_1, E_2):
    return wg_Smat(lambda_0) * np.matrix([[E_1], [E_2]])

# Transmission
def measure_pow(lambda_0, port):

    power_lin = abs(TMM(lambda_0, 1.,0.)[port])**2 # Linear power [W]
    if power_lin <= 0:
        power_dB = -100.
    else:
        power_dB = 10 * math.log10(power_lin)
    return power_dB

def measure_phase(lambda_0, port):
    return cmath.phase(TMM(lambda_0, 1.,0.)[port])

def sweep(lambda_min, lambda_max, lambda_points):
    """"""
    lambda_0 = np.linspace(lambda_min, lambda_max, lambda_points)
    T_t = np.zeros(len(lambda_0))
    for ii in range(len(lambda_0)):
        T_t[ii] = measure_pow(lambda_0[ii], 1)
    plot_transmission(lambda_0, T_t)

def sweep_phase(lambda_min, lambda_max, lambda_points):
    """"""
    lambda_0 = np.linspace(lambda_min, lambda_max, lambda_points)
    p_t = np.zeros(len(lambda_0))
    for ii in range(len(lambda_0)):
        p_t[ii] = measure_phase(lambda_0[ii], 1)
    plot_phase(lambda_0, p_t)

def plot_transmission(lambda_0, T_t):
    plt.plot(lambda_0 * 1e9, T_t)
    plt.xlabel('Wavelength (nm)', Fontsize=14)
    plt.ylabel('Power transmission (dB)', Fontsize=14)
    plt.ylim([-100, 0])
    plt.show()

def plot_phase(lambda_0, p_t):
    plt.plot(lambda_0 * 1e9, np.unwrap(p_t))
    plt.xlabel('Wavelength (nm)', Fontsize=14)
    plt.ylabel('Phase (rad)', Fontsize=14)
    plt.show()

#sweep(1500e-9,1600e-9,1000)

sweep_phase(1500e-9,1600e-9,1000)