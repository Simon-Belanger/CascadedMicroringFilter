import math
from misc.utils import *
import cmath

# Waveguide attributes
IL = 3 # Insertion loss [dB]
C = 0.5  # power cross-coupling [lin]

def get_k(C, IL):
    return math.sqrt(10**(-IL/10)) * math.sqrt(C) * 1j

def get_t(C,IL):
    return math.sqrt(10**(-IL/10)) * math.sqrt(1 - C)

def wg_Smat(C, IL):
    k = get_k(C, IL)
    t = get_t(C, IL)
    return np.matrix([[t, k],
                      [k, t]])

# TODO add wavelength dependance for the coupling coefficient (power)

# TODO use complex kappa S-param


def TMM(E_1, E_2):
    return wg_Smat(C, IL) * np.matrix([[E_1], [E_2]])

def measure_pow(port):

    power_lin = abs(TMM(1.,0.)[port])**2 # Linear power [W]
    if power_lin <= 0:
        power_dB = -100.
    else:
        power_dB = 10 * math.log10(power_lin)
    return power_dB

def measure_phase(port):
    return cmath.phase(TMM(1.,0.)[port])

print("Thru pow = " + str(measure_pow(0)))
print("Cross pow = " + str(measure_pow(1)))

print("Thru phase = " + str(measure_phase(0)))
print("Cross phase = " + str(measure_phase(1)))