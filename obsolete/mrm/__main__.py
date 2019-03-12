# Code to test the response of a microring resonator when being phase modulated
from math import cos,sin,atan,pi,exp,sqrt
import numpy as np
import cmath
import matplotlib.pyplot as plt

class MRM(object):

    def __init__(self, a, r):
        self.a = a
        self.r = r

    def transmission(self, phi):
        """ Compute the field transmission for a given detuning. """
        return (self.r - self.a * cmath.exp(1j * phi)) / (1 - self.a * self.r * cmath.exp(1j * phi))

    def phase(self, phi):
        """ Compute the effective phase shift for a given detuning. """
        return pi + phi + atan((self.r * sin(phi))/(self.a - self.r * cos(phi))) + atan((self.r * self.a * sin(phi))/(1 - self.r * self.a * cos(phi)))

    def sweep_detuning(self, phi):
        """ Sweep the detuning from -pi to pi. """
        mod, phs = [],[]
        for i in phi:
            mod.append(abs(self.transmission(i))**2)
            phs.append(self.phase(i))
        plt.plot(phi/pi, mod)
        plt.plot(phi/pi, phs)
        plt.xlabel("Detuning $\phi$")
        plt.show()

    def sweep_complex(self, phi):
        """ Plot the symbol constellation on the complex plane."""
        trace = []
        # Trace
        for i in phi:
            mod = abs(self.transmission(i)) ** 2 / abs(self.transmission(phi[0])) ** 2
            phs = self.phase(i) - self.phase(phi[0])
            trace.append(cmath.rect(mod, phs))

        states = self.get_states([min(phi),min(phi)/1.5,max(phi)/1.5,max(phi)],phi[0])
        print(states)

        plt.plot(np.asarray(trace).real, np.asarray(trace).imag)
        plt.scatter(np.asarray(states).real, np.asarray(states).imag)
        plt.xlim([-1,1]); plt.ylim([-1, 1])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('Re(E)')
        plt.ylabel('Im(E)')
        plt.grid()
        plt.show()
        print(str((self.phase(phi[-1])-self.phase(phi[0]))/pi))
        return trace, states

    def get_states(self, phi_states, phi_0):
        pts = []
        for i in phi_states:
            mod = abs(self.transmission(i)) ** 2 / abs(self.transmission(phi_0)) ** 2
            phs = self.phase(i) - self.phase(phi_0)
            pts.append(cmath.rect(mod, phs))
        return pts

def IQ(I=[-1, -1/3, 1/3, 1]):
    """ Using the In-phase signal of one BPSK, this will produce the constellation of the IQ modulator. """

    ref_sig = np.asarray([-1, -1/3, 1/3, 1])
    ref_const = np.add.outer(ref_sig, ref_sig*cmath.exp(1j*pi/2))

    # In phase
    I = np.asarray(I)
    plt.scatter(I.real,I.imag)
    plt.xlim([-1,1]); plt.ylim([-1, 1])
    plt.title("In-phase")
    plt.show()

    # In quadrature
    Q = I * cmath.exp(1j*pi/2)
    plt.scatter(Q.real,Q.imag)
    plt.xlim([-1,1]); plt.ylim([-1, 1])
    plt.title("Quadrature")
    plt.show()

    # Total
    M = np.add.outer(I,Q)
    plt.scatter(ref_const.real, ref_const.imag)
    plt.scatter(M.real, M.imag)
    plt.title("In-phase + Quadrature")
    plt.grid()
    plt.show()

if __name__ == "__main__":
    MR = MRM(a=0.855, r=0.85) #a=0.855, r=0.85
    halfspan = pi/73 #pi/73
    phi = np.linspace(-halfspan, halfspan, 1000)
    MR.sweep_detuning(phi)
    trace, states = MR.sweep_complex(phi)


    IQ(states)





