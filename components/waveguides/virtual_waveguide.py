from misc.utils import *
import math
import cmath as cm
import numpy as np

class virtual_WG(object):
    """ Virtual Waveguide that can be used as a template for other subtypes of waveguides. """

    def __init__(self):
        pass

    def get_transfer_matrix(self, phi):
        """ Transfer matrix of the simple waveguide [2x2]. """
        return np.matrix([[0, cm.exp(-1j*phi)],
                          [cm.exp(1j*phi), 0]])

    def get_scattering_matrix(self, phi):
        """ Scattering matrix of the simple waveguide [2x2]. """
        return np.matrix([[0, cm.exp(-1j*phi)],
                          [cm.exp(-1j*phi), 0]])


if __name__ == '__main__':

    wg = virtual_WG()

    print(wg.get_transfer_matrix(1))




