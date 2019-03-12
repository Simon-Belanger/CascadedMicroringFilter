import math
from misc.utils import *

class DC_bidirectionnal(object):

    def __init__(self, C=0.5, IL=3):
        # IL = Insertion loss [dB]
        # C = power cross-coupling [lin]

        # Build the Transfer matrix
        self.T_mat(C, IL)


    @staticmethod
    def get_k(C, IL):
        return math.sqrt(10**(-IL/10)) * math.sqrt(C) * 1j

    @staticmethod
    def get_t(C,IL):
        return math.sqrt(10**(-IL/10)) * math.sqrt(1 - C)

    def S_mat(self, C, IL):
        k = self.get_k(C, IL)
        t = self.get_t(C, IL)
        return np.matrix([[t, k],
                          [k, t]])

    def T_mat(self, C, IL):
        Tc = 
        self.T = 1/k * np.matrix([[t, -Tc],
                                  [1, -t]])

if __name__ == '__main__':

    dc = DC_bidirectionnal(0.5, 3)