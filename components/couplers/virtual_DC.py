from misc.utils import *
import math

class virtual_DC(object):
    """"""

    def __init__(self, k_amp, Tc):

        # Build the Transfer matrix
        self.T_mat(self.get_k(k_amp), self.get_t(k_amp, Tc), Tc)

    @staticmethod
    def get_k(k_amp):
        return k_amp * 1j

    @staticmethod
    def get_t(k_amp, Tc):
        return math.sqrt(Tc - k_amp ** 2)

    def T_mat(self, k, t, Tc):
        self.T = 1 / k * np.matrix([[-Tc, 0, 0, t],
                                    [0, -Tc, t, 0],
                                    [0, -t, 1, 0],
                                    [-t, 0, 0, 1]])

    def coupling_matrix_complete(self, k, t, r=0.0, a=0.0):
        """
        Transfer matrix (TM) for the directional coupler obtained from the S matrix definition.

        :param k: field cross-coupling coefficient.
        :param t: field thru-coupling coefficient.
        :param r: field back-reflection coefficient. (default = 0.0)
        :param a: field add-coupling coefficient. (default = 0.0)
        :return: Coupling Transfer Matrix.
        """
        return s2t(np.matrix([[r, t, k, a],
                              [t, r, a, k],
                              [k, a, r, t],
                              [a, k, t, r]]))


if __name__ == '__main__':

    # Instantiate object
    DC = virtual_DC(0.3, 1.)

    print(DC.T)

