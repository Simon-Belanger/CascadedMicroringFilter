from math import sqrt
import numpy as np
import matplotlib.pyplot as plt

def flatTopCoupling(order, powerCoupling):
    """ Apodization profiles for the directional couplers. 
        B.E.Little Maximally Flat response for a N rings filter where N is the order 

    :param k_pow: Power cross-coupling coefficient.
    :type k_pow: float
    :return: Fiels cross-coupling coefficient list for all N+1 couplers.
    :rtype: list
    """
    if order == 1:
        k = sqrt(powerCoupling)
        return [k, k]
    if order == 2:
        k1, k2 = sqrt(powerCoupling), sqrt(0.250) * powerCoupling
        return [k1, k2, k1]
    elif order == 3:
        k1, k2 = sqrt(powerCoupling), sqrt(0.125) * powerCoupling
        return [k1, k2, k2, k1]
    elif order == 4:
        k1, k2, k3 = sqrt(powerCoupling), sqrt(0.100) * powerCoupling, sqrt(0.040) * powerCoupling
        return [k1, k2, k3, k2, k1]
    elif order == 5:
        k1, k2, k3 = sqrt(powerCoupling), sqrt(0.0955) * powerCoupling, sqrt(0.0295) * powerCoupling
        return [k1, k2, k3, k3, k2, k1]
    elif order == 6:
        k1, k2, k3, k4 = sqrt(powerCoupling), sqrt(0.0915) * powerCoupling, sqrt(0.0245) * powerCoupling, sqrt(0.0179) * powerCoupling
        return [k1, k2, k3, k4, k3, k2, k1]
    else:
        print('Order not supported.')

def displayApodizationCurves_couplingFixed(powerCoupling, ax=None):
    """[summary]

    Args:
        order ([type]): [description]
    """
    axNone=False
    if ax==None:
        axNone = True
        fig, ax = plt.subplots()

    # Main 
    for order in [2,3,4,5,6]:
        posVec = list(range(0, order+1))
        print(posVec)
        couplingCoeffs = flatTopCoupling(order, powerCoupling)
        ax.plot(np.asarray(posVec), 10*np.log(np.asarray(couplingCoeffs)**2/np.asarray(couplingCoeffs[0])**2))
        ax.set_xlabel('Index of the coupler relative to central coupler')
        ax.set_ylabel(r'Normalized power coupling coefficient $\frac{\kappa_i^2}{\kappa_1^2}$ [dB]')

    if axNone:
        plt.show()

def displayApodizationCurves_orderFixed(order, couplingList, ax=None):
    """[summary]

    Args:
        order ([type]): [description]
    """
    axNone=False

    if ax==None:
        axNone=True
        fig, ax = plt.subplots()

    # Main 
    for powerCoupling in couplingList:
        posVec = list(range(1, order+2))
        couplingCoeffs = flatTopCoupling(order, powerCoupling)
        ax.plot(np.asarray(posVec), 10*np.log(np.asarray(couplingCoeffs)**2/np.asarray(couplingCoeffs[0])**2), label='{:.2f}'.format(powerCoupling))
        ax.set_xlabel('Index of the coupler relative to central coupler')
        ax.set_ylabel(r'Normalized power coupling coefficient $\frac{\kappa_i^2}{\kappa_1^2}$ [dB]')
        ax.set_xlim([min(posVec), max(posVec)])
        ax.legend()

    if axNone:
        plt.show()