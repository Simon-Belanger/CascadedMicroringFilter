import numpy as np
""" Apodization profiles for the directional couplers."""

def flattop2(k_pow):
    """B.E.Little Maximally Flat response for a 2 rings filter

    :param k_pow: Power cross-coupling coefficient.
    :type k_pow: float
    :return: Fiels cross-coupling coefficient list for all 6 couplers.
    :rtype: list
    """
    k1, k2 = np.sqrt(k_pow), np.sqrt(0.250) * k_pow
    k3 = k1
    return [k1, k2, k3]

def flattop3(k_pow):
    """B.E.Little Maximally Flat response for a 3 rings filter

    :param k_pow: Power cross-coupling coefficient.
    :type k_pow: float
    :return: Fiels cross-coupling coefficient list for all 6 couplers.
    :rtype: list
    """
    k1, k4 = np.sqrt(k_pow), np.sqrt(k_pow)
    k2, k3 = np.sqrt(0.125) * k_pow, np.sqrt(0.125) * k_pow
    return [k1, k2, k3, k4]

def flattop4(k_pow):
    """B.E.Little Maximally Flat response for a 4 rings filter

    :param k_pow: Power cross-coupling coefficient.
    :type k_pow: float
    :return: Fiels cross-coupling coefficient list for all 6 couplers.
    :rtype: list
    """
    k1, k5 = np.sqrt(k_pow), np.sqrt(k_pow)
    k2, k4 = np.sqrt(0.100) * k_pow, np.sqrt(0.100) * k_pow
    k3 = np.sqrt(0.040) * k_pow
    return [k1, k2, k3, k4, k5]

def flattop5(k_pow):
    """B.E.Little Maximally Flat response for a 5 rings filter

    :param k_pow: Power cross-coupling coefficient.
    :type k_pow: float
    :return: Fiels cross-coupling coefficient list for all 6 couplers.
    :rtype: list
    """
    k1, k2, k3 = np.sqrt(k_pow), np.sqrt(0.0955) * k_pow, np.sqrt(0.0295) * k_pow
    k4, k5, k6 = k3, k2, k1
    return [k1, k2, k3, k4, k5, k6]

def flattop6(k_pow):
    """B.E.Little Maximally Flat response for a 6 rings filter

    :param k_pow: Power cross-coupling coefficient.
    :type k_pow: float
    :return: Fiels cross-coupling coefficient list for all 6 couplers.
    :rtype: list
    """
    k1, k7 = np.sqrt(k_pow), np.sqrt(k_pow)
    k2, k6 = np.sqrt(0.0915) * k_pow, np.sqrt(0.0915) * k_pow
    k3, k5 = np.sqrt(0.0245) * k_pow, np.sqrt(0.0245) * k_pow
    k4 = np.sqrt(0.0179) * k_pow
    return [k1, k2, k3, k4, k5, k6, k7]

def flattopN(k_pow):
    """B.E.Little Maximally Flat response for a N rings filter

    :param k_pow: Power cross-coupling coefficient.
    :type k_pow: float
    :return: Fiels cross-coupling coefficient list for all 6 couplers.
    :rtype: list
    """
    pass


if __name__ == "__main__":


    print([i ** 2 for i in flattop5(0.3)])