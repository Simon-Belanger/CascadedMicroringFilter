import numpy as np
from scipy.sparse import spdiags


def thermal_crosstalk(num_rings, kappa_1, kappa_2, kappa_3):
    """Creates the crosstalk matrix according to parameters
        num_rings : number of actuators
        kappa_1   : direct heating
        kappa_2   : indirect heating, first order
        kappa_3   : indirect heating, second order
    """
    data = np.array([[kappa_1] * num_rings, [kappa_2] * num_rings, [kappa_2] * num_rings, [kappa_3] * num_rings,
                     [kappa_3] * num_rings])
    diags = np.array([0, -1, 1, -2, 2])
    return np.asmatrix(spdiags(data, diags, num_rings, num_rings).toarray())


def s2t(s):
    """Convert a 4x4 S-matrix to a 4x4 T-matrix."""

    if len(s) != 4 and len(s[0]) != 4:
        print('Error: Input S-matrix must have dimensions 4x4.')
        return

    s11, s21, s31, s41 = s[0, 0], s[1, 0], s[2, 0], s[3, 0]
    s12, s22, s32, s42 = s[0, 1], s[1, 1], s[2, 1], s[3, 1]
    s13, s23, s33, s43 = s[0, 2], s[1, 2], s[2, 2], s[3, 2]
    s14, s24, s34, s44 = s[0, 3], s[1, 3], s[2, 3], s[3, 3]

    t11 = (s11 * s32 * s43 - s11 * s33 * s42 - s12 * s31 * s43 + s12 * s41 * s33 + s13 * s31 * s42 - s13 * s41 *
           s32)/(s31 * s42 - s41 * s32)
    t12 = (s11 * s32 * s44 - s11 * s34 * s42 - s12 * s31 * s44 + s12 * s41 * s34 + s14 * s31 * s42 - s14 * s41 *
           s32)/(s31 * s42 - s41 * s32)
    t13 = (s11 * s42 - s12 * s41) / (s31 * s42 - s41 * s32)
    t14 = (s12 * s31 - s11 * s32) / (s31 * s42 - s41 * s32)
    t21 = (s21 * s32 * s43 - s21 * s33 * s42 - s22 * s31 * s43 + s22 * s41 * s33 + s23 * s31 * s42 - s23 * s41 *
           s32)/(s31 * s42 - s41 * s32)
    t22 = (s21 * s32 * s44 - s21 * s34 * s42 - s22 * s31 * s44 + s22 * s41 * s34 + s24 * s31 * s42 - s24 * s41 *
           s32)/(s31 * s42 - s41 * s32)
    t23 = (s21 * s42 - s22 * s41) / (s31 * s42 - s41 * s32)
    t24 = (s22 * s31 - s21 * s32) / (s31 * s42 - s41 * s32)
    t31 = (s32 * s43 - s33 * s42) / (s31 * s42 - s41 * s32)
    t32 = (s32 * s44 - s34 * s42) / (s31 * s42 - s41 * s32)
    t33 = s42 / (s31 * s42 - s41 * s32)
    t34 = -s32 / (s31 * s42 - s41 * s32)
    t41 = (s41 * s33 - s31 * s43) / (s31 * s42 - s41 * s32)
    t42 = (s41 * s34 - s31 * s44) / (s31 * s42 - s41 * s32)
    t43 = -s41 / (s31 * s42 - s41 * s32)
    t44 = s31 / (s31 * s42 - s41 * s32)

    return np.matrix([[t11, t12, t13, t14], [t21, t22, t23, t24], [t31, t32, t33, t34], [t41, t42, t43, t44]])


def t2s(t):
    """Convert a 4x4 T-matrix to a 4x4 S-matrix."""

    if len(t) != 4 and len(t[0]) != 4:
        print('Error: Input T-matrix must have dimensions 4x4.')
        return

    t11, t21, t31, t41 = t[0, 0], t[1, 0], t[2, 0], t[3, 0]
    t12, t22, t32, t42 = t[0, 1], t[1, 1], t[2, 1], t[3, 1]
    t13, t23, t33, t43 = t[0, 2], t[1, 2], t[2, 2], t[3, 2]
    t14, t24, t34, t44 = t[0, 3], t[1, 3], t[2, 3], t[3, 3]

    s11 = (t13 * t44 - t14 * t43) / (t33 * t44 - t34 * t43)
    s21 = (t23 * t44 - t24 * t43) / (t33 * t44 - t34 * t43)
    s31 = t44 / (t33 * t44 - t34 * t43)
    s41 = -t43 / (t33 * t44 - t34 * t43)
    s12 = (t14 * t33 - t13 * t34) / (t33 * t44 - t34 * t43)
    s22 = (t24 * t33 - t23 * t34) / (t33 * t44 - t34 * t43)
    s32 = -t34 / (t33 * t44 - t34 * t43)
    s42 = t33 / (t33 * t44 - t34 * t43)
    s13 = (t11 * t33 * t44 - t11 * t34 * t43 - t13 * t44 * t31 + t13 * t34 * t41 + t14 * t43 * t31 - t14 * t33 *
           t41)/(t33 * t44 - t34 * t43)
    s23 = (t21 * t33 * t44 - t21 * t34 * t43 - t23 * t44 * t31 + t23 * t34 * t41 + t24 * t43 * t31 - t24 * t33 *
           t41)/(t33 * t44 - t34 * t43)
    s33 = (t34 * t41 - t44 * t31) / (t33 * t44 - t34 * t43)
    s43 = (t43 * t31 - t33 * t41) / (t33 * t44 - t34 * t43)
    s14 = (t12 * t33 * t44 - t12 * t34 * t43 - t13 * t44 * t32 + t13 * t34 * t42 + t14 * t43 * t32 - t14 * t33 *
           t42)/(t33 * t44 - t34 * t43)
    s24 = (t22 * t33 * t44 - t22 * t34 * t43 - t23 * t44 * t32 + t23 * t34 * t42 + t24 * t43 * t32 - t24 * t33 *
           t42)/(t33 * t44 - t34 * t43)
    s34 = (t34 * t42 - t44 * t32) / (t33 * t44 - t34 * t43)
    s44 = (t43 * t32 - t33 * t42) / (t33 * t44 - t34 * t43)

    return np.matrix([[s11, s12, s13, s14], [s21, s22, s23, s24], [s31, s32, s33, s34], [s41, s42, s43, s44]])

def listmat_multiply(listmat):
    """Returns the product of subsequent muplitplications of matrices in a list."""
    H = listmat[0]
    for el in listmat[1::]:
        H = H * el
    return H