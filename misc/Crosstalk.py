import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags

def build_crosstalk_matrix(crosstalk_terms, matrix_size=None):
    """Creates a NxN symmetric crosstalk matrix, normalised.

        crosstalk_terms[0]   : self coupling
        crosstalk_terms[1]   : first neighbor coupling
        crosstalk_terms[2]   : second neighbor coupling
        ...
        crosstalk_terms[N]   : Nth neighbor coupling
    """

    if matrix_size==None:
        matrix_size = len(crosstalk_terms)

    diagonal_terms = crosstalk_terms[::-1][0:-1]+crosstalk_terms
    diagonal_indices = list(range(0,-len(crosstalk_terms),-1))[::-1][0:-1]+list(range(len(crosstalk_terms)))

    return diags(diagonal_terms, diagonal_indices, shape=(matrix_size, matrix_size)).toarray()

def exponential_crosstalk(num_rings, coupling_strength, plot_profile=False):
    """
    Returns the crosstalk coefficients that model an exponential crosstalk characteristic

        num_rings : Number of rings in the filter [int].
        coupling_strength :  The coupling strenght gives the ammount of crosstalk between the resonators [int].
    outputs:
        coeff: Crosstalk coefficients [Nx1].

    """
    x = np.linspace(0, num_rings - 1, num_rings)
    try:
        coeff = np.exp(-1/coupling_strength * x).tolist()
    except ZeroDivisionError:
        coeff = [1.] + ([0.]*(num_rings-1))

    if plot_profile==True:
        plt.stem(x,coeff)
        plt.xlabel('Neighbor')
        plt.ylabel('Crosstalk coefficient')
        plt.xticks(x)
        plt.show()

    return coeff

def linear_crosstalk():
    pass

if __name__ == '__main__':

    # Unit test
    assert exponential_crosstalk(5, 0) == [1.,0.,0.,0.,0.]

    # Build a crosstalk matrix of size 5x5 with exponential decay
    print(build_crosstalk_matrix(exponential_crosstalk(5, 0.1)))

    # Plot the profile
    exponential_crosstalk(5, 1, True)
