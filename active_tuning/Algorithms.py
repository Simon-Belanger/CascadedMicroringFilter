import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

# Coordinates descent Algorithm
def CoordsDescent(MRF, order_cycle, number_iter=2, plot_maps=False):
    """
    Coarse tuning of a MRF object using the coordinates descent algorithm.

    Inputs:
        MRF         : Microring Filter instance on which to perform the tuning
        order_cycle : array with the order in which to perform the tuning for the n rings (e.g. [0 1 2] or [2 1 0] or [1 2 0])

    example: CoordsDescent(MRF, [0, 1, 2, 3, 4], 2)
    """

    # TODO : Set the min and max bias value according to realistic values. Set the resolution according to the DC source.

    # Possible bias values
    bias_min, bias_max, bias_points  = 0, 5, 100
    bias_testpoints = np.linspace(bias_min, bias_max, bias_points).tolist()

    # Turn on the laser and set the wavelength
    phase, power = [],[]
    for i in range(1, number_iter + 1):  # For each iteration
        for j in range(len(MRF.Rings)):  # For each ring
            drop_list, thru_list = [],[]
            for k in bias_testpoints:  # For each bias value
                MRF.apply_bias(order_cycle[j], k)
                drop_list.append(MRF.measure_power(MRF.target_wavelength)[0])
                thru_list.append(MRF.measure_power(MRF.target_wavelength)[1])
                label = 'Iteration #' + str(i) + ' Ring #' + str(order_cycle[j])
            if plot_maps==True:
                plotsweep(bias_testpoints, drop_list, thru_list, label)
            MRF.apply_bias(order_cycle[j], bias_testpoints[drop_list.index(max(drop_list))])
            phase.append(MRF.get_total_phase(MRF.target_wavelength))
            power.append(max(drop_list))
    return phase, power

# Gradient descent Algorithm
def GradientDescent(MRF):
    pass

# Nelder Mead simplex algorithm
def NelderMead(MRF, options=None):
    """Fine tuning of a MRF object using the Nelder Mead simplex algorithm.


        Minimization of scalar function of one or more variables using the Nelder-Mead algorithm.
    """
    # Nelder-Mead algorithm parameters
    if options==None:
        options = {'maxiter':       None,       # Maximum allowed number of iterations
                    'maxfev':       None,       # Maximum allowed number offunction evaluations
                                                # Will default to N*200, where N is the number of variables, if neither maxiter or maxfev is set.
                                                # If both maxiter and maxfev are set, minimization will stop at the first reached.
                    'xatol':        0.001,      # Absolute error in xopt between iterations that is acceptable for convergence.
                    'fatol':        0.001,      # Absolute error in func(xopt) between iterations that is acceptable for convergence.
                    'adaptive':     False,      # Adapt algorithm parameters to dimensionality of problem. Useful for high-dimensional minimization [1].
                    'disp':         True,       # Set to True to print convergence messages.
                    'initial_simplex': None}    # Initial simplex. If given, overrides x0. initial_simplex[j,:] should contain the coordinates
                                                # of the j-th vertex of the N+1 vertices in the simplex, where N is the dimension.

    # Optimization
    MRF.NM_phase = []
    MRF.NM_power = []
    res = optimize.minimize(fun=MRF.minimize_MRF,
                            x0=MRF.bias,            # Initial guess. Array of real elements of size (n,), where ‘n’ is the number of independent variables.
                            args=(),                # Extra arguments passed to the objective function and its derivatives (fun, jac and hess functions).
                            method='Nelder-Mead',
                            tol=1e-6,               # Tolerance for termination. For detailed control, use solver-specific options.
                            callback=callbackNM,          # Function called after each iteration.
                            options=options)        # A dictionary of solver options.
    return res

Nfeval = 1
def callbackNM(Xi):
    global Nfeval
    #self.NM_phase.append(self.get_total_phase())
    #self.NM_power.append(self.measure_power(lambda_0)[0])
    print('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}   {4: 3.6f}   {5: 3.6f}'.format(Nfeval, Xi[0], Xi[1], Xi[2], Xi[3], Xi[4]))
    Nfeval += 1

def tuneMRF(MRF, wavelength):
    """Tune/stabilise the MRF object using coarse + fine algorithms"""

    MRF.target_wavelength = wavelength

    CoordsDescent(MRF, [0,1,2,3,4], 2)  # Coordinates descent
    NelderMead(MRF)  # Nelder Mead

def plotsweep(bias_testpoints, drop_list, thru_list, label):
    """"""
    plt.plot(bias_testpoints, drop_list, label='Drop port')
    plt.plot(bias_testpoints, thru_list, label='Thru port')
    plt.xlabel('Voltage [V]')
    plt.ylabel('Transmission [dB]')
    plt.title(label)
    plt.legend()
    plt.show()

def plot_convergence(phase, power):
    """ Plot the convergence of the algorithm.

        Inputs:
            phase   :   Convergence of the total phase at the drop port vs iteration#.
            power   :   Convergence of the measured power at the drop port vs iteration#.

        Outputs:
            (Convergence plot is displayed on screen)
    """
    fig, ax1 = plt.subplots()
    ax1.set_title('Convergence plot')

    ax1.set_xlabel('Iteration')
    ax1.set_ylabel(r'$\log_{10}\left|\left|\phi\right|\right|$', color='b')
    ax1.tick_params('y', colors='b')
    ax2 = ax1.twinx()
    ax2.set_ylabel('Drop port transmission [dB]', color='r')
    ax2.tick_params('y', colors='r')

    ax1.plot(phase, 'b')
    ax2.plot(power, 'r')
    fig.show()
