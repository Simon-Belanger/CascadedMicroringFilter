import csv
from components import MRF
import numpy as np
from misc.Filter import Analyze_results

"""

To generate MRF objects from a parameter file in csv format. Analyze them and save the data.

"""



def list_from_csv(csv_filename):
    """
    Read the csv file and extract the different variations as well as their parameters. Store the information for 
    each device in a dictionnary. Put all the different dicts in a list.


    The output list should look like this:

        variation_list = [var1_dict, var2_dict, ...]

    While each dict in the list should look like this:

        varX_dict = {'DeviceID': 'Num1', 'param1': X, 'param2': Y, ...}

        where 'DeviceID' is the name of the device (variation) and param1-N are the different design parameters

    example:

        variation_list = list_from_csv('file.csv'):
    """

    variation_list = []
    with open(csv_filename, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            variation_list.append(dict(row))
    return variation_list


def condition_parameters(variation_list):
    """
    Do some data analysis depending on the type of the input data.

    Convert strings to float

    :return:
    """

    for var in variation_list:

        var['R [um]'] = float(var['R [um]'])
        var['wR [nm]'] = float(var['wR [nm]'])
        var['wB [nm]'] = float(var['wB [nm]'])
        var['CA [deg]'] = float(var['CA [deg]'])
        try:
            var['g1 [nm]'] = float(var['g1 [nm]'])
        except ValueError:
            var['g1 [nm]'] = []
        try:
            var['g2 [nm]'] = float(var['g2 [nm]'])
        except ValueError:
            var['g2 [nm]'] = []
        try:
            var['g3 [nm]'] = float(var['g3 [nm]'])
        except ValueError:
            var['g3 [nm]'] = []
        try:
            var['g4 [nm]'] = float(var['g4 [nm]'])
        except ValueError:
            var['g4 [nm]'] = []
        try:
            var['g5 [nm]'] = float(var['g5 [nm]'])
        except ValueError:
            var['g5 [nm]'] = []
        try:
            var['g6 [nm]'] = float(var['g6 [nm]'])
        except ValueError:
            var['g6 [nm]'] = []

    return variation_list


def MRF_obj_from_file(dict_variation):
    """

    generate an object from the parameters in the dict to use this object for simulations/analysis.



    :return:
    """

    dict = {'DeviceID': '3R_V1', 'R [um]': 2.5, 'wR [nm]': 450.0, 'wB [nm]': 340.0, 'CA [deg]': 100.0, 'g1 [nm]': 200.0,'g2 [nm]': 210.0, 'g3 [nm]': 300.0, 'g4 [nm]': 200.0, 'g5 [nm]': [], 'g6 [nm]': []}

    # Extract the different parameters from the dictionnary
    name = dict['DeviceID']
    R = dict['R [um]'] * 1e-6
    wR = dict['wR [nm]'] * 1e-9
    wB = dict['wB [nm]'] * 1e-9
    CA = dict['CA [deg]']
    gap = [dict['g1 [nm]'], dict['g2 [nm]'], dict['g3 [nm]'], dict['g4 [nm]'], dict['g5 [nm]'], dict['g6 [nm]']]
    gap2 = (np.asarray([x for x in gap if x!= []])*1e-9).tolist()
    num_rings = 5

    # TODO diectional coupler class between gap and kappa->MRF
    # TODO for now, use the S parameters data files that have been previously generated

    # Create the MRF object with the proper parameters and return the object
    mrf = MRF.MRF(name, num_rings, R, k_amp=MRF.flattop5(0.3), Tc=[1., 1., 1., 1., 1., 1.], alpha_wg=3., crosstalk_coeff=[1, 0.75, 0.1])
    return mrf


def sweep_variations(variation_list):
    """
    Sweep the dict and simulate the circuit response for all devices.

    :return:
    """

    # Sweep variations and create circuit models for each
    for var in variation_list:

        # Create the circuit model
        lambda_0 = np.linspace(1530, 1560, 1000) * 1e-9
        mrf = MRF_obj_from_file(var)
        E_drop, E_thru = mrf.sweep(lambda_0, E_in=1, E_add=0, plot_results=False)
        P_drop, P_thru = 10 * np.log10(abs(E_drop) ** 2), 10 * np.log10(abs(E_thru) ** 2)

        # Analyze the results
        dat = Analyze_results(lambda_0, P_drop, P_thru, Max_width=3, FSR_min=10)
        print(dat)

        # Save the results to files
        #save_files(mrf)

    # TODO


def save_files():
    """
    save all the data files after a simulation.

    :return:
    """

    # TODO

    pass

def montecarlo():
    """

    From the parameters dict,

    :return:
    """
    # TODO

    pass

if __name__ == '__main__':

    """"""

    # Import variations in a list
    var_list = list_from_csv('DomHould_params.csv')

    # format parameters to be floats
    var_list_c = condition_parameters(var_list)

    # Sweep the different variations in the list
    MRFlist = sweep_variations(var_list_c)





