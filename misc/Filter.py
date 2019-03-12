from components import MRF
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, argrelextrema

"""
Functions for MRF filter analysis. 
"""


def Analyze_results(lambda_0, P_drop, P_thru, Max_width=3, FSR_min=10):
    """

    From a dataset of wavelength, Drop port and thru port power, return the different peaks and their properties.
    As well as the fsr

    General:
        - FSR : Free spectral range

    Drop Port:
        - IL_D : Insertion loss [dB]
        - BW_3 : 3 dB bandwidth [nm]
        - BW_20 : 20 dB bandwidth [nm]
        - OBRR : Out of band rejection ratio [dB]
        - Ripples_D : Ripples amplitude [dB]
        # -Roll-off (dB/decade) to implement later

    Thru Port:
        - IL_T : Insertion loss [dB]
        - ER : Extinction ratio [dB]
        - Ripples_T : Ripples amplitude [dB]

    returns a dict with all the values and their key.

    :param lambda_0: Wavelength array [m x 1].
    :type lambda_0: ndarray
    :param P_drop: Drop port power values array [m x 1].
    :type P_drop: ndarray
    :param P_thru: Thru port power values array [m x 1].
    :type P_thru: ndarray
    :param Max_width: Estimated width of the peak [nm].
    :type Max_width: float
    :param FSR_min: Estimated FSR, peaks that are closer than this won't be considered [nm].
    :type FSR_min: float
    :return: List containing peak info, FSR
    :rtype: list, list
    """

    # Find the peaks
    peaks, _ = find_peaks(P_drop, height=None, threshold=None, distance=None, prominence=1, width=FSR_min, wlen=None, rel_height=0.5)

    # Find the FSR from calculated peaks
    fsr = FSR(lambda_0, peaks)

    # Sweep through each peak and returns data on it's shape
    peak_data = []
    for p in peaks:
        # Initialize new dict to store information about the peak
        peak_dict = {'IL_D': None, 'IL_T': None, 'ER': None, 'BW_3': None, 'BW_20': None, 'Ripples_D': None,
                     'Ripples_T': None, 'OBRR': None}

        # Noise floor
        NoiseFloor = min(P_drop)

        # TODO Interpolate properly to get a precise bandwidth value

        # Center the peak to analyze it more closely
        lambda_0_c, P_drop_c, P_thru_c = center_peak(p, lambda_0, P_drop, P_thru, Max_width)
        # plt.plot(lambda_0_c * 1e9, P_drop_c, label='Drop', marker='x')
        # plt.plot(lambda_0_c * 1e9, P_thru_c, label='Thru', marker='o')

        # Calculate the Insertion Loss (IL) for both ports
        peak_dict['IL_D'], peak_dict['IL_T'] = -max(P_drop_c), -max(P_thru_c)

        # Calculate the bandwidth for the cropped dataset
        peak_dict['BW_3'], peak_dict['BW_20'] = bandwidth(3, lambda_0_c, P_drop_c), bandwidth(20, lambda_0_c, P_drop_c)

        # Calculate the ripples in the passband
        max_val = argrelextrema(P_drop_c, np.greater_equal, order=1)
        min_val = argrelextrema(P_drop_c, np.less_equal, order=1)
        peak_dict['Ripples_D'] = max(P_drop_c[max_val[0]]) - min(P_drop_c[min_val[0][1:-1]])
        # plt.hlines(max(P_drop_c[max_val[0]]), min(lambda_0_c)*1e9, max(lambda_0_c)*1e9, colors='k', linestyles='solid', label='', data=None)
        # plt.hlines(min(P_drop_c[min_val[0][1:-1]]), min(lambda_0_c) * 1e9, max(lambda_0_c) * 1e9, colors='k',linestyles='solid', label='', data=None)

        # Calculate the ripples in the thru port
        max_val = argrelextrema(P_thru_c, np.greater_equal, order=1)
        min_val = argrelextrema(P_thru_c, np.less_equal, order=1)
        peak_dict['Ripples_T'] = max(P_thru_c[max_val[0][1:-1]]) - min(P_thru_c[min_val[0]])
        # plt.hlines(max(P_thru_c[max_val[0][1:-1]]), min(lambda_0_c)*1e9, max(lambda_0_c)*1e9, colors='k', linestyles='solid', label='', data=None)
        # plt.hlines(min(P_thru_c[min_val[0]]), min(lambda_0_c) * 1e9, max(lambda_0_c) * 1e9, colors='k',linestyles='solid', label='', data=None)

        # Calculate the Extinction ratio (ER)
        peak_dict['ER'] = max(P_thru_c) - min(P_thru_c)

        # Calculate the Out of Band rejection ratio (OBRR)
        peak_dict['OBRR'] = peak_dict['IL_D'] - NoiseFloor

        # TODO Calculate the ripples in the power floor
        # TODO Calculate the roll-off

        peak_data.append(peak_dict)

        plt.plot(lambda_0 * 1e9, P_drop, label='Drop')
        plt.plot(lambda_0 * 1e9, P_thru, label='Thru')
        plt.plot(lambda_0[peaks] * 1e9, P_drop[peaks], label='peaks', marker='o', linestyle='', color='blue')
        plt.xlabel('wavelength (nm)', Fontsize=14)
        plt.ylabel('Power transmission (dB)', Fontsize=14)
        plt.legend(loc='lower right')
        plt.show()


    return peak_data, fsr

def FSR(lambda_0, peaks):
    """
    Returns the distance between adjacent peaks in unit given by lambda_0.

    :param lambda_0: Wavelength array [m x 1].
    :type lambda_0: ndarray
    :param peaks: Peaks array [n x 1].
    :type peaks: ndarray
    :return: Array containing the FSR between adjacent peaks [(n-1) x 1].
    :rtype: ndarray
    """
    if len(peaks)>1:
        FSR = [x - lambda_0[peaks][i - 1] for i, x in enumerate(lambda_0[peaks])][1:]
        return FSR
    else:
        print("""Warning: Not enough information about the spectrum was given. 
                At least 2 resonance peaks must be shown for the FSR to be calculated.""")
        return None

def center_peak(peak, lambda_0, P_drop, P_thru, Max_width):
    """
    Centers the spectrum on a single peak given by index peak. Returns cropped data.

    :return:
    :param peak: Index corresponding to the peak of the interest in the spectrum.
    :type peak: int
    :param lambda_0: Wavelength array [m x 1].
    :type lambda_0: ndarray
    :param P_drop: Drop port power values array [m x 1].
    :type P_drop: ndarray
    :param P_thru: Thru port power values array [m x 1].
    :type P_thru: ndarray
    :param Max_width: Width of the peak (estimation) [nm].
    :type Max_width: float
    :return: lambda_c, P_drop_c, P_thru_c
    :rtype: ndarray, ndarray, ndarray
    """
    C1 = lambda_0 >= (lambda_0[peak] - Max_width * 1e-9 / 2)  # Condition 1 : Minimal wavelength (peak isolation)
    C2 = lambda_0 <= (lambda_0[peak] + Max_width * 1e-9 / 2)  # Condition 2 : Maximal wavelength (peak isolation)
    new_range = np.logical_and(C1, C2) # Combined condition

    return lambda_0[new_range], P_drop[new_range], P_thru[new_range]

def bandwidth(BW_excursion, lambda_0, P_drop):
    """
    This function returns the linewidth of the the passband filter.

    :param BW_excursion: Power attenatuation at which the bandwidth is measured (e.g. 3 dB).
    :type BW_excursion: float
    :param lambda_0: Wavelength array [m x 1].
    :type lambda_0:
    :param P_drop: Drop port power values array [m x 1].
    :type P_drop: ndarray
    :return: Width of the peak.
    :rtype: float
    """
    thres = max(P_drop) - BW_excursion # Compute the theshold

    return max(lambda_0[P_drop >= thres]) - min(lambda_0[P_drop >= thres])


if __name__ == '__main__':

    # Wavelength Array
    lambda_0 = np.linspace(1550, 1580, 1000) * 1e-9

    # Build the mrf and get the data
    mrf = MRF.MRF(name='1', num_rings=5, R=2.5e-6, k_amp=MRF.flattop5(0.3), Tc=[1., 1., 1., 1., 1., 1.], alpha_wg=3., crosstalk_coeff=[1, 0.75, 0.1])
    E_drop, E_thru = mrf.sweep(lambda_0, E_in=1, E_add=0, plot_results=False)
    P_drop, P_thru = 10*np.log10(abs(E_drop)**2), 10*np.log10(abs(E_thru)**2)

    # Analyze the filter response
    dat = Analyze_results(lambda_0, P_drop, P_thru, 5, 10)
    print(dat)
