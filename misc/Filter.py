"""
Methods used to characterize the performances of an add-drop filter using the through/drop spectrum results.

Parameters & performance metrics extracted from the spectrum:

    General:
        - longitudinal modes : resonant wavelenghts [m] or frequencies [Hz]
        - longitudinal modes order (m) : integers corresponding to the order of the resonance (m = 1 corresponds to lambda = OPL)
        - FSR : Free spectral range [nm] (if showing more thant two peaks)

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

Main code -> ringAnalysis class


Author      : Simon Belanger-de Villers
Date        : 2018
Last edited : October 20th 2020
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths, argrelextrema
import statistics
import sys, math
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/')
from misc.utils import powerLinear2powerdB, field2PowerdB
from misc import constants

# Development
#TODO Mean of minimal values == Noise floor
#TODO Calculate the ripples in the power floor
#TODO Calculate the roll-off
#TODO Interpolate properly to get a precise bandwidth value

class ringAnalyisGeneral(object):
    """
    Class used to analyze the filter response of a general add/drop filter.
    """
    def __init__(self, wavelength, transmissionDrop=None, transmissionThrough=None, peakFindingOptions={'estimateHeight': True, 'estimateProminence':True}):
        
        # Input Data
        self.wavelength             = wavelength            # Wavelength array                                  [m]
        self.transmissionDrop       = transmissionDrop      # Transmission measured at the drop port array      [dB]
        self.transmissionThrough    = transmissionThrough   # Transmission measured at the through port array   [dB]

        # Options for peak finding
        self.estimateHeight         = peakFindingOptions['estimateHeight']
        self.estimateProminence     = peakFindingOptions['estimateProminence']

        # Find the transmission peaks 
        if transmissionDrop is not None:
            self.findTransmissionPeaksDrop()
        elif  transmissionThrough is not None:
            self.findTransmissionPeaksThrough()
        else:
            print(' The drop port transmission array required to proceed. ')

    def findTransmissionPeaksDrop(self):
        " use find_peaks to find the transmission peaks at the drop port of the filter [dB]. "

        # Find the peaks
        if self.estimateHeight:
            estimatedHeight = np.real(max(self.transmissionDrop))
            height      = (min(0.5*estimatedHeight,2.*estimatedHeight), max(0.5*estimatedHeight, 2.*estimatedHeight))          # (-10, 10) Required height of peaks (min, max)
            print('\nFinding peaks with an estimated height of = {} dB ...\n'.format(height))
        else:
            height=None
        threshold   = None
        distance    = None          # estimateFSRSampleSize Estimate the FSR ... a threshold based approach would be better
        if self.estimateProminence:
            estimatedProminence = np.abs(max(self.transmissionDrop) - min(self.transmissionDrop))
            prominence  = (0.7*estimatedProminence, 1.25*estimatedProminence)      # Estimated height of peaks vs noise floor [dB]
            print('\nFinding peaks with an estimated prominence of = {} dB ...\n'.format(prominence))
        else:
            prominence=None
        width       = None 
        wlen        = None
        rel_height  = 0.5
        plateau_size= None
        self.peaks, self.properties = find_peaks(self.transmissionDrop, height=height, threshold=threshold, distance=distance, prominence=prominence,
                                                                            width=width, wlen=wlen, rel_height=rel_height, plateau_size=plateau_size)

    def findTransmissionPeaksThrough(self):
        " use find_peaks to find the transmission peaks at the through port of the filter [dB]. "
        print(' Peak finding at the through port is not yet implemented. ')
        pass

    # Attributes
    @property
    def freeSpectralRange(self):
        'Free spectral range [m] of the resonator'
        if len(self.peaks)>1:
            FSR = [x - self.wavelength[self.peaks][i - 1] for i, x in enumerate(self.wavelength[self.peaks])][1:]
            return FSR
        else:
            print("Warning: Not enough information about the spectrum was given. At least 2 resonance peaks must be shown for the FSR to be calculated.")
            return None

    @property
    def insertionLossDrop(self):
        'Drop port insertion loss [dB]'
        return self.properties['peak_heights']

    @property
    def insertionLossThrough(self):
        'Through port insertion loss [dB]'
        return statistics.median(self.transmissionThrough)

    @property
    def extinctionRatio(self):
        ' Extinction ration of the through port [dB] '
        return self.insertionLossThrough - 10 * self.transmissionThrough[self.peaks]

    @property
    def noiseFLoorLevel(self):
        ' Average value of the noise floor [dB]'
        return statistics.median(self.transmissionDrop)

    @property
    def outOfBandRejectionRatio(self):
        ' Out of band rejection ratio [dB]'
        return self.properties['peak_heights'] - self.noiseFLoorLevel

    # Methods
    def measureRollOff(self):
        ' Plot vs frequency in db for Roll Off. '
        normedFrequency = constants.c/self.wavelength-self.measureLongitudinalModes('frequency')[0]
        powerLog        = self.transmissionDrop

        # Measure the Roll-off in dB/decade
        #plt.plot(np.diff(powerLog)/np.diff(normedFrequency))
        #plt.show()
        #plt.plot(normedFrequency, powerLog)
        #plt.xscale('log');plt.grid()
        #plt.show()
        return normedFrequency, powerLog
        

    def measureBandwidth(self, transmissionLevel=-3):
        ' Measure the bandwidth of the passband (drop port) at the given transmission level e.g. -3 dB'
        rel_height = 10**(transmissionLevel/10)
        results_half = peak_widths(self.transmissionDrop, self.peaks, rel_height=1-rel_height)
        return sampleSizeToMetres(results_half[0][:], self.wavelength)

    def measureLongitudinalModes(self, type='wavelength'):
        ' Return an array containing the resonant wavelengths [m] or frequencies [Hz] of the resonator. '
        if type == 'wavelength':
            return self.wavelength[self.peaks]
        elif type == 'frequency':
            return constants.c / self.wavelength[self.peaks]
        else: 
            print('Error : Bad type, exiting'); exit()

    def measureLongitudinalModesOrder(self, effectiveIndex, radius):
        ' Return the mode order (m) of the longitudinal modes using the theoretical effective index and ring radius.'
        print('Warning: This method is not very accurate')
        resonantWavelengths = self.measureLongitudinalModes('wavelength')
        return np.round(2*np.pi*effectiveIndex*radius/resonantWavelengths)

    def annotateSpectrum(self):
        ' Place the different metrics on a spectrum plot. s'
        pass

    def displayProps(self):
        # For the drop port

        print('wvl modes = ' +str(self.measureLongitudinalModes('wavelength')))
        print('freq modes = ' + str(self.measureLongitudinalModes('frequency')))
        print('FSR = ' + str(self.freeSpectralRange))
        print('IL_D =' + str(self.insertionLossDrop))
        print('Noise Floor = ' + str(self.noiseFLoorLevel))
        print('OBRR = ' + str(self.outOfBandRejectionRatio))
        print('BW3 = ' + str(self.measureBandwidth(-3)))
        print('BW20 = ' + str(self.measureBandwidth(-20)))
        print('BW40 = ' + str(self.measureBandwidth(-40)))

        # For the through port
        if self.transmissionThrough is not None:
            print('IL_T = ' + str(self.insertionLossThrough))
            print('ER_T = ' + str(self.extinctionRatio))

    def gui(self):
        # TODO : Implement a GUI
        pass

# Exceptions
class dataError(Exception): 
    pass

#    # Find the FSR from calculated peaks
#    fsr = FSR(wavelength, peaks)
#
#    # Sweep through each peak and returns data on it's shape
#    peak_data = []
#    for p in peaks:
#        # Initialize new dict to store information about the peak
#        peak_dict = {'IL_D': None, 'IL_T': None, 'ER': None, 'BW_3': None, 'BW_20': None, 'Ripples_D': None,
#                     'Ripples_T': None, 'OBRR': None}
#
#        # Noise floor
#        NoiseFloor = min(P_drop)
#
#        # Center the peak to analyze it more closely
#        lambda_0_c, P_drop_c, P_thru_c = center_peak(p, wavelength, P_drop, P_thru, Max_width)
#        # plt.plot(lambda_0_c * 1e9, P_drop_c, label='Drop', marker='x')
#        # plt.plot(lambda_0_c * 1e9, P_thru_c, label='Thru', marker='o')
#
#        # Calculate the Insertion Loss (IL) for both ports
#        peak_dict['IL_D'], peak_dict['IL_T'] = -max(P_drop_c), -max(P_thru_c)
#
#        # Calculate the bandwidth for the cropped dataset
#        peak_dict['BW_3'], peak_dict['BW_20'] = bandwidth(3, lambda_0_c, P_drop_c), bandwidth(20, lambda_0_c, P_drop_c)
#
#        # Calculate the ripples in the passband
#        max_val = argrelextrema(P_drop_c, np.greater_equal, order=1)
#        min_val = argrelextrema(P_drop_c, np.less_equal, order=1)
#        peak_dict['Ripples_D'] = max(P_drop_c[max_val[0]]) - min(P_drop_c[min_val[0][1:-1]])
#        # plt.hlines(max(P_drop_c[max_val[0]]), min(lambda_0_c)*1e9, max(lambda_0_c)*1e9, colors='k', linestyles='solid', label='', data=None)
#        # plt.hlines(min(P_drop_c[min_val[0][1:-1]]), min(lambda_0_c) * 1e9, max(lambda_0_c) * 1e9, colors='k',linestyles='solid', label='', data=None)
#
#        # Calculate the ripples in the thru port
#        max_val = argrelextrema(P_thru_c, np.greater_equal, order=1)
#        min_val = argrelextrema(P_thru_c, np.less_equal, order=1)
#        peak_dict['Ripples_T'] = max(P_thru_c[max_val[0][1:-1]]) - min(P_thru_c[min_val[0]])
#        # plt.hlines(max(P_thru_c[max_val[0][1:-1]]), min(lambda_0_c)*1e9, max(lambda_0_c)*1e9, colors='k', linestyles='solid', label='', data=None)
#        # plt.hlines(min(P_thru_c[min_val[0]]), min(lambda_0_c) * 1e9, max(lambda_0_c) * 1e9, colors='k',linestyles='solid', label='', data=None)
#
#        # Calculate the Extinction ratio (ER)
#        peak_dict['ER'] = max(P_thru_c) - min(P_thru_c)
#
#        # Calculate the Out of Band rejection ratio (OBRR)
#        peak_dict['OBRR'] = peak_dict['IL_D'] - NoiseFloor
#
#        # TODO Calculate the ripples in the power floor
#        # TODO Calculate the roll-off
#
#        peak_data.append(peak_dict)
#
#        plt.plot(wavelength * 1e9, P_drop, label='Drop')
#        plt.plot(wavelength * 1e9, P_thru, label='Thru')
#        plt.plot(wavelength[peaks] * 1e9, P_drop[peaks], label='peaks', marker='o', linestyle='', color='blue')
#        plt.xlabel('wavelength (nm)', Fontsize=14)
#        plt.ylabel('Power transmission (dB)', Fontsize=14)
#        plt.legend(loc='lower right')
#        plt.show()

#    return peak_data, fsr

def measureFSR(wavelength, peaks):
    " Returns the distance between adjacent peaks in wavelength units [m]. "
    if len(peaks)>1:
        FSR = [x - wavelength[peaks][i - 1] for i, x in enumerate(wavelength[peaks])][1:]
        return FSR
    else:
        print("Warning: Not enough information about the spectrum was given. At least 2 resonance peaks must be shown for the FSR to be calculated.")
        return None

def center_peak(peak, wavelength, P_drop, P_thru, Max_width):
    " Centers the spectrum on a single peak given by index peak. Returns cropped data. "
    C1 = wavelength >= (wavelength[peak] - Max_width * 1e-9 / 2)  # Condition 1 : Minimal wavelength (peak isolation)
    C2 = wavelength <= (wavelength[peak] + Max_width * 1e-9 / 2)  # Condition 2 : Maximal wavelength (peak isolation)
    new_range = np.logical_and(C1, C2) # Combined condition

    return wavelength[new_range], P_drop[new_range], P_thru[new_range]

def sampleSizeToMetres(sampleSize, wavelength):
    ' Convert the bandwidth in samples to the units of wavelength [m]'
    return sampleSize * (max(wavelength)-min(wavelength)) / len(wavelength)

def bandwidth(BW_excursion, wavelength, P_drop):
    """ This function returns the linewidth of the the passband filter. """
    thres = max(P_drop) - BW_excursion # Compute the theshold

    return max(wavelength[P_drop >= thres]) - min(wavelength[P_drop >= thres])

def measureBandwidth(wavelength, powerDropLinear, peaks, transmissionLevel=-3):
    ' Measure the bandwidth of the passband (drop port) at the given transmission level e.g. -3 dB'

    rel_height = 10**(transmissionLevel/10)
    results_half = peak_widths(powerDropLinear, peaks, rel_height=1-rel_height)
    return sampleSizeToMetres(results_half[0][:], wavelength)*1e9

# Misc functions 
def estimateFSRSampleSize(wavelength, estimatedFSR=35e-9):
    return len(wavelength)/(max(wavelength)-min(wavelength)) * estimatedFSR