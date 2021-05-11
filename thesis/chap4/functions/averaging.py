import numpy as np
import scipy.constants as cte

# Methods
def apply_smoothing_filter(wavelength, time, quantity, filter_bw):
    """ applies a smoothing filter to a measurement quantity . From J. Cauchon. """

    FWHM = 1e3*(wavelength[int(wavelength.shape[0]/2)])**2 / filter_bw / cte.c

    # conversion from FWHM to sigma
    sigma = FWHM/(np.sqrt(8*np.log(2)))
    gaussian = np.exp(-((time - time[int(time.shape[0]/2)])/sigma)**2)
    
    filtered_time_domain = gaussian * np.fft.fftshift(np.fft.fft(quantity))
    filtered_quantity = np.fft.ifft(np.fft.fftshift(filtered_time_domain))

    return filtered_quantity