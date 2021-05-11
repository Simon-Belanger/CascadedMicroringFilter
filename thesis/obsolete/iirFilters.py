from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

# Filter parameters
filterOrder = 10
cutoffFrequencies = [100, 200]
filterType = 'bandpass'


b, a = signal.butter(filterOrder, cutoffFrequencies, filterType, analog=True)
w, h = signal.freqs(b, a)
plt.plot(w, 20 * np.log10(abs(h)))
plt.xscale('log')
plt.title('Butterworth filter frequency response')
plt.xlabel('Frequency [radians / second]')
plt.ylabel('Amplitude [dB]')
plt.margins(0, 0.1)
plt.grid(which='both', axis='both')
plt.axvline(100, color='green') # cutoff frequency
plt.axvline(200, color='green') # cutoff frequency
plt.show()