import numpy as np
import matplotlib.pyplot as plt 
import scipy.constants as cte

class LunaMeasurement:
    """ Luna OVA5000 measurement treatment """

    def __init__(self, file, filter_bw=30):

        self.file = file
        self.filter_bw = filter_bw
        self.load()

        # self.insertion_loss = self.apply_smoothing_filter(self.insertion_loss)
        # self.min_loss = self.apply_smoothing_filter(self.min_loss)
        # self.max_loss = self.apply_smoothing_filter(self.max_loss)
        # self.group_delay = self.apply_smoothing_filter(self.group_delay)


    @property
    def wavelength(self):
        return self.meas[:,0]

    @property
    def frequency(self):
        return self.meas[:,1]
    
    @property
    def insertion_loss(self):
        return self.apply_smoothing_filter(self.meas[:,2])

    @property
    def group_delay(self):
        return self.apply_smoothing_filter(self.meas[:,3])



    def load(self):
        meas = np.loadtxt(self.file, skiprows=9)
        self.meas = meas
        assert meas.shape[-1] == 25, "Measurement file is incomplete"


        # self.group_delay = meas[:,3]
        self.chromatic_dispersion = meas[:,4]
        self.polarization_dependent_loss = meas[:,5]
        self.polarization_mode_dispersion = meas[:,6]
        self.linear_phase_deviation = meas[:,7]
        self.quadratic_phase_deviation = meas[:,8]
        self.jones_a_amplitude = meas[:,9]
        self.jones_b_amplitude = meas[:,10]
        self.jones_c_amplitude = meas[:,11]
        self.jones_d_amplitude = meas[:,12]
        self.jones_a_phase   = meas[:,13]
        self.jones_b_phase = meas[:,14] 
        self.jones_c_phase = meas[:,15] 
        self.jones_d_phase = meas[:,16] 
        self.time = meas[:,17]  
        self.time_domain_amplitude = meas[:,18] 
        self.time_domain_wavelength = meas[:,19]    
        self.min_loss = meas[:,20]  
        self.max_loss = meas[:,21]  
        self.second_order_pmd = meas[:,22]
        self.phase_ripple_linear = meas[:,23]
        self.phase_ripple_quadratic = meas[:,24]

        self.jones = np.zeros((self.wavelength.shape[0], 2, 2), dtype=complex)
        self.jones[:,0,0] = self.jones_a_amplitude * np.exp(1j*self.jones_a_phase)
        self.jones[:,0,1] = self.jones_b_amplitude * np.exp(1j*self.jones_b_phase)
        self.jones[:,1,0] = self.jones_c_amplitude * np.exp(1j*self.jones_c_phase)
        self.jones[:,1,1] = self.jones_d_amplitude * np.exp(1j*self.jones_d_phase)

        # il = 10*np.log10(1/2*(np.abs(self.jones[:,0,0])**2 + np.abs(self.jones[:,0,1])**2 + 
                # np.abs(self.jones[:,1,0])**2 + np.abs(self.jones[:,1,1])**2))


        # plt.plot(self.wavelength, il)

        return self


    def apply_smoothing_filter(self, quantity):
        """ applies a smoothing filter to a measurement quantity """

        FWHM = 1e3*(self.wavelength[int(self.wavelength.shape[0]/2)])**2 / self.filter_bw / cte.c

        # conversion from FWHM to sigma
        sigma = FWHM/(np.sqrt(8*np.log(2)))
        gaussian = np.exp(-((self.time - self.time[int(self.time.shape[0]/2)])/sigma)**2)
        
        filtered_time_domain = gaussian * np.fft.fftshift(np.fft.fft(quantity))
        filtered_quantity = np.fft.ifft(np.fft.fftshift(filtered_time_domain))

        return filtered_quantity






# meas = LunaMeasurement("tanh1100/1580_1594_drop.txt", 500)

# plt.figure()
# plt.plot(meas.wavelength, meas.group_delay)
# plt.xlim((1575,1610))
# # plt.figure()
# # plt.plot(luna.time_domain_amplitude)
# plt.show()
