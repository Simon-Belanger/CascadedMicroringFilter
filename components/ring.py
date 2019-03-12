# Ring Class, generates a model for a single ring and it's parameters
class Ring(object):

    def __init__(self, radius, neff, alpha_wg):
        " Constructor for the ring object. "
        self.radius     = radius        # Radius of the ring [m]
        self.neff       = neff          # Effective index of the waveguide [-]
        self.alpha_wg   = alpha_wg      # Propagation losses [dB/cm]

        self.L_rt       = self.get_roundtrip_length() # Roundtrip length of the ring resonator [m]
        self.alpha      = self.get_loss_coefficient() # Loss coefficient of the ring resonator [m-1]

    def get_roundtrip_length(self):
        " Getter for the roundtrip length attribute. "
        return 2*pi * self.radius

    def get_loss_coefficient(self):
        " Getter for the loss coefficient attribute. "
        return 0.23 / 2 * 100 * self.alpha_wg

    def get_complex_phase(self, wavelength):
        " Getter for the complex phase attribute. "

        # TODO : Consider dispersion (at least first order)
        # TODO : make a waveguide object and associate it to the ring object.

        # MODE data
        # neff_data   = np.array([2.4992, 2.4725, 2.4449, 2.4163, 2.3868])  # calculated using MODE
        # lambda_data = np.array([1500, 1525, 1550, 1575, 1600]) * 1e-9  # wavelength points in the MODE simulation
        # neff = np.interp(lambda_0, lambda_data, neff_data)

        return (2*pi * self.neff / wavelength - 1j * self.alpha) * self.L_rt

    def get_roundtrip_phase(self, wavelength):
        """ Setter for the roundtrip phase attribute.
        An input wavelength that is on resonance will yield a 0 (mod 2pi) roundtrip phase.
        """
        return (self.get_complex_phase(wavelength).real)%(2*pi)

    def set_phase_deviation(self):
        " Setter for the phase deviation attribute. "
        pass

    def set_tuning_phase(self):
        " Setter for the tuning phase attribute. "
        pass

    def get_total_phase(self):
        " Getter for the total phase attribute. "
        pass

    def get_propagation_matrix(self, wavelength):
        " Getter for the propagation matrix attribute. "