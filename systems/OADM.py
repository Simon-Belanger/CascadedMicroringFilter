class OADM(object):
	""" 
	Optical Add-Drop Multiplexer (OADM) semi-abstract class used as a template for 
	simulations. 
	"""

	Input_port = 0
	Through_port = 0
	Drop_ports = [1,2,3,4]
	Add_ports = [1,2,3,4]

	def __init__(self, num_channels):

		self.num_channels = num_channels # Number of independant channels in the OADM
		self.filters      = [ChanDropFilter(...) for i in range(num_channels)]  # Ring resonator list

		pass

class ROADM(OADM):
	""" 
	Reconfigurable Optical Add-Drop Multiplexer (ROADM) semi-abstract class used as a 
	template for simulations. 
	"""

	def __init__(self, num_channels):
		super().__init__(num_channels)

class ChanDropFilter(object):
	"""
	Channel dropping filter that can be either tunable or not and which implements 
	functions for using it as part of a system.

	Channel dropping filters in their more general expression are passband filters centered 
	at a specific central wavelength and with a given bandwidth.

	"""

	central_wavelength = ...
	bandwidth = ...


	def __init__(object):
		pass

	def measure_drop(self, Input_port):
		""" Measure the dropped signal. """
		drop = []
		for wavelength, power in zip(Input_port.wavelength, Input_port.power):
			if (wavelength > central_wavelength-bandwidth/2) and (wavelength < central_wavelength+bandwidth/2) :
				drop.append(power)
			else:
				drop.append(0)
		return drop

	def measure_through(self, Input_port):
		""" Measure the signal going through the bus. """
		pass


if __name__ == "__main__":

	from signals import constant_spectrum_signal
	import matplotlib.pyplot as plt

	# 1) Generate an input signal
	input_sig = constant_spectrum_signal(100, 1500e-9, 1600e-9, 1)
	plt.plot(input_sig.wavelength, input_sig.power)
	plt.show()

	# 2) Make a filter and filter the input signal
	fil = ChanDropFilter(1550e-9, 5e-9)
	drop_sig = fil.measure_drop(input_sig)
	#through_sig = fil.measure_through(input_sig)

	#3) 
	plt.plot(drop_sig.wavelength, input_sig.power)
	plt.show()

