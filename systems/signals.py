import numpy as np
import unittest

class spectrum_signal(object):
	""" 
	Basic spectral signal consisting of a wavelength field and a constant power field (data field).

	Inputs :
		wavelength_list : Wavelength list 	[N x 1]
		power_list		: Power list 		[N x 1]
	Attributes:
		spectrum_signal.wavelength 	: Wavelength list 	[N x 1]
		spectrum_signal.power 		: Power list 		[N x 1]

	
	* Note that the power field could be pretty much anything. To be general it should be an S parameter 
		(complex amplitude and phase)

	Simon Bélanger-de Villers
	June 2019
	"""
	def __init__(self, wavelength_list, power_list):
		self.wavelength = wavelength_list
		self.power 		= power_list
		if len(self.wavelength) != len(self.power):
			raise Exception('The number of wavelength points should match the number of power points.')

class constant_spectrum_signal(spectrum_signal):
	""" 
	Power transmission spectrum object.

	Basic for now
	constant signal

	"""
	def __init__(self, number_points, min_wavelength, max_wavelength, power):
		wavelength_list = np.linspace(1500e-9,1600e-9,100).tolist()
		power_list 		= [power]*number_points
		super().__init__(wavelength_list, power_list)


# Unit Testing
class MyTestCase(unittest.TestCase):

	def test_equal(self):
		pass

	def test_notequal(self):
		""" Check that you get an error if the length of the two lists is not the same. """
		with self.assertRaises(Exception) as context:
			spectrum_signal([1,2,3], [1,3])
		self.assertTrue('The number of wavelength points should match the number of power points.' in str(context.exception))

if __name__ == "__main__":

	unittest.main()
