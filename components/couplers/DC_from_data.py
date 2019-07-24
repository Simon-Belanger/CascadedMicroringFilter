import numpy as np
from misc.utils import *
import math

class DC_from_data(object):
	"""
	Directional coupler object from coupling data.

	Very basic wavelength dependance.

	Inputs :
		wavelength_vec 	: Wavelength data [m] (list)
		kappa_vec 		: Power cross-coupling data [m] (list)	
		t_vec 		 	: Power thru-coupling data [m] (list)	

	"""

	def __init__(self, wavelength_vec, kappa_vec, t_vec):

		self.wavelength_vec = wavelength_vec
		self.kappa_vec 		= kappa_vec
		self.t_vec 			= t_vec

	def get_transfer_matrix(self, wavelength):

		k = self.get_closestvalue_from_dataset(self.wavelength_vec, self.kappa_vec, wavelength)
		t = self.get_closestvalue_from_dataset(self.wavelength_vec, self.kappa_vec, wavelength)
		r = 0.0
		a = 0.0

		return s2t(np.matrix([[r, t, a, k],
								[t, r, k, a],
								[a, k, r, t],
								[k, a, t, r]]))

	@staticmethod
	def get_closestvalue_from_dataset(wavelength_dataset, value_dataset, target_wavelength):
		"""
		Find the coefficient that is closest to the target wavelength.
		"""
		array = np.asarray(wavelength_dataset)
		idx = (np.abs(array - target_wavelength)).argmin()
		return value_dataset[idx]

if __name__ == "__main__":

	dc = DC_from_data([1,2,3], [4,5,6], [7,8,9])
	print(dc.get_closestvalue_from_dataset([1,2,3], [4,5,6], 4))



