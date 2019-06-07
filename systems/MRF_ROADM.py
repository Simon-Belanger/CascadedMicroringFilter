from OADM import ROADM
from signals import constant_spectrum_signal

# Components required by the system in order to make it work

class PSR(object):
	"""
	Polarization Splitter Rotator (PSR) 

	Component used to split the TE0 and TM0 into the TE0 in two distinct branches.

	Simon Bélanger-de Villers
	Juin 2019
	"""

	def __init__(self):
		pass

class EC(object):
	"""
	Edge coupler (EC)

	Component used to split the TE0 and TM0 into the TE0 in two distinct branches.

	Simon Bélanger-de Villers
	Juin 2019
	"""

	def __init__(self):
		pass


class PowTaps(object):
	"""
	Power taps for measuring the optical power at the drop port of the channel.

	Simon Bélanger-de Villers
	Juin 2019
	"""

	def __init__(self):
		pass

# The actual system

class MRF_ROADM(ROADM):
	"""
	Reconfigurable Optical Add-Drop Multiplexer (ROADM) made with cascaded microring filters (MRFs)

	Simon Bélanger-de Villers
	Juin 2019
	"""

	def __init__(self):
		super().__init__()
		pass

if __name__ == "__main__":
	pass