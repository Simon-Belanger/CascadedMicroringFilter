import math

import numpy as np
from misc import lumerical


# THis file implements a ring object which has attributes and methods

# Future : make a waveguide object to be children of ring

class ring(object):

    def __init__(self, radius=2.5e-6, wg_width=450.e-9, wg_height=220.e-9, slab_thickness=0.):
        self.radius = radius
        self.wg_width = wg_width
        self.wg_height = wg_height
        self.slab_thickness = slab_thickness

    def fsr(self, wavelength):
        """"""
        L = 2 * math.pi * self.radius
        ng = self.ng(wavelength)
        return wavelength**2/(L * ng)*1e9

    def neff(self, wavelength):
        """ Get effective index using Lumerical MODE. """

        # Open a MODE instance
        mode1 = lumerical.mode()

        # Calculate ng for the given set of parameters
        mode1.inst.load("waveguide.lms")
        mode1.inst.select("WG")
        mode1.inst.set("Ridge width", self.wg_width)
        mode1.inst.set("Ridge thickness", self.wg_height)
        mode1.inst.set("Slab thickness", self.slab_thickness)
        mode1.inst.setnamed("FDE", "bend radius", 2.5e-6)
        mode1.inst.setnamed("FDE", "wavelength", wavelength)
        mode1.inst.findmodes()
        neff = np.squeeze(np.real(mode1.inst.getresult("FDE::data::mode1", "neff")))
        # Close session
        mode1.close_session()

        return neff

    def ng(self, wavelength):
        """ Get group index using Lumerical MODE. """
        # Open a MODE instance
        mode1 = lumerical.mode()

        # Calculate ng for the given set of parameters
        mode1.inst.load("waveguide.lms")
        mode1.inst.select("WG")
        mode1.inst.set("Ridge width", self.wg_width)
        mode1.inst.set("Ridge thickness", self.wg_height)
        mode1.inst.set("Slab thickness", self.slab_thickness)
        mode1.inst.setnamed("FDE", "bend radius", 2.5e-6)
        mode1.inst.setnamed("FDE", "wavelength", wavelength)
        mode1.inst.findmodes()
        ng = np.squeeze(np.real(mode1.inst.getresult("FDE::data::mode1", "ng")))
        # Close session
        mode1.close_session()

        return ng


r1 = ring()
print(r1.fsr(1550e-9))




