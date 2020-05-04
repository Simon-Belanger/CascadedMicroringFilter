"""Test script to perform circuit simulations of Microring Modulators.
"""
from MRM import MRM_Static
import os
from misc.dataIO import saveData, loadData
from components.couplers.mrmCoupler import mrmDC
from components.modulators.pnJunction import pnJunction

# Parameters
radius          = 10e-6     # [m]
loss_per_cm     = 2         # [dB/cm]
n_eff           = 2.5       # [-]
n_g             = 3.8       # [-]
gap             = 'CC'      # [m] Value or 'CC': Critically coupled, 'OC': Over Coupled, 'UC': Under Coupled
pnCoverage      = 70        # [%]
phaseshifter    = pnJunction(os.getcwd()+'/data/pnJunction/test1')
directionalCoupler = mrmDC(os.getcwd()+'/data/directionalCoupler/coupler1')

# Microring
mrm1 = MRM_Static(radius, loss_per_cm, n_eff, n_g, gap, pnCoverage, directionalCoupler, phaseshifter)

# Measuring results
directionalCoupler.plotCouplingCoefficient()
phaseshifter.plotPhaseShifterEfficiency()
phaseshifter.plotPropagationLosses()
mrm1.measureTransmissionMRM()