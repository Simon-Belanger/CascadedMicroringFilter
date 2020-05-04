""" Code used to implement PN Junction simulation from the Python API directly into Lumerical. No Lumerical script files.
"""


from mrmfuncs import *
import time

# DEVICE materials 
silicon     = {'name': 'silicon', 'EM': 'Si (Silicon) - Palik', 'CT': 'Si (Silicon)', 'HT': 'Si (Silicon)'}
oxide       = {'name': 'oxide', 'EM': 'SiO2 (Glass) - Palik', 'CT': 'SiO2 (Glass) - SZE', 'HT': 'SiO2 (Glass) - SZE'}
aluminium    = {'name': 'aluminium', 'EM': 'Al (Aluminium) - Palik', 'CT': 'Al (Aluminium) - CRC', 'HT': 'Al (Aluminium) - CRC'}

## Parameters
waveguideOptions = {'ridgeWidth':       500e-9, # Ridge waveguide width                                                 [m]
                    'ridgeThickness':   220e-9, # Thickness of the ridge                                                [m]
                    'slabThickness':    90e-9,  # Thickness of the slab                                                 [m]
                    'boxThickness':     2e-6,   # Thickness of the BOX                                                  [m]
                    'claddingThickness':2e-6}   # Thickness of the cladding                                             [m]

pnJunctionOptions = {'pConc':   3.6e17,         # Light P doping Carrier concentration                                  [cm-3]
                     'ppConc':  1.8e17,         # Int. P doping Carrier concentration                                   [cm-3]
                     'pppConc': 4.4e20,         # Heavy P doping Carrier concentration                                  [cm-3]
                     'nConc':   2.1e17,         # Light N doping Carrier concentration                                  [cm-3]
                     'npConc':  2.9e18,         # Int. N doping Carrier concentration                                   [cm-3]
                     'nppConc': 4.4e20,         # Heavy N doping Carrier concentration                                  [cm-3]
                     
                     'pDmin':   0e-6,           # P implant start distance measured from the center of the ridge        [m]
                     'pDmax':   0e-6,           # P implant stop distance measured from the center of the ridge         [m]
                     'ppDmin':  0e-6,           # P+ implant start distance measured from the center of the ridge       [m]
                     'ppDmax':  2.5e-6,         # P+ implant stop distance measured from the center of the ridge        [m]
                     'pppDmin': 0.75e-6,        # P++ implant start distance measured from the center of the ridge      [m]
                     'pppDmax': 6.5e-6,         # P++ implant stop distance measured from the center of the ridge       [m]
                     'nDmin':   0e-6,           # N implant start distance measured from the center of the ridge        [m]
                     'nDmax':   0e-6,           # N implant stop distance measured from the center of the ridge         [m]
                     'npDmin':  0e-6,           # N+ implant start distance measured from the center of the ridge       [m]
                     'npDmax':  2.5e-6,         # N+ implant stop distance measured from the center of the ridge        [m]
                     'nppDmin': 0.75e-6,        # N++ implant start distance measured from the center of the ridge      [m]
                     'nppDmax': 6.5e-6}         # N++ implant stop distance measured from the center of the ridge       [m]

# Other parameters
simZ            = 8e-6      # Simulation depth for the junction [m]
contactWidth    = 500e-9    # Width of the aluminum contacts [m]


## Functions
def addMaterialToDevice(deviceApp, material):
    'Add a material to DEVICE using elements from a dict.'
    deviceApp.addmodelmaterial()
    deviceApp.set("name", material['name'])
    deviceApp.addmaterialproperties("EM", material['EM'])
    deviceApp.select("materials::" + material['name'])
    deviceApp.addmaterialproperties("CT", material['CT'])
    deviceApp.select("materials::" + material['name'])
    deviceApp.addmaterialproperties("HT", material['HT'])

## Build the junction using Python Code instead of Lumerical Code
with lumapi.DEVICE(hide=False) as device:

    # Add materials to DEVICE
    addMaterialToDevice(device, silicon)
    addMaterialToDevice(device, oxide)
    addMaterialToDevice(device, aluminium)

    # Add structures in DEVICE
    # TODO : Make a doped waveguide object out of this and make a function buildJunctionGeometryDevice
    device.addrect(name='cladding', x_min=-pnJunctionOptions['pppDmax'], x_max=pnJunctionOptions['nppDmax'], y_min=0, y_max=waveguideOptions['claddingThickness'], z=0, z_span=simZ, material=silicon['name'])
    device.addrect(name='slab', x_min=-pnJunctionOptions['pppDmax'], x_max=pnJunctionOptions['nppDmax'], y_min=0, y_max=waveguideOptions['slabThickness'], z=0, z_span=simZ, material=silicon['name'])
    device.addrect(name='ridge', x=0, x_span=waveguideOptions['ridgeWidth'], y_min=0, y_max=waveguideOptions['ridgeThickness'], z=0, z_span=simZ, material=silicon['name'])
    device.addrect(name='BOX', x_min=-pnJunctionOptions['pppDmax'], x_max=pnJunctionOptions['nppDmax'], y_min=-waveguideOptions['boxThickness'], y_max=0, z=0, z_span=simZ, material=oxide['name'])
    device.addrect(name='anode', x_min=-pnJunctionOptions['pppDmax'], x_max=-pnJunctionOptions['pppDmax']+contactWidth, y_min=waveguideOptions['slabThickness'], y_max=waveguideOptions['claddingThickness'], z=0, z_span=simZ, material=aluminium['name'])
    device.addrect(name='cathode', x_min=pnJunctionOptions['nppDmax']-contactWidth, x_max=pnJunctionOptions['nppDmax'], y_min=waveguideOptions['slabThickness'], y_max=waveguideOptions['claddingThickness'], z=0, z_span=simZ, material=aluminium['name'])

    # Add ground for R Slab measurements
    device.addrect(name='ground', x=0, x_span=waveguideOptions['ridgeWidth'], y_min=0, y_max=waveguideOptions['ridgeThickness'], z=0, z_span=simZ, material=aluminium['name'], enabled=0)

    # Add Simulation region + CHARGE solver + Boundary conditions
    device.select("simulation region")
    device.set('dimension', 3)              #  1 = 2D x-normal, 2 =  2D y-normal, 3 = 2D z-normal, 4 = 3D
    device.set('x min', -pnJunctionOptions['pppDmax'])
    device.set('x max', pnJunctionOptions['nppDmax'])
    device.set('y', waveguideOptions['ridgeThickness']/2)
    device.set('y span', 2e-6)
    device.set('z', 0)
    device.set('z_span', simZ)

    # Add FEEM solver + boundary conditions + mesh

    time.sleep(100)

