import sys
import os
sys.path.append('/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter')
from misc.lumerical import lumapi
import pickle

class DC_FDTD(object):
    """
    Directional coupler where the S parameters are extracted from FDTD simulations.
    
    Author      : Simon Bélanger-de Villers
    Created     : July 23 2019 
    Last edited : July 24 2019
    """

    meshAccuracy    = 1             # Mesh accuracy setting for the simulation
    startWavelength = 1530e-9       # Start wavelength for the simulation 
    stopWavelength  = 1565e-9       # Stop wavelength for the simulation

    def __init__(self, radius, waveguideWidth, ribHeight, gap):
        """
        Defines the parameters for the general structure
        e.g. simulation time, mesh accuracy
        """
        self.radius         = radius
        self.waveguideWidth = waveguideWidth
        self.ribHeight      = ribHeight
        self.gap            = gap

        pass

    def bindScript(self, lsfFilename):
        pass

    def measure(self):
        """
        Opens a Lumerical FDTD instance, set the simulation parameters and run the 
        simulation. Extract the data from the instance and close it.
        """

        # Open a new FDTD session
        self.handle = lumapi.open('fdtd')

        # Put variables in the workspace
        lumapi.putDouble(self.handle, 'R', self.radius)
        lumapi.putDouble(self.handle, 'W', self.waveguideWidth)
        lumapi.putDouble(self.handle, 'H', self.ribHeight)
        lumapi.putDouble(self.handle, 'g', self.gap)
        lumapi.putDouble(self.handle, 'mesh', self.meshAccuracy)
        lumapi.putDouble(self.handle, 'lambda_min', 1530e-9)
        lumapi.putDouble(self.handle, 'lambda_max', 1565e-9)

        # Load the simulation file 
        pathToSimFile = "/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/components/couplers/dependancies/simFiles/DC_seriesrings.fsp"
        lumapi.evalScript(self.handle, "load('" + pathToSimFile + "');")

        # Execute the script to run the simulation 
        pathToScriptFile = "/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/components/couplers/dependancies/scripts/DC_seriesrings.lsf"
        lumapi.evalScript(self.handle, "addpath('" +  os.path.dirname(os.path.abspath(pathToScriptFile)) + "/');")
        lumapi.evalScript(self.handle, os.path.basename(os.path.abspath(pathToScriptFile)).split('.')[0] + ";")

        # Get back variables
        self.data = {'S11': lumapi.getVar(self.handle, "S11"), 
                'S21': lumapi.getVar(self.handle, "S21"), 
                'S31': lumapi.getVar(self.handle, "S31"),
                'S41': lumapi.getVar(self.handle, "S41")}

        # Close the FDTD session
        lumapi.close(self.handle)

        self.saveDataToFile()

    def saveDataToFile(self):
        " Save the data stored in a dict to file using pickle. "

        pathToDataFile = "/Users/simonbelanger/Documents/UL/Silicon_Photonics/Python/CascadedMicroringFilter/components/couplers/dependancies/data/"
        dataFilename = "DCSR_R=" + str(self.radius*1e6).replace('.','_') + "um_" \
            + "wgW=" + str(self.waveguideWidth*1e9).replace('.','_') + "nm_" \
            + "wgH=" + str(self.ribHeight*1e9).replace('.','_') + "nm_" \
            + "g=" + str(self.gap*1e9).replace('.','_') + "nm_" \
            + "MSH=" + str(self.meshAccuracy) + ".pickle"
        outfile = open(pathToDataFile + dataFilename, 'wb')
        pickle.dump(self.data, outfile)
        outfile.close()


    def loadDataFromFile(self, datafile):
        """
        Load data from a datafile. 

        datafile = complete path to datafile with extension.
        """
        self.data = pickle.load( open( datafile, "rb" ))

    def plotData(self):
        """ Plot the data stored in the object. """
        pass

    def dataExists(self):
        """ 
        Check if data exists for the given properties. Important to do versionning 
        on the properties, each version must have a given set of properties.

        Example : if (!dataExists)
                    self.measure()
                    self.storeData()
        """
        pass

    def storeData(self):
        """
        Store the data in a proper file e.g. .JSON, .mat, .pckl etc
        Add the results to a a lookup table to check the parameters in dataExists().
        """
        pass


if __name__ == "__main__":
    
    dc = DC_FDTD(2.5e-6, 500e-9, 220e-9, 200e-9)
    dc.measure()