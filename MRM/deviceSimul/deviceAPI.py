"""
Lumerical DEVICE junction simulation using the Python API

Functions to perform EO characterization of the process simulated junction in 
order to extract the compact model parameters

Author 		: Simon Belanger-de Villers (simon.belanger-de-villers.1@ulaval.ca)
Created 	: August 9 2019
Last edited : August 29 2019
"""
import os, subprocess, time, datetime
import numpy as np
import matplotlib.pyplot as plt
import platform, imp

# Import Lumapi
if platform.system() == "Linux":
    lumapi = imp.load_source("lumapi", str("/opt/lumerical/fdtd/api/python/lumapi.py"))
elif platform.system() == "Windows":
    lumapi = imp.load_source("lumapi", "C:\\Program Files\\Lumerical\\FDTD\\api\\python\\lumapi.py")
elif platform.system() == "Darwin":
    lumapi = imp.load_source("lumapi", "/Applications/Lumerical 2020a.app/Contents/API/Python/lumapi.py")

# Default parameters for neff sweep
neffSweepDefaults = {'edge_per_wvl': 3, 'poly_order': 3, 'refinement': 20}
BWSweepDefaults = {}

# Misc parameters
hideGUI = True

def buildJunction(device):
    """ Build the junction in Lumerical DEVICE using the constructor script.

        This function creates a new .ldev file which sometimes gets corrupted.
    """
    t = time.time()
    deviceProjectFileLocation = os.getcwd() + "/deviceSimul/"
    device.putv("deviceProjectFileLocation", str(deviceProjectFileLocation))
    device.feval(deviceProjectFileLocation + "lumericalScriptFiles/buildJunction.lsf")
    print("{0} - The junction has been built successfully in {1:.2f} seconds.\n".format(datetime.datetime.now(), time.time() - t))

def importTDR(device, TDRDatafile):
    """"
    Import the data from a .tdr file into a device simulation file.

    deviceSimul     : Handle to an opened instance of device.
    TDRDatafile     : Path to the TDR datafile to be imported. (e.g. /gel/usr/<idul>/MRMProcessSimulations/recipe1_dir/recipe1_fps.tdr)
    """

    # Create a copy of .tdr as .h5
    h5FilePath = os.getcwd() + "/deviceSimul/" + os.path.basename(TDRDatafile).split('.')[0] +".h5"
    subprocess.check_call(['cp', TDRDatafile, h5FilePath])
    
    # Put the name of the file as a variable and run the script
    device.putv("filename", str(h5FilePath))
    device.feval(os.getcwd() + "/deviceSimul/lumericalScriptFiles/" + "importDopingProfile.lsf")
    
    # Remove the copy of the file 
    subprocess.check_call(['rm', h5FilePath])

    print("{0} - Data from {1} has been successfully imported in the simulation file.\n".format(datetime.datetime.now(), TDRDatafile + ".tdr"))

def neffSweep(device, wavelength, bias, plotResults=True):
    """ Returns the modulator's effective index vs reverse bias voltage. """

    # Pass sweep parameters to DEVICE instance 
    device.putv("V", bias)
    device.putv("lambda_0", wavelength)
    device.putv("edge_per_wvl", neffSweepDefaults['edge_per_wvl'])
    device.putv("poly_order", neffSweepDefaults['poly_order'])
    device.putv("refinement", neffSweepDefaults['refinement'])

    # Run the script in DEVICE
    device.feval(os.getcwd() + "/deviceSimul/lumericalScriptFiles/" + "getneff.lsf")

    # Retrieve the results from DEVICE
    neff = np.squeeze(device.getv("neff"))

    # Delete the log files, .mat file and the sweep dir
    subprocess.check_call(['rm', '-r', os.getcwd() + '/deviceSimul/junction_voltage'])
    subprocess.check_call(['rm', os.getcwd() + '/deviceSimul/junction_p0.log'])
    subprocess.check_call(['rm', os.getcwd() + '/deviceSimul/wg_charge.mat'])
    subprocess.check_call(['rm', os.getcwd() + '/deviceSimul/junction.ldev'])

    if plotResults==True:
        # Plot the results for d_neff
        plt.plot(bias, np.real(neff))
        plt.xlabel("Bias [V]")
        plt.ylabel("Effective index change")
        plt.title("Re(n_eff) vs bias for lambda={0:f} nm".format(wavelength))
        plt.show()

        # Plot the results for d_alpha
        # TODO
    

    return neff

def bandwidthSweep(device, wavelength, bias, plotResults=True):
    """
    Return the modulator RC bandwidth versus bias voltage.
    """

    # Pass parameters
    device.putv("V", format_bias(bias))
    # Mesh
    device.putv("min_edge", 0.025)
    device.putv("max_edge", 0.05)
    # Transient Solver
    device.putv("min_step", 0.05)
    device.putv("max_step", 4)
    device.putv("abs_lim", 0.001)
    device.putv("rel_lim", 0.001)
    # Convergence criterion
    device.putv("conv_criteria", "update")
    device.putv("up_abs_tol", 0.0001)
    device.putv("up_rel_tol", 1e-6)
    device.putv("res_abs_tol", 0.0001)
    # Cathode parameters
    device.putv("t_ramp", 1)
    device.putv("t_hold", 5)
    device.putv("squarepulse", False)
    # Save data to matlab file
    device.putv("Datafilename", "")

    BW_arr = []
    for i in range(1, len(format_bias(bias)) + 1, 1):
        print("{0} - RC bandwidth simulation progress {1:f}.\n".format(datetime.datetime.now(),i/(len(format_bias(bias)) + 1)))


        device.putv("min_step", 0.2 / i)
        device.putv("max_step", 0.2 / i)
        V_step = format_bias(bias)[i - 1]
        device.putv("V_step", V_step)
        device.putv("matfilename",
        device.getv("Datafilename") + "V=" + str.replace(str(V_step), ".", "_"))

        # Run the script
        device.feval(os.getcwd() + "/deviceSimul/lumericalScriptFiles/" + "getBW.lsf")
        BW_arr.append(device.getv("BW"))
        if BW_arr[-1] >= 200e9:
            print("Warning : The 3dB bandwidth measured at {0:f} is measured to be {1:f} which is unrealistically large.".format(BW_arr[-1], V_step))

    return np.squeeze(BW_arr)

def format_bias(bias):
    " Format the bias so it gets accepted in Lumerical DEVICE."
    return bias.astype(float)

def measureIV():
    " Measure the IV curve obtained from device. Obsolete"

    # Read through the file
    V, I = [], []
    f = open(os.getcwd()+'/deviceSimul/IVcurve.txt','r')
    for line in f.read().split('\n')[0:-1]:
        V.append(float(line.split('\t')[0]))
        I.append(float(line.split('\t')[1]))
    f.close()

    # Plot the IV curve
    plt.plot(V, I)
    plt.show()

    # Plot the static resistance curve 
    plt.plot(np.array(V)/np.array(I))
    plt.show()

    # Plot the dynamic resistance curve
    dvdi= np.diff(V) / np.diff(I)
    V2 = (np.array(V)[:-1] + np.array(V)[1:]) / 2
    plt.plot(V2, dvdi)
    plt.show()

def neffDatafile(bias, neff, filename):
    """
    ex: neffDatafile([1,2,3,4], [5+1j,6+2j,7+3j,8+4j], os.getcwd() + '/circuitSimul/effective_index.txt')
    """
    with open(filename, "w+") as f:
        #f.writelines("(voltage(a.u.) real(neff)    imag(neff))\n")
        for V, n in zip(bias, neff):
            f.writelines(" {0:f}\t{1:f}\t{2:f}\n".format(V, np.real(n), np.imag(n)))
        print("{0} - Effective index data has been succesfully save in file {1}.\n".format(datetime.datetime.now(), str(filename)))

def RCbandwidthDatafile(bias, bandwidth, filename):
    """
    ex: neffDatafile([0, 1], [100e9, 100e9], os.getcwd() + '/circuitSimul/RCbandwidth.txt')
    """
    with open(filename, "w+") as f:
        #f.writelines("("Voltage","Corresponding 3DB Cut-off Bandwidth")\n")
        for V, BW in zip(bias, bandwidth):
            f.writelines(" {0:f}\t{1:f}\n".format(V, BW))
        print("{0} - RC bandwidth data has been succesfully save in file {1}.\n".format(datetime.datetime.now(), str(filename)))

def writeCSV(bias, resistance, capacitance, filename):
    " Formats the output data in a CSV file. " 

    with open(filename, 'w+') as f:
        f.writelines("reverse-bias voltage [V], series resistance [ohm/m], junction capacitance [F/m]\n" )
        for V, R, C in zip(bias, resistance, capacitance):
            f.writelines("{0:f},{1:f},{2:f}\n".format(V, R, C))
        print("{0} - output data has been saved in a CSV file.".format(datetime.datetime.now()))

def characterizeEO(TDRDatafile):
    ""

    # Open an instance of DEVICE
    # filename = os.getcwd() + "/deviceSimul/" + "junction.ldev"
    with lumapi.DEVICE(hide=hideGUI) as device:
        print("{0} - Opening an instance of Lumerical DEVICE.\n".format(datetime.datetime.now()))

        # Build the junction using the constructor script
        buildJunction(device)
    
        # Import the doping profile
        importTDR(device, TDRDatafile)

        # Measure the series resistance and junction capacitance for the reverse bias applied
        bias = np.linspace(0, 4, 10)
       	resistance = bias
        capacitance = bias

        # Save the data in a CSV file
        csvfilename = os.getcwd() + "/deviceSimul/" + "resistance_capacitance.csv"
        writeCSV(bias, resistance, capacitance, csvfilename)




        # Run the simulation to acquire d_neff and d_alpha and save to file
        #bias = np.linspace(0,4,10)
        #neff = neffSweep(device, 1550e-9, np.linspace(0,4,10), plotResults=True)
        #neffDatafile(bias, neff, os.getcwd() + '/circuitSimul/effective_index.txt')
        
        # Measure bandwidth and put in a textfile
        #bias = np.linspace(0,4,10)
        #bandwidth = bandwidthSweep(device, 1550e-9, bias, plotResults=True)
        #RCbandwidthDatafile(bias, bandwidth, os.getcwd() + '/circuitSimul/RCbandwidth.txt')

        # Measure capacitance
        # Measure resitance
        # Save neff, R and C in a CSV file

    # Close the instance of DEVICE
    device.close()
    print("{} - The opened instance of Lumerical DEVICE has been closed succesfully.\n".format(datetime.datetime.now()))