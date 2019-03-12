import numpy as np
import matplotlib.pyplot as plt
import math


# Important constants [global]
c = 299792458 # Velocity of light in vacuum [m/s]

def obtenirdata():
    """

    Open the datafile and extract the different arrays to a dict.

    :return:
    """

    file = open('WG_RW=500nm_RH=220nm_SH=0nm_R=2,5um.txt', 'r')

    # Read the first line and extract keys
    keys = file.readline().replace(' \n','').split(', ')

    # Create a dict from the different keys
    data = {keys[0]: [], keys[1]: [], keys[2]: [], keys[3]: [], keys[4]: [], keys[5]: [], keys[6]: [], keys[7]: [], keys[8]: []}

    # Read the remaining of the file and extract data, put data in different lists

    for dat_line in file.read().split(' \n'):
        linesplit = dat_line.split(', ')

        if (len(linesplit) == 9):
            for i in range (0,len(linesplit)):
                data[keys[i]].append(float(linesplit[i]))

    # Remove the extra values for f_ng, ng, f_D and D
    data['f_ng'].pop()
    data['ng'].pop()
    data['f_D'].pop();data['f_D'].pop()
    data['D'].pop();data['D'].pop()

    file.close()
    return data


def extractcompactmodel(data):

    """
    Extract parameters for compact model, ng, D.


    lambda_c = central wavelength
    """


    # Find the value of n_eff, n_g and D at lambda_c
    lambda_c = c/data['f'][5]*1e9
    n_eff = data['real(neff)'][5]
    n_g = data['ng'][5]
    D = data['D'][5]



    return lambda_c, n_eff, n_g, D

def ang_freq(lambda_i):
    return (2 * math.pi * c) / lambda_i

def buildcompactmodel(n_eff, n_g, D, lambda_c, lambda_0):

    omega = ang_freq(lambda_0)
    omega_c = ang_freq(lambda_c)
    d_omega = omega - omega_c
    beta_0 = (omega_c * n_eff) / c
    beta_1 = n_g / c
    beta_2 = -(D * lambda_c**2)/(2 * math.pi * c)
    beta = beta_0 + beta_1 * d_omega + 1/2 * beta_2 * d_omega**2
    return beta

def plotCM(data, neff, ng, D, lambda_c):

    f = data['f']
    beta_CM = []
    for i in range(0,len(f)):
        lambda_0 = c/f[i]
        beta_CM.append(buildcompactmodel(neff, ng, D, lambda_c, lambda_0))

    plt.plot(f, beta_CM)
    #plt.show()

def plotbeta(data):

    plt.plot(data['f'], data['real(beta)'], label='Re(beta)')
    #plt.plot(data['f'], data['imag(beta)'], label='Im(beta)')
    #plt.legend()
    plt.show()


def comparebeta2compactmodel():

    pass


plotCM(obtenirdata(), 2.44685, 4.06793, 0.00103741, 1553.4253839616972e-9)
plotbeta(obtenirdata())


