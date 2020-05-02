"""
Directional coupler model that implements a dependancy to the gap.

Author      : Simon Bélanger-de Villers (simon.belanger-de-villers.1@ulaval.ca)
Created     : 2018
Last edited : October 17th 2019
"""
import math
import matplotlib.pyplot as plt
import numpy as np

class DCMeasured(object):
    """
        This class models a directional coupler (point coupling) where the only parameter is
        the gap between the 2 adjacent waveguides. As such     
    
            kappa(gap) = a * exp( 1/b * gap)

        where 
            gap     : is the gap between the 2 waveguides [nm].
            kappa   : field cross-coupling coefficient [].


            a : field cross-coupling coefficient for a gap of 0 nm [].
            b : rate of decay of the cross coupling coefficient vs gap [1/nm]
    """
    a = 0.3  # Between 0 and 1
    b = 900  # Between 0 and inf
    
    def __init__(self, gap):
        self.gap = gap

    def getCouplingCoefficient(self):
        return self.measureCouplingCoefficient(self.a, self.b, self.gap)

    @staticmethod
    def measureCouplingCoefficient(a, b, gap):
        return a * math.exp( -1/b * gap)

    def plotKappa(self):
        " Plot the relationship between the coupling coefficient and the gap. "

        gap, kappa     = np.linspace(200,1000,100), [] # Arbitrary but realistic
        for g in gap:
            kappa.append(self.measureCouplingCoefficient(self.a, self.b, g))

        plt.plot(gap, kappa)
        plt.xlabel('Gap [nm]')
        plt.ylabel('Field cross-coupling coefficient')
        plt.show()

if __name__ == '__main__':

    dc = DCMeasured(300)
    dc.b = 100
    dc.plotKappa()
    print(dc.getCouplingCoefficient())
