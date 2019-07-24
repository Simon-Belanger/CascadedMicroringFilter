from math import pi, exp, sin, sqrt

class DC_various(object):
    """
    DC made with 

    Design Space Exploration of Microring Resonators in 
    Silicon Photonic Interconnects : Impact of the ring curvature
    Maisam Bahadori
    JLT, 2018
    """

    def __init__(self):

        # Define parameters for a 400x220 nm waveguide
        self.aE = 0.242422
        self.aO = 0.077526
        self.gE = 0.010687
        self.gO = 0.006129

        self.d  = 200e-9

        self.L = 10e-6

        pass

    def kappa(self, wavelength):
        """ Cross-coupling coefficient of the DC. """
        
        # Obtain curvature function and simplify the term by using B [ B(x) ]
        B = self.curvatureFunction("Straight")

        # Obtain x parameters 
        xE = self.gE * self.L 
        xO = self.gO * self.L

        # Compute the coupling coefficient
        kappa = sin(pi/wavelength * (self.aE/self.gE * exp(-self.gE * self.d) * B(xE) + self.aO/self.gO * exp(-self.gO * self.d) * B(xO)))
        return kappa

    def curvatureFunction(self, geometry):
        """ Returns the curvature function for the given geometry 
            which is an anonymous function that depends on x. 
            """

        if (geometry in "Straight"):
            return lambda x: x
        elif (geometry in "Straight + Bends"):
            pass
        elif (geometry in "Bus - Ring"):
            return lambda x: sqrt(2 * pi * x)
        elif (geometry in "Bus - Racetrack"):
            pass
        elif (geometry in "Ring - Ring"):
            pass
        elif (geometry in "Bent"):
            pass
        else:
            print("Error : Geometry not implemented")

    def curvatureParameter(self, geometry):
        """ Returns the expression of x depending on the parameters. 
            """

if __name__ == '__main__':

    k = DC_various()
    print(k.kappa(1500e-9))