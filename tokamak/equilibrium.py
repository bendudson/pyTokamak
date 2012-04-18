
class FluxSurface:
    """ Describes a single flux surface 
    and provides some useful operations on it
    
    Data members
    ------------
    
    f    = R * Bt
    
    """

    def __init__(self, f):
        self.f = f
        self._theta = [0,1] # Theta values
        return

    def interpolate(self, var, theta):
        """ Interpolate a variable at a given angle
        
        """
        theta = theta % (2.*pi)  # Get between 0 and 2pi
        
        # Find if a close to a theta grid point
        
        
    def integral(self, var):
        """ Flux surface integral
        
        result = int var dl
        
        where dl is the poloidal arc length
        """
        
    def deriv(self, var):
        """ Flux surface derivative
        
        result = d/dl(var)
        
        result(theta) 

        where dl is the poloidal arc length
        
        """
        
        # Find closest theta below this one
        t0 = theta - 0.1
        # Closest theta above this
        t1 = theta + 0.1
        
        # Arc length between them
        dl = 0.1 ### REPLACE
        
        # Define a 'closure' function to return
        def dvardl(theta):
            return ( var(t1) - var(t0) ) / dl
        
        return dvardl
            

    def average(self, var):
        """ Flux-surface average a variable
        
        Input can be either an array, or a function
        which returns the value as a function of poloidal angle
        
        """
        
        return integral(lambda x: var(x) / Bp(x)) / integral(lambda x: 1./Bp(x))
    
    def R(self, theta):
        """ Major radius in meters as function of angle
        """
        return interpolate(self._R, theta)

    def Z(self, theta):
        """ Height in meters as function of angle
        around flux surface
        """
        return interpolate(self._Z, theta)
        
    def B(self, theta):
        """ B field as function of angle
        """
        return interpolate(self._B, theta)
        
    def Bp(self, theta):
        """ Poloidal field as function of angle
        """
        return interpolate(self._Bp, theta)
    
    def Bt(self, theta):
        """ Toroidal field as function of angle
        """
        return f / R(theta)
    
    ###### Useful functions
    
    def Bsqav(self):
        """ Calculate < B^2 > 
        """
        return average(lambda x: B(x)**2)

    
    def circumference(self):
        """ Circumference of the flux surface in meters
        """
        
