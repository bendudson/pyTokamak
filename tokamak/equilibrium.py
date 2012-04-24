#
# 
# Requires the NumPy library
#

from numpy import size, copy, sqrt, linspace, pi, zeros
from utils import integrate, interpPeriodic

class FluxSurface:
    """ Describes a single flux surface 
    and provides some useful operations on it
    
    Public data members
    -------------------
    
    f    = R * Bt
    
    Public functions
    ----------------
    
    These are created by __init__
    
    R(theta)    Major radius in meters
    Z(theta)    Height in meters
    B(theta)    Magnetic field in Teslas
    Bp(theta)   Poloidal magnetic field
    Bt(theta)   Toroidal magnetic field
    
    Private data members
    --------------------
    
    Implementation details which may change
    
    _R, _Z    = arrays of (R,Z) coordinates  [m]
    
    _Bp       = Poloidal field [T]
    _Bt       = Toroidal field [T]
    _B        = Field magnitude [T]
    
    _ntheta   = Number of points in the arrays
    _theta    = array of theta values on each point
    
    """

    def __init__(self, f, r, z, Bp, theta=None):
        """ Initialise
        """
        self.f = f
        
        # Check sizes of input arrays
        self._ntheta = size(r)
        if (size(z) != self._ntheta) or (size(Bp) != self._ntheta):
            raise ValueError("r, z, and Bp arrays must have the same size")
        
        # Copy input arrays to prevent modification from outside
        self._R = copy(r)
        self._Z = copy(z)
        self._Bp = copy(Bp)
        
        # Calculate other components of B
        self._Bt = f / r
        self._B   = sqrt(self._Bt**2 + self._Bp**2)
        
        if theta == None:
            # Define theta to be equally spaced (arbitrary) in range ]0,2pi]
            self._theta = linspace(0.0, 2.*pi, num=self._ntheta, endpoint=False)
        else:
            if size(theta) != self._ntheta:
                raise ValueError("theta must have the same number of points as r,z and Bp")
            self._theta = theta
        
        self._dldt = zeros(self._ntheta)
        dtheta = 2.*pi / self._ntheta
        for i in range(self._ntheta):
            drdt = (r[ (i+1) % self._ntheta ] - r[ (i-1) % self._ntheta ]) / (2.*dtheta)
            dzdt = (z[ (i+1) % self._ntheta ] - z[ (i-1) % self._ntheta ]) / (2.*dtheta)
            self._dldt[i] = sqrt(drdt**2 + dzdt**2)

        # Create functions for external access to values
        self.R = interpPeriodic(self._theta, self._R, copy=False)
        self.Z = interpPeriodic(self._theta, self._Z, copy=False)
        self.B = interpPeriodic(self._theta, self._B, copy=False)
        self.Bp = interpPeriodic(self._theta, self._Bp, copy=False)
        self.Bt = interpPeriodic(self._theta, self._Bt, copy=False)
        
        self._dldtheta = interpPeriodic(self._theta, self._dldt, copy=False)

    def max(self, var):
        """ Maximum value of a variable over theta 
        
        Parameters
        ----------
        
        var(theta) is a callable object
        
        """
        # Only sample at the grid points
        return max([var(theta) for theta in self._theta])
    
    def integral(self, var):
        """ Flux surface integral
        
        result = int var dl
        
        where dl is the poloidal arc length
        """
        
        # Integrate  ( var * dl/dtheta ) dtheta
        return integrate( lambda x: var(x) * self._dldtheta(x), 0, 2*pi )
        
    def deriv(self, var, dtheta = 0.01):
        """ Flux surface derivative
        
        result = d/dl(var)
        
        result(theta) 

        where dl is the poloidal arc length
        
        """
        
        # Define a function which returns derivative
        # of var at theta, using central differencing
        def dvardl(theta):
            dvdtheta = (var(theta + dtheta) - var(theta - dtheta))/(2.*dtheta)
            return dvdtheta / self._dldtheta(theta)
        
        # Return this new function 
        return dvardl
    
    def average(self, var):
        """ Flux-surface average a variable
        
        Input can be either an array, or a function
        which returns the value as a function of poloidal angle
        
        """
        
        return self.integral(lambda x: var(x) / self.Bp(x)) / self.integral(lambda x: 1./self.Bp(x))
        
    ###### Useful functions
    
    def Bsqav(self):
        """ Calculate < B^2 > 
        """
        return self.average(lambda x: self.B(x)**2)
    

from copy import deepcopy
from numpy import searchsorted

class Equilibrium:
    """ Represents an axisymmetric tokamak equilibrium 
    
    """
    def __init__(self, fix=None):
        """ Construct an Equilibrium object 
        
        Parameters
        ----------
        
        If no inputs are given, an empty equilibrium is created
        which can be added to afterwards
        
        fix = Dictionary of values describing fixed-boundary solution
              on flux surfaces. Assumed to contain the following keys:
              
            'npsi'   Number of points in poloidal flux psi
            'npol'   Number of poloidal points
            
            'psi'    Poloidal flux
        """
        
        if fix != None:
            # Check that fix has the required keys
            required = ['npsi', 'npol', 'psi', 'f(psi)', 'R', 'Z', 'Bp']
            
            for r in required:
                if not fix.has_key(r):
                    raise ValueError("Missing key: " + r)
            
            if not fix.has_key("psinorm"):
                # Add normalised psi, assuming core to edge
                psi = fix['psi']
                fix['psinorm'] = (psi - psi[0]) / (psi[-1] - psi[0])

            # Make a deep copy of the data so it can't be modified afterwards
            # in unpredictable ways
            self.fix = deepcopy(fix);
        else:
            self.fix = None
        
    def getFluxSurface(self, psi):
        """ Return a FluxSurface object at a given psi 
        
        Will be created by interpolation if needed
        
        Paramaters
        ----------
        
        psi     =  Normalised poloidal flux
        
        """
        
        if self.fix != None:
            # Find the psi index this value comes before
            ind = searchsorted(self.fix['psinorm'], psi)
            psiarr = self.fix['psinorm']
            
            if (ind == 0) or (ind == self.fix['npsi']):
                raise ValueError("normalised psi value %e out of range %e to %e" % (psi, psiarr[0], psiarr[-1]))
            
            
            # Indices
            im = ind-1
            ip = ind
            
            # Weights for interpolation
            wp = (psi - psiarr[im]) / (psiarr[ip] - psiarr[im])
            wm = 1. - wp
            
            # Interpolate
            def interp1d(var):
                return wp*var[ip] + wm*var[im]
            def interp2d(var):
                return wp*var[ip,:] + wm*var[im,:]
            
            f = interp1d(self.fix['f(psi)'])
            R = interp2d(self.fix['R'])
            Z = interp2d(self.fix['Z'])
            Bp = interp2d(self.fix['Bp'])
            
            return FluxSurface(f, R, Z, Bp)
        
        # No data
        return None
    
    
    def surfaces(self, psimin=0.0, psimax=1.0, n=None):
        """ Iterate over flux surfaces 
        
        Keywords
        --------
        
        psimin (0.0)   = Minimum normalised poloidal flux. 0 = core
        psimax (1.0)   = Maximum normalised poloidal flux. 1 = edge
        
        Example
        -------
        
        for s in eq.surfaces():
           <code>
        """
        
        if self.fix != None:
            # Got flux-surface data
            
            ind0 = searchsorted(self.fix['psinorm'], psimin)
            ind1 = searchsorted(self.fix['psinorm'], psimax)

            # Create a generator which will iterate over surfaces
            def surfgen():
                for i in range(ind0, ind1):
                    f = (self.fix['f(psi)'])[i]
                    R = (self.fix['R'])[i,:]
                    Z = (self.fix['Z'])[i,:]
                    Bp = (self.fix['Bp'])[i,:]
                    yield FluxSurface(f, R, Z, Bp)
            
            return surfgen()
        
        # No data
        return None
