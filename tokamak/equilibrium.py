#
# 
# Requires the NumPy library
#

from numpy import size, copy, sqrt, linspace, pi, zeros
from utils import integrate, interpPeriodic, interp1d

import species

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
    
    
    Public functions 
    ----------------
    
    pressure(psi)    Pressure [Pa]
    dens(psi)        Density [m^-3]
    temp(psi)        Temperature [eV]
    
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
            required = ['npsi', 'npol', 'psi', 'f(psi)', 'p', 'R', 'Z', 'Bp']
            
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
            
            # Create a function to return the pressure
            self.pressure = interp1d(self.fix['psinorm'], self.fix['p'], copy=False)

            # See if density, temperature profiles also set
            if fix.has_key("ne"):
                self.setDensity(self.fix["ne"])
        else:
            self.fix = None
    
    def ddpsi(self, f, dpsi=0.01):
        """ Calculate psi derivative 
        
        """
        psi = self.fix['psi']
        pnorm = psi[-1] - psi[0]  # normalised psi = psi / pnorm

        def dfdpsi(psi):
            fp = psi + dpsi
            if fp > 1.0:
                fp = 1.0
            fm = psi - dpsi
            if fp < 0.:
                fp = 0.
            return ( f(fp) - f(fm) ) / (fp - fm) / pnorm
        return dfdpsi
    
    def setDensity(self, dens, psi=None):
        """ Sets the density profile
        
        Parameters
        ----------
        
        dens - Density in m^-3
               Can be:

               1. A function taking normalised psi and returning density

                  ne = dens(psinorm)

               2. An array of values on psi grid (uniform between 0 and 1 if not)
               
               3. A constant

        psi  - Optional array of normalised psi values
        
        """
        
        # Check if we can call it to get a value
        try:
            val = dens(0.5);
            densfunc = dens; # Use this new function
        except:
            # Not a function type
            # Check if it's an array type
            try:
                val = dens[0]
                p = psi
                if p == None:
                    p = linspace(0.0, 1.0, num=len(dens), endpoint=False)
                densfunc = interp1d(p, dens)
            except:
                # Should be a constant
                val = dens
                densfunc = lambda x: dens
                
        # val should now contain a number
        if size(val) > 1:
            raise ValueError("dens argument doesn't yield a single value")
        try:
            test = 2.*val + val*val        
        except:
            raise ValueError("dens argument must give a numerical value")
        
        # Checks passed, so can set the density function
        self.dens = densfunc

        # Set temperature function to use the density, assuming Ti = Te
        self.temp = lambda x: self.pressure(x) / (2. * 1.602e-19 * self.dens(x) )
        
    def setTemperature(self, T, psi=None):
        """ Sets the temperature profile
        
        Parameters
        ----------
        
        T    - temperature in eV
               Can be:

               1. A function taking normalised psi and returning density

                  ne = T(psinorm)

               2. An array of values on psi grid (uniform between 0 and 1 if not)
               
               3. A constant

        psi  - Optional array of normalised psi values
        
        """
        
        # Check if we can call it to get a value
        try:
            val = T(0.5);
            Tfunc = T; # Use this new function
        except:
            # Not a function type
            # Check if it's an array type
            try:
                val = T[0]
                p = psi
                if p == None:
                    p = linspace(0.0, 1.0, num=len(T), endpoint=False)
                Tfunc = interp1d(p, T)
            except:
                # Should be a constant
                val = T
                Tfunc = lambda x: T
                
        # val should now contain a number
        if size(val) > 1:
            raise ValueError("T argument doesn't yield a single value")
        try:
            test = 2.*val + val*val
        except:
            raise ValueError("T argument must give a numerical value")
        
        # Checks passed, so can set the density function
        self.temp = Tfunc

        # Set density function to use the temperature
        self.dens = lambda x: self.pressure(x) / (2. * 1.602e-19 * self.temp(x) )

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
            
            if ind == 0:
                ind = 1
                
            if (ind == self.fix['npsi']):
                raise ValueError("normalised psi value %e out of range %e to %e" % (psi, psiarr[0], psiarr[-1]))
            
            
            # Indices
            im = ind-1
            ip = ind
            
            # Weights for interpolation
            wp = (psi - psiarr[im]) / (psiarr[ip] - psiarr[im])
            wm = 1. - wp
            
            # Interpolate
            def inter1d(var):
                return wp*var[ip] + wm*var[im]
            def inter2d(var):
                return wp*var[ip,:] + wm*var[im,:]
            
            f = inter1d(self.fix['f(psi)'])
            R = inter2d(self.fix['R'])
            Z = inter2d(self.fix['Z'])
            Bp = inter2d(self.fix['Bp'])
            
            # Create the flux surface
            f = FluxSurface(f, R, Z, Bp)
        
            f.psinorm = psi

            # Create some species
            try:
                # Try to get profiles
                d = self.dens(psi)
                T = self.temp(psi)
                dndpsi = self.ddpsi(self.dens)(psi)
                dTdpsi = self.ddpsi(self.temp)(psi)
                f.species = species.genSpecies(T, d, AA=[None, 2], 
                                               dTdpsi=dTdpsi, dndpsi=dndpsi)
            except:
                print "Warning: No species information"
            return f
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
                    f = FluxSurface(f, R, Z, Bp)
                    try:
                        # Try to get profiles
                        psi = (self.fix['psinorm'])[i]
                        d = self.dens(psi)
                        T = self.temp(psi)
                        dndpsi = self.ddpsi(self.dens)(psi)
                        dTdpsi = self.ddpsi(self.temp)(psi)
                        f.species = species.genSpecies(T, d, AA=[None, 2], 
                                                       dTdpsi=dTdpsi, dndpsi=dndpsi)
                        f.psinorm = psi
                    except:
                        raise
                    yield f
            
            return surfgen()
        
        # No data
        return None
