#
# Some useful maths routines
#
# If available, imports functions from standard library or SciPy.
# Otherwise defines a functional alternative
# 
# Defines the following functions:
# 
#   y = erf(x)     Error function
#   y = integrate( f, a, b )  Integrate a function f from a to b
#   
#   f = interpPeriod( x, y, copy=True, period=2.*pi)  Return a callable object

try:
    # erf function available from Python 3.2
    from math import erf
except ImportError:
    # Try SciPy
    try:
        from scipy.special import erf
    except ImportError:
        # No erf function, so define
        def erf(x):
            # save the sign of x
            sign = 1 if x >= 0 else -1
            x = abs(x)

            # constants
            a1 =  0.254829592
            a2 = -0.284496736
            a3 =  1.421413741
            a4 = -1.453152027
            a5 =  1.061405429
            p  =  0.3275911

            # A&S formula 7.1.26
            t = 1.0/(1.0 + p*x)
            y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x)
            return sign*y # erf(-x) = -erf(x)

try:
    import scipy.integrate.quad as integrate
    #print("Using SciPy integrate.quad")
except ImportError:
    # Define simple integrate function
    def integrate( f, a, b, n=100 ):
        """Approximate the definite integral of f from a to b by Simpson's rule"""
        
        if n % 2 != 0:
            raise ValueError("n must be even!")
        
        h  = (float(b) - a)/n
        si = 0.0
        sp = 0.0
    
        for i in range(1, n, 2):
            xk = a + i*h
            si += f(xk)
            
        for i in range(2, n, 2):
            xk = a + i*h
            sp += f(xk)
            
        s = 2*sp + 4*si + f(a) + f(b)
        return (h/3)*s
    #print("Using Simpson's rule for integration. Install SciPy for better method")


from numpy import pi, copy, array, searchsorted, size, where

class interpPeriodic:
    def __init__(self, x, y, copy=True, period=2.*pi):
        """
        """
        
        # Store x & y, either just references, or a copy of the data
        self._x = array(x, copy=copy)
        self._y = array(y, copy=copy)
        self._period = period
        
    def __call__(self, x):
        """ Interpolate 
        
        """
        
        new_x = x % self._period

        # Find the index this value comes before
        new_ind = searchsorted(self._x, new_x) % len(self._x)
        
        # Simple linear interpolation for now
        
        im  = (new_ind - 1) % len(self._x)
        ip  = new_ind
        
        xm = self._x[im]
        xp = self._x[ip]

        xm    = where(xm < xp, xm, xm - self._period)
        new_x = where(new_x > xp, new_x - self._period, new_x)
        
        ym = self._y[im]
        yp = self._y[ip]
        
        slope = (yp - ym) / (xp - xm)
        
        return ym + slope * (new_x - xm)
