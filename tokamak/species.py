from numpy import size, sqrt

class Species:
    """Class to represent a particle species
    
    Members
    -------

    T          = Temperature in eV
    density    = Density in m^-3
    atomicMass = Mass in atomic mass units (proton = 1)
    unitCharge = Charge in elementary units (proton = 1, electron = -1)
    mass       = mass in kg
    charge     = charge in Coulombs
    
    vth        = Thermal speed (m/s)
    """
    def __init__(self, T, n, AA=None, ZZ=1):
        # Set basic quantities
        self.T = T
        self.density = n
        if AA == None:
            self.atomicMass = 1./1860.  # Electron mass in atomic mass units
            self.unitCharge = -1      # Charge in elementary charge units
        else:
            self.atomicMass = AA
            self.unitCharge = ZZ
        
        # Derived quantities
        self.mass   = self.atomicMass * 1.67262158e-27
        self.charge = self.unitCharge * 1.602e-19
        self.Vth    = sqrt(2.*1.602e-19 * self.T / self.mass)

def genSpecies(T, n, AA=None, ZZ=1):
    """Create a list of species classes
    
    Inputs
    ------

    T  =  Temperature (eV, scalar or array)
    n  =  density   (m^-3, scalar or array)
    AA  = atomic number (scalar or array). None means electron
    ZZ  = charge (scalar or array)
    """
    
    # Get number of species
    ns = max(size(T), size(n), size(AA), size(ZZ)) # Number of ion species
    
    # Check sizes of arrays
    if (size(T) != 1) and (size(T) != ns):
        raise ValueError("first argument T must have either 1 or "+str(ns)+" elements")
    if (size(n) != 1) and (size(n) != ns):
        raise ValueError("second argument n must have either 1 or "+str(ns)+" elements")
    if (size(AA) != 1) and (size(AA) != ns):
        raise ValueError("Keyword AA must have either 1 or "+str(ns)+" elements")
    if (size(ZZ) != 1) and (size(ZZ) != ns):
        raise ValueError("Keyword AA must have either 1 or "+str(ns)+" elements")
    
    s = []
    for i in range(ns):
        Ti = T[i] if size(T) != 1 else T
        ni = n[i] if size(n) != 1 else n
        Ai = AA[i] if size(AA) != 1 else AA
        Zi = ZZ[i] if size(ZZ) != 1 else ZZ
        s.append(Species(Ti, ni, Ai, Zi))
    return s
