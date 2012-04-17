#
# Neoclassical calculations of viscosities and bootstrap current
#

from species import *
from numpy import zeros, sqrt, pi

def collisionTime(si, sj):
    """ Collision time between two species
    
    """
    eps0   = 8.8542e-12
    coolog = 20. # Coulomb Logarithm
    
    rat  = (4.*pi*eps0)**2 / si.charge
    chdn = sj.density * si.charge
    fac  = si.mass / sj.charge
    return (3./(16.*sqrt(pi)))*(fac**2)*(si.Vth**3)*rat/(chdn*coolog)

def collisionTimes(spec):
    """Calculate collision time for a list of species

    Inputs
    ------

    spec = a list of Species classes
    
    For single ion species Deuterium plasma with Ti=Te, could use
    
    tau = collisionTime(genSpecies(Te, Ne, AA=[None, 2]))

    Return
    ------
    
    Returns collision frequencies
 
    tau_ij
    
    the collision time for species i colliding with species j
    
    In the example above, tau contains

    tau[0,0] = tau_ee
    tau[1,0] = tau_ie
    tau[0,1] = tau_ei
    tau[1,1] = tau_ii
    
    """

    ns = len(spec) # Number of species
    
    # If only one species, collision time with itself
    if ns == 1:
        return collisionTime(spec, spec)
    
    # Array of collision times
    tau = zeros([ns, ns])

    for i in range(ns):
        for j in range(ns):
            tau[i,j] = collisionTime(spec[i], spec[j])
    return tau


def viscosity(spec, conl, ft, vmax=5):
    """
    Neoclassical viscosity calculation
    
    Inputs
    ======

    spec  = list of Species
    conl  = connection length (meters)
    ft    = Trapped fraction
    
    Optional inputs
    ===============
    
    vmax (=5)  
    """
    
    for s in spec:
        vtha = s.Vth # Thermal speed
        omta = vtha / conl
        vsttau = 8./(3.*pi) * ft*omta
    

def bootstrapHS(ft):
    """ Bootstrap current calculation using Hirshman-Sigmar formalism
    from Nucl. Fusion 21 (1981) 1079
    
    Inputs
    ======

    ft   = Trapped fraction
    
    """
    
    if (ft < 0.) or (ft > 1.):
        raise ValueError("Trapped fraction ft must be between 0 and 1")
    
    
