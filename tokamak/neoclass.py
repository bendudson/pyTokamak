#
# Neoclassical calculations of viscosities and bootstrap current
#
# Adapted from H.Wilson's SCENE code. Parts just translated straight from
# Fortran, so could do with some Pythonifying
#
# Ben Dudson, University of York
#

from species import Species, genSpecies
from numpy import zeros, sqrt, pi, linalg, size
from math import exp

from utils import erf, integrate

def collisionTime(si, sj):
    """ Collision time between two species
    
    Parameters
    ----------
    
    si and sj should have the following members:

    s.charge   = charge in Coulombs
    s.mass     = mass in kg
    s.density  = density in m^-3
    s.Vth      = thermal speed in m/s
    
    Returns
    -------
    
    Returns collision time in seconds
    
    """
    eps0   = 8.8542e-12
    coolog = 20. # Coulomb Logarithm
    
    #coolog = 24. - log(sqrt(ne*1.0d-6)*bk/te)

    rat  = (4.*pi*eps0)**2 / si.charge
    chdn = sj.density * si.charge
    fac  = si.mass / sj.charge
    return (3./(16.*sqrt(pi)))*(fac**2)*(si.Vth**3)*rat/(chdn*coolog)

def collisionTimes(spec):
    """Calculate collision time for a list of species

    Parameters
    ----------

    spec = a list of Species classes, assumed to have the following members:
        .charge   = charge in Coulombs
        .mass     = mass in kg
        .density  = density in m^-3
        .Vth      = thermal speed in m/s

    Returns
    -------
    
    Returns collision times in seconds
 
    tau_ij
    
    the collision time for species i colliding with species j
    
    Example
    -------

    For single ion species Deuterium plasma with Ti=Te, could use
    
    tau = collisionTime(genSpecies(Te, Ne, AA=[None, 2]))

    tau[0,0] = tau_ee
    tau[1,0] = tau_ie
    tau[0,1] = tau_ei
    tau[1,1] = tau_ii
    
    """

    ns = size(spec) # Number of species
    
    # If only one species, collision time with itself
    if ns == 1:
        return collisionTime(spec, spec)
    
    # Array of collision times
    tau = zeros([ns, ns])

    for i, si in enumerate(spec):
        for j, sj in enumerate(spec):
            tau[i,j] = collisionTime(si, sj)
    return tau
    
def friction(spec):
    """ Calculate friction coefficients
    
    Parameters
    ----------
    
    spec  = list of Species with members:
        .charge   = charge in Coulombs
        .mass     = mass in kg
        .density  = density in m^-3
        .Vth      = thermal speed in m/s
    
    Returns
    -------
    
    Returns a class with data members
    
    rl11
    rl12
    rl21
    rl22
    
    Example
    -------
    
    For single ion species Deuterium plasma with Ti=Te, try
    
    from tokamak import neoclass
    r = neoclass.friction(neoclass.genSpecies(1000, 1e19, AA=[None, 2]))
    
    """
    
    ns = size(spec) # Number of species
    
    sum00 = zeros(ns)
    sum10 = zeros(ns)
    sum01 = zeros(ns)
    sum11 = zeros(ns)
    
    for i, si in enumerate(spec):
        for j, sj in enumerate(spec):
            tau_ij = collisionTime(si, sj) # Collision time
            xab = sj.Vth/si.Vth            # Ratio of thermal speeds xab(i,j)
            
            rmat00 = -(1. + si.mass / sj.mass) / (1. + xab**2)**1.5
            sum00[i] += si.density*si.mass*rmat00/tau_ij
            
            rmat01 = -1.5*(1. + si.mass/sj.mass)/(1. + xab**2)**2.5
            sum01[i] += si.density*si.mass*rmat01/tau_ij
            
            rmat11 = -(3.25 + 4.*xab**2 + 7.5*xab**4) / (1. + xab**2)**2.5
            sum11[i] += si.density*si.mass * rmat11/tau_ij
        
        sum10[i] = sum10[i]
    
    rl11 = zeros([ns, ns])
    rl12 = zeros([ns, ns])
    rl21 = zeros([ns, ns])
    rl22 = zeros([ns, ns])

    for i, si in enumerate(spec):
        for j, sj in enumerate(spec):
            
            tau_ij = collisionTime(si, sj) # Calculate collision time
            xab = sj.Vth/si.Vth  # Ratio of thermal speeds xab(i,j)
            
            rnat00 = (1. + si.mass/sj.mass)/(1. + xab**2)**1.5
            rl11[i,j] = si.density*si.mass*rnat00 / tau_ij
            if i == j:
                rl11[i,j] += sum00[i]
                
            rnat10 = 1.5*(1. + si.mass/sj.mass)/(1. + xab**2)**2.5
            rl21[i,j] = si.density*si.mass * rnat10/tau_ij
            if i == j:
                rl21[i,j] += sum10[i]
                
            rnat01 = (si.T / (sj.T * xab))*1.5*(1.+sj.mass/si.mass)/(1.+xab**2)**2.5
            rl12[i,j] = si.density*si.mass*rnat01/tau_ij
            if i == j:
                rl12[i,j] += sum01[i]
                
            rnat11 = 6.75*(si.T/sj.T)*xab**2/(1.+xab**2)**2.5
            rl22[i,j] = si.density*si.mass*rnat11/tau_ij
            if i ==j:
                rl22[i,j] += sum11[i]
            
    # Create a structure to hold the results
    class FrictionData:
        """ Neoclassical friction coefficients
        
        Members
        -------
        
        rl11
        rl12
        rl21
        rl22
        
        """
        pass
    
    result = FrictionData()
    result.rl11 = rl11
    result.rl21 = rl21
    result.rl12 = rl12
    result.rl22 = rl22
    
    return result

def chandrasekhar(y):
    """ Chandrasekhar function for dynamical friction
    
    """
    
    if abs(y) < 1e-3:
        return  (2./sqrt(pi)) * y/3.
    
    return (erf(y) - (2./sqrt(pi))*y*exp(-y*y)) / (2.*y*y)
   

def phig(y):
    """
    Modified error function
    
    erf(y) - chandrasekhar(y)
    
    """
    return erf(y) - chandrasekhar(y)


def trappedFraction(surf):
    """ Calculate the trapped particle fraction
    
    Parameters
    ----------
    
    surf  = Flux surface object
        .max( f(theta) )          Maximum value over theta
        .B(theta)                 Magnetic field strength [T]
        .Bsqav()                  < B^2 >
        .average( func(theta) )   Flux surface average
    """
    try:
        # See if surface already has a trapped fraction value
        ft = surf.trappedFraction
        # If so, return it
        return ft
    except:
        pass # Just catch the error
    
    try:
        if surf.psinorm < 1e-2:
            surf.trappedFraction = 0.
            return 0.
    except:
        pass

    bigint =  integrate( lambda rla: 
                         rla / surf.average( lambda x: sqrt(1. - rla*surf.B(x)) ),
                         0.0,   # Lower limit
                         1./surf.max(surf.B))  # Upper limit
    
    # Set a variable in surface object so we don't have to calculate all that again
    surf.trappedFraction = 1. - 3.*surf.Bsqav() * bigint / 4.
    
    # Sanity check the value
    if (surf.trappedFraction < 0.) or (surf.trappedFraction > 1.):
        raise ValueError("Trapped fraction ft must be between 0 and 1")
    
    return surf.trappedFraction
    
def connectionLength(surf):
    """ Flux surface connection length, the 
    distance along a field-line to complete one poloidal turn.
    In the limit of large aspect-ratio, this becomes qR
    
    Calculates poloidal loop integral
    
    cl = integral  B / Bp  dl
    
    Parameters
    ----------
    
    surf   = Flux surface class
        .B(theta)    Total field [T]
        .Bp(theta)   Poloidal field [T]
        .integral( func(theta) )   Poloidal integral
    
    An additional field
        .connectionLength

    is checked and returned if it exists. Once connection
    length is calculated, this field is set. This means
    that connectionLength is only calculated once per surf object
    
    Returns
    -------
    
    Connection length in meters
    
    """
    try:
         # See if surface already has a value
        conl = surf.connectionLength
        return conl
    except:
        pass
    
    surf.connectionLength = surf.integral(lambda x: surf.B(x) / surf.Bp(x)) / (2.*pi)
    
    return surf.connectionLength

def viscosity(surf, spec, vmax=5., nv=500, colldamping=0.):
    """
    Neoclassical viscosity calculation
    
    Parameters
    ----------

    surf  = Flux surface object (see equilibrium.py)
        .aspectRatio    = aspect ratio a/R
        .Bsqav() # < B^2 > flux-surface average
        .average( func(theta) )
        .deriv( func(theta) )
        .B(theta)
        .Bp(theta)
        
    spec  = list of Species
    
    Optional inputs
    ---------------
    
    vmax (=5.)  Maximum velocity to consider (units of vth)
    nv (=500)   Number of velocity grid points
    
    colldamping (=0.) determines collisional suppression:
       ]0,1[   Collisional suppression
       None    Shiang's result that a=1 experiences no collisional suppression
    
    """
    
    conl = connectionLength(surf)
    ft   = trappedFraction(surf)

    Bsq = surf.Bsqav() # < B^2 > flux-surface average
    
    # < (dB/dl)^2 Bp^2 / B^2 >
    dotav = surf.average(lambda x: ( surf.deriv(surf.B)(x) * surf.Bp(x) / surf.B(x) )**2 )
    
    ns = size(spec) # Number of species

    nv = 500 # Number of velocity points
    dv=vmax/(nv-1.)
    
    vd   = zeros([ns, nv])
    vt   = zeros([ns, nv])
    vtot = zeros([ns, nv])
    
    # Calculate collision times
    coltau = collisionTimes(spec)

    # Populate velocity mesh values
    # This code really needs some love
    
    for i, si in enumerate(spec):
        vtha = si.Vth # Thermal speed
        omta = vtha / conl
        vsttau = 8./(3.*pi) * ft*omta*Bsq / (dotav * vtha**2)
        
        for n in range(nv):
            vv = (n+0.5)*dv
            
            for j, sj in enumerate(spec):
                y    = vv * si.Vth / sj.Vth
                
                vdij = phig(y) / (coltau[i,j]*vv**3)
                vd[i,n] += vdij
                
                vtij = (((erf(y) - 3.*chandrasekhar(y)) / vv**3) + 4.*(si.T/sj.T + (si.Vth/sj.Vth)**2)*chandrasekhar(y)/vv)/coltau[i,j]
                vt[i,n] += vtij
            
            vd[i,n] *= 3.*sqrt(pi)/4.
            vt[i,n] *= 3.*sqrt(pi)/4.
            
            if colldamping == None:
                # Shiang's result
                fac1=1.+(1.- surf.aspectRatio**2)*vsttau*vd[i,n]/vv
                fac2=1.+5.*pi*vt[i,n]/(8.*vv*omta)
            else:
                fac1=1.+colldamping*vsttau*vd[i,n]/vv
                fac2=1.+colldamping*5.*pi*vt[i,n]/(8.*vv*omta)
            vtot[i,n] = vd[i,n] / (fac1 * fac2)
    
    # perform velocity integrals to calculate k's
    # No need to store vd,vt,vtot arrays
    rk11 = zeros(ns)
    rk12 = zeros(ns)
    rk22 = zeros(ns)

    vis = zeros([ns, 3])
    
    for i, si in enumerate(spec):
        for n in range(nv):
            vv = (n+0.5)*dv
            arg=(vv**4)*exp(-vv**2)*vtot[i,n]*dv
            rk11[i] += arg
            rk12[i] += arg*vv**2
            rk22[i] += arg*vv**4
        
        rk11[i] *= ft*8./(3.*(1.-ft)*sqrt(pi))
        rk12[i] *= ft*8./(3.*(1.-ft)*sqrt(pi))
        rk22[i] *= ft*8./(3.*(1.-ft)*sqrt(pi))

        print("viscosity: ", rk11[i], rk12[i], rk22[i])

        vis[i,0] = rk11[i]*si.density*si.mass
        vis[i,1] = (rk12[i] - 2.5*rk11[i])*si.density*si.mass
        vis[i,2] = (rk22[i] - 5.*rk12[i] + 6.25*rk11[i])*si.density*si.mass
        
    return vis

def collisionality(surf, spec=None):
    """ Calculate collisionality for a collection of species 
    on a flux surface.
    
    """
    if spec == None:
        try:
            spec = surf.species
        except:
            raise ValueError("Must specify species information")
        
    # Get collision frequencies
    tau = collisionTimes(spec)

    ns = len(spec) # Number of species

    # Calculate collision frequencies for each species 
    nu = zeros(ns)
    for i in range(ns):
        nu[i] = sum(1./tau[i,:]) # Sum collision frequency with each species
    
    # Calculate banana bounce frequency
    length = connectionLength(surf)
    eps = surf.eps()
    wb = zeros(ns)
    for i in range(ns):
        wb[i] = sqrt(2.*eps)*spec[i].Vth / length
    
    return nu / (eps * wb)

def bootstrapHS(surf, spec=None):
    """ Bootstrap current calculation using Hirshman-Sigmar formalism
    from Nucl. Fusion 21 (1981) 1079
    
    Parameters
    ----------
    
    surf = flux surface
        .f         = R*Bt
        .Bsqav()   = <B^2> flux-surface average
    
    spec = List of species objects
        .T       = Temperature in eV
        .dTdpsi  = Derivative of T w.r.t poloidal flux
        .density = Density in m^-3
        .dndpsi  = Derivative of density w.r.t poloidal flux
    
    Returns
    -------
    
    Flux-surface average parallel bootstrap current < Jbs . B >
    (scalar)

    """
    
    # If no species specified, get from surface
    if spec == None:
        try:
            spec = surf.species
        except:
            raise ValueError("Must specify species information")
    
    # Calculate viscosity components
    vis = viscosity(surf, spec)

    # Calculate friction coefficients
    fric = friction(spec)
    
    ns = size(spec) # Number of species

    # Load matrices
    vvec = zeros(2*ns)
    tmat = zeros([2*ns, 2*ns])
    rlmat = zeros([2*ns, 2*ns])
    
    bk = 1.6022e-19
    
    ## This loop could be tidied up
    for i in range(0,2*ns):
        di = i % ns
        si = spec[di]
        
        pda  =  bk * ( si.T * si.dndpsi + si.dTdpsi * si.density )
        v1a  = -surf.f * pda / (si.charge * si.density)
        v2a  = -surf.f * bk * si.dTdpsi / si.charge
        
        if i < ns:
            vvec[i] = vis[di,0]*v1a + vis[di,1]*v2a
        else:
            vvec[i] = vis[di,1]*v1a + vis[di,2]*v2a
        
        for j in range(0,2*ns):
            if (i < ns) and (j < ns):
                rlmat[i,j] = -fric.rl11[i,j]
                if i == j:
                    tmat[i,j] = vis[i,0]
            elif (i >= ns) and (j < ns):
                rlmat[i,j] = fric.rl21[i-ns, j]
                if i-ns == j:
                    tmat[i,j] = vis[i-ns,1]
            elif (i < ns) and (j >= ns):
                rlmat[i,j] = fric.rl12[i,j-ns]
                if i == j-ns:
                    tmat[i,j] = vis[i,1]
            else:
                rlmat[i,j] = -fric.rl22[i-ns,j-ns]
                if i == j:
                    tmat[i,j] = vis[i-ns,2]
    
    # Solve using LU decomposition (LAPACK routine) to get parallel velocity
    try:
        upar = linalg.solve(rlmat + tmat, vvec)
    except:
        return 0.

    print("upar = ", upar)

    # Sum contribution from all species
    bstrap = sum([ s.charge*s.density*upar[i] for i,s in enumerate(spec) ])
    bstrap /= sqrt(surf.Bsqav())
    
    return bstrap


def bootstrapSimple(surf, spec=None):
    """ Calculates the Bootstrap current in large aspect ratio, collisionless
    approximation using formula 4.9.2 from Tokamaks 2nd Ed. by Wesson p. 173
    
    Returns
    -------

    < Jb >
    
    note: no factor of B
    
    """
    if spec == None:
        try:
            spec = surf.species
        except:
            raise ValueError("Must specify species information")
    
    # Average Bp * R to convert between d/dpsi and d/dr
    BpR = surf.average(lambda x: surf.Bp(x) * surf.R(x))
    Bp = surf.average(surf.Bp)
    
    for s in spec:
        dTdr = s.dTdpsi * BpR
        
        if s.atomicMass < 0.1:
            # Electron
            dTedr = dTdr
            Te = s.T
            
            n = s.density
            dndr = s.dndpsi * BpR
        else:
            # ion species
            dTidr = dTdr
            Ti = s.T
            
    Jb = 1.602e-19*(sqrt(surf.eps()) * n  / Bp) * ( 2.44*(Te + Ti)*dndr/n + 0.69*dTedr - 0.42*dTidr )
    
    return Jb

def bootstrapWesson(surf, species=None, trapfrac=None):
    """ Bootstrap current calculation
    spanning range of collisionality and aspect ratio
    
    From Tokamaks 2nd Ed by Wesson. Appendix 14.12 
    
    Parameters
    ----------

    
    Returns
    -------
    
    < J . B > flux-surface averaged bootstrap current
    
    Raises
    ------
    
    
    """
    if species == None:
        try:
            species = surf.species
        except:
            raise ValueError("Must specify species information")
        
    # Find which species is the electron and main ion species
    eind = None
    iind = None
    for i, s in enumerate(species):
        if s.atomicMass < 0.1:
            eind = i
        else:
            if iind != None:
                if s.density > species[iind].density:
                    iind = i
            else:
                iind = i
    if eind == None:
        raise ValueError("No electron species given")
    if iind == None:
        raise ValueError("No ion species given")

    # Get electron and ion temperatures in eV
    Te = species[eind].T
    Teprime = species[eind].dTdpsi
    Ti = species[iind].T
    Tiprime = species[iind].dTdpsi
    
    # Calculate pressures in Pascals and pressure gradients
    charge = 1.602e-19
    pe = charge * Te * species[eind].density
    peprime = charge * (Te * species[eind].dndpsi + species[eind].dTdpsi * species[eind].density)
    pi = charge * Ti * species[iind].density
    piprime = charge * (Ti * species[iind].dndpsi + species[iind].dTdpsi * species[iind].density)
    
    if trapfrac == None:
        # Calculate trapped fraction
        trapfrac = trappedFraction(surf)
    x = trapfrac / (1 - trapfrac) # Ratio of trapped to circulating particles
    
    # Get collisionalities
    nustar = collisionality(surf)
    
    nustar_e = nustar[eind] # Electron collisionality
    nustar_i = nustar[iind] # Dominant ion species collisionality

    eps = surf.eps() # Inverse aspect ratio
    
    # Calculate coefficients
    c1 = (4. + 2.6*x) / (  (1. + 1.02*sqrt(nustar_e) + 1.07*nustar_e)*(1. + 1.07*eps**1.5 * nustar_e)  )
    
    c2 = (Ti / Te) * c1
    
    c3 = (7. + 6.5*x) / (  (1. + 0.57*sqrt(nustar_e) + 0.61*nustar_e)*(1. + 0.61*eps**1.5*nustar_e)  )  - 2.5*c1
    
    d = -1.17 / (1. + 0.46*x)

    c4 = (  (d + 0.35*sqrt(nustar_i))/(1.+0.7*sqrt(nustar_i)) + 2.1*eps**3*nustar_i**2  )* c2 / (1. - eps**3*nustar_i**2) / (1. + eps**3*nustar_e**2)

    # Combine into expression for < J . B >
    
    Jb = ((surf.f * x * pe) / (2.4 + 5.4*x + 2.6*x**2)) * ( c1*peprime/pe + c2*piprime/pi + c3*Teprime/Te + c4*Tiprime/Ti )
    
    return Jb
