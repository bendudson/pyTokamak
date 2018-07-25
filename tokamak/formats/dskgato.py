import numpy as np
from utils import file_numbers

def read(f):
    """ Reads a DSKGATO ('t') file
    
    Parameters
    ----------
    
    infile = Input file. Can either be a file-like object,
             or a string. If a string, then treated as a file name
             and opened.
    
    Returns
    -------
    
    A dictionary with the following keys:

    npsi                Number of points in psi
    npol                Number of poloidal points

    psi      [npsi]     psi values
    f(psi)   [npsi]     = R * Bt
    p        [npsi]     Pressure

    R        [npsi, npol]  Major radius [m]
    Z        [npsi, npol]  Height [m]

    Bp       [npsi, npol]  Poloidal magnetic field [T]
    Bt       [npsi, npol]  Toroidal magnetic field [T]
    
    q        [npsi]    Safety factor

    pprime   [npsi]    Pressure gradient dP/dpsi
    ffprime  [npsi]    f * df/dpsi
    
    may contain the following additional keys:
    
    Br       [npsi, npol]   Radial magnetic field [T]
    Bz       [npsi, npol]   Vertical magnetic field [T]
    
    """
    
    if isinstance(f, basestring):
        # If the input is a string, treat as file name
        with open(f) as fh: # Ensure file is closed
            return read(fh) # Call again with file object
    
    # First line contains the date
    date = f.readline()
    if not date:
        raise IOError("Cannot read from input file "+str(filename))
    
    # Second is description
    desc = f.readline()
                
    token = file_numbers(f)
    
    # Third contains number of mesh points
    try:
        npsi = int(token.next())
        ntheta = int(token.next())
        isym = int(token.next())
    except StopIteration:
        raise IOError("Unexpected end of file while reading grid size")
    except ValueError:
        raise IOError("Third line should contain npsi, ntheta and isym")
    
    # Check values
    if (isym < 0) or (isym > 1):
        raise IOError("isym must be either 0 or 1")
    if (npsi < 1) or (ntheta < 1):
        raise IOError("Invalid npsi="+str(npsi)+" or ntheta=" + str(ntheta))
    
    # Read normalisation factors

    try:
        rcnt = float(token.next())
        xma  = float(token.next())
        zma  = float(token.next())
        btor = float(token.next())
        curtot = float(token.next())
        eaxe   = float(token.next())
        dnorm  = float(token.next())
    except:
        raise IOError("Couldn't read normalisation factors")
    
    def read_array(n, name="Unknown"):
        data = np.zeros([n])
        try:
            for i in np.arange(n):
                data[i] = float(token.next())
        except:
            raise IOError("Failed reading array '"+name+"' of size ", n)
        return data

    def read_2d(nx, ny, name="Unknown"):
        data = np.zeros([nx, ny])
        for i in np.arange(nx):
            data[i,:] = read_array(ny, name+"["+str(i)+"]")
        return data

    # Read 1D arrays
    psiflux = read_array(npsi, "psiflux")
    fnorm   = read_array(npsi, "fnorm")
    ffpnorm = read_array(npsi, "ffpnorm")
    ponly   = read_array(npsi, "ponly")
    pponly  = read_array(npsi, "pponly")
    qsf     = read_array(npsi, "qsf")
    d       = read_array(npsi, "d")
    
    dpdz = read_array(ntheta, "dpdz")
    dpdr = read_array(ntheta, "dpdr")
    
    # 2D arrays
    
    xnorm = read_2d(ntheta, npsi, "xnorm")
    znorm = read_2d(ntheta, npsi, "znorm")
    
    # Try to read Br and Bz (may be present)
    try:
        Br = read_2d(ntheta, npsi, "Br")
        Bz = read_2d(ntheta, npsi, "Bz")
    except:
        Br = Bz = None
    
    ny = ntheta

    if isym == 1:
        # Fill in values for up-down symmetric case
        print("Grid is up-down symmetric. Reflecting grid about midplane")
        ny = tsize = 2*(ntheta - 1) + 1
    
        def reflect(data, mapfunc = lambda x:x):
            """ Reflect a variable about midplane
            Optionally supply a mapping function"""
            data2 = np.zeros([tsize, npsi])
            # Copy the original data
            for i in np.arange(ntheta):
                data2[i,:] = data[i,:]
                # Now fill in the remainder
            for i in np.arange(ntheta, tsize):
                t0 = tsize - 1 - i
                data2[i,:] = mapfunc(data[t0,:])
            return data2
    
        xnorm = reflect(xnorm)
        znorm = reflect(znorm, lambda x: 2.*zma - x) # Reflect about zma
        if Br != None:
            Br = reflect(Br, lambda x:-x) # Br reverses
        if Bz != None:
            Bz = reflect(Bz) # Bz remains the same
        theta = tsize

    # Make sure we have Br, Bz and Bpol

    if (Br == None) or (Bz == None):
        # Calculate Bpol from psi then Br and Bz from Bpol
        # Use dpsi = R*Bp dx (for now)
        Bpol = np.zeros([ny, npsi])
        
        def deriv(f):
            n = np.size(f)
            dfdi = np.zeros(n)
            dfdi[1:-1] = (f[2:n] - f[0:-2])/2.  # Central difference in the middle
            dfdi[0] = f[1] - f[0]
            dfdi[-1] = f[-1] - f[-2]
            return dfdi
            
        for i in np.arange(ntheta):
            drdi = deriv(xnorm[i, :])
            dzdi = deriv(znorm[i, :])
            dldi = sqrt(drdi**2 + dzdi**2) # Arc length
            dpsidi = deriv(psiflux)
            
            Bpol[i, :] = dpsidi / (dldi * xnorm[i,:])
    else:
        Bpol = np.sqrt(Br**2 + Bz**2)
    
    # Calculate toroidal field
    Btor = fnorm / xnorm
    
    #########################################
    # Create a dictionary of values to return
    # 
    # Need to transpose 2D arrays to [psi, theta] 
    # to be consistent with elite inputs
    
    var = {"npsi":npsi, "npol":ny,   # Sizes
           
           "psi":psiflux,
           "f(psi)":fnorm,
           "p":ponly,
           
           "R": np.transpose(xnorm),
           "Z": np.transpose(znorm),

           "Bp":np.transpose(Bpol),
           "Bt":np.transpose(Btor),

           "q":qsf,

           "ffprime":ffpnorm,
           "pprime":pponly}

    if Br != None:
        var['Br'] = np.transpose(Br)
    if Bz != None:
        var['Bz'] = np.transpose(Bz)
           
    return var
