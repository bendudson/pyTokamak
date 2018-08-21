"""
Read and write 'G' formatted equilibria. This is an R-Z free boundary
format.

Format of G-EQDSK file is specified here:
  https://fusion.gat.com/THEORY/efit/g_eqdsk.html
"""

import numpy as np
from .utils import file_numbers, writef

try:
  basestring
except NameError:
  basestring = str

def read(f):
    """ Reads a G-EQDSK file

    Parameters
    ----------

    f = Input file. Can either be a file-like object,
        or a string. If a string, then treated as a file name
        and opened.

    Returns
    -------

    """

    if isinstance(f, basestring):
        # If the input is a string, treat as file name
        with open(f) as fh: # Ensure file is closed
            return read(fh) # Call again with file object

    # Read the first line, which should contain the mesh sizes
    desc = f.readline()
    if not desc:
        raise IOError("Cannot read from input file")

    s = desc.split() # Split by whitespace
    if len(s) < 3:
        raise IOError("First line must contain at least 3 numbers")

    idum = int(s[-3])
    nxefit = int(s[-2])
    nyefit = int(s[-1])

    # Use a generator to read numbers
    token = file_numbers(f)

    try:
        xdim   = float(token.next())
        zdim   = float(token.next())
        rcentr = float(token.next())
        rgrid1 = float(token.next())
        zmid   = float(token.next())

        rmagx  = float(token.next())
        zmagx  = float(token.next())
        simagx = float(token.next())
        sibdry = float(token.next())
        bcentr = float(token.next())

        cpasma = float(token.next())
        simagx = float(token.next())
        xdum   = float(token.next())
        rmagx  = float(token.next())
        xdum   = float(token.next())

        zmagx  = float(token.next())
        xdum   = float(token.next())
        sibdry = float(token.next())
        xdum   = float(token.next())
        xdum   = float(token.next())
    except:
        xdim = float(next(token))
        zdim = float(next(token))
        rcentr = float(next(token))
        rgrid1 = float(next(token))
        zmid = float(next(token))

        rmagx = float(next(token))
        zmagx = float(next(token))
        simagx = float(next(token))
        sibdry = float(next(token))
        bcentr = float(next(token))

        cpasma = float(next(token))
        simagx = float(next(token))
        xdum = float(next(token))
        rmagx = float(next(token))
        xdum = float(next(token))

        zmagx = float(next(token))
        xdum = float(next(token))
        sibdry = float(next(token))
        xdum = float(next(token))
        xdum = float(next(token))

    # Read arrays
    def read_array(n, name="Unknown"):
        data = np.zeros([n])
        try:
            for i in np.arange(n):
                try:
                    data[i] = float(token.next())
                except:
                    data[i] = float(next(token))
        except:
            raise IOError("Failed reading array '"+name+"' of size ", n)
        return data

    def read_2d(nx, ny, name="Unknown"):
        data = np.zeros([nx, ny])
        for i in np.arange(nx):
            data[i,:] = read_array(ny, name+"["+str(i)+"]")
        return data

    fpol   = read_array(nxefit, "fpol")
    pres   = read_array(nxefit, "pres")
    workk1 = read_array(nxefit, "workk1")
    workk2 = read_array(nxefit, "workk2")
    psi    = read_2d(nxefit, nyefit, "psi")
    qpsi   = read_array(nxefit, "qpsi")

    # Read boundary and limiters, if present
    try:
        nbdry = int(token.next())
        nlim  = int(token.next())
    except:
        nbdry = int(next(token))
        nlim = int(next(token))

    if nbdry > 0:
        rbdry = np.zeros([nbdry])
        zbdry = np.zeros([nbdry])
        for i in range(nbdry):
            try:
                rbdry[i] = float(token.next())
                zbdry[i] = float(token.next())
            except:
                rbdry[i] = float(next(token))
                zbdry[i] = float(next(token))
    else:
        rbdry = [0]
        zbdry = [0]

    if nlim > 0:
        xlim = np.zeros([nlim])
        ylim = np.zeros([nlim])
        for i in range(nlim):
            try:
                xlim[i] = float(token.next())
                ylim[i] = float(token.next())
            except:
                xlim[i] = float(next(token))
                ylim[i] = float(next(token))
    else:
        xlim = [0]
        ylim = [0]

    # Construct R-Z mesh
    r = np.zeros([nxefit, nyefit])
    z = r.copy()
    for i in range(nxefit):
        r[i,:] = rgrid1 + xdim*i/float(nxefit-1)
    for j in range(nyefit):
        z[:,j] = (zmid-0.5*zdim) + zdim*j/float(nyefit-1)

    # Create dictionary of values to return
    result = {'nx': nxefit, 'ny':nyefit,        # Number of horizontal and vertical points
              'r':r, 'z':z,                     # Location of the grid-poinst
              'rdim':xdim, 'zdim':zdim,         # Size of the domain in meters
              'rcentr':rcentr, 'bcentr':bcentr, # Reference vacuum toroidal field (m, T)
              'rgrid1':rgrid1,                  # R of left side of domain
              'zmid':zmid,                      # Z at the middle of the domain
              'rmagx':rmagx, 'zmagx':zmagx,     # Location of magnetic axis
              'simagx':simagx, # Poloidal flux at the axis (Weber / rad)
              'sibdry':sibdry, # Poloidal flux at plasma boundary (Weber / rad)
              'cpasma':cpasma,
              'psi':psi,    # Poloidal flux in Weber/rad on grid points
              'fpol':fpol,  # Poloidal current function on uniform flux grid
              'pressure':pres,  # Plasma pressure in nt/m^2 on uniform flux grid
              'qpsi':qpsi,  # q values on uniform flux grid
              'nbdry':nbdry, 'rbdry':rbdry, 'zbdry':zbdry, # Plasma boundary
              'nlim':nlim, 'xlim':xlim, 'ylim':ylim} # Wall boundary

    return result

def write(f, data):
    """ Write a G-EQDSK file

    """

    if isinstance(f, basestring):
        # If the input is a string, treat as file name
        with open(f, "w") as fh: # Ensure file is closed
            return write(fh, data) # Call again with file object

    nx = int(data['nx'])
    ny = int(data['ny'])

    # Write description
    f.write("geqdsk 0" + str(nx) + str(ny))

    writef(f, data['rdim'])
    writef(f, data['zdim'])
    writef(f, data['rcentr'])
    writef(f, data['rgrid1'])
    writef(f, data['zmid'])

    writef(f, data['rmagx'])
    writef(f, data['zmagx'])
    writef(f, data['simagx'])
    writef(f, data['sibdry'])
    writef(f, data['bcentr'])

    writef(f, data['cpasma'])
    writef(f, data['simagx'])
    writef(f, 0.0)
    writef(f, data['rmagx'])
    writef(f, 0.0)

    writef(f, data['zmagx'])
    writef(f, 0.0)
    writef(f, data['sibdry'])
    writef(f, 0.0)
    writef(f, 0.0)

    # Write the arrays
    writef(f, data['fpol'])
    writef(f, data['pressure'])
    writef(f, np.zeros(nx))  # Workk1
    writef(f, np.zeros(nx))  # Workk2
    writef(f, data['psi'])
    writef(f, data['qpsi'])
