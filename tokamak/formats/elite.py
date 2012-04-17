
from utils import file_tokens

def read(filename):
    """ Reads an ELITE .eqin format file
    """
    f = open(infile)
    
    desc = f.readline()
    if not desc:
        raise IOError("Cannot read from input file")
    
    # Define a generator to get the next token from the file
    token = file_tokens(f)
    
    try:
        npsi = int(token.next())
        npol = int(token.next())
    except StopIteration:
        raise IOError("Unexpected end of file while reading grid size")
    except ValueError:
        raise IOError("Second line should contain Npsi and Npol")
    
    # Result will be a dictionary of variables
    var=dict()
    
    while True:
        try:
            varname = token.next()
        except:
            break
        
        try: 
            if (varname == 'R') or (varname == 'Z') or (varname[0] == 'B'):
                # A 2D variable
            
                data = np.zeros([npsi, npol])
                for j in np.arange(npol):
                    for i in np.arange(npsi):
                        data[i,j] = float(token.next())
            else:
                # Assume it's a 1D (psi) variable
                data = np.zeros([npsi])
                for i in np.arange(npsi):
                    data[i] = float(token.next())
        
        except StopIteration:
            raise IOError("Unexpected end of file while reading " + varname)
        except ValueError:
            raise IOError("Expecting float while reading " + varname)
        except:
            raise IOError("Out of cheese. Variable " + varname)
        
        # Add this variable to the dictionary
        var[varname] = data
    
    # Finished reading data
    f.close()
    
    # Calculate toroidal field

    try:
        Bt = np.zeros([npsi, npol])
        f = var['f(psi)']
        R = var['R']
        for i in np.arange(npsi):
            Bt[i,:] = f[i] * R[i,:]
        var['Bt'] = Bt
    except KeyError:
        raise IOError("Need f(psi) and R to calculate Bt")

    # Temperatures: If only have one, set equal

    if var.has_key('Te') and not var.has_key('Ti'):
        var['Ti'] = var['Te']

    if var.has_key('Ti') and not var.has_key('Te'):
        var['Te'] = var['Ti']

    # Calculate pressure

    if not var.has_key('p'):
        if var.has_key('ne') and var.has_key('Te'):
            # Calculate using Ni, Te and Ti
            var['p'] = var['ne'] * (var['Te'] + var['Ti']) * 1.602e-19 * (4.0*3.14159*1e-7)
            # Could check against pprime
        else:
            # No plasma quantities to use, so integrate pprime
            raise IOError("Can't calculate pressure")
    
    return var
