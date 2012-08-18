# 
# Grad-Shafranov elliptic operator
# 

try:
    from scipy.sparse import lil_matrix
    from numpy import ones
except:
    print("Couldn't load necessary libraries")
    raise

# Map from (r,z) location to row/column index
def indx(r,z):
    return r*nz + z

# Create a sparse matrix given arrays of R and Z values
def matrix(R, Z):
    nr = R.size
    nz = Z.size
    
    # Create a sparse matrix
    N = nr * nz
    A = lil_matrix((N, N))
    
    # Fixed boundaries, so set to identity matrix
    A.setdiag(ones(N))
    
    # Loop over the core of the domain
    for i in range(1,nr-1):
        for j in range(1,nz-1):
            rp = 0.5*(R[i+1] + R[i]) # R_{i+1/2}
            rm = 0.5*(R[i] + R[i-1]) # R_{i-1/2}
            
            # Mesh spacing
            dr = 0.5*(R[i+1] - R[i-1])
            dz = 0.5*(Z[j+1] - Z[j-1])
            
            ind = indx(i,j)
            
            A[ind, ind]         = -(R[i]/dr**2)*(1./rp + 1./rm) - 2./dz**2
            A[ind, indx(i+1,j)] =  (R[i]/dr**2)/rp
            A[ind, indx(i-1,j)] =  (R[i]/dr**2)/rm
            A[ind, indx(i,j+1)] = 1./dz**2
            A[ind, indx(i,j-1)] = 1./dz**2
    # Done. Convert to CSR format
    return A.tocsr()

if __name__ == "__main__":
    # Test case
    
    from numpy import linspace, zeros, arange, transpose
    from scipy.sparse.linalg import spsolve
    
    nr = 50
    nz = 50
    
    # Define R and Z arrays
    r = linspace(0.1, 3, nr)
    z = linspace(-1,1, nz)
    
    # Create the elliptic operator sparse matrix
    A = matrix(r, z)
    
    # Create a vector for b, containing all zeros
    N = nr*nz
    b = zeros(N)
    
    # Set the values on  boundaries to 1
    #b[indx(0, arange(nz))] = 1.
    #b[indx(nr-1, arange(nz))] = 1.
    
    

    # Direct solve sparse matrix
    x = spsolve(A, b)
    
    # View X as a 2D array
    xv = x.view()
    xv.shape = (nr, nz)
    
    # Make a filled contour plot
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    im = plt.imshow(transpose(xv), interpolation='bilinear', origin='lower',
                    cmap=cm.hot, extent=(r[0],r[-1],z[0],z[-1]))
    plt.colorbar()
    plt.show()
    
