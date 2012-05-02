#
# Example of reading an input file, extracting flux surface data
# 


#########################################
# Read an equilibrium from a DSKGATO 't' formatted file

from tokamak.formats import dskgato

print "Reading DSKGATO input file"
data = dskgato.read("test.dskgato")

#########################################
# Create an equilibrium object from the data

from tokamak.equilibrium import Equilibrium

equil = Equilibrium(data)
equil.setDensity(1e20) # Set a constant density

#########################################
# Get a single flux-surface (interpolating)
# and calculate the bootstrap current

from tokamak import neoclass
from numpy import linspace, zeros

npoints = 20

psi = linspace(0.5, 0.8, npoints, endpoint=True)
bs = zeros(npoints)
bs2 = zeros(npoints)
for i, p in enumerate(psi):
    f = equil.getFluxSurface(p)
    bs2[i] = neoclass.bootstrapHS(f)
    bs[i] = neoclass.bootstrapSimple(f)
    print p, bs[i], bs2[i]

