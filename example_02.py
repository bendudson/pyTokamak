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
equil.setDensity(1e18) # Set a constant density

#########################################
# Get a single flux-surface (interpolating)
# and calculate the bootstrap current

from tokamak import neoclass
from numpy import linspace, zeros

npoints = 20

psi = linspace(0.5, 0.8, npoints, endpoint=True)
bs_s = zeros(npoints)
bs_w = zeros(npoints)
bs_hs = zeros(npoints)

for i, p in enumerate(psi):
    f = equil.getFluxSurface(p)
    bs_s[i] = neoclass.bootstrapSimple(f) * f.average(f.B)
    bs_w[i] = neoclass.bootstrapWesson(f)
    bs_hs[i] = neoclass.bootstrapHS(f)
    print p, bs_s[i], bs_w[i], bs_hs[i]


