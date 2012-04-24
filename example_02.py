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

#########################################
# Get a single flux-surface (interpolating)
# and calculate the connection length for it

from tokamak import neoclass

f = equil.getFluxSurface(0.8)   # Normalised psi = 0.8

print neoclass.connectionLength(f)


#########################################
# Iterate over flux surfaces

for f in equil.surfaces():
    print neoclass.connectionLength(f)
