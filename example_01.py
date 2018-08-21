#
# Some tests of the equilibrium and neoclass routines
# 

from tokamak.equilibrium import FluxSurface
from tokamak import neoclass

from numpy import linspace, pi, zeros, sin, cos, sqrt

Rmaj = 5.   # Major radius
amin = 1.   # Minor radius
Bt0 = 10.
Bp0 = 1.

nth = 20 
theta = linspace(0,2*pi, nth, endpoint=False)

r = Rmaj + amin * cos(theta)
z = sin(theta)
Bp = zeros(nth) + Bp0

# Create a flux surface object
f = FluxSurface(Rmaj*Bt0, r, z, Bp)

q = Bt0*amin / (Bp0 * Rmaj) # Large aspect-ratio approximation
eps = amin / Rmaj

print("Inverse aspect ratio eps = ", eps)

print("")
print("Large aspect ratio q = ", q)
print("q = ", f.integral(lambda x: f.Bt(x) / (f.Bp(x) * f.R(x))) / (2 * pi))

print("")
print("Large aspect ratio qR = ", q * Rmaj)
print("Connection length = ", neoclass.connectionLength(f))

print("")
print("Large aspect-ratio trapped = ", 1.46 * sqrt(eps))
print("Trapped fraction  = ", neoclass.trappedFraction(f))
