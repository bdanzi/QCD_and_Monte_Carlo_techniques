#
#   Testing Monte Carlo integration.
#   We have a test function:  g0=1
#   we calculate the integral \int_0^x\int_0^1dxdy
#
#   We calculate the result and its error.
#
#   
#   Authors: H. Jung, 
#            A. Bermudez Martinez, L.I. Estevez Banos, J. Lidrych, M. Mendizabal Morentin, 
#            S. Taheri Monfared, P. L.S. Connor, Q. Wang, H. Yang, R. Zlebcik 
#
# Import what is needed

from ROOT import gRandom
from math import sqrt

# Loop over random points in the 2D space

npoints = 2000000

# initialise random number generator: 
gRandom.SetSeed(32767)

nhit = 0

for n in range(npoints):
    x,y = gRandom.Uniform(), gRandom.Uniform()
    # insert here the values for x,y and the rejection condition   
    # accept if y < x
    if y < x:
        nhit += 1

# Calculate the integral

Int = float(nhit)/npoints

# the uncertainty is calcualted from a Binominal distribution

sigma2 = (1. - Int)/nhit
error = sqrt(sigma2)*Int

# Show the results

print (" integral for  is: ", Int, "+/-", error)
print (" true integral  is : 0.5 ")
