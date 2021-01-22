#   Testing Monte Carlo integration.
#   We have a test function:  g0=3x**2
#   We weight the histo with the function value.
#   It is integrated over the region from 0 to 1.
#      
#   We calculate the result and its error.
#
#   The random number generator is RANLUX by M. Luescher.
#   M. Luescher, A portable high-quality random number generator 
#   for lattice field theory simulations, 
#   Computer Physics Communications 79 (1994) 100
#   
#   Authors: H. Jung, 
#            A. Bermudez Martinez, L.I. Estevez Banos, J. Lidrych, M. Mendizabal Morentin, 
#            S. Taheri Monfared, P. L.S. Connor, Q. Wang, H. Yang, R. Zlebcik 
#
# Define function which we want to integrate

from ROOT import gRandom
from math import sqrt

def g0(z):
    return 3*z*z; 

# initialise random number generator: 
gRandom.SetSeed(32767)

# Calculate the sum and sum2 of the function values at the random points
xg0 = xg00 = 0
npoints = 1000000
for n in range(npoints):
    x0 = gRandom.Uniform()
    f  = g0(x0) 
    xg0 += f
    xg00+= f**2

# Calculate average and average to the squared values

avg  = xg0  / npoints
avg2 = xg00 / npoints

# Get the average deviation squared

sigma2 = avg2 - avg*avg

# The actual error behaves as 1/sqrt(npoints)

error = sqrt(sigma2/npoints) ;

# Finally the result

print (" integral for 3x**2 is: ", avg, "+/-", error)
print (" true integral for 3x**2 is : 1.0 ")
