#
#   Testing Monte Carlo integration.
#   We have a test function:  g0=(1-x)**pow1/x
#   We weight the histo with the function value.
#   It is integrated over the region from xmin to 1.
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

# Import what will be needed

from ROOT import gRandom, TMath
from math import sqrt, log

# Define function which we want to integrate

def g0(z):
    return (1 - z)**5 / z


npoints = 1000000
# The lower limit of the integration, the upper is 1
xmin = 1e-4

# initialise random number generator: 
gRandom.SetSeed(32767)

# Generate npoints randomly with importance sampling

xg0 = xg00 = 0
for n1 in range(npoints):
    # here do the calculation with importance sampling
    x0 = xmin**gRandom.Uniform()
    weight = x0*log(1/xmin)
    # here do the calculation using linear sampling
    # x0 = xmin+(1-xmin)*gRandom.Uniform()
    # weight = 1-xmin
    f  = g0(x0) 
    ff = weight*f
    xg0  +=  ff
    xg00 +=  ff**2

# Calculate the MC integral

xg0  /= npoints
xg00 /= npoints
sigma2 = xg00 - xg0*xg0
error  = sqrt(sigma2/npoints)
print (" integral for g(x) = (1-x)**5/x is: ",xg0,"+/-", error)

# Get the exact value using incomplete beta function [https://en.wikipedia.org/wiki/Beta_function]

x = 1-xmin
a = 5
b = -0.9999999999
print ("Exact value ", TMath.BetaIncomplete(x, a+1, b+1) * TMath.Beta(a+1, b+1))
