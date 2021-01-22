#!/usr/bin/env  python3
#  
#   Constructing a random number generator
#
#   Authors: H. Jung, 
#            A. Bermudez Martinez, L.I. Estevez Banos, J. Lidrych, M. Mendizabal Morentin, 
#            S. Taheri Monfared, P. L.S. Connor, Q. Wang, H. Yang, R. Zlebcik 
#   

# First import from ROOT libraries which are needed

from  ROOT import TH2D, TCanvas, gRandom, TRandom3

# Define the 2D histograms to display the correlation between random numbers

hConguent = TH2D("Congruent", "congruental random numbers ",100, 0, 1, 100, 0, 1)
hRANLUX = TH2D("RANLUX", "RANLUX",100, 0, 1, 100, 0, 1)

# Define our first congruent generator

def randCon():
    # construct random numbers according the formula in exercise 1
    im = 139968
    ia = 205
    ic = 29573
    last = 4711
    randCon.r = (ia*randCon.r+ic) % im #get next random number
    return float(randCon.r) / im #return normalized one

ic = 29573

randCon.r = ic #Starting value (seed)

# Number of points to be generated

npoints = 100000

# Generate npoint pairs of random numbers with congruent generator and fill it to 2D histogram

#myRandom = TRandom3
#myRandom = gRandom.TRandom
# Initialize random number generator.
gRandom.SetSeed()
ranlux = TRandom3

for n in range(npoints):
    # fill the 2 random numbers into a scatter plot to compare the correlation
    xC = randCon()
    yC = randCon()
    hConguent.Fill(xC, yC)

    # next use the ranlux generator and compare also the 2 random numbers

    xRL = gRandom.Uniform()
    yRL = gRandom.Uniform()
    #xRL = ranlux
    #yRL = ranlux
    hRANLUX.Fill(xRL, yRL)

# Plotting the results

c = TCanvas()
hConguent.Draw()
c.Draw()

d = TCanvas()
hRANLUX.Draw()
d.Draw()
