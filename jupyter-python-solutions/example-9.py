# Calculation of the Higgs cross section
#   Monte Carlo integration.
#   We calculate the x-section for p p -> Higgs  + X
#   and plot the distribution in pt and rapidity
#
#   We use our own evolved pdf.
#
# ![gg2H](https://upload.wikimedia.org/wikipedia/commons/0/0e/Higgs-gluon-fusion.svg)  
#   
#   Authors: H. Jung, 
#            A. Bermudez Martinez, L.I. Estevez Banos, J. Lidrych, M. Mendizabal Morentin, 
#            S. Taheri Monfared, P. L.S. Connor, Q. Wang, H. Yang, R. Zlebcik 
#   
# Import what is necessary

from math import pi, sin, cos, log, log10, sqrt, exp
from ROOT import gRandom, gStyle, TLorentzVector, TCanvas, TH1D


# The cross section of the Higgs

def sigma(m2):
    # calculate partonic x-section
    #  
    # calculate sigma0 (qq -> gamma^*)
    # double aem=1./137.;
    # result = 4.*pi*pi*aem/3.;
    #
    # calculate sigma0(gg->H) 
    # sigma0 = as^2/pi /576 * GF * sqrt(2)*\Delta(1.-tau) 
    aS = 0.1      # alphaS
    GF = 1.166e-5 # Fermi constant
    result = aS**2 * sqrt(2.) * GF /pi / 576. 
    return result 

# here we use the evolution code from example-7

# Define function which returns 2D point according to the Gaussian distribution

def gauss2D(sigma):
    kT  = sigma * sqrt(-2*log(gRandom.Uniform()));
    phi = 2*pi * gRandom.Uniform()
    kx, ky = kT*cos(phi), kT*sin(phi)
    return (kx, ky)

# Get PDF at the scale $q^2$ (currently no dependence on $q^2$)

def get_starting_pdf(q2):
    xmin = 0.0001
    xmax = 0.999
    # get x value according to g(x)\sim 1/x            
    x = xmin * (xmax/xmin)**gRandom.Uniform()
    weightx = x*log(xmax/xmin) 
    # use: xg(x) = 3 (1-x)**5
    pdf = 3.*pow((1-x),5)/x

    # now generate instrinsic kt according to a gauss distribution  
    kx, ky = gauss2D(0.7)
    pVec = TLorentzVector()
    pVec.SetXYZM(kx, ky, x*Eb, 0.)
    return pVec, weightx * pdf

# Calculate Sudakov form factor

def sudakov(t0):
    #   here we calculate  from the sudakov form factor the next t > t0
    epsilon = 0.1
    as2pi = 0.1/(2.*pi)
    Ca = 3.
    # for Pgg use fact = 2*Ca
    fac = 2.*Ca 
    # use fixed alphas and only the 1/(1-z) term of the splitting fct

    r1 = gRandom.Uniform()
    Pint=log((1.-epsilon)/epsilon) # for 1/(1-z) splitting fct 
    t2 = -log(r1)/fac/as2pi/Pint
    t2 = t0 * exp(t2)
    assert(t2 >= t0)
    return t2

# The splitting function which is needed to get the z on the branching point

def splitting():
    epsilon = 0.1

    as2pi = 0.1/2./pi

    #	here we calculate  the splitting variable z for 1/z and  1/(1-z)
    #	use Pgg=6(1/z + 1/(1-z))  // for large z we use z -> 1-z

    g0 = 6.*as2pi * log((1.-epsilon)/epsilon)
    g1 = 6.*as2pi * log((1.-epsilon)/epsilon)
    gtot = g0 + g1 

    zmin = epsilon
    zmax = 1.-epsilon

    r1 = gRandom.Uniform()
    r2 = gRandom.Uniform()

    z = zmin * (zmax/zmin)**r2
    if r1 > g0/gtot:
        z = 1. - z
    weightz = 1.
    return z

# Evolve the PDF between scales $q^2_0$ and $q^2$, the kinematics of the parton in the beginning is described by four-vector p0

def evolve_pdf(q20, q2, p0):
    x = p0.Pz()/Eb
    kx = p0.Px()
    ky = p0.Py()
    weightx = 1.
    t1 = q20
    tcut = q2
    while t1 < tcut:
        # here we do now the evolution
        # first we generate t
        t0 = t1
        t1 = sudakov(t0) 
        
        # now we generate z
        z = splitting()
        #   since the sudakov involves only the 1/(1-z) part 
        #   of the splitting fct, we need to weight each branching
        #   by the ratio of the integral of the full and 
        #   approximate splitting fct

        ratio_splitting = 2 # for using Pgg

        if  t1 < tcut:
            x = x*z
            weightx = weightx *ratio_splitting
            # 
            # use here the angular ordering condition: sqrt(t1) = qt/(1-z) 
            # and apply this also to the propagator gluons
            #
            phi = 2*pi*gRandom.Uniform()
            kx +=  sqrt(t1)*cos(phi)*(1.-z)
            ky +=  sqrt(t1)*sin(phi)*(1.-z)                     
            #   kx += sqrt(t1)*cos(phi)
            #   ky += sqrt(t1)*sin(phi)                     
    k = TLorentzVector()
    k.SetXYZM(kx, ky, x*Eb, 0.)
    return k, weightx

# get parton distribution

def getpdf(q2):
    # generate starting distribution in x and kt
    q20 = 1
    k0, weightx0 = get_starting_pdf(q20)
    x0 = k0.Pz()/Eb
    
    # now do the evolution	
    k, weightxf = evolve_pdf(q20, q2, k0)  
    xf = k.Pz()/Eb
     
    weightx = weightx0 * weightxf
    momsum0 = x0*weightx0
    momsum = xf*weightx
    
    return k, weightx 
        
# end evolution code from example-7
        

# Book the histograms

histo1  =  TH1D("x1","x1",100, -4, 0.)
histo2  =  TH1D("x2","x2",100, -4, 0.)
histo3  =  TH1D("kt1 ","kt1 ",100, 0, 100.)
histo4  =  TH1D("kt2 ","kt2 ",100, 0, 100.)
histo5  =  TH1D("pt ","pt ",50, 0, 200.)
histo6  =  TH1D("eta ","eta",50, -8, 8.)
histo7  =  TH1D("Mass ","mass",50, 60., 160.)

s=4*3500.*3500 # center of mass energy
Eb = sqrt(s)/2 # Beam energy

print( " sqrt(s) = " , sqrt(s) )

# Scale of the process and the mass window

q2 = 10000 # this gives reasonable results... remember we use simplifications
mass_min = 0.
mass_max = 200.

m_higgs = 125.
Gamma = 0.4

# Loop over MC events

sum0 = sum00 = 0
nacc = 0
#npoints = 100000
npoints = 10000000



for n1 in range(npoints):
    # generate p4 of incoming parton 1
    pA, weightx1 = getpdf(q2)
    # generate p4 of incoming parton 2
    pB, weightx2 = getpdf(q2)
    pB.SetPz(-pB.Pz())

    x1 = pA.Pz() / Eb
    x2 =-pB.Pz() / Eb
    # plot dxg(x)/dlogx *Jacobian, Jacobian dlogx/dx = 1/x
    histo1.Fill(log10(x1),weightx1/log(10))
    histo2.Fill(log10(x2),weightx2/log(10))
    histo3.Fill(pA.Pt(), weightx1)
    histo4.Fill(pB.Pt(), weightx2)

    # total p4
    pH = pA + pB
    mass = pH.M()

    # calculate rapidity of Higgs
    # rapidity = 0.5 * log(x1/x2) 
    if mass < mass_min or mass > mass_max:
        continue
    # note, pdfs are already included in weightx  
    
    Hx1x2 = weightx1*weightx2* sigma(mass**2)
    # multiply with a Breit Wigner resonance
    Hx1x2 = Hx1x2 * Gamma/(pow((mass-m_higgs),2) + pow(Gamma,2)/4)/2./pi
             
    # Change units to nb
    gev2nb = 0.3893857E+6
    ff = Hx1x2 * gev2nb

    nacc += 1
    sum0  +=  ff
    sum00 +=  ff**2
    # weighting with 1/x0:
    # plot dxg(x)/dlogx *Jacobian, Jacobian dlogx/dx = 1/x

    histo5.Fill(pH.Pt(), ff)
    histo6.Fill(pH.Rapidity(), ff)
    histo7.Fill(mass,  ff)

# Normalization and evaluation of the error

sum0  /= npoints
sum00 /= npoints
sigma2 = sum00 - sum0*sum0 
error = sqrt(sigma2/npoints) 
print (" nr of events accepted: ", nacc )
print (" integral for Higgs xsection is [pb]: " , sum0 * 1000., " +/- " , error*1000.)

# Plotting of x1, x2 and kT1 and kT2

gStyle.SetOptStat(0) # get rid of statistics box
c = TCanvas()
# divide the canvas in 1 parts in x and 1 in y
c.Divide(2,2)
c.cd(1)
histo1.Scale(1./npoints, "width")
histo1.Draw()
c.cd(2)
histo2.Scale(1./npoints, "width")
histo2.Draw()
c.cd(3)
histo3.Draw()
c.cd(4)
histo4.Draw()
c.Draw()

# Draw the pT, eta and mass of the Higgs boson

c5 = TCanvas()
histo5.Draw()
c5.Draw()
c6 = TCanvas()
histo6.Draw()
c6.Draw()
c7 = TCanvas()
histo7.Draw()
c7.Draw()


