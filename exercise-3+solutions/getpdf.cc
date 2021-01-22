#include "courselib.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

/*
	Parton evolution:
      1. generate starting distribution in x and add intrisnic transverse mom.
      2. evolve sarting distribution to large values of t
            
	The random number generator is RANLUX by M. Luescher.
      M. Lüscher, A portable high-quality random number generator 
      for lattice field theory simulations, 
      Computer Physics Communications 79 (1994) 100  
	http://luscher.web.cern.ch/luscher/ranlux/index.html      

*/

void gauss2D (double sigma, double &kx, double &ky)
{
    double kT  = sigma * sqrt(-2*log(Rand()));
    double phi = 2*M_PI * Rand();
    kx = kT*cos(phi);
    ky = kT*sin(phi);
}



void get_starting_pdf (double xmin, double q2, double& weightx, double& x, double& kx, double& ky)
{
    double xmax = 0.999;

    // get x value according to g(x)\sim 1/x            
    x = xmin*pow(xmax/xmin,Rand());
    weightx =  x*log(xmax/xmin) ;
    // use gluon like parton density: xg(x) = 3 (1-x)**5
    double pdf = 3.* pow((1.-x),5.)/x ;
    weightx = weightx * pdf;

    // now generate instrinsic kt according to a gauss distribution     
    gauss2D(0.7, kx, ky);
}

void sudakov (double t0, double& t)
{
    // here we calculate  from the sudakov form factor the next t > t0
    const double epsilon=0.1;
    const double as2pi=0.1/2./M_PI;
    double Pint;
    // for Pgg use fact = 2*Ca
    static double Ca = 3. ;
    static double fac=2*Ca; 

    // for Pqq use fact = Cf
    //static double Cf = 4./3.;
    //static double fac=Cf; 
    double r1,t2;
    /* 
       frist do the simple case with fixed alphas and 
       only the 1/(1-z) term of the splitting fct
     */
    r1 = Rand();
    Pint=log((1.-epsilon)/epsilon);  	/* for 1/(1-z) splitting fct */
    t2 = -log(r1)/fac/as2pi/Pint;
    t2 = t0 * exp(t2);
    t = t2;   
    if ( t2 < t0 ) {        
        cout << "FATAL t0 > t:  t0= "<<t0<< "t= "<<t2<< "epsilon= "<<epsilon <<endl;
    }
}

void splitting (double& z, double& weightz)
{
    double epsilon=0.1;
    double as2pi=0.1/2./M_PI;
    double r1,r2;
    double g0,g1,gtot;
    double zmin, zmax;

    // here we calculate  the splitting variable z for 1/z and  1/(1-z)
    // use Pgg=6 (1/z + 1/(1-z))
    g0 = 6.*as2pi * log((1.-epsilon)/epsilon) ;
    g1 = 6.*as2pi * log((1.-epsilon)/epsilon) ;
    gtot = g0 + g1 ;
    zmin = epsilon;
    zmax = 1.-epsilon ;
    r1 = Rand();
    r2 = Rand();
    z = zmin*pow((zmax/zmin),r2);
    if ( r1 > g0/gtot ) { z = 1. - z;}
    weightz = 1.;   

    // cout << "splitting z= "<<z << " weightz = " << weightz <<endl;

}

void evolve_pdf (double xmin, double q2, double x0, double kx0, double ky0, double& weightx, double& x, double& kx, double& ky)
{
    double t0=1.,t1, t=0.,tcut,phi;
    double z=0, weightz;
    double ratio_splitting;
    x = x0;
    kx = kx0;
    ky = ky0;
    weightx = 1. ;
    t1 = t = t0 ;
    tcut = q2;
    while (t1 < tcut ) {
        // here we do now the evolution
        // first we generate t
        t0 = t1;
        sudakov(t0, t1) ;
        // then we generate z
        splitting(z, weightz);
        /*
           since the sudakov involves only the 1/(1-z) part of the 
           splitting fct, we need to weight each branching by the 
           ratio of the integral of the full and approxiamte 
           splitting fct
         */
        ratio_splitting= 2; /* for using Pgg */


        if ( t1 < tcut ) { 
            t = t1;
            x = x*z;
            weightx = weightx *weightz*ratio_splitting;
            phi = 2*M_PI*Rand();  
            /* 
               use here the qt ordering condition:  sqrt(t1) = qt
               and apply this also to the propagator gluons
             */    
            kx +=  sqrt(t1)*cos(phi)*(1.-z);
            ky +=  sqrt(t1)*sin(phi)*(1.-z);
        }
    }
    //            cout << " evolve end  t= "<< t<< " q2 " <<q2 << " x0 = " <<x0 << "x =" << x <<endl;
}


		double static momsum0=0., momsum=0.;


void getpdf (double xmin, double q2,double& weightx, double& x, double& kx, double& ky)
{
    double weightx0, x0, kx0, ky0;
    double weightxf, xf, kxf, kyf;
    weightx0= weightxf = 0;
    // generate starting distribution in x and kt
    get_starting_pdf(xmin, q2, weightx0, x0, kx0, ky0);
    x = x0;
    kx = kx0;
    ky = ky0;
    weightx = weightx0;


    // now do the evolution	
    evolve_pdf(xmin, q2, x0, kx0, ky0, weightxf, xf, kxf, kyf); 
    weightx = weightx0 *weightxf ; 
    momsum0 = momsum0+x0*weightx0 ;
    momsum = momsum+xf*weightx ;

    x = xf;
    kx = kxf;
    ky = kyf;

    // printf (" in getpdf mom sum:  %20.4f     %20.4f \n",momsum0,momsum );
}
void getnorm (double& norm)
{
		norm = momsum0/momsum ;
//            printf (" in getnorm  mom sum:  %20.4f  %20.4f  %20.4f \n",momsum0,momsum,norm );
}
