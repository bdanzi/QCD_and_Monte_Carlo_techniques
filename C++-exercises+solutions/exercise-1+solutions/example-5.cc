x0#include "courselib.h"
#include <cmath>
#include <TMath.h>
#include <iostream>
using namespace std;

/*
   Testing Monte Carlo integration.
   We have a test function:  g0=(1-x)**pow1/x
   We weight the histo with the function value.
   It is integrated over the region from xmin to 1.

   We calculate the result and its error.

   The random number generator is RANLUX by M. Luescher.
   M. Lüscher, A portable high-quality random number generator 
   for lattice field theory simulations, 
   Computer Physics Communications 79 (1994) 100  

   Authors: H. Jung, R. Zlebcik 
   
*/


double g0(double z)
{	
    double result=pow(1 - z, 5) / z;
    return result; 
}

int main(int argc,char **argv)
{
    const int npoints = 1000000;
    const double xmin = 1e-4;

    // initialise random number generator: rlxd_init( luxory level, seed )
    rlxd_init(2, 32767);

    double xg0 = 0, xg00 = 0;

    for (int n1 = 0; n1 < npoints; ++n1) {

        // here do the calcualtion with importance sampling
        double x0 = pow(xmin, Rand());
        //genero random number che seguono con pdf 1/x/log(xmax/xmin)
        double weight = x0*log(1/xmin) ;
        // here do the calcualtion using linear sampling
        // x0 = xmin+(1-xmin)*Rand();
        // weight = 1-xmin ;
        double f  = g0(x0); 
        double ff = weight*f;
        xg0  +=  ff;
        xg00 +=  ff*ff; 
    }

    xg0  /= npoints;
    xg00 /= npoints;
    double sigma2 = xg00 - xg0*xg0 ;
    double error  = sqrt(sigma2/npoints) ;
    cout<<" integral for g(x) = (1-x)**5/x is: "<<xg0<<"+/-"<< error<<endl;

    double x = 1-xmin;
    double b = -0.9999999999, a = 5;
    cout << "Exact value "<< TMath::BetaIncomplete(x, a+1, b+1) * TMath::Beta(a+1, b+1) << endl;

    return EXIT_SUCCESS;
}

