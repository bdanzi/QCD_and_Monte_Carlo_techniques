#include "courselib.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

/*
   Testing Monte Carlo integration.
   We have a test function:  g0=3x**2
   We weight the histo with the function value.
   It is integrated over the region from 0 to 1.
      
   We calculate the result and its error.

   The random number generator is RANLUX by M. Luescher.
   M. Lüscher, A portable high-quality random number generator 
   for lattice field theory simulations, 
   Computer Physics Communications 79 (1994) 100
   
   Authors: H. Jung, R. Zlebcik 
     
*/


double g0(double z)
{	
    return 3*z*z; 
}

int main(int argc,char **argv)
{
    const int npoints = 10000000;

    // initialise random number generator: rlxd_init( luxory level, seed )
    rlxd_init(2,32767);

    double xg0 = 0, xg00 = 0;
    // using uniformly random numbers to obtain the integral (b-1)=1-0
    for (int n1 = 0; n1 < npoints; ++n1) {
        double x0 = Rand();;
        double f  = g0(x0); 
        double ff = f;
        xg0  +=  ff;
        xg00 +=  ff*ff; //per calcolare l'errore
    }

    xg0  /= npoints;
    xg00 /= npoints;
    double sigma2 = xg00 - xg0*xg0 ;
    double error = sqrt(sigma2/npoints) ; //deviazione standard
    cout<<" integral for 3x**2 is: "<<xg0<<"+/-"<< error<<endl;
    cout<<" true integral for 3x**2 is : 1.0 "<< endl;

    return EXIT_SUCCESS;
}

