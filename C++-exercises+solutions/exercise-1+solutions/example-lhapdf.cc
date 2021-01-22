#include "courselib.h"
#include <cmath>
#include "TH1.h"
#include "TFile.h" 
#include <iostream>
#include <string>
#include "LHAPDF/LHAPDF.h"
#include <iostream>
#include <cstdlib>

using namespace std;

/*
   We calculate the momentum sum rule of a parton density function by performing explicitly the integral.
   We use the LHAPDF library to access the different PDF parameterisations.
   We check the momentum sum rule for a standard set, and for the LO* parameterisations.

   We use LHAPDF 6.2.0

   The random number generator is RANLUX by M. Luescher.
   M. Lüscher, A portable high-quality random number generator 
   for lattice field theory simulations, 
   Computer Physics Communications 79 (1994) 100  
   http://luscher.web.cern.ch/luscher/ranlux/index.html   
   
   Authors: H. Jung, R. Zlebcik 
      
*/

int main()
{
    const int SUBSET = 0;
// 	const string NAME = "MRST2004nlo";
      const string NAME = "PB-TMDNLO-HERAI+II-2018-set1";
//      const string NAME = "PB-TMDNLO-HERAI+II-2018-set2";
// 	const string NAME = "MRST2007lomod";		// LO * pdf //
//    const string NAME = "MRSTMCal";			// LO** pdf //
    LHAPDF::initPDFSet(NAME, LHAPDF::LHGRID, SUBSET);

    const int npoints = 100000;


//  initialise random number generator
//  rlxd_init( luxory level, seed )
    rlxd_init(2,32767);

    const double Q = 10.0, mz = 91.2;
    cout << "alphas(mz) = " << LHAPDF::alphasPDF(mz) << endl;
    cout << "qcdlam4    = " << LHAPDF::getLam4(SUBSET) << endl;
    cout << "qcdlam5    = " << LHAPDF::getLam5(SUBSET) << endl;
    cout << "orderPDF   = " << LHAPDF::getOrderPDF() << endl;
    cout << "xmin       = " << LHAPDF::getXmin(SUBSET) << endl;
    cout << "xmax       = " << LHAPDF::getXmax(SUBSET) << endl;
    cout << "q2min      = " << LHAPDF::getQ2min(SUBSET) << endl;
    cout << "q2max      = " << LHAPDF::getQ2max(SUBSET) << endl;
    cout << "orderalfas = " << LHAPDF::getOrderAlphaS() << endl;
    cout << "num flav   = " << LHAPDF::getNf() << endl;
    cout << "name       = " << NAME << endl;
    cout << "number     = " << LHAPDF::numberPDF() << endl;
    cout << endl;

    //const int NUMBER = numberPDF();

    LHAPDF::initPDF(0);

    const double xmin = LHAPDF::getXmin(0);
    const double xmax = LHAPDF::getXmax(0);

    //      xmax = 0.05 ;
    //      xmin = 0.0001 ;

    double sum0=0, sum00=0;

    for (int n1 = 0; n1 < npoints; ++n1) {
        // for simple integration
        // x = xmin + (xmax-xmin)*Rand();
        
        
        // for importance sampling

        double x = xmin * pow((xmax/xmin), Rand());
        double f=0;
        int flavor;
        // the pdf from LHAPDFLIB is called via: xfx(x,Q,flavor) 
        //		with x=fractional momentum
        //		     Q = sqrt(q2), the sqrt of the scale
        //		     flavor = -6,.. ,6 the flavor code of the parton, 
        //		     flavor: 0=gluon, 1=down, 2=up, 3=strange, 4=charm, 5=bottom 6=top           

        //  sum over all flavors for mom sum rule
        for ( flavor = -6; flavor <= 6; flavor++)
            f += LHAPDF::xfx(x,Q,flavor);

        //  take only flavor 1 (2) for flavor sum rule
        flavor = 1;
        //f = (LHAPDF::xfx(x,Q,flavor) - LHAPDF::xfx(x,Q,-flavor))/x;

        //  take only gluon for gluon momentum fraction
        flavor = 0;
        //f = LHAPDF::xfx(x,Q,flavor);

        // for simple integration
        // ff = f*(xmax-xmin);
        
        // for importance sampling
        //      divide f(x) with g(x) = 1/x since we generate x according to g(x). 
        double ff = f*x*log(xmax/xmin);
        sum0  +=  ff;
        sum00 +=  ff*ff; 
        //                                  
    }

    sum0  /= npoints;
    sum00 /= npoints;
    double sigma2 = sum00 - sum0*sum0 ;
    double error = sqrt(sigma2/npoints) ;
    cout<<" mom sum rule is: "<<sum0<<"+/-"<< error<<endl;

	return EXIT_SUCCESS;

}
