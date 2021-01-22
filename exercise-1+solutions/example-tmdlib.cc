#include "courselib.h"
#include <cmath>
#include "TH1.h"
#include "TFile.h" 
#include "TCanvas.h"
#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>
#include <string>
#include <iostream>
#include <cstdlib>
#include "TMDlib.h"

using namespace std;
using namespace TMDlib;
TMD::TMD(void) {
   cout << "TMD object is being created "  << endl;    
    }


/*
   We calculate the integral ove kt of a TMD parton distribution

   We use TMDlib 2.0.1

   The random number generator is RANLUX by M. Luescher.
   M. Lüscher, A portable high-quality random number generator 
   for lattice field theory simulations, 
   Computer Physics Communications 79 (1994) 100  
   http://luscher.web.cern.ch/luscher/ranlux/index.html    
   
   Authors: H. Jung, R. Zlebcik 
  
*/

int main (int argc,char **argv)
{
	TH1::SetDefaultSumw2();
	TApplication* gMyRootApp = new TApplication("My ROOT Application", &argc, argv);
	int nx = 100;
	double xmin = 0.0001 ;
	double xmax = 0.9;
    
	double logxmin = log10(xmin) ;
	double logxmax = log10(xmax);
    
	TH1D *histo1 = new TH1D("TMDint","TMDint",nx, logxmin, logxmax);

	const int npoints = 100000;
	double up, ubar, down, dbar, strange, sbar, charm, cbar, bottom, bbar, gluon;

      string name;
      //name="PB-NLO-HERAI+II-2018-set1";
      name="PB-NLO-HERAI+II-2018-set2";
      TMD TMDmy;

      TMDmy.TMDinit(name);
      cout << " TMDSet Description: " << TMDmy.TMDgetDesc() << endl;
      cout << " TMDlib: Order alphas = " << TMDmy.TMDgetOrderAlphaS() << endl;
      cout << " TMDlib: order PDF    = " << TMDmy.TMDgetOrderPDF() << endl;
      cout << " TMDlib: xmin         = " << TMDmy.TMDgetXmin() << endl;
      cout << " TMDlib: xmax         = " << TMDmy.TMDgetXmax() << endl;
      cout << " TMDlib: q2min        = " << TMDmy.TMDgetQ2min() << endl;
      cout << " TMDlib: q2max        = " << TMDmy.TMDgetQ2max() << endl;
      cout << " TMDlib: number       = " << TMDmy.TMDnumberPDF(name) << endl;
      cout << endl;


//  initialise random number generator
//  rlxd_init( luxory level, seed )
    rlxd_init(2,32767);


	double sum0=0, sum00=0;
    
    
	double kt2min = 0.01;
      double kt2max= 100.;    
	double x = 0.01;
	double xbar = 0.;
	double q = 10.;
      
      //cout << " " << logxmin << " " << logxmax << endl;
      
      for (int i1= 1; i1 <= nx; ++i1) {
          
		double logx = logxmin + i1 * (logxmax-logxmin)/nx;
            // cout << logx << endl ;
            x = pow(10,logx) ; //inverse formula for obtaining correct x
            sum0 = 0.;
            sum00 = 0.;
            double sigma2 = 0.;
            double error = 0.;
          
		for (int n1 = 0; n1 < npoints; ++n1) {
		// calculate the integral over kt
            
		// for importance sampling

		double kt2 = kt2min * pow((kt2max/kt2min), Rand());
            double weightkt2 = kt2 * log(kt2max/kt2min);
            double kt = sqrt(kt2) ;
		// get TMD with the correct x,kt,q
       	TMDmy.TMDpdf(x, xbar, kt, q, up, ubar, down, dbar, strange, sbar,  charm, cbar, bottom, bbar, gluon);

		double ff = gluon*weightkt2;
            //cout << " x = " << x<<" gluon " << gluon << " " << weightkt2 << endl; 
		sum0  +=  ff;
		sum00 +=  ff*ff; 
		//                                  
		}
		sum0  /= npoints;
		sum00 /= npoints;
		sigma2 = sum00 - sum0*sum0 ;
		error = sqrt(sigma2/npoints) ;
		// cout<<" integral at " << x <<" is: "<<sum0<<"+/-"<< error<<endl;
        	histo1->SetBinContent(i1, sum0);
        	histo1->SetBinError(i1, error);
	}
	gStyle->SetPadTickY(1); // ticks at right side
	gStyle->SetOptStat(0); // get rid of statistics box

	TCanvas *c = new TCanvas("ctest", "", 0, 0, 500, 500);
	gPad->SetLogy();
	histo1->Draw();
      c->Draw();
	c->Print("example-TMDlib.pdf");

	// write histogramm out to file
	TFile file("output-example-TMDlib.root","RECREATE");
	histo1->Write();
	file.Close();
	gMyRootApp->Run();
	return EXIT_SUCCESS;

}
