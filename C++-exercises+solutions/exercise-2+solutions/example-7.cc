#include "courselib.h"
#include <cmath>
#include "TH1.h"
#include "TFile.h" 
#include "TCanvas.h"
#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TF1.h>
#include <iostream>

using namespace std;

/*
   Parton evolution:
     1. generate starting distribution in x and add intrisnic transverse mom.
     2. evolve sarting distribution to large values of t

   Plot starting distribution and intrinsic kt
   Plot evolved distribution as function of x and kt

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
    const double xmax = 1;

    // get x value according to g(x)\sim 1/x            
    x = xmin*pow(xmax/xmin, Rand());
    // use: xg(x) = 3 (1-x)**5

    double pdf     = pow(1-x, 5) * 3./x;
    weightx = x*log(xmax/xmin) * pdf ;

    // now generate instrinsic kt according to a gauss distribution  
    gauss2D(0.7, kx, ky);
    // cout << "in getpdf:  x " << x << " kx " << kx <<" ky " << ky <<endl; 
}

void sudakov (double t0, double& t)
{
//  here we calculate  from the sudakov form factor the next t > t0
    double epsilon = 0.1;
    double as2pi = 0.1/2./M_PI, Ca=3.;
    double Pint;
// for Pgg use fact = 2*Ca

    double fac=2.*Ca; 
    double t2;
    // use fixed alphas and only the 1/(1-z) term of the splitting fct

    double r1 = Rand();
    Pint=log((1.-epsilon)/epsilon); /* for 1/(1-z) splitting fct */
    t2 = -log(r1)/fac/as2pi/Pint;
    t2 = t0 * exp(t2);
    t = t2;   
    // cout << " t0= "<<t0<< " t= "<<t2<< " epsilon= "<<epsilon <<endl;
    if (t2 < t0) {        
        cout << "FATAL t0 > t:  t0= "<<t0<< "t= "<<t2<< "epsilon= "<<epsilon <<endl;
    }
}

void splitting(double& z, double& weightz)
{
    const double epsilon = 0.1;

    double as2pi = 0.1/2./M_PI;

//		here we calculate  the splitting variable z for 1/z and  1/(1-z)
//		use Pgg=6 (1./z + 1./(1.-z))

    double g0 = 6.*as2pi * log((1.-epsilon)/epsilon);
    double g1 = 6.*as2pi * log((1.-epsilon)/epsilon);
    double gtot = g0 + g1 ;

    double zmin = epsilon;
    double zmax = 1.-epsilon;

    double r1 = Rand();
    double r2 = Rand();

    z = zmin*pow((zmax/zmin),r2);
    if (r1 > g0/gtot)
        z = 1. - z;
//
    weightz = 1.;		
}

void evolve_pdf(double xmin, double q2, double x0, double kx0, double ky0, double& weightx, double& x, double& kx, double& ky)
{
    double t0=1.,t1, tcut;
    double z=0, weightz;
    double ratio_splitting;
    x = x0;
    kx = kx0;
    ky = ky0;
    weightx = 1. ;
    t1 = t0 ;
    tcut = q2;
    while (t1 < tcut ) {
        // here we do now the evolution
        // first we generate t
        t0 = t1;
        sudakov(t0, t1) ;
        // now we generate z
        splitting(z, weightz);
        /* 
           since the sudakov involves only the 1/(1-z) part 
           of the splitting fct, we need to weight each branching
           by the ratio of the intgral of the full and 
           approximate splitting fct
         */
        ratio_splitting = 2; /* for using Pgg*/

        if ( t1 < tcut ) { 
            x = x*z;
            weightx = weightx *weightz*ratio_splitting;
            /* 
               use here the angular ordering condition: sqrt(t1) = qt/(1-z) 
               and apply this also to the propagator gluons
             */
            double phi = 2*M_PI*Rand();  
            kx +=  sqrt(t1)*cos(phi)*(1.-z);
            ky +=  sqrt(t1)*sin(phi)*(1.-z);                     
            //                     kx = kx + sqrt(t1)*cos(phi);
            //                     ky = ky + sqrt(t1)*sin(phi);                     
        }
    }
    // cout << " evolve end  t= "<< t<< " q2 " <<q2 << " x0 = " <<x0 << "x =" << x <<endl;
}



int main(int argc,char **argv)
{
    TH1::SetDefaultSumw2();
    TApplication* gMyRootApp = new TApplication("My ROOT Application", &argc, argv);
    const int npoints = 1e6;

    const double xmin = 0.00001, q2=100*100;

    // initialise random number generator: rlxd_init( luxory level, seed )
    rlxd_init(2,32767);

    // book histograms: TH1D("label","title",nr of bins, xlow,xhigh )
    TH1D *histo1  = new TH1D("x0","x0",100, -5, 0.);
    TH1D *histo2  = new TH1D("kt0 ","kt0 ",100, 0, 10.);
    TH1D *histo3  = new TH1D("x","x",100, -5, 0.);
    TH1D *histo4  = new TH1D("kt ","kt ",1000, 0, 1000.);




    for (int n1 = 0; n1 < npoints; n1++ ) {

        double x , kx, ky, weightx ;
        double x0, kx0, ky0, weightx0 ;

        // generate starting distribution in x and kt
        get_starting_pdf(xmin, q2, weightx0, x0, kx0, ky0);
        // now do the evolution	
        evolve_pdf(xmin, q2, x0, kx0, ky0, weightx, x, kx, ky);   
        weightx = weightx0 * weightx;         

        if (x < 1) {
            // weighting with 1/x0:
            // plot dxg(x)/dlogx *Jacobian, Jacobian dlogx/dx = 1/x
            // log(x) = 2.3026 log10(x)
            double kt0 = sqrt(kx0*kx0 + ky0*ky0);
            double kt  = sqrt(kx*kx + ky*ky);

            histo1->Fill(log10(x0), weightx0/log(10));
            histo2->Fill(kt0,weightx0);
            histo3->Fill(log10(x),  weightx/log(10));
            histo4->Fill(kt, weightx);
        }
    }          
//                                  

    gStyle->SetPadTickY(1); // ticks at right side
    gStyle->SetOptStat(0); // get rid of statistics box
    TCanvas *c = new TCanvas("ctest", "", 0, 0, 500, 500);
    // divide the canvas in 2 parts in x and 2 in y
    c -> Divide(2,2);

    c->cd(1);
    histo1->Scale(1./npoints, "width");
    // gPad->SetLogy();
    histo1->Draw("hist e");

    c->cd(2);
    histo3->Scale(1./npoints, "width");
    histo3->Draw("hist e");


    c->cd(3);
    histo2->Scale(1./npoints, "width");
    gPad->SetLogy();
    histo2->Draw("hist e");


    c->cd(4);
    histo4->Scale(1./npoints, "width");
    gPad->SetLogy();
    histo4->Draw("hist e");

    c->Draw();
    c->Print("example7.pdf");

    // write histogramm out to file
    TFile file("output-example-7.root","RECREATE");
    histo1->Write();
    histo2->Write();
    histo3->Write();
    histo4->Write();
    file.Close();

    gMyRootApp->Run();
    return EXIT_SUCCESS;
}
