#include "ranlxd.h"
#include <cmath>
#include "TH1.h"
#include "TFile.h" 
#include "TCanvas.h"
#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>
#include <cstdlib>

using namespace std;

/*
   Monte Carlo integration.
   We calculate the x-section for p p -> Higgs  + X
   and plot the distribution in pt and rapidity

   We use our own evolved pdf.

   The random number generator is RANLUX by M. Luescher.
   M. Lüscher, A portable high-quality random number generator 
   for lattice field theory simulations, 
   Computer Physics Communications 79 (1994) 100  
   http://luscher.web.cern.ch/luscher/ranlux/index.html      

*/

typedef struct {
    double e;
    double px,py,pz;
    double m;
} PVEC;
#define DOTPR(p1,p2) (p1.e*p2.e-p1.px*p2.px-p1.py*p2.py-p1.pz*p2.pz)

void getpdf (double xmin, double q2,double& weightx, double& x, double& kx, double& ky);

double sigma(double m2, double sh, double th)
{	
	static const double pi = 3.14159265358979323848;
	double result;
      result = 0;
      // calculate sigma0 (qq -> gamma)
	// double aem=1./137.;
	// result = 4.*pi*pi*aem/3.;
      
      // calculate sigma0(gg->H) 
      // sigma0 = as^2/pi /576 * GF * sqrt(2)*\Delta(1.-tau) 
      double const as=0.1, GF=1.166E-5;
      result = as * as * sqrt(2.) * GF /pi / 576.;
	return(result); 
}

void getnorm (double& norm);

int main(int argc,char **argv)
{
    TApplication* gMyRootApp = new TApplication("My ROOT Application", &argc, argv);
    int n1, npoints = 1.E7, nacc=0;
    double x1,x2;
    double x1min=0.0001, x2min=0.0001;
    double Hx1x2, mass,sh=0,th=0,q2;
    double mass_min, mass_max;
    double m_higgs = 125., Gamma=0.4; /* Width from PDG 2014 */
    static const double pi = 3.14159265358979323848;
    double s=4*3500.*3500.,stest;
    double pt2, rapidity;
    double weightx1, weightx2;
    double sum0, sum00, ff, sigma2, error;
    double kx,ky;
    double a,b;


    static const double gev2nb=0.3893857E+6;

    PVEC pin1,pin2, p_a, p_b, pout;

    // book the histogram TH1F("label","title",nr of bins, xlow,xhigh )
    TH1F *histo1 = new TH1F("x1","x1",100, -4, 0.);
    TH1F *histo10 = new TH1F(*histo1); 
    double binwidth1=4./100.;
    TH1F *histo2 = new TH1F("x2","x2",100, -4, 0.);
    TH1F *histo20 = new TH1F(*histo2); 
    TH1F *histo3 = new TH1F("kt1 ","kt1 ",100, 0, 100.);
    TH1F *histo4 = new TH1F("kt2 ","kt2 ",100, 0, 100.);
    TH1F *histo5 = new TH1F("pt ","pt ",50, 0, 200.);
    TH1F *histo6 = new TH1F("eta ","eta",50, -8, 8.);
    TH1F *histo7 = new TH1F("Mass ","mass",50, 60., 160.);


    // initialise random number generator: rlxd_init( luxory level, seed )
    rlxd_init(2,132144);

    sum0 = sum00 = 0;
    // define incoming protons      
    pin1.px = 0;
    pin1.py = 0;
    pin1.pz = sqrt(s)/2;
    pin1.e = sqrt(s)/2;

    pin2.px = 0;
    pin2.py = 0;
    pin2.pz = -sqrt(s)/2;
    pin2.e = sqrt(s)/2;

    stest = sqrt(2*DOTPR(pin1,pin2));

    cout << " sqrt(s) = " << stest << endl;

    q2 = 10000;
    mass_min = 0.;
    mass_max = 200.;

    for ( n1 = 0; n1 < npoints; n1++ ) {
        // this is using importance sampling
        // generate starting distribution in x and kt
        // define incoming partons: parton 1
        getpdf(x1min, q2, weightx1, x1, kx, ky);
        p_a.px = kx;
        p_a.py = ky;
        p_a.pz = sqrt(s)/2.*x1;
        p_a.e = sqrt(p_a.px*p_a.px + p_a.py*p_a.py + p_a.pz*p_a.pz);
        double kt21 = p_a.px*p_a.px+p_a.py*p_a.py;
        getpdf(x2min, q2, weightx2, x2, kx, ky);

        // define incoming partons: parton 2
        p_b.px = kx;
        p_b.py = ky;
        p_b.pz = -sqrt(s)/2.*x2;
        p_b.e = sqrt(p_b.px*p_b.px + p_b.py*p_b.py + p_b.pz*p_b.pz);
        double kt22 = p_b.px*p_b.px+p_b.py*p_b.py;

        // plot dxg(x)/dlogx *Jacobian, Jacobian dlogx/dx = 1/x
        histo1->Fill(log10(x1),weightx1/2.3026);
        histo2->Fill(log10(x2),weightx2/2.3026);
        histo3->Fill(sqrt(kt21),weightx1);
        histo4->Fill(sqrt(kt22),weightx2);


        // calculate 4-vector of DY   
        pout.px = p_a.px + p_b.px;
        pout.py = p_a.py + p_b.py;
        pout.pz = p_a.pz + p_b.pz;
        pout.e = p_a.e + p_b.e;
        pt2 = pout.px*pout.px + pout.py*pout.py;
        mass = sqrt(DOTPR(pout,pout));

        // calculate rapidity of Higgs
        rapidity = 0.5 * log(x1/x2) ;
        //         if (mass > mass_min && mass < mass_max && abs(rapidity) < 3 ) {
        if (mass > mass_min && mass < mass_max ) {

            // 		note, pdfs are already included in weightx            
            Hx1x2 = weightx1*weightx2;
            Hx1x2 = Hx1x2*sigma(mass*mass, sh, th) ;
            Hx1x2 = Hx1x2 * Gamma/(pow((mass-m_higgs),2) + pow(Gamma,2)/4)/2./pi;
            Hx1x2 = Hx1x2*gev2nb;         
            if ( x1 < 1) {
                nacc = nacc + 1;
                ff=Hx1x2; 
                sum0 = sum0 + ff;
                sum00 = sum00 + ff*ff; 
                // weighting with 1/x0:
                // plot dxg(x)/dlogx *Jacobian, Jacobian dlogx/dx = 1/x

                histo5->Fill(sqrt(pt2),ff);
                histo6->Fill(rapidity,ff);
                histo7->Fill(mass,ff);
            }
        }          
    }
    //             
    double norm;                     
    getnorm(norm);
    cout << " norm " << norm << endl;


    sum0 = sum0/npoints;
    sum00 = sum00/npoints;
    sigma2 = sum00 - sum0*sum0 ;
    error = sqrt(sigma2/npoints) ;
    cout << " nr of events accepted: "<< nacc << endl;
    cout<< " integral for Higgs xsection is [pb]: " << sum0 * 1000.<< " +/- " << error*1000.<< endl;
    // write histogramm out to file
    TFile file("output-example9.root","RECREATE");
    //general root settings 
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetPadTickY(1); // ticks at right side
    gStyle->SetOptStat(0); // get rid of statistics box
    TCanvas *c = new TCanvas("ctest", "" ,0, 0, 500, 500);
    // divide the canvas in 1 parts in x and 1 in y
    c -> Divide(3,3);
    c -> cd(1);
    a=1./npoints/binwidth1; 
    b=0; 
    histo10->Add(histo1,histo1,a,b); 
    histo10 -> Draw();
    c -> cd(2);
    a=1./npoints/binwidth1; 
    b=0; 
    histo20->Add(histo2,histo2,a,b); 
    histo20 -> Draw();
    c -> cd(3);
    histo3 -> Draw();
    c -> cd(4);
    histo4 -> Draw();
    c -> cd(5);
    histo5 -> Draw();
    c -> cd(6);
    histo6 -> Draw();
    c -> cd(7);
    histo7 -> Draw();
    c-> Draw();
    c-> Print("example9.pdf");

    histo1->Write();
    histo2->Write();
    histo3->Write();
    histo4->Write();
    histo5->Write();
    histo6->Write();
    histo7->Write();
    file.Close();

    gMyRootApp->Run();
    return EXIT_SUCCESS;
    }
