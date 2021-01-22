#include "courselib.h"
#include <cmath>
#include "TH1.h"
#include "TFile.h" 
#include "TCanvas.h"
#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>
#include <cstdlib>
#include "TLorentzVector.h"

using namespace std;

/*
   Monte Carlo integration.
   We calculate the x-section for p p -> Higgs  + X
   and plot the distribution in pt and eta

   The random number generator is RANLUX by M. Luescher.
   M. Lüscher, A portable high-quality random number generator 
   for lattice field theory simulations, 
   Computer Physics Communications 79 (1994) 100  
   http://luscher.web.cern.ch/luscher/ranlux/index.html      

*/


double sigma (double m2, double sh, double th)
{	
    // calculate sigma0 (qq -> gamma^*)
    // double aem=1./137.;
    // result = 4.*pi*pi*aem/3.;

    // calculate sigma0(gg->H) 
    // sigma0 = as^2/pi /576 * GF * sqrt(2)*\Delta(1.-tau) 
    const double as = 0.1;   //alphaS
    const double GF = 1.166E-5; //Fermi constant
    double result = as * as * sqrt(2.) * GF /M_PI / 576. ;
    //      cout << " M = " << sqrt(m2) << endl;
    return result; 
}


void gauss2D (double sigma, double &kx, double &ky)
{
    double kT  = sigma * sqrt(-2*log(Rand()));
    double phi = 2*M_PI * Rand();
    kx = kT*cos(phi);
    ky = kT*sin(phi);
}


void getpdf (double xmin, double q2, double& weightx, double& x, double& kx, double& ky)
{
    double xmax = 0.999;

    x = xmin*pow(xmax/xmin, Rand());
    weightx =  x*log(xmax/xmin) ;
    // this is for the simple case            
    double pdf = 3.*pow((1-x),5)/x;

    weightx = weightx * pdf;

    gauss2D(0.7, kx, ky);
}

int main (int argc,char **argv)
{
    TH1::SetDefaultSumw2();
    TApplication* gMyRootApp = new TApplication("My ROOT Application", &argc, argv);
    const int npoints = 1.0e6;
    const double x1min = 0.00001, x2min = 0.00001;
    const double m_higgs = 125., Gamma = 0.4; /* Width from PDG 2014 */
    const double s = 4*3500.*3500.;

    const double gev2nb = 0.3893857E+6;

    // book the histogram TH1D("label","title",nr of bins, xlow,xhigh )
    TH1D *histo1  = new TH1D("x1","x1",100, -5, 0.);
    TH1D *histo2  = new TH1D("x2","x2",100, -5, 0.);
    TH1D *histo3  = new TH1D("kt1 ","kt1 ",100, 0, 10.);
    TH1D *histo4  = new TH1D("kt2 ","kt2 ",100, 0, 10.);
    TH1D *histo5  = new TH1D("pt ","pt ",50, 0, 10.);
    TH1D *histo6  = new TH1D("eta ","eta",50, -8, 8.);
    TH1D *histo7  = new TH1D("Mass ","mass",50, 60., 160.);

    // initialise random number generator: rlxd_init( luxory level, seed )
    rlxd_init(2,32767);

    const double q2 = 10000;
    const double mass_min = 0.;
    const double mass_max = 200.;

    double sum0 = 0, sum00 = 0;
    int nacc = 0;

    for (int n1 = 0; n1 < npoints; ++n1) {
        // this is using importance sampling
        // generate starting distribution in x and kt
        double kx, ky, weightx1, weightx2;
        double x1, x2;

        // define incoming partons: parton 1
        getpdf(x1min, q2, weightx1, x1, kx, ky);
        TLorentzVector pA;
        pA.SetXYZM(kx, ky, sqrt(s)/2.*x1, 0.);


        // define incoming partons: parton 2
        getpdf(x2min, q2, weightx2, x2, kx, ky);
        TLorentzVector pB;
        pB.SetXYZM(kx, ky, -sqrt(s)/2.*x2, 0.);

        // plot dxg(x)/dlogx *Jacobian, Jacobian dlogx/dx = 1/x
        histo1->Fill(log10(x1),weightx1/log(10));
        histo2->Fill(log10(x2),weightx2/log(10));
        histo3->Fill(pA.Pt(), weightx1);
        histo4->Fill(pB.Pt(), weightx2);


        TLorentzVector pH = pA + pB;
        double mass = pH.M();

        // calculate rapidity of Higgs
        //double rapidity = 0.5 * log(x1/x2) ;
        //         if (mass > mass_min && mass < mass_max && abs(rapidity) < 3 ) {
        if (mass > mass_min && mass < mass_max ) {

            double sh=0., th=0.;
            // 		note, pdfs are already included in weightx            
            double Hx1x2 = weightx1*weightx2;
            Hx1x2 = Hx1x2*sigma(mass*mass, sh, th) ;
            // multiply with a Breit Wigner resonance
            Hx1x2 = Hx1x2 * Gamma/(pow((mass-m_higgs),2) + pow(Gamma,2)/4)/2./M_PI;
            Hx1x2 = Hx1x2*gev2nb;         
            if (x1 < 1) {
                ++nacc;
                double ff=Hx1x2; 
                sum0  +=  ff;
                sum00 +=  ff*ff; 
                // weighting with 1/x0:
                // plot dxg(x)/dlogx *Jacobian, Jacobian dlogx/dx = 1/x

                histo5->Fill(pH.Pt(), ff);
                histo6->Fill(pH.Rapidity(), ff);
                histo7->Fill(mass,  ff);
            }
        }          
    }
    //                                  
    sum0  /= npoints;
    sum00 /= npoints;
    double sigma2 = sum00 - sum0*sum0 ;
    double error = sqrt(sigma2/npoints) ;
    cout << " nr of events accepted: "<< nacc << endl;
    cout<< " integral for Higgs xsection is [pb]: " << sum0 * 1000.<< " +/- " << error*1000.<< endl;

    // write histogramm out to file
    TFile file("output-example8.root","RECREATE");
    //general root settings 
    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gStyle->SetPadTickY(1); // ticks at right side
    gStyle->SetOptStat(0); // get rid of statistics box
    TCanvas *c = new TCanvas("ctest", "" ,0, 0, 500, 500);
    // divide the canvas in 1 parts in x and 1 in y
    c->Divide(3,3);
    c->cd(1);
    histo1->Scale(1./npoints, "width");
    histo1->Draw();
    c->cd(2);
    histo2->Scale(1./npoints, "width");
    histo2->Draw();
    c->cd(3);
    histo3->Draw();
    c->cd(4);
    histo4->Draw();
    c->cd(5);
    histo5->Draw();
    c->cd(6);
    histo6->Draw();
    c->cd(7);
    histo7->Draw();
    c-> Draw();
    c-> Print("example8.pdf");

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
