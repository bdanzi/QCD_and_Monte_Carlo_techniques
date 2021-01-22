#include "courselib.h"
#include <cmath>
#include "TH1.h"
#include "TFile.h" 
#include "TCanvas.h"
#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>

#include <iostream>

using namespace std;

/*
   Monte Carlo integration.
   We calculate the Sudakov form factor
   and plot it for different upper scales as function of pt

   The random number generator is RANLUX by M. Luescher.
   M. Lüscher, A portable high-quality random number generator 
   for lattice field theory simulations, 
   Computer Physics Communications 79 (1994) 100  
   http://luscher.web.cern.ch/luscher/ranlux/index.html      


*/


//Alpha strong
double alphas (double q)
{
    double QCDlam = 0.2; //Lambda QCD for 3 flavours
    double Qlam0  = 1;    //scale freezing
    int    nf     = 3;    //number of flavours
    double beta0  = (33 - 2*nf) / 6; 

    double Qval = max(Qlam0, q);
    double d1   = Qval/QCDlam;
    double as   = M_PI / (beta0*log(d1));
    //      cout << " alphas = " << as << " q = " << q << endl;
    return as;
}

// gluon -> glouon or quark -> quark splitting
double Splitting (double z )
{
    double temp;
    // Pgg
    //temp = 6*(1/z -2 +z*(1-z) + 1/(1-z));
    // Pqq
    temp = 4./3.*((1.+z*z)/(1.-z));
    if (z > 0.99999) {
        cout <<" Splitting large z = " << z<< " P = "<< temp<< endl;
        temp = 0;
    }
    if (z < 0.000001) {
        cout <<" Splitting small z = " << z<< " P = "<< temp<< endl;
        temp = 0;
    }
    return temp;
}


double suda (double t1, double t2, double x, double y)
{
    double q0 = 0.1;
    t1 = max(q0, t1);

    if ( x >= 1 ){ cout << " Suda: argument 1 out of range: "<< x << endl;}
    if ( x <= 0 ){ cout << " Suda: argument 1 out of range: "<< x << endl;}
    if ( y >= 1 ){ cout << " Suda: argument 2 out of range: "<< y << endl;}
    if ( y <= 0 ){ cout << " Suda: argument 2 out of range: "<< y << endl;}
    //     	 cout << "Suda:  t2 = "<< t2 << " t1 = "<< t1 << endl;

    double d1 = t2/t1;
    double q2 = t1*pow(d1,x);
    double q  = sqrt(q2);
    // we generate here z1 = 1-z, 
    // because we have a pole in the splitting functions \sim 1/(1-z)       
    double z1min = 0.01;
    double z1max = 0.99;
    double d2 = z1max/z1min;
    double z1= z1min*pow(d2,y);

    // cout << " z1 "<<z1<< "zmin " <<z1min << "zmax " <<z1max<< endl;
    // cout << " d2 "<<d2<< " y  " << y << endl;
    if (z1max <= z1min) { z1 = 0;}

    double result = 0;
    if(z1 > 0) {
        double z  = 1. - z1;
        double d3 = q;
        double temp = alphas(d3)/2/M_PI * Splitting(z) /q2;	 
        result = temp*q2*log(t2/t1)*z1* log(z1max/z1min);
        // cout << " result :" << result <<endl;
    }
    // cout << " z1 " << z1 << " t " << q2 << " result " << result<< endl; 
    result = max(0., result);

    return result; 
}


int main (int argc,char **argv)
{
    TH1::SetDefaultSumw2();
    TApplication* gMyRootApp = new TApplication("My ROOT Application", &argc, argv);
    const int npoints = 100000;
    const int ntmax = 20;
    const double tmin = 1., tmax = 500.; 

    // initialise random number generator: rlxd_init( luxory level, seed )
    rlxd_init(2,32767);


    // book the histogram: TH1D("label","title",nr of bins, xlow,xhigh )
    TH1D *histo1 = new TH1D("sudakov","sudakov",ntmax, tmin, tmax);

    // loop over t1
    const double delta = (tmax-tmin)/ntmax;

    
    for (int  nt = 0; nt < ntmax; ++nt) {  
        double sum0 = 0, sum00 = 0;
        double t1 = tmin + delta*(nt+0.5);
        double t2 = tmax; // select here the upper scale t2 = tmax
        // cout << " tmax = "<< t2 << " t1 = "<< t1 << " delta " << delta<< endl;
        for (int n1 = 0; n1 < npoints; ++n1) {
            double x1 = Rand();
            double y1 = Rand();
            double ff = suda(t1, t2, x1, y1);
            sum0  +=  ff;
            sum00 +=  ff*ff; 
        }
        //                                  
        sum0  /= npoints;
        sum00 /= npoints;
        double sigma2 = sum00 - sum0*sum0;
        double error = sqrt(sigma2/npoints);

        double sudakov = exp(-sum0);
        double sudError = sudakov*error; //Error of the sudakov
        cout << " t2 = "<< t2 << " t1 = "<< t1 << " Delta_S = " << sudakov << " +-" << sudError << endl;
        //histo1->Fill(t1+0.001*delta,sudakov);
        histo1->SetBinContent(nt+1, sudakov);
        histo1->SetBinError(nt+1, sudError);
    }

    gStyle->SetPadTickY(1); // ticks at right side
    gStyle->SetOptStat(0); // get rid of statistics box

    TCanvas *c = new TCanvas("ctest", "", 0, 0, 500, 500);
    gPad->SetLogy();
    histo1->Draw();
    c->Draw();
    c->Print("example6.pdf");

    // write histogramm out to file
    TFile file("output-example6.root","RECREATE");
    histo1->Write();
    file.Close();
    gMyRootApp->Run();
    return EXIT_SUCCESS;
}
