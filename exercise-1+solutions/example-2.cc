//#include "ranlxd.h"
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
   Cosntructing a gauss random number generator, by adding 
   several random numbers.

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
    const int npoints = 1000000;

    //Init histograms for given number of points in sum
    int pointsA[] = {1, 2, 3, 6, 12, 20, 40};
    vector<int> points(pointsA, pointsA + sizeof(pointsA)/sizeof(pointsA[0]));
    vector<TH1D*> histos(points.size());
    for(unsigned i = 0; i < points.size(); ++i) {
        TString n = TString::Format("r%d",points[i]);
        histos[i] = new TH1D(n, n, 100, -6, 6);
    }


    // initialise random number generator: rlxd_init( luxory level, seed )
    rlxd_init(2,32767);
	for (int n1 = 0; n1 < npoints; n1++ ) {

        // here calculate the sum over the random numbers and fill them into the histos
        // calcualte it for the sum over 1,3,6,10,12,20,40 random numbers
        // in root then fit the distribution and determine the mean and sigma

        for(unsigned iH = 0; iH < histos.size(); ++iH) {
            int nP = points[iH];
            double x = 0.;
            for(int n2 = 0; n2 < nP; ++n2)
                x += Rand();
            histos[iH]->Fill((x - nP/2.)/sqrt(nP/12.));
        }

    }

    cout<<" Gauss Random number generator "<< endl;
    TCanvas *c = new TCanvas("ctest", "" , 0, 0, 1000, 500);
    // divide the canvas in 1 parts in x and 1 in y
    c->Divide(2,3);
    //enter the first part of the canvas (upper left)

    for(unsigned iH = 0; iH < points.size(); ++iH) {
        c->cd(iH+1);
        histos[iH]->SetLineColor(iH+1);
        histos[iH]->DrawCopy("hist e");
    }
    c->Print("example2.pdf(");

    //Plot normalized into single Graph
    TCanvas *d = new TCanvas("dtest", "" , 0, 0, 1000, 500);
    for(unsigned iH = 0; iH < points.size(); ++iH) {
        histos[iH]->Scale(1./histos[iH]->Integral());
        histos[iH]->Draw(iH == 0 ? "hist e": "hist e same");
    }
    histos[0]->SetMaximum(0.05);
    d->BuildLegend();

    d->Print("example2.pdf)");


    // write histogramm out to file
    TFile file("output-example2.root","RECREATE");
    for(unsigned iH = 0; iH < points.size(); ++iH)
        histos[iH]->Write();

    file.Close();
    gMyRootApp->Run();
    return EXIT_SUCCESS;
}
