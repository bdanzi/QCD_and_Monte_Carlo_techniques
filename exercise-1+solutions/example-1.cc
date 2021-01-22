#include "courselib.h"
#include "TH2.h"
#include "TFile.h" 
#include "TCanvas.h"
#include <TApplication.h>
#include <TStyle.h>
#include <TROOT.h>
#include <iostream>
using namespace std;

/*
   Cosntructing a random number generator

   Authors: H. Jung, R. Zlebcik 
   
*/


int main (int argc,char **argv)
{
    TApplication* gMyRootApp = new TApplication("My ROOT Application", &argc, argv);
    const int npoints = 100000;
    const long long ia=205, ic=29573, im=139968;
    long long last=4711;

    const double xlow=0.,xup=1.,ylow=0.,yup=1;
    //  TH2D("label", "title", nBinsX, xlow, xhigh,  nBinsY, ylow, yhigh )
    TH2D *histo1 = new TH2D("congruental random numbers","congruental random numbers ",100, xlow, xup, 100, ylow, yup);
    TH2D *histo2 = new TH2D("RANLUX","RANLUX",100, xlow, xup, 100, ylow, yup);

    // initialise random number generator: rlxd_init( luxory level, seed )
    rlxd_init(2,32767);
    for (int n1 = 0; n1 < npoints; ++n1 ) {
        //   construct random numbers according the formula in exercise 1
        //   construct here two random numbers according to the congruential method
        last = (ia*last+ic) % im;
        //cout << ia*last+ic << endl;
        // cout << " i1 = "<< ia*last+ic << " im = " << im << " result =" << last << endl; 
        //            if ( last == 0. ){ last = fmod(ia*last+ic,im);}
        double xC = double(last)/im;
        //last = fmod(ia*last+ic,im);
        last = (ia*last+ic) % im;
        // if ( last == 0. ){ last = fmod(ia*last+ic,im);}
        double yC = double(last)/im;
        // fill the 2 random numbers into a scatter plot to compare the correlation
        histo1->Fill(xC, yC);
        // cout << " x,y " << x << " " <<  y  << endl;             
        // next use the ranlux generator and compare also the 2 random numbers

        double xRL = Rand();
        double yRL = Rand();
        histo2->Fill(xRL, yRL);
    }

    cout<<" Congruential Random number generator "<< endl;
    TCanvas *c = new TCanvas("ctest", "" ,0, 0, 1000, 500);
    // divide the canvas in 1 parts in x and 1 in y
    c->Divide(2,1);
    //enter the first part of the canvas (upper left)
    c->cd(1);
    histo1->Draw();

    c->cd(2);
    histo2->Draw();

    c->Print("example1.pdf");

    // write histogramm out to file
    TFile file("output-example1.root","RECREATE");
    histo1->Write();
    histo2->Write();
    file.Close();

    gMyRootApp->Run();
    return EXIT_SUCCESS;
}

