//Generate points distributed with some errors around a circle
//Fit a circle through the points and draw
//To run the script, do, eg
//   root > .x fitCircle.C   (10000 points by default)
//   root > .x fitCircle.C(100);  (with only 100 points
//   root > .x fitCircle.C(100000);  with ACLIC
//
//Author: Rene Brun

#include "TCanvas.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TMath.h"
#include "TArc.h"
#include "TVirtualFitter.h"

TGraph *gr;

//____________________________________________________________________
void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
   //minimisation function computing the sum of squares of residuals
   Int_t np = gr->GetN();
   f = 0;
   Double_t *x = gr->GetX();
   Double_t *y = gr->GetY();
   for (Int_t i=0;i<np;i++) {
      Double_t u = x[i] - par[0];
      Double_t v = y[i] - par[1];
      Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
      f += dr*dr;
   }
}

//____________________________________________________________________
void fitCircle(Int_t n=10000) {
   //generates n points around a circle and fit them
   TCanvas *c1 = new TCanvas("c1","c1",600,600);
   c1->SetGrid();
   gr = new TGraph();
   if (n> 999) gr->SetMarkerStyle(1);
   else        gr->SetMarkerStyle(3);
   gr->SetMarkerStyle( 21);
   TRandom3 r;
   for( int iloop = 0; iloop < 100; iloop++){
     gr->Set(0);
     Double_t x,y;
     for (Int_t i=0;i<4;i++) {
       r.Circle(x,y,r.Gaus(4,0.1));
       gr->SetPoint(i,x,y);
       std::cout<< i << "\t" << x << "\t" << y << std::endl;
     }
     c1->DrawFrame(-5,-5,5,5);
     gr->Draw("p");

     //Fit a circle to the graph points
     TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
     TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
     fitter->SetFCN(myfcn);

     fitter->SetParameter(0, "x0",   0, 0.1, 0,0);
     fitter->SetParameter(1, "y0",   0, 0.1, 0,0);
     fitter->SetParameter(2, "R",    1, 0.1, 0,0);

     Double_t arglist[1] = {0};
     fitter->ExecuteCommand("MIGRAD", arglist, 0);
     std::cout<< fitter->GetParameter(0) << "\t" << fitter->GetParameter(1) << "\t" << fitter->GetParameter(2) << std::endl;
     //Draw the circle on top of the points
     TArc *arc = new TArc(fitter->GetParameter(0),
			  fitter->GetParameter(1),
			  fitter->GetParameter(2));
     arc->SetLineColor(kRed);
     arc->SetLineWidth(4);
     arc->SetFillStyle(0);
     arc->Draw();
     c1->Update();
     c1->Modified();
     getchar();
     gSystem->ProcessEvents();
   }
}
