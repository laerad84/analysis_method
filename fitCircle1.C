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
Double_t RCut = 200;

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
      if( TMath::Abs(dr) < RCut ){
	f += dr*dr;
      }
   }
}

//____________________________________________________________________
void fitCircle1(Int_t n=300) {
   //generates n points around a circle and fit them
  TCanvas *c1 = new TCanvas("c1","c1",1800,1800);
  c1->Draw();
  c1->Divide(3,3);
  for( int i = 0; i< 9; i++){
    c1->cd(i+1);
    gPad->SetGrid();
    gPad->DrawFrame(-300,-300,300,300);
  }
  c1->Update();
  c1->Modified();
  gSystem->ProcessEvents();

  TGraph* grHit[10];
  for( int i = 0; i< 10;i++){
    grHit[i] = new TGraph();
    grHit[i]->SetMarkerColor(1);
    grHit[i]->SetMarkerStyle(3);
  }

  c1->Update();
  c1->Modified();
  //gSystem->ProcessEvents();

  TRandom3 r;
  Double_t x,y;
  Double_t x1,y1;
  Double_t r1=400;
  Double_t r2=300;
  for (Int_t i=0;i<n*2;i++) {
    r.Circle(x,y,r.Gaus(r1,3));
    x1 = x + 200;
    y1 = y + 80;
    //if( x1*x1+y1*y1 < 400){ continue; }
    //if( y1 > 0 ){ continue; }
    if( TMath::Abs(x1)>300 || TMath::Abs(y1)>300){ continue;}
    grHit[0]->SetPoint(grHit[0]->GetN(),x1,y1);
  }


  for( int i = 0; i< n*2; i++){
    r.Circle(x,y,r.Gaus(r2,3));
    x1 = x - r2;
    y1 = y;
    //if( x1*x1+y1*y1 < 10000){ continue; }
    //if( y1 > 0 ){ continue; }
    if( TMath::Abs(x1)>300 || TMath::Abs(y1)>300){ continue;}
    grHit[0]->SetPoint(grHit[0]->GetN(),x1,y1);
  }
  for( int i = 0; i< 200; i++){
    r.Circle(x,y,r.Gaus(r2,3));
    x1 = x -r2;
    y1 = y +r2;
    if( x1*x1+y1*y1 < 10000){ continue; }
    if( TMath::Abs(x1)>300 || TMath::Abs(y1)>300){ continue;}
    grHit[0]->SetPoint(grHit[0]->GetN(),x1,y1);
  }

  Double_t pseudo_fit[10][3]={{0}};
  pseudo_fit[0][0] = 0;
  pseudo_fit[0][1] = -270;
  pseudo_fit[0][2] = 310;

  pseudo_fit[1][0] = 100;
  pseudo_fit[1][1] = 200;
  pseudo_fit[1][2] = 500;

  pseudo_fit[2][0] = -200;
  pseudo_fit[2][1] = 200;
  pseudo_fit[2][2] = 200;


  for( int itrk = 0; itrk <9; itrk++){
    if( grHit[itrk]->GetN() < 3 ){ break; }
    gr = grHit[itrk];

    c1->cd(itrk+1);
    grHit[itrk]->Draw("p");
    //can->cd();
    //grHit[itrk]->Draw("p");
    const Int_t NFit = 10;
    Double_t fitPar[3]={pseudo_fit[itrk][0],pseudo_fit[itrk][1],pseudo_fit[itrk][2]};
    Double_t result[NFit][3]={{0}};
    Double_t RCutArr[NFit]={0};
    Double_t RInit = 500;
    for( int i = 0; i< NFit; i++){
      RCutArr[i] = RInit*TMath::Power(0.7,i);
      if( RCutArr[i] < 30){ RCutArr[i] = 30;}
    }
    Double_t RMSArr[NFit]={10,10,10,10,10,5,5,5,5,5};

    //Fit a circle to the graph points
    TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
    fitter->SetFCN(myfcn);

    for( int i = 0; i< NFit; i++){
      Double_t arglist[1] = {0};
      RCut = RCutArr[i];
      fitter->SetParameter(0, "x0",   fitPar[0], RMSArr[i], -1000,1000);
      fitter->SetParameter(1, "y0",   fitPar[1], RMSArr[i], -1000,1000);
      fitter->SetParameter(2, "R",    fitPar[2], RMSArr[i],200,1000);
      fitter->ExecuteCommand("MIGRAD", arglist, 0);

      for( int j = 0; j< 3; j++){
	fitPar[j] = fitter->GetParameter(j);
	result[i][j] = fitter->GetParameter(j);
      }
    }

    for( int i = 0; i< NFit; i++){
      std::cout<< result[i][0] << "\t" << result[i][1] << "\t" << result[i][2] << std::endl;
    }
    TArc* arc[NFit];
    for( int i = 0; i< NFit; i++){
      arc[i] = new TArc(result[i][0],result[i][1],result[i][2]);
      arc[i]->SetLineColor(i+1);
      arc[i]->SetFillStyle(0);
      arc[i]->SetLineWidth(2);
      arc[i]->Draw("same");
    }

    std::cout<< itrk << std::endl;
    c1->Update();
    c1->Modified();
    gSystem->ProcessEvents();
    std::cout<< "SetPoint to next " << std::endl;
    for( Int_t ip = 0; ip < grHit[itrk]->GetN();ip++){
      //std::cout<< ip << std::endl;
      Double_t dist = TMath::Sqrt(TMath::Power(grHit[itrk]->GetX()[ip]-result[NFit-1][0],2) +
				  TMath::Power(grHit[itrk]->GetY()[ip]-result[NFit-1][1],2));
      //std::cout<< ip << "\t" << dist << std::endl;
      if( TMath::Abs(dist-result[NFit-1][2]) > 30){
	grHit[itrk+1]->SetPoint(grHit[itrk+1]->GetN(),
				grHit[itrk]->GetX()[ip],
				grHit[itrk]->GetY()[ip]);}
    }

    std::cout<< "End SetPoint :" << itrk <<"\t" << grHit[itrk]->GetN() << "\t" << grHit[itrk+1]->GetN()<< std::endl;
    std::cout<< "Fit result   :"
	     << result[NFit-1][0] << "\t"
	     << result[NFit-1][1] << "\t"
	     << result[NFit-1][2] << "\t"
	     << std::endl;

    getchar();
  }
  gSystem->ProcessEvents();

}
