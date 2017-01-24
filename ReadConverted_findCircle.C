#include "TCanvas.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TMath.h"
#include "TArc.h"
#include "TVirtualFitter.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include <vector>
TGraph* gr;
TGraphErrors* gre;
Double_t RCut = 200;

void myfcn( Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t ){
  Int_t np = gr->GetN();
  f = 0;
  Double_t *x = gr->GetX();
  Double_t *y = gr->GetY();
  for( Int_t i = 0; i < np; i++){
    Double_t u = x[i] - par[0];
    Double_t v = y[i] - par[1];
    Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
    if( TMath::Abs(dr) < RCut){
      f += dr*dr;
    }
  }
}
void myfcnWeight( Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t ){
  Int_t np = gre->GetN();
  f = 0;
  Double_t *x = gre->GetX();
  Double_t *y = gre->GetY();
  Double_t *w = gre->GetEX();

  for( Int_t i = 0; i< np; i++){
    Double_t u = x[i] - par[0];
    Double_t v = y[i] - par[1];
    Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
    if( TMath::Abs(dr) < RCut ){
      f += dr*dr*w[i];
    }
  }
}
void myfcnNHit( Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t ){
  Int_t np = gr->GetN();
  f = 0;
  Double_t *x = gr->GetX();
  Double_t *y = gr->GetY();
  for( int i = 0; i < np; i++){
    Double_t u = x[i] - par[0];
    Double_t v = y[i] - par[1];
    Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
    if( TMath::Abs(dr) < RCut ){
      f += -1;
    }
  }
}

void ReadConverted_findCircle(){
  std::cout <<"ReadConverted" << std::endl;
  gSystem->Load("lib/libanalysis_method.dylib");
  TFile* tf = new TFile("TestHit.root");
  TTree* tr = (TTree*)tf->Get("padHit");
  TClonesArray* arr = new TClonesArray("TPCPadHit");
  tr->SetBranchAddress("HitArr",&arr);

  TCanvas* can = new TCanvas("can","can",1800,1800);
  can->Divide(3,3);
  for( int i = 0; i< 9; i++){
    can->cd(i+1);
    gPad->DrawFrame(-300,-300,300,300);
    gPad->SetGrid();
  }
  can->Update();
  can->Modified();
  gSystem->ProcessEvents();
  std::cout<< "Prepared" << std::endl;
  const Int_t NFit = 10;
  TGraph* grHit[10];
  TGraph* grCluster[10];
  for( int ifit = 0; ifit < 10; ifit++){
    grHit[ifit] = new TGraph();
    grHit[ifit]->SetMarkerColor(1);
    grHit[ifit]->SetMarkerStyle(3);
    grCluster[ifit] = new TGraph();
    grCluster[ifit]->SetMarkerColor(2);
    grCluster[ifit]->SetMarkerStyle(4);
  }

  std::cout<< "Starting Loop" << std::endl;
  for( int ievt = 0; ievt < tr->GetEntries(); ievt++){
    if( ievt > 40 ){ break; }
    tr->GetEntry(ievt);
    for( int i = 0; i< 10; i++){
      grHit[i]->Set(0);
      grCluster[i]->Set(0);
    }


    /// Hit version ///
    TH1D* hisF = new TH1D("hisF","hisF",1000,0,10000);
    std::vector<int> MCTrackIDList;
    for( int ientry  = 0; ientry < arr->GetEntries(); ientry++){
      //std::cout<< ientry << "/"  << arr->GetEntries() << std::endl;
      TPCPadHit* hit = (TPCPadHit*)arr->At(ientry);
      Bool_t bNoSame = kTRUE;
      if( MCTrackIDList.size() == 0 ){
        MCTrackIDList.push_back(hit->MCTrackID());
      }else{
        for( int i = 0; i < MCTrackIDList.size(); i++){
	  if( MCTrackIDList.at(i) == hit->MCTrackID()){
	    bNoSame = kFALSE;
	    break;
	  }
        }
        if( bNoSame ){ MCTrackIDList.push_back(hit->MCTrackID());}
      }
    }

    std::cout<< "MCTrackIDList : " << MCTrackIDList.size() <<  std::endl;
    for( int i = 0; i< MCTrackIDList.size(); i++){
      std::cout<< i << "\t" << MCTrackIDList.at(i) << std::endl;
    }

    TVector3 InnerPos(0,0,0);
    TVector3 OuterPos(0,0,0);

    /// Get Initial information ///
    Int_t minimumRow=99;
    Int_t maximumRow=-1;
    std::vector<TPCPadHit> hitArr;
    for( int ientry = 0; ientry < arr->GetEntries(); ientry++){
      TPCPadHit* hit = (TPCPadHit*)arr->At(ientry);
      hitArr.push_back(*hit);
      //std::cout<< hit->Row() << std::endl;
      //hit->Position().Print();
      if( hit->Row() < 5 ){ continue; }
      if( hit->Row() < minimumRow ){ minimumRow = hit->Row(); }
      if( hit->Row() > maximumRow ){ maximumRow = hit->Row(); }
      //if( hit->MCTrackID() != MCTrackIDList.at(0)){ continue; }

      grHit[0]->SetPoint(grHit[0]->GetN(),hit->Position().Z(), hit->Position().X());
      //grHit[0]->SetPointError(grHit[0]->GetN()-1,hit->Energy(),0);
    }
    std::vector<TPCPadHitCluster> clvec = HitClustering( hitArr );
    for( int ientry = 0; ientry < clvec.size(); ientry++){
      grCluster[0]->SetPoint(grCluster[0]->GetN(), clvec.at(ientry).GetClusterPos().Z(),clvec.at(ientry).GetClusterPos().X());
    }

    std::cout<< clvec.size() << std::endl;
    Double_t pseudo_fit[NFit][3]={{0}};
    for( int itrk = 0; itrk < 1; itrk++){

      can->cd(itrk+1);
      grHit[itrk]->Draw("p");
      grCluster[itrk]->Draw("p");
      gr = grHit[itrk];
      //gr = grCluster[itrk];
      TVirtualFitter::SetDefaultFitter("Minuit");
      TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
      //fitter->SetFCN(myfcn);
      fitter->SetFCN(myfcnNHit);

      //Double_t fitPar[3]={pseudo_fit[itrk][0],pseudo_fit[itrk][1],300};
      Double_t fitPar[3]={-200,300,300};
      Double_t RCutArr[NFit]={0};
      Double_t RInit = 60;
      Int_t    Result[NFit] = {0};
      for( int i = 0; i< NFit; i++){
	RCutArr[i] = RInit*TMath::Power(0.8,i);
	if( RCutArr[i] < 30){ RCutArr[i] = 30;}
      }
      Double_t RMSArr[NFit]={40,30,20,10,10,5,5,5,5,5};
      Double_t result[NFit][3]={{0}};
      for( int i = 0; i< NFit; i++){
        std::cout<< "Fit" << std::endl;
        Double_t arglist[1] = {0};
	RCut = RCutArr[i];
	fitter->SetParameter(0, "x0",   fitPar[0], RMSArr[i], -2000,2000);
	fitter->SetParameter(1, "y0",   fitPar[1], RMSArr[i], -2000,2000);
	fitter->SetParameter(2, "R",    fitPar[2], RMSArr[i],200,4000);
	fitter->ExecuteCommand("MIGRAD", arglist, 0);
	std::cout<< "FitEnd" << std::endl;
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
      can->Update();
      can->Modified();
      gSystem->ProcessEvents();
      getchar();
      for( int i = 0; i< NFit; i++){
	arc[i]->Delete();
      }
    }
  }
  TCanvas* can1 = new TCanvas("can1","can1",800,800);
  hisF->Draw();
}
