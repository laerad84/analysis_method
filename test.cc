#include "TCanvas.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TMath.h"
#include "TArc.h"
#include "TVirtualFitter.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TPCData.h"
#include "ConverterFunction.h"

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
int main( int argc, char** argv ){
//void ReadConverted_findCircle_1(){
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
    grCluster[ifit]->SetMarkerColor(ifit+1);
    grCluster[ifit]->SetMarkerStyle(ifit+20);
  }
  TH1D* hisF = new TH1D("hisF","hisF",1000,0,10000);

  std::cout<< "Starting Loop" << std::endl;
  for( int ievt = 0; ievt < tr->GetEntries(); ievt++){
    if( ievt > 40 ){ break; }
    tr->GetEntry(ievt);
    for( int i = 0; i< 10; i++){
      grHit[i]->Set(0);
      grCluster[i]->Set(0);
    }


    /// Hit version ///
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
    std::vector<TPCPadHit> hitCopy = hitArr;
    std::cout << hitArr.size() << std::endl;
    hitArr.erase( hitArr.begin() );
    Int_t clusterNumbers = 0;
    TPCPadHitCluster clusterList[300];
    Int_t nloop=0;
    std::vector<TPCPadHitCluster> clvec;
    HitClusteringB( hitArr, clvec );
    for( int i = 0; i< clvec.size(); i++){
      std::cout<< i << "\t" << clvec.at(i).HitArr().size() << std::endl;
    }

    getchar();
  }

  TCanvas* can1 = new TCanvas("can1","can1",800,800);
  hisF->Draw();
}
