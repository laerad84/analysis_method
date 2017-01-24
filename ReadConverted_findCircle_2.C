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

void ReadConverted_findCircle_2(){
  std::cout <<"ReadConverted" << std::endl;
  gSystem->Load("lib/libanalysis_method.dylib");
  TFile* tf = new TFile("TestHit.root");
  TTree* tr = (TTree*)tf->Get("padHit");
  TClonesArray* arr = new TClonesArray("TPCPadHit");
  tr->SetBranchAddress("HitArr",&arr);


  TCanvas* can = new TCanvas("can","can",1100,1100);
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
  const Int_t NFit = 30;
  TGraph* grHit[30];
  TGraph* grCluster[30];
  std::vector<TVector3> HitPointArr[30];
  for( int ifit = 0; ifit < 30; ifit++){
    grHit[ifit] = new TGraph();
    grHit[ifit]->SetMarkerColor(ifit%5+1);
    grHit[ifit]->SetMarkerStyle(3);
    grCluster[ifit] = new TGraph();
    grCluster[ifit]->SetMarkerColor(ifit%5+1);
    grCluster[ifit]->SetMarkerStyle(ifit%6+20);
  }
  TH1D* hisF = new TH1D("hisF","hisF",1000,0,10000);

  std::cout<< "Starting Loop" << std::endl;
  for( int ievt = 0; ievt < tr->GetEntries(); ievt++){
    if( ievt > 100 ){ break; }
    tr->GetEntry(ievt);
    for( int i = 0; i< 30; i++){
      grHit[i]->Set(0);
      grCluster[i]->Set(0);
      HitPointArr[i].clear();
    }


    /// Hit version ///
    std::vector<int> MCTrackIDList;
    for( int ientry  = 0; ientry < arr->GetEntries(); ientry++){
      std::cout<< ientry << "/"  << arr->GetEntries() << std::endl;
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
    std::cout<< "Copy Hit Array " << std::endl;
    std::vector<TPCPadHit> hitArr;
    for( int ientry = 0; ientry < arr->GetEntries(); ientry++){
      TPCPadHit* hit = (TPCPadHit*)arr->At(ientry);
      hitArr.push_back(*hit);
      if( hit->Row() < 5 ){ continue; }
      if( hit->Row() < minimumRow ){ minimumRow = hit->Row(); }
      if( hit->Row() > maximumRow ){ maximumRow = hit->Row(); }
      //if( hit->MCTrackID() != MCTrackIDList.at(0)){ continue; }
      grHit[0]->SetPoint(grHit[0]->GetN(),hit->Position().Z(), hit->Position().X());
    }
    std::cout<< "End Copy " << std::endl;
    std::cout<< "Clustering"  << std::endl;
    if( hitArr.size() < 5 ){ continue; }
    std::vector<TPCCluster> clvec;
    std::cout<< "hitArrsize: " << hitArr.size() << std::endl;
    for( int id = 0; id < hitArr.size(); id++){
      std::cout<< id << "\t"
	       << hitArr.at(id).m_PadID << "\t"
	       << hitArr.at(id).m_Phi << std::endl;
    }

    HitClustering( hitArr, clvec );
    std::cout<< "clvecsize : " << clvec.size() << std::endl;

    for( int id = 0; id< clvec.size(); id++){
      /*
      std::cout<< id << "\t"
	       << (int)clvec.at(id).NHit << "\t"
	       << (int)clvec.at(id).ID   << "\t"
	       << (double)clvec.at(id).Energy   << "\t"
	       << (int)clvec.at(id).Row  << "\t"
	       << (double)clvec.at(id).Col  << "\t"
	       << (double)clvec.at(id).PhiMax << "\t"
	       << (double)clvec.at(id).PhiMin << "\t"
	       << clvec.at(id).Position.X() << "\t"
	       << clvec.at(id).Position.Y() << "\t"
	       << clvec.at(id).Position.Z() << std::endl;
      */
      grHit[1]->SetPoint(grHit[1]->GetN(),clvec.at(id).Position.Z(),clvec.at(id).Position.X());
    }


    /// Blocking Test ///
    std::cout<< clvec.size() << std::endl;
    if( clvec.size() > 200 ){ continue; }

    std::vector<TPCCluster> clVec0 = clvec;
    //std::vector<TPCCluster> block = TPCClusterBlocker(clvec,0);
    //std::cout<< block.size() << std::endl;

    std::vector<TPCCluster> Block[10];
    std::vector<TPCCluster> TrackCand[30];
    HCircle circle[30];
    HSprial sprial[30];
    TArc*    arc[30];

    Int_t nGraph = 0;
    Int_t nBlock = 0;
    Int_t nCand = 0;
    Int_t nCircle = 0;
    for( int n = 0; n< 9; n++){
      if( clvec.size() == 0 ){ break; }
      Block[nBlock] = TPCClusterBlocker(clvec,0);
      nBlock++;
      if( nBlock == 10 ){ break; }
    }

    for( int n = 0; n < nBlock;n++){
      while( Block[n].size() > 3 ){
	Int_t CandSize = BlockDivider( Block[n], TrackCand[nCand],0);
	if( TrackCand[nCand].size() == 0 ){
	  CandSize = BlockDivider( Block[n], TrackCand[nCand],1);
	}
	std::cout <<Block[n].size() << "\t" <<  TrackCand[nCand].size() << std::endl;
	if( CandSize == 0 ){
	  break;
	}else{
	  nCand++;
	}
      }
      if( nCand >= 30 ){ break; }
    }

    std::cout<< nBlock << "\t" << nCand << std::endl;
    for( int n = 0; n < nCand; n++){
      for( int p = 0; p < TrackCand[n].size(); p++){
	grCluster[n]->SetPoint(grCluster[n]->GetN(),
			       TrackCand[n].at(p).Position.Z(),
			       TrackCand[n].at(p).Position.X());
      }
      if( TrackCand[n].size() >=5 ){
	Int_t firstIndex  = 0;
	Int_t middleIndex = (int)((TrackCand[n].size()-1)/2.);
	Int_t lastIndex   = TrackCand[n].size()-1;
	if( TrackCand[n].size() >= 6 ){
	  lastIndex  = TrackCand[n].size()-2;
	}
	if(TrackCand[n].size() >= 7){
	  firstIndex = 1;
	}
	if(TrackCand[n].size() >= 8 ){
	  lastIndex = TrackCand[n].size()-3;
	}
	sprial[nCircle] = GenerateSprial(TrackCand[n].at(firstIndex).Position,
					 TrackCand[n].at(middleIndex).Position,
					 TrackCand[n].at(lastIndex).Position);
	circle[nCircle] = ThreePointCircle( TrackCand[n].at(firstIndex).Position,
					    TrackCand[n].at(middleIndex).Position,
					    TrackCand[n].at(lastIndex).Position);
	circle[nCircle].ID = nCircle;
	sprial[nCircle].ID = nCircle;
	sprial[nCircle].Print();
	arc[nCircle] = new TArc(circle[nCircle].Z,circle[nCircle].X,circle[nCircle].R);
	nCircle++;
      }
    }
    ///Calculate Circles
    for( int i = 0; i< nCircle; i++){
      arc[i]->SetLineColor(i%6+1);
      arc[i]->SetFillStyle(0);
      arc[i]->SetLineWidth(2);
      arc[i]->SetLineStyle( i%5 +1 );
    }
    std::cout<< "Sprial" << std::endl;
    for( int i = 0; i< nCircle; i++){
      sprial[i].Print();
    }


    can->cd(1);
    grHit[0]->Draw("P");
    can->cd(3);
    for( int i = 0; i < nCand; i++){
      grCluster[i]->Draw("P");
    }
    for( int i = 0; i < nCircle; i++){
      arc[i]->Draw();
    }
    can->Update();
    can->Modified();
    gSystem->ProcessEvents();
    getchar();
    for( int i = 0; i< nCircle; i++){
      arc[i]->Delete();
    }

  }

  //TCanvas* can1 = new TCanvas("can1","can1",800,800);
  //hisF->Draw();
}
