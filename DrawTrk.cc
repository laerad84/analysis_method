/// Read track Data and Helix fitting
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "Data.h"
#include "CalculationFunction.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TVirtualFitter.h"
#include <vector>
#include <list>
const double pimass = 139.57018;//MeV
const double pmass  = 938.272046;//MeV
const double kpmass = 493.667;//MeV
TGraph2D* grTrackSample;

void HelixCircle( Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t ){
  Int_t np = grTrackSample->GetN();
  f = 0;
  Double_t *x = grTrackSample->GetX();
  Double_t *y = grTrackSample->GetY();
  Double_t *z = grTrackSample->GetZ();
  Int_t nHit = 0;
  for( Int_t i = 0; i< np; i++){
    //par[0] : x center
    //par[1] : y Yoffset at theta == 0;
    //par[2] : z center
    //par[3] : r
    //par[4] : dY/dTheta

    Double_t u = x[i] - par[0];
    Double_t t = z[i] - par[1];
    Double_t dr = par[2] - TMath::Sqrt( u*u + t*t );
    if( dr > 20 ){ continue; }
    f += dr*dr;
    nHit++;
  }
}

int main( int argc, char** argv){
//void DrawTrk(){
  /*
  TString mSystemName = gSystem->GetFromPipe("uname");
  std::cout<< mSystemName << std::endl;
  if( mSystemName == "Darwin"){
    gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method.dylib");
    gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  }else if( mSystemName == "Linux"){
    gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method.so");
    gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  }
  */
  TFile* itf = new TFile("Proton_trk.root");
  TTree* itr = (TTree*)itf->Get("convTree");

  TClonesArray* TrkArr = new TClonesArray("TGraph2D");
  itr->SetBranchAddress("TrkArr",&TrkArr);

  TFile* otf = new TFile("TestDist.root","recreate");
  TTree* otr = new TTree("Output","");
  Int_t    nComb;
  Int_t    nTrk;
  Double_t dY[400];
  Double_t dR[400];
  Int_t    nMerged;
  Int_t    nTrackBefore;
  Int_t    nTrackAfter;
  TClonesArray* MergedTrkArr = new TClonesArray("TGraph2D");
  TClonesArray* MergedSprialArr= new TClonesArray("HSprial");
  otr->Branch("nComb",&nComb,"nComb/I");
  otr->Branch("nTrk",&nTrk,"nTrk/I");
  otr->Branch("dY",dY,"dY[nComb]/D");
  otr->Branch("dR",dR,"dR[nComb]/D");
  otr->Branch("nMerged");
  otr->Branch("MergedTrkArr",&MergedTrkArr,256000,99);
  otr->Branch("MergedSprialArr",&MergedSprialArr,256000,99);

  otr->Branch("nTrackBefore",&nTrackBefore,"nTrackBefore/I");
  otr->Branch("nTrackAfter",&nTrackAfter,"nTrackAfter/I");

  TGraph* grX[20];
  for( int i = 0; i< 20; i++ ){
    grX[i] = new TGraph();
    grX[i]->SetMarkerStyle(20+i%4);
    grX[i]->SetMarkerColor(1+i%5);
  }

  TCanvas* can = new TCanvas("can","can",800,800);
  can->Divide(2,2);
  can->cd(1);
  gPad->DrawFrame(-300,-300,300,300);
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter *circleFitter = TVirtualFitter::Fitter(0,3);
  circleFitter->SetFCN(HelixCircle);

  for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
  //for( int ievt = 0; ievt < 100; ievt++){
    itr->GetEntry(ievt);
    std::cout<< ievt << std::endl;
    for( int i = 0; i< 20; i++){
      grX[i]->Set(0);
    }
    MergedTrkArr->Clear();
    MergedSprialArr->Clear();

    if( TrkArr->GetEntries() > 20 ){ continue; }
    nTrk = TrkArr->GetEntries();
    TClonesArray* sprialArr = new TClonesArray("HSprial");
    TClonesArray* TrackArr = new TClonesArray("TGraph2D");
    std::cout<< "trk" << std::endl;
    Int_t nSizeTrack = 0;
    for( int itrk = 0; itrk < TrkArr->GetEntries(); itrk++){
      std::cout<<"itrk "<<  itrk << std::endl;
      TGraph2D* itrack = (TGraph2D*)TrkArr->At(itrk);
      grTrackSample = itrack;
      if( itrack->GetN() < 5){
	continue;
      }

      for( int ip = 0; ip < itrack->GetN(); ip++){
	grX[itrk]->SetPoint( ip, itrack->GetZ()[ip],itrack->GetX()[ip]);
	//grTrackSample->SetPoint( ip, itrack->GetZ()[ip],itrack->GetX()[ip]);
      }
      std::cout << "Index " << std::endl;
      Int_t Index[3]={0,(int)((itrack->GetN()-1)/2),itrack->GetN() -1};
      TVector3 vec[3];

      for( int i = 0; i< 3; i++){
	vec[i].SetXYZ(itrack->GetX()[Index[i]],
		       itrack->GetY()[Index[i]],
		       itrack->GetZ()[Index[i]]);
      }
      std::cout<< "Sprial" << std::endl;
      HSprial sprial = GenerateSprial(vec[0],vec[1],vec[2]);
      sprial.ID = itrk;
      circleFitter->SetParameter(0,"X", sprial.X, 10, sprial.X-50, sprial.X+50 );
      circleFitter->SetParameter(1,"Z", sprial.Z, 10, sprial.Z-50, sprial.Z+50);
      circleFitter->SetParameter(2,"R", sprial.R, 10, sprial.R-50, sprial.R+50);
      Double_t arglist[1] = {0};
      circleFitter->ExecuteCommand("MIGRAD",arglist, 0 );
      sprial.X = circleFitter->GetParameter(0);
      sprial.Z = circleFitter->GetParameter(1);
      sprial.R = circleFitter->GetParameter(2);
      sprial.EX = circleFitter->GetParError(0);
      sprial.EZ = circleFitter->GetParError(1);
      sprial.ER = circleFitter->GetParError(2);
      sprial.InitPos = TVector3(itrack->GetX()[0],itrack->GetY()[0],itrack->GetZ()[0]);
      sprial.FinalPos= TVector3(itrack->GetX()[Index[2]],itrack->GetY()[Index[2]],itrack->GetZ()[Index[2]]);
      new((*sprialArr)[nSizeTrack]) HSprial(sprial);
      new((*TrackArr)[nSizeTrack]) TGraph2D(*itrack);
      nSizeTrack++;
    }
    std::cout<< "Combine Tracks" << std::endl;
    std::cout<< ievt << "\t" << TrkArr->GetEntries() <<"\t" << TrackArr->GetEntries() << std::endl;
    nComb = 0;


    std::list<int> tmpIDList;
    std::list<int> tmpIDList_bck;
    std::vector<std::vector<int> > clusterArr;
    for( int itrk = 0; itrk < sprialArr->GetEntries(); itrk++){
      std::cout<< itrk << std::endl;
      TGraph2D* itrack = (TGraph2D*)TrackArr->At(itrk);
      if( itrack->GetN() < 5 ){ continue; }
      tmpIDList.push_back( itrk );
      tmpIDList_bck.push_back( itrk );
    }
    std::cout<< "tmpIDList.size() " << tmpIDList.size() << std::endl;
    for( std::list<int>::iterator i = tmpIDList.begin();
	 i != tmpIDList.end();
	 i++){
      std::cout << *i << "\t";
    }
    std::cout<< std::endl;

    while( tmpIDList.size() > 0){
      std::vector<int> indexList;
      //indexList.clear();
      int firstIndex = tmpIDList.front();
      indexList.push_back(firstIndex);
      tmpIDList.pop_front();
      for(int i = 0 ; i < indexList.size(); i++){
	std::cout<< indexList.at(i) << std::endl;
	TGraph2D* itrack = (TGraph2D*)TrackArr->At(indexList.at(i));
	HSprial*  isprial= (HSprial*)sprialArr->At(indexList.at(i));
	for( std::list<int>::iterator j = tmpIDList.begin();
	     j != tmpIDList.end();
	     j++){

	  TGraph2D* jtrack = (TGraph2D*)TrackArr->At(*j);
	  HSprial* jsprial = (HSprial*)sprialArr->At(*j);
	  Double_t r0,y0;
	  Double_t r1,y1;
	  isprial->CalculateDist(jtrack,r0,y0);
	  jsprial->CalculateDist(itrack,r1,y1);
	  if( (r0 < 10 && y0 < 20) ||
	      (r1 < 10 && y1 < 20) ){
	    std::cout<< "ADD" << std::endl;
	    std::cout<< indexList.size() << "\t" <<indexList.at(i) << "\t" <<  *j << std::endl;
	    Int_t k = (int)(*j);
	    std::cout<< k << std::endl;
	    indexList.push_back(k);
	    std::cout<< indexList.size() << std::endl;
	    std::list<int>::iterator dj = j--;
	    tmpIDList.erase(dj);
	  }
	}
      }
      clusterArr.push_back(indexList);
    }
    nTrackBefore = tmpIDList_bck.size();
    nTrackAfter  = clusterArr.size();
    std::cout<< "//////////////////////////////////////////////////" << std::endl;
    std::cout<< "Cluster Size : " << clusterArr.size() << std::endl;
    std::cout<< "//////////////////////////////////////////////////" << std::endl;

    for( std::vector<std::vector<int> >::iterator i = clusterArr.begin();
	 i  != clusterArr.end();
	 i++ ){
      std::cout<< "Cluster" << std::endl;
      for( std::vector<int>::iterator j = (*i).begin();
	   j != (*i).end();
	   j++){
	std::cout<< *j << "\t";
      }
      std::cout<< std::endl;
    }
    //getchar();
    otr->Fill();
  }
  otr->Write();
  otf->Close();
}
