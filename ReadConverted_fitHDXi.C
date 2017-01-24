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

const double pimass = 139.57018;//MeV
const double pmass  = 938.272046;//MeV
const double kpmass = 493.667;//MeV
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

void ReadConverted_fitHDXi(){
  std::cout <<"ReadConverted" << std::endl;
  gSystem->Load("lib/libanalysis_method.dylib");
  TFile* itf  = new TFile("HDXiE42_conv.root");
  TTree* itr = (TTree*)itf->Get("convTree");
  TClonesArray* HitArr = new TClonesArray("TPCPadHit");
  GsimGenParticleData*   particle = new GsimGenParticleData();
  GsimDetectorEventData* tpcData  = new GsimDetectorEventData();
  GsimDetectorEventData* ftofData = new GsimDetectorEventData();
  GsimDetectorEventData* DC1Data  = new GsimDetectorEventData();
  GsimDetectorEventData* DC2Data  = new GsimDetectorEventData();
  GsimDetectorEventData* DC3Data  = new GsimDetectorEventData();
  GsimDetectorEventData* NBARData = new GsimDetectorEventData();
  itr->SetBranchAddress("HitArr",&HitArr);
  itr->SetBranchAddress("GenParticle.",&particle);
  //itr->SetBranchAddress("TPC.",&tpcData);
  itr->SetBranchAddress("FTOF.",&ftofData);
  itr->SetBranchAddress("DC1.",&DC1Data);
  itr->SetBranchAddress("DC2.",&DC2Data);
  itr->SetBranchAddress("DC3.",&DC3Data);
  itr->SetBranchAddress("NBAR.",&NBARData);

  TFile* otf = new TFile("circleHDXi.root","recreate");
  TTree* otr = new TTree("analysisTree","circle fitting result");

  Int_t iEvent;
  Int_t iMCEvent;
  Int_t ABORT;
  Int_t NCl;
  Int_t NSprial;
  Int_t NSprialM;
  Int_t NSprialP;

  TClonesArray* SprialArr = new TClonesArray("HSprial");
  TClonesArray* SprialMArr = new TClonesArray("HSprial");
  TClonesArray* SprialPArr = new TClonesArray("HSprial");

  TClonesArray* ClusterArr= new TClonesArray("TPCCluster");
  TClonesArray* Particle = new TClonesArray("TLorentzVector");
  TClonesArray* vertex   = new TClonesArray("TVector3");
  Int_t NParticle=0;
  Int_t Index[30];
  otr->Branch("iEvent",&iEvent,"iEvent/I");
  otr->Branch("iMCEvent",&iMCEvent,"iMCEvent/I");
  otr->Branch("ABORT",&ABORT,"ABORT/I");

  otr->Branch("SprialArr",&SprialArr,256000,99);
  otr->Branch("SprialMArr",&SprialMArr,256000,99);
  otr->Branch("SprialPArr",&SprialPArr,256000,99);
  otr->Branch("ClusterArr",&ClusterArr,256000,99);
  otr->Branch("NCl",&NCl,"NCl/I");
  otr->Branch("NSprial",&NSprial,"NSprial/I");
  otr->Branch("NSprialM",&NSprialM,"NSprialM/I");
  otr->Branch("NSprialP",&NSprialP,"NSprialP/I");
  otr->Branch("NParticle",&NParticle,"NParticle/I");
  otr->Branch("Index",Index,"Index[NParticle]/I");
  otr->Branch("Particle",&Particle, 256000,99);
  otr->Branch("Vertex",&vertex,256000,99);
  const Int_t  NFit = 30;
  TGraph* grHit[NFit];
  TGraph* grCluster[NFit];
  std::vector<TVector3> HitPointArr[NFit];
  for( int ifit = 0; ifit < NFit; ifit++){
    grHit[ifit] = new TGraph();
    grHit[ifit]->SetMarkerColor(ifit%5+1);
    grHit[ifit]->SetMarkerStyle(3);
    grCluster[ifit] = new TGraph();
    grCluster[ifit]->SetMarkerColor(ifit%5+1);
    grCluster[ifit]->SetMarkerStyle(ifit%6+20);
  }

  std::cout<< "Starting Loop" << std::endl;
  iEvent = 0;
  iMCEvent = 0;
  //for( int ievt = 0; ievt < 100; ievt++){
  for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
    //if( ievt > 100 ){ break; }
    itr->GetEntry(ievt);
    if( ievt % 100 == 0){ std::cout<< ievt << "/" << itr->GetEntries() << std::endl;}
    /// initialize
    iMCEvent = ievt;
    for( int i = 0; i< NFit; i++){
      grHit[i]->Set(0);
      grCluster[i]->Set(0);
      HitPointArr[i].clear();
    }
    ClusterArr->Clear();
    SprialArr->Clear();
    Particle->Clear();
    vertex->Clear();
    NParticle = 0;
    NCl = 0;
    NSprial = 0;
    NSprialM = 0;
    NSprialP = 0;
    /// Hit version ///
    std::vector<int> MCTrackIDList;
    for( int ientry  = 0; ientry < HitArr->GetEntries(); ientry++){
      TPCPadHit* hit = (TPCPadHit*)HitArr->At(ientry);
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

    /*
      std::cout<< "MCTrackIDList : " << MCTrackIDList.size() <<  std::endl;
      for( int i = 0; i< MCTrackIDList.size(); i++){
      std::cout<< i << "\t" << MCTrackIDList.at(i) << std::endl;
      }
    */

    /// Get Initial information ///
    Int_t minimumRow=99;
    Int_t maximumRow=-1;
    //std::cout<< "Copy Hit Array " << std::endl;

    std::vector<TPCPadHit> hitArr;
    for( int ientry = 0; ientry < HitArr->GetEntries(); ientry++){
      TPCPadHit* hit = (TPCPadHit*)HitArr->At(ientry);
      hitArr.push_back(*hit);
      if( hit->Row() < 5 ){ continue; }
      if( hit->Row() < minimumRow ){ minimumRow = hit->Row(); }
      if( hit->Row() > maximumRow ){ maximumRow = hit->Row(); }
    }
    //std::cout<< "End Copy " << std::endl;
    if( hitArr.size() < 5 ){ continue; }
    std::cout<< "Clustering"  << std::endl;
    std::vector<TPCCluster> clvec;
    HitClustering( hitArr, clvec );
    if( clvec.size() > 200 ){ continue; }
    /// SetClusterData
    for( int i = 0; i< clvec.size(); i++){
      new((*ClusterArr)[i]) TPCCluster(clvec.at(i));
    }

    /// Blocking Test ///
    std::cout<< clvec.size() << std::endl;
    if( clvec.size() > 200 ){ continue; }
    NCl = clvec.size();

    std::vector<TPCCluster> clVec0 = clvec;
    //std::vector<TPCCluster> block = TPCClusterBlocker(clvec,0);
    //std::cout<< block.size() << std::endl;
    std::vector<TPCCluster> Block[50];
    std::vector<TPCCluster> TrackCand[50];

    HCircle circle[30];
    HSprial sprial[30];
    TArc*    arc[30];

    Int_t nGraph = 0;
    Int_t nBlock = 0;
    Int_t nCand = 0;
    Int_t nCircle = 0;
    for( int n = 0; n< NFit; n++){
      if( clvec.size() == 0 ){ break; }
      Block[nBlock] = TPCClusterBlocker(clvec,0);
      nBlock++;
      if( nBlock == 10 ){ break; }
    }
    std::cout<< "Block : " << nBlock << std::endl;
    for( int n = 0; n < nBlock;n++){
      std::cout<< n << std::endl;
      while( Block[n].size() > 3 ){
	Int_t CandSize = BlockDivider( Block[n], TrackCand[nCand],0);
	//std::cout<< "Cand : " << CandSize << std::endl;
	if( TrackCand[nCand].size() == 0 ){
	  //std::cout<< "InnerCand" << std::endl;
	  CandSize = BlockDivider( Block[n], TrackCand[nCand],1);
	  //std::cout<< "InnerCand : " << CandSize << std::endl;
	}
	//std::cout <<Block[n].size() << "\t" <<  TrackCand[nCand].size() << std::endl;
	if( CandSize == 0 ){
	  break;
	}else{
	  nCand++;
	}
      }

      if( nCand >= 30 ){ break; }
    }
    std::cout<< nCand << std::endl;
    if( nCand >= 30 ){ continue;}
    //std::cout<< nBlock << "\t" << nCand << std::endl;

    std::cout<< "Fit Circle : " << nCand << std::endl;
    for( int n = 0; n < nCand; n++){
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
	new((*SprialArr)[nCircle]) HSprial(sprial[nCircle]);
	sprial[nCircle].Print();
	nCircle++;
      }
    }
    NSprial = nCircle;
    //if( nCircle != 5 ){ continue; }


    std::vector<HSprial> plusTrack;
    std::vector<HSprial> minusTrack;
    for( int i = 0; i< nCircle; i++){
      if( sprial[i].RL > 0 ){
	new((*SprialPArr)[NSprialP]) HSprial(sprial[i]);
	plusTrack.push_back( sprial[i]);
	NSprialP++;
      }else{
	new((*SprialMArr)[NSprialM]) HSprial(sprial[i]);
	minusTrack.push_back(sprial[i]);
	NSprialM++;
      }
    }

    /// Find minimum
    TVector3 vertexPoint[2];
    TVector3 Momentum[2][2];
    Double_t DistSquare = 1e9;
    Double_t DistMin[2] = {1e9,1e9};
    Int_t    TrackIndex[2]={-1,-1};
    Bool_t   bFound[3][2]={{false}};
    Double_t Dist[3][2]={{1e9}};
    TVector3 pos[3][2];
    Int_t    index[3][2]={{-1}};
    Int_t    KPIndex = -1;
    if( NSprialP == 3 && NSprialM == 2 ){
      std::cout<< "Sprial : " << plusTrack.size() << "\t" << minusTrack.size() << std::endl;
      for( int ip = 0; ip < plusTrack.size();ip++){
	for( int im = 0; im < minusTrack.size(); im++){
	  //Double_t dist;
	  //TVector3 position(0,0,0);
	  //bFound[ip][im] = CalculateCrossing(plusTrack.at(ip),minusTrack.at(im),dist, position);
	  //TVector3 position[2];
	  //CalculateCircleCrossing( plusTrack.at(ip),minusTrack.at(im),position[0],position[1]);
	  //position[0].Print();
	  //position[1].Print();
	  //getchar();
	  //bFound[ip][im] = CalculateCrossing(plusTrack.at(ip),minusTrack.at(im),Dist[ip][im],position);
	  bFound[ip][im] = CalculateCrossing(plusTrack.at(ip),minusTrack.at(im),Dist[ip][im],pos[ip][im]);
	  std::cout<< ip << "\t" << im << "\t" << bFound[ip][im] << std::endl;
	  pos[ip][im].Print();
	  index[ip][im]= ip*minusTrack.size() + im;
	}
      }
      std::cout<< "Calculate "<< std::endl;
      for( int ip = 0; ip < plusTrack.size(); ip++){
	for( int im = 0; im < minusTrack.size(); im++){
	  if( Dist[ip][im] < DistMin[im]){
	    DistMin[im]    = Dist[ip][im];
	    TrackIndex[im] = ip;
	  }
	}
      }
      std::cout<< "Momentum " << std::endl;
      if( TrackIndex[0] != TrackIndex[1] && TrackIndex[0] != -1 && TrackIndex[1] != -1){
	/// Calculate momentum ///
	if( TrackIndex[0] == 0 ){
	  if( TrackIndex[1] == 1 ){
	    KPIndex = 2;
	  }else{
	    KPIndex = 1;
	  }
	}else if( TrackIndex[0] ==1 ){
	  if( TrackIndex[1] == 2 ){
	    KPIndex = 0;
	  }else{
	    KPIndex = 2;
	  }
	}else{
	  if( TrackIndex[1] == 0 ){
	    KPIndex = 1;
	  }else{
	    KPIndex = 0;
	  }
	}



	TVector3 tan[2][2];
	vertexPoint[0] = pos[TrackIndex[0]][0];
	vertexPoint[1] = pos[TrackIndex[1]][1];
	vertexPoint[0].Print();
	vertexPoint[1].Print();
	new((*vertex)[0]) TVector3(vertexPoint[0]);
	new((*vertex)[1]) TVector3(vertexPoint[1]);

	std::cout << "momentum" << std::endl;
	TLorentzVector s[5];
	Double_t energy;
	CalculateCircleTangent( plusTrack.at(TrackIndex[0]),
				pos[TrackIndex[0]][0],tan[0][0]);
	Momentum[0][0] = tan[0][0]*CalculateMomentum( plusTrack.at(TrackIndex[0]),1.);
	energy = TMath::Sqrt( Momentum[0][0].Mag2() + pmass*pmass );
	s[0].SetPxPyPzE(Momentum[0][0].X(),
			Momentum[0][0].Y(),
			Momentum[0][0].Z(),
			energy);

	std::cout<< "Check" << std::endl;
	CalculateCircleTangent( minusTrack.at(0),
				pos[TrackIndex[0]][0],tan[0][1]);
	Momentum[0][1] = tan[0][1]*CalculateMomentum( minusTrack.at(0),1.);
	energy = TMath::Sqrt( Momentum[0][1].Mag2() + pimass*pimass );
	s[1].SetPxPyPzE(Momentum[0][1].X(),
			Momentum[0][1].Y(),
			Momentum[0][1].Z(),
			energy);

	std::cout<< "Check" << std::endl;
	CalculateCircleTangent( plusTrack.at(TrackIndex[1]),
				pos[TrackIndex[1]][1],tan[1][0]);
	Momentum[1][0] = tan[1][0]*CalculateMomentum( plusTrack.at(TrackIndex[1]),1.);
	energy = TMath::Sqrt( Momentum[1][0].Mag2() + pmass*pmass );
	s[2].SetPxPyPzE(Momentum[1][0].X(),
			Momentum[1][0].Y(),
			Momentum[1][0].Z(),
			energy);

	std::cout<< "Check" << std::endl;
	CalculateCircleTangent( minusTrack.at(1),
				pos[TrackIndex[1]][1],tan[1][1]);
	Momentum[1][1] = tan[1][1]*CalculateMomentum( minusTrack.at(1),1.);
	energy = TMath::Sqrt( Momentum[1][1].Mag2() + pimass*pimass );
	s[3].SetPxPyPzE(Momentum[1][1].X(),
			Momentum[1][1].Y(),
			Momentum[1][1].Z(),
			energy);
	TVector3 kpMom;
	CalculateCircleTangent( plusTrack.at(KPIndex),
				TVector3(0,0,-143),kpMom);
	energy = TMath::Sqrt( kpMom.Mag2() + kpmass*kpmass);
	s[4].SetPxPyPzE(kpMom.X(),
		       kpMom.Y(),
		       kpMom.Z(),
		       energy);
	for( int il = 0; il < 5; il++){
	  new((*Particle)[il]) TLorentzVector(s[il]);
	}
	NParticle = 5;

	std::cout<< "Check" << std::endl;
      }
    }
    std::cout<< "Circle : " << plusTrack.size() << "\t" << minusTrack.size() << std::endl;
    std::cout<< "Fill" << std::endl;
    otr->Fill();
    std::cout<< "End Fill" << std::endl;
  }
  otr->Write();
  otf->Close();
}
