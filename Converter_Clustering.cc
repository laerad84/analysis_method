#include "Data.h"
#include "TClonesArray.h"
#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "CalculationFunction.h"
#include "ConverterFunction.h"
#include "TSystem.h"
#include "TVirtualFitter.h"

int main( int argc, char** argv ){



//void Converter_Clustering(){
  /*
  gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  gSystem->Load("./lib/libanalysis_method");
  */
  std::string ifilename = argv[1];
  //std::string ifilename = "Lambda.root";
  //std::string ifilename = "Proton.root";
  //std::string ifilename = "HEXiE42.root";
  std::string ofilename = ifilename.substr(0,ifilename.rfind(".")) + "_trk.root";
  std::cout<< ifilename << std::endl;
  std::cout<< ofilename << std::endl;

  TFile* itf = new TFile(ifilename.c_str());
  TTree* itr = (TTree*)itf->Get("eventTree00");
  GsimGenParticleData*   particle = new GsimGenParticleData();
  GsimDetectorEventData* tpcData = new GsimDetectorEventData();
  GsimDetectorEventData* ftofData = new GsimDetectorEventData();
  GsimDetectorEventData* DC1Data  = new GsimDetectorEventData();
  GsimDetectorEventData* DC2Data  = new GsimDetectorEventData();
  GsimDetectorEventData* DC3Data  = new GsimDetectorEventData();
  GsimDetectorEventData* NBARData = new GsimDetectorEventData();
  itr->SetBranchAddress("GenParticle.",&particle);
  itr->SetBranchAddress("TPC.",&tpcData);
  itr->SetBranchAddress("FTOF.",&ftofData);
  itr->SetBranchAddress("DC1.",&DC1Data);
  itr->SetBranchAddress("DC2.",&DC2Data);
  itr->SetBranchAddress("DC3.",&DC3Data);
  itr->SetBranchAddress("NBAR.",&NBARData);

  TFile* otf = new TFile(ofilename.c_str(),"recreate");
  TTree* otr = new TTree("convTree","TPCTrack Data");
  Int_t EventID;
  TClonesArray* HitArr   = new TClonesArray("TPCPadHit");
  TClonesArray* ClArr    = new TClonesArray("TPCCluster");
  TClonesArray* TrkArr   = new TClonesArray("TGraph2D");
  Int_t nHit;
  Int_t nCl;
  Int_t nTrk;
  otr->Branch("EventID",&EventID,"EventID/I");
  otr->Branch("GenParticle.",&particle,16000,99);
  otr->Branch("nHit",&nHit,"nHit/I");
  otr->Branch("nCl",&nCl,"nCl/I");
  otr->Branch("nTrk",&nTrk,"nTrk/I");

  otr->Branch("HitArr",&HitArr,256000,99);
  otr->Branch("ClArr",&ClArr,256000,99);
  otr->Branch("TrkArr",&TrkArr,256000,99);

  otr->Branch("DC1.",&DC1Data,256000,99);
  otr->Branch("DC2.",&DC2Data,256000,99);
  otr->Branch("DC3.",&DC3Data,256000,99);
  otr->Branch("NBAR.",&NBARData,256000,99);
  otr->Branch("FTOF.",&ftofData,256000,99);

  //otr->Branch("clArr",&clArr,2560000,99);

  std::cout<< "Event Loop" << std::endl;
  for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
    //for( int ievt = 0; ievt < 100; ievt++){
    itr->GetEntry(ievt);
    std::cout<< ievt << std::endl;
    HitArr->Clear();
    ClArr->Clear();
    TrkArr->Clear();

    EventID = ievt;
    TClonesArray* trackArr = particle->briefTracks;
    TClonesArray* tpcHitArr= tpcData->hits;

    Int_t nTrack = trackArr->GetEntries();
    //Int_t nHit   = tpcHitArr->GetEntries();
    std::vector<TPCPadHit>        padHitArr = ConvertTPC(tpcData, particle);
    ///Save HitArr
    for( int i = 0; i< padHitArr.size(); i++){
      new((*HitArr)[i]) TPCPadHit(padHitArr.at(i));
    }
    Bool_t bCluster = HitClustering( HitArr, ClArr );
    nHit = HitArr->GetEntries();
    nCl  = ClArr->GetEntries();

    if( bCluster ){
      std::vector<Int_t> RootIDList = GetListOfTrackRoot( ClArr );
      for( Int_t arrIndex = 0; arrIndex < RootIDList.size(); arrIndex++){
	TClonesArray* subClArr = new TClonesArray("TPCCluster");
	Bool_t bResult = ClusterBlocker( ClArr, subClArr,RootIDList, RootIDList.at(arrIndex));
	if( !bResult){ continue; }
	TGraph2D* gTrack = new TGraph2D();
	gTrack->SetUniqueID(arrIndex);
	for( Int_t index = 0; index < subClArr->GetEntries(); index++){
	  TPCCluster* subCl = (TPCCluster*)subClArr->At(index);
	  gTrack->SetPoint( index, subCl->Position.X(), subCl->Position.Y(), subCl->Position.Z());
	}
	new((*TrkArr)[arrIndex]) TGraph2D(*gTrack);
      }
      nTrk = TrkArr->GetEntries();
    }else{
      /// some trouble in clustering ///
      ;
    }
    otr->Fill();
  }
  std::cout<< "Write" << std::endl;
  otr->Write();
  HitArr->Clear();
  std::cout<< "Close" << std::endl;
  otf->Close();
}
