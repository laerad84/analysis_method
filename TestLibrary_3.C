#include "TVirtualFitter.h"

const double pimass = 139.57018;//MeV
const double pmass  = 938.272046;//MeV
const double kpmass = 493.667;//MeV
TGraph2D* grTrackSample;
void HelixY( Int_t &, Double_t *, Double_t &f, Double_t *par , Int_t ){
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
    Double_t t = z[i] - par[2];
    Double_t theta = TMath::ATan2( u,t);
    Double_t expY = theta*par[4] + par[1];
    Double_t v = y[i] - expY;
    Double_t dr = par[3] - TMath::Sqrt( u*u + t*t );
    //if( dr > 20 ){ continue; }
    f += v*v;
    nHit++;
  }
}


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
    //if( dr > 20 ){ continue; }
    f += dr*dr;
    nHit++;
  }
}


void TestLibrary_3(){/// Save cluster data
  TString mSystemName = gSystem->GetFromPipe("uname");
  std::cout<< mSystemName << std::endl;
  if( mSystemName == "Darwin"){
    gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method.dylib");
  }else if( mSystemName == "Linux"){
    gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method.so");
  }

  //GenFitter* fitter = new GenFitter("Data/TPC.root");
  //TVirtualFitter::SetDefaultFitter("Minuit");
  grTrackSample = new TGraph2D();
  //fitter->SetField("~/local/hep/E42/E42/data/field/SC_FieldMap.root",TVector3(0,0,0),TVector3(0,0,0));
  //TFile* itf = new TFile("HDE42_conv.root");
  //TTree* itr = (TTree*)itf->Get("convTree");
  //TFile* itf = new TFile("TestHit.root");
  //TTree* itr = (TTree*)itf->Get("padHit");

  TChain* itr = new TChain("convTree");
  //itr->Add("MCData/HDE42_conv_*.root");
  //itr->Add("HD_0_conv.root");
  itr->Add("lambda_conv.root");
  TClonesArray* HitArr = new TClonesArray("TPCPadHit");
  GsimGenParticleData*   particle = new GsimGenParticleData();

  grTrackSample = new TGraph2D();

  /*
  GsimDetectorEventData* tpcData  = new GsimDetectorEventData();
  GsimDetectorEventData* ftofData = new GsimDetectorEventData();
  GsimDetectorEventData* DC1Data  = new GsimDetectorEventData();
  GsimDetectorEventData* DC2Data  = new GsimDetectorEventData();
  GsimDetectorEventData* DC3Data  = new GsimDetectorEventData();
  GsimDetectorEventData* NBARData = new GsimDetectorEventData();
  */
  itr->SetBranchAddress("HitArr",&HitArr);
  itr->SetBranchAddress("GenParticle.",&particle);
  //itr->SetBranchAddress("TPC.",&tpcData);
  /*
  itr->SetBranchAddress("FTOF.",&ftofData);
  itr->SetBranchAddress("DC1.",&DC1Data);
  itr->SetBranchAddress("DC2.",&DC2Data);
  itr->SetBranchAddress("DC3.",&DC3Data);
  itr->SetBranchAddress("NBAR.",&NBARData);
  */
  TFile* otf = new TFile("TestFit.root","recreate");
  TTree* otr = new TTree("anaTree","Test");
  TClonesArray* pSprialArr = new TClonesArray("HSprial");
  TClonesArray* mSprialArr = new TClonesArray("HSprial");
  TClonesArray* pMomentumArr = new TClonesArray("TLorentzVector");
  TClonesArray* mMomentumArr = new TClonesArray("TLorentzVector");
  TClonesArray* lMomentumArr = new TClonesArray("TLorentzVector");
  TClonesArray* PositionArr  = new TClonesArray("TVector3");
  TClonesArray* ClArr = new TClonesArray("TPCCluster");
  Int_t         EventNumber;
  Int_t         nPTrack;
  Int_t         nMTrack;
  Int_t         nLambda;
  Double_t      Distance[100];
  otr->Branch("ClArr",&ClArr,256000,99);
  otr->Branch("EventNumber",&EventNumber,"EventNumber/I");
  TH1D* hislMass = new TH1D("hislMass","hislMass",3000,500,3500);
  TH2D* hisClusterSizeEne  = new TH2D("hisClustersizeEne","hisCLustersizeEne",30,0,30,100,0,1000);
  TH1D* hisDist = new TH1D("hisDist","hisDist",400,0,400);
  TH2D* hisCross= new TH2D("hisCross","hisCross",300,-300,300,300,-300,300);
  TH3D* hisBox = new TH3D("hisBox","hisBox;X;Y;Z",30,-300,300,30,-300,300,30,-300,300);
  TH2D* hisMag = new TH2D("hisMag","hisMag",100,0,2000,100,0.5,1.5);
  TCanvas* can = new TCanvas("can","can",1200,1200);
  //TCanvas* can = new TCanvas("can","can",1200,600);
  //can->Divide(2,1);
  //can->cd(1);
  //gPad->DrawFrame(-300,-300,300,300);
  //can->cd(2);
  can->Divide(2,2);
  can->cd(1);
  gPad->DrawFrame(-300,-300,300,300);
  can->cd(2);
  hisBox->Draw();
  can->cd(3);
  hisBox->Draw();

  const Int_t maxNBlock = 50;
  const Int_t maxNTrackCand    = 50;
  TGraph* grCandidates[maxNBlock];
  for( int i = 0; i< maxNBlock; i++){
    grCandidates[i] = new TGraph();
    grCandidates[i]->SetMarkerStyle( 20+i%8);
    grCandidates[i]->SetMarkerColor( i/7 +1);
  }


  //for( int ievt = 0; ievt <itr->GetEntries() ;ievt++){
  for( int ievt = 0; ievt < 1000;ievt++){
    EventNumber = ievt;
    ClArr->Clear();
    itr->GetEntry( ievt );
    TPolyMarker3D* total = new TPolyMarker3D();
    total->SetMarkerStyle(20);
    TPolyMarker3D* track[maxNBlock];
    for( int i = 0; i< maxNBlock; i++){
      track[i] = new TPolyMarker3D();
      track[i]->SetMarkerStyle( 20+ i%8 );
      track[i]->SetMarkerColor( i/7 +1 );
    }
    std::cout<< "Clustering" << std::endl;
    for( int igr =0; igr < maxNBlock; igr++){
      grCandidates[igr]->Set(0);
    }

    Bool_t bCluster = HitClustering( HitArr, ClArr );
    if( !bCluster ){ continue; }
    if( ClArr->GetEntries() > 300 ){ continue; }
    otr->Fill();
  }
  can->cd(4);
  hislMass->Draw();
  otf->cd();
  otr->Write();
  otf->Close();
}
