
#include "ConstField.h"
#include "Exception.h"
#include "FieldManager.h"
#include "KalmanFitterRefTrack.h"
#include "StateOnPlane.h"
#include "Track.h"
#include "TrackPoint.h"
#include "MaterialEffects.h"
#include "RKTrackRep.h"
#include "EventDisplay.h"
#include "SpacepointMeasurement.h"
#include <vector>
#include <algorithm>

#include "TSystem.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2Poly.h"
#include "TMath.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TEveManager.h"
#include "TGeoManager.h"
#include "TGeoMaterialInterface.h"

#include "FieldManager.h"
#include "VariableField.h"

#include "GsimData/GsimGenParticleData.h"
#include "GsimData/GsimDetectorEventData.h"
#include "GsimData/GsimDetectorHitData.h"
#include "GsimData/GsimDigiData.h"
#include "GsimData/GsimTimeData.h"
#include "GsimData/GsimTrackData.h"

#include "DataReader.h"
#include "DataIO.h"
#include "DataIOAll.h"


//void TestRead(){
int main( int argc, char** argv){
  new TGeoManager("TPCGeom", "E42 GeoMetry" );
  gGeoManager->Import("TPC.root");
  TVector3 FieldCenter(0,0,-143);
  TVector3 FieldRotation(0,0,0);

  //genfit::VariableField* field =   new genfit::VariableField("/Users/jwlee/local/hep/E42/E42/data/field/SC_FieldMap.root",FieldCenter, FieldRotation);
  genfit::FieldManager::getInstance()->init( new genfit::VariableField("/Users/jwlee/local/hep/E42/E42/data/field/SC_FieldMap.root",FieldCenter, FieldRotation));
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

  TFile* itf = new TFile("TrackOutput.root");
  TTree* itr = (TTree*)itf->Get("TrackTree");

  Int_t EventID;
  Int_t TrackID;
  Int_t TrackIDMC;
  Int_t PID;
  Double_t Charge;
  TVector3* InitMom;
  TVector3* InitPos;
  Int_t nPoint;
  Double_t eArr[120];
  Double_t tArr[120];
  Double_t xArr[120];
  Double_t yArr[120];
  Double_t zArr[120];

  itr->SetBranchAddress("EventID",&EventID);
  itr->SetBranchAddress("TrackID",&TrackID);
  itr->SetBranchAddress("TrackIDMC",&TrackIDMC);
  itr->SetBranchAddress("PID",&PID);
  itr->SetBranchAddress("Charge",&Charge);
  itr->SetBranchAddress("InitMom",&InitMom);
  itr->SetBranchAddress("InitPos",&InitPos);
  itr->SetBranchAddress("nPoint",&nPoint);
  itr->SetBranchAddress("eArr",eArr);
  itr->SetBranchAddress("tArr",tArr);
  itr->SetBranchAddress("xArr",xArr);
  itr->SetBranchAddress("yArr",yArr);
  itr->SetBranchAddress("zArr",zArr);


  TCanvas* can = new TCanvas("can","can",800,800);

  for( int ievt = 0; ievt <  10; ievt++){
    //for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
    itr->GetEntry( ievt );
    //std::cout<< ievt << "\t" << EventID <<  "\t" << nPoint << "\t" << eArr[0] << std::endl;
    std::cout<< ievt << std::endl;
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep( PID );
    genfit::Track fitTrack( rep, InitPos, initMom);
    TGraph* gr = new TGraph();
    for( int ipoint  = 0; ipoint < nPoint; ipoint++){
      gr->SetPoint( ipoint, xArr[ipoint], zArr[ipoint] );
      TVectorD hitCoords(3);
      hitCoords[0] = xArr[ipoint];
      hitCoords[1] = yArr[ipoint];
      hitCoords[2] = zArr[ipoint];
      TMatrixDSym hitCov(3);
      hitCov.UnitMatrix();
      genfit::AbsMeasurement* measurement = new genfit::SpacepointMeasurement( hitCoords, hitCov, 0, ipoint, NULL);
      fitTrack.insertPoint( new genfit::TrackPoint( measurement, &fitTrack));
    }
    if( ievt < 20 ){
      gr->Draw("AP");
      can->Update();
      can->Modifiled();
    }

    assert(fitTrack.checkConsistency());
    fitter->processTrack(&fitTrack);
    fitTrack.getFittedState().Print();
    assert( fitTrack.checkConsistency());
  }
}
