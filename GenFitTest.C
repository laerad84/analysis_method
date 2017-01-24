#include "AbsFinitePlane.h"
#include "AbsFitterInfo.h"
#include "AbsMeasurement.h"
#include "AbsTrackRep.h"
#include "DetPlane.h"

#include "KalmanFittedStateOnPlane.h"
#include "AbsKalmanFitter.h"
#include "KalmanFitter.h"
#include "KalmanFitterRefTrack.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitStatus.h"

#include "DAF.h"
#include "GFGbl.h"

#include "MeasuredStateOnPlane.h"
#include "MeasurementOnPlane.h"
#include "FullMeasurement.h"
#include "PlanarMeasurement.h"
#include "ProlateSpacepointMeasurement.h"
#include "mySpacepointMeasurement.h"
#include "mySpacepointDetectorHit.h"
#include "MeasurementProducer.h"
#include "MeasurementFactory.h"


#include "RectangularFinitePlane.h"
#include "ReferenceStateOnPlane.h"
#include "SharedPlanePtr.h"

#include "SpacepointMeasurement.h"
#include "StateOnPlane.h"

#include "Tools.h"
#include "TrackCand.h"
#include "TrackCandHit.h"
#include "RKTools.h"
#include "RKTrackRep.h"
#include "StepLimits.h"

#include "FieldManager.h"
#include "ConstField.h"
#include "VariableField.h"

#include "Exception.h"
#include "StateOnPlane.h"
#include "Track.h"
#include "TrackPoint.h"
#include "MaterialEffects.h"
#include "EventDisplay.h"

#include <vector>
#include <algorithm>

#include "TSystem.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH2Poly.h"
#include "TMath.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TEveManager.h"
#include "TGeoManager.h"
#include "TGeoMaterialInterface.h"
#include "TRandom.h"

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
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  TVector3 FieldCenter(0,0,0);
  TVector3 FieldRotation(0,0,0);

  //TFile* itf = new TFile("TrackOutputProton.root");
  TFile* itf = new TFile("protonTest_track.root");
  TTree* itr = (TTree*)itf->Get("TrackTree");

  Int_t EventID;
  Int_t TrackID;
  Int_t TrackIDMC;
  Int_t PID;
  Double_t Charge;
  Int_t    nPoint;
  Int_t    DetID[120];
  Double_t eArr[120];
  Double_t tArr[120];
  Double_t xArr[120];
  Double_t yArr[120];
  Double_t zArr[120];
  TVector3* InitMom=NULL;
  TVector3* InitPos=NULL;

  itr->SetBranchAddress("EventID"  ,&EventID);
  itr->SetBranchAddress("TrackID"  ,&TrackID);
  itr->SetBranchAddress("TrackIDMC",&TrackIDMC);
  itr->SetBranchAddress("PID"      ,&PID);
  itr->SetBranchAddress("Charge"   ,&Charge);
  itr->SetBranchAddress("nPoint"   ,&nPoint);
  itr->SetBranchAddress("DetID"    ,DetID);
  itr->SetBranchAddress("eArr"     ,eArr);
  itr->SetBranchAddress("xArr"     ,xArr);
  itr->SetBranchAddress("yArr"     ,yArr);
  itr->SetBranchAddress("zArr"     ,zArr);
  itr->SetBranchAddress("tArr"     ,tArr);
  itr->SetBranchAddress("InitMom"  ,&InitMom);
  itr->SetBranchAddress("InitPos"  ,&InitPos);
  //itr->Print();
  //itr->Show(1);
  //itr->GetEntry(0);
  //InitPos->Print();
  TFile* otf = new TFile("hisDistribution.root","recreate");
  TH2D*  hisResMom = new TH2D("hisResMom","hisResMom",1000,-0.5, 0.5, 1000, 0, 0.5 );/// pm 50 %  and upto 2GeV
  TTree* otr = new TTree("FitResult","fitting result tree");
  TVector3 posInit;
  TVector3 momInit;
  TVector3 posFit;
  TVector3 momFit;
  otr->Branch("posInit",&posInit);
  otr->Branch("posFit" ,&posFit);
  otr->Branch("momInit",&momInit);
  otr->Branch("momFit" ,&momFit);

  //genfit::VariableField* field =   new genfit::VariableField("/Users/jwlee/local/hep/E42/E42/data/field/SC_FieldMap.root",FieldCenter, FieldRotation);
  genfit::FieldManager::getInstance()->init( new genfit::VariableField("/Users/jwlee/local/hep/E42/E42/data/field/SC_FieldMap.root",FieldCenter, FieldRotation));
  //genfit::EventDisplay* display = genfit::EventDisplay::getInstance();

  otf->cd();
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Initial point generator setting
  std::cout<< "Initial pointer generator setting " << std::endl;
  /// Preparation for data handler ///
  //- fitter with KalmanFitter
  //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter();
  //- hit array
  // detector ID - Should be changed by channel?
  // MeasurementFactory : Create Absmeasurment from data
  // MeasurementProducer : Create Measurements for one specific detector types

  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
  TClonesArray myDetectorHitArray("genfit::mySpacepointDetectorHit");
  int myDetId(1);
  genfit::MeasurementFactory<genfit::AbsMeasurement> factory;
  genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::mySpacepointMeasurement> myProducer(& myDetectorHitArray);
  // attach producer to the factory
  //factory.addProducer(1, &myProducer);

  factory.addProducer(1,&myProducer);
  /*
  for( int i = 1; i< 33; i++){
    factory.addProducer(i, &myProducer);
  }
  */

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* can = new TCanvas("can","can", 800,600);

  //TCanvas* can = new TCanvas("can","can",800,800);
  std::cout<< itr->GetEntries() << std::endl;
  //for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
  for( int ievt = 0; ievt < 1000; ievt++){ /// Short test
    //std::cout<< "GetEntry" << std::endl;
    itr->GetEntry( ievt );
    std::cout<< ievt << std::endl;
    if( nPoint < 5 ){ continue; }

    genfit::SharedPlanePtr measurementPlane(new genfit::DetPlane(TVector3(0,0,15), TVector3(0,0,1)));

    // initialization
    if( myDetectorHitArray.GetEntries() != 0 ){
      myDetectorHitArray.Clear();
    }
    genfit::TrackCand myCand;
    // Setting initial position //
    TVector3 initialPos( InitPos->X()*0.1, InitPos->Y()*0.1, InitPos->Z()*0.1);/// mm to cm
    TVector3 initialMom( InitMom->X()*0.001*gRandom->Gaus(1,0.01), InitMom->Y()*0.001*gRandom->Gaus(1,0.01), InitMom->Z()*0.001*gRandom->Gaus(1,0.01));//// MeV to GeV

    //std::cout<< "Set Point" << std::endl;
    double resolution = 1;// cm
    for( int ipoint = 0; ipoint < nPoint; ipoint++){
      TVector3 currentPos( xArr[ipoint]*0.1, yArr[ipoint]*0.1, zArr[ipoint]*0.1);
      TMatrixDSym cov(3);
      for (int i = 0; i < 3; ++i){
	cov(i,i) = resolution*resolution;
      }
      new(myDetectorHitArray[ipoint]) genfit::mySpacepointDetectorHit( currentPos, cov);
      //myCand.addHit(DetID[ipoint]+1, ipoint );
      myCand.addHit(1, ipoint );
    }

    ///insert Initial guess ///



    TVector3 initialGuessV( initialPos.X(), initialPos.Y(), initialPos.Z() );
    TVector3 initialGuessP( initialMom.X(), initialMom.Y(), initialMom.Z() );
    TMatrixDSym covSeed(6);
    /// Position resolution ???
    for( int i = 0; i< 3; i++){
      covSeed(i,i) = resolution* resolution;
    }
    /// Momentum resolution ??? Other resolution
    for( int i = 3; i< 6; i++){
      covSeed(i,i) = pow(resolution/nPoint/sqrt(3),2);
      //covSeed(i,i) = pow( 5, 2 );
    }

    //////
    //InitPos->Print();
    //InitMom->Print();

    genfit::RKTrackRep* rep = new genfit::RKTrackRep( PID );
    myCand.setPosMomSeedAndPdgCode( initialGuessV, initialGuessP, PID);
    myCand.setCovSeed( covSeed );
    genfit::Track fitTrack( myCand, factory, rep );
    // Keep initial state
    genfit::MeasuredStateOnPlane stateRef(rep);
    rep->setPosMomCov(stateRef,initialGuessV,initialGuessP,covSeed);
    const genfit::StateOnPlane   origState( stateRef );
    //std::cout<< "Try fitting" << std::endl;
    try{
      fitter->processTrackWithRep(&fitTrack, rep);
      std::cout<< "TryFit" << std::endl;
      //fitter->processTrack(&fitTrack, rep);
      std::cout<< "End" << std::endl;
    }
    catch( genfit::Exception& e ){
      std::cerr << e.what();
      std::cerr << "Exception, next Track" << std::endl;
      continue;
    }

    if( !(fitTrack.getFitStatus()->isFitted()) ){
      std::cout<< "Get rid of bad fitting." << std::endl;
      continue;
    }

    std::cout<< "Catch" << std::endl;
    //assert( fitTrack.checkConsistency());
    //std::cout<< "Add Fit track" << std::endl;

    //std::cout << "Get TrackPoint" << std::endl;
    genfit::TrackPoint* tp = fitTrack.getPointWithFitterInfo(0, rep );
    if( tp == NULL ){ continue; }
    //std::cout<< "KalmanFittedSateOnPlane" << std::endl;
    //static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate();
    //std::cout<< "KalmanFittedStateOnPlane" << std::endl;
    genfit::StateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
    genfit::MeasuredStateOnPlane fitMeasurement((static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getFittedState()));


    std::cout<< "Position" << std::endl;
    fitMeasurement.getPos().Print();
    std::cout<< "Momentum" << std::endl;
    fitMeasurement.getMom().Print();
    std::cout<< "OriginalState" << std::endl;
    origState.getPlane()->getO().Print();
    ///// Get raw hit point /////
    /*
    std::vector<genfit::TrackPoint*> pointArray = fitTrack.getPointsWithMeasurement();
    std::cout<<"RawMeasurement:" <<  pointArray.at(0)->getNumRawMeasurements() << std::endl;
    for( int ipoint  = 0; ipoint < pointArray.size(); ipoint++){
      genfit::AbsMeasurement* mes = pointArray.at(ipoint)->getRawMeasurement();
      //TVectorD point = mes->getRawHitCoords();
      //std::cout<< point[0] << "\t" << point[1] << "\t" << point[2] << std::endl;
      (mes->getRawHitCoords()).Print();//TVectorD
    }
    */
    genfit::StateOnPlane stateplane(rep);

    //display->addEvent( &fitTrack );

    //std::cout<< "Try extrapolate" << std::endl;
    try{
      //std::cout<< "extrapolate" << std::endl;
      //rep->extrapolateToPlane( kfsop, origState.getPlane());
      rep->extrapolateToPlane( kfsop, measurementPlane);

      //std::cout<< "extrapolate" << std::endl;
    }
    catch(genfit::Exception& e ){
      std::cerr << "Exception, next track" << std::endl;
      std::cerr << e.what();
      continue;
    }

    std::cout<< "Get Plane" << std::endl;
    genfit::SharedPlanePtr p(kfsop.getPlane());
    p->getO().Print();
    //std::cout << "Get position and momentum" << std::endl;

    origState.getPosMom(posInit, momInit);
    kfsop.getPosMom(posFit, momFit);
    momFit.Print();
    otr->Fill();

  }
  myDetectorHitArray.Clear();
  //display->open();
  //std::cout<< "End" << std::endl;
  otr->Write();
  otf->Close();
}
