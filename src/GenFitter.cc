#include "GenFitter.h"

ClassImp(GenFitter)

GenFitter::GenFitter( char* geometryFilename ){

  //new TGeoManager("TPCGeom","E42 Geometry");
  //gGeoManager->Import(geometryFilename);
  //genfit::MaterialEffects::getInstance()->init( new genfit::TGeoMaterialInterface() );
  Init(geometryFilename);
  /*
  genfit::Fieldmanager::GetInstance()->init( new genfit::VariableField("/Users/jwlee/local/hep/E42/E42/data/field/SC_FieldMap.root", FieldCenter), FieldRotation);
  fitter = new genfit::KalmanFitterRefTrack();
  factory.addProducer(1, &myProducer);
  myDetectorHitArray = new TClonesArray("genfit::mySpacepointDetectorHit");
  */

}

GenFitter::~GenFitter(){
  delete myProducer;
  delete myDetectorHitArray;
  delete factory;
  delete fitter;
}

void GenFitter::Init(char* geometryFilename ){
  new TGeoManager("TPCGeom", "E42 Geometry");
  gGeoManager->Import( geometryFilename );
  genfit::MaterialEffects::getInstance()->init( new genfit::TGeoMaterialInterface() );
  fitter             = new genfit::KalmanFitterRefTrack();
  factory            = new genfit::MeasurementFactory<genfit::AbsMeasurement>();
  myDetectorHitArray = new TClonesArray("genfit::mySpacepointDetectorHit");
  myProducer         = new genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::mySpacepointMeasurement>( myDetectorHitArray);
  factory->addProducer(1,myProducer);
}

void GenFitter::SetField(char* fieldFilename, TVector3 fieldCenter, TVector3 fieldRotation){
  genfit::FieldManager::getInstance()->init( new genfit::VariableField(fieldFilename, fieldCenter,fieldRotation));
}

Bool_t GenFitter::Fit( TClonesArray* clarr, HSprial sprial , Int_t PID, TVector3& posFit, TVector3& momFit ){
  if ( PID != 2212 &&
       PID != 321  &&
       PID != 211  &&
       PID != -211 ){
    std::cerr << "PID is not in the list " << std::endl;
  }


  std::cout<< __PRETTY_FUNCTION__ << ":" << __LINE__ << std::endl;

  Bool_t bFit = false;
  //// set point
  std::cout<< __PRETTY_FUNCTION__ << ":" << __LINE__ << std::endl;

  genfit::SharedPlanePtr measurementPlane( new genfit::DetPlane(TVector3(0,0,0),TVector3(0,0,1)));

  std::cout<< __PRETTY_FUNCTION__ << ":" << __LINE__ << std::endl;


  //if( myDetectorHitArray->GetEntries() != 0 ){ myDetectorHitArray->Clear(); }
  std::cout<< __PRETTY_FUNCTION__ << ":" << __LINE__ << std::endl;

  genfit::TrackCand myCand;
  std::cout<< __PRETTY_FUNCTION__ << ":" << __LINE__ << std::endl;

  Double_t resolution = 1;
  Int_t number_of_point = clarr->GetEntries();
  for( int index = 0; index < clarr->GetEntries(); index++ ){
    std::cout<< index << std::endl;
    TPCCluster* cl = (TPCCluster*)clarr->At(index);
    std::cout<< index << std::endl;
    TVector3 HitPoint( cl->Position.X()*0.1,
		       cl->Position.Y()*0.1,
		       cl->Position.Z()*0.1);
    TMatrixDSym cov(3);
    for( int i = 0; i< 3; i++){
      cov(i,i) = resolution*resolution;
    }
    std::cout<< index << std::endl;
    new((*myDetectorHitArray)[index]) genfit::mySpacepointDetectorHit( HitPoint, cov );
    std::cout<< index << std::endl;
    myCand.addHit(1,index);
  }
  std::cout<< __PRETTY_FUNCTION__ << ":" << __LINE__ << std::endl;

  TMatrixDSym covSeed(6);
  for( int i = 0; i< 3; i++){
    covSeed(i,i) = resolution*resolution;
  }
  for( int i = 3; i< 6; i++){
    covSeed(i,i) = TMath::Power( resolution / number_of_point / sqrt(3), 2 );
  }
  std::cout<< __PRETTY_FUNCTION__ << ":" << __LINE__ << std::endl;

  genfit::RKTrackRep* rep = new genfit::RKTrackRep( PID );
  myCand.setPosMomSeedAndPdgCode( sprial.InitPos, sprial.InitMom, PID);
  myCand.setCovSeed( covSeed );
  genfit::Track fitTrack( myCand, *factory, (genfit::AbsTrackRep*)rep );
  genfit::MeasuredStateOnPlane stateRef( (genfit::AbsTrackRep*)rep );
  rep->setPosMomCov( stateRef, sprial.InitPos, sprial.InitMom, covSeed );
  const genfit::StateOnPlane stateref( (genfit::AbsTrackRep*)rep );
  std::cout<< __PRETTY_FUNCTION__ << ":" << __LINE__ << std::endl;

  try{
    fitter->processTrackWithRep(&fitTrack, (genfit::AbsTrackRep*)rep);
  }catch( genfit::Exception& e ){
    std::cerr << e.what();
    std::cerr << "Exception, next Track" << std::endl;
    return bFit;
  }

  if( !(fitTrack.getFitStatus())->isFitted()){
    std::cerr << "Bad Error" << std::endl;
    return bFit;
  }

  genfit::TrackPoint* tp = fitTrack.getPointWithFitterInfo( 0, (genfit::AbsTrackRep*)rep );
  if( tp == NULL ){
    return bFit;
  }else{
    bFit = true;
  }
  std::cout<< __PRETTY_FUNCTION__ << ":" << __LINE__ << std::endl;

  genfit::StateOnPlane kfsop(*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo((genfit::AbsTrackRep*)rep))->getBackwardUpdate()));
  genfit::MeasuredStateOnPlane fitMeasurement((static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo((genfit::AbsTrackRep*)rep))->getFittedState()));
  kfsop.getPosMom( posFit, momFit );
  std::cout<< __PRETTY_FUNCTION__ << ":" << __LINE__ << std::endl;

  myDetectorHitArray->Clear();
  return bFit;
}
