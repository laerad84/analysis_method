#include "AbsFinitePlane.h"
#include "AbsFitterInfo.h"
#include "AbsMeasurement.h"
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
#include "AbsTrackRep.h"
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
#ifndef GENFITTER__H__
#define GENFITTER__H__

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
#include "TClonesArray.h"
#include "GsimData/GsimGenParticleData.h"
#include "GsimData/GsimDetectorEventData.h"
#include "GsimData/GsimDetectorHitData.h"
#include "GsimData/GsimDigiData.h"
#include "GsimData/GsimTimeData.h"
#include "GsimData/GsimTrackData.h"
#include <iostream>

#include "Data.h"
#include "TPCData.h"


class GenFitter {
 public:

  GenFitter();//( char* geometryFilename );
  GenFitter( char* geometryFilename );
  ~GenFitter();


  void   Init(char* geometryFilename);
  void   SetField( char* fieldFilename, TVector3 fieldCenter, TVector3 fieldRotation);
  Bool_t Fit( TClonesArray* clarr, HHelix sprial , Int_t PID, TVector3& posFit, TVector3& momFit );

  genfit::AbsKalmanFitter* fitter;
  genfit::MeasurementFactory<genfit::AbsMeasurement>* factory;
  genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::mySpacepointMeasurement> *myProducer;
  TClonesArray* myDetectorHitArray;//("genfit::mySpacepointDetectorHit");

  ClassDef( GenFitter, 1 )
};

R__EXTERN GenFitter* gGenFitter;

#endif //GENFITTER__H__
