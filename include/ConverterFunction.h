#ifndef ConverterFunction__H__
#define ConverterFunction__H__

//#include "TPCData.h"
#include "TPCGlobals.h"
#include "Data.h"
#include "Fitter.h"

#include "GsimData/GsimDetectorEventData.h"
#include "GsimData/GsimGenParticleData.h"
#include "GsimData/GsimTrackData.h"
#include "GsimData/GsimDetectorHitData.h"

#include "TClonesArray.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TClass.h"
#include <vector>
#include <iostream>


std::vector<TPCHit > ConvertTPC(GsimDetectorEventData* detData,GsimGenParticleData* particle);/// Convert GsimData to PadHit.

/// Clusterer
Bool_t HitClustering( TClonesArray* HitArr, TClonesArray* clusterArr );/// Clustering hits
TPCTrack ConvertToTrack( TClonesArray* clusterArr);/// Tracking cluster

Bool_t AdjCluster( TPCCluster c_0, TPCCluster c_1, Int_t DRow);/// Check adjacent cluster
Bool_t AdjCluster( TPCCluster* c_0, TPCCluster* c_1, Int_t DRow);/// Check adjacent cluster
std::vector<Int_t> GetListOfTrackRoot( TClonesArray* clArr );
Bool_t ClusterBlocker( TClonesArray* clrr, TClonesArray* blockArr, std::vector<Int_t> BlockRoot, Int_t BlockIndex );///
////////////////////////////////////////////////////////////
/// TrackHandler : function class for handling track data
////////////////////////////////////////////////////////////
class TrackHandler : public TObject {
 public:
  TrackHandler();
  virtual    ~TrackHandler();
  Bool_t     ConversionToTrack( TClonesArray* clusterArr , TPCTrack& track);
  void       EvalTrackInit( TPCTrack& track, TVector3 initPoint, Double_t bfield = 1.0 );//// Set track parameters.
  TPCTrack*  MergingTrack( TPCTrack* trk0, TPCTrack* trk1 );
 public:
  ClassDef( TrackHandler, 0)
};

R__EXTERN TrackHandler* gTrackHandler;



#endif //ConverterFunction__H__
