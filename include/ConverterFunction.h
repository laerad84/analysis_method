#ifndef ConverterFunction__H__
#define ConverterFunction__H__

#include "TPCData.h"
#include "TPCGlobals.h"
#include "Data.h"

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

std::vector<TPCPadHit >       ConvertTPC(GsimDetectorEventData* detData,GsimGenParticleData* particle);
//std::vector<TPCPadHitCluster> HitClustering(std::vector<TPCPadHit> &hitArr);
/// Clusterer
Bool_t HitClusteringB( std::vector<TPCPadHit> &hitArr, std::vector<TPCPadHitCluster> &clusterArr );//will be delete
Bool_t HitClustering( std::vector<TPCPadHit> &hitArr, std::vector<TPCCluster> & clusterArr );///will be delete
Bool_t HitClustering( TClonesArray* HitArr, std::vector<TPCCluster>& clusterArr );///will be delete
Bool_t HitClustering( TClonesArray* HitArr, TClonesArray* clusterArr );
/// Blocker
std::vector<TPCCluster> TPCClusterBlocker( std::vector<TPCCluster> &clArr, Int_t Direction);///will be delete
std::vector<TPCCluster> TPCClusterSingleBlocker( std::vector<TPCCluster> &clArr);/// will be delete

/// Divide block
Int_t  BlockDivider( std::vector<TPCCluster> &clArr, std::vector<TPCCluster> &trackCand, Int_t Mode = 0);/// will be delete

///
Bool_t AdjCluster( TPCCluster c_0, TPCCluster c_1, Int_t DRow);
Bool_t AdjCluster( TPCCluster* c_0, TPCCluster* c_1, Int_t DRow);
void   ConnectionFinder(std::vector<TPCCluster>& clArr);/// will be delete
std::vector< std::vector<TPCCluster> > ConnectionBlocker( std::vector<TPCCluster> &clArr );///will bedelete
std::vector<Int_t> GetListOfTrackRoot( TClonesArray* clArr );
Bool_t ClusterBlocker( TClonesArray* clrr, TClonesArray* blockArr, std::vector<Int_t> BlockRoot, Int_t BlockIndex );

#endif //ConverterFunction__H__
