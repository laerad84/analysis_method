#ifndef TPCMCDataConverter__H__
#define TPCMCDataConverter__H__
#include "GsimData/GsimDetectorEventData.h"
#include "GsimData/GsimGenParticleData.h"
#include "GsimData/GsimTrackData.h"
#include "GsimData/GsimDetectorHitData.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TTree.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include <vector>

#include "TPCData.h"
#include "TPCGlobals.h"

class TPCDataConverter {
 public:
  TPCDataConverter(TTree* tr);
  virtual ~TPCDataConverter();
  GsimDetectorEventData* GetDetectorEventData();
  //std::vector<TPCPadHit> Convert( int ievent );
  void Convert( int ievent );
  void Clustering();
  std::vector<TPCPadHit>        GetPadHitArr() const { return m_PadHitArr; }
  std::vector<TPCPadHitCluster> GetClusterArr() const { return m_ClusterArr; }

  void SetData();
  void ClearSubData();


  TGraph2D* grHit;
  TGraph2D* grPadHit;
  TGraph2D* grCluster;
  TGraph2D* grTrack[16];
  TGraph*   grPadHitXZ;
  TGraph*   grClusterXZ;
  TGraph*   grPadHitYZ;
  TGraph*   grClusterYZ;

  TClonesArray* GetDetDigi( int id );
  Int_t         AddDetector( char* detName );

 private:
  std::vector<TPCPadHit> m_PadHitArr;
  std::vector<TPCPadHitCluster> m_ClusterArr;

  TClonesArray* m_TrackArr;
  TClonesArray* m_HitArr;
  GsimDetectorEventData* m_det;
  GsimDetectorEventData* m_digiDet[32];
  GsimGenParticleData*   m_particle;
  TTree*                 m_tr;
  TPCPoly*               m_poly;
  TPCRowYHistArray*      m_RowYHistArr;

  Int_t                  m_nDet;

 public:
  ClassDef(TPCDataConverter,1)
};


#endif //TPCMCDataConverter__H__
