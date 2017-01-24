#ifndef TPCData__H__
#define TPCData__H__
#include "TObject.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TPCIDManager.h"
#include "TPCMinimizer.h"
#include "TPCGlobals.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <vector>

//class TPCPadHit {
class TPCPadHit : public TObject {
 public:
  TPCPadHit();
  TPCPadHit(const TPCPadHit& right);
  TPCPadHit(int padID, TVector3 pos, double energy);
  virtual ~TPCPadHit();
  virtual void  Clear(Option_t* opt="");

  void SetTrackID(int trackID){ m_TrackID = trackID; }
  void SetClusterID( int clusterID){ m_ClusterID = clusterID; }
  void SetHitID(int hitID){ m_HitID = hitID; }
  void SetMCTrackID( int trackID){ m_MCTrackID = trackID; }
  void SetPadID( int padID ){ m_PadID = padID; }

  Int_t TrackID()   const { return m_TrackID; }
  Int_t ClusterID() const { return m_ClusterID; }
  Int_t HitID()     const { return m_HitID; }
  Int_t MCTrackID() const { return m_MCTrackID; }
  Int_t PadID()     const { return m_PadID;}

  void SetPosition( TVector3 pos ){ m_pos = pos; }
  void SetPosError( TVector3 posErr ){ m_posErr = posErr; }
  void SetEnergy( double energy ){ m_energy = energy; }
  void SetRowCol( int row, int col ){ m_Row = row; m_Col = col; }
  void SetRad( double rad ){ m_Rad = rad;}
  void SetPhi( double phi ){ m_Phi = phi;}

  TVector3 Position() const { return m_pos; }
  TVector3 PosError() const { return m_posErr; }
  Double_t Energy()   const { return m_energy; }
  Int_t    Row()      const { return m_Row; }
  Int_t    Col()      const { return m_Col; }
  Double_t Phi()      const { return -180 + TPCGlobals::sTPC_Pad_Parameter[m_Row][4] + (m_Col+0.5)*TPCGlobals::sTPC_Pad_Parameter[m_Row][3];}
  Double_t Rad()      const { return TPCGlobals::sTPC_Pad_Parameter[m_Row][2]; }

  Int_t    m_TrackID;
  Int_t    m_ClusterID;
  Int_t    m_MCTrackID;
  Int_t    m_HitID;

  Int_t    m_PadID;
  Int_t    m_Row;
  Int_t    m_Col;
  Double_t m_Phi;
  Double_t m_Rad;

  TVector3 m_pos;
  TVector3 m_posErr;

  Double_t m_energy;
  Double_t m_signal;
  Double_t m_rawTiming;
  Double_t m_timing;

  Bool_t   IsClustered;
  Bool_t   IsTracked;
 public:
  ClassDef( TPCPadHit, 1)
};


class TPCPadHitCluster : public TObject {
 public:
  TPCPadHitCluster();
  TPCPadHitCluster(int rowID);
  TPCPadHitCluster(const TPCPadHitCluster& right);
  virtual ~TPCPadHitCluster();
  virtual void  Clear(Option_t* opt="");


  Int_t    m_TrackID;
  Int_t    m_ClusterID;
  Int_t    m_MCTrackID;
  Int_t    m_MotherID;
  std::vector<Int_t>     m_DaughterID;
  std::vector<TPCPadHit> m_HitArr;

  Int_t    m_Row;
  Double_t m_Col;
  Int_t    m_ColMin;
  Int_t    m_ColMax;
  Double_t m_YMin;
  Double_t m_YMax;
  Double_t m_PhiMin;
  Double_t m_PhiMax;
  std::vector<Int_t>     m_ColArr;
  std::vector<Double_t>  m_EArr;
  Int_t    m_nHit;
  TVector3 m_ClusterPos;
  TVector3 m_ClusterPosErr;
  Double_t m_TotalEnergy;

  Bool_t   IsOuterCluster;// Checking Cluster was most outer row of TPC
  Bool_t   IsTracked;// Checking Cluster was tracked
  Bool_t   IsUpdated;// Checking Cluster data updated

  void SetTrackID(int trackID){ m_TrackID = trackID; IsTracked = true;}
  void SetClusterID( int clusterID){ m_ClusterID = clusterID;}
  void SetMCTrackID( int trackID){ m_MCTrackID = trackID; }

  Int_t TrackID()   const { return m_TrackID; }
  Int_t MCTrackID() const { return m_MCTrackID; }
  Int_t ClusterID() const { return m_ClusterID; }
  std::vector<Int_t>     ColIDArr() const { return m_ColArr; }
  std::vector<TPCPadHit> HitArr()   const { return m_HitArr; }

  Double_t    Energy() const { return m_TotalEnergy;}
  Int_t       RowID()  const { return m_Row; }
  Double_t    ColID()  const { return m_Col; }
  Int_t       ColMax() const { return m_ColMax; }
  Int_t       ColMin() const { return m_ColMin; }
  Double_t    YMax()   const { return m_YMax; }
  Double_t    YMin()   const { return m_YMin; }
  Double_t    PhiMin() const { return m_PhiMin; }
  Double_t    PhiMax() const { return m_PhiMax; }
  Int_t       nHit()          const { return m_nHit; }
  TVector3    GetClusterPos() const { return m_ClusterPos; }
  Int_t       Mother()        const { return m_MotherID; }
  std::vector<Int_t> Daughter() const { return m_DaughterID; }
  void        AddPadHit(TPCPadHit hit);
  void        Evaluate();

  void        SetMother( Int_t motherID ){ m_MotherID = motherID; }
  Int_t       AddDaughter( Int_t daughterID ){ m_DaughterID.push_back( daughterID );return m_DaughterID.size();}
  bool        operator < ( const TPCPadHitCluster& cl ) const { return (RowID() < cl.RowID()); }
  //TPCPadHitCluster& operator = ( const TPCPadHitCluster& right );
  Bool_t      IsSortable() const { return kTRUE;}

  Int_t       Compare( TObject* obj){
    int this_row = this->RowID();
    int that_row = ((TPCPadHitCluster*)obj)->RowID();
    if( this_row > that_row ){ return 1;}
    else if( this_row == that_row){ return 0;}
    else{ return -1;}
  }
 public:
  ClassDef(TPCPadHitCluster, 0)
};

TPCPadHitCluster  TPCClusterer( int rowID, std::vector<TPCPadHit> &hitArr);
TPCPadHitCluster  TPCClusterer( int rowID, std::vector<TPCPadHit> &hitArr, Int_t MCTrackID );
class TPCTrack : public TObject {
 public:
  TPCTrack();
  ~TPCTrack();

  double FitCircle();
  double getCircX()     const {return m_circX; }
  double getCircZ()     const {return m_circZ; }
  double getCircR()     const {return m_circR; }
  double getCircXErr()  const {return m_circXErr; }
  double getCircZErr()  const {return m_circZErr; }
  double getCircRErr()  const {return m_circRErr; }
  double getCircChisq() const {return m_circFitChisq; }

 private:
  Int_t  m_TrackID;
  Int_t  m_nCluster;
  std::vector<TPCPadHitCluster> m_ClusterArr;

  Double_t m_circX;
  Double_t m_circZ;
  Double_t m_circR;
  Double_t m_circXErr;
  Double_t m_circZErr;
  Double_t m_circRErr;
  Double_t m_circFitChisq;

  Bool_t  m_circFit;

 public:
  ClassDef( TPCTrack, 0 )
};

std::vector<TPCPadHitCluster> TPCTrackBlocking(std::vector<TPCPadHitCluster> &clusterArr,Int_t Direction=0);//0 inner->outer 1 outer->inner
std::vector<TPCPadHitCluster> TPCRootClusters(std::vector<TPCPadHitCluster> clusterArr);
std::vector<TPCPadHitCluster> TPCEdgeClusters(std::vector<TPCPadHitCluster> clusterArr);
std::vector<TPCPadHitCluster> TPCFindEdgeBlock( std::vector<TPCPadHitCluster> &clusterArr );
std::vector<TPCPadHitCluster> TPCFindBlock(  std::vector<TPCPadHitCluster> &clusterArr, int method = 0);
std::vector<TPCPadHitCluster> TPCFindBlockCircle( std::vector<TPCPadHitCluster> &clusterArr, int Change = 4);
TGraph* DrawClusterArrXZ( std::vector<TPCPadHitCluster> clusterArr);
TGraph* DrawClusterArrYZ( std::vector<TPCPadHitCluster> clusterArr);
Bool_t  AdjacentCluster( TPCPadHitCluster c_0, TPCPadHitCluster c_1, int drow );
void    ClusteringCalculator( TPCPadHitCluster c_0, TPCPadHitCluster c_1, double& length, double& rho, double& phi, double& dydl);
bool    CircleCrossing( TVector3 rxy0, TVector3 rxy1, TVector2& CrossingVec1, TVector2& CrossingVec2);
std::vector<TVector2>  CircleTangent( TVector3 rxy0, TVector3 rxy1);
bool    yDistCalculator(TPCCircleFitResult rst0, TPCCircleFitResult rst1, Double_t& dist, TVector3& Pos);

/*
class TPCPad : public TObject {
 public:
  TPCPad();
  TPCPad(int padID, double depE, double y );
  TPCPad( const TPCPad* right);
  virtual ~TPCPad();

  void Init();
  void Evaluation();

  Int_t id;
  Double_t energy;
  Double_t hitY;

  Int_t row;
  Int_t col;
};

class TPCCluster : public TObject {
  TPCCluster();
  TPCCluster(int rowID);
  TPCCluster( const TPCCluster& right);
  virtual ~TPCCluster();

  Int_t row;
  Double_t col;

};

*/
#endif //TPCData__H__
