#ifndef Data_H__
#define Data_H__
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TArc.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include <iostream>
#include "TPCIDManager.h"

class HCircle : public TObject {
 public:
  HCircle();
  HCircle(Int_t id, Double_t x, Double_t z, Double_t r);
  HCircle(const HCircle& right);
  virtual ~HCircle();

  Int_t    ID;
  Double_t X;
  Double_t Z;
  Double_t R;
  //virtual void Draw();
  void Print();
 public:
  ClassDef( HCircle, 1 )
};

class HHelix : public TObject {
 public:
  HHelix();
  HHelix(Int_t id, Double_t x, Double_t y, Double_t z, Double_t r, Double_t dY, Double_t dTheta, Int_t rl);
  HHelix(const HHelix& right);
  virtual ~HHelix();

  Int_t    ID;        /// Helix ID
  Double_t X;         /// Helix center X
  Double_t Y;         /// Helix starting Y
  Double_t Z;         /// Helix center Z
  Double_t R;         /// Helix Radius
  Double_t EX;        /// Helix X error
  Double_t EZ;        /// Helix Z error
  Double_t EY;        /// Helix Y error
  Double_t ER;        /// Helix R error
  Double_t EDYDTheta; /// DYDTheta error
  Double_t DTheta;    /// Theta changes from initial position to final position
  Double_t DY;        /// Y changes from initial position to final position
  Int_t    RL;        /// Right curved and left curved
  TVector3 InitPos;   /// initial position
  TVector3 InitMom;   /// estimated momentum
  TVector3 FinalPos;  /// final position
  TVector3 FitMom;    /// fitted momentum

  ///// parameters will be removed
  Int_t    nCluster;    /// number of clusters in the helix
  Double_t TrackLength; /// calculated track length
  Double_t DepE;        /// energy deposition

  void ResetAll();
  void Print();
  TPolyLine3D* GenerateHelix();
  TPolyLine*   GenerateArc();
  TF1*         GenerateXY();
  TF1*         GenerateZY();
  void         CalculateDist(TVector3 hitpos, Double_t& dR, Double_t& dY);
  void         CalculateDist(TGraph2D* grTrack, Double_t& dR, Double_t& dY);
  //virtual void DrawXY();
  //virtual void DrawYZ();
  //virtual void DrawZX();
 public:
  ClassDef( HHelix, 1 )
};

class HTrack : public TObject {
 public:
  HTrack();
  HTrack(const HTrack& right);
  virtual ~HTrack();
  Int_t     ID;
  TGraph2D  track;
  TVector3  Momentum;
  Double_t  R;
  Double_t  DYDL;
  Double_t  DYDT;
  Double_t  Length;
  Double_t  DepE;
 public:
  ClassDef( HTrack, 1)
};

class HLine : public TObject {
 public:
  HLine();
  HLine(Int_t id, TVector3 slope, TVector3 offset);
  HLine(const HLine& right);
  virtual ~HLine();
  // ax+by+cz+offset= 0;
  Int_t ID;
  TVector3 Direction;
  TVector3 ClosestPoint;
  TVector3 Offset;
  Bool_t   bNorm;

  void Normalize();
  void Print();
  //virtual void DrawXY();
  //virtual void DrawYZ();
  //virtual void DrawZX();

 public:
  ClassDef( HLine, 1)
};
/*
class TPCHit : public TObject {
 public:
  TPCHit();
  TPCHit( int id, Double_t y, Double_t energy );
  TPCHit( const TPCHit& right );
  virtual ~TPCHit();
  Int_t ID;
  Int_t Row;
  Int_t Col;
  Double_t Y;
  Double_t Energy;
  TVector3 Position;
  TVector3 PositionError;

  Int_t PadID;
  Int_t ClusterID;

 public:
  ClassDef( TPCHit, 0 )
};
*/

////////////////////////////////////////////////////////////////////////////////////////////////////
/// TPCHit : Hit class
////////////////////////////////////////////////////////////////////////////////////////////////////
class TPCHit : public TObject {
 public:
  TPCHit();
  TPCHit(const TPCHit& right);
  TPCHit(int padID, TVector3 pos, double energy);
  virtual ~TPCHit();
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
  ClassDef( TPCHit, 1)
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// TPCCluster : Cluster class
////////////////////////////////////////////////////////////////////////////////////////////////////
class TPCCluster : public TObject {
 public:
  TPCCluster();
  TPCCluster( int id, Double_t y, Double_t energy);
  TPCCluster( const TPCCluster& right );
  virtual ~TPCCluster();

  Int_t ID;                     /// Cluster parameters
  ////// Cluster parameters //////
  Int_t Row;
  std::vector<Int_t>   ColArr;
  std::vector<Double_t> EArr;
  Double_t Energy;
  Double_t Col;
  Double_t ColMin;
  Double_t ColMax;
  TVector3 Position;
  TVector3 PositionErr;
  Double_t YMin;
  Double_t YMax;
  Double_t PhiMax;
  Double_t PhiMin;

  Int_t    NHit;                /// Number of hits in cluster

  ////// Mother & Daughter parameters //////
  Int_t    Mother;
  Int_t    Daughter;
  Double_t MinMotherDist;
  Double_t MinDaughterDist;
  Double_t MotherY;
  Double_t DaughterY;

  Int_t    nMother;              /// Number of Mothers
  Int_t    nDaughter;            /// Number of Daughters
  std::vector<Int_t> MotherID;   /// Array of mother cluster's ID
  std::vector<Int_t> DaughterID; /// Array of Daughter cluster's ID
  Bool_t   Blocked;              /// Blocked Indicator
  void     Print();
 public:
  ClassDef( TPCCluster, 1 )
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// TPCTrack : track data for fitting
////////////////////////////////////////////////////////////////////////////////////////////////////
class TPCTrack : public TObject {
 public:
  TPCTrack();
  TPCTrack(const TPCTrack& right);
  virtual ~TPCTrack();
  Int_t     ID;              /// Track ID
  Int_t     InitialRowID;    /// Row ID of the starting cluster
  Int_t     FinalRowID;      /// Row ID of the final cluster
  Int_t     PID;             /// Estimated particle id
  Int_t     nCluster;        /// Number of clusters
  TGraph2D* track;           /// 3d track points
  TGraph2D* trackErr;        /// Not used currently
  std::vector<Double_t> EArr;/// Energy array
  HHelix    helix;           /// Helix parameter
  TVector3  Momentum;        /// Measured momentum
  TVector3  MomentumErr;     /// Measure momentum Error from fitting
  Double_t  Length;          /// path length
  Double_t  DepE;            /// deposit energy
  Double_t  ChisqY;          /// Chisquare Y
  Double_t  ChisqR;          /// Chisquare R
  void Init();
  void Print();
 public:
  ClassDef( TPCTrack, 1)
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// HPiM : anti pion
////////////////////////////////////////////////////////////////////////////////////////////////////
class HPiM : public TObject {
 public:
  HPiM();
  HPiM( const HPiM& right );
  virtual  ~HPiM();
  Int_t          ID;
  Double_t       RecMass;           /// PDG value
  Double_t       Energy;
  TLorentzVector Momentum;          /// Momentum
  TVector3       Position;

  /// Track;
  TPCTrack track;
  void Init();
 public:
  ClassDef( HPiM, 1 )
};
////////////////////////////////////////////////////////////////////////////////////////////////////
/// HKaonP : kaon+
////////////////////////////////////////////////////////////////////////////////////////////////////
class HKaonP : public TObject {
 public:
  HKaonP();
  HKaonP( const HKaonP& right );
  virtual ~HKaonP();
  Int_t          ID;
  Double_t       RecMass;
  Double_t       Energy;
  TLorentzVector Momentum;
  TVector3       Position;
  /// Track;
  TPCTrack  track;
  void Init();
 public:
  ClassDef( HKaonP, 1 )
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// HProton : proton
////////////////////////////////////////////////////////////////////////////////////////////////////
class HProton : public TObject {
 public:
  HProton();
  HProton( const HProton& right );
  virtual ~HProton();
  Int_t          ID;
  Double_t       RecMass;
  Double_t       Energy;
  TLorentzVector Momentum;
  TVector3       Position;
  /// Track;
  TPCTrack       track;
  void Init();
 public:
  ClassDef( HProton, 1 )
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// HLambda : Lambda
////////////////////////////////////////////////////////////////////////////////////////////////////
class HLambda : public TObject{
 public:
  HLambda();
  HLambda( const HLambda& right );
  virtual ~HLambda();
  Int_t          ID;          /// lambda ID
  Double_t       RecMass;     /// reconstructed lambda mass with two tracks
  Double_t       Energy;
  TLorentzVector Momentum;    /// momentum of lambda, energy with PDG lambda mass
  TVector3       Position;    /// Decay Position of lambda.
  HProton*       proton;
  HPiM*          pion;

  void Init();
 public:
  ClassDef( HLambda, 1 )
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// HCascade : Cascade
////////////////////////////////////////////////////////////////////////////////////////////////////
class HCascade : public TObject{
 public:
  HCascade();
  HCascade( const HCascade& right );
  virtual ~HCascade();
  Int_t          ID;          /// lambda ID
  Double_t       RecMass;     /// reconstructed lambda mass with two tracks
  Double_t       Energy;
  TLorentzVector Momentum;    /// momentum of lambda, energy with PDG lambda mass
  TVector3       Position;    /// Decay Position of lambda.
  HLambda*       lambda;      /// lambda
  HPiM*          pion;        /// pion

  void           Init();
 public:
  ClassDef( HCascade, 1 )
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// HDibaryonLL : Dibaryon reconstructed with 2 lambda
////////////////////////////////////////////////////////////////////////////////////////////////////
class HDibaryonLL : public TObject {
 public:
  HDibaryonLL();
  HDibaryonLL( const HDibaryonLL& right );
  virtual ~HDibaryonLL();
  Int_t          ID;        /// dibaryon ID
  Double_t       RecMass;   /// reconstructed HD mass
  Double_t       Energy;
  TLorentzVector Momentum;  /// momentum of HD.
  TVector3       Position;  /// decay vertex of HD
  HLambda*       lambda0;   /// lambda
  HLambda*       lambda1;   /// lambda

  void Init();
 public:
  ClassDef( HDibaryonLL, 1 )
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// HDibaryonPC : Dibaryon reconstructed with proton and cascade
////////////////////////////////////////////////////////////////////////////////////////////////////
class HDibaryonPC : public TObject {
 public:
  HDibaryonPC();
  HDibaryonPC( const HDibaryonPC& right );
  virtual  ~HDibaryonPC();
  Int_t          ID;
  Double_t       RecMass;
  Double_t       Energy;
  TLorentzVector Momentum;
  TVector3       Position;
  HCascade*      cascade;
  HPiM*          pion;

  void Init();
 public:
  ClassDef( HDibaryonPC, 1 )
};

#endif
