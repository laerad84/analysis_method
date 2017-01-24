#ifndef Data_H__
#define Data_H__
#include "TVector3.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TArc.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TMath.h"
#include "TF1.h"
#include <iostream>

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
  ClassDef( HCircle, 0 )
};

class HSprial : public TObject {
 public:
  HSprial();
  HSprial(Int_t id, Double_t x, Double_t y, Double_t z, Double_t r, Double_t dY, Double_t dTheta, Int_t rl);
  HSprial(const HSprial& right);
  virtual ~HSprial();

  Int_t    ID;
  Double_t X;
  Double_t Y;
  Double_t Z;
  Double_t R;
  Double_t DTheta;
  Double_t DY;
  Int_t    RL;
  TVector3 InitPos;
  TVector3 InitMom;
  TVector3 FinalPos;
  TVector3 FitMom;

  Int_t    nCluster;
  Double_t TrackLength;
  Double_t DepE;

  void Print();
  TPolyLine3D* GenerateHelix();
  TPolyLine*   GenerateArc();
  TF1*         GenerateXY();
  TF1*         GenerateZY();
  //virtual void DrawXY();
  //virtual void DrawYZ();
  //virtual void DrawZX();
 public:
  ClassDef( HSprial, 1 )
};

class HLine : public TObject {
 public:
  HLine();
  HLine(Int_t id, Double_t slopeX, Double_t slopeY, Double_t slopeZ, Double_t offset);
  HLine(const HLine& right);
  virtual ~HLine();
  // ax+by+cz+offset= 0;
  Int_t ID;
  Double_t SlopeX;
  Double_t SlopeY;
  Double_t SlopeZ;
  Double_t Offset;

  void Print();
  //virtual void DrawXY();
  //virtual void DrawYZ();
  //virtual void DrawZX();

 public:
  ClassDef( HLine, 0)
};

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
 public:
  ClassDef( TPCHit, 0 )
};

class TPCCluster : public TObject {
 public:
  TPCCluster();
  TPCCluster( int id, Double_t y, Double_t energy);
  TPCCluster( const TPCCluster& right );
  virtual ~TPCCluster();
  Int_t ID;
  Int_t Row;
  Double_t Energy;

  Double_t Col;
  Double_t ColMin;
  Double_t ColMax;

  TVector3 Position;
  Double_t YMin;
  Double_t YMax;
  Double_t PhiMax;
  Double_t PhiMin;

  Int_t    NHit;

  Int_t    Mother;
  Int_t    Daughter;
  Double_t MinMotherDist;
  Double_t MinDaughterDist;
  Double_t MotherY;
  Double_t DaughterY;

  Int_t    nMother;
  Int_t    nDaughter;
  std::vector<Int_t> MotherID;
  std::vector<Int_t> DaughterID;
  Bool_t   Blocked;
  void     Print();
 public:
  ClassDef( TPCCluster, 1 )
};


class TempCluster : public TObject{
 public:
  TempCluster();
  TempCluster( const TempCluster& right );
  virtual ~TempCluster();
  void Init();
  Int_t ID;
  Int_t nMother;
  Int_t Mother[16];
  Int_t nDaughter;
  Int_t Daughter[16];
  std::vector<Int_t> idList;
 public:
  ClassDef( TempCluster, 1)
};


#endif
