#ifndef TPCSprial__H__
#define TPCSprial__H__
#include "TPbject.h"
#include "TArc.h"
#include "TMath.h"
#include "CalculationFunction.h"

class TPCSprial : public TObject {
 public:
  TPCSprial();
  TPCSprial(const TPCSprial& right);
  ~TPCSprial();

  Double_t x;
  Double_t y;
  Double_t z;
  Double_t DThetaDy;
  Double_t R;

  TVector3 InitialPoint;
  TVector3 MiddlePoint;
  TVector3 FinalPoint;

  Int_t    RL;//RightHanded : -1, LeftHanded : 1

  void GenerateSprial(TVector3 init, TVector3 mid, TVector finl);

  void DrawXY();
  void DrawYZ();
  void DrawZX();

 public:
  ClassDef( TPCSprial, 0)
};

TVector2   ThreePointCircle(TVector3 init, TVector3 middle, TVector3 final);
TPCSprial* generateSprial(TVector3 initial, TVector3 middle, TVector3 final);
#endif

TPCSprial::TPCSprial(): x(0),y(0),z(0),DThetaDy(0),R(0),InitialPoint(0,0,0),MiddlePoint(0,0,0),FinalPoint(0,0,0){
  ;
}
TPCSprial::TPCSprial( const TPCSprial& right){
  x = right.x;
  y = right.y;
  z = right.z;
  DThetaDy = right.DThetaDy;
  R = right.R;
  InitialPoint(right.InitialPoint);
  MiddlePoint(right.MiddlePoint);
  FinalPoint(right.FinalPoint);
}
TPCSprial::~TPCSprial(){
  ;
}

void TPCSprial::GenerateSprial(TVector3 init, TVector3 mid, TVector3 final){
  ;
}
