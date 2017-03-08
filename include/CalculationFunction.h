#ifndef CalculationFunction__H__
#define CalculationFunction__H__
#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TArc.h"
#include "TPolyMarker3D.h"
#include "TPolyLine.h"
#include "TClonesArray.h"
#include "TObject.h"
#include "TClass.h"

#include "Data.h"

HCircle  ThreePointCircle(TVector3 init, TVector3 middle, TVector3 final);//return z,x,r
HCircle  ThreePointCircle(std::vector<double> xarr, std::vector<double> zarr);//return z,x,r
TVector3 DyDtheta(TVector3 init, TVector3 middle, TVector3 final);//return Dy, Dtheta, RL
HHelix  GenerateSprial( TVector3 init, TVector3 middle, TVector3 final );
TPolyMarker3D* GenerateHitView( TClonesArray* arr );
TPolyMarker3D* GenerateClusterView( TClonesArray* arr );
TPolyMarker3D* GenerateClusterView( std::vector<TPCCluster> arr);
Double_t CalculateDeltaR( HHelix sprial, TVector3 vec );
Double_t CalculateDeltaY( HHelix sprial, TVector3 vec );
Bool_t   CalculateCrossing( HHelix sprial0, HHelix sprial1, Double_t& dist, TVector3& pos );
Double_t CalculateY(HHelix sprial, TVector3 vec );
Bool_t   CalculateCircleTangent( HHelix sprial0, TVector3 pos, TVector3& vec);
Bool_t   CalculateCircleCrossing( HHelix sprial0, HHelix sprial1, TVector3& vec0, TVector3& vec1);
Double_t CalculateMomentumR( HHelix sprial, double dtesla =1.);
Double_t CalculateMomentumY( HHelix sprial, double dtesla = 1. );
Double_t CalculateMomentum( HHelix sprial, double dtesla = 1. );

TPolyMarker3D* GenerateCrossingPoints( std::vector<HHelix> PSprial, std::vector<HHelix> MSprial );

class HyperCalculator : public TObject{
 public:
  HyperCalculator();
  virtual ~HyperCalculator();
  // List of calculation funciton

  static TVector3 CalculateMomentum( HHelix sprial, TVector3 pos, Double_t dtesla = 1.);
  static Bool_t   CCCrossing( HHelix sprial0, HHelix sprial1, TVector3& vec, double& dist );//// Calculate circle - circle crossing
  static Bool_t   CLCrossing( HHelix sprial, HLine line, TVector3& vec, double& dist );//// Calculate circle - line crossing
  static Bool_t   LLCrossing( HLine line0, HLine line1, TVector3& vec, double& dist );//// Calculate line -line crossing

  
 public:
  ClassDef( HyperCalculator, 0 )
};

#endif
