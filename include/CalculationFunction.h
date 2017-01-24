#ifndef CalculationFunction__H__
#define CalculationFunction__H__
#include <vector>
#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "TMath.h"
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
HSprial  GenerateSprial( TVector3 init, TVector3 middle, TVector3 final );
TPolyMarker3D* GenerateHitView( TClonesArray* arr );
TPolyMarker3D* GenerateClusterView( TClonesArray* arr );
TPolyMarker3D* GenerateClusterView( std::vector<TPCCluster> arr);
Double_t CalculateDeltaR( HSprial sprial, TVector3 vec );
Double_t CalculateDeltaY( HSprial sprial, TVector3 vec );
Bool_t   CalculateCrossing( HSprial sprial0, HSprial sprial1, Double_t& dist, TVector3& pos );
Double_t CalculateY(HSprial sprial, TVector3 vec );
Bool_t   CalculateCircleTangent( HSprial sprial0, TVector3 pos, TVector3& vec);
Bool_t   CalculateCircleCrossing( HSprial sprial0, HSprial sprial1, TVector3& vec0, TVector3& vec1);
Double_t CalculateMomentumR( HSprial sprial, double dtesla =1.);
Double_t CalculateMomentumY( HSprial sprial, double dtesla = 1. );
Double_t CalculateMomentum( HSprial sprial, double dtesla = 1. );

TPolyMarker3D* GenerateCrossingPoints( std::vector<HSprial> PSprial, std::vector<HSprial> MSprial );

#endif
