#ifndef FITTER__H__
#define FITTER__H__

#include <iostream>
#include <vector>
#include "TObject.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"

#include "Data.h"
#include "CalculationFunction.h"
//#include "ConverterFunction.h"

class HelixFitter : public TObject {
 public :

  HelixFitter();
  ~HelixFitter();
  static Double_t RFitCut;//limits on R diff
  static Double_t YFitCut;//linits on Y diff
  HHelix  FitResult;
  static Bool_t  bWeight;
  static Bool_t  bEnergy;
  void FitData( TPCTrack* track );
  void Fit(TGraph2D* track, TGraph2D* weights = NULL );
  void CircleFit();
  void HelixFit();
  static TGraph2D*     fitTrack;
  static TGraph2D*     fitWeight;
  static TGraphErrors* fitYTheta;
  static std::vector<Double_t> EArr;
  static void          Chi2MinimizerCircle( Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t );
  //static void          Chi2MinimizerHelix( Int_t & Double_t *, Double_t &f, Double_t *par, Int_t );
  HHelix               GetResult() const { return FitResult;}

  ClassDef( HelixFitter, 0 )
};

R__EXTERN HelixFitter* gHelixFitter;

#endif //FITTER__H__
