#ifndef RECONSTRUCTION__H__
#define RECONSTRUCTION__H__
#include "Data.h"
#include "CalculationFunction.h"
#include "Particles.h"
#include "TObject.h"
#include "TDatabasePDG.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class HRec : public TObject {
 public:
  HRec();
  virtual ~HRec();
  //// reconstruct basic particles from track
  static HPiM        RecPiM( TPCTrack* track );
  static HProton     RecProton( TPCTrack* track );
  static HKaonP      RecKaonP( TPCTrack* track );

  //// reconstruct decay particles from basic track;
  static HLambda     RecLambda(HProton p, HPiM pi);
  static HCascade    RecCascade( HLambda l, HPiM pi );
  static HDibaryonLL RecHDibaryonLL( HLambda l0, HLambda l1);
  static HDibaryonPC RecHDibaryonPC( HCascade cc, HProton p );


 public:
  ClassDef( HRec, 0)
};



#endif //RECONSTRUCTION
