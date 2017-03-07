#ifndef RECONSTRUCTION__H__
#define RECONSTRUCTIONI__H__
#include "Data.h"
#include "CalculationFunction.h"
#include "Particles.h"
#include "TObject.h"
#include "TDatabasePDG.h"
class HRec : public TObject {
 public:
  HRec();
  virtual ~HRec();
  //// reconstruct basic particles from track
  HPiM        RecPiM( TPCTrack* track );
  HProton     RecProton( TPCTrack* track );
  HKaonP      RecKaonP( TPCTrack* track );


  //// reconstruct decay particles from basic track;
  HLambda     RecLambda(HProton p, HPiM pi);
  HCascade    RecCascade( HLambda l, HPiM pi );
  HDibaryonLL RecHDibaryonLL( HLambda l0, HLambda l1);
  HDibaryonPC RecHDibaryonPC( HCascade cc, HProton p );


 public:
  ClassDef( HRec, 0)
};



#endif //RECONSTRUCTION
