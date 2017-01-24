#ifndef TPCGlobals__H__
#define TPCGlobals__H__
#include "TObject.h"
class TPCGlobals {
 public :
  TPCGlobals(){;};
  virtual ~TPCGlobals(){;};
  static const Double_t sTPC_Pad_Parameter[32][6];
  static const Int_t    sTPC_PadRowMaxID[32];
  static const Int_t    sTPC_NMaximum;
  static const Int_t    sTPC_NEperEnergy;
  static const Double_t sTPC_zOffset;
  static const Int_t    sTPC_GEM_DefusionR;

  ClassDef( TPCGlobals, 1 )
};


#endif //TPCGlobals__H__
