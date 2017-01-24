#ifndef TPCFitter__H__
#define TPCFitter__H__
#include "TGraph.h"

#include "AbsFinitePlane.h"
#include "DetPlane.h"
#include "RectangularFinitePlane.h"
#include "ReferenceStateOnPlane.h"
#include "StateOnPlane.h"
#include "SharedPlanePtr.h"

#include "MeasuredStateOnPlane.h"
#include "MeasurementOnPlane.h"
#include "FullMeasurement.h"
#include "PlanarMeasurement.h"
#include "ProlateSpacepointMeasurement.h"
#include "SpacepointMeasurement.h"
#include "AbsMeasurement.h"
#include "muSpacepointMeasurement.h"

#include "AbsFitterInfo.h"
#include "AbsTrackRep.h"

#include "KalmanFitter.h"
#include "KalmanFitterRefTrack.h"
#include "KalmanFitterInfo.h"
#include "KalmanFitStatus.h"
#include "AbsKamanFitter.h"

#include "Tools.h"
#include "TrackCand.h"
#include "TrackCandHit.h"
#include "RKTools.h"
#include "RKTrackRep.h"
#include "StepLimits.h"

#include "FieldManager.h"


class TPCCircleFitter : public TObject {
 public:


  static TGraph* gr;
  //TPCCircleFitter();
  //~TPCCircleFitter();

  //void Init();


 private:

 public:
  ClassDef( TPCCircleFitter, 0 )
}

#endif
