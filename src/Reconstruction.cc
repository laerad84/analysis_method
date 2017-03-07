#include "Reconstruction.h"

ClassImp( HRec )

HRec::HRec(){
  ;
}
HRec::~HRec(){
  ;
}
HPiM        HRec::RecPiM( TPCTrack* track ){
  HPiM pim;
  pim.track    = *track;
  pim.RecMass  = gPiM->Mass()*0.001;
  Double_t energy = TMath::Sqrt(pim.RecMass*pim.RecMass + track->Momentum.Mag2());
  pim.Momentum.SetPxPyPzE( track->Momentum.X(), track->Momentum.Y(), track->Momentum.Z(), energy);
  pim.Position = track->helix.InitPos;
  pim.Energy = energy;
  return pim;
}
HProton     HRec::RecProton( TPCTrack* track ){
  HProton proton;
  proton.track    = *track;
  proton.RecMass  = gProton->Mass()*0.001;
  Double_t energy = TMath::Sqrt(proton.RecMass*proton.RecMass + track->Momentum.Mag2());
  proton.Momentum.SetPxPyPzE( track->Momentum.X(), track->Momentum.Y(), track->Momentum.Z(), energy);
  proton.Position = track->helix.InitPos;
  proton.Energy   = energy;
  return proton;
}
HKaonP      HRec::RecKaonP( TPCTrack* track ){
  HKaonP kaonp;
  kaonp.track   = *track;
  kaonp.RecMass = gProton->Mass()*0.001;
  Double_t energy = TMath::Sqrt(kaonp.RecMass*kaonp.RecMass + track->Momentum.Mag2());
  kaonp.Momentum.SetPxPyPzE( track->Momentum.X(), track->Momentum.Y(), track->Momentum.Z(), energy);
  kaonp.Position = track->helix.InitPos;
  kaonp.Energy   = energy;
  return kaonp;
}
HLambda     HRec::RecLambda( HProton p, HPiM pi){
  HLambda lambda;
  lambda.RecMass = gLambda->Mass()*0.001;
  //// lambda reconstruction code
  //// Circle - circle fitting
  Double_t dist;
  TVector3 recPos;
  HyperCalculator::CCCrossing( p.track.helix, pi.track.helix, recPos, dist );
  //// calculate momentum
  lambda.Energy = p.Energy + pi.Energy;

  return lambda;
}
HCascade    HRec::RecCascade( HLambda l, HPiM pi ){
  HCascade cascade;
  cascade.RecMass = gCascade->Mass()*0.001;
  //// cascade reconstruction code
  //// circle - line fitting

  return cascade;

}
HDibaryonLL HRec::RecHDibaryonLL( HLambda l0, HLambda l1){
  HDibaryonLL HDLL;
  HDLL.RecMass = gHDLL->Mass()*0.001;
  //// HDLL reconstruction code
  //// line line fitting

  return HDLL;
}
HDibaryonPC HRec::RecHDibaryonPC( HCascade cc, HProton p ){
  HDibaryonPC HDPC;
  HDPC.RecMass = gHDPC->Mass()*0.001;
  //// HDPC reconstruction code;
  //// circle - line fitting

  return HDPC;
}
