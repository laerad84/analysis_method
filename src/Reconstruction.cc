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
  std::cout<< "pim : " << pim.track.helix.R << std::endl;
  pim.track.Momentum = HyperCalculator::CalculateMomentum( track->helix, track->helix.InitPos, 1);
  pim.RecMass  = gPiM->Mass()*1000.;
  Double_t energy = TMath::Sqrt(pim.RecMass*pim.RecMass + pim.track.Momentum.Mag2());
  pim.Momentum.SetPxPyPzE( pim.track.Momentum.X(),
			   pim.track.Momentum.Y(),
			   pim.track.Momentum.Z(),
			   energy);
  pim.Position = track->helix.InitPos;
  pim.Energy = energy;
  return pim;
}
HProton     HRec::RecProton( TPCTrack* track ){
  HProton proton;
  proton.track          = *track;
  proton.track.Momentum = HyperCalculator::CalculateMomentum( track->helix, track->helix.InitPos, 1);
  proton.RecMass  = gProton->Mass()*1000.;
  Double_t energy = TMath::Sqrt(proton.RecMass*proton.RecMass + proton.track.Momentum.Mag2());
  proton.Momentum.SetPxPyPzE( proton.track.Momentum.X(),
			      proton.track.Momentum.Y(),
			      proton.track.Momentum.Z(),
			      energy);
  proton.Position = track->helix.InitPos;
  proton.Energy   = energy;
  return proton;
}
HKaonP      HRec::RecKaonP( TPCTrack* track ){
  HKaonP kaonp;
  kaonp.track   = *track;
  kaonp.track.Momentum = HyperCalculator::CalculateMomentum( track->helix, track->helix.InitPos, 1 );
  kaonp.RecMass = gProton->Mass()*1000.;
  Double_t energy = TMath::Sqrt(kaonp.RecMass*kaonp.RecMass + kaonp.track.Momentum.Mag2());
  kaonp.Momentum.SetPxPyPzE( kaonp.track.Momentum.X(),
			     kaonp.track.Momentum.Y(),
			     kaonp.track.Momentum.Z(), 
			     energy);
  kaonp.Position = track->helix.InitPos;
  kaonp.Energy   = energy;
  return kaonp;
}
HLambda     HRec::RecLambda( HProton p, HPiM pi){
  HLambda lambda;
  lambda.RecMass = gLambda->Mass()*1000.;
  //// lambda reconstruction code
  //// Circle - circle fitting
  Double_t dist= 0;
  TVector3 recPos;
  //std::cout<< "Calculate Crossing " << std::endl;
  //std::cout<< p.track.helix.R << "\t" <<  p.track.helix.X << "\t" << p.track.helix.Z << std::endl;
  HyperCalculator::CCCrossing( p.track.helix, pi.track.helix, recPos, dist );
  std::cout<< "Lambda rec position : " << std::endl;
  std::cout<< "Dist : " << dist << std::endl;
  recPos.Print(); 
  //std::cout<< "Calculate Momentum " << std::endl;
  TVector3 momentumP = HyperCalculator::CalculateMomentum( p.track.helix, recPos );
  TVector3 momentumPi= HyperCalculator::CalculateMomentum( pi.track.helix, recPos );

  Double_t       ep  = TMath::Sqrt( TMath::Power(gProton->Mass()*1000,2) + momentumP.Mag2());
  Double_t       epi = TMath::Sqrt( TMath::Power(gPiM->Mass()*1000,2) + momentumPi.Mag2());

  std::cout << "Energy : " << ep << "\t" << epi << std::endl;
  TLorentzVector lvp; 
  lvp.SetPxPyPzE(momentumP.X(), momentumP.Y(), momentumP.Z(),ep);
  TLorentzVector lvpi;
  lvpi.SetPxPyPzE(momentumPi.X(), momentumPi.Y(), momentumPi.Z(), epi );
  p.Momentum         = lvp;
  pi.Momentum        = lvpi;
  p.Energy           = lvp.E();
  pi.Energy          = lvpi.E();
  //// calculate momentum  
  lambda.Energy   = p.Energy + pi.Energy;
  lambda.Momentum = p.Momentum + pi.Momentum;
  lambda.RecMass  = lambda.Momentum.M();
  lambda.Position = recPos;
  lambda.Dist     = TMath::Abs(dist);
  lambda.proton   = new HProton(p);
  lambda.pion     = new HPiM(pi);
  return lambda;
}
HCascade    HRec::RecCascade( HLambda l, HPiM pi ){
  HCascade cascade;  
  cascade.RecMass = gCascade->Mass()*1000.;
  cascade.Energy  = l.Energy + pi.Energy;
  cascade.Momentum= l.Momentum + pi.Momentum;

  //// cascade reconstruction code
  //// circle - line fitting

  return cascade;
}
HDibaryonLL HRec::RecHDibaryonLL( HLambda l0, HLambda l1){
  HDibaryonLL HDLL;
  //HDLL.RecMass = gHDLL->Mass()*1000.;
  //// HDLL reconstruction code
  //// line line fitting
  HLine line0 = HyperCalculator::GenerateLine( l0.Momentum, l0.Position );
  HLine line1 = HyperCalculator::GenerateLine( l1.Momentum, l1.Position );
  Double_t dist = 0;
  TVector3 recPos;
  HyperCalculator::LLCrossing( line0, line1 , recPos, dist );  
  HDLL.Momentum  = l0.Momentum + l1.Momentum;
  HDLL.RecMass   = HDLL.Momentum.M();
  HDLL.Energy    = l0.Momentum.E() + l1.Momentum.E();
  HDLL.Position  = recPos;
  HDLL.Dist      = TMath::Abs(dist);
  HDLL.lambda0   = new HLambda( l0 );
  HDLL.lambda1   = new HLambda( l1 );
  return HDLL;
}

HDibaryonPC HRec::RecHDibaryonPC( HCascade cc, HProton p ){
  HDibaryonPC HDPC;
  HDPC.RecMass = gHDPC->Mass()*1000.;
  //// HDPC reconstruction code;
  //// circle - line fitting

  return HDPC;
}
