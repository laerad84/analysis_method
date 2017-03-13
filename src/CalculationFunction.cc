#include "CalculationFunction.h"

HCircle ThreePointCircle(std::vector<double> xarr, std::vector<double> zarr){
  Int_t iIndex = 0;
  Int_t mIndex = (int)((xarr.size() -1)/2.);
  Int_t fIndex = xarr.size()-1;
  if( xarr.size() > 3 ){
    fIndex = xarr.size() -2;
  }
  if( xarr.size() > 5 ){
    iIndex = 1;
  }
  std::vector<Int_t> IndexArr;
  IndexArr.push_back(iIndex);
  IndexArr.push_back(mIndex);
  IndexArr.push_back(fIndex);

  double r,x,z;

  Double_t A = (xarr.at(iIndex)*xarr.at(iIndex) - xarr.at(mIndex)*xarr.at(mIndex) + zarr.at(iIndex)*zarr.at(iIndex) - zarr.at(mIndex)*zarr.at(mIndex))/2.;
  Double_t B = (xarr.at(mIndex)*xarr.at(mIndex) - xarr.at(fIndex)*xarr.at(fIndex) + zarr.at(mIndex)*zarr.at(mIndex) - zarr.at(fIndex)*zarr.at(fIndex))/2.;

  x = ((zarr.at(mIndex) - zarr.at(iIndex))*B - (zarr.at(fIndex) - zarr.at(mIndex))*A) /
    ((zarr.at(fIndex) - zarr.at(mIndex))*(xarr.at(mIndex) - xarr.at(iIndex)) - (xarr.at(fIndex) - xarr.at(mIndex))*(zarr.at(mIndex) - zarr.at(iIndex)));
  z = ((xarr.at(mIndex) - xarr.at(iIndex))*B - (xarr.at(fIndex) - xarr.at(mIndex))*A) /
    ((xarr.at(fIndex) - xarr.at(mIndex))*(zarr.at(mIndex) - zarr.at(iIndex)) - (zarr.at(fIndex) - zarr.at(mIndex))*(xarr.at(mIndex) - xarr.at(iIndex)));

  double rsq = 0;
  for( int i = 0; i< 3; i++){
    std::cout<< (xarr.at(IndexArr.at(i)) - x)*(xarr.at(IndexArr.at(i)) - x) + (zarr.at(IndexArr.at(i)) -z)*(zarr.at(IndexArr.at(i)) - z) << std::endl;
    rsq = rsq + (xarr.at(IndexArr.at(i)) - x)*(xarr.at(IndexArr.at(i)) - x) + (zarr.at(IndexArr.at(i)) -z)*(zarr.at(IndexArr.at(i)) - z);
  }
  rsq = rsq/3.;
  r = TMath::Sqrt(rsq);
  //TVector3 result(z,x,r);
  HCircle result(-1,x,z,r);
  return result;
}

HCircle ThreePointCircle( TVector3 init, TVector3 middle, TVector3 final){
  Int_t iIndex = 0;
  Int_t mIndex = 1;
  Int_t fIndex = 2;
  std::vector<double> xarr;
  std::vector<double> zarr;
  std::vector<Int_t> IndexArr;
  IndexArr.push_back(0);
  IndexArr.push_back(1);
  IndexArr.push_back(2);
  xarr.push_back( init.X());
  xarr.push_back( middle.X());
  xarr.push_back( final.X());
  zarr.push_back( init.Z());
  zarr.push_back( middle.Z());
  zarr.push_back( final.Z());

  double r,x,z;
  Double_t A = (xarr.at(iIndex)*xarr.at(iIndex) - xarr.at(mIndex)*xarr.at(mIndex) + zarr.at(iIndex)*zarr.at(iIndex) - zarr.at(mIndex)*zarr.at(mIndex))/2.;
  Double_t B = (xarr.at(mIndex)*xarr.at(mIndex) - xarr.at(fIndex)*xarr.at(fIndex) + zarr.at(mIndex)*zarr.at(mIndex) - zarr.at(fIndex)*zarr.at(fIndex))/2.;

  x = ((zarr.at(mIndex) - zarr.at(iIndex))*B - (zarr.at(fIndex) - zarr.at(mIndex))*A) /
    ((zarr.at(fIndex) - zarr.at(mIndex))*(xarr.at(mIndex) - xarr.at(iIndex)) - (xarr.at(fIndex) - xarr.at(mIndex))*(zarr.at(mIndex) - zarr.at(iIndex)));
  z = ((xarr.at(mIndex) - xarr.at(iIndex))*B - (xarr.at(fIndex) - xarr.at(mIndex))*A) /
    ((xarr.at(fIndex) - xarr.at(mIndex))*(zarr.at(mIndex) - zarr.at(iIndex)) - (zarr.at(fIndex) - zarr.at(mIndex))*(xarr.at(mIndex) - xarr.at(iIndex)));

  double rsq = 0;
  for( int i = 0; i< 3; i++){
    std::cout<< (xarr.at(IndexArr.at(i)) - x)*(xarr.at(IndexArr.at(i)) - x) + (zarr.at(IndexArr.at(i)) -z)*(zarr.at(IndexArr.at(i)) - z) << std::endl;
    rsq = rsq + (xarr.at(IndexArr.at(i)) - x)*(xarr.at(IndexArr.at(i)) - x) + (zarr.at(IndexArr.at(i)) -z)*(zarr.at(IndexArr.at(i)) - z);
  }
  rsq = rsq/3.;
  r = TMath::Sqrt(rsq);
  //TVector3 result(z,x,r);
  HCircle  result(-1,x,z,r);
  return result;
}

TVector3 DyDtheta(TVector3 init, TVector3 middle, TVector3 final ){
  HCircle  circle = ThreePointCircle(init,middle,final);

  Double_t cx = circle.X;
  Double_t cz = circle.Z;
  Double_t cr = circle.R;
  Double_t initTheta   = TMath::ATan2(init.X()-cx,init.Z()-cz);
  Double_t middleTheta = TMath::ATan2(middle.X()-cx,middle.Z()-cz);
  Double_t finalTheta  = TMath::ATan2(final.X()-cx,final.Z()-cz);
  Double_t centerTheta = TMath::ATan2(cx, cz);

  Double_t Dy = final.Y() - init.Y();
  Double_t Dtheta(0);
  Double_t RL=0;

  //Calculate Dl, Dtheta, RL
  if( (initTheta < middleTheta && middleTheta < finalTheta ) ||
      (finalTheta < initTheta  && initTheta < middleTheta  ) ||
      (middleTheta < finalTheta && finalTheta < initTheta )){
    //Normal order
    RL = 1;
    if( initTheta < finalTheta ){
      Dtheta = finalTheta - initTheta;
    }else{
      Dtheta = 2*TMath::Pi() + (finalTheta - initTheta);
    }
  }else{
    //Inverse order
    RL = -1;
    if( initTheta > finalTheta ){
      Dtheta = finalTheta - initTheta;
    }else{
      Dtheta = -2*TMath::Pi() + (finalTheta - initTheta );
    }
  }

  TVector3 result(Dy, Dtheta, RL);
  return result;
}
HHelix GenerateSprial( TVector3 init, TVector3 middle, TVector3 final){

  HCircle circle = ThreePointCircle( init, middle, final);
  TVector3 parameter = DyDtheta( init, middle, final );
  HHelix sprial_result( circle.ID, circle.X, (init.Y()+final.Y())/2., circle.Z, circle.R, parameter.X(),parameter.Y(),(int)(parameter.Z()));
  /// initial Y should be recalculated ///
  sprial_result.InitPos = init;
  sprial_result.FinalPos = final;
  std::cout<< sprial_result.ID << "\t" << sprial_result.X << std::endl;
  return sprial_result;

}
TPolyMarker3D* GenerateHitView( TClonesArray* arr ){
  TClass* cl = arr->GetClass();
  if( strcmp( "TPCHit", cl->GetName()) != 0 ){
    return NULL;
  }
  TPolyMarker3D* marker = new TPolyMarker3D();
  for( int iarr = 0; iarr < arr->GetEntries(); iarr++){
    TPCHit* hit = (TPCHit*)arr->At(iarr);
    marker->SetPoint( iarr,
		      hit->Position().X(),
		      hit->Position().Y(),
		      hit->Position().Z());
  }
  return marker;
}
TPolyMarker3D* GenerateClusterView( TClonesArray* arr ){
  TClass* cl  = arr->GetClass();
  if( strcmp( "TPCCluster", cl->GetName() ) != 0){
    return NULL;
  }
  TPolyMarker3D* marker = new TPolyMarker3D();
  for( int iarr = 0; iarr < arr->GetEntries(); iarr++){
    TPCCluster* cl = (TPCCluster*)arr->At(iarr);
    marker->SetPoint(iarr,cl->Position.X(),cl->Position.Y(),cl->Position.Z());
  }
  return marker;
}
TPolyMarker3D* GenerateClusterView( std::vector<TPCCluster> arr){
  TPolyMarker3D* marker = new TPolyMarker3D();
  for( int iarr = 0; iarr < arr.size(); iarr++){
    marker->SetPoint(iarr,
		     arr.at(iarr).Position.X(),
		     arr.at(iarr).Position.Y(),
		     arr.at(iarr).Position.Z());
  }
  return marker;
}
Double_t CalculateDeltaR( HHelix sprial, TVector3 vec){
  Double_t dx = vec.X() - sprial.X;
  Double_t dz = vec.Z() - sprial.Z;
  Double_t dr = TMath::Sqrt( dx*dx + dz*dz ) - sprial.R;
  return dr;
}
Double_t CalculateDeltaY( HHelix sprial, TVector3 vec ){
  Double_t initTheta = TMath::ATan2(sprial.InitPos.X()-sprial.X,sprial.InitPos.Z()-sprial.Z);
  Double_t theta     = TMath::ATan2( vec.X() - sprial.X, vec.Z() - sprial.Z);
  Double_t dTheta    = 0;
  dTheta  = theta - initTheta;
  if( dTheta >= 2*TMath::Pi() ){ dTheta = dTheta - 2*TMath::Pi();}
  if( dTheta <= -2*TMath::Pi()){ dTheta = dTheta - 2*TMath::Pi();}

  Double_t dy    = vec.Y() - dTheta*sprial.DY/sprial.DTheta - sprial.Y;
  return dy;
}
Double_t CalculateY( HHelix sprial, TVector3 vec ){
  Double_t initTheta = TMath::ATan2(sprial.InitPos.X()-sprial.X,sprial.InitPos.Z()-sprial.Z);
  Double_t theta     = TMath::ATan2( vec.X() - sprial.X, vec.Z() - sprial.Z);
  Double_t dTheta    = theta - initTheta;
  //if( dTheta >= 2*TMath::Pi() ){ dTheta = dTheta - 2*TMath::Pi();}
  //if( dTheta <= -2*TMath::Pi()){ dTheta = dTheta - 2*TMath::Pi();}
  Double_t y    = dTheta*sprial.DY/sprial.DTheta + sprial.InitPos.Y();
  return y;
}
Bool_t CalculateCircleCrossing( HHelix sprial0, HHelix sprial1, TVector3& vec0, TVector3& vec1 ){
  Double_t r0 = sprial0.R;
  Double_t r1 = sprial1.R;
  TVector2 v0(sprial0.X,sprial1.Z);
  TVector2 v1(sprial1.X,sprial1.Z);
  Double_t dist = (v1 - v0).Mod();

  if( dist > r0 +r1 ||
      dist < TMath::Abs(r0 -r1 )){
    Double_t y0 = CalculateY( sprial0, TVector3( sprial1.X,0,sprial1.Z));
    Double_t y1 = CalculateY( sprial1, TVector3( sprial0.X,0,sprial0.Z));
    TVector3 pVec(sprial0.X - sprial1.X, 0, sprial0.Z - sprial1.Z);
    TVector3 nVec = pVec.Unit();
    TVector3 vv0 = TVector3( sprial0.X, y0, sprial0.Z ) + sprial0.R*nVec;
    TVector3 vv1 = TVector3( sprial1.X, y1, sprial1.Z ) - sprial1.R*nVec;
    vec0 = vv0;
    vec1 = vv1;
    return false;
  }
  Double_t phi = (v1-v0).Phi();
  TVector2 unit = (v1-v0).Unit();
  double  a = (v1-v0).X();
  double  b = (v1-v0).Y();
  double  c = ( v1.Mod2() - v0.Mod2() + r0*r0 -r1*r1)/2.;

  double distPoint = TMath::Abs( a*v0.X() + b*v0.Y() -c )/TMath::Sqrt(a*a +b*b);
  double crossingLength =  TMath::Sqrt(r0*r0-distPoint*distPoint);

  double dTheta    = TMath::ACos( distPoint/r0);
  double Phi       = unit.Phi();
  TVector2 vecCent;
  vecCent.SetMagPhi(distPoint, Phi);
  TVector2 targetCenter( 0,-143);
  TVector3 vecCentral( vecCent.X(), 0, vecCent.Y() );
  TVector3 vecCrossing( vecCent.X(),0, vecCent.Y());
  vecCrossing.RotateY(90*TMath::DegToRad());
  vecCrossing.SetMag(crossingLength);
  TVector3 vecOffset(sprial0.X,0,sprial0.Z);
  vec0 = vecCentral + vecCrossing + vecOffset;
  vec1 = vecCentral - vecCrossing + vecOffset;

  return true;
}
Bool_t CalculateCircleTangent( HHelix sprial0, TVector3 pos, TVector3& vec ){
  Double_t r0 = sprial0.R;
  TVector3 vecCenter = TVector3(pos.X()-sprial0.X, 0,pos.Z() - sprial0.Z);
  if( sprial0.RL > 0 ){
    vecCenter.RotateY(90*TMath::DegToRad());
  }else{
    vecCenter.RotateY(-90*TMath::DegToRad());
  }
  vec = vecCenter.Unit();
  vec.SetY( sprial0.DY/sprial0.DTheta/sprial0.R );
  return true;
}

Bool_t CalculateCrossing( HHelix sprial0, HHelix sprial1, Double_t& dist, TVector3& pos ){
  TVector3 v0;
  TVector3 v1;
  if( !CalculateCircleCrossing( sprial0, sprial1, v0,v1) ){// not connected
    dist = (v0 - v1).Mag();
    pos = TVector3( (v0.X() + v1.X())/2,(v0.Y() + v1.Y())/2,(v0.Z() + v1.Z())/2 );
    return false;
  }
  Double_t y00 = CalculateY( sprial0, v0);
  Double_t y01 = CalculateY( sprial0, v1);
  Double_t y10 = CalculateY( sprial1, v0);
  Double_t y11 = CalculateY( sprial1, v1);
  Double_t yCenter= 0;
  std::cout<< "Y : " << y00 << "\t" << y01 <<"\t" << y10 << "\t" << y11 << std::endl;
  if( TMath::Abs( y00 - y10 ) < TMath::Abs( y01 - y11)){
    yCenter = (y00 + y10 )/2.;
    dist = TMath::Abs(y00-y10);
    pos = TVector3( v0.X(), yCenter, v0.Z());
  }else{
    yCenter = (y01 + y11 )/2.;
    dist = TMath::Abs( y01 - y11 );
    pos = TVector3( v1.X(), yCenter, v1.Z());
  }
  return true;
}
double CalculateMomentumR( HHelix sprial, double dtesla){
  //std::cout<< __PRETTY_FUNCTION__ << std::endl;
  return sprial.R * dtesla * 0.3;
}
double CalculateMomentumY( HHelix sprial, double dtesla){
  //std::cout<< __PRETTY_FUNCTION__ << std::endl;
  double momentumR = CalculateMomentumR( sprial, dtesla );
  return momentumR * sprial.DY / TMath::Abs( sprial.R * sprial.DTheta );
}

double CalculateMomentum( HHelix sprial, double dtesla ){
  //std::cout<< __PRETTY_FUNCTION__ << std::endl;
  double pR = CalculateMomentumR( sprial, dtesla );
  double pY = CalculateMomentumY( sprial, dtesla );
  double totalMomentum = TMath::Sqrt( pR*pR + pY*pY );
  return totalMomentum;
}

TPolyMarker3D* GenerateCrossingPoints( std::vector<HHelix> PSprial, std::vector<HHelix> MSprial ){
  /// Calculation crossing point combination between particle +/-
  TPolyMarker3D* marker = new TPolyMarker3D();
  for( int ip = 0; ip < PSprial.size(); ip++){
    for( int im = 0; im < MSprial.size(); im++){
      Double_t dist = 0;
      TVector3 pos(0,0,0);
      Bool_t cross = CalculateCrossing( PSprial.at(ip), MSprial.at(im), dist, pos);
      if( !cross ){ continue; }
      if( dist > 40 ){ continue; }
      marker->SetNextPoint( pos.X(), pos.Y(), pos.Z());
    }
  }
  return marker;
}

ClassImp( HyperCalculator )
HyperCalculator::HyperCalculator(){
  ;
}
HyperCalculator::~HyperCalculator(){
  ;
}
HLine HyperCalculator::GenerateLine( TLorentzVector mom, TVector3 pos ){
  HLine line;
  TVector3 vec( mom.Px(), mom.Py(), mom.Pz());
  line.Direction = vec.Unit();
  line.Offset    = pos;
  line.bNorm     = kTRUE;
  return line;
}
HLine HyperCalculator::GenerateLine( TVector3 mom, TVector3 pos ){
  HLine line;
  line.Direction = mom.Unit();
  line.Offset    = pos;
  line.bNorm     = kTRUE;
  return line;  
}


TVector3 HyperCalculator::CalculateMomentum( HHelix sprial, TVector3 pos, Double_t dtesla ){
  TVector3 momentum(0,0,0);
  // Calculate center to point direction
  TVector3 vecNorm( pos.X() - sprial.X,0, pos.Z() - sprial.Z);
  if( sprial.RL > 0 ){
    // RL > 0 : positive
    vecNorm.RotateY( TMath::DegToRad() * 90 );
  }else if( sprial.RL < 0 ){
    // RL < 0 : negative
    vecNorm.RotateY(-TMath::DegToRad() * 90 );
  }else{
    // bad sprial 
    return momentum;
  }
  // Calculate momentum
  Double_t px = vecNorm.Unit().X()*sprial.R*0.3*dtesla;//MeV/c
  Double_t pz = vecNorm.Unit().Z()*sprial.R*0.3*dtesla;//MeV/c
  Double_t py = sprial.DY/TMath::Abs(sprial.DTheta)*0.3*dtesla;//MeV/c
  momentum.SetXYZ( px, py, pz);
  return momentum; 
}

Bool_t HyperCalculator::CCCrossing( HHelix sprial0, HHelix sprial1, TVector3& vec, double& dist ){
  Double_t centerDist = TMath::Sqrt(TMath::Power( sprial0.X - sprial1.X, 2)+
				    TMath::Power( sprial0.Z - sprial1.Z, 2));
  if( centerDist > sprial0.R + sprial1.R ||
      centerDist < TMath::Abs( sprial0.R - sprial1.R )){
    std::cout<< "Notconnected : "<< centerDist << "\t" << sprial0.R << "\t" << sprial1.R  << std::endl;

    /// Not connected 
    /// Calculate closet point 
    TVector3 c0( sprial0.X, 0, sprial0.Z);
    TVector3 c1( sprial1.X, 0, sprial1.Z); 
    TVector3 vD =( c1 - c0 ).Unit();
    TVector3 cc0 = c0 + sprial0.R*vD;
    TVector3 cc1 = c1 - sprial1.R*vD;
    
    /// Calculate y0;
    Double_t iTheta0 = TMath::ATan2( sprial0.InitPos.X() - sprial0.X,
				     sprial0.InitPos.Z() - sprial0.Z);
    Double_t theta0  = TMath::ATan2( cc0.X() - sprial0.X,
				     cc0.Z() - sprial0.Z);
    Double_t y0      = (theta0 - iTheta0)*sprial0.DY/sprial0.DTheta - sprial0.InitPos.Y();

    /// Calculate y1;
    Double_t iTheta1 = TMath::ATan2( sprial1.InitPos.X() - sprial1.X,
				     sprial1.InitPos.Z() - sprial1.Z);
    Double_t theta1  = TMath::ATan2( cc1.X() - sprial1.X,
				     cc1.Z() - sprial1.Z);
    Double_t y1      = (theta1 - iTheta1)*sprial1.DY/sprial1.DTheta - sprial1.InitPos.Y();

    cc0.SetY(y0);
    cc1.SetY(y1);
    
    TVector3 vMiddle( (cc0.X() + cc1.X())/2.,
		      (cc0.X() + cc1.Y())/2.,
		      (cc0.Z() + cc1.Z())/2.);
    dist = (cc0 - cc1).Mag();
    vec = vMiddle;
    return false; // not connected
  }else{ 
    /// Connected 
    std::cout<< "Connected" << std::endl;
    TVector3 c0( sprial0.X, 0, sprial0.Z);
    TVector3 c1( sprial1.X, 0, sprial1.Z);
    TVector3 RUnit = (c1 - c0).Unit();
    TVector3 aUnit = RUnit;
    aUnit.RotateY(90*TMath::DegToRad());

    Double_t d0 = (centerDist*centerDist + sprial0.R*sprial0.R - sprial1.R*sprial1.R)/2./centerDist;
    Double_t d1 = (centerDist*centerDist - sprial0.R*sprial0.R + sprial1.R*sprial1.R)/2./centerDist;
    Double_t a  = 0;
    TVector3 cc0(0,0,0);
    TVector3 cc1(0,0,0);
    if( sprial0.R > sprial1.R ){
      a   = TMath::Sqrt( sprial0.R*sprial0.R - d0*d0 );
      cc0 = c0 + d0*RUnit + a*aUnit;
      cc1 = c0 + d0*RUnit - a*aUnit;
    }else{
      a   = TMath::Sqrt( sprial1.R*sprial1.R - d1*d1 );
      cc0 = c1 - d1*RUnit + a*aUnit;
      cc1 = c1 - d1*RUnit - a*aUnit;
    }

    Double_t y00 = CalculateY( sprial0, cc0 );
    Double_t y01 = CalculateY( sprial0, cc1 );
    Double_t y10 = CalculateY( sprial1, cc0 );
    Double_t y11 = CalculateY( sprial1, cc1 );
    Double_t yC  = 0;
    if( TMath::Abs( y00 - y10 ) < TMath::Abs( y01 - y11 )){
      // cc0 is cross point 
      yC   = (y00 + y10 )/2.;
      dist = TMath::Abs( y00 - y10 );
      vec.SetXYZ(cc0.X(),yC,cc0.Z());
    }else{
      // cc1 is cross point 
      yC   = (y10 + y11);
      dist = TMath::Abs( y10 - y11 );
      vec.SetXYZ(cc1.X(),yC,cc1.Z());
    }
    std::cout<< "Circle cross point" << std::endl;
    std::cout<< dist << std::endl;
    vec.Print();

    return true; // connected 
  }
}

Bool_t HyperCalculator::CLCrossing( HHelix sprial, HLine line, TVector3& vec, double& dist ){
  Double_t mag = TMath::Sqrt( line.Direction.Z()*line.Direction.Z() + line.Direction.X()*line.Direction.X() );
  TVector2 vec_cc( sprial.Z, sprial.X);
  TVector2 vec_direction(line.Direction.Z()/mag, line.Direction.X()/mag);
  TVector2 vec_normal(line.Direction.X()/mag, -line.Direction.Z()/mag);
  TVector2 vec_delta( sprial.Z - line.Offset.Z(), sprial.X - line.Offset.X());
  Double_t distance = vec_normal.X()*vec_delta.X() + vec_normal.Y()*vec_delta.Y();
  TVector2 vec_point( vec_normal.X()*distance, vec_normal.Y()*distance );
  TVector3 vec0;
  TVector3 vec1;
  TVector3 meanHitPosition;
  if( sprial.R >= TMath::Abs(distance) ){// crossed
    /// Calculate crossing point
    Double_t d  = TMath::Sqrt( sprial.R*sprial.R + distance*distance );
    vec0.SetX(vec_cc.X() - vec_point.X() + vec_direction.X()*d);
    vec0.SetY(vec_cc.Y() - vec_point.Y() + vec_direction.Y()*d);

    vec1.SetX(vec_cc.X() - vec_point.X() - vec_direction.X()*d);
    vec1.SetY(vec_cc.Y() - vec_point.Y() - vec_direction.Y()*d);

    /// Calculate helix Y
    Double_t y0  = CalculateY( sprial, TVector3( vec0.Y(), 0, vec0.Y()));
    Double_t y1  = CalculateY( sprial, TVector3( vec1.Y(), 0, vec1.Y()));
    Double_t yl0 = 0;
    Double_t yl1 = 0;
    /// Calculate line Y
    if( vec_direction.X() != 0 ){/// z component
      yl0 = ((vec0.X() - line.Offset.Z())/line.Direction.Z() * line.Direction.Y()) + line.Offset.Y();
      yl1 = ((vec1.X() - line.Offset.Z())/line.Direction.Z() * line.Direction.Y()) + line.Offset.Y();
    }else{/// x component
      yl0 = ((vec0.Y() - line.Offset.X())/line.Direction.X() * line.Direction.Y()) + line.Offset.Y();
      yl1 = ((vec1.Y() - line.Offset.X())/line.Direction.X() * line.Direction.Y()) + line.Offset.Y();
    }
    if( TMath::Abs( y0 - yl0 ) < TMath::Abs( y1 - yl1 )){
      vec.SetXYZ( vec0.Z(), (y0+yl0)/2., vec0.X());
      dist = y0-yl0;
      return true;
    }else{
      vec.SetXYZ( vec1.Z(), (y1+yl1)/2., vec1.X());
      dist = y1-yl1;
      return true;
    }
    return true;
  }else{ // not crossed
    vec0.SetX(vec_cc.X() -sprial.R*vec_normal.X() );// closest point on circle
    vec0.SetY(vec_cc.Y() -sprial.R*vec_normal.Y() );// closest point on circle
    vec1.SetX(vec_cc.X() -vec_point.X() );
    vec1.SetY(vec_cc.Y() -vec_point.Y() );

    Double_t y0 = CalculateY( sprial, TVector3( vec0.Y(), 0, vec0.X()));
    Double_t y1;
    if( vec_direction.X() != 0 ){/// z component
      y1 = ((vec1.X() - line.Offset.Z())/line.Direction.Z() * line.Direction.Y()) + line.Offset.Y();
    }else{/// x component
      y1 = ((vec1.Y() - line.Offset.X())/line.Direction.X() * line.Direction.Y()) + line.Offset.Y();
    }
    vec.SetXYZ( (vec0.Z()+vec1.Z())/2., (y0+y1)/2., (vec0.X()+vec0.X())/2. );
    TVector3 delta_line( vec0.Z() - vec1.Z(), y0-y1, vec0.X() - vec1.X());
    dist = delta_line.Mag();
    return false;
  }
}
Bool_t HyperCalculator::LLCrossing( HLine line0, HLine line1, TVector3& vec, double& dist ){

  // Get Unit vector of the line
  Double_t mag0 = line0.Direction.Mag();
  Double_t mag1 = line1.Direction.Mag(); 
  TVector3 v0(line0.Direction.X()/mag0,line0.Direction.Y()/mag0,line0.Direction.Z()/mag0);
  TVector3 v1(line1.Direction.X()/mag1,line1.Direction.Y()/mag1,line1.Direction.Z()/mag1);
  // Get shortest vector l0 -> l1 
  TVector3 vCross = (v0.Cross(v1)).Unit();
  Double_t d0     = line0.Offset.Dot(vCross);
  Double_t d1     = line1.Offset.Dot(vCross);
  TVector3 vn     = (d1 - d0)*vCross;

  // Calculate displacement of offsets
  TVector3 vl     = line1.Offset - line0.Offset - vn;
  TVector3 vls0   =  vl - ( vl.Dot( v0 ))*v0;
  TVector3 vls1   =  vl - ( vl.Dot( v1 ))*v1;
  TVector3 vD0    =  vls0.Mag()*(vls0.Dot(v1))*v0;
  TVector3 vD1    =  vls1.Mag()*(vls1.Dot(v0))*v1;
  
  vD0.Print();
  vD1.Print();

  // Calculate closest points on two lines
  TVector3 vC0    = line0.Offset + vD0;
  TVector3 vC1    = line1.Offset - vD1;

  /// Center Position
  vec.SetXYZ((vC0.X() + vC1.X())/2.,
	     (vC0.Y() + vC1.Y())/2.,
	     (vC0.Z() + vC1.Z())/2.);  
  /// distance between two lines
  dist = TMath::Abs( d0 - d1 );
  return true;  
}
