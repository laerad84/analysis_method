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
HSprial GenerateSprial( TVector3 init, TVector3 middle, TVector3 final){

  HCircle circle = ThreePointCircle( init, middle, final);
  TVector3 parameter = DyDtheta( init, middle, final );
  HSprial sprial_result( circle.ID, circle.X, (init.Y()+final.Y())/2., circle.Z, circle.R, parameter.X(),parameter.Y(),(int)(parameter.Z()));
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
    marker->SetPoint( iarr, hit->Position.X(), hit->Position.Y(), hit->Position.Z());
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
Double_t CalculateDeltaR( HSprial sprial, TVector3 vec){
  Double_t dx = vec.X() - sprial.X;
  Double_t dz = vec.Z() - sprial.Z;
  Double_t dr = TMath::Sqrt( dx*dx + dz*dz ) - sprial.R;
  return dr;
}
Double_t CalculateDeltaY( HSprial sprial, TVector3 vec ){
  Double_t initTheta = TMath::ATan2(sprial.InitPos.X()-sprial.X,sprial.InitPos.Z()-sprial.Z);
  Double_t theta     = TMath::ATan2( vec.X() - sprial.X, vec.Z() - sprial.Z);
  Double_t dTheta    = 0;
  dTheta  = theta - initTheta;
  if( dTheta >= 2*TMath::Pi() ){ dTheta = dTheta - 2*TMath::Pi();}
  if( dTheta <= -2*TMath::Pi()){ dTheta = dTheta - 2*TMath::Pi();}

  Double_t dy    = vec.Y() - dTheta*sprial.DY/sprial.DTheta - sprial.Y;
  return dy;
}
Double_t CalculateY( HSprial sprial, TVector3 vec ){
  Double_t initTheta = TMath::ATan2(sprial.InitPos.X()-sprial.X,sprial.InitPos.Z()-sprial.Z);
  Double_t theta     = TMath::ATan2( vec.X() - sprial.X, vec.Z() - sprial.Z);
  Double_t dTheta    = theta - initTheta;
  //if( dTheta >= 2*TMath::Pi() ){ dTheta = dTheta - 2*TMath::Pi();}
  //if( dTheta <= -2*TMath::Pi()){ dTheta = dTheta - 2*TMath::Pi();}
  Double_t y    = dTheta*sprial.DY/sprial.DTheta + sprial.InitPos.Y();
  return y;
}
Bool_t CalculateCircleCrossing( HSprial sprial0, HSprial sprial1, TVector3& vec0, TVector3& vec1 ){
  Double_t r0 = sprial0.R;
  Double_t r1 = sprial1.R;
  TVector2 v0(sprial0.X,sprial1.Z);
  TVector2 v1(sprial1.X,sprial1.Z);
  Double_t dist = (v1 - v0).Mod();
  if( dist > r0 +r1 ||
      dist < TMath::Abs(r0 -r1 )){
    vec0=TVector3(0,0,0);
    vec1=TVector3(0,0,0);
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
Bool_t CalculateCircleTangent( HSprial sprial0, TVector3 pos, TVector3& vec ){
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

Bool_t CalculateCrossing( HSprial sprial0, HSprial sprial1, Double_t& dist, TVector3& pos ){
  TVector3 v0;
  TVector3 v1;
  if( !CalculateCircleCrossing( sprial0, sprial1, v0,v1) ){
    dist = 1e8;
    pos= TVector3(0,0,0);
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
double CalculateMomentumR( HSprial sprial, double dtesla){
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
  return sprial.R * dtesla / 3.36;
}
double CalculateMomentumY( HSprial sprial, double dtesla){
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
  double momentumR = CalculateMomentumR( sprial, dtesla );
  return momentumR * sprial.DY / TMath::Abs( sprial.R * sprial.DTheta );
}

double CalculateMomentum( HSprial sprial, double dtesla ){
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
  double pR = CalculateMomentumR( sprial, dtesla );
  double pY = CalculateMomentumY( sprial, dtesla );
  double totalMomentum = TMath::Sqrt( pR*pR + pY*pY );
  return totalMomentum;
}

TPolyMarker3D* GenerateCrossingPoints( std::vector<HSprial> PSprial, std::vector<HSprial> MSprial ){
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
