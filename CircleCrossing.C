TVector3 CircleCrossing(TVector3 rxy0, TVector3 rxy1){
  gROOT->ProcessLine("#include <vector>");
  gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  gSystem->Load("./lib/libanalysis_method");


  Double_t r0 = rxy0.X();
  Double-t r1 = rxy1.X();

  TVector2 v0( rxy0.Y(), rxy0.Z());
  TVector2 v1( rxy1.Y(), rxy1.Z());

  Double_t dist = (v1 - v0).Mod();
  Double_t phi  = (v1 - v0).Phi();
  TVector2 unit = (v1 - v0).Unit();

  double a = (v1.X() - v0.X() );
  double b = (v1.Y() - v0.Y() );
  double c = (x1*x1 - x0*x0 + y1*y1 - y0*y0 + r0*r0 - r1*r1)/2.;

  double distPoint = TMath::Abs( a*v0.X() + b*v0.Y() - c)/TMath::Sqrt(a*a + b*b);
  double dTheta    = TMath::ACos( distPoint/r0);
  double Phi       = unit.Phi();

  double slide = a/(-b);
  double offset = c/(b);


  TF1* func = new TF1("func","[0]+[1]*x",-100,100);
  func->SetParameter(0,offset);
  func->SetParameter(1,slide);
  std::cout<< distPoint << "\t" << r0 << "\t" << dTheta << "\t" << Phi << std::endl;
  TVector2 vecCent;
  vecCent.SetMagPhi( distPoint, Phi);
  TVector2 vec[2];
  vec[0].SetMagPhi(r0,Phi + dTheta );
  vec[1].SetMagPhi(r0,Phi - dTheta );
  vec[0] = vec[0] + v0;
  vec[1] = vec[1] + v0;


  if( vec[0].Mag() < vec[1].Mag() ){
    return vec[0];
  }else{
    return vec[1];6

  }

}
