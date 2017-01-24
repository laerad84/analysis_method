void CalculatePoint(){


}

  void CalCross(){
  gROOT->ProcessLine("#include <vector>");
  gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  gSystem->Load("./lib/libanalysis_method");

  Double_t r0 = 20;
  Double_t x0 = 5;
  Double_t y0 = 5;
  Double_t r1 = 15;
  Double_t x1 = 15;
  Double_t y1 = 15;

  TVector2 v0(5,5);
  TVector2 v1(15,15);
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
  TMarker* markCent = new TMarker(vecCent.X(), vecCent.Y(), 21);
  TMarker* mark[2];
  for( int i = 0; i< 2; i++){
    mark[i] = new TMarker(vec[i].X(),vec[i].Y(),20);
  }

  std::cout<< dist << std::endl;
  unit.Print();

  TCanvas* can = new TCanvas("can","can",800,800);
  gPad->DrawFrame(-100,-100,100,100);
  TArc* arc0 = new TArc(x0,y0,r0);
  TArc* arc1 = new TArc(x1,y1,r1);

  arc0->SetFillColor(0);
  arc1->SetFillColor(0);
  arc0->SetFillStyle(0);
  arc1->SetFillStyle(0);


  arc0->Draw();
  arc1->Draw();
  mark[0]->Draw();
  mark[1]->Draw();
  markCent->Draw();
  func->Draw("same");
  }
