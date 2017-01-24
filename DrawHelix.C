void DrawHelix(){
  TGraph2D* gr = new TGraph2D();
  for( int i = 0; i< 1000; i++){
    Double_t x = TMath::Cos(i*TMath::DegToRad());
    Double_t y = TMath::Sin(i*TMath::DegToRad());
    Double_t z = i/20.;
    gr->SetPoint( i,x,y,z );
  }
  gr->SetMarkerColor(2);
  gr->Draw("APL");
}
