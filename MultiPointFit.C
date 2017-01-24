
TGraphErrors* gre;
void chisqCal( Int_t &, Double_t * Double_t &f, Double_t *par, Int_t ){
  Int_t np = gre->GetN();
  f = 0;
  Double_t *x = gre->GetX();
  Double_t *y = gre->GetY();
  Double_t *ex = gre->GetEX();
  Double_t *ey = gre->GetEY();
  for( int ip  =0; ip < np; ip++){
    Double_t u = (x[ip] - par[0])/ex[ip];
    Double_t v = (y[ip] - par[1])/ey[ip];
    Double_t dx = (x[ip] - par[0]);
    Double_t dy = (y[ip] - par[1]);
    Double_t dr = par[2] - TMath::Sqrt(dx*dx + dy*dy);
    f+= dr*dr;
  }
}

void MultiPointFit(){
  gre = new TGraphErrors();
  TCanvas* can = new TCanvas("can","can",800,800);
  Int_t nPoint = 5;
  can->DrawFrame(-5,-5,5,5);

  TRandom3 r;
  Double_t x,y;
  for( int i = 0; i< nPoint; i++){
    r.Circle(x,y,r.Gaus(4,0.3));
    gre->SetPoint( i, x, y );
  }

  gre->Draw("P");

  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter * fitter = TVirtualFitter::Fitter(0,3);
  fitter->SetFCN(chisqCal);
  fitter->SetParameter( 0, "x0", 0,0.1,0,0);
  fitter->SetParameter( 1, "y0", 0,0.1,0,0);
  fitter->SetParameter( 2, "R", 1, 0.1, 0,0);

  Double_t arglist[1] = {0};
  fitter->ExecuteCommand("MIGRAD",arglist, 0);

  TArc* arc = new TArc(fitter->GetParameter(0),
		       fitter->GetParameter(1),
		       fitter->GetParameter(2));
  arc->SetLineColor(kRed);
  arc->SetLineWidth(4);
  arc->SetFillStyle(0);
  arc->Draw();

}
