TSpline3* spl;

double func( double* x, double* par ){
  double value = par[0]+ par[1]*( spl->Eval( x[0] - par[2]));
  if( TMath::Abs(x[0] - par[2] ) > 10 ){ return par[0];}
  return value;
}

void SPL(){

  TH1D* his = new TH1D("his","his",100,-10,10);
  for( int i = 0; i< 1000; i++){
    his->Fill(gRandom->Gaus(0,3));
  }
  spl = new TSpline3(his);
  TF1* f = new TF1("f",func,-10,10,3);
  f->SetParameters( 10, 2, 5);

  TCanvas* can = new TCanvas("can","can",1600,800);
  can->Divide(2,1);
  can->cd(1);
  spl->Draw();

  can->cd(2);
  f->Draw();
}
