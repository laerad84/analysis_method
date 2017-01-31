void DrawMomentumResolution(){
  gStyle->SetOptStat(0);

  TH2D* hisPRes[2];
  char* name[2] ={"Proton","Pion"};
  char* filename[2] ={"Proton10k","Pion"};
  char* parameter[2]={"pMom","piMom"};
  for( int i = 0; i< 2; i++){
    hisPRes[i]= new TH2D(Form("hisPRes%d",i),Form("Momentum Resolution of %s;Momentum [MeV];P_{rec}/P_{MC}",name[i]),20,0,800,200,0,2);
  }
  for( int i = 0; i< 2; i++){
    TChain* chProton = new TChain("Output");
    chProton->Add(Form("%s_*_trk_ana.root",filename[i]));
    chProton->Project(hisPRes[i]->GetName(),Form("MergedSprialArr.R/3.36/%s.Mag():%s.Mag()",parameter[i],parameter[i]),"nTrackAfter==1");
    chProton->Draw(Form("MergedSprialArr.R/3.36/%s.Mag():%s.Mag()",parameter[i],parameter[i]),"nTrackAfter==1");
  }
  TF1* func = new TF1("func","gaus",0.9,1.1);
  TH1D* hisRes[2];
  hisPRes[0]->FitSlicesY(func);
  hisRes[0] = (TH1D*)gROOT->FindObject(Form("%s_2",hisPRes[0]->GetName()));
  hisRes[0]->SetTitle("Momentum Resolution;Momentum [MeV];#Delta p/p");
  //hisRes[0]->Draw();
  hisPRes[1]->FitSlicesY(func);
  hisRes[1] = (TH1D*)gROOT->FindObject(Form("%s_2",hisPRes[1]->GetName()));
  hisRes[1]->SetTitle("Momentum Resolution;Momentum [MeV];#Delta p/p");
  //hisRes[1]->Draw();

  hisRes[0]->SetMarkerStyle(20);
  hisRes[1]->SetMarkerStyle(21);
  hisRes[0]->SetLineColor(1);
  hisRes[0]->SetMarkerColor(1);
  hisRes[1]->SetLineColor(2);
  hisRes[1]->SetMarkerColor(2);

  TCanvas* can = new TCanvas("can","can",2100,700);
  can->Divide(3,1);
  can->cd(1);
  hisPRes[0]->Draw("colz");
  can->cd(2);
  hisPRes[1]->Draw("colz");
  can->cd(3);
  hisRes[0]->Draw();
  hisRes[1]->Draw("same");
  TLegend* leg = new TLegend(0.1,0.6,0.4,0.9);
  leg->AddEntry(hisRes[0],"Proton");
  leg->AddEntry(hisRes[1],"Pion");
  leg->Draw();
}
