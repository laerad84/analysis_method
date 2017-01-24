void Analysis_HDXi(){
  gSystem->Load("lib/libanalysis_method.dylib");

  TFile* itf = new TFile("circleHDXi.root");
  TTree* itr = (TTree*)itf->Get("analysisTree");

  Int_t iEvent;
  Int_t iMCEvent;
  Int_t ABORT;
  Int_t NCl;
  Int_t NSprial;
  Int_t NSprialM;
  Int_t NSprialP;

  TClonesArray* SprialArr = new TClonesArray("HSprial");
  TClonesArray* SprialMArr = new TClonesArray("HSprial");
  TClonesArray* SprialPArr = new TClonesArray("HSprial");

  TClonesArray* ClusterArr= new TClonesArray("TPCCluster");
  TClonesArray* Particle = new TClonesArray("TLorentzVector");
  TClonesArray* Vertex   = new TClonesArray("TVector3");
  Int_t NParticle=0;
  Int_t Index[30];
  itr->SetBranchAddress("iEvent",&iEvent);
  itr->SetBranchAddress("iMCEvent",&iMCEvent);
  itr->SetBranchAddress("ABORT",&ABORT);

  itr->SetBranchAddress("SprialArr",&SprialArr);
  itr->SetBranchAddress("SprialMArr",&SprialMArr);
  itr->SetBranchAddress("SprialPArr",&SprialPArr);
  itr->SetBranchAddress("ClusterArr",&ClusterArr);
  itr->SetBranchAddress("NCl",&NCl);
  itr->SetBranchAddress("NSprial",&NSprial);
  itr->SetBranchAddress("NSprialM",&NSprialM);
  itr->SetBranchAddress("NSprialP",&NSprialP);
  itr->SetBranchAddress("NParticle",&NParticle);
  itr->SetBranchAddress("Index",Index);
  itr->SetBranchAddress("Particle",&Particle);
  itr->SetBranchAddress("Vertex",&Vertex);

  TH1D* hislambdaMass  = new TH1D("hislambdaMass","lambda Mass",125,1000,2500);
  TH1D* hisLambdaM[2];
  for( int i = 0; i<2; i++){
    hisLambdaM[i] = new TH1D(Form("hisLambdaM_%d",i),Form("hisLambda_%d",i),125,1000,2500);
  }
  for( int ievt =0; ievt < itr->GetEntries(); ievt++){
    itr->GetEntry(ievt);
    if(!(NSprialM == 2 && NSprialP == 3)){
      continue;
    }
    std::cout<< ievt << std::endl;
    std::cout<< Particle->GetEntries() << std::endl;
    if( Particle->GetEntries() != 5){ continue;}
    TLorentzVector* lv[5];
    for( int i = 0; i< 5; i++){
      lv[i] = new TLorentzVector();
      lv[i] = (TLorentzVector*)Particle->At(i);
    }
    std::cout<< "lambda" << std::endl;
    std::cout<< Vertex->GetEntries() << std::endl;
    TVector3* vec[2];
    TVector3* position[2];
    Double_t dist[2];
    dist[0] = 0;
    dist[1] = 1;

    for( int i = 0; i< 2; i++){
      vec[i] = new TVector3();
      vec[i] = (TVector3*)Vertex->At(i);
      position[i] = new TVector3(vec[i].X(),vec[i].Y(),vec[i].Z()+143);
      dist[i] = TMath::Sqrt( vec[i]->X()*vec[i]->X() + (vec[i].Z()+143)*(vec[i].Z()+143));
    }


    TLorentzVector* xiLv = new TLorentzVector();
    TLorentzVector* lambda[2];
    for( int i = 0; i< 2; i++){
      std::cout<< i << std::endl;
      lambda[i] = new TLorentzVector();
      lambda[i]->SetPxPyPzE( lv[i*2]->Px() + lv[i*2+1]->Px(),
			     lv[i*2]->Py() + lv[i*2+1]->Py(),
			     lv[i*2]->Pz() + lv[i*2+1]->Pz(),
			     lv[i*2]->E() + lv[i*2+1]->E());


      lambda[i]->Print();
    }



    std::cout<< "Fill" << std::endl;
    hislambdaMass->Fill(lambda[0]->M());
    hislambdaMass->Fill(lambda[1]->M());
    if( dist[0] < dist[1]){
      hisLambdaM[0]->Fill(lambda[0]->M());
      hisLambdaM[1]->Fill(lambda[1]->M());
    }else{
      hisLambdaM[0]->Fill(lambda[1]->M());
      hisLambdaM[1]->Fill(lambda[0]->M());
    }
  }
  hislambdaMass->Draw();
  hisLambdaM[0]->SetLineColor(2);
  hisLambdaM[1]->SetLineColor(3);
  hisLambdaM[0]->Draw("same");
  hisLambdaM[1]->Draw("same");
}
