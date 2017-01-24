void DrawHD(){
  gSystem->Load("lib/libanalysis_method.dylib");

  TFile* itf = new TFile("circleHD.root");
  TTree* itr = (TTree*)itf->Get("analysisTree");

  TCanvas* can = new TCanvas("can","can",800,800);
  can->DrawFrame(-300,-300,300,300);
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


  TH1D* hislambdaMass  = new TH1D("hislambdaMass","lambda Mass",1000,1000,2000);
  TH1D* hisHDMass = new TH1D("hisHDMass","hisHDMass",1000,1000,4000);
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
    TLorentzVector* lambda[2];
    TLorentzVector* HDlv = new TLorentzVector();
    for( int i = 0; i< 2; i++){
      std::cout<< i << std::endl;
      lambda[i] = new TLorentzVector();
      lambda[i]->SetPxPyPzE( lv[i*2]->Px() + lv[i*2+1]->Px(),
			     lv[i*2]->Py() + lv[i*2+1]->Py(),
			     lv[i*2]->Pz() + lv[i*2+1]->Pz(),
			     lv[i*2]->E() + lv[i*2+1]->E());


      lambda[i]->Print();
    }
    HDlv->SetPxPyPzE(lambda[0]->Px()+lambda[1]->Px(),
		     lambda[0]->Py()+lambda[1]->Py(),
		     lambda[0]->Pz()+lambda[1]->Pz(),
		     lambda[0]->E()+lambda[1]->E());
    hisHDMass->Fill(HDlv->M());
    std::cout<< "Fill" << std::endl;
    hislambdaMass->Fill(lambda[0]->M());
    hislambdaMass->Fill(lambda[1]->M());

    TArc* arc[30];
    for(int i = 0; i< SprialArr->GetEntries(); i++){
      HSprial* sprial = (HSprial*)SprialArr->At(i);
      arc[i] = new TArc(sprial.Z,sprial.X,sprial.R);
      arc[i]->SetFillStyle(0);
      arc[i]->Draw();
      if(sprial.RL > 0){
	arc[i]->SetLineColor(2);
      }else{
	arc[i]->SetLineColor(3);
      }
    }
    can->Update();
    can->Modified();
    gSystem->ProcessEvents();
    getchar();
    for( int i = 0; i< SprialArr->GetEntries(); i++){
      delete arc[i];
    }
  }
  hislambdaMass->Draw();
  //hisHDMass->Draw();
}
