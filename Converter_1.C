void Converter_1(){
  gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  gSystem->Load("./lib/libanalysis_method");

  std::string ifilename = "HDTest_2lambda.root";
  std::string ofilename = "TestHit.root";
  TFile* itf = new TFile(ifilename.c_str());
  TTree* itr = (TTree*)itf->Get("eventTree00");
  GsimGenParticleData*   particle = new GsimGenParticleData();
  GsimDetectorEventData* tpcData = new GsimDetectorEventData();

  itr->SetBranchAddress("GenParticle.",&particle);
  itr->SetBranchAddress("TPC.",&tpcData);

  TFile* otf = new TFile(ofilename.c_str(),"recreate");
  TTree* otr = new TTree("padHit","TPCPadData");
  Int_t EventID;
  TClonesArray* HitArr = new TClonesArray("TPCPadHit");
  TClonesArray* clArr  = new TClonesArray("TPCPadHitCluster");

  std::cout<<"Branch" << std::endl;

  otr->Branch("EventID",&EventID,"EventID/I");
  otr->Branch("GenParticle.",&particle,16000,99);
  otr->Branch("HitArr",&HitArr,256000,99);
  otr->Branch("clArr",&clArr,2560000,99);

  std::cout<< "Event Loop" << std::endl;
  for( int ievt = 0; ievt < 100; ievt++){
    itr->GetEntry(ievt);
    std::cout<< ievt << std::endl;
    HitArr->Clear();
    EventID = ievt;
    TClonesArray* trackArr = particle->briefTracks;
    TClonesArray* tpcHitArr= tpcData->hits;
    Int_t nTrack = trackArr->GetEntries();
    Int_t nHit   = tpcHitArr->GetEntries();
    std::vector<TPCPadHit>        padHitArr = ConvertTPC(tpcData, particle);
    std::vector<TPCPadHitCluster> padClArr  = HitClustering( padHitArr );

    ///Save HitArr
    for( int i = 0; i< padHitArr.size(); i++){
      new((*HitArr)[i]) TPCPadHit(padHitArr.at(i));
    }
    ///Save ClusterArr
    for( int i = 0; i< padClArr.size(); i++){
      std::cout<< padClArr.at(i).RowID() << std::endl;
      new((*clArr)[i]) TPCPadHitCluster(padClArr.at(i));
    }
    otr->Fill();
  }

  std::cout<< "Write" << std::endl;
  otr->Write();
  HitArr->Clear();
  std::cout<< "Close" << std::endl;
  otf->Close();
}
