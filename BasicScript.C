void BasicScript(){
  gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  gSystem->Load("./lib/libanalysis_method");

  std::string ifilename = "HDTest_2lambda.root";
  TFile* itf = new TFile(ifilename.c_str());
  TTree* itr = (TTree*)itf->Get("eventTree00");
  GsimGenParticleData*   particle = new GsimGenParticleData();
  GsimDetectorEventData* tpcData = new GsimDetectorEventData();

  itr->SetBranchAddress("GenParticle.",&particle);
  itr->SetBranchAddress("TPC.",&tpcData);
  std::vector<TPCPadHit> padHitArr;

  for( int ievt = 0; ievt < 10; ievt++){
    itr->GetEntry(ievt);

    TClonesArray* trackArr = particle->briefTracks;
    TClonesArray* tpcHitArr= tpcData->hits;

    Int_t nTrack = trackArr->GetEntries();
    Int_t nHit   = tpcHitArr->GetEntries();

    std::cout<< nTrack << "\t" << nHit << std::endl;
    padHitArr.clear();
    padHitArr = ConvertTPC(tpcData, particle);
    std::cout<< padHitArr.size() << std::endl;
    for( int i = 0; i< padHitArr.size(); i++){

      std::cout<< i << "\t" << padHitArr.at(i).Col() << "\t" << padHitArr.at(i).Row() << std::endl;
    }
  }



}
