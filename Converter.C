void Converter(){
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
  otr->Branch("EventID",&EventID,"EventID/I");
  otr->Branch("GenParticle.",&particle,16000,99);
  DataHandler* datahandler = new DataHandler(otr);
  datahandler->BranchOfConverted();


  //for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
  for( int ievt = 0; ievt < 20; ievt++){
    itr->GetEntry(ievt);
    EventID = ievt;
    TClonesArray* trackArr = particle->briefTracks;
    TClonesArray* tpcHitArr= tpcData->hits;

    Int_t nTrack = trackArr->GetEntries();
    Int_t nHit   = tpcHitArr->GetEntries();

    //std::cout<< nTrack << "\t" << nHit << std::endl;

    //padHitArr.clear();
    std::vector<TPCPadHit> padHitArr;
    std::cout<< padHitArr.size() << std::endl;
    padHitArr = ConvertTPC(tpcData, particle);
    std::cout<< padHitArr.size() << std::endl;
    datahandler->SetData(padHitArr);

    /*
    std::cout<< padHitArr.size() << std::endl;
    for( int i = 0; i< padHitArr.size(); i++){

      std::cout<< i                           << "\t"
	       << padHitArr.at(i).Col()       << "\t"
	       << padHitArr.at(i).Row()       << "\t"
	       << padHitArr.at(i).Energy()    << "\t"
	       << padHitArr.at(i).MCTrackID() << "\t"
	       << padHitArr.at(i).TrackID()   << "\t"
	       << std::endl;

    }
    */
    otr->Fill();
  }
  otr->Write();
  otf->Close();
}
