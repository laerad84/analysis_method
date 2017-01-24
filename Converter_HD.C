void Converter_HD(){
  gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  gSystem->Load("./lib/libanalysis_method");


  std::string ifilename = "lambda.root";
  //std::string ifilename = "Proton.root";
  //std::string ifilename = "HEXiE42.root";
  std::string ofilename = ifilename.substr(0,ifilename.rfind(".")) + "_Aconv.root";
  std::cout<< ifilename << std::endl;
  std::cout<< ofilename << std::endl;

  TFile* itf = new TFile(ifilename.c_str());
  TTree* itr = (TTree*)itf->Get("eventTree00");
  GsimGenParticleData*   particle = new GsimGenParticleData();
  GsimDetectorEventData* tpcData = new GsimDetectorEventData();
  GsimDetectorEventData* ftofData = new GsimDetectorEventData();
  GsimDetectorEventData* DC1Data  = new GsimDetectorEventData();
  GsimDetectorEventData* DC2Data  = new GsimDetectorEventData();
  GsimDetectorEventData* DC3Data  = new GsimDetectorEventData();
  GsimDetectorEventData* NBARData = new GsimDetectorEventData();
  itr->SetBranchAddress("GenParticle.",&particle);
  itr->SetBranchAddress("TPC.",&tpcData);
  itr->SetBranchAddress("FTOF.",&ftofData);
  itr->SetBranchAddress("DC1.",&DC1Data);
  itr->SetBranchAddress("DC2.",&DC2Data);
  itr->SetBranchAddress("DC3.",&DC3Data);
  itr->SetBranchAddress("NBAR.",&NBARData);

  TFile* otf = new TFile(ofilename.c_str(),"recreate");
  TTree* otr = new TTree("convTree","TPCPadData");
  Int_t EventID;
  TClonesArray* HitArr = new TClonesArray("TPCPadHit");
  TClonesArray* clArr  = new TClonesArray("TPCPadHitCluster");

  otr->Branch("EventID",&EventID,"EventID/I");
  otr->Branch("GenParticle.",&particle,16000,99);
  otr->Branch("HitArr",&HitArr,256000,99);
  otr->Branch("DC1.",&DC1Data,256000,99);
  otr->Branch("DC2.",&DC2Data,256000,99);
  otr->Branch("DC3.",&DC3Data,256000,99);
  otr->Branch("NBAR.",&NBARData,256000,99);
  otr->Branch("FTOF.",&ftofData,256000,99);

  //otr->Branch("clArr",&clArr,2560000,99);

  std::cout<< "Event Loop" << std::endl;
  for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
  //for( int ievt = 0; ievt < 1000; ievt++){
    itr->GetEntry(ievt);
    std::cout<< ievt << std::endl;
    HitArr->Clear();
    EventID = ievt;
    TClonesArray* trackArr = particle->briefTracks;
    TClonesArray* tpcHitArr= tpcData->hits;
    Int_t nTrack = trackArr->GetEntries();
    Int_t nHit   = tpcHitArr->GetEntries();
    std::vector<TPCPadHit>        padHitArr = ConvertTPC(tpcData, particle);

    ///Save HitArr
    for( int i = 0; i< padHitArr.size(); i++){
      new((*HitArr)[i]) TPCPadHit(padHitArr.at(i));
    }

    /*
    ///Save ClusterArr
    for( int i = 0; i< padClArr.size(); i++){
      std::cout<< padClArr.at(i).RowID() << std::endl;
      new((*clArr)[i]) TPCPadHitCluster(padClArr.at(i));
    }
    */
    otr->Fill();
  }

  std::cout<< "Write" << std::endl;
  otr->Write();
  HitArr->Clear();
  std::cout<< "Close" << std::endl;
  otf->Close();
}
