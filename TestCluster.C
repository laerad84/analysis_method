void TestCluster(){
  gSystem->Load("lib/libanalysis_method.dylib");
  TFile* tf = new TFile("TestCluster.root","recreate");
  TTree* tr = new TTree("test","test");
  Int_t nArr = 0;
  TClonesArray* arr = new TClonesArray("TempCluster");
  tr->Branch("nArr",&nArr,"nArr/I");
  tr->Branch("arr",&arr,256000,99);

  for( int i = 0; i< 20; i++){
    nArr = 0;
    for( int j = 0; j< i+1; j++){
      TempCluster cl;
      cl.Init();
      cl.ID = j;
      for( int im = 0; im < i; im++){
	if( im >= 15 ){ continue; }
	cl.Mother[im] = im;
	cl.nMother = im+1;
	cl.idList.push_back(im);
      }

      new((*arr)[j]) TempCluster(cl);
      nArr++;
    }
    tr->Fill();
  }
  tr->Write();
  tf->Close();
}
