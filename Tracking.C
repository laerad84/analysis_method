void Tracking(){


  gSystem->Load("~/local/hep/E42/e42/lib/so/libGsimData");
  gSystem->Load("~/work/E42_ana/TestField/lib/libTestField");
  gSystem->Load("~/work/E42_ana/TPCDataAnalyzer/lib/libTPCDataAnalyzer");
  gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method");

  TFile* irawtf = new TFile("protonTest.root");
  TTree* irawtr = (TTree*)irawtf->Get("eventTree00");
  GsimDetectorEventData* det = new GsimDetectorEventData();
  irawtr->SetBrachAddress("TPCPad.",&det);

  TFile* otf = new TFile("Track.root","recreate");
  TTree* otr = new TTree("TrackTree","TrackEvent");


  for( int ievent = 0; ievent < 1; ievent++){
    irawtr->GetEntry( ievent );
    TClonesArray* arr = det->hits;
    Int_t nHits = arr->GetEntries();

    for( int ihit = 0; ihit < nHits; ihit++){


    }




  }

}
