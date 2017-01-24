void DataConvTester(){
  gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  gSystem->Load("./lib/libanalysis_method");

  TFile* itf  = new TFile("HDTest_2lambda.root");
  TTree* itr  = (TTree*)itf->Get("eventTree00");
  TPCDataConverter* converter = new TPCDataConverter(itr);

  TCanvas* can = new TCanvas("can","can",1600,800);
  can->Divide(2,1);
  TH3D* his    = new TH3D("his","his;X;Y;Z",150,-300,300,150,-300,300,150,-300,300);
  TH2D* hisHit = new TH2D("hisHit","his;X;Z",150,-300,300,150,-300,300);
  //his->SetMarkerStyle(20);
  //hisHit->SetMarkerStyle(21);
  hisHit->SetMarkerColor(2);

  for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
    //for( int ievt = 0; ievt < 1; ievt++){
    std::vector<TPCPadHit> HitArr = converter->Convert(ievt);
    std::cout<< HitArr.size() << std::endl;
    his->Reset();
    hisHit->Reset();
    TGraph* gr = new TGraph();
    gr->SetMarkerStyle(21);
    for( int i = 0; i< HitArr.size(); i++){
      gr->SetPoint( i, HitArr.at(i).Row(), HitArr.at(i).Col());
      //std::cout<< HitArr.at(i).TrackID() << std::endl;
      his->Fill( HitArr.at(i).Position().X(),
		 HitArr.at(i).Position().Y(),
		 HitArr.at(i).Position().Z()
		 );
    }

    /*
    GsimDetectorEventData* det    = converter->GetDetectorEventData();
    TClonesArray* arr = det->hits;
    for( int i = 0; i< arr->GetEntries(); i++){
      GsimDetectorHitData* hit = (GsimDetectorHitData*)arr->At(i);
      hisHit->Fill( hit->r.Z(),hit->r.X());
    }
    */
    for( int i = 0; i< 3000; i++){
      Int_t row = -1, col=-1;
      TVector2 vec = gTPCIDHandler->GetXZ( i );
      bool f = gTPCIDHandler->GetRowCol( vec.X(), vec.Y(), row, col );
      if( !f ){continue;}
      TVector2 vec1 = gTPCIDHandler->GetXZ( row, col );
      hisHit->Fill(vec1.Y(), vec1.X());
    }

    //gr->Draw("AP");
    can->cd(1);
    his->Draw();
    can->cd(2);
    hisHit->Draw();
    can->Update();
    can->Modified();
    getchar();
    gSystem->ProcessEvents();
  }
}
