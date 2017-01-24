void DataClusteringTest(){
  gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  gSystem->Load("./lib/libanalysis_method");

  TFile* itf  = new TFile("HDTest_2lambda.root");
  TTree* itr  = (TTree*)itf->Get("eventTree00");
  TPCDataConverter* converter = new TPCDataConverter(itr);

  TCanvas* can = new TCanvas("can","can",1600,800);
  can->Divide(2,1);
  can->Draw();
  TH3D* his    = new TH3D("his","his;X;Y;Z",150,-300,300,150,-300,300,150,-300,300);
  TH2D* hisHit = new TH2D("hisHit","his;X;Z",150,-300,300,150,-300,300);
  TPCPoly* polyHit = new TPCPoly("polyHit","test");
  TPCPoly* polyCluster = new TPCPoly("polyCluster","test");

  //his->SetMarkerStyle(20);
  //hisHit->SetMarkerStyle(21);
  hisHit->SetMarkerColor(2);

  for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
    //for( int ievt = 0; ievt < 1; ievt++){

    //std::vector<TPCPadHit> HitArr = converter->Convert(ievt);
    std::cout<< "Convert" << std::endl;
    converter->Convert(ievt);
    std::cout<< "Clustering " << std::endl;
    converter->Clustering();
    std::cout<< "End clustering " << std::endl;
    std::vector<TPCPadHitCluster> clusterArr = converter->GetClusterArr();
    std::cout<< "GetCluster" << std::endl;
    std::vector<TPCPadHitCluster> cl[4];
    TGraph* grX[4];
    TGraph* grY[4];
    int nCluster = 0;
    for( int i = 0; i< 4; i++){
      cl[i] = TPCFindBlock( clusterArr );
      if( cl[i].size() < 3){
	cl[i] = TPCFindBlock( clusterArr );
      }
      if( cl[i].size() != 0 ){
	nCluster++;
	grX[i] = DrawClusterArrXZ(cl[i]);
	grY[i] = DrawClusterArrYZ(cl[i]);
	grX[i]->SetMarkerColor(i+3);
	grY[i]->SetMarkerColor(i+3);
	grX[i]->SetMarkerStyle(20);
	grY[i]->SetMarkerStyle(20);
      }
    }

    if( nCluster < 4 ){ continue; }

    his->Reset();
    hisHit->Reset();
    polyHit->Reset();
    polyCluster->Reset();

    /*
    TGraph* gr = new TGraph();
    TGraph* grCluster = new TGraph();

    gr->SetMarkerStyle(21);
    grCluster->SetMarkerStyle(21);
    grCluster->SetMarkerColor(2);
    for( int i = 0; i< HitArr.size(); i++){
      gr->SetPoint( i, HitArr.at(i).Position().Z(), HitArr.at(i).Position().X());
      //std::cout<< HitArr.at(i).TrackID() << std::endl;
      his->Fill( HitArr.at(i).Position().X(),
		 HitArr.at(i).Position().Y(),
		 HitArr.at(i).Position().Z()
		 );
      polyHit->Fill(HitArr.at(i).Position().X(), HitArr.at(i).Position().Z());
    }


    std::cout<< "Clustering" << std::endl;

    std::vector<TPCPadHitCluster> clusterArr;

    for( int irow = 0; irow < 32; irow++){
      Bool_t  bRun = kTRUE;

      while( HitArr.size() != 0 ){
	Int_t ilength  = HitArr.size();
	TPCPadHitCluster cluster = TPCClusterer( irow, HitArr );
	Int_t flength  = HitArr.size();
	Int_t clength  = cluster.HitArr().size();
	std::cout<< "Row: " << irow << "\t" << ilength << "\t" << flength << "\t" << clength << std::endl;
	if( cluster.HitArr().size() != 0 ){
	  clusterArr.push_back(cluster);
	  //std::cout<< cluster.RowID() << "\t" << cluster.ColID() << "\t" << cluster.GetClusterPos().X() <<"\t" << cluster.GetClusterPos().Z() << std::endl;
	  polyCluster->Fill(cluster.GetClusterPos().X(),cluster.GetClusterPos().Z());
	  grCluster->SetPoint( grCluster->GetN(), cluster.GetClusterPos().Z(),cluster.GetClusterPos().X() );
	}else{
	  break;
	}
      }
    }
    */

    converter->grPadHitXZ->SetMarkerStyle(21);
    converter->grClusterXZ->SetMarkerStyle(22);
    converter->grClusterXZ->SetMarkerColor(2);
    converter->grPadHitYZ->SetMarkerStyle(21);
    converter->grClusterYZ->SetMarkerStyle(22);
    converter->grClusterYZ->SetMarkerColor(2);

    can->cd(1);
    gPad->DrawFrame(-300,-300,300,300);
    converter->grPadHitXZ->Draw("P");
    converter->grClusterXZ->Draw("P");
    for( int i = 0; i< 4; i++){
      grX[i]->Draw("P");
    }
    can->cd(2);
    gPad->DrawFrame(-300,-300,300,300);
    converter->grPadHitYZ->Draw("P");
    converter->grClusterYZ->Draw("P");
    for( int i = 0; i< 4; i++){
      grY[i]->Draw("P");
    }
    can->Update();
    can->Modified();
    gSystem->ProcessEvents();
    getchar();

  }
}
