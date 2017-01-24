void DataClusteringTest_Circle(){
  gROOT->ProcessLine("#include <vector>");
  gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  gSystem->Load("./lib/libanalysis_method");

  TFile* itf  = new TFile("HDTest_2lambda.root");
  TTree* itr  = (TTree*)itf->Get("eventTree00");
  TPCDataConverter* converter = new TPCDataConverter(itr);

  TCanvas* can = new TCanvas("can","can",2400,1600);
  can->Divide(3,2);
  can->Draw();
  TH3D* his    = new TH3D("his","his;X;Y;Z",150,-300,300,150,-300,300,150,-300,300);
  TH2D* hisHit = new TH2D("hisHit","his;X;Z",150,-300,300,150,-300,300);
  TPCPoly* polyHit = new TPCPoly("polyHit","test");
  TPCPoly* polyCluster = new TPCPoly("polyCluster","test");

  TH1D* hisChi2Temp = new TH1D("hisChi2Temp","hisChi2Temp",1000,0,100);
  TH1D* hisDistTemp = new TH1D("hisDistTemp","hisDistTemp",1000,0,1000);
  TH1D* hisDistCircle = new TH1D("hisDistCircle","hisDistCircle", 1000,0,100);
  TH1D* hisDyCircle   = new TH1D("hisDyCircle","hisDyCircle",1000,-50,50);
  TH2D* hisDDyCircle  = new TH2D("hisDDyCircle","hisDDyCircle",100,0,100,100,-50,50);
  hisChi2Temp->SetLineColor(1);
  hisDistTemp->SetLineColor(2);
  //his->SetMarkerStyle(20);
  //hisHit->SetMarkerStyle(21);
  hisHit->SetMarkerColor(2);
  TGraph* grTemp = new TGraph();
  grTemp->SetMarkerStyle(20);
  grTemp->SetMarkerColor(3);
  TArc* arc = new TArc(0,0,1);
  arc->SetLineColor(2);
  arc->SetFillStyle(0);
  arc->SetLineWidth(3);
  TArc* arcThree = new TArc(0,0,1);
  arcThree->SetLineColor(4);
  arcThree->SetFillStyle(0);
  arcThree->SetLineWidth(3);
  TPCCircleFitter* fitter  = new TPCCircleFitter();
  for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
    //for( int ievt = 0; ievt < 1; ievt++){
    //hisChi2Temp->Reset();
    grTemp->Set(0);
    //std::vector<TPCPadHit> HitArr = converter->Convert(ievt);
    std::cout<< "Convert" << std::endl;
    converter->Convert(ievt);
    std::cout<< "Clustering " << std::endl;
    converter->Clustering();
    std::cout<< "End clustering " << std::endl;
    std::vector<TPCPadHitCluster> clusterArr = converter->GetClusterArr();
    if ( clusterArr.size() < 3 ){ continue; }
    std::cout<< "GetCluster" << std::endl;
    std::vector<TPCPadHitCluster> cl[4];
    TGraph* grX[4];
    TGraph* grY[4];
    int nCluster = 0;

    //std::vector<TPCPadHitCluster> tempCluster = TPCFindEdgeBlock( clusterArr );
    std::vector<TPCPadHitCluster> tempCluster = TPCFindBlocKCircle( clusterArr );
    std::cout<< tempCluster.size() << std::endl;
    if( tempCluster.size() < 3 ) {continue; }

    /// Circle fitter
    //cl[0] = TPCFindEdgeBlock( clusterArr );
    cl[0] = tempCluster;
    std::cout<< "EDGE: " << cl[0].size() << std::endl;
    grX[0] = DrawClusterArrXZ(cl[0]);
    grY[0] = DrawClusterArrYZ(cl[0]);
    grX[0]->SetMarkerColor(0+3);
    grY[0]->SetMarkerColor(0+3);
    grX[0]->SetMarkerStyle(20+3);
    grY[0]->SetMarkerStyle(20+3);

    std::vector<double> xarr;
    std::vector<double> zarr;
    for( int icl = 0; icl < tempCluster.size(); icl++){
      xarr.push_back( tempCluster.at(icl).GetClusterPos().X());
      zarr.push_back( tempCluster.at(icl).GetClusterPos().Z());
    }

    TPCCircleFitResult result = fitter->Fit(xarr,zarr);

    result.GetFitPar().Print();
    result.GetFitParErr().Print();
    arc->SetX1( result.GetFitPar().Z());
    arc->SetY1( result.GetFitPar().Y());
    arc->SetR1( result.GetFitPar().X());
    arc->SetR2( result.GetFitPar().X());

    TPCCircleFitResult resultCircle = fitter->ThreePointCircleFit( xarr, zarr);
    arcThree->SetX1( resultCircle.GetFitPar().Z() );
    arcThree->SetY1( resultCircle.GetFitPar().Y() );
    arcThree->SetR1( resultCircle.GetFitPar().X() );
    arcThree->SetR2( resultCircle.GetFitPar().X() );

    //// Calculate dydl
    Double_t rad = resultCircle.GetFitPar().X();
    Double_t theta[3];
    Double_t dtheta[3];
    Double_t l[3];
    Double_t dl[3];
    Double_t y[3];
    Double_t dy[3];
    for( int i = 0; i< 3; i++){
      theta[i] = TMath::ATan2(xarr.at(i)- resultCircle.GetFitPar().Y(), zarr.at(i)- resultCircle.GetFitPar().Z());
      l[i]     = rad*theta[i];
      y[i]     = tempCluster.at(i).GetClusterPos().Y();
    }
    for( int i = 0; i< 3; i++){
      dtheta[i] = dtheta[(i+1)%3] - dtheta[i%3];
      dl[i]     = TMath::Abs(l[(i+1)%3] - l[i%3]);
      dy[i]     = y[(i+1)%3] - y[i];
    }

    Double_t avgDyDl = (dy[0]/dl[0] + dy[1]/dl[1])/2;

    converter->grPadHitXZ->SetMarkerStyle(21);
    converter->grClusterXZ->SetMarkerStyle(22);
    converter->grClusterXZ->SetMarkerColor(2);
    converter->grPadHitYZ->SetMarkerStyle(21);
    converter->grClusterYZ->SetMarkerStyle(22);
    converter->grClusterYZ->SetMarkerColor(2);
    his->Reset();
    hisHit->Reset();
    polyHit->Reset();
    polyCluster->Reset();



    for( int icl = 0; icl < clusterArr.size(); icl++){
      double chi2 = result.CalculateChi2( clusterArr.at(icl).GetClusterPos().X(),
					   clusterArr.at(icl).GetClusterPos().Z());
      double dist = result.CalculateDist( clusterArr.at(icl).GetClusterPos().X(),
					  clusterArr.at(icl).GetClusterPos().Z());

      double distCircle = resultCircle.CalculateDist( clusterArr.at(icl).GetClusterPos().X(),
						      clusterArr.at(icl).GetClusterPos().Z());
      double tempTheta  = TMath::ATan2( clusterArr.at(icl).GetClusterPos().X() - resultCircle.GetFitPar().Y(),
					clusterArr.at(icl).GetClusterPos().Z() - resultCircle.GetFitPar().Z());
      double tempDTheta = tempTheta - theta[2];
      double tempDl     = rad*tempDTheta;

      double tempDy     = avgDyDl * tempDl;
      double DeltaY     =  y[2] + tempDy - clusterArr.at(icl).GetClusterPos().Y();


      hisChi2Temp->Fill(chi2);
      hisDistTemp->Fill(dist);

      hisDistCircle->Fill( distCircle );
      hisDyCircle->Fill( DeltaY );
      hisDDyCircle->Fill( distCircle, DeltaY);
      //if( chi2 < 2 && dist < 10 ){
      if( TMath::Abs(DeltaY) < 15 && distCircle < 10 ){
	grTemp->SetPoint( grTemp->GetN(),
			  clusterArr.at(icl).GetClusterPos().Z(),
			  clusterArr.at(icl).GetClusterPos().X());


      }
    }

    //if( nCluster < 4 ){ continue; }


    can->cd(1);
    gPad->DrawFrame(-300,-300,300,300);
    converter->grPadHitXZ->Draw("P");
    converter->grClusterXZ->Draw("P");
    grX[0]->Draw("P");
    grTemp->Draw("P");
    arc->Draw();
    arcThree->Draw();

    /*
    for( int i = 0; i< 4; i++){
      grX[i]->Draw("P");
    }
    */

    can->cd(2);
    gPad->DrawFrame(-300,-300,300,300);
    converter->grPadHitYZ->Draw("P");
    converter->grClusterYZ->Draw("P");
    /*
    for( int i = 0; i< 4; i++){
      grY[i]->Draw("P");
    }
    */
    grY[0]->Draw("P");

    can->cd(3);
    hisChi2Temp->DrawNormalized();
    hisDistTemp->DrawNormalized("same");

    can->cd(4);
    hisDistCircle->Draw();

    can->cd(5);
    hisDyCircle->Draw();

    can->cd(6);
    hisDDyCircle->Draw("colz");
    can->Update();
    can->Modified();
    gSystem->ProcessEvents();
    getchar();

  }
}
