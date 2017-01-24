void DrawSamples(){
  gROOT->ProcessLine("#include <vector>");
  gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  gSystem->Load("./lib/libanalysis_method");
  gStyle->SetOptStat(0);
  TFile* itf  = new TFile("HDTest_2lambda.root");
  TTree* itr  = (TTree*)itf->Get("eventTree00");
  TPCDataConverter* converter = new TPCDataConverter(itr);

  TCanvas* can = new TCanvas("can","can",800,1200);
  can->Divide(2,3);
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
  TArc* arcThree = new TArc(0,0,1);
  arcThree->SetLineColor(4);
  arcThree->SetFillStyle(0);
  arcThree->SetLineWidth(3);
  TPCCircleFitter* fitter  = new TPCCircleFitter();

  Int_t nTrack;
  Int_t eventList[3] = {2,4,6};

  TH2D* HisXZ;
  TH2D* hisYZ;
  //for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
  for( int ievt = 0; ievt < 3; ievt++){
    //for( int ievt = 0; ievt < 1; ievt++){
    //hisChi2Temp->Reset();
    nTrack = 0;
    grTemp->Set(0);
    //std::vector<TPCPadHit> HitArr = converter->Convert(ievt);
    std::cout<< "Convert" << std::endl;
    converter->Convert(eventList[ievt]);
    std::cout<< "Clustering " << std::endl;
    converter->Clustering();
    std::cout<< "End clustering " << std::endl;
    std::vector<TPCPadHitCluster> clusterArr = converter->GetClusterArr();
    if ( clusterArr.size() < 3 ){ continue; }
    std::cout<< "GetCluster" << std::endl;
    hisXZ = new TH2D(Form("hisXZ_%d",eventList[ievt]), Form("XZ distribution (Event #%d) ; Z [mm]; X [mm]",eventList[ievt]), 600,-300,300,600,-300,300);
    hisYZ = new TH2D(Form("hisXZ_%d",eventList[ievt]), Form("YZ distribution (Event #%d) ; Z [mm]; Y [mm]",eventList[ievt]), 600,-300,300,600,-300,300);
    std::vector<TPCPadHitCluster> cl[10];
    TGraph* grX[10];
    TGraph* grY[10];
    std::vector<double> xarr[10];
    std::vector<double> zarr[10];
    std::vector<double> yarr[10];
    std::vector<double> larr[10];
    TGraph* grXZ[10];
    TGraph* grYZ[10];
    TGraph* grYL[10];
    TGraph* grPadHitXZ;
    TGraph* grPadHitYZ;
    TPCCircleFitResult result[10];
    TArc* arc[10];


    //std::vector<TPCPadHitCluster> tempCluster = TPCFindEdgeBlock( clusterArr );
    Int_t nCenterCluster = 0;
    for( int icl = 0; icl < clusterArr.size(); icl++){
      if( clusterArr.at(icl).RowID() < 5 ){
	nCenterCluster++;
      }
    }
    std::vector<TPCPadHitCluster> tempCluster[10];
    Int_t nCluster = 0;
    //tempCluster[nCluster] = TPCFindBlockCircle( clusterArr,30 );
    //std::cout<< tempCluster[nCluster].size() << std::endl;

    while( clusterArr.size() > nCenterCluster && nCluster< 10 && clusterArr.size() > 3 ){
      Int_t nloop = 0;
      do{
	tempCluster[nCluster] = TPCFindBlockCircle( clusterArr, 30 );
	std::cout<<"Loop : " <<  nCluster << "\t" <<  tempCluster[nCluster].size() << std::endl;
	nloop = nloop +1;
	if( nloop > 100){ break; }
      }while( tempCluster[nCluster].size() <= 6 && tempCluster[nCluster].size() > 0 && clusterArr.size() > 3);
      if( nloop > 100 ){ break; }
      if( tempCluster[nCluster].size() != 0 ){
	nCluster = nCluster + 1;
      }
    }
    Int_t nClusterCheck = 0;
    for( int i = 0; i< nCluster; i++){
      if( tempCluster[i].size() > 6 ){ nClusterCheck++; }
    }
    nCluster = nClusterCheck;
    std::cout<< "Cluster Number : " << nCluster << std::endl;
    for( int i = 0; i< nCluster; i++){
      cl[i] = tempCluster[i];
      std::cout<< i << "\t" <<  tempCluster[i].size() << std::endl;
    }

    if( nCluster == 0 || nCluster == 10 ){ continue; }
    std::cout<< "Cluster Setting " << std::endl;
    std::vector<TPCPadHitCluster> NormalCluster = converter->GetClusterArr();

    for( int i = 0; i< nCluster; i++){
      grX[i] = DrawClusterArrXZ(cl[i]);
      grY[i] = DrawClusterArrYZ(cl[i]);
      grXZ[i] = new TGraph();
      grYZ[i] = new TGraph();
      grYL[i] = new TGraph();

      grX[i]->SetMarkerColor(i+1);
      grY[i]->SetMarkerColor(i+1);
      grX[i]->SetMarkerStyle(20+i);
      grY[i]->SetMarkerStyle(20+i);
      grXZ[i]->SetMarkerColor(i+1);
      grYZ[i]->SetMarkerColor(i+1);
      grXZ[i]->SetMarkerStyle(20+i);
      grYZ[i]->SetMarkerStyle(20+i);
      arc[i] = new TArc(0,0,1);
      arc[i]->SetLineColor(i+1);
      arc[i]->SetFillStyle(0);
      arc[i]->SetLineStyle(2);
      arc[i]->SetLineWidth(2);

      for( int icl = 0; icl < cl[i].size(); icl++){
	xarr[i].push_back( cl[i].at(icl).GetClusterPos().X());
	zarr[i].push_back( cl[i].at(icl).GetClusterPos().Z());
	yarr[i].push_back( cl[i].at(icl).GetClusterPos().Y());
	grXZ[i]->SetPoint( grXZ[i]->GetN(), cl[i].at(icl).GetClusterPos().Z(),cl[i].at(icl).GetClusterPos().X());
	grYZ[i]->SetPoint( grYZ[i]->GetN(), cl[i].at(icl).GetClusterPos().Z(),cl[i].at(icl).GetClusterPos().Y());
      }
      result[i] = fitter->ThreePointCircleFit(xarr[i],zarr[i]);
      std::cout << "Cluster : " << i << " status " << std::endl;
      for( int icl = 0; icl< cl[i].size(); icl++){
	std::cout<< icl << "\t" <<  cl[i].at(icl).RowID() << "\t" << cl[i].at(icl).ColID()  << "\t" << result[i].CalculateDist(cl[i].at(icl).GetClusterPos().X(), cl[i].at(icl).GetClusterPos().Z())<< std::endl;
      }
    }
    /// Print ///
    std::cout<< "Print" << std::endl;
    for( int i = 0; i< nCluster; i++){
      std::cout<< i << std::endl;
      std::cout<< "YL" << std::endl;
      for( int icl  = 0; icl < cl[i].size(); icl++){
	double theta = TMath::ATan2( xarr[i].at(icl) - result[i].GetFitPar().Y(),
				     zarr[i].at(icl) - result[i].GetFitPar().Z());
	double l     = theta * result[i].GetFitPar().X();
	larr[i].push_back( l );
	grYL[i]->SetPoint( grYL[i]->GetN(), l, yarr[i].at(icl));
      }
      grYL[i]->Fit("pol1");
      std::cout<< "FUNC" << std::endl;
      TF1* func = grYL[i]->GetFunction("pol1");
      double offset = func->GetParameter(0);
      double sloff  = func->GetParameter(1);
      result[i].SetDyDl( sloff );
      result[i].SetDyOff(offset);

      ///////// get remain clusters ///////////
      std::cout<< "Get remain cluster" << std::endl;
      for( int icl = 0; icl < NormalCluster.size(); icl++ ){
	double dist = result[i].CalculateDist( NormalCluster.at(icl).GetClusterPos().X(),
					       NormalCluster.at(icl).GetClusterPos().Z());
	double tempTheta = TMath::ATan2( NormalCluster.at(icl).GetClusterPos().X() - result[i].GetFitPar().Y(),
					 NormalCluster.at(icl).GetClusterPos().Z() - result[i].GetFitPar().Z());
	double l         = tempTheta * result[i].GetFitPar().X();
	double dy        = NormalCluster.at(icl).GetClusterPos().Y() - (result[i].GetDyDl() * l + result[i].GetDyOff());
	if( TMath::Abs( dist ) < 5  &&
	    TMath::Abs( dy )   < 10  ){
	  grXZ[i]->SetPoint(grXZ[i]->GetN(), NormalCluster.at(icl).GetClusterPos().Z(), NormalCluster.at(icl).GetClusterPos().X());
	  grYZ[i]->SetPoint(grYZ[i]->GetN(), NormalCluster.at(icl).GetClusterPos().Z(), NormalCluster.at(icl).GetClusterPos().Y());
	}
      }
      //result.GetFitPar().Print();
      //result.GetFitParErr().Print();
      arc[i]->SetX1( result[i].GetFitPar().Z());
      arc[i]->SetY1( result[i].GetFitPar().Y());
      arc[i]->SetR1( result[i].GetFitPar().X());
      arc[i]->SetR2( result[i].GetFitPar().X());
    }


    //if( nCluster < 4 ){ continue; }


    can->cd(1+ievt*2);

    //gPad->DrawFrame(-300,-300,300,300);
    hisXZ->GetYaxis()->SetTitleSize(0.05);
    hisXZ->GetYaxis()->SetLabelSize(0.04);
    hisXZ->GetYaxis()->SetTitleOffset(1.2);
    hisXZ->GetXaxis()->SetTitleSize(0.05);
    hisXZ->GetXaxis()->SetLabelSize(0.04);


    hisXZ->Draw();

    grPadHitXZ = (TGraph*)(converter->grPadHitXZ->Clone());
    grPadHitXZ->Draw("P");
    //converter->grClusterXZ->Draw("P");
    for( int i = 0; i< nCluster; i++){
      arc[i]->Draw();
    }
    for( int i = 0; i< nCluster; i++){
      grXZ[i]->Draw("P");
    }
    for( int i = 0; i< nCluster; i++){
      //grX[i]->Draw("P");
    }

    can->cd(2+ievt*2);
    hisYZ->GetYaxis()->SetTitleSize(0.05);
    hisYZ->GetYaxis()->SetLabelSize(0.04);
    hisYZ->GetYaxis()->SetTitleOffset(1.2);
    hisYZ->GetXaxis()->SetTitleSize(0.05);
    hisYZ->GetXaxis()->SetLabelSize(0.04);

    hisYZ->Draw();
    //gPad->DrawFrame(-300,-300,300,300);
    grPadHitYZ = (TGraph*)(converter->grPadHitYZ->Clone());
    grPadHitYZ->Draw("P");
    for( int i = 0; i< nCluster; i++){
      grYZ[i]->Draw("P");
      //grY[i]->Draw("P");
    }

    //converter->grClusterYZ->Draw("P");

    /*
    for( int i = 0; i< 4; i++){
      grY[i]->Draw("P");
    }
    */

    /*
    can->cd(3);
    //hisChi2Temp->DrawNormalized();
    //hisDistTemp->DrawNormalized("same");
    can->cd(4);
    hisDistCircle->Draw();

    can->cd(5);
    hisDyCircle->Draw();

    can->cd(6);
    hisDDyCircle->Draw("colz");
    can->Update();
    can->Modified();
    gSystem->ProcessEvents();
    std::cout<< ievt << std::endl;
    getchar();
    */
  }
}
