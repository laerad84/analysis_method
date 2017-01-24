void DataClusteringTest_Circle_1(){
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine("#include <algorithm>");
  gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  gSystem->Load("./lib/libanalysis_method");

  const double pimass = 139.57018;//MeV
  const double pmass  = 938.272046;//MeV


  TFile* itf  = new TFile("HDTest_2lambda.root");
  TTree* itr  = (TTree*)itf->Get("eventTree00");
  TPCDataConverter* converter = new TPCDataConverter(itr);
  converter->AddDetector("NBAR.");

  TCanvas* can = new TCanvas("can","can",2400,1600);
  can->Divide(3,2);
  can->Draw();


  TH1D* hisFit         = new TH1D("hisFit","hisFit",20,0,20);
  TH3D* his            = new TH3D("his","his;X;Y;Z",150,-300,300,150,-300,300,150,-300,300);
  TH2D* hisHit         = new TH2D("hisHit","his;X;Z",150,-300,300,150,-300,300);
  TPCPoly* polyHit     = new TPCPoly("polyHit","test");
  TPCPoly* polyCluster = new TPCPoly("polyCluster","test");
  TH1D* hisChi2Temp    = new TH1D("hisChi2Temp","hisChi2Temp",1000,0,100);
  TH1D* hisDistTemp    = new TH1D("hisDistTemp","hisDistTemp",1000,0,1000);
  TH1D* hisDistCircle  = new TH1D("hisDistCircle","hisDistCircle", 1000,0,100);
  TH1D* hisDyCircle    = new TH1D("hisDyCircle","hisDyCircle",1000,-50,50);
  TH2D* hisDDyCircle   = new TH2D("hisDDyCircle","hisDDyCircle",100,0,100,100,-50,50);
  TH1D* hisDR          = new TH1D("hisDR","hisDR",100,-50,50);
  TH1D* hisDY          = new TH1D("hisDY","hisDY",100,-50,50);
  TH2D* hisDYSlope     = new TH2D("hisDYSlope","hisDYSlope;dy;Slope", 100,-10,10,100,-5,5);
  TH1D* hisYDist       = new TH1D("hisYDist","hisYDist",400,-200,200);
  TH1D* hisMDist       = new TH1D("hisMDist","hisMDist",500,1000,1500);

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


  Int_t status = 0;
  Int_t nTrack;
  //for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
  for( int ievt = 0; ievt < 100; ievt++){
    status = 0;
    //for( int ievt = 0; ievt < 1; ievt++){
    //hisChi2Temp->Reset();
    nTrack = 0;
    grTemp->Set(0);
    //std::vector<TPCPadHit> HitArr = converter->Convert(ievt);
    std::cout<< "Convert" << std::endl;
    converter->Convert(ievt);
    std::cout<< "Clustering " << std::endl;
    converter->Clustering();
    std::cout<< "End clustering " << std::endl;

    std::cout<< "GetCluster" << std::endl;
    std::vector<TPCPadHitCluster> clusterArr = converter->GetClusterArr();
    std::cout<< "GetNBAR digi Data " << std::endl;
    TClonesArray* NBARArr = converter->GetDetDigi(0);
    std::cout<< "NBARSize: " << NBARArr->GetEntries() << std::endl;

    if( NBARArr->GetEntries() < 4 ){
      hisFit->Fill(status);
      continue;
    }else{ status++; }
    if( clusterArr.size() < 3){
      hisFit->Fill(status);
      continue;
    }else{ status++;}

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
    TPCCircleFitResult result[10];
    TArc* arc[10];

    /// Cluster finding
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

    if( nCluster == 0 || nCluster == 10 ){
      hisFit->Fill(status);
      continue;
    }else{
      status++;
    }

    std::cout<< "Cluster Setting " << std::endl;
    std::vector<TPCPadHitCluster> NormalCluster = converter->GetClusterArr();
    std::cout<< (*NormalCluster.begin()).RowID() << std::endl;


    ///// Sorting test /////
    /*
    std::cout<< "RowID array" << std::endl;
    for( int i = 0; i< NormalCluster.size(); i++){
      std::cout<< i << "\t" << NormalCluster.at(i).RowID() << std::endl;
    }
    std::cout<< "Sorting" << std::endl;
    std::sort(NormalCluster.begin(), NormalCluster.end());
    for( int i = 0; i< NormalCluster.size(); i++){
      std::cout<< i << "\t" << NormalCluster.at(i).RowID() << std::endl;
    }
    getchar();
    continue;
    */

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
	std::cout<< icl << "\t"
		 << cl[i].at(icl).RowID() << "\t"
		 << cl[i].at(icl).ColID() << "\t"
		 << result[i].CalculateDist(cl[i].at(icl).GetClusterPos().X(), cl[i].at(icl).GetClusterPos().Z())<< std::endl;
      }
    }

    Int_t nFitTrack = 0;

    for( int i = 0; i< nCluster; i++){
      if( grXZ[i]->GetN() >= 6 ){
	nFitTrack++;
      }
    }

    if( nFitTrack != 4){
      hisFit->Fill(status);
    }else{
      status++;
      hisFit->Fill( status);
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
      double slope  = func->GetParameter(1);
      result[i].SetDyDl(slope);
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
	  hisDR->Fill(dist);
	  hisDY->Fill(dy);
	  hisDYSlope->Fill(dy,slope);
	}
      }
      //result.GetFitPar().Print();
      //result.GetFitParErr().Print();
      arc[i]->SetX1( result[i].GetFitPar().Z());
      arc[i]->SetY1( result[i].GetFitPar().Y());
      arc[i]->SetR1( result[i].GetFitPar().X());
      arc[i]->SetR2( result[i].GetFitPar().X());
    }



    for( int i = 0; i< nCluster; i++){
      std::cout<< "Clustering result : " << i << std::endl;
      result[i].GetFitPar().Print();
      std::cout<< "Momentum : " << result[i].GetFitPar().Print()*1.0/3.36 << " MeV" << std::endl;//MeV
      result[i].SetID(i);
      result[i].SetMinimumY(10000);
      result[i].SetCoupleID(-1);
    }

    std::cout<< "Crossing Point " << std::endl;
    for( int i = 0; i< nCluster; i++){
      double dist = 0;
      result[i].SetMinimumY( 100000 );
      for( int j = i+1; j< nCluster; j++){
	TVector3 pos;
	bool rst = yDistCalculator( result[i], result[j], dist, pos );
	if( !rst ){ continue; }
	if( TMath::Abs(dist) > 50 ){ continue; }
	std::cout<< dist << std::endl;
	if( TMath::Abs(dist) < result[i].GetMinimumY() ){
	  result[i].SetCoupleID(j);
	  result[i].SetMinimumY(TMath::Abs(dist));
	  result[i].SetStartingPoint( pos );
	  result[i].CalculateStartingMomentum(1.5);
	  result[j].SetCoupleID(i);
	  result[j].SetMinimumY(TMath::Abs(dist));
	  result[j].SetStartingPoint( pos );
	  result[j].CalculateStartingMomentum(1.5);

	}
	pos.Print();
	hisYDist->Fill(dist);
      }
    }

    TMarker* mark[100];
    Int_t    markerIndex = 0;
    //// Old crossing point ////
    /*
    for( int i = 0; i< nCluster; i++){
      for( int j = i+1; j < nCluster; j++){
	TVector2 CrossingXZ[2];
	bool bCrossing  = CircleCrossing( result[i].GetFitPar(), result[j].GetFitPar(), CrossingXZ[0], CrossingXZ[1]);
	if( !bCrossing ){ continue; }else{
	  mark[markerIndex] = new TMarker( CrossingXZ[0].Y(), CrossingXZ[0].X(),22);
	  mark[markerIndex]->SetMarkerSize(3);
	  markerIndex = markerIndex +1;
	  mark[markerIndex] = new TMarker( CrossingXZ[1].Y(), CrossingXZ[1].X(),22);
	  mark[markerIndex]->SetMarkerSize(3);
	  markerIndex = markerIndex +1;
	}
      }
      }*/

    /////////////////////////////
    ///// New crossing point ////
    for( int i = 0; i< nCluster; i++){
      if( result[i].GetCoupleID() < 0){ continue; }
      if( result[i].GetID() == result[result[i].GetCoupleID()].GetCoupleID()){
	mark[markerIndex] = new TMarker(result[i].GetStartingPoint().Z(),result[i].GetStartingPoint().X(),22);
	mark[markerIndex]->SetMarkerSize(3);
	markerIndex++;


	TVector3 vec0[2];
	TVector3 vec1[2];
	/// i = pi j = proton
	result[i].GetStartingMomentum(vec0[0],vec1[1]);
	result[j].GetStartingMomentum(vec1[0],vec1[1]);
	/// combination -charge
	/// (0[0], 1[1]), (0[1],1[0]);
	//// pion - proton

	TLorentzVector lv[2];
	Double_t e[2] ={0};
	Double_t mass = 0;
	TLorentzVector s;
	e[0] = TMath::Sqrt( vec0[0].Mag2() + pimass*pimass);
	e[1] = TMath::Sqrt( vec1[1].Mag2() + pmass*pmass);
	lv[0].SetPxPyPzE(vec0[0].X(),vec0[0].Y(),vec0[0].Z(),e[0]);
	lv[1].SetPxPyPzE(vec1[1].X(),vec1[1].Y(),vec1[1].Z(),e[1]);
	s = (lv[0]+lv[1]);
	mass = s.Mag();
	hisMDist->Fill(mass);

	e[0] = TMath::Sqrt( vec0[0].Mag2() + pmass*pmass);
	e[1] = TMath::Sqrt( vec1[1].Mag2() + pimass*pimass);
	lv[0].SetPxPyPzE(vec0[0].X(),vec0[0].Y(),vec0[0].Z(),e[0]);
	lv[1].SetPxPyPzE(vec1[1].X(),vec1[1].Y(),vec1[1].Z(),e[1]);
	s = (lv[0]+lv[1]);
	mass = s.Mag();
	hisMDist->Fill(mass);

	e[0] = TMath::Sqrt( vec0[1].Mag2() + pimass*pimass);
	e[1] = TMath::Sqrt( vec1[0].Mag2() + pmass*pmass);
	lv[0].SetPxPyPzE(vec0[1].X(),vec0[1].Y(),vec0[1].Z(),e[0]);
	lv[1].SetPxPyPzE(vec1[0].X(),vec1[0].Y(),vec1[0].Z(),e[1]);
	s = (lv[0]+lv[1]);
	mass = s.Mag();
	hisMDist->Fill(mass);

	e[0] = TMath::Sqrt( vec0[1].Mag2() + pmass*pmass);
	e[1] = TMath::Sqrt( vec1[0].Mag2() + pimass*pimass);
	lv[0].SetPxPyPzE(vec0[1].X(),vec0[1].Y(),vec0[1].Z(),e[0]);
	lv[1].SetPxPyPzE(vec1[0].X(),vec1[0].Y(),vec1[0].Z(),e[1]);
	s = (lv[0]+lv[1]);
	mass = s.Mag();
	hisMDist->Fill(mass);

      }
    }




    //if( nCluster < 4 ){ continue; }


    can->cd(1);
    gPad->DrawFrame(-300,-300,300,300);
    converter->grPadHitXZ->Draw("P");
    converter->grClusterXZ->Draw("P");
    for( int i = 0; i< nCluster; i++){
      arc[i]->Draw();
    }
    for( int i = 0; i< nCluster; i++){
      grXZ[i]->Draw("P");
    }
    for( int i = 0; i< nCluster; i++){
      grX[i]->Draw("P");
    }

    for( int i = 0; i< markerIndex; i++){
      mark[i]->Draw();
    }

    can->cd(2);
    gPad->DrawFrame(-300,-300,300,300);
    converter->grPadHitYZ->Draw("P");
    converter->grClusterYZ->Draw("P");
    for( int i = 0; i< nCluster; i++){
      grYZ[i]->Draw("P");
      grY[i]->Draw("P");
    }
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

  }
  can->cd(1);
  hisFit->Draw();
  can->cd(2);
  hisDR->Draw();
  can->cd(3);
  hisDY->Draw();
  can->cd(4);
  hisDYSlope->Draw("colz");
  can->cd(5);
  hisYDist->Draw();
  can->cd(6);
  hisMDist->Draw();
}
