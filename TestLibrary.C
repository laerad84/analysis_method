void TestLibrary(){
  TString mSystemName = gSystem->GetFromPipe("uname");
  std::cout<< mSystemName << std::endl;
  if( mSystemName == "Darwin"){
    gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method.dylib");
  }else if( mSystemName == "Linux"){
    gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method.so");
  }

  TFile* itf = new TFile("HDE42_All.root");
  TTree* itr = (TTree*)itf->Get("convTree");
  TClonesArray* HitArr = new TClonesArray("TPCPadHit");
  GsimGenParticleData*   particle = new GsimGenParticleData();
  GsimDetectorEventData* tpcData  = new GsimDetectorEventData();
  GsimDetectorEventData* ftofData = new GsimDetectorEventData();
  GsimDetectorEventData* DC1Data  = new GsimDetectorEventData();
  GsimDetectorEventData* DC2Data  = new GsimDetectorEventData();
  GsimDetectorEventData* DC3Data  = new GsimDetectorEventData();
  GsimDetectorEventData* NBARData = new GsimDetectorEventData();
  itr->SetBranchAddress("HitArr",&HitArr);
  itr->SetBranchAddress("GenParticle.",&particle);
  //itr->SetBranchAddress("TPC.",&tpcData);
  itr->SetBranchAddress("FTOF.",&ftofData);
  itr->SetBranchAddress("DC1.",&DC1Data);
  itr->SetBranchAddress("DC2.",&DC2Data);
  itr->SetBranchAddress("DC3.",&DC3Data);
  itr->SetBranchAddress("NBAR.",&NBARData);

  TH1D* hisDist = new TH1D("hisDist","hisDist",400,0,400);
  TH2D* hisCross= new TH2D("hisCross","hisCross",300,-300,300,300,-300,300);
  TH3D* hisBox = new TH3D("hisBox","hisBox;X;Y;Z",30,-300,300,30,-300,300,30,-300,300);
  TCanvas* can = new TCanvas("can","can",2000,1000);
  //TCanvas* can = new TCanvas("can","can",1200,600);
  //can->Divide(2,1);
  //can->cd(1);
  //gPad->DrawFrame(-300,-300,300,300);
  //can->cd(2);
  can->Divide(2,1);
  can->cd(1);
  gPad->DrawFrame(-300,-300,300,300);
  can->cd(2);
  hisBox->Draw();

  const Int_t maxNClusterBlock = 50;
  const Int_t maxNTrackCand    = 50;
  TGraph* gr = new TGraph();
  gr->SetMarkerStyle(21);
  for( int ievt = 0; ievt < itr->GetEntries() ;ievt++){
    itr->GetEntry( ievt );
    std::cout<< "Clustering" << std::endl;
    gr->Set(0);

    //std::vector<TPCCluster> ClArr;
    TClonesArray* ClArr = new TClonesArray("TPCCluster");
    HitClustering( HitArr, ClArr );
    for( int i = 0; i < ClArr.size(); i++){
      std::cout<< ClArr.at(i).ID << "\t"
	       << ClArr.at(i).Mother << "\t"
	       << ClArr.at(i).Daughter << "\t"
	       << ClArr.at(i).nMother << "\t"
	       << ClArr.at(i).nDaughter << std::endl;
    }


    for( int i = 0; i< ClArr.size(); i++){
      gr->SetPoint( i, ClArr.at(i).Position.Z(),ClArr.at(i).Position.X());
    }
    std::cout << "blocking" << std::endl;
    //ConnectionFinder( ClArr );

    std::vector< std::vector<TPCCluster> > clBlockArr = ConnectionBlocker( ClArr );

    /// Connection Blocker
    std::cout<< clBlockArr.size() << std::endl;


    getchar();
    std::cout<<"Number of Hit     : " << HitArr.GetEntries() << "\n"
	     <<"Number of Cluster : " << ClArr.size() << "\n"
	     << "End" << std::endl;

    if( ClArr.size() > 200 ){ continue; }
    TPolyMarker3D* ClEvent = GenerateClusterView( ClArr );


    std::cout << "Divide into Block" << std::endl;
    std::vector<TPCCluster> Block[maxNClusterBlock];
    std::vector<TPCCluster> TrackCand[maxNTrackCand];
    /// divide into block ///
    Int_t BlockIndex     = 0;
    Int_t TrackCandIndex = 0;
    while( ClArr.size() != 0){
      Block[BlockIndex] = TPCClusterBlocker(ClArr,0);
      BlockIndex++;
      if( BlockIndex >= maxNClusterBlock){
	break;
      }
    }
    if( BlockIndex == maxNClusterBlock){
      std::cout<< "Number of Block reached to the limit, 50 " << std::endl;
    }
    for( int iBlock =0; iBlock < BlockIndex; iBlock++){
      while( Block[iBlock].size() > 3 ){
	Int_t CandSize = BlockDivider( Block[iBlock], TrackCand[TrackCandIndex],0);//root to Center
	if( CandSize == 0){
	  CandSize = BlockDivider( Block[iBlock], TrackCand[TrackCandIndex],1);//Branch to center
	}
	if( CandSize == 0 ){ break; }
	TrackCandIndex++;
	if( TrackCandIndex >= maxNTrackCand ){ break; }
      }
      if( TrackCandIndex >= maxNTrackCand ){ break; }
    }
    if( TrackCandIndex >= maxNTrackCand){
      std::cout<< "Number of Block reached to the limit, 50" << std::endl;
    }
    std::cout<< "Number of Block     : " << BlockIndex << "\n"
	     << "Number of TrackCand : " << TrackCandIndex << "\n"
	     << std::endl;

    std::cout << "Sprial Fitting " << std::endl;
    Int_t HelixIndex = 0;
    //std::vector<HSprial> SprialArr;
    std::vector<HSprial>      PSprialArr;
    std::vector<HSprial>      MSprialArr;
    std::vector<TPolyLine3D*> helixArr;
    std::vector<TArc*>        arcArr;
    /*
    for( int iCand = 0; iCand < TrackCandIndex; iCand++){
      if( TrackCand[iCand].size() < 5 ){ continue; }
      Int_t firstIndex  = 0;
      Int_t middleIndex = (int)((TrackCand[iCand].size() -1)/2.);
      Int_t lastIndex   = TrackCand[iCand].size() -1;
      if( TrackCand[iCand].size() >= 6 ){ lastIndex = TrackCand[iCand].size()-2;}
      if( TrackCand[iCand].size() >= 7 ){ firstIndex = 1; }
      if( TrackCand[iCand].size() >= 8 ){ lastIndex = TrackCand[iCand].size()-3;}
      HSprial sprial;
      if( TrackCand[iCand].at(firstIndex).Row < TrackCand[iCand].at(lastIndex).Row ){
	sprial = GenerateSprial( TrackCand[iCand].at(firstIndex).Position,
				 TrackCand[iCand].at(middleIndex).Position,
				 TrackCand[iCand].at(lastIndex).Position);
      }else{
	sprial = GenerateSprial( TrackCand[iCand].at(lastIndex).Position,
				 TrackCand[iCand].at(middleIndex).Position,
				 TrackCand[iCand].at(firstIndex).Position);
      }
      sprial.ID = HelixIndex;
      TPolyLine3D* helix  = sprial.GenerateHelix();
      TArc*        arc    = sprial.GenerateArc();
      //SprialArr.push_back(sprial);
      if( sprial.RL > 0 ){
	PSprialArr.push_back( sprial );
      }else{
	MSprialArr.push_back( sprial );
      }
      helixArr.push_back(helix);
      arcArr.push_back(arc);
      HelixIndex++;

    }
    */
    for( int ip = 0; ip < PSprialArr.size(); ip++){
      for( int im = 0; im < MSprialArr.size(); im++){
	Double_t dist;
	TVector3 pos(0,0,0);
	Bool_t cross = CalculateCrossing( PSprialArr.at(ip), MSprialArr.at(im), dist, pos);
	if( !cross ){ continue; }
	hisDist->Fill( dist );
	hisCross->Fill(pos.Z(),pos.X());
      }
    }

    can->cd(1);
    gr->Draw("P");
    for( Int_t iarc = 0; iarc < arcArr.size(); iarc++){
      TArc* arc = arcArr.at(iarc);
      arc->SetFillStyle(0);
      arc->SetLineColor(iarc+1);
      arc->Draw();
    }

    can->cd(2);
    TPolyMarker3D* crossingPoints = GenerateCrossingPoints(PSprialArr, MSprialArr);
    crossingPoints->SetMarkerStyle(22);
    crossingPoints->SetMarkerColor(2);
    ClEvent->SetMarkerStyle(21);

    ClEvent->Draw();
    crossingPoints->Draw();
    for( Int_t isprial = 0; isprial < helixArr.size(); isprial++ ){
      TPolyLine3D* line = helixArr.at(isprial);
      line->SetLineColor( isprial+1);
      line->SetLineWidth(2);
      line->Draw();
    }

    can->Update();
    can->Modified();
    gSystem->ProcessEvents();
    getchar();
    delete ClEvent;
    delete crossingPoints;
    for( Int_t isprial = 0; isprial < helixArr.size(); isprial++ ){
      TPolyLine3D* line = helixArr.at(isprial);;
      TArc* arc = arcArr.at(isprial);
      delete line;
      delete arc;
    }

  }
  hisCross->Draw("colz");
}
