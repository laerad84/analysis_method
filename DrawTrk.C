/// Read track Data and Helix fitting

const double pimass = 139.57018;//MeV
const double pmass  = 938.272046;//MeV
const double kpmass = 493.667;//MeV
TGraph2D* grTrackSample;

void HelixCircle( Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t ){
  Int_t np = grTrackSample->GetN();
  f = 0;
  Double_t *x = grTrackSample->GetX();
  Double_t *y = grTrackSample->GetY();
  Double_t *z = grTrackSample->GetZ();
  Int_t nHit = 0;
  for( Int_t i = 0; i< np; i++){
    //par[0] : x center
    //par[1] : y Yoffset at theta == 0;
    //par[2] : z center
    //par[3] : r
    //par[4] : dY/dTheta

    Double_t u = x[i] - par[0];
    Double_t t = z[i] - par[1];
    Double_t dr = par[2] - TMath::Sqrt( u*u + t*t );
    if( dr > 20 ){ continue; }
    f += dr*dr;
    nHit++;
  }
}


void DrawTrk(){
  TString mSystemName = gSystem->GetFromPipe("uname");
  std::cout<< mSystemName << std::endl;
  if( mSystemName == "Darwin"){
    gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method.dylib");
    gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  }else if( mSystemName == "Linux"){
    gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method.so");
    gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  }

  TFile* itf = new TFile("Proton_trk.root");
  TTree* itr = (TTree*)itf->Get("convTree");

  TClonesArray* TrkArr = new TClonesArray("TGraph2D");
  itr->SetBranchAddress("TrkArr",&TrkArr);

  TFile* otf = new TFile("TestDist.root","recreate");
  TTree* otr = new TTree("Output","");
  Int_t    nComb;
  Int_t    nTrk;
  Double_t dY[400];
  Double_t dR[400];
  Int_t    nMerged;
  TClonesArray* MergedTrkArr = new TClonesArray("TGraph2D");
  TClonesArray* MergedSprialArr= new TClonesArray("HSprial");
  otr->Branch("nComb",&nComb,"nComb/I");
  otr->Branch("nTrk",&nTrk,"nTrk/I");
  otr->Branch("dY",dY,"dY[nComb]/D");
  otr->Branch("dR",dR,"dR[nComb]/D");
  otr->Branch("nMerged");
  otr->Branch("MergedTrkArr",&MergedTrkArr,256000,99);
  otr->Branch("MergedSprialArr",&MergedSprialArr,256000,99);
  TGraph* grX[20];
  for( int i = 0; i< 20; i++ ){
    grX[i] = new TGraph();
    grX[i]->SetMarkerStyle(20+i%4);
    grX[i]->SetMarkerColor(1+i%5);
  }

  TCanvas* can = new TCanvas("can","can",800,800);
  can->Divide(2,2);
  can->cd(1);
  gPad->DrawFrame(-300,-300,300,300);
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter *circleFitter = TVirtualFitter::Fitter(0,3);
  circleFitter->SetFCN(HelixCircle);

  for( int ievt = 0; ievt < itr->GetEntries(); ievt++){
  //for( int ievt = 0; ievt < 100; ievt++){
    itr->GetEntry(ievt);
    std::cout<< ievt << std::endl;
    for( int i = 0; i< 20; i++){
      grX[i]->Set(0);
    }
    MergedTrkArr->Clear();
    MergedSprialArr->Clear();

    if( TrkArr->GetEntries() > 20 ){ continue; }
    nTrk = TrkArr->GetEntries();
    TClonesArray* sprialArr = new TClonesArray("HSprial");
    std::cout<< "trk" << std::endl;
    for( int itrk = 0; itrk < TrkArr->GetEntries(); itrk++){
      std::cout<<"itrk "<<  itrk << std::endl;
      TGraph2D* itrack = (TGraph2D*)TrkArr->At(itrk);
      grTrackSample = itrack;
      for( int ip = 0; ip < itrack->GetN(); ip++){
	grX[itrk]->SetPoint( ip, itrack->GetZ()[ip],itrack->GetX()[ip]);
	//grTrackSample->SetPoint( ip, itrack->GetZ()[ip],itrack->GetX()[ip]);
      }
      std::cout << "Index " << std::endl;
      Int_t Index[3]={0,(int)((itrack->GetN()-1)/2),itrack->GetN() -1};
      TVector3 vec[3];

      for( int i = 0; i< 3; i++){
	vec[i].SetXYZ(itrack->GetX()[Index[i]],
		       itrack->GetY()[Index[i]],
		       itrack->GetZ()[Index[i]]);
      }
      std::cout<< "Sprial" << std::endl;
      HSprial sprial = GenerateSprial(vec[0],vec[1],vec[2]);
      sprial.ID = itrk;
      circleFitter->SetParameter(0,"X", sprial.X, 10, sprial.X-50, sprial.X+50 );
      circleFitter->SetParameter(1,"Z", sprial.Z, 10, sprial.Z-50, sprial.Z+50);
      circleFitter->SetParameter(2,"R", sprial.R, 10, sprial.R-50, sprial.R+50);
      Double_t arglist[1] = {0};
      circleFitter->ExecuteCommand("MIGRAD",arglist, 0 );
      sprial.X = circleFitter->GetParameter(0);
      sprial.Z = circleFitter->GetParameter(1);
      sprial.R = circleFitter->GetParameter(2);
      sprial.EX = circleFitter->GetParError(0);
      sprial.EZ = circleFitter->GetParError(1);
      sprial.ER = circleFitter->GetParError(2);
      sprial.InitPos = TVector3(itrack->GetX()[0],itrack->GetY()[0],itrack->GetZ()[0]);
      sprial.FinalPos= TVector3(itrack->GetX()[Index[2]],itrack->GetY()[Index[2]],itrack->GetZ()[Index[2]]);
      new((*sprialArr)[itrk]) HSprial(sprial);
    }
    std::cout<< "Combine Tracks" << std::endl;
    std::cout<< ievt << "\t" << TrkArr->GetEntries() << std::endl;
    nComb = 0;

    std::vector<int> MergedIndexArr;

    for( int itrk = 0; itrk < sprialArr->GetEntries(); itrk++){
      std::cout<< itrk << std::endl;
      Bool_t bMerged = false;
      for( int is = 0; is < MergedIndexArr.size(); is++){
	if( itrk == MergedIndexArr.at(is)){
	  bMerged = true;
	}
      }
      if( bMerged ){ continue; }

      Bool_t bMergeWork = false;
      HSprial* isprial = (HSprial*)sprialArr->At(itrk);
      TGraph2D* itrack = (TGraph2D*)TrkArr->At(itrk);
      if( itrack == NULL ){ std::cout<< "NULL" << std::endl; continue; }
      if( itrack->GetN() < 3 ){ continue; }
      TGraph2D* grTemp = new TGraph2D(*itrack);
      TGraph2D* grTest = new TGraph2D();
      grTest->SetPoint(0, itrack->GetX()[0],itrack->GetY()[0],itrack->GetZ()[0]);
      Bool_t bFind = false;
      Int_t index  = itrk+1;
      std::cout<< grTemp->GetN() << std::endl;
      std::vector<Int_t> connectedList;
      connectedList.push_back(itrk);
      for( int jtrk = itrk+1; jtrk <TrkArr->GetEntries(); jtrk++){
	TGraph2D* jtrack  = (TGraph2D*)TrkArr->At(jtrk);
	HSprial*  jsprial = (HSprial*)sprialArr->At(jtrk);
	if( jtrack->GetN() < 3 ){ continue; }
	Double_t r0,y0;
	Double_t r1,y1;
	isprial->CalculateDist(jtrack,r0,y0);
	jsprial->CalculateDist(itrack,r1,y1);
	if( (r0 < 10 && y0 < 20) ||
	    (r1 < 10 && y1 < 20) ){
	  std::cout<< "Merged" << std::endl;

	  connectedList.push_back(jtrk);
	  MergedIndexArr.push_back(jtrk);

	  Int_t initsize = grTemp->GetN();
	  for( int ip = 0; ip < jtrack->GetN(); ip++){
	    Int_t point= initsize + ip;
	    grTemp->SetPoint(point, jtrack->GetX()[ip],jtrack->GetY()[ip],jtrack->GetZ()[ip]);
	  }
	}
	nComb++;
      }




      std::cout<< "iend" << std::endl;
    }


    /*
    std::cout<< "Write track " << MergedIndexArr.size() <<  std::endl;
    nMerged = 0;
    for( int itrk = 0; itrk < TrkArr->GetEntries(); itrk++){
      Bool_t bMerged = false;
      for( int is = 0; is < MergedIndexArr.size(); is++){
	if( itrk == MergedIndexArr.at(is)){
	  bMerged = true;
	}
      }
      if( bMerged ){ continue; }
      std::cout<< itrk << "\t" << TrkArr->GetEntries() << std::endl;
      TGraph2D* itrack = (TGraph2D*)TrkArr->At(itrk);
      std::cout<< itrack << std::endl;
      if( itrack == NULL ){ continue; }
      if( itrack->GetN() < 3 ){ continue; }
      new((*MergedTrkArr)[nMerged]) TGraph2D(*itrack);
      nMerged++;
    }
    */

    /*
    can->cd(1);
    for( int itrk = 0; itrk < TrkArr->GetEntries(); itrk++){
      grX[itrk]->Draw("P");
    }
    can->Update();
    can->Modified();
    gSystem->ProcessEvents();
    getchar();
    */
    otr->Fill();
  }
  otr->Write();
  otf->Close();
}
