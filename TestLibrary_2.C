#include "TVirtualFitter.h"

const double pimass = 139.57018;//MeV
const double pmass  = 938.272046;//MeV
const double kpmass = 493.667;//MeV
TGraph2D* grTrackSample;
void HelixY( Int_t &, Double_t *, Double_t &f, Double_t *par , Int_t ){
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
    Double_t t = z[i] - par[2];
    Double_t theta = TMath::ATan2( u,t);
    Double_t expY = theta*par[4] + par[1];
    Double_t v = y[i] - expY;
    Double_t dr = par[3] - TMath::Sqrt( u*u + t*t );
    //if( dr > 20 ){ continue; }
    f += v*v;
    nHit++;
  }
}


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
    //if( dr > 20 ){ continue; }
    f += dr*dr;
    nHit++;
  }
}


void TestLibrary_2(){/// reconstruction ///
  TString mSystemName = gSystem->GetFromPipe("uname");
  std::cout<< mSystemName << std::endl;
  if( mSystemName == "Darwin"){
    gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method.dylib");
    gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  }else if( mSystemName == "Linux"){
    gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method.so");
    gSystem->Load("~/local/hep/E42/E42/lib/so/libGsimData.dylib");
  }

  //GenFitter* fitter = new GenFitter("Data/TPC.root");
  //TVirtualFitter::SetDefaultFitter("Minuit");
  grTrackSample = new TGraph2D();
  //fitter->SetField("~/local/hep/E42/E42/data/field/SC_FieldMap.root",TVector3(0,0,0),TVector3(0,0,0));
  //TFile* itf = new TFile("HDE42_conv.root");
  //TTree* itr = (TTree*)itf->Get("convTree");
  //TFile* itf = new TFile("TestHit.root");
  //TTree* itr = (TTree*)itf->Get("padHit");

  TChain* itr = new TChain("convTree");
  //itr->Add("MCData/HDE42_conv_*.root");
  //itr->Add("HD_0_conv.root");
  itr->Add("lambda_Aconv.root");
  //itr->Add("Proton_conv.root");
  TClonesArray* HitArr = new TClonesArray("TPCPadHit");
  GsimGenParticleData*   particle = new GsimGenParticleData();

  grTrackSample = new TGraph2D();

  /*
  GsimDetectorEventData* tpcData  = new GsimDetectorEventData();
  GsimDetectorEventData* ftofData = new GsimDetectorEventData();
  GsimDetectorEventData* DC1Data  = new GsimDetectorEventData();
  GsimDetectorEventData* DC2Data  = new GsimDetectorEventData();
  GsimDetectorEventData* DC3Data  = new GsimDetectorEventData();
  GsimDetectorEventData* NBARData = new GsimDetectorEventData();
  */
  itr->SetBranchAddress("HitArr",&HitArr);
  itr->SetBranchAddress("GenParticle.",&particle);
  //itr->SetBranchAddress("TPC.",&tpcData);
  /*
  itr->SetBranchAddress("FTOF.",&ftofData);
  itr->SetBranchAddress("DC1.",&DC1Data);
  itr->SetBranchAddress("DC2.",&DC2Data);
  itr->SetBranchAddress("DC3.",&DC3Data);
  itr->SetBranchAddress("NBAR.",&NBARData);
  */
  TFile* otf = new TFile("TestFit_Alambda.root","recreate");
  TTree* otr = new TTree("anaTree","Test");
  TClonesArray* pSprialArr = new TClonesArray("HSprial");
  TClonesArray* mSprialArr = new TClonesArray("HSprial");
  TClonesArray* pMomentumArr = new TClonesArray("TLorentzVector");
  TClonesArray* mMomentumArr = new TClonesArray("TLorentzVector");
  TClonesArray* lMomentumArr = new TClonesArray("TLorentzVector");
  TClonesArray* PositionArr  = new TClonesArray("TVector3");
  TClonesArray* ClArr = new TClonesArray("TPCCluster");
  Int_t         EventNumber;
  Int_t         nPTrack;
  Int_t         nMTrack;
  Int_t         nLambda;
  Double_t      Distance[100];
  TVector3      pMom;
  TVector3      piMom;
  TVector3      lMom;
  TVector3      lPos;
  //otr->Branch("GenParticle.",&particle,25600,99);
  otr->Branch("HitArr",&HitArr,25600,99);
  otr->Branch("pSprialArr",&pSprialArr,25600,99);
  otr->Branch("mSprialArr",&mSprialArr,25600,99);
  otr->Branch("EventNumber",&EventNumber,"EventNumber/I");
  otr->Branch("GenParticle.",&particle);
  otr->Branch("pMomentumArr",&pMomentumArr,256000,99);
  otr->Branch("mMomentumArr",&mMomentumArr,256000,99);
  otr->Branch("lMomentumArr",&lMomentumArr,256000,99);
  otr->Branch("PositionArr",&PositionArr,256000,99);
  otr->Branch("nPTrack",&nPTrack,"nPTrack/I");
  otr->Branch("nMTrack",&nMTrack,"nMTrack/I");
  otr->Branch("ClArr",&ClArr,256000,99);
  otr->Branch("lMom",&lMom);
  otr->Branch("lPos",&lPos);
  otr->Branch("pMom",&pMom);
  otr->Branch("piMom",&piMom);
  otr->Branch("nLambda",&nLambda,"nLambda/I");
  otr->Branch("Distance",Distance,"Distance[nLambda]/D");


  TH1D* hislMass = new TH1D("hislMass","hislMass",3000,500,3500);
  TH2D* hisClusterSizeEne  = new TH2D("hisClustersizeEne","hisCLustersizeEne",30,0,30,100,0,1000);
  TH1D* hisDist = new TH1D("hisDist","hisDist",400,0,400);
  TH2D* hisCross= new TH2D("hisCross","hisCross",300,-300,300,300,-300,300);
  TH3D* hisBox = new TH3D("hisBox","hisBox;X;Y;Z",30,-300,300,30,-300,300,30,-300,300);
  TH2D* hisMag = new TH2D("hisMag","hisMag",100,0,2000,100,0.5,1.5);
  TCanvas* can = new TCanvas("can","can",1200,1200);
  //TCanvas* can = new TCanvas("can","can",1200,600);
  //can->Divide(2,1);
  //can->cd(1);
  //gPad->DrawFrame(-300,-300,300,300);
  //can->cd(2);
  can->Divide(2,2);
  can->cd(1);
  gPad->DrawFrame(-300,-300,300,300);
  can->cd(2);
  hisBox->Draw();
  can->cd(3);
  hisBox->Draw();
  //TFile* otf = new TFile("Output.root","recreate");

  const Int_t maxNBlock = 50;
  const Int_t maxNTrackCand    = 50;
  TGraph* grCandidates[maxNBlock];
  for( int i = 0; i< maxNBlock; i++){
    grCandidates[i] = new TGraph();
    grCandidates[i]->SetMarkerStyle( 20+i%8);
    grCandidates[i]->SetMarkerColor( i/7 +1);
  }

  pMom  = TVector3(0,0,0);
  piMom = TVector3(0,0,0);
  lMom  = TVector3(0,0,0);
  lPos  = TVector3(0,0,0);

  for( int ievt = 0; ievt <itr->GetEntries() ;ievt++){
  //for( int ievt = 0; ievt < 1000;ievt++){
    if( ievt  == 781 ){ continue;}
    std::cout<< "*************************************************"<< std::endl;
    std::cout<< "EventNumber : " << ievt << std::endl;
    std::cout<< "*************************************************"<< std::endl;
    EventNumber = ievt;
    ClArr->Clear();
    itr->GetEntry( ievt );

    Int_t iLambdaTrack=-1;
    TClonesArray* trackArr = particle->briefTracks;
    for( int itrack = 0; itrack < trackArr->GetEntries(); itrack++){
      GsimTrackData* gsimtrack  = (GsimTrackData*)trackArr->At(itrack);
      if( gsimtrack->mother == -1 && gsimtrack->pid == 3122){
	iLambdaTrack = gsimtrack->track;
	lMom.SetXYZ(gsimtrack->p.X(),gsimtrack->p.Y(),gsimtrack->p.Z());
	lPos.SetXYZ(gsimtrack->end_v.X(),gsimtrack->end_v.Y(),gsimtrack->end_v.Z());
      }
    }
    for( int itrack = 0; itrack < trackArr->GetEntries(); itrack++){
      GsimTrackData* gsimtrack  = (GsimTrackData*)trackArr->At(itrack);
      if( gsimtrack->mother == iLambdaTrack){
	if ( gsimtrack->pid == 2212){
	  pMom.SetXYZ(gsimtrack->p.X(),gsimtrack->p.Y(),gsimtrack->p.Z());
	}else if( gsimtrack->pid == -211 ){
	  piMom.SetXYZ(gsimtrack->p.X(),gsimtrack->p.Y(),gsimtrack->p.Z());
	}
      }
    }


    TPolyMarker3D* total = new TPolyMarker3D();
    total->SetMarkerStyle(20);
    TPolyMarker3D* track[maxNBlock];
    for( int i = 0; i< maxNBlock; i++){
      track[i] = new TPolyMarker3D();
      track[i]->SetMarkerStyle( 20+ i%8 );
      track[i]->SetMarkerColor( i/7 +1 );
    }
    std::cout<< "Clustering" << std::endl;
    for( int igr =0; igr < maxNBlock; igr++){
      grCandidates[igr]->Set(0);
    }

    Bool_t bCluster = HitClustering( HitArr, ClArr );
    if( !bCluster ){
      std::cout<< "*************************" << std::endl;
      std::cout<< "Clustering : No Cluster "  << std::endl;
      std::cout<< "*************************" << std::endl;
      continue;
    }
    if( ClArr->GetEntries() > 200 ){
      std::cout<< "*************************" << std::endl;
      std::cout<< "Clustering : Too many entries"  << std::endl;
      std::cout<< "*************************" << std::endl;
      continue;
    }

    std::cout<< "ClArr : " << ClArr->GetEntries() << std::endl;

    for( int i = 0; i < ClArr->GetEntries(); i++){
      TPCCluster* cl = (TPCCluster*)(ClArr->At(i));
      hisClusterSizeEne->Fill(cl->NHit, cl->Energy);
      //cl->Print();
      /*
      std::cout<< i <<"\t" <<  cl->ID << "\t"
	       << cl->MotherID.size() << "\t" << cl->DaughterID.size() << "\t"
	       << cl->Mother          << "\t" << cl->Daughter << "\t"
	       << std::endl;
      */
      total->SetNextPoint( cl->Position.X(), cl->Position.Y(), cl->Position.Z());
    }


    //std::cout<< "GetList " << std::endl;
    std::vector<Int_t> RootIDList = GetListOfTrackRoot( ClArr );
    //std::cout<<"Root ID list size : " << RootIDList.size() << std::endl;
    if( RootIDList.size() > maxNBlock ){ std::cout << "too many blocks " << std::endl; }
    Int_t nCand = 0;
    Int_t nHelix = 0;
    pSprialArr->Clear();
    mSprialArr->Clear();
    pMomentumArr->Clear();
    mMomentumArr->Clear();
    lMomentumArr->Clear();

    std::vector<TPolyLine3D*> HelixArr;
    std::vector<TPolyLine*>   ArcArr;

    Int_t nPsp = 0;
    Int_t nMsp = 0;
    Int_t nMTrack = 0;
    Int_t nPTrack = 0;
    for( Int_t arrIndex = 0; arrIndex < RootIDList.size(); arrIndex++){
      TClonesArray* subClArr = new TClonesArray("TPCCluster");
      Bool_t bResult = ClusterBlocker( ClArr, subClArr,RootIDList, RootIDList.at(arrIndex));
      if( !bResult){ continue; }
      for( Int_t index = 0; index < subClArr->GetEntries(); index++){
	TPCCluster* subCl = (TPCCluster*)subClArr->At(index);
	/*
	std::cout<< subCl->Mother       << "\t" << subCl->ID << "\t" << subCl->Daughter << "\t"
		 << subCl->nMother      << "\t" << subCl->nDaughter << "\t"
		 << subCl->Row          << "\t" << subCl->Col << "\t"
		 << subCl->Position.X() << "\t" << subCl->Position.Z() << "\t"
		 << std::endl;
	*/
	grCandidates[nCand]->SetPoint( index, subCl->Position.Z(), subCl->Position.X());
	track[nCand]->SetNextPoint(subCl->Position.X(), subCl->Position.Y(), subCl->Position.Z());
      }

      if( subClArr->GetEntries() >= 4 ){
	grTrackSample->Set(0);
	for( Int_t ip = 0; ip < subClArr->GetEntries(); ip++){
	  TPCCluster* clref = (TPCCluster*)subClArr->At(ip);
	  grTrackSample->SetPoint( ip ,clref->Position.X(),
				   clref->Position.Y(),
				   clref->Position.Z());
	}


	Int_t firstIndex  = 0;
	Int_t middleIndex = (int)((subClArr->GetEntries()-1)/2.);
	Int_t finalIndex  =  subClArr->GetEntries() -1;
	if( ((TPCCluster*)(subClArr->At(finalIndex)))->Row == 31 ){
	  finalIndex = subClArr->GetEntries() -2;
	}
	HSprial sprial = GenerateSprial(((TPCCluster*)(subClArr->At(firstIndex)))->Position,
					((TPCCluster*)(subClArr->At(middleIndex)))->Position,
					((TPCCluster*)(subClArr->At(finalIndex)))->Position);
	/*
	TVirtualFitter::SetDefaultFitter("Minuit");
	TVirtualFitter *mfitter = TVirtualFitter::Fitter(0,5);
	mfitter->SetFCN(HelixFunc);
	std::cout<< "test" << std::endl;
	mfitter->SetParameter(0,"X", sprial.X, 20, 0,0 );
	Double_t InitTheta= TMath::ATan2( sprial.InitPos.X() - sprial.X, sprial.InitPos.Z() - sprial.Z );
	mfitter->SetParameter(1,"Y", sprial.InitPos.Y() - InitTheta*sprial.DY/sprial.DTheta, 10, 0,0);
	mfitter->SetParameter(2,"Z", sprial.Z, 20, 0, 0);
	mfitter->SetParameter(3,"R", sprial.R, 20, 0, 0);
	mfitter->SetParameter(4,"DYDT", sprial.DY/sprial.DTheta, 20, 0, 0);
	std::cout<< "test" << std::endl;
	Double_t arglist[1] = {0};
	mfitter->ExecuteCommand("MIGRAD",arglist, 0 );
	std::cout<< nHelix << "\t"
		 << mfitter->GetParameter(0) << "\t"
		 << mfitter->GetParameter(1) << "\t"
		 << mfitter->GetParameter(2) << "\t"
		 << mfitter->GetParameter(3) << "\t"
		 << mfitter->GetParameter(4) << "\t"
		 << std::endl;
	*/



	TVirtualFitter::SetDefaultFitter("Minuit");
	TVirtualFitter *circleFitter = TVirtualFitter::Fitter(0,3);
	circleFitter->SetFCN(HelixCircle);
	circleFitter->SetParameter(0,"X", sprial.X, 10, sprial.X-50, sprial.X+50 );
	circleFitter->SetParameter(1,"Z", sprial.Z, 10, sprial.Z-50, sprial.Z+50);
	circleFitter->SetParameter(2,"R", sprial.R, 10, sprial.R-50, sprial.R+50);
	Double_t arglist[1] = {0};
	circleFitter->ExecuteCommand("MIGRAD",arglist, 0 );

	std::cout<< "**************************************" << std::endl;
	std::cout<< sprial.R << "\t" << circleFitter->GetParameter(2) << std::endl;
	std::cout<< sprial.X << "\t" << circleFitter->GetParameter(0) << std::endl;
	std::cout<< sprial.Z << "\t" << circleFitter->GetParameter(1) << std::endl;
	std::cout<< "**************************************" << std::endl;
	sprial.X = circleFitter->GetParameter(0);
	sprial.Z = circleFitter->GetParameter(1);
	sprial.R = circleFitter->GetParameter(2);


	Double_t AbsoluteMomentum = CalculateMomentum( sprial, 1.);
	TVector3 tangentialVec(0,0,0);
	CalculateCircleTangent(sprial,((TPCCluster*)(subClArr->At(firstIndex)))->Position,tangentialVec);
	TVector3 circleTangent(tangentialVec.X(),
			       sprial.DY/(sprial.R*sprial.DTheta),
			       tangentialVec.Z());
	TVector3 direction = circleTangent.Unit();
	TVector3 momentum = AbsoluteMomentum*direction;

	/*
	sprial.InitMom = momentum;
	TVector3 MomentumFit(0,0,0);
	TVector3 PositionFit(0,0,0);
	Int_t PID = 0;
	if( sprial.RL > 0 ){
	  PID = 2212; /// proton
	}else{
	  PID = -211;
	}
	fitter->Fit(subClArr , sprial, PID, PositionFit, MomentumFit  );

	std::cout<<
	momentum.Print();
	MomentumFit.Print();
	hisMag->Fill(momentum.Mag(), MomentumFit.Mag()/momentum.Mag());
	*/

	if( sprial.RL > 0 ){
	  new((*pSprialArr)[nPsp]) HSprial(sprial);
	  nPsp++;
	}else{
	  new((*mSprialArr)[nMsp]) HSprial(sprial);
	  nMsp++;
	}

	HelixArr.push_back(sprial.GenerateHelix());
	ArcArr.push_back(sprial.GenerateArc());
	nHelix++;
	/*
	can->cd(1);
	total->Draw();
	(sprial.GenerateArc())->Draw();

	can->Update();
	can->Modified();
	gSystem->ProcessEvents();
	getchar();
	*/
      }
      nCand++;
    }
    std::cout<< "Sprial End" << std::endl;
    nLambda = 0;
    nMTrack = mSprialArr->GetEntries();
    nPTrack = pSprialArr->GetEntries();


    if( pSprialArr->GetEntries() >= 0 &&
	pSprialArr->GetEntries() < 5 &&
	mSprialArr->GetEntries() >= 0 &&
	mSprialArr->GetEntries() <5 ){
      Int_t pIndex = 0;

      for( int ip = 0; ip < pSprialArr->GetEntries() ; ip++){
	for( int im = 0; im < mSprialArr->GetEntries(); im++){
	  Double_t dist=0;
	  TVector3 position(0,0,0);
	  TVector3 pTangent(0,0,0);
	  TVector3 mTangent(0,0,0);
	  TVector3 pMomentum(0,0,0);
	  TVector3 mMomentum(0,0,0);
	  TLorentzVector lvp(0,0,0,0);
	  TLorentzVector lvm(0,0,0,0);
	  TLorentzVector lvl(0,0,0,0);

	  Bool_t rst = CalculateCrossing( *(HSprial*)(pSprialArr->At(ip)),
					  *(HSprial*)(mSprialArr->At(im)),
					  dist,
					  position);
	  if( !rst ){ continue; }

	  CalculateCircleTangent(*(HSprial*)(pSprialArr->At(ip)),
				 //(HSprial*)(mSprialArr->At(im)),
				 position,
				 pTangent);

	  CalculateCircleTangent(*(HSprial*)(mSprialArr->At(im)),
				 position,
				 mTangent);

	  Double_t pTm = CalculateMomentum( *(HSprial*)(pSprialArr->At(ip)), 1.);
	  Double_t mTm = CalculateMomentum( *(HSprial*)(mSprialArr->At(im)), 1.);
	  pTangent.SetMag(pTm);
	  mTangent.SetMag(mTm);
	  pMomentum = pTangent;
	  mMomentum = mTangent;
	  Double_t pE = TMath::Sqrt( pMomentum.Mag2() + pmass*pmass);
	  Double_t mE = TMath::Sqrt( mMomentum.Mag2() + pimass*pimass);
	  lvp.SetPxPyPzE( pMomentum.X(),
			  pMomentum.Y(),
			  pMomentum.Z(),
			  pE);
	  lvm.SetPxPyPzE( mMomentum.X(),
			  mMomentum.Y(),
			  mMomentum.Z(),
			  mE);
	  lvl = lvp + lvm;
	  lvl.Print();
	  new((*pMomentumArr)[pIndex]) TLorentzVector(lvp);
	  new((*mMomentumArr)[pIndex]) TLorentzVector(lvm);
	  new((*lMomentumArr)[pIndex]) TLorentzVector(lvl);
	  new((*PositionArr)[pIndex])  TVector3(position);
	  hislMass->Fill(lvl.M());
	  Distance[nLambda] = dist;

	  nLambda++;
	  std::cout<< "END" << std::endl;
	  pIndex++;
	}
      }
    }
    otr->Fill();

  }
  can->cd(4);
  hislMass->Draw();
  otf->cd();
  otr->Write();
  otf->Close();
}
