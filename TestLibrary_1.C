#include "TVirtualFitter.h"

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


void TestLibrary_1(){
  TString mSystemName = gSystem->GetFromPipe("uname");
  std::cout<< mSystemName << std::endl;
  if( mSystemName == "Darwin"){
    gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method.dylib");
  }else if( mSystemName == "Linux"){
    gSystem->Load("~/work/E42_ana/analysis_method/lib/libanalysis_method.so");
  }
  GenFitter* fitter = new GenFitter("Data/TPC.root");
  //TVirtualFitter::SetDefaultFitter("Minuit");
  grTrackSample = new TGraph2D();
  fitter->SetField("~/local/hep/E42/E42/data/field/SC_FieldMap.root",TVector3(0,0,0),TVector3(0,0,0));
  //TFile* itf = new TFile("HDE42_All.root");
  //TTree* itr = (TTree*)itf->Get("convTree");
  //TFile* itf = new TFile("TestHit.root");
  //TTree* itr = (TTree*)itf->Get("padHit");
  //TFile* itf = new TFile("HD_0_conv.root");
  TFile* itf = new TFile("ProtonTest_conv.root");
  TTree* itr = (TTree*)itf->Get("convTree");

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
  TClonesArray* ClArr = new TClonesArray("TPCCluster");
  //TTree* otr = new TTree("testTre","test");
  //otr->Branch("ClArr",&ClArr,256000,-1);

  const Int_t maxNBlock = 50;
  const Int_t maxNTrackCand    = 50;
  TGraph* grCandidates[maxNBlock];
  for( int i = 0; i< maxNBlock; i++){
    grCandidates[i] = new TGraph();
    grCandidates[i]->SetMarkerStyle( 20+i%8);
    grCandidates[i]->SetMarkerColor( i/7 +1);
  }


  for( int ievt = 0; ievt <100 ;ievt++){
    ClArr->Clear();
    itr->GetEntry( ievt );
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
    if( !bCluster ){ continue; }
    std::cout<< "ClArr : " << ClArr->GetEntries() << std::endl;
    for( int i = 0; i < ClArr->GetEntries(); i++){
      TPCCluster* cl = (TPCCluster*)(ClArr->At(i));
      hisClusterSizeEne->Fill(cl->NHit, cl->Energy);
      //cl->Print();
      std::cout<< i <<"\t" <<  cl->ID << "\t"
	       << cl->MotherID.size() << "\t" << cl->DaughterID.size() << "\t"
	       << cl->Mother          << "\t" << cl->Daughter << "\t"
	       << std::endl;
      total->SetNextPoint( cl->Position.X(), cl->Position.Y(), cl->Position.Z());
    }


    std::cout<< "GetList " << std::endl;
    std::vector<Int_t> RootIDList = GetListOfTrackRoot( ClArr );
    std::cout<<"Root ID list size : " << RootIDList.size() << std::endl;
    if( RootIDList.size() > maxNBlock ){ std::cout << "too many blocks " << std::endl; }
    Int_t nCand = 0;
    Int_t nHelix = 0;
    std::vector<HSprial>      SprialArr;
    std::vector<TPolyLine3D*> HelixArr;
    std::vector<TPolyLine*>   ArcArr;

    for( Int_t arrIndex = 0; arrIndex < RootIDList.size(); arrIndex++){
      TClonesArray* subClArr = new TClonesArray("TPCCluster");
      Bool_t bResult = ClusterBlocker( ClArr, subClArr,RootIDList, RootIDList.at(arrIndex));
      if( !bResult){ continue; }
      for( Int_t index = 0; index < subClArr->GetEntries(); index++){
	TPCCluster* subCl = (TPCCluster*)subClArr->At(index);
	std::cout<< subCl->Mother       << "\t" << subCl->ID << "\t" << subCl->Daughter << "\t"
		 << subCl->nMother      << "\t" << subCl->nDaughter << "\t"
		 << subCl->Row          << "\t" << subCl->Col << "\t"
		 << subCl->Position.X() << "\t" << subCl->Position.Z() << "\t"
		 << std::endl;
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
	circleFitter->SetParameter(0,"X", sprial.X, 10, 0,0 );
	circleFitter->SetParameter(1,"Z", sprial.Z, 10, 0, 0);
	circleFitter->SetParameter(2,"R", sprial.R, 10, 0, 0);
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
	/*
	///// Fit Dy/Dtheta
	TGraph* graphYL = new TGraph();
	Double_t StartAngle = TMath::ATan2( sprial.InitPos.X() - sprial.X,
					    sprial.InitPos.Z() - sprial.Z);
	for( Int_t ip = 0; ip < subClArr->GetEntries(); ip++){
	  TPCCluster* clref = (TPCCluster*)subClArr->At(ip);
	  /// calculate length
	  Double_t CurrentAngle = TMath::ATan2( clref->Position.X() - sprial.X,
					       clref->Position.Z() - sprial.Z);
	  Double_t deltaAngle = CurrentAngle - StartAngle;
	  if( deltaAngle*sprial.RL < 0 ){ /// crossed pi ///
	    if( deltaAngle > 0 ){
	      deltaAngle = 2*TMath::Pi()-deltaAngle;
	    }else{
	      deltaAngle = 2*TMath::Pi()+deltaAngle;
	    }
	  }
	  double l = deltaAngle*sprial.R;
	  graphYL->SetPoint( ip,l,clref->Position.Y());
	}
	graphYL->Fit("pol1");
	TF1* func = graphYL->GetFunction("pol1");
	std::cout<< func->GetParameter(0) << "\t"
		 << func->GetParameter(1) << std::endl;
	double YOffset = func->GetParameter(0);
	double DTheta = TMath::Pi()/2.;
	double DY     = func->GetParameter(1)/sprial.R*TMath::Pi()/2.;
	sprial.DTheta = DTheta;
	sprial.Y      = func->GetParameter(0);
	sprial.DY     = DY;
	*/

	/*
	TVirtualFitter *yFitter      = TVirtualFitter::Fitter(0,2);
	yFitter->SetFCN(HelixY);
	Double_t InitTheta= TMath::ATan2( sprial.InitPos.X() - sprial.X, sprial.InitPos.Z() - sprial.Z );
	yFitter->SetParameter(0,"Y", sprial.InitPos.Y() - InitTheta*sprial.DY/sprial.DTheta, 10, 0,0);
	yFitter->SetParameter(1,"DYDT", sprial.DY/sprial.DTheta, 20, 0, 0);

	Double_t arglist1[1] = {0};
	yFitter->ExecuteCommand("MIGRAD", arglist1, 0 );

	/// Make Sprial ///
	HSprial fitSprial;
	fitSprial.X      = circleFitter->GetParameter(0);
	fitSprial.Y      = yFitter->GetParameter(0);
	fitSprial.Z      = circleFitter->GetParameter(1);
	fitSprial.R      = circleFitter->GetParameter(2);
	fitSprial.DTheta = TMath::Pi()/2;
	fitSprial.DY     = yFitter->GetParameter(2)*TMath::Pi()/2;
	fitSprial.RL     = sprial.RL;
	fitSprial.InitPos= TVector3( 0,fitSprial.Y, fitSprial.R);

	std::cout<< "**************************************" << std::endl;
	std::cout<< nHelix << "\t"
		 << circleFitter->GetParameter(0) << "\t"
		 << circleFitter->GetParameter(1) << "\t"
		 << circleFitter->GetParameter(2) << "\t"
		 << yFitter->GetParameter(0) << "\t"
		 << yFitter->GetParameter(1) << "\t"
		 << std::endl;
	std::cout<< "**************************************" << std::endl;
	std::cout<< "test" << std::endl;
	*/


	///Genfit tester ///

	Double_t AbsoluteMomentum = CalculateMomentum( sprial, 1.);
	TVector3 tangentialVec(0,0,0);
	CalculateCircleTangent(sprial,((TPCCluster*)(subClArr->At(firstIndex)))->Position,tangentialVec);
	TVector3 circleTangent(tangentialVec.X(),
			       sprial.DY/(sprial.R*sprial.DTheta),
			       tangentialVec.Z());
	TVector3 direction = circleTangent.Unit();
	TVector3 momentum = AbsoluteMomentum*direction;
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

	SprialArr.push_back(sprial);
	HelixArr.push_back(sprial.GenerateHelix());
	ArcArr.push_back(sprial.GenerateArc());
	nHelix++;
      }
      nCand++;
    }




    for( int iHelix = 0; iHelix < nHelix; iHelix++){
      ArcArr.at(iHelix).SetLineColor( iHelix%7 +1);
      HelixArr.at(iHelix).SetLineColor( iHelix%7 +1);
      ArcArr.at(iHelix).SetLineWidth( 3 );
      HelixArr.at(iHelix).SetLineWidth( 3 );
    }
    std::cout<< "Draw" << std::endl;
    can->cd(1);
    for( int iCand  = 0; iCand < nCand; iCand++){
      grCandidates[iCand]->Draw("P");
    }
    for( int iHelix = 0; iHelix < nHelix; iHelix++){
      ArcArr.at(iHelix).Draw();
    }
    can->cd(2);
    hisBox->Draw();
    for( int iCand = 0; iCand < nCand; iCand++){
      track[iCand]->Draw("P");
    }
    for( int iHelix = 0; iHelix < nHelix; iHelix++){
      HelixArr.at(iHelix).Draw();
    }



    can->cd(3);
    hisBox->Draw();
    total->Draw("P");
    can->Update();
    can->Modified();
    gSystem->ProcessEvents();
    getchar();
    for( int iCand = 0; iCand < maxNBlock; iCand++){
      track[iCand]->Delete();
    }
    for( int iHelix = 0; iHelix < nHelix; iHelix++){
      TPolyLine3D* line = HelixArr.at(iHelix);
      TPolyLine*   arc  = ArcArr.at(iHelix);
      delete arc;
      delete line;
    }

    //can->cd(2)->Clear();
  }
  //hisClusterSizeEne->Write();
  can->cd(4);
  hisMag->Draw();
}
