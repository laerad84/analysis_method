
#include "Fitter.h"

ClassImp( HelixFitter )
Bool_t HelixFitter::bWeight = kFALSE;
Bool_t HelixFitter::bEnergy = kFALSE;
Double_t HelixFitter::RFitCut = 50;
Double_t HelixFitter::YFitCut = 20;
TGraph2D* HelixFitter::fitTrack= new TGraph2D();
TGraph2D* HelixFitter::fitWeight = new TGraph2D();
TGraphErrors* HelixFitter::fitYTheta = new TGraphErrors();

std::vector<Double_t> HelixFitter::EArr;

HelixFitter* gHelixFitter = new HelixFitter();

HelixFitter::HelixFitter(){
  ;
}
HelixFitter::~HelixFitter(){
  ;
}
void HelixFitter::FitData( TPCTrack* track ){
  EArr.clear();
  EArr = track->EArr;
  /*
  for( int i = 0; i< track->EArr.size(); i++){
    EArr.push_back( track->EArr.at(i));
  }
  */
  HelixFitter::Fit( track->track );/// currently no consideration for cluster error
  track->helix = FitResult;
  track->Length = track->helix.TrackLength;
  FitResult.Print();
}
void HelixFitter::Fit( TGraph2D* track, TGraph2D* weights ){
  FitResult.ResetAll();
  //fitTrack = track;
  fitTrack->Set(0);
  for( int i = 0; i< track->GetN(); i++){
    fitTrack->SetPoint(i, track->GetX()[i], track->GetY()[i], track->GetZ()[i]);
    //std::cout << i << "\t" <<  fitTrack->GetX()[i] << std::endl;
  }

  if( weights != NULL ){
    HelixFitter::bWeight = kTRUE;
  }else{
    HelixFitter::bWeight = kFALSE;
  }
  CircleFit();
  HelixFit();
}
void HelixFitter::CircleFit( ){
  std::cout << __PRETTY_FUNCTION__ << " : " << __LINE__ << std::endl;

  Int_t nPoint = HelixFitter::fitTrack->GetN();
  /// Three point circle fit ///
  Int_t    pIndex[3];
  pIndex[0] = 0;
  pIndex[2] = HelixFitter::fitTrack->GetN()-1;
  pIndex[1] = (int)(pIndex[2]/2);
  if( pIndex[0] == pIndex[1] ||
      pIndex[1] == pIndex[2] ){
    std::cerr << "Cannot fit circle." << std::endl;
    return;
  }
  TVector3 vec[3];
  for( int i = 0; i< 3; i++){
    vec[i].SetXYZ( HelixFitter::fitTrack->GetX()[pIndex[i]],
		   HelixFitter::fitTrack->GetY()[pIndex[i]],
		   HelixFitter::fitTrack->GetZ()[pIndex[i]] );
  }
  FitResult = GenerateSprial( vec[0], vec[1], vec[2] );
  TVirtualFitter *fitter = TVirtualFitter::Fitter(0,3);
  TVirtualFitter::SetDefaultFitter("Minuit");
  fitter->SetFCN(HelixFitter::Chi2MinimizerCircle);
  fitter->SetParameter(0, "X",FitResult.X, 20, 0, 0);
  fitter->SetParameter(1, "Z",FitResult.Z, 20, 0, 0);
  fitter->SetParameter(2, "R",FitResult.R, 20, 0, 0);
  std::cout<< "PreFit" << std::endl;
  Double_t arglist[1] = {0};
  fitter->ExecuteCommand("MIGRAD",arglist, 0 );
  std::cout<< "Fit" << std::endl;
  FitResult.X=fitter->GetParameter(0);
  FitResult.Z=fitter->GetParameter(1);
  FitResult.R=fitter->GetParameter(2);
  FitResult.EX=fitter->GetParError(0);
  FitResult.EZ=fitter->GetParError(1);
  FitResult.ER=fitter->GetParError(2);
  std::cout << "//////////////////////////////////////////////////" << std::endl;
  FitResult.Print();
  std::cout << "//////////////////////////////////////////////////" << std::endl;
}
void HelixFitter::HelixFit( ){
  std::cout << __PRETTY_FUNCTION__ << " : " << __LINE__ << std::endl;

  Int_t nPoint = HelixFitter::fitTrack->GetN();
  fitYTheta->Set(0);
  TVector3 initVec( HelixFitter::fitTrack->GetX()[0]-FitResult.X,
		    0,
		    HelixFitter::fitTrack->GetZ()[0]-FitResult.Z);/// No Y position
  TVector3 lastVec( HelixFitter::fitTrack->GetX()[nPoint-1]-FitResult.X,
		    0,
		    HelixFitter::fitTrack->GetZ()[nPoint-1]-FitResult.Z);
  Double_t commonRL = 0;
  Double_t lastInternal = lastVec.Dot( initVec );
  Double_t lastcosVal   = lastInternal/initVec.Mag()/lastVec.Mag();
  Double_t lastTheta    = TMath::ACos( lastcosVal );

  for( Int_t ip = 1; ip < nPoint; ip++){
    TVector3 curVec(HelixFitter::fitTrack->GetX()[ip]-FitResult.X,
		    0,
		    HelixFitter::fitTrack->GetZ()[ip]-FitResult.Z);
    TVector3 External = initVec.Cross( curVec );// init X cur
    if( External.Y() != 0){
      commonRL += External.Y()/TMath::Abs(External.Y());
    }
  }

  if( nPoint -1  == 0 ){
    std::cerr <<__PRETTY_FUNCTION__ << ":"
	      << __LINE__ << " : "
	      << "Wrong enetries" << std::endl;
    return;
  }
  if( commonRL == 0 ){
    std::cerr <<__PRETTY_FUNCTION__ << ":"
	      << __LINE__ << " : "
	      << "Wrong value for commonRL" << std::endl;
    return;
  }

  commonRL = commonRL /( nPoint -1 );
  // commonRL < 0 : right handed circle
  // commonRL > 0 : left handed circle
  if( commonRL < 0 ){ FitResult.RL = -1;}
  else{ FitResult.RL = 1; }

  for( Int_t ip = 0; ip < nPoint; ip++){
    TVector3 curVec(HelixFitter::fitTrack->GetX()[ip]-FitResult.X,
		    0,
		    HelixFitter::fitTrack->GetZ()[ip]-FitResult.Z);
    Double_t Internal = initVec.Dot( curVec );// init . cur
    Double_t cosVal   = Internal/initVec.Mag()/curVec.Mag();//
    Double_t theta    = TMath::ACos( cosVal );// 0 <= theta <= pi
    TVector3 External = initVec.Cross( curVec );// init X cur
    Double_t thetaSign = 0;
    if( FitResult.RL * External.Y()  < 0 ){ /// point crosses boundary ///
      if( FitResult.RL > 0 ){
	thetaSign = 2*TMath::Pi() - theta;
      }else{
	thetaSign = -2*TMath::Pi() + theta;
      }
    }else{
      thetaSign = FitResult.RL*theta;
    }
    if( ip == nPoint -1 ){ lastTheta = thetaSign;}
    fitYTheta->SetPoint( ip, thetaSign, HelixFitter::fitTrack->GetY()[ip]);
  }

  fitYTheta->Fit("pol1","","",0.0,lastTheta*1.1);
  TF1*     rstFunc   = fitYTheta->GetFunction("pol1");
  Double_t StartingY = rstFunc->Eval(0);
  Double_t EndY      = rstFunc->Eval(lastTheta);
  Double_t dTheta    = lastTheta;
  Double_t dY        = EndY - StartingY;
  Double_t EDYDTheta = rstFunc->GetParError(1);
  Double_t EY        = rstFunc->GetParError(0);
  FitResult.DTheta = dTheta;
  FitResult.DY     = dY;
  FitResult.Y      = StartingY;
  FitResult.EY     = EY;
  FitResult.EDYDTheta= EDYDTheta;
  FitResult.InitPos= TVector3((initVec.Unit()*FitResult.R).X()+ FitResult.X,
			      StartingY,
			      (initVec.Unit()*FitResult.R).Z()+ FitResult.Z);
  FitResult.FinalPos=TVector3((lastVec.Unit()*FitResult.R).X()+ FitResult.X,
			      EndY,
			      (lastVec.Unit()*FitResult.R).Z()+ FitResult.Z);
  FitResult.TrackLength = TMath::Abs(FitResult.R * FitResult.DTheta);

  /*
  FitResult.X=fitter->GetParameter(0);
  FitResult.Z=fitter->GetParameter(1);
  FitResult.R=fitter->GetParameter(2);
  FitResult.EX=fitter->GetParError(0);
  FitResult.EZ=fitter->GetParError(1);
  FitResult.ER=fitter->GetParError(2);
  */

}
void HelixFitter::Chi2MinimizerCircle( Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t ){
  std::cout<< __PRETTY_FUNCTION__ << " : " << __LINE__ << std::endl;
  Int_t nPoint = fitTrack->GetN();
  f = 0;
  Double_t *x = HelixFitter::fitTrack->GetX();
  Double_t *y = HelixFitter::fitTrack->GetY();
  Double_t *z = HelixFitter::fitTrack->GetZ();
  Double_t *EX;
  Double_t *EZ;

  if( HelixFitter::bWeight ){ // Weight fitting
    EX = HelixFitter::fitWeight->GetX();
    EZ = HelixFitter::fitWeight->GetZ();
  }else{ // Normal Fitting
    ;
  }
  HelixFitter::bWeight = kFALSE;
  std::cout<< nPoint << std::endl;
  for( int ip = 0; ip < nPoint ; ip++){
    std::cout << ip << std::endl;
    std::cout << ip << "\t"
	      << HelixFitter::fitTrack->GetX()[ip] << "\t"
	      << HelixFitter::fitTrack->GetY()[ip] << "\t"
	      << HelixFitter::fitTrack->GetZ()[ip] << "\t"
	      << std::endl;
    std::cout << ip << "\t" << x[ip] << "\t" << z[ip] << std::endl;
    Double_t u = x[ip] - par[0];
    Double_t v = z[ip] - par[1];
    Double_t energy = EArr.at(ip);
    std::cout<< u << "\t" << v <<  std::endl;
    if( HelixFitter::bWeight ){ // Weight fitting
      Double_t dr = par[2] - TMath::Sqrt( u*u + v*v );
      Double_t er = TMath::Sqrt(EX[ip]*EX[ip]+EZ[ip]*EZ[ip]);
      if( dr < RFitCut ){
	f += dr*dr/(er*er)*energy;
      }
    }else{ // Normal Fitting
      Double_t dr = par[2] - TMath::Sqrt( u*u +v*v );
      if( dr < RFitCut ){
	f += dr*dr*energy;
      }
    }
  }
}
