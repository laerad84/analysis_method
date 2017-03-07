#include "ConverterFunction.h"
//#include "TPCData.h"
#include "Data.h"
#include "GsimData/GsimDetectorHitData.h"
#include "GsimData/GsimDetectorEventData.h"
#include "GsimData/GsimGenParticleData.h"
#include "GsimData/GsimTrackData.h"
#include "TRandom.h"
#include "TProfile.h"


std::vector<TPCHit> ConvertTPC(GsimDetectorEventData* detData, GsimGenParticleData* particle){

  std::vector<TPCHit> PadHitArr;
  TClonesArray* m_HitArr = detData->hits;
  TClonesArray* m_TrackArr = particle->briefTracks;

  /// global object for TPC hist array ///
  gTPCConvHistArray->Reset();

  std::vector<std::vector<Int_t> > hitTrackIDList;
  std::vector<Int_t> trackList;

  for( int ihit = 0; ihit < m_HitArr->GetEntries(); ihit++){
    GsimDetectorHitData* hit = (GsimDetectorHitData*)m_HitArr->At(ihit);

    //// set track ID
    Int_t trackID = 0;
    Bool_t bTrackFind = kFALSE;
    for( int itid = 0; itid < trackList.size(); itid++){
      if( trackList.at(itid) == hit->track ){
	trackID = itid;
	bTrackFind = kTRUE;
	hitTrackIDList.at(trackID).push_back( ihit );
	break;
      }
    }
    if( !bTrackFind ){
      std::vector<int> list;
      hitTrackIDList.push_back( list );
      trackList.push_back( hit->track );
      bTrackFind = kTRUE;
      trackID = trackList.size()-1;
      hitTrackIDList.at(trackID).push_back(ihit);
    }
  }

  //// generate average hit point by pad // TPCHit
  //std::vector<TPCHit> m_PadHitArr;
  Int_t iTrack = 0;
  Int_t iHit   = 0;
  Int_t indexhit = 0;

  for( int iarr = 0; iarr < hitTrackIDList.size(); iarr++){
    gTPCConvHistArray->Reset();

    // make TPCHit array by track //
    Int_t nFill = 0;
    for( int jarr = 0; jarr < hitTrackIDList.at(iarr).size(); jarr++){
      GsimDetectorHitData* hit = (GsimDetectorHitData*)m_HitArr->At(hitTrackIDList.at(iarr).at(jarr));
      //Double_t diffusionR      = 0.18*TMath::Sqrt((hit->r.Y()+300)/10);/// diffusion radius by R, drift
      Double_t diffusionR      = 2.3*TMath::Sqrt((hit->r.Y()+300)/150.);/// diffusion radius by R, drift
      Double_t diffusionG      = 0.1; // diffusion radius by GEM, ignored in the current version
      Int_t    nElectrons      = (int)(hit->edep*TPCGlobals::sTPC_NEperEnergy*4);//edep:MeV NEperEnergy : 2/0.00015 ->4 electrons per 150eV. Change to real value;
      if( TMath::Abs(hit->r.Y()) >= 300 ){ continue; }
      for( int ielec = 0; ielec < nElectrons; ielec++ ){
	double xhit = gRandom->Gaus(hit->r.X(),diffusionR);
	double zhit = gRandom->Gaus(hit->r.Z(),diffusionR);
	double yhit = gRandom->Gaus(hit->r.Y(),0.5);
	if( TMath::Abs(yhit) >= 300  ){ continue; }
	//std::cout<< xhit << "\t" << yhit << "\t" << zhit << std::endl;
	Int_t rowID  = -1;
	Int_t colID  = -1;
	Bool_t fTPCHit = gTPCIDHandler->GetRowCol(xhit,zhit,rowID,colID);
	if( !fTPCHit ){ continue; }
	if( rowID < 0 || rowID >= 32){ continue; }
	gTPCConvHistArray->hisRowY[rowID]->Fill( colID, yhit );
      }

      int rid, cid;
      gTPCIDHandler->GetRowCol(hit->r.X(),hit->r.Z(), rid, cid );
      //std::cout<< "RCID " << rid << "\t" << cid  << std::endl;
      nFill++;
    }

    /// Gather PadHit information from hists ///
    for( int irow = 0; irow < 32; irow++){
      TProfile* prof = gTPCConvHistArray->hisRowY[irow]->ProfileX(); /// Profile for calculating y position
      TH1D*     prjt = gTPCConvHistArray->hisRowY[irow]->ProjectionX();
      for( int icol = 0; icol < TPCGlobals::sTPC_Pad_Parameter[irow][1]; icol++ ){
	double ymean = prof->GetBinContent(icol+1);
	Int_t  padID  = gTPCIDHandler->GetID( irow, icol );
	double   yerr = prof->GetBinError( icol+1 );
	double   ene  = prjt->GetBinContent( icol+1 );
	if( ene < 10 ){ continue; }//75eV
	TPCHit padHit;
	TVector2 TPCXZ = gTPCIDHandler->GetXZ( irow, icol );
	TVector3 vec(TPCXZ.X(), ymean, TPCXZ.Y());
	padHit.SetPadID(padID);
	padHit.SetHitID( iHit );
	padHit.SetTrackID( iTrack);
	padHit.SetMCTrackID( trackList.at(iarr));
	padHit.SetPosition(vec);
	padHit.SetRowCol( irow, icol );
	padHit.SetEnergy( ene );
	padHit.m_signal=0;
	padHit.SetRad( gTPCIDHandler->GetRad(irow,icol));
	padHit.SetPhi( gTPCIDHandler->GetPhi(irow,icol));
	PadHitArr.push_back(padHit);
	iHit++;
	/*
	std::cout<< padID << "\t"
		 << padHit.PadID() << "\t"
		 << PadHitArr.at(PadHitArr.size()-1).PadID() << "\t"
		 << PadHitArr.size() << std::endl;
	*/
      }
      prof->Delete();
      prjt->Delete();
    }
    iTrack++;
  }

  return PadHitArr;
}
Bool_t HitClustering( TClonesArray* HitArr, TClonesArray* clusterArr ){
  //clusterArr.clear();
  clusterArr->Clear();

  std::vector<TPCHit> hitArr;
  for( int i = 0; i< HitArr->GetEntries(); i++){
    TPCHit* hit = (TPCHit*)HitArr->At(i);
    hitArr.push_back(*hit);
  }

  TPCCluster clArr[200];

  Int_t ClusterID = 0;
  Int_t clIndex = 0;
  for( int irow = 0; irow < 32; irow++){
    //std::cout<< "irow " << irow << std::endl;
    Int_t nRow = 0;
    for( int i = 0; i< hitArr.size(); i++){
      if( hitArr.at(i).Row() == irow){ nRow++; }
    }
    while(nRow > 0 ){
      std::vector<TPCHit> sub_hitArr;
      Int_t clusterHitSize = 0;
      ///Check first hit ///
      for( int ihit = 0; ihit < hitArr.size(); ihit++){
	if( hitArr.at(ihit).Row() == irow ){
	  if( sub_hitArr.size() == 0){
	    sub_hitArr.push_back( hitArr.at(ihit));
	    hitArr.erase(hitArr.begin()+ihit);
	    break;
	  }
	}
      }
      if( sub_hitArr.size() == 0 ){ break; }

      /// Add adjacent hits to sub_hitArr
      for( int ihit = 0; ihit < sub_hitArr.size(); ihit++){
	for( int jhit = 0; jhit < hitArr.size();){
	  if( TMath::Abs(hitArr.at(jhit).Col() - sub_hitArr.at(ihit).Col()) < 2 &&
	      TMath::Abs(hitArr.at(jhit).Position().Y() - sub_hitArr.at(ihit).Position().Y()) < 10 &&
	      hitArr.at(jhit).Row() == irow ){
	    sub_hitArr.push_back(hitArr.at(jhit));
	    hitArr.erase( hitArr.begin() + jhit);
	  }else{
	    jhit++;
	  }
	}
      }

      if( sub_hitArr.size() == 0 ){ break; }
      Int_t RowID(-1);
      Double_t ColID(0);
      Double_t ColRMS(0);
      Double_t X(0),Y(0),Z(0);
      Double_t YMin(9999);
      Double_t YMax(-9999);
      Double_t ColMin(9999);
      Double_t ColMax(-9999);
      Double_t PhiMin(9999);
      Double_t PhiMax(-9999);
      Double_t TotalEnergy(0);

      RowID = irow;
      for( int ihit = 0; ihit < sub_hitArr.size(); ihit++){
	TotalEnergy += sub_hitArr.at(ihit).Energy();
	ColID       += sub_hitArr.at(ihit).Energy()*sub_hitArr.at(ihit).Col();
	ColRMS      += sub_hitArr.at(ihit).Energy()*TMath::Power(sub_hitArr.at(ihit).Col(),2);
	Y           += sub_hitArr.at(ihit).Energy()*sub_hitArr.at(ihit).Position().Y();

	if( YMin > sub_hitArr.at(ihit).Position().Y() ){ YMin = sub_hitArr.at(ihit).Position().Y();}
	if( YMax < sub_hitArr.at(ihit).Position().Y() ){ YMax = sub_hitArr.at(ihit).Position().Y();}
	if( ColMin > sub_hitArr.at(ihit).Position().Y()){ ColMin = sub_hitArr.at(ihit).Col();}
	if( ColMax < sub_hitArr.at(ihit).Position().Y()){ ColMax = sub_hitArr.at(ihit).Col();}
	if( PhiMin > sub_hitArr.at(ihit).Phi() ){ PhiMin = sub_hitArr.at(ihit).Phi();}
	if( PhiMax < sub_hitArr.at(ihit).Phi() ){ PhiMax = sub_hitArr.at(ihit).Phi();}
      }
      ColID = ColID / TotalEnergy;
      ColRMS= ColRMS/ TotalEnergy - ColID*ColID;

      Double_t Rad = TPCGlobals::sTPC_Pad_Parameter[RowID][2];
      Double_t ColThetaDeg = (-180 + TPCGlobals::sTPC_Pad_Parameter[RowID][4] + (ColID+0.5)*TPCGlobals::sTPC_Pad_Parameter[RowID][3]);
      X = Rad * TMath::Sin( ColThetaDeg * TMath::DegToRad() );
      Y = Y / TotalEnergy;
      Z = TPCGlobals::sTPC_zOffset + Rad * TMath::Cos( ColThetaDeg * TMath::DegToRad() );

      /// Set Cluster information using sub_hitArr
      clArr[ClusterID].ID       = ClusterID;
      clArr[ClusterID].Energy   = TotalEnergy;
      clArr[ClusterID].Row      = irow;
      clArr[ClusterID].YMin     = YMin;
      clArr[ClusterID].YMax     = YMax;
      clArr[ClusterID].Col      = ColID;
      clArr[ClusterID].ColMin   = ColMin;
      clArr[ClusterID].ColMax   = ColMax;
      clArr[ClusterID].PhiMin   = PhiMin;
      clArr[ClusterID].PhiMax   = PhiMax;
      clArr[ClusterID].Position = TVector3(X,Y,Z);
      clArr[ClusterID].NHit     = sub_hitArr.size();
      clArr[ClusterID].Mother   = -1;
      clArr[ClusterID].Daughter = -1;
      clArr[ClusterID].MinMotherDist = -1;
      clArr[ClusterID].MinDaughterDist = -1;
      clArr[ClusterID].nDaughter= 0;
      clArr[ClusterID].nMother  = 0;
      clArr[ClusterID].MotherID.clear();
      clArr[ClusterID].DaughterID.clear();
      clArr[ClusterID].Blocked = false;
      //new((*clusterArr)[ClusterID]) TPCCluster(cluster);

      ClusterID++;
      nRow = 0;
      for( int i = 0; i< hitArr.size(); i++){
	if( hitArr.at(i).Row() == irow){ nRow++; }
      }
      if( ClusterID >= 200 ){ return false; }
    }
  }

  ////// connection
  Int_t ClusterEntries = ClusterID;
  for( Int_t icl = 0; icl < ClusterEntries; icl++){
    Int_t cRow = clArr[icl].Row;
    Int_t Mother = clArr[icl].ID;

    //// Find Daughter at next row
    Bool_t bFindNext = false;
    for( Int_t iloop = 0; iloop < 3; iloop++){
      for( Int_t jcl =0; jcl < ClusterEntries; jcl++){
	if( icl == jcl ){ continue; }
	Int_t dRow  = clArr[jcl].Row;
	if( cRow != dRow - (iloop+1) ){ continue; }
	if( !AdjCluster(clArr[icl],clArr[jcl],iloop+1) ){
	  /// not connected
	  continue;
	}else{
	  /// connected
	  if( clArr[icl].MotherY < -500 ){/// mother is not set
	    if( iloop > 0 ){ continue; }
	    TVector3 dp(clArr[icl].Position.X() - clArr[jcl].Position.X(),
			clArr[icl].Position.Y() - clArr[jcl].Position.Y(),
			clArr[icl].Position.Z() - clArr[jcl].Position.Z());
	    Double_t dist = dp.Mag();
	    clArr[icl].DaughterID.push_back( clArr[jcl].ID);
	    clArr[jcl].MotherID.push_back(clArr[icl].ID);
	    clArr[icl].nDaughter = clArr[icl].DaughterID.size();
	    clArr[jcl].nMother   = clArr[jcl].MotherID.size();
	    if(clArr[icl].MinDaughterDist <= 0 ){
	      clArr[icl].Daughter = clArr[jcl].ID;
	      clArr[icl].MinDaughterDist = dist;
	      clArr[icl].DaughterY = clArr[jcl].Position.Y();
	    }else{
	      if( clArr[icl].MinDaughterDist < dist ){
		clArr[icl].Daughter = clArr[jcl].ID;
		clArr[icl].MinDaughterDist = dist;
		clArr[icl].DaughterY = clArr[jcl].Position.Y();
	      }
	    }

	    if( clArr[jcl].MinMotherDist <= 0 ){
	      clArr[jcl].Mother = clArr[icl].ID;
	      clArr[jcl].MinMotherDist = dist;
	      clArr[jcl].MotherY   = clArr[icl].Position.Y();
	    }else{
	      if( clArr[jcl].MinMotherDist < dist ){
		clArr[jcl].Mother = clArr[icl].ID;
		clArr[jcl].MinMotherDist = dist;
		clArr[jcl].MotherY   = clArr[icl].Position.Y();
	      }
	    }
	    bFindNext = true;
	  }else{//// mother is exist
	    /// Mother exist ///
	    //Constraints on Y

	    TVector3 dp(clArr[icl].Position.X() - clArr[jcl].Position.X(),
			clArr[icl].Position.Y() - clArr[jcl].Position.Y(),
			clArr[icl].Position.Z() - clArr[jcl].Position.Z());
	    ///// Should be edited for dist to  xz dist
	    Double_t dist       = dp.Mag();
	    Double_t DYDLMother = (clArr[icl].Position.Y()-clArr[icl].MotherY)/clArr[icl].MinMotherDist;
	    Double_t YExpect    = DYDLMother * dist + clArr[icl].Position.Y();
	    Double_t DeltaY     = clArr[jcl].Position.Y() - YExpect;

	    if( TMath::Abs(DeltaY) > 20 ){  continue; }

	    clArr[icl].DaughterID.push_back( clArr[jcl].ID);
	    clArr[jcl].MotherID.push_back(clArr[icl].ID);
	    clArr[icl].nDaughter = clArr[icl].DaughterID.size();
	    clArr[jcl].nMother   = clArr[jcl].MotherID.size();
	    if(clArr[icl].MinDaughterDist <= 0 ){
	      clArr[icl].Daughter = clArr[jcl].ID;
	      clArr[icl].MinDaughterDist = dist;
	      clArr[icl].DaughterY = clArr[jcl].Position.Y();
	    }else{
	      if( clArr[icl].MinDaughterDist < dist ){
		clArr[icl].Daughter = clArr[jcl].ID;
		clArr[icl].MinDaughterDist = dist;
		clArr[icl].DaughterY = clArr[jcl].Position.Y();
	      }
	    }

	    if( clArr[jcl].MinMotherDist <= 0 ){
	      clArr[jcl].Mother = clArr[icl].ID;
	      clArr[jcl].MinMotherDist = dist;
	      clArr[jcl].MotherY   = clArr[icl].Position.Y();
	    }else{
	      if( clArr[jcl].MinMotherDist < dist ){
		clArr[jcl].Mother = clArr[icl].ID;
		clArr[jcl].MinMotherDist = dist;
		clArr[jcl].MotherY   = clArr[icl].Position.Y();
	      }
	    }
	    bFindNext = true;
	  }
	}
      }
      if( bFindNext ){ break; }
    }
  }
  for( int i = 0; i< ClusterEntries; i++){
    new((*clusterArr)[i]) TPCCluster(clArr[i]);
  }
  return true;
}
TPCTrack ConvertToTrack( TClonesArray* clusterArr ){
  TPCTrack track;
  track.Init();
  Int_t Entries = clusterArr->GetEntries();
  Int_t initialRow = 32;
  Int_t finalRow = -1;
  for( int icl = 0; icl < Entries; icl++){
    TPCCluster* cl = (TPCCluster*)clusterArr->At(icl);
    if( cl->Row > finalRow ){ finalRow = cl->Row;}
    if( cl->Row < initialRow ){ initialRow = cl->Row;}
    track.track->SetPoint( icl, cl->Position.X(), cl->Position.Y(), cl->Position.Z());
    track.trackErr->SetPoint( icl, cl->PositionErr.X(), cl->PositionErr.Y(), cl->PositionErr.Z());
    track.EArr.push_back( cl->Energy );
    track.DepE += cl->Energy;
  }
  track.InitialRowID = initialRow;
  track.FinalRowID   = finalRow;

  return track;
}
Bool_t AdjCluster( TPCCluster c_0, TPCCluster c_1, Int_t DRow ){
  /*
  Double_t c0Min = c_0.PhiMin;
  Double_t c0Max = c_0.PhiMax;
  Double_t c1Min = c_1.PhiMin;
  Double_t c1Max = c_1.PhiMax;
  */

  Double_t c0Min = c_0.PhiMin - TPCGlobals::sTPC_Pad_Parameter[c_0.Row][3]/2.*TMath::DegToRad();
  Double_t c0Max = c_0.PhiMax + TPCGlobals::sTPC_Pad_Parameter[c_0.Row][3]/2.*TMath::DegToRad();
  Double_t c1Min = c_1.PhiMin - TPCGlobals::sTPC_Pad_Parameter[c_1.Row][3]/2.*TMath::DegToRad();
  Double_t c1Max = c_1.PhiMax + TPCGlobals::sTPC_Pad_Parameter[c_1.Row][3]/2.*TMath::DegToRad();

  if(( c_1.Row - c_0.Row == DRow) &&
     ( c0Min <= c1Max && c0Max >= c1Min ) &&
     ( TMath::Abs( c_1.Position.Y() - c_0.Position.Y()) < 40 )){
    return kTRUE;
  }else{
    /*
    if( c_1.Row - c_0.Row == DRow ){
      std::cout<<"Debug: " << DRow << "\t"
	       << c_0.Row << "\t" << c_0.PhiMin << "\t" <<  c_0.PhiMax << "\t"
	       << c_1.Row << "\t" << c_1.PhiMin << "\t" <<  c_1.PhiMax << std::endl;
    }
    */
    return kFALSE;
  }
}
Bool_t AdjCluster( TPCCluster* c_0, TPCCluster* c_1, Int_t DRow){

  Double_t c0Min = c_0->PhiMin;// - TPCGlobals::sTPC_Pad_Parameter[c_0->Row][3]/2.*TMath::DegToRad();
  Double_t c0Max = c_0->PhiMax;// + TPCGlobals::sTPC_Pad_Parameter[c_0->Row][3]/2.*TMath::DegToRad();
  Double_t c1Min = c_1->PhiMin;// - TPCGlobals::sTPC_Pad_Parameter[c_1->Row][3]/2.*TMath::DegToRad();
  Double_t c1Max = c_1->PhiMax;// + TPCGlobals::sTPC_Pad_Parameter[c_1->Row][3]/2.*TMath::DegToRad();

  if(( c_1->Row - c_0->Row == DRow) &&
     ( c0Min <= c1Max && c0Max >= c1Min ) ){
    //&&
    //( TMath::Abs( c_1->Position.Y() - c_0->Position.Y()) < 40 )){
    return kTRUE;

  }else{
    /*
    if( c_1->Row - c_0->Row == DRow ){
      std::cout<<"Debug: " << DRow << "\t"
	       << c_0->Row << "\t" << c_0->PhiMin << "\t" <<  c_0->PhiMax << "\t"
	       << c_1->Row << "\t" << c_1->PhiMin << "\t" <<  c_1->PhiMax << std::endl;
    }
    */
    return kFALSE;
  }

}
std::vector<Int_t> GetListOfTrackRoot(TClonesArray* clArr){// Searching track's root points.
  std::vector<Int_t> BlockIndexList;
  for( Int_t i = 0; i< clArr->GetEntries(); i++){
    TPCCluster* cl = (TPCCluster*)clArr->At(i);
    if( cl->nDaughter != 1){ continue; }
    if( cl->nMother == 0 ){
	BlockIndexList.push_back(i);
    }else if( cl->nMother == 1 ){
      Int_t motherID = cl->Mother;
      TPCCluster* mcl = (TPCCluster*)clArr->At(motherID);
      if( mcl->nDaughter > 1 ){
	BlockIndexList.push_back( i );
      }
    }else if( cl->nMother > 1 ){
      BlockIndexList.push_back(i);
    }
  }
  return BlockIndexList;
}
Bool_t ClusterBlocker( TClonesArray* clArr, TClonesArray* blockArr ,std::vector<Int_t> BlockRoot, Int_t BlockIndex){/// TClonesArray("TPCCluster");
  if( BlockIndex >= clArr->GetEntries() ){ return false;}
  if( blockArr->GetEntries() != 0 ){ blockArr->Clear();}

  Bool_t bFind = true;
  Int_t clIndex = BlockIndex;
  Int_t iCL = 0;
  while( bFind ){
    TPCCluster* cl = (TPCCluster*)clArr->At(clIndex);
    if( clIndex != cl->ID ){ std::cerr << "Hmm..... "<< __PRETTY_FUNCTION__ << std::endl;
      bFind = false;
      return false;
    }
    new((*blockArr)[iCL]) TPCCluster(*cl);
    iCL++;
    if( cl->nDaughter == 1 ){
      bool existInList = false;
      for( int isize = 0; isize < BlockRoot.size(); isize++){
	if( cl->Daughter == BlockRoot.at(isize)){
	  existInList = true;
	}
      }
      if( existInList ){ bFind = false; break; }
      clIndex = cl->Daughter;
    }else{
      bFind = false;
    }
  }
  return true;
}

ClassImp( TrackHandler )
TrackHandler* gTrackHandler = new TrackHandler();

TrackHandler::TrackHandler(){;}
TrackHandler::~TrackHandler(){;}
Bool_t TrackHandler::ConversionToTrack( TClonesArray* clusterArr, TPCTrack& track){
  track.Init();
  Int_t Entries = clusterArr->GetEntries();
  Int_t initialRow = 32;
  Int_t finalRow   = -1;
  for( int icl = 0; icl < Entries; icl++){
    TPCCluster* cl = (TPCCluster*)clusterArr->At(icl);
    if( cl->Row > finalRow ){ finalRow = cl->Row;}
    if( cl->Row < initialRow ){ initialRow = cl->Row;}
    track.track->SetPoint( icl, cl->Position.X(), cl->Position.Y(), cl->Position.Z());
    track.trackErr->SetPoint( icl, cl->PositionErr.X(), cl->PositionErr.Y(), cl->PositionErr.Z());
    track.EArr.push_back( cl->Energy );
    track.DepE += cl->Energy;
  }
  track.InitialRowID = initialRow;
  track.FinalRowID   = finalRow;
  /// Fitting Track with helix ///
  gHelixFitter->FitData( &track );
  return true;
}
void TrackHandler::EvalTrackInit(TPCTrack& track, TVector3 initPoint, Double_t bField ){//Setting initial track information
  /// Calculate momentum with initial point
}

//// returns TPCTrack. if condition were not admitted, return zerosize track.
TPCTrack* TrackHandler::MergingTrack( TPCTrack* trk0, TPCTrack* trk1 ){
  TPCTrack* track = new TPCTrack();

  /// Check track parameters
  if( trk0->helix.RL != trk1->helix.RL ){ return track; }
  if( trk0->helix.R > 6000 || trk0->helix.R < 100 ){ return track; }
  if( trk1->helix.R > 6000 || trk1->helix.R < 100 ){ return track; }

  Int_t direction = 0;/// 1 : trk0 -- trk1,  -1: trk1 -- trk0
  if( trk0->InitialRowID > trk1->FinalRowID ){
    direction  = -1;
  }
  if( trk0->FinalRowID < trk1->InitialRowID ){
    direction = 1;
  }
  if( direction == 0 ){ return track; }

  Double_t DistY0=0;
  Double_t DistY1=0;
  Double_t DistR0=0;
  Double_t DistR1=0;
  trk0->helix.CalculateDist( trk1->track, DistR0, DistY0 );
  trk1->helix.CalculateDist( trk0->track, DistR1, DistY1 );

  //// Check Connection
  if(( DistY0 < 10 && DistR0 < 10 ) ||
     ( DistY1 < 10 && DistR0 < 10 )){
    Int_t nPoint0 = trk0->track->GetN();
    Int_t nPoint1 = trk1->track->GetN();
    if( direction == 1){/// common order
      for( int ip = 0; ip < nPoint0; ip++){
	track->track->SetPoint( ip,
				trk0->track->GetX()[ip],
				trk0->track->GetY()[ip],
				trk0->track->GetZ()[ip]);
	track->EArr.push_back( trk0->EArr.at(ip));
	track->DepE += trk0->EArr.at(ip);

      }
      for( int ip = 0; ip < nPoint1; ip++){
	track->track->SetPoint( ip + nPoint0,
				trk1->track->GetX()[ip],
				trk1->track->GetY()[ip],
				trk1->track->GetZ()[ip]);
	track->EArr.push_back( trk1->EArr.at(ip));
	track->DepE += trk1->EArr.at(ip);
      }
      track->InitialRowID = trk0->InitialRowID;
      track->FinalRowID   = trk1->FinalRowID;

    }else if( direction == -1 ){/// reverse order
      for( int ip = 0; ip < nPoint1; ip++){
	track->track->SetPoint( ip,
				trk1->track->GetX()[ip],
				trk1->track->GetY()[ip],
				trk1->track->GetZ()[ip]);
	track->EArr.push_back( trk1->EArr.at(ip));
	track->DepE += trk1->EArr.at(ip);
      }
      for( int ip = 0; ip < nPoint0; ip++){
	track->track->SetPoint( ip + nPoint1,
				trk0->track->GetX()[ip],
				trk0->track->GetY()[ip],
				trk0->track->GetZ()[ip]);
	track->EArr.push_back( trk0->EArr.at(ip));
	track->DepE += trk0->EArr.at(ip);
      }
      track->InitialRowID = trk1->InitialRowID;
      track->FinalRowID   = trk0->FinalRowID;
    }
  }else{ return track; }
  return track;
}
