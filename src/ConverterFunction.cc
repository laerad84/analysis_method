#include "ConverterFunction.h"
#include "TPCData.h"
#include "GsimData/GsimDetectorHitData.h"
#include "GsimData/GsimDetectorEventData.h"
#include "GsimData/GsimGenParticleData.h"
#include "GsimData/GsimTrackData.h"
#include "TRandom.h"
#include "TProfile.h"

std::vector<TPCPadHit> ConvertTPC(GsimDetectorEventData* detData, GsimGenParticleData* particle){

  std::vector<TPCPadHit> PadHitArr;
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

  //// generate average hit point by pad // TPCPadHit
  //std::vector<TPCPadHit> m_PadHitArr;
  Int_t iTrack = 0;
  Int_t iHit   = 0;
  Int_t indexhit = 0;

  for( int iarr = 0; iarr < hitTrackIDList.size(); iarr++){
    gTPCConvHistArray->Reset();

    // make TPCPadHit array by track //
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
	TPCPadHit padHit;
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

/*
std::vector<TPCPadHitCluster> HitClustering(std::vector<TPCPadHit> &hitArr){
  std::vector<TPCPadHitCluster> clusterArr;
  Int_t clusterID = 0;
  for( int irow = 0; irow < 32; irow++){
    while( hitArr.size() != 0 ){
      TPCPadHitCluster cluster  = TPCClusterer( irow, hitArr );
      if( cluster.HitArr().size() != 0 ){
	cluster.SetClusterID(clusterID);
	cluster.m_nHit = cluster.HitArr().size();
	cluster.Evaluate();
	std::cout<<"HitClustering"
		 << cluster.HitArr().size() << "\t"
		 << cluster.nHit()          << "\t"
		 << cluster.Energy()        << "\t"
		 << cluster.YMin()          << "\t"
		 << cluster.YMax()          << std::endl;
	clusterArr.push_back( cluster );
	clusterID++;
      }else{ break; }
    }
  }
  return clusterArr;
}*/
Bool_t HitClusteringB( std::vector<TPCPadHit> &hitArr, std::vector<TPCPadHitCluster> &clusterArr ){
  clusterArr.clear();
  Int_t clusterID = 0;
  Int_t clIndex   = 0;
  //for( int irow = 0; irow < 32; irow++){
  for( int irow = 0; irow < 12; irow++){
    while( hitArr.size() != 0 ){
      TPCPadHitCluster cluster  = TPCClusterer( irow, hitArr );
      if( cluster.HitArr().size() != 0 ){

	clIndex = clusterArr.size();
	cluster.SetClusterID(clusterID);
	cluster.m_nHit = cluster.HitArr().size();
	cluster.Evaluate();
	clusterArr.push_back( cluster );
	clusterArr.at(clIndex).m_HitArr.clear();
	for( int ihit = 0; ihit < cluster.HitArr().size(); ihit++){
	  clusterArr.at(clIndex).AddPadHit( cluster.HitArr().at(ihit));
	  //std::cout<< "cl : "<< clusterArr.at(clIndex).HitArr().size() << std::endl;
	}
	clusterArr.at(clIndex).Evaluate();
	/*
	std::cout<< clIndex  << "\t"
		 << cluster.HitArr().size() << "\t"
		 << cluster.nHit()          << "\t"
		 << cluster.Energy()        << "\t"
		 << cluster.YMin()          << "\t"
		 << cluster.YMax()          << "\t"
		 << clusterArr.at(clIndex).HitArr().size() << "\t"
		 << clusterArr.at(clIndex).nHit() << std::endl;
	*/
	//std::cout<< clusterArr.at(clIndex).HitArr().size() << std::endl;
	clusterID++;
      }else{ break; }
    }
    /*
    for ( int ihit = 0; ihit < clusterArr.size(); ihit++){
      std::cout << ihit << "\t"
		<< clusterArr.at(ihit).HitArr().size() << "\t"
		<< clusterArr.at(ihit).nHit() << std::endl;
    }
    */

    if( hitArr.size() == 0){ break; }
  }
  //return clusterArr;
  return true;
}

Bool_t HitClustering( TClonesArray* HitArr, std::vector<TPCCluster> &clusterArr ){
  TClass* cl = HitArr->GetClass();
  if( strcmp( "TPCPadHit", cl->GetName()) != 0 ){
    return false;
  }

  std::vector<TPCPadHit> hitArr;
  for( int i = 0; i< HitArr->GetEntries(); i++){
    TPCPadHit* hit = (TPCPadHit*)HitArr->At(i);
    hitArr.push_back(*hit);
  }
  HitClustering( hitArr, clusterArr );
  return true;
}
Bool_t HitClustering( std::vector<TPCPadHit>& hitArr, std::vector<TPCCluster> &clusterArr ){
  clusterArr.clear();
  Int_t ClusterID = 0;
  Int_t clIndex = 0;
  for( int irow = 4; irow < 32; irow++){
    //std::cout<< "irow " << irow << std::endl;
    Int_t nRow = 0;
    for( int i = 0; i< hitArr.size(); i++){
      if( hitArr.at(i).Row() == irow){ nRow++; }
    }
    while(nRow > 0 ){
      std::vector<TPCPadHit> sub_hitArr;
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
      TPCCluster cluster(ClusterID, Y, TotalEnergy);
      cluster.ID       = ClusterID;
      cluster.Row      = irow;
      cluster.YMin     = YMin;
      cluster.YMax     = YMax;
      cluster.Col      = ColID;
      cluster.ColMin   = ColMin;
      cluster.ColMax   = ColMax;
      cluster.PhiMin   = PhiMin;
      cluster.PhiMax   = PhiMax;
      cluster.Position = TVector3(X,Y,Z);
      cluster.NHit     = sub_hitArr.size();
      cluster.Mother   = -1;
      cluster.Daughter = -1;
      cluster.MinMotherDist = -1;
      cluster.MinDaughterDist = -1;
      cluster.nDaughter= 0;
      cluster.nMother  = 0;
      cluster.MotherID.clear();
      cluster.DaughterID.clear();
      cluster.Blocked = false;
      //std::cout<< PhiMin << "\t" << cluster.PhiMin << "\t" << PhiMax << "\t" << cluster.PhiMax << std::endl;
      /*
      std::cout<< cluster.ID << "\t"
	       << cluster.Mother << "\t"
	       << cluster.Daughter << "\t"
	       << cluster.nMother << "\t"
	       << cluster.nDaughter << "\t"
	       << cluster.MotherID.size() << "\t"
	       << cluster.DaughterID.size() << "\t"
	       << cluster.Energy << "\t"
	       << cluster.PhiMin << "\t"
	       << cluster.PhiMax << "\t"
	       << cluster.NHit << std::endl;
      */
      clusterArr.push_back( cluster );
      ClusterID++;

      nRow = 0;
      for( int i = 0; i< hitArr.size(); i++){
	if( hitArr.at(i).Row() == irow){ nRow++; }
      }
      if( ClusterID >= 200 ){ return false; }
    }
  }



  ////// connection
  for( Int_t icl = 0; icl < clusterArr.size(); icl++){
    Int_t cRow = clusterArr.at( icl ).Row;
    Int_t Mother = clusterArr.at(icl).ID;

    //// Find Daughter at next row
    Bool_t bFindNext = false;
    for( Int_t iloop = 0; iloop < 3; iloop++){
      for( Int_t jcl =0; jcl < clusterArr.size(); jcl++){
	if( icl != jcl ){ continue; }
	Int_t dRow  = clusterArr.at(jcl).Row;
	if( cRow != dRow - 1 ){ continue; }
	if( !AdjCluster( clusterArr.at(icl), clusterArr.at(jcl),1) ){
	  /// not connected
	continue;
	}else{
	  /// connected
	  TVector3 dp  = clusterArr.at(icl).Position - clusterArr.at(jcl).Position;
	  Double_t dist = dp.Mag();
	  clusterArr.at(icl).DaughterID.push_back( clusterArr.at(jcl).ID);
	  clusterArr.at(jcl).MotherID.push_back( clusterArr.at(icl).ID);
	  clusterArr.at(icl).nDaughter = clusterArr.at(icl).DaughterID.size();
	  clusterArr.at(jcl).nMother   = clusterArr.at(jcl).MotherID.size();
	  if( clusterArr.at(icl).MinDaughterDist > 0 &&
	      clusterArr.at(icl).MinDaughterDist > dist ){
	    clusterArr.at(icl).MinDaughterDist = dist;
	  }
	  if( clusterArr.at(jcl).MinMotherDist > 0 &&
	      clusterArr.at(jcl).MinMotherDist > dist ){
	    clusterArr.at(jcl).MinMotherDist = dist;
	  }
	  bFindNext = true;
	}
      }
      if( bFindNext ){ break; }
    }
  }
  /*
  for( int i = 0; i< clusterArr.size(); i++){
    std::cout<< clusterArr.at(i).ID << "\t"
	     << clusterArr.at(i).Mother << "\t"
	     << clusterArr.at(i).Daughter << "\t"
	     << clusterArr.at(i).nMother << "\t"
	     << clusterArr.at(i).nDaughter << "\t"
	     << clusterArr.at(i).MotherID.size() << "\t"
	     << clusterArr.at(i).DaughterID.size() << "\t"
	     << clusterArr.at(i).Energy << "\t"
	     << clusterArr.at(i).PhiMin << "\t"
	     << clusterArr.at(i).PhiMax << "\t"
	     << clusterArr.at(i).NHit << std::endl;
  }
  */
  return true;
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
     ( TMath::Abs( c_1.Position.Y() - c_0.Position.Y()) < 60 )){
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


void ConnectionFinder( std::vector<TPCCluster> &clArr ){

  for( Int_t icl = 0; icl < clArr.size(); icl++){
    Int_t cRow = clArr.at( icl ).Row;
    Int_t Mother = clArr.at(icl).ID;

    //// Find Daughter at next row
    Bool_t bFindNext = false;
    for( Int_t iloop = 0; iloop < 3; iloop++){
      for( Int_t jcl =0; jcl < clArr.size(); jcl++){
	if( icl != jcl ){ continue; }
	Int_t dRow  = clArr.at(jcl).Row;
	if( cRow != dRow - 1 ){ continue; }
	if( !AdjCluster( clArr.at(icl), clArr.at(jcl),1) ){
	  /// not connected
	continue;
	}else{
	  /// connected
	  TVector3 dp  = clArr.at(icl).Position - clArr.at(jcl).Position;
	  Double_t dist = dp.Mag();
	  clArr.at(icl).DaughterID.push_back( clArr.at(jcl).ID);
	  clArr.at(jcl).MotherID.push_back( clArr.at(icl).ID);
	  clArr.at(icl).nDaughter = clArr.at(icl).DaughterID.size();
	  clArr.at(jcl).nMother   = clArr.at(jcl).MotherID.size();
	  if( clArr.at(icl).MinDaughterDist > 0 &&
	      clArr.at(icl).MinDaughterDist > dist ){
	    clArr.at(icl).MinDaughterDist = dist;
	  }
	  if( clArr.at(jcl).MinMotherDist > 0 &&
	      clArr.at(jcl).MinMotherDist > dist ){
	    clArr.at(jcl).MinMotherDist = dist;
	  }
	  bFindNext = true;
	}
      }
      if( bFindNext ){ break; }
    }
  }
}
std::vector< std::vector<TPCCluster> > ConnectionBlocker( std::vector<TPCCluster> &clArr ){
  std::vector< std::vector<TPCCluster> > BlockList;
  if( clArr.size() > 200 ) { return BlockList; }

  //std::cout<< "Collect root events " << std::endl;
  std::vector<Int_t> sidList;
  for( int i = 0; i< clArr.size(); i++){
    if( clArr.at(i).nDaughter > 1 ){
      //std::cout<< clArr.at(i).DaughterID.size() << std::endl;
      for( int IDsize = 0; IDsize < clArr.at(i).DaughterID.size(); IDsize++){
	sidList.push_back( clArr.at(i).DaughterID.at(IDsize));
      }
    }else if( clArr.at(i).nMother == 0){
      sidList.push_back( clArr.at(i).ID );
    }
  }
  //std::cout<< "Blocking " << std::endl;
  for( int id = 0; id < sidList.size(); id++){
    std::vector<TPCCluster> clBlock;
    Int_t currentID = sidList.at(id);
    while( clArr.at(currentID).nDaughter < 2 ){
      clBlock.push_back( clArr.at(currentID) );
      if( clArr.at(currentID).nDaughter == 0 ){ continue; }
      currentID = clArr.at(currentID).Daughter;
    }
    if( clBlock.size() == 0 ){ continue; }
    BlockList.push_back( clBlock );
  }
  //std::cout<< "End Blocking" << std::endl;
  return BlockList;
}

Bool_t HitClustering( TClonesArray* HitArr, TClonesArray* clusterArr ){
  //clusterArr.clear();
  clusterArr->Clear();

  std::vector<TPCPadHit> hitArr;
  for( int i = 0; i< HitArr->GetEntries(); i++){
    TPCPadHit* hit = (TPCPadHit*)HitArr->At(i);
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
      std::vector<TPCPadHit> sub_hitArr;
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


std::vector<TPCCluster> TPCClusterBlocker( std::vector<TPCCluster> &clArr, Int_t Direction ){
  std::vector<TPCCluster> clBlock;
  Int_t ID=-1;
  Int_t minimumID=32;
  Int_t MaximumID= -1;
  std::vector<int> BlockIndexArr;

  switch( Direction ){
  case 0:
    // Blocking inner to outer

    for( Int_t i = 0; i< clArr.size(); i++){
      if( clArr.at(i).Row < minimumID ){
	minimumID = clArr.at(i).Row;
	ID = i;
      }
    }
    clBlock.push_back( clArr.at(ID) );
    clArr.erase( clArr.begin() + ID );
    for( int icl = 0; icl < clBlock.size(); icl++){
      Bool_t bFind = false;
      if( bFind ){ break; }
      for( Int_t jcl = 0; jcl< clArr.size();){
	if( AdjCluster( clBlock.at(icl), clArr.at(jcl),1)){
	  clArr.at(jcl).Mother = clBlock.at(icl).ID;
	  clArr.at(jcl).MotherID.push_back( clBlock.at(icl).ID);
	  clBlock.at(icl).nDaughter = clBlock.at(icl).nDaughter +1;
	  clBlock.at(icl).DaughterID.push_back( clArr.at(jcl).ID);
	  clBlock.push_back(clArr.at(jcl));
	  clArr.erase(clArr.begin() + jcl );
	  bFind = true;
	}else{
	  jcl++;
	}
      }
    }
    break;
    /*
      case 1:
    // Blocking outer to inner
    for( Int_t i = 0; i< clArr.size(); i++){
      if( clArr.at(i).Row > MaximumID ){
	MaximumID = clArr.at(i).Row;
	ID = i;
      }
    }

    clBlock.push_back( clArr.at(ID) );
    clArr.erase( clArr.begin() + ID );

    for( int icl = 0; icl < clBlock.size(); icl++){
      Bool_t bFindFirst  = false;
      Bool_t bFindSecond = false;
      Bool_t bFindThird  = false;
      Bool_t bFind = false;
      for( Int_t iRepeat = 1; iRepeat <= 3; iRepeat++){
	if( bFind ){ break; }
	for( Int_t jcl = 0; jcl< clArr.size();){
	  if( AdjCluster( clArr.at(icl), clBlock.at(jcl),iRepeat)){
	    clArr.at(jcl).Mother = clBlock.at(icl).ID;
	    clArr.at(jcl).MotherID.push_back( clBlock.at(icl).ID);
	    clBlock.at(icl).nDaughter = clBlock.at(icl).nDaughter +1;
	    clBlock.at(icl).DaughterID.push_back( clArr.at(jcl).ID);
	    clBlock.push_back(clArr.at(jcl));
	    clArr.erase(clArr.begin() + jcl );
	    bFind = true;
	  }else{
	    jcl++;
	  }
	}
      }
    }
    break;
    */

  case 2:

    for( Int_t i = 0; i< clArr.size(); i++){
      if( clArr.at(i).Row < minimumID ){
	minimumID = clArr.at(i).Row;
	ID = i;
      }
    }
    clBlock.push_back( clArr.at(ID) );
    clArr.erase( clArr.begin() + ID );
    for( int Row = minimumID+1; Row < 32; Row++){
      //std::cout<< Row << std::endl;
      for( int jcl = 0; jcl < clArr.size();){
	if( clArr.at(jcl).Row != Row ){ jcl++; continue; }
	Int_t dRow = 1;
	Bool_t bFind = false;
	for( int icl = 0; icl < clBlock.size(); icl++){
	  if( clBlock.at(icl).Row != clArr.at(jcl).Row - dRow ){ continue; }
	  //std::cout << icl << "/" << clBlock.size() << std::endl;
	  if( AdjCluster( clBlock.at(icl), clArr.at(jcl), dRow) ){
	    clArr.at(jcl).Mother = clBlock.at(icl).ID;
	    clArr.at(jcl).MotherID.push_back( clBlock.at(icl).ID);
	    clBlock.at(icl).nDaughter = clBlock.at(icl).nDaughter +1;
	    clBlock.at(icl).DaughterID.push_back( clArr.at(jcl).ID);
	    clBlock.push_back(clArr.at(jcl));
	    clArr.erase( clArr.begin()+jcl);
	    bFind = true;
	    break;
	  }
	}

	if( bFind ){ continue; }

	dRow = 2;
	for( int icl = 0; icl < clBlock.size(); icl++){
	  if( clBlock.at(icl).Row != clArr.at(jcl).Row - dRow ){ continue; }
	  if( AdjCluster( clBlock.at(icl), clArr.at(jcl), dRow) ){
	    clArr.at(jcl).Mother = clBlock.at(icl).ID;
	    clArr.at(jcl).MotherID.push_back( clBlock.at(icl).ID);
	    clBlock.at(icl).nDaughter = clBlock.at(icl).nDaughter +1;
	    clBlock.at(icl).DaughterID.push_back( clArr.at(jcl).ID);
	    clBlock.push_back(clArr.at(jcl));
	    clArr.erase( clArr.begin()+jcl);
	    bFind = true;
	    break;
	  }
	}
	if( bFind ){ continue; }

	for( int icl = 0; icl < clBlock.size(); icl++){
	  if( clBlock.at(icl).Row != clArr.at(jcl).Row - dRow ){ continue; }
	  if( AdjCluster( clBlock.at(icl), clArr.at(jcl), dRow) ){
	    clArr.at(jcl).Mother = clBlock.at(icl).ID;
	    clArr.at(jcl).MotherID.push_back( clBlock.at(icl).ID);
	    clBlock.at(icl).nDaughter = clBlock.at(icl).nDaughter +1;
	    clBlock.at(icl).DaughterID.push_back( clArr.at(jcl).ID);
	    clBlock.push_back(clArr.at(jcl));
	    clArr.erase( clArr.begin()+jcl);
	    bFind = true;
	    break;
	  }
	}
	if( bFind ){ continue; }
	jcl++;
      }
    }

    for( int jcl = 0; jcl < clArr.size();){
      for( int icl = 0; icl < clBlock.size(); icl++){

      }
    }

  default:
    break;
  }
  return clBlock;
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





std::vector<TPCCluster> TPCClusterSingleBlocker( std::vector<TPCCluster> &clArr){
  std::vector<TPCCluster> clBlock;
  Int_t ID=-1;
  Int_t minimumID=32;
  Int_t MaximumID= -1;
  std::vector<int> BlockIndexArr;

  for( Int_t i = 0; i< clArr.size(); i++){
    if( clArr.at(i).Row < minimumID ){
      minimumID = clArr.at(i).Row;
      ID = i;
    }
  }

  clBlock.push_back( clArr.at(ID) );
  clArr.erase( clArr.begin() + ID );

  for( int icl = 0; icl < clBlock.size(); icl++){
    Bool_t bFind = false;
    // find next clusters //
    Int_t nFind = 0;
    for( Int_t iRepeat = 1; iRepeat <= 3; iRepeat++){
      if( bFind ){break;}
      nFind = 0;
      for( Int_t jcl = 0; jcl < clArr.size(); jcl++){
	if( AdjCluster( clBlock.at(icl), clArr.at(jcl),iRepeat)){
	  bFind = true;
	  nFind++;
	}
      }
    }
    if( nFind > 1 ){// branched track
      break;
    }
    bFind = false;
    for( Int_t iRepeat = 1; iRepeat <= 3; iRepeat++){
      if( bFind ){ break; }
      for( Int_t jcl = 0; jcl< clArr.size();){
	if( AdjCluster( clBlock.at(icl), clArr.at(jcl),iRepeat)){
	  clArr.at(jcl).Mother = clBlock.at(icl).ID;
	  clArr.at(jcl).MotherID.push_back(clBlock.at(icl).ID);
	  clBlock.at(icl).nDaughter = clBlock.at(icl).nDaughter +1;
	  clBlock.at(icl).DaughterID.push_back(clArr.at(jcl).ID);
	  clBlock.push_back(clArr.at(jcl));
	  clArr.erase(clArr.begin() + jcl );
	  bFind = true;
	}else{
	  jcl++;
	}
      }
    }
  }

  return clBlock;
}

Int_t BlockDivider( std::vector<TPCCluster> &clArr, std::vector<TPCCluster> &trackCand, Int_t Mode ){
  Int_t trackSize=0;
  Int_t MotherID = -1;
  Bool_t bBranch = false;
  switch( Mode ){
  case 0 : // Edge to center stop at branching point
    //// find edge ////
    for( int icl = 0; icl < clArr.size(); icl++ ){
      if( clArr.at(icl).nDaughter == 0){
	MotherID = clArr.at(icl).Mother;
	trackCand.push_back( clArr.at(icl) );
	clArr.erase( clArr.begin() + icl);
	trackSize++;
	break;
      }
    }
    if( trackSize == 0 ){ break; }
    /// tracking mother ///
    bBranch = false;
    for( int icl = 0; icl < trackCand.size(); icl++){
      if( trackCand.at(icl).Mother == -1 ){ break; }
      if( bBranch ){ break; }
      for( int jcl = 0; jcl< clArr.size(); ){
	if( clArr.at(jcl).ID == MotherID){
	  MotherID = clArr.at( jcl ).Mother;
	  trackCand.push_back( clArr.at(jcl) );
	  trackSize++;
	  if( clArr.at( jcl ).nDaughter != 1 ){
	    bBranch = true;
	    jcl++;
	  }else{
	    clArr.erase( clArr.begin() + jcl );
	  }
	  break;
	}else{
	  jcl++;
	}
      }
    }
    break;
  case 1 : // branchingPoint to Center
    //Searching branching Point ( non ROOT point )
    for( int icl = 0; icl < clArr.size(); icl++){
      if(clArr.at(icl).nDaughter > 1 && clArr.at(icl).Mother != -1 ){
	MotherID = clArr.at(icl).Mother;
	trackCand.push_back( clArr.at(icl));
	clArr.erase( clArr.begin() + icl );
	trackSize++;
	break;
      }
    }
    if( trackSize == 0){ break; }

    // tracking Mother
    bBranch = false;
    for( int icl = 0; icl < trackCand.size(); icl++){
      if( trackCand.at(icl).Mother == -1 ){ break; }
      if( bBranch ){ break; }
      for( int jcl = 0; jcl< clArr.size(); ){
	if( clArr.at(jcl).ID == MotherID){
	  MotherID = clArr.at( jcl ).Mother;
	  trackCand.push_back( clArr.at(jcl) );
	  trackSize++;
	  if( clArr.at( jcl ).nDaughter != 1 ){
	    bBranch = true;
	    jcl++;
	  }else{
	    clArr.erase( clArr.begin() + jcl );
	  }
	  break;
	}else{
	  jcl++;
	}
      }
    }
    break;
  default:
    break;
  }
  return trackSize;
}
