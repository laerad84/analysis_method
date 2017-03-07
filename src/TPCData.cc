#include "TPCData.h"
#include "TVector3.h"

ClassImp( TPCPadHit )

TPCPadHit::TPCPadHit(){
  m_TrackID   = -1;
  m_ClusterID = -1;
  m_HitID     = -1;
  m_PadID     = -1;
  m_MCTrackID = -1;
  m_Row       = -1;
  m_Col       = -1;
  m_pos = TVector3(0,0,0);
  m_posErr = TVector3(0,0,0);
  m_Phi=0.;
  m_Rad=0.;
  m_energy = 0;
  m_signal = 0;
  m_rawTiming = 0;
  m_timing = 0;

  IsClustered = false;
  IsTracked   = false;
}
TPCPadHit::TPCPadHit(const TPCPadHit& right):
  m_TrackID(right.m_TrackID),
  m_MCTrackID(right.m_MCTrackID),
  m_ClusterID(right.m_ClusterID),
  m_HitID(right.m_HitID),
  m_PadID(right.m_PadID),
  m_pos(right.m_pos),
  m_posErr(right.m_posErr),
  m_Row(right.m_Row),
  m_Col(right.m_Col),
  m_Phi(right.m_Phi),
  m_Rad(right.m_Rad),
  m_energy(right.m_energy),
  m_signal(right.m_signal),
  m_rawTiming(right.m_rawTiming),
  m_timing(right.m_timing),
  IsClustered(right.IsClustered),
  IsTracked(right.IsTracked)
{
  ;
}
TPCPadHit::TPCPadHit(int padID, TVector3 pos, double energy ) :
  m_TrackID(-1),
  m_ClusterID(-1),
  m_HitID(-1),
  m_PadID(padID),
  m_MCTrackID(-1),
  m_Row(gTPCIDHandler->GetRow(padID)),
  m_Col(gTPCIDHandler->GetCol(padID)),
  m_pos(pos),
  m_posErr(TVector3(0,0,0)),
  m_energy(energy),
  m_signal(0.),
  m_Phi(0.),
  m_Rad(0.),
  m_rawTiming(0),
  m_timing(0),
  IsClustered(false),
  IsTracked(false){
  /*
    m_TrackID = -1;
    m_ClusterID = -1;
    m_HitID = -1;
    m_PadID = padID;
    m_MCTrackID = -1;

    m_Row = gTPCIDHandler->GetRow(padID);
    m_Col  = gTPCIDHandler->GetCol(padID);
    m_pos =  pos;
    m_posErr = TVector3(0,0,0);

    m_energy = energy;
    m_signal = 0;
    m_rawTiming = 0;
    m_timing = 0;

    IsClustered = false;
    IsTracked   = false;
  */

}
TPCPadHit::~TPCPadHit(){
  Clear();
}
void TPCPadHit::Clear(Option_t* opt){

  m_TrackID   = -1;
  m_ClusterID = -1;
  m_HitID     = -1;
  m_PadID     = -1;
  m_MCTrackID = -1;
  m_Row       = -1;
  m_Col       = -1;
  m_pos = TVector3(0,0,0);
  m_posErr = TVector3(0,0,0);

  m_energy = 0;
  m_signal = 0;
  m_rawTiming = 0;
  m_timing = 0;

  IsClustered = false;
  IsTracked   = false;

}

ClassImp(TPCPadHitCluster)
TPCPadHitCluster::TPCPadHitCluster(){
  m_TrackID = -1;
  m_ClusterID = -1;
  m_Row = -1;
  m_YMin      = 999;
  m_YMax      = -999;
  m_ColMax    = -1;
  m_ColMin    = 999;
  m_nHit = 0;
  m_TotalEnergy =0;
  m_ClusterPos=TVector3(0,0,0);
  m_ClusterPosErr=TVector3(0,0,0);
  m_MotherID = -1;

  IsOuterCluster = false;
  IsTracked = false;
  IsUpdated = false;
}
TPCPadHitCluster::TPCPadHitCluster(int rowID):
  m_TrackID(-1),
  m_ClusterID(-1),
  m_Row(rowID),
  m_YMin     (999),
  m_YMax     (-999),
  m_ColMax   (-1),
  m_ColMin   (999),
  m_nHit(0),
  m_TotalEnergy(0),
  m_ClusterPos(TVector3(0,0,0)),
  m_ClusterPosErr(TVector3(0,0,0)),
  m_MotherID(-1),
  IsOuterCluster(false),
  IsTracked(false),
  IsUpdated(false)
{
  ;
}
TPCPadHitCluster::TPCPadHitCluster( const TPCPadHitCluster& right):
  m_TrackID(right.m_TrackID),
  m_ClusterID(right.m_ClusterID),
  m_Row(right.m_Row),
  m_Col(right.m_Col),
  m_ColMin(right.m_ColMin),
  m_ColMax(right.m_ColMax),
  m_MotherID(right.m_MotherID),
  m_nHit(right.m_MotherID),
  m_ClusterPos(right.m_ClusterPos),
  m_ClusterPosErr(right.m_ClusterPosErr),
  m_TotalEnergy(right.m_TotalEnergy),
  IsOuterCluster(right.IsOuterCluster),
  IsTracked(right.IsOuterCluster),
  IsUpdated(right.IsUpdated)
{
  ;
}


void TPCPadHitCluster::Clear(Option_t* opt){
  m_TrackID = -1;
  m_ClusterID = -1;
  m_Row = -1;
  m_YMin      = 999;
  m_YMax      = -999;
  m_ColMax    = -1;
  m_ColMin    = 999;
  m_nHit = 0;
  m_TotalEnergy =0;
  m_ClusterPos=TVector3(0,0,0);
  m_ClusterPosErr=TVector3(0,0,0);
  m_MotherID = -1;
  IsOuterCluster = false;
  IsTracked = false;
  IsUpdated = false;
}
TPCPadHitCluster::~TPCPadHitCluster(){
  m_HitArr.clear();
}
void TPCPadHitCluster::AddPadHit( TPCPadHit hit ){
  m_nHit++;
  m_HitArr.push_back( hit );
  m_ColArr.push_back( hit.Col() );
  m_EArr.push_back( hit.Energy() );
}
void TPCPadHitCluster::Evaluate(){
  //m_nHit = m_HitArr.size();
  m_Row  = m_HitArr.at(0).Row();
  Double_t meanColID = 0;
  Double_t rmsColID  = 0;
  Double_t meanY     = 0;
  Double_t rmsY      = 0;
  m_TotalEnergy      = 0;
  m_ColMin = 999;
  m_ColMax = -1;
  m_YMin   = 999;
  m_YMax   = -999;
  m_PhiMax = -999;
  m_PhiMin = 999;

  /// Check size
  if( m_HitArr.size() != m_ColArr.size() ||
      m_HitArr.size() != m_EArr.size() ){
    std::cerr << "Array size dose not Matched!!!" << std::endl;
    std::cerr << __PRETTY_FUNCTION__ << " : " << __LINE__ << std::endl;
  }

  //m_nHit = m_HitArr.size();
  for( int ihit = 0; ihit < m_HitArr.size(); ihit++){
    meanColID     += m_HitArr.at(ihit).Col()*m_HitArr.at(ihit).Energy();
    m_TotalEnergy += m_HitArr.at(ihit).Energy();
    meanY         += m_HitArr.at(ihit).Position().Y()*m_HitArr.at(ihit).Energy();
    if( m_YMin > m_HitArr.at(ihit).Position().Y() ){ m_YMin = m_HitArr.at(ihit).Position().Y(); }
    if( m_YMax < m_HitArr.at(ihit).Position().Y() ){ m_YMax = m_HitArr.at(ihit).Position().Y(); }
    if( m_ColMin > m_HitArr.at(ihit).Col() ){ m_ColMin = m_HitArr.at(ihit).Col(); }
    if( m_ColMax < m_HitArr.at(ihit).Col() ){ m_ColMax = m_HitArr.at(ihit).Col(); }
    if( m_PhiMin > m_HitArr.at(ihit).Phi() ){ m_PhiMin = m_HitArr.at(ihit).Phi(); }
    if( m_PhiMax < m_HitArr.at(ihit).Phi() ){ m_PhiMax = m_HitArr.at(ihit).Phi(); }
  }

  //std::cout<< "cluster Phi: "<< m_Row << "\t" << m_PhiMin << "\t" << m_PhiMax << std::endl;

  meanColID = meanColID / m_TotalEnergy;
  meanY     = meanY / m_TotalEnergy;

  for( int ihit = 0; ihit< m_HitArr.size(); ihit++){
    rmsColID += TMath::Power((m_HitArr.at(ihit).Col() - meanColID),2)*m_EArr.at(ihit);
  }
  rmsColID = TMath::Sqrt(rmsColID/m_TotalEnergy);

  ///// cluster position calculation /////
  m_Col = meanColID;
  double ColthetaDeg = (-180 + TPCGlobals::sTPC_Pad_Parameter[m_Row][4] + (meanColID+0.5)*TPCGlobals::sTPC_Pad_Parameter[m_Row][3]);
  double rad         = TPCGlobals::sTPC_Pad_Parameter[m_Row][2];
  double xPos        = rad * TMath::Sin( ColthetaDeg * TMath::DegToRad() );
  double zPos        = TPCGlobals::sTPC_zOffset + rad * TMath::Cos( ColthetaDeg * TMath::DegToRad() );
  double yPos        = meanY;

  ///// cluster position error calculation /////
  double padHalfLength = TPCGlobals::sTPC_Pad_Parameter[m_Row][5]/2.;
  double thetaErr      = rmsColID*TPCGlobals::sTPC_Pad_Parameter[m_Row][3]*TMath::DegToRad()*rad;
  double rErr          = padHalfLength;
  double theta         = ColthetaDeg*TMath::DegToRad();

  double xErr          = TMath::Sqrt(TMath::Power(rErr*TMath::Sin(theta), 2) + TMath::Power(thetaErr*TMath::Cos(theta),2));
  double yErr          = rmsY;
  double zErr          = TMath::Sqrt(TMath::Power(rErr*TMath::Cos(theta), 2) + TMath::Power(thetaErr*TMath::Sin(theta),2));

  /*
    std::cout<< "Hit Size: " << m_HitArr.size() << "Col : " << meanColID << std::endl;
    for( int ihit = 0; ihit <m_HitArr.size(); ihit++){
    std::cout<<m_HitArr.at(ihit).Row() << "\t" <<  m_HitArr.at(ihit).Col() << "\t" << m_HitArr.at(ihit).Energy()<< std::endl;
    }
  */
  m_ClusterPos    = TVector3( xPos, yPos, zPos );
  m_ClusterPosErr = TVector3( xErr, yErr, zErr );

}

TPCPadHitCluster TPCClusterer( int rowID, std::vector<TPCPadHit> &hitArr){
  //std::cout<< __PRETTY_FUNCTION__ << std::endl;
  bool isGood = false;
  TPCPadHitCluster cluster;
  //std::vector<TPCPadHit> subHitArr;
  /// Simple version : only check RowID;
  Int_t clusterHitSize = 0;
  for( int ihit = 0; ihit < hitArr.size(); ihit++){
    if( hitArr.at(ihit).Row() == rowID ){
      if( cluster.nHit() == 0 ){
	//subHitArr.push_back(hitArr.at(ihit));
	cluster.AddPadHit( hitArr.at(ihit) );
	hitArr.erase( hitArr.begin()+ihit);
	break;
      }
    }
  }
  //if( subHitArr.size() == 0 ){ return cluster; }
  if( cluster.HitArr().size()==0 ){ return cluster; }

  for( int iclhit  = 0; iclhit < cluster.HitArr().size(); iclhit++){
    for( int ihit = 0; ihit < hitArr.size(); ){
      if( TMath::Abs(hitArr.at(ihit).Col() - cluster.HitArr().at(iclhit).Col()) < 2 &&
	  TMath::Abs(hitArr.at(ihit).Position().Y() - cluster.HitArr().at(iclhit).Position().Y()) < 10 &&
	  hitArr.at(ihit).Row() == rowID ){
	/// Add Y direction decision
	//subHitArr.push_back(hitArr.at(ihit));
	cluster.AddPadHit(hitArr.at(ihit));
	hitArr.erase( hitArr.begin() + ihit);
      }else{
	ihit++;
      }
    }
  }

  //std::cout<< "clusterSize :" << "\t" << cluster.HitArr().size() << std::endl;
  //std::cout<< "clusterSizeCK : "  << cluster.HitArr().size() << std::endl;
  cluster.Evaluate();

  return cluster;
}
TPCPadHitCluster TPCClusterer( int rowID, std::vector<TPCPadHit> &hitArr, Int_t MCTrackID){
  bool isGood = false;
  TPCPadHitCluster cluster( rowID );
  for( int ihit = 0; ihit < hitArr.size(); ihit++){
    if( hitArr.at(ihit).Row() == rowID &&
	hitArr.at(ihit).MCTrackID() == MCTrackID ){
      cluster.AddPadHit( hitArr.at(ihit) );
      hitArr.erase( hitArr.begin()+ihit);
      break;
    }
  }

  if( cluster.HitArr().size() == 0 ) { return cluster; }
  for( int iclhit = 0; iclhit < cluster.ColIDArr().size(); iclhit++){
    for( int ihit = 0; ihit < hitArr.size(); ){
      if( TMath::Abs(hitArr.at(ihit).Col() - cluster.HitArr().at(iclhit).Col()) < 2 && /// adjected pad
	  TMath::Abs(hitArr.at(ihit).Position().Y() - cluster.HitArr().at(iclhit).Position().Y() ) < 20 &&
	  hitArr.at(ihit).Row()       == rowID &&
	  hitArr.at(ihit).MCTrackID() == MCTrackID ){/// 20mm in y dist
	cluster.AddPadHit( hitArr.at( ihit ));
	hitArr.erase( hitArr.begin() + ihit );
      }else{
	ihit++;
      }
    }
  }

  cluster.Evaluate();
  return cluster;
}

/*
ClassImp(TPCTrack)
TPCTrack::TPCTrack(){
  m_TrackID = -1;
  m_nCluster = 0;
  m_circX = 0;
  m_circZ = 0;
  m_circR = 0;
  m_circXErr = 0;
  m_circZErr = 0;
  m_circRErr = 0;
  m_circFitChisq = 0;
  m_circFit = kFALSE;
  m_ClusterArr.clear();
}
TPCTrack::~TPCTrack(){
  m_ClusterArr.clear();
}
double TPCTrack::FitCircle(){

  double chi;
  return chi;
}
*/

std::vector<TPCPadHitCluster> TPCTrackBlocking( std::vector<TPCPadHitCluster> &clusterArr, Int_t Direction){
  //find Minimum Row clusterArr
  std::cout<< "TPCTrackBlocking" << std::endl;
  std::vector<TPCPadHitCluster> clusterBlock;
  if( Direction == 0 ){// Inner -> Outer
    std::cout<< "Inner" << std::endl;
    Int_t minimumID = 32;
    Int_t ID=-1;
    for( Int_t i = 0; i< clusterArr.size(); i++){
      if( clusterArr.at(i).RowID() < minimumID ){
        minimumID = clusterArr.at(i).RowID();
        ID = i;
      }
    }
    std::vector<int> BlockIndexArr;
    clusterBlock.push_back(clusterArr.at(ID));
    clusterArr.erase(clusterArr.begin() + ID);
    // find cluster in XZ
    for( Int_t icl = 0; icl < clusterBlock.size(); icl++){
      // find Cluster in next Row
      std::cout<< icl << "\t" << clusterBlock.size() << std::endl;
      Bool_t bFindFirst  = false;
      Bool_t bFindSecond = false;
      Bool_t bFindThird  = false;
      for( Int_t i = 0; i< clusterArr.size();){
	//std::cout<< clusterBlock.at(icl).PhiMin() << "\t" << clusterBlock.at(icl).PhiMax()
	//<< clusterArr.at(icl).PhiMin() << "\t" << clusterArr.at(icl).PhiMax() << std::endl;

	if( AdjacentCluster(clusterBlock.at(icl), clusterArr.at(i), 1)){
	    clusterArr.at(i).SetMother(clusterBlock.at(icl).ClusterID());
	    clusterBlock.at(icl).AddDaughter(clusterArr.at(i).ClusterID());
	    clusterBlock.push_back(clusterArr.at(i));
	    clusterArr.erase(clusterArr.begin() + i );
	    bFindFirst = true;
	}else{
	  i++;
	}
      }
      if( !bFindFirst ){
	for( Int_t i = 0; i< clusterArr.size();){
	  if( AdjacentCluster(clusterBlock.at(icl), clusterArr.at(i), 2)){
	      clusterArr.at(i).SetMother(clusterBlock.at(icl).ClusterID());
	      clusterBlock.at(icl).AddDaughter(clusterArr.at(i).ClusterID());
	      clusterBlock.push_back(clusterArr.at(i));
	      clusterArr.erase(clusterArr.begin() + i );
	      bFindFirst = true;
	    }else{
	      i++;
	  }
	}
	if( !bFindSecond ){
	  for( Int_t i = 0; i< clusterArr.size();){
	    if( AdjacentCluster(clusterBlock.at(icl), clusterArr.at(i), 3)){
	      clusterArr.at(i).SetMother(clusterBlock.at(icl).ClusterID());
	      clusterBlock.at(icl).AddDaughter(clusterArr.at(i).ClusterID());
	      clusterBlock.push_back(clusterArr.at(i));
	      clusterArr.erase(clusterArr.begin() + i );
	      bFindFirst = true;
	    }else{
	      i++;
	    }
	  }
	}
      }
    }
  }else if( Direction == 1 ){// Outer -> Inner
    Int_t maximumID = -1;
    Int_t ID=-1;
    std::vector<TPCPadHitCluster> clusterBlock;
    for( Int_t i = 0; i< clusterArr.size(); i++){
      if( clusterArr.at(i).RowID() > maximumID ){
        maximumID = clusterArr.at(i).RowID();
        ID = i;
      }
    }
    std::vector<int> BlockIndexArr;
    clusterBlock.push_back(clusterArr.at(ID));
    clusterArr.erase(clusterArr.begin() + ID);
    // find cluster in XZ
    for( Int_t icl = 0; icl < clusterBlock.size(); icl++){
      // find Cluster in next Row
      Bool_t bFindFirst  = false;
      Bool_t bFindSecond = false;
      Bool_t bFindThird  = false;
      for( Int_t i = 0; i< clusterArr.size();){
	std::cout<< icl << "\t" << clusterBlock.size() << std::endl;
	if( AdjacentCluster(clusterArr.at(icl), clusterBlock.at(i), 1)){
	  clusterArr.at(i).SetMother(clusterBlock.at(icl).ClusterID());
	  clusterBlock.at(icl).AddDaughter(clusterArr.at(i).ClusterID());
	  clusterBlock.push_back(clusterArr.at(i));
	  clusterArr.erase(clusterArr.begin() + i );
	  bFindFirst = true;
	}else{
	  i++;
	}
      }
      if( !bFindFirst ){
	for( Int_t i = 0; i< clusterArr.size();){
	  if( AdjacentCluster(clusterArr.at(icl), clusterBlock.at(i), 2)){
	    clusterArr.at(i).SetMother(clusterBlock.at(icl).ClusterID());
	    clusterBlock.at(icl).AddDaughter(clusterArr.at(i).ClusterID());
	    clusterBlock.push_back(clusterArr.at(i));
	    clusterArr.erase(clusterArr.begin() + i );
	    bFindFirst = true;
	}else{
	    i++;
	  }
	}
	if( !bFindSecond ){
	  for( Int_t i = 0; i< clusterArr.size();){
	    if( AdjacentCluster(clusterArr.at(icl), clusterBlock.at(i), 3)){
	      clusterArr.at(i).SetMother(clusterBlock.at(icl).ClusterID());
	      clusterBlock.at(icl).AddDaughter(clusterArr.at(i).ClusterID());
	      clusterBlock.push_back(clusterArr.at(i));
	      clusterArr.erase(clusterArr.begin() + i );
	      bFindFirst = true;
	    }else{
	      i++;
	    }
	  }
	}
      }
    }


  }
  return clusterBlock;
}

std::vector<TPCPadHitCluster> TPCEdgeClusters(std::vector<TPCPadHitCluster> clusterArr){
  std::vector<TPCPadHitCluster> semiEdgeCluster;
  std::vector<TPCPadHitCluster> EdgeCluster;
  for( int i = 0; i< clusterArr.size(); i++){
    if( TMath::Abs(clusterArr.at(i).GetClusterPos().Y()) < 250 &&
        TMath::Sqrt(TMath::Power(clusterArr.at(i).GetClusterPos().X(),2) +
		    TMath::Power(clusterArr.at(i).GetClusterPos().Z(),2)) < 220 ){
      ;
    }else{
      semiEdgeCluster.push_back(clusterArr.at(i));
    }
  }
  for( int i = 0; i< semiEdgeCluster.size(); i++){
    ///Searching outer ////
    Int_t nOuter     = 0;
    Double_t meanPhi = (semiEdgeCluster.at(i).PhiMax() + semiEdgeCluster.at(i).PhiMin())/2.;
    Double_t widthPhi= TMath::Abs(semiEdgeCluster.at(i).PhiMax() - semiEdgeCluster.at(i).PhiMin())/2.;
    Double_t meanY   =  semiEdgeCluster.at(i).GetClusterPos().Y();
    Double_t widthY  =  TMath::Abs(semiEdgeCluster.at(i).YMax() - semiEdgeCluster.at(i).YMin())/2.;

    for( int j = 0; j< semiEdgeCluster.size(); j++){
      if( i==j ){ continue; }
      Int_t RowDelta = semiEdgeCluster.at(j).RowID() - semiEdgeCluster.at(i).RowID();
      if( RowDelta <= 0 ){ continue; }
      TVector3 delta = semiEdgeCluster.at(j).GetClusterPos() - semiEdgeCluster.at(i).GetClusterPos();
      if( delta.Mag() > 70){ continue;}

      Double_t meanPhiSub = (semiEdgeCluster.at(j).PhiMax() + semiEdgeCluster.at(j).PhiMin())/2.;
      Double_t widthPhiSub= TMath::Abs(semiEdgeCluster.at(j).PhiMax() - semiEdgeCluster.at(j).PhiMin())/2.;
      Double_t meanYSub   =  semiEdgeCluster.at(j).GetClusterPos().Y();
      Double_t widthYSub  =  TMath::Abs(semiEdgeCluster.at(j).YMax() - semiEdgeCluster.at(j).YMin())/2.;

      if( TMath::Abs(meanPhiSub - meanPhi) > (widthPhi*(RowDelta))){ continue;}
      if( TMath::Abs(meanYSub - meanY) > RowDelta*70){continue;}
      nOuter++;
    }
    if( nOuter == 0 ){
      EdgeCluster.push_back(semiEdgeCluster.at(i));
    }
  }
  return EdgeCluster;
}
std::vector<TPCPadHitCluster> TPCRootClusters(std::vector<TPCPadHitCluster> clusterArr){
  std::vector<TPCPadHitCluster> RootCluster;
  for( int i = 0; i< clusterArr.size(); i++){
    if( clusterArr.at(i).RowID() < 5 ){ continue; }
    Int_t nRoot = 0;

    Double_t meanPhi = (clusterArr.at(i).PhiMax() + clusterArr.at(i).PhiMin())/2.;
    Double_t widthPhi= TMath::Abs(clusterArr.at(i).PhiMax() - clusterArr.at(i).PhiMin())/2.;
    Double_t meanY   =  clusterArr.at(i).GetClusterPos().Y();
    Double_t widthY  =  TMath::Abs(clusterArr.at(i).YMax() - clusterArr.at(i).YMin())/2.;

    for( int j = 0; j< clusterArr.size(); j++){
      if( clusterArr.at(j).RowID() < 5 ){ continue; }
      if( i == j ){ continue; }

      Int_t RowDelta = clusterArr.at(j).RowID() - clusterArr.at(i).RowID();
      if( RowDelta > 0 ){ continue; }
      RowDelta = TMath::Abs(RowDelta);
      TVector3 delta = clusterArr.at(j).GetClusterPos() - clusterArr.at(i).GetClusterPos();
      if( delta.Mag() > 70){ continue;}

      Double_t meanPhiSub = (clusterArr.at(j).PhiMax() + clusterArr.at(j).PhiMin())/2.;
      Double_t widthPhiSub= TMath::Abs(clusterArr.at(j).PhiMax() - clusterArr.at(j).PhiMin())/2.;
      Double_t meanYSub   =  clusterArr.at(j).GetClusterPos().Y();
      Double_t widthYSub  =  TMath::Abs(clusterArr.at(j).YMax() - clusterArr.at(j).YMin())/2.;

      if( TMath::Abs(meanPhiSub - meanPhi) > (widthPhi*(RowDelta))){ continue;}
      nRoot++;
    }
    if( nRoot == 0 ){
      RootCluster.push_back(clusterArr.at(i));
    }
  }
  return RootCluster;
}
std::vector<TPCPadHitCluster> TPCFindBlock( std::vector<TPCPadHitCluster> &clusterArr, int method ){

  std::vector<TPCPadHitCluster> tempCluster;
  if( clusterArr.size() < 3 ){ return tempCluster;} /// analyze events including more than 3 cluster
  Int_t maximumRow = -1;
  Int_t maxRowClusterID;

  /// find maximum row ID
  for( int icl = 0; icl < clusterArr.size(); icl++){
    if( clusterArr.at(icl).RowID() > maximumRow ){
      maximumRow = clusterArr.at(icl).RowID();
      maxRowClusterID = icl;
    }
  }

  /// start block finding
  tempCluster.push_back(clusterArr.at(maxRowClusterID));
  std::cout<< "starting "
	   << tempCluster.at(0).RowID() << "\t"
	   << tempCluster.at(0).PhiMin() << "\t"
	   << tempCluster.at(0).PhiMax() << std::endl;
  clusterArr.erase( clusterArr.begin() + maxRowClusterID );
  std::cout<< "Inner:0" << std::endl;

  //// find inner block
  for( int icl = 0; icl < clusterArr.size(); ){
    if( AdjacentCluster( clusterArr.at(icl), tempCluster.at(0),1)){
      tempCluster.push_back(clusterArr.at(icl));
      clusterArr.erase( clusterArr.begin() + icl );
    }else{
      icl++;
    }
  }
  std::cout<< tempCluster.size() << std::endl;
  if( tempCluster.size() != 2 ){
    std::cout << "Y division is required " << std::endl;
    return tempCluster;
  }
  std::cout<< "Inner:1" << std::endl;
  //find inner cluster
  for( int icl = 0; icl < clusterArr.size(); ){
    if( AdjacentCluster( clusterArr.at(icl), tempCluster.at(1),1)){
      tempCluster.push_back(clusterArr.at(icl));
      clusterArr.erase( clusterArr.begin() + icl );
    }else{
      icl++;
    }
  }
  if( tempCluster.size() != 3 ){
    std::cout << tempCluster.size() << std::endl;
    std::cout << "Y division is required " << std::endl;
    return tempCluster;
  }
  std::cout<< "Calculate : " << tempCluster.size() << std::endl;

  double length = 0;
  double rho    = 0;
  double dydl   = 0;
  double phi    = 0;
  double dphidl = 0;
  double ddydl  = 0;

  std::vector<double> vlength;
  std::vector<double> vrho;
  std::vector<double> vdydl;
  std::vector<double> vphi;

  std::vector<double> vddydl;
  std::vector<double> vdphidl;
  std::cout<< "ClusteringCal" << std::endl;
  ClusteringCalculator( tempCluster.at(0), tempCluster.at(1), length, rho, phi, dydl);
  vlength.push_back(length);
  vrho.push_back(rho);
  vdydl.push_back(dydl);
  vphi.push_back(phi);
  std::cout<< "ClusteringCal" << std::endl;
  ClusteringCalculator( tempCluster.at(1), tempCluster.at(2), length, rho, phi, dydl);

  vlength.push_back(length);
  vrho.push_back(rho);
  vdydl.push_back(dydl);
  vphi.push_back(phi);

  dphidl = (vphi.at(0) - vphi.at(1))/(vrho.at(0)+vrho.at(1))*2;
  ddydl  = (vdydl.at(0) - vdydl.at(1));
  vddydl.push_back(ddydl);
  vdphidl.push_back(dphidl);

  std::cout<< "EDGE" << std::endl;
  //Find next clusters

  Int_t  cTempID        = 3;
  Bool_t bfind_one_next = true;
  Bool_t bfind_two_next = true;
  Double_t th_phi = 2;
  Double_t th_dy  = 2;
  while( bfind_one_next || bfind_two_next){
    bfind_one_next = false;
    bfind_two_next = false;
    std::cout<< "Clustering" << std::endl;
    for( int icl = 0; icl < clusterArr.size(); ){
      //if( TMath::Abs(clusterArr.at(icl).RowID() - tempCluster.at(2).RowID()) > 1 ){ icl++; }
      if( !(AdjacentCluster(clusterArr.at(icl),tempCluster.at(cTempID-1), 1))){ icl++; }
      else{
	std::cout<< "passed" << std::endl;
	double length_sub;
	double rho_sub;
	double dydl_sub;
	double phi_sub;
	double dphidl_sub;
	double ddydl_sub;
	ClusteringCalculator( tempCluster.at(cTempID-1), clusterArr.at(icl), length_sub, rho_sub, phi_sub, dydl_sub);
	std::cout << "ClCal" << std::endl;
	dphidl_sub = (vphi.at(cTempID-2) - phi_sub)/(vrho.at(cTempID-2)+rho_sub)*2;
	ddydl_sub  = (vdydl.at(cTempID-2) - dydl_sub);
	std::cout << "ClCal" << std::endl;

	std::cout<< dphidl_sub <<"\t" << ddydl_sub << std::endl;
	if( dphidl_sub*TMath::RadToDeg() < th_phi && TMath::Abs(ddydl_sub) < th_dy ){
	  vlength.push_back( length_sub );
	  vrho.push_back( rho_sub );
	  vdydl.push_back( dydl_sub );
	  vphi.push_back( phi_sub );
	  vdphidl.push_back( dphidl_sub);
	  vddydl.push_back( ddydl_sub);
	  tempCluster.push_back(clusterArr.at(icl));
	  clusterArr.erase( clusterArr.begin() + icl);
	  cTempID++;
	  bfind_one_next = true;
	  break;
	}else{
	  icl++;
	}
      }
    }

    if( !bfind_one_next ){
      for( int icl = 0; icl < clusterArr.size(); ){
	//if( TMath::Abs(clusterArr.at(icl).RowID() - tempCluster.at(2).RowID()) > 1 ){ icl++; }
	if( !(AdjacentCluster(clusterArr.at(icl), tempCluster.at(cTempID-1), 2))){ icl++; }
	else{
	  double length_sub;
	  double rho_sub;
	  double dydl_sub;
	  double phi_sub;
	  double dphidl_sub;
	  double ddydl_sub;
	  ClusteringCalculator( tempCluster.at(cTempID-1), clusterArr.at(icl), length_sub, rho_sub, phi_sub, dydl_sub);
	  dphidl_sub = (vphi.at(cTempID-2) - phi_sub)/(vrho.at(cTempID-2)+rho_sub)*2;
	  ddydl_sub  = (vdydl.at(cTempID-2) - dydl_sub);

	  //dphidl_sub = (vphi.at(cTempID-2) -vphi.at(cTempID-1))/(vrho.at(cTempID-2)+vrho.at(cTempID-1))*2;
	  //ddydl_sub= (vdydl.at(cTempID-2) -vdydl.at(cTempID-1));
	  std::cout<< dphidl_sub <<"\t" << ddydl_sub << std::endl;
	  if( dphidl_sub*TMath::RadToDeg() < th_phi && TMath::Abs(ddydl_sub) < th_dy ){
	    vlength.push_back( length_sub );
	    vrho.push_back( rho_sub );
	    vdydl.push_back( dydl_sub );
	    vphi.push_back( phi_sub );
	    vdphidl.push_back( dphidl_sub);
	    vddydl.push_back( ddydl_sub);
	    tempCluster.push_back(clusterArr.at(icl));
	    clusterArr.erase( clusterArr.begin() + icl);
	    cTempID++;
	    bfind_two_next = true;
	    break;
	  }else{
	    icl++;
	  }
	}
      }
    }
  }
  std::cout<<"Clustered : " <<  tempCluster.size() << std::endl;
  return tempCluster;
}
std::vector<TPCPadHitCluster> TPCFindEdgeBlock( std::vector<TPCPadHitCluster> &clusterArr){
  /*
    std::vector<TPCPadHitCluster> tempCluster;
    if( clusterArr.size() < 3 ) { return tempCluster; }

    Int_t maximumRow = -1;
    Int_t maxRowClusterID = -1;

    for( int icl  = 0; icl < clusterArr.size(); icl++){
    if( clusterArr.at( icl).RowID() > maximumRow ){
    maximumRow = clusterArr.at( icl).RowID();
    maxRowClusterID = icl;
    }
    }
    if( maximumRow > 31 || maximumRow < 0 ){ return tempCluster;}
    if( maxRowClusterID < 0 || clusterArr.size() ){ return tempCluster; }

    tempCluster.push_back( clusterArr.at( maxRowClusterID ));
    clusterArr.erase( clusterArr.begin() + maxRowClusterID );

    /// Find inner connected clusters
    for(int icl = 0; icl < clusterArr.size();){
    if( AdjacentCluster( clusterArr.at( icl ), tempCluster.at(0), 1)){
    tempCluster.push_back( clusterArr.at(icl));
    clusterArr.erase( clusterArr.begin() + icl );
    }else{
    icl++;
    }
    }

    std::cout<< "EdgeCluster 1: " << tempCluster.size() << "\t" << clusterArr.size() << std::endl;

    for( int icl = 0; icl < clusterArr.size();){
    for( int jcl = 0; jcl < tempCluster.size() -1; jcl++ ){
    if( TMath::Abs(tempCluster.at(jcl).RowID() - tempCluster.at(0).RowID() ) > 1 ){ continue; }
    if( AdjacentCluster( clusterArr.at( icl ), tempCluster.at(jcl+1), 1)){
    tempCluster.push_back( clusterArr.at( icl ));
    clusterArr.erase(clusterArr.begin() + icl );
    }else{
    icl++;
    }
    }
    }
    std::cout<< "EdgeCluster 2: " << tempCluster.size() << "\t" << clusterArr.size() << std::endl;
    return tempCluster;
  */

  std::vector<TPCPadHitCluster> tempCluster;
  if( clusterArr.size() < 3 ){ return tempCluster;} /// analyze events including more than 3 cluster
  Int_t maximumRow = -1;
  Int_t maxRowClusterID;

  /// find maximum row ID
  for( int icl = 0; icl < clusterArr.size(); icl++){
    if( clusterArr.at(icl).RowID() > maximumRow ){
      maximumRow = clusterArr.at(icl).RowID();
      maxRowClusterID = icl;
    }
  }

  /// start block finding
  tempCluster.push_back(clusterArr.at(maxRowClusterID));
  std::cout<< "starting "
	   << tempCluster.at(0).RowID() << "\t"
	   << tempCluster.at(0).PhiMin() << "\t"
	   << tempCluster.at(0).PhiMax() << std::endl;
  clusterArr.erase( clusterArr.begin() + maxRowClusterID );
  std::cout<< "Inner:0" << std::endl;

  //// find inner block
  for( int icl = 0; icl < clusterArr.size(); ){
    if( AdjacentCluster( clusterArr.at(icl), tempCluster.at(0),1) &&
	TMath::Abs(clusterArr.at(icl).GetClusterPos().Y() - tempCluster.at(0).GetClusterPos().Y()) < 20 ){
      tempCluster.push_back(clusterArr.at(icl));
      clusterArr.erase( clusterArr.begin() + icl );
    }else{
      icl++;
    }
  }
  std::cout<< tempCluster.size() << std::endl;
  if( tempCluster.size() != 2 ){
    std::cout << "Y division is required " << std::endl;
    return tempCluster;
  }
  std::cout<< "Inner:1" << std::endl;
  //find inner cluster
  for( int icl = 0; icl < clusterArr.size(); ){
    if( AdjacentCluster( clusterArr.at(icl), tempCluster.at(1),1) &&
	TMath::Abs(clusterArr.at(icl).GetClusterPos().Y() - (tempCluster.at(1).GetClusterPos().Y() + (tempCluster.at(1).GetClusterPos().Y() - tempCluster.at(0).GetClusterPos().Y()))) < 5 ){
      tempCluster.push_back(clusterArr.at(icl));
      clusterArr.erase( clusterArr.begin() + icl );
    }else{
      icl++;
    }
  }
  if( tempCluster.size() != 3 ){
    std::cout << tempCluster.size() << std::endl;
    std::cout << "Y division is required " << std::endl;
    return tempCluster;
  }
  return tempCluster;

}
std::vector<TPCPadHitCluster> TPCFindBlockCircle( std::vector<TPCPadHitCluster> &clusterArr, int Change){
  std::vector<TPCPadHitCluster> tempCluster = TPCFindEdgeBlock( clusterArr );
  if( tempCluster.size() != 3 ){ return tempCluster; }
  std::cout<< "Starting Blocking" << std::endl;
  Int_t algorithmChange = Change;// clusterArr.size();
  TPCCircleFitter* fitter = new TPCCircleFitter();
  Int_t    nBlocking  = 0;
  Int_t    minimumRow = tempCluster.at(2).RowID();
  Double_t rad        = 0;
  Double_t avgDyDl    = 0;
  std::vector<Double_t> vtheta;
  std::vector<Double_t> vl;
  std::vector<Double_t> vy;
  std::vector<Double_t> vDist;

  std::vector<Double_t> vDtheta;
  std::vector<Double_t> vDl;
  std::vector<Double_t> vDy;
  std::vector<Double_t> xarr;
  std::vector<Double_t> zarr;
  std::vector<Double_t> yarr;

  std::cout<< "Starting Loop" << std::endl;
  Bool_t bFindNext = true;
  while( bFindNext ){
    bFindNext = false;
    std::cout<<"nBlocking: " << nBlocking << std::endl;
    vtheta.clear();
    vl.clear();
    vy.clear();
    vDist.clear();

    vDtheta.clear();
    vDl.clear();
    vDy.clear();
    xarr.clear();
    yarr.clear();
    zarr.clear();

    minimumRow = 32;
    for( int icl = 0; icl < tempCluster.size(); icl++ ){
      xarr.push_back( tempCluster.at(icl).GetClusterPos().X());
      yarr.push_back( tempCluster.at(icl).GetClusterPos().Y());
      zarr.push_back( tempCluster.at(icl).GetClusterPos().Z());
      if( tempCluster.at(icl).RowID() < minimumRow){ minimumRow = tempCluster.at(icl).RowID();}
    }
    std::cout<< "Fitting" << std::endl;
    TPCCircleFitResult circleResult;
    if( nBlocking < algorithmChange ){
      circleResult = fitter->ThreePointCircleFit( xarr,zarr );
    }else{
      circleResult = fitter->Fit(xarr,zarr);
    }
    std::cout<< "Parameter calculation" << std::endl;
    rad = circleResult.GetFitPar().X();
    avgDyDl = 0;
    for( int icl = 0; icl < tempCluster.size(); icl++){
      double theta = TMath::ATan2( xarr.at(icl) - circleResult.GetFitPar().Y(),
				   zarr.at(icl) - circleResult.GetFitPar().Z() );
      double l     = rad*theta;
      double y     = tempCluster.at(icl).GetClusterPos().Y();
      double dist  = circleResult.CalculateDist(tempCluster.at(icl).GetClusterPos().X(),
						tempCluster.at(icl).GetClusterPos().Z());
      vtheta.push_back( theta );
      vl.push_back(l);
      vy.push_back(y);
      vDist.push_back(TMath::Abs( dist) );
    }
    for( int icl = 0; icl < tempCluster.size() -1; icl++){
      vDtheta.push_back( vtheta.at(icl+1) - vtheta.at(icl));
      vDl.push_back( vl.at(icl+1) - vl.at(icl) );
      vDy.push_back( vy.at(icl+1) - vy.at(icl) );
      avgDyDl = avgDyDl + (vDy.at(icl) / vDl.at(icl));
    }
    //avgDyDl = avgDyDl/(tempCluster.size() - 1 );
    avgDyDl = (vy.at(tempCluster.size()-1) - vy.at(0))/(vl.at(tempCluster.size()-1) - vl.at(0));
    //////// find next cluster //////
    std::cout<< "Find next cluster" << std::endl;
    std::vector<double> rArr;
    std::vector<int>    idArr;
    /// Make candidate points
    for( int icl = 0; icl < clusterArr.size(); ){
      //if( TMath::Abs(clusterArr.at(icl).RowID() - minimumRow) > 4 ){
      if( minimumRow - clusterArr.at(icl).RowID() < 0 ||
	  minimumRow - clusterArr.at(icl).RowID() > 3 ){
	icl++;
	continue;
      }
      if( clusterArr.at(icl).RowID() < 5 ){
	icl++;
	continue;
      }
      /*
	if( !(AdjacentCluster( clusterArr.at(icl),tempCluster.at(tempCluster.size()-1), 1) ||
	AdjacentCluster( clusterArr.at(icl),tempCluster.at(tempCluster.size()-1), 2 )) ){
	icl++;
	continue;
	}
      */
      double chi2 = circleResult.CalculateChi2( clusterArr.at(icl).GetClusterPos().X(),
						clusterArr.at(icl).GetClusterPos().Z());
      double dist = circleResult.CalculateDist( clusterArr.at(icl).GetClusterPos().X(),
						clusterArr.at(icl).GetClusterPos().Z());
      double tempTheta = TMath::ATan2( clusterArr.at(icl).GetClusterPos().X() - circleResult.GetFitPar().Y(),
				       clusterArr.at(icl).GetClusterPos().Z() - circleResult.GetFitPar().Z());
      double tempDTheta =tempTheta - vtheta.at(vtheta.size() -1 );
      double tempDl     = rad*tempDTheta;
      double tempDy     = avgDyDl * tempDl;
      double DeltaY     = vy.at( vy.size() -1 ) + tempDy - clusterArr.at(icl).GetClusterPos().Y();
      double yCut = 10;
      double rCut = 15;

      if( nBlocking < algorithmChange ){
	if( TMath::Abs( DeltaY ) < yCut && dist < rCut && dist < rCut*2/TMath::Sqrt(tempCluster.size()) && TMath::Abs( DeltaY) < yCut*2/TMath::Sqrt(tempCluster.size())){//&& TMath::Abs(tempDl) < 100 ){
	  std::cout<< "DL:" << tempDl << std::endl;
	  bFindNext = true;
	  std::cout<< clusterArr.at(icl).RowID() << std::endl;
	  rArr.push_back(dist);
	  idArr.push_back(icl);
	  //tempCluster.push_back( clusterArr.at(icl) );
	  //clusterArr.erase( clusterArr.begin() + icl );
	  icl++;
	}else{
	  icl++;
	}
      }else{
	if( TMath::Abs( DeltaY ) < yCut && dist < rCut && dist < rCut*2/TMath::Sqrt(tempCluster.size()) && TMath::Abs( DeltaY) < yCut*2/TMath::Sqrt(tempCluster.size()) &&  chi2 < 3){//&& TMath::Abs(tempDl) < 100 ){// && chi2 < 5 ){
	  std::cout<< "DL:" << tempDl << std::endl;
	  bFindNext = true;
	  std::cout<< clusterArr.at(icl).RowID() << std::endl;
	  rArr.push_back(dist);
	  idArr.push_back(icl);
	  //tempCluster.push_back( clusterArr.at(icl) );
	  //clusterArr.erase( clusterArr.begin() + icl );
	  icl++;
	}else{
	  icl++;
	}
      }
    }
    /// Good point decision
    if( bFindNext ){
      Double_t minimumR = 1000;
      Int_t    candID;
      for( int icd = 0; icd < rArr.size(); icd++ ){
	if( minimumR > rArr.at(icd) ){
	  minimumR = rArr.at(icd);
	  candID = idArr.at(icd);
	}
      }
      tempCluster.push_back( clusterArr.at(candID));
      clusterArr.erase(clusterArr.begin() + candID);
    }
    nBlocking = nBlocking +1;
  }
  /// PrintOut clusterBlock
  for( int i = 0; i< tempCluster.size()-1; i++){
    std::cout<< i << "\t" << tempCluster.at(i).RowID() << "\t" << tempCluster.at(i).ColID() << "\t" << vl.at(i) << "\t" << vtheta.at(i) << "\t" << vy.at(i) << "\t" << vDist.at(i) << "\t" << vDy.at(i)/vDl.at(i) << std::endl;
  }

  return tempCluster;
}
Bool_t AdjacentCluster( TPCPadHitCluster c_0, TPCPadHitCluster c_1, int drow ){
  Double_t c0Min = c_0.PhiMin() - TPCGlobals::sTPC_Pad_Parameter[c_0.RowID()][3]/2.*TMath::DegToRad();
  Double_t c0Max = c_0.PhiMax() + TPCGlobals::sTPC_Pad_Parameter[c_0.RowID()][3]/2.*TMath::DegToRad();
  Double_t c1Min = c_1.PhiMin() - TPCGlobals::sTPC_Pad_Parameter[c_1.RowID()][3]/2.*TMath::DegToRad();
  Double_t c1Max = c_1.PhiMax() + TPCGlobals::sTPC_Pad_Parameter[c_1.RowID()][3]/2.*TMath::DegToRad();
  if( (c_1.RowID() - c_0.RowID()) == drow &&
      ( c0Min < c1Max && c0Max > c1Min )
      //((c0Min < c1Min && c0Max > c1Min ) ||
      //(c0Min >= c1Min && c0Min < c1Max ))// &&
      //((c0Min > c1Min && c0Min < c1Max ) ||
      //(c0Max > c1Min && c0Max < c1Max ) ||
      //(c0Min < c1Min && c0Max > c1Max )) &&
      //TMath::Abs( c_1.GetClusterPos().Y() - c_0.GetClusterPos().Y()) < 40 ){
      ){
    return kTRUE;
  }else{
    if( (c_1.RowID() - c_0.RowID()) == drow ){
      std::cout<<"Debug: " << drow << "\t"
	       << c_0.RowID() << "\t" << c_0.PhiMin() << "\t" <<  c_0.PhiMax() << "\t"
	       << c_1.RowID() << "\t" << c_1.PhiMin() << "\t" <<  c_1.PhiMax() << std::endl;
    }
    return kFALSE;
  }
}
void ClusteringCalculator( TPCPadHitCluster c_0, TPCPadHitCluster c_1, double& length, double& rho ,double& phi, double& dydl){
  double dl
    = TMath::Power( c_0.GetClusterPos().X() - c_1.GetClusterPos().X(), 2)
    + TMath::Power( c_0.GetClusterPos().Z() - c_1.GetClusterPos().Z(), 2);
  dl = TMath::Sqrt( dl );
  double   dy = c_0.GetClusterPos().Y() - c_1.GetClusterPos().Y();
  TVector3 dv = TVector3( c_0.GetClusterPos().Z() - c_1.GetClusterPos().Z(),
			  c_0.GetClusterPos().X() - c_1.GetClusterPos().X(),
			  c_0.GetClusterPos().Y() - c_1.GetClusterPos().Y());
  length = (c_0.GetClusterPos() - c_1.GetClusterPos()).Mag();
  rho    = dl;
  dydl   = dy/dl;
  phi = dv.Phi();
}
TGraph* DrawClusterArrXZ( std::vector<TPCPadHitCluster> clusterArr){
  TGraph* clusterEventXZ = new TGraph();
  for( int i = 0; i< clusterArr.size(); i++){
    clusterEventXZ->SetPoint(i, clusterArr.at( i ).GetClusterPos().Z(), clusterArr.at(i).GetClusterPos().X());
  }
  return clusterEventXZ;
}
TGraph* DrawClusterArrYZ( std::vector<TPCPadHitCluster> clusterArr){
  TGraph* clusterEventYZ = new TGraph();
  for( int i = 0; i< clusterArr.size(); i++){
    clusterEventYZ->SetPoint(i, clusterArr.at( i ).GetClusterPos().Z(), clusterArr.at(i).GetClusterPos().Y());
  }
  return clusterEventYZ;
}
bool CircleCrossing( TVector3 rxy0, TVector3 rxy1, TVector2& CrossingVec1, TVector2& CrossingVec2){

  Double_t r0 = rxy0.X();
  Double_t r1 = rxy1.X();
  TVector2 v0( rxy0.Y(), rxy0.Z());
  TVector2 v1( rxy1.Y(), rxy1.Z());

  Double_t dist = (v1 - v0).Mod();
  if( dist > r0 + r1 ||
      dist < TMath::Abs( r0 - r1)){
    CrossingVec1 = TVector2(0,0);
    CrossingVec2 = TVector2(0,0);
    return false;
  }

  Double_t phi  = (v1 - v0).Phi();
  TVector2 unit = (v1 - v0).Unit();

  double  a = (v1-v0).X();
  double  b = (v1-v0).Y();
  double  c = ( v1.Mod2() - v0.Mod2() + r0*r0 -r1*r1)/2.;

  double distPoint = TMath::Abs( a*v0.X() + b*v0.Y() -c )/TMath::Sqrt(a*a +b*b);
  double dTheta    = TMath::ACos( distPoint/r0);
  double Phi       = unit.Phi();

  TVector2 vecCent;
  vecCent.SetMagPhi(distPoint, Phi);
  TVector2 targetCenter( 0,-143);

  TVector2 v[2];
  v[0].SetMagPhi( r0, Phi - dTheta );
  v[1].SetMagPhi( r0, Phi + dTheta );
  //CrossingVec1 = v[0] + v0;
  //CrossingVec2 = v[1] + v0;
  v[0] = v[0] + v0;
  v[1] = v[1] + v0;
  if( (v[0] - targetCenter).Mod() < (v[1] - targetCenter).Mod() ){
    CrossingVec1 = v[0];
    CrossingVec2 = v[1];
  }else{
    CrossingVec2 = v[0];
    CrossingVec1 = v[1];
  }
  return true;
}
std::vector<TVector2> CircleTangent( TVector3 rxy0, TVector3 rxy1){

  std::vector<TVector2> tangentialVec;


  Double_t r0 = rxy0.X();
  Double_t r1 = rxy1.X();
  TVector2 v0( rxy0.Y(), rxy0.Z());
  TVector2 v1( rxy1.Y(), rxy1.Z());
  TVector2 CrossingVec1;
  TVector2 CrossingVec2;
  Double_t dist = (v1 - v0).Mod();
  if( dist > r0 + r1 ||
      dist < TMath::Abs( r0 - r1)){
    CrossingVec1 = TVector2(0,0);
    CrossingVec2 = TVector2(0,0);
    return tangentialVec;
  }

  Double_t phi  = (v1 - v0).Phi();
  TVector2 unit = (v1 - v0).Unit();

  TVector2 unit1 = (v0 - v1).Unit();
  Double_t phi1  = (v0 - v0).Phi();

  double  a = (v1-v0).X();
  double  b = (v1-v0).Y();
  double  c = ( v1.Mod2() - v0.Mod2() + r0*r0 -r1*r1)/2.;

  double distPoint = TMath::Abs( a*v0.X() + b*v0.Y() -c )/TMath::Sqrt(a*a +b*b);
  double dTheta    = TMath::ACos( distPoint/r0);
  double Phi       = unit.Phi();

  double distPoint1 = TMath::Abs( a*v1.X() + b*v1.Y() - c)/TMath::Sqrt(a*a +b*b);
  double dTheta1    = TMath::ACos( distPoint/r1);
  double Phi1       = unit1.Phi();

  TVector2 vecCent;
  vecCent.SetMagPhi(distPoint, Phi);
  TVector2 targetCenter( 0,-143);

  Double_t Momentum0 = r0*1.0/3.36; // MeV 1T
  Double_t Momentum1 = r1*1.0/3.36; // MeV 1T


  TVector2 v[4];
  v[0].SetMagPhi( r0, Phi - dTheta );
  v[1].SetMagPhi( r0, Phi + dTheta );
  v[2].SetMagPhi( r1, Phi1 - dTheta1 );
  v[3].SetMagPhi( r1, Phi1 + dTheta1 );

  TVector2 vTan[8];
  vTan[0].SetMagPhi( Momentum0, Phi - dTheta - TMath::Pi()/2.);
  vTan[1].SetMagPhi( Momentum0, Phi - dTheta + TMath::Pi()/2.);
  vTan[2].SetMagPhi( Momentum0, Phi + dTheta - TMath::Pi()/2.);
  vTan[3].SetMagPhi( Momentum0, Phi + dTheta + TMath::Pi()/2.);

  vTan[4].SetMagPhi( Momentum1, Phi1 - dTheta1 - TMath::Pi()/2.);
  vTan[5].SetMagPhi( Momentum1, Phi1 - dTheta1 + TMath::Pi()/2.);
  vTan[6].SetMagPhi( Momentum1, Phi1 + dTheta1 - TMath::Pi()/2.);
  vTan[7].SetMagPhi( Momentum1, Phi1 + dTheta1 + TMath::Pi()/2.);



  for( int i = 0; i< 8; i++){
    tangentialVec.push_back(vTan[i]);
  }
  ///// Calculate tangential vector /////
  return tangentialVec;
}
bool    yDistCalculator(TPCCircleFitResult rst0, TPCCircleFitResult rst1, Double_t& dist, TVector3& Pos){
  TVector2 crossingCandidate[2];
  bool result = CircleCrossing(rst0.GetFitPar(),rst1.GetFitPar(),crossingCandidate[0],crossingCandidate[1]);
  if( !result ){
    dist = 0;
    Pos.SetXYZ(0,0,0);
    return result;
  }

  double y00 = rst0.CalculateDeltaY(crossingCandidate[0]);
  double y01 = rst0.CalculateDeltaY(crossingCandidate[1]);
  double y10 = rst1.CalculateDeltaY(crossingCandidate[0]);
  double y11 = rst1.CalculateDeltaY(crossingCandidate[1]);
  std::cout<< TMath::Abs( y00 - y10 ) << "\t" << TMath::Abs( y01 - y11 ) << std::endl;
  Double_t yCenter = 0;
  if( TMath::Abs( y00 - y10 ) < TMath::Abs( y01 -y11 ) ){
    yCenter = (y00 + y10)/2.;
    dist = y00 - y10;
    Pos.SetXYZ(crossingCandidate[0].X(),yCenter,crossingCandidate[0].Y());
  }else{
    yCenter = (y01 + y11)/2.;
    dist = y01 - y11;
    Pos.SetXYZ(crossingCandidate[1].X(),yCenter,crossingCandidate[1].Y());
  }
  return true;
}
