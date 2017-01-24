#include "TPCMCDataConverter.h"

ClassImp( TPCDataConverter )
TPCDataConverter::TPCDataConverter(TTree* tr){
	m_tr       = tr;
	m_det      = new GsimDetectorEventData();
	/*
  for( int i = 0; i< 32 ; i++){
	 m_digiDet[i] = new GsimDetectorEventData();
  }
  */
	m_particle = new GsimGenParticleData();

	m_poly        = new TPCPoly("TPCDataConverter_poly","TPCDataConverter_poly");
	m_RowYHistArr = new TPCRowYHistArray("TPCDataConverter_RowY",0);

	grHit     = new TGraph2D();
	grPadHit  = new TGraph2D();
	grCluster = new TGraph2D();
	for( int i = 0; i< 16; i++){
		grTrack[i]   = new TGraph2D();
	}

	grPadHitXZ  = new TGraph();
	grClusterXZ = new TGraph();
	grPadHitYZ  = new TGraph();
	grClusterYZ = new TGraph();

	m_nDet = 0;
	SetData();
}

TPCDataConverter::~TPCDataConverter(){
	;
}


void TPCDataConverter::SetData(){
	if( m_tr == NULL ){ std::cerr << "Error: No Tree" << std::endl;}
	m_tr->SetBranchAddress("GenParticle.",&m_particle);
	m_tr->SetBranchAddress("TPC.",&m_det);
}

TClonesArray* TPCDataConverter::GetDetDigi( int id ){
	if( id < 0 || id >= m_nDet ){ return NULL; }
	return m_digiDet[id]->digi;
}
Int_t TPCDataConverter::AddDetector( char* detName ){
	m_digiDet[m_nDet] = new GsimDetectorEventData();
	Int_t rtn = m_tr->SetBranchAddress(detName,&m_digiDet[m_nDet]);
	if( rtn != 0 ){
		std::cout<< "Wrong detector Name : " << detName  << std::endl;
		delete m_digiDet[m_nDet];
		return m_nDet;
	}else{
		std::cout << "Add Detector : " << detName << std::endl;
		m_nDet++;
		return m_nDet;


	}
}


GsimDetectorEventData* TPCDataConverter::GetDetectorEventData(){
	return m_det;
}

//std::vector<TPCPadHit> TPCDataConverter::Convert( int ievent ){
void TPCDataConverter::Convert( int ievent ){
	m_tr->GetEntry( ievent );
	m_TrackArr = m_particle->briefTracks;
	m_HitArr   = m_det->hits;
	Double_t depEArr[TPCGlobals::sTPC_NMaximum];
	Double_t depYArr[TPCGlobals::sTPC_NMaximum];


	//// initialize data handlers
	grHit->Set(0);
	grPadHit->Set(0);
	grPadHitXZ->Set(0);
	grPadHitYZ->Set(0);

	m_PadHitArr.clear();
	m_ClusterArr.clear();

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
	//std::cout<< trackList.size() << "\t" << hitTrackIDList.size() << std::endl;
	//// generate average hit point by pad // TPCPadHit
	//std::vector<TPCPadHit> m_PadHitArr;
	Int_t iTrack = 0;
	for( int iarr = 0; iarr < hitTrackIDList.size(); iarr++){
		//m_poly->Reset();
		m_RowYHistArr->Reset();
		// make TPCPadHit array by track //
		//std::cout<< "digitize : " << trackList.at(iarr) << std::endl;
		Int_t nFill = 0;
		for( int jarr = 0; jarr < hitTrackIDList.at(iarr).size(); jarr++){
			GsimDetectorHitData* hit = (GsimDetectorHitData*)m_HitArr->At(hitTrackIDList.at(iarr).at(jarr));
			//Double_t diffusionR      = 0.18*TMath::Sqrt((hit->r.Y()+300)/10);/// diffusion radius by R, drift
			Double_t diffusionR      = 15*TMath::Sqrt((hit->r.Y()+300)/150);/// diffusion radius by R, drift
			Double_t diffusionG      = 0.1; // diffusion radius by GEM, ignored in the current version
			Int_t    nElectrons      = (int)(hit->edep*TPCGlobals::sTPC_NEperEnergy*20);
			grHit->SetPoint( grHit->GetN(), hit->r.X(), hit->r.Y(), hit->r.Z());

			/*
			 if( hit->track == trackList.at(iarr)){
			 std::cout<< hit->track << " : " << nElectrons << std::endl;
			 }
			 */
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
				/*
				 if( colID < 0 || colID >= TPCGlobals::sTPC_Pad_Parameter[rowID][1]){
				 continue;
				 }
				 */
				//std::cout<< xhit <<"\t" << yhit << "\t" << zhit << std::endl;
				m_RowYHistArr->hisRowY[rowID]->Fill( colID, yhit );
			}

			int rid, cid;
			gTPCIDHandler->GetRowCol(hit->r.X(),hit->r.Z(), rid, cid );
			//std::cout<< "RCID " << rid << "\t" << cid  << std::endl;
			nFill++;
		}

		//std::cout<< "nFill " << nFill << std::endl;
		/// Get PadHit information from hists
		for( int irow = 0; irow < 32; irow++){
			TProfile* prof = m_RowYHistArr->hisRowY[irow]->ProfileX();
			TH1D*     prjt = m_RowYHistArr->hisRowY[irow]->ProjectionX();
			for( int icol = 0; icol < TPCGlobals::sTPC_Pad_Parameter[irow][1]; icol++ ){
				double ymean = prof->GetBinContent(icol+1);
				//std::cout<< ymean << std::endl;
				//if( ymean < 2 ){ continue; }
				//std::cout <<trackList.at(iarr) << "\t" <<  icol << "\t" << irow << std::endl;
				Int_t  padID  = gTPCIDHandler->GetID( irow, icol );
				double   yerr = prof->GetBinError( icol+1 );
				double   ene  = prjt->GetBinContent( icol+1 );
				if( ene < 10 ){ continue; }
				TPCPadHit padHit;
				TVector2 TPCXZ = gTPCIDHandler->GetXZ( irow, icol );
				TVector3 vec(TPCXZ.X(), ymean, TPCXZ.Y());
				padHit.SetTrackID( iTrack);
				padHit.SetMCTrackID( trackList.at(iarr));
				padHit.SetPosition(vec);
				padHit.SetRowCol( irow, icol );
				padHit.SetEnergy( ene );//number of hit electron

				m_PadHitArr.push_back(padHit);
				//std::cout<< m_PadHitArr.at(m_PadHitArr.size() -1 ).Col() << std::endl;
				//std::cout<<"Arr Size: " <<  m_PadHitArr.size() << std::endl;
				grPadHit->SetPoint(grPadHit->GetN(), padHit.Position().X(), padHit.Position().Y(), padHit.Position().Z());
				grPadHitXZ->SetPoint( grPadHitXZ->GetN(), padHit.Position().Z(), padHit.Position().X());
				grPadHitYZ->SetPoint( grPadHitYZ->GetN(), padHit.Position().Z(), padHit.Position().Y());
				iTrack++;
			}
			prof->Delete();
			prjt->Delete();
		}
	}
	//std::cout<< "HitArrSize : " << m_PadHitArr.size() << std::endl;
	//return m_PadHitArr;
}
void TPCDataConverter::Clustering(){
	m_ClusterArr.clear();
	grCluster->Set(0);
	grClusterXZ->Set(0);
	grClusterYZ->Set(0);
	std::vector<TPCPadHit> HitArr = m_PadHitArr;
	Int_t clusterID = 0;
	for( int irow = 0; irow < 32; irow++){
		while( m_PadHitArr.size() != 0 ){
			TPCPadHitCluster cluster  = TPCClusterer( irow, HitArr );
			if( cluster.HitArr().size() != 0 ){
				cluster.SetClusterID(clusterID);
				m_ClusterArr.push_back( cluster );
				grCluster->SetPoint( grCluster->GetN(),
									cluster.GetClusterPos().X(),
									cluster.GetClusterPos().Y(),
									cluster.GetClusterPos().Z());
				grClusterXZ->SetPoint( grClusterXZ->GetN(),
									  cluster.GetClusterPos().Z(),
									  cluster.GetClusterPos().X());
				grClusterYZ->SetPoint( grClusterYZ->GetN(),
									  cluster.GetClusterPos().Z(),
									  cluster.GetClusterPos().Y());
			}else{ break; }
		}
	}
}

void TPCDataConverter::ClearSubData(){
	for( int i = 0; i< 32; i++){

	}
}
