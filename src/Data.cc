#include "Data.h"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class HCircle
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(HCircle)
HCircle::HCircle(): ID(-1), X(0),Z(0),R(0){
  ;
}
HCircle::HCircle(Int_t id, Double_t x, Double_t z, Double_t r): ID(id),X(x),Z(z),R(r){
  ;
}
HCircle::HCircle(const HCircle& right){
  ID = right.ID;
  X  = right.X;
  Z  = right.Z;
  R  = right.R;
}

HCircle::~HCircle(){
  ;
}

void HCircle::Print(){
  std::cout<< "ID  : " << ID << "\n"
	   << "XZR : " << X << ", " << Z << ", " << R << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class HHelix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(HHelix)
HHelix::HHelix(): ID(-1),X(0),Y(0),Z(0),R(0),DTheta(0),DY(0),RL(0),InitPos(0,0,0),InitMom(0,0,0),FinalPos(0,0,0),FitMom(0,0,0),nCluster(0),TrackLength(0),DepE(0){
  EX=0;
  EZ=0;
  ER=0;
}
HHelix::HHelix(Int_t id, Double_t x, Double_t y, Double_t z, Double_t r, Double_t dY, Double_t dTheta, Int_t rl):
  ID(id),X(x),Y(y),Z(z),R(r),DTheta(dTheta),DY(dY),RL(rl),InitPos(0,0,0),InitMom(0,0,0),FinalPos(0,0,0),FitMom(0,0,0),nCluster(0),TrackLength(0),DepE(0) {
  EX=0;
  EZ=0;
  ER=0;
  EY=0;
}
HHelix::HHelix(const HHelix& right){
  ID = right.ID;
  X  = right.X;
  Y  = right.Y;
  Z  = right.Z;
  R  = right.R;
  EX = right.EX;
  EZ = right.EZ;
  EY = right.EY;
  ER = right.ER;
  DTheta = right.DTheta;
  DY     = right.DY;
  EDYDTheta= right.EDYDTheta;
  RL     = right.RL;

  InitPos    = right.InitPos;
  InitMom    = right.InitMom;
  FinalPos   = right.FinalPos;
  FitMom     = right.FitMom;
  nCluster   = right.nCluster;
  TrackLength= right.TrackLength;
  DepE       = right.DepE;
}
HHelix::~HHelix(){
  ;
}
void HHelix::ResetAll(){
  ID=-1;
  X=0;
  Y=0;
  Z=0;
  R=0;
  DTheta      = 0;
  DY          = 0;
  RL          = 0;
  InitPos     = TVector3(0,0,0);
  InitMom     = TVector3(0,0,0);
  FinalPos    = TVector3(0,0,0);
  FitMom      = TVector3(0,0,0);
  nCluster    = 0;
  TrackLength = 0;
  DepE = 0;
  EX   = 0;
  EZ   = 0;
  EY   = 0;
  ER   = 0;
  EDYDTheta = 0;
}
void HHelix::Print(){
  std::cout<< "Helix ID  : " << ID << "\n"
	   << "XYZR      : " << X << ", " << Y  << ", " << Z  << ", " << R << "\n"
	   << "EXZR      : " << EX <<", " << EZ << ", " << ER << "\n"
	   << "DTheta    : " << DTheta << "\n"
	   << "DY        : " << DY << "\n"
	   << "RL        : " << RL << "\n"
	   << "nCluster  : " << nCluster << "\n"
	   << std::endl;
  InitPos.Print();
  FinalPos.Print();

}
TPolyLine3D* HHelix::GenerateHelix(){
  TPolyLine3D* tmpLine = new TPolyLine3D();
  Double_t c         = this->DY/( this->DTheta );
  Double_t thetaOff  = TMath::ATan2( InitPos.X() - X, InitPos.Z() - Z);
  Double_t lastTheta = TMath::ATan2( FinalPos.X() - X, FinalPos.Z() -Z );
  Double_t midTheta  = ( thetaOff + lastTheta )/2.;
  Double_t deltaTheta = this->DTheta;
  //std::cout<< DY << "\t" << R << "\t" << DTheta << "\t" << c <<  std::endl;
  /*
  for( int i = -180; i <= 180; i++){
    Double_t tmpX = this->R * TMath::Sin((double)(i)*TMath::DegToRad() + thetaOff ) + X;
    Double_t tmpZ = this->R * TMath::Cos((double)(i)*TMath::DegToRad() + thetaOff ) + Z;
    Double_t tmpY = this->DY/( this->DTheta )*((double)(i)*TMath::DegToRad()) + InitPos.Y();
    tmpLine->SetNextPoint( tmpX, tmpY, tmpZ );
  }
  */
  for( int i = -50; i <=150; i++){
    Double_t theta = ((double)(i))*deltaTheta/100.;
    Double_t tmpX = this->R * TMath::Sin(theta + thetaOff ) + X;
    Double_t tmpZ = this->R * TMath::Cos(theta + thetaOff ) + Z;
    Double_t tmpY = this->DY/( this->DTheta )*(theta) + InitPos.Y();
    if( TMath::Abs(tmpX) > 300 ||
	TMath::Abs(tmpY) > 300 ||
	TMath::Abs(tmpZ) > 300 ){
      continue;
    }
    tmpLine->SetNextPoint( tmpX, tmpY, tmpZ );
  }
  return tmpLine;
}
void HHelix::CalculateDist( TVector3 hitpos, Double_t& dR, Double_t& dY ){
  Double_t theta = TMath::ATan2(hitpos.X() - X, hitpos.Z() - Z );
  Double_t thetaOff= TMath::ATan2( InitPos.X() - X, InitPos.Z() - Z);
  Double_t tmpX  = R * TMath::Sin( theta ) + X;
  Double_t tmpZ  = R * TMath::Cos( theta ) + Z;
  Double_t tmpY  = DY/DTheta*(theta-thetaOff) + InitPos.Y();
  Double_t tmpR  = TMath::Sqrt( (hitpos.X()-X)*(hitpos.X()-X) + (hitpos.Z()-Z)*(hitpos.Z()-Z));
  Double_t pitch = 2*DY/DTheta*TMath::Pi();
  dY = hitpos.Y() - tmpY;
  if( TMath::Abs(dY) > pitch/2 ){
    if( dY > 0 ){ dY= pitch -dY;}
    else{ dY = dY + pitch;}
  }
  dR = tmpR - R;
  //std::cout<< dY << "\t" << dR << std::endl;
}
void HHelix::CalculateDist( TGraph2D* grTrack, Double_t& dR, Double_t& dY ){
  Int_t nPoint=0;
  Double_t RSum=0;
  Double_t YSum=0;
  for( int i = 0; i< grTrack->GetN(); i++){
    TVector3 pos(grTrack->GetX()[i], grTrack->GetY()[i],grTrack->GetZ()[i]);
    Double_t dy,dr;
    CalculateDist( pos, dr, dy );
    RSum += TMath::Abs(dr);
    YSum += TMath::Abs(dy);
    nPoint++;
  }
  dR = RSum / nPoint;
  dY = YSum / nPoint;
}
TPolyLine* HHelix::GenerateArc(){
  TPolyLine* tmpLine = new TPolyLine();

  Double_t c         = this->DY/( this->DTheta );
  Double_t thetaOff  = TMath::ATan2( InitPos.X() - X, InitPos.Z() - Z);
  Double_t lastTheta = TMath::ATan2( FinalPos.X() - X, FinalPos.Z() -Z );
  Double_t midTheta  = ( thetaOff + lastTheta )/2.;
  Double_t deltaTheta = this->DTheta;

  for( int i = -50; i <=150; i++){
    Double_t theta = ((double)(i))*deltaTheta/100.;
    Double_t tmpX = this->R * TMath::Sin(theta + thetaOff ) + X;
    Double_t tmpZ = this->R * TMath::Cos(theta + thetaOff ) + Z;
    Double_t tmpY = this->DY/( this->DTheta )*(theta) + InitPos.Y();
    if( TMath::Abs(tmpX) > 300 ||
	TMath::Abs(tmpY) > 300 ||
	TMath::Abs(tmpZ) > 300 ){
      continue;
    }
    tmpLine->SetNextPoint( tmpZ, tmpX );
  }
  return tmpLine;
}
TF1* HHelix::GenerateXY(){
  Double_t thetaOffset = TMath::ATan2( InitPos.X() -X ,InitPos.Z() - Z);
  Double_t c = this->DY/this->DTheta;
  Double_t yOff = InitPos.Y();
  TF1* func = new TF1("func","[0]*sin([1]*(x-[2])+[3])+[4]");
  func->SetParameter(0,this->R);
  func->SetParameter(1,1./c);
  func->SetParameter(2,yOff);
  func->SetParameter(3,thetaOffset);
  func->SetParameter(4,X);
  func->SetRange(-300,300);
  return func;
}
TF1* HHelix::GenerateZY(){
  Double_t thetaOffset = TMath::ATan2( InitPos.X() -X ,InitPos.Z() - Z);
  Double_t c = this->DY/this->DTheta;
  Double_t yOff = InitPos.Y();
  TF1* func = new TF1("func","[0]*cos([1]*(x-[2])+[3])+[4]");
  func->SetParameter(0,this->R);
  func->SetParameter(1,1./c);
  func->SetParameter(2,yOff);
  func->SetParameter(3,thetaOffset);
  func->SetParameter(4,Z);
  func->SetRange(-300,300);
  return func;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class HTrack
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(HTrack)
HTrack::HTrack():ID(-1),R(0),DYDL(0),DYDT(0),Length(0),DepE(0){
  Momentum = TVector3(0,0,0);
}
HTrack::HTrack(const HTrack& right){
  ID = right.ID;
  R  = right.R;
  DYDL=right.DYDL;
  DYDT=right.DYDT;
  Length=right.Length;
  DepE=right.DepE;
  Momentum = right.Momentum;
}
HTrack::~HTrack(){
  ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class HLine
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(HLine)
HLine::HLine():ID(-1){
  Direction    = TVector3(0,0,0);
  ClosestPoint = TVector3(0,0,0);
  Offset       = TVector3(0,0,0);
  bNorm        = kFALSE;
}
HLine::HLine(Int_t id, TVector3 slope, TVector3 offset){
  ID        = id;
  Direction = slope;
  Offset    = offset;
  bNorm     = false;
}

HLine::HLine(const HLine& right){
  ID        = right.ID;
  Direction = right.Direction;
  Offset =   right.Offset;
  bNorm     = bNorm;
}

HLine::~HLine(){
  ;
}
void HLine::Print(){
  std::cout<< "Line ID  : " << ID << "\n"
	   << std::endl;
  std::cout<< "Direction " << std::endl;
  Direction.Print();
  std::cout<< "Offset" << std::endl;
  Offset.Print();
  std::cout<< "ClosestPoint" << std::endl;
  ClosestPoint.Print();
}
void HLine::Normalize(){
  Double_t NormDirection = Direction.Mag();
  Direction = Direction.Unit();
  bNorm = kTRUE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class TPCHit
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
  TPCHit::TPCHit():ID(-1),Row(-1),Col(-1),Y(0),Energy(0){
  ;
  }
  TPCHit::TPCHit( int id, Double_t y, Double_t energy): ID( id ),Row(-1),Col(-1),Y(y),Energy(energy){
  ;
  }
  TPCHit::TPCHit(const TPCHit& right ){
  ID     = right.ID;
  Row    = right.Row;
  Col    = right.Col;
  Y      = right.Y;
  Energy = right.Energy;
  }
  TPCHit::~TPCHit(){
  ;
  }
*/

ClassImp( TPCHit )

TPCHit::TPCHit(){
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
TPCHit::TPCHit(const TPCHit& right):
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
TPCHit::TPCHit(int padID, TVector3 pos, double energy ) :
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
TPCHit::~TPCHit(){
  Clear();
}
void TPCHit::Clear(Option_t* opt){

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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class TPCCluster
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

ClassImp(TPCCluster)
TPCCluster::TPCCluster():ID(-1),Energy(0),Row(-1),Col(-1),ColMin(0),ColMax(0),Position(0,0,0),YMin(0),YMax(0),NHit(0),Mother(-1),Daughter(-1),nMother(0),nDaughter(0),MinMotherDist(-1),MinDaughterDist(-1),MotherY(-999),DaughterY(-999),Blocked(false){
  MotherID.clear();
  DaughterID.clear();
}
TPCCluster::TPCCluster(int id, Double_t y, Double_t energy):ID(id),Energy(energy),Row(-1),Col(-1),ColMin(0),ColMax(0),Position(0,y,0),YMin(0),YMax(0),NHit(0),Mother(-1),Daughter(-1),nMother(0),nDaughter(0),MinMotherDist(-1),MinDaughterDist(-1),MotherY(-999),DaughterY(-999),Blocked(false){
  MotherID.clear();
  DaughterID.clear();
}
TPCCluster::TPCCluster(const TPCCluster& right){
  ID     = right.ID;
  Energy = right.Energy;
  Row    = right.Row;
  Col    = right.Col;
  ColMin = right.ColMin;
  ColMax = right.ColMax;
  Position   = right.Position;
  YMin       = right.YMin;
  YMax       = right.YMax;
  NHit       = right.NHit;
  PhiMax     = right.PhiMax;
  PhiMin     = right.PhiMin;
  Mother     = right.Mother;
  Daughter   = right.Daughter;
  MinMotherDist   = right.MinMotherDist;
  MinDaughterDist = right.MinDaughterDist;
  MotherY    = right.MotherY;
  DaughterY  = right.DaughterY;

  nMother    = right.nMother;
  nDaughter  = right.nDaughter;
  //MotherID   = right.MotherID;
  //DaughterID = right.DaughterID;
  Blocked    = right.Blocked;

  MotherID.clear();
  DaughterID.clear();

  for( int isize = 0; isize < right.MotherID.size(); isize++){
    MotherID.push_back( right.MotherID.at(isize));
  }
  for( int isize = 0; isize < right.DaughterID.size(); isize++){
    DaughterID.push_back( right.DaughterID.at(isize));
  }
  //Print();
}
TPCCluster::~TPCCluster(){
  ;
}
void TPCCluster::Print(){
  std::cout<< "Print" << std::endl;
  std::cout<< ID << "\t"
	   << Mother << "\t" << Daughter << "\t"
	   << nMother << "\t" << nDaughter << "\t"
	   << MotherID.size() << "\t"   << DaughterID.size() << "\t"
	   << Energy << "\t" << NHit<< "\t"
	   << PhiMin << "\t" << PhiMax << "\t"
	   << MinMotherDist << "\t" << MinDaughterDist << "\t"
	   << MotherY << "\t" << DaughterY << "\t"
	   << std::endl;
  std::cout<< "Components : " << std::endl;
  std::cout<< "Mother  : " << "\t";
  for( int i = 0; i< MotherID.size(); i++){
    std::cout<< MotherID.at(i) << "\t";
  }
  std::cout<< std::endl;
  std::cout<< "Daughter : " << "\t";
  for( int j = 0; j << DaughterID.size();j++){
    std::cout<< DaughterID.at(j) << "\t";
  }
  std::cout<< std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class TPCTrack
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClassImp(TPCTrack)
TPCTrack::TPCTrack(){
  ID           = -1;
  PID          = 0;
  InitialRowID = -1;
  FinalRowID   = -1;
  nCluster     = 0;
  track        = new TGraph2D();
  trackErr     = new TGraph2D();
  EArr.clear();
  Momentum     = TVector3(0,0,0);
  Length       = 0;
  DepE         = 0;
  ChisqY       = 0;
  ChisqR       = 0;
}

TPCTrack::TPCTrack( const TPCTrack& right){
  ID           = right.ID;
  PID          = right.PID;
  InitialRowID = right.InitialRowID;
  FinalRowID   = right.FinalRowID;
  track        = right.track;
  trackErr     = right.trackErr;
  EArr.clear();
  for( int i = 0; i< right.EArr.size(); i++){
    EArr.push_back( right.EArr.at(i));
  }
  Momentum = right.Momentum;
  Length   = right.Length;
  DepE     = right.DepE;
  ChisqR   = right.ChisqR;
  ChisqY   = right.ChisqY;
  helix    = right.helix;
}

TPCTrack::~TPCTrack(){
  ;
}

void TPCTrack::Init(){
  ID           = -1;
  PID          = 0;
  InitialRowID = -1;
  FinalRowID   = -1;
  nCluster     = 0;
  track->Set(0);
  trackErr->Set(0);
  EArr.clear();

  Length = 0;
  DepE = 0;
  ChisqR = 0;
  ChisqY = 0;
  Momentum.SetXYZ(0,0,0);
  helix.ResetAll();
}

void TPCTrack::Print(){
  std::cout << "Track ID : " << ID << "\n"
	    << "RowID    : " << InitialRowID << "\t" << FinalRowID << "\n"
	    << "PID      : " << PID    << "\n"
	    << "Length   : " << Length << "\t"
	    << "DepositE : " << DepE   << "\n";
  helix.Print();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class HPiM
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClassImp( HPiM )
HPiM::HPiM(){
  ID          = -1;
  RecMass     = 0;
  Energy      = 0;
  Momentum    = TLorentzVector(0,0,0,0);
  Position    = TVector3(0,0,0);
}

HPiM::HPiM( const HPiM& right ){
  ID        = right.ID;
  RecMass   = right.RecMass;
  Energy    = right.Energy;
  Momentum  = right.Momentum;
  Position  = right.Position;
}

HPiM::~HPiM(){
  ;
}

void HPiM::Init(){
  ID = -1;
  RecMass = 0;
  Energy      = 0;
  Momentum = TLorentzVector(0,0,0,0);
  Position = TVector3(0,0,0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class HKaonP
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClassImp( HKaonP )
HKaonP::HKaonP(){
  ID          = -1;
  RecMass     = 0;
  Energy      = 0;
  Momentum    = TLorentzVector(0,0,0,0);
  Position    = TVector3(0,0,0);
}

HKaonP::HKaonP( const HKaonP& right ){
  ID        = right.ID;
  RecMass   = right.RecMass;
  Energy    = right.Energy;
  Momentum  = right.Momentum;
  Position  = right.Position;
}

HKaonP::~HKaonP(){
  ;
}

void HKaonP::Init(){
  ID = -1;
  RecMass = 0;
  Energy      = 0;
  Momentum = TLorentzVector(0,0,0,0);
  Position = TVector3(0,0,0);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class HProton
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClassImp( HProton )
HProton::HProton(){
  ID          = -1;
  RecMass     = 0;
  Energy      = 0;
  Momentum    = TLorentzVector(0,0,0,0);
  Position    = TVector3(0,0,0);
}

HProton::HProton( const HProton& right ){
  ID        = right.ID;
  RecMass   = right.RecMass;
  Energy    = right.Energy;
  Momentum  = right.Momentum;
  Position  = right.Position;
}

HProton::~HProton(){
  ;
}

void HProton::Init(){
  ID = -1;
  RecMass = 0;
  Energy      = 0;
  Momentum = TLorentzVector(0,0,0,0);
  Position = TVector3(0,0,0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class HLambda
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClassImp( HLambda )
HLambda::HLambda(){
  ID=-1;
  RecMass     = 0;
  Energy      = 0;
  Momentum    = TLorentzVector(0,0,0,0);
  Position    = TVector3(0,0,0);
  proton      = new HProton();
  pion        = new HPiM();
}
HLambda::HLambda( const HLambda& right ){
  ID        = right.ID;
  RecMass   = right.RecMass;
  Energy    = right.Energy;
  Momentum  = right.Momentum;
  Position  = right.Position;
  proton    = right.proton;
  pion      = right.pion;
}
HLambda::~HLambda(){
  ;
}
void HLambda::Init(){
  ID = -1;
  RecMass = 0;
  Energy      = 0;
  Momentum = TLorentzVector(0,0,0,0);
  Position = TVector3(0,0,0);
  proton->Init();
  pion->Init();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class HCascade
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClassImp( HCascade )
HCascade::HCascade(){
  ID=-1;
  RecMass     = 0;
  Energy      = 0;
  Momentum    = TLorentzVector(0,0,0,0);
  Position    = TVector3(0,0,0);
  lambda      = new HLambda();
  pion        = new HPiM();
}

HCascade::HCascade( const HCascade& right ){
  ID        = right.ID;
  RecMass   = right.RecMass;
  Energy    = right.Energy;
  Momentum  = right.Momentum;
  Position  = right.Position;
  lambda    = right.lambda;
  pion      = right.pion;
}
HCascade::~HCascade(){
  ;
}
void HCascade::Init(){
  ID = -1;
  RecMass = 0;
  Energy      = 0;
  Momentum = TLorentzVector(0,0,0,0);
  Position = TVector3(0,0,0);
  lambda->Init();
  pion->Init();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class HDibaryonLL
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClassImp( HDibaryonLL )
HDibaryonLL::HDibaryonLL(){
  ID=-1;
  RecMass     = 0;
  Momentum    = TLorentzVector(0,0,0,0);
  Position    = TVector3(0,0,0);
  lambda0     = new HLambda();
  lambda1     = new HLambda();
}
HDibaryonLL::HDibaryonLL( const HDibaryonLL& right ){
  ID        = right.ID;
  RecMass   = right.RecMass;
  Energy    = right.Energy;
  Momentum  = right.Momentum;
  Position  = right.Position;
  lambda0   = right.lambda0;
  lambda1   = right.lambda1;
}
HDibaryonLL::~HDibaryonLL(){
  ;
}
void HDibaryonLL::Init(){
  ID = -1;
  RecMass = 0;
  Energy      = 0;
  Momentum = TLorentzVector(0,0,0,0);
  Position = TVector3(0,0,0);
  lambda0->Init();
  lambda1->Init();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Class HDibaryonPC
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClassImp( HDibaryonPC )
HDibaryonPC::HDibaryonPC(){
  ID=-1;
  RecMass     = 0;
  Energy      = 0;
  Momentum    = TLorentzVector(0,0,0,0);
  Position    = TVector3(0,0,0);
  cascade     = new HCascade();
  pion        = new HPiM();
}
HDibaryonPC::HDibaryonPC( const HDibaryonPC& right ){
  ID        = right.ID;
  RecMass   = right.RecMass;
  Energy    = right.Energy;
  Momentum  = right.Momentum;
  Position  = right.Position;
  cascade   = right.cascade;
  pion      = right.pion;
}
HDibaryonPC::~HDibaryonPC(){
  ;
}
void HDibaryonPC::Init(){
  ID = -1;
  RecMass = 0;
  Energy      = 0;
  Momentum = TLorentzVector(0,0,0,0);
  Position = TVector3(0,0,0);
  cascade->Init();
  pion->Init();
}
