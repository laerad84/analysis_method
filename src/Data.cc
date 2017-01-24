#include "Data.h"
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

ClassImp(HSprial)
HSprial::HSprial(): ID(-1),X(0),Y(0),Z(0),R(0),DTheta(0),DY(0),RL(0),InitPos(0,0,0),InitMom(0,0,0),FinalPos(0,0,0),FitMom(0,0,0),nCluster(0),TrackLength(0),DepE(0){
  ;
}
HSprial::HSprial(Int_t id, Double_t x, Double_t y, Double_t z, Double_t r, Double_t dY, Double_t dTheta, Int_t rl):
  ID(id),X(x),Y(y),Z(z),R(r),DTheta(dTheta),DY(dY),RL(rl),InitPos(0,0,0),InitMom(0,0,0),FinalPos(0,0,0),FitMom(0,0,0),nCluster(0),TrackLength(0),DepE(0) {
  ;
}
HSprial::HSprial(const HSprial& right){
  ID = right.ID;
  X  = right.X;
  Y  = right.Y;
  Z  = right.Z;
  R  = right.R;
  DTheta = right.DTheta;
  DY     = right.DY;
  RL     = right.RL;
  InitPos    = right.InitPos;
  InitMom    = right.InitMom;
  FinalPos   = right.FinalPos;
  FitMom     = right.FitMom;
  nCluster   = right.nCluster;
  TrackLength= right.TrackLength;
  DepE       = right.DepE;
}
HSprial::~HSprial(){
  ;
}
void HSprial::Print(){
  std::cout<< "Sprial ID : " << ID << "\n"
	   << "XYZR      : " << X << ", " << Y << ", " << Z << ", " << R << "\n"
	   << "DTheta    : " << DTheta << "\n"
	   << "DY        : " << DY << "\n"
	   << "RL        : " << RL << "\n"
	   << "nCluster  : " << nCluster << "\n"
	   << std::endl;
  InitPos.Print();
  FinalPos.Print();

}
TPolyLine3D* HSprial::GenerateHelix(){
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

TPolyLine* HSprial::GenerateArc(){
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

TF1* HSprial::GenerateXY(){
  TF1* func = new TF1("func","sin(x)");
  return func;
}
TF1* HSprial::GenerateZY(){
  TF1* func = new TF1("func","sin(x)");
  return func;
}

ClassImp(HLine)
HLine::HLine():ID(-1),SlopeX(0),SlopeY(0),SlopeZ(0),Offset(0){
  ;
}
HLine::HLine(Int_t id, Double_t slopeX,Double_t slopeY, Double_t slopeZ, Double_t offset):
  ID(id),SlopeX(slopeX),SlopeY(slopeY),SlopeZ(slopeX),Offset(offset)
{
  ;
}
HLine::HLine(const HLine& right){
  ID = right.ID;
  SlopeX = right.SlopeX;
  SlopeY = right.SlopeY;
  SlopeZ = right.SlopeZ;
  Offset = right.Offset;
}
HLine::~HLine(){
  ;
}
void HLine::Print(){
  std::cout<< "Line ID  : " << ID << "\n"
	   << "SlopeXYZ : " << SlopeX << ", " << SlopeY << ", " << SlopeZ << "\n"
	   << "Offset   : " << Offset
	   << std::endl;
}
ClassImp(TPCHit)
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
ClassImp(TempCluster)
TempCluster::TempCluster(){
  Init();
}
TempCluster::TempCluster( const TempCluster& right ){
  Init();
  ID = right.ID;
  nMother = right.nMother;
  nDaughter = right.nDaughter;

  for( int i = 0; i< nMother; i++){
    Mother[i] = right.Mother[i];
  }
  for( int i = 0; i< nDaughter;i++){
    Daughter[i] = right.Daughter[i];
  }

}
TempCluster::~TempCluster(){
  ;
}
void TempCluster::Init(){
  ID = -1;
  nMother = 0;
  nDaughter = 0;
  idList.clear();

  for( int i = 0; i< 16; i++){
    Mother[i] = -1;
    Daughter[i] = -1;
  }

}
