#include "TPCIDManager.h"
#include <iostream>

ClassImp(TPCIDHandler)

TPCIDHandler* gTPCIDHandler = new TPCIDHandler();

TPCIDHandler::TPCIDHandler(double offset){
  m_tpcOffset = offset;
}

TPCIDHandler::~TPCIDHandler(){
  ;
}

int TPCIDHandler::GetRow( int ID ){
  int tRow = 0;
  for( int i = 0; i< 32; i++){
    if( ID < TPCGlobals::sTPC_PadRowMaxID[i] ){
      tRow = i;
      break;
    }
  }
  return tRow;
}

int TPCIDHandler::GetCol( int ID ){
  int tRow = this->GetRow( ID );
  int tCol = 0;
  if( tRow  == 0 ){
    tCol = ID;
  }else{
    tCol = ID - TPCGlobals::sTPC_PadRowMaxID[tRow -1];
  }
  return tCol;
}

TVector2 TPCIDHandler::GetXZ(int row, int col){
  Double_t rad   = TPCGlobals::sTPC_Pad_Parameter[row][2];
  Double_t sTheta= -180+TPCGlobals::sTPC_Pad_Parameter[row][4];
  Double_t theta = (-180+TPCGlobals::sTPC_Pad_Parameter[row][4] + (col+0.5)*TPCGlobals::sTPC_Pad_Parameter[row][3] )*TMath::DegToRad();
  Double_t z = rad*TMath::Cos(theta)+m_tpcOffset;
  Double_t x = rad*TMath::Sin(theta);
  TVector2 vec(x,z);
  return vec;
}

TVector2 TPCIDHandler::GetXZ(int ID){
  int tRow = this->GetRow(ID);
  int tCol = this->GetCol(ID);
  return this->GetXZ( tRow, tCol );
}

double TPCIDHandler::GetWidth(int row, int col ){
  double dtheta = TPCGlobals::sTPC_Pad_Parameter[row][3]*TMath::DegToRad();
  double rad    = TPCGlobals::sTPC_Pad_Parameter[row][2];
  double width = rad*dtheta/2.;
  return width;
}
double TPCIDHandler::GetWidth(int ID){
  int row = this->GetRow( ID );
  int col = this->GetCol( ID );
  return this->GetWidth( row, col );
}
double TPCIDHandler::GetLength( int row, int col ){
  return TPCGlobals::sTPC_Pad_Parameter[row][5]/2.;
}

double TPCIDHandler::GetLength( int ID ){
  int row = this->GetRow(ID);
  int col = this->GetCol(ID);

  return this->GetLength( row, col );
}

bool TPCIDHandler::GetRowCol( double x, double z, int& row, int& col ){
  TVector2 vec(z-m_tpcOffset,x);
  Double_t rad = vec.Mod();
  Double_t phi = vec.Phi();//rad
  if( phi > TMath::Pi()){ phi = phi - 2*TMath::Pi();}
  bool hit = false;

  if( rad < TPCGlobals::sTPC_Pad_Parameter[0][2]-TPCGlobals::sTPC_Pad_Parameter[0][5]/2. ||
      rad > TPCGlobals::sTPC_Pad_Parameter[31][2]+TPCGlobals::sTPC_Pad_Parameter[31][5]/2.){
    row = -1;
    col = -1;
    return hit;
  }

  row = -1;
  for( int i = 0; i< 32; i++){
    if( rad < TPCGlobals::sTPC_Pad_Parameter[i][2]+TPCGlobals::sTPC_Pad_Parameter[i][5]/2.){
      row = i;
      break;
    }
  }
  col = -1;
  double starting_phi = (-180 + TPCGlobals::sTPC_Pad_Parameter[row][4])*TMath::DegToRad();
  double delta_phi    = TPCGlobals::sTPC_Pad_Parameter[row][3]*TMath::DegToRad();
  double starting_deg = (-180 + TPCGlobals::sTPC_Pad_Parameter[row][4]);
  double delta_deg    = TPCGlobals::sTPC_Pad_Parameter[row][3];
  //std::cout<< starting_phi << "\t" << delta_phi << "\t" << phi << "\t" << row << "\t" << (int)((phi - starting_phi)/delta_phi) << std::endl;

  if( TMath::Abs(phi) > TMath::Abs(starting_phi) ){
    col = -2;
    return hit;
  }

  /*
  std::cout << starting_deg << "\t" << delta_deg << "\t"
	    << phi*TMath::RadToDeg() << "\t"
	    << phi*TMath::RadToDeg() - starting_deg << "\t"
	    << x << "\t" << z
	    << std::endl;
  */


  col = (int)((phi - starting_phi)/delta_phi);
  hit = true;
  return hit;
}

bool TPCIDHandler::GetRowCol( int ID, int& row, int& col ){
  row = this->GetRow(ID);
  col = this->GetCol(ID);
  return true;
}

int  TPCIDHandler::GetID( int row, int col ){
  int ID = (int)(TPCGlobals::sTPC_PadRowMaxID[row]-TPCGlobals::sTPC_Pad_Parameter[row][1]+col);
  return ID;
}

int  TPCIDHandler::GetID( double x, double z ){
  int row;
  int col;
  bool hit = this->GetRowCol( x, z, row, col );
  if( !hit ){ return -1; }
  return this->GetID( row, col );
}
void TPCIDHandler::DrawPoints(){
  TGraph* gr = new TGraph();
  TVector2 pos;
  for( int i = 0; i< TPCGlobals::sTPC_NMaximum; i++){
    pos = this->GetXZ( i );
    gr->SetPoint( i, pos.Y(), pos.X() );
  }
  gr->Draw("AP");
}

double TPCIDHandler::GetPhi( int row, int col ){
  // center angle
  double phi = (-180 +TPCGlobals::sTPC_Pad_Parameter[row][4]) + (col+0.5)*TPCGlobals::sTPC_Pad_Parameter[row][3];
  return phi;
}
double TPCIDHandler::GetPhi( int ID ){
  int row;
  int col;
  this->GetRowCol( ID, row, col );
  return GetPhi( row, col );
}
double TPCIDHandler::GetRad( int row, int col ){
  double rad=TPCGlobals::sTPC_Pad_Parameter[row][2];
  return rad;
}
double TPCIDHandler::GetRad( int ID ){
  int row;
  int col;
  this->GetRowCol( ID, row, col );
  return GetRad( row, col );
}

/*
class TPCPoly : public TH2Poly {
 public:
  TPCPoly();
  TPCPoly(char* name, char* title );
  virtual ~TPCPoly();
  virtual void Reset{ TH2Poly::Reset("");}
  void Init();
  virtual int Fill( int ibin, double w = 1);
  virtual int Fill( double x, double z , double w = 1);
  virtual int FillY(double x, double z , double y, double e = 1);

 private:
  double m_Offset;
  double m_pLength;
  double m_sTheta;
  double m_dTheta;
  double m_cRad;
  int    m_nPad;
  double m_x[5];
  double m_y[5];
  double m_xPad;
  double m_yPad;

 public:
  ClassDef( TPCPoly, 0 )
};
*/
ClassImp(TPCPoly);
TPCPoly::TPCPoly(){
  m_Offset = -143;//mm
  this->Initialize(-300,300,-300,300,25,25);
  this->SetName("NoName");
  this->SetTitle("NoTitle");
  this->SetFloat();
  this->Init();
}
TPCPoly::TPCPoly( char* name, char* title ){
  m_Offset = -143;
  this->Initialize( -300, 300, -300, 300, 25, 25);
  this->SetName( name );
  this->SetTitle( title );
  this->SetFloat( kFALSE );
  this->Init();
}
TPCPoly::~TPCPoly(){
  TH2Poly::~TH2Poly();
}
void TPCPoly::Init(){
  if( gTPCIDHandler == NULL ){
    gTPCIDHandler = new TPCIDHandler();
  }
  for( int i = 0; i< 32; i++){
    m_pLength = TPCGlobals::sTPC_Pad_Parameter[i][5];
    m_sTheta  = (-180+TPCGlobals::sTPC_Pad_Parameter[i][4])*TMath::DegToRad();
    m_dTheta  = TPCGlobals::sTPC_Pad_Parameter[i][3]*TMath::DegToRad();
    m_cRad    = TPCGlobals::sTPC_Pad_Parameter[i][2];
    m_nPad    = TPCGlobals::sTPC_Pad_Parameter[i][1];
    for( int j = 0; j< m_nPad; j++){
      //for( int j = 0; j< 1; j++){
      m_x[0] = (m_cRad+(m_pLength/2.))*TMath::Cos(j*m_dTheta+m_sTheta)     + m_Offset;
      m_x[1] = (m_cRad+(m_pLength/2.))*TMath::Cos((j+1)*m_dTheta+m_sTheta) + m_Offset;
      m_x[2] = (m_cRad-(m_pLength/2.))*TMath::Cos((j+1)*m_dTheta+m_sTheta) + m_Offset;
      m_x[3] = (m_cRad-(m_pLength/2.))*TMath::Cos(j*m_dTheta+m_sTheta)     + m_Offset;
      m_x[4] = m_x[0];

      m_y[0] = (m_cRad+(m_pLength/2.))*TMath::Sin(j*m_dTheta+m_sTheta);
      m_y[1] = (m_cRad+(m_pLength/2.))*TMath::Sin((j+1)*m_dTheta+m_sTheta);
      m_y[2] = (m_cRad-(m_pLength/2.))*TMath::Sin((j+1)*m_dTheta+m_sTheta);
      m_y[3] = (m_cRad-(m_pLength/2.))*TMath::Sin(j*m_dTheta+m_sTheta);
      m_y[4] = m_y[0];
      this->AddBin(5, m_x, m_y);
      //this->AddBin(4, m_x, m_y);
    }

  }
}
int TPCPoly::Fill( int ibin, double w){
  TVector2 position = gTPCIDHandler->GetXZ( ibin );
  return TH2Poly::Fill( position.Y(), position.X(), w);
}
int TPCPoly::Fill( double x, double z, double w){
  return TH2Poly::Fill( z, x, w);
}
int TPCPoly::FillY( double x, double z, double y, double e){
  return 0;
}

/*
class TPCPolyRTheta : public TH2Poly {
 public:
  TPCPolyRTheta();
  TPCPolyRTheta(char* name, char* title);
  virtual ~TPCPolyRTheta();
  virtual void Reset{ TH2Poly::Reset("");}
  void Init();
  virtual int Fill( int ibin, double w = 1 );
  virtual int Fill( double x, double z, double w = 1);
  virtual int FillY( double x, double z, double y, double e = 1 );
 private:
  double m_pLength;
  double m_sTheta;
  double m_dTheta;
  double m_cRad;
  int    m_nPad;
  double m_x[5];
  double m_y[5];
  double m_xPad;
  double m_yPad;

 public:
  ClassDef( TPCPolyRTheta, 0 )
};
*/
TPCPolyRTheta::TPCPolyRTheta(){
  this->Initialize( -180,180, 0, 450, 25, 25 );
  this->SetName("NoName");
  this->SetTitle("NoTitle");
  this->SetFloat();
  this->Init();
}
TPCPolyRTheta::TPCPolyRTheta(char* name, char* title ){
  this->Initialize( -180,180, 0, 450, 25, 25 );
  this->SetName(name);
  this->SetTitle(title);
  this->SetFloat(kFALSE);
  this->Init();
}

TPCPolyRTheta::~TPCPolyRTheta(){
  ;
  //TH2Poly::~TH2Poly();
}

void TPCPolyRTheta::Init(){
  if( gTPCIDHandler == NULL ){
    gTPCIDHandler = new TPCIDHandler();
  }
  for( int i = 0; i<32; i++){
    double startingAngle = -180+TPCGlobals::sTPC_Pad_Parameter[i][4];
    int    ndiv          = TPCGlobals::sTPC_Pad_Parameter[i][1];
    double centerR       = TPCGlobals::sTPC_Pad_Parameter[i][2];
    double padLength     = TPCGlobals::sTPC_Pad_Parameter[i][5];
    for( int j = 0; j< ndiv; j++){
      m_x[0] = startingAngle + 2*TMath::Abs(startingAngle)*j/ndiv;
      m_x[1] = startingAngle + 2*TMath::Abs(startingAngle)*j/ndiv;
      m_x[2] = startingAngle + 2*TMath::Abs(startingAngle)*(j+1)/ndiv;
      m_x[3] = startingAngle + 2*TMath::Abs(startingAngle)*(j+1)/ndiv;
      m_x[4] = m_x[0];
      m_y[0] = centerR - padLength/2.;
      m_y[1] = centerR + padLength/2.;
      m_y[2] = centerR + padLength/2.;
      m_y[3] = centerR - padLength/2.;
      m_y[4] = m_y[0];
      this->AddBin(5,m_x,m_y);
    }
  }
}

int TPCPolyRTheta::Fill( int ibin, double w ){
  double r(0);
  double theta(0);
  int row(-1);
  int col(-1);
  gTPCIDHandler->GetRowCol(ibin, row, col);
  if( row < 0  || row >= 32 ){return -1;}

  r = TPCGlobals::sTPC_Pad_Parameter[row][2];
  double stheta = -180 + TPCGlobals::sTPC_Pad_Parameter[row][3];
  double dtheta = 2*TMath::Abs(stheta)/TPCGlobals::sTPC_Pad_Parameter[row][1];
  theta = stheta + dtheta*col;
  return TH2Poly::Fill( theta, r , w);
}

int TPCPolyRTheta::Fill( double x, double z, double w ){
  TVector2 vec( z - 143, x);
  double r = vec.Mod();
  double theta = vec.Phi()*TMath::RadToDeg();
  return TH2Poly::Fill( theta, r, w );
}
int TPCPolyRTheta::FillY( double x, double z, double y, double e){
  return 0;
}

ClassImp(TPCRowYHistArray)
TPCRowYHistArray* gTPCConvHistArray = new TPCRowYHistArray("gTPCConvHistArray",0);

TPCRowYHistArray::TPCRowYHistArray(char* hisName, double thres){
  m_hisName   = hisName;
  m_Threshold = thres;
  Init();
}

TPCRowYHistArray::~TPCRowYHistArray(){
  ;
}

void TPCRowYHistArray::Init(){
  for( int i = 0; i< 32; i++){
    hisRowY[i] = new TH2D(Form("%s_%d",m_hisName.c_str(),i),
			  Form("%s:%d;ColID;Y[mm]",m_hisName.c_str(),i),
			  (int)TPCGlobals::sTPC_Pad_Parameter[i][1],0,(int)TPCGlobals::sTPC_Pad_Parameter[i][1],
			  300,-300,300);
  }
}
void TPCRowYHistArray::Fill( int row, int col, double y, double e ){
  hisRowY[row]->Fill(col,y,e);
}
void TPCRowYHistArray::Reset(){
  for( int i = 0; i< 32; i++){
    hisRowY[i]->Reset();
  }
}

void TPCRowYHistArray::Clustering(){
  Double_t ydist;

}
void TPCRowYHistArray::ExportData(){
  ;
}
