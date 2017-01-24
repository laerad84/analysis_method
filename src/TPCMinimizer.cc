#include "TPCMinimizer.h"
ClassImp(TPCCircleFitter)
int TPCCircleFitter::value = 0;
std::vector<double> TPCCircleFitter::EArr;
TGraphErrors* TPCCircleFitter::gr = new TGraphErrors();

TPCCircleFitter::TPCCircleFitter() : m_fit_r(0), m_fit_x(0), m_fit_z(0), m_fitErr_x(1), m_fitErr_z(1), m_fitErr_r(1) {
  std::cout<< "Initialize" << std::endl;
}

TPCCircleFitter::~TPCCircleFitter(){

  std::cout<<"Destruct" << std::endl;
}

void TPCCircleFitter::PrintAll(){
  std::cout<< "PrintAll" << std::endl;
}

void TPCCircleFitter::TestFunction(){
  std::cout<< "Static" << std::endl;

}
void TPCCircleFitter::Chi2Minimizer( Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t ){

  Int_t nPoint = TPCCircleFitter::gr->GetN();
  f = 0;
  Double_t* x = TPCCircleFitter::gr->GetY();
  Double_t* z = TPCCircleFitter::gr->GetX();
  for( Int_t ipoint = 0; ipoint < nPoint; ipoint++){
    Double_t u = x[ipoint] - par[0];
    Double_t v = z[ipoint] - par[1];
    Double_t dr= par[2] - TMath::Sqrt( u*u + v*v );
    f += dr*dr;
  }
}

void TPCCircleFitter::Chi2MinimizerWithWeight( Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t  ){
  Int_t nPoint = TPCCircleFitter::gr->GetN();
  f = 0;
  Double_t* x = TPCCircleFitter::gr->GetY();
  Double_t* z = TPCCircleFitter::gr->GetX();
  Double_t* ex = TPCCircleFitter::gr->GetEY();
  Double_t* ez = TPCCircleFitter::gr->GetEX();
  for( Int_t ipoint = 0; ipoint < nPoint; ipoint++){
    Double_t u = x[ipoint] - par[0];
    Double_t v = z[ipoint] - par[1];
    Double_t theta = TMath::ATan2( u, v );
    Double_t Err_r = TMath::Sqrt(TMath::Power(ex[ipoint]*TMath::Sin(theta),2) + TMath::Power(ez[ipoint]*TMath::Cos(theta),2));
    Double_t dr = (par[2] - TMath::Sqrt( u*u + v*v ))/Err_r;
    f += dr*dr;
  }
}

void TPCCircleFitter::ProcessFit( double initial_r, double initial_err_r,
				  double initial_x, double initial_err_x,
				  double initial_z, double initial_err_z ){

  TVirtualFitter *fitter = TVirtualFitter::Fitter(0,3);
  TVirtualFitter::SetDefaultFitter("Minuit");
  fitter->SetFCN(TPCCircleFitter::Chi2Minimizer);
  fitter->SetParameter(0, "x0", initial_x,initial_err_x,0,0 );
  fitter->SetParameter(1, "z0", initial_z,initial_err_z,0,0 );
  fitter->SetParameter(2, "R" , initial_r,initial_err_r,0,0 );

  Double_t arglist[1] = {0};
  fitter->ExecuteCommand("MIGRAD", arglist, 0);
  std::cout<< fitter->GetParameter(0) << "\t"
	   << fitter->GetParameter(1) << "\t"
	   << fitter->GetParameter(2) << "\n";

  m_fit_r = fitter->GetParameter(2);
  m_fit_x = fitter->GetParameter(0);
  m_fit_z = fitter->GetParameter(1);
  m_fitErr_r = fitter->GetParError(2);
  m_fitErr_x = fitter->GetParError(0);
  m_fitErr_z = fitter->GetParError(1);

}

TPCCircleFitResult TPCCircleFitter::Fit(std::vector<double> xarr, std::vector<double> zarr){
  /// initial three point circle fitting ///
  std::vector<double> tmpxarr;
  std::vector<double> tmpzarr;
  tmpxarr.push_back( xarr.at(0) );
  tmpzarr.push_back( zarr.at(0) );
  Int_t middlePoint = (int)(xarr.size()/2);
  tmpxarr.push_back( xarr.at(middlePoint) );
  tmpzarr.push_back( zarr.at(middlePoint) );
  tmpxarr.push_back( xarr.at(xarr.size()-1));
  tmpzarr.push_back( zarr.at(zarr.size()-1));
  TPCCircleFitResult tmpResult = this->ThreePointCircleFit( tmpxarr, tmpzarr);

  this->SetPoints(xarr, zarr);
  this->ProcessFit(tmpResult.GetFitPar().X(),10, tmpResult.GetFitPar().Y(),10,tmpResult.GetFitPar().Z(), 10);

  TVector3 par(m_fit_r, m_fit_x, m_fit_z);
  TVector3 parErr(m_fitErr_r, m_fitErr_x, m_fitErr_z);
  TPCCircleFitResult result(par, parErr);
  return result;
}

/*
int  TPCCircleFitter::SetPoint( double x, double z){
  TPCCircleFitter::gr->SetPoint( TPCCircleFitter::gr->GetN(), z, x );
  return TPCCircleFitter::gr->GetN();
}

int TPCCircleFitter::SetPoint( double x, double z, double e){
  TPCCircleFitter::gr->SetPoint( TPCCircleFitter::gr->GetN(), z, x );
  EArr.push_back(e);
  return TPCCircleFitter::gr->GetN();
}
*/


void TPCCircleFitter::SetPoints( int n, double* xarr, double* zarr){
  for( int i = 0; i< n; i++){
    TPCCircleFitter::gr->SetPoint(i, zarr[i], xarr[i]);
  }
}
void TPCCircleFitter::SetPoints( int n, double* xarr, double* zarr, double* earr){
  TPCCircleFitter::EArr.clear();
  for( int i = 0; i< n; i++){
    TPCCircleFitter::gr->SetPoint(i, zarr[i], xarr[i]);
    TPCCircleFitter::EArr.push_back( earr[i]);
  }
}
void TPCCircleFitter::SetPoints( int n, double* xarr, double* zarr, double* earr, double* xarrErr, double* zarrErr){
  TPCCircleFitter::EArr.clear();
  for( int i = 0; i< n; i++){
    TPCCircleFitter::gr->SetPoint(i, zarr[i], xarr[i]);
    TPCCircleFitter::gr->SetPointError( i, zarrErr[i], xarrErr[i] );
    TPCCircleFitter::EArr.push_back( earr[i]);
  }
}
void TPCCircleFitter::SetPoints( std::vector<double> xarr, std::vector<double> zarr){
  Int_t n = xarr.size();
  if( n > zarr.size() ){
    n = zarr.size();
  }
  for( int i = 0; i< n; i++){
    TPCCircleFitter::gr->SetPoint(i, zarr.at(i), xarr.at(i));
  }
}
void TPCCircleFitter::PrintFitResult(){
  std::cout<< "Fit Result" << std::endl;
  std::cout<< "Chisq : " << m_fitChisq << std::endl;
  std::cout<< "r : " << m_fit_r << "\t" << m_fitErr_r << std::endl;
  std::cout<< "x : " << m_fit_x << "\t" << m_fitErr_x << std::endl;
  std::cout<< "z : " << m_fit_z << "\t" << m_fitErr_z << std::endl;
}
void TPCCircleFitter::Init(){
  TPCCircleFitter::gr->Set(0);
  EArr.clear();
  m_fit_r = 0;
  m_fit_x = 0;
  m_fit_z = 0;
  m_fitErr_r = 1;
  m_fitErr_x = 1;
  m_fitErr_z = 1;
  m_fitChisq = 0;
}

TPCCircleFitResult TPCCircleFitter::ThreePointCircleFit(std::vector<double> xarr, std::vector<double> zarr){
  Int_t iIndex = 0;
  Int_t mIndex = (int)((xarr.size() -1)/2.);
  Int_t fIndex = xarr.size()-1;
  if( xarr.size() > 3 ){
    fIndex = xarr.size() -2;
  }
  if( xarr.size() > 5 ){
    iIndex = 1;
  }
  std::vector<Int_t> IndexArr;
  IndexArr.push_back(iIndex);
  IndexArr.push_back(mIndex);
  IndexArr.push_back(fIndex);

  double r,x,z;

  Double_t A = (xarr.at(iIndex)*xarr.at(iIndex) - xarr.at(mIndex)*xarr.at(mIndex) + zarr.at(iIndex)*zarr.at(iIndex) - zarr.at(mIndex)*zarr.at(mIndex))/2.;
  Double_t B = (xarr.at(mIndex)*xarr.at(mIndex) - xarr.at(fIndex)*xarr.at(fIndex) + zarr.at(mIndex)*zarr.at(mIndex) - zarr.at(fIndex)*zarr.at(fIndex))/2.;

  x = ((zarr.at(mIndex) - zarr.at(iIndex))*B - (zarr.at(fIndex) - zarr.at(mIndex))*A) /
    ((zarr.at(fIndex) - zarr.at(mIndex))*(xarr.at(mIndex) - xarr.at(iIndex)) - (xarr.at(fIndex) - xarr.at(mIndex))*(zarr.at(mIndex) - zarr.at(iIndex)));
  z = ((xarr.at(mIndex) - xarr.at(iIndex))*B - (xarr.at(fIndex) - xarr.at(mIndex))*A) /
    ((xarr.at(fIndex) - xarr.at(mIndex))*(zarr.at(mIndex) - zarr.at(iIndex)) - (zarr.at(fIndex) - zarr.at(mIndex))*(xarr.at(mIndex) - xarr.at(iIndex)));

  double rsq = 0;
  for( int i = 0; i< 3; i++){
    std::cout<< (xarr.at(IndexArr.at(i)) - x)*(xarr.at(IndexArr.at(i)) - x) + (zarr.at(IndexArr.at(i)) -z)*(zarr.at(IndexArr.at(i)) - z) << std::endl;
    rsq = rsq + (xarr.at(IndexArr.at(i)) - x)*(xarr.at(IndexArr.at(i)) - x) + (zarr.at(IndexArr.at(i)) -z)*(zarr.at(IndexArr.at(i)) - z);
  }
  rsq = rsq/3.;
  r = TMath::Sqrt(rsq);
  TVector3 vrst( r, x, z);
  TVector3 vrstErr(0,0,0);
  TPCCircleFitResult result(vrst, vrstErr);
  return result;
}


ClassImp( TPCCircleFitResult)
TPCCircleFitResult::TPCCircleFitResult() : m_DyDl(0), m_DyOff(0){;}
TPCCircleFitResult::TPCCircleFitResult( TVector3 par, TVector3 parErr) : m_DyDl(0), m_DyOff(0){
  m_fitpar    = par;
  m_fitparErr = parErr;
}
TPCCircleFitResult::~TPCCircleFitResult(){
  ;
}

Double_t TPCCircleFitResult::CalculateChi2(double x, double z){
  Double_t chi2 = 0;
  if( m_fitparErr[0] != 0 &&
      m_fitparErr[1] != 0 &&
      m_fitparErr[2] != 0 ){
    chi2 = TMath::Power(TMath::Sqrt(( x - m_fitpar[1] )*( x - m_fitpar[1] ) + ( z - m_fitpar[2] )*( z - m_fitpar[2] )) - m_fitpar[0],2)/m_fitparErr[0]/m_fitparErr[0];
  }else{
    chi2 = TMath::Power(TMath::Sqrt(( x - m_fitpar[1] )*( x - m_fitpar[1] ) + ( z - m_fitpar[2] )*( z - m_fitpar[2] )) - m_fitpar[0], 2);
  }
  return chi2;
}
Double_t TPCCircleFitResult::CalculateChi2( TVector2 vec ){
  Double_t chi2 = 0;
  Double_t x = vec.X();
  Double_t z = vec.Y();
  if( m_fitparErr[0] != 0 &&
      m_fitparErr[1] != 0 &&
      m_fitparErr[2] != 0 ){
    chi2 = TMath::Power(TMath::Sqrt(( x - m_fitpar[1] )*( x - m_fitpar[1] ) + ( z - m_fitpar[2] )*( z - m_fitpar[2] )) - m_fitpar[0],2)/m_fitparErr[0]/m_fitparErr[0];
  }else{
    chi2 = TMath::Power(TMath::Sqrt(( x - m_fitpar[1] )*( x - m_fitpar[1] ) + ( z - m_fitpar[2] )*( z - m_fitpar[2] )) - m_fitpar[0], 2);
  }
  return chi2;
}

Double_t TPCCircleFitResult::CalculateDist( double x, double z ){
  Double_t dist = TMath::Abs( TMath::Sqrt(( x - m_fitpar[1] ) * ( x - m_fitpar[1] ) + (z - m_fitpar[2] )* ( z - m_fitpar[2])) -m_fitpar[0] );
  return dist;
}

Double_t TPCCircleFitResult::CalculateDeltaY( TVector2 pos ){
  Double_t theta = TMath::ATan2( pos.X() - m_fitpar[1], pos.Y() - m_fitpar[2]);
  Double_t l     = theta*m_fitpar[0];
  Double_t y     = l*m_DyDl + m_DyOff;
  return y;
}

void TPCCircleFitResult::CalculateStartingMomentum(Double_t tesla){
  TVector3 vec0(0,0,0);
  TVector3 vec1(0,0,0);
  Double_t circleMomentum = m_fitpar[0]*tesla/3.36;// MeV tesla
  Double_t tangentialY    = m_DyDl;
  TVector2 originPoint( m_fitpar[2], m_fitpar[1]);//z,x
  TVector2 startPoint( m_StartingPoint.Z(),m_StartingPoint.X());
  TVector2 pointingVec = startPoint - originPoint;
  Double_t angle = pointingVec.Phi();
  TVector2 v[2];
  Double_t py[2]={0};
  v[0].SetMagPhi(circleMomentum, angle + TMath::Pi()/2.); // left circle
  py[0] = circleMomentum*tangentialY;
  v[1].SetMagPhi(circleMomentum, angle - TMath::Pi()/2.); // right circle
  py[1] = circleMomentum*tangentialY;

  m_StartingMomentum[0].SetXYZ( v[0].Y(), py[0], v[0].X());
  m_StartingMomentum[1].SetXYZ( v[1].Y(), py[1], v[1].X());
}
