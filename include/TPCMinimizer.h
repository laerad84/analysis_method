#ifndef TPCMinimizer__H__
#define TPCMinimizer__H__
#include <iostream>
#include <vector>

#include "TObject.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TVirtualFitter.h"
#include "TArc.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
class TPCCircleFitResult : public TObject {
 public:
  TPCCircleFitResult( );
  TPCCircleFitResult(TVector3 par, TVector3 parErr);
  virtual ~TPCCircleFitResult();

  TVector3 GetFitPar() const { return m_fitpar; }
  TVector3 GetFitParErr() const { return m_fitparErr; }
  Double_t CalculateChi2( double x, double z);
  Double_t CalculateChi2( TVector2 vec);
  Double_t CalculateDist( double x, double z);
  Double_t CalculateDeltaY( TVector2 pos );
  void     CalculateStartingMomentum(double tesla = 1.0);
  void     SetID( int id ){ m_id = id; }
  void     SetCoupleID( Int_t id ){ m_coupleID = id;}
  void     SetMinimumY( double y ){ m_minimumY = y;}
  void     SetDyDl( double dydl){ m_DyDl = dydl; }
  void     SetDyOff( double dyOff){ m_DyOff = dyOff;}
  void     SetStartingPoint(TVector3 vec){ m_StartingPoint = vec;}
  void     SetStartingMomentum( TVector3 vec0, TVector3 vec1){
    m_StartingMomentum[0] = vec0;
    m_StartingMomentum[1] = vec1;
  }
  Int_t    GetID()            const { return m_id;}
  Int_t    GetCoupleID()      const { return m_coupleID;}
  Double_t GetMinimumY()      const { return m_minimumY;}
  TVector3 GetStartingPoint() const { return m_StartingPoint; }
  void     GetStartingMomentum( TVector3& vec0, TVector3& vec1){ vec0 = m_StartingMomentum[0]; vec1 = m_StartingMomentum[1];}
  Double_t GetDyDl()  const { return m_DyDl;}
  Double_t GetDyOff() const { return m_DyOff;}

 private:
  Int_t    m_id;
  TVector3 m_fitpar;
  TVector3 m_fitparErr;
  double   m_Chi;
  double   m_DyDl;
  double   m_DyOff;

  Int_t    m_coupleID;
  Double_t m_minimumY;
  TVector3 m_StartingPoint;
  TVector3 m_StartingMomentum[2];


 public:
  ClassDef( TPCCircleFitResult, 1 )
};


////// new TPCCircleFitResult  under construction
/*
class TPCCircleFitResult_new : public TObject {
public:
	TPCCircleFitResult_new();
	TPCCircleFitResult_new(TVector3 par, TVector3 parErr);
	TPCCircleFitResult_new(Int_t ID, TVector3 par, TVector3 parErr );
	virtual TPCCircleFitResult_new();
	
	Int_t ID()				const { return m_id; }
	TVector3 FitPar()		const { return m_fitpar; }
	TVector3 FitParErr()	const { return m_fitparErr; }
	Double_t Chi()			const { return m_Chi; }
	Double_t DyDl()			const { return m_DyDl; }
	Double_t DyOff()		const { return m_DyOff; }
	Double_t TrackLength()	const { return m_trackLength; }
	Double_t totalCharge()	const { return m_trackLength; }

	Int_t    CoupleID()		const { return m_coupleID; }
	Double_t MinimumY()		const { return m_minimumY; }
	TVector3 StartingPoint()const { return m_StartingPoint; }
	void     StartingMomentum( TVector3& vec0, TVector3& vec1 ){ vec0 = m_StartingMomentum[0]; vec1 = m_StartingMomentum[1]; }
	
	void SetID( Int_t id ){ m_ID = id; }
	void SetParameters( TVector3 par, TVector3 parErr ){ m_fitpar = par; m_fitparErr = parErr; }
	void SetChi( Double_t chi ){ m_Chi = chi; }
	void SetDyDl( Double_t DyDl ){ m_DyDl = DyDl; }
	void SetDyOff( Double_t DyOff ){ m_DyOff = DyOff; }
	void SetTrackLength( Double_t trackLength ){ m_trackLength = trackLength; }
	void SetTotalCharge( Double_t totalCharge ){ m_totalCharge = totalCharge; }
	void SetCoupleID( Int_t coupleID ){ m_coupleID = coupleID; }
	void SetStartingPoint( TVector3 startingPoint ){ m_StartingPoint = startingPoint; }
	void SetStartingMomentum( TVector3 startingMomentum0, TVector3 startingMomentum1){ m_StartingMomentum[0] = startingMomentum0; m_StartingMomentum[1] = startingMomentum1; }
	
	Double_t CalculateDist( double x , double z);
	Double_t CalculateY( double x, double z );
	Double_t CalculateY( TVector2 pos );
	Double_t CalculateStartingMomentum(double tesla = 1.0 );
	
	
private:
	/// basic parameters
	Int_t m_id;           //Circle fit id, same with track id
	TVector3 m_fitpar;    // fitting result x : radius, y : x pos on TPC z : z pos on TPC
	TVector3 m_fitparErr; // fitting result error
	Double_t m_Chi;       // chisquare
	Double_t m_DyDl;      // Y changes by length. length : length of arc from z axis
	Double_t m_DyOff;     // Y offset at length = 0;
	Double_t m_trackLength; // track length ( point to point )
	Double_t m_totalCharge; // total charge of track
	std::vector<TVector3> clusterArr; // Cluster position information
	std::vector<Double_t> eArr;       // Cluster energy deposit information
	
	/// information of coupled result
	Int_t    m_coupleID;            //Coupled circle's ID;
	Double_t m_minimumY;            // y dist at crossing point
	TVector3 m_StartingPoint;       // crossing point
	TVector3 m_StartingMomentum[2]; // momentum at crossing point
	
	Bool_t   m_IsCoupled;   // Flag for coupling
	Bool_t   m_IsEvaluated; // Flag for evaluation
	
public:
	ClassDef( TPCCircleFitResult_new, 1 )
};
*/




class TPCCircleFitter : public TObject{
 public:
  TPCCircleFitter();
  ~TPCCircleFitter();
  void PrintAll();
  void SetPoints(int n, double* xarr, double* zarr);
  void SetPoints(int n, double* xarr, double* zarr, double* earr);
  void SetPoints(int n, double* xarr, double* zarr, double* earr, double* xarrErr, double* zarrErr);

  void SetPoints(std::vector<double> xarr, std::vector<double> zarr);


  /*
  int  SetPoint( double x, double z );
  int  SetPoint( double x, double z , double e);
  int  SetPoint( double x, double z , double e, double xerr, double zerr);
  */

  void ProcessFit( double initial_r = 300,double initial_err_r = 10,  double initial_x = 0, double initial_err_x = 50, double initial_y = 0 , double initial_err_y = 50 );
  void PrintFitResult();
  void Init();

  TPCCircleFitResult Fit( std::vector<double> xarr, std::vector<double> zarr);
  TPCCircleFitResult ThreePointCircleFit(std::vector<double> xarr, std::vector<double> zarr);

  static void    Chi2Minimizer( Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t );
  static void    Chi2MinimizerWithWeight( Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t );
  static         std::vector<Double_t> EArr;
  static int     value;
  static         TGraphErrors* gr;
  static void    TestFunction();

 private:
  ////// fit parameters //////

  ////// Fit result //////
  Double_t m_fit_r;
  Double_t m_fit_x;
  Double_t m_fit_z;
  Double_t m_fitErr_r;
  Double_t m_fitErr_x;
  Double_t m_fitErr_z;
  Double_t m_fitChisq;

  /*{
    std::cout<< "Static" << std::endl;
    }*/
  ClassDef( TPCCircleFitter, 1)
};

/*
class TPCCircleFitter : public TObject{
 public :
  TPCCircleFitter();
  TPCCircleFitter(double magneticField);
  ~TPCCircleFitter();

  static TGraphErrors* sTrack;
  static void chi2minimizer( Int_t &, Double_t *, Double_t &f,Double_t *par, Int_t )

};
*/


class TPCCircleTrack : public TObject {
 public:
  TPCCircleTrack();
  ~TPCCircleTrack();
  void GetCirCle(std::vector<double> x, std::vector<double> z, std::vector<double> fitResult);




 private:

 public:
  ClassDef( TPCCircleTrack, 1)
};

#endif
