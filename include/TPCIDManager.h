#ifndef TPCIDManager__H__
#define TPCIDManager__H__
#include "TPCGlobals.h"
#include "TObject.h"
#include "TVector2.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH2Poly.h"
#include "TH2.h"
#include "TH1.h"
#include "TMath.h"
class TPCIDHandler {
 public:
  TPCIDHandler(double offset = -143);
  virtual ~TPCIDHandler();

  TVector2 GetXZ(int row, int col );
  TVector2 GetXZ(int ID);

  double GetWidth( int row, int col );
  double GetWidth( int ID );
  double GetLength( int row, int col );
  double GetLength( int ID );
  bool   GetRowCol( int id, int& row, int& col );
  bool   GetRowCol( double x, double z, int& row, int& col );
  int    GetID( double x, double z );
  int    GetID( int row, int col );
  int    GetRow( int ID );
  int    GetCol( int ID );
  double GetPhi( int row, int col );
  double GetPhi( int ID );
  double GetRad( int row, int col );
  double GetRad( int ID );

  /// Map reading version will be updated ///


  void   DrawPoints();

 private :
  double m_tpcOffset;
 public :
  ClassDef( TPCIDHandler, 0 )
};

R__EXTERN TPCIDHandler *gTPCIDHandler;
/*
class TPCPadHit {
 public:
  TPCPadHit(int padID);
  ~TPCPadHit();
  Int_t m_PadID;
  std::vector<double> m_yArr;
  std::vector<double> m_tArr;
  std::vector<double> m_eArr;
  int   AddHit(double y, double t, double e );
  void  Set(double y, double t, double e);
  void  Set(std::vector<double> y, std::vector<double> t, std::vector<double> e );
};
*/
class TPCPoly : public TH2Poly {
 public:
  TPCPoly();
  TPCPoly(char* name, char* title );
  virtual ~TPCPoly();
  virtual void Reset(){ TH2Poly::Reset("");}
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

class TPCPolyRTheta : public TH2Poly {
 public:
  TPCPolyRTheta();
  TPCPolyRTheta(char* name, char* title);
  virtual ~TPCPolyRTheta();
  virtual void Reset(){ TH2Poly::Reset("");}
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

class TPCRowYHistArray {
 public:
  TPCRowYHistArray(char* hisName, double thres);
  ~TPCRowYHistArray();

  TH2D* hisRowY[32];

  void Fill(int row, int col, double y, double e =1.);
  void Reset();
  void Clustering();
  void ExportData();
 private:
  void Init();
  std::string m_hisName;
  Double_t    m_Threshold;
 public:

  ClassDef( TPCRowYHistArray, 0 )
};
R__EXTERN TPCRowYHistArray *gTPCConvHistArray;

/*
class TPCPadViewer {
 public :
  TH2D* m_his_yt[32];//
  TPCPolyRTheta* m_his_rtheta;
  TPCPoly* m_his_xz;
 private :

 public :
  ClassDef( TPCPadViwer, 0 )
};

*/

#endif //TPCIDManager__H__
