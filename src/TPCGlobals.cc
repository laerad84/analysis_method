#include "TPCGlobals.h"
ClassImp( TPCGlobals )
const Double_t TPCGlobals::sTPC_Pad_Parameter[32][6] = {
  {0,	48,	14.5,	7.5	,0.	,9.},
  {1,	48,	24.,	7.5	,0.	,9.},
  {2,	72,	33.5,	5.	,0.	,9.},
  {3,	96,	43.,	3.75	,0.	,9.},
  {4,	120,	52.5,	3.	,0.	,9.},
  {5,	144,	62.,	2.5	,0.	,9.},
  {6,	168,	71.5,	2.14286	,0.	,9.},
  {7,	192,	81.,	1.875	,0.	,9.},
  {8,	216,	90.5,	1.66667	,0.	,9.},
  {9,	240,	100.,	1.5	,0.	,9.},
  {10,	208,	111.25,	1.49375	,24.65	,12.5},
  {11,	218,	124.25,	1.32844	,35.2	,12.5},
  {12,	230,	137.25,	1.2	,42	,12.5},
  {13,	214,	150.25,	1.09093	,63.27	,12.5},
  {14,	212,	163.25,	1.      ,74	,12.5},
  {15,	214,	176.25,	0.923084,	81.23	,12.5},
  {16,	220,	189.25,	0.857182,	85.71	,12.5},
  {17,	224,	202.25,	0.801786,	90.2	,12.5},
  {18,	232,	215.25,	0.751552,	92.82	,12.5},
  {19,	238,	228.25,	0.707227,	95.84	,12.5},
  {20,	244,	241.25,	0.667869,	98.52	,12.5},
  {21,	232,	254.25,	0.632672,	106.61	,12.5},
  {22,	218,	267.25,	0.60101	,	114.49	,12.5},
  {23,	210,	280.25,	0.573238,	119.81	,12.5},
  {24,	206,	293.25,	0.547111,	123.648	,12.5},
  {25,	202,	306.25,	0.523267,	127.15	,12.5},
  {26,	200,	319.25,	0.5014  ,	129.86	,12.5},
  {27,	196,	332.25,	0.481327,	132.83	,12.5},
  {28,	178,	345.25,	0.463371,	138.76	,12.5},
  {29,	130,	358.25,	0.446154,	151	,12.5},
  {30,	108,	371.25,	0.430185,	156.77	,12.5},
  {31,	90,	384.25,	0.415333,	161.31	,12.5}};

const Int_t    TPCGlobals::sTPC_PadRowMaxID[32]={
  48  ,96  ,168 ,264 ,384 ,528 ,696 ,888 ,
  1104,1344,1552,1770,2000,2214,2426,2640,
  2860,3084,3316,3554,3798,4030,4248,4458,
  4664,4866,5066,5262,5440,5570,5678,5768
};
const Int_t    TPCGlobals::sTPC_NMaximum      = 5768;
const Int_t    TPCGlobals::sTPC_NEperEnergy   = 2/0.00015;// nE/MeV
const Double_t TPCGlobals::sTPC_zOffset       = -143;//mm
const Int_t    TPCGlobals::sTPC_GEM_DefusionR = 20;//mm
