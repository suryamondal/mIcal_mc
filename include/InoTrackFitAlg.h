#ifndef INOTRACKFITALG_H
#define INOTRACKFITALG_H
//AlgFitTrackCam

#include <vector>
using std::vector;
#include "TVector3.h"
#include <TRandom3.h>
#include <InoMuRange.h>
#include "TGeoManager.h"
#include "vect_manager.h"
#include "micalDetectorParameterDef.hh"
#include "MultiSimAnalysis.hh"
#include "micalFieldPropagator.hh"
#include "SwimSwimmer.h"
//#include "TGeoManager.h"

//#include "G4ParticleDefinition.hh"
//#include "G4MuonPlus.hh"
//#include "G4MuonMinus.hh"
//#include "G4Material.hh"
//#include "G4VEmFluctuationModel.hh"
//#include "G4EnergyLossForExtrapolator.hh"

#include <string>
#include <cstdlib>

const unsigned int doubleLa=500;
const unsigned int shiftLa=250;
class InoTrackCand;

//typedef struct {
//  InoHit* csh;
//}HitStruct;

typedef struct {
  InoCluster* csh;
}ClustStruct;


typedef struct{
  bool   Straight;
  double XPos;
  double YPos;
  double ZPos;
  int    PlaneView;
  double XPosErrSq;
  double YPosErrSq;
  int    numInList;
  double cltime;
}TrkDataStruct;

typedef struct{
  double x_k0;
  double x_k1;
  double x_k2;
  double x_k3;
  double x_k4;
  int    x_k5; //withcls;
  bool   x_k6; //true for cluster in linear part, false : bended back
}FiltDataStruct;


class InoTrackFitAlg 
{
 public:
  
  InoTrackFitAlg();
  virtual ~InoTrackFitAlg();
  
  virtual void RunAlg();
  void InitialFramework();
  void RunTheFitter();
  
  void StoreFilteredData(const int NewPlane);
  void StoreFilteredData_sr(const int NewPlane, double*, bool);
  void FillGapsInTrack();
	
  void GetFitData(int& Plane1, int& Plane2);
  
  void ShowerStrips();
  void RemoveTrkHitsInShw();
  void ShowerSwim();
  
  void GoBackwards(const bool first);
  void GoForwards(const bool first);
  
  void GetPropagator(double *istate, double Bx, double By, double dBxbx, double dBxdy, double dBydx, double dBydy, double dz, TGeoMaterial* material);
  
  bool Swim(double* StateVector, double* Output, const int Plane, const int NewPlane, const bool GoForward, double* dS=0, double* Range=0, double* dE=0);
  bool Swim(double* StateVector, double* Output, const double zbeg, const int NewPlane, const bool GoForward, double* dS=0, double* Range=0, double* dE=0);
  bool Swim(double* StateVector, double* Output, const int Plane, const double zend, const bool GoForward, double* dS=0, double* Range=0, double* dE=0);
  bool Swim_new(double* StateVector, double* Output, const int Plane,  int& NewPlane, const bool GoForward, double* dS=0, double* Range=0, double* dE=0);
  
  void TrackElementMerging(double *Tr1, double TargetZ, double *Tr2);
  void GetInitialCovarianceMatrix(const bool FirstIteration);
  bool PredictedStateCov(double* StateVector, const int Plane, int& NewPlane, const bool GoForward, double *ax__minus, int isHalf=0, double* dS=0, double* Range=0); //, double* dE=0);
  double GetEnergyLoss(double* istate, double dz, double &axi, double &aT_max, double &aI, TGeoMaterial* material);
  void GetMultipleScattering(double* mstate,double Bx,double By, double dz, /*double axi,*/ double aT_max, double aI, TGeoMaterial* material);
  void ExtrapCovMatrix();  //Put it in C_k_intermediate
  void ExtrapCovMatrixall(); //Change C_k_minus[lm][kl] also
  
  void CalcKalmanGain(double *x__minus,const int NewPlane);
  //void UpdateStateVector(const int Plane, const int NewPlane, const bool GoForward);
  
  //void UpdateStateVector_new(const int Plane, const int NewPlane, double* Output, const bool GoForward);
  void KalmanFilterStateVector(double *x__minus, const int Plane, const bool GoForward, double *xk);
  
  void UpdateCovMatrix(const int NewPlane);
  void MoveArrays(const int NewPlane, const bool GoForward);
  //void CheckValues(double* Input, const int NewPlane);
  
  void SetTrackProperties( double* Input); //CandFitTrackCamHandle &cth);
  /*
    void SetPropertiesFromFinderTrack(InotTrackCand &cth);
    void SpectrometerSwim(CandFitTrackCamHandle &cth);
  */
  //void SetRangeAnddS(); //CandFitTrackCamHandle& cth);
  void TimingFit(); //CandFitTrackCamHandle &cth);
  
  //  bool NDPlaneIsActive(int plane, float u, float v, float projErr);
  virtual void Trace(const char *c) const;
  
  void ResetCovarianceMatrix();
  
  void SetT( );
  void CalculateTrace(){;};
  
  bool DirectionFromFinderHits(InoTrack *trk, double& FinderPathLength, double& FinderDistance);
  bool DirectionFromFinderHitsOldFunc(InoTrack *trk, double& FinderPathLength, double& FinderDistance); 
  //Abhijit's Work ADB 2015/05/06
  bool CheckFCPC(double *x_k, bool GoForward);
  int CheckFCPCUpOrDn(double *x_k, bool DirExtraPol, int MaxMinPlane, bool GoDir);

 private:
  
  vector<ClustStruct> SlcClustData[doubleLa]; //GMA put a very large value, but need to be put from database 
  vector<ClustStruct> InitTrkClustData[doubleLa]; //Only Finder track cluster
  
  vector<TrkDataStruct> TrkClustsData[doubleLa]; //TrkHitsData[doubleLa]; // TrkStripData[150];
  vector<FiltDataStruct> FilteredData[doubleLa];
  
  double ZPosLayer[doubleLa];
  double IcalX;
  double IcalY;

  int checkfcorpc;
  int FCorPC; //=0;
  int FCorPCForward; // FCorPCUp; //=0;
  int FCorPCBackward; //FCorPCDn; // = 0;
  /* bool CheckMatUp;//   =0; */
  /* bool CheckMatDn;//   =0; */
  /* bool LocalPosUp;//   =0; */
  /* bool LocalPosDn;//   =0; */
  int FCPC;

  double ShiftInX;
  double ShiftInY;
  double ShiftInZ;

  Int_t nbfield;
  Double_t bave;
  Double_t t_bave;
  Bool_t EndofRange;
  Int_t EndofRangePlane;
  Bool_t LastIteration;
  double x_k4_biased;
  int UseGeoSwimmer;
  double x_k[6]; // x_k5 for tagging used hit or not 5];
  double x_k_minus[6]; // 5];
  double C_k[5][5];
  double C_k_minus[5][5];
  double C_k_intermediate[5][5];
  double F_k[5][5];
  double F_k_minus[5][5];
  double Q_k[5][5];
  double Q_k_minus[5][5];
  double K_k[5][2];
  double StateIter[5];
  
  int H_k[2][5];
  int Identity[5][5];
  
  double VtxCov[5];
  double EndCov[5];
  double prevstate[5];
  double prevpredn[5];
  
  int OtLStrip;
  int MaxPlane;
  int MinPlane;
  //  unsigned nhits;
  double DeltaZ;
  double DeltaPlane;  
  
  InoTrackCand* fTrackCand;
  
  bool debug_fit;
  
  bool ZIncreasesWithTime;
  bool FirstIteration;
  bool PassTrack;
  bool fMT;
  
  double xxin;
  double yyin;
  double txin;
  double tyin;
  double B_in;
  
  double ds;
  double drange;
  double I;
  double xi;
  double T_max;
  double BetheBloch;
  
  //  double L;
  //  double Lz;
  double GPL;
  double RNG;
  double ChiSquare;
  double MagicRatio;
  double EndState[5];
  double MaxPlaneData[6];
  double MinPlaneData[6];

  bool SaveData;
  bool SwimThroughShower;

  int ShowerEntryPlane;
  
  int NIter;
  int TotalNSwimFail;

  int NumFinderStrips;
  double MeanTrackTime;
  
  double StripListTime;

  InoFittedTrack_Manager* inoFittedTrack_pointer; 
  InoTrackCand_Manager* inoTrackCand_pointer;
  InoHit_Manager* inoHit_pointer; 
  InoCluster_Manager *InoCluster_pointer;
  MultiSimAnalysis *pAnalysis; 
  micalFieldPropagator *pFieldMap;
  const InoTrack* fFinderTrack;

  double StripXWidth;
  double StripYWidth;
  double LayerThickness;
  unsigned int    nLayer;
  int nHit;

  TVector3 shiftvector; //Shift in position, while track bend back to the same layer  
  
  TGeoManager* icalGeometry;
  TGeoManager* abc;
  TGeoMaterial *localmat;
  TRandom3 *PoissonRn;
  InoMuRange *IcalRange;

  double CorrTimeError;
  double UnCorrTimeError;
  double TimeError;

  double pargasxyz[3];
};
#endif   // ALGFITTRACKCAM_H
