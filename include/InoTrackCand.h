#ifndef INOTRACKCAND_H
#define INOTRACKCAND_H
//CandTrackHandle
#include "TObject.h"
#include <map>
class InoTrack;
class InoShowerCand;
class InoVertex;
class InoHit;
class InoTrackCand //: public CandRecoHandle
{

public:
  InoTrackCand();
  InoTrackCand(const InoTrackCand &cdh);
  InoTrackCand(InoTrack *trk, bool forward);
  virtual ~InoTrackCand();

  virtual void Trace(const char *c = "") const;

  void SetU(Int_t,Float_t);
  void SetV(Int_t,Float_t);
  void SetdS(Int_t,Float_t);
  void SetTrackPointXError(Int_t,Float_t);
  void SetTrackPointYError(Int_t,Float_t);
  void SetRange(Int_t plane,Float_t g_cm2);

  void Set2dS(Int_t,Float_t);
  void Set2Range(Int_t plane,Float_t g_cm2);

  int GetFCPC() const;
  void SetFCPC(int);


  void SetT(Int_t,Double_t);

  void SetVtxTrace(Double_t);
  void SetVtxTraceZ(Double_t);
  void SetVtxnActiveUpstream(Int_t);
  void SetEndTrace(Double_t);
  void SetEndTraceZ(Double_t);
  void SetEndnActiveDownstream(Int_t);
  void SetVtxDistToEdge(Double_t);
  void SetEndDistToEdge(Double_t);

  virtual void ClearMaps(); // clears all STL maps, including fdS and fRange
  virtual void ClearUVT(); // same as ClearMap, left for backward compatibility

  //  Bool_t IsXPosValid(Int_t) const; // returns true if U and V positions have been set for this plane
  //  Bool_t IsXPosValid(Int_t) const; // returns true if U and V positions have been set for this plane
  Float_t GetU(Int_t) const; // U position at specified plane
  Float_t GetV(Int_t) const; // V position at specified plane
  Float_t GetZ(Int_t) const; // Z position at specified plane

  //  Int_t GetNTrackPlane(PlaneView::PlaneView_t = PlaneView::kUnknown) const;

  Float_t GetTrackPointXError(Int_t) const; // X-position error at specified plane
  Float_t GetTrackPointYError(Int_t) const; // Y-position error at specified plane
  Double_t GetT(Int_t) const; // time at specified plane
  //  Double_t GetT(Int_t,StripEnd::StripEnd_t) const; // time at specified plane
  //  Double_t GetT(StripEnd::StripEnd_t,Int_t) const; // time at specified plane
  Float_t GetdS(Int_t = -1) const; // travel distance from vertex
  Float_t GetRange(Int_t = -1) const; // g/cm**2 from vertex

  Float_t Get2dS(Int_t = -1) const; // travel distance from vertex
  Float_t Get2Range(Int_t = -1) const; // g/cm**2 from vertex



  Double_t GetVtxTrace() const;
  Double_t GetVtxTraceZ() const;
  Int_t GetVtxnActiveUpstream() const;
  Double_t GetEndTrace() const;
  Double_t GetEndTraceZ() const;
  Int_t GetEndnActiveDownstream() const;
  Double_t GetVtxDistToEdge() const;
  Double_t GetEndDistToEdge() const;

  virtual Double_t GetScore() const; // larger is better

  //  Int_t IsInShower(CandStripHandle *) const;
  //  void SetInShower(CandStripHandle *, Int_t);

  Double_t GetMomentum() const;
  void SetMomentum(Double_t);
 
  Double_t GetTheta() const;
  void SetTheta(Double_t);
  
  Double_t GetThErr() const;
  void SetThErr(Double_t);  

  Double_t GetPhi() const;
  void SetPhi(Double_t);
  
  Double_t GetPhErr() const;
  void SetPhErr(Double_t);

  Double_t GetdSExtra() const;
  void SetdSExtra(Double_t);

  Double_t GetRangeExtra() const;
  void SetRangeExtra(Double_t);

  Int_t GetFitType() const {return fType;}
  void SetFitType(Int_t typ) {fType = typ;}

  

  //  Bool_t IsUnphysical(Float_t trkFrac=0.660,Float_t asymCut=0.500,
  //                      Float_t xtalkFrac=0.660,Float_t xtalkCut=2.0){return false;};
  Bool_t IsContained();

  void SetNTrackStrip(Int_t);
  void SetNTrackDigit(Int_t);
  void SetNTimeFitDigit(Int_t);
  void SetTimeFitChi2(Double_t);
  void SetTimeForwardFitRMS(Double_t);
  void SetTimeForwardFitNDOF(Int_t);
  void SetTimeBackwardFitRMS(Double_t);
  void SetTimeBackwardFitNDOF(Int_t);

  Int_t GetNTrackStrip() const;
  Int_t GetNTrackDigit() const;
  Int_t GetNTimeFitDigit() const;
  Double_t GetTimeFitChi2() const;
  Double_t GetTimeForwardFitRMS() const;
  Int_t GetTimeForwardFitNDOF() const;
  Double_t GetTimeBackwardFitRMS() const;
  Int_t GetTimeBackwardFitNDOF() const;

  mutable map<const InoHit*, Int_t> fInShower;
  mutable map<Int_t,Float_t> fUPos;
  mutable map<Int_t,Float_t> fVPos;
  mutable map<Int_t,Float_t> fdS;                           // in meters
  mutable map<Int_t,Float_t> fRange;                       // in g/cm**2
  mutable map<Int_t,Float_t> fXPosError;
  mutable map<Int_t,Float_t> fYPosError;  
  mutable map<Int_t,Double_t> fTime[2];

  // in meters for being points, which is not used to 
  //update track parameter, but included in tracks
  mutable map<Int_t,Float_t> f2dS;
  mutable map<Int_t,Float_t> f2Range;    // in g/cm**2

  Double_t fVtxTrace;
  Double_t fVtxTraceZ;
  Double_t fEndTrace;
  Double_t fEndTraceZ;
  Double_t fVtxDistToEdge;
  Double_t fEndDistToEdge;
  Int_t fVtxnActiveUpstream;
  Int_t fEndnActiveDownstream;
  
  Int_t fNTrackStrip;               // # of strips that have InShower<=1
  Int_t fNTrackDigit;               // # of digits that have InShower<=1
  Int_t fNTimeFitDigit;          // # of digits used to determine timing
  Double_t fTimeFitChi2;
  Double_t fTimeForwardFitRMS;
  Int_t fTimeForwardFitNDOF;
  Double_t fTimeBackwardFitRMS;
  Int_t fTimeBackwardFitNDOF;

  Double_t fTimeSlope;
  Double_t fTimeOffset;
  Double_t fMomentum;
  Double_t fTheta;
  Double_t fPhi;
  Double_t fErrTh;
  Double_t fErrPh;

  Double_t fdSExtra;
  Double_t fRangeExtra; 


  InoTrack*     fTrack; 
  InoVertex*    fVertex;
  InoVertex*    fTerm;

  
  public : 
  Int_t GetNDaughters() const; 
  //////////////////////////////////////////////////
  //From  CandRecoHandle
  Int_t GetNStrip(Int_t i) const; 

  // intersection between two candreco objects
  //  Int_t GetNStrip(const CandRecoHandle *, 
  //                  PlaneView::PlaneView_t = PlaneView::kUnknown) const;

  Int_t GetNDigit(Int_t) const;
  
  Int_t GetNPlane(Int_t) const;
  
  Int_t GetBegPlane(Int_t) const;
  Int_t GetEndPlane(Int_t) const;   

  void SetVtxU(Double_t);
  Double_t GetVtxU() const;
  
  void SetVtxV(Double_t);
  Double_t GetVtxV() const;

 void SetVtxZ(Double_t);
  Double_t GetVtxZ() const;

  void SetVtxT(Double_t);
  Double_t GetVtxT() const;

  void SetVtxPlane(Int_t);
  Int_t GetVtxPlane() const;

  void SetEndU(Double_t);
  Double_t GetEndU() const;
  
  void SetEndV(Double_t);
  Double_t GetEndV() const;

  void SetEndZ(Double_t);
  Double_t GetEndZ() const;

  void SetEndT(Double_t);
  Double_t GetEndT() const;

  void SetEndPlane(Int_t);
  Int_t GetEndPlane() const;

  void SetVtxDirCosU(Double_t);
  Double_t GetVtxDirCosU() const;

  void SetVtxDirCosV(Double_t);
  Double_t GetVtxDirCosV() const;

  void SetVtxDirCosZ(Double_t);
  Double_t GetVtxDirCosZ() const;

  void SetEndDirCosU(Double_t);
  Double_t GetEndDirCosU() const;

  void SetEndDirCosV(Double_t);
  Double_t GetEndDirCosV() const;

  void SetEndDirCosZ(Double_t);
  Double_t GetEndDirCosZ() const;

  void SetDirCosU(Double_t);
  Double_t GetDirCosU() const;

  void SetDirCosV(Double_t);
  Double_t GetDirCosV() const;

  void SetDirCosZ(Double_t);
  Double_t GetDirCosZ() const;
 
  void SetTimeSlope(Double_t);
  Double_t GetTimeSlope() const;

  void SetTimeOffset(Double_t);
  Double_t GetTimeOffset() const;

  //  CalTimeType::CalTimeType_t GetCalTimeType() const;

  //  void CalibrateSigMapped(UInt_t fEncoded,Float_t);
  //  void CalibrateMIP(UInt_t fEncoded,Float_t);

  Double_t GetPulse() const; // CalStripType::CalStripType_t = CalStripType::kMIP) const;
  Double_t GetPlanePulse(Int_t iplane) const; //, CalStripType::CalStripType_t = CalStripType::kMIP) const;

  //InoFitterTrackHandle

  void SetFinderMomentum(double);
  double GetFinderMomentum() const;

  void SetMomentumdS(double);
  double GetMomentumdS() const;

  void SetMomentumRange(double);
  double GetMomentumRange() const;
  
  void SetMomentumCurve(double);
  double GetMomentumCurve() const;

  void SetEndMomentumCurve(double);
  double GetEndMomentumCurve() const;

  double GetEMCharge() const;
  void SetEMCharge(double);

  double GetVtxQPError() const;
  void SetVtxQPError(double);

  double GetVtxUError() const;
  void SetVtxUError(double);
  
  double GetVtxVError() const;
  void SetVtxVError(double);
  ////////////////////////////////
  double GetVtxXX() const;
  void SetVtxXX(double);
  
  double GetVtxYY() const;
  void SetVtxYY(double);  
  
  double GetVtxTX() const;
  void SetVtxTX(double);
  
  double GetVtxTY() const;
  void SetVtxTY(double);    
  ////////////////////////////////
  double GetVtxdU() const;
  void SetVtxdU(double);
  
  double GetVtxdV() const;
  void SetVtxdV(double);
  ////////////////////////////////
  double GetVtxdUError() const;
  void SetVtxdUError(double);
  
  double GetVtxdVError() const;
  void SetVtxdVError(double);  

  void SetBave(Double_t);
  Double_t GetBave() const;

  void SetEndQP(Double_t);
  void SetPlaneChi2(Int_t,Double_t);
  void SetPlaneQP(Int_t,Double_t);
  void SetEndUError(Double_t);
  void SetEndVError(Double_t);
  void SetEnddUError(Double_t);
  void SetEnddVError(Double_t);
  void SetEndQPError(Double_t);
  void SetNSwimFail(Int_t);

  Double_t GetEndQP() const;
  Float_t GetPlaneChi2(Int_t) const;
  Float_t GetPlaneQP(Int_t) const;
  Double_t GetEndUError() const;
  Double_t GetEndVError() const;
  Double_t GetEnddUError() const;
  Double_t GetEnddVError() const;
  Double_t GetEndQPError() const;
  Int_t GetNSwimFail() const;

  Double_t GetChi2() const;
  Double_t Getcval() const {return mclight;}
  void SetChi2(Double_t);
  void Setcval(Double_t c) {mclight= c;}

  Int_t GetNDOF() const;
  void SetNDOF(Int_t);

  Double_t GetRangeBiasedQP() const;
  void SetRangeBiasedQP(Double_t qp);

  Int_t GetNIterate() const;
  void SetNIterate(Int_t nit);

  vector<InoCluster*>ClustsInTrack;
  unsigned int GetClusterEntries() const {return ClustsInTrack.size();}
  
  void SetEndRPCmod(Int_t ivar);
  void SetVtxRPCmod(Int_t ivar);
  Int_t GetVtxRPCmod();
  Int_t GetEndRPCmod();
 
 private :
    
  double mFinderMomentum;  
  double mMomentumdS;
  double mMomentumRange;
  double mMomentumCurve;
  double mEndMomentumCurve;
  double mEMCharge;
  double mVtxQPError;
  int FCPC;
  int checkfcpc;
  Int_t mNDOF;
  
  double mVtxXX;
  double mVtxYY;  
  double mVtxTX;
  double mVtxTY;  
  double mVtxdU;
  double mVtxdV;
  
  double errth;
  double mVtxUError;
  double mVtxVError;
  double mVtxdUError;
  double mVtxdVError;

  Double_t mEndQP;
  Double_t mBave;
  Float_t mPlaneChi2[500];
  Float_t mPlaneQP[500];
  Double_t mEndUError;
  Double_t mEndVError;
  Double_t mEnddUError;
  Double_t mEnddVError;
  Double_t mEndQPError;
  Int_t mNSwimFail;
  Double_t mChi2;
  Double_t mclight;
  Double_t fQP_rangebiased;
  Int_t fNIterate;
  Int_t fType;
  
};

#endif                                              // INOTRACKCAND_H
