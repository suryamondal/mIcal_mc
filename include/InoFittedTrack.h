#ifndef CANDFITTRACK_H
#define CANDFITTRACK_H
//CandFitTrack
#include <iostream>
#include "InoTrack.h"

using namespace std;
#include <map>
#include "InoTrackCand.h"
//class InoFittedTrackHandle;

class InoFittedTrack // : public CandTrack
{
  //  friend class InoFittedTrackHandle;

 public:
  //  static InoFittedTrackHandle MakeCandidate(AlgHandle &ah,
  //                                                       CandContext &cx);
  InoFittedTrack();

  std::ostream& FormatToOStream(std::ostream& os,
                                        Option_t *option="") const;

  InoFittedTrack(InoTrack &ah);
  InoFittedTrack(const InoFittedTrack &rhs);
   ~InoFittedTrack();
   void CreateLocalHandle();
protected:
  //GMA open this  virtual InoFittedTrack *Dup() const;
  //  virtual Bool_t IsEquivalent(const TObject *rhs) const;

  Double_t fEMCharge;                     // in units of positron charge
  Double_t fChi2;
  Double_t fMomentumCurve;
  Double_t fMomentumRange;
  Bool_t fPass;                               // = 1 if successfully fit
  Double_t fBave;
  Int_t fNDOF;
  Double_t fCPUTime;              // time spent to create this candidate
  Double_t fVtxQPError;
  Int_t fNIterate;
  Double_t fVtxUError;
  Double_t fVtxVError;
  Double_t fVtxdUError;
  Double_t fVtxdVError;

  Double_t fEndQP;
  Double_t fEndQPError;
  Double_t fEndUError;
  Double_t fEndVError;
  Double_t fEnddUError;
  Double_t fEnddVError;
  Int_t fNSwimFail;

  mutable map<Int_t,Float_t> fPlaneQP;     // filtered q/p at this plane
  mutable map<Int_t,Float_t> fPlaneChi2; // filtered chi2 at this pln

  //  CandTrackHandle * fFinderTrack;

  InoTrackCand* fTrackCand;
  InoTrack*     fTrack;

};

#endif                                                 // CANDFITTRACK_H
