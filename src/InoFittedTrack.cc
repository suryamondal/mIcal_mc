#include "InoFittedTrack.h"
//#include "InoFittedTrackHandle.h"
#include "TString.h"

InoFittedTrack::InoFittedTrack() :
  fEMCharge(0.),
  fChi2(0.),
  fMomentumCurve(0.),
  fMomentumRange(0.),
  fPass(0),
  fBave(0),
  fNDOF(0),
  fCPUTime(0.),
  fVtxQPError(0.),
  fNIterate(0),
  fVtxUError(0),
  fVtxVError(0),
  fVtxdUError(0),
  fVtxdVError(0),
  fEndQP(0.),
  fEndQPError(0.),
  fEndUError(0.),
  fEndVError(0.),
  fEnddUError(0.),
  fEnddVError(0.),
  fNSwimFail(0),
  fTrackCand(0),
  fTrack(0)
{
  cout<< "Begin InoFittedTrack::InoFittedTrack() ctor: " << endl;
    //      << "UidInt = " << GetUidInt()
    //      << ", ArchUidInt " << GetArchUidInt() << endl
    //      << "No. of links = " << GetNLinks() << endl
    //      << "End InoFittedTrack::InoFittedTrack() ctor." << endl;
}

/*
//______________________________________________________________________
InoFittedTrack::InoFittedTrack(AlgHandle &ah) :
  CandTrack(ah),     // Should be the next class up on inheritance chain
  fEMCharge(0.),
  fChi2(0.),
  fMomentumCurve(0.),
  fMomentumRange(0.),
  fPass(0),
  fBave(0),
  fNDOF(0),
  fCPUTime(0.),
  fVtxQPError(0.),
  fNIterate(0),
  fVtxUError(0),
  fVtxVError(0),
  fVtxdUError(0),
  fVtxdVError(0),
  fEndQP(0.),
  fEndQPError(0.),
  fEndUError(0.),
  fEndVError(0.),
  fEnddUError(0.),
  fEnddVError(0.),
  fNSwimFail(0),  
  fTrackcand(0),
  fTrack(0)
{

// The sole purpose of this constructor is to transmit the AlgHandle
// up the inheritance chain to CandBase without having to invoke the
// full constructor of an intermediate Candidate type which the highest
// level Candidate might inherit from.  One only wants to create the
// LocalHandle and invoke the RunAlg() method in the lowest level class.
}

//______________________________________________________________________
InoFittedTrack::InoFittedTrack(AlgHandle &ah, CandHandle &ch,
                                                      CandContext &cx) :
  CandTrack(ah),     // Should be the next class up on inheritance chain
  fEMCharge(0.),
  fChi2(0.),
  fMomentumCurve(0.),
  fMomentumRange(0.),
  fPass(0),
  fBave(0),
  fNDOF(0),
  fCPUTime(0.),
  fVtxQPError(0.),
  fNIterate(0),
  fVtxUError(0),
  fVtxVError(0),
  fVtxdUError(0),
  fVtxdVError(0),
  fEndQP(0.),
  fEndQPError(0.),
  fEndUError(0.),
  fEndVError(0.),
  fEnddUError(0.),
  fEnddVError(0.),
  fNSwimFail(0),
  fTrackCand(0),
  fTrack(0)
{
  CreateLocalHandle();
  MSG("Cand", Msg::kDebug)
    << "Begin InoFittedTrack::InoFittedTrack(AlgHandle &, CandHandle &, "
    << "CandContext &) ctor: " << endl
    << "UidInt = " << GetUidInt()
    << ", ArchUidInt " << GetArchUidInt() << endl
    << "No. of links = " << GetNLinks() << endl
    << "End InoFittedTrack::InoFittedTrack(AlgHandle &, CandHandle &, "
    << "CandContext &) ctor." << endl;
  
  // Run Algorithm to construct Candidate
  {                                                   // Start of scope.
    InoFittedTrackHandle csh(this);            // csh will go out of scope
    ch = csh;                                       // after setting ch.
  }                                                     // End of scope.
  ah.RunAlg(ch, cx);
}

*/

//______________________________________________________________________
InoFittedTrack::InoFittedTrack(const InoFittedTrack &rhs) :
  //  CandTrack(rhs),    // Should be the next class up on inheritance chain
  fEMCharge(rhs.fEMCharge),
  fChi2(rhs.fChi2),
  fMomentumCurve(rhs.fMomentumCurve),
  fMomentumRange(rhs.fMomentumRange),
  fPass(rhs.fPass),
  fBave(rhs.fBave),
  fNDOF(rhs.fNDOF),
  fCPUTime(rhs.fCPUTime),
  fVtxQPError(rhs.fVtxQPError),
  fNIterate(rhs.fNIterate),
  fVtxUError(rhs.fVtxUError),
  fVtxVError(rhs.fVtxVError),
  fVtxdUError(rhs.fVtxdUError),
  fVtxdVError(rhs.fVtxdVError),
  fEndQP(rhs.fEndQP),
  fEndQPError(rhs.fEndQPError),
  fEndUError(rhs.fEndUError),
  fEndVError(rhs.fEndVError),
  fEnddUError(rhs.fEnddUError),
  fEnddVError(rhs.fEnddVError),
  fNSwimFail(rhs.fNSwimFail),
  fTrackCand(rhs.fTrackCand),
  fTrack(rhs.fTrack)
{
  map<Int_t,Float_t>::iterator fPlaneChi2Iter;
  map<Int_t,Float_t>::iterator fPlaneQPIter;

  for (fPlaneChi2Iter = rhs.fPlaneChi2.begin();
       fPlaneChi2Iter!=rhs.fPlaneChi2.end(); fPlaneChi2Iter++) {
    fPlaneChi2[fPlaneChi2Iter->first] = fPlaneChi2Iter->second;
  }
  for (fPlaneQPIter = rhs.fPlaneQP.begin();
       fPlaneQPIter!=rhs.fPlaneQP.end(); fPlaneQPIter++) {
    fPlaneQP[fPlaneQPIter->first] = fPlaneQPIter->second;
  }
  
  //CreateLocalHandle(); // Moved to Dup function following copy-ctor call
  cout<< "Begin InoFittedTrack::InoFittedTrack(const InoFittedTrack &rhs) ctor:"
    //      << endl << "UidInt = " << GetUidInt()
    //      << ", ArchUidInt " << GetArchUidInt() << endl
    //      << "No. of links = " << GetNLinks() << endl
    //      << "End InoFittedTrack::InoFittedTrack(const InoFittedTrack &rhs) ctor."
                                                                << endl;
}

//______________________________________________________________________
InoFittedTrack::~InoFittedTrack() {
  cout << "Begin InoFittedTrack::~InoFittedTrack() dtor: " << endl;
  //       << "UidInt = " << GetUidInt()
  //       << ", ArchUidInt " << GetArchUidInt() << endl
  //       << "No. of links = " << GetNLinks() << endl
  //       << "End InoFittedTrack::~InoFittedTrack() dtor." << endl;
}

//______________________________________________________________________
void InoFittedTrack::CreateLocalHandle() {
  //GMA open it  SetLocalHandle(new InoFittedTrackHandle(this));
}

/*
//______________________________________________________________________
InoFittedTrack *InoFittedTrack::Dup() const
{

// Base copy ctor dups owned pointers, but defers copying Daughter List.
// Daughter List copy is made in the derived class Dup() function.
// This is because base class copy constructor hasn't yet created
// fLocalHandle with a CandHandle* of the full derived type.
  InoFittedTrack *cb = new InoFittedTrack(*this);     // Copy-ctor dups ptrs
  cb->CreateLocalHandle();   // Initializes fLocalHandle after copy-ctor
  TIter iterdau = GetDaughterIterator();
  CandHandle *dau;
  while ((dau=(CandHandle *) iterdau())) cb->AddDaughterLink(*dau);
  return cb;
}
*/

/*
//______________________________________________________________________
Bool_t InoFittedTrack::IsEquivalent(const TObject *rhs) const
{
  Bool_t result = true;
  //  if (!CandTrack::IsEquivalent(rhs)) result = false;  // superclass test
  TestDisplayCandBanner("InoFittedTrack");
  const InoFittedTrack* rCnd = dynamic_cast<const InoFittedTrack*>(rhs);
  if (rCnd == NULL) return false;

  result = TestEquality("fEMCharge",         this->fEMCharge,
                                             rCnd->fEMCharge) && result;
  result = TestEquality("fChi2",                 this->fChi2,
                                                 rCnd->fChi2) && result;
  result = TestEquality("fMomentumCurve",
                                        this->fMomentumCurve,
                                        rCnd->fMomentumCurve) && result;
  result = TestEquality("fPass",        this->fPass,
                                        rCnd->fPass)          && result;
  result = TestEquality("fBave",        this->fBave,
                                        rCnd->fBave)          && result;
  result = TestEquality("fNDOF",         this->fNDOF, 
                                         rCnd->fNDOF)         && result;
  result = TestEquality("fCPUTime",      this->fCPUTime, 
                                         rCnd->fCPUTime)      && result;
  result = TestEquality("fNIterate",
                            this->fNIterate,
                            rCnd->fNIterate)                  && result;
  result = TestEquality("fVtxQPError",   this->fVtxQPError, 
                                         rCnd->fVtxQPError)   && result; 
  result = TestEquality("fVtxUError",   this->fVtxUError, 
                                         rCnd->fVtxUError)   && result; 
  result = TestEquality("fVtxVError",   this->fVtxVError, 
                                         rCnd->fVtxVError)   && result; 
  result = TestEquality("fVtxdUError",   this->fVtxdUError, 
                                         rCnd->fVtxdUError)   && result; 
  result = TestEquality("fVtxdVError",   this->fVtxdVError, 
                                         rCnd->fVtxdVError)   && result; 
  result = TestEquality("fFinderTrack",   this->fFinderTrack, 
                                         rCnd->fFinderTrack)   && result;
  result = TestEquality("fEndQP",
                        this->fEndQP,
                        rCnd->fEndQP)                     && result;
  result = TestEquality("fPlaneChi2",    this->fPlaneChi2, 
                                         rCnd->fPlaneChi2)    && result;
  result = TestEquality("fPlaneQP",      this->fPlaneQP, 
                                         rCnd->fPlaneQP)      && result;
  result = TestEquality("fEndUError",    this->fEndUError, 
                                         rCnd->fEndUError)    && result;
  result = TestEquality("fEndVError",    this->fEndVError, 
                                         rCnd->fEndVError)    && result;
  result = TestEquality("fEnddUError",   this->fEnddUError, 
                                         rCnd->fEnddUError)   && result;
  result = TestEquality("fEnddVError",   this->fEnddVError, 
                                         rCnd->fEnddVError)   && result;
  result = TestEquality("fEndQPError",   this->fEndQPError, 
                                         rCnd->fEndQPError)   && result;
  result = TestEquality("fNSwimFail",    this->fNSwimFail, 
                                         rCnd->fNSwimFail)    && result;

  return result;
}
*/

//______________________________________________________________________
//InoFittedTrackHandle InoFittedTrack::MakeCandidate(AlgHandle &ah,
//                                                        CandContext &cx)
//{
//  InoFittedTrackHandle csh;
//  new InoFittedTrack(ah, csh, cx);        // csh owns the new InoFittedTrack
//  return csh;
//}

//______________________________________________________________________
std::ostream& InoFittedTrack::FormatToOStream(std::ostream& os,
					      Option_t *option) const {
  //  CandTrack::FormatToOStream(os,option);
  
  TString opt(option);
  if (!opt.Contains("v0")) { // v0 means suppress the data values
    //    const TString& indent = GetIndentString();
    
    //    os << indent << GetDataIndent()
    os   << "EMCharge " << fEMCharge
	 << " Chi2 " << fChi2
	 << " MomentumCurve " << fMomentumCurve
	 << " NIterate " << fNIterate
	 << "NDOF " << fNDOF 
	 << " " << (fPass?"Pass":"Fail")
	 << " CPUTime " << fCPUTime
	 << " VtxUError " << fVtxUError
	 << " VtxVError " << fVtxVError
	 << " VtxdUError " << fVtxdUError
	 << " VtxdVError " << fVtxdVError
	 << endl;
  }
  return os;
  
}
