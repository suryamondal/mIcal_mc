#include "TMath.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include "TVector3.h"

//#include "RecoBase/CandSliceHandle.h"
//#include "CandStrip.h"
#include "vect_manager.h"
#include "InoVertex.h"
#include "InoTrackCand.h"
#include "InoShowerCand.h"

//GMA need start and end points for both forward and backward fit.
//    Thus, need four inoVertex paramter for this, now we have only two
InoTrackCand::InoTrackCand(const InoTrackCand &rhs) {
  //GMA Initialise all parameters here 

  fInShower=   rhs.fInShower;    
  fUPos=       rhs.fUPos;      
  fVPos=       rhs.fVPos;      
  fdS=         rhs.fdS;                            
  fRange=      rhs.fRange;                        
  fXPosError=  rhs.fXPosError; 
  fYPosError=  rhs.fYPosError; 
  fTime[0]=    rhs.fTime[0];  
  fTime[1]=    rhs.fTime[1];  

  fTrack =                 rhs.fTrack;                 
  fVtxTrace=		   rhs.fVtxTrace;              
  fVtxTraceZ=		   rhs.fVtxTraceZ;             
  fEndTrace=		   rhs.fEndTrace;              
  fEndTraceZ=		   rhs.fEndTraceZ;             
  fVtxDistToEdge=	   rhs.fVtxDistToEdge;         
  fEndDistToEdge=	   rhs.fEndDistToEdge;         
  fVtxnActiveUpstream=	   rhs.fVtxnActiveUpstream;    
  fEndnActiveDownstream=   rhs.fEndnActiveDownstream;  
  			                               
  fNTrackStrip=            rhs.fNTrackStrip;             
  fNTrackDigit=            rhs.fNTrackDigit;             
  fNTimeFitDigit=          rhs.fNTimeFitDigit;           
  fTimeFitChi2=		   rhs.fTimeFitChi2;           
  fTimeForwardFitRMS=	   rhs.fTimeForwardFitRMS;     
  fTimeForwardFitNDOF=	   rhs.fTimeForwardFitNDOF;    
  fTimeBackwardFitRMS=	   rhs.fTimeBackwardFitRMS;    
  fTimeBackwardFitNDOF=	   rhs.fTimeBackwardFitNDOF;   
  			                               
  fTimeSlope=		   rhs.fTimeSlope;             
  fTimeOffset=		   rhs.fTimeOffset;           
  fVertex =                rhs.fVertex;
  fTerm =                  rhs.fTerm;
  ClustsInTrack=           rhs.ClustsInTrack;
}

InoTrackCand::InoTrackCand( ) //  :  fTrack(trk)
{
  //GMA Initialise all parameters here  
  /*
  fInShower=   0;
  fUPos=       0;
  fVPos=       0;
  fdS=         0;
  fRange=      0;
  fXPosError=  0;
  fYPosError=  0;
  fTime[0]=    0;
  fTime[1]=    0;
  */ 

  fTrack = 0;
  fVtxTrace=-999;
  fVtxTraceZ=-999;
  fEndTrace=-999;
  fEndTraceZ=-999;
  fVtxDistToEdge=-999;
  fEndDistToEdge=-999;
  fVtxnActiveUpstream=-999;
  fEndnActiveDownstream=-999;
  
  fNTrackStrip=-999;            
  fNTrackDigit=-999;            
  fNTimeFitDigit=-999;          
  fTimeFitChi2=-999;
  fTimeForwardFitRMS=-999;
  fTimeForwardFitNDOF=-999;
  fTimeBackwardFitRMS=-999;
  fTimeBackwardFitNDOF=-999;
  
  fTimeSlope=-999;
  fTimeOffset=-999;
  
  ClustsInTrack.clear();
}

InoTrackCand::InoTrackCand(InoTrack *trk, bool forward) :  fTrack(trk) {
  //  InoTrackCand();
  //GMA Initialise all parameters here  
  fTrack = trk;
  
  double dxdz=trk->GetBegXDir();
  double dydz=trk->GetBegYDir();
  
  double phi = atan2(dydz, dxdz);
  double zdir= 1./pow (1+dxdz*dxdz+dydz*dydz,0.5);

  double xdir = sin(acos(zdir))*cos(phi); 
  double ydir = sin(acos(zdir))*sin(phi);

  TVector3 tmp3v(xdir, ydir, zdir);
  
  fVertex = new InoVertex();
  fVertex->SetU(trk->GetBegXPos());
  fVertex->SetV(trk->GetBegYPos());
  fVertex->SetZ(trk->GetBegZPos());
  fVertex->SetPlane(trk->GetBegZPlane());

  fVertex->SetDirCosine(tmp3v);

  dxdz=trk->GetEndXDir();
  dydz=trk->GetEndYDir();
  
  phi = atan2(dydz, dxdz);
  zdir= 1./pow(1+dxdz*dxdz+dydz*dydz,0.5);

  xdir = sin(acos(zdir))*cos(phi); 
  ydir = sin(acos(zdir))*sin(phi);

  TVector3  tmp3vm(xdir, ydir, zdir);
  
  fTerm = new InoVertex();
  fTerm->SetU(trk->GetEndXPos());
  fTerm->SetV(trk->GetEndYPos());
  fTerm->SetZ(trk->GetEndZPos());
  fTerm->SetPlane(trk->GetEndZPlane());

  fTerm->SetDirCosine(tmp3vm);

  mFinderMomentum = pow(pow(trk->GetEndXPos()-trk->GetBegXPos(), 2.) + 
			pow(trk->GetEndYPos()-trk->GetBegYPos(), 2.) + 
			pow(trk->GetEndZPos()-trk->GetBegZPos(), 2.) , 0.5);
  
  ClustsInTrack.clear();
  int ncls = trk->GetEntries();
  for (unsigned ij=0; ij<trk->GetEntries(); ij++) {
    if (((forward && ij<=ncls/2.) || (!forward && ij>=ncls/2.)) && (!trk->ClustsInTrack[ij]->GetStraight())) continue;
    
    ClustsInTrack.push_back(trk->ClustsInTrack[ij]);
  }
}

//______________________________________________________________________
InoTrackCand::~InoTrackCand()
{
  //  ClustsInTrack.clear();
}

//InoTrackCand *InoTrackCand::DupHandle() const
//{
//   return (new InoTrackCand(*this));
//}


//______________________________________________________________________
void InoTrackCand::Trace(const char *c) const
{
  cout <<" InorackCand : **********Begin InoTrackCand::Trace(\"" << c << "\")" << endl
       << "Information from InoTrackCand's CandHandle: " << endl;
  //  CandHandle::Trace(c);
  cout <<" InoTrackCand : ********End InoTrackCand::Trace(\"" << c << "\")" << endl;
}


//NavKey InoTrackCand::KeyFromSlice(const InoTrackCand *reco) {
//  if (reco->GetCandSlice()) {
//    return static_cast<Int_t>(reco->GetCandSlice()->GetUidInt());
//  }
//  return 0;
//}

void InoTrackCand::SetTrackPointXError(Int_t plane, Float_t tpos) {
  fXPosError[plane] = tpos;
}

void InoTrackCand::SetTrackPointYError(Int_t plane, Float_t tpos) {
  fYPosError[plane] = tpos;
}

void InoTrackCand::SetU(Int_t plane, Float_t tpos)  {
  fUPos[plane] = tpos;
}

void InoTrackCand::SetV(Int_t plane, Float_t tpos)  {
  fVPos[plane] = tpos;
}

void InoTrackCand::SetdS(Int_t plane, Float_t ds) {
  fdS[plane] = ds;
}

void InoTrackCand::SetRange(Int_t plane, Float_t range) {
  fRange[plane] = range;
}

void InoTrackCand::Set2dS(Int_t plane, Float_t ds) {
  f2dS[plane] = ds;
}

void InoTrackCand::Set2Range(Int_t plane, Float_t range) {
  f2Range[plane] = range;
}

void InoTrackCand::SetT(Int_t plane, Double_t time)
{
  /*
  switch (stripend_t) {
  case StripEnd::kNegative:
    track->fTime[0][plane] = time;
    break;
  case StripEnd::kPositive:
    track->fTime[1][plane] = time;
    break;
  default:
  track->fTime[0][plane] = time;
  break;
  }
  */
  fTime[0][plane] = time;
}

void InoTrackCand::ClearMaps()
{
  fUPos.clear();
  fVPos.clear();
  fTime[0].clear();
  fTime[1].clear();
  fdS.clear();
  fRange.clear();
}

void InoTrackCand::ClearUVT()
{
  ClearMaps();
}

/*
Bool_t InoTrackCand::IsXPosValid(Int_t plane) const
{
if (UPos.count(plane) && fVPos.count(plane)) {
    return kTRUE;
  }
  return kFALSE;
}
*/
Float_t InoTrackCand::GetTrackPointXError(Int_t plane) const {
  if(fXPosError.count(plane)){
    return fXPosError[plane];
  }
  return -99999.;
}

Float_t InoTrackCand::GetTrackPointYError(Int_t plane) const {
  if(fYPosError.count(plane)) {
    return fYPosError[plane];
  }
  return -99999.;
}

Float_t InoTrackCand::GetU(Int_t plane) const {
  if (fUPos.count(plane)) {
    return fUPos[plane];
  }
  return -99999.;
}

Float_t InoTrackCand::GetV(Int_t plane) const {
  if (fVPos.count(plane)) {
    return fVPos[plane];
  }
  return -99999.;
}

Float_t InoTrackCand::GetZ(Int_t plane) const {
  InoHit_Manager *pinohit = InoHit_Manager::APointer;
  for (unsigned i=0; i<pinohit->InoHit_list.size() ; i++) {
    if (pinohit->InoHit_list[i]->GetZPlane()==plane) {
      return pinohit->InoHit_list[i]->GetZPos();
    }
  }
  return -1.;
}

/*
Double_t InoTrackCand::GetT(Int_t plane) const // ,StripEnd::StripEnd_t stripend_t) const {
  if (fTime[0].count(plane)) return fTime[0][plane];

  const CandTrack *track = dynamic_cast<const CandTrack*>(GetCandBase());
  if (stripend_t==StripEnd::kNegative) {
    if ((track->fTime[0]).count(plane)) {
      return track->fTime[0][plane];
    }
    return -99999.;
  }
  if (stripend_t==StripEnd::kPositive) {
    if ((track->fTime[1]).count(plane)) {
      return track->fTime[1][plane];
    }
    return -99999.;
  }
  if ((track->fTime[0]).count(plane) && (track->fTime[1]).count(plane)) {
    return min(track->fTime[0][plane],track->fTime[1][plane]);
  }
  else if ((track->fTime[0]).count(plane)) {
    return track->fTime[0][plane];
  }
  else if ((track->fTime[1]).count(plane)) {
    return track->fTime[1][plane];
  }

  return -99999.;
}

Double_t InoTrackCand::GetT(StripEnd::StripEnd_t stripend_t,Int_t plane) const {
  return GetT(plane,stripend_t);
}
*/

Double_t InoTrackCand::GetT(Int_t plane) const {
  if (fTime[0].count(plane)) return fTime[0][plane];
  return -99999;
}


//GMA add +- half lenght before(after) first(last) layer
Float_t InoTrackCand::GetdS(Int_t plane) const {
  if (plane<0) {plane = GetVtxPlane(); }
  if (fdS.count(plane)) { return fdS[plane]; }
  return -1.;
}

Float_t InoTrackCand::GetRange(Int_t plane) const {
  if (plane<0) { plane = GetVtxPlane();  }
  if (fRange.count(plane)) {return fRange[plane]; }
  return -1.;
}

Float_t InoTrackCand::Get2dS(Int_t plane) const {
  if (plane<0) {plane = GetVtxPlane(); }
  if (f2dS.count(plane)) { return f2dS[plane]; }
  return -1.;
}

Float_t InoTrackCand::Get2Range(Int_t plane) const {
  if (plane<0) { plane = GetVtxPlane();  }
  if (f2Range.count(plane)) {return f2Range[plane]; }
  return -1.;
}


Double_t InoTrackCand::GetScore() const {
  //GMA set plane in three categories
  // 0 : only X-axis hit
  // 1 : only Y-axis hit
  // 2 : Both X and Y-axis hit 

  Double_t score = 0.;
  score = (Double_t)(GetNDaughters());
  Int_t dbegpln = GetBegPlane(0)-GetBegPlane(1);
  Int_t dendpln = GetEndPlane(0)-GetEndPlane(1);
  Double_t dbegpln2 = (Double_t)(dbegpln*dbegpln);
  Double_t dendpln2 = (Double_t)(dendpln*dendpln);
  score -= (dbegpln2+dendpln2);
  return score;
}

/*
Int_t InoTrackCand::GetNTrackPlane(PlaneView::PlaneView_t planeview_t) const
{
  CandStripHandleItr stripItr(GetDaughterIterator());
  CandStripHandleKeyFunc *stripKf = stripItr.CreateKeyFunc();
  stripKf->SetFun(CandStripHandle::KeyFromPlane);
  stripItr.GetSet()->AdoptSortKeyFunc(stripKf);
  stripKf = 0;

  Int_t plane=0;
  Int_t oldplane=0;
  CandStripHandle *strip;
  while ((strip = dynamic_cast<CandStripHandle*>(stripItr()))) {
    if (IsInShower(strip)<=1) {
      PlaneView::PlaneView_t planeview = strip->GetPlaneView();
      if (planeview!=PlaneView::kA && planeview!=PlaneView::kB &&
          (planeview_t==PlaneView::kUnknown || planeview==planeview_t)) {
        if (!plane || strip->GetPlane()!=oldplane) {
          plane++;
        }
        oldplane = strip->GetPlane();
      }
    }
  }
  return plane;
}
*/

Double_t InoTrackCand::GetMomentum() const {
  return fMomentum;
}

void InoTrackCand::SetMomentum(Double_t momentum) {
  fMomentum = momentum;
}

int InoTrackCand::GetFCPC() const {
  return FCPC;
}

void InoTrackCand::SetFCPC(int val) {
  FCPC = val;
}


//------------------------------------------


//------------------------------------------
Double_t InoTrackCand::GetTheta() const {
  return fTheta;
}

void InoTrackCand::SetTheta(Double_t theta) {
  fTheta = theta;
}
//------------------------------------------
void InoTrackCand::SetThErr(Double_t aerrth) {
  fErrTh = aerrth;
}
Double_t InoTrackCand::GetThErr() const {
  return fErrTh;
}
//------------------------------------------
Double_t InoTrackCand::GetPhi() const {
  return fPhi;
}

void InoTrackCand::SetPhi(Double_t phi) {
  fPhi = phi;
}
//------------------------------------------
void InoTrackCand::SetPhErr(Double_t errph) {
  fErrPh = errph;
}
Double_t InoTrackCand::GetPhErr() const {
  return fErrPh;
}
//------------------------------------------
Double_t InoTrackCand::GetdSExtra() const {
  return fdSExtra;
}

void InoTrackCand::SetdSExtra(Double_t dsextra) {
  fdSExtra = dsextra;
}

Double_t InoTrackCand::GetRangeExtra() const {
  return fRangeExtra;
}

void InoTrackCand::SetRangeExtra(Double_t rangeextra) {
  fRangeExtra = rangeextra;
}

/*
Int_t InoTrackCand::IsInShower(CandStripHandle *striphandle) const {

  const CandTrack *track = dynamic_cast<const CandTrack*>(GetCandBase());
  const CandStrip *strip = dynamic_cast<const CandStrip*>(striphandle->GetCandBase());
  if (track->fInShower.count(strip)>0) {
    return track->fInShower[strip];
  }
  return 0;
}

void InoTrackCand::SetInShower(CandStripHandle *striphandle, Int_t ival) {
  CandTrack *track = dynamic_cast<CandTrack*>(GetOwnedCandBase());
  const CandStrip *strip = dynamic_cast<const CandStrip*>(striphandle->GetCandBase());
  track->fInShower[strip] = ival;
}
*/

void InoTrackCand::SetVtxTrace(Double_t dvar) {
  fVtxTrace = dvar;
}

void InoTrackCand::SetVtxTraceZ(Double_t dvar) {
  fVtxTraceZ = dvar;
}

void InoTrackCand::SetVtxnActiveUpstream(int ivar) {
  fVtxnActiveUpstream = ivar;
}

void InoTrackCand::SetEndTrace(Double_t dvar) {
  fEndTrace = dvar;
}

void InoTrackCand::SetEndTraceZ(Double_t dvar) {
  fEndTraceZ = dvar;
}

void InoTrackCand::SetEndnActiveDownstream(int ivar) {
  fEndnActiveDownstream = ivar;
}

void InoTrackCand::SetVtxDistToEdge(Double_t dvar) {
  fVtxDistToEdge = dvar;
}

void InoTrackCand::SetEndDistToEdge(Double_t dvar) {
  fEndDistToEdge = dvar;
}

Double_t InoTrackCand::GetVtxTrace() const {
  return fVtxTrace;
}

Double_t InoTrackCand::GetVtxTraceZ() const {
  return fVtxTraceZ;
}

Int_t InoTrackCand::GetVtxnActiveUpstream() const {
  return fVtxnActiveUpstream;
}

Double_t InoTrackCand::GetEndTrace() const {
  return fEndTrace;
}

Double_t InoTrackCand::GetEndTraceZ() const {
  return fEndTraceZ;
}

Int_t InoTrackCand::GetEndnActiveDownstream() const {
  return fEndnActiveDownstream;
}

Double_t InoTrackCand::GetVtxDistToEdge() const {
  return fVtxDistToEdge;
}

Double_t InoTrackCand::GetEndDistToEdge() const {
  return fEndDistToEdge;
}

/* 
//_____________________________________________________________________
Bool_t InoTrackCand::BelongsWithTrack(InoTrackCand * trk, 
				      Double_t tolTPos2, Double_t tolZPos, Double_t tolTime){
  
  if(!trk)return false;
  
  Int_t pSMPlaneLast = ac.GetInt("SMPlaneLast");
  Int_t pSMPlaneFirst = ac.GetInt("SMPlaneFirst");
  VldContext vldc = *vldcptr;
  Double_t zGapSM=0;
  UgliGeomHandle ugh(vldc); 
  if(vldc.GetDetector() == Detector::kFar){  
    // calculate Z gap between SM for later use
    PlexPlaneId scintlastid(vldc.GetDetector(),pSMPlaneLast,kFALSE);
    PlexPlaneId scintfirstid(vldc.GetDetector(),pSMPlaneFirst,kFALSE);
    UgliScintPlnHandle scintlast = ugh.GetScintPlnHandle(scintlastid);
    UgliScintPlnHandle scintfirst = ugh.GetScintPlnHandle(scintfirstid);
    
    if (scintlast.IsValid() && scintfirst.IsValid()) {  
      zGapSM=scintfirst.GetZ0()-scintlast.GetZ0()-0.0594;
    } 
  } 
  Double_t dz =this->GetVtxZ()-trk->GetVtxZ();
  Double_t du= this->GetVtxU()-trk->GetVtxU();
  Double_t dv= this->GetVtxV()-trk->GetVtxV();  
  Double_t dt= this->GetVtxT()-trk->GetVtxT();
  
  //compensate dz for SM gap if necessary
  if (vldcptr->GetDetector()==Detector::kFar &&
      this->GetVtxPlane()<=pSMPlaneLast &&
      trk->GetVtxPlane()>=pSMPlaneFirst) { 
    dz+=zGapSM;
  } else if (vldcptr->GetDetector()==Detector::kFar &&
	     this->GetVtxPlane()>=pSMPlaneLast &&
	     trk->GetVtxPlane()<=pSMPlaneFirst) {
    dz-=zGapSM;
  }
  if (du*du+dv*dv<tolTPos2 &&
      fabs(dz)<tolZPos && 
      fabs(dt)<tolTime) return true;
  
  return false;
}
//_____________________________________________________________________
Bool_t InoTrackCand::BelongsWithShower(CandShowerHandle * shw, 
				       Double_t tolTPos2, Double_t tolZPos, Double_t tolTime){
  if(!shw)return false;
  
  Int_t pSMPlaneLast = ac.GetInt("SMPlaneLast");
  Int_t pSMPlaneFirst = ac.GetInt("SMPlaneFirst");
  Double_t zGapSM=0;
  VldContext vldc = *vldcptr;
  UgliGeomHandle ugh(vldc); 
  // calculate Z gap between SM for later use
  if(vldc.GetDetector() == Detector::kFar){     
    PlexPlaneId scintlastid(vldc.GetDetector(),pSMPlaneLast,kFALSE);
    PlexPlaneId scintfirstid(vldc.GetDetector(),pSMPlaneFirst,kFALSE);
    UgliScintPlnHandle scintlast = ugh.GetScintPlnHandle(scintlastid);  
    UgliScintPlnHandle scintfirst = ugh.GetScintPlnHandle(scintfirstid);
    if (scintlast.IsValid() && scintfirst.IsValid()) {  
      zGapSM=scintfirst.GetZ0()-scintlast.GetZ0()-0.0594;
    } 
  }
  
  Float_t tolTPos=TMath::Sqrt(tolTPos2);
  Int_t fwdTrkPlane=min(this->GetVtxPlane(),this->GetEndPlane());
  Int_t bckTrkPlane=max(this->GetVtxPlane(),this->GetEndPlane());
  Int_t fwdShwPlane=min(shw->GetVtxPlane(),shw->GetEndPlane());
  Int_t bckShwPlane=max(shw->GetVtxPlane(),shw->GetEndPlane());
  Double_t dt= this->GetVtxT()-shw->GetVtxT();
  if(this->IsTPosValid(shw->GetVtxPlane())){
    dt = shw->GetVtxT()-this->GetT(shw->GetVtxPlane());
  }
  if( bckShwPlane>=fwdTrkPlane && fwdShwPlane<=bckTrkPlane && fabs(dt)<tolTime){   // trk overlaps in Z and Time
    MSG("RecoBase",Msg::kDebug) << " trk and shower overlap in Z and time " << endl;
    Bool_t matchu=false;
    Bool_t matchv=false;
    for(Int_t iloop=0;iloop<2;iloop++){
      
      // if on second loop through and have at least one matching view, loosen tolerance
      if(iloop >0 && (matchu || matchv))tolTPos=tolTPos*2;
      // if on second loop and no matches so far, give up
      if(iloop >0 && !matchu && !matchv) break;
      
      for (Int_t iplane=fwdShwPlane; iplane<=bckShwPlane;iplane++){  // look for trk passing through shower
        Double_t trkU = 0, trkV = 0;
        //if(matchu || matchv)tolTPos = tolTPos*2;
        PlexPlaneId plnid(GetVldContext()->GetDetector(),iplane,false);
        if(plnid.GetPlaneView()==PlaneView::kU){
          trkU=this->GetU(iplane);
          if(!this->IsTPosValid(iplane)){
            if(abs(iplane-this->GetVtxPlane())<abs(iplane-this->GetEndPlane())){
              UgliScintPlnHandle scintpln = ugh.GetScintPlnHandle(plnid);
              Double_t dz=scintpln.GetZ0()-this->GetVtxZ();
              trkU = (this->GetVtxU() + dz*this->GetVtxDirCosU()/this->GetVtxDirCosZ());
            } else {
              UgliScintPlnHandle scintpln = ugh.GetScintPlnHandle(plnid);
              Double_t dz=scintpln.GetZ0()-this->GetEndZ();
              trkU = (this->GetEndU() + dz*this->GetEndDirCosU()/this->GetEndDirCosZ());
            }
          }
          if(shw->GetNStrips(iplane)>0 && 
             trkU>=shw->GetMinU(iplane)-tolTPos && 
             trkU<=shw->GetMaxU(iplane)+tolTPos) {
            MSG("RecoBase", Msg::kDebug)<< " u match! (details follow)" << endl;
            matchu=true;  
          }
        } else if(plnid.GetPlaneView()==PlaneView::kV) {
          trkV=this->GetV(iplane);
          if(!this->IsTPosValid(iplane)){
            if(abs(iplane-this->GetVtxPlane())<abs(iplane-this->GetEndPlane())){
              UgliScintPlnHandle scintpln = ugh.GetScintPlnHandle(plnid);
              Double_t dz=scintpln.GetZ0()-this->GetVtxZ();
              trkV = (this->GetVtxV() + dz*this->GetDirCosV()/this->GetDirCosZ());
            } else {
              UgliScintPlnHandle scintpln = ugh.GetScintPlnHandle(plnid);
              Double_t dz=scintpln.GetZ0()-this->GetEndZ();
              trkV= (this->GetEndV() + dz*this->GetEndDirCosV()/this->GetEndDirCosZ());
            }
          }

          if(shw->GetNStrips(iplane)>0 && 
             trkV>=shw->GetMinV(iplane)-tolTPos && 
             trkV<=shw->GetMaxV(iplane)+tolTPos) {
            MSG("RecoBase", Msg::kDebug)<< " v match! ((details follow)" << endl;
            matchv=true;
          }
        }

        MSG("RecoBase",Msg::kDebug) << " plane  " << iplane 
                                    << " trk u " << trkU 
                                    << " min/max U " 
                                    << shw->GetMinU(iplane) << "/" 
                                    << shw->GetMaxU(iplane)
                                    << " trk v " << trkV 
                                    << " min/max V " 
                                    << shw->GetMinV(iplane) << "/" 
                                    << shw->GetMaxV(iplane) << endl; 
	
        if(matchu && matchv) {
          MSG("RecoBase", Msg::kDebug)<< " 3D  match! " << endl;
          return true;
        }
      }
    }
  }
  
  if( shw->GetVtxZ()>min(this->GetVtxZ(),this->GetEndZ())-tolZPos  &&
      shw->GetVtxZ()<max(this->GetVtxZ(),this->GetEndZ())+tolZPos) {
    Double_t du= this->GetVtxU()-shw->GetVtxU();
    Double_t dv= this->GetVtxV()-shw->GetVtxV();        
    dt= this->GetVtxT()-shw->GetVtxT();
    if(this->IsTPosValid(shw->GetVtxPlane())){
      du = shw->GetVtxU()-this->GetU(shw->GetVtxPlane());
      dv = shw->GetVtxV()-this->GetV(shw->GetVtxPlane());
      dt = shw->GetVtxT()-this->GetT(shw->GetVtxPlane());
    }
    MSG("RecoBase",Msg::kDebug)
      << "    at plane " << shw->GetVtxPlane() << " dt2  " << du*du+dv*dv << "/" << tolTPos2 << " dt " << dt*1.e9 << "/" << tolTime*1e9 <<  "\n";
    
    if(dt<tolTime && 
       du*du+dv*dv<tolTPos2) return true;   // check whether vertex within minimum distance from trk 
  }

  // finally, check distance between vertices 
  Double_t trkZ=this->GetVtxZ();
  Double_t trkU=this->GetVtxU();
  Double_t trkV=this->GetVtxV();
  Double_t trkT=this->GetVtxT();
  Int_t trkPlane=this->GetVtxPlane();
  Double_t trkDirU=this->GetVtxDirCosU();
  Double_t trkDirV=this->GetVtxDirCosV();
  Double_t trkDirZ=this->GetVtxDirCosZ();

  if(fabs(shw->GetVtxZ()-this->GetEndZ())<fabs(shw->GetVtxZ()-this->GetVtxZ())){
    trkZ=this->GetEndZ();
    trkU=this->GetEndU();
    trkV=this->GetEndV();
    trkT=this->GetEndT();
    trkDirU=this->GetEndDirCosU();
    trkDirV=this->GetEndDirCosV();
    trkDirZ=this->GetEndDirCosZ();
    trkPlane=this->GetEndPlane();
  }
  
  Double_t dz=trkZ-shw->GetVtxZ();
  dt= trkT-shw->GetVtxT(); 
  //compensate dz for SM gap if necessary
  if (vldcptr->GetDetector()==Detector::kFar &&
      trkPlane<=pSMPlaneLast &&
      shw->GetVtxPlane()>=pSMPlaneFirst) { 
    dz+=zGapSM;
  } else if (vldcptr->GetDetector()==Detector::kFar &&
           trkPlane>=pSMPlaneLast &&
	     shw->GetVtxPlane()<=pSMPlaneFirst) { 
    dz-=zGapSM;
  }
  
  Double_t du = (trkU - dz*trkDirU/trkDirZ) - shw->GetVtxU();
  Double_t dv = (trkV - dz*trkDirV/trkDirZ) - shw->GetVtxV();
  MSG("RecoBase",Msg::kDebug)
    << "    dvertex shower/track " << du
    << " " << dv << " " << dz <<" "  << dt*1.e9 << "\n";
  if (du*du+dv*dv<tolTPos2 &&
      fabs(dz)<tolZPos && 
      fabs(dt)<tolTime) return true;

  // now try shower End
  dz=trkZ-shw->GetEndZ();
  dt= trkT-shw->GetVtxT(); 
  //compensate dz for SM gap if necessary
  if (vldcptr->GetDetector()==Detector::kFar &&
      trkPlane<=pSMPlaneLast &&
      shw->GetVtxPlane()>=pSMPlaneFirst) { 
    dz+=zGapSM;
  } else if (vldcptr->GetDetector()==Detector::kFar &&
           trkPlane>=pSMPlaneLast &&
	     shw->GetEndPlane()<=pSMPlaneFirst) {
    dz-=zGapSM;
  }
  du = (trkU - dz*trkDirU/trkDirZ) - shw->GetEndU();
  dv = (trkV - dz*trkDirV/trkDirZ) - shw->GetEndV();
  
  MSG("RecoBase",Msg::kDebug)
    << "    dvertex shower/track " << du
    << " " << dv << " " << dz <<" "  << dt*1.e9 << "\n";
  if (du*du+dv*dv<tolTPos2 &&
      fabs(dz)<tolZPos && 
      fabs(dt)<tolTime) return true;
  
  return false;
}

//_____________________________________________________________________

Bool_t InoTrackCand::IsUnphysical(Float_t trkFrac,Float_t asymCut,
                                     Float_t xtalkFrac,Float_t xtalkCut) {
  // loop through track planes and calculate #gaps/#planes 
  // also consider balance of U/V hit planes
  // (for ND consider only area covered by partial plane)
  
  Float_t totTrkPlanes = 0;   // total number of hit planes in track
  Float_t nTrkPlanesU  = 0;   // total number of valid planes 
  Float_t nTrkPlanesV  = 0;   // between beg/end in U and V
  Float_t nHitPlanesU  = 0;   // total number of valid hit planes 
  Float_t nHitPlanesV  = 0;   // between beg/end in U and V
  Float_t nXTalkPlanes = 0;   // total number of xtalk-like track planes
  
  UgliGeomHandle ugh(*this->GetVldContext());
  for(int ipln=this->GetBegPlane();ipln<=this->GetEndPlane();ipln++){
    MSG("RecoBase",Msg::kDebug) << " plane = " << ipln << endl;
    
    //add up total number of track planes containing hits:
    if(this->IsTPosValid(ipln)) totTrkPlanes += 1;
    
    PlexPlaneId scintid(this->GetVldContext()->GetDetector(),ipln,kFALSE);
    UgliScintPlnHandle scintpln = ugh.GetScintPlnHandle(scintid);
    if(scintpln.IsValid()){
      MSG("RecoBase",Msg::kDebug) << "Valid plane (has scintillator)" << endl;
      //add up total number of valid track planes:
      if(this->GetVldContext()->GetDetector()==Detector::kFar) {        
        if(scintid.GetPlaneView()==PlaneView::kU) { 
	  nTrkPlanesU += 1;
        } else if(scintid.GetPlaneView()==PlaneView::kV) {
	  nTrkPlanesV += 1;
	}
      }  else if(this->GetVldContext()->GetDetector()==Detector::kNear){  //for neardet, only consider the partial plane region:
        Float_t z = scintpln.GetZ0();
        Float_t u = this->GetVtxU();
        Float_t v = this->GetVtxV();
        if(this->IsTPosValid(ipln)) {
          u = this->GetU(ipln); 
          v = this->GetV(ipln);
        } else {
          Int_t count = 1;
          while(!this->IsTPosValid(ipln+count) && 
                ipln+count<this->GetEndPlane()) count +=1;
          PlexPlaneId scintid_upp(this->GetVldContext()->GetDetector(),
                                  ipln+count,kFALSE);
          UgliScintPlnHandle scintpln_upp = ugh.GetScintPlnHandle(scintid_upp);
          Float_t upp_u = this->GetU(ipln+count);
          Float_t upp_v = this->GetV(ipln+count);
          Float_t upp_z = scintpln_upp.GetZ0();
          count = 1;
          while(!this->IsTPosValid(ipln-count) && 
                ipln-count>this->GetBegPlane()) count +=1;
          PlexPlaneId scintid_low(this->GetVldContext()->GetDetector(),
                                  ipln-count,kFALSE);
          UgliScintPlnHandle scintpln_low = ugh.GetScintPlnHandle(scintid_low);
          Float_t low_u = this->GetU(ipln-count);
          Float_t low_v = this->GetV(ipln-count);
          Float_t low_z = scintpln_low.GetZ0();
          u = (z - low_z)*(upp_u - low_u)/(upp_z - low_z);
          v = (z - low_z)*(upp_v - low_v)/(upp_z - low_z);
        }
        
        //define a region for the ND to do this test in:
        const Float_t part_u_min = -0.2 * Munits::m;
        const Float_t part_v_min = -2.4 * Munits::m;
        const Float_t part_u_max = 2.4 * Munits::m;
        const Float_t part_v_max = 0.2 * Munits::m;
        if( u<part_u_min || u>part_u_max || //partial planes
            v<part_v_min || v>part_v_max ) continue;
        
        Float_t r2 = u*u+v*v;
        const Float_t coil_r     =  0.5 * Munits::m;
        const Float_t detect_r   =  2.8 * Munits::m;
        const Float_t coil_r2    = coil_r * coil_r;
        const Float_t detect_r2  = detect_r * detect_r; 
        if( r2<=coil_r2 || r2>=detect_r2 ) continue;  //coil hole + detector edge   
	
        if(scintid.GetPlaneView()==PlaneView::kU) { 
	  nTrkPlanesU += 1;
	} else if(scintid.GetPlaneView()==PlaneView::kV) { 
	  nTrkPlanesV += 1;
        }
      }
      if(this->IsTPosValid(ipln)) {
        //add upp number of valid planes that have hits:
        if(scintid.GetPlaneView()==PlaneView::kU) { 
	  nHitPlanesU += 1;
        } else if(scintid.GetPlaneView()==PlaneView::kV) {
	  nHitPlanesV += 1;
	}
        //also count number of planes that have only xtalk-like hits
        if(this->GetPlaneCharge(ipln,CalStripType::kPE)<xtalkCut) {
          nXTalkPlanes+=1;
	}
      }
    }
  }
  
  MSG("RecoBase",Msg::kDebug) 
    << endl << "Total Number of track planes:" 
    << ( this->GetEndPlane() - this->GetBegPlane() + 1) << endl
    << "Total Number of hit track planes:" << totTrkPlanes << endl
    << "Number of (valid) track planes U: " << nTrkPlanesU 
    << " V: " << nTrkPlanesV << endl
    << "Number of (valid) hit track planes U: " << nHitPlanesU 
    << " V: " << nHitPlanesV << endl
    << "Number of xtalk-like (valid) hit track planes:  " 
    << nXTalkPlanes << endl;
  
  //if track hits are mainly in spectrometer in ND do not perform tests
  if(this->GetVldContext()->GetDetector()==Detector::kNear){
    if(nHitPlanesU+nHitPlanesV<totTrkPlanes-nHitPlanesU-nHitPlanesV) 
      return false;
  }

  //check for nonsensical values:
  if(nTrkPlanesU+nTrkPlanesV<=0) return true; 
  if(nHitPlanesU+nHitPlanesV<=0) return true;

  //now check for gap planes, xtalk and asymmetry:
  if((nHitPlanesU+nHitPlanesV)/(nTrkPlanesU+nTrkPlanesV)<trkFrac) {
    MSG("RecoBase",Msg::kDebug) 
      << "IsUnphysical because trkFrac = "
      << (nHitPlanesU+nHitPlanesV)/(nTrkPlanesU+nTrkPlanesV)
      << " and cut-off is = " << trkFrac << endl;
    return true;
  }
  
  if(nXTalkPlanes/(nHitPlanesU+nHitPlanesV)>xtalkFrac) {
    MSG("RecoBase",Msg::kDebug) 
      << "IsUnphysical because xtalkFrac = "
      << nXTalkPlanes/(nHitPlanesU+nHitPlanesV)
      << " and cut-off is = " << xtalkFrac << endl;
    return true;
  }
  
  if(TMath::Abs(nHitPlanesU-nHitPlanesV)/
     TMath::Max(nHitPlanesU,nHitPlanesV)>asymCut) {
    MSG("RecoBase",Msg::kDebug) 
      << "IsUnphysical because asym = "
      << ( TMath::Abs(nHitPlanesU-nHitPlanesV) / 
           TMath::Max(nHitPlanesU,nHitPlanesV) )
      << " and cut-off is = " << asymCut << endl;
    return true;
  }

  return false;
}
*/
void InoTrackCand::SetNTrackStrip(Int_t n) {
  fNTrackStrip = n;
}

void InoTrackCand::SetNTrackDigit(Int_t n) {
  fNTrackDigit = n;
}

void InoTrackCand::SetNTimeFitDigit(Int_t n) {
  fNTimeFitDigit = n;
}

void InoTrackCand::SetTimeFitChi2(Double_t x) {
  fTimeFitChi2 = x;
}

void InoTrackCand::SetTimeForwardFitRMS(Double_t x) {
  fTimeForwardFitRMS = x;
}
void InoTrackCand::SetTimeForwardFitNDOF(Int_t x) {
  fTimeForwardFitNDOF = x;
}

void InoTrackCand::SetTimeBackwardFitRMS(Double_t x) {
  fTimeBackwardFitRMS = x;
}

void InoTrackCand::SetTimeBackwardFitNDOF(Int_t x) {
  fTimeBackwardFitNDOF = x;
}

Int_t InoTrackCand::GetNTrackStrip() const {
  return fNTrackStrip;
}

Int_t InoTrackCand::GetNTrackDigit() const {
  return fNTrackDigit;
}

Int_t InoTrackCand::GetNTimeFitDigit() const {
  return fNTimeFitDigit;
}

Double_t InoTrackCand::GetTimeFitChi2() const {
  return fTimeFitChi2;
}

Double_t InoTrackCand::GetTimeForwardFitRMS() const {
  return fTimeForwardFitRMS;
}

Int_t InoTrackCand::GetTimeForwardFitNDOF() const {
  return fTimeForwardFitNDOF;
}

Double_t InoTrackCand::GetTimeBackwardFitRMS() const {
  return fTimeBackwardFitRMS;
}

Int_t InoTrackCand::GetTimeBackwardFitNDOF() const {
  return fTimeBackwardFitNDOF;
}

bool InoTrackCand::IsContained() {
  
  bool contained = (GetVtxTrace()>0.1 && GetEndTrace()>0.1 &&
                    GetVtxnActiveUpstream()>2 && GetEndnActiveDownstream()>2 &&
                    GetVtxDistToEdge()>0.1 && GetEndDistToEdge()>0.1 );
  return contained;
}

Int_t InoTrackCand::GetNDaughters() const {
  return ClustsInTrack.size();
}

//______________________________________________________________________
Int_t InoTrackCand::GetNStrip(Int_t iuv) const {
  Int_t n=0;
  for (unsigned i=0; i<ClustsInTrack.size(); i++) {
    InoCluster* inoclust = ClustsInTrack[i];
    for (unsigned j=0; j<inoclust->HitsInCluster.size(); j++) {
      InoHit* inohit = inoclust->GetHit(j);
      if (inohit->GetView()==2 || inohit->GetView()==iuv) {n++;}
    }
    
  }
  return n;
  
  /*
  if (iuv==2) { 
    return fTrack->HitsInTrack.size(); 
  } else {
    Int_t n=0;
    for (unsigned i=0; i<fTrack->HitsInTrack.size(); i++) {
      InoHit* inohit = fTrack->HitsInTrack[i];
      if (inohit->GetView()==2 || inohit->GetView()==iuv) {n++;}
    }
    return n;
  }
  return -999;
  */
}

Double_t InoTrackCand::GetPlanePulse(Int_t plane) const {
  for (unsigned i=0; i<ClustsInTrack.size(); i++) {
    InoCluster* inoclust = ClustsInTrack[i];
    if (inoclust->GetZPlane()==plane) return inoclust->GetPulse(); 
  }
  return -999;
  
  /* 08/02/09
  for (unsigned i=0; i<fTrack->HitsInTrack.size(); i++) {
    InoHit* inohit = fTrack->HitsInTrack[i];
    if (inohit->GetZPlane()==plane) return inohit->GetPulse(); 
  }
  return -999;
  */
}

Double_t InoTrackCand::GetPulse() const {
  Double_t pulse = 0;
  
  for (unsigned i=0; i<ClustsInTrack.size(); i++) {
    InoCluster* inoclust = ClustsInTrack[i];
    for (unsigned j=0; j<inoclust->HitsInCluster.size(); j++) {
      InoHit* inohit = inoclust->GetHit(j);
      pulse +=inohit->GetPulse();
    }
  }
  
  /* 08/02/09
  for (unsigned i=0; i<fTrack->HitsInTrack.size(); i++) {
    InoHit* inohit = fTrack->HitsInTrack[i];
    pulse +=inohit->GetPulse();
  }
  */

  return pulse;
}


Int_t InoTrackCand::GetNDigit(Int_t) const {
  //GMA put it properly

  return ClustsInTrack.size();
  // 08/02/09  return fTrack->HitsInTrack.size();
  /*
  Int_t n=0;
  TIter stripItr(GetDaughterIterator());
  CandStripHandle *strip;
  while ((strip = dynamic_cast<CandStripHandle*>(stripItr()))) {
    n += strip->GetNDigit(stripend_t);
  }
  return n;
  */
}

Int_t InoTrackCand::GetBegPlane(Int_t iuv) const  {
  if (iuv==2) { 
    return ClustsInTrack[0]->GetZPlane();
  } else {
    for (unsigned i=0; i<ClustsInTrack.size(); i++) {
      if (ClustsInTrack[i]->GetView()==iuv) {
	return ClustsInTrack[i]->GetZPlane();
      }
    }
  }
  return 5000;
}
						    
Int_t InoTrackCand::GetEndPlane(Int_t iuv) const {
  int nsize = ClustsInTrack.size();
  if (iuv==2) { 
    return ClustsInTrack[nsize-1]->GetZPlane();
  } else {
    for (int i=nsize-1; i>=0; i--) {
      if (ClustsInTrack[i]->GetView()==iuv) {
	return ClustsInTrack[i]->GetZPlane();
      }
    }
  }
  return -20;
}

//______________________________________________________________________
Int_t InoTrackCand::GetNPlane(Int_t iuv) const {
  //  if (planeview_t==PlaneView::kUnknown) {return GetNDaughters();  }
  //GMA : Modify this for more than one hits in a layer
  //  cout <<" iuv "<< iuv<<endl;
  if (iuv==2) { 
    return ClustsInTrack.size(); 
  } else {
    Int_t n=0;
    //    cout <<"iuv "<< iuv<<" "<<ClustsInTrack.size()<<endl;
    for (unsigned i=0; i<ClustsInTrack.size(); i++) {
      InoCluster* inoclust = ClustsInTrack[i];
      //      cout <<"inoclust->GetView() "<< inoclust->GetView()<<endl;
      if (inoclust->GetView()==2 || inoclust->GetView()==iuv) {n++;}
    }
    return n;
  }
}


//NavKey InoTrackCand::KeyFromSlice(const InoTrackCand *reco) {
//  if (reco->GetCandSlice()) {return static_cast<Int_t>(reco->GetCandSlice()->GetUidInt()); }
//  return 0;
//}

//______________________________________________________________________
void InoTrackCand::SetVtxU(Double_t dvar) {
  fVertex->SetU(dvar);
}

Double_t InoTrackCand::GetVtxU() const {
  return fVertex->GetU();
}

//______________________________________________________________________

void InoTrackCand::SetVtxV(Double_t dvar) {
  fVertex->SetV(dvar);
}

Double_t InoTrackCand::GetVtxV() const {
  return fVertex->GetV();
}

//______________________________________________________________________

void InoTrackCand::SetVtxZ(Double_t dvar) {
  fVertex->SetZ(dvar);
}

Double_t InoTrackCand::GetVtxZ() const {
  return fVertex->GetZ();
}

//______________________________________________________________________

void InoTrackCand::SetVtxT(Double_t dvar) {
  fVertex->SetT(dvar);
}

Double_t InoTrackCand::GetVtxT() const {
  return fVertex->GetT();
}

//______________________________________________________________________

void InoTrackCand::SetVtxPlane(Int_t ivar) {
  fVertex->SetPlane(ivar);
}

Int_t InoTrackCand::GetVtxPlane() const {
  return fVertex->GetPlane();
}

void InoTrackCand::SetVtxRPCmod(Int_t ivar) {
  for(unsigned int ij=0; ij<ClustsInTrack.size(); ij++) {
    if(ClustsInTrack[ij]->GetZPlane() == ivar) {
      fVertex->SetRPCmod(ClustsInTrack[ij]->GetRPCmod());
      break;
    }
  }
}

Int_t InoTrackCand::GetVtxRPCmod() {
  return fVertex->GetRPCmod();
}

//______________________________________________________________________

void InoTrackCand::SetEndU(Double_t dvar) {
  fTerm->SetU(dvar);
}

Double_t InoTrackCand::GetEndU() const {
  return fTerm->GetU();
}

//______________________________________________________________________
void InoTrackCand::SetEndV(Double_t dvar) {
  fTerm->SetV(dvar);
}


Double_t InoTrackCand::GetEndV() const {
  return fTerm->GetV();
}

//______________________________________________________________________
void InoTrackCand::SetEndZ(Double_t dvar) {
  fTerm->SetZ(dvar);
}

Double_t InoTrackCand::GetEndZ() const {
  return fTerm->GetZ();
}

//______________________________________________________________________

void InoTrackCand::SetEndT(Double_t dvar) {
  fTerm->SetT(dvar);
}

Double_t InoTrackCand::GetEndT() const {
  return fTerm->GetT();
}

//______________________________________________________________________

void InoTrackCand::SetEndPlane(Int_t ivar) {
  fTerm->SetPlane(ivar);
}

Int_t InoTrackCand::GetEndPlane() const {
  return fTerm->GetPlane();
}

void InoTrackCand::SetEndRPCmod(Int_t ivar) {
  for(unsigned int ij=0; ij<ClustsInTrack.size(); ij++) {
    if(ClustsInTrack[ij]->GetZPlane() == ivar) {
      fTerm->SetRPCmod(ClustsInTrack[ij]->GetRPCmod());
      break;
    }
  }
}

Int_t InoTrackCand::GetEndRPCmod() {
  return fTerm->GetRPCmod();
}
//______________________________________________________________________
void InoTrackCand::SetDirCosU(Double_t dvar) {
  fVertex->GetDirCosine().SetX(dvar);
}

Double_t InoTrackCand::GetDirCosU() const {
  return fVertex->GetDirCosine().X();
}

//______________________________________________________________________
void InoTrackCand::SetDirCosV(Double_t dvar) {
 fVertex->GetDirCosine().SetY(dvar);
}

Double_t InoTrackCand::GetDirCosV() const {
  return fVertex->GetDirCosine().Y();
}

//______________________________________________________________________
void InoTrackCand::SetDirCosZ(Double_t dvar) {
  fVertex->GetDirCosine().SetZ(dvar);
}

Double_t InoTrackCand::GetDirCosZ() const {
  return fVertex->GetDirCosine().Z();
}

//______________________________________________________________________
void InoTrackCand::SetVtxDirCosU(Double_t dvar) {
  fVertex->GetDirCosine().SetX(dvar);
}

Double_t InoTrackCand::GetVtxDirCosU() const {
  return fVertex->GetDirCosine().X();
}

//______________________________________________________________________
void InoTrackCand::SetVtxDirCosV(Double_t dvar) {
  fVertex->GetDirCosine().SetY(dvar);
}

Double_t InoTrackCand::GetVtxDirCosV() const {
  return fVertex->GetDirCosine().Y();
}

//______________________________________________________________________
void InoTrackCand::SetVtxDirCosZ(Double_t dvar) {
  fVertex->GetDirCosine().SetZ(dvar);
}

Double_t InoTrackCand::GetVtxDirCosZ() const {
  return fVertex->GetDirCosine().Z();
}

//______________________________________________________________________
void InoTrackCand::SetEndDirCosU(Double_t dvar) {
  fTerm->GetDirCosine().SetX(dvar);
}

Double_t InoTrackCand::GetEndDirCosU() const {
  return  fTerm->GetDirCosine().X();
}

//______________________________________________________________________
void InoTrackCand::SetEndDirCosV(Double_t dvar) {
  fTerm->GetDirCosine().SetX(dvar);
}

Double_t InoTrackCand::GetEndDirCosV() const {
 return  fTerm->GetDirCosine().Y();
}

//______________________________________________________________________
void InoTrackCand::SetEndDirCosZ(Double_t dvar) {
  fTerm->GetDirCosine().SetZ(dvar);
}

Double_t InoTrackCand::GetEndDirCosZ() const {
  return  fTerm->GetDirCosine().Z();
}

//______________________________________________________________________
void InoTrackCand::SetTimeSlope(Double_t dvar) {
  fTimeSlope = dvar;
}

Double_t InoTrackCand::GetTimeSlope() const {
  return fTimeSlope;
}

//______________________________________________________________________
void InoTrackCand::SetTimeOffset(Double_t dvar) {
  fTimeOffset = dvar;
}

Double_t InoTrackCand::GetTimeOffset() const {
  return fTimeOffset;
}

//InoFittedTrackHandle

double InoTrackCand::GetMomentumRange() const {
  return mMomentumRange;
}

void InoTrackCand::SetMomentumRange(double momentum) {
  mMomentumRange = momentum;
}

double InoTrackCand::GetMomentumdS() const {
  return mMomentumdS;
}

void InoTrackCand::SetMomentumdS(double momentum) {
  mMomentumdS = momentum;
}

double InoTrackCand::GetFinderMomentum() const {
  return mFinderMomentum;
}

void InoTrackCand::SetFinderMomentum(double momentum) {
  mFinderMomentum = momentum;
}

double InoTrackCand::GetMomentumCurve() const {
  return mMomentumCurve;
}

void InoTrackCand::SetMomentumCurve(double momentum) {
  mMomentumCurve = momentum;
}

double InoTrackCand::GetEndMomentumCurve() const {
  return mEndMomentumCurve;
}

void InoTrackCand::SetEndMomentumCurve(double momentum) {
  mEndMomentumCurve = momentum;
}

double InoTrackCand::GetEMCharge() const {
  return mEMCharge;
}

void InoTrackCand::SetEMCharge(double emcharge) {
  mEMCharge = emcharge;
}

//----------------------------------------------------------------------
void InoTrackCand::SetVtxQPError(double error) {
  mVtxQPError = error;
}

//----------------------------------------------------------------------
double InoTrackCand::GetVtxQPError() const {
  return mVtxQPError;
}

//----------------------------------------------------------------------
void InoTrackCand::SetVtxUError(double error) {
  mVtxUError = error;
}

//----------------------------------------------------------------------
double InoTrackCand::GetVtxUError() const {
  return mVtxUError;
}

//----------------------------------------------------------------------
void InoTrackCand::SetVtxVError(double error) {
  mVtxVError = error;
}

//----------------------------------------------------------------------
double InoTrackCand::GetVtxVError() const {
  return mVtxVError;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void InoTrackCand::SetVtxXX(double XX) {
	mVtxXX = XX;
}
//----------------------------------------------------------------------
double InoTrackCand::GetVtxXX() const {
	return mVtxXX;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void InoTrackCand::SetVtxYY(double YY) {
	mVtxYY = YY;
}
//----------------------------------------------------------------------
double InoTrackCand::GetVtxYY() const {
	return mVtxYY;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void InoTrackCand::SetVtxTX(double TX) {
	mVtxTX = TX;
}
//----------------------------------------------------------------------
double InoTrackCand::GetVtxTX() const {
	return mVtxTX;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void InoTrackCand::SetVtxTY(double TY) {
	mVtxTY = TY;
}
//----------------------------------------------------------------------
double InoTrackCand::GetVtxTY() const {
	return mVtxTY;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void InoTrackCand::SetVtxdU(double tx) {
  mVtxdU = tx;
}

//----------------------------------------------------------------------
double InoTrackCand::GetVtxdU() const {
  return mVtxdU;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void InoTrackCand::SetVtxdV(double ty) {
  mVtxdV = ty;
}

//----------------------------------------------------------------------
double InoTrackCand::GetVtxdV() const {
  return mVtxdV;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//----------------------------------------------------------------------
void InoTrackCand::SetVtxdUError(double error) {
  mVtxdUError = error;
}

//----------------------------------------------------------------------
double InoTrackCand::GetVtxdUError() const {
  return mVtxdUError;
}

//----------------------------------------------------------------------
void InoTrackCand::SetVtxdVError(double error) {
  mVtxdVError = error;
}

//----------------------------------------------------------------------
double InoTrackCand::GetVtxdVError() const {
  return mVtxdVError;
}

//----------------------------------------------------------------------
Double_t InoTrackCand::GetEndQP() const {
  return mEndQP;
}

//----------------------------------------------------------------------
void InoTrackCand::SetEndQP(Double_t qp0) {
  mEndQP = qp0;
}

//----------------------------------------------------------------------
Double_t InoTrackCand::GetBave() const {
  return mBave;
}

//----------------------------------------------------------------------
void InoTrackCand::SetBave(Double_t bave) {
  mBave = bave;
}

//----------------------------------------------------------------------
Float_t InoTrackCand::GetPlaneChi2(Int_t iplane) const
{
  /*
  map<Int_t,Float_t>::iterator iter;
  iter = (dynamic_cast<const InoFittedTrack *>(GetCandBase())->fPlaneChi2).find(iplane);
  if (iter!=(dynamic_cast<const InoFittedTrack *>(GetCandBase())->fPlaneChi2).end()) {
    return dynamic_cast<const InoFittedTrack *>(GetCandBase())->fPlaneChi2[iplane];
  }
  return -1.;
  */
  return mPlaneChi2[iplane];
}

//----------------------------------------------------------------------
Float_t InoTrackCand::GetPlaneQP(Int_t iplane) const
{
  /*
  map<Int_t,Float_t>::iterator iter;
  iter = (dynamic_cast<const InoFittedTrack *>(GetCandBase())->fPlaneQP).find(iplane);
  if (iter!=(dynamic_cast<const InoFittedTrack *>(GetCandBase())->fPlaneQP).end()) {
    return dynamic_cast<const InoFittedTrack *>(GetCandBase())->fPlaneQP[iplane];
  }
  return -1.;
  */
  return mPlaneQP[iplane];
}

//----------------------------------------------------------------------
void InoTrackCand::SetPlaneQP(Int_t iplane, Double_t qp) {
  mPlaneQP[iplane] = (Float_t)qp;
}

//----------------------------------------------------------------------
void InoTrackCand::SetPlaneChi2(Int_t iplane, Double_t chi2) {
  mPlaneChi2[iplane] = (Float_t)chi2;
}

//----------------------------------------------------------------------
void InoTrackCand::SetEndUError(Double_t error) {
  mEndUError = error;
}

//----------------------------------------------------------------------
void InoTrackCand::SetEndVError(Double_t error) {
  mEndVError = error;
}

//----------------------------------------------------------------------
void InoTrackCand::SetEnddUError(Double_t error) {
  mEnddUError = error;
}

//----------------------------------------------------------------------
void InoTrackCand::SetEnddVError(Double_t error) {
  mEnddVError = error;
}

//----------------------------------------------------------------------
void InoTrackCand::SetEndQPError(Double_t error) {
  mEndQPError = error;
}

//----------------------------------------------------------------------
Double_t InoTrackCand::GetEndUError() const {
  return mEndUError;
}

//----------------------------------------------------------------------
Double_t InoTrackCand::GetEndVError() const {
  return mEndVError;
}

//----------------------------------------------------------------------
Double_t InoTrackCand::GetEnddUError() const {
  return mEnddUError;
}

//----------------------------------------------------------------------
Double_t InoTrackCand::GetEnddVError() const {
  return mEnddVError;
}

//----------------------------------------------------------------------
Double_t InoTrackCand::GetEndQPError() const {
  return mEndQPError;
}

//----------------------------------------------------------------------
void InoTrackCand::SetNSwimFail(Int_t nswimfail) {
  mNSwimFail = nswimfail;
}

//----------------------------------------------------------------------
Int_t InoTrackCand::GetNSwimFail() const {
  return mNSwimFail;
}

Double_t InoTrackCand::GetChi2() const {
  return mChi2;
}

void InoTrackCand::SetChi2(Double_t chi2) {
  mChi2 = chi2;
}

//----------------------------------------------------------------------
void InoTrackCand::SetNDOF(Int_t ndof) {
  mNDOF = ndof;
}

//----------------------------------------------------------------------
Int_t InoTrackCand::GetNDOF() const {
  return mNDOF;
}

//----------------------------------------------------------------------
Double_t InoTrackCand::GetRangeBiasedQP() const {
  return fQP_rangebiased;
}

//----------------------------------------------------------------------
void InoTrackCand::SetRangeBiasedQP(Double_t qp0) {
  fQP_rangebiased = qp0;
}

//----------------------------------------------------------------------
Int_t InoTrackCand::GetNIterate() const {
  return fNIterate;
}

//----------------------------------------------------------------------
void InoTrackCand::SetNIterate(Int_t nit) {
  fNIterate = nit;
}

/*
CalTimeType::CalTimeType_t InoTrackCand::GetCalTimeType() const {
  TIter stripItr(GetDaughterIterator());
  CandStripHandle *strip = dynamic_cast<CandStripHandle*>(stripItr());
  if (strip) {
    return strip->GetCalTimeType();
  }
  else {
    return CalTimeType::kNone;
  }
}

//______________________________________________________________________
void InoTrackCand::CalibrateSigMapped(UInt_t encoded, Float_t ph) {
  CandReco *candreco = dynamic_cast<CandReco *>(GetOwnedCandBase());
  candreco->fSigMapped[encoded] = ph;
}

//______________________________________________________________________
void InoTrackCand::CalibrateMIP(UInt_t encoded, Float_t ph) {
  CandReco *candreco = dynamic_cast<CandReco *>(GetOwnedCandBase());
  candreco->fMIP[encoded] = ph;
}
*/

