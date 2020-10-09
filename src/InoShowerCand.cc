#include "TMath.h"

#include <cassert>
#include <cmath>
#include <iostream>
//#include <map>

#include "InoCluster.h"
#include "InoHit.h"
#include "InoTrackCand.h"
//#include "InoFitterTrackHandle.h"
#include "InoShowerCand.h"

/*
//______________________________________________________________________
InoShowerCand::InoShowerCand()
{
}


//______________________________________________________________________
InoShowerCand::InoShowerCand(const InoShowerCand &cdh) :
  CandRecoHandle(cdh)
{
}

//______________________________________________________________________
InoShowerCand::InoShowerCand(CandShower *cd) :
  CandRecoHandle(cd)
{
}

//______________________________________________________________________
InoShowerCand::~InoShowerCand()
{
}

//______________________________________________________________________
InoShowerCand *InoShowerCand::DupHandle() const
{
   return (new InoShowerCand(*this));
}

//______________________________________________________________________
void InoShowerCand::Trace(const char *c) const
{
  MSG("Cand", Msg::kDebug)
    << "**********Begin InoShowerCand::Trace(\"" << c << "\")" << endl
           << "Information from InoShowerCand's CandHandle: " << endl;
  CandHandle::Trace(c);
  MSG("Cand", Msg::kDebug)
     << "**********End InoShowerCand::Trace(\"" << c << "\")" << endl;
}






//______________________________________________________________________
void InoShowerCand::SetU(Int_t plane, Float_t tpos)  {
  CandShower *shower = dynamic_cast<CandShower*>(GetOwnedCandBase());
  shower->fUPos[plane] = tpos;
}

Float_t InoShowerCand::GetU(Int_t plane) const {
  const CandShower *shower = dynamic_cast<const CandShower*>(GetCandBase());
  if ((shower->fUPos).count(plane)) return shower->fUPos[plane];
  return -99999.;
}

//______________________________________________________________________

void InoShowerCand::SetV(Int_t plane, Float_t tpos)  {
  CandShower *shower = dynamic_cast<CandShower*>(GetOwnedCandBase());
  shower->fVPos[plane] = tpos;
}

Float_t InoShowerCand::GetV(Int_t plane) const {
  const CandShower *shower = dynamic_cast<const CandShower*>(GetCandBase());
  if ((shower->fVPos).count(plane)) return shower->fVPos[plane];
  return -99999.;
}

//______________________________________________________________________

Float_t InoShowerCand::GetZ(Int_t plane) const {
  TIter stripItr(GetDaughterIterator());
  while (CandStripHandle *strip =
	 dynamic_cast<CandStripHandle*>(stripItr())) {
    if (strip->GetPlane()==plane) return strip->GetZPos();
  }
  return -1.;
}

//______________________________________________________________________

void InoShowerCand::SetT(Int_t plane, StripEnd::StripEnd_t stripend_t, Double_t time) {
  CandShower *shower = dynamic_cast<CandShower*>(GetOwnedCandBase());
  switch (stripend_t) {
  case StripEnd::kNegative:
    shower->fTime[0][plane] = time;
    break;
  case StripEnd::kPositive:
    shower->fTime[1][plane] = time;
    break;
  default:
    shower->fTime[0][plane] = time;
    break;
  }
}

//______________________________________________________________________

void InoShowerCand::ClearUVT() {
  CandShower *shower = dynamic_cast<CandShower*>(GetOwnedCandBase());
  shower->fUPos.clear();
  shower->fVPos.clear();
  shower->fTime[0].clear();
  shower->fTime[1].clear();
}

//_____________________________________________________________________

Float_t InoShowerCand::GetMinU(Int_t iplane, Double_t minPE) const{
  CandStripHandleItr sItr(GetDaughterIterator());
  CandStripHandle * begstrip = sItr();
  if(!begstrip)return -999;
  const VldContext * vld = begstrip->GetVldContext();
  if(!vld)return -999;
  PlexPlaneId plnid(vld->GetDetector(),iplane,false);
  if(plnid.GetPlaneView()!=PlaneView::kU)return -999.;
  CandStripHandleItr stripItr(GetDaughterIterator());
  CandStripHandleKeyFunc *stripKf = stripItr.CreateKeyFunc();
  stripKf->SetFun(CandStripHandle::KeyFromPlane);
  stripItr.GetSet()->AdoptSortKeyFunc(stripKf);
  stripKf = 0;
  stripItr.GetSet()->Slice(iplane);

  Float_t minU=999.;
  while (const CandStripHandle *strip =
         dynamic_cast<const CandStripHandle*>(stripItr())) {
    
    if(iplane<GetBegPlane() || iplane > GetEndPlane()) return 0.;
    
    if(strip->GetTPos() < minU && strip->GetCharge(CalDigitType::kPE) > minPE )
       minU=strip->GetTPos();
  }
  if(minU==999)minU=0;
  return minU;
}
//_____________________________________________________________________

Float_t InoShowerCand::GetMinV(Int_t iplane, Double_t minPE) const{
  CandStripHandleItr sItr(GetDaughterIterator());
  CandStripHandle * begstrip = sItr();
  if(!begstrip)return -999;
  const VldContext * vld = begstrip->GetVldContext();
  if(!vld)return -999;
  PlexPlaneId plnid(vld->GetDetector(),iplane,false);
  if(plnid.GetPlaneView()!=PlaneView::kV)return -999.;
  if(iplane<GetBegPlane() || iplane > GetEndPlane()) return 0.;
  CandStripHandleItr stripItr(GetDaughterIterator());
  CandStripHandleKeyFunc *stripKf = stripItr.CreateKeyFunc();
  stripKf->SetFun(CandStripHandle::KeyFromPlane);
  stripItr.GetSet()->AdoptSortKeyFunc(stripKf);
  stripKf = 0;
  stripItr.GetSet()->Slice(iplane);
 
  Float_t minV=999.;
  while (const CandStripHandle *strip =
         dynamic_cast<const CandStripHandle*>(stripItr())) {
    if(strip->GetTPos()<minV && strip->GetCharge(CalDigitType::kPE) > minPE )
       minV=strip->GetTPos();
  }
  if(minV==999.)minV=0;
  return minV;
}

//_____________________________________________________________________

Float_t InoShowerCand::GetMaxU(Int_t iplane, Double_t minPE) const{
  CandStripHandleItr sItr(GetDaughterIterator());
  CandStripHandle * begstrip = sItr();
  if(!begstrip)return 999;
  const VldContext * vld = begstrip->GetVldContext();
  if(!vld)return 999;
  PlexPlaneId plnid(vld->GetDetector(),iplane,false);
  if(plnid.GetPlaneView()!=PlaneView::kU)return 999.;
  if(iplane<GetBegPlane() || iplane > GetEndPlane()) return 0.;

  CandStripHandleItr stripItr(GetDaughterIterator());
  CandStripHandleKeyFunc *stripKf = stripItr.CreateKeyFunc();
  stripKf->SetFun(CandStripHandle::KeyFromPlane);
  stripItr.GetSet()->AdoptSortKeyFunc(stripKf);
  stripKf = 0;
  stripItr.GetSet()->Slice(iplane);
 
  Float_t maxU=-999.;
  while (const CandStripHandle *strip =
         dynamic_cast<const CandStripHandle*>(stripItr())) {
    if(strip->GetTPos() > maxU&& strip->GetCharge(CalDigitType::kPE) > minPE ) maxU=strip->GetTPos();
  }
  if(maxU==-999.)maxU=0;
  return maxU;
}

//_____________________________________________________________________

Float_t InoShowerCand::GetMaxV(Int_t iplane, Double_t minPE) const {
  CandStripHandleItr sItr(GetDaughterIterator());
  CandStripHandle * begstrip = sItr();
  if(!begstrip)return 999;
  const VldContext * vld = begstrip->GetVldContext();
  if(!vld)return 999;
  PlexPlaneId plnid(vld->GetDetector(),iplane,false);
  if(plnid.GetPlaneView()!=PlaneView::kV)return 999.;
  if(iplane<GetBegPlane() || iplane > GetEndPlane()) return 0.;

  CandStripHandleItr stripItr(GetDaughterIterator());
  CandStripHandleKeyFunc *stripKf = stripItr.CreateKeyFunc();
  stripKf->SetFun(CandStripHandle::KeyFromPlane);
  stripItr.GetSet()->AdoptSortKeyFunc(stripKf);
  stripKf = 0;
  stripItr.GetSet()->Slice(iplane);
  Float_t maxV=-999.;
  while (const CandStripHandle *strip =
         dynamic_cast<const CandStripHandle*>(stripItr())) {
    if(strip->GetTPos()>maxV && strip->GetCharge(CalDigitType::kPE)  > minPE)
       maxV=strip->GetTPos();
  }
  if(maxV==-999)maxV=0;
  return maxV;
}
//_____________________________________________________________________

Int_t InoShowerCand::GetNStrips(Int_t iplane){
  Int_t nstrips=0;
  CandStripHandleItr stripItr(GetDaughterIterator());
  while (const CandStripHandle *strip =
         dynamic_cast<const CandStripHandle*>(stripItr())) {
    if(strip->GetPlane()==iplane)nstrips++;
  }
  return nstrips;
}

//_____________________________________________________________________

Double_t InoShowerCand::GetT(Int_t plane,StripEnd::StripEnd_t stripend_t) const {
  const CandShower *shower = dynamic_cast<const CandShower*>(GetCandBase());
  if (stripend_t==StripEnd::kNegative) {
    if ((shower->fTime[0]).count(plane)) {
      return shower->fTime[0][plane];
    }
    return -99999.;
  }
  if (stripend_t==StripEnd::kPositive) {
    if ((shower->fTime[1]).count(plane)) {
      return shower->fTime[1][plane];
    }
    return -99999.;
  }
  if ((shower->fTime[0]).count(plane) && (shower->fTime[1]).count(plane)) {
    return min(shower->fTime[0][plane],shower->fTime[1][plane]);
  }
  else if ((shower->fTime[0]).count(plane)) {
    return shower->fTime[0][plane];
  }
  else if ((shower->fTime[1]).count(plane)) {
    return shower->fTime[1][plane];
  }
  return -99999.;
}

Double_t InoShowerCand::GetT(StripEnd::StripEnd_t stripend_t,Int_t plane) const {
  return this->GetT(plane,stripend_t);
}

Double_t InoShowerCand::GetT(Int_t plane) const {
  return this->GetT(plane,StripEnd::kWhole);
}

//__________________________________________________________________
Bool_t InoShowerCand::IsTPosValid(Int_t plane) const {
  const CandShower *shower = dynamic_cast<const CandShower*>(GetCandBase());
  if ((shower->fUPos).count(plane) && (shower->fVPos).count(plane)) {
    return kTRUE;
  }
  return kFALSE;
}
//_____________________________________________________________________
Bool_t InoShowerCand::BelongsWithTrack(CandTrackHandle * trk, 
				       AlgConfig & ac, 
				       const VldContext * vldcptr, 
				       Double_t tolTPos2, Double_t tolZPos, Double_t tolTime){
  
  if(!trk)return false;
  
  Int_t pSMPlaneLast = ac.GetInt("SMPlaneLast");
  Int_t pSMPlaneFirst = ac.GetInt("SMPlaneFirst");
  VldContext vldc = *vldcptr;
  UgliGeomHandle ugh(vldc); 
  Double_t zGapSM=0;
  // "Fixed": MAK, Feb 8, 2005
  // GetScintPlnHandle was being called in the ND for planes which didn't
  // have any scint. This fix only stops the block below from being called
  // in the near detector.  It doesn't require that 
  // pSMPlaneLast, pSMPlaneFirst actually refer to planes with scintillator!!
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
  
  MSG("RecoBase",Msg::kDebug) 
    << " comparing shower at " 
    << this->GetVtxPlane() 
    << " " << this->GetVtxU() 
    << " " << this->GetVtxV() 
    << " with track at " << trk->GetVtxPlane() 
    << " " << trk->GetVtxU() 
    << " " << trk->GetVtxV() << endl;
  Float_t tolTPos = TMath::Sqrt(tolTPos2); 
  Int_t fwdTrkPlane=min(trk->GetVtxPlane(),trk->GetEndPlane());
  Int_t bckTrkPlane=max(trk->GetVtxPlane(),trk->GetEndPlane());
  Int_t fwdShwPlane=min(this->GetVtxPlane(),this->GetEndPlane());
  Int_t bckShwPlane=max(this->GetVtxPlane(),this->GetEndPlane());
  
  Double_t dt = trk->GetVtxT()-this->GetVtxT();
  //  if(trk->IsTPosValid(this->GetVtxPlane())){
  //   dt = this->GetVtxT()-trk->GetT(this->GetVtxPlane());
  //   }
  
  if( bckShwPlane>=fwdTrkPlane && fwdShwPlane<=bckTrkPlane && fabs(dt)<tolTime){   // trk overlaps in Z and Time
    MSG("RecoBase",Msg::kDebug) << " trk and shower overlap in Z and time " << endl;
    Bool_t matchu=false;
    Bool_t matchv=false;
    for(Int_t iloop=0;iloop<2;iloop++){      
      if(iloop >0 && (matchu || matchv))tolTPos=tolTPos*2;
      // if on second loop and no matches so far, give up
      if(iloop >0 && !matchu && !matchv) break;
      for (Int_t iplane=fwdShwPlane; iplane<=bckShwPlane;iplane++){  // look for trk passing through shower
	
        Double_t trkU(0.),trkV(0.);
        PlexPlaneId plnid(GetVldContext()->GetDetector(),iplane,false);
        if(plnid.GetPlaneView()==PlaneView::kU){
          trkU=trk->GetU(iplane);
          if(!trk->IsTPosValid(iplane)){
            if(abs(iplane-trk->GetVtxPlane())<abs(iplane-trk->GetEndPlane())){
              UgliScintPlnHandle scintpln = ugh.GetScintPlnHandle(plnid);
              Double_t dz=scintpln.GetZ0()-trk->GetVtxZ();
              trkU = trk->GetVtxU() + dz*trk->GetVtxDirCosU()/trk->GetVtxDirCosZ();
            } else {
              UgliScintPlnHandle scintpln = ugh.GetScintPlnHandle(plnid);
              Double_t dz=scintpln.GetZ0()-trk->GetEndZ();
              trkU = (trk->GetEndU() + dz*trk->GetEndDirCosU()/trk->GetEndDirCosZ());
            }
          }
          if(trkU>=this->GetMinU(iplane)-tolTPos && trkU<=this->GetMaxU(iplane)+tolTPos  && this->GetNStrips(iplane)>0){
            MSG("RecoBase", Msg::kDebug)<< " u match! (details follow)" << endl;
            matchu=true;
          }
        } else if(plnid.GetPlaneView()==PlaneView::kV) {
          trkV=trk->GetV(iplane);
          if(!trk->IsTPosValid(iplane)){
            if(abs(iplane-trk->GetVtxPlane())<abs(iplane-trk->GetEndPlane())){
              UgliScintPlnHandle scintpln = ugh.GetScintPlnHandle(plnid);
              Double_t dz=scintpln.GetZ0()-trk->GetVtxZ();
              trkV = (trk->GetVtxV() + dz*trk->GetDirCosV()/trk->GetDirCosZ());
            } else {
              UgliScintPlnHandle scintpln = ugh.GetScintPlnHandle(plnid);
              Double_t dz=scintpln.GetZ0()-trk->GetEndZ();
              trkV= (trk->GetEndV() + dz*trk->GetEndDirCosV()/trk->GetEndDirCosZ());
            }
          }
          if(trkV>=this->GetMinV(iplane)-tolTPos && trkV<=this->GetMaxV(iplane)+tolTPos && this->GetNStrips(iplane)>0){
            MSG("RecoBase", Msg::kDebug)<< " v match! (details follow)" << endl;
            matchv=true;
          }
        }
        MSG("RecoBase",Msg::kDebug) << " plane  " << iplane 
                                    << " trk u " << trkU 
                                    << " min/max U " << this->GetMinU(iplane) 
                                    << "/" << this->GetMaxU(iplane) 
                                    << " trk v " << trkV 
                                    << " min/max V " << this->GetMinV(iplane) 
                                    << "/" << this->GetMaxV(iplane) << endl; 
        if(matchu && matchv) {
          MSG("RecoBase", Msg::kDebug)<< " 3D  match! " << endl;
          return true;
        }
      }
    }
  }
  
  if( this->GetVtxZ()>min(trk->GetVtxZ(),trk->GetEndZ())-tolZPos  &&
      this->GetVtxZ()<max(trk->GetVtxZ(),trk->GetEndZ())+tolZPos) {
    Double_t du= trk->GetVtxU()-this->GetVtxU();
    Double_t dv= trk->GetVtxV()-this->GetVtxV();        
    if(trk->IsTPosValid(this->GetVtxPlane())){
      MSG("RecoBase",Msg::kDebug) << " TPos is valid @ shower vtx " 
                                  << " track u/v " << trk->GetU(this->GetVtxPlane()) 
                                  << " " << trk->GetV(this->GetVtxPlane()) << endl;
      
      du = this->GetVtxU()-trk->GetU(this->GetVtxPlane());
      dv = this->GetVtxV()-trk->GetV(this->GetVtxPlane());
    }
    MSG("RecoBase",Msg::kDebug)
      << "    at plane " << this->GetVtxPlane() << " dt2  " 
      << du*du+dv*dv << "/" << tolTPos2
      << " dt " << dt*1.e9 
      << "/" << tolTime*1e9 <<  "\n";
    if(dt<tolTime && 
       du*du+dv*dv<tolTPos2) return true;   // check whether vertex within minimum distance from trk 
  }
  // finally, check distance between vertices 
  Double_t trkZ=trk->GetVtxZ();
  Double_t trkU=trk->GetVtxU();
  Double_t trkV=trk->GetVtxV();
  Double_t trkT=trk->GetVtxT();
  Double_t trkDirU=trk->GetVtxDirCosU();
  Double_t trkDirV=trk->GetVtxDirCosV();
  Double_t trkDirZ=trk->GetVtxDirCosZ();
  Int_t trkPlane=trk->GetVtxPlane();
  if(fabs(this->GetVtxZ()-trk->GetEndZ())<fabs(this->GetVtxZ()-trk->GetVtxZ())){
    trkZ=trk->GetEndZ();
    trkU=trk->GetEndU();
    trkV=trk->GetEndV();
    trkT=trk->GetEndT();
    trkPlane=trk->GetEndPlane();
    trkDirU=trk->GetEndDirCosU();
    trkDirV=trk->GetEndDirCosV();
    trkDirZ=trk->GetEndDirCosZ();
  }
  
  Double_t dz = this->GetVtxZ()-trkZ;
  dt = this->GetVtxT()-trkT;
  //compensate dz for SM gap if necessary
  if (vldcptr->GetDetector()==Detector::kFar &&
      this->GetVtxPlane()<=pSMPlaneLast &&
      trkPlane>=pSMPlaneFirst) {
    dz+=zGapSM;
  } else if (vldcptr->GetDetector()==Detector::kFar &&
           this->GetVtxPlane()>=pSMPlaneLast &&
	     trkPlane<=pSMPlaneFirst) {
    dz-=zGapSM;
  }

  Double_t du = (trkU + dz*trkDirU/trkDirZ) - this->GetVtxU();
  Double_t dv = (trkV + dz*trkDirV/trkDirZ) - this->GetVtxV();
  
  MSG("RecoBase",Msg::kDebug)
    << "    dvertex shower/track " << du
    << " " << dv << " " << dz <<" "  << dt*1.e9 << "\n";
  if (du*du+dv*dv<tolTPos2 &&
      fabs(dz)<tolZPos && 
      fabs(dt)<tolTime) return true;
  
  dz = this->GetEndZ()-trkZ;
  dt = this->GetVtxT()-trkT;
  //compensate dz for SM gap if necessary
  if (vldcptr->GetDetector()==Detector::kFar &&
      this->GetVtxPlane()<=pSMPlaneLast &&
      trkPlane>=pSMPlaneFirst) { 
    dz+=zGapSM;
  } else if (vldcptr->GetDetector()==Detector::kFar &&
           this->GetVtxPlane()>=pSMPlaneLast &&
	     trkPlane<=pSMPlaneFirst) { 
    dz-=zGapSM;
  }
  du = (trkU + dz*trkDirU/trkDirZ) - this->GetEndU();
  dv = (trkV + dz*trkDirV/trkDirZ) - this->GetEndV();
  
  MSG("RecoBase",Msg::kDebug)
    << "    dvertex shower/track " << du
    << " " << dv << " " << dz <<" "  << dt*1.e9 << "\n";
  if (du*du+dv*dv<tolTPos2 &&
      fabs(dz)<tolZPos && 
      fabs(dt)<tolTime) return true;

  return false;
}
//_____________________________________________________________________
Bool_t InoShowerCand::BelongsWithShower(InoShowerCand * shw, 
                                          AlgConfig & ac, 
                                          const VldContext * vldcptr, 
                                          Double_t tolTPos2, Double_t tolZPos, Double_t tolTime){

  if(!shw)return false;

  Int_t pSMPlaneLast = ac.GetInt("SMPlaneLast");
  Int_t pSMPlaneFirst = ac.GetInt("SMPlaneFirst");
  VldContext vldc = *vldcptr;
  UgliGeomHandle ugh(vldc); 
  Double_t zGapSM=0;

  // "Fixed": MAK, Feb 8, 2005
  // GetScintPlnHandle was being called in the ND for planes which didn't
  // have any scint. This fix only stops the block below from being called
  // in the near detector.  It doesn't require that 
  // pSMPlaneLast, pSMPlaneFirst actually refer to planes with scintillator!!
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
  
  MSG("RecoBase",Msg::kDebug) << " shower extents " << this->GetVtxZ() << " " << this->GetEndZ() << " " << shw->GetVtxZ() << " " << shw->GetEndZ() << endl;
  Double_t dz = this->GetVtxZ()-shw->GetVtxZ();
  Double_t dzend= this->GetEndZ()-shw->GetVtxZ();
  Double_t dzend2= this->GetVtxZ()-shw->GetEndZ();
  Double_t dzend3= this->GetEndZ()-shw->GetEndZ();
  Double_t du= this->GetVtxU()-shw->GetVtxU();
  Double_t dv= this->GetVtxV()-shw->GetVtxV();  
  Double_t dt= this->GetVtxT()-shw->GetVtxT();
 
  //compensate dz for SM gap if necessary
  if (vldcptr->GetDetector()==Detector::kFar &&
      this->GetVtxPlane()<=pSMPlaneLast &&
      shw->GetVtxPlane()>=pSMPlaneFirst) { 
    dz+=zGapSM;
  } else if (vldcptr->GetDetector()==Detector::kFar &&
           this->GetVtxPlane()>=pSMPlaneLast &&
	     shw->GetVtxPlane()<=pSMPlaneFirst) { 
    dz-=zGapSM;
  }
  
  if (vldcptr->GetDetector()==Detector::kFar &&
      this->GetEndPlane()<=pSMPlaneLast &&
      shw->GetVtxPlane()>=pSMPlaneFirst) { 
    dzend+=zGapSM;
  } else if (vldcptr->GetDetector()==Detector::kFar &&
           this->GetEndPlane()>=pSMPlaneLast &&
	     shw->GetVtxPlane()<=pSMPlaneFirst) { 
    dzend-=zGapSM;
  }
  if(fabs(dz)>fabs(dzend)){
    du= this->GetEndU()-shw->GetVtxU();
    dv= this->GetEndV()-shw->GetVtxV(); 
    dz=dzend;
  }
  if (vldcptr->GetDetector()==Detector::kFar &&
      this->GetVtxPlane()<=pSMPlaneLast &&
      shw->GetEndPlane()>=pSMPlaneFirst) { 
    dzend2+=zGapSM;
  } else if (vldcptr->GetDetector()==Detector::kFar &&
           this->GetVtxPlane()>=pSMPlaneLast &&
	     shw->GetEndPlane()<=pSMPlaneFirst) { 
    dzend2-=zGapSM;
  }
  if(fabs(dz)>fabs(dzend2)){
    du= this->GetVtxU()-shw->GetEndU();
    dv= this->GetVtxV()-shw->GetEndV(); 
    dz=dzend2;
  }
  if (vldcptr->GetDetector()==Detector::kFar &&
      this->GetEndPlane()<=pSMPlaneLast &&
      shw->GetEndPlane()>=pSMPlaneFirst) { 
    dzend3+=zGapSM;
  } else if (vldcptr->GetDetector()==Detector::kFar &&
           this->GetEndPlane()>=pSMPlaneLast &&
	     shw->GetEndPlane()<=pSMPlaneFirst) { 
    dzend3-=zGapSM;
  }
  if(fabs(dz)>fabs(dzend3)){
    du= this->GetEndU()-shw->GetEndU();
    dv= this->GetEndV()-shw->GetEndV(); 
    dz=dzend3;
  }

 Double_t tolTp=tolTPos2;
  Double_t tolZ=tolZPos;
  if(this->GetEnergy()>15.0 || shw->GetEnergy()>15.0){
    tolTp=tolTPos2*4;
  }
  MSG("RecoBase",Msg::kDebug)
    << "    dvertex shower/shower " << du
    << " " << dv << " " << dz <<" "  << dt*1.e9 << "\n";
  if(( (du*du+dv*dv<tolTp &&fabs(dz)<tolZ) ||  
       (du*du<tolTp/4. && dv*dv<tolTp*4. && fabs(dz)<tolZ/2.) ||
       (du*du<tolTp*4 && dv*dv<tolTp/4. && fabs(dz)<tolZ/2.))
     &&  fabs(dt)<tolTime) return true;
  
  return false;
}

//______________________________________________________________________
void InoShowerCand::AddCluster(CandClusterHandle *cluster) {
  const TObjArray &clusterlist =
    (dynamic_cast<CandShower*>(GetOwnedCandBase()))->fClusterList;
  Bool_t found(0);
  for (Int_t i=0; i<=clusterlist.GetLast() && !found; i++) {
    CandClusterHandle *target =
      dynamic_cast<CandClusterHandle*>(clusterlist.At(i));
    if (*cluster == *target) found = 1;
  }
  if (!found) {                 // don't want to duplicate object in list
    CandClusterHandle * cl=cluster->DupHandle();
    (dynamic_cast<CandShower*>(GetOwnedCandBase()))->fClusterList.
      Add(cl);
  }
  return;
}
//______________________________________________________________________
void InoShowerCand::RemoveCluster(CandClusterHandle *cluster) {
  const TObjArray &clusterlist =
    (dynamic_cast<CandShower*>(GetOwnedCandBase()))->fClusterList;
  Bool_t found(0);
  for (Int_t i=0; i<=clusterlist.GetLast() && !found; i++) {
    CandClusterHandle *target =
      dynamic_cast<CandClusterHandle*>(clusterlist.At(i));
    if (*cluster == *target){
      (dynamic_cast<CandShower*>(GetOwnedCandBase()))->fClusterList.
        RemoveAt(i);
      delete target;
      (dynamic_cast<CandShower*>(GetOwnedCandBase()))->fClusterList.
	Compress();
      return;
    }
  }
  return;
}

//_____________________________________________________________________
bool InoShowerCand::IsContained(){
  // determine whether shower is contained.

  CandStripHandleItr sItr(GetDaughterIterator());
  CandStripHandle * begstrip = sItr();
  if(!begstrip)return true;
  const VldContext * vld = begstrip->GetVldContext();
  if(!vld)return true;

  Bool_t contained =true;
  double totcharge = 0;
  double containedcharge = 0;
  UgliGeomHandle ugh(*vld); 
  PlaneOutline pl;
  float detzmin=0;
  float detzmax=999;
  ugh.GetZExtent(detzmin,detzmax,-1);
  float u=0,v=0;
  float xedge,yedge,dist;

  while(CandStripHandle * strip = sItr()){
    float z = strip->GetZPos();
    if(z<detzmin+0.1*Munits::m ||
       z>detzmax-0.1*Munits::m) contained=false;
    PlexPlaneId plnid(vld->GetDetector(),strip->GetPlane(),false);
    if(plnid.GetPlaneView()==PlaneView::kV){
      v=strip->GetTPos();
      u=GetU(strip->GetPlane());
    }
    else if(plnid.GetPlaneView()==PlaneView::kU){
      u=strip->GetTPos();
      v=GetV(strip->GetPlane());
    }
    float x = 0.707*(u-v);
    float y = 0.707*(u+v);
    pl.DistanceToEdge(x, y,
                      plnid.GetPlaneView(),
                      plnid.GetPlaneCoverage(),
                      dist, xedge, yedge);
    bool isInside = pl.IsInside( x, y,
                                 plnid.GetPlaneView(),
                                 plnid.GetPlaneCoverage());
    
    totcharge += strip->GetCharge();
    isInside &= dist>0.1*Munits::m;
    if(isInside) {
      containedcharge +=strip->GetCharge();
    }
  }
  if(totcharge>0) contained = containedcharge/totcharge>0.9;
  return contained;
}

//______________________________________________________________________

const CandClusterHandle *InoShowerCand::GetCluster(Int_t i) const {
  const TObjArray *fClusterList = &(dynamic_cast<const CandShower*>
				    (GetCandBase())->fClusterList);
  if (i>fClusterList->GetLast()) {
    return 0;
  }
  return dynamic_cast<const CandClusterHandle*>(fClusterList->At(i));
}

//______________________________________________________________________

Int_t InoShowerCand::GetLastCluster() const {
  return dynamic_cast<const CandShower*>(GetCandBase())->fClusterList.GetLast();
}

//______________________________________________________________________

Bool_t InoShowerCand::IsUnphysical(Float_t xtalkFrac,Float_t xtalkCut) {
  if(this->GetNStrip()<=0) return true;
  Float_t nxtalk = 0;
  CandStripHandleItr stripItr(GetDaughterIterator());
  while (const CandStripHandle *strip =
         dynamic_cast<const CandStripHandle*>(stripItr())) {
    if(strip->GetCharge(CalDigitType::kPE) < xtalkCut) nxtalk+=1;      
  }
  if(nxtalk/Float_t(this->GetNStrip())>xtalkFrac) return true;
  return false;
}

//______________________________________________________________________

void InoShowerCand::SetEnergy(Double_t energy, 
                                 InoShowerCand::ShowerType_t showertype) {
  
  switch (showertype){
  case  kWtCC:
    dynamic_cast<CandShower*>(GetOwnedCandBase())->fEnergy = energy;
    dynamic_cast<CandShower*>(GetOwnedCandBase())->fEnergy_wtCC = energy;
    break;
  case kCC:
    dynamic_cast<CandShower*>(GetOwnedCandBase())->fEnergy_CC = energy;
    break;
  case kWtNC:
    dynamic_cast<CandShower*>(GetOwnedCandBase())->fEnergy_wtNC = energy;
    break;
  case kNC:
    dynamic_cast<CandShower*>(GetOwnedCandBase())->fEnergy_NC = energy;
    break;
  case kEM:
    dynamic_cast<CandShower*>(GetOwnedCandBase())->fEnergy_EM = energy;
    break;
  }
}

Double_t InoShowerCand::GetEnergy(InoShowerCand::ShowerType_t showertype) const {
  switch (showertype){
  case kWtCC:
    return dynamic_cast<const CandShower*>(GetCandBase())->fEnergy_wtCC;
  case kCC:
    return dynamic_cast<const CandShower*>(GetCandBase())->fEnergy_CC;
  case kWtNC:
    return dynamic_cast<const CandShower*>(GetCandBase())->fEnergy_wtNC;
  case kNC:
    return dynamic_cast<const CandShower*>(GetCandBase())->fEnergy_NC;
  case kEM:
    return dynamic_cast<const CandShower*>(GetCandBase())->fEnergy_EM;
  default:  
    MSG("RecoBase",Msg::kWarning) << " requested shower energy of unknown type " << endl;
    return 0;
  }
}

void InoShowerCand::CalibrateEnergy(CandTrackHandle * trk,  AlgConfig & ac){
  
  CandStripHandleItr sItr(GetDaughterIterator());
  CandStripHandle * begstrip = sItr();
  const VldContext * vld = begstrip->GetVldContext();

  Double_t fCCWtLowScale = ac.GetDouble("CCWtLowScale");
  Double_t fCCWtLowC1 = ac.GetDouble("CCWtLowC1");
  Double_t fCCWtLowC2 = ac.GetDouble("CCWtLowC2");
  Double_t fCCWtLowC3 = ac.GetDouble("CCWtLowC3");
  Double_t fCCWtLowC4 = ac.GetDouble("CCWtLowC4");
  Double_t fCCWtHiC0 = ac.GetDouble("CCWtHiC0");
  Double_t fCCWtHiC1 = ac.GetDouble("CCWtHiC1");

  Double_t fCCLinLowScale = ac.GetDouble("CCLinLowScale");
  Double_t fCCLinLowC1 = ac.GetDouble("CCLinLowC1");
  Double_t fCCLinLowC2 = ac.GetDouble("CCLinLowC2");
  Double_t fCCLinHiC0 = ac.GetDouble("CCLinHiC0");
  Double_t fCCLinHiC1 = ac.GetDouble("CCLinHiC1");

  Double_t fNCWtLowScale = ac.GetDouble("NCWtLowScale");
  Double_t fNCWtMidScale = ac.GetDouble("NCWtMidScale");
  Double_t fNCWtLowC1 = ac.GetDouble("NCWtLowC1");
  Double_t fNCWtLowC2 = ac.GetDouble("NCWtLowC2");
  Double_t fNCWtLowC3 = ac.GetDouble("NCWtLowC3");
  Double_t fNCWtMidC0 = ac.GetDouble("NCWtMidC0");
  Double_t fNCWtMidC1 = ac.GetDouble("NCWtMidC1");
  Double_t fNCWtMidC2 = ac.GetDouble("NCWtMidC2");
  Double_t fNCWtMidC3 = ac.GetDouble("NCWtMidC3");
  Double_t fNCWtHiC0 = ac.GetDouble("NCWtHiC0");
  Double_t fNCWtHiC1 = ac.GetDouble("NCWtHiC1");

  Double_t fNCLinLowScale = ac.GetDouble("NCLinLowScale");
  Double_t fNCLinLowC1 = ac.GetDouble("NCLinLowC1");
  Double_t fNCLinLowC2 = ac.GetDouble("NCLinLowC2");
  Double_t fNCLinHiC0 = ac.GetDouble("NCLinHiC0");
  Double_t fNCLinHiC1 = ac.GetDouble("NCLinHiC1");

  Double_t fEMLowScale = ac.GetDouble("EMLowScale");
  Double_t fEMLowC1 = ac.GetDouble("EMLowC1");
  Double_t fEMLowC2 = ac.GetDouble("EMLowC2");
  Double_t fEMHiC0 = ac.GetDouble("EMHiC0");
  Double_t fEMHiC1 = ac.GetDouble("EMHiC1");

  Double_t shw_linmipsum=GetCharge(CalStripType::kMIP);
  Double_t shw_wtmipsum=0;
  Double_t totshw_wtmipsum=0;
  Double_t totshw_linmipsum = GetCharge(CalStripType::kMIP);

  // loop through strips once, getting linear version of shower sum, used
  // to obtain deweighting factors as well

  Int_t maxpln = 0;
  Int_t minpln = 999;
  Int_t shared_planes = 0;
  Int_t shared_strips = 0;
  Int_t nshwstp = 0;
  Double_t dedx_1 = 0.;
  Double_t dedx_2 = 0.;
  Double_t avg_dedx = 0.;
  Double_t reco_emu = 0.;
  Double_t reco_dircosz = 0.;
  Double_t reco_emu_endshw=0.;

  if (trk) {
    reco_dircosz = fabs(trk->GetDirCosZ());
    CandFitTrackHandle* fittrk = dynamic_cast<CandFitTrackHandle*>(trk); // if there is a fitted track then use it for improved tracking + direction in the shower
    reco_emu = sqrt(trk->GetMomentum()*trk->GetMomentum()+ 0.10566*0.10566);
    if(fittrk){
      reco_dircosz = fabs(fittrk->GetDirCosZ());
      if(fittrk->GetMomentum()>0.) {reco_emu = sqrt(fittrk->GetMomentum()*fittrk->GetMomentum()+ 0.10566*0.10566);}
    }
    
    CandStripHandleItr stripItr(GetDaughterIterator());
    // check all strips in shower,and see determine length of track through shower
    while (CandStripHandle *strip = dynamic_cast<CandStripHandle*>(stripItr())) {
      if (strip) {
        nshwstp++;
        if((!fittrk && trk->FindDaughter(strip)) || (fittrk && fittrk->FindDaughter(strip))){     
          shared_strips++;
          if(strip->GetPlane()>maxpln) maxpln = strip->GetPlane();
          if(strip->GetPlane()<minpln) minpln = strip->GetPlane(); 
        }
      }
    }
    shared_planes = maxpln - minpln + 1;
    
    reco_emu_endshw = reco_emu-shared_planes/(30.*reco_dircosz);
    
    if(fittrk && fittrk->GetPlaneQP(maxpln)!=0.){
      Double_t pln_qp = fittrk->GetPlaneQP(maxpln)*(1.013*fabs(fittrk->GetPlaneQP(maxpln)));
      reco_emu_endshw = sqrt(1./pln_qp*1./pln_qp + 0.10566*0.10566);
    }

    if(reco_emu_endshw>reco_emu){reco_emu_endshw = reco_emu-shared_planes/(30.*reco_dircosz);}// get best estimate for muon energy at the end of the shw using fittrk info if present

    dedx_1 = DeDx(reco_emu*1000)/1.985; // get dedx at vtx, normalised to MIPs
    if(reco_emu_endshw>0.106){ // stop non-physical situations where more energy is lost than initially present
      dedx_2 = DeDx(reco_emu_endshw*1000)/1.985; // get approx dedx at end of shw
      avg_dedx = (dedx_1 > dedx_2)? 0.5*(dedx_1+dedx_2):dedx_1; // stop muons which stop in the shw from skewing the dedx factor
    }
    else{avg_dedx = dedx_1;}
    shw_linmipsum -= avg_dedx*shared_strips/reco_dircosz; // remove muon chg from shw
  }
  else{
    shw_linmipsum=totshw_linmipsum;
  }
  
  Double_t CCEnergyLin = (shw_linmipsum<fCCLinLowScale) ? shw_linmipsum*fCCLinLowC1 + shw_linmipsum*shw_linmipsum*fCCLinLowC2: fCCLinHiC0 +  shw_linmipsum*fCCLinHiC1;
  if(shw_linmipsum<=0) CCEnergyLin=0;
        
  Double_t NCEnergyLin = (totshw_linmipsum<fNCLinLowScale) ? totshw_linmipsum*fNCLinLowC1 + totshw_linmipsum*totshw_linmipsum*fNCLinLowC2: fNCLinHiC0 +  totshw_linmipsum*fNCLinHiC1;
  if(totshw_linmipsum<=0) NCEnergyLin=0;
        
  Double_t EMEnergy = (totshw_linmipsum<fEMLowScale) ? totshw_linmipsum*fEMLowC1  + totshw_linmipsum*totshw_linmipsum*fEMLowC2: fEMHiC0 +  totshw_linmipsum*fEMHiC1;
  if(totshw_linmipsum<=0) EMEnergy=0;

  // for old software, exponential deweighting function is used.  For SS, use cubic

  Double_t deweightfactorCC = 1.0;
  Double_t deweightfactorNC = 1.0;
 
  Double_t fdeweightLowScale = ac.GetDouble("deweightLowScale");
  Double_t fdeweightC0 = ac.GetDouble("deweightC0");
  Double_t fdeweightC1 = ac.GetDouble("deweightC1");
  Double_t fdeweightC2 = ac.GetDouble("deweightC2");
  Double_t fdeweightC3 = ac.GetDouble("deweightC3");

  //  if(CCEnergyLin < fdeweightLowScale){
  //    deweightfactorCC =  fdeweightC0 + 
  //      fdeweightC1 * CCEnergyLin +
  //      fdeweightC2 * CCEnergyLin * CCEnergyLin +
  //      fdeweightC3 * CCEnergyLin * CCEnergyLin * CCEnergyLin;
  //  }

  // Find final deweight power
  double SigHalfPoint = 8.0;
  double SigGrad = 0.5;      
  deweightfactorCC = 0.25 + 0.75/(1. + exp(-(CCEnergyLin-SigHalfPoint)*SigGrad));
 
  if(NCEnergyLin < fdeweightLowScale){
    deweightfactorNC =  fdeweightC0 + 
      fdeweightC1 * NCEnergyLin +
      fdeweightC2 * NCEnergyLin * NCEnergyLin +
      fdeweightC3 * NCEnergyLin * NCEnergyLin * NCEnergyLin;
  }
  
  // now that deweighting factors are known, calculated de-weighted sum 
    CandStripHandleItr stripItr(GetDaughterIterator());
  // check all strips in shower,and see whether they are shared with the track 
  while (CandStripHandle *strip = dynamic_cast<CandStripHandle*>(stripItr())) {
    if (strip) {
      totshw_wtmipsum+=pow(GetStripCharge(strip,StripEnd::kWhole),deweightfactorNC);
      Double_t trackEloss = 0;  
      // if track crosses strip, subtract chg dependent on dedx factor (calc earlier) and direction
      if(trk){
        CandFitTrackHandle* fittrk = dynamic_cast<CandFitTrackHandle*>(trk);
        if((!fittrk && trk->FindDaughter(strip)) || (fittrk && fittrk->FindDaughter(strip))){     
          trackEloss=avg_dedx/reco_dircosz;
        }
      }
      if(GetStripCharge(strip,StripEnd::kWhole) - trackEloss > 0){
        shw_wtmipsum+=pow(GetStripCharge(strip,StripEnd::kWhole)-trackEloss,deweightfactorCC);
      }
    }
  }
  
  Double_t CCEnergyWt =  fCCWtHiC0 +  shw_wtmipsum*fCCWtHiC1;
  if( shw_wtmipsum<fCCWtLowScale){
    CCEnergyWt =  shw_wtmipsum*fCCWtLowC1 + 
      shw_wtmipsum*shw_wtmipsum*fCCWtLowC2 + 
      shw_wtmipsum*shw_wtmipsum*shw_wtmipsum*fCCWtLowC3 + 
      shw_wtmipsum*shw_wtmipsum*shw_wtmipsum*shw_wtmipsum*fCCWtLowC4;
  }
  if(shw_wtmipsum<=0) CCEnergyWt=0;
  
  Double_t NCEnergyWt = fNCWtHiC0 + totshw_wtmipsum*fNCWtHiC1;
  if(totshw_wtmipsum<fNCWtLowScale){ 
    NCEnergyWt = totshw_wtmipsum*fNCWtLowC1 + 
      totshw_wtmipsum*totshw_wtmipsum*fNCWtLowC2 + 
      totshw_wtmipsum*totshw_wtmipsum*totshw_wtmipsum*fNCWtLowC3; 
  }

  if(totshw_wtmipsum>=fNCWtLowScale && totshw_wtmipsum<fNCWtMidScale){
    NCEnergyWt =  fNCWtMidC0 +  totshw_wtmipsum*fNCWtMidC1 + totshw_wtmipsum*totshw_wtmipsum*fNCWtMidC2 + totshw_wtmipsum*totshw_wtmipsum*totshw_wtmipsum*fNCWtMidC3;
  }

  if(totshw_wtmipsum<=0) NCEnergyWt=0;

  // cout << " CC deweightmips = " << shw_wtmipsum << " NC deweight_mips = " <<totshw_wtmipsum  << endl;
//  cout << " CC deweight E = " << CCEnergyWt << " NC dwewightE = " << NCEnergyWt<< endl;
  
  // for near detector, compensate for near/far shower completeness difference
  // and observed difference between near/far mip scale on muon tracks
  
  if (vld && vld->GetDetector()==Detector::kNear){
    CCEnergyWt *= (.97 + exp(-(CCEnergyWt+11.6)/9.3)); 
    CCEnergyLin *= (.99 + exp(-(CCEnergyLin+28.2)/11.8)); 
    NCEnergyWt *= (1.01 + exp(-(NCEnergyWt+4.33)/5.)); 
    NCEnergyLin *= (1. + exp(-(NCEnergyLin+10.1)/5.)); 
    EMEnergy *= 1.06; 
  }
  
  // now apply energy corrections, corresponding to EnergyCorrections::MasakiMay17thCGScaled.
  
  if(vld){
    float etemplin = CCEnergyLin;
    CCEnergyLin = EnergyCorrections::ShowerEnergyConversionDogwood(etemplin,*vld);
  }
  
  
  SetEnergy(CCEnergyWt);
  SetEnergy(CCEnergyWt,kWtCC);
  SetEnergy(CCEnergyLin,kCC);
  SetEnergy(NCEnergyWt,kWtNC);
  SetEnergy(NCEnergyLin,kNC);
  SetEnergy(EMEnergy,kEM);
}

//______________________________________________________________________

Double_t InoShowerCand::DeDx(Double_t emu){
  
  Double_t dedx = 0.; 
  Double_t mm = 105.658389;
  Double_t mm_2 = mm*mm;

  if(emu*emu-mm_2<0.){return 0.;}

  Double_t a_2     = 1./(137.036*137.036);
  Double_t pi = 3.141;
  Double_t Na = 6.023;
  Double_t lamda_2 = 3.8616*3.8616;
  Double_t Z_A = 0.5377;
  Double_t me  = 0.51099906;
  Double_t me_2 = me*me;  
  Double_t beta = sqrt(emu*emu-mm_2)/emu;
  Double_t beta_2 = beta*beta;
  Double_t gamma = emu/105.658389;
  Double_t gamma_2 = gamma*gamma;
  Double_t I = 68.7e-6;
  Double_t I_2 = I*I;
  Double_t p = sqrt(emu*emu-mm_2);
  Double_t p_2 = p*p;
  Double_t Em = 2 * me * p_2 / ( me_2 + mm_2 + 2*me*emu );
  Double_t Em_2= Em*Em;
  Double_t emu_2 = emu*emu;
  Double_t X0 = 0.165;
  Double_t X1 = 2.503;
  Double_t X = log10(beta*gamma);
  Double_t a = 0.165;
  Double_t C = -3.3;
  Double_t m = 3.222;
  Double_t d = 0;
  if(X0 < X && X < X1){d = 4.6052 * X + a * pow(X1-X,m) + C;}
  if(X > X1){d = 4.6052 * X + C;}

  dedx =  a_2 * 2*pi * Na * lamda_2 * Z_A * (me / beta_2) *( log( 2*me*beta_2*gamma_2*Em/I_2 ) - 2*beta_2 + 0.25*(Em_2/emu_2) - d );
  
  return 10*dedx;
  
}


//______________________________________________________________________
// Navigation Helpers
//______________________________________________________________________

NavKey InoShowerCand::KeyFromSlice(const InoShowerCand *reco) {
  if (reco->GetCandSlice()) {
    return static_cast<Int_t>(reco->GetCandSlice()->GetUidInt());
  }
  return 0;

}
*/
