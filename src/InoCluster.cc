#include "InoCluster.h"
#include "InoHit.h"
#include "TMath.h"
#include <iostream>
InoCluster::InoCluster(InoHit* hit) :
  fZPlane(-1),fRPCmod(-1),fBegXStrip(-1), fEndXStrip(-1),fBegYStrip(-1), fEndYStrip(-1), 
  fBegTime(0.), fEndTime(0.),
  fBegXPos(-999.), fEndXPos(999.),
  fBegYPos(-999.), fEndYPos(999.),
  fXPos(-999.), fYPos(-999.), fZPos(0.), fXPulse(0.), fYPulse(0.), 
  fTrkFlag(0), fShwFlag(0),
  fTrkPlnFlag(0), fShwPlnFlag(0),
  //  fPlaneView(-1),
  fDigits(0),fNDFlag(1),       //asmQ what is the difference between fDigit and fView
  fXPosErr(999.),
  fYPosErr(999.),
  //  fXPosErr(0.02/sqrt(12.)),
  //  fYPosErr(0.02/sqrt(12.)),
  fView (2),
  InTrack(false), isStraight(true),
  StripXWidth(0.0196),
  StripYWidth(0.0196)
{
  paradef = micalDetectorParameterDef::AnPointer; //AAR: 
  StripXWidth = paradef->GetXStrwd()/1000;
  StripYWidth = paradef->GetYStrwd()/1000;
  //LayerThickness = (1/m)*2*(paradef->GetParlay(2)+paradef->GetParirlay(2)); //(1/m)*2*paradef->GetParlay(2);

  this->AddHit(hit);
}


InoCluster::~InoCluster()
{
  HitsInCluster.clear();
}


void InoCluster::AddHit(InoHit* hit) {
  if(HitsInCluster.size()==0) {
    HitsInCluster.push_back(hit);
    fZPlane=hit->GetZPlane();
    fRPCmod = hit->GetRPCmod();
    if (hit->GetXPosErr()<100) {
      fBegXStrip=hit->GetXStrip()->GetStrip();
      fEndXStrip=hit->GetXStrip()->GetStrip();
      fBegXPos=hit->GetXPos();
      fEndXPos=hit->GetXPos();
    }
    if (hit->GetYPosErr()<100) {
      fBegYStrip=hit->GetYStrip()->GetStrip();
      fEndYStrip=hit->GetYStrip()->GetStrip();
      fBegYPos=hit->GetYPos();
      fEndYPos=hit->GetYPos();
    }

    fBegTime=hit->GetTime();
    fEndTime=hit->GetTime();
    fZPos=hit->GetZPos();
    fView = hit->GetView();
    fXPosErr = hit->GetXPosErr();
    fYPosErr = hit->GetYPosErr();
    //    fPlaneView=hit->GetPlaneView();            //asmQ what is fView for why is it set to 2
  } else {
    if(this->ContainsHit(hit)==true) {return;}
    HitsInCluster.push_back(hit);   
    
    if (hit->GetXPosErr()<100) {
      if(hit->GetXStrip()->GetStrip()<fBegXStrip) fBegXStrip=hit->GetXStrip()->GetStrip();
      if(hit->GetXStrip()->GetStrip()>fEndXStrip) fEndXStrip=hit->GetXStrip()->GetStrip();
      if(hit->GetXPos()<fBegXPos) fBegXPos=hit->GetXPos();
      if(hit->GetXPos()>fEndXPos) fEndXPos=hit->GetXPos();
      if (fView ==1) fView = 2;
    }
    
    if (hit->GetYPosErr()<100) {
      if(hit->GetYStrip()->GetStrip()<fBegYStrip) fBegYStrip=hit->GetYStrip()->GetStrip();
      if(hit->GetYStrip()->GetStrip()>fEndYStrip) fEndYStrip=hit->GetYStrip()->GetStrip();
      if(hit->GetYPos()<fBegYPos) fBegYPos=hit->GetYPos();
      if(hit->GetYPos()>fEndYPos) fEndYPos=hit->GetYPos();
      if (fView ==0) fView = 2;
    }
    if(hit->GetTime()<fBegTime) fBegTime=hit->GetTime();
    if(hit->GetTime()>fEndTime) fEndTime=hit->GetTime();
  }
  
  //fDigit gives total # of X+Y  strips  //nXstrip -> total # of x strip in this cluster
  if (hit->GetXPosErr()<100) { 
    fDigits += 1; // hit->GetCandStripHandle()->GetNDaughters();
    unsigned int nxstrip = GetXEntries();
    fXPos = (fXPos*(nxstrip-1)+hit->GetXPos())/(1.0*nxstrip);
    if (nxstrip==1) {  
      fXPosErr = hit->GetXPosErr();
    } else {
      fXPosErr = hit->GetXPosErr()/2; //GMA 09/02/09 Put these separately from data   //asm why so?
    } 
    fXPulse += hit->GetXPulse();
  }
  
  if (hit->GetYPosErr()<100) { 
    fDigits += 1;       
    unsigned int nystrip = GetYEntries();
    fYPos = (fYPos*(nystrip-1)+hit->GetYPos())/(1.0*nystrip);
    if (nystrip==1) {
      fYPosErr = hit->GetYPosErr();
    } else {
      fYPosErr = hit->GetYPosErr()/2; //GMA 09/02/09 Put these separately from data
    } 
    fYPulse += hit->GetYPulse();
  }
  return;
}

bool InoCluster::ContainsHit(InoHit* hit) {
  for(unsigned int ij=0; ij<HitsInCluster.size(); ++ij) {
    if(hit==HitsInCluster[ij]) {return true;}
  }
  return false;
}

unsigned int InoCluster::GetXEntries() {
  unsigned int nxhit = 0;
  for(unsigned int ij=0; ij<HitsInCluster.size(); ++ij) {
    if (HitsInCluster[ij]->GetXPosErr() < 100) nxhit++;
  }
  return nxhit;
}

unsigned int InoCluster::GetYEntries() {
  unsigned int nyhit = 0;
  for(unsigned int ij=0; ij<HitsInCluster.size(); ++ij) {
    if (HitsInCluster[ij]->GetYPosErr() < 100) nyhit++;
  }
  return nyhit;
}

unsigned int InoCluster::GetXProjEntries() {
  unsigned int nxhit = 0; 
   
  for(unsigned int ij=0; ij<HitsInCluster.size(); ++ij) {
    bool is_identicl=0;
    if (HitsInCluster[ij]->GetXPosErr() < 100){ 
      for ( unsigned int jk=0; jk<ij;jk++){
	if( HitsInCluster[ij]->GetXStripNum() == HitsInCluster[ij]->GetXStripNum()){ 
	  is_identicl=1; break;
	}
      }
      if(!is_identicl)nxhit++; 
    }
  }
  return nxhit;
}
unsigned int InoCluster::GetYProjEntries() {
  unsigned int nyhit = 0;
  
  for(unsigned int ij=0; ij<HitsInCluster.size(); ++ij) {
    bool is_identicl=0;
    if (HitsInCluster[ij]->GetYPosErr() < 100){
      for ( unsigned int jk=0; jk<ij;jk++){
	if( HitsInCluster[ij]->GetYStripNum() == HitsInCluster[ij]->GetYStripNum()) { 
	  is_identicl=1; break;
	}
      }
      if(!is_identicl)nyhit++;
    }
  }
  return nyhit;
}

int InoCluster::IsHitAssoc(InoHit* hit) const {
  double TimeWindow = 9999.9; //GMA 290709 Optimise this
  
  if ((hit->GetZPlane()==fZPlane) &&
      (hit->GetXPosErr()>100 || (hit->GetXPos()>(fBegXPos-2*StripXWidth) && hit->GetXPos()<(fEndXPos+2*StripXWidth))) &&
      (hit->GetYPosErr()>100 || (hit->GetYPos()>(fBegYPos-2*StripYWidth) && hit->GetYPos()<(fEndYPos+2*StripYWidth))) &&
      (hit->GetTime()>(fBegTime-TimeWindow) && hit->GetTime()<(fEndTime+TimeWindow))) {
    return 1; 
  }  else {
    return 0;
  }
}

int InoCluster::IsShwAssoc(InoCluster* clust) const {
  //Need to check for signle plane hits
  double TimeWindow = 99.9; int ShwAssocNum = 0;
  
  if( (fEndTime-clust->GetBegTime())>-TimeWindow || (clust->GetEndTime()-fBegTime)>-TimeWindow ) {
    
    if( abs(clust->GetZPlane()-fZPlane)<5 &&
	(clust->GetXPos()<-100 ||
	 ((clust->GetEndXPos()-fBegXPos)>-6*StripXWidth && (fEndXPos-clust->GetBegXPos())>-6*StripXWidth)) &&
	(clust->GetYPos()<-100 ||
	 ((clust->GetEndYPos()-fBegYPos)>-6*StripYWidth && (fEndYPos-clust->GetBegYPos())>-6*StripYWidth))) {
      
      if( ( abs(clust->GetZPlane()-fZPlane)<3  &&
	    (clust->GetXPos()<-100 || 
	     ((clust->GetEndXPos()-fBegXPos)>-StripXWidth && (fEndXPos-clust->GetBegXPos())>-StripXWidth)) &&
	    (clust->GetYPos()<-100 || 
	     ((clust->GetEndYPos()-fBegYPos)>-StripYWidth && (fEndYPos-clust->GetBegYPos())>-StripYWidth)))
	  
          || ( clust->GetZPlane()==fZPlane &&
	       (clust->GetXPos()<-100 || 
	        ((clust->GetEndXPos()-fBegXPos)>-3*StripXWidth && (fEndXPos-clust->GetBegXPos())>-3*StripXWidth)) &&
	       (clust->GetYPos()<-100 || 
               ((clust->GetEndYPos()-fBegYPos)>-3*StripYWidth && (fEndYPos-clust->GetBegYPos())>-3*StripYWidth ))) ) {
	ShwAssocNum=2;
      } else {
	ShwAssocNum=1;
      }
    }
  }
  
  if(ShwAssocNum==2 && this->GetHitEntries()<3 && clust->GetHitEntries()<2) {ShwAssocNum=1;}
  
  return ShwAssocNum;
}


int InoCluster::IsTrkAssoc(InoCluster* clustm, InoCluster* clustp) const {
  double TimeWindow = 99.9; int TrkAssocNum = 0;
  double NDScale=2; //GMAA 1;
  //  double fact; //asm
  //  double min=0.2; double max=0.8;
  
  // Configure for correct detector instrumentation
  int PlaneGap = 1;
  
  // Check timing proximity
  
  if (fXPos <-100 && clustm->GetXPos()<-100 &&  clustp->GetXPos()<-100) return TrkAssocNum;
  if (fYPos <-100 && clustm->GetYPos()<-100 &&  clustp->GetYPos()<-100) return TrkAssocNum;

  if(( (fEndTime-clustm->GetBegTime())>-TimeWindow && (clustm->GetEndTime()-fBegTime)>-TimeWindow) &&
     ( (fEndTime-clustp->GetBegTime())>-TimeWindow && (clustp->GetEndTime()-fBegTime)>-TimeWindow) ) {
    
    // If more than two planes away, scale back width of cluster
    // and then treat as if only two planes away
    
    double mXPos=clustm->GetXPos();
    double mYPos=clustm->GetYPos();
    
    if((fZPlane-clustm->GetZPlane())>=PlaneGap) {
      double mScale = double(PlaneGap)/double(fZPlane-clustm->GetZPlane());
      
      mXPos=fXPos+mScale*(mXPos-fXPos);
      mYPos=fYPos+mScale*(mYPos-fYPos);
    }
    
    double pXPos=clustp->GetXPos();
    double pYPos=clustp->GetYPos();
    
    if((clustp->GetZPlane()-fZPlane)>=PlaneGap) {
      double pScale = double(PlaneGap)/double(clustp->GetZPlane()-fZPlane);
      
      pXPos=fXPos+pScale*(pXPos-fXPos);
      pYPos=fYPos+pScale*(pYPos-fYPos);
    }

    //asm[][][][][][][][][][][][][]This function is changed to get tracks with gracing angle[][][][]
    
    if (fXPos >-100 && clustm->GetXPos()>-100 &&  clustp->GetXPos()>-100 &&
	// TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*fact*StripXWidth &&
	TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*1.1*StripXWidth &&
	TMath::Abs(mXPos-fXPos) < 30*StripXWidth &&
	TMath::Abs(pXPos-fXPos) < 30*StripXWidth) {
      TrkAssocNum = 1;
    } else if (fXPos >-100 && clustm->GetXPos()>-100) {
      if (TMath::Abs(mXPos-fXPos) < 5*StripXWidth) TrkAssocNum = 1;
    } else if (fXPos >-100 && clustp->GetXPos()>-100) {
      if (TMath::Abs(pXPos-fXPos) < 5*StripXWidth) TrkAssocNum = 1;
    } else if (clustm->GetXPos() >-100 && clustp->GetXPos()>-100) {
      if (TMath::Abs(pXPos-mXPos) < 5*StripXWidth) TrkAssocNum = 1;
    } else {
      TrkAssocNum = 1;
    }

    if (fYPos >-100 && clustm->GetYPos()>-100 &&  clustp->GetYPos()>-100 &&
	// TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*fact*StripYWidth &&
	TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*1.1*StripYWidth  &&
	TMath::Abs(mYPos-fYPos) < 30*StripYWidth &&
	TMath::Abs(pYPos-fYPos) < 30*StripYWidth)  {
      TrkAssocNum += 2;
    } else if (fYPos >-100 && clustm->GetYPos()>-100) {
      if (TMath::Abs(mYPos-fYPos) < 5*StripYWidth) TrkAssocNum += 2;
    } else if (fYPos >-100 && clustp->GetYPos()>-100) {
      if (TMath::Abs(pYPos-fYPos) < 5*StripYWidth) TrkAssocNum += 2;
    } else if (clustm->GetYPos() >-100 && clustp->GetYPos()>-100) {
      if (TMath::Abs(pYPos-mYPos) < 5*StripYWidth) TrkAssocNum += 2;
    } else {
      TrkAssocNum += 2;
    }
  
    if (TrkAssocNum!=3){
      NDScale= 3;    
      if(TrkAssocNum==1 || TrkAssocNum==0){
	if (TrkAssocNum==0) NDScale= 2;  
	if( fYPos >-100 && clustm->GetYPos()>-100 &&  clustp->GetYPos()>-100 &&
            TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*3.1*StripYWidth &&
	    ((TMath::Abs(mYPos-fYPos) >= 25*StripYWidth && TMath::Abs(mYPos-fYPos) < 50*StripYWidth  &&
	      TMath::Abs(pYPos-fYPos) >= 25*StripYWidth && TMath::Abs(pYPos-fYPos) < 50*StripYWidth)||
	     (TrkAssocNum==0&&TMath::Abs(mYPos-fYPos) < 50*StripYWidth && TMath::Abs(pYPos-fYPos)<50*StripYWidth)))
	  { 
	    TrkAssocNum += 2; 
	  } 
      }
      if (TrkAssocNum==2 ){                 
	if (fXPos >-100 && clustm->GetXPos()>-100 &&  clustp->GetXPos()>-100 &&
	    TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*3.1*StripXWidth &&
	    ((TMath::Abs(mXPos-fXPos) >=25*StripXWidth && TMath::Abs(mXPos-fXPos) < 50*StripXWidth &&
	      TMath::Abs(pXPos-fXPos) >=25*StripXWidth && TMath::Abs(pXPos-fXPos) < 50*StripYWidth )||
	     (NDScale==2&&TMath::Abs(mYPos-fYPos) < 50*StripYWidth&&TMath::Abs(pYPos-fYPos)<50*StripYWidth))) 
	  {
	    TrkAssocNum=3; 
	  } 
      }
      NDScale=2;  
    }
    
    if (TrkAssocNum<3){
      TrkAssocNum = 0;
      if((TMath::Abs(mYPos-fYPos) >=50*StripYWidth && TMath::Abs(mYPos-fYPos) < 100*StripYWidth  &&
	  TMath::Abs(pYPos-fYPos) >=50*StripYWidth && TMath::Abs(pYPos-fYPos) < 100*StripYWidth  &&
	  TMath::Abs(pXPos-fXPos) < 100*StripYWidth  ) ||
	 (TMath::Abs(mXPos-fXPos) >=50*StripXWidth && TMath::Abs(mXPos-fXPos)< 100*StripXWidth &&
	  TMath::Abs(pXPos-fXPos) >=50*StripXWidth && TMath::Abs(pXPos-fXPos)< 100*StripXWidth &&
	  TMath::Abs(pYPos-fYPos) < 100*StripYWidth )){
	
	if (fYPos >-100 && clustm->GetYPos()>-100 &&  clustp->GetYPos()>-100 &&
	    TMath::Abs(mYPos+pYPos-2*fYPos) > NDScale*3.1*StripYWidth &&
	    TMath::Abs(mYPos+pYPos-2*fYPos) < NDScale*8.1*StripYWidth) {
	  TrkAssocNum = 1; 
	}
        if (fXPos >-100 && clustm->GetXPos()>-100 &&  clustp->GetXPos()>-100 &&
	    TMath::Abs(mXPos+pXPos-2*fXPos) > NDScale*3.1*StripXWidth && 
	    TMath::Abs(mXPos+pXPos-2*fXPos) < NDScale*8.1*StripXWidth) {
	  TrkAssocNum += 2; 
	}
      }
    }
    //asm[][][][][][][][][][][][][][][][][][][][]Over[][][][][][][][][][][]
    
    /*    
    // If more than two planes away, scale back width of cluster
    // and then treat as if only two planes away
    double mBegXPos=clustm->GetBegXPos();
    double mEndXPos=clustm->GetEndXPos();
    
    double mBegYPos=clustm->GetBegYPos();
    double mEndYPos=clustm->GetEndYPos();
    
    if((fZPlane-clustm->GetZPlane())>=PlaneGap)  {
    double mScale = double(PlaneGap)/double(fZPlane-clustm->GetZPlane());
    
    mBegXPos=fBegXPos+mScale*(mBegXPos-fBegXPos);
    mEndXPos=fEndXPos+mScale*(mEndXPos-fEndXPos);
    
    mBegYPos=fBegYPos+mScale*(mBegYPos-fBegYPos);
    mEndYPos=fEndYPos+mScale*(mEndYPos-fEndYPos);
    
    //    cout <<"InoCluster::trkassom "<< mScale<<" "<<mBegXPos<<" "<<mEndXPos<<" "<<mBegYPos<<" "<<mEndYPos<<endl;
    }
    
    double pBegXPos=clustp->GetBegXPos();
    double pEndXPos=clustp->GetEndXPos();
 
    double pBegYPos=clustp->GetBegYPos();
    double pEndYPos=clustp->GetEndYPos();
   
    if((clustp->GetZPlane()-fZPlane)>=PlaneGap) {
      double pScale = double(PlaneGap)/double(clustp->GetZPlane()-fZPlane);
      
      pBegXPos=fBegXPos+pScale*(pBegXPos-fBegXPos);
      pEndXPos=fEndXPos+pScale*(pEndXPos-fEndXPos);

      pBegYPos=fBegYPos+pScale*(pBegYPos-fBegYPos);
      pEndYPos=fEndYPos+pScale*(pEndYPos-fEndYPos);
      //    cout <<"InoCluster::trkassop "<< pScale<<" "<<pBegXPos<<" "<<pEndXPos<<" "<<pBegYPos<<" "<<pEndYPos<<endl;    
    } 
    

    // Scale tolerance for matching clusters for cases where there are gaps
    double k0 = 0.5*(clustp->GetZPlane()-clustm->GetZPlane()-(2*PlaneGap));

    min = min + 0.1*k0; 
    max = max - 0.1*k0; 
        
    // Determine associations
    // For tracks with +ve dtpos/dz
    //    cout <<"InoCluster::IsTrkAssoc "<< NDScale<<" "<<StriXpWidth<<" " <<fEndXPos-clustm->GetBegXPos()<<" "<<clustp->GetEndXPos()-fBegXPos<<" " <<fEndYPos-clustm->GetBegYPos()<<" " <<clustp->GetEndYPos()-fBegYPos<<endl;


    if( ((fEndXPos <-100 || clustm->GetBegXPos()<-100 || fEndXPos-clustm->GetBegXPos())>-(NDScale*1.1*StripXWidth)) && 
	((fBegXPos <-100 || clustp->GetEndXPos()<-100 || clustp->GetEndXPos()-fBegXPos)>-(NDScale*1.1*StripXWidth)) &&
	((fEndYPos <-100 || clustm->GetBegYPos()<-100 || fEndYPos-clustm->GetBegYPos())>-(NDScale*1.1*StripYWidth)) && 
	((fBegYPos <-100 || clustp->GetEndYPos()<-100 || clustp->GetEndYPos()-fBegYPos)>-(NDScale*1.1*StripYWidth))) {

      // Clusters don't overlap
      if( (fBegXPos <-100 || clustm->GetEndXPos()<-100 || (fBegXPos-clustm->GetEndXPos())>-(NDScale*0.1*StripXWidth)) ||
	  (fEndXPos <-100 || clustp->GetBegXPos()<-100 || (clustp->GetBegXPos()-fEndXPos)>-(NDScale*0.1*StripXWidth)) ||
	  (fBegYPos <-100 || clustm->GetEndYPos()<-100 || (fBegYPos-clustm->GetEndYPos())>-(NDScale*0.1*StripYWidth)) ||
	  (fEndYPos <-100 || clustp->GetBegYPos()<-100 || (clustp->GetBegYPos()-fEndYPos)>-(NDScale*0.1*StripYWidth))) {
        if( TMath::Abs( (clustp->GetBegXPos()-fEndXPos)-(fBegXPos-clustm->GetEndXPos()) )<(NDScale*2.1*StripXWidth) ||
	    ( ((min*mEndXPos)+(max*pEndXPos))>(fBegXPos-(NDScale*0.5*StripXWidth)) && ((max*mBegXPos)+(min*pBegXPos))<(fEndXPos+(NDScale*0.5*StripXWidth)) ) ||
	    TMath::Abs( (clustp->GetBegYPos()-fEndYPos)-(fBegYPos-clustm->GetEndYPos()) )<(NDScale*2.1*StripYWidth) ||
	    ( ((min*mEndYPos)+(max*pEndYPos))>(fBegYPos-(NDScale*0.5*StripYWidth)) && ((max*mBegYPos)+(min*pBegYPos))<(fEndYPos+(NDScale*0.5*StripYWidth)) ) //GMA May need criteria for single hit
	    )
          {TrkAssocNum=2;} 
      }
      
      // Overlapping clusters
      if( (clustm->GetEndXPos()<-100 || fBegXPos<-100 || (clustm->GetEndXPos()-fBegXPos)>-(NDScale*1.1*StripXWidth)) && 
	  (fEndXPos<-100 || clustp->GetBegXPos()<-100 || (fEndXPos-clustp->GetBegXPos())>-(NDScale*1.1*StripXWidth)) &&
	  (clustm->GetEndYPos()<-100 || fBegYPos<-100 || (clustm->GetEndYPos()-fBegYPos)>-(NDScale*1.1*StripYWidth)) && 
	  (fEndYPos<-100 || clustp->GetBegYPos()<-100 || (fEndYPos-clustp->GetBegYPos())>-(NDScale*1.1*StripYWidth)) ) {
        if(TrkAssocNum<1) TrkAssocNum=1;
      }
    }
    
    
    // For tracks with -ve dtpos/dz
    if( (fBegXPos <-100 ||clustm->GetEndXPos()<-100 || (fBegXPos-clustm->GetEndXPos())<(NDScale*1.1*StripXWidth)) && 
	(clustp->GetBegXPos()<-100 || fEndXPos<-100 || (clustp->GetBegXPos()-fEndXPos)<(NDScale*1.1*StripXWidth)) &&
	(fBegYPos<-100 || clustm->GetEndYPos()<-100 || (fBegYPos-clustm->GetEndYPos())<(NDScale*1.1*StripYWidth)) && 
	(clustp->GetBegYPos()<-100 || fEndYPos<-100 || (clustp->GetBegYPos()-fEndYPos)<(NDScale*1.1*StripYWidth)) ) {
      
      // Clusters don't overlap
      if( (fEndXPos<-100 || clustm->GetBegXPos()<-100 || (fEndXPos-clustm->GetBegXPos())<(NDScale*0.1*StripXWidth)) || 
	  (clustp->GetEndXPos()<-100 || fBegXPos<-100 || (clustp->GetEndXPos()-fBegXPos)<(NDScale*0.1*StripXWidth)) || 
	  (fEndYPos<-100 || clustm->GetBegYPos()<-100 || (fEndYPos-clustm->GetBegYPos())<(NDScale*0.1*StripYWidth)) || 
	  (clustp->GetEndYPos()<-100 || fBegYPos<-100 || (clustp->GetEndYPos()-fBegYPos)<(NDScale*0.1*StripYWidth))) {
        if( TMath::Abs( (clustp->GetEndXPos()-fBegXPos)-(fEndXPos-clustm->GetBegXPos()) )<(NDScale*2.1*StripXWidth) ||
	    ( ((min*pEndXPos)+(max*mEndXPos))>(fBegXPos-(NDScale*0.5*StripXWidth)) && ((max*pBegXPos)+(min*mBegXPos))<(fEndXPos+(NDScale*0.5*StripYWidth)) ) ||
	    TMath::Abs( (clustp->GetEndYPos()-fBegYPos)-(fEndYPos-clustm->GetBegYPos()) )<(NDScale*2.1*StripYWidth) ||
	    ( ((min*pEndYPos)+(max*mEndYPos))>(fBegYPos-(NDScale*0.5*StripYWidth)) && ((max*pBegYPos)+(min*mBegYPos))<(fEndYPos+(NDScale*0.5*StripYWidth)) ))
          {TrkAssocNum=2;}
      }
      
      // Overlapping clusters
      if( (clustm->GetBegXPos()<-100 || fEndXPos<-100 || (clustm->GetBegXPos()-fEndXPos)<(NDScale*1.1*StripXWidth)) && 
	  (fBegXPos<-100 || clustp->GetEndXPos()<-100 || (fBegXPos-clustp->GetEndXPos())<(NDScale*1.1*StripXWidth)) && 
	  (clustm->GetBegYPos()<-100 || fEndYPos<-100 || (clustm->GetBegYPos()-fEndYPos)<(NDScale*1.1*StripYWidth)) && 
	  (fBegYPos<-100 || clustp->GetEndYPos()<-100 || (fBegYPos-clustp->GetEndYPos())<(NDScale*1.1*StripYWidth))) {
        if(TrkAssocNum<1) TrkAssocNum=1;
      }
    }
    */
  }

  //  cout <<"xxpostyyyy "<< mXPos<<" "<< fXPos<<" "<<pXPos<<" "<<mYPos<<" "<< mYPos<<" "<<pYPos<<endl;
  
  return TrkAssocNum;
}

InoHit* InoCluster::GetHit(unsigned int ij) const {
  if(ij<HitsInCluster.size()) {return HitsInCluster[ij];}
  else {return 0;}
}

int InoCluster::IsDiffuseShwAssoc(InoCluster* clr) const {
  double win = 99.9;
  int assoc = 0;
  if(fEndTime-clr->GetBegTime()>-win || clr->GetEndTime()-fBegTime>-win){
    if( clr->GetZPlane()-fZPlane<9 && clr->GetZPlane()-fZPlane>-9 && 
	clr->GetEndXStrip()-fBegXStrip>-21 && fEndXStrip-clr->GetBegXStrip()>-21 &&
	clr->GetEndYStrip()-fBegYStrip>-21 && fEndYStrip-clr->GetBegYStrip()>-21 ) {
      if( ( clr->GetZPlane()-fZPlane<5 && clr->GetZPlane()-fZPlane>-5 && 
	    clr->GetEndXStrip()-fBegXStrip>-11 && fEndXStrip-clr->GetBegXStrip()>-11 && 
	    clr->GetEndYStrip()-fBegYStrip>-11 && fEndYStrip-clr->GetBegYStrip()>-11) ) assoc=2;
      else assoc=1;
    }
  }
  return assoc;
}

bool InoCluster::isIdentical(InoCluster* icls) {
  if ((GetHitEntries()) == (icls->GetHitEntries())) {
    for (unsigned ij=0; ij<GetHitEntries(); ij++) {
      if ( !GetHit(ij)->isIdentical( icls->GetHit(ij))) return false;
    }
  } else { return false ;}
  return true;
}

void InoCluster::Print() {
  cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
  cout<<"Cluster Info :::"<<endl;
  cout<<"Plane of Cluster = "<<fZPlane<<" RPCmod "<<fRPCmod<<endl;
  cout<<"fView = "<<fView<<" Cluster Size "<<HitsInCluster.size()<<endl;
  cout<<"Begin X Strip = "<<fBegXStrip<<" End X Strip = "<<fEndXStrip<<endl;
  cout<<"Begin Y Strip = "<<fBegYStrip<<" End Y Strip = "<<fEndYStrip<<endl;
  cout<<"Begin Time = "<<fBegTime<<" End Time = "<<fEndTime<<endl;
  cout<<"Begin X Pos = ("<<fBegXPos<<","<<fBegYPos<<")"<<endl;
  cout<<"End X Pos = ("<<fEndXPos<<","<<fEndYPos<<")"<<endl;
  for(unsigned int ij=0; ij<HitsInCluster.size(); ij++) {
    HitsInCluster[ij]->Print();
  }
  cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
}
