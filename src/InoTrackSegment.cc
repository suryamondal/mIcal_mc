#include "InoTrackSegment.h"
#include "InoCluster.h"
#include "InoHit.h"
#include "math.h"
#include <iostream>

InoTrackSegment::InoTrackSegment(InoCluster* clustm, InoCluster* clust0, InoCluster* clustp):
  fSeedSegment(0), fPartner(0), fUID(0), fBegZPlane(999), fEndZPlane(-999), fBegVtxZ(999.), 
  fEndVtxZ(-999.), fTrkFlag(0), fTmpTrkFlag(0), fNPlanes(0), 
  fBegTime(0.), fEndTime(0.), StripWidth(4.108e-2)
{
  //  fPlaneView=clust0->GetPlaneView();
  fBegTime=clust0->GetBegTime();
  fEndTime=clust0->GetEndTime();

  this->AddCluster(clustm); this->AddCluster(clust0); this->AddCluster(clustp);
}

InoTrackSegment::~InoTrackSegment() {
  ClustersInSegment.clear(); 
  fBegAssociatedSegList.clear();
  fEndAssociatedSegList.clear();
  fBegPreferredSegList.clear();
  fEndPreferredSegList.clear();
  fBegMatchedSegList.clear();
  fEndMatchedSegList.clear();
  //  delete fSeedSegment;
  //  delete fPartner;
}

void InoTrackSegment::AddCluster(InoCluster* clust) {
  if(clust) {
    if(this->ContainsCluster(clust)==false) {
      ClustersInSegment.push_back(clust);
      if( fBegZPlane > clust->GetZPlane()){ fBegZPlane = clust->GetZPlane(); fBegVtxZ =  clust->GetZPos(); }
      if( fEndZPlane < clust->GetZPlane()){ fEndZPlane = clust->GetZPlane(); fEndVtxZ =  clust->GetZPos(); }
      
      if( clust->GetBegTime()<fBegTime ) { fBegTime=clust->GetBegTime(); }
      if( clust->GetEndTime()>fEndTime ) { fEndTime=clust->GetEndTime(); }
    }
  }  
  return;
}

bool InoTrackSegment::ContainsCluster(InoCluster* clust) {
  for(unsigned int ij=0; ij<ClustersInSegment.size(); ++ij) {
    if(clust==ClustersInSegment[ij]) {return true;}
  }
  return false;
}

InoCluster* InoTrackSegment::GetCluster(unsigned int ij) {
  if(ij<ClustersInSegment.size()) {
    return ClustersInSegment[ij];
  } else {return 0;}
}

void InoTrackSegment::AddSegment(InoTrackSegment* segment) {
  for(unsigned int ij=0; ij<segment->GetEntries(); ++ij) {
    this->AddCluster(segment->GetCluster(ij));    
  }
}

unsigned int InoTrackSegment::GetEntries() const {
  return ClustersInSegment.size();
}

int InoTrackSegment::GetBegZPlane() const {
  return fBegZPlane;
}


int InoTrackSegment::GetEndZPlane() const
{
  return fEndZPlane;
}


bool InoTrackSegment::IsAssoc(InoTrackSegment* segment) {
  unsigned int ij;
  bool assoc=false;
  bool flag=false; 
  InoTrackSegment* Seg1 = this;
  InoTrackSegment* Seg2 = segment;
  //------------------------------
  //If the two clusters overlap at all
  if(Seg1->GetEndZPlane()>=Seg2->GetBegZPlane()){
    flag=true;
    //All clusters in Seg1 in planes overlapping Seg2 should also be in Seg2
    for(ij=0;ij<Seg1->GetEntries();ij++){
      InoCluster* clr = Seg1->GetCluster(ij);
      if(clr->GetZPlane()>=Seg2->GetBegZPlane()){ 
        if(!(Seg2->ContainsCluster(clr))){
          flag=false; 
          break; //one mismatch is too many
        }
      }
    }

    if(flag) { //only bother with 2nd check if 1st was ok
      //All clusters in Seg2 in planes overlapping Seg1 should also be in Seg1
      for(ij=0;ij<Seg2->GetEntries();ij++){
	InoCluster* clr = Seg2->GetCluster(ij);
	if(clr->GetZPlane()<=Seg1->GetEndZPlane()){ 
	  if(!(Seg1->ContainsCluster(clr))){
	    flag=false; 
	    break; //one mismatch is too many 
	  }
	}
      }
    }
    //If this is true and they overlap by more than one plane then they are associated
    if(Seg1->GetEndZPlane()>Seg2->GetBegZPlane()){ if(flag) assoc=true;}
  }
  //----------------------------
  //If the segments overlap by 1 or fewer planes need to do some more work...    
  if(Seg1->GetEndZPlane()<=Seg2->GetBegZPlane()){
    int idb=0,idb0=0;
    int bpln=Seg2->GetEndZPlane()+1;
    //find the first (lowest plane no) two clusters in the second segment
    for(ij=0;ij<Seg2->GetEntries();ij++) {
      InoCluster* clr = Seg2->GetCluster(ij);
      if(clr->GetZPlane()<bpln && clr->GetZPlane()>Seg2->GetBegZPlane()+2){
        bpln=clr->GetZPlane(); idb=ij;
      } else { 
	if(clr->GetZPlane()==Seg2->GetBegZPlane()) idb0=ij; 
      }
    }
    InoCluster* clrb = Seg2->GetCluster(idb); InoCluster* clrb0 = Seg2->GetCluster(idb0);
    double bBegXPos = clrb->GetBegXPos(); double bEndXPos = clrb->GetEndXPos();
    double b0BegXPos = clrb0->GetBegXPos(); double b0EndXPos = clrb0->GetEndXPos();

    double bBegYPos = clrb->GetBegYPos(); double bEndYPos = clrb->GetEndYPos();
    double b0BegYPos = clrb0->GetBegYPos(); double b0EndYPos = clrb0->GetEndYPos();

    int ide=0,ide0=0;
    int epln=Seg1->GetBegZPlane()-1;
    //find the last (highest plane no) two clusters in the first segment
    for(ij=0;ij<Seg1->GetEntries();ij++){
      InoCluster* clr = Seg1->GetCluster(ij);
      if(clr->GetZPlane()>epln && clr->GetZPlane()<Seg1->GetEndZPlane()-2){
        epln=clr->GetZPlane(); ide=ij;
      } else { 
	if(clr->GetZPlane()==Seg1->GetEndZPlane()) ide0=ij; 
      }
    }
    InoCluster* clre = Seg1->GetCluster(ide); InoCluster* clre0 = Seg1->GetCluster(ide0);
    double eBegXPos = clre->GetBegXPos(); double eEndXPos = clre->GetEndXPos();
    double e0BegXPos = clre0->GetBegXPos(); double e0EndXPos = clre0->GetEndXPos();

    double eBegYPos = clre->GetBegYPos(); double eEndYPos = clre->GetEndYPos();
    double e0BegYPos = clre0->GetBegYPos(); double e0EndYPos = clre0->GetEndYPos();

    //Look at how these clusters overlap stripwise    
    //GMA need optimisation 01/07/2008 (move from 2* to 3*)
    if ( ( ( bEndXPos-b0BegXPos>-2*StripWidth && b0EndXPos-bBegXPos>-2*StripWidth && 
	     eEndXPos-e0BegXPos>-2*StripWidth && e0EndXPos-eBegXPos>-2*StripWidth ) || 
	   ( ( bBegXPos-b0BegXPos>-2*StripWidth && e0BegXPos-eBegXPos>-2*StripWidth ) || 
	     ( bEndXPos-b0EndXPos>-2*StripWidth && e0EndXPos-eEndXPos>-2*StripWidth ) ) || 
	   ( ( bBegXPos-b0BegXPos<2*StripWidth && e0BegXPos-eBegXPos<2*StripWidth ) || 
	     ( bEndXPos-b0EndXPos<2*StripWidth && e0EndXPos-eEndXPos<2*StripWidth ) ) ) &&
	 
	 ( ( bEndYPos-b0BegYPos>-2*StripWidth && b0EndYPos-bBegYPos>-2*StripWidth && 
	     eEndYPos-e0BegYPos>-2*StripWidth && e0EndYPos-eBegYPos>-2*StripWidth ) || 
	   ( ( bBegYPos-b0BegYPos>-2*StripWidth && e0BegYPos-eBegYPos>-2*StripWidth ) || 
	     ( bEndYPos-b0EndYPos>-2*StripWidth && e0EndYPos-eEndYPos>-2*StripWidth ) ) || 
	   ( ( bBegYPos-b0BegYPos<2*StripWidth && e0BegYPos-eBegYPos<2*StripWidth ) || 
	     ( bEndYPos-b0EndYPos<2*StripWidth && e0EndYPos-eEndYPos<2*StripWidth ) ) )) {
      //if the segments overlap by 1 plane then we can now say if they are associated or not
      if(Seg2->GetBegZPlane()==Seg1->GetEndZPlane()){ if(flag) assoc=true; 
      } else { 
        double z1,z2,x1,x2,y1,y2,dx1,dx2, dy1,dy2,dirtmp1,dirtmp2,win;
	double dir1x, dir1y, dir2x, dir2y;

        z1=Seg1->GetEndZPos();  z2=Seg2->GetBegZPos();

        x1=Seg1->GetEndXPos();  x2=Seg2->GetBegXPos();
	y1=Seg1->GetEndYPos();  y2=Seg2->GetBegYPos();

	dir1x = Seg1->GetEndXDir(); dir1y = Seg1->GetEndYDir();
	//        dir1=Seg1->GetEndDir(); 
	dx1=x2-x1-dir1x*(z2-z1); dy1=y2-y1-dir1y*(z2-z1);    
	
        dir2x=Seg2->GetBegXDir(); dir2y=Seg2->GetBegYDir();
	//        dir2=Seg2->GetBegDir(); 
	dx2=x2-x1-dir2x*(z2-z1); dy2=y2-y1-dir2y*(z2-z1);

        win=0.1+0.35*(z2-z1); //GMA need correction   // we need to set this window in order to take care of bending of track 
	dirtmp1=(1+dir1x*dir2x)/(sqrt(1+dir1x*dir1x)*sqrt(1+dir2x*dir2x));
        dirtmp2=(1+dir1y*dir2y)/(sqrt(1+dir1y*dir1y)*sqrt(1+dir2y*dir2y));
	
        if( ( (dx1>-win&&dx1<win)||(dx2>-win&&dx2<win) ) &&
	    ( (dy1>-win&&dy1<win)||(dy2>-win&&dy2<win) ) &&
	    ( dirtmp1>0.65 ) && ( dirtmp2>0.65 ) ) {assoc=true; 
        }
      }
    }
  }
  //-------------------------------
  return assoc;
}


double InoTrackSegment::GetBegXPos() {
  double tot=0.0,begt=0.0;
  //Loop over clusters in segment
  unsigned int nclusters = this->GetEntries();
  unsigned int nhits = 0;
  InoHit* hit=0;
  int ifirst = 0;
  while ( tot==0 && ifirst <4) {
    for( unsigned int ij=0; ij<nclusters; ++ij) {
      //find the clusters on the first plane in the segment
      //    int ifirst = 0;
      //    while ( tot==0 && ifirst <4) {
      if(ClustersInSegment[ij]->GetZPlane()==fBegZPlane+ifirst) {
	nhits = ClustersInSegment[ij]->GetHitEntries();
	
	//loop over hits in cluster
	for(unsigned int jk=0;jk<nhits;++jk){
	  hit = ClustersInSegment[ij]->GetHit(jk);
	  if(hit) {
	    if (hit->GetXPosErr() <100) {
	      begt+=hit->GetXPos();
	      tot+=1.0;
	    }
	  }
	}
      }
    }
    if(tot>0) return (begt/tot); 
    ifirst++;
  }
  return 0;
}

double InoTrackSegment::GetEndXPos() {
  double tot=0.0,endt=0.0;
  //Loop over clusters in segment
  unsigned int nclusters = this->GetEntries();
  unsigned int nhits = 0;
  InoHit* hit=0;
  int ifirst = 0;
  while ( tot==0 && ifirst <4) {
  for( unsigned int ij=0; ij<nclusters; ++ij) {
    //find the clusters on the first plane in the segment
    //    int ifirst = 0;
    //    while ( tot==0 && ifirst <4) {
    if(ClustersInSegment[ij]->GetZPlane()==fEndZPlane-ifirst) {
      nhits = ClustersInSegment[ij]->GetHitEntries();
      
      //loop over hits in cluster
      for(unsigned int jk=0;jk<nhits;++jk) {
	hit = ClustersInSegment[ij]->GetHit(jk);
	if (hit) {
	  if (hit->GetXPosErr() <100) {
	    endt+=hit->GetXPos();
	    tot+=1.0;
	  }
	}
      }
    }
  }
  if(tot>0) return (endt/tot); 
  ifirst++;
  }
  return 0;
}

double InoTrackSegment::GetBegYPos() {
  double tot=0.0,begt=0.0;
  //Loop over clusters in segment
  unsigned int nclusters = this->GetEntries();
  unsigned int nhits = 0;
  InoHit* hit=0;
  int ifirst = 0;
  while ( tot==0 && ifirst <4) {
    for( unsigned int ij=0; ij<nclusters; ++ij) {
      //find the clusters on the first plane in the segment
      //      int ifirst = 0;
      //    while ( tot==0 && ifirst <4) {
      if(ClustersInSegment[ij]->GetZPlane()==fBegZPlane+ifirst) {
	nhits = ClustersInSegment[ij]->GetHitEntries();
	
	//loop over hits in cluster
	for(unsigned int jk=0;jk<nhits;++jk){
	  hit = ClustersInSegment[ij]->GetHit(jk);
	  if(hit) {
	    if (hit->GetYPosErr() <100) {
	      begt+=hit->GetYPos();
	      tot+=1.0;
	    }
	  }
	}
      }
    }
    if(tot>0) return (begt/tot); 
    ifirst++;
  }
  return 0;
}

double InoTrackSegment::GetEndYPos() {
  double tot=0.0,endt=0.0;
  //Loop over clusters in segment
  unsigned int nclusters = this->GetEntries();
  unsigned int nhits = 0;
  InoHit* hit=0;
  int ifirst = 0;
  while ( tot==0 && ifirst <4) { 
    for( unsigned int ij=0; ij<nclusters; ++ij) {
      //find the clusters on the first plane in the segment
      //      int ifirst = 0;
      //      while ( tot==0 && ifirst <4) {
      if(ClustersInSegment[ij]->GetZPlane()==fEndZPlane-ifirst) {
	nhits = ClustersInSegment[ij]->GetHitEntries();
	
	//loop over hits in cluster
	for(unsigned int jk=0;jk<nhits;++jk) {
	  hit = ClustersInSegment[ij]->GetHit(jk);
	  if (hit) {
	    if (hit->GetYPosErr() <100) {
	      endt+=hit->GetYPos();
	      tot+=1.0;
	    }
	  }
	}
      }
    }
    if(tot>0) return (endt/tot); 
    ifirst++;
  }
  return 0;
}



/* 
double* InoTrackSegment::GetBegDir() {
  double z=0.0,t=0.0;
  double sw=0.0,swx=0.0,swx2=0.0;
  double swy=0.0,swyx=0.0;

  double zyy=0.0,tyy=0.0;
  double swyy=0.0,swxyy=0.0,swx2yy=0.0;
  double swyyy=0.0,swyxyy=0.0;

  unsigned int nclusters = this->GetEntries();
  unsigned int nhits = 0;
  InoHit* hit=0;
  double* angle=0;
  for( unsigned int ij=0; ij<nclusters; ++ij) {
    //find the clusters on the first plane in the segment
    if(ClustersInSegment[ij]->GetZPlane()<fBegZPlane+10) {
      nhits = ClustersInSegment[ij]->GetEntries();
      
      //loop over hits in cluster
      for(unsigned int jk=0;jk<nhits;++jk) {
	hit = ClustersInSegment[ij]->GetHit(jk);
	if(hit) {
	  if (hit->GetXPosErr()<100) {
	    z=hit->GetZPos(); t=hit->GetXPos();
	    sw+=1.0; swx+=z; swx2+=z*z; 
	    swy+=t; swyx+=t*z;
	  }
	  
	  if (hit->GetYPosErr()<100) {
	    zyy=hit->GetZPos(); tyy=hit->GetYPos();
	    swyy+=1.0; swxyy+=zyy; swx2yy+=zyy*zyy; 
	    swyyy+=tyy; swyxyy+=tyy*zyy;
	  } 
	}
      }
    }
  }
  if((swx*swx-sw*swx2)!=0) {angle[0]= (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
  if((swxyy*swxyy-swyy*swx2yy)!=0) {angle[1]= (swxyy*swyyy-swyy*swyxyy)/(swxyy*swxyy-swyy*swx2yy);} 

  return angle;
}
*/

double InoTrackSegment::GetBegXDir() {
  double z=0.0,t=0.0;
  double sw=0.0,swx=0.0,swx2=0.0;
  double swy=0.0,swyx=0.0;

  unsigned int nclusters = this->GetEntries();
  unsigned int nhits = 0;
  InoHit* hit=0;
  double angle=0;
  for( unsigned int ij=0; ij<nclusters; ++ij) {
    //find the clusters on the first plane in the segment
    if(ClustersInSegment[ij]->GetZPlane()<fBegZPlane+5) {  //31Aug  10->5
      nhits = ClustersInSegment[ij]->GetHitEntries();
      
      //loop over hits in cluster
      for(unsigned int jk=0;jk<nhits;++jk) {
	hit = ClustersInSegment[ij]->GetHit(jk);
	if(hit) {
	  if (hit->GetXPosErr()<100) {
	    z=hit->GetZPos(); t=hit->GetXPos();
	    sw+=1.0; swx+=z; swx2+=z*z; 
	    swy+=t; swyx+=t*z;
	  }
	  
	}
      }
    }
  }
  if((swx*swx-sw*swx2)!=0) {angle= (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
  return angle;
}


double InoTrackSegment::GetBegYDir() {
  double zyy=0.0,tyy=0.0;
  double swyy=0.0,swxyy=0.0,swx2yy=0.0;
  double swyyy=0.0,swyxyy=0.0;

  unsigned int nclusters = this->GetEntries();
  unsigned int nhits = 0;
  InoHit* hit=0;
  double angle=0;
  for( unsigned int ij=0; ij<nclusters; ++ij) {
    //find the clusters on the first plane in the segment
    if(ClustersInSegment[ij]->GetZPlane()<fBegZPlane+5) {   //31Aug 10->5
      nhits = ClustersInSegment[ij]->GetHitEntries();
      
      //loop over hits in cluster
      for(unsigned int jk=0;jk<nhits;++jk) {
	hit = ClustersInSegment[ij]->GetHit(jk);
	if(hit) {
	  if (hit->GetYPosErr()<100) {
	    zyy=hit->GetZPos(); tyy=hit->GetYPos();
	    swyy+=1.0; swxyy+=zyy; swx2yy+=zyy*zyy; 
	    swyyy+=tyy; swyxyy+=tyy*zyy;
	  } 
	}
      }
    }
  }
  if((swxyy*swxyy-swyy*swx2yy)!=0) {angle= (swxyy*swyyy-swyy*swyxyy)/(swxyy*swxyy-swyy*swx2yy);} 

  return angle;
}

/* 
double* InoTrackSegment::GetEndDir() {
  double z=0.0,t=0.0;
  double sw=0.0,swx=0.0,swx2=0.0;
  double swy=0.0,swyx=0.0;

  double zyy=0.0,tyy=0.0;
  double swyy=0.0,swxyy=0.0,swx2yy=0.0;
  double swyyy=0.0,swyxyy=0.0;

  unsigned int nclusters = this->GetEntries();
  unsigned int nhits = 0;
  InoHit* hit=0;
  double* angle=0;
  for( unsigned int ij=0; ij<nclusters; ++ij) {
    //find the clusters on the first plane in the segment
    if(ClustersInSegment[ij]->GetZPlane()>fEndZPlane-10) {
      nhits = ClustersInSegment[ij]->GetEntries();
      
      //loop over hits in cluster
      for(unsigned int jk=0;jk<nhits;++jk) {
	hit = ClustersInSegment[ij]->GetHit(jk);
	if(hit) {
	  if (hit->GetXPosErr()<100) {
	    z=hit->GetZPos(); t=hit->GetXPos();
	    sw+=1.0; swx+=z; swx2+=z*z; 
	    swy+=t; swyx+=t*z;
	  }
	  
	  if (hit->GetYPosErr()<100) {
	    zyy=hit->GetZPos(); tyy=hit->GetYPos();
	    swyy+=1.0; swxyy+=zyy; swx2yy+=zyy*zyy; 
	    swyyy+=tyy; swyxyy+=tyy*zyy;
	  } 
	}
      }
    }
  }
  if((swx*swx-sw*swx2)!=0) {angle[0]= (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
  if((swxyy*swxyy-swyy*swx2yy)!=0) {angle[1]= (swxyy*swyyy-swyy*swyxyy)/(swxyy*swxyy-swyy*swx2yy);} 

  return angle;
}
*/

double InoTrackSegment::GetEndXDir()
{
  double z=0.0,t=0.0;
  double sw=0.0,swx=0.0,swx2=0.0;
  double swy=0.0,swyx=0.0;

  unsigned int nclusters = this->GetEntries();
  unsigned int nhits = 0;
  InoHit* hit=0;
  double angle(0);
  for( unsigned int ij=0; ij<nclusters; ++ij) {
    //find the clusters on the first plane in the segment
    if(ClustersInSegment[ij]->GetZPlane()>fEndZPlane-5) {   //31Aug  10->5
      nhits = ClustersInSegment[ij]->GetHitEntries();
      
      //loop over hits in cluster
      for(unsigned int jk=0;jk<nhits;++jk) {
	hit = ClustersInSegment[ij]->GetHit(jk);
	if(hit) {
	  if (hit->GetXPosErr()<100) {
	    z=hit->GetZPos(); t=hit->GetXPos();
	    sw+=1.0; swx+=z; swx2+=z*z; 
	    swy+=t; swyx+=t*z;
	  }
	}
      }
    }
  }
  if((swx*swx-sw*swx2)!=0) {angle= (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}

  return angle;
}

double InoTrackSegment::GetEndYDir()
{

  double zyy=0.0,tyy=0.0;
  double swyy=0.0,swxyy=0.0,swx2yy=0.0;
  double swyyy=0.0,swyxyy=0.0;

  unsigned int nclusters = this->GetEntries();
  unsigned int nhits = 0;
  InoHit* hit=0;
  double angle(0);
  for( unsigned int ij=0; ij<nclusters; ++ij) {
    //find the clusters on the first plane in the segment
    if(ClustersInSegment[ij]->GetZPlane()>fEndZPlane-5) { // 31Aug 10->5
      nhits = ClustersInSegment[ij]->GetHitEntries();
      
      //loop over hits in cluster
      for(unsigned int jk=0;jk<nhits;++jk) {
	hit = ClustersInSegment[ij]->GetHit(jk);
	if(hit) {
	  
	  if (hit->GetYPosErr()<100) {
	    zyy=hit->GetZPos(); tyy=hit->GetYPos();
	    swyy+=1.0; swxyy+=zyy; swx2yy+=zyy*zyy; 
	    swyyy+=tyy; swyxyy+=tyy*zyy;
	  } 
	}
      }
    }
  }
  if((swxyy*swxyy-swyy*swx2yy)!=0) {angle= (swxyy*swyyy-swyy*swyxyy)/(swxyy*swxyy-swyy*swx2yy);} 

  return angle;
}


double InoTrackSegment::GetBegZPos() const {
  return fBegVtxZ;
}

double InoTrackSegment::GetEndZPos() const {
  return fEndVtxZ;
}

void InoTrackSegment::AddAssocSegToBeg(InoTrackSegment* seg)  { //void AddTriToBeg(InoTrackSegment* seg); 
  fBegAssociatedSegList.push_back(seg);
  return;
}

void InoTrackSegment::AddAssocSegToEnd(InoTrackSegment* seg) { //void AddTriToEnd(InoTrackSegment* seg);
  fEndAssociatedSegList.push_back(seg);
  return;
}

InoTrackSegment* InoTrackSegment::GetAssocSegBeg(unsigned int ij)  { //InoTrackSegment* GetTriBeg(unsigned int ij);
  if(ij<fBegAssociatedSegList.size()) {return fBegAssociatedSegList[ij];}
  else return 0;
}

InoTrackSegment* InoTrackSegment::GetAssocSegEnd(unsigned int ij) { //InoTrackSegment* GetTriEnd(unsigned int ij);
  if(ij<fEndAssociatedSegList.size()) {return fEndAssociatedSegList[ij];}
  else return 0;
}

void InoTrackSegment::AddPrefSegToBeg(InoTrackSegment* seg) { //void AddAssocTriToBeg(TrackSegmentAtNu* segment);
  fBegPreferredSegList.push_back(seg);
  return;
}

void InoTrackSegment::AddPrefSegToEnd(InoTrackSegment* seg) { //void AddAssocTriToEnd(TrackSegmentAtNu* segment);
  fEndPreferredSegList.push_back(seg);
  return;
}

InoTrackSegment* InoTrackSegment::GetPrefSegBeg(unsigned int ij) { //TrackSegmentAtNu* GetAssocTriBegAt(Int_t ij);
  if(ij<fBegPreferredSegList.size()) {return fBegPreferredSegList[ij];}
  else return 0;
}

InoTrackSegment* InoTrackSegment::GetPrefSegEnd(unsigned int ij) { //TrackSegmentAtNu* GetAssocTriEndAt(Int_t ij);
  if(ij<fEndPreferredSegList.size())  {return fEndPreferredSegList[ij];}
  else return 0;
}

void InoTrackSegment::AddMatchSegToBeg(InoTrackSegment* seg) { //void AddSegToBeg(TrackSegmentAtNu* segment);
  fBegMatchedSegList.push_back(seg);
  return;
}


void InoTrackSegment::AddMatchSegToEnd(InoTrackSegment* seg) { //void AddSegToEnd(TrackSegmentAtNu* segment)
  fEndMatchedSegList.push_back(seg);
  return;
}

InoTrackSegment* InoTrackSegment::GetMatchSegBeg(unsigned int ij) { //TrackSegmentAtNu* GetSegBegAt(Int_t ij)
  if(ij<fBegMatchedSegList.size()) {return fBegMatchedSegList[ij];}
  else return 0;
}


InoTrackSegment* InoTrackSegment::GetMatchSegEnd(unsigned int ij) { //TrackSegmentAtNu* GetSegEndAt(Int_t ij)
  if(ij<fEndMatchedSegList.size()) {return fEndMatchedSegList[ij];}
  else return 0;
}


//int InoTrackSegment::GetPlaneView() const{  return fPlaneView;}

//GMA  Check this straightness (Here only X-axis)
double InoTrackSegment::GetScore(vector<InoTrackSegment*> *BegSegBank, vector<InoTrackSegment*> *EndSegBank) {
  int begplane2, begplane1, endplane1, endplane2;
  int nplane, plane;
  
  vector<InoCluster*> TempContainer;
  bool IsInTemp;
  
  //Store the beginning and end of the segment locally
  begplane1=this->GetBegZPlane();
  endplane1=this->GetEndZPlane();
  begplane2=this->GetBegZPlane();
  endplane2=this->GetEndZPlane();
  //
  //  begplane2============begplane1===========endplane1============endplane2
  //            BegSegBank             this              EndSegBank
  //

  // Store pointers to all clusters
  unsigned int nclusters = this->GetEntries();  
  for(unsigned int ij=0; ij<nclusters; ++ij) {
    TempContainer.push_back(this->GetCluster(ij));
  }

  unsigned int nsegments;
  unsigned int nclustersintemp;
  int tempplane;
  if(BegSegBank) {
    nsegments = (*BegSegBank).size();
    for(unsigned int jk=0; jk<nsegments; ++jk) {
      //looking for the earliest plane in the Beginning Segments
      tempplane = (*BegSegBank)[jk]->GetBegZPlane();
      if(tempplane<begplane2) begplane2=tempplane;
      //Loop over clusters in this segment and add any we don't have
      nclusters=(*BegSegBank)[jk]->GetEntries();
      for(unsigned int kl=0; kl<nclusters; ++kl) {
        InoCluster* clust = (*BegSegBank)[jk]->GetCluster(kl);
        
        IsInTemp=false;
        nclustersintemp = TempContainer.size();
        for(unsigned int lm=0; lm<nclustersintemp; ++lm) {
          if(TempContainer[lm]==clust) {IsInTemp=true; break;}
        }
        if(IsInTemp==false) {TempContainer.push_back(clust);}
      }
    }
  }

  if(EndSegBank) {
    nsegments = (*EndSegBank).size();
    for(unsigned int jk=0; jk<nsegments; ++jk) {
      //looking for the latest plane in the End Segments
      tempplane = (*EndSegBank)[jk]->GetEndZPlane();
      if(tempplane>endplane2) endplane2=tempplane;
      //Loop over clusters in this segment and add any we don't have
      nclusters = (*EndSegBank)[jk]->GetEntries();
      for(unsigned int kl=0; kl<nclusters; ++kl) {
        InoCluster* clust = (*EndSegBank)[jk]->GetCluster(kl);
        
        IsInTemp=false;
        nclustersintemp = TempContainer.size();
        for(unsigned int lm=0; lm<nclustersintemp; ++lm) {
          if(TempContainer[lm]==clust) {IsInTemp=true; break;}
        }
        if(IsInTemp==false) {TempContainer.push_back(clust);}
      }
    }
  }

  //Convert planes to a single view co-ordinate system where begplane2 = 0;
  begplane1-=begplane2; begplane1/=2;
  endplane1-=begplane2; endplane1/=2;
  endplane2-=begplane2; endplane2/=2;

  //calculate the number of planes we will be considering
  nplane=1+endplane2;
  double* TX = new double[nplane];
  double* ZX = new double[nplane];
  double* WX = new double[nplane];
  for(int ij=0; ij<nplane; ++ij) {TX[ij]=0.; ZX[ij]=0.; WX[ij]=0.;}

  int kl, kp;  
  double am,c;
  double dt2,sn;
  double score, dstraightness, straightness, expected;
  double sw, swz, swt, swzt, swzz;
  unsigned int nhits;

  nclusters = TempContainer.size();

  for(unsigned int ij=0; ij<nclusters; ++ij) {
    InoCluster* clust = TempContainer[ij];
    //Convert clust->GetZPlane() to a single view co-ordinate system where begplane2 = 0;
    plane=(clust->GetZPlane()-begplane2)/2;

    if(!(plane<0 || plane>=nplane)) {
      sw=0.; swz=0.; swt=0.;
      nhits = clust->GetHitEntries();
      for(unsigned int k1=0; k1<nhits; ++k1) {
        InoHit* hit = clust->GetHit(k1);
        
        swz+=hit->GetPulse()*hit->GetZPos();
        swt+=hit->GetPulse()*hit->GetXPos();
        sw+=hit->GetPulse();
      }
      
      if(sw>0.){
        ZX[plane]=swz/sw; TX[plane]=swt/sw;
        // Weight segments on planes spanned by seed segment. 
        // Deweight segments on other planes.
        if( plane+1>begplane1 && plane-1<endplane1 ) {WX[plane]=5.;} 
        else {WX[plane]=0.5;}
      }
    }
  }
    
  score=0.; straightness=0.; expected=0.;
  
  for(int ij=0; ij<nplane; ++ij) {
    
    if(WX[ij]>0.){
      swz=0.; swt=0.; swzz=0.; swzt=0.; sw=0.; sn=0.;
      dstraightness=0.;
      
      kl=ij-5; kp=ij+5;
            
      if(kl<0) {kl=0;}
      if(kp>nplane-1) {kp=nplane-1;}
      

      // Fit this section
      for(int k1=kl; k1<kp+1; ++k1){
        if(WX[k1]>0.) {
          swz+=WX[k1]*ZX[k1]; swt+=WX[k1]*TX[k1]; 
          swzz+=WX[k1]*ZX[k1]*ZX[k1]; swzt+=WX[k1]*ZX[k1]*TX[k1]; 
          sw+=WX[k1]; sn+=1.;
        }
      }
      
      // Calculate deviation from fit at this plane
      if(sn>2.){
        am=(sw*swzt-swz*swt)/(sw*swzz-swz*swz);
        c=(swt*swzz-swz*swzt)/(sw*swzz-swz*swz);

        dt2=pow(TX[ij]-(am*ZX[ij]+c),2);

        if(dt2<1.e-5) {straightness+=1;}
        else {straightness+=dt2/1.e-5;}
        expected+=1;
      }
    }
  }

  // Protect against divide by zero
  if(expected==0) {expected=1; straightness=1;}

  score = 1.e4 + double(TempContainer.size()) - pow(straightness/expected,0.5);

  if(score<0.) {score=0.;}
 
  if(TX) {delete [] TX;}
  if(ZX) {delete [] ZX;}
  if(WX) {delete [] WX;}

  //For Y

  double* TY = new double[nplane];
  double* ZY = new double[nplane];
  double* WY = new double[nplane];
  for(int ij=0; ij<nplane; ++ij) {TY[ij]=0.; ZY[ij]=0.; WY[ij]=0.;}

  for(unsigned int ij=0; ij<nclusters; ++ij) {
    InoCluster* clust = TempContainer[ij];
    //Convert clust->GetZPlane() to a single view co-ordinate system where begplane2 = 0;
    plane=(clust->GetZPlane()-begplane2)/2;

    if(!(plane<0 || plane>=nplane)) {
      sw=0.; swz=0.; swt=0.;
      nhits = clust->GetHitEntries();
      for(unsigned int k1=0; k1<nhits; ++k1) {
        InoHit* hit = clust->GetHit(k1);
        
        swz+=hit->GetPulse()*hit->GetZPos();
        swt+=hit->GetPulse()*hit->GetYPos();
        sw+=hit->GetPulse();
      }
      
      if(sw>0.){
        ZY[plane]=swz/sw; TY[plane]=swt/sw;
        // Weight segments on planes spanned by seed segment. 
        // Deweight segments on other planes.
        if( plane+1>begplane1 && plane-1<endplane1 ) {WY[plane]=5.;} 
        else {WY[plane]=0.5;}
      }
    }
  }
  
  straightness=0.; expected=0.;

  for(int ij=0; ij<nplane; ++ij) {

    if(WY[ij]>0.){
      swz=0.; swt=0.; swzz=0.; swzt=0.; sw=0.; sn=0.;
      dstraightness=0.;

      kl=ij-5; kp=ij+5;
            
      if(kl<0) {kl=0;}
      if(kp>nplane-1) {kp=nplane-1;}
      

      // Fit this section
      for(int k1=kl; k1<kp+1; ++k1){
        if(WY[k1]>0.) {
          swz+=WY[k1]*ZY[k1]; swt+=WY[k1]*TY[k1]; 
          swzz+=WY[k1]*ZY[k1]*ZY[k1]; swzt+=WY[k1]*ZY[k1]*TY[k1]; 
          sw+=WY[k1]; sn+=1.;
        }
      }
      
      // Calculate deviation from fit at this plane
      if(sn>2.){
        am=(sw*swzt-swz*swt)/(sw*swzz-swz*swz);
        c=(swt*swzz-swz*swzt)/(sw*swzz-swz*swz);

        dt2=pow(TY[ij]-(am*ZY[ij]+c),2);

        if(dt2<1.e-5) {straightness+=1;}
        else {straightness+=dt2/1.e-5;}
        expected+=1;
      }
    }
  }

  // Protect against divide by zero
  if(expected==0) {expected=1; straightness=1;}
  double score1 = 1.e4 + double(TempContainer.size()) - pow(straightness/expected,0.5);
  if (score1 <0.) {score1 = 0.;}
  score += score1;

  if(score<0.) {score=0.;}
 
  if(TY) {delete [] TY;}
  if(ZY) {delete [] ZY;}
  if(WY) {delete [] WY;}

  return score;
}
