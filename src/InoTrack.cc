#include "InoTrack.h"
#include "InoCluster.h"

#include "TMath.h"
#include "InoTrackSegment.h"
#include <iostream>

InoTrack::InoTrack() :
  fSegment(0), fBegZPlane(999), fEndZPlane(-999), fBegZ(999.), fEndZ(-999.),
  fUID(0), fUsed(-1)
  //, fPlaneView(-1)
{
  fSegment = 0;
}

InoTrack::InoTrack(InoTrackSegment* segment) :
  fSegment(0), fBegZPlane(999), fEndZPlane(-999), fBegZ(999.), fEndZ(-999.),
  fUID(0), fUsed(-1)
  //, fPlaneView(-1)
{
  fSegment = segment;
}


InoTrack::~InoTrack() {
  ClustsInTrack.clear();
}

void InoTrack::AddCluster(InoCluster* clust) {
  //  if(this->GetEntries()==0) { fPlaneView=clust->GetPlaneView(); }
  
  if(this->ContainsClust(clust)==true) {return;}
  ClustsInTrack.push_back(clust);   
  if( fBegZPlane > clust->GetZPlane()){ fBegZPlane = clust->GetZPlane(); fBegZ = clust->GetZPos();}
  if( fEndZPlane < clust->GetZPlane()){ fEndZPlane = clust->GetZPlane(); fEndZ = clust->GetZPos();}
  return;
}

//AAR this function is for joining 2 track 071211//Kolahal add this function into code
void InoTrack::AddTrack(InoTrack* trk)
{
  //  if(this->GetEntries()==0) { fPlaneView=clust->GetPlaneView(); }
  for(unsigned int ii=0;ii<trk->GetEntries();ii++){
    InoCluster * clust = trk->GetCluster(ii);
    this->AddCluster(clust);
  }
  return;
}

bool InoTrack::ContainsClust(InoCluster* clust) const {
  for(unsigned int i=0; i<ClustsInTrack.size(); ++i) {
    if(clust==ClustsInTrack[i]) {return true;}
  }
  return false;
}


InoCluster* InoTrack::GetCluster(unsigned int i) const {
  if(i<ClustsInTrack.size()) {return ClustsInTrack[i];}
  else {return 0;}
}

double InoTrack::GetBegXPos() {
  //  double tryx = 0.0;
  double tot=0.0,begt=0.0;
  unsigned int nclusts = ClustsInTrack.size();
  
  //     for(unsigned int kl = 0; kl < nclusts; kl++) {
  //       InoCluster* TestClust = ClustsInTrack[kl];
  //       if (!TestClust->GetStraight()) continue;
  //       if (kl==0)
  //       tryx = TestClust->GetXPos();
  //     }

  
  int ifirst=0;
  while(tot==0 && ifirst<int(nclusts-3)) {
    for(unsigned int ij=0; ij<nclusts; ++ij) {
      //find the clusts on the first plane in the track
      InoCluster* clust = ClustsInTrack[ij];
      if (!clust->GetStraight()) continue;
      
      if(clust->GetZPlane()==fBegZPlane+ifirst) {
	if (clust->GetXPosErr()<100) {
	  begt+=clust->GetXPos();
	  tot+=1.0;
	}
      }
    }
    if(tot>0) return (begt/tot);
    ifirst++;
  }
  
  return 0;
}

double InoTrack::GetEndXPos() {
  double tot=0.0,endt=0.0;
  unsigned int nclusts = ClustsInTrack.size();
  
  int ifirst=0;
  while(tot==0 && ifirst<int(nclusts-3)) {
    for(unsigned int ij=0; ij<nclusts; ++ij) {
      //find the clusts on the first plane in the track
      InoCluster* clust = ClustsInTrack[ij];
      if (!clust->GetStraight()) continue;
      if(clust->GetZPlane()==fEndZPlane-ifirst) {
	if (clust->GetXPosErr()<100) {
	  endt+=clust->GetXPos();
	  tot+=1.0;
	}
      }
    }
    if(tot>0) return (endt/tot); 
    ifirst++;
  }
  return 0;
}


double InoTrack::GetBegYPos() {
  double tot=0.0,begt=0.0;
  unsigned int nclusts = ClustsInTrack.size();
  
  int ifirst=0;
  while(tot==0 && ifirst<int(nclusts-3)) {
    for(unsigned int ij=0; ij<nclusts; ++ij) {
      //find the clusts on the first plane in the track
      InoCluster* clust = ClustsInTrack[ij];
      if (!clust->GetStraight()) continue;
      if(clust->GetZPlane()==fBegZPlane+ifirst) {
	if (clust->GetYPosErr()<100){
	  begt+=clust->GetYPos();
	  tot+=1.0;
	}
      }
    }
    if(tot>0) return (begt/tot); 
    ifirst++;
  }
  return 0;
}

double InoTrack::GetEndYPos() {
  double tot=0.0,endt=0.0;
  unsigned int nclusts = ClustsInTrack.size();
  
  int ifirst=0;
  while(tot==0 && ifirst<int(nclusts-3)) {
    for(unsigned int ij=0; ij<nclusts; ++ij) {
      //find the clusts on the first plane in the track
      InoCluster* clust = ClustsInTrack[ij];
      if (!clust->GetStraight()) continue;
      if(clust->GetZPlane()==fEndZPlane-ifirst) {
	if (clust->GetYPosErr()<100){
	  endt+=clust->GetYPos();
	  tot+=1.0;
	}
      }
    }
    if(tot>0) return (endt/tot); 
    ifirst++;
  }
  return 0;
}


/* 
double* InoTrack::GetDir(int Plane1, int Plane2) {
  
  unsigned int nclusts = ClustsInTrack.size();
  
  double* angle=0;
  
  for( unsigned int ij=0; ij<nclusts; ++ij) {
    InoCluster* clust = ClustsInTrack[ij];
    //find the clusts in the end region of the track
    int ifirst = 0;
    while (angle[0]==0 && ifirst <4) {
      double z=0.0,t=0.0;
      double sw=0.0,swx=0.0,swx2=0.0;
      double swy=0.0,swyx=0.0;
      
      if(clust->GetZPlane()>=Plane1+int(ifirst/2) && clust->GetZPlane()<=Plane2-int((ifirst+1)/2)) {
	if (clust->GetXPosErr()<100){
	  z=clust->GetZPos(); 
	  t=clust->GetXPos();
	  sw+=1.0; 
	  swx+=z; 
	  swx2+=z*z; 
	  swy+=t; 
	  swyx+=t*z;
	}
      }
      if(TMath::Abs(swx*swx-sw*swx2)>1e-10) {angle[0]= (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
      ifirst++;
    }

    ifirst = 0;
    while (angle[1]==0 && ifirst <4) {
      double z=0.0,t=0.0;
      double sw=0.0,swx=0.0,swx2=0.0;
      double swy=0.0,swyx=0.0;
      
      if(clust->GetZPlane()>=Plane1+int(ifirst/2) && clust->GetZPlane()<=Plane2-int((ifirst+1)/2)) {
	if (clust->GetYPosErr()<100){
	  z=clust->GetZPos(); 
	  t=clust->GetYPos();
	  sw+=1.0; 
	  swx+=z; 
	  swx2+=z*z; 
	  swy+=t; 
	  swyx+=t*z;
	}
      }
      if(TMath::Abs(swx*swx-sw*swx2)>1e-10) {angle[1]= (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
      ifirst++;
    }
  }

  return angle;
}
*/

double InoTrack::GetXDir(int Plane1, int Plane2) {
  unsigned int nclusts = ClustsInTrack.size();
  
  double sum_iX = 0.0;
  double sum_iZ = 0.0;
  double sum_iXZ= 0.0;
  double sum_i2Z= 0.0;
  
  unsigned int NPlane = (unsigned int) abs(Plane2-Plane1);
  
  if (nclusts < 12) NPlane = 3;
  
  //if (nclusts < NPlane) NPlane = nclusts;
  
  for(unsigned int ij = 0; ij < NPlane; ij++)
    {
      InoCluster* TestClust = ClustsInTrack[ij];
      sum_iX	+= TestClust->GetXPos();
      sum_iZ	+= TestClust->GetZPos();
      sum_iXZ	+= TestClust->GetXPos()*TestClust->GetZPos();
      sum_i2Z += TestClust->GetZPos()*TestClust->GetZPos();
    }
  
  double slope_mx = (sum_iXZ - (sum_iX * sum_iZ)/NPlane)/(sum_i2Z-(sum_iZ*sum_iZ)/NPlane);
  
  /*
    double angle(0);
    double z=0.0,t=0.0;
    double sw=0.0,swx=0.0,swx2=0.0;
    double swy=0.0,swyx=0.0;
    
    for (unsigned int ij=0; ij<nclusts; ++ij) {
      InoCluster* clust = ClustsInTrack[ij];
      if (!clust->GetStraight()) continue;
      
      // find the clusts in the end region of the track
      // cout <<"Plane ij "<<ij<<" "<<clust->GetZPlane()<<" "<<clust->GetXPosErr()<<" "<<clust->GetXPos()<<" "<<clust->GetZPos()<<" "<<Plane1<<" "<<Plane2<<endl;
      if(clust->GetZPlane()>=Plane1 && clust->GetZPlane()<=Plane2) {
	if (clust->GetXPosErr()<100) {
	  z=clust->GetZPos();
	  t=clust->GetXPos();
	  sw+=1.0;
	  swx+=z;
	  swx2+=z*z;
	  swy+=t;
	  swyx+=t*z;
	  //	cout <<"sxy "<<sw<<" "<<swx<<" "<<swy<<endl;
	}
      }
      if(TMath::Abs(swx*swx-sw*swx2) >1e-10) {angle= (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
      //    cout <<"xx "<<swx*swx-sw*swx2<<" "<<Plane1<<" "<<Plane2<<" "<<fBegZPlane<<" "<<fEndZPlane<<" "<<angle<<endl;
      if (TMath::Abs(angle)>1e-10) return angle;
    }
    // cout <<"InoTrack::GetXDir============ "<<angle<<endl;
    */
  return slope_mx;
}

double InoTrack::GetYDir(int Plane1, int Plane2) {
  unsigned int nclusts = ClustsInTrack.size();
  
  double sum_iY = 0.0;
  double sum_iZ = 0.0;
  double sum_iYZ= 0.0;
  double sum_i2Z= 0.0;
  
  unsigned int NPlane = (unsigned int) abs(Plane2-Plane1);
  if (nclusts < NPlane) NPlane = nclusts;
  
  for(unsigned int ij = 0; ij < NPlane; ij++) {
    InoCluster* TestClust = ClustsInTrack[ij];
    sum_iY	+= TestClust->GetYPos();
    sum_iZ	+= TestClust->GetZPos();
    sum_iYZ	+= TestClust->GetYPos()*TestClust->GetZPos();
    sum_i2Z += TestClust->GetZPos()*TestClust->GetZPos();
  }
  
  double slope_my = (sum_iYZ - (sum_iY * sum_iZ)/NPlane)/(sum_i2Z-(sum_iZ*sum_iZ)/NPlane);
  
  /*
  double angle(0);
  double z=0.0,t=0.0;
  double sw=0.0,swx=0.0,swx2=0.0;
  double swy=0.0,swyx=0.0;
  
  for( unsigned int ij=0; ij<nclusts; ++ij) {
    InoCluster* clust = ClustsInTrack[ij];
    if (!clust->GetStraight()) continue;
    //find the clusts in the end region of the track
    //     cout <<"Planeyy ij "<<ij<<" "<<clust->GetZPlane()<<" "<<clust->GetYPosErr()<<" "<<clust->GetYPos()<<" "<<clust->GetZPos()<<" "<<Plane1<<" "<<Plane2<<endl;
    if(clust->GetZPlane()>=Plane1 && clust->GetZPlane()<=Plane2) {
      if (clust->GetYPosErr()<100) {
	z=clust->GetZPos();
	t=clust->GetYPos();
	sw+=1.0;
	swx+=z;
	swx2+=z*z;
	swy+=t;
	swyx+=t*z;
      }
    }
    if(TMath::Abs(swx*swx-sw*swx2)>1e-10) {angle = (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
    
    //    cout <<"yy "<<swx*swx-sw*swx2<<" "<<Plane1<<" "<<Plane2<<" "<<fBegZPlane<<" "<<fEndZPlane<<" "<<angle<<endl;
    if (TMath::Abs(angle)>1e-10) return angle;
  }
  */
  return slope_my;
}


/* 
double InoTrack::GetXDir(int Plane1, int Plane2) {
  
  unsigned int nclusts = ClustsInTrack.size();
  
  double angle(0);
  
  for( unsigned int ij=0; ij<nclusts; ++ij) {
    InoCluster* clust = ClustsInTrack[ij];
    if (!clust->GetStraight()) continue;
    //find the clusts in the end region of the track
    int ifirst = 0;
    while (angle==0 && ifirst <4) {
      double z=0.0,t=0.0;
      double sw=0.0,swx=0.0,swx2=0.0;
      double swy=0.0,swyx=0.0;
      //GMAGMA Why enter only once ??????===============
      //      cout <<"Plane ij "<<ij<<" "<<clust->GetZPlane()<<" "<<clust->GetXPosErr()<<" "<<clust->GetXPos()<<" "<<clust->GetZPos()<<" "<<Plane1<<" "<<int(ifirst/2)<<" "<<Plane2<<endl;
      if(clust->GetZPlane()>=Plane1+int(ifirst/2) && clust->GetZPlane()<=Plane2-int((ifirst+1)/2)) {
	if (clust->GetXPosErr()<100){
	  z=clust->GetZPos(); 
	  t=clust->GetXPos();
	  sw+=1.0; 
	  swx+=z; 
	  swx2+=z*z; 
	  swy+=t; 
	  swyx+=t*z;
//	  cout <<"sxy "<<sw<<" "<<swx<<" "<<swy<<endl;
	}
      }
      if(TMath::Abs(swx*swx-sw*swx2) >1e-10) {angle= (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
      //      cout <<"xx "<<swx*swx-sw*swx2<<" "<<angle<<" "<<ifirst<<endl;
      ifirst++;
      if (TMath::Abs(angle)>1e-10) return angle;
    }
  }
  
  //  cout <<"InoTrack::GetXDir============ "<<angle<<endl;
  
  return angle;
}

double InoTrack::GetYDir(int Plane1, int Plane2) {
  
  unsigned int nclusts = ClustsInTrack.size();
  
  double angle(0);

  for( unsigned int ij=0; ij<nclusts; ++ij) {
    InoCluster* clust = ClustsInTrack[ij];
    if (!clust->GetStraight()) continue;
    //find the clusts in the end region of the track
    int ifirst = 0;
    while (angle==0 && ifirst <4) {
      double z=0.0,t=0.0;
      double sw=0.0,swx=0.0,swx2=0.0;
      double swy=0.0,swyx=0.0;
      
      if(clust->GetZPlane()>=Plane1+int(ifirst/2) && clust->GetZPlane()<=Plane2-int((ifirst+1)/2)) {
	if (clust->GetYPosErr()<100){
	  z=clust->GetZPos(); 
	  t=clust->GetYPos();
	  sw+=1.0; 
	  swx+=z; 
	  swx2+=z*z; 
	  swy+=t; 
	  swyx+=t*z;
	}
      }
      if(TMath::Abs(swx*swx-sw*swx2)>1e-10) {angle = (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
      ifirst++;
    }
  }
  return angle;
}

*/


/*
double* InoTrack::GetBegDir() {
  return this->GetDir(fBegZPlane,fBegZPlane+10);
}
*/

double InoTrack::GetBegXDir() {
  return this->GetXDir(fBegZPlane,min(fBegZPlane+6,fEndZPlane)); //10->5
}

double InoTrack::GetBegYDir() {
  return this->GetYDir(fBegZPlane,min(fBegZPlane+6,fEndZPlane)); //10->5
}

/*
double* InoTrack::GetEndDir() {
  return this->GetDir(fEndZPlane-10,fEndZPlane);
}
*/

double InoTrack::GetEndXDir() {
  return this->GetXDir(max(fEndZPlane-5,fBegZPlane), fEndZPlane); //10->5
}

double InoTrack::GetEndYDir() {
  return this->GetYDir(max(fEndZPlane-5,fBegZPlane), fEndZPlane);  //10->5
}

void InoTrack::InsertCluster(vector<InoCluster*>::iterator it ,InoCluster* clust) {
  //  if(this->GetEntries()==0) { fPlaneView=clust->GetPlaneView(); }
  
  if(this->ContainsClust(clust)==true) {return;}
  ClustsInTrack.insert(it,clust);
  
  if( fBegZPlane > clust->GetZPlane()){ fBegZPlane = clust->GetZPlane(); fBegZ = clust->GetZPos();}
  if( fEndZPlane < clust->GetZPlane()){ fEndZPlane = clust->GetZPlane(); fEndZ = clust->GetZPos();}
  
  return;
}


void InoTrack::SetStraight(){
  int ifor=0;
 
  int ncls = GetEntries();

  for (int jk=0; jk<=int(ncls/2.); jk++) {
    if (GetCluster(jk+1)->GetZPlane() > GetCluster(jk)->GetZPlane()) {
      ifor++; GetCluster(jk)->SetStraight(true);
    } else {
      GetCluster(jk)->SetStraight(false); 
    }
  }
  for (int jk=int(ncls/2.); jk<ncls-1; jk++) {
    if (GetCluster(jk+1)->GetZPlane() > GetCluster(jk)->GetZPlane()) {
      ifor++; GetCluster(jk+1)->SetStraight(true);
    } else {
      GetCluster(jk+1)->SetStraight(false); 
    }
  }

  /*
  int ncls = GetEntries()-1;
  for (int jk=0; jk<ncls; jk++) {
    if (GetCluster(jk+1)->GetZPlane() > GetCluster(jk)->GetZPlane()) {
      ifor++; GetCluster(jk+1)->SetStraight(true);
      if (jk==0) GetCluster(jk)->SetStraight(true);
    } else {
      GetCluster(jk+1)->SetStraight(false); 
      if (jk==0) GetCluster(jk)->SetStraight(false);
    }
  }

   
  if (ifor<ncls/2.) {
    for (int jk=0; jk<ncls; jk++) {
      if (GetCluster(jk+1)->GetZPlane() > GetCluster(jk)->GetZPlane()) {
	GetCluster(jk+1)->SetStraight(false); 
	if (jk==0) GetCluster(jk)->SetStraight(false);
      } else {
	GetCluster(jk+1)->SetStraight(true);
	if (jk==0) GetCluster(jk)->SetStrainght(true);
      }
    }
  }
  */
}

