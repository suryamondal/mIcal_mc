#include <cassert>
#include <iostream>
#include "InoStrip.h"

//______________________________________________________________________

//______________________________________________________________________
InoStrip::InoStrip() {
  fView = -1; 
  fStrip =-1000;
  fPlane =-1000; 
  fXYPos =-1000; 
  fZPos = -1000; 
  iTrueTime = -1000; 
  iSmrTime = -1000; 
  fPulse = 1.; //Put it ONE, such that no effect on digital signal

  fMomentum=-1000.;
  fTheta = -1000.;
  fPhi =-1000.;

  fXgen = -1000.;
  fYgen = -1000.;
  fZgen = -1000.;

  fId    = -1;  
  fRPCmod =-1;
}

//______________________________________________________________________
InoStrip::InoStrip(InoStrip* cd) {
  fView = cd->fView;  
  fStrip = cd->fStrip;
  fPlane = cd->fPlane; 
  fXYPos = cd->fXYPos; 
  fZPos = cd->fZPos;  
  iTrueTime =  cd->iTrueTime; 
  iSmrTime =  cd->iSmrTime; 
  fId   =  cd->fId;
  fRPCmod =  cd->fRPCmod;
  iSmrTime =cd->iSmrTime;
  iTrueTime = cd->iTrueTime;
  

}

//______________________________________________________________________
InoStrip::~InoStrip() {
}

//______________________________________________________________________
InoStrip *InoStrip::DupHandle() const {
   return (new InoStrip(*this));
}


//______________________________________________________________________
void InoStrip::Trace(const char *c) const {
  std::cout<<"**********Begin InoStrip::Trace(\"" << c << "\")" << std::endl
	   << "Information from InoStrip's CandHandle: " << std::endl;
  //  CandHandle::Trace(c);
  std::cout
    << "**********End InoStrip::Trace(\"" << c << "\")" << std::endl;
}

/*
int GetRPCmodX() {
 
  unsigned int idet = fId;
  idet>>=15; //8 are dummy and 7 for strip
  return idet%65536; // 2**16
}

void GetXYPosXYZ(double xypos, double zpos) {
  unsigned int idet = fId;
  idet>>=8;
  int ninxy = idet%128;
  idet>>=7;
  int ninch = idet%8;
  idet>>=3;
  int ninmo = idet%8;  
  idet>>=3;
  int ninla = idet%256;
  idet>>=8;
  int nindt = idet%4;
  idet>>=2;
  
  if (idet==0) {
    xypos = (1/m)*( (nindt-1)*(2*parino[0]+gapino) - parlay[0]  + (2*ninmo+1)*parmod[0] -pargas[0] + xstrwd*(ninx+0.5)); //gma use global variables (for all these three co-ordinates)
    zpos = (1/m)*(-parino[2] + 2*parirnlay[2] + (2*ninla+1)*parlay[2]);
  } else {
    xypos = (1/m)*(- parmod[1]  + (2*ninch+1)*parchm[1] -pargas[1] + ystrwd*(niny+0.5));
    zpos = (1/m)*(-parino[2] + 2*parirnlay[2] + (2*ninla+1)*parlay[2]);
  }
}


int GetPlaneViewX() {
 
  unsigned int idet = fId;
  idet>>=31;
  return (int)idet;
}
*/
