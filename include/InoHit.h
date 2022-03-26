#ifndef INOHITATNU_H
#define INOHITATNU_H
//HitCam
//GMA use proper fView value to calcualte X/Y-direction(length) peoperly

//#include "MessageService/MsgService.h"
#include "InoStrip.h"
#include "micalDetectorParameterDef.hh"
#include "MultiSimAnalysis.hh"

using namespace std;

class InoHit
{

public:
  InoHit();
  InoHit(InoStrip* fx, InoStrip* fy);
  InoHit(InoStrip* fx);
  InoHit(InoHit* hit);
  ~InoHit();

  InoStrip* GetXStrip() const {return fXStrip;};
  InoStrip* GetYStrip() const {return fYStrip;};
  void Print();

  double GetPulse() const { 
    if (fView==2) {
      return 0.5*(fXPulse+fYPulse);
    } else if (fView==1) {
      return fYPulse;
    } else {
      return fXPulse;
    }
  };
  
  double GetXPulse() const {return fXPulse;};  
  void SetXPulse(double q) {fXPulse=q;};

  double GetYPulse() const {return fYPulse;};  
  void SetYPulse(double q) {fYPulse=q;};

  int GetXPlane() const     {return fXPlane;};
  void SetXPlane(int Plane) {fXPlane = Plane;};

  int GetYPlane() const     {return fYPlane;};
  void SetYPlane(int Plane) {fYPlane = Plane;};

  int GetXStripNum() const     {return fXStripNum;};
  void SetXStripNum(int strip) {fXStripNum = strip;};

  int GetYStripNum() const     {return fYStripNum;};
  void SetYStripNum(int strip) {fYStripNum = strip;};

  int GetRPCmod() const;

double GetTime() const; // Returns Avg of X and Y time with the delay correction (in ns).
double GetXTime() const {return DigiToTimeConv*fXTime;}; // Returns Smeared Time of X strip (in ns).
double GetYTime() const {return DigiToTimeConv*fYTime;}; // Returns Smeared Time of Y strip (in ns).
double GetXTimeCorr() const; // Returns Smeared Time of X strip with delay correctino (in ns).
double GetYTimeCorr() const; // Returns Smeared Time of Y strip with delay correctino (in ns).
double GetXTrueTime() const {return DigiToTimeConv*fXTrueTime;}; // Returns True Time of X strip (in ns).
double GetYTrueTime() const {return DigiToTimeConv*fYTrueTime;}; // Returns True Time of X strip (in ns).
int GetXpdgId() const {return fXpdgId;};
int GetYpdgId() const {return fYpdgId;};

double GetXPos() const    {return fXPos;};
void SetXPos(double tpos) {fXPos = tpos;};
double GetXPosErr() const    {return fXPosErr;};

double GetYPos() const    {return fYPos;};
void SetYPos(double tpos) {fYPos = tpos;};
double GetYPosErr() const    {return fYPosErr;};

int GetTrkFlag() const    {return fTrackFlag;};
void SetTrkFlag(int flag) {fTrackFlag=flag;};

int GetShwFlag() const    { return fShowerFlag; };
void SetShwFlag(int flag) { fShowerFlag=flag; };

  int GetUID() const   {return fUid;};
  void SetUID(int uid) {fUid=uid;};

  double GetZPos() const    {return fZPos;};
  void SetZPos(double zpos) {fZPos = zpos;};

  int GetZPlane() const    {return fZPlane;};
  void SetZPlane(int zpl) {fZPlane = zpl;};

  int GetView() const    {return fView;};
  void SetView(int zpl) {fView = zpl;};

  int IsDiffuseShwAssoc(InoHit* hit) const; //new
  int IsShwAssoc(InoHit* hit) const; //new

  void SetMomentum(double f) {fMomentum = f;};
  void SetTheta(double f) {fTheta = f;};
  void SetPhi(double f) {fPhi = f;};

  double GetMomentum() {return fMomentum ;};
  double GetTheta() {return fTheta ;};
  double GetPhi() {return fPhi;};

  bool isIdentical(InoHit* hit);

private:
  micalDetectorParameterDef* paradef;
  MultiSimAnalysis* pAnalysis;
  double StripXWidth;
  double StripYWidth;
  InoStrip* fXStrip;   //Strip for X-axis 
  InoStrip* fYStrip;   //Strip for Y-axis
  int fUid;            //User ID for track/cluster
  //  int fPlaneView;

  //Hit parameters may not be exactly the strips paramters, thus 
  //for the time being copy those, lateron modified those as equired

  int fXStripNum;      //StripID for X-axis strips 
  double fXTime;      // Timing of X-strip (ns)
  double fXTrueTime;      // True Timing of X-strip (ns)
  int fYStripNum;     //StripID for Y-axis strips  
  double fYTime;     // Timing of Y-strip (ns)
  double fYTrueTime;      // True Timing of X-strip (ns)

  int fXpdgId; // Particle ID which produces hit in X-strip
  int fYpdgId; // Particle ID which produces hit in Y-strip

  int fTrackFlag;   //Tag for track association
  int fShowerFlag;    //Tag for shower association
  double fXPulse;  //Not useful for RPC, but let it be
  int fXPlane;     // =0 indicating strip measured X-axis
  double fYPulse; // Not useful for RPC, but let it be
  int fYPlane;     // =1 indicating strip measured X-axis

  double fXPos;    //X-position (m)
  double fXPosErr; //Error in X-position (m)
  double fYPos;    //Y-position (m)
  double fYPosErr; //Error in Y-position (m)
  int    fZPlane;  //Z-plane
  double fZPos;    //Z-postion (m)
  int    fView;    // 
                   // 0 : only X(U)-axis hit
                   // 1 : only Y(V)-axis hit
                   // 2 : Both X and Y-axis hit 

  double fMomentum; //Momentum of track which is behind this hit
  double fTheta;  //Theta
  double fPhi;    // Phi

double DigiToTimeConv;
double SignalSpeed;

};

#endif
