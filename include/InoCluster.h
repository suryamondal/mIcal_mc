#ifndef INOCLUSTER_H
#define INOCLUSTER_H
// ClusterCam

#include <vector>
#include "micalDetectorParameterDef.hh"
using std::vector;

class InoHit;
class InoCluster
{
 public:
  InoCluster(InoHit* hit);
  virtual ~InoCluster();
 
  void Print();
  void AddHit(InoHit* hit);
  bool ContainsHit(InoHit* hit);

  int IsHitAssoc(InoHit* hit) const;
  int IsShwAssoc(InoCluster* clust) const;
  int IsTrkAssoc(InoCluster* clustm, InoCluster* clustp) const;
  int IsDiffuseShwAssoc(InoCluster* clr) const;

  unsigned int GetHitEntries() const {return HitsInCluster.size();}
  unsigned int GetXEntries();
  unsigned int GetYEntries();
  unsigned int GetXProjEntries();
  unsigned int GetYProjEntries();

  InoHit* GetHit(unsigned int i) const;
  int GetDigits() const { return fDigits; }

  int GetZPlane() const {return fZPlane;}
  int GetRPCmod() const {return fRPCmod;}
  int GetBegXStrip() const {return fBegXStrip;}
  int GetEndXStrip() const {return fEndXStrip;}
  int GetBegYStrip() const {return fBegYStrip;}
  int GetEndYStrip() const {return fEndYStrip;}

  int GetView() const {return fView;}

  double GetZPos() {return fZPos;}
  double GetXPos() {return fXPos;}
  double GetYPos() {return fYPos;}
  

  double GetPulse() const { 
    if (fXPulse >0 && fYPulse >0) {
      return 0.5*(fXPulse+fYPulse);
    } else if (fXPulse >0) {
      return fXPulse;
    } else {
      return fYPulse;
    }
  };  

  double GetXPulse() {return fXPulse;}
  double GetYPulse() {return fYPulse;}

  double GetTime() const {return 0.5*(fBegTime+fEndTime);}
  double GetBegTime() const {return fBegTime;}
  double GetEndTime() const {return fEndTime;}

  double GetBegXPos() const {return fBegXPos;}
  double GetEndXPos() const {return fEndXPos;}
  double GetBegYPos() const {return fBegYPos;}
  double GetEndYPos() const {return fEndYPos;}

  double GetXPosErr() const    {return fXPosErr;}
  double GetYPosErr() const    {return fYPosErr;}

  void SetXPosErr(double err) {fXPosErr = err;}
  void SetYPosErr(double err) {fYPosErr = err;}

  int GetTrkFlag() const {return fTrkFlag;}
  int GetShwFlag() const {return fShwFlag;}
  int GetTrkPlnFlag() const {return fTrkPlnFlag;}
  int GetShwPlnFlag() const {return fShwPlnFlag;}

  void SetTrkFlag(int flag) {fTrkFlag=flag;}
  void SetShwFlag(int flag) {fShwFlag=flag;}
  void SetTrkPlnFlag(int flag) {fTrkPlnFlag=flag;}
  void SetShwPlnFlag(int flag) {fShwPlnFlag=flag;}
  
  void SetNDFlag(int flag) {fNDFlag=flag;}
  int GetNDFlag() const {return fNDFlag;}
  void SetInTrack (bool tag) { InTrack=tag;}
  bool GetInTrack () { return InTrack;}
  void SetInShower (bool tag) { InShower=tag;}
  bool GetInShower () { return InShower;}

  void SetStraight (bool tag) { isStraight=tag;}
  bool GetStraight () { return isStraight;}

  int fZPlane;
  int fRPCmod;
  int fBegXStrip;
  int fEndXStrip;
  int fBegYStrip;
  int fEndYStrip;

  double fBegTime;
  double fEndTime;
  double fBegXPos;
  double fEndXPos;
  double fBegYPos;
  double fEndYPos;
  double fXPos;
  double fYPos;
  double fZPos;
  double fXPulse;
  double fYPulse;
  int fTrkFlag;
  int fShwFlag;
  int fTrkPlnFlag;
  int fShwPlnFlag;

  int fDigits;
  int fNDFlag;

  double fXPosErr; //Error in X-position (m)
  double fYPosErr; //Error in Y-position (m
  int    fView;    // 
                   // 0 : only X(U)-axis hit
                   // 1 : only Y(V)-axis hit
                   // 2 : Both X and Y-axis hit 

  //GMA put in public place to get hits size in InoTrackFinder.cc file for informations of visualisation2
  vector<InoHit*> HitsInCluster;
  bool isIdentical(InoCluster* icls);
 private:
  micalDetectorParameterDef *paradef;
  bool   InTrack;
  bool   InShower;
  bool   isStraight; //Is it in straight section of track
  double StripXWidth;
  double StripYWidth;

};

#endif
