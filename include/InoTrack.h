#ifndef INOTRACK_H
#define INOTRACK_H
//TrackCam

#include "TObject.h"
#include <vector>
#include<cstdlib>
#include <cmath>
using namespace std;
using std::vector;
class inoTrackCand;
class InoTrackSegment;
//class InoHit;
class InoCluster;

class InoTrack //: public TObject
{
  friend class InotTrackCand;
 public:
  InoTrack();
  InoTrack(InoTrackSegment* segment);
  ~InoTrack();
  
  InoTrackSegment* GetInoTrackSegment() const {return fSegment;};

  int GetBegZPlane() const {return fBegZPlane;};
  int GetEndZPlane() const {return fEndZPlane;};

  //  void AddHit(InoHit* hit);
  //  bool ContainsHit(InoHit* hit) const;
  //  InoHit* GetHit(unsigned int i) const;

  void AddCluster(InoCluster* clust);
  void AddTrack(InoTrack* trk);  
  
  bool ContainsClust(InoCluster* clust) const;
  InoCluster* GetCluster(unsigned int i) const;
  void InsertCluster(vector<InoCluster*>::iterator it,InoCluster* cls);
  vector<InoCluster*>::iterator begin(){return ClustsInTrack.begin();};
  vector<InoCluster*>::iterator end(){return ClustsInTrack.end();};

  double GetBegXPos();
  double GetEndXPos();

  double GetBegYPos();
  double GetEndYPos();

  double GetXDir(int Plane1, int Plane2);
  double GetBegXDir();
  double GetEndXDir();

  double GetYDir(int Plane1, int Plane2);
  double GetBegYDir();
  double GetEndYDir();

  double GetBegZPos() const {return fBegZ;};
  double GetEndZPos() const {return fEndZ;};

  unsigned int GetEntries() const {return ClustsInTrack.size();} // HitsInTrack.size();};

  int GetUID() const {return fUID;};
  void SetUID(int UIDNum) {fUID=UIDNum;};

  int GetUsed() const {return fUsed;};
  void SetUsed(int UIDNum) {fUsed=UIDNum;};

  void SetStraight(); //11Nov09 sortout clusters according toslope in Z

  //  vector<InoHit*>HitsInTrack;
  vector<InoCluster*>ClustsInTrack;

 private: 
  InoTrackSegment* fSegment;
  int fBegZPlane;
  int fEndZPlane;
  double fBegZ;
  double fEndZ;
  int fUID;
  int fUsed;

};

#endif
