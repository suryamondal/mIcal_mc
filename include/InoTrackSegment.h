#ifndef INOTRACKSEGMENT_H
#define INOTRACKSEGMENT_H
//TrackSegmentCam
//GMA need to put different Z-plane corresponding to X and Y-strip separately

#include <vector>
using std::vector;

class InoCluster;

class InoTrackSegment
{
 public:
  InoTrackSegment(InoCluster* clustm, InoCluster* clust0, InoCluster* clustp);
  virtual ~InoTrackSegment();
  
  void AddCluster(InoCluster* clust);
  bool ContainsCluster(InoCluster* clust);
  InoCluster* GetCluster(unsigned int i);
  unsigned int GetEntries() const;
  int GetBegZPlane() const;
  int GetEndZPlane() const;

  //GMA put it properly in .cc file 
  int GetBegZXPlane() const {return  GetBegZPlane();};
  int GetEndZXPlane() const {return  GetEndZPlane();};
  int GetBegZYPlane() const {return  GetBegZPlane();};
  int GetEndZYPlane() const {return  GetEndZPlane();};

  double GetBegXDir();
  double GetBegYDir();
  double GetBegTPos();
  double GetBegXPos();
  double GetBegYPos();
  double GetBegZPos() const;

  double GetEndXDir();
  double GetEndYDir();
  double GetEndXPos();
  double GetEndYPos();
  double GetEndZPos() const;
  void AddSegment(InoTrackSegment* segment);
  bool IsAssoc(InoTrackSegment* segment);


  //Associated Segments (Lowest level of association)
  void AddAssocSegToBeg(InoTrackSegment* seg);   
  void AddAssocSegToEnd(InoTrackSegment* seg);   
  InoTrackSegment* GetAssocSegBeg(unsigned int i);
  InoTrackSegment* GetAssocSegEnd(unsigned int i);
  unsigned int GetNAssocSegBeg() const{return fBegAssociatedSegList.size();};  
  unsigned int GetNAssocSegEnd() const{return fEndAssociatedSegList.size();};  

  //Preferred Segments (Second level of association)
  void AddPrefSegToBeg(InoTrackSegment* seg);          
  void AddPrefSegToEnd(InoTrackSegment* seg);          
  InoTrackSegment* GetPrefSegBeg(unsigned int i);      
  InoTrackSegment* GetPrefSegEnd(unsigned int i);      
  unsigned int GetNPrefSegBeg() const{return fBegPreferredSegList.size();};
  unsigned int GetNPrefSegEnd() const{return fEndPreferredSegList.size();};

  //Matched Segments (Third level of association)
  void AddMatchSegToBeg(InoTrackSegment* seg);
  void AddMatchSegToEnd(InoTrackSegment* seg);
  InoTrackSegment* GetMatchSegBeg(unsigned int i); 
  InoTrackSegment* GetMatchSegEnd(unsigned int i); 
  unsigned int GetNMatchSegBeg() const{return fBegMatchedSegList.size();};
  unsigned int GetNMatchSegEnd() const{return fEndMatchedSegList.size();};


  void SetTmpTrkFlag(int flag){ fTmpTrkFlag=flag; };
  int GetTmpTrkFlag() const { return fTmpTrkFlag; };

  void SetUID(int uid){ fUID=uid; };
  int GetUID() const { return fUID; };

  void SetTrkFlag(int flag){ fTrkFlag=flag; };
  int GetTrkFlag() const { return fTrkFlag; };

  void SetSeedSegment(InoTrackSegment* segment) {fSeedSegment = segment;};
  InoTrackSegment* GetSeedSegment() {return (InoTrackSegment*)(fSeedSegment);};

  void SetNPlanes(int nplanes){ fNPlanes=nplanes; };
  int GetNPlanes() const { return fNPlanes; };

  double GetScore(vector<InoTrackSegment*> *BegSegBank=0, vector<InoTrackSegment*> *EndSegBank=0);

  double GetBegTime() const { return fBegTime; };
  double GetEndTime() const { return fEndTime; };

  void SetPartner(InoTrackSegment* segment) {fPartner = segment;};
  InoTrackSegment* GetPartner() {return (InoTrackSegment*)(fPartner);};


  vector<InoCluster*> ClustersInSegment;
 private:
  vector<InoTrackSegment*> fBegAssociatedSegList;
  vector<InoTrackSegment*> fEndAssociatedSegList;
  vector<InoTrackSegment*> fBegPreferredSegList;
  vector<InoTrackSegment*> fEndPreferredSegList;
  vector<InoTrackSegment*> fBegMatchedSegList;
  vector<InoTrackSegment*> fEndMatchedSegList;

  InoTrackSegment* fSeedSegment; // Last segment in chain of segments to which this segment belongs

  InoTrackSegment* fPartner;

  int fUID;
  int fBegZPlane;
  int fEndZPlane;
  double fBegVtxZ;
  double fEndVtxZ;
  int fTrkFlag;
  int fTmpTrkFlag;
  //  int fPlaneView;
  int fNPlanes;
  double fBegTime;
  double fEndTime;
  const double StripWidth;

};

#endif
