#ifndef INOTRACKFINDER_H
#define INOTRACKFINDER_H
//AlgTrackCamList

#include <vector>
#include "vect_manager.h"
#include "InoTrackSegment.h"
#include "micalDetectorParameterDef.hh"
#include "MultiSimAnalysis.hh"
// #include "G4SystemOfUnits.hh"
using std::vector;

class InoTrackFinder // : public AlgBase
{
  
 public:
  InoTrackFinder();
  ~InoTrackFinder();
  void RunAlg() ;
  void Trace(const char *c) const;

  void RunTheFinder(); //CandSliceHandle* slice);
  void FormTheHits(); //CandSliceHandle* slice);
  void FormTheClusters();
  void IDTrkAndShwClusters();
  void FormTriplets();
  void FindAllAssociations();
  void FindPreferredJoins();
  void FindMatchedJoins();

  void FirstComparison();
  void FormTracks();
  void JoinTracks();
  void FormFinalTracks(); //CandSliceHandle* slice);
  void LookForHitsAcrossGap(InoTrack* Trk);
  void JoinCurvedTrack();
  void ExtendTrack(InoTrack* Trk);
  void FillGapsInTrack(InoTrack* Trk);
  void CleanAndFilled();
  void ClearUp();

  InoCluster_Manager* inoCluster_pointer; 
  MultiSimAnalysis *pAnalysis;
  micalDetectorParameterDef* paradef;
  
 private:
  //  vector<InoHit*> AllHitBank[500]; //GMA should be given though database of static constant

  vector<InoHit*> HitBank[500];  //GMA Same

  //  vector<InoHit*> HitsInTracks[2];

  vector<InoCluster*> ClusterBank[500];
  //  vector<InoCluster*> AllClusterBank[500];

  vector<InoCluster*> ClusterList[500];

  vector<InoCluster*> ClustsInTracks[2];

  vector<InoTrackSegment*> SegmentBank[500];
  vector<InoTrackSegment*> NDSegmentBank[500];
  
  vector<InoTrackSegment*> ViewSegBank[2];
  vector<InoTrackSegment*> TempTrack[2];

  vector<InoTrackSegment*> PossibleJoins[2];

  vector<InoTrack*> FinalTrackBank[2];

  vector<InoCluster*> ClustList[500];
  vector<InoTrack*> FinalTrackTempHolder[2];

  int NumModules;
  int ModuleType;  //2=near detector 1=fardetector
  //  const double StripWidth;

  double PECut, PECut2;
  InoHit_Manager* inoHit_pointer; 
  InoRPCStrip_Manager* inoRPC_pointer;
  //  InoCluster_Manager* inoCluster_pointer; 
  InoTrack_Manager* inoTrack_pointer; 
  unsigned long long int NoisyStrpX[3][150][8][8];
  unsigned long long int NoisyStrpY[3][150][8][8];
  double StripXWidth;
  double StripYWidth;
  int nXStrips;
  int nYStrips;
  double GasChmX;
  double GasChmY;
  double LayerThickness;
  int PlanesInModule;

  double DigiToTimeConv;
  double SignalSpeed;
  
    
};

#endif // INOTRACKFINDER_H
