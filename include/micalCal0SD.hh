#ifndef micalcal0SD_h
#define micalcal0SD_h 1
#include "G4VSensitiveDetector.hh"
#include "micalCal0Hit.hh"
#include "MultiSimAnalysis.hh"
#include "micalCal0SDMessenger.hh"
#include "micalDetectorParameterDef.hh"
#include "InoHit.h"
#include <vector>
#include "vect_manager.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class micalcal0SDMessenger;



class micalcal0SD : public G4VSensitiveDetector
{
  
public:
  micalcal0SD(G4String name);
  ~micalcal0SD();
  
  void Initialize(G4HCofThisEvent*HCE);
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*HCE);
  void clear();
  void DrawAll();
  void PrintAll();
  void SetCorrTimeSmear(G4double val);
  void SetUnCorrTimeSmear(G4double val);
  void SetCorrInefficiency(G4double val);
  void SetUnCorrXInefficiency(G4double val);
  void SetUnCorrYInefficiency(G4double val);
  void SetTimeToDigiConv(G4double val);
  void SetSignalSpeed(G4double val);
  void SetCorrNoise1(G4double val);
  void SetCorrNoise2(G4double val);
  void SetCorrNoise3(G4double val);
  void SetRandomNoise(G4int val);
  void SetRootRandom(G4int val);
  int GetRandomXY(double& GapX, TH2D* tmphistx);
  
private:
  micalcal0HitsCollection *cal0Collection;
  micalcal0SDMessenger* cal0SDMessenger;
  MultiSimAnalysis *pAnalysis;
  micalDetectorParameterDef* paradef;
  /*
    const G4int numberInMO;
    const G4int numberInCH;
    const G4int numberInLA;
  */
  G4int numberInINO;
  G4int numberInMO;
  G4int numberInCH;
  G4int numberInLA;
  
  // double ShiftInX;
  // double ShiftInY;
  // double ShiftInZ;
  int RootRandom;
  //const G4int numberInX;
  //const G4int numberInY;
  G4int numberInX;  //AAR if stripwidth is variable then total number of strips will also change.
  G4int numberInY;  //AAR
  
  const G4int numberInT;
  const G4int numberInCell;
  G4int InCell;

  int NewMultiplicity;
  double TimeCorrSmr;
  double TimeUnCorrSmr;
  double CorrIneffiPar;
  double UnCorrXIneffiPar;
  double UnCorrYIneffiPar;
  double TimeToDigiConv;
  double SignalSpeed;
  double CorrNoisePar1;
  double CorrNoisePar2;
  double CorrNoisePar3;
  int RandomNoisePar;
  float ResolutionOfRPCAlign[3][65536];
  
  float gapino;
  float parino[3]; //={1606.01*cm, 706.01*cm, 598.0*cm};
  float parlay[3]; //={1600.0*cm, 700.02*cm, 4.26*cm};
  float parmod[3]; //={100.01*cm, 700.01*cm, 4.25*cm};
  float parchm[3]; //={100.0*cm,  100.0*cm, 4.24*cm};
  //  float parair[3]; //={98.005*cm,   98.005*cm, 1.25*cm};
  //  float parirnlay[3]; //={1606.0*cm,  706.0*cm, 1.5*cm};
  float parcup[3];//={98.004*cm,  98.004*cm, 0.42*cm};
  // float parg10[3];//={98.003*cm,  98.003*cm, 0.40*cm};
  float parqurz[3];//={98.002*cm,  98.002*cm, 0.30*cm};
  float pargas[3];//={98.001*cm,  98.001*cm, 0.10*cm};

  float parirlay[3]; //thickness of iron layer
  float parhcoil[3];  //space for coild
  float parcoilsupport[3]; // Coild support

  float Xstrwd;// = 2.0*cm;
  float Ystrwd;// = 2.0*cm;
  int   nINODet;
  unsigned long long int DeadStrpX[3][150][8][8];
  unsigned long long int DeadStrpY[3][150][8][8];
  unsigned int CellDetID[20000];
  //  unsigned long CellDetID[20000];
  vector<InoHit> InoHit_list;
  
  //     InoHit_Manager* inoHit_pointer;
  InoStripX_Manager* inoStripX_pointer;
  InoStripY_Manager* inoStripY_pointer;
  public:
  float histxmn, histxmx, histymn, histymx, histzmn, histzmx;
  unsigned int twopow31;
};




#endif

