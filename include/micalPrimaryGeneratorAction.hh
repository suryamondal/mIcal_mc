#ifndef micalPrimaryGeneratorAction_h
#define micalPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "MultiSimAnalysis.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
//#include <iostream.h>
//#include <fstream.h>
#include <iostream>
#include <fstream>
#include "micalDetectorParameterDef.hh"
#include "TFile.h"
class G4ParticleGun;
class G4Event;
class micalDetectorConstruction;
class micalPrimaryGeneratorMessenger;

class micalPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  micalPrimaryGeneratorAction(micalDetectorConstruction*, MultiSimAnalysis *panalysis);    
  ~micalPrimaryGeneratorAction();
  void OpenFileCORSIKA();
  void CloseFileCORSIKA();
  void OpenFileFLUX();
  void CloseFileFLUX();
  static micalPrimaryGeneratorAction* AnPointer;
  void GeneratePrimaries(G4Event*);
  void SetRunNumber(G4int p) { RunNumber = p;}
  void SetInputFlag(G4int p);// { InputFlag = p;}
  void SetFirstEvt(G4int p) { FirstEvt = p;}
  void SetRndmFlag(G4String p) { rndmFlag = p;}
  void SetPartId(G4int p) { partId = p;}
  void SetIncEnergySmr(G4double p) {incEnergySmr = p;}
  void SetIncEnergy(G4double p) {incEnergy = p;}
  void SetIncDirection(G4ThreeVector p) {incDirection = p;}  
  void SetIncThetaSmr(G4double p) {incThetaSmr =p;}
  void SetIncPhiSmr(G4double p) {incPhiSmr =p;}
  void SetIncPosition(G4ThreeVector p);// {incPosition = p;} 
  void SetIncVxSmr(G4double p) {incVxSmr =p;}
  void SetIncVySmr(G4double p) {incVySmr =p;}
  void SetIncVzSmr(G4double p) {incVzSmr =p;}
  void SetCorsikaFileName(G4String p) {CorsikaFileName=p;}
  void SetCorsikaFileDir(G4String p) {CorsikaFileDir=p;}
  void SetFluxFileName(G4String p) {FluxFileName=p;}
  //For inputflag 1
  void SetPowerCosTheta(G4double p) {PowCosTheta=p;}
  void SetDetectorThetaCover(G4double p) {DetThetaCov=p*pivalGA/180.;}
  void SetPowerCosmicEnergy(G4double p) {PowCosmicEnr=p;}
  void SetELowLimit(G4double p) {ELowLim=p*GeV;}
  void SetEUpLimit(G4double p) {EUpLim=p*GeV;}
  void SetTopTrgLay(G4int p) {toptrgly=p;}
  void SetBottomTrgLay(G4int p) {bottomtrgly=p;}
  double GetCosmicEnergy(double ELimLow, double ELimUp);
  int LinePlaneInt(double* Line, double* Plane, double* Point);
  double energy_func(double x1);
  
  TFile *FileFLUX;
  TFile *FileCORSIKA;
  TTree *TreeCORSIKA;
  TH3F *muFlux, *mupFlux, *munFlux;
  TH3D *corsikaFlux;
  TH3D *corsikaFluxMuM;
  TH3D *corsikaFluxMuP;
  G4int InputFlag; 
  G4int FirstEvt;  
  G4ParticleGun* particleGun;
  micalDetectorConstruction* micalDetector;
  micalDetectorParameterDef* paradef;    
  MultiSimAnalysis *pAnalysis;    
  micalPrimaryGeneratorMessenger* gunMessenger;
  G4int RunNumber;
  G4String rndmFlag;
  G4int partId;
  G4double incEnergy;
  G4double incEnergySmr;
  G4ThreeVector incDirection;
  double MuPProb;
  double MuMProb;
  double TotalProb;
  G4double PowCosTheta;
  G4double DetThetaCov;
  G4double PowCosmicEnr;
  G4double norm1;
  G4double norm2;
  G4double ENorm;
  G4double ELowLim;
  G4double EUpLim;
  double pivalGA; // = 3.14159265;
  double pargas[3];
  double RPCLayerPosZ[12];
  double IRONLayerPosZ[13];
  double LayerZdim[12];
  double IronLayerZdim[13];
  double StackPosInWorld[3];
  double WorldXDim;
  double WorldYDim;
  double WorldZDim;
  int initialise;
  int initialiseCor;
  int g_nevt;
  G4ThreeVector incPosition;
  G4double incVxSmr;
  G4double incVySmr;
  G4double incVzSmr;
  G4double incThetaSmr;
  G4double incPhiSmr;

  int bottomtrgly;
  int toptrgly;
  
  G4double fHalfLayerThickness;
  G4int nLayer;
  G4int nIRLayer;
  
  int nINODet;
  float parino[3];
  float gapino;

  double enerin[20];
  
  static const int iC_npartmax = 100;
  int           iC_nevt;
  int           iC_npart;
  int           iC_cpid[iC_npartmax];   //[iC_npart]
  float         iC_cvx[iC_npartmax];   //[iC_npart]
  float         iC_cvy[iC_npartmax];   //[iC_npart]
  float         iC_cpx[iC_npartmax];   //[iC_npart]
  float         iC_cpy[iC_npartmax];   //[iC_npart]
  float         iC_cpz[iC_npartmax];   //[iC_npart]
  float         iC_eventweight;
  
  G4String CorsikaFileDir;
  G4String CorsikaFileName;
  G4String FluxFileName;
};

#endif
