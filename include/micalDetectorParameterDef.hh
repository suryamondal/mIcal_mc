#ifndef micalDetectorParameterDef_h
#define micalDetectorParameterDef_h 1
//#include "G4SIunits.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include "G4ios.hh"
#include <G4String.hh>
// #include "TFile.h"
// #include "TNtuple.h"
// #include <libpq-fe.h>
#include <strings.h>

class micalDetectorParameterDef {
  public:
  micalDetectorParameterDef();
  ~micalDetectorParameterDef(){};
  static micalDetectorParameterDef* AnPointer;
  void UpdateDetectorParameterDef();

  double GetParworld(int i) { return parworld[i];}
  double GetParroom(int i) { return parroom[i];}
  double GetParairroom(int i) { return parairroom[i];} 
  double GetParairroom2(int i) { return parairroom2[i];} 
  double GetParstaircaseair(int i) { return parstaircaseair[i];} 
  double GetParstaircasel(int i) { return parstaircasel[i];} 
  double GetParstaircase(int i) { return parstaircase[i];} 
  double GetParino(int i) { return parino[i];}
  double GetParmagnet(int i) { return parmagnet[i];}
  double GetShiftInX()       { return ShiftInX;} 
  double GetShiftInY()       { return ShiftInY;} 
  double GetShiftInZ(int i)       { return ShiftInZ[i];} 
  double GetAlShiftInAirBox() {return AlShiftInAirBox;}
  double GetAirBoxShiftInFRP() {return AirBoxShiftInFRP;}
  double GetFRPshiftInLayer(int i) {return FRPshiftInLayer[i];}

  void SetINOMPositionGlobalX(double xpos1){INOMPositionGlobalX = xpos1;};
  void SetINOMPositionGlobalY(double xpos1){INOMPositionGlobalY = xpos1;};
  void SetINOMPositionGlobalZ(double xpos1){INOMPositionGlobalZ = xpos1;};
  void SetLayerPosInStack(int ilay, double zpos1) {LayerPosInStack[ilay] = zpos1;}
  double GetLayerPosInStack(int ilay) {return LayerPosInStack[ilay];}
  double GetINOMPositionGlobalX() {return INOMPositionGlobalX;}
  double GetINOMPositionGlobalY() {return INOMPositionGlobalY;}
  double GetINOMPositionGlobalZ() {return INOMPositionGlobalZ;}
  double GetParcoilspacerpc(int i) { return parcoilspacerpc[i];}
  double GetParcoilspaceiron(int i) { return parcoilspaceiron[i];}
  double GetParairgap1(int i) { return parairgap1[i];}
  double GetParairgap2(int i) { return parairgap2[i];}
  double GetParairgap3(int i) { return parairgap3[i];}
  double GetParairgap4(int i) { return parairgap4[i];}
  double GetParlay(int i) { return parlay[i];}
  double GetParchm(int i) { return parchm[i];}
  double GetParirlay(int i) { return parirlay[i];}
  double GetParspacerA(int i) { return parspacerA[i];}
  double GetParspacerB(int i) { return parspacerB[i];}
  double GetParspacerC(int i) { return parspacerC[i];}
  double GetParspacerD(int i) { return parspacerD[i];}
  double GetParFrpBox(int i) { return parfrpbox[i];}
  double GetParAirBox(int i) { return parairbox[i];}
  double GetParG10Trap1(int i) { return parg10Trap1[i];}
  double GetParG10Trap2(int i) { return parg10Trap2[i];}
  // double GetParg10(int i) { return parg10[i];}
  double GetParal(int i) { return paral[i];}
  double GetParALCutBig(int i) { return paralCutBig[i];}
  double GetParALCutSmall(int i) { return paralCutSmall[i];}
  double GetParhoneycomb(int i) { return parhoneycomb[i];}
  double GetParHoneyCombCutBig(int i) { return parhoneycombCutBig[i];}
  double GetParHoneyCombCutSmall(int i) { return parhoneycombCutSmall[i];}
  double GetParcup(int i) { return parcup[i];}
  double GetParCupCutBig(int i) { return parcupCutBig[i];}
  double GetParCupCutSmall(int i) { return parcupCutSmall[i];}
  double GetParmylar(int i) { return parmylar[i];}
  double GetParMylarCutBig(int i) { return parmylarCutBig[i];}
  double GetParMylarCutSmall(int i) { return parmylarCutSmall[i];}
  double GetParcoat(int i) { return parcoat[i];}
  double GetParCoatCutBig(int i) { return parcoatCutBig[i];}
  double GetParCoatCutSmall(int i) { return parcoatCutSmall[i];}
  double GetParqurz(int i) { return parqurz[i];}
  double GetParQurzCutBig(int i) { return parqurzCutBig[i];}
  double GetParQurzCutSmall(int i) { return parqurzCutSmall[i];}
  double GetPargas(int i) { return pargas[i];}
  double GetParGasCutBig(int i) { return pargasCutBig[i];}
  double GetParGasCutSmall(int i) { return pargasCutSmall[i];}
  double GetParvcoil(int i) { return parvcoil[i];}
  double GetParhcoil(int i) { return parhcoil[i];}
  double GetParcurvedcoil(int i) { return parcurvedcoil[i];}
  double GetParcoilsupport(int i) { return parcoilsupport[i];}
  double GetXStrwd() { return Xstrwd;}
  double GetYStrwd() { return Ystrwd;}
  //double GetFeThickness(){ return nFeThickness;}
  //double GetAirGap_btwn_FeLayers(){ return nAirGap;}
  int GetnXStrip() { return nXStrip;} //AAR:
  int GetnYStrip() { return nYStrip;} //AAR:
  int GetnStack() { return nStack;}
  int GetnLayer() { return nLayer;}
  int GetnModule() { return nModule;}
  int GetnChamber() { return nChamber;} 
  int GetnIRLayer(){return  nIRLayer;}
  int GetnSpacerA(){return nSpacerA;}
  int GetnSpacerB(){return nSpacerB;}
  int GetnSpacerC(){return nSpacerC;}
  int GetnSpacerD(){return nSpacerD;}
  int GetnCoil(){return nCoil;}
  int GetNumino(){return numino;}
  double GetINOroomPos(int j){return INOroomPos[j];}
  double GetStackPosInRoom(int j){return StackPosInRoom[j];} 
  double GetRPCLayerPosZ(int j){return RPCLayerPosZ[j];}
  double GetIRONLayerPosZ(int j){return IRONLayerPosZ[j];}
  double GetRoomWallThicknessZ(){return RoomWallThicknessZ;}
  double GetRoomWallThickness(){return RoomWallThickness;}
  double GetLayerZdim(int j) {return LayerZdim[j];}
  double GetIronLayerZdim(int j) {return IronLayerZdim[j];}
  int GetnScintInUnit() {return nScintInUnit;}
  int GetnUnitTop() {return nUnitTop;}
  int GetnUnitWall() {return nUnitWall;}
  int GetnScintLayer() {return nScintLayer;}
  double GetPartopscint(int j) {return partopscint[j];}
  double GetParwallscint(int j) {return parwallscint[j];}
  double GetScintFromBottom() {return ScintFromBottom;}
  double GetScintUnitX() {return ScintUnitX;}
  double GetScintUnitY() {return ScintUnitY;}
  double GetScintUnitZ() {return ScintUnitZ;}
  double GetAirGapScintTop() {return AirGapScintTop;}
  double GetAirGapScintWall() {return AirGapScintWall;}

private : 
  int numino;
  double nFeThickness;
  float nAirGap ;

  int nXStrip;
  int nYStrip;

  double StackPosInRoom[3];
  double INOroomPos[3];
  double RPCLayerPosZ[10];
  double IRONLayerPosZ[11];

  // # of RPC Stackts
  int nStack;
  // # of layer
  int nLayer;
  // # of chamber
  int nChamber;
  // # ofcolumns of  module 
  int nModule; //8
  // # of iron layer
  int nIRLayer; //151 //will be calculated from nLayer

  int nSpacerA; //4
  int nSpacerB; //8
  int nSpacerC; //6
  int nSpacerD; //21
  //# of coils per module
  int nCoil; //4
  
  // INO detector size   
  float parworld[3];
  float parroom[3];
  float parairroom[3];
  float parairroom2[3];
  float parstaircaseair[3];
  float parstaircasel[3];
  float parstaircase[3];

  float parmagnet[3];
  float parino[3]; //={1606.01*cm, 706.01*cm, 598.0*cm};
  //space for coil in RPC layer
  float parcoilspacerpc[3];
  //space for coil in IRON layer
  float parcoilspaceiron[3];
  //space for horizontal air gap in IRON layer
  float parairgap1[3];
  //space for vertical air gap in IRON layer
  float parairgap2[3];
  //space for air gap between coils in IRON layer
  float parairgap3[3];
  //space for air gap between coils in IRON layer
  float parairgap4[3];
  // layer
  float parlay[3]; //={1600.0*cm, 700.02*cm, 4.26*cm};
  // chamber size
  float parchm[3]; //={100.0*cm,  100.0*cm, 4.24*cm};
  //iron layer
  float parirlay[3];

  //spacer of type A
  float parspacerA[3];
  //spacer of type B
  float parspacerB[3];
  //spacer of type C
  float parspacerC[3];
  //spacer of type D
  float parspacerD[3];

  //  FRP Tray which includes g10 electronics and RPC
  float parfrpbox[3];
  // Inner Walls of FRP Tray.
  float parairbox[3];
  // Triangle Electronics Board for DFE
  float parg10Trap1[8];
  // Triangle Electronics Board for H.V. supply
  float parg10Trap2[8];

  //Dimensions in (mm)

  //aluminium containing honeycomb
  float paral[3]; //={870.06, 917.56, 9.38}; 
  float paralCutBig[4];//={95, 136.35, 136.35, 1};
  float paralCutSmall[4];//={22.3, 33.537, 33.537, 1};
  
  //honeycomb containing air volume
  float parhoneycomb[3];//={870.05, 917.55, 9.23};
  float parhoneycombCutBig[4];//={95, 136.35, 136.35, 9.23};
  float parhoneycombCutSmall[4];//={22.3, 33.537, 33.537, 9.23};
  
  // cu plate which includes cu, g10, qrz and gas
  float parcup[3];//={870.04, 917.54, 4.23};
  float parcupCutBig[4];//={95, 136.35, 136.35, 4.23};
  float parcupCutSmall[4];//={22.3, 33.537, 33.537, 4.23};
  
  //mylar on graphite
  float parmylar[3];//={870.03, 917.53, 4.13};
  float parmylarCutBig[4];//={95, 136.35, 136.35, 4.13};
  float parmylarCutSmall[4];//={22.3, 33.537, 33.537, 4.13};
  
  // graphite coating on glass
  float parcoat[3];//={870.02, 917.52, 4.03};
  float parcoatCutBig[4];//={95, 136.35, 136.35, 4.03};
  float parcoatCutSmall[4];//={22.3, 33.537, 33.537, 4.03};
  
  // rpc glass which includes qrz and gas
  float parqurz[3];//={870.01, 917.51, 4};
  float parqurzCutBig[4];//={95, 136.35, 136.35, 4};
  float parqurzCutSmall[4];//={22.3, 33.537, 33.537, 4};
  
  // rpc gas module
  float pargas[3];//={870, 917.5, 1};
  float pargasCutBig[4];//={95, 136.35, 136.35, 1};
  float pargasCutSmall[4];//={22.3, 33.537, 33.537, 1};
  
  //coil in vertical direction
  float parvcoil[3];
  
  //coil in horizontal direction
  float parhcoil[3];
  
  //curved portion of coil
  double parcurvedcoil[5];
  
  //coil support
  float parcoilsupport[3];
  
  //  strip width
  float Xstrwd;// = 2.0*cm;
  float Ystrwd;// = 2.0*cm;
  float XYstrwd;

  float delta;

  float SpacerWidthAB; // = 20*mm;
  float SpacerWidthCD; // = 40*mm;
  float SpacerALength; // = 125*mm;
  float SpacerBLength; // = 250*mm;
  float SpacerCLength; // = 252.5*mm;
  float SpacerDLength; // = 252.5*mm;

  float RoomWallThickness;
  float RoomWallThicknessZ;
  float RoomXX;
  float RoomYY;
  float RoomZZ;
  float RoomXX2;
  float RoomYY2;
  float RoomZZ2;
  float StairCaseAirX;
  float StairCaseAirY;
  float StairCaseAirZ;
  float StairCaseLX;
  float StairCaseLY;
  float StairCaseLZ;
  float StairCaseX;
  float StairCaseY;
  float StairCaseZ;

  float ChamberSize; // = 1000*mm;

  float CutSmall; // = 919.99*mm;
  float CutBig; // = 919.99*mm;
  float ChamberLength; // = 1000*mm;

  float CoilHLength; // = 3782*mm;
  float CoilVLength; // = 7250*mm;
  float CoilThickness; // = 40*mm;
  float CoilWidth; // = 312.5*mm;
  float CurvedCoilInRadii; // = 178*mm;
  float CurvedCoilOutRadii; // = 258*mm;
  float CurvedCoilPhiMin; // = 0*rad;
  float CurvedCoilPhiMax; // = M_PI/2*rad;
  float CoilSupportLength; // = 100*mm;

  float AirGapIronLayer; // = 2.5*mm;

  //RPC
  float GasChamberSizeX; // = 870*mm; 
  float GasChamberSizeY; // = 917.5*mm;

  float GasGapThickness; // = 1*mm;
  float QurzThickness; // = 3*mm;
  float CoatThickness; // = 0.03*mm;
  float MylarThickness; // = 0.1*mm;
  float CopperThickness; // = 0.1*mm;
  float HoneyCombThickness; // = 5*mm;
  float AluminiumThickness; // = 0.150*mm;

  float BigCutRPC; // = 95*mm;
  float SmallCutRPC; // = 22.3*mm;
  float DeltaChm; // = 0.01*mm
  float DeltaCut; // = 2*mm;
  float FRPBoxXDim; // = 955*mm;
  float FRPBoxYDim; // = 970*mm;
  float FRPBoxZDim; // = 17*mm;
  float FRPThickness; // = 5*mm;
  float FRPZThickness;

  float G10Thickness; // = 2.5*mm;
  float G10TrapDX1; // = 0.01*mm;
  float G10Triangle1Length; // = 125.0*mm;
  float G10Triangle2Length; // = 100.0*mm;
  float G10TrapDY; // = 2.5*mm;

  // Scint Dimensions
  int nScintInUnit;// = 8;
  int nUnitTop;// = 11;
  int nUnitWall;// = 5;
  int nScintLayer;// = 3;
  float ScintUnitX;// = 5.0*cm;
  float ScintUnitY;// = 4.6*m;
  float ScintUnitZ;// = 1*cm;
  float AirGapScintTop;
  float AirGapScintWall;
  float ScintFromBottom;// = 300*mm;

  float partopscint[3];
  float parwallscint[3];

  //To calculate Shift of RPC from FRP Box Center
  float RPCShiftX; // = 119*mm;
  float RPCShiftY; // = 86*mm;
  float RPCShiftZ; // = 9*mm;
  double ShiftInRPCX;
  double ShiftInRPCY;
  double ShiftInRPCZ;

  double ShiftInX;
  double ShiftInY;
  double ShiftInZ[10];

  double AlShiftInAirBox;
  double AirBoxShiftInFRP;
  double FRPshiftInLayer[10];

  //To calculate pos of G10Trap1 in FRPBox
  float G10ShiftX1; // = 92*mm;
  float G10ShiftY1; // = 60*mm;

  //To calculate pos of G10Trap2 in FRPBox
  float G10ShiftX2; // = 13*mm;
  float G10ShiftY2; // = 1*mm;

  //Ellipsoid World SemiAxis;
  float WorldXX;
  float WorldYY;
  float WorldZZ;

  double LayerZdim[10];
  double IronLayerZdim[11];

  double LayerPosInStack[20];
  double INOMPositionGlobalX;
  double INOMPositionGlobalY;
  double INOMPositionGlobalZ;

};
#endif
