#include "micalDetectorParameterDef.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "micalParameterMessenger.hh"
// #include "DatabaseManager.hh"

using namespace std;

micalDetectorParameterDef *micalDetectorParameterDef::AnPointer;


micalDetectorParameterDef::micalDetectorParameterDef(){

  AnPointer = this;

  micalParameterMessenger* detectorConfig = new micalParameterMessenger();

  nFeThickness= 56*mm;
  nAirGap = 45*mm;
  XYstrwd = 30*mm;
  nLayer = 10; 
  nIRLayer = nLayer+1; 

  nChamber = 2;
  nModule =  2;
  nSpacerA = 4;
  nSpacerB = 2;
  nSpacerC = 4;
  nSpacerD = 2;
  
  nCoil = 2;
  nStack = 1;
  CoilThickness = 30*mm/2;
  CoilWidth = 744*mm/2;
  CoilHLength = 1600*mm/2;
  CoilVLength = 1700*mm/2;
  CurvedCoilInRadii = 185*mm;
  CurvedCoilOutRadii = 215*mm;
  CurvedCoilPhiMin = 0*rad;
  CurvedCoilPhiMax = M_PI/2*rad;
  CoilSupportLength = 100*mm;

  AirGapIronLayer = 2.5*mm;

  SpacerWidthAB = 40*mm/2;
  SpacerWidthCD = 80*mm/2;
  SpacerALength = 400*mm/2;
  SpacerBLength = 702*mm/2;
  SpacerCLength = 750.0*mm/2;
  SpacerDLength = 700.0*mm/2;

  ChamberSize = 1000*mm;
  GasChamberSizeX = 870*mm; 
  GasChamberSizeY = 917.5*mm;

  nScintInUnit = 8;
  nUnitTop = 11;
  nUnitWall = 5;
  nScintLayer = 3;
  ScintUnitX = 5.0*cm;
  ScintUnitY = 4.6*m;
  ScintUnitZ = 1*cm;
  AirGapScintTop = 10*mm;;
  AirGapScintWall = 10*mm;
  ScintFromBottom = 300*mm;
  
  GasGapThickness = 1*mm;
  QurzThickness = 3*mm;
  CoatThickness = 0.03*mm;
  MylarThickness = 0.1*mm;
  CopperThickness = 0.1*mm;
  HoneyCombThickness = 5*mm;
  AluminiumThickness = 0.150*mm;

  BigCutRPC = 95*mm;
  SmallCutRPC = 22.3*mm;
  DeltaChm = 0.01*mm;
  DeltaCut = 2*mm;

  FRPBoxXDim = 955*mm;
  FRPBoxYDim = 970*mm;
  FRPBoxZDim = 17*mm;
  FRPThickness = 5*mm;

  G10Thickness = 2.5*mm;
  G10TrapDX1 = 0.01*mm;
  G10Triangle1Length = 125.0*mm;
  G10Triangle2Length = 100.0*mm;
  G10TrapDY = 2.5*mm;

  //To calculate Shift of RPC from FRP Box Center
  RPCShiftX = 119*mm;
  RPCShiftY = 86*mm;
  RPCShiftZ = 9*mm;

  //To calculate pos of G10Trap1 in FRPBox
  G10ShiftX1 = 92*mm;
  G10ShiftY1 = 60*mm;

  //To calculate pos of G10Trap2 in FRPBox
  G10ShiftX2 = 13*mm;
  G10ShiftY2 = 1*mm;

  WorldXX = 1*km;
  WorldYY = 1*km;
  WorldZZ = 1*km;

  RoomWallThickness = 230*mm;
  RoomWallThicknessZ = 220*mm;
  
  RoomXX = 6*m;
  RoomYY = 12*m;
  RoomZZ = 3*m;

  RoomXX2 = 6*m;
  RoomYY2 = 18*m;
  RoomZZ2 = 3*m;

  StairCaseAirX = 6*m;
  StairCaseAirY = 3*m;
  StairCaseAirZ = RoomZZ;// + RoomWallThicknessZ;

  StairCaseX = 3.5*m;
  StairCaseY = StairCaseAirY/2 - 10*cm;
  StairCaseZ = 15*cm;

  StairCaseLX = 1.5*m;
  StairCaseLY = 3*m - 2*mm;;
  StairCaseLZ = RoomWallThicknessZ - 70*mm;
  
  LayerZdim[0] = 45*mm; LayerZdim[1] = 45*mm;
  LayerZdim[2] = 45*mm; LayerZdim[3] = 45*mm;
  LayerZdim[4] = 45*mm; LayerZdim[5] = 45*mm;
  LayerZdim[6] = 45*mm; LayerZdim[7] = 45*mm;
  LayerZdim[8] = 45*mm; LayerZdim[9] = 45*mm;
  
  IronLayerZdim[0] = 56*mm; IronLayerZdim[1] = 56*mm; 
  IronLayerZdim[2] = 56*mm; IronLayerZdim[3] = 56*mm; 
  IronLayerZdim[4] = 56*mm; IronLayerZdim[5] = 56*mm; 
  IronLayerZdim[6] = 56*mm; IronLayerZdim[7] = 56*mm; 
  IronLayerZdim[8] = 56*mm; IronLayerZdim[9] = 56*mm; 
  IronLayerZdim[10] = 56*mm; 
  
  IRONLayerPosZ[0] = -505*mm;

  for(int lyx=0; lyx<nLayer; lyx++) {
    RPCLayerPosZ[lyx] = IRONLayerPosZ[lyx] + IronLayerZdim[lyx]/2 + LayerZdim[lyx]/2;
    IRONLayerPosZ[lyx+1] = RPCLayerPosZ[lyx] + LayerZdim[lyx]/2 + IronLayerZdim[lyx+1]/2;
    cout<<"layerpos["<<lyx<<"] = "<<RPCLayerPosZ[lyx]<<" "<<IRONLayerPosZ[lyx+1]<<endl;
  }

  
  // ShiftInRPCStackX = 0.0*mm;
  // ShiftInRPCStackY = 540*mm;
  // ShiftInRPCStackZ = 475*mm;
  UpdateDetectorParameterDef();
}

void micalDetectorParameterDef::UpdateDetectorParameterDef(){
  // chamber size
  parchm[0] = ChamberSize; //100.0*cm;
  parchm[1] = ChamberSize; //100.0*cm;
  parchm[2] =nAirGap*0.5; //3.24*cm;
  
  // rpc gas rectangle(2mm)
  pargas[0] = GasChamberSizeX; //870*mm; //919.99*mm;
  pargas[1] = GasChamberSizeY; //917.5*mm;
  pargas[2] = GasGapThickness;

  // rpc gas Big Cut (2mm)
  pargasCutBig[0] = BigCutRPC;
  pargasCutBig[1] = sqrt(2)*pargasCutBig[0] + DeltaCut;
  pargasCutBig[2] = sqrt(2)*pargasCutBig[0] + DeltaCut;
  pargasCutBig[3] = GasGapThickness;

  // rpc gas Small Cut (2mm)
  pargasCutSmall[0] = SmallCutRPC;
  pargasCutSmall[1] = sqrt(2)*pargasCutSmall[0] + DeltaCut;
  pargasCutSmall[2] = sqrt(2)*pargasCutSmall[0] + DeltaCut;
  pargasCutSmall[3] = GasGapThickness;

  // rpc glass which includes qrz and gas (3mm)
  parqurz[0] =pargas[0] + DeltaChm;
  parqurz[1] =pargas[1] + DeltaChm;
  parqurz[2] =pargas[2] + QurzThickness; //0.4*cm;

  // qurz Big Cut 1 (3mm)
  parqurzCutBig[0] = BigCutRPC;
  parqurzCutBig[1] = sqrt(2)*parqurzCutBig[0] + DeltaCut; //HalfLength Along X
  parqurzCutBig[2] = sqrt(2)*parqurzCutBig[0] + DeltaCut; //HalfLength Along Y
  parqurzCutBig[3] = pargasCutBig[3] + QurzThickness; //HalfLength Along Z

  // qurz Small Cut 2 (3mm)
  parqurzCutSmall[0] = SmallCutRPC;
  parqurzCutSmall[1] = sqrt(2)*parqurzCutSmall[0] + DeltaCut; //HalfLength Along X
  parqurzCutSmall[2] = sqrt(2)*parqurzCutSmall[0] + DeltaCut; //HalfLength Along Y
  parqurzCutSmall[3] = pargasCutSmall[3] + QurzThickness; //HalfLength Along Z

  //graphite coat on glass (30um)
  parcoat[0] =pargas[0] + DeltaChm;
  parcoat[1] =pargas[1] + DeltaChm;
  parcoat[2] =parqurz[2] + CoatThickness; // 0.4030*cm;

  // graphite Big Cut 1 (30um)
  parcoatCutBig[0] = BigCutRPC;
  parcoatCutBig[1] = sqrt(2)*parcoatCutBig[0] + DeltaCut; //HalfLength Along X
  parcoatCutBig[2] = sqrt(2)*parcoatCutBig[0] + DeltaCut; //HalfLength Along Y
  parcoatCutBig[3] = parqurzCutBig[3] + CoatThickness; //HalfLength Along Z

  // graphite Small Cut 2 (30um)
  parcoatCutSmall[0] = SmallCutRPC;
  parcoatCutSmall[1] = sqrt(2)*parcoatCutSmall[0] + DeltaCut; //HalfLength Along X
  parcoatCutSmall[2] = sqrt(2)*parcoatCutSmall[0] + DeltaCut; //HalfLength Along Y
  parcoatCutSmall[3] = parqurzCutSmall[3] + CoatThickness; //HalfLength Along Z

  //mylar on graphite (100um)
  parmylar[0] =parcoat[0] + DeltaChm;
  parmylar[1] =parcoat[1] + DeltaChm;
  parmylar[2] =parcoat[2] + MylarThickness; // 0.4130*cm;

  // mylar Big Cut 1 (100um)
  parmylarCutBig[0] = BigCutRPC;
  parmylarCutBig[1] = sqrt(2)*parmylarCutBig[0] + DeltaCut; //HalfLength Along X
  parmylarCutBig[2] = sqrt(2)*parmylarCutBig[0] + DeltaCut; //HalfLength Along Y
  parmylarCutBig[3] = parcoatCutBig[3] + MylarThickness; //HalfLength Along Z

  // mylar Small Cut 2 (100um)
  parmylarCutSmall[0] = SmallCutRPC;
  parmylarCutSmall[1] = sqrt(2)*parmylarCutSmall[0] + DeltaCut; //HalfLength Along X
  parmylarCutSmall[2] = sqrt(2)*parmylarCutSmall[0] + DeltaCut; //HalfLength Along Y
  parmylarCutSmall[3] = parcoatCutSmall[3] + MylarThickness; //HalfLength Along Z

  // copper strip (100um)
  parcup[0]=parmylar[0] + DeltaChm;
  parcup[1]=parmylar[1] + DeltaChm;
  parcup[2]=parmylar[2] + CopperThickness; //0.4230*cm;
  
  // copper Big Cut 1 (100um)
  parcupCutBig[0] = BigCutRPC;
  parcupCutBig[1] = sqrt(2)*parcupCutBig[0] + DeltaCut; //HalfLength Along X
  parcupCutBig[2] = sqrt(2)*parcupCutBig[0] + DeltaCut; //HalfLength Along Y
  parcupCutBig[3] = pargasCutBig[3] + CopperThickness; //HalfLength Along Z

  // copper Small Cut 2 (100um)
  parcupCutSmall[0] = SmallCutRPC;
  parcupCutSmall[1] = sqrt(2)*parcupCutSmall[0] + DeltaCut; //HalfLength Along X
  parcupCutSmall[2] = sqrt(2)*parcupCutSmall[0] + DeltaCut; //HalfLength Along Y
  parcupCutSmall[3] = parmylarCutSmall[3] + CopperThickness; //HalfLength Along Z

  // honeycomb containing air volume (5mm)
  parhoneycomb[0]=parcup[0] + DeltaChm;
  parhoneycomb[1]=parcup[1] + DeltaChm;
  parhoneycomb[2]=parcup[2] + HoneyCombThickness; //0.9230*cm;
  
  // honeycomb Big Cut 1 (5mm)
  parhoneycombCutBig[0] = BigCutRPC;
  parhoneycombCutBig[1] = sqrt(2)*parhoneycombCutBig[0] + DeltaCut; //HalfLength Along X
  parhoneycombCutBig[2] = sqrt(2)*parhoneycombCutBig[0] + DeltaCut; //HalfLength Along Y
  parhoneycombCutBig[3] = parcupCutBig[3] + HoneyCombThickness; //HalfLength Along Z

  // honeycomb Small Cut 2 (5mm)
  parhoneycombCutSmall[0] = SmallCutRPC;
  parhoneycombCutSmall[1] = sqrt(2)*parhoneycombCutSmall[0] + DeltaCut; //HalfLength Along X
  parhoneycombCutSmall[2] = sqrt(2)*parhoneycombCutSmall[0] + DeltaCut; //HalfLength Along Y
  parhoneycombCutSmall[3] = parcupCutSmall[3] + HoneyCombThickness; //HalfLength Along Z

  // aluminium containing honeycomb (150um)
  paral[0]=parhoneycomb[0] + DeltaChm;
  paral[1]=parhoneycomb[1] + DeltaChm;
  paral[2]=parhoneycomb[2] + AluminiumThickness; //0.9380*cm;
  
  // aluminium Big Cut 1 (150um)
  paralCutBig[0] = BigCutRPC;
  paralCutBig[1] = sqrt(2)*paralCutBig[0] + DeltaCut; //HalfLength Along X
  paralCutBig[2] = sqrt(2)*paralCutBig[0] + DeltaCut; //HalfLength Along Y
  paralCutBig[3] = parhoneycombCutBig[3] + AluminiumThickness; //HalfLength Along Z

  // aluminium Small Cut 2 (150um)
  paralCutSmall[0] = SmallCutRPC;
  paralCutSmall[1] = sqrt(2)*paralCutSmall[0] + DeltaCut; //HalfLength Along X
  paralCutSmall[2] = sqrt(2)*paralCutSmall[0] + DeltaCut; //HalfLength Along Y
  paralCutSmall[3] = parhoneycombCutSmall[3] + AluminiumThickness; //HalfLength Along Z

  // FRP Tray Dimensions
  parfrpbox[0] = FRPBoxXDim; //955*mm; 
  parfrpbox[1] = FRPBoxYDim; //970*mm; 
  parfrpbox[2] = FRPBoxZDim; //17*mm;

  //  Domensions of InnerWall of FRP Tray. Airbox in which G10 Trapezoid and Aluminium are placed
  parairbox[0] = parfrpbox[0] - FRPThickness;
  parairbox[1] = parfrpbox[1] - FRPThickness;
  parairbox[2] = parfrpbox[2] - FRPThickness + 0.5*mm;

  // The position of RPC w.r.t. the inner walls of FRP Tray. Shift to be added in the data while reconstucting the position of hit.
  ShiftInX = /*119*mm*/ RPCShiftX - parairbox[0] + paral[0];
  ShiftInY = /*86*mm*/ RPCShiftY - parairbox[1] + paral[1];
  
  AlShiftInAirBox = -parairbox[2]+paral[2];
  AirBoxShiftInFRP = parfrpbox[2] - parairbox[2];
  for(int ij=0; ij<nLayer; ij++) {
      FRPshiftInLayer[ij] = -(LayerZdim[ij]/2) + parfrpbox[2];
      ShiftInZ[ij] = AlShiftInAirBox + AirBoxShiftInFRP + FRPshiftInLayer[ij];
  }

  // ShiftInZ = /*9*mm*/ RPCShiftZ - parfrpbox[2] + parairbox[2];

  //g10 trapzoid 1 (Triangle board for DFE)
  parg10Trap1[0] = G10Triangle1Length; //125*mm; //Side Length of Triangle
  parg10Trap1[1] = G10TrapDX1; //0.01*mm; //(Half-length along x at the surface positioned at -dz)
  parg10Trap1[2] = sqrt(2)*parg10Trap1[0]/2; //Hypotenuse of triangle (Half-length along x at the surface positioned at +dz)
  parg10Trap1[3] = G10Thickness; //2.5*mm; // Thickness of the board (Half-length along y)
  parg10Trap1[4] = (sqrt((parg10Trap1[0]*parg10Trap1[0])-(parg10Trap1[0]*parg10Trap1[0]/2.))/2.); //(Half-length along z)
  parg10Trap1[5] = /*92*mm*/ G10ShiftX1 - parairbox[0] + parg10Trap1[0]/2; // xpos wrt center of airbox
  parg10Trap1[6] = /*60*mm*/ G10ShiftY1 - parairbox[1] + parg10Trap1[0]/4; // ypos wrt center of airbox
  parg10Trap1[7] = -parairbox[2] + parg10Trap1[3]; // zpos wrt center of airbox
  
  //g10 trapzoid 2 (Triangle board for HV power supply)
  parg10Trap2[0] = G10Triangle2Length; //100*mm; //Side Length of Triangle
  parg10Trap2[1] = G10TrapDX1; //0.01*mm; //(Half-length along x at the surface positioned at -dz)
  parg10Trap2[2] = sqrt(2)*parg10Trap2[0]/2; //Hypotenuse of triangle (Half-length along x at the surface positioned at +dz)
  parg10Trap2[3] = G10Thickness; //2.5*mm; // Thickness of the board (Half-length along y)
  parg10Trap2[4] = (sqrt((parg10Trap2[0]*parg10Trap2[0])-(parg10Trap2[0]*parg10Trap2[0]/2.))/2.); //(Half-length along z)
  parg10Trap2[5] = parairbox[0] - parg10Trap2[0]/2 - G10ShiftX2; //13*mm; // xpos wrt center of airbox
  parg10Trap2[6] = parairbox[0] - parg10Trap2[0]/4 - G10ShiftY2; //1*mm; // ypos wrt center of airbox
  parg10Trap2[7] = -parairbox[2] + parg10Trap2[3]; // zpos wrt center of airbox
  
  //spacerA
  parspacerA[0]=SpacerWidthAB; //2*cm;
  parspacerA[1]=SpacerALength; //12.5*cm;
  parspacerA[2]=nAirGap*0.5;
  
  //spacerB
  parspacerB[0]=SpacerWidthAB; //2*cm;
  parspacerB[1]=SpacerBLength; //25*cm;
  parspacerB[2]=nAirGap*0.5;
  
  //spacerC
  parspacerC[0]=SpacerWidthCD; //2*cm;
  parspacerC[1]=SpacerCLength; //25.25*cm;
  parspacerC[2]=nAirGap*0.5;
  
  //spacerD
  parspacerD[0]=SpacerWidthCD; //4*cm;
  parspacerD[1]=SpacerDLength; //25.25*cm;
  parspacerD[2]=nAirGap*0.5;
  
  // AIR in ROOM Dimensions;
  parairroom[0] = RoomXX/2;
  parairroom[1] = RoomYY/2;
  parairroom[2] = RoomZZ/2;

  // AIR in ROOM Dimensions;
  parairroom2[0] = RoomXX2/2;
  parairroom2[1] = RoomYY2/2;
  parairroom2[2] = RoomZZ2/2;
 
  // StairCase Air
  parstaircaseair[0] = StairCaseAirX/2;
  parstaircaseair[1] = StairCaseAirY/2;
  parstaircaseair[2] = StairCaseAirZ/2;

  // StairCase Landing
  parstaircasel[0] = StairCaseLX/2;
  parstaircasel[1] = StairCaseLY/2;
  parstaircasel[2] = StairCaseLZ/2;

  // StairCase 
  parstaircase[0] = StairCaseX/2;
  parstaircase[1] = StairCaseY/2;
  parstaircase[2] = StairCaseZ/2;

  // ROOM Dimensions;
  parroom[0] = parairroom[0] + RoomWallThickness;
  parroom[1] = parairroom[1] + parairroom2[1] + 3*parstaircaseair[1] + 3*RoomWallThickness;
  parroom[2] = 2*parairroom[2] + 1.5*RoomWallThicknessZ;
  
  //Cuboid World Size
  parworld[0] = WorldXX;
  parworld[1] = WorldYY;
  parworld[2] = parroom[2]*1.01;
 
  //coil in vertical direction
  parvcoil[0] = CoilThickness; //4*cm;
  parvcoil[1] = CoilWidth; //31.25*cm;
  parvcoil[2] = (nLayer*nAirGap*0.5)+(nIRLayer*nFeThickness*0.5) + 100*mm; //CoilVLength; //725*cm ;
  
  //coil in horizontal direction
  parhcoil[0]=CoilHLength; //378.2*cm;
  parhcoil[1]=CoilWidth; //31.25*cm;
  parhcoil[2]=CoilThickness; //4*cm ;
  
  //curved portion of coil
  parcurvedcoil[0]= CurvedCoilInRadii; //17.8*cm;
  parcurvedcoil[1]= CurvedCoilOutRadii; //25.8*cm;
  parcurvedcoil[2]= CoilWidth; //31.25*cm ;
  parcurvedcoil[3]= CurvedCoilPhiMin; //0*rad ;
  parcurvedcoil[4]= CurvedCoilPhiMax; //M_PI/2*rad ;

  //Top Scintillator
  partopscint[0] = nUnitTop*nScintInUnit*ScintUnitX/2;
  partopscint[1] = ScintUnitY/2;
  partopscint[2] = ScintUnitZ/2;
  
  //Wall Scintillator
  parwallscint[0] = nUnitWall*nScintInUnit*ScintUnitX/2;
  parwallscint[1] = ScintUnitY/2;
  parwallscint[2] = ScintUnitZ/2;
  
  //Magnet
  parmagnet[0] = parchm[0]*2 + 2*mm; // 1606.01*cm;
  parmagnet[1] = parchm[1]*2 + 2*mm; //706.01*cm;
  parmagnet[2] = parvcoil[2] + CurvedCoilOutRadii + 5*mm;//parairroom[2];

  //INO
  parino[0] = 1.1*(partopscint[0]+ScintUnitZ); // 1606.01*cm;
  parino[1] = 1.1*(partopscint[1]+ScintUnitZ); //706.01*cm;
  parino[2] = parairroom[2];//1.01*(parwallscint[0] + ScintFromBottom/2);
  
  // RPC layer
  parlay[0]=parchm[0]*2; //800.0*cm; //1600.0*cm;
  parlay[1]=parchm[1]*2; // 700.02*cm;
  parlay[2]=nAirGap*0.5;// 3.26*cm;
  
  //space for coil in RPC layer
  parcoilspacerpc[0] = CoilThickness + 25*mm;
  parcoilspacerpc[1] = CoilWidth + 25*mm;
  parcoilspacerpc[2] = parlay[2];
  
  // iron layer
  parirlay[0]=parlay[0];
  parirlay[1]=parlay[1];
  parirlay[2]=nFeThickness*0.5;
  
  //space for coil in iron layer
  parcoilspaceiron[0] = CoilThickness + 25*mm;
  parcoilspaceiron[1] = CoilWidth + 25*mm;
  parcoilspaceiron[2] = parirlay[2];
  
  //space for vertical air gap in iron layer
  parairgap1[0] = AirGapIronLayer;
  parairgap1[1] = 803*mm - parcoilspaceiron[1] - 1*mm;
  parairgap1[2] = parirlay[2];
  
  //space for vertical air gap in iron layer
  parairgap2[0] = AirGapIronLayer;
  parairgap2[1] = (parlay[1] - 2*parcoilspaceiron[1] - parairgap1[1])*0.5 - 1*mm;
  parairgap2[2] = parirlay[2];
  
  //space for air gap between coil in iron layer
  parairgap3[0] = (parchm[0] - AirGapIronLayer)/2 - 1*mm;
  parairgap3[1] = AirGapIronLayer;
  parairgap3[2] = parirlay[2];
  
  //space for air gap between coil in iron layer
  parairgap4[0] = parchm[0] - parcoilspaceiron[0] - 1*mm;
  parairgap4[1] = AirGapIronLayer;
  parairgap4[2] = parirlay[2];
  
  //coil support
  parcoilsupport[0]=CoilSupportLength; //10*cm;
  parcoilsupport[1]=CoilWidth; //31.25*cm;
  parcoilsupport[2]=CoilSupportLength; //10*cm ;
  
  //  strip width
  Xstrwd = XYstrwd*mm; //100 strips
  Ystrwd = XYstrwd*mm; // 100 strips
  nXStrip =  int(1.999*pargas[0]/Xstrwd)+1; //parchm[0]*2/Xstrwd; //GMA14
  nYStrip =  int(1.999*pargas[1]/Ystrwd)+1; //parchm[1]*2/Xstrwd;

  
  // for(int ij=0; ij<nLayer; ij++) {
  //   RPCLayerPosZ[ij] = -((int(0.5*nIRLayer)-ij)*2*(parirlay[2]+ parlay[2]))+(parirlay[2]+ parlay[2]);
  // }

  // for(int ij=0; ij<nIRLayer; ij++) {
  //   IRONLayerPosZ[ij] = -( (int(0.5*nIRLayer)-ij)*2*(parirlay[2]+ parlay[2]) );
  // }

  StackPosInRoom[0] = 0.0*mm;
  StackPosInRoom[1] = parairroom[1] - parino[1] - 2000*mm;
  StackPosInRoom[2] = 0.0;//-parairroom[2] + parino[2] + 5*mm;

  INOroomPos[0] = 0.0*mm;
  INOroomPos[1] = -parroom[1] + 2*RoomWallThickness + parairroom[1] + 2*parstaircaseair[1];
  INOroomPos[2] = -parroom[2] + RoomWallThicknessZ + parairroom[2];
  
  // cout<<"StackPosInRoom "<<StackPosInRoom[0]<<" "<<StackPosInRoom[1]<<" "<<StackPosInRoom[2]<<endl;
  // cout<<"INOroomPos "<<INOroomPos[0]<<" "<<INOroomPos[1]<<" "<<INOroomPos[2]<<endl;
  // for(int ij=0; ij<nLayer; ij++) {
  //   cout<<"Layer["<<ij<<"] = "<<RPCLayerPosZ[ij]<<endl;
  // }
  // for(int ij=0; ij<nIRLayer; ij++) {
  //   cout<<"Layer["<<ij<<"] = "<<IRONLayerPosZ[ij]<<endl;
  // }

}

