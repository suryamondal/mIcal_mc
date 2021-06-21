#ifndef micalDetectorConstruction_h
#define micalDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4Ellipsoid.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "micalCal0SD.hh"
#include "G4SDManager.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4MagIntegratorDriver.hh"

class G4Box;
class G4Trd;
class G4Ellipsoid;
//class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class micalDetectorMessenger;
class G4MaterialCutsCouple;
class micalElectroMagneticField;
class G4EqMagElectricField;
class G4MagInt_Driver;
class G4MagIntegratorStepper;
class G4ChordFinder;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class micalDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  micalDetectorConstruction();
  ~micalDetectorConstruction();

public:
     
  void SetUniformMagField(G4double);
     
  G4VPhysicalVolume* Construct();

  void UpdateGeometry();
  G4SubtractionSolid* ConstructRPCBox(float* parBox, float* parCutBig, float* parCutSmall);

     
public:
  
  void PrintCalorParameters(); 
                    
  G4Material* GetActiveMaterial()       {return ActiveMaterial;};

  const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
  G4VPhysicalVolume* Getphysi_World() {return physiWorld;};           
  const G4VPhysicalVolume* Getcal0()       {return physiGASR;};
  const G4VPhysicalVolume* Getcal1()        {return physiGASR;};
  const G4VPhysicalVolume* Getcal2()        {return physiGASR;};
  G4double GetAbsThickness(){return nFeThickness;}
  G4double GetGapThickness(){return nAirGap;}
  G4double GetStripWidth(){return Xstrwd;}
  //G4FieldManager*  GetLocalFieldManager() { return fLocalFieldManager ;}
  const G4Box*           InoINOM() {return solidINOM;};

private:
  
  G4Material*        ActiveMaterial;
  G4Material*        Iron;
  G4Material*        Copper;
  G4Material*        Aluminium;
  G4Material*        CoatMaterial;
  G4Material*        MylarMaterial;
  G4Material*        HoneyCombMaterial;
  G4Material*        G10;
  G4Material*        Air;
  G4Material*        SiO2;
  G4Material*        CarbonFRP;
  G4Material*        SF6;
  G4Material*        ConcreteMaterial;
  G4Material*        WorldMaterial;


  G4Box* solidWorld;
  G4LogicalVolume* logicWorld;
  G4VPhysicalVolume* physiWorld;
  
  G4Box* solidBuilding;
  G4LogicalVolume* logicBuilding;
  G4VPhysicalVolume* physiBuilding;

  G4Box* solidAirRoom;
  G4LogicalVolume* logicAirRoom;
  G4VPhysicalVolume* physiAirRoom;

  G4Box* solidAirRoomUp;
  G4LogicalVolume* logicAirRoomUp;
  G4VPhysicalVolume* physiAirRoomUp;

  G4Box* solidAirRoom2;
  G4LogicalVolume* logicAirRoom2;
  G4VPhysicalVolume* physiAirRoom2;

  G4Box* solidAirRoom2Up;
  G4LogicalVolume* logicAirRoom2Up;
  G4VPhysicalVolume* physiAirRoom2Up;

  G4Box* solidAirStairCase;
  G4LogicalVolume* logicAirStairCase;
  G4VPhysicalVolume* physiAirStairCase;

  G4Box* solidStairCase;
  G4LogicalVolume* logicStairCase;
  G4VPhysicalVolume* physiStairCase;

  G4Box* solidStairCaseL;
  G4LogicalVolume* logicStairCaseL;
  G4VPhysicalVolume* physiStairCaseL;

  G4Box* solidAirRoom3;
  G4LogicalVolume* logicAirRoom3;
  G4VPhysicalVolume* physiAirRoom3;

  G4Box* solidTopScint[3];
  G4LogicalVolume* logicTopScint[3];
  G4VPhysicalVolume* physiTopScint[3];
  
  G4Box* solidWallScintFront[3];
  G4LogicalVolume* logicWallScintFront[3];
  G4VPhysicalVolume* physiWallScintFront[3];

  G4Box* solidWallScintBack[3];
  G4LogicalVolume* logicWallScintBack[3];
  G4VPhysicalVolume* physiWallScintBack[3];

  G4Box* solidWallScintRight[3];
  G4LogicalVolume* logicWallScintRight[3];
  G4VPhysicalVolume* physiWallScintRight[3];

  G4Box* solidWallScintLeft[3];
  G4LogicalVolume* logicWallScintLeft[3];
  G4VPhysicalVolume* physiWallScintLeft[3];

  G4Box* solidINOM;
  G4LogicalVolume* logicINOM;
  G4VPhysicalVolume* physiINOM;
  
  G4Box* solidMagnet;
  G4LogicalVolume* logicMagnet;
  G4VPhysicalVolume* physiMagnet;
  
  G4SubtractionSolid* solidLAYE[10];
  G4LogicalVolume* logicLAYE[10];
  G4VPhysicalVolume* physiLAYE[10];
  
  G4SubtractionSolid* solidIRLAYE[11];
  G4LogicalVolume* logicIRLAYE[11]; 
  G4VPhysicalVolume* physiIRLAYE[11];
  
  G4Box* solidVCOIL;
  G4LogicalVolume* logicVCOIL;
  G4VPhysicalVolume* physiVCOIL;
  
  G4Box* solidHCOIL;
  G4LogicalVolume* logicHCOIL;
  G4VPhysicalVolume* physiHCOIL;
  
  G4Tubs* solidCurvedCOIL;
  G4LogicalVolume* logicCurvedCOIL;
  G4VPhysicalVolume* physiCurvedCOIL;

  G4Box* solidCOILSupport;
  G4LogicalVolume* logicCOILSupport;
  G4VPhysicalVolume* physiCOILSupport;
  
  G4Box* solidSpacerA[10];
  G4LogicalVolume* logicSpacerA[10];
  G4VPhysicalVolume* physiSpacerA[10];
  
  G4Box* solidSpacerB[10];
  G4LogicalVolume* logicSpacerB[10];
  G4VPhysicalVolume* physiSpacerB[10];
  
  G4Box* solidSpacerC[10];
  G4LogicalVolume* logicSpacerC[10];
  G4VPhysicalVolume* physiSpacerC[10];
  
  G4Box* solidSpacerD[10];
  G4LogicalVolume* logicSpacerD[10];
  G4VPhysicalVolume* physiSpacerD[10];
  
  G4Box* solidFRPBox;
  G4LogicalVolume* logicFRPBox;
  G4VPhysicalVolume* physiFRPBox;
  
  G4Trd* solidG10Trap1;
  G4LogicalVolume* logicG10Trap1;
  G4VPhysicalVolume* physiG10Trap1;
  
  G4Trd* solidG10Trap2;
  G4LogicalVolume*   logicG10Trap2;
  G4VPhysicalVolume* physiG10Trap2;

  G4Box* solidAirBox;
  G4LogicalVolume* logicAirBox;
  G4VPhysicalVolume* physiAirBox;

  G4SubtractionSolid* solidAL;
  G4LogicalVolume* logicAL;
  G4VPhysicalVolume* physiAL;

  G4SubtractionSolid* solidHoneyComb;
  G4LogicalVolume* logicHoneyComb;
  G4VPhysicalVolume* physiHoneyComb;

  G4SubtractionSolid* solidCUPL;
  G4LogicalVolume* logicCUPL;
  G4VPhysicalVolume* physiCUPL;
  
  G4SubtractionSolid* solidMYLAR;
  G4LogicalVolume* logicMYLAR;
  G4VPhysicalVolume* physiMYLAR;
  
  G4SubtractionSolid* solidCOAT;
  G4LogicalVolume* logicCOAT;
  G4VPhysicalVolume* physiCOAT;
  
  G4SubtractionSolid* solidQURZ;
  G4LogicalVolume* logicQURZ;
  G4VPhysicalVolume* physiQURZ;

  G4SubtractionSolid* solidGASR;
  G4LogicalVolume* logicGASR;
  G4VPhysicalVolume* physiGASR;

  G4UniformMagField* magField;
  G4SDManager* SDman;
  micalcal0SD * cal0SD;
  G4String  ecalName;
  G4Region* aRegion0;
  
  G4VisAttributes* visWhite;
  G4VisAttributes* visYellow;
  G4VisAttributes* visGray;
  G4VisAttributes* visCyan;
  G4VisAttributes* visBlue;
  G4VisAttributes* visRed;
  G4VisAttributes* visGreen;
  G4VisAttributes* visMagenta;
  G4VisAttributes* visNull;

  micalElectroMagneticField* inoical0Field;
  G4EqMagElectricField* fEquation;
  G4MagInt_Driver* pIntgrDriver;
  G4MagIntegratorStepper* pStepper;
  G4ChordFinder* pChordFinder;

protected:
  G4FieldManager*         fieldMgr;
  //  G4FieldManager*         fFieldManager ;
  //  G4FieldManager*         fLocalFieldManager ;
  
  
private:
  double nFeThickness, nAirGap,Xstrwd,B;
  void DefineMaterials();
  G4VPhysicalVolume* ConstructCalorimeter();     
  void ConstructFieldMap();
  
  //Variable to store the parameters
  float  parworld[3];
  float  parroom[3];
  float  parairroom[3];
  float parairroom2[3];
  float parstaircaseair[3];
  float parstaircasel[3];
  float parstaircase[3];
  float  parino[3];
  float parmagnet[3];
  float  parlay[3];
  float  parchm[3];
  float  parirlay[3];
  float  parcoilspacerpc[3];
  float  parcoilspaceiron[3];
  float  parairgap1[3];
  float  parairgap2[3];
  float  parairgap3[3];
  float  parairgap4[3];
  float  parspacerA[3];
  float  parspacerB[3];
  float  parspacerC[3];
  float  parspacerD[3];
  float parfrpbox[3];
  float parairbox[3];
  float parg10Trap1[8];
  float parg10Trap2[8];
  float paral[3];  
  float paralCutBig[4];
  float paralCutSmall[4];
  float parhoneycomb[3];
  float parhoneycombCutBig[4];
  float parhoneycombCutSmall[4];
  float parcup[3];
  float parcupCutBig[4];
  float parcupCutSmall[4];
  float parmylar[3];
  float parmylarCutBig[4];
  float parmylarCutSmall[4];
  float parcoat[3];
  float parcoatCutBig[4];
  float parcoatCutSmall[4];
  float parqurz[3];
  float parqurzCutBig[4];
  float parqurzCutSmall[4];
  float pargas[3];
  float pargasCutBig[4];
  float pargasCutSmall[4];
  float parvcoil[3];
  float parhcoil[3];
  float parcurvedcoil[5];
  float parcoilsupport[3];
  float RoomWallThickness;
  float RoomWallThicknessZ;

  double LayerZdim[10];
  double IronLayerZdim[11];

  double RPCLayerPosZ[10];
  double IRONLayerPosZ[11];
  
  double ShiftInX;
  double ShiftInY;
  double ShiftInZ[10];
  double AlShiftInAirBox;
  double AirBoxShiftInFRP;
  double FRPshiftInLayer[10];
  
  // Scint Dimensions
  int nScintInUnit;// = 8;
  int nUnitTop;// = 11;
  int nUnitWall;// = 5;
  int nScintLayer;// = 3;
  float ScintUnitX;// = 5.0*cm;
  float ScintUnitY;// = 4.6*m;
  float ScintUnitZ;// = 1*cm;
  float ScintFromBottom;// = 300*mm;
  float partopscint[3];
  float parwallscint[3];
  float postopcoil;
  float posmagnet;
  float AirGapScintTop;
  float AirGapScintWall;

  int nStack;
  int nLayer;
  int nChamber;
  int nIRLayer;
  int nSpacerA;
  int nSpacerB;
  int nSpacerC;
  int nSpacerD;
  int nCoil;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

