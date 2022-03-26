//#define uniformField
#define arbField 

#include "micalDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Trd.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
//#include "G4PVReplica.hh"
#include "G4PVDivision.hh"

#ifdef uniformField
#include "G4UniformMagField.hh"
#endif

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4RotationMatrix.hh"

#include "G4UniformMagField.hh"
//#include "micalMagneticField.hh"
#include "micalElectroMagneticField.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4PropagatorInField.hh"
#include "G4TransportationManager.hh"

#include "G4ClassicalRK4.hh"
#include "G4SimpleRunge.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleHeum.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4GDMLParser.hh"
#define debug

#ifdef debug
#include "G4Timer.hh"
#endif

#include "G4ProductionCuts.hh"
#include "G4MaterialCutsCouple.hh"
#include "micalDetectorParameterDef.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalDetectorConstruction::micalDetectorConstruction()
  :ActiveMaterial(0),
   solidWorld(0),logicWorld(0),physiWorld(0),
   solidBuilding(0),logicBuilding(0),physiBuilding(0),
  solidAirRoom(0),logicAirRoom(0),physiAirRoom(0),
  solidAirRoomUp(0),logicAirRoomUp(0),physiAirRoomUp(0),
  solidAirRoom2(0),logicAirRoom2(0),physiAirRoom2(0),
  solidAirRoom2Up(0),logicAirRoom2Up(0),physiAirRoom2Up(0),
  solidAirStairCase(0), logicAirStairCase(0), physiAirStairCase(0), 
  solidStairCase(0), logicStairCase(0), physiStairCase(0), 
   solidStairCaseL(0), logicStairCaseL(0), physiStairCaseL(0), 
  solidAirRoom3(0),logicAirRoom3(0),physiAirRoom3(0),
  solidINOM(0),logicINOM(0), physiINOM(0),
   solidVCOIL(0), logicVCOIL(0), physiVCOIL(0),
   solidHCOIL(0), logicHCOIL(0), physiHCOIL(0),
   solidCurvedCOIL(0), logicCurvedCOIL(0), physiCurvedCOIL(0),
   solidCOILSupport(0), logicCOILSupport(0), physiCOILSupport(0),
   solidFRPBox(0),logicFRPBox(0),physiFRPBox(0),
   solidG10Trap1(0),logicG10Trap1(0),physiG10Trap1(0),
   solidG10Trap2(0),logicG10Trap2(0),physiG10Trap2(0),
   solidAirBox(0),logicAirBox(0),physiAirBox(0),
   solidAL(0), logicAL(0), physiAL(0),
   solidHoneyComb(0), logicHoneyComb(0), physiHoneyComb(0),
   solidCUPL(0), logicCUPL(0), physiCUPL(0),
   solidMYLAR(0), logicMYLAR(0), physiMYLAR(0),
   solidCOAT(0), logicCOAT(0), physiCOAT(0),
   solidQURZ(0), logicQURZ(0), physiQURZ(0),
   solidGASR(0), logicGASR(0), physiGASR(0),
   magField(0),
  visWhite(0),
  visYellow(0),
  visGray(0),
  visCyan(0),
  visBlue(0),
  visRed(0),
  visGreen(0),
  visMagenta(0),
  visNull(0),
   inoical0Field(0),
   fEquation(0),
   pIntgrDriver(0),
   pStepper(0),
   pChordFinder(0) {
  
  //#include "micalDetectorParameterDef.hh"
  
  // materials
  DefineMaterials();
  
  SDman = G4SDManager::GetSDMpointer();
  ecalName = "mical0cal0";
  cal0SD = new micalcal0SD(ecalName);
  SDman->AddNewDetector(cal0SD);
  
  aRegion0 = new G4Region("Calor_EBlock");
  
  // create commands for interactive definition of the calorimeter
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalDetectorConstruction::~micalDetectorConstruction() {
  delete   SF6; SF6=0;
  delete G10 ; G10=0;
  delete ActiveMaterial;  ActiveMaterial=0;
  delete visWhite;      visWhite=0;
  delete visYellow;      visYellow=0;
  delete visGray;      visGray=0;
  delete visCyan;      visCyan=0;
  delete visBlue;      visBlue=0;
  delete visRed;      visRed=0;
  delete visGreen;      visGreen=0;
  delete visMagenta;      visMagenta=0;
  delete visNull;       visNull=0;
  delete inoical0Field; inoical0Field=0;
  delete fEquation;      fEquation=0;
  delete pStepper;      pStepper=0;   
  delete pChordFinder; pChordFinder=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* micalDetectorConstruction::Construct() {
  
#ifdef uniformField
  SetUniformMagField(0*tesla);
#else
  
#ifdef arbField
  ConstructFieldMap();
#endif
#endif
  
  //micalDetectorParameterDef* paradef = micalDetectorParameterDef::AnPointer;
  G4VPhysicalVolume* physiworld;
  
  physiworld = ConstructCalorimeter();
  return physiworld;
  
  //return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalDetectorConstruction::DefineMaterials() { 
  
  
  //This function illustrates the possible ways to define materials
  
  G4String symbol;             //a=mass of a mole;
  G4double density;      //z=mean number of protons;  
  // n=number of nucleons in an isotope;
  
  G4int ncomponents,natoms;
  G4double fractionmass;
  
  G4NistManager* mat = G4NistManager::Instance();  
  mat->SetVerbose(1);
  // define Elements
  G4Element* C = mat->FindOrBuildElement("C");
  G4Element* S = mat->FindOrBuildElement("S");
  G4Element* F = mat->FindOrBuildElement("F");
  G4Element* H = mat->FindOrBuildElement("H");
  G4Element* O = mat->FindOrBuildElement("O");
  G4Element* Si = mat->FindOrBuildElement("Si");
  
  SF6 =  new G4Material("SF6",density = 0.006164*g/cm3, ncomponents=2);
  SF6->AddElement(S, natoms=1);
  SF6->AddElement(F, natoms=6);
 
  CarbonFRP = new G4Material("FRPCarbon",density = 3.*g/cm3,ncomponents=1);
  CarbonFRP->AddElement(C, fractionmass=100.*perCent);

  G4Material* Butane =   mat->FindOrBuildMaterial("G4_BUTANE");
  G4Material* Freon =   mat->FindOrBuildMaterial("G4_FREON-13B1");
 
  // G4Material* tmpFE = mat->FindOrBuildMaterial("G4_Fe");
  // Iron =  new G4Material("G4_Fe",0.0001*g/cm3,tmpFE);
  Iron = mat->FindOrBuildMaterial("G4_Fe");
  
  Copper =  mat->FindOrBuildMaterial("G4_Cu");
  Aluminium =  mat->FindOrBuildMaterial("G4_Al");
  CoatMaterial = mat->FindOrBuildMaterial("G4_GRAPHITE");
  MylarMaterial = mat->FindOrBuildMaterial("G4_MYLAR");
  HoneyCombMaterial =  mat->FindOrBuildMaterial("G4_POLYETHYLENE");
  Air = mat->FindOrBuildMaterial("G4_AIR");

  ConcreteMaterial = mat->FindOrBuildMaterial("G4_CONCRETE");

  // G4Material* tmpSiO2 = mat->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  // SiO2 = new G4Material("G4_SILICON_DIOXIDE",0.00001*g/cm3,tmpSiO2);
  SiO2 = mat->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

  //define materials by fractional mass
  
  G10 = new G4Material("G10", density= 1.09*g/cm3, ncomponents=4); // NemaG10
  // G10 = new G4Material("G10", density= 1.700*g/cm3, ncomponents=4); // NemaG10
  G10->AddElement(Si, natoms=1);
  G10->AddElement(O , natoms=2);
  G10->AddElement(C , natoms=3);
  G10->AddElement(H , natoms=3);
  
  // G10 = new G4Material("G10",  z=8.455, a=16.91*g/mole, density=1.09*g/cm3);

  ActiveMaterial = new G4Material("rpcgas",density= 0.00418*g/cm3, ncomponents=3);
  ActiveMaterial->AddMaterial(Freon, fractionmass=95.15*perCent);
  ActiveMaterial->AddMaterial(Butane, fractionmass=4.51*perCent);
  ActiveMaterial->AddMaterial(SF6, fractionmass=0.34*perCent);

  WorldMaterial = new G4Material("EarthSiO2",density=2.65*g/cm3, ncomponents=2);
  WorldMaterial->AddElement(Si, natoms=1);
  WorldMaterial->AddElement(O , natoms=2);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* micalDetectorConstruction::ConstructCalorimeter() {
  // cout <<"G4VPhysicalVolume* micalDetectorConstruction::ConstructCalorimeter()"<<endl;
  // Clean old geometry, if any//
  
  cout<<" ConstructCalorimeter  "<<endl;
  system("free");

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  // Get detector definitions from micalDetectorParameterDef;
  micalDetectorParameterDef* paradef = micalDetectorParameterDef::AnPointer;
  
  for (int ij=0; ij<3; ij++) {
    parworld[ij] = paradef->GetParworld(ij);
    parroom[ij] = paradef->GetParroom(ij);
    parairroom[ij] = paradef->GetParairroom(ij); 
    parairroom2[ij] = paradef->GetParairroom2(ij); 
    parstaircaseair[ij] = paradef->GetParstaircaseair(ij); 
    parstaircasel[ij] = paradef->GetParstaircasel(ij); 
    parstaircase[ij] = paradef->GetParstaircase(ij); 
    parino[ij] = paradef->GetParino(ij);
    parmagnet[ij] = paradef->GetParmagnet(ij);
    parlay[ij] = paradef->GetParlay(ij);
    parchm[ij] = paradef->GetParchm(ij);
    parirlay[ij] = paradef->GetParirlay(ij);
    parcoilspacerpc[ij]= paradef->GetParcoilspacerpc(ij);
    parcoilspaceiron[ij]=paradef->GetParcoilspaceiron(ij);
    parairgap1[ij] = paradef->GetParairgap1(ij);
    parairgap2[ij] = paradef->GetParairgap2(ij);
    parairgap3[ij] = paradef->GetParairgap3(ij);
    parairgap4[ij] = paradef->GetParairgap4(ij);
    parspacerA[ij] = paradef->GetParspacerA(ij);
    parspacerB[ij] = paradef->GetParspacerB(ij);
    parspacerC[ij] = paradef->GetParspacerC(ij);
    parspacerD[ij] = paradef->GetParspacerD(ij);
    parfrpbox[ij] = paradef->GetParFrpBox(ij);
    parairbox[ij] = paradef->GetParAirBox(ij);
    parvcoil[ij] = paradef->GetParvcoil(ij);
    parhcoil[ij] = paradef->GetParhcoil(ij);
    parcoilsupport[ij] = paradef->GetParcoilsupport(ij);
    paral[ij] = paradef->GetParal(ij);
    parhoneycomb[ij] = paradef->GetParhoneycomb(ij);
    parcup[ij] = paradef->GetParcup(ij);
    parmylar[ij] = paradef->GetParmylar(ij);
    parcoat[ij] = paradef->GetParcoat(ij);
    parqurz[ij] = paradef->GetParqurz(ij);
    pargas[ij] = paradef->GetPargas(ij);
    partopscint[ij] = paradef->GetPartopscint(ij);
    parwallscint[ij] = paradef->GetParwallscint(ij);
  }

  for (int ij=0; ij<4; ij++) {
    paralCutBig[ij] = paradef->GetParALCutBig(ij);
    paralCutSmall[ij] = paradef->GetParALCutSmall(ij);
    parhoneycombCutBig[ij] = paradef->GetParHoneyCombCutBig(ij);
    parhoneycombCutSmall[ij] = paradef->GetParHoneyCombCutSmall(ij);
    parcupCutBig[ij] = paradef->GetParCupCutBig(ij);  
    parcupCutSmall[ij] = paradef->GetParCupCutSmall(ij);
    parmylarCutBig[ij] = paradef->GetParMylarCutBig(ij);
    parmylarCutSmall[ij] = paradef->GetParMylarCutSmall(ij);
    parcoatCutBig[ij] = paradef->GetParCoatCutBig(ij);
    parcoatCutSmall[ij] = paradef->GetParCoatCutSmall(ij);
    parqurzCutBig[ij] = paradef->GetParQurzCutBig(ij);
    parqurzCutSmall[ij] = paradef->GetParQurzCutSmall(ij);
    pargasCutBig[ij] = paradef->GetParGasCutBig(ij);
    pargasCutSmall[ij] = paradef->GetParGasCutSmall(ij);
  }

  for (int ij=0; ij<5; ij++) {parcurvedcoil[ij] = paradef->GetParcurvedcoil(ij);}
  for (int ij=0; ij<8; ij++) {parg10Trap1[ij] = paradef->GetParG10Trap1(ij);}
  for (int ij=0; ij<8; ij++) {parg10Trap2[ij] = paradef->GetParG10Trap2(ij);}
  
  nLayer = paradef->GetnLayer();
  nChamber = paradef->GetnChamber();
  nStack = paradef->GetnStack();
  nIRLayer = paradef->GetnIRLayer();
  nSpacerA = paradef->GetnSpacerA();
  nSpacerB = paradef->GetnSpacerB();
  nSpacerC = paradef->GetnSpacerC();
  nSpacerD = paradef->GetnSpacerD();
  nCoil = paradef->GetnCoil();
  RoomWallThicknessZ = paradef->GetRoomWallThicknessZ();
  RoomWallThickness = paradef->GetRoomWallThickness();

  nScintInUnit = paradef->GetnScintInUnit();
  nUnitTop = paradef->GetnUnitTop();
  nUnitWall = paradef->GetnUnitWall();
  nScintLayer = paradef->GetnScintLayer();
  ScintUnitX = paradef->GetScintUnitX();
  ScintUnitY = paradef->GetScintUnitY();
  ScintUnitZ = paradef->GetScintUnitZ();
  AirGapScintTop = paradef->GetAirGapScintTop();
  AirGapScintWall = paradef->GetAirGapScintWall();
  ScintFromBottom = paradef->GetScintFromBottom();

  ShiftInX = paradef->GetShiftInX();
  ShiftInY = paradef->GetShiftInY();
  AlShiftInAirBox = paradef->GetAlShiftInAirBox();
  AirBoxShiftInFRP = paradef->GetAirBoxShiftInFRP();
  for(int ij=0; ij<nIRLayer; ij++) {
    IronLayerZdim[ij] = paradef->GetIronLayerZdim(ij);
    IRONLayerPosZ[ij] = paradef->GetIRONLayerPosZ(ij);
  }

  for(int ij=0; ij<nLayer; ij++) {
    LayerZdim[ij] = paradef->GetLayerZdim(ij);
    RPCLayerPosZ[ij] = paradef->GetRPCLayerPosZ(ij);
    ShiftInZ[ij] = paradef->GetShiftInZ(ij);
    FRPshiftInLayer[ij] = paradef->GetFRPshiftInLayer(ij);
  }

  cout<< "-----------------------------------"<<endl;
  cout<< "Set iron thickness     :"<<2* parirlay[2]<<"mm"<<endl;
  cout<< "Set thickness of airGap:"<<2* parlay[2]<<"mm"<<endl;
  cout<< "Set stripwidth   X:Y   :"<<paradef->GetXStrwd()<<"mm"<<" : "<<paradef->GetXStrwd()<<"mm"<<endl;
  cout<< "Set no. of strips X:Y   :"<<paradef->GetnXStrip()<<" : "<<paradef->GetnYStrip()<<endl;
  cout<< "Number of iron layers  :"<<nIRLayer<<endl;
  cout<< "Number of RPC layers  :"<< nLayer<<endl;
  cout<< "Number of RPC stacks  :"<< nStack<<endl;
  cout<< "Number of RoomWallThicknessZ : "<< RoomWallThicknessZ <<endl;
  cout<< "Number of RoomWallThickness : "<< RoomWallThickness <<endl;
  cout<<"ShiftInX = "<< ShiftInX <<endl;
  cout<<"ShiftInY = "<< ShiftInY <<endl;
  cout<< "-----------------------------------"<<endl;

  PrintCalorParameters();

  double xpos = 0;
  double ypos = 0;
  double zpos = 0;

  //     
  // World
  //
  //  solidWorld = 0; logicWorld=0; physiWorld=0;
  solidWorld = new G4Box("WorldSolid",parworld[0],parworld[1],parworld[2]);
  
  logicWorld = new G4LogicalVolume(solidWorld,
				   Air,
				   "WorldLogic",0);
  physiWorld = new G4PVPlacement(0,
				 G4ThreeVector(0,0,0),
				 logicWorld,
				 "WorldPhysi",
				 0,
				 false,
				 0);
  
  // Building
  solidBuilding = new G4Box("BuildingSolid",parroom[0],parroom[1],parroom[2]);
  logicBuilding = new G4LogicalVolume(solidBuilding,
				      ConcreteMaterial,
				      "BuildingLogic",0);
  
  // Air Room
  solidAirRoom = new G4Box("AirRoomSolid",parairroom[0],parairroom[1],parairroom[2]);
  logicAirRoom = new G4LogicalVolume(solidAirRoom,
				     Air,
				     "AirRoomLogic",0);
  
  // Air Room Up stairs
  solidAirRoomUp = new G4Box("AirRoomSolidUp",parairroom[0],parairroom[1],parairroom[2]);
  logicAirRoomUp = new G4LogicalVolume(solidAirRoomUp,
				      Air,
				       "AirRoomUpLogic",0);
  
  solidAirRoom3 = new G4Box("AirRoom3Solid",parstaircaseair[0],parstaircaseair[1],parairroom[2]);
  logicAirRoom3 = new G4LogicalVolume(solidAirRoom3,
				      Air,
				      "AirRoom3Logic",0);
  
  solidAirStairCase = new G4Box("AirStairCaseSolid",parstaircaseair[0],parstaircaseair[1],parairroom[2]);
  logicAirStairCase = new G4LogicalVolume(solidAirStairCase,
					  Air,
					  "AirStairCaseLogic",0);
  
  solidStairCaseL = new G4Box("StairCaseLSolid",parstaircasel[0],parstaircasel[1],parstaircasel[2]);
  logicStairCaseL = new G4LogicalVolume(solidStairCaseL,
					ConcreteMaterial,
					"StairCaseLLogic",0);

  solidStairCase = new G4Box("StairCaseSolid",parstaircase[0],parstaircase[1],parstaircase[2]);
  logicStairCase = new G4LogicalVolume(solidStairCase,
					ConcreteMaterial,
					"StairCaseLogic",0);
  // Air Room
  solidAirRoom2 = new G4Box("AirRoom2Solid",parairroom2[0],parairroom2[1],parairroom2[2]);
  logicAirRoom2 = new G4LogicalVolume(solidAirRoom2,
				      Air,
				      "AirRoom2Logic",0);
  
  // Air Room Up stairs
  solidAirRoom2Up = new G4Box("AirRoom2SolidUp",parairroom2[0],parairroom2[1],parairroom2[2]);
  logicAirRoom2Up = new G4LogicalVolume(solidAirRoom2Up,
					Air,
					"AirRoom2UpLogic",0);
  
  //miniICAL main structure
  solidINOM = new G4Box("INOM", parino[0], parino[1], parino[2]);
  logicINOM = new G4LogicalVolume(solidINOM,
				  Air,
				  "INOMlog",0);

  //Magnet main structure
  solidMagnet = new G4Box("MagnetSolid", parmagnet[0], parmagnet[1], parmagnet[2]);
  logicMagnet = new G4LogicalVolume(solidMagnet,
				  Air,
				  "MagnetLog",0);
  
  G4RotationMatrix Rot;
  Rot.rotateY(0*rad);
  
  // creating space for coil in RPC layer
  for(int lx=0; lx<nLayer; lx++) {
    // cout<<lx<<" "<<IronLayerZdim[lx]<<" "<<IronLayerZdim[lx]/2<<endl;
    G4Box* solidNoSpaceRPCLayer = new G4Box("NoSpaceRPCLayer", parlay[0], parlay[1], LayerZdim[lx]/2);
    G4Box* solidCoilSpaceRPC = new G4Box("CoilSpaceRPC", parcoilspacerpc[0], parcoilspacerpc[1], LayerZdim[lx]/2);
    // RPC Layers in INO structuregunMessenger = new micalPrimaryGeneratorMessenger(this);
    xpos = -parchm[0];
    for(int jk=0;jk<2;jk++) {
      ypos = -parlay[1] + parcoilspacerpc[1] + 800*mm;
      for(int ij=0;ij<nCoil;ij++) {
	G4ThreeVector trans(xpos,ypos,0);
	G4Transform3D transform(Rot,trans);
	if((ij==0)&&(jk==0)) {
	  solidLAYE[lx] = new G4SubtractionSolid("LAYE", solidNoSpaceRPCLayer, solidCoilSpaceRPC,transform); 
	} else {
	  solidLAYE[lx] = new G4SubtractionSolid("LAYE", solidLAYE[lx], solidCoilSpaceRPC,transform); 
	}
	ypos = -ypos;
      }
      xpos = -xpos;
    }
    logicLAYE[lx] = new G4LogicalVolume(solidLAYE[lx],
					Air,
					"LAYElog");
    
    //spacer type A
    solidSpacerA[lx] = new G4Box("SpacerA",parspacerA[0],parspacerA[1],LayerZdim[lx]/2);
    logicSpacerA[lx] = new G4LogicalVolume(solidSpacerA[lx],
				       Iron,
				       "SpacerAlog");
    xpos = -parlay[0] + parspacerA[0];
    ypos = -parlay[1] + parspacerA[1];
    for(int jk=0;jk<nSpacerA/2;jk++) {
      for (int ij=0; ij<nSpacerA/2; ij++) {
	physiSpacerA[lx] = new  G4PVPlacement(0,
					      G4ThreeVector((ij==0)?xpos:-xpos,ypos,0),
					      logicSpacerA[lx],
					      "SpacerAphy",
					      logicLAYE[lx],
					      false,
					      ij+(nSpacerA/2)*jk);
      }
      ypos=-ypos;
    }
    
    //spacer type B
    solidSpacerB[lx] = new G4Box("Spacer2",parspacerB[0],parspacerB[1],LayerZdim[lx]/2);
    logicSpacerB[lx] = new G4LogicalVolume(solidSpacerB[lx],
					   Iron,
					   "SpacerBlog");
    xpos = - parlay[0] + parspacerB[0];
    for (int ij=0; ij<nSpacerB; ij++) {
      physiSpacerB[lx] = new  G4PVPlacement(0,
					    G4ThreeVector(xpos,0,0),
					    logicSpacerB[lx],
					    "SpacerBphy",
					    logicLAYE[lx],
					    false,
					    ij);
      xpos=-xpos;
    }
    
    //spacer type C
    solidSpacerC[lx] = new G4Box("Spacer3",parspacerC[0],parspacerC[1],LayerZdim[lx]/2);
    logicSpacerC[lx] = new G4LogicalVolume(solidSpacerC[lx],
				       Iron,
				       "SpacerClog");
    xpos = -parchm[0];
    ypos = -parlay[1] + parspacerC[1];
    for(int jk=0;jk<nSpacerC/2;jk++) {
      for (int ij=0; ij<nSpacerC/2; ij++) {
	physiSpacerC[lx] = new  G4PVPlacement(0,
					  G4ThreeVector((ij==0)?xpos:-xpos,ypos,0),
					  logicSpacerC[lx],
					  "SpacerCphy",
					  logicLAYE[lx],
					  false,
					  ij+(nSpacerC/2)*jk);
      }
      ypos=-ypos;
    }
    
    //spacer type D
    solidSpacerD[lx] = new G4Box("Spacer4",parspacerD[0],parspacerD[1],LayerZdim[lx]/2);
    logicSpacerD[lx] = new G4LogicalVolume(solidSpacerD[lx],
				       Iron,
				       "SpacerDlog");
    xpos = - parchm[1];
    for (int ij=0; ij<nSpacerD; ij++) {
      physiSpacerD[lx] = new  G4PVPlacement(0,
					G4ThreeVector(xpos,0,0),
					logicSpacerD[lx],
					"SpacerDphy",
					logicLAYE[lx],
					false,
					ij);
      xpos=-xpos;
    }
    
  } // for(int lx=0; lx<nLayer; lx++) {
  
  //creating space for coil in iron layer
  for(int irlx=0; irlx<nIRLayer; irlx++) {
    // cout<<irlx<<" "<<IronLayerZdim[irlx]<<" "<<IronLayerZdim[irlx]/2<<endl;
    G4Box* solidNoSpaceIRONLayer = new G4Box("NoSpaceIRONlayer", parirlay[0], parirlay[1], IronLayerZdim[irlx]/2);
    G4Box* solidCoilSpaceIRON = new G4Box("CoilSpaceRPC", parcoilspaceiron[0], parcoilspaceiron[1], IronLayerZdim[irlx]/2);
    G4SubtractionSolid* solidNoGapIRONLayer = 0;
    xpos = -parchm[0];
    for(int jk=0;jk<2;jk++) {
      ypos = -parlay[1] + parcoilspacerpc[1] + 800*mm;
      for(int ij=0;ij<nCoil;ij++) {
	G4ThreeVector  trans(xpos,ypos,0);
	G4Transform3D transform(Rot,trans);
	if((ij==0)&&(jk==0)) {
	  solidNoGapIRONLayer = new G4SubtractionSolid("IRLAYE", solidNoSpaceIRONLayer, solidCoilSpaceIRON,transform);
	} else {
	  solidNoGapIRONLayer = new G4SubtractionSolid("IRLAYE", solidNoGapIRONLayer, solidCoilSpaceIRON,transform);
	}
	ypos = -ypos;
      }
      xpos = -xpos;
    }
    
    //creating space for horizontal air gap in iron layer
    G4Box* solidAirGap1 = new G4Box("AirGap1", parairgap1[0], parairgap1[1], IronLayerZdim[irlx]/2);
    G4SubtractionSolid* solidVGapIRONLayer1 = 0;
    ypos = 0.0*cm;
    xpos = -parchm[0];
    for(int ij=0;ij<2;ij++) {
      G4ThreeVector  trans(xpos,ypos,0);
      G4Transform3D transform(Rot,trans);
      if(ij==0) {
	solidVGapIRONLayer1 = new G4SubtractionSolid("IRLAYE", solidNoGapIRONLayer, solidAirGap1,transform);
	// cout<<"solidVGapIRONLayer1 (ij,xpos,ypos) = ("<<ij<<","<<xpos<<","<<ypos<<")"<<endl;
      } else {
	solidVGapIRONLayer1 = new G4SubtractionSolid("IRLAYE", solidVGapIRONLayer1, solidAirGap1,transform);
	// cout<<"solidVGapIRONLayer1 (ij,xpos,ypos) = ("<<ij<<","<<xpos<<","<<ypos<<")"<<endl;
      }
      xpos = -xpos;
    }
  
    //creating space for vertical air gap in iron layer
    G4Box* solidAirGap2 = new G4Box("AirGap2", parairgap2[0], parairgap2[1], IronLayerZdim[irlx]/2);
    G4SubtractionSolid* solidVGapIRONLayer2 = 0;
    xpos = -parchm[0];
    for(int jk=0;jk<2;jk++) {
      ypos = -parlay[1] + parairgap2[1];
      for(int ij=0;ij<2;ij++) {
	G4ThreeVector trans(xpos,ypos,0);
	G4Transform3D transform(Rot,trans);
	if(ij==0 && jk==0) {
	  solidVGapIRONLayer2 = new G4SubtractionSolid("IRLAYE", solidVGapIRONLayer1, solidAirGap2,transform);
	  // cout<<"solidVGapIRONLayer2 (jk,ij,xpos,ypos) = ("<<jk<<","<<ij<<","<<xpos<<","<<ypos<<")"<<endl;
	} else { 
	  solidVGapIRONLayer2 = new G4SubtractionSolid("IRLAYE", solidVGapIRONLayer2, solidAirGap2,transform);
	  // cout<<"solidVGapIRONLayer2 (jk,ij,xpos,ypos) = ("<<jk<<","<<ij<<","<<xpos<<","<<ypos<<")"<<endl;
	}
	ypos = -ypos;
      }
      xpos = -xpos;
    }
  
    G4Box* solidAirGap3 = new G4Box("AirGap3", parairgap3[0], parairgap3[1], IronLayerZdim[irlx]/2);
    G4SubtractionSolid* solidHGapIRONLayer1 = 0;
    ypos = 0;
    xpos = parirlay[0] - parairgap3[0] - 1*mm;
    for(int ij=0;ij<2;ij++) {
      G4ThreeVector trans(xpos,ypos,0);
      G4Transform3D transform(Rot,trans);
      if(ij==0) {
	solidHGapIRONLayer1 = new G4SubtractionSolid("IRLAYE", solidVGapIRONLayer2, solidAirGap3,transform); 
	// cout<<"solidHGapIRONLayer1 (ij,xpos,ypos) = ("<<ij<<","<<xpos<<","<<ypos<<")"<<endl;
      } else {
	solidHGapIRONLayer1 = new G4SubtractionSolid("IRLAYE", solidHGapIRONLayer1, solidAirGap3,transform);
	// cout<<"solidHGapIRONLayer1 (ij,xpos,ypos) = ("<<ij<<","<<xpos<<","<<ypos<<")"<<endl;
      }
      xpos = -xpos;
    }

    G4Box* solidAirGap4 = new G4Box("AirGap4", parairgap4[0], parairgap4[1], IronLayerZdim[irlx]/2);
    xpos = 0;
    ypos = -803*mm;
    for(int ij=0;ij<2;ij++) {
      G4ThreeVector trans(xpos,ypos,0);
      G4Transform3D transform(Rot,trans);
      if(ij==0) {
	solidIRLAYE[irlx] = new G4SubtractionSolid("IRLAYE", solidHGapIRONLayer1, solidAirGap4,transform); 
	// cout<<"solidIRLAYE (ij,xpos,ypos) = ("<<","<<ij<<","<<xpos<<","<<ypos<<")"<<endl;
      } else {
	solidIRLAYE[irlx] = new G4SubtractionSolid("IRLAYE", solidIRLAYE[irlx], solidAirGap4,transform);
	// cout<<"solidIRLAYE (ij,xpos,ypos) = ("<<","<<ij<<","<<xpos<<","<<ypos<<")"<<endl;
      }
      ypos = -ypos;
    }
    
    logicIRLAYE[irlx] = new G4LogicalVolume(solidIRLAYE[irlx],
					    Iron,
					    // Air,
					    "IRLAYElog");
    
    G4bool allLocal = true ;
    logicIRLAYE[irlx]->SetFieldManager( fieldMgr,allLocal );
    
  } // for(int irlx=0; irlx<nIRLayer; irlx++) {
    
  solidFRPBox = new G4Box("FRPBox",parfrpbox[0],parfrpbox[1],parfrpbox[2]);
  logicFRPBox = new G4LogicalVolume(solidFRPBox,
				    CarbonFRP,
				    "FRPBoxlog");

  solidAirBox = new G4Box("AirBox",parairbox[0],parairbox[1],parairbox[2]);
  logicAirBox = new G4LogicalVolume(solidAirBox,
				    Air,
				    "AIRBOXlog");

  // Triangle Electronics Board for DFE
  solidG10Trap1 = new G4Trd("G10Trap1",parg10Trap1[1],parg10Trap1[2],parg10Trap1[3],parg10Trap1[3],parg10Trap1[4]);
  logicG10Trap1 = new G4LogicalVolume(solidG10Trap1,
				      G10,
				      "G10Trap1log");
  
  solidG10Trap2 = new G4Trd("G10Trap2",parg10Trap2[1],parg10Trap2[2],parg10Trap2[3],parg10Trap2[3],parg10Trap2[4]);
  logicG10Trap2 = new G4LogicalVolume(solidG10Trap2,
				      G10,
				      "G10Trap2log");
  
  // Aluminium
  solidAL = ConstructRPCBox(paral,paralCutBig,paralCutSmall);
  solidAL->SetName("AL");
  logicAL = new G4LogicalVolume(solidAL,
 				Aluminium,
   				"ALlog");
  
  // HoneyComb
  solidHoneyComb = ConstructRPCBox(parhoneycomb,parhoneycombCutBig,parhoneycombCutSmall);
  solidHoneyComb->SetName("HneyComb");
  logicHoneyComb = new G4LogicalVolume(solidHoneyComb,
   				       HoneyCombMaterial,
   				       "HoneyComblog");

  // Copper
  solidCUPL = ConstructRPCBox(parcup,parcupCutBig,parcupCutSmall);
  solidCUPL->SetName("CUPL");
  logicCUPL = new G4LogicalVolume(solidCUPL,
   				  Copper,
   				  "CUPLlog");
  
  // Mylar
  solidMYLAR = ConstructRPCBox(parmylar,parmylarCutBig,parmylarCutSmall);
  solidMYLAR->SetName("MYLAR");
  logicMYLAR = new G4LogicalVolume(solidMYLAR,
   				   MylarMaterial,
				   "MYLARlog");

  // Graphite Coating
  solidCOAT = ConstructRPCBox(parcoat,parcoatCutBig,parcoatCutSmall);
  solidCOAT->SetName("COAT");
  logicCOAT = new G4LogicalVolume(solidCOAT,
   				  CoatMaterial,
   				  "COATlog");
  
  // Quartz
  solidQURZ = ConstructRPCBox(parqurz,parqurzCutBig,parqurzCutSmall);
  solidQURZ->SetName("QURZ");
  logicQURZ = new G4LogicalVolume(solidQURZ,
   				  SiO2,
   				  "QURZlog");

  // RPC GAS
  solidGASR = ConstructRPCBox(pargas,pargasCutBig,pargasCutSmall);
  solidGASR->SetName("GASR");
  logicGASR = new G4LogicalVolume(solidGASR,
   				  ActiveMaterial,
   				  "GASRlog");
  
  // coil in vertical direction
  solidVCOIL = new G4Box("VCOIL", parvcoil[0], parvcoil[1], parvcoil[2]); 
  logicVCOIL = new G4LogicalVolume(solidVCOIL,
				   Copper,//Iron,
				   "VCOILlog");
  
  // coil in horizontal direction
  solidHCOIL = new G4Box("HCOIL", parhcoil[0], parhcoil[1], parhcoil[2]); 
  logicHCOIL = new G4LogicalVolume(solidHCOIL,
				   Copper,//Iron,
				   "HCOILlog");

  // curved portion of coil
  solidCurvedCOIL = new G4Tubs("CurvedCOIL", parcurvedcoil[0], parcurvedcoil[1], parcurvedcoil[2], parcurvedcoil[3],  parcurvedcoil[4] ); 
  logicCurvedCOIL = new G4LogicalVolume(solidCurvedCOIL,
					Copper,//Iron,
					"CurvedCOILlog");
  
  //coil support
  solidCOILSupport = new G4Box("COILSupport", parcoilsupport[0], parcoilsupport[1], parcoilsupport[2]); 
  logicCOILSupport = new G4LogicalVolume(solidCOILSupport,
					 Iron,
					 "COILSupportlog");
  
  logicGASR->SetSensitiveDetector(cal0SD);
  logicGASR->SetRegion(aRegion0);
  aRegion0->AddRootLogicalVolume(logicGASR);
  
  // RPC gas inside RPC glass(Quartz)
  xpos = ypos = 0.0*cm;
  physiGASR = new  G4PVPlacement(0,
				 G4ThreeVector(xpos,ypos,0),
				 logicGASR,
				 "GASRphy",
				 logicQURZ,
				 false,
				 0);
  
  //Quartz inside Graphite
  physiQURZ = new  G4PVPlacement(0,
				 G4ThreeVector(0,0,0),
				 logicQURZ,
				 "QURZphy",
				 logicCOAT,
				 false,
				 0);
  
  // Positioning graphite inside mylar
  physiCOAT = new  G4PVPlacement(0,
				 G4ThreeVector(0,0,0),
				 logicCOAT,
				 "COATphy",
				 logicMYLAR,
				 false,
				 0);
  
  // Positioning mylar inside copper
  physiMYLAR = new  G4PVPlacement(0,
				  G4ThreeVector(0,0,0),
				  logicMYLAR,
				  "MYLARphy",
				  logicCUPL,
				  false,
				  0);
  
  // Copper inside honeycomb
  physiCUPL = new  G4PVPlacement(0,
				 G4ThreeVector(0,0,0),
				 logicCUPL,
				 "CUPLphy",
				 logicHoneyComb,
				 false,
				 0);
  
  
  // honeycomb inside Al
  physiHoneyComb = new  G4PVPlacement(0,
				      G4ThreeVector(0,0,0),
				      logicHoneyComb,
				      "HoneyCombphy",
				      logicAL,
				      false,
				      0);
  
  // cout<<"xpos = "<< ShiftInX <<endl;
  // cout<<"ypos = "<< ShiftInY <<endl;

  xpos = ShiftInX;
  ypos = ShiftInY;
  zpos = AlShiftInAirBox; //-(parairbox[2]-paral[2]);
  // cout<<"al zpos "<<zpos<<endl;
  // Al inside g10
  physiAL = new  G4PVPlacement(0,
			       G4ThreeVector(xpos,ypos,zpos),
			       logicAL,
			       "ALhy",
			       logicAirBox,
			       false,
			       0);
  
  //G10 trap position in logicAirBox
  G4RotationMatrix rotmG10Trap1  = G4RotationMatrix();
  rotmG10Trap1.rotateY(45*deg);
  rotmG10Trap1.rotateX(-90*deg);
  xpos = parg10Trap1[5];
  ypos = parg10Trap1[6];
  zpos = parg10Trap1[7];
  G4ThreeVector posG10Trap1 = G4ThreeVector(xpos,ypos,zpos);
  G4Transform3D transformG10Trap1 = G4Transform3D(rotmG10Trap1,posG10Trap1);
  physiG10Trap1 = new G4PVPlacement(transformG10Trap1,
  				    logicG10Trap1,
  				    "G10Trap1phy",
  				    logicAirBox,
  				    false,
  				    0); 
  
  //G10 trap position in logicAirBox
  G4RotationMatrix rotmG10Trap2 = G4RotationMatrix();
  rotmG10Trap2.rotateY(-45*deg);
  rotmG10Trap2.rotateX(90*deg);
  xpos = parg10Trap2[5];
  ypos = parg10Trap2[6];
  zpos = parg10Trap2[7];
  G4ThreeVector posG10Trap2 = G4ThreeVector(xpos,ypos,zpos);
  G4Transform3D transformG10Trap2 = G4Transform3D(rotmG10Trap2,posG10Trap2);

  physiG10Trap2 = new G4PVPlacement(transformG10Trap2,
				    logicG10Trap2,
				    "G10Trap2phy",
				    logicAirBox,
				    false,
				    0); 

  // cout<<"zpos = "<< ShiftInZ <<endl;
  xpos = 0;
  ypos = 0;
  zpos = AirBoxShiftInFRP; //ShiftInZ;

  physiAirBox = new  G4PVPlacement(0,
				   G4ThreeVector(xpos,ypos,zpos),
				   logicAirBox,
				   "AIRBOXphy",
				   logicFRPBox,
				   false,
				   0);
  
  // Positioning Air Box(containing RPC + G10 Trapezoid) in layer
  if(nStack>2) {
    cout<<"nStacks > 2 shouldn't happen."<<endl;
    nChamber = 2;
  }
  
  for(int lx=0; lx<nLayer; lx++) {
    zpos = FRPshiftInLayer[lx]; //-((LayerZdim[lx]/2)-parfrpbox[2]);
    // cout<<"lx "<<lx<<" "<<zpos<<endl;
    for (int ij=0; ij<nStack; ij++) {
      if(nStack==1) {
	ypos = 0;
      } else {
	ypos = (ij) ? -1*parchm[1] : parchm[1];
      }
      physiFRPBox = new  G4PVPlacement(0,
				       G4ThreeVector(0,ypos,zpos),
				       logicFRPBox,
				       "FRPBOXphy",
				       logicLAYE[lx],
				       false,
				       ij);
    }
  }
  
  //iron Layers in INO structure
  for (int ij=0; ij<nIRLayer; ij++) {
    zpos = IRONLayerPosZ[ij];
    cout<<"IRON "<<ij<<" "<<zpos<<endl;
    //-( (int(0.5*nIRLayer)-ij)*2*(parirlay[2]+ parlay[2]) );
    physiIRLAYE[ij] = new  G4PVPlacement(0,
					 G4ThreeVector(0,0,zpos),
					 logicIRLAYE[ij],
					 "IRLAYEphy",
					 logicMagnet,
					 false,
					 ij);
  }
  
  //RPC Layers in INO structure
  for (int ij=0; ij<nLayer; ij++) {
    zpos = RPCLayerPosZ[ij];
    cout<<"RPC "<<ij<<" "<<zpos<<endl;
    //-((int(0.5*nIRLayer)-ij)*2*(parirlay[2]+ parlay[2]))+(parirlay[2]+ parlay[2]);
    physiLAYE[ij] = new  G4PVPlacement(0,
				       G4ThreeVector(0,0,zpos),
				       logicLAYE[ij],
				       "LAYEphy",
				       logicMagnet,
				       false,
				       ij);
  }
  
  //coil in vertical direction in INO module
  xpos = -parchm[0];
  for(int jk=0;jk<2;jk++) {
    ypos = -parlay[1] + parcoilspacerpc[1] + 800*mm;
    for(int ij=0;ij<nCoil;ij++) {
      physiVCOIL = new  G4PVPlacement(0,
				      G4ThreeVector(xpos,ypos,0),
				      logicVCOIL,
 				      "VCOILphy",
				      logicMagnet,
				      false,
				      ij+nCoil*jk);
      ypos = -ypos;
    }
    xpos = -xpos;
  }

  // coil in horizontal direction
  // zpos = -parino[2] + parhcoil[2];
  zpos = -parvcoil[2] - (parcurvedcoil[0] + parcurvedcoil[1])/2;
  postopcoil = -(-parvcoil[2] - (parcurvedcoil[0] + parcurvedcoil[1])/2);//-zpos;
  for(int jk=0;jk<2;jk++) {
    ypos = -parlay[1] + parcoilspacerpc[1] + 800*mm;
    for(int ij=0;ij<nCoil;ij++) {
      // cout<<"jk,ij "<<jk<<","<<ij<<" "<<ypos<<" "<<zpos<<" "<<endl;
      physiHCOIL = new  G4PVPlacement(0,
  				      G4ThreeVector(0,ypos,zpos),
  				      logicHCOIL,
  				      "HCOILphy",
  				      logicMagnet,
  				      false,
  				      ij+nCoil*jk);
      ypos = -ypos;
    }
    zpos = -zpos;
  }
  
  //curved portion of coil at 4 corners
  G4RotationMatrix xRot,yRot,zRot,totalRot ;
  G4double theta=M_PI/2;
  yRot.rotateY(2*theta);
  xRot.rotateX(-theta);
  zRot.rotateZ(0*theta);
  totalRot=xRot*yRot*zRot;
  
  xpos = -parchm[0] + (parcurvedcoil[0] + parcurvedcoil[1])/2;
  zpos = -parvcoil[2];
  for(int kl=0;kl<2;kl++) {
    for(int jk=0;jk<2;jk++) {
      ypos =  -parlay[1] + parcoilspacerpc[1] + 800*mm;
      for(int ij=0;ij<nCoil;ij++) {
	physiCurvedCOIL = new  G4PVPlacement(G4Transform3D(totalRot,G4ThreeVector(xpos,ypos,zpos)),
					     logicCurvedCOIL,
					     "CurvedCOILphy",
					     logicMagnet,
					     false,
					     ij+nCoil*jk+2*nCoil*kl);
	ypos = -ypos;
      }
      xpos = -xpos;
      if(kl==0)
	yRot.rotateY(-2*theta);
      else
	yRot.rotateY(2*theta);
      totalRot=xRot*yRot*zRot;
    }
    zpos = -zpos;
    xRot.rotateX(-2*theta);
    totalRot=xRot*yRot*zRot;
  }
  
  // //coil support
  // zpos = -parino[2] + 2*parhcoil[2] + parcoilsupport[2];
  // for(int kl=0;kl<2;kl++) {
  //   for(int jk=0;jk<nCoil;jk++) {
  //     ypos =  -parlay[1]+(2*jk+5)*parchm[1];
  //     for(int ij=0;ij<nCoilSupport;ij++) {
  // 	xpos =  -parlay[0] + (3.5*ij + 4.5)*parchm[0];
  // 	physiCOILSupport = new  G4PVPlacement(0,		  //no rotation
  // 					      G4ThreeVector(xpos,ypos,zpos), //position
  // 					      logicCOILSupport,	   //its logical volume
  // 					      "COILSupportphy",		//its name
  // 					      logicINOM,	   //its mother  volume
  // 					      false,		//no boolean operation
  // 					      ij+nCoilSupport*jk+nCoil*nCoilSupport*kl);	                                                                 //copy number
  //     }
  //   }
  //   zpos = -zpos;
  // }

  xpos = 0.0;
  ypos = 0.0;
  zpos = -parino[2] + parmagnet[2];
  posmagnet = zpos;
  physiMagnet = new  G4PVPlacement(0,
				 G4ThreeVector(xpos,ypos,zpos),
				 logicMagnet,
				 "MagnetPhy",
				 logicINOM,
				 false,
				 0);
  
  // xpos = 0.0;
  // ypos = 0.0;
  // zpos = posmagnet + parmagnet[2] + 180 + partopscint[2];
  // char namexs[100];
  // char namexl[100];
  // char namexp[100];
  // for(int nsc=0; nsc<nScintLayer; nsc++) {
  //   sprintf(namexs,"TopScintSolid_%i",nsc);
  //   sprintf(namexl,"TopScintLogic_%i",nsc);
  //   sprintf(namexp,"TopScintPhysi_%i",nsc);
  //   solidTopScint[nsc] = new G4Box(namexs,partopscint[0],partopscint[1],partopscint[2]);
  //   logicTopScint[nsc] = new G4LogicalVolume(solidTopScint[nsc],
  // 					     Air,
  // 					     namexl);
  //   physiTopScint[nsc] = new  G4PVPlacement(0,
  // 					    G4ThreeVector(xpos,ypos,zpos),
  // 					    logicTopScint[nsc],
  // 					    namexp,
  // 					    logicINOM,
  // 					    false,
  // 					    nsc);
  //   zpos += AirGapScintTop + ScintUnitZ;
  // }
  
  // xpos = 0.0;
  // ypos = -partopscint[1] - 1*mm - parwallscint[2];
  // zpos = -parino[2] + ScintFromBottom + parwallscint[0];
  // for(int nsc=0; nsc<nScintLayer; nsc++) {
  //   sprintf(namexs,"WallScintFrontSolid_%i",nsc);
  //   sprintf(namexl,"WallScintFrontLogic_%i",nsc);
  //   sprintf(namexp,"WallScintFrontPhysi_%i",nsc);
  //   solidWallScintFront[nsc] = new G4Box(namexs,parwallscint[1],parwallscint[2],parwallscint[0]);
  //   logicWallScintFront[nsc] = new G4LogicalVolume(solidWallScintFront[nsc],
  // 					     Air,
  // 					     namexl);
  //   physiWallScintFront[nsc] = new  G4PVPlacement(0,
  // 					    G4ThreeVector(xpos,ypos,zpos),
  // 					    logicWallScintFront[nsc],
  // 					    namexp,
  // 					    logicINOM,
  // 					    false,
  // 					    nsc);
  //   ypos += -AirGapScintWall - ScintUnitZ;
  // }
  
  // xpos = 0.0;
  // ypos = partopscint[1] + 1*mm + parwallscint[2];
  // zpos = -parino[2] + ScintFromBottom + parwallscint[0];
  // for(int nsc=0; nsc<nScintLayer; nsc++) {
  //   sprintf(namexs,"WallScintBackSolid_%i",nsc);
  //   sprintf(namexl,"WallScintBackLogic_%i",nsc);
  //   sprintf(namexp,"WallScintBackPhysi_%i",nsc);
  //   solidWallScintBack[nsc] = new G4Box(namexs,parwallscint[1],parwallscint[2],parwallscint[0]);
  //   logicWallScintBack[nsc] = new G4LogicalVolume(solidWallScintBack[nsc],
  // 					     Air,
  // 					     namexl);
  //   physiWallScintBack[nsc] = new  G4PVPlacement(0,
  // 					    G4ThreeVector(xpos,ypos,zpos),
  // 					    logicWallScintBack[nsc],
  // 					    namexp,
  // 					    logicINOM,
  // 					    false,
  // 					    nsc);
  //   ypos += AirGapScintWall + ScintUnitZ;
  // }

  // xpos = partopscint[0] + 1*mm + parwallscint[2];
  // ypos = 0.0;
  // zpos = -parino[2] + ScintFromBottom + parwallscint[0];
  // for(int nsc=0; nsc<nScintLayer; nsc++) {
  //   sprintf(namexs,"WallScintRightSolid_%i",nsc);
  //   sprintf(namexl,"WallScintRightLogic_%i",nsc);
  //   sprintf(namexp,"WallScintRightPhysi_%i",nsc);
  //   solidWallScintRight[nsc] = new G4Box(namexs,parwallscint[2],parwallscint[1],parwallscint[0]);
  //   logicWallScintRight[nsc] = new G4LogicalVolume(solidWallScintRight[nsc],
  // 					     Air,
  // 					     namexl);
  //   physiWallScintRight[nsc] = new  G4PVPlacement(0,
  // 					    G4ThreeVector(xpos,ypos,zpos),
  // 					    logicWallScintRight[nsc],
  // 					    namexp,
  // 					    logicINOM,
  // 					    false,
  // 					    nsc);
  //   xpos += AirGapScintWall + ScintUnitZ;
  // }

  // xpos = -partopscint[0] - 1*mm - parwallscint[2];
  // ypos = 0.0;
  // zpos = -parino[2] + ScintFromBottom + parwallscint[0];
  // for(int nsc=0; nsc<nScintLayer; nsc++) {
  //   sprintf(namexs,"WallScintLeftSolid_%i",nsc);
  //   sprintf(namexl,"WallScintLeftLogic_%i",nsc);
  //   sprintf(namexp,"WallScintLeftPhysi_%i",nsc);
  //   solidWallScintLeft[nsc] = new G4Box(namexs,parwallscint[2],parwallscint[1],parwallscint[0]);
  //   logicWallScintLeft[nsc] = new G4LogicalVolume(solidWallScintLeft[nsc],
  // 					     Air,
  // 					     namexl);
  //   physiWallScintLeft[nsc] = new  G4PVPlacement(0,
  // 					    G4ThreeVector(xpos,ypos,zpos),
  // 					    logicWallScintLeft[nsc],
  // 					    namexp,
  // 					    logicINOM,
  // 					    false,
  // 					    nsc);
  //   xpos += -AirGapScintWall - ScintUnitZ;
  // }

  xpos = paradef->GetStackPosInRoom(0); 
  ypos = paradef->GetStackPosInRoom(1); 
  zpos = paradef->GetStackPosInRoom(2); 
  physiINOM = new  G4PVPlacement(0,
				 G4ThreeVector(xpos,ypos,zpos),
				 logicINOM,
				 "INOMphy",
				 logicAirRoom,
				 false,
				 0);
  
  xpos = paradef->GetINOroomPos(0); 
  ypos = paradef->GetINOroomPos(1); 
  zpos = paradef->GetINOroomPos(2); 
  physiAirRoom = new G4PVPlacement(0,
				   G4ThreeVector(xpos,ypos,zpos),
				   logicAirRoom,
				   "AirRoom",
				   logicBuilding,
				   false,
				   0);

  xpos = paradef->GetINOroomPos(0); 
  ypos = paradef->GetINOroomPos(1); 
  zpos = -paradef->GetINOroomPos(2); 
  physiAirRoomUp = new G4PVPlacement(0,
  				     G4ThreeVector(xpos,ypos,zpos),
  				     logicAirRoomUp,
  				     "AirRoomUp",
  				     logicBuilding,
  				     false,
  				     0);
  
  ypos = -parroom[1] + RoomWallThickness + parstaircaseair[1];
  xpos = 0.0*mm;
  zpos = paradef->GetINOroomPos(2);
  physiAirRoom3 = new G4PVPlacement(0,
  				    G4ThreeVector(xpos,ypos,zpos),
  				   logicAirRoom3,
  				    "AirRoom3_1",
  				    logicBuilding,
  				    false,
  				    0);
  
  ypos = -parroom[1] + RoomWallThickness + parstaircaseair[1];
  xpos = 0.0*mm;
  zpos = -paradef->GetINOroomPos(2);
  physiAirRoom3 = new G4PVPlacement(0,
  				    G4ThreeVector(xpos,ypos,zpos),
  				    logicAirRoom3,
  				    "AirRoom3_2",
  				    logicBuilding,
  				    false,
  				    1);
  
  ypos = parroom[1] - RoomWallThickness - parstaircaseair[1];
  xpos = 0.0*mm;
  zpos = paradef->GetINOroomPos(2);
  physiAirRoom3 = new G4PVPlacement(0,
  				    G4ThreeVector(xpos,ypos,zpos),
  				    logicAirRoom3,
  				    "AirRoom3_3",
  				    logicBuilding,
  				    false,
  				    2);
  
  ypos = parroom[1] - RoomWallThickness - parstaircaseair[1];
  xpos = 0.0*mm;
  zpos = -paradef->GetINOroomPos(2);
  physiAirRoom3 = new G4PVPlacement(0,
  				    G4ThreeVector(xpos,ypos,zpos),
  				    logicAirRoom3,
  				    "AirRoom3_4",
  				    logicBuilding,
  				    false,
  				    3);
  
  xpos = -parstaircaseair[0] + 1*mm + parstaircasel[0];
  ypos = 0.0*mm;
  zpos = 0.0*mm;
  physiStairCaseL = new G4PVPlacement(0,
				      G4ThreeVector(xpos,ypos,zpos),
				      logicStairCaseL,
				      "StairCaseL",
				      logicAirStairCase,
				      false,
				      0);
  // cout<<"staircaseL "<<xpos<<" "<<ypos<<" "<<zpos<<endl;
    
  xpos = -parstaircaseair[0] + 2*parstaircasel[0] + parstaircase[0] - 10*cm;
  ypos = -parstaircaseair[1] + 1*mm + parstaircase[1];//0.0*mm;
  zpos = -parstaircaseair[2]/2 + 1*mm;// + parstaircase[0]*sin(45*deg);//0.0*mm;
  G4RotationMatrix rotmstrcs1;
  rotmstrcs1.rotateY(22*deg);
  physiStairCase = new G4PVPlacement(
				      G4Transform3D(rotmstrcs1,G4ThreeVector(xpos,ypos,zpos)),
				      logicStairCase,
				      "StairCase_1",
				      logicAirStairCase,
				      false,
				      0);
  // cout<<"staircase1 "<<xpos<<" "<<ypos<<" "<<zpos<<endl;

  xpos = -parstaircaseair[0] + 2*mm + 2*parstaircasel[0] + parstaircase[0] - 10*cm;
  ypos = parstaircaseair[1] - 1*mm - parstaircase[1];//0.0*mm;
  zpos = parstaircaseair[2]/2 - 1*mm;//0*mm;//-parstaircaseair[2]/2 + 1*mm;//0.0*mm;
  G4RotationMatrix rotmstrcs2;
  rotmstrcs2.rotateY(-22*deg);
  physiStairCase = new G4PVPlacement(
				      G4Transform3D(rotmstrcs2,G4ThreeVector(xpos,ypos,zpos)),
				      logicStairCase,
				      "StairCase_2",
				      logicAirStairCase,
				      false,
				      1);
  // cout<<"staircase2 "<<xpos<<" "<<ypos<<" "<<zpos<<endl;

  ypos = -parroom[1] + 3*RoomWallThickness + 3*parstaircaseair[1] + 2*parairroom[1];
  xpos = 0.0*mm;
  zpos = paradef->GetINOroomPos(2);;//0.0*mm;//paradef->GetINOroomPos(2);
  physiAirStairCase = new G4PVPlacement(0,
  					G4ThreeVector(xpos,ypos,zpos),
  					logicAirStairCase,
  					"AirStairCase",
  					logicBuilding,
  					false,
  					0);
  
  ypos = -parroom[1] + 3*RoomWallThickness + 3*parstaircaseair[1] + 2*parairroom[1];
  xpos = 0.0*mm;
  zpos = -paradef->GetINOroomPos(2);;//0.0*mm;//paradef->GetINOroomPos(2);
  physiAirStairCase = new G4PVPlacement(0,
  					G4ThreeVector(xpos,ypos,zpos),
  					logicAirStairCase,
  					"AirStairCase",
  					logicBuilding,
  					false,
  					1);
  
  ypos = parroom[1] - 2*RoomWallThickness - 2*parstaircaseair[1] - parairroom2[1];
  xpos = 0.0*mm;
  zpos = paradef->GetINOroomPos(2);
  physiAirRoom2 = new G4PVPlacement(0,
  				   G4ThreeVector(xpos,ypos,zpos),
  				    logicAirRoom2,
  				   "AirRoom2",
  				   logicBuilding,
  				   false,
  				   0);
  // cout<<"xpos "<<xpos<<" "<<ypos<<" "<<zpos<<endl;

  ypos = parroom[1] - 2*RoomWallThickness - 2*parstaircaseair[1] - parairroom2[1];
  xpos = 0.0*mm;
  zpos = -paradef->GetINOroomPos(2);
  physiAirRoom2Up = new G4PVPlacement(0,
  				   G4ThreeVector(xpos,ypos,zpos),
  				    logicAirRoom2Up,
  				   "AirRoom2Up",
  				   logicBuilding,
  				   false,
  				   1);

  // cout<<"xpos "<<xpos<<" "<<ypos<<" "<<zpos<<endl;
  physiBuilding = new G4PVPlacement(0,
				    G4ThreeVector(0,0,0),
				    logicBuilding,
				    "Building",
				    logicWorld,
				    false,
				    0);

  // G4GDMLParser parser;
  // parser.Write("mical_world.gdml",logicWorld);
  // system("rm geo.gdml");
  // system("cp detector_world.gdml geo.gdml");
  // system("rm detector_world.gdml");

  
  // Visualization attributes
  if(!visWhite){visWhite=new G4VisAttributes(true,G4Colour(1.0,1.0,1.0));}
  if(!visYellow){ visYellow=new G4VisAttributes(true,G4Colour(1.0,1.0,0.0));}
  if(!visGray){visGray=new G4VisAttributes(true,G4Colour(0.5,0.5,0.5));}
  if(!visCyan){visCyan=new G4VisAttributes(true,G4Colour(0.0,1.0,1.0));}
  if(!visBlue){visBlue=new G4VisAttributes(true, G4Colour(0.0,0.0,1.0));}
  if(!visRed){visRed=new G4VisAttributes(true,G4Colour(1.0,0.0,0.0));}
  if(!visGreen){visGreen=new G4VisAttributes(true,G4Colour(0.0,1.0,0.0));}
  if(!visMagenta){visMagenta=new G4VisAttributes(true,G4Colour(1.0,0.0,1.0));visMagenta->SetForceSolid(true);}
  if(!visNull){visNull=new G4VisAttributes(G4VisAttributes::Invisible);}
  
  logicWorld->SetVisAttributes(visNull);
  logicBuilding->SetVisAttributes(visNull);
  logicAirRoom->SetVisAttributes(visRed);
  logicAirRoom2->SetVisAttributes(visNull);
  logicAirRoom3->SetVisAttributes(visNull);
  logicAirRoomUp->SetVisAttributes(visNull);
  logicAirRoom2Up->SetVisAttributes(visNull);
  logicINOM->SetVisAttributes(visWhite);
  logicMagnet->SetVisAttributes(visCyan);

  logicAirStairCase->SetVisAttributes(visNull);
  logicStairCaseL->SetVisAttributes(visNull);
  logicStairCase->SetVisAttributes(visNull);

  for(int irlx=0; irlx<nIRLayer; irlx++) {
    logicIRLAYE[irlx]->SetVisAttributes(visBlue);
  }
  for(int lx=0; lx<nLayer; lx++) {
    logicLAYE[lx]->SetVisAttributes(visNull);
    logicSpacerA[lx]->SetVisAttributes(visNull);
    logicSpacerB[lx]->SetVisAttributes(visNull);
    logicSpacerC[lx]->SetVisAttributes(visNull);
    logicSpacerD[lx]->SetVisAttributes(visNull);
  }
  
  logicFRPBox->SetVisAttributes(visNull);
  logicG10Trap1->SetVisAttributes(visNull);
  logicG10Trap2->SetVisAttributes(visNull);
  logicAirBox->SetVisAttributes(visNull);
  logicAL->SetVisAttributes(visNull);
  logicHoneyComb->SetVisAttributes(visNull);
  logicCUPL->SetVisAttributes(visNull);
  logicMYLAR->SetVisAttributes(visNull);
  logicCOAT->SetVisAttributes(visNull);
  logicQURZ->SetVisAttributes(visNull);
  logicGASR->SetVisAttributes(visNull);
  // for(int nsc=0; nsc<nScintLayer; nsc++) {
  //   logicTopScint[nsc]->SetVisAttributes(visGreen);
  //   logicWallScintFront[nsc]->SetVisAttributes(visNull);
  //   logicWallScintBack[nsc]->SetVisAttributes(visNull);
  //   logicWallScintRight[nsc]->SetVisAttributes(visNull);
  //   logicWallScintLeft[nsc]->SetVisAttributes(visNull);
  // }
  logicVCOIL->SetVisAttributes(visNull);
  logicHCOIL->SetVisAttributes(visNull);
  logicCurvedCOIL->SetVisAttributes(visNull);
  logicCOILSupport->SetVisAttributes(visNull);
  
  cout<<" ConstructCalorimeter end  "<<endl;
  system("free");
  
  return physiWorld;
}

void micalDetectorConstruction::PrintCalorParameters() {
  cout << "------------------------------------------------------------"<<endl;
  cout <<"Detector Parameters"<<endl; 
  cout <<"parworld "<<parworld[0]<<"*mm, "<<parworld[1]<<"*mm, "<<parworld[2]<<"*mm"<<endl;
  cout <<"parroom "<<parroom[0]<<"*mm, "<<parroom[1]<<"*mm, "<<parroom[2]<<"*mm"<<endl;
  cout <<"parairroom "<<parairroom[0]<<"*mm, "<<parairroom[1]<<"*mm, "<<parairroom[2]<<"*mm"<<endl; 
  cout <<"parairroom2 "<<parairroom2[0]<<"*mm, "<<parairroom2[1]<<"*mm, "<<parairroom2[2]<<"*mm"<<endl;
  cout <<"parstaircaseair "<<parstaircaseair[0]<<"*mm, "<<parstaircaseair[1]<<"*mm, "<<parstaircaseair[2]<<"*mm"<<endl;
  cout <<"parstaircasel "<<parstaircasel[0]<<"*mm, "<<parstaircasel[1]<<"*mm, "<<parstaircasel[2]<<"*mm"<<endl;
  cout <<"parstaircase "<<parstaircase[0]<<"*mm, "<<parstaircase[1]<<"*mm, "<<parstaircase[2]<<"*mm"<<endl;
  // cout <<" "<<[0]<<"*mm, "<<[1]<<"*mm, "<<[2]<<"*mm"<<endl;
  cout <<"parino "<<parino[0]<<"*mm, "<<parino[1]<<"*mm, "<<parino[2]<<"*mm"<<endl;
  cout <<"parlay "<<parlay[0]<<"*mm, "<<parlay[1]<<"*mm, "<<parlay[2]<<"*mm"<<endl;
  cout <<"parchm "<<parchm[0]<<"*mm, "<<parchm[1]<<"*mm, "<<parchm[2]<<"*mm"<<endl;
  cout <<"parirlay "<<parirlay[0]<<"*mm, "<<parirlay[1]<<"*mm, "<<parirlay[2]<<"*mm"<<endl;
  cout <<"parcoilspacerpc "<<parcoilspacerpc[0]<<"*mm, "<<parcoilspacerpc[1]<<"*mm, "<<parcoilspacerpc[2]<<"*mm"<<endl;
  cout <<"parcoilspaceiron "<<parcoilspaceiron[0]<<"*mm, "<<parcoilspaceiron[1]<<"*mm, "<<parcoilspaceiron[2]<<"*mm"<<endl;
  cout <<"parairgap1 "<<parairgap1[0]<<"*mm, "<<parairgap1[1]<<"*mm, "<<parairgap1[2]<<"*mm"<<endl;
  cout <<"parairgap2 "<<parairgap2[0]<<"*mm, "<<parairgap2[1]<<"*mm, "<<parairgap2[2]<<"*mm"<<endl;
  cout <<"parairgap3 "<<parairgap3[0]<<"*mm, "<<parairgap3[1]<<"*mm, "<<parairgap3[2]<<"*mm"<<endl;
  cout <<"parairgap4 "<<parairgap4[0]<<"*mm, "<<parairgap4[1]<<"*mm, "<<parairgap4[2]<<"*mm"<<endl;
  cout <<"parspacerA "<<parspacerA[0]<<"*mm, "<<parspacerA[1]<<"*mm, "<<parspacerA[2]<<"*mm"<<endl;
  cout <<"parspacerB "<<parspacerB[0]<<"*mm, "<<parspacerB[1]<<"*mm, "<<parspacerB[2]<<"*mm"<<endl;
  cout <<"parspacerC "<<parspacerC[0]<<"*mm, "<<parspacerC[1]<<"*mm, "<<parspacerC[2]<<"*mm"<<endl;
  cout <<"parspacerD "<<parspacerD[0]<<"*mm, "<<parspacerD[1]<<"*mm, "<<parspacerD[2]<<"*mm"<<endl;
  // cout <<"parg10 "<<parg10[0]<<"*mm, "<<parg10[1]<<"*mm, "<<parg10[2]<<"*mm"<<endl;
  cout <<"parfrpbox "<<parfrpbox[0]<<"*mm, "<<parfrpbox[1]<<"*mm, "<<parfrpbox[2]<<"*mm"<<endl;
  cout <<"parairbox "<<parairbox[0]<<"*mm, "<<parairbox[1]<<"*mm, "<<parairbox[2]<<"*mm"<<endl;
  cout <<"parg10Trap1 "<<parg10Trap1[0]<<"*mm, "<<parg10Trap1[1]<<"*mm, "<<parg10Trap1[2]<<"*mm, "<<parg10Trap1[3]<<"*mm, "<<parg10Trap1[4]<<"*mm, "<<parg10Trap1[5]<<"*mm, "<<parg10Trap1[6]<<"*mm, "<<parg10Trap1[7]<<"*mm"<<endl;
  cout <<"parg10Trap2 "<<parg10Trap2[0]<<"*mm, "<<parg10Trap2[1]<<"*mm, "<<parg10Trap2[2]<<"*mm, "<<parg10Trap2[3]<<"*mm, "<<parg10Trap2[4]<<"*mm, "<<parg10Trap2[5]<<"*mm, "<<parg10Trap2[6]<<"*mm, "<<parg10Trap2[7]<<"*mm"<<endl;
  cout <<"paral "<<paral[0]<<"*mm, "<<paral[1]<<"*mm, "<<paral[2]<<"*mm"<<endl;
  cout <<"paralCutBig "<<paralCutBig[0]<<"*mm, "<<paralCutBig[1]<<"*mm, "<<paralCutBig[2]<<"*mm, "<<paralCutBig[3]<<"*mm"<<endl;
  cout <<"paralCutSmall "<<paralCutSmall[0]<<"*mm, "<<paralCutSmall[1]<<"*mm, "<<paralCutSmall[2]<<"*mm, "<<paralCutSmall[3]<<"*mm"<<endl;
  cout <<"parhoneycomb "<<parhoneycomb[0]<<"*mm, "<<parhoneycomb[1]<<"*mm, "<<parhoneycomb[2]<<"*mm"<<endl;
  cout <<"parhoneycombCutBig "<<parhoneycombCutBig[0]<<"*mm, "<<parhoneycombCutBig[1]<<"*mm, "<<parhoneycombCutBig[2]<<"*mm, "<<parhoneycombCutBig[3]<<"*mm"<<endl;
  cout <<"parhoneycombCutSmall "<<parhoneycombCutSmall[0]<<"*mm, "<<parhoneycombCutSmall[1]<<"*mm, "<<parhoneycombCutSmall[2]<<"*mm, "<<parhoneycombCutSmall[3]<<"*mm"<<endl;
  cout <<"parcup "<<parcup[0]<<"*mm, "<<parcup[1]<<"*mm, "<<parcup[2]<<"*mm"<<endl;
  cout <<"parcupCutBig "<<parcupCutBig[0]<<"*mm, "<<parcupCutBig[1]<<"*mm, "<<parcupCutBig[2]<<"*mm, "<<parcupCutBig[3]<<"*mm"<<endl;
  cout <<"parcupCutSmall "<<parcupCutSmall[0]<<"*mm, "<<parcupCutSmall[1]<<"*mm, "<<parcupCutSmall[2]<<"*mm, "<<parcupCutSmall[3]<<"*mm"<<endl;
  cout <<"parmylar "<<parmylar[0]<<"*mm, "<<parmylar[1]<<"*mm, "<<parmylar[2]<<"*mm"<<endl;
  cout <<"parmylarCutBig "<<parmylarCutBig[0]<<"*mm, "<<parmylarCutBig[1]<<"*mm, "<<parmylarCutBig[2]<<"*mm, "<<parmylarCutBig[3]<<"*mm"<<endl;
  cout <<"parmylarCutSmall "<<parmylarCutSmall[0]<<"*mm, "<<parmylarCutSmall[1]<<"*mm, "<<parmylarCutSmall[2]<<"*mm, "<<parmylarCutSmall[3]<<"*mm"<<endl;
  cout <<"parcoat "<<parcoat[0]<<"*mm, "<<parcoat[1]<<"*mm, "<<parcoat[2]<<"*mm"<<endl;
  cout <<"parcoatCutBig "<<parcoatCutBig[0]<<"*mm, "<<parcoatCutBig[1]<<"*mm, "<<parcoatCutBig[2]<<"*mm, "<<parcoatCutBig[3]<<"*mm"<<endl;
  cout <<"parcoatCutSmall "<<parcoatCutSmall[0]<<"*mm, "<<parcoatCutSmall[1]<<"*mm, "<<parcoatCutSmall[2]<<"*mm, "<<parcoatCutSmall[3]<<"*mm"<<endl;
  cout <<"parqurz "<<parqurz[0]<<"*mm, "<<parqurz[1]<<"*mm, "<<parqurz[2]<<"*mm"<<endl;
  cout <<"parqurzCutBig "<<parqurzCutBig[0]<<"*mm, "<<parqurzCutBig[1]<<"*mm, "<<parqurzCutBig[2]<<"*mm, "<<parqurzCutBig[3]<<"*mm"<<endl;
  cout <<"parqurzCutSmall "<<parqurzCutSmall[0]<<"*mm, "<<parqurzCutSmall[1]<<"*mm, "<<parqurzCutSmall[2]<<"*mm, "<<parqurzCutSmall[3]<<"*mm"<<endl;
  cout <<"pargas "<<pargas[0]<<"*mm, "<<pargas[1]<<"*mm, "<<pargas[2]<<"*mm"<<endl;
  cout <<"pargasCutBig "<<pargasCutBig[0]<<"*mm, "<<pargasCutBig[1]<<"*mm, "<<pargasCutBig[2]<<"*mm, "<<pargasCutBig[3]<<"*mm"<<endl;
  cout <<"pargasCutSmall "<<pargasCutSmall[0]<<"*mm, "<<pargasCutSmall[1]<<"*mm, "<<pargasCutSmall[2]<<"*mm, "<<pargasCutSmall[3]<<"*mm"<<endl;
  cout <<"parvcoil "<<parvcoil[0]<<"*mm, "<<parvcoil[1]<<"*mm, "<<parvcoil[2]<<"*mm"<<endl;
  cout <<"parhcoil "<<parhcoil[0]<<"*mm, "<<parhcoil[1]<<"*mm, "<<parhcoil[2]<<"*mm"<<endl;
  cout <<"parcurvedcoil "<<parcurvedcoil[0]<<"*mm, "<<parcurvedcoil[1]<<"*mm, "<<parcurvedcoil[2]<<"*mm, "<<parcurvedcoil[3]<<"*mm, "<<parcurvedcoil[4]<<"*mm"<<endl;
  cout <<"parcoilsupport "<<parcoilsupport[0]<<"*mm, "<<parcoilsupport[1]<<"*mm, "<<parcoilsupport[2]<<"*mm"<<endl;
  cout << "------------------------------------------------------------"<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void micalDetectorConstruction::SetUniformMagField(G4double fieldValue) {
  //apply a global uniform magnetic field along Y axis
  //  G4FieldManager* 
  fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  
  if(magField) { delete magField;}		//delete the existing magn field
  
  if(fieldValue!=0.0) {			// create a new one if non nul
    magField = new G4UniformMagField(G4ThreeVector(0.0, fieldValue, 0.)); //Field along y-axis
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  } else {
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }
  // fEquation =  new G4EqMagElectricField(magField);
  // pStepper = new G4SimpleRunge(fEquation);
  // G4ChordFinder *pChordFinder = new G4ChordFinder((G4UniformMagField*)magField,1.e-1*mm, pStepper);
  // G4bool fieldLocal1 = true;
  // fieldMgr = new G4FieldManager(magField,pChordFinder,fieldLocal1);
  // G4TransportationManager::GetTransportationManager()->SetFieldManager(fieldMgr);

}

void micalDetectorConstruction::ConstructFieldMap() {
  static G4bool fieldIsInitialized = false;
  //And finally that it was not initialized previously
  if (!fieldIsInitialized) {
    inoical0Field=new micalElectroMagneticField(); //GMA14 ("file_magField.dat");
    
    G4double field = inoical0Field->GetConstantFieldvalue();
    if (field == 0) {
      inoical0Field = NULL;
      cout << "***************************" << endl;
      cout << "*                         *" << endl;
      cout << "*  Magnetic Field is off  *" << endl;
      cout << "*                         *" << endl;
      cout << "***************************" << endl;
    } else {
      cout << "***************************" << endl;
      cout << "*                         *" << endl;
      cout << "*  Magnetic Field is on   *" << endl;
      cout << "*                         *" << endl;
      cout << "***************************" << endl;
      cout << endl;
      cout << " Field Value " << field << endl;
    }
    
    fieldMgr= G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(inoical0Field);
    
    // G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs((G4MagneticField*)inoical0Field);
    fEquation =  new G4EqMagElectricField(inoical0Field);
    
    //pStepper = new G4ExplicitEuler(fEquation);
    //pStepper = new G4ImplicitEuler(fEquation);
    pStepper = new G4SimpleRunge(fEquation);
    //pStepper = new G4ClassicalRK4(fEquation);
    
    //pStepper = new G4SimpleHeum( fEquation );         
    //pStepper = new G4HelixExplicitEuler( fEquation ); 
    //pStepper = new G4HelixImplicitEuler( fEquation ); 
    //pStepper = new G4HelixSimpleRunge( fEquation );   
    //pStepper = new G4CashKarpRKF45( fEquation );      
    //pStepper = new G4RKG3_Stepper( fEquation );       
    
    //    G4ChordFinder *pChordFinder = new G4ChordFinder((G4MagneticField*)inoical0Field,
    //                                                    1.e-1*mm, pStepper);
    
    pIntgrDriver = new G4MagInt_Driver(0.000001*mm,pStepper,pStepper->GetNumberOfVariables() );
    pChordFinder = new G4ChordFinder(pIntgrDriver);
    
    pChordFinder->SetDeltaChord(1.0e-3*mm);
    fieldMgr->SetChordFinder(pChordFinder);
    fieldMgr->SetDeltaOneStep(1.0e-3*mm);
    fieldMgr->SetDeltaIntersection(1.0e-4*mm);
    G4PropagatorInField* fieldPropagator= G4TransportationManager::GetTransportationManager()->GetPropagatorInField();
    fieldPropagator->SetMinimumEpsilonStep(1.e-5*mm);
    fieldPropagator->SetMaximumEpsilonStep(1.e-2*mm);
    fieldIsInitialized = true;
  }
  
#ifdef debug
  G4cout	<< "micalDetectorConstruction: Starting timer!!!"<< G4endl;
  G4Timer timer;
  timer.Start();
#endif
}

G4SubtractionSolid* micalDetectorConstruction::ConstructRPCBox(float* parBox, float* parCutBig, float* parCutSmall) {
  // Rotation for Cutezoids
  G4RotationMatrix rotmBigCut1 = G4RotationMatrix();
  rotmBigCut1.rotateZ(-45*deg);
  G4RotationMatrix rotmBigCut2 = G4RotationMatrix();
  rotmBigCut2.rotateZ(45*deg);
  G4RotationMatrix rotmSmallCut1 = G4RotationMatrix();
  rotmSmallCut1.rotateZ(-45*deg);
  G4RotationMatrix rotmSmallCut2 = G4RotationMatrix();
  rotmSmallCut2.rotateZ(45*deg);
  
  
  // RPC GAS
  G4Box* solidRect = new G4Box("Rect", parBox[0], parBox[1], parBox[2]);
  // Cut1 Big
  G4Box* solidCutBig1 = new G4Box("CutBig1",parCutBig[1],parCutBig[2],parCutBig[3]);
  G4Box* solidCutBig2 = new G4Box("CutBig2",parCutBig[1],parCutBig[2],parCutBig[3]);
  // Cut2 Small
  G4Box* solidCutSmall1 = new G4Box("CutSmall1",parCutSmall[1],parCutSmall[2],parCutSmall[3]);
  G4Box* solidCutSmall2 = new G4Box("CutSmall2",parCutSmall[1],parCutSmall[2],parCutSmall[3]);

  // Big Gas Cut Transform
  G4ThreeVector positionBigCut1 = G4ThreeVector(-parBox[0], -parBox[1],0);
  G4Transform3D transformBigCut1 = G4Transform3D(rotmBigCut1,positionBigCut1);
  G4ThreeVector positionBigCut2 = G4ThreeVector(parBox[0], parBox[1],0);
  G4Transform3D transformBigCut2 = G4Transform3D(rotmBigCut2,positionBigCut2);			

  // Small Gas Cut Transform
  G4ThreeVector positionSmallCut1 = G4ThreeVector(-parBox[0],parBox[1],0);
  G4Transform3D transformSmallCut1 = G4Transform3D(rotmSmallCut1,positionSmallCut1);			
  G4ThreeVector positionSmallCut2 = G4ThreeVector(parBox[0],-parBox[1],0);
  G4Transform3D transformSmallCut2 = G4Transform3D(rotmSmallCut2,positionSmallCut2);			

  G4SubtractionSolid* solidSub1 =  new G4SubtractionSolid("SubSolid1", solidRect, solidCutBig1,transformBigCut1);

  G4SubtractionSolid* solidSub2 =  new G4SubtractionSolid("SubSolid2", solidSub1, solidCutBig2,transformBigCut2);

  G4SubtractionSolid* solidSub3 =  new G4SubtractionSolid("SubSolid3", solidSub2, solidCutSmall1,transformSmallCut1);
  
  G4SubtractionSolid* solidRPCBox = new G4SubtractionSolid("FinalSubtraction", solidSub3, solidCutSmall2,transformSmallCut2);

  // return solidSub1;
  // return solidSub2;
  // return solidSub3;
  return solidRPCBox;

}



#include "G4RunManager.hh"

void micalDetectorConstruction::UpdateGeometry() {
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

