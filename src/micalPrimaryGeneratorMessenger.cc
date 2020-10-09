
// $Id: micalPrimaryGeneratorMessenger.cc,v 1.8 2002/12/16 16:37:27 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "micalPrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalPrimaryGeneratorMessenger::micalPrimaryGeneratorMessenger(
                                          micalPrimaryGeneratorAction* micalGun)
:micalAction(micalGun)
{
  
  mIcalDir = new G4UIdirectory("/mical/");
  mIcalDir->SetGuidance("UI commands of this example");
  
  gunDir = new G4UIdirectory("/mical/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");

  RunNumberCmd = new G4UIcmdWithAnInteger("/mical/gun/runNumber",this);
  RunNumberCmd->SetGuidance("Run number for these events");
  RunNumberCmd->SetParameterName("Run_number",true, true);
  RunNumberCmd->SetDefaultValue(1);
  
  GeneratorCmd = new G4UIcmdWithAnInteger("/mical/gun/Gen",this);
  GeneratorCmd->SetGuidance(" Is it Nu-Gen input or PID input ?");
  GeneratorCmd->SetGuidance(" Choice 0(G4_MC), 1(Exp. mu Flux), 2(Corsika 3D hist flux), 3(Corsika Event File)");
  GeneratorCmd->SetParameterName("inputId",true, true);
  GeneratorCmd->SetDefaultValue(0);

  FirstEvtCmd = new G4UIcmdWithAnInteger("/mical/gun/firstEvt",this);
  FirstEvtCmd->SetGuidance("Starting points at nuance output file");
  FirstEvtCmd->SetParameterName("first_evt",true, true);
  FirstEvtCmd->SetDefaultValue(1);

  RndmCmd = new G4UIcmdWithAString("/mical/gun/rndm",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on(default), off");
  RndmCmd->SetParameterName("choice",true, true);
  RndmCmd->SetDefaultValue("on");
  RndmCmd->SetCandidates("on off");

  CorsikaFileDirCmd = new G4UIcmdWithAString("/mical/gun/corsika_file_dir",this);
  CorsikaFileDirCmd->SetGuidance(" Corsika directory name");
  CorsikaFileDirCmd->SetParameterName("corsika_dir",true, true);
  CorsikaFileDirCmd->SetDefaultValue("./");
  
  CorsikaFileNameCmd = new G4UIcmdWithAString("/mical/gun/corsika_file",this);
  CorsikaFileNameCmd->SetGuidance(" CorsikaFile name");
  CorsikaFileNameCmd->SetParameterName("corsika_file",true, true);
  CorsikaFileNameCmd->SetDefaultValue("corsika_file");

  partIdCmd = new G4UIcmdWithAnInteger("/mical/gun/pid",this);
  partIdCmd->SetGuidance("Set incident energy of particle ID (PDG)");
  partIdCmd->SetParameterName("partID",true, true);
  partIdCmd->SetDefaultValue(13);

  incEnergyCmd = new G4UIcmdWithADoubleAndUnit("/mical/gun/energy",this);
  incEnergyCmd->SetGuidance("Set incident energy of particle.");
  incEnergyCmd->SetParameterName("Energy",true, true);
  incEnergyCmd->SetDefaultUnit("GeV");
  incEnergyCmd->SetDefaultValue(15.);
  incEnergyCmd->SetUnitCandidates("eV keV MeV GeV TeV PeV");

  incEnergySmrCmd = new G4UIcmdWithADoubleAndUnit("/mical/gun/ensmear",this);
  incEnergySmrCmd->SetGuidance("Set incident energy of particle.");
  incEnergySmrCmd->SetParameterName("EnergySmr",true, true);
  incEnergySmrCmd->SetDefaultUnit("GeV");
  incEnergySmrCmd->SetDefaultValue(0.1);
  incEnergySmrCmd->SetUnitCandidates("eV keV MeV GeV TeV");

  incDirectionCmd = new G4UIcmdWith3Vector("/mical/gun/incdir",this);
  incDirectionCmd->SetGuidance("Set the incident direction ");
  incDirectionCmd->SetParameterName("cosx","cosy","cosz",true, true);
  incDirectionCmd->SetDefaultValue(G4ThreeVector(0.,0.,1.));

  incThetaSmrCmd = new G4UIcmdWithADoubleAndUnit("/mical/gun/thsmear",this);
  incThetaSmrCmd->SetGuidance("Set incident energy of particle.");
  incThetaSmrCmd->SetParameterName("ThetaSmr",true, true);
  incThetaSmrCmd->SetDefaultUnit("mrad");
  incThetaSmrCmd->SetDefaultValue(0.0);
  incThetaSmrCmd->SetUnitCandidates("degree rad mrad");

  incPhiSmrCmd = new G4UIcmdWithADoubleAndUnit("/mical/gun/phsmear",this);
  incPhiSmrCmd->SetGuidance("Set incident energy of particle.");
  incPhiSmrCmd->SetParameterName("PhiSmr",true, true);
  incPhiSmrCmd->SetDefaultUnit("mrad");
  incPhiSmrCmd->SetDefaultValue(1.0);
  incPhiSmrCmd->SetUnitCandidates("degree rad mrad");

  incPositionCmd = new G4UIcmdWith3VectorAndUnit("/mical/gun/incpos",this);
  incPositionCmd->SetGuidance("Set the incident position ");
  incPositionCmd->SetParameterName("x","y","z",true, true);
  incPositionCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
  incPositionCmd->SetDefaultUnit("cm");
  incPositionCmd->SetUnitCandidates("microm mm cm m km");

  incVxSmrCmd = new G4UIcmdWithADoubleAndUnit("/mical/gun/vxsmear",this);
  incVxSmrCmd->SetGuidance("Set uniform smearing in X-position");
  incVxSmrCmd->SetParameterName("VxSmr",true, true);
  incVxSmrCmd->SetDefaultUnit("cm");
  incVxSmrCmd->SetDefaultValue(10);

  incVySmrCmd = new G4UIcmdWithADoubleAndUnit("/mical/gun/vysmear",this);
  incVySmrCmd->SetGuidance("Set uniform smearing in Y-position");
  incVySmrCmd->SetParameterName("VySmr",true, true);
  incVySmrCmd->SetDefaultUnit("cm");
  incVySmrCmd->SetDefaultValue(10);

  incVzSmrCmd = new G4UIcmdWithADoubleAndUnit("/mical/gun/vzsmear",this);
  incVzSmrCmd->SetGuidance("Set uniform smearing in Z-position");
  incVzSmrCmd->SetParameterName("VzSmr",true, true);
  incVzSmrCmd->SetDefaultUnit("cm");
  incVzSmrCmd->SetDefaultValue(10);

  PowerCosThetaCmd = new G4UIcmdWithADouble("/mical/gun/PowerCosTheta",this);
  PowerCosThetaCmd->SetGuidance("Power of cos theta in cosmic muon spectrum");
  PowerCosThetaCmd->SetParameterName("PowerCosTheta",true, true);
  PowerCosThetaCmd->SetDefaultValue(2.15);

  DetectorThetaCoverCmd = new G4UIcmdWithADouble("/mical/gun/DetectorThetaCover",this);
  DetectorThetaCoverCmd->SetGuidance("Detector theta coverage in degrees");
  DetectorThetaCoverCmd->SetParameterName("DetectorThetaCover",true, true);
  DetectorThetaCoverCmd->SetDefaultValue(50);
  
  PowerCosmicEnergyCmd = new G4UIcmdWithADouble("/mical/gun/PowerCosmicEnergy",this);
  PowerCosmicEnergyCmd->SetGuidance("Power of energy dependence in cosmic muon spectrum");
  PowerCosmicEnergyCmd->SetParameterName("PowerCosmicEnergy",true, true);
  PowerCosmicEnergyCmd->SetDefaultValue(2.2);
  
  ELowLimitCmd = new G4UIcmdWithADouble("/mical/gun/ELowLimit",this);
  ELowLimitCmd->SetGuidance("Lower Limit of energy for cosmic muon energy in (GeV)");
  ELowLimitCmd->SetParameterName("ELowLimit",true, true);
  ELowLimitCmd->SetDefaultValue(0.05);
  
  EUpLimitCmd = new G4UIcmdWithADouble("/mical/gun/EUpLimit",this);
  EUpLimitCmd->SetGuidance("Upper Limit of energy for cosmic muon energy in (GeV)");
  EUpLimitCmd->SetParameterName("EUpLimit",true, true);
  EUpLimitCmd->SetDefaultValue(0.05);
  
  TopTrgLayCmd = new G4UIcmdWithAnInteger("/mical/gun/TopTrgLay",this);
  TopTrgLayCmd->SetGuidance("Top Trigger layer of detector");
  TopTrgLayCmd->SetParameterName("TopTrgLay",true, true);
  TopTrgLayCmd->SetDefaultValue(9);

  BottomTrgLayCmd = new G4UIcmdWithAnInteger("/mical/gun/BottomTrgLay",this);
  BottomTrgLayCmd->SetGuidance("Bottom Trigger layer of detector");
  BottomTrgLayCmd->SetParameterName("BottomTrgLay",true, true);
  BottomTrgLayCmd->SetDefaultValue(2);
  

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalPrimaryGeneratorMessenger::~micalPrimaryGeneratorMessenger() {
  if (mIcalDir){delete mIcalDir;}
  if (gunDir){delete gunDir;}

  if (RunNumberCmd){delete RunNumberCmd;}
  if (GeneratorCmd){delete GeneratorCmd;}
  if (FirstEvtCmd){delete FirstEvtCmd;}
  if (RndmCmd){delete RndmCmd;}
  if (CorsikaFileDirCmd){delete CorsikaFileDirCmd;}
  if (CorsikaFileNameCmd){delete CorsikaFileNameCmd;}
  if (partIdCmd){delete  partIdCmd;}
  if (incEnergyCmd){delete incEnergyCmd;}
  if (incEnergySmrCmd){delete incEnergySmrCmd;}
  if (incDirectionCmd){delete  incDirectionCmd;}
  if (incThetaSmrCmd){delete  incThetaSmrCmd;}
  if (incPhiSmrCmd){delete  incPhiSmrCmd;}
  if (incPositionCmd){delete  incPositionCmd;}
  if (incVxSmrCmd){delete incVxSmrCmd;}
  if (incVySmrCmd){delete incVySmrCmd;}
  if (incVzSmrCmd){delete incVzSmrCmd;}
  if (PowerCosThetaCmd){delete PowerCosThetaCmd;}
  if (DetectorThetaCoverCmd){delete DetectorThetaCoverCmd;}
  if (PowerCosmicEnergyCmd){delete PowerCosmicEnergyCmd;}
  if (ELowLimitCmd){delete ELowLimitCmd;}
  if (EUpLimitCmd){delete EUpLimitCmd;}
  if (TopTrgLayCmd){delete TopTrgLayCmd;}
  if (BottomTrgLayCmd){delete BottomTrgLayCmd;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) { 
  
  if( command == RunNumberCmd )
    { micalAction->SetRunNumber(RunNumberCmd->GetNewIntValue(newValue));}
  if( command == GeneratorCmd )
    { micalAction->SetInputFlag(GeneratorCmd->GetNewIntValue(newValue));}
  if( command == FirstEvtCmd )
    { micalAction->SetFirstEvt(FirstEvtCmd->GetNewIntValue(newValue));}
  if( command == RndmCmd )
    { micalAction->SetRndmFlag(newValue);}
  if( command == CorsikaFileDirCmd )
    { micalAction->SetCorsikaFileDir(newValue);}
  if( command == CorsikaFileNameCmd )
    { micalAction->SetCorsikaFileName(newValue);}
  if( command == partIdCmd )
    { micalAction->SetPartId(partIdCmd->GetNewIntValue(newValue));}
  if( command == incEnergyCmd )
    { micalAction->SetIncEnergy(incEnergyCmd->GetNewDoubleValue(newValue));}
  if( command == incEnergySmrCmd )
    { micalAction->SetIncEnergySmr(incEnergySmrCmd->GetNewDoubleValue(newValue));}
  if( command == incDirectionCmd )
    { micalAction->SetIncDirection(incDirectionCmd->GetNew3VectorValue(newValue));}
  if( command == incThetaSmrCmd )
    { micalAction->SetIncThetaSmr(incThetaSmrCmd->GetNewDoubleValue(newValue));}
  if( command == incPhiSmrCmd )
    { micalAction->SetIncPhiSmr(incPhiSmrCmd->GetNewDoubleValue(newValue));}
  if( command == incPositionCmd )
    { micalAction->SetIncPosition(incPositionCmd->GetNew3VectorValue(newValue));}
  if( command == incVxSmrCmd )
    { micalAction->SetIncVxSmr(incVxSmrCmd->GetNewDoubleValue(newValue));}
  if( command == incVySmrCmd )
    { micalAction->SetIncVySmr(incVySmrCmd->GetNewDoubleValue(newValue));}
  if( command == incVzSmrCmd )
    { micalAction->SetIncVzSmr(incVzSmrCmd->GetNewDoubleValue(newValue));}
  if ( command == PowerCosThetaCmd)
    {micalAction->SetPowerCosTheta(PowerCosThetaCmd->GetNewDoubleValue(newValue));}
  if ( command == DetectorThetaCoverCmd)
    {micalAction->SetDetectorThetaCover(DetectorThetaCoverCmd->GetNewDoubleValue(newValue));}
  if ( command == PowerCosmicEnergyCmd)
    {micalAction->SetPowerCosmicEnergy(PowerCosmicEnergyCmd->GetNewDoubleValue(newValue));}
  if ( command == ELowLimitCmd)
    {micalAction->SetELowLimit(ELowLimitCmd->GetNewDoubleValue(newValue));}
  if ( command == EUpLimitCmd)
    {micalAction->SetEUpLimit(EUpLimitCmd->GetNewDoubleValue(newValue));}
  if ( command == TopTrgLayCmd)
    {micalAction->SetTopTrgLay(TopTrgLayCmd->GetNewIntValue(newValue));}
  if ( command == BottomTrgLayCmd)
    {micalAction->SetBottomTrgLay(BottomTrgLayCmd->GetNewIntValue(newValue));}
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

