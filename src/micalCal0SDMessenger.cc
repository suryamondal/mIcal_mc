// $Id: micalDetectorMessenger.cc,v 1.9 2003/09/15 15:38:18 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "micalCal0SDMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalcal0SDMessenger::micalcal0SDMessenger(
                                           micalcal0SD* cal0SD)
:micalcal0SDptr(cal0SD) { 
  micalDir = new G4UIdirectory("/mical/");
  micalDir->SetGuidance("UI commands of this example");
  
  cal0SDDir = new G4UIdirectory("/mical/cal0SD/");
  cal0SDDir->SetGuidance("digi control");

  CorrTimeSmrCmd = new G4UIcmdWithADoubleAndUnit("/mical/cal0SD/CorrTimeSmr",this);
  CorrTimeSmrCmd->SetGuidance("Set Correlated Time Smear");
  CorrTimeSmrCmd->SetParameterName("CorrTimeSmr",true, true);
  CorrTimeSmrCmd->SetDefaultUnit("ns");
  CorrTimeSmrCmd->SetDefaultValue(1.0);
  CorrTimeSmrCmd->SetUnitCategory("Time");

  UnCorrTimeSmrCmd = new G4UIcmdWithADoubleAndUnit("/mical/cal0SD/UnCorrTimeSmr",this);
  UnCorrTimeSmrCmd->SetGuidance("Set UnCorrelated Time Smear");
  UnCorrTimeSmrCmd->SetParameterName("UnCorrTimeSmr",true, true);
  UnCorrTimeSmrCmd->SetDefaultUnit("ns");
  UnCorrTimeSmrCmd->SetDefaultValue(1.0);
  UnCorrTimeSmrCmd->SetUnitCategory("Time");

  TimeToDigiConvCmd = new G4UIcmdWithADouble("/mical/cal0SD/TimeToDigiConv",this);
  TimeToDigiConvCmd->SetGuidance("Set time to digi Conv Assuming Minimum scale of timing ~100 ps = 0.1 ns");
  TimeToDigiConvCmd->SetParameterName("TimeToDigiConv",true, true);
  TimeToDigiConvCmd->SetDefaultValue(0.1);

  SigSpeedCmd = new G4UIcmdWithADouble("/mical/cal0SD/SigSpeed",this);
  SigSpeedCmd->SetGuidance("Set signal speed in the strip in the units of ns/strip");
  SigSpeedCmd->SetParameterName("SigSpeed",true, true);
  SigSpeedCmd->SetDefaultValue(0.15);

  RootRandomCmd = new G4UIcmdWithAnInteger("/mical/cal0SD/RootRandom",this);
  RootRandomCmd->SetGuidance("To switch root random on or off");
  RootRandomCmd->SetParameterName("RootRandom",true,true);
  RootRandomCmd->SetDefaultValue(0);

  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalcal0SDMessenger::~micalcal0SDMessenger() {
  delete micalDir;
  delete cal0SDDir;
  delete CorrTimeSmrCmd;    delete UnCorrTimeSmrCmd;
  delete SigSpeedCmd;
  delete TimeToDigiConvCmd;
  delete RootRandomCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalcal0SDMessenger::SetNewValue(G4UIcommand* command,G4String newValue) { 
  if( command == CorrTimeSmrCmd )
    { micalcal0SDptr->SetCorrTimeSmear(CorrTimeSmrCmd->GetNewDoubleValue(newValue));}
  if( command == UnCorrTimeSmrCmd )
    { micalcal0SDptr->SetUnCorrTimeSmear(UnCorrTimeSmrCmd->GetNewDoubleValue(newValue));}
  if( command == TimeToDigiConvCmd )
    { micalcal0SDptr->SetTimeToDigiConv(TimeToDigiConvCmd->GetNewDoubleValue(newValue));}
  if( command == SigSpeedCmd )
    { micalcal0SDptr->SetSignalSpeed(SigSpeedCmd->GetNewDoubleValue(newValue));}
  if( command == RootRandomCmd )
    { micalcal0SDptr->SetRootRandom(RootRandomCmd->GetNewIntValue(newValue));}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
