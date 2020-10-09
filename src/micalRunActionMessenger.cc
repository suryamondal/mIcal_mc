#include "micalRunActionMessenger.hh"
#include "micalRunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include <sstream>

micalRunActionMessenger::micalRunActionMessenger(micalRunAction* aRunAction)
:theRunAction(aRunAction)
{
  runDirectory = new G4UIdirectory("/mical/");
  runDirectory->SetGuidance("UI commands of this example");

  runDir = new G4UIdirectory("/mical/run/");
  runDir->SetGuidance("RunAction control");


  //runDirectory = new G4UIdirectory("/mysetrun/");
  //runDirectory->SetGuidance("My set commands.");

  runIDCmd = new G4UIcommand("/mical/run/SetRunID",this);
  runIDCmd->SetGuidance("Set run ID");
  runIDCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  G4UIparameter* p1 = new G4UIparameter("runID",'i',true);
  p1->SetDefaultValue(1000);
  runIDCmd->SetParameter(p1);

  InputDirCmd = new G4UIcmdWithAString("/mical/run/input_dir",this);
  InputDirCmd->SetGuidance(" Input Directory name");
  InputDirCmd->SetParameterName("input_dir",true, true);
  InputDirCmd->SetDefaultValue("./");

  OutputDirCmd = new G4UIcmdWithAString("/mical/run/output_dir",this);
  OutputDirCmd->SetGuidance(" OutputFile name");
  OutputDirCmd->SetParameterName("output_dir",true, true);
  OutputDirCmd->SetDefaultValue("./");

  InputFileCmd = new G4UIcmdWithAString("/mical/run/input_file",this);
  InputFileCmd->SetGuidance(" InputFile name");
  InputFileCmd->SetParameterName("input_file",true, true);
  InputFileCmd->SetDefaultValue("input_file");

  OutputFileCmd = new G4UIcmdWithAString("/mical/run/output_file",this);
  OutputFileCmd->SetGuidance(" OutputFile name");
  OutputFileCmd->SetParameterName("output_file",true, true);
  OutputFileCmd->SetDefaultValue("output_file");

  CollatedFileCmd = new G4UIcmdWithAString("/mical/run/collated_input_file",this);
  CollatedFileCmd->SetGuidance("Collated InputFile name");
  CollatedFileCmd->SetParameterName("collated_input_file",true, true);
  CollatedFileCmd->SetDefaultValue("collated_input_file");

  FirstEvtCmd = new G4UIcmdWithAnInteger("/mical/run/firstEvt",this);
  FirstEvtCmd->SetGuidance("Starting points at nuance output file");
  FirstEvtCmd->SetParameterName("first_evt",true, true);
  FirstEvtCmd->SetDefaultValue(1);

  InputOutputCmd = new G4UIcmdWithAnInteger("/mical/run/inout",this);
  InputOutputCmd->SetGuidance("Input/output option, 0:GEN->RECO, 1:GEN->DIGI, 2:GEN->SIM, 3: SIM -> RECO, 4: SIM -> DIFI, 5 : DIGI -> RECO");
  InputOutputCmd->SetParameterName("InputOutput",true, true);
  InputOutputCmd->SetDefaultValue(0);

  isVisOutCmd = new G4UIcmdWithAnInteger("/mical/run/isVis",this);
  isVisOutCmd->SetGuidance("Visualisation option, 0:off, 1:ascii file, 2:histogram, 3:histogram+Tgraph");
  isVisOutCmd->SetParameterName("isVisOut",true, true);
  isVisOutCmd->SetDefaultValue(0);

  isXtermOutCmd = new G4UIcmdWithAnInteger("/mical/run/isXterm",this);
  isXtermOutCmd->SetGuidance("Verbose option, 0:normal, 1:ellaborate, 2:debug ");
  isXtermOutCmd->SetParameterName("isXtermOut",true, true);
  isXtermOutCmd->SetDefaultValue(0);

  CollatedCmd = new G4UIcmdWithAnInteger("/mical/run/collated_input",this);
  CollatedCmd->SetGuidance("Collated InputFile name");
  CollatedCmd->SetParameterName("collated_input_file",true, true);
  CollatedCmd->SetDefaultValue(0);

}

micalRunActionMessenger::~micalRunActionMessenger() {
  delete runIDCmd;
  delete runDirectory;
}

void micalRunActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue) {
  const char* nv = (const char*)newValue;
  if( command==runIDCmd ) {
    G4int id;
    std::istringstream is(nv);
    is >> id;
    
    theRunAction->SetRunID(id);
  }
  if( command == InputDirCmd )
    { theRunAction->SetInputDirectory(newValue);}
  
  if( command == OutputDirCmd )
    { theRunAction->SetOutputDirectory(newValue);}
  
  if( command == InputFileCmd )
    { theRunAction->SetInputFile(newValue);}
  
  if( command == OutputFileCmd )
    { theRunAction->SetOutputFile(newValue);}

  if( command == CollatedFileCmd )
    { theRunAction->SetCollatedFile(newValue);} 
  
  if( command == FirstEvtCmd )
    { theRunAction->SetFirstEvt(FirstEvtCmd->GetNewIntValue(newValue));}
  
  if( command == InputOutputCmd )
    { theRunAction->SetInputOutput(InputOutputCmd->GetNewIntValue(newValue));}
  if( command == isVisOutCmd )
    { theRunAction->SetisVisOut(isVisOutCmd->GetNewIntValue(newValue));}
  if( command == isXtermOutCmd )
    { theRunAction->SetisXtermOut(isXtermOutCmd->GetNewIntValue(newValue));}
  if( command == CollatedCmd )
    { theRunAction->SetCollatedIn(CollatedCmd->GetNewIntValue(newValue));}
  
}

