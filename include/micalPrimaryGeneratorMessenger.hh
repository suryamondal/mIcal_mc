#ifndef micalPrimaryGeneratorMessenger_h
#define micalPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include "micalPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class micalPrimaryGeneratorMessenger: public G4UImessenger {
public:
  micalPrimaryGeneratorMessenger(micalPrimaryGeneratorAction*);
  ~micalPrimaryGeneratorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  micalPrimaryGeneratorAction* micalAction;
  
  G4UIdirectory*               mIcalDir;
  G4UIdirectory*               gunDir; 

  G4UIcmdWithAnInteger*        RunNumberCmd;
  G4UIcmdWithAnInteger*        GeneratorCmd;
  G4UIcmdWithABool*            inputCmd;
  
  G4UIcmdWithAnInteger*        FirstEvtCmd;
  G4UIcmdWithAString*          RndmCmd;
  G4UIcmdWithAString*          CorsikaFileDirCmd;
  G4UIcmdWithAString*          CorsikaFileNameCmd;
  G4UIcmdWithAString*          FluxFileNameCmd;
  // InputFlag == 0 commands
  G4UIcmdWithAnInteger*        partIdCmd;
  G4UIcmdWithADoubleAndUnit*   incEnergyCmd;
  G4UIcmdWithADoubleAndUnit*   incEnergySmrCmd;
  G4UIcmdWith3Vector*          incDirectionCmd;
  G4UIcmdWithADoubleAndUnit*   incThetaSmrCmd;
  G4UIcmdWithADoubleAndUnit*   incPhiSmrCmd;
  G4UIcmdWith3VectorAndUnit*   incPositionCmd;
  G4UIcmdWithADoubleAndUnit*   incVxSmrCmd;
  G4UIcmdWithADoubleAndUnit*   incVySmrCmd;
  G4UIcmdWithADoubleAndUnit*   incVzSmrCmd;
  
  // Input Flag == 1 commands
  G4UIcmdWithADouble* PowerCosThetaCmd;
  G4UIcmdWithADouble* DetectorThetaCoverCmd;
  G4UIcmdWithADouble* PowerCosmicEnergyCmd;
  G4UIcmdWithADouble* ELowLimitCmd;
  G4UIcmdWithADouble* EUpLimitCmd;
  G4UIcmdWithAnInteger* TopTrgLayCmd;
  G4UIcmdWithAnInteger* BottomTrgLayCmd;


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

