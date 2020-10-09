
// $Id: micalDetectorMessenger.hh,v 1.6 2002/12/16 16:37:26 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef micalCal0SDMessenger_h
#define micalCal0SDMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "micalCal0SD.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithADouble.hh"

class micalcal0SD;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class micalcal0SDMessenger: public G4UImessenger
{
  public:
    micalcal0SDMessenger(micalcal0SD* );
   ~micalcal0SDMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    micalcal0SD* micalcal0SDptr;
    
    G4UIdirectory*             micalDir;
    G4UIdirectory*             cal0SDDir;
    G4UIcmdWithADoubleAndUnit* CorrTimeSmrCmd;
    G4UIcmdWithADoubleAndUnit* UnCorrTimeSmrCmd;
    G4UIcmdWithADouble* CorrIneffiCmd;
    G4UIcmdWithADouble* UnCorrXIneffiCmd;
    G4UIcmdWithADouble* UnCorrYIneffiCmd;
    G4UIcmdWithADouble* TimeToDigiConvCmd;
    G4UIcmdWithADouble* SigSpeedCmd;
    G4UIcmdWithADouble* CorrNoiseCmd1;
    G4UIcmdWithADouble* CorrNoiseCmd2;
    G4UIcmdWithADouble* CorrNoiseCmd3;
    G4UIcmdWithAnInteger* RandomNoiseCmd;
    G4UIcmdWithAnInteger* RootRandomCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

