//
// $Id: micalEventActionMessenger.hh,v 1.7 2002/12/16 16:37:26 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef micalEventActionMessenger_h
#define micalEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class micalEventAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class micalEventActionMessenger: public G4UImessenger
{
  public:
    micalEventActionMessenger(micalEventAction*);
   ~micalEventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    micalEventAction*     eventAction;
    G4UIdirectory*        eventDir;   
    G4UIcmdWithAString*   DrawCmd;
    G4UIcmdWithAnInteger* PrintCmd;    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
