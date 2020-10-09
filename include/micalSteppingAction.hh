
// $Id: micalSteppingAction.hh,v 1.8 2003/09/15 15:38:14 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef micalSteppingAction_h
#define micalSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "micalElectroMagneticField.hh"

class micalDetectorConstruction;
class micalEventAction;
class MultiSimAnalysis;
class micalElectroMagneticField;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class micalSteppingAction : public G4UserSteppingAction
{
  public:
    micalSteppingAction(micalDetectorConstruction*, micalEventAction*, MultiSimAnalysis* panalysis);
   ~micalSteppingAction();

    void UserSteppingAction(const G4Step*);
    
  private:
    micalDetectorConstruction* detector;
    micalEventAction*          eventaction;  
    MultiSimAnalysis*          pAnalysis;  
  micalElectroMagneticField* ical0Field;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
