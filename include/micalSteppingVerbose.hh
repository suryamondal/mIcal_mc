// $Id: micalSteppingVerbose.hh,v 1.8 2003/09/15 15:38:15 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class micalSteppingVerbose;

#ifndef micalSteppingVerbose_h
#define micalSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class micalSteppingVerbose : public G4SteppingVerbose
{
 public:   

   micalSteppingVerbose();
  ~micalSteppingVerbose();

   void StepInfo();
   void TrackingStarted();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
