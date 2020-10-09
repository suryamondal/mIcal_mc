#ifndef micalPhysicsList_h
#define micalPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"
#include "micalDetectorConstruction.hh"
#include "G4ProductionCutsTable.hh"


class micalDetectorConstruction;
class G4VPhysicsConstructor;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class micalPhysicsList: public G4VModularPhysicsList
{
public:
    micalPhysicsList(micalDetectorConstruction* adet);
   ~micalPhysicsList();
   
    void micalAddPhysicsList(const G4String& name);
    
   
   public:
  // SetCuts() 
  virtual void SetCuts();
  void SetRegionCut(G4double);
  void AddPhysicsList(const G4String& name);
  void GetRange(G4double);
  
   private:
 // enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };
   micalDetectorConstruction* pDet;
   G4ProductionCuts* cuts;
};


#endif
