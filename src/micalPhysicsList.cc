#include "micalPhysicsList.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4HadronElasticPhysics.hh"

// #include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"


#include "G4ProcessTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "micalDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalPhysicsList::micalPhysicsList(micalDetectorConstruction* adet):G4VModularPhysicsList(), 
 
 pDet(0)
 
 
  {
	  G4LossTableManager::Instance();
      pDet = adet;
	   
	   // default cut value  (1.0mm) 
  // defaultCutValue = 1.0*mm;
    // G4DataQuestionaire it(photon, neutron, no, no, no, neutronxs);
    G4cout << "<<< Geant4 Physics List: micalPhysicsList " <<G4endl;
    G4cout <<G4endl;
    defaultCutValue = 0.7*mm;
    G4int ver = 1;
    SetVerboseLevel(ver);

    // EM Physics
    RegisterPhysics(new G4EmStandardPhysics(ver));

    // Synchroton Radiation & GN Physics
    RegisterPhysics(new G4EmExtraPhysics(ver));
    // Decays
    RegisterPhysics(new G4DecayPhysics(ver));

    // Hadron physics
    RegisterPhysics(new G4HadronElasticPhysics(ver) ); 
    RegisterPhysics(new G4HadronPhysicsQGSP_BERT_HP(ver));
    
    // Ion Physics
    RegisterPhysics(new G4IonPhysics(ver));

  }
  
  micalPhysicsList::~micalPhysicsList()
  
  {
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  
  void micalPhysicsList::SetCuts() {
  if (verboseLevel >1){
    G4cout << "micalPhysicsList::SetCuts: default cut length : "
      << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  
    SetCutsWithDefault(); 
    
    G4Region* region;
    
    G4String regName;
    
    G4ProductionCuts* cutsb;
    
    regName = "Calor_EBlock";
    
    region = G4RegionStore::GetInstance()->GetRegion(regName);
  cutsb = new G4ProductionCuts;
  cutsb->SetProductionCut(0.01*mm); // same cuts for gamma, e- and e+
  cutsb->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("gamma"));
  
  
  region->SetProductionCuts(cutsb);
  
//}

//Set cut values for absorber materials

//  G4MaterialCutsCouple* cutsAbsorber = new G4MaterialCutsCouple(pDet->GetAbsorberMaterial(), cuts);


//  const G4MaterialCutsCouple* xxx = pDet->GetAbsorbMaterialCut();	
//  G4cout <<" material name "<< xxx->GetMaterial()->GetName()<<G4endl;


   G4ProductionCutsTable* production = G4ProductionCutsTable::GetProductionCutsTable();
   production->UpdateCoupleTable(pDet->Getphysi_World());

 // cuts->SetProductionCut(2*mm);
 // G4MaterialCutsCouple* cutsAbsorber = new G4MaterialCutsCouple(pDet->GetAbsorberMaterial(), cuts);
 // cutsAbsorber->PhysicsTableUpdated();
//  cutsAbsorber->SetProductionCuts(cuts);
 // pDet->GetlogicAbsorber()->SetMaterialCutsCouple(cutsAbsorber);

 // production->UpdateCoupleTable();

//  const G4MaterialCutsCouple* cutsAbsorber =  production->GetMaterialCutsCouple(1);
  	
//  G4cout <<"material "<< cutsAbsorber->GetMaterial()->GetName()<<G4endl;
 // cutsAbsorber->SetProductionCuts((G4ProductionCuts*)cuts);

//  production->ScanAndSetCouple(pDet->GetlogicAbsorber(), (G4MaterialCutsCouple*)cutsAbsorber, allregion->GetRegion("DefaultRegionForTheWorld"));
	

	G4cout <<"couple proct@ SetsCut "<<production->GetTableSize()<<G4endl;
	G4cout <<"couple proct@ SetsCut "<<production->GetTableSize()<<G4endl;
	G4cout <<"couple proct@ SetsCut "<<production->GetTableSize()<<G4endl;
	
   for (unsigned i=0; i<production->GetTableSize(); i++) {
     const G4MaterialCutsCouple* getcouple = production->GetMaterialCutsCouple(i);
     G4cout <<"couple proct "<<i<<" "<<getcouple->GetMaterial()->GetName()<<G4endl; //" "<<getcouple->
  	}	
  		

//  if (this->verboseLevel >0)
  
//  G4VUserPhysicsList::DumpCutValuesTable();  



	
}


