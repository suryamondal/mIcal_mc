/*
Change in official code
/usr/local/physics/geant4.10.00.p01-install/bin/../include/Geant4/G4NeutronHPVector.hh: In member function -F¡void G4NeutronHPVector::IntegrateAndNormalise()¢:-A
/usr/local/physics/geant4.10.00.p01-install/bin/../include/Geant4/G4NeutronHPVector.hh:403: error: -F¡__isnan¢ is not a member of ¡std¢-A


*/
// $Id: exampleN03.cc,v 1.22 2005/05/26 12:21:05 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <fstream>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4UI_USE_WIN32
#include "G4UIWin32.hh"
#endif

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#include "G4VisManager.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "micalDetectorConstruction.hh"
#include "micalDetectorParameterDef.hh"

//#include "QGSP_BERT.hh"
//#include "QGSP_BERT_HP.hh"
// #include "LHEP.hh"
// failed #include "LBE.hh"
//#include "QBBC.hh"
//#include "QGS_BIC.hh"
// #include "G4GenericPhysicsList.hh"
//#include "QGSP_BIC.hh"


//#include "MyLHEP.hh"
#include "micalPhysicsList.hh"
#include "micalPrimaryGeneratorAction.hh"
#include "micalRunAction.hh"
#include "micalEventAction.hh"
#include "micalSteppingAction.hh"
#include "micalSteppingVerbose.hh"
#include "Randomize.hh"
#include <time.h>
#include "MultiSimAnalysis.hh"
#include "micalFieldPropagator.hh"
#include "G4GDMLParser.hh" //GMA14q
#include "vect_manager.h"
#include "G4VEmFluctuationModel.hh"
#include "G4EnergyLossForExtrapolator.hh"

#include "TTimeStamp.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"

std::ofstream outFile;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// bin/Linux-g++/inoical0 vis.mac 1
// vis.mac -> input macro file
// 1 : after excute through .mac file it will be ready for futher input
//     suitable for interactive mode, 
//     otherwise (0) release pad, suitable for batch mode       

int main(int argc,char** argv) {
  G4cout <<"argc "<<argc<<" "<<argv[0]<<" "<<argv[1]<<" "<<argv[2]<<G4endl;
  // choose the Random engine
  // HepRandom::setTheEngine(new RanecuEngine);
  long seeds[2]={1327511442, 1202219559};// {12334457,1239075};
  // GMAA HepRandom::setTheSeeds(seeds);
  //  G4long myseed = 345354;
  //  HepRandom::setTheSeed(myseed);
  
  G4GDMLParser parser;   //GMA14
  // long seeds[2]={12334457,1239075};
  TTimeStamp tttxx;
  int systime = Long64_t(tttxx.AsDouble()*1e3)%Long64_t(1e9);
  seeds[0] = systime;
  seeds[1] = (systime*G4UniformRand());
  
  // time_t systime = time(NULL);
  // seeds[0] = (long) systime;
  // seeds[1] = (long) (systime*G4UniformRand());
  
  cout << "seed1: " << seeds[0] << "; seed2: " << seeds[1] << endl;
  CLHEP::HepRandom::setTheSeeds(seeds);
  CLHEP::HepRandom::showEngineStatus();
  //InoMuRange_Manager* MuRangeManager = new InoMuRange_Manager();
  
  TTimeStamp ttxx;		// initialising with system time
  gRandom->SetSeed(Long64_t(ttxx.AsDouble()*1e3)%Long64_t(1e9)+seeds[1]*G4UniformRand());
  // gRandom->SetSeed(clock());
  // cout << "clock "<<clock()<<endl;
  cout<<"SetSeed main "<<ttxx.AsDouble()<<" "<<(Long64_t(ttxx.AsDouble()*1000)%Long64_t(1e9))<<" "<<gRandom->GetSeed()<<endl;
  
  
  // my Verbose output class
  G4VSteppingVerbose::SetInstance(new micalSteppingVerbose);
  
  // ROOT and Histogramming Environment
  G4String prefix("mysimul");
  
  //if (argc>=2) {prefix=argv[1];}
  micalDetectorParameterDef* paradef = new micalDetectorParameterDef();
  MultiSimAnalysis *panalysis = new MultiSimAnalysis(prefix);


  InoGeometry_Manager* geoManager = new InoGeometry_Manager(); //GMA14
  cout <<"geometry reading complete..."<<endl;
  
  micalFieldPropagator   *pfield= new  micalFieldPropagator() ;
  
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;
  G4RunManager::GetRunManager()->SetRunIDCounter(1);
  // set mandatory initialization classes
  
  //  micalDetectorParameterDef* paradef = new micalDetectorParameterDef(); //GMA14
  
  micalDetectorConstruction* detector = new micalDetectorConstruction;
  // cout <<"xxx "<<endl;
  runManager->SetUserInitialization(detector);
  // cout <<"1xxx "<<endl;
  runManager->SetUserInitialization(new micalPhysicsList (detector));
  //  runManager->SetUserInitialization(new LHEP);
  //  runManager->SetUserInitialization(new QGSP_BIC);
  //  runManager->SetUserInitialization(new MyLHEP(detector));
  // cout <<"2xxx "<<endl;
  
  // cout <<"3xxx "<<endl;
  // cout <<"4xxx "<<endl;
  // set user action classes
  runManager->SetUserAction(new micalPrimaryGeneratorAction(detector, panalysis));
  // cout <<"5xxx "<<endl;
  runManager->SetUserAction(new micalRunAction);
  // cout <<"6xxx "<<endl;
  micalEventAction* eventaction = new micalEventAction;
  // cout <<"7xxx "<<endl;
  runManager->SetUserAction(eventaction);
  // cout <<"8xxx "<<endl;
  //  runManager->SetUserAction(new micalSteppingAction(detector, eventaction));
  runManager->SetUserAction(new micalSteppingAction(detector, eventaction, panalysis));
  // cout <<"9xxx "<<endl;
  //Initialize G4 kernel
  runManager->Initialize();
  cout <<"10xxx "<<endl;
  
  // //  Write the GDML files
  // parser.Write("detector_world.gdml", G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
  // cout <<"11xxx "<<endl;
  // system("rm geo_mical_world.gdml");
  // cout <<"12xxx "<<endl;
  // system("cp detector_world.gdml geo_mical_world.gdml");
  // cout <<"13xxx "<<endl;
  // system("rm detector_world.gdml");
  // cout <<"14xxx "<<endl;

  //  Write the GDML files
  //parser.Write("detector_world.gdml", G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
  //  parser.Write("detector_log.gdml", G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume()->GetLogicalVolume());  

#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif
  
  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  if (argc!=1)   // batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  else
    {  // interactive mode : define UI session
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute init_vis.mac");
#else
      UImanager->ApplyCommand("/control/execute init.mac");
#endif
      if (ui->IsGUI())
	UImanager->ApplyCommand("/control/execute gui.mac");
      ui->SessionStart();
      delete ui;
#endif
    }
  
  if(panalysis->InputOutput==0|| panalysis->InputOutput==4|| panalysis->InputOutput==5) {
    cout<<" output files created:" <<endl;
  }
  
  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete geoManager;
  delete runManager;
  delete panalysis;
  delete pfield; 
  return 0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
// yum install xerces-c-devel
// yum install libconfig-devel.x86_64

// cmake -DCMAKE_INSTALL_PREFIX=/usr/local/physics/geant4.10.00.p01-install -DGEANT4_USE_SYSTEM_CLHEP=ON -DCLHEP_ROOT_DIR=/usr/local -DGEANT4_USE_GDML=ON -DGEANT4_USE_QT=ON -DGEANT4_USE_OPENGL_X11=ON -DGEANT4_USE_RAYTRACER_X11=ON /usr/local/physics/geant4.10.00.p01

// cmake -DCMAKE_INSTALL_PREFIX=/usr/local/physics/geant4.10.00.p01-install -DGEANT4_USE_SYSTEM_CLHEP=ON -DCLHEP_ROOT_DIR=/usr/local -DGEANT4_USE_GDML=ON -DGEANT4_USE_QT=ON -DGEANT4_USE_OPENGL_X11=ON /usr/local/physics/geant4.10.00.p01

// cmake -DCMAKE_INSTALL_PREFIX=/usr/local/physics/geant4.10.00.p01-install -DGEANT4_USE_SYSTEM_CLHEP=ON -DCLHEP_ROOT_DIR=/usr/local  /usr/local/physics/geant4.10.00.p01

*/
