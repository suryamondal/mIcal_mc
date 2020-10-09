//
// $Id: micalRunAction.cc,v 1.15 2003/11/25 16:50:13 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "micalRunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalRunAction::micalRunAction()
{//MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
  // cout<<"micalRunAction::micalRunAction()..."<<endl;
theRunActMessenger = new micalRunActionMessenger(this);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalRunAction::~micalRunAction() {
  delete theRunActMessenger; theRunActMessenger=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalRunAction::BeginOfRunAction(const G4Run* aRun) {
  // cout<<"In micalRunAction::BeginOfRunAction(const G4Run* aRun) "<<endl;
  MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
  pAnalysis->InputOutput=InputOutput;
  pAnalysis->isVisOut=isVisOut;
  pAnalysis->isXtermOut=isXtermOut;
  pAnalysis->FirstEvt=FirstEvt;
  pAnalysis->collatedIn=collatedIn;
  G4cout << "### RunID "<< aRun->GetRunID() <<" Nevt " <<aRun->GetNumberOfEvent()<<" I/O " <<InputOutput<< G4endl;
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  char runno[40];
  int ij=aRun->GetRunID();
  //G4String filename="nuance_nc_9gev_";  //body of the file name
  sprintf(runno,"%d",ij*10); //body of the file name
  //itoa(ij,runno,10);
  
  //Directory of the input an output file:
  if (pAnalysis->InputOutput==0 ||pAnalysis->InputOutput==1 || pAnalysis->InputOutput==2) {
    input_title = a_dir_title ;//"../../../Detector_Resolution/Detresolution_010110/Nuance_Energy_9gev/NC1/";
    output_title= r_dir_title; //"../Analysis_NC/";
  } else if (pAnalysis->InputOutput==3 ||pAnalysis->InputOutput==4) {
    input_title = a_dir_title; //"../Analysis_NC/";
    output_title= r_dir_title; //"../Analysis_NC/";
  } else if (pAnalysis->InputOutput==5) {
    input_title	= a_dir_title; //"../Analysis_NC/";
    output_title= r_dir_title; //"../Analysis_NC/";
  }

  collated_title = a_dir_title;
  
  input_title.append(a_file_title);
  //input_title.append(runno);
  
  output_title.append(r_file_title);
  //output_title.append(runno);
  
  if(collatedIn) {
    collated_title.append(c_file_title);
  }
  
  //cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
  if (pAnalysis->InputOutput==0 ||pAnalysis->InputOutput==1 || pAnalysis->InputOutput==2) { 
    /*input_title.append(".dat");*/
  } else if (pAnalysis->InputOutput==3 ||pAnalysis->InputOutput==4) {
    input_title.append("_sim.root");
  } else if (pAnalysis->InputOutput==5) {
    input_title.append("_digi.root");
  }
  if (pAnalysis->InputOutput==0 ||pAnalysis->InputOutput==3 || pAnalysis->InputOutput==5) {
    output_title.append("_reco");
  } else if (pAnalysis->InputOutput==1 ||pAnalysis->InputOutput==4) {
    output_title.append("_digi");
  } else if (pAnalysis->InputOutput==2)	{
    output_title.append("_sim");
  }
  pAnalysis->OpenRootfiles(input_title,output_title,collated_title); //VALGRIND
  cout <<"input_title = "<< input_title<<endl;
  cout <<"output_title = "<<output_title<<endl;
  //cout<<"#############################"<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalRunAction::EndOfRunAction(const G4Run* ) {
  // cout<<" In micalRunAction::EndOfRunAction(const G4Run* ) ..."<<endl;
  MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
  // cout<<"MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;"<<endl;
  pAnalysis->CloseRootfiles();
  // cout<<"pAnalysis->CloseRootfiles();"<<endl;
  //micalPrimaryGeneratorAction *pgPointer = micalPrimaryGeneratorAction::AnPointer;
  //if(r_file){r_file->cd();r_file->Close();}
  //if(a_file) a_file.close();
  //pgPointer->nuance_input.close();
  //delete pAnalysis;
}

void micalRunAction::SetRunID(G4int run_id) {
  run_ID = run_id ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
