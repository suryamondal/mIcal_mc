// $Id: micalRunAction.hh,v 1.8 2003/09/15 15:38:14 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef micalRunAction_h
#define micalRunAction_h 1
#include "TFile.h"
#include "TH1.h"
#include "G4UserRunAction.hh"
#include "micalRunActionMessenger.hh"
#include "MultiSimAnalysis.hh"
#include "globals.hh"
#include "iostream"
#include <fstream>
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;

class micalRunAction : public G4UserRunAction
{
  public:
    micalRunAction();
   ~micalRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
    void SetRunID(G4int run_id);
    void SetInputDirectory(G4String p) {  a_dir_title = p;}
    void SetOutputDirectory(G4String p) { r_dir_title = p;}

    void SetInputFile(G4String p) {  a_file_title = p;}
    void SetOutputFile(G4String p) { r_file_title = p;}

  void SetCollatedIn(G4int p) { collatedIn = p; }
  void SetCollatedFile(G4String p) { c_file_title = p; }
  
    void SetFirstEvt(G4int p) { FirstEvt = p;}
    void SetInputOutput(G4int p) {InputOutput=p ;}
    void SetisVisOut(G4int p) {isVisOut=p ;}
    void SetisXtermOut(G4int p) {isXtermOut=p ;}
    TFile *r_file;
    ofstream a_file;
        
  //    MultiSimAnalysis *pAnalysis;
    micalRunActionMessenger *theRunActMessenger;
    private:
    G4String input_title;
    G4String output_title;
  G4String collated_title;
    G4String a_dir_title;
    G4String r_dir_title;
    G4String a_file_title;
    G4String r_file_title;
  G4String c_file_title;

    G4int run_no;
    G4int evt_no;
    G4int run_ID;

    G4int FirstEvt;
    G4int InputOutput;
  G4int collatedIn;
    G4int isVisOut;
    G4int isXtermOut;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

