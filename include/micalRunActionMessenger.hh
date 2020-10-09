#ifndef micalRunActionMessenger_h
#define micalRunActionMessenger_h 1

class G4UIdirectory;
class micalRunAction;
#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

class micalRunActionMessenger: public G4UImessenger
{
public:
    micalRunActionMessenger(micalRunAction*);
    ~micalRunActionMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);

  private:
    micalRunAction * theRunAction;
    
  private: //commands
    G4UIdirectory *             runDirectory;
    G4UIdirectory*              runDir;
    G4UIcommand *               runIDCmd ;

    G4UIcmdWithAString*          InputDirCmd;
    G4UIcmdWithAString*          OutputDirCmd;
    G4UIcmdWithAString*          InputFileCmd;
    G4UIcmdWithAString*          OutputFileCmd;
  G4UIcmdWithAString*          CollatedFileCmd;
  
    G4UIcmdWithAnInteger*        FirstEvtCmd;
    G4UIcmdWithAnInteger*        InputOutputCmd;
    G4UIcmdWithAnInteger*        isVisOutCmd;
    G4UIcmdWithAnInteger*        isXtermOutCmd;
   G4UIcmdWithAnInteger*        CollatedCmd;
};

#endif

