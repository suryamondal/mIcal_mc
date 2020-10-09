
#ifndef micalElectroMagneticField_H
#define micalElectroMagneticField_H
#include "TGeoManager.h"
#include "G4ElectroMagneticField.hh"
#include "micalDetectorParameterDef.hh"
#include "MultiSimAnalysis.hh"
#include "micalFieldPropagator.hh"
//#include <TFile.h>
//#include <TTree.h>
#include "G4ThreeVector.hh"
#include <fstream>
using namespace std;


class G4FieldManager;

class micalElectroMagneticField: public G4ElectroMagneticField {
public:
  micalElectroMagneticField();
  //  micalElectroMagneticField(const G4String &name); //GMA14
  ~micalElectroMagneticField();  
  
  static micalElectroMagneticField* EMFPointer;
  // Access functions
  void ElectroMagneticField(const double Point[3], double Bfield[6]) const;
  G4ThreeVector MagneticField(const G4ThreeVector Point) const;
  G4ThreeVector ElectricField(const G4ThreeVector Point) const;
  virtual void GetFieldValue(const double Point[3], double* Bfield) const;
  G4double GetConstantFieldvalue() const {return fval;}
  // void GetConstantFieldvalue(double Point[3],double* Bfield ){
  //                 GetFieldValue(Point,Bfield);}
  G4bool   DoesFieldChangeEnergy() const { return false; }
  MultiSimAnalysis *pAnalysis;
  micalFieldPropagator *pFieldMap;
  micalDetectorParameterDef* paradef;
  
protected:
  // Find the global Field Manager
  G4FieldManager* GetGlobalFieldManager(); 


 int temp;
private:
  G4double   fval;             // Field value
  G4int      npts;             // Number of poinst
  G4double   xoff;             // Offset
  G4double*  pos;              // Position
  G4double*  slope;            // Slope
  G4double*  intercept;        // Intercept
  G4String   filename;         // field map file

  double xpos;
  double ypos;
  double zpos;
  int dofinput;
  int ndata;
  int stepSize;

  double irlayZ; //nFeThickness/2;
  double rpclayZ; //nAirGap/2;
  double ironrpcZ; // nFeThickness + nAirGap;
  int nLayer;
  int nIRLayer;
  double posINOMworldX;
  double posINOMworldY;
  double posINOMworldZ;
  double PosLayerTop;
  double PosLayerBot;
};

#endif
