///////////////////////////////////////////////////////////////////////////////
// File: micalElectroMagneticField.cc
// Description: User Field class implementation.
///////////////////////////////////////////////////////////////////////////////

#include <fstream>
#include "micalutils.hh"
#include "G4FieldManager.hh"
#include "micalElectroMagneticField.hh"
#include "micalDetectorParameterDef.hh"
#define polint 0
#define interpolate 1

//Constructor and destructor:
//micalElectroMagneticField *micalElectroMagneticField::FdPointer;
//GMA14 micalElectroMagneticField::micalElectroMagneticField(const G4String &file_magField):fval(1), pos(0), slope(0), intercept(0) {
micalElectroMagneticField::micalElectroMagneticField():fval(1), pos(0), slope(0), intercept(0) {
  //  EMFPointer = this;
  paradef = micalDetectorParameterDef::AnPointer;
  pFieldMap	= micalFieldPropagator::FdPointer;
  pAnalysis = MultiSimAnalysis::AnPointer;
  nIRLayer = paradef->GetnIRLayer();
  posINOMworldX = paradef->GetINOMPositionGlobalX();
  posINOMworldY = paradef->GetINOMPositionGlobalY();
  posINOMworldZ = paradef->GetINOMPositionGlobalZ();
  irlayZ = paradef->GetParirlay(2);
  rpclayZ = paradef->GetParlay(2);
  ironrpcZ = irlayZ + rpclayZ;
  PosLayerTop = paradef->GetLayerPosInStack(nIRLayer-1);
   PosLayerBot = paradef->GetLayerPosInStack(0);
  
  
}

micalElectroMagneticField::~micalElectroMagneticField() {
  if (pos)		{delete[] pos;			pos=0;}
  if (slope)		{delete[] slope;		slope=0;}
  if (intercept)	{delete[] intercept;	intercept=0;}
}

// Member functions

void micalElectroMagneticField::ElectroMagneticField(const double x[3], double B[6]) const {
  cout<< "x " <<x[0]<<" "<<x[1]<<endl;
  B[0] = 0*tesla;
  B[1] = 0*tesla;
  B[2] = 0*tesla;
}

G4ThreeVector micalElectroMagneticField::ElectricField(const G4ThreeVector point) const {
  
  G4double x[3],B[6];
  G4ThreeVector v;
  
  x[0] = point.x();
  x[1] = point.y();
  x[2] = point.z();
  micalElectroMagneticField::ElectroMagneticField(x, B);
  v.setX(B[3]);
  v.setY(B[4]);
  v.setZ(B[5]);
  return v;
}

G4ThreeVector micalElectroMagneticField::MagneticField(const G4ThreeVector point) const {
  
  G4double x[3],B[6];
  G4ThreeVector v;
  
  x[0] = point.x();
  x[1] = point.y();
  x[2] = point.z();
  micalElectroMagneticField::ElectroMagneticField(x, B);
  v.setX(B[0]);
  v.setY(B[1]);
  v.setZ(B[2]);
  return v;
}

void micalElectroMagneticField::GetFieldValue(const double x[3], double* B) const {
  //AAR: in this routine, position coordinates are in mm and Magnetic field is in tesla (T).
  
  double Bx= 0.0;
  double By= 0.0;
  
  pFieldMap->ElectroMagneticField(x, Bx, By, 1);
  
  B[0] = Bx;		//GMA-magnetic field
  B[1] = By;		//GMA-magnetic field
  B[2] = 0.0*tesla;	//GMA-magnetic field
  B[3] = 0.0*MeV/cm;
  B[4] = 0.0*MeV/cm;
  B[5] = 0.0*MeV/cm;
  

}
