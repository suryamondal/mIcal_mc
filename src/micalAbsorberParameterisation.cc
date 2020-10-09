//
// $Id: ExN02ChamberParameterisation.cc,v 1.8 2003/05/28 09:54:09 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "micalAbsorberParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Trd.hh"
//#include "G4Box.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalAbsorberParameterisation::micalAbsorberParameterisation( G4int    NoChambers, 
							      G4double startZ, 
							      G4double halfWidthZ,
							      G4double spacing,
							      G4double polAngle,
							      G4double starthalfX,
							      G4double halfY) {
  fNoChambers =  NoChambers; 
  fStartZ     =  startZ; 
  fHalfWidthZ =  halfWidthZ;
  fSpacing    =  spacing;
  fPolAngle   =  polAngle;
  fStartHalfX =  starthalfX;
  fHalfY      =  halfY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalAbsorberParameterisation::~micalAbsorberParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalAbsorberParameterisation::ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const {
  G4double      Zposition= fStartZ + copyNo * fSpacing;
  G4ThreeVector origin(0,0,Zposition);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalAbsorberParameterisation::ComputeDimensions (G4Trd& trackerChamber, const G4int copyNo, const G4VPhysicalVolume*) const {
  G4double x1 = fStartHalfX + copyNo*fSpacing*tan(fPolAngle);
  G4double x2 = fStartHalfX + (copyNo*fSpacing+2*fHalfWidthZ)*tan(fPolAngle);
  
  trackerChamber.SetXHalfLength1(x1);
  trackerChamber.SetXHalfLength2(x2);
  trackerChamber.SetYHalfLength1(fHalfY);
  trackerChamber.SetYHalfLength2(fHalfY);
  trackerChamber.SetZHalfLength(fHalfWidthZ);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


