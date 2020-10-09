//
// $Id: ExN02ChamberParameterisation.cc,v 1.8 2003/05/28 09:54:09 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "micalActiveParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
//#include "G4Trd.hh"
#include "G4Box.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalActiveParameterisation::micalActiveParameterisation( G4int    NoChambers, 
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

micalActiveParameterisation::~micalActiveParameterisation()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalActiveParameterisation::ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol) const {
  G4double      Zposition= fStartZ + copyNo * fSpacing;
  G4ThreeVector origin(0,0,Zposition);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalActiveParameterisation::ComputeDimensions
(G4Box& trackerChamber, const G4int copyNo, const G4VPhysicalVolume*) const {
  G4double x1 = fStartHalfX + copyNo*fSpacing*tan(fPolAngle);
  //GMA  G4double x2 = fStartHalfX + (copyNo*fSpacing+2*fHalfWidthZ)*tan(fPolAngle);
  
  trackerChamber.SetXHalfLength(x1);
  //  trackerChamber.SetXHalfLength2(x2);
  trackerChamber.SetYHalfLength(fHalfY);
  //  trackerChamber.SetYHalfLength2(fHalfY);
  trackerChamber.SetZHalfLength(fHalfWidthZ);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


