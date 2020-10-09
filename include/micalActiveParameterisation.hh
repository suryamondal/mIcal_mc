// $Id: ExN02ChamberParameterisation.hh,v 1.9 2003/11/10 14:29:53 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//
//  A parameterisation that describes a series of boxes along Z
//    The boxes have equal width, & their lengths are a linear equation.
//    They are spaced an equal distance apart, starting from given location.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef micalActiveParameterisation_H
#define micalActiveParameterisation_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"

class G4VPhysicalVolume;
class G4Box;

// Dummy declarations to get rid of warnings ...
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;
class G4Ellipsoid;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class micalActiveParameterisation : public G4VPVParameterisation
{ 
  public:
  
    micalActiveParameterisation(G4int    NoChambers, 
				  G4double startZ, 
				  G4double halfWidthZ, 
				  G4double spacing,
				  G4double polAngle,
				  G4double starthalfX,
				  G4double halfY);
    virtual				 
   ~micalActiveParameterisation();
   
    void ComputeTransformation (const G4int copyNo,
                                G4VPhysicalVolume* physVol) const;
    
    void ComputeDimensions (G4Box & trackerLayer, const G4int copyNo,
                            const G4VPhysicalVolume* physVol) const;

  private:  // Dummy declarations to get rid of warnings ...

    void ComputeDimensions (G4Trd&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Trap&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Cons&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Sphere&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Orb&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Torus&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Para&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Hype&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Tubs&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polycone&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Polyhedra&,const G4int,const G4VPhysicalVolume*) const {}
    void ComputeDimensions (G4Ellipsoid&,const G4int,const G4VPhysicalVolume*) const {}
  private:

    G4int    fNoChambers;   
    G4double fStartZ;     // center Z-co-ordinate for first cell
    G4double fHalfWidthZ;   // HalfWidth of absorber
    G4double fSpacing;      //  The distance between the chambers' center
    G4double fPolAngle;       // Angle of trapizoid
    G4double fStartHalfX;      // Halflength along X-axis for 1st chamber
    G4double fHalfY;            //Half lenght along y-Axix

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


