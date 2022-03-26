// This class contains all the data of a particle:
// initial position, current position, momentum, mass & charge.
// Momentum is always forward in time.
//
// Units are: position  : meters
//            momentum  : GeV/c
//            mass      : GeV/c^2
//
#ifndef SWIMPARTICLE_H
#define SWIMPARTICLE_H

#include "TVector3.h"           // ROOT TVector class
#include "TMinuit.h"

class SwimParticle

//......................................................................

{
public:
  SwimParticle() {;}
  virtual ~SwimParticle() {;}
  SwimParticle(const TVector3 position,
  const TVector3 momentum,
  double         mass = 0.105658357, //*TMinuit::GeV,
  double         charge = -1.0);
  const TVector3 GetInitPosition()    const;
  const TVector3 GetPosition()        const;
  const TVector3 GetMomentum()        const;
  TVector3       GetDirection()       const;
  double         GetMomentumModulus() const;
  double         GetEnergy()          const;
  double         GetMass()            const;
  double         GetCharge()          const;
  double         GetS()               const;
  double         GetRange()           const;
  double         GetVxB()         const;

  void SetPosition(const TVector3 position);
  void SetMomentum(const TVector3 momentum);
  void SetMass(double mass);
  void SetCharge(double charge);
  void AddS(double s);
  void AddRange(double range);
  void AddVxB(double VxB);

private:
  TVector3 fInitPosition; // Initial Position of a particle
  TVector3 fPosition;     // Position of a particle
  TVector3 fMomentum;     // Momentum of a particle
  double   fMass;         // Mass of a particle
  double   fCharge;       // Charge of a particle
  double   fS;            // Path
  double   fRange;        // Range
  double   fVxB;      // integrated v x B
};

#endif // SWIMPARTICLE_H
