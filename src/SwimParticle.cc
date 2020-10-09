//
// This class contains all the data of a particle:
// initial position, current position, momentum, mass & charge.
//
#include "SwimParticle.h"
#include "TMath.h"
//......................................................................

SwimParticle::SwimParticle(const TVector3 position,
                           const TVector3 momentum,
                           double         mass,
                           double         charge) :
  fInitPosition(position),
  fPosition(position),
  fMomentum(momentum),
  fMass(mass),
  fCharge(charge),
  fS(0),
  fRange(0),
  fVxB(0)
{}

//......................................................................

const TVector3 SwimParticle::GetInitPosition() const {
  return fInitPosition;
} 

//......................................................................

const TVector3 SwimParticle::GetPosition() const {
  return fPosition;
} 

//......................................................................

const TVector3 SwimParticle::GetMomentum() const {
  return fMomentum;
} 

//......................................................................

TVector3 SwimParticle::GetDirection() const {
  double modulus = SwimParticle::GetMomentumModulus();  
  if (modulus!=0.0)
    return TVector3(fMomentum.X()/modulus,
                    fMomentum.Y()/modulus,
                    fMomentum.Z()/modulus);
  else return TVector3(0.0,0.0,0.0);
}

//......................................................................

double SwimParticle::GetMomentumModulus() const {
  return TMath::Sqrt(fMomentum.X()*fMomentum.X()
		     + fMomentum.Y()*fMomentum.Y()
		     + fMomentum.Z()*fMomentum.Z());
} 

//......................................................................

double SwimParticle::GetEnergy() const {
  double momentumSqr = fMomentum.X()*fMomentum.X()
    + fMomentum.Y()*fMomentum.Y()
    + fMomentum.Z()*fMomentum.Z();
  return TMath::Sqrt(fMass*fMass+momentumSqr);
}

//......................................................................

double SwimParticle::GetMass() const  { 
  return fMass; 
} 

//......................................................................

double SwimParticle::GetCharge() const { 
  return fCharge; 
} 

//......................................................................

double SwimParticle::GetS() const { 
  return fS; 
} 

//......................................................................

double SwimParticle::GetRange() const { 
  return fRange; 
} 

//......................................................................

double SwimParticle::GetVxB() const {
  return fVxB;
}

//......................................................................

void SwimParticle::SetPosition(const TVector3 position) { 
  fPosition = position;
}

//......................................................................

void SwimParticle::SetMomentum(const TVector3 momentum) {
  fMomentum = momentum;
}

//......................................................................

void SwimParticle::SetMass(double mass) {
  fMass = mass;
}

//......................................................................

void SwimParticle::SetCharge(double charge) {
  fCharge = charge;
}

//......................................................................

void SwimParticle::AddS(double s) {
  fS += s;
}

//......................................................................

void SwimParticle::AddRange(double range) {
  fRange += range;
}

void SwimParticle::AddVxB(double VxB) {
  fVxB += VxB;
}

