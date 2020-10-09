#include <iostream>
#include "InoVertex.h"

// for Form()
#include "TString.h"
#include "TMinuit.h"
#include "TVector3.h"

//______________________________________________________________________
InoVertex::InoVertex() {
  fU = 0.;
  fV = 0.;
  fZ = 0.;
  fT = 0.;
  fPlane = -1;
  fDirCosine = TVector3(0,0,0);
  fQbyP = 0;
}

//______________________________________________________________________
InoVertex::InoVertex(const InoVertex &rhs) :
  TObject(rhs),
  fU(rhs.fU),
  fV(rhs.fV),
  fZ(rhs.fZ),
  fT(rhs.fT),
  fPlane(rhs.fPlane),
  fDirCosine(rhs.fDirCosine),
  fQbyP(rhs.fQbyP)
{
}

//______________________________________________________________________
InoVertex::~InoVertex()
{
}

//______________________________________________________________________
Bool_t InoVertex::operator==(const InoVertex& rhs) const {
  if (this->fU == rhs.fU &&
      this->fV == rhs.fV &&
      this->fZ == rhs.fZ &&
      this->fT == rhs.fT &&
      this->fPlane == rhs.fPlane &&
      this->fDirCosine == rhs.fDirCosine &&
      this->fQbyP == rhs.fQbyP) {
    return true;
  }
  return false;
}

//______________________________________________________________________
Int_t InoVertex::GetPlane() const {
  return fPlane;
}

//______________________________________________________________________
Int_t InoVertex::GetRPCmod() const {
  return fRPCmod;
}

//______________________________________________________________________
Double_t InoVertex::GetT() const {
  return fT;
}

//______________________________________________________________________
Double_t InoVertex::GetU() const {
  return fU;
}

//______________________________________________________________________
Double_t InoVertex::GetV() const {
  return fV;
}

//______________________________________________________________________
Double_t InoVertex::GetZ() const {
  return fZ;
}

//______________________________________________________________________
void InoVertex::SetPlane(Int_t ivar) {
  fPlane = ivar;
}

//______________________________________________________________________
void InoVertex::SetRPCmod(Int_t ivar) {
  fRPCmod = ivar;
}

//______________________________________________________________________
void InoVertex::SetT(Double_t dvar) {
  fT = dvar;
}

//______________________________________________________________________
void InoVertex::SetU(Double_t dvar) {
  fU = dvar;
}

//______________________________________________________________________
void InoVertex::SetV(Double_t dvar) {
  fV = dvar;
}

//______________________________________________________________________
void InoVertex::SetZ(Double_t dvar) {
  fZ = dvar;
}

//______________________________________________________________________
const char* InoVertex::AsString(Option_t* /* option */ ) const {
  return Form(" pln=%4d U=%f V=%f Z=%f T(ns)=%f",
              fPlane,fU,fV,fZ,fT); ///Munits::ns);
}

void InoVertex::SetDirCosine(TVector3 dvar) {
  fDirCosine = dvar;
}

TVector3 InoVertex::GetDirCosine() const {
  return fDirCosine;
}


void InoVertex::SetQbyP(Double_t dvar) {
  fQbyP = dvar;
}

Double_t InoVertex::GetQbyP() const {
  return fQbyP;
}
