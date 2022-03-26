#ifndef INOVERTEX_H
#define INOVERTEX_H

#include "TObject.h"
#include "TVector3.h"
class InoVertex : public TObject
{

public:
  InoVertex();
  InoVertex(const InoVertex &rhs);
  virtual ~InoVertex();

  Bool_t operator==(const InoVertex& rhs) const;

  Int_t GetPlane() const;
  Int_t GetRPCmod() const;
  Double_t GetT() const;
  Double_t GetU() const;
  Double_t GetV() const;
  Double_t GetZ() const;
  TVector3 GetDirCosine() const;
  Double_t GetQbyP() const;

  void SetPlane(Int_t);
  void SetRPCmod(Int_t);
  void SetT(Double_t);
  void SetU(Double_t);
  void SetV(Double_t);
  void SetZ(Double_t);
  void SetDirCosine(TVector3 t3);
  void SetQbyP(Double_t);

  const char* AsString(Option_t *option="") const;

private:
  Double_t fU;
  Double_t fV;
  Double_t fZ;
  Double_t fT;
  Int_t fPlane;
  Int_t fRPCmod;
  TVector3 fDirCosine;
  Double_t fQbyP;

};

#endif                                                       // INOVERTEX_H
