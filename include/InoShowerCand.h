#ifndef INOSHOWERCAND_H
#define INOSHOWERCAND_H
//CandShowerHandle


class InoHit;
class InoCluster; //InoCluster;
class InoTrackCand; // CandTrackHandle;

class InoShowerCand
{

public:
  InoShowerCand(){;};
  //  InoShowerCand(const InoShowerCand &cdh){;};
  //  InoShowerCand(InoShowerCand *cd){;};
  virtual ~InoShowerCand(){;};
  //  virtual InoShowerCand *DupHandle() const;
  //  virtual void Trace(const char *c = "") const;
  
  typedef enum EShowerType{
    kCC=0,
    kWtCC=1,
    kNC=2,
    kWtNC=3,
    kEM=4
  }ShowerType_t;

  /*
  void SetU(Int_t,Float_t);
  Float_t GetU(Int_t) const; // U position at specified plane

  void SetV(Int_t,Float_t);
  Float_t GetV(Int_t) const; // V position at specified plane

  Float_t GetZ(Int_t) const; // Z position at specified plane

  virtual void ClearUVT(); // clears all STL maps
  Bool_t IsTPosValid(Int_t) const; // returns true if U and V positions have been set for this plane  
  Float_t GetMinU(Int_t plane, Double_t minPE=0) const;
  Float_t GetMaxU(Int_t plane, Double_t minPE=0) const;
  Float_t GetMinV(Int_t plane, Double_t minPE=0) const;
  Float_t GetMaxV(Int_t plane, Double_t minPE=0) const;
  Int_t GetNStrips(Int_t);

  void SetT(Int_t,StripEnd::StripEnd_t,Double_t);
  Double_t GetT(Int_t) const; // time at specified plane
  Double_t GetT(Int_t,StripEnd::StripEnd_t) const; // time at specified plane
  Double_t GetT(StripEnd::StripEnd_t,Int_t) const; // time at specified plane

  virtual Bool_t BelongsWithTrack(CandTrackHandle * trk, 
                                  AlgConfig & ac, 
                                  const VldContext * vldcptr, 
                                  Double_t tolTPos2, Double_t tolZPos, Double_t tolTime);
  virtual Bool_t BelongsWithShower(InoShowerCand * shw, 
                                   AlgConfig & ac, 
                                   const VldContext * vldcptr, 
                                   Double_t tolTPos2, Double_t tolZPos, Double_t tolTime);
  void AddCluster(InoCluster *);
  void RemoveCluster(InoCluster *);
  const InoCluster *GetCluster(Int_t) const;
  Int_t GetLastCluster() const;
  bool IsContained();
  virtual Bool_t IsUnphysical(Float_t xtalkFrac=0.660,
                              Float_t xtalkCut=2.0);

  void SetEnergy(Double_t, InoShowerCand::EShowerType=kCC);
  Double_t GetEnergy(InoShowerCand::EShowerType=kCC) const;
  void CalibrateEnergy(CandTrackHandle *associatedtrk,   AlgConfig & ac);
  static NavKey KeyFromSlice(const InoShowerCand *);

  Double_t DeDx(Double_t);

  private:

  Double_t fEnergy;
  Double_t fEnergy_wtCC;
  Double_t fEnergy_CC;
  Double_t fEnergy_wtNC;
  Double_t fEnergy_NC;
  Double_t fEnergy_EM;
  TObjArray fClusterList;         // Components owned from CandShower v2
  
  mutable map<Int_t,Float_t> fUPos;
  mutable map<Int_t,Float_t> fVPos;
  mutable map<Int_t,Double_t> fTime[2];
  */

};

#endif                                              // INOSHOWERCAND_H
