//

#ifndef micalcal0Hit_h
#define micalcal0Hit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class micalcal0Hit : public G4VHit
{
  public:

      micalcal0Hit();
      ~micalcal0Hit();
      micalcal0Hit(const micalcal0Hit &right);
      const micalcal0Hit& operator=(const micalcal0Hit &right);
      G4int operator==(const micalcal0Hit &right) const;

      inline void *operator new(size_t);
      inline void operator delete(void *aHit);

      void Draw();
      void Print();

  private:
      G4double edep;
      G4int    pdgid;  //Particle ID
      G4ThreeVector pos;
      G4double localx;
      G4double localy;
      G4ThreeVector mom; //Momentum of track at earliest energy deposite, has meaning only for muon track, not usefull at all for hadronic shower
      G4double toff;
      G4double tofx;
      G4double tofy;
      unsigned long    HitId;
  public:
      inline void SetpdgId(G4int id)
      { pdgid = id; }
      inline void SetEdep(G4double de)
      { edep = de; }
      inline void AddEdep(G4double de) 
      { edep +=de; }
      inline G4double GetEdep()
      { return edep; }
      inline void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      inline G4ThreeVector GetPos()
      { return pos; }
      inline void SetMom(G4ThreeVector xyz)
      { mom = xyz; }
      inline G4ThreeVector GetMom()
      { return mom; }
      inline void SetTime(G4double tf)
      { toff = tf; }
      inline G4double GetTime()
      { return toff; }

      inline void SetLocalXPos(G4double xyz)
     { localx = xyz; }
      inline G4double GetLocalXPos()
      { return localx; }
  
      inline void SetLocalYPos(G4double xyz)
      { localy = xyz; }
      inline G4double GetLocalYPos()
      { return localy; }
  

  //      inline void SetTimeX(G4double tf)
  //      { tofx = tf; }
  //      inline G4double GetTimeX()
  //      { return tofx; }
  //      inline void SetTimeY(G4double tf)
  //      { tofy = tf; }
  //     inline G4double GetTimeY()
  //      { return tofy; }
      inline void SetHitId (unsigned long id)
      { HitId = id; }
      inline unsigned long GetHitId()
      { return HitId; }
      inline G4int GetpdgId()
      { return pdgid; }
};

typedef G4THitsCollection<micalcal0Hit> micalcal0HitsCollection;

extern G4Allocator<micalcal0Hit> micalcal0HitAllocator;

inline void* micalcal0Hit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) micalcal0HitAllocator.MallocSingle();
  return aHit;
}

inline void micalcal0Hit::operator delete(void *aHit)
{
  micalcal0HitAllocator.FreeSingle((micalcal0Hit*) aHit);
}

#endif
