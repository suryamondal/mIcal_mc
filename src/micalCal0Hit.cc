
#include "micalCal0Hit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
using namespace std;

G4Allocator<micalcal0Hit> micalcal0HitAllocator;

micalcal0Hit::micalcal0Hit() {
  pdgid=-25;
  edep = 0;
  toff = 1000000;
  
}

micalcal0Hit::~micalcal0Hit()
{;}

micalcal0Hit::micalcal0Hit(const micalcal0Hit &right)
  : G4VHit()
{
  pdgid  = right.pdgid;
  edep = right.edep;
  pos = right.pos;
  mom = right.mom;
  toff = right.toff;
  HitId = right.HitId;
}

const micalcal0Hit& micalcal0Hit::operator=(const micalcal0Hit &right) {
  pdgid = right.pdgid;
  edep = right.edep;
  pos = right.pos;
  mom = right.mom;
  toff = right.toff;
  HitId = right.HitId;

  return *this;
}

G4int micalcal0Hit::operator==(const micalcal0Hit &right) const {
  return (this==&right) ? 1 : 0;
}

void micalcal0Hit::Draw() {
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(0.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

void micalcal0Hit::Print() {
  cout<<"hit "<<HitId<<" "<<pos<<endl;
}


