#ifndef INOSTRIP_H
#define INOSTRIP_H
//CandStripHandle

class InoStrip
{

public:
  InoStrip();
  InoStrip(InoStrip* cd);
  ~InoStrip();
  InoStrip *DupHandle() const;

  void Trace(const char *c = "") const;

  int    GetPlaneView() const { return fView;};
  int    GetStrip() const { return fStrip;};
  double GetXYPos() const {return fXYPos;};

  int    GetPlane() const { return fPlane; };
  double GetZPos() const {return fZPos;};
  int GetTrueTime() const {return iTrueTime;};
  int GetSmrTime() const {return iSmrTime;};
  //  double GetTimeY() const {return fTimeY;};
  double GetPulse() const {return fPulse;}; 

  void   SetPlaneView( int f) {fView=f;};
  void   SetStrip(int f) { fStrip=f;};
  void   SetXYPos(double f) {fXYPos=f;};

  void   SetPlane( int f) {fPlane=f; };
  void   SetZPos( double f) { fZPos=f;};
  void   SetTrueTime( int ia) { iTrueTime=ia;};
  void   SetSmrTime( int ia) { iSmrTime=ia;};
  //  void   SetTimeY( double f) { fTimeY=f;};
  void   SetPulse( double f) { fPulse=f;}; 
  void   AddPulse(double f) {fPulse +=f;};
  void   SetRPCmod( int i) { fRPCmod = i;}; 
  int    GetRPCmod() {return fRPCmod;};

  void SetMomentum(double f) {fMomentum = f;};
  void SetTheta(double f) {fTheta = f;};
  void SetPhi(double f) {fPhi = f;};

  double GetMomentum() {return fMomentum ;};
  double GetTheta() {return fTheta ;};
  double GetPhi() {return fPhi;};

  void SetGenPosX(double f) {fXgen = f;}
  void SetGenPosY(double f) {fYgen = f;}
  void SetGenPosZ(double f) {fZgen = f;}


  double GetGenPosX() {return fXgen;}
  double GetGenPosY() {return fYgen;}
  double GetGenPosZ() {return fZgen;}

  void   SetId(int id) {fId= id;}
  int    GetId() { return fId;};

  void   SetpdgId(int id) {pdgid= id;}
  int    GetpdgId() { return pdgid;};

  int    fRPCmod; //RPC module ID
  int    fView;   // 0/1 for X/Y-axis
  int    fStrip;   //StripID
  double fXYPos;  //< Transverse position of strip (m).
  int    fPlane;  //Z-Plane, start from -ve Z of INOICAL
  double fZPos; //< Z position of strip (m)
  int    pdgid;  //pdgId of the particle which produce this hit
  int    iTrueTime; //time in ns
  int    iSmrTime; //time in ns
  //  double fTimeY; //time in ns
  double fPulse; //Pulse height in fC; //For INOCAL we are not using it
  //for MC information
  double fMomentum; //Momentum of track which is behind this hit
  double fTheta;  //Theta
  double fPhi;    // Phi

  double fXgen; //Generated track x-coordinate
  double fYgen; //Generated track y-coordinate
  double fZgen; //Generated track z-coordinate

  int    fId; // id, will replace all other individual ids

};
#endif                                              // INOSTRIP_H
