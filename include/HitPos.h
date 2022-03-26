
// ********************************************************************
// * License and Disclaimer                                           *
// *This  code  implementation is the result of  the  scientific and  *
// *technical work of the India Based Neutrino Observatory (INO)      *
// *group of Prof.Dr.Naba Mondal at the Tata Institute of Fundamental *
// *Research (TIFR),Mumbai.                                           *
// *                                                                  *
// *Neither the authors of this software system, nor their employing  *
// *institutes,nor the agencies providing financial support for this  *
// *work  make  any representation or  warranty, express or implied,  *
// *regarding  this  software system or assume any liability for its  *
// * use.                                                             *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications. Copyrights belong  *
// * to the respective author(s).                                     *
// *  Author (s) : Deepak Samuel     contact: samuel@tifr.res.in      *
// *  Last Edited on: 5.04.2009                                       *
// ********************************************************************
#ifndef ROOT_HitPos
#define ROOT_HitPos

#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"

class HitPos;

class HitPos : public TObject {
public:
   int TrackType;// Track Type: -1: hits, -101 : Primary hits -2: clulster, -3: triplet, -4: track,-44: track /fitter?
   int HitNum;// Hit Number
   int PrimHitNum; //Primary Number 
   int CluNum;// Cluster Number
   int TriNum;// Triplet Number
   int FitUpNum;// Fit  Number
   int FitDownNum;// Fit  Number
   int FindNum; // Finder Number
   int ShowerHitNum; //Shower Hit Number
   int ParCode;// Particle Code
   //int FSetNum;// Finder Track Set Number
   int Fup;// Finder set up  even 0 odd 1
   float XX, YY, ZZ; // X, Y, Z positions in case of tracks , vertex in case of particle info
   float pmag, pt, pp; // in case of particle info, pmagnitude, ptheta and pphi
  // int NEveNum;
   HitPos();
   //   HitPos(const HitPos& orig); //GMA14
   HitPos(Float_t random);
   virtual ~HitPos();
   ClassDef(HitPos,2)  //A track segment
};
#endif
