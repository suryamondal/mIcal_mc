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


#ifndef ROOT_Hits
#define ROOT_Hits



#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TH1.h"
#include "TBits.h"
#include "TMath.h"
#include "HitPos.h"

class HitPos;
class Hits;
class Hits : public TObject {
 public:
  int ENum; // Event Number
  int NEle;// Number of HitPos Objects (i.e NHits+NClus+NTrips+Ntracs)
  int NHits; // Number of Hits not sets!
  int NPrimHits; // Number of Primary hits sets //max 10  added on Dec7
  int NClus; // Number of Clusters sets
  int NTrips;// Number of Triplets sets
  int NFitUp;// Number of Fit up sets
  int NFitDown;// Number of Fit up sets
  int NFinders;// Number of Finder sets   
 // int NFTracs;// Number of Finder Tracks -4
  int NParticles;// Number of Particles
  int NRecTracks; //Number of Reconstructed Tracks 
  int NShowerHits; //Number of Shower Hits
  TRef fLastHit;         //reference pointer to last Hit/Cluster/Trip or track
  TClonesArray  *fHits;            //->array with all Hits
  static TClonesArray *fgHits;
  Bool_t         fIsValid;  
       //
public:
   Hits();
   virtual ~Hits();
 HitPos    *AddHits(Float_t random, Float_t ptmin=1);
 void ClearTracks();// clears previous track objects..
 TClonesArray *GetHits() const {return fHits;}
     ClassDef(Hits,1)  //Event structure 
 };

#endif
 




