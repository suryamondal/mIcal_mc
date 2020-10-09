
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

#include "Hits.h"

TClonesArray *Hits::fgHits = 0;

Hits::Hits()  {
  // Create an Hit object.
  if (!fgHits) fgHits = new TClonesArray("HitPos", 1000);
  fHits = fgHits;
  NEle=0;
  //  Ntracks = 0;
  //Nevents++; 
  //int Event::NEvents=0;
}

//______________________________________________________________________________
Hits::~Hits()
{
   
}

//______________________________________________________________________________
HitPos * Hits::AddHits(Float_t random, Float_t ptmin) {
  // Add a new track to the list of tracks for this event.
  ptmin=0;// variable not used in INO
  TClonesArray &hitp = *fHits;
  HitPos *hitpos = new(hitp[NEle++]) HitPos(random);
  //Save reference to last Track in the collection of Tracks
  fLastHit= hitpos;//
  return hitpos;
}


void Hits::ClearTracks() {
  if (fHits) {
    fHits->Clear("C"); // will also call Track::Clear
    NEle = 0; // !!! reset counter !!!
  }
}



