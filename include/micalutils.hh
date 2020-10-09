///////////////////////////////////////////////////////////////////////////////
// File: micalutils.hh
// Description: General utilities.
///////////////////////////////////////////////////////////////////////////////
#ifndef micalutils_hh
#define micalutils_hh 1

#include <iostream>
#include <fstream>
#include "globals.hh"


G4String operator+(const G4String&, const int);
G4String operator+(const G4String&, const double);
// "number " + i = "number i"

std::ifstream& readName(std::ifstream&, G4String&);
// It reads a name into G4String between quotes and skips lines 
// beginning with #. and if found *ENDDO returns.

std::ifstream& findDO(std::ifstream&, const G4String&);
// It reads until a *DO str is found.

std::ostream& tab(std::ostream&);
// It add a tab.

std::istream& jump(std::istream&);
// It ignores character until the end of line.

bool openGeomFile(std::ifstream& is, const G4String& pathname, 
		  const G4String& filename);
// It opens the geometry file, either locally (if it exists) or "remotely".


#endif
