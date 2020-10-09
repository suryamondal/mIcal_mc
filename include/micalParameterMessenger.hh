#ifndef ICAL0PARAMETERMESSENGER_H
#define ICAL0PARAMETERMESSENGER_H

#include "globals.hh"

class micalParameterMessenger
{
  public:
    micalParameterMessenger();
    virtual ~micalParameterMessenger();

    static micalParameterMessenger* AnPointer;

    //Getters
    G4String GetParameterLocation() {
      return parameterLocation;
    }

    G4String GetGeometryLocation() {
      return geometryLocation;
    }

    G4String GetStripInfoLocation() {
      return stripInfoLocation;
    }

    G4double GetnFeThickness() {
      return nFeThickness;
    }

    G4double GetnAirGap() {
      return nAirGap;
    }

    G4int GetnLayer() {
      return nLayer;
    }

    G4double GetXYstrwd() {
      return XYstrwd;
    }

    G4String GetFileVersion() {
      return fileVersion;
    }

    //Setters
    void SetParameterLocation(G4String value) {
      parameterLocation = value;
    }

    void SetGeometryLocation(G4String value) {
      geometryLocation = value;
    }

    void SetStripInfoLocation(G4String value) {
      stripInfoLocation = value;
    }

    void SetnFeThickness(G4double value) {
      nFeThickness = value;
    }

    void SetnAirGap(G4double value) {
      nAirGap = value;
    }

    void SetnLayer(G4int value) {
      nLayer = value;
    }
  
    void SetXYstrwd(G4double value) {
      XYstrwd = value;
    }

    void SetFileVersion(G4String value) {
      fileVersion = value;
    }

  protected:

private:
  G4String parameterLocation;
  G4String geometryLocation;
  G4String stripInfoLocation;
  G4String fileVersion;
  G4int nLayer;
  G4double nFeThickness;
  G4double nAirGap;
  G4double XYstrwd;
};

#endif // ICAL0PARAMETERMESSENGER_H
