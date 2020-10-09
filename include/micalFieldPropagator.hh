#ifndef MICALFIELDPROPAGATOR_HH
#define MICALFIELDPROPAGATOR_HH
#include "micalDetectorParameterDef.hh"
#include "MultiSimAnalysis.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "vect_manager.h"
#include "TH1.h"
#include "TFile.h"
#include "TH2.h"
#include <fstream>

class micalFieldPropagator {
public:
  micalFieldPropagator();
  ~micalFieldPropagator();
  void ElectroMagneticField(const double Point[3], double& B1, double& B2, int ftype);
  void PrintFieldMap();
  void F2int(int* ag_f2i, double& bx_f2i, double& by_f2i);
  double bilinearInterpolation(double* f_bli,double* arg_bli,int* ag_bli);
public:
 static micalFieldPropagator* FdPointer;
  // Access functions
 
  micalDetectorParameterDef *paradef;
  MultiSimAnalysis *pAnalysis;

protected:
  // Find the global Field Manager
  int temp;
private:
  
  TFile* pMagFile;
  int slotsize;
  TH2D* fieldxin;
  TH2D* fieldyin;

  TGeoManager* geometryIn;
  
  double irlayZ; //nFeThickness/2;
  double rpclayZ; //nAirGap/2;
  double ironrpcZ; // nFeThickness + nAirGap;
  int nLayer;
  int nIRLayer;
  double parino[3];
  double gapino;
  double parirlay[3];
  double IRONLayerPosZ[300];
  double RPCLayerPosZ[300];
  double StackPosInRoom[3];
  double INOroomPos[3];
  double PosLayerTop;
  double PosLayerBot;
  double StackShift[3];
  
};

#endif // MICALFIELDPROPAGATOR_HH
