#include <fstream>
#include <iostream>
#include <InoMuRange.h>

InoMuRange::InoMuRange() {
  int M = 141;
  int ijk = 0;
  
  double muP = 0.0;		double MuMom[141];
  //  double muIP= 0.0;		
  double MuInP[141];
  double muGr= 0.0;		double CGrBB[141];
  double muAl= 0.0;		double AluBB[141];
  double muGl= 0.0;		double GlsBB[141];
  double muFe= 0.0;		double IrnBB[141];
  double muCu= 0.0;		double CopBB[141];
  
  ifstream RangeMomentumLists;
  RangeMomentumLists.open("Range.txt");
  // RangeMomentumLists.open("/products/GEANT4.10/ICALDOCS/Range.txt");
  while(!RangeMomentumLists.eof()) {
    if (ijk > 140) break;
    RangeMomentumLists>>muP>>muFe>>muCu>>muAl>>muGr>>muGl;
    //cout<<"muP "<<muP<<"     muFe "<<muFe<<"     ijk "<<ijk<<endl;
    
    //muP		MeV/c
    //MuMom		GeV/c
    //MuInP		GeV/c ]^-1
    
    MuMom[ijk] = 1.e-3*muP;
    MuInP[ijk] = 1.0/(1.e-3*muP);
    CGrBB[ijk] = muGr;
    AluBB[ijk] = muAl;
    GlsBB[ijk] = muGl;
    IrnBB[ijk] = muFe;
    CopBB[ijk] = muCu;
    
    ijk++;
  }
  
  Spline1D_CGr = new TSpline3("Cubic Spline", MuMom, CGrBB, M, "b2e2", 0, 0);
  Spline1D_Alu = new TSpline3("Cubic Spline", MuMom, AluBB, M, "b2e2", 0, 0);
  Spline1D_Gls = new TSpline3("Cubic Spline", MuMom, GlsBB, M, "b2e2", 0, 0);
  Spline1D_Irn = new TSpline3("Cubic Spline", MuMom, IrnBB, M, "b2e2", 0, 0);
  Spline1D_Cop = new TSpline3("Cubic Spline", MuMom, CopBB, M, "b2e2", 0, 0);
  
  RGrInP = new TSpline3("Cubic Spline", CGrBB, MuInP, M, "b2e2", 0, 0);
  RAlInP = new TSpline3("Cubic Spline", AluBB, MuInP, M, "b2e2", 0, 0);
  RGlInP = new TSpline3("Cubic Spline", GlsBB, MuInP, M, "b2e2", 0, 0);
  RIrInP = new TSpline3("Cubic Spline", IrnBB, MuInP, M, "b2e2", 0, 0);
  RCuInP = new TSpline3("Cubic Spline", CopBB, MuInP, M, "b2e2", 0, 0);
}

InoMuRange::~InoMuRange() {
  delete Spline1D_CGr;		delete RGrInP;
  delete Spline1D_Alu;		delete RAlInP;
  delete Spline1D_Gls;		delete RGlInP;
  delete Spline1D_Irn;		delete RIrInP;
  delete Spline1D_Cop;		delete RCuInP;
}

double InoMuRange::MaterialMuRange(double Z, double P) {
  double MuRangeZ= 0.0;
  
  if (Z == 6.0)	{
    MuRangeZ = Spline1D_CGr->Eval(P);
  } else if (Z == 10.8046) {
    MuRangeZ = Spline1D_Gls->Eval(P);
  }  else if (Z == 13.0) {
    MuRangeZ = Spline1D_Alu->Eval(P);
  } else if (Z == 26.0) {
    MuRangeZ = Spline1D_Irn->Eval(P);
  } else if (Z == 29.0) {
    MuRangeZ = Spline1D_Cop->Eval(P);
  } else {
    cout<<"information not available..."<<endl;
  }
  return MuRangeZ;
}

double InoMuRange::FirstDerivative(double R, double Z) {
  double d1fdr1	= 0.0;
  double dr		= 1.0;
  if (Z == 6.0) {
    d1fdr1 = -(RGrInP->Eval(R+4*dr)-40*RGrInP->Eval(R+2*dr)+256*RGrInP->Eval(R+dr)-256*RGrInP->Eval(R-dr)+40*RGrInP->Eval(R-2*dr)-RGrInP->Eval(R-4*dr))/(360*dr);
  } else if (Z == 10.8046) {
    d1fdr1 = -(RGlInP->Eval(R+4*dr)-40*RGrInP->Eval(R+2*dr)+256*RGrInP->Eval(R+dr)-256*RGrInP->Eval(R-dr)+40*RGrInP->Eval(R-2*dr)-RGrInP->Eval(R-4*dr))/(360*dr);
  } else if (Z == 13.0) {
    d1fdr1 = -(RAlInP->Eval(R+4*dr)-40*RAlInP->Eval(R+2*dr)+256*RAlInP->Eval(R+dr)-256*RAlInP->Eval(R-dr)+40*RAlInP->Eval(R-2*dr)-RAlInP->Eval(R-4*dr))/(360*dr);
  } else if (Z == 26.0) {
    d1fdr1 = -(RIrInP->Eval(R+4*dr)-40*RIrInP->Eval(R+2*dr)+256*RIrInP->Eval(R+dr)-256*RIrInP->Eval(R-dr)+40*RIrInP->Eval(R-2*dr)-RIrInP->Eval(R-4*dr))/(360*dr);
  } else if (Z == 29.0) {
    d1fdr1 = -(RCuInP->Eval(R+4*dr)-40*RCuInP->Eval(R+2*dr)+256*RCuInP->Eval(R+dr)-256*RCuInP->Eval(R-dr)+40*RCuInP->Eval(R-2*dr)-RCuInP->Eval(R-4*dr))/(360*dr);
  }
  return d1fdr1;
}

double InoMuRange::SecndDerivative(double R, double Z) {
  double d2fdr2	= 0.0;
  double dr		= 1.0;
  if (Z == 6.0)	{ //GMA14, what about less than 6
    d2fdr2 = (-RGrInP->Eval(R+2*dr) + 16*RGrInP->Eval(R+dr) - 30*RGrInP->Eval(R) + 16*RGrInP->Eval(R-dr) - RGrInP->Eval(R-2*dr))/(12*dr*dr);
  } else if (Z == 10.8046) {
    d2fdr2 = (-RGlInP->Eval(R+2*dr) + 16*RGlInP->Eval(R+dr) - 30*RGlInP->Eval(R) + 16*RGlInP->Eval(R-dr) - RGlInP->Eval(R-2*dr))/(12*dr*dr);
  } else if (Z == 13.0)	{
    d2fdr2 = (-RAlInP->Eval(R+2*dr) + 16*RAlInP->Eval(R+dr) - 30*RAlInP->Eval(R) + 16*RAlInP->Eval(R-dr) - RAlInP->Eval(R-2*dr))/(12*dr*dr);
  } else if (Z == 26.0)	{
    d2fdr2 = (-RIrInP->Eval(R+2*dr) + 16*RIrInP->Eval(R+dr) - 30*RIrInP->Eval(R) + 16*RIrInP->Eval(R-dr) - RIrInP->Eval(R-2*dr))/(12*dr*dr);
  } else if (Z == 29.0)	{
    d2fdr2 = (-RCuInP->Eval(R+2*dr) + 16*RCuInP->Eval(R+dr) - 30*RCuInP->Eval(R) + 16*RCuInP->Eval(R-dr) - RCuInP->Eval(R-2*dr))/(12*dr*dr);
  }
  return d2fdr2;
}

double InoMuRange::ThirdDerivative(double R, double Z) {
  double d3fdr3	= 0.0;
  double dr		= 1.0;
  if (Z == 6.0)	{
    d3fdr3 = (-RGrInP->Eval(R+3*dr) + 8*RGrInP->Eval(R+2*dr) - 13*RGrInP->Eval(R+dr) + 13*RGrInP->Eval(R-dr) - 8*RGrInP->Eval(R-2*dr) + RGrInP->Eval(R-3*dr))/(8*dr*dr*dr);
  } else if (Z == 10.8046) {
    d3fdr3 = (-RGlInP->Eval(R+3*dr) + 8*RGlInP->Eval(R+2*dr) - 13*RGlInP->Eval(R+dr) + 13*RGlInP->Eval(R-dr) - 8*RGlInP->Eval(R-2*dr) + RGlInP->Eval(R-3*dr))/(8*dr*dr*dr);
  } else if (Z == 13.0)	{
    d3fdr3 = (-RAlInP->Eval(R+3*dr) + 8*RAlInP->Eval(R+2*dr) - 13*RAlInP->Eval(R+dr) + 13*RAlInP->Eval(R-dr) - 8*RAlInP->Eval(R-2*dr) + RAlInP->Eval(R-3*dr))/(8*dr*dr*dr);
  } else if (Z == 26.0) {
    d3fdr3 = (-RIrInP->Eval(R+3*dr) + 8*RIrInP->Eval(R+2*dr) - 13*RIrInP->Eval(R+dr) + 13*RIrInP->Eval(R-dr) - 8*RIrInP->Eval(R-2*dr) + RIrInP->Eval(R-3*dr))/(8*dr*dr*dr);
  } else if (Z == 29.0)	{
    d3fdr3 = (-RCuInP->Eval(R+3*dr) + 8*RCuInP->Eval(R+2*dr) - 13*RCuInP->Eval(R+dr) + 13*RCuInP->Eval(R-dr) - 8*RCuInP->Eval(R-2*dr) + RCuInP->Eval(R-3*dr))/(8*dr*dr*dr);
  }
  return d3fdr3;
}

double InoMuRange::FourthDerivative(double R, double Z) {
  double d4fdr4	= 0.0;
  double dr		= 1.0;
  if (Z == 6.0) {
    d4fdr4 = (-RGrInP->Eval(R+3*dr) + 12*RGrInP->Eval(R+2*dr) -39*RGrInP->Eval(R+dr) + 56*RGrInP->Eval(R) - 39*RGrInP->Eval(R-dr) + 12*RGrInP->Eval(R-2*dr) - RGrInP->Eval(R-3*dr))/(6*dr*dr*dr*dr);
  } else if (Z == 10.8046) {
    d4fdr4 = (-RGlInP->Eval(R+3*dr) + 12*RGlInP->Eval(R+2*dr) -39*RGlInP->Eval(R+dr) + 56*RGlInP->Eval(R) - 39*RGlInP->Eval(R-dr) + 12*RGlInP->Eval(R-2*dr) - RGlInP->Eval(R-3*dr))/(6*dr*dr*dr*dr);
  } else if (Z == 13.0)	{
    d4fdr4 = (-RAlInP->Eval(R+3*dr) + 12*RAlInP->Eval(R+2*dr) -39*RAlInP->Eval(R+dr) + 56*RAlInP->Eval(R) - 39*RAlInP->Eval(R-dr) + 12*RAlInP->Eval(R-2*dr) - RAlInP->Eval(R-3*dr))/(6*dr*dr*dr*dr);
  } else if (Z == 26.0)	{
    d4fdr4 = (-RIrInP->Eval(R+3*dr) + 12*RIrInP->Eval(R+2*dr) -39*RIrInP->Eval(R+dr) + 56*RIrInP->Eval(R) - 39*RIrInP->Eval(R-dr) + 12*RIrInP->Eval(R-2*dr) - RIrInP->Eval(R-3*dr))/(6*dr*dr*dr*dr);
  } else if (Z == 29.0)	{
    d4fdr4 = (-RCuInP->Eval(R+3*dr) + 12*RCuInP->Eval(R+2*dr) -39*RCuInP->Eval(R+dr) + 56*RCuInP->Eval(R) - 39*RCuInP->Eval(R-dr) + 12*RCuInP->Eval(R-2*dr) - RCuInP->Eval(R-3*dr))/(6*dr*dr*dr*dr);
  }
  
  return d4fdr4;
}
