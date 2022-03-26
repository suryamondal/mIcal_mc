#include <vect_manager.h>
#include <fstream>
#include "TSpline.h"

class InoMuRange
{
 private:
  
 public:
  InoMuRange();
  virtual ~InoMuRange();
  double MaterialMuRange(double Z, double P);
  double FirstDerivative(double RangeZ, double Z);
  double SecndDerivative(double RangeZ, double Z);
  double ThirdDerivative(double RangeZ, double Z);
  double FourthDerivative(double RangeZ, double Z);
  
  TSpline3 *Spline1D_CGr;		TSpline3 *RGrInP;
  TSpline3 *Spline1D_Alu;		TSpline3 *RAlInP;
  TSpline3 *Spline1D_Gls;		TSpline3 *RGlInP;
  TSpline3 *Spline1D_Irn;		TSpline3 *RIrInP;
  TSpline3 *Spline1D_Cop;		TSpline3 *RCuInP;
};
