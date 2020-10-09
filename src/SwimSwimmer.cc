//
// Swimming through a particle, forward or backward
//
#include "SwimSwimmer.h"
#include "SwimParticle.h"

#include "TMath.h"
#include <string>
#include <cassert>
#include <iostream>
using namespace std;
//......................................................................

SwimSwimmer::SwimSwimmer(double dist, double halfgap ) : //const VldContext& vldc) :
  //  fMagField(0),
  //  fStepMax(0.0852), //254*Munits::m),
  fStepMin(1.0e-6), //(*Munits::m),
  fAcc(1.0e-4),
  fNmaxStep(1000),
  fStepSize(0.005),
  fSPI(-1)
{
  paradef = micalDetectorParameterDef::AnPointer;//AAR:these variables added to include variable airgap
  fDistance = dist;
  fHalfLayerThickness = halfgap;
  fHalfAirGap=paradef->GetParlay(2)/1000;
  fStepMax = 2*halfgap;
}


SwimSwimmer::SwimSwimmer(int pln, double dist, double halfgap) : //const VldContext& vldc) :
  //  fMagField(0),
  fStepMin(1.0e-6), //(*Munits::m),
  fAcc(1.0e-4),
  fNmaxStep(1000),
  fStepSize(0.005),        // AAR0.02->0.005
  fSPI(-1)
{
  //paradef = micalDetectorParameterDef::AnPointer;//AAR:these variables added to include variable airgap
  //fDistance = dist;
  //fHalfLayerThickness = halfgap;
  //fHalfAirGap=paradef->GetParirlay(2)/1000;
  //fStepMax = 2*halfgap;
  
  
  fDistance = dist;
  fHalfLayerThickness = halfgap;
  fStepMax = dist; // 2*halfgap;
  lastCrossingShift = TVector3(0,0,0);
  startposXYZ = TVector3(0,0,0);
  Plane = pln;
}


SwimSwimmer::~SwimSwimmer() {
  //  this->DeleteBfield();
  //  this->DeleteSwimGeo();
  //  this->DeleteStepper();
  
  //  if (fStepData) {
  //    delete fStepData;
  //    fStepData = 0;
  //  }
}

//......................................................................

bool SwimSwimmer::SwimForward(SwimParticle& particle, double& b_ave) { //, (SwimCondition& c)
  
  // Arrow of time, true for dt>0, false for dt<0
  //  fStepData->SetIsForward(true);
  fIsForward = true;
  nbfield=0;
  b_ave=0;
  b_ave=Swim(particle); //, c) 
  
  if(b_ave!=0){ return true;
  } else {return false;}
}

//......................................................................

bool SwimSwimmer::SwimBackward(SwimParticle& particle, double& b_ave) { // , (SwimCondition& c)
  // Arrow of time, true for dt>0, false for dt<0
  //  fStepData->SetIsForward(false);
  
  fIsForward = false;
  nbfield=0;
  b_ave=0;
  
  b_ave=Swim(particle); //, c);
  
  if(b_ave!=0) { return true;
  } else { return false;}
}

//......................................................................

double SwimSwimmer::Swim(SwimParticle& particle) { // (SwimCondition& c)

  pFieldMap = micalFieldPropagator::FdPointer;
  
  //  SwimGeo::SwimMaterial_t material;
  double              fStep;
  double              distToNextPlane = 0;
  double              pThres = 0.05; //*Munits::GeV;
  //  bool                satisfied = true; //false;
  bool                zDirInitial = false, zDir;   // z-momentum direction
  double     fbave=0;
  double  fnbfield=0;
  
  
  // Check if particle has a momentum < 50MeV
  if (particle.GetMomentumModulus()<pThres) {
    //    satisfied = c.Satisfied(particle);
    //    return satisfied;
    return false;
  }
  
  // initial layer is undefined (gets set in GetSwimMaterial)
  //  fStepData->SetSPI(-1);
  SetSPI(-1);
  
  for (int ij=0; ij<fNmaxStep; ++ij) {
    const TVector3 xyz = particle.GetPosition();
    const TVector3 direction = particle.GetDirection();
    Double_t momZ= particle.GetMomentum().Z();
    
    // particle traveling in position/negative z-direction
    if ((fIsForward && momZ>0.0) ||
	(!(fIsForward) && momZ<0.0)) {
      zDir = true;
    } else {
      zDir = false;
    }
    
    if (ij==0) {  zDirInitial = zDir;}
    
    //<< "ZDir" << zDir<<" " << fIsForward <<" "<< momZ<< endl;
    if (zDir!=zDirInitial) { return false ;
      break; }         // particle changed swimming direction. Exit for loop
    
    // Calculate the distance to the next plane
    // GMA Stepsize along Z-axis
    if (distToNextPlane+fStepSize <fDistance) {
      fStep = fStepSize;
      distToNextPlane += fStepSize; // fSwimGeo->DistToNextPlane(xyz, direction);
    } else {
      fStep = fDistance-distToNextPlane;
      distToNextPlane = fDistance;
    }
    //cout<<" fStep "<< fStep << " "<<ij<<endl;
    
    //    // Far detector limit: 11.6 meters, which always works for near detector
    if (distToNextPlane==0.0 || distToNextPlane>=11.6 || distToNextPlane >fDistance) // *Munits::m)
      break;           // Outside the detector. Exit for loop
    //AAR: distToNextPlane>11.6 is this condition required: though it will never be satisfied in out case
    //======================================================================
    // Step Action: StepOnce, dEdx
    //======================================================================
    
    //    fStepper->Action(particle,fStepData);
    //  SwimGeo::SwimMaterial_t material = stepData->GetSwimMaterial();
    if (particle.GetMomentumModulus()!=0.0) {
      // Used to calculate path, range traveled in this step
      double density; //  = 7.847; //SwimGeo::GetSwimMaterialDensity(material);
      TVector3 startpos, startmom;
      startpos = particle.GetPosition();
      startmom = particle.GetMomentum();
      
      //GMA Change momentum too (convert old fortran code to C++)
      // 1. Transform co-ordisnate system such that magnetic field is along Z' axis
      // 2. Get distance to the crossing point of heliz and plane
      // 3. Get the track paramters at the crossing point
      // 4. Return back to INO co-ordinate system
      
      double pos[3], newpos[3];
      double mom[3], newmom[3];
      double dir[3];
      double projpos[3], newprojpos[3];
      double track[6], trackout[6];
      
      pos[0] = startpos.X();
      pos[1] = startpos.Y();
      pos[2] = startpos.Z();
      
      int ilay = (int)fabs((distToNextPlane)/(2*fHalfLayerThickness));
      double localZ = fabs(distToNextPlane - 2*ilay*fHalfLayerThickness) -fHalfLayerThickness;
      
      double B[6];
      double Bx,By;
      double buvz[3];
      double pos1[3];
      pos1[0] =pos[0]*1000;
      pos1[1] =pos[1]*1000;
      pos1[2] =pos[2]*1000;
      
      if (abs(localZ) < fHalfLayerThickness-fHalfAirGap){ // 125) { //AAR:these variables added to include variable airgap
	pFieldMap->ElectroMagneticField( pos1, Bx,By,1);
	B[0] =Bx*1000;B[1]= By*1000; B[2]=0;
	density = 7.847; //iron
      } else if (abs(localZ) < fHalfLayerThickness-0.004) {
	B[0]=0;B[1]=0; B[2]=0;
	density = 0.00001; // //Foam (neglecting strips)
      } else if (abs(localZ) < fHalfLayerThickness-0.003) {
	B[0]=0;B[1]=0; B[2]=0;
	density = 1.09;// G10;
      } else if (abs(localZ) < fHalfLayerThickness-0.001) {
	B[0]=0;B[1]=0; B[2]=0;
	density = 2.32;// quartz;
      } else {
	B[0]=0;B[1]=0; B[2]=0;
	density = 0.00137; //RPC gas
      }
      if(TMath::Sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2])!=0){ 
	fbave += TMath::Sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);
	fnbfield++;
      }
      
      buvz[0] = 0.3*B[0]; buvz[1] = 0.3*B[1]; buvz[2]= 0.3*B[2];
      
      mom[0] = startmom.X();
      mom[1] = startmom.Y();
      mom[2] = startmom.Z();
      
      double mommag = startmom.Mag();
      //cout<<"mommag = "<<mommag<<endl;
      
      for (int jk=0; jk<3; jk++) {dir[jk] = mom[jk]/mommag;/*cout <<mom[jk];*/}
      //Projected position;
      double dist = fabs(fStep/dir[2]); //startmom.CosTheta(); //GMA does it need absolute value ?
      
      
      // If material is Air or Scint
      // straight line approximation
      
      //double buvz[3] ={0.3, 0.3, 0.0}; //GMA-magnetic field
      //double buvz[3] ={0.45, 0.0, 0.0}; //GMA-magnetic field
      // GMA BX=By=1Tesla, Bz=0 and 0.3 factor comes from
      // the conversion factor, pt=0.3Br;
      
      double charge = particle.GetCharge();
      if (fIsForward) { charge *=-1.;}
      //GMA 06/02/2009 Chage -ve why ???????????????
      
      if (!fIsForward) { for (int jk=0; jk<3; jk++) {dir[jk] = -dir[jk];};}
      
      for (int jk=0; jk<3; jk++) {projpos[jk] = pos[jk] + dist*dir[jk];}
      
      //include bending in x-y plane
      double momxy = pow(max(1.e-12, mom[0]*mom[0] + mom[1]*mom[1]), 0.5);
      double momyz = pow(max(1.e-12, mom[1]*mom[1] + mom[2]*mom[2]), 0.5);
      double momzx = pow(max(1.e-12, mom[2]*mom[2] + mom[0]*mom[0]), 0.5);
      
      //GMA 06/02/2009 Looks like we swapped charge convention in the remaining coding
      projpos[0] -=0.5*charge*fStep*fStep * (buvz[2]/momxy - buvz[1]/momzx);
      projpos[1] -=0.5*charge*fStep*fStep * (buvz[0]/momyz - buvz[2]/momxy);
      
      // GMA 06/02/2009 Calculate average momentum in those two points (start and end)
      
      double rot[9]; //rot[3][3];
      // Calculate rotation matrix to rotate co-oridnate system such that
      // Z' is the direction of local Magnetic field
      anal_getnrot(buvz, rot);
      //      cout<<"buvz "<<buvz[0]<<" "<<buvz[1]<<" "<<buvz[2]<<" "<<rot[0]<<" "<<rot[1]<<" "<<rot[2]<<" "<<rot[3]<<" "<<rot[4]<<" "<<rot[5]<<" "<<rot[6]<<" "<<rot[7]<<" "<<rot[8]<<endl;
      
      //rotate local postion in new frame
      
      //      if ((posz>0.7&& posz<0.95) ||(posz >4.9 && posz<5.1)) cout<<"posmom1 "<<projpos[0]<<" "<< projpos[1]<<" "<< projpos[2]<<" "<<pos[0]<<" "<< pos[1]<<" "<< pos[2]<<" "<<mom[0]<<" "<< mom[1]<<" "<< mom[2]<<endl;
      anal_rotme(pos,rot,newpos);
      anal_rotme(mom,rot,newmom);
      anal_rotme(projpos,rot,newprojpos);
      
      //      if ((posz>0.7&& posz<0.95) ||(posz >4.9 && posz<5.1)) cout<<"posmom2 "<<newprojpos[0]<<" "<< newprojpos[1]<<" "<< newprojpos[2]<<" "<<newpos[0]<<" "<< newpos[1]<<" "<< newpos[2]<<" "<<newmom[0]<<" "<< newmom[1]<<" "<< newmom[2]<<endl;
      
      for (int jk=0; jk<3; jk++) {
	track[jk] = newpos[jk];
	
	if (fIsForward) {
	  track[jk+3] = newmom[jk];
	} else {
	  track[jk+3] = -newmom[jk];
	}
      }
      
      //Construct a plane whose co-oridnate is (0,0,pos[2]) and (0,0,mom[2]) and
      //rotate to new co-ordinate system

      double distep = 0; // total track length
      track_move_pt_align(buvz,track, charge, newprojpos, trackout, distep);//AAR: 1st argument added to pass the Field array.
      
      for (int jk=0; jk<3; jk++) {
	newpos[jk] = trackout[jk];
	
	if (fIsForward) {
	  newmom[jk] = trackout[jk+3];
	} else {
	  newmom[jk] = -trackout[jk+3];
	}
      }
      
      //return back to original co-ordinate system
      anal_rotmet(newpos,rot,pos);
      anal_rotmet(newmom,rot,mom);
      
      //Forward backward part is taken care off by the trace_track_planef ??

      //===============================================================
      // Set position, momentum of the particle
      //===============================================================
      
      TVector3 position(pos[0], pos[1], pos[2]);
      particle.SetPosition(position);
      
      //GMA use momentum dependent energy loss factor
      double tmpmom = mommag;
      double elos = (14.9+0.96*fabs(log(tmpmom*2.8))+0.033*tmpmom*(1.0-pow(tmpmom,-0.33)))*1e-2*1.1;
      
      //      cout <<"elos =============="<<tmpmom<<" "<<elos<<endl;
      if (fIsForward) {
	tmpmom -=density*distep*elos; // 0.175;
      } else {
	tmpmom +=density*distep*elos; // 0.175;
      }
      
      for (int jk=0; jk<3; jk++) {mom[jk] = tmpmom*mom[jk]/mommag;/*cout <<mom[jk];*/}
      //cout <<"new"<<mommag<< endl;
      TVector3 momentum(mom[0], mom[1], mom[2]);
      
      particle.SetMomentum(momentum);
      
      // Calculate path, range in this step, add to particle totals
      
      particle.AddS(distep);
      particle.AddRange(density*distep); // /(Munits::g/Munits::cm2));
      
      if (distToNextPlane >=fDistance) break;
    }
    
  }
  
  fbave /= fnbfield; 
  return fbave;
}

//......................................................................

bool SwimSwimmer::SwimForward(SwimParticle& particle, int & nextp, double& b_ave) {  //, (SwimCondition& c)

  // Arrow of time, true for dt>0, false for dt<0
  //  fStepData->SetIsForward(true);
  fIsForward = true;
  nbfield=0;
  b_ave=0;
  b_ave=Swim(particle, nextp); //, c);
  if(b_ave!=0) { return true;
  } else { return false;}
}

//......................................................................

bool SwimSwimmer::SwimBackward(SwimParticle& particle, int & nextp, double& b_ave) { // , (SwimCondition& c)

  // Arrow of time, true for dt>0, false for dt<0
  //  fStepData->SetIsForward(false);
  
  fIsForward = false;
  
  nbfield=0;
  b_ave=0;
  b_ave=Swim(particle, nextp); //, c);
  if(b_ave!=0) { return true;
  } else { return false;} 
}

//......................................................................

double SwimSwimmer::Swim(SwimParticle& particle, int& nextplane) {  // (SwimCondition& c)

  pFieldMap = micalFieldPropagator::FdPointer;
  
  double              fStep=.1;
  double              distToNextPlane = 0;
  double              pThres = 0.05; //*Munits::GeV;
  double     fbave=0; 
  double  fnbfield=0;
  //  bool                satisfied = true; //false;
  //  bool                zDirInitial = false, zDir;   // z-momentum direction
  
  // Check if particle has a momentum < 50MeV
  if (particle.GetMomentumModulus()<pThres) {
    return false;
  }
  
  SetSPI(-1);
  nextplane=-99;
  
  for (int ij=0; ij<fNmaxStep; ++ij) {
    
    const TVector3 xyz = particle.GetPosition();
    if (ij==0) { startposXYZ = xyz;}
    
    const TVector3 direction = particle.GetDirection();
    Double_t posZ= xyz.z();
    
    if (particle.GetMomentumModulus()<pThres/2.) return false;
    
    // GMA Stepsize along particle direction
    double invdircos = 10000.;
    if (direction[2] >0) {
      invdircos = 1./max(0.00001,direction[2]);
    } else if (direction[2] <0)  {
      invdircos = 1./min(-0.00001,direction[2]);
    }
    
    
    if (fabs(posZ+fStepSize*direction[2]-startposXYZ.z())<fDistance) { //*direction[2]) {
      fStep = fStepSize; // /max(0.001,direction[2]);
    } else {
      fStep = (fDistance-fabs(posZ-startposXYZ.z()))*fabs(invdircos); // /max(0.001,direction[2]);
      if (abs(fStep)<0.000001) {
	// <<"fDistance "<<fStep<<" "<< fDistance<<" "<< posZ<<" "<<posZ-startposXYZ.z()<<" "<<  invdircos<<endl;
      }
    }
    
    distToNextPlane +=fStep;
    
    //======================================================================
    // Step Action: StepOnce, dEdx
    //======================================================================
    
    if (particle.GetMomentumModulus()!=0.0) {
      // Used to calculate path, range traveled in this step
      double density; //  = 7.847; //SwimGeo::GetSwimMaterialDensity(material);
      TVector3 startpos, startmom;
      startpos = particle.GetPosition();
      startmom = particle.GetMomentum();
      
      //GMA Change momentum too (convert old fortran code to C++)
      // 1. Transform co-ordisnate system such that magnetic field is along Z' axis
      // 2. Get distance to the crossing point of heliz and plane
      // 3. Get the track paramters at the crossing point
      // 4. Return back to INO co-ordinate system
      
      double pos[3], newpos[3];
      double mom[3], newmom[3];
      double dir[3];
      double projpos[3], newprojpos[3];
      double track[6], trackout[6];
      
      pos[0] = startpos.X();
      pos[1] = startpos.Y();
      pos[2] = startpos.Z();
      
      //      int ilay = (int)fabs((distToNextPlane)/(2*fHalfLayerThickness));
      //      double localZ = fabs(distToNextPlane - 2*ilay*fHalfLayerThickness) -fHalfLayerThickness;
      double localZ = startpos.z() - startposXYZ.z();
      
      double B[6];
      double buvz[3];
      double pos1[3];
      double Bx,By;
      pos1[0]=pos[0]*1000;
      pos1[1]=pos[1]*1000;
      pos1[2]=pos[2]*1000;
      
      if (abs(localZ) < fHalfLayerThickness-fHalfAirGap){ // 125) { //AAR:these variables added to include variable airgap
	if(pFieldMap==NULL)cout<<"ERRRRR"<<endl;
	pFieldMap->ElectroMagneticField( pos1, Bx,By,1);
	B[0] = Bx*1000; B[1] =By*1000; B[2]=0;
	density = 7.847; //iron
      } else if (abs(localZ) < fHalfLayerThickness-0.004) {
	B[0] = 0.0; B[1] = 0.0; B[2]= 0.0;
	density = 0.00001; // //Foam (neglecting strips)
      } else if (abs(localZ) < fHalfLayerThickness-0.003) {
	B[0] = 0.0; B[1] = 0.0; B[2]= 0.0;
	density = 1.09;// G10;
      } else if (abs(localZ) < fHalfLayerThickness-0.001) {
	B[0] = 0.0; B[1] = 0.0; B[2]= 0.0;
	density = 2.32;// quartz;
      } else {
	B[0] = 0.0; B[1] = 0.0; B[2]= 0.0;
	density = 0.00137; //RPC gas
      }
      if( TMath::Sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2])!=0){
	fbave += TMath::Sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);                     //AAR
	fnbfield++;
      }
      
      buvz[0] = 0.3*B[0]; buvz[1] = 0.3*B[1]; buvz[2]= 0.3*B[2];                //AAR
      mom[0] = startmom.X();
      mom[1] = startmom.Y();
      mom[2] = startmom.Z();
      
      double mommag = startmom.Mag();
      for (int jk=0; jk<3; jk++) {dir[jk] = mom[jk]/mommag;}
      
      //Projected position;
      //      double dist = fabs(fStep/dir[2]); //startmom.CosTheta(); //GMA does it need absolute value ?
      
      double dist = fStep; // *direction[2];
      
      // If material is Air or Scint
      // straight line approximation
      
      //double buvz[3] ={0.3, 0.3, 0.0}; //GMA-magnetic field
      //	double buvz[3] ={0.45, 0.0, 0.0}; //GMA-magnetic field
      // GMA BX=By=1Tesla, Bz=0 and 0.3 factor comes from
      // the conversion factor, pt=0.3Br;
      
      double charge = particle.GetCharge();
      if (fIsForward) { charge *=-1.;}
      //GMA 06/02/2009 Chage -ve why ???????????????
      
      if (!fIsForward) { for (int jk=0; jk<3; jk++) {dir[jk] = -dir[jk];};}
      
      for (int jk=0; jk<3; jk++) {projpos[jk] = pos[jk] + dist*dir[jk];}
      
      //include bending in x-y plane
      double momxy = pow(max(1.e-12, mom[0]*mom[0] + mom[1]*mom[1]), 0.5);
      double momyz = pow(max(1.e-12, mom[1]*mom[1] + mom[2]*mom[2]), 0.5);
      double momzx = pow(max(1.e-12, mom[2]*mom[2] + mom[0]*mom[0]), 0.5);
      
      //GMA 06/02/2009 Looks like we swapped charge convention in the remaining coding
      projpos[0] -=0.5*charge*(fStep*fStep*dir[2]*dir[2]) * (buvz[2]/momxy - buvz[1]/momzx);
      projpos[1] -=0.5*charge*(fStep*fStep*dir[2]*dir[2]) * (buvz[0]/momyz - buvz[2]/momxy);
      
      // GMA 06/02/2009 Calculate average momentum in those two points (start and end)
      
      double rot[9]; //rot[3][3];
      // Calculate rotation matrix to rotate co-oridnate system such that
      // Z' is the direction of local Magnetic field
      anal_getnrot(buvz, rot);
      
      //rotate local postion in new frame
      
      anal_rotme(pos,rot,newpos);
      anal_rotme(mom,rot,newmom);
      anal_rotme(projpos,rot,newprojpos);
      
      for (int jk=0; jk<3; jk++) {
	track[jk] = newpos[jk];
	
	if (fIsForward) {
	  track[jk+3] = newmom[jk];
	} else {
	  track[jk+3] = -newmom[jk];
	}
      }
      
      double distep = 0; // total track length
      //if((buvz[0]*buvz[0]+buvz[1]*buvz[1]+buvz[2]*buvz[2])>0.0001)
      track_move_pt_align(buvz,track, charge, newprojpos, trackout, distep);//AAR: 1st argument added to pass the Field array.
      //else{ for(int ii=0;ii<6;ii++){ trackout[ii]=track[ii];} }
      for (int jk=0; jk<3; jk++) {
	newpos[jk] = trackout[jk];
	
	if (fIsForward) {
	  newmom[jk] = trackout[jk+3];
	} else {
	  newmom[jk] = -trackout[jk+3];
	}
      }
      
      //return back to original co-ordinate system
      anal_rotmet(newpos,rot,pos);
      anal_rotmet(newmom,rot,mom);
      
      //Forward backward part is taken care off by the trace_track_planef ??
      
      //===============================================================
      // Set position, momentum of the particle
      //===============================================================
      
      TVector3 position(pos[0], pos[1], pos[2]);
      particle.SetPosition(position);
      
      //GMA use momentum dependent energy loss factor
      double tmpmom = mommag;
      double elos = (14.9+0.96*fabs(log(tmpmom*2.8))+0.033*tmpmom*(1.0-pow(tmpmom,-0.33)))*1e-2*1.1;
      
      //       <<"elos =============="<<tmpmom<<" "<<elos<<endl;
      if (fIsForward) {
	tmpmom -=density*distep*elos; // 0.175;
      } else {
	tmpmom +=density*distep*elos; // 0.175;
      }
      
      for (int jk=0; jk<3; jk++) {mom[jk] = tmpmom*mom[jk]/mommag;}
      TVector3 momentum(mom[0], mom[1], mom[2]);
      particle.SetMomentum(momentum);
      
      // Calculate path, range in this step, add to particle totals
      
      particle.AddS(distep);
      particle.AddRange(density*distep); // /(Munits::g/Munits::cm2));
      
      //      <<"posmom4 "<<ij<<" "<<pos[0]<<" "<< pos[1]<<" "<< pos[2]<<" "<<mom[0]<<" "<< mom[1]<<" "<< mom[2]<<" "<<fStep<<" "<<invdircos<<" "<<dir[2]<<" "<<direction[2]<<endl;
      
      if (fabs(pos[2]-startposXYZ.z())<0.5*fStepSize*dir[2] && ij!=0) {
	//Crossing boundary again
	if (lastCrossingShift.Mag() !=0) {
	  cout<<"Crossing same layer twice, what is wrong ???????????? "<<endl;
	  cout<<"start "<<ij<<" "<< pos[2]<<" "<<startposXYZ.z()<<" "<<fStepSize<<" "<<  particle.GetMomentum().Mag()<<" "<< particle.GetMomentum().Theta()<<" "<<particle.GetMomentum().Phi()<<endl;
	}
	lastCrossingShift = position - startposXYZ;
      } else if (fabs(pos[2]-startposXYZ.z())>=fDistance || fStep<1.e-5) {
	if (pos[2]-startposXYZ.z()>0) {
	  nextplane = Plane + 1; //(pos[2]-startposXYZ.z()>0) ? 1: -1;
	} else {
	  nextplane = Plane - 1;
	}
	//	<<" next "<<Plane<<" "<<nextplane <<" "<<pos[2]<<" "<<startposXYZ.z()<<" "<<pos[2]-startposXYZ.z()<<" "<<fDistance<<endl;
	break;
      }
      
    }
  }
  
  //  return satisfied;
  fbave /=fnbfield;
  if (nextplane !=-99) {
    return fbave;
    //return true;
  } else {
    return 0;
  }
}

void SwimSwimmer::anal_getnrot( double* n, double* rot ) {
  //C             get rotation matrix to align Z-axis
  //C             with direction given by a  vector
  //  real*8 n(3),rot(*)
  
  double thd,phd,nn,nl[3], initrot[9];
  
  nn = sqrt( n[0]*n[0]+ n[1]*n[1]+ n[2]*n[2]);
  if( nn ==0.0 ) nn=1.00;
  nl[0]=n[0]/nn;
  nl[1]=n[1]/nn;
  nl[2]=n[2]/nn;
  
  thd=acos(max( -1.0, min( 1.0, nl[2] )));
  if( nl[0] != 0.0 || nl[1] !=0.0 ) {
    phd = atan2( nl[1],nl[0] );
  } else {
    phd = 0.0;
  }
  anal_getarot(thd,phd,initrot);
  for (int ij=0; ij<9; ij++) {rot[ij] = initrot[ij];}
  
}


void SwimSwimmer::anal_getarot(double thd, double phd, double* rot) {
  //  C             get rotation matrix to align Z-axis
  // C             with direction given by theta, phi
  
  double cthd,sthd,cphd,sphd;
  
  cthd=cos(thd);
  sthd=sin(thd);
  cphd=cos(phd);
  sphd=sin(phd);
  
  rot[0] = pow(sphd,2.)+cthd*pow(cphd,2.);
  rot[1] = (cthd-1)*sphd*cphd;
  rot[2] = - sthd*cphd;
  
  rot[3] = (cthd-1)*sphd*cphd;
  rot[4] = pow(cphd,2.)+cthd*pow(sphd,2.);
  rot[5] = - sthd*sphd;
  
  rot[6] = sthd*cphd;
  rot[7] = sthd*sphd;
  rot[8] = cthd;
  
  /*
    rot[0][0] = pow(sphd,2.)+cthd*pow(cphd,2.);
    rot[1][0] = (cthd-1)*sphd*cphd;
    rot[2][0] = - sthd*cphd;
    
    rot[0][1] = (cthd-1)*sphd*cphd;
    rot[2][1] = pow(cphd,2.)+cthd*pow(sphd,2.);
    rot[2][1] = - sthd*sphd;
    
    rot[0][2] = sthd*cphd;
    rot[1][2] = sthd*sphd;
    rot[2][2] = cthd;
  */
}


void SwimSwimmer::track_move_pt_align( double* b, double* tin, double q, double* x, double* tout, double& dl) {//AAR: 1st argument added to pass the Field array.
/*
C -------------------------------------------------------------------
C
C  move trajectory to the point of the closest approach in x-y to a point
C
C Input:
C       tin(6) - initial trajectory parameters  1-3 position 4-6 momentum
C       q   - charge
C       x(2)   - point to move to
C Output:
C       tout(6) - trajectory parameters at the closest approach
C       dl   - pathlength
C
C
C  $Id: rico_track_move_point.F,v 1.1.1.1 1999/09/16 03:10:50 ts Exp $
C
C  $Log: rico_track_move_point.F,v $
C  Revision 1.1.1.1  1999/09/16 03:10:50  ts
C  first release
C
C
C Tomasz Skwarnicki 7/8/98
C adopted from mcfast's move_wtk_point_bz by Paul Avery (see below)
c  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c      integer function move_wtk_point_bz(w1, x, bf, w2, s3d)
c
c  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c  Takes a helix expressed in w form and calculates the new parameters
c  nearest (nearest in the x-y plane) the point x(1-3) in a solenoidal
c  B field.
c
c  w1       w track structure (read)
c           Initial track parameters
c
c  x(3)     DFLOAT array (read)
c           Point to project track to
c
c  bf       B field structure (read)
c           B field information
c
c *w2       w track structure (write)
c           Track parameters at position closest to x(1-3)
c
c *s3d      DFLOAT variable (write)
c           3-D arc length the track moved between points
c
c return
c           0 ==> all OK
c           1 ==> track cannot be projected
c  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c   The equations of motion are
c
c     Px = Px0*cos(rho*s) - Py0*sin(rho*s)
c     Py = Py0*cos(rho*s) + Px0*sin(rho*s)
c     Pz = Pz0
c
c      x = x0 + {Px0*sin(rho*s) - Py0*[1-cos(rho*s)]} / a = x0 + (Py-Py0)/a
c      y = y0 + {Py0*sin(rho*s) + Px0*[1-cos(rho*s)]} / a = y0 - (Px-Px0)/a
c      z = z0 + ct*s
c
c   where s = arc length in r-phi plane
c         a = c_b * Bfield * q (c_b is defined in const.inc)
c       rho = a / Pt
c
c   We have to find the point closest to (xc,yc,zc). The solution is
c
c        cos(rho*s) = (1. - rho*(dx*Py0 - dy*Px0) / Pt) / norm
c        sin(rho*s) = (-rho*(dx*Px0 + dy*Py0) / Pt) / norm
c
c        norm = sqrt(1+2*rho*(dy*Px0 - dx*Py0)/Pt + rho**2*(dx**2 + dy**2))
c          dx = x0 - xc
c          dy = y0 - yc
c  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      double precision tin(6),q,x(2),tout(6),dl
C................................................................
#include "RichTracksProd/rich_magfield.inc"
c     local variables
      DOUBLE PRECISION rho, delx, dely, sinps, cosps, dcosps
      DOUBLE PRECISION alpha, a, ainv, sovp, ptinv, px, py
c  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      DOUBLE PRECISION pt
c     quit if no pt
*/
  
  double rho, delx, dely, sinps, cosps, dcosps;
  double alpha, a, ainv, sovp, ptinv, px, py;
  
  //GMA use these from database
  double c_b = 1.0;
  
  double brich = pow(b[0]*b[0]+b[1]*b[1]+b[2]*b[2],0.5); //GMA-magnetic field
  
  if(brich<0.00000001)brich=0.00001;
  //double brich = pow(0.0*0.0+0.0*0.0+0.0*0.0,0.5); //GMA-magnetic field
  
  double pt = pow(pow(tin[3],2.)+pow(tin[4],2.),0.5);
  if(pt != 0.0 && q != 0.0) {
    //    a = - c_b * brich * q; //GMA
    
    a = c_b * brich * q;
    ptinv = 1. / pt;
    rho = a * ptinv;
    delx = tin[0] - x[0];
    dely = tin[1] - x[1];
    ainv = 1. / a;
    px = tin[3];
    py = tin[4];
    cosps =  1. - rho*(delx*py - dely*px) * ptinv;
    sinps = -rho*(delx*px + dely*py) * ptinv;
    sovp = atan2(sinps, cosps) * ainv;
    
    alpha =  1. / sqrt(cosps*cosps + sinps*sinps);
    cosps = cosps * alpha;
    sinps = sinps * alpha;
    dcosps = 1. - cosps;
    tout[3] = px*cosps - py*sinps;
    tout[4] = py*cosps + px*sinps;
    tout[5] = tin[5];
    tout[0] = tin[0] + (px*sinps - py*dcosps) * ainv;
    tout[1] = tin[1] + (py*sinps + px*dcosps) * ainv;
    tout[2] = tin[2] + abs(sovp) * tin[5];
    dl = fabs(sovp) * pow(pow(tin[3],2.0)+pow(tin[4],2.)+pow(tin[5],2.),0.5);
  } else {
    for (int ij=0; ij<6; ij++) {tout[ij] = tin[ij];}
    dl = 0.;
  }
}


void SwimSwimmer::trace_track_planef(double * b, double* track, double charge, double* pl, double& dlstep) { //AAR: 1st argument added to pass the Field array.
/* Intersect track with a plane

    Inputs/Output:    track(6) - track parameters:
                            1-3 position  4-6 momentum
    InPuts       charge   - track charge
	         pl(6) - plane definitions 1-3 origion 4-6 normal
    Output       dlstep - pathlenght;
C.
C.......................................................................
C.
C.  rico_trace_track_planes - Intersect track with ordered planes
C.
C.
C.  Inputs:    track(6) - track parameters:
C.                            1-3 position  4-6 momentum
C.             charge   - track charge
C.
C.             npl - number of planes
C.            *pl(6,npl) - plane definitions 1-3 origion 4-6 normal
C.
C.  Outputs:   nc - number of planes crossed
C.             xc(6,nc) - intersection points and momenta
C
C.            *pl will be modified if normal vectors are not pointing
C.                towards the origin of coordinates
C.
C.
C.      track must be closer to the coordinate origin than a plane
C.      it is also not allowed to move away from the plane before
C.          intersecting it
C.      tracing stops if these conditions are not satisified
C.
C.
C.  Tomasz Skwarnicki 7/3/98
C.                    8/25/99 fix bug with no checking on number of plains crossed 
C.
C
C  $Id: rico_trace_track_planes.F,v 1.1.1.1 1999/09/16 03:10:50 ts Exp $
C
C  $Log: rico_trace_track_planes.F,v $
C  Revision 1.1.1.1  1999/09/16 03:10:50  ts
C  first release
C
C.......................................................................
      double precision track(6),charge
      integer npl
      double precision  pl(6,*)
C .....
      integer nc
      double precision xc(6,*)
C ............................................................

#include "RichTracksProd/rich_magfield.inc"
*      double precision brich
*      double precision c_b
*      parameter (brich = -14.9279) 
*      parameter (c_b = 0.299792458D0 * 0.1D0 / 100.0D0 )
*/      
  double brich;
  double c_b = 1.0;
  double trk[6],trkn[6];
  //  double dlstep,dlstep_old;
  //  double dist,dist_new,dist_old;
  //  double step,step_acc;
  //C step size when away from bounderies (5mm)
  //  step=0.005D0;
  //C requested accuracy in step when crossing bounderies 0.1mm
  //  step_acc=0.0001D0;
  
  bool nohelix;
  double  bf,chargel,ptot,pperp,rh,wh,xh,yh,phih,dt,x;
  rh = wh = xh = yh = 0;
  
  //C..............................................................
  //  nc = 0;
  //C...............................................................
  //C make sure all normal vectors point into the origin of coordinates
  //  C flip if necessary
  //  dist = pl[0]*pl[3] + pl[1]*pl[4] + pl[2]*pl[5];
  //  if( dist >0 ) {
  //    pl[3]= -pl[3];
  //    pl[4]= -pl[4];
  //    pl[5]= -pl[5];
  //  }
  
  // -------------------- ini from rico_track_move -----------------
  for (int ij=0; ij<6; ij++) {
    trk[ij]= track[ij];
  }
  
  //c ------------------------------------------------
  bf = brich =pow(b[0]*b[0]+b[1]*b[1]+b[2]*b[2],0.5);
  chargel = charge;
  
  // total momentum
  ptot = pow( pow(trk[3],2.)+pow(trk[4],2.)+pow(trk[5],2.), 0.5);
  // momentum perpendicular to the magnetic field direction
  pperp=pow(pow(trk[3],2.)+pow(trk[4],2.), 0.5);
  // see if need to do helix or straight line
  nohelix = (bf==0.0 || chargel ==0.0 || pperp <1.0e-4);
  if(!nohelix ) {
    // work with positive field, invert charge if field is negative
    if( bf <0.0 ) chargel = - chargel;
    // change units to tesla
    //      bf = abs(bf*0.1)
    // radius of helix in m
    //c      rh = pperp /(0.29979*bf*abs(chargel)) * 100.0
    rh = pperp /(c_b*abs(bf*chargel));
    // frequency
    wh = - pperp / rh * chargel;
    // center of the helix
    xh = trk[0] + rh * trk[4]/pperp * chargel;
    yh = trk[1] - rh * trk[3]/pperp * chargel;
  }
  // ------------------------------------------------------------
  dlstep = (trk[0]-pl[0]*pl[3]) + (trk[1]-pl[1]*pl[4]) + (trk[2]-pl[2]*pl[5]);
  
  if( nohelix ) {
    dt = dlstep; // / ptot;
    trkn[0]=trk[0]+trk[3]*dt;
    trkn[1]=trk[1]+trk[4]*dt;
    trkn[2]=trk[2]+trk[5]*dt;
    trkn[3]=trk[3];
    trkn[4]=trk[4];
    trkn[5]=trk[5];
  } else {
    // ---------------------- from rico_track_move ---------------
    dt = dlstep / ptot;
    // initial phase
    phih = atan2( trk[1]-yh,trk[0]-xh );
    // final phase
    x = phih + wh*dt;
    // calculate new position and momentum vector
    trkn[0] = xh + rh*cos(x);
    trkn[1] = yh + rh*sin(x);
    trkn[2] = trk[2] + trk[5]*dt;
    
    trkn[3] = - rh*wh * sin(x);
    trkn[4] =   rh*wh * cos(x);
    trkn[5] =   trk[5];
    // -----------------------------------------------------------
  }
  
  for (int ij=0; ij<6; ij++) {
    track[ij]= trkn[ij];
  }
}

void SwimSwimmer::anal_rotme(double* xx,double* rot,double* yy) {
  //C                 rotation of the vector
  double xt[3];
  for (int ij=0; ij<3; ij++) { xt[ij] = xx[ij];}
  
  for (int ij=0; ij<3; ij++) {
    yy[ij]=0;
    for (int jk=0; jk<3; jk++) {
      yy[ij] += rot[jk+3*ij]*xt[jk]; //rot[jk][ij]*xt[jk];
    }
  }
}

void SwimSwimmer::anal_rotmet(double* xx, double* rot,double* yy) {
  //C             inverse rotation of the vector
  double xt[3];
  
  for (int ij=0; ij<3; ij++) {xt[ij] = xx[ij]; }
  for (int ij=0; ij<3; ij++) {
    yy[ij]=0;
    for (int jk=0; jk<3; jk++) {
      yy[ij] +=rot[ij+3*jk]*xt[jk];// rot[ij][jk]*xt[jk];
    }
  }
}

