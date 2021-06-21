#include "micalPrimaryGeneratorAction.hh"

#include "micalDetectorConstruction.hh"
#include "micalPrimaryGeneratorMessenger.hh"
#include "micalDetectorParameterDef.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4Box.hh"

#include "math.h"
#include "CLHEP/Random/RandGauss.h"

//using namespace std;
//#include "micalDetectorParameterDef.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
micalPrimaryGeneratorAction *micalPrimaryGeneratorAction::AnPointer;
micalPrimaryGeneratorAction::micalPrimaryGeneratorAction(
							 micalDetectorConstruction* micalDC, MultiSimAnalysis *panalysis)
  :micalDetector(micalDC), pAnalysis(panalysis) { 
  AnPointer =this;
  // G4cout<<" Initialized micalPrimaryGeneratorAction Constructor"<<endl;
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  //Default settings :
  SetRunNumber(0);
  SetInputFlag(0);
  SetFirstEvt(1);
  SetRndmFlag("on");
  SetPartId(13);
  SetIncEnergy(2.0*GeV);
  SetIncEnergySmr(0*MeV);
  SetIncDirection(G4ThreeVector(0.0,0.0,-1.0));
  SetIncThetaSmr(0*mrad);
  SetIncPhiSmr(0*mrad);
  SetIncPosition(G4ThreeVector(0.0*cm,0.0*cm,0*cm));
  SetIncVxSmr(0*cm);
  SetIncVySmr(0*cm);
  SetIncVzSmr(0*cm);
  
  enerin[0] = 0.3;
  enerin[1] = 0.4;
  enerin[2] = 0.5;
  enerin[3] = 0.6;
  enerin[4] = 0.7;
  enerin[5] = 0.8;
  enerin[6] = 0.9;
  enerin[7] = 1.0;
  enerin[8] = 1.1;
  enerin[9] = 1.2;
  enerin[10] = 1.3;
  enerin[11] = 1.4;
  enerin[12] = 1.5;
  enerin[13] = 1.8;
  enerin[14] = 2.0;
  enerin[15] = 2.5;
  enerin[16] = 3.0;
  enerin[17] = 4.0;
  enerin[18] = 5.0;
  enerin[19] = 6.0;
  
  pivalGA = acos(-1.0);
  initialise = 0;
  initialiseCor = 0;
  g_nevt=-1;
  
  //create a messenger for this class
  gunMessenger = new micalPrimaryGeneratorMessenger(this);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalPrimaryGeneratorAction::~micalPrimaryGeneratorAction() {
  // G4cout<<"micalPrimaryGeneratorAction Distructor"<<endl;
  if (particleGun)	{delete particleGun;}
  if (gunMessenger) {delete gunMessenger;}
}

void micalPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  // cout<<"micalPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {..."<<endl;
  // int nevt, npart, ipart;
  // G4int      pid;
  G4double   vx=0, vy=0, vz=0;
  // G4double   px=0, py=0, pz=0;
  // G4double   atime;
  
  if (initialise==0) {
    // gunMessenger = new micalPrimaryGeneratorMessenger(this);
    paradef = micalDetectorParameterDef::AnPointer;
    for(int ij=0; ij<3; ij++) {
      parino[ij] = paradef->GetParino(ij);
      pargas[ij] = paradef->GetPargas(ij);
      StackPosInWorld[ij] = paradef->GetStackPosInRoom(ij) + paradef->GetINOroomPos(ij);
      cout<<"GetStackPosInRoom "<<ij<<" "<<paradef->GetStackPosInRoom(ij)<<" "<<paradef->GetINOroomPos(ij)<<endl;
    }
    cout<<"StackPosInWorld[3] "<<StackPosInWorld[0]<<" "<<StackPosInWorld[1]<<" "<<StackPosInWorld[2]<<endl;
    
    WorldXDim = paradef->GetParworld(0);
    WorldYDim = paradef->GetParworld(1);
    WorldZDim = paradef->GetParworld(2);
    nINODet = paradef->GetNumino();
    nLayer=paradef->GetnLayer();
    nIRLayer=paradef->GetnIRLayer();
    for(int ij=0; ij<nLayer; ij++) {
      RPCLayerPosZ[ij] = paradef->GetRPCLayerPosZ(ij);
      LayerZdim[ij] = paradef->GetLayerZdim(ij);
    }
    for(int ij=0; ij<nIRLayer; ij++) {
      IRONLayerPosZ[ij] = paradef->GetIRONLayerPosZ(ij);
      IronLayerZdim[ij] = paradef->GetIronLayerZdim(ij);
    }
    g_nevt=-1;
    initialise = 1;
    FirstEvt=pAnalysis->FirstEvt;
  }

  //this function is called at the begining of event
  // default particle kinematic
  
  pAnalysis->irun = RunNumber;  //Keep an option that in a file, one may have more than two run number
  g_nevt++;
  // cout<<"I am here..."<<InputFlag<<endl;
  //inputFlag = 0 for G4MC & =0 for Nuance 2 :GINIE
  if (InputFlag==0) { // G4MC case
    pAnalysis->ievt=g_nevt;
    pAnalysis->ngent = ((unsigned)particleGun->GetNumberOfParticles() <=pAnalysis->ngenmx) ?
      particleGun->GetNumberOfParticles() : pAnalysis->ngenmx;
    
    for (int ij=0; ij<particleGun->GetNumberOfParticles(); ij++)	{
      //Option to have multiple particle
      G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
      // partId = 13;
      G4ParticleDefinition* particle = particleTable->FindParticle(partId);
      particleGun->SetParticleDefinition(particle);
      
      G4ThreeVector ini_Dir(incDirection);
      double in_Energy = incEnergy*MeV;
      // vx = incPosition.x()*mm;
      // vy = incPosition.y()*mm;
      // vz = incPosition.z()*mm;
      
      vx = 0.0;
      vy = 0.0;
      vz = 0.0;
      // vz = IRONLayerPosZ[nIRLayer-1] + 2*mm + IronLayerZdim[nIRLayer-1]/2;
      // vz = IRONLayerPosZ[0] - 2*mm - IronLayerZdim[0]/2;

      // cout<<"PGA before = "<<vx<<" "<<vy<<" "<<vz<<endl;

      
      vx += StackPosInWorld[0];
      vy += StackPosInWorld[1];
      vz += StackPosInWorld[2];
      
      // cout<<"PGA after = "<<vx<<" "<<vy<<" "<<vz<<endl;

      if (rndmFlag=="on") {
	G4double theta=0;
	if(incThetaSmr>=0    && incDirection.theta() !=0) {
	  theta = G4RandGauss::shoot(0,incThetaSmr);
	} else if(incThetaSmr<0 || incDirection.theta()==0) {
	  theta = incThetaSmr*(2*G4UniformRand()-1);
	}
	
	G4double phi=0;
	if(incPhiSmr>=0    && incDirection.theta() !=0) {
	  phi = G4RandGauss::shoot(0,incPhiSmr);
	} else if(incPhiSmr<0 || incDirection.theta()==0) {
	  phi = incPhiSmr*(2*G4UniformRand()-1);
	}
	
	// ini_Dir.setTheta(ini_Dir.theta()+theta);
	// ini_Dir.setPhi(ini_Dir.phi()+phi);

	theta = ini_Dir.theta();
	bool tmpDD = true;//false;
	while(tmpDD) {
	  double xx11 = G4UniformRand();// - 1;
	  if(xx11<1.00 && xx11>0.5) {
	    theta = acos(-xx11);
	    break;
	  }
	}
	theta = acos(-1.0);//-0.95);
	phi = pivalGA*(2*G4UniformRand()-1);
	
	// theta = acos(-1);
	// phi = 0;

	ini_Dir.setTheta(theta);
	ini_Dir.setPhi(phi);
	
	// if(incEnergySmr>=0) {
	//   in_Energy += G4RandGauss::shoot(0,incEnergySmr);
	//   if (in_Energy <1*MeV) in_Energy=1*MeV;
	// } else if(incEnergySmr<0) {
	//   in_Energy += incEnergySmr*(2*G4UniformRand()-1);
	//   if (in_Energy <1*MeV) in_Energy=1*MeV;
	// }
	
	in_Energy = 1*GeV;//(3.5 + 3.0*(2*G4UniformRand()-1))*GeV;
	// in_Energy = 0.7*GeV;
	// in_Energy = incEnergy*MeV;
	// in_Energy = enerin[g_nevt/5000]*GeV;
	// cout<<vx<<","<< vy<<","<< vz<<endl;
	// 	if(incVxSmr>=0&&incVySmr>=0&&incVzSmr>=0) {
	// 	  vx += 0*mm;
	// 	  vy += 0*mm;
	// 	  vz += 0*mm;
	// 	} else {
	// 	  vx += incVxSmr*(2*G4UniformRand()-1.)*mm;
	// 	  vy += incVySmr*(2*G4UniformRand()-1.)*mm;
	// 	  vz += incVzSmr*(2*G4UniformRand()-1.)*mm;
	// 	}
      }
      
      particleGun->SetParticleMomentumDirection(ini_Dir);
      particleGun->SetParticleMomentum(in_Energy);
      particleGun->SetParticlePosition(G4ThreeVector(vx, vy, vz));
      particleGun->GeneratePrimaryVertex(anEvent);
      if (ij < (int)pAnalysis->ngenmx) {
	pAnalysis->pidin[ij] = particleGun->GetParticleDefinition()->GetPDGEncoding();
	pAnalysis->posxin[ij]= particleGun->GetParticlePosition().x()/m;
	pAnalysis->posyin[ij]= particleGun->GetParticlePosition().y()/m;
	pAnalysis->poszin[ij]= particleGun->GetParticlePosition().z()/m;
	double momentum = particleGun->GetParticleMomentum()/GeV; 
	if(particle->GetPDGCharge()==0) {
	  pAnalysis->momin[ij] = momentum;
	} else {
	  pAnalysis->momin[ij] = momentum*(particle->GetPDGCharge());
	}
	pAnalysis->thein[ij] = particleGun->GetParticleMomentumDirection().theta();
	pAnalysis->phiin[ij] = particleGun->GetParticleMomentumDirection().phi();
	
	// cout<<"PGA = "<<ij<<" "<<particleGun->GetParticlePosition()/m<<" "<<pAnalysis->momin[ij]<<" "<<particleGun->GetParticleMomentum()/GeV<<" "<<pAnalysis->thein[ij]*180/pivalGA<<" "<<pAnalysis->phiin[ij]*180/pivalGA<<endl;
      } // if (ij < (int)pAnalysis->ngenmx)
    } // for (int ij=0; ij<particleGun->GetNumberOfParticles(); ij++)
  } else if (InputFlag==1) {
    // Cosmic Flux from data
    pAnalysis->ievt=g_nevt;
    pAnalysis->ngent = ((unsigned)particleGun->GetNumberOfParticles() <=pAnalysis->ngenmx) ?
      particleGun->GetNumberOfParticles() : pAnalysis->ngenmx;
    for (int ij=0; ij<particleGun->GetNumberOfParticles(); ij++) {
      G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
      G4ParticleDefinition* particle = particleTable->FindParticle(partId);
      particleGun->SetParticleDefinition(particle);
      
      G4ThreeVector ini_Dir(incDirection);
      double rand_theta;
      double theta;
      double theta_gen;
      double phi_gen;
      double phi;
      int brkpt = 1;
      double Point2[3];
      // double energy;
      double vertexX;
      double vertexY;
      double vertexZ;
      double Ini_Theta = 0;
      double Ini_Phi = 0;

      while(brkpt) {
	vx = 2.4*pargas[0]*(G4UniformRand()-0.5);
	vy = 2.4*pargas[1]*(G4UniformRand()-0.5);
	vz = RPCLayerPosZ[toptrgly];
	//   for Standalone Theta generation
	rand_theta = G4UniformRand();
	theta = acos(pow((1-rand_theta*norm1),(1./(PowCosTheta+1.0))))*rad;
	phi = pivalGA*(2*G4UniformRand()-1)*rad;

	double Line1[6];
	double Plane1[6];
	double Point1[3] = {-100000.,-100000.,-100000.};
	Line1[0] = vx;
	Line1[1] = vy;
	Line1[2] = vz;
	Line1[3] = -sin(theta)*cos(phi);
	Line1[4] = -sin(theta)*sin(phi);
	Line1[5] = -cos(theta);
	Plane1[0] =  0;
	Plane1[1] =  0;
	Plane1[2] =  RPCLayerPosZ[bottomtrgly];;
	Plane1[3] = 0;
	Plane1[4] = 0;
	Plane1[5] = 1;
	
	int trgCheck = LinePlaneInt(Line1,Plane1,Point1);
	if(trgCheck == 1) {
	  if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	    double Line2[6];
	    double Plane2[6];
	    for(int xxi=0;xxi<3;xxi++) {Point2[xxi] = -100000000.;}	    
	    for(int lmn=0;lmn<3;lmn++) {Point2[lmn] = Point1[lmn];}
	    Line2[0] = StackPosInWorld[0] + vx;
	    Line2[1] = StackPosInWorld[1] + vy;
	    Line2[2] = StackPosInWorld[2] + vz;
	    Line2[3] = -sin(theta)*cos(phi);
	    Line2[4] = -sin(theta)*sin(phi);
	    Line2[5] = -cos(theta);

	    Plane2[0] =  0;
	    Plane2[1] =  0;
	    Plane2[2] =  WorldZDim - 1*mm;
	    Plane2[3] = 0;
	    Plane2[4] = 0;
	    Plane2[5] = 1;

	    int TopPlane = LinePlaneInt(Line2,Plane2,Point2);
	    if(TopPlane ==1) {
	      if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
		vertexX = Point2[0];
		vertexY = Point2[1];
		vertexZ = Point2[2];
		theta_gen = theta;
		phi_gen = phi;
		Ini_Theta = theta + pivalGA*rad;
		Ini_Phi = phi - pivalGA*rad;
		if(Ini_Phi < -pivalGA) {
		  Ini_Phi = Ini_Phi + 2*pivalGA;
		} else if(Ini_Phi > pivalGA) {
		  Ini_Phi = Ini_Phi - 2*pivalGA;
		}
		brkpt = 0;
	      } // if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
	    } // if(TopPlane ==1) {
	  } // if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	} // if(trgCheck == 1) {
      } // while(brkpt) {

      double in_Energy = GetCosmicEnergy(ELowLim,EUpLim);
      
      ini_Dir.setTheta(Ini_Theta);
      ini_Dir.setPhi(Ini_Phi);
      particleGun->SetParticleMomentumDirection(ini_Dir);
      particleGun->SetParticleEnergy(in_Energy);
      particleGun->SetParticlePosition(G4ThreeVector(vertexX, vertexY, vertexZ));
      particleGun->GeneratePrimaryVertex(anEvent);
      if (ij < (int)pAnalysis->ngenmx) {
	pAnalysis->pidin[ij] = particleGun->GetParticleDefinition()->GetPDGEncoding();
	pAnalysis->posxin[ij]= particleGun->GetParticlePosition().x();
	pAnalysis->posyin[ij]= particleGun->GetParticlePosition().y();
	pAnalysis->poszin[ij]= particleGun->GetParticlePosition().z();
	if(particle->GetPDGCharge()==0) {
	  pAnalysis->momin[ij] = particleGun->GetParticleMomentum()/GeV;
	} else {
	  pAnalysis->momin[ij] = (particleGun->GetParticleMomentum()/GeV)*(particle->GetPDGCharge());
	}
	pAnalysis->thein[ij] = particleGun->GetParticleMomentumDirection().theta();
	pAnalysis->phiin[ij] = particleGun->GetParticleMomentumDirection().phi();
      } // if (ij < (int)pAnalysis->ngenmx) {    
    } // for (int ij=0; ij<particleGun->GetNumberOfParticles(); ij++) {
  } else if (InputFlag==2) {
    // Cosmic Flux 3D hist Corsika
    if(initialiseCor==0) {
      FileCORSIKA->cd();
      corsikaFlux = (TH3D*)FileCORSIKA->Get("flux");
      initialiseCor = 1;
    }
    pAnalysis->ievt=g_nevt;
    pAnalysis->ngent = ((unsigned)particleGun->GetNumberOfParticles() <=pAnalysis->ngenmx) ? particleGun->GetNumberOfParticles() : pAnalysis->ngenmx;
    for (int ij=0; ij<particleGun->GetNumberOfParticles(); ij++) {
      G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
      G4ParticleDefinition* particle = particleTable->FindParticle(partId);
      particleGun->SetParticleDefinition(particle);
      G4ThreeVector ini_Dir(incDirection);
      double vertexX;
      double vertexY;
      double vertexZ;
      double Ini_Theta = 0;
      double Ini_Phi = 0;
      double Ini_Enrgy = 0;
      double costheta;
      double phi;
      double logEnrgy;
      double theta;
      double enrgy;
      int brkpt = 1;
      double Point2[3];
      while(brkpt) {
	vx = 2.4*pargas[0]*(G4UniformRand()-0.5);
	vy = 2.4*pargas[1]*(G4UniformRand()-0.5);
	vz = RPCLayerPosZ[toptrgly];
	
	corsikaFlux->GetRandom3(logEnrgy,costheta,phi);
	phi = phi*pivalGA/180;
	phi = phi*rad;
	theta = acos(costheta);
	enrgy = pow(10,logEnrgy);
	double Line1[6];
	double Plane1[6];
	double Point1[3] = {-100000.,-100000.,-100000.};
	Line1[0] = vx;
	Line1[1] = vy;
	Line1[2] = vz;
	Line1[3] = -sin(theta)*cos(phi);
	Line1[4] = -sin(theta)*sin(phi);
	Line1[5] = -cos(theta);
	Plane1[0] =  0;
	Plane1[1] =  0;
	Plane1[2] =  RPCLayerPosZ[bottomtrgly];;
	Plane1[3] = 0;
	Plane1[4] = 0;
	Plane1[5] = 1;
	
	int trgCheck = 1; //LinePlaneInt(Line1,Plane1,Point1);
	if(trgCheck == 1) {
	  // if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	  {
	    double Line2[6];
	    double Plane2[6];
	    for(int xxi=0;xxi<3;xxi++) {Point2[xxi] = -100000000.;}	    
	    for(int lmn=0;lmn<3;lmn++) {Point2[lmn] = Point1[lmn];}
	    Line2[0] = StackPosInWorld[0] + vx;
	    Line2[1] = StackPosInWorld[1] + vy;
	    Line2[2] = StackPosInWorld[2] + vz;
	    Line2[3] = -sin(theta)*cos(phi);
	    Line2[4] = -sin(theta)*sin(phi);
	    Line2[5] = -cos(theta);

	    Plane2[0] =  0;
	    Plane2[1] =  0;
	    Plane2[2] =  WorldZDim - 1*mm;
	    Plane2[3] = 0;
	    Plane2[4] = 0;
	    Plane2[5] = 1;
	    int TopPlane = 1; //LinePlaneInt(Line2,Plane2,Point2);
	    if(TopPlane ==1) {
	      // if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
	      {
		vertexX = 0;//Point2[0];
		vertexY = 0;//Point2[1];
		vertexZ = 0;//Point2[2];
		Ini_Theta = theta;
		Ini_Phi = phi;
		if(Ini_Phi < -pivalGA) {
		  Ini_Phi = Ini_Phi + 2*pivalGA;
		} else if(Ini_Phi > pivalGA) {
		  Ini_Phi = Ini_Phi - 2*pivalGA;
		}
		Ini_Enrgy = enrgy;
		brkpt = 0;
	      } // if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
	    } // if(TopPlane ==1) {
	  } // if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	} // // if(trgCheck == 1) {
      } // while(brkpt) {	
      
      ini_Dir.setTheta(Ini_Theta);
      ini_Dir.setPhi(Ini_Phi);
      particleGun->SetParticleMomentumDirection(ini_Dir);
      particleGun->SetParticleEnergy(Ini_Enrgy);
      particleGun->SetParticlePosition(G4ThreeVector(vertexX, vertexY, vertexZ));
      particleGun->GeneratePrimaryVertex(anEvent);
      if (ij < (int)pAnalysis->ngenmx) {
	pAnalysis->pidin[ij] = particleGun->GetParticleDefinition()->GetPDGEncoding();
	pAnalysis->posxin[ij]= particleGun->GetParticlePosition().x();
	pAnalysis->posyin[ij]= particleGun->GetParticlePosition().y();
	pAnalysis->poszin[ij]= particleGun->GetParticlePosition().z();
	if(particle->GetPDGCharge()==0) {
	  pAnalysis->momin[ij] = particleGun->GetParticleMomentum()/GeV;
	} else {
	  pAnalysis->momin[ij] = (particleGun->GetParticleMomentum()/GeV)*(particle->GetPDGCharge());
	}
	pAnalysis->thein[ij] = particleGun->GetParticleMomentumDirection().theta();
	pAnalysis->phiin[ij] = particleGun->GetParticleMomentumDirection().phi();
      } // if (ij < (int)pAnalysis->ngenmx) {
    } // for (int ij=0; ij<particleGun->GetNumberOfParticles(); ij++) {
  } else if (InputFlag==3) {
    if (initialiseCor==0) {
      TreeCORSIKA = (TTree*)FileCORSIKA->Get("corsikaTreeAll");
      TreeCORSIKA->SetBranchAddress("iC_nevt",&iC_nevt);
      TreeCORSIKA->SetBranchAddress("iC_eventweight",&iC_eventweight);
      TreeCORSIKA->SetBranchAddress("iC_npart",&iC_npart);
      TreeCORSIKA->SetBranchAddress("iC_cpid",iC_cpid);
      TreeCORSIKA->SetBranchAddress("iC_cvx",iC_cvx);
      TreeCORSIKA->SetBranchAddress("iC_cvy",iC_cvy);
      TreeCORSIKA->SetBranchAddress("iC_cpx",iC_cpx);
      TreeCORSIKA->SetBranchAddress("iC_cpy",iC_cpy);
      TreeCORSIKA->SetBranchAddress("iC_cpz",iC_cpz);
      initialiseCor=1;
    } // if (initialiseCor==0) {
    if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
      TreeCORSIKA->GetEntry(FirstEvt+g_nevt);
      pAnalysis->ievt		= iC_nevt;
      pAnalysis->ievt_wt	= iC_eventweight;
      pAnalysis->ngent = (iC_npart < (int)pAnalysis->ngenmx) ? iC_npart : pAnalysis->ngenmx;
      for (unsigned int count=0; count<pAnalysis->ngent; count++) {
	if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {
	  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	  G4ParticleDefinition* particle = particleTable->FindParticle(iC_cpid[count]);
	  particleGun->SetParticleDefinition(particle);
	  particleGun->SetParticlePosition(G4ThreeVector(0,0,0));

	  /*** Check this properly ---> MeV or GeV ***/
	  G4ThreeVector tmp3v(iC_cpx[count]*MeV,iC_cpy[count]*MeV,iC_cpz[count]*MeV);
	  // G4ThreeVector tmp3v(iC_cpx[count]*GeV,iC_cpy[count]*GeV,iC_cpz[count]*GeV);
	  
	  double Point2[3];
	  // double energy;
	  double vertexX;
	  double vertexY;
	  double vertexZ;
	  double Ini_Theta = 0;
	  double Ini_Phi = 0;
	  G4ThreeVector ini_Dir(tmp3v.unit());
	  double theta, phi;
	  int brkpt = 1;
	  int brkcnt = 0;

	  pAnalysis->ngenerated = 0;
	  pAnalysis->naperture = 0;

	  while(brkpt) {
	    vx = 1.2*pargas[0]*(2*G4UniformRand()-1.0);
	    vy = 1.2*pargas[1]*(2*G4UniformRand()-1.0);
	    vz = RPCLayerPosZ[toptrgly];
	    
	    if(abs(vx)<pargas[0] && abs(vy)<pargas[1]) {
	      pAnalysis->ngenerated++;
	    }
	    phi = ini_Dir.phi();
	    theta = ini_Dir.theta();
	    double Line1[6];
	    double Plane1[6];
	    double Point1[3] = {-100000.,-100000.,-100000.};
	    Line1[0] = vx;
	    Line1[1] = vy;
	    Line1[2] = vz;
	    Line1[3] = -sin(theta)*cos(phi);
	    Line1[4] = -sin(theta)*sin(phi);
	    Line1[5] = -cos(theta);
	    Plane1[0] =  0;
	    Plane1[1] =  0;
	    Plane1[2] =  RPCLayerPosZ[bottomtrgly];;
	    Plane1[3] = 0;
	    Plane1[4] = 0;
	    Plane1[5] = 1;
	    int trgCheck = LinePlaneInt(Line1,Plane1,Point1);
	    if(trgCheck == 1) {
	      if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
		pAnalysis->naperture++;
		double Line2[6];
		double Plane2[6];
		for(int lmn=0;lmn<3;lmn++) {Point2[lmn] = Point1[lmn];}
		Line2[0] = StackPosInWorld[0] + vx;
		Line2[1] = StackPosInWorld[1] + vy;
		Line2[2] = StackPosInWorld[2] + vz;
		Line2[3] = -sin(theta)*cos(phi);
		Line2[4] = -sin(theta)*sin(phi);
		Line2[5] = -cos(theta);
		
		Plane2[0] =  0;
		Plane2[1] =  0;
		Plane2[2] =  WorldZDim - 1*mm;
		Plane2[3] = 0;
		Plane2[4] = 0;
		Plane2[5] = 1;
		int TopPlane = LinePlaneInt(Line2,Plane2,Point2);
		if(TopPlane ==1) {
		  if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
		    vertexX = Point2[0];
		    vertexY = Point2[1];
		    vertexZ = Point2[2];
		    Ini_Theta = theta;
		    Ini_Phi = phi;
		    if(Ini_Phi < -pivalGA) {
		      Ini_Phi = Ini_Phi + 2*pivalGA;
		    } else if(Ini_Phi > pivalGA) {
		      Ini_Phi = Ini_Phi - 2*pivalGA;
		    }
		    brkpt = 0;

		  } // if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
		} // if(TopPlane ==1) {
	      } // if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	    }// if(trgCheck == 1) {
	    brkcnt++;
	    if(brkcnt>1000) {
	      cout<<"brkcnt "<<brkcnt<<endl;
	      vertexX = 0;
	      vertexY = 0; //Point2[1];
	      vertexZ = WorldZDim - 1*mm; //Point2[2];
	      Ini_Theta = theta;
	      Ini_Phi = phi;
	      if(Ini_Phi < -pivalGA) {
		Ini_Phi = Ini_Phi + 2*pivalGA;
	      } else if(Ini_Phi > pivalGA) {
		Ini_Phi = Ini_Phi - 2*pivalGA;
	      }
	      brkpt = 0;
	      cout<<"brkt "<<brkcnt<<" "<<Ini_Theta*180/pivalGA <<" "<<Ini_Phi*180/pivalGA<<endl;
	    }
	  } // while(brkpt) {	
	  // cout << " vertexX " << vertexX << " vertexY " << vertexY << " vertexZ " << vertexZ << endl;
	  
	  // Ini_Enrgy = enrgy;
	  
	  ini_Dir.setTheta(Ini_Theta);
	  ini_Dir.setPhi(Ini_Phi);
	  particleGun->SetParticleMomentumDirection(ini_Dir);
	  particleGun->SetParticleMomentum(tmp3v.mag());
	  particleGun->SetParticlePosition(G4ThreeVector(vertexX, vertexY, vertexZ));
	  particleGun->GeneratePrimaryVertex(anEvent);
	  if (count < (int)pAnalysis->ngenmx) {
	    pAnalysis->pidin [count]= particleGun->GetParticleDefinition()->GetPDGEncoding();//_cpid[count];
	    pAnalysis->posxin[count]= particleGun->GetParticlePosition().x();
	    pAnalysis->posyin[count]= particleGun->GetParticlePosition().y();
	    pAnalysis->poszin[count]= particleGun->GetParticlePosition().z();
	    if(particle->GetPDGCharge()==0){
	      pAnalysis->momin[count] = particleGun->GetParticleMomentum()/GeV;
	    } else {
	      pAnalysis->momin[count] = (particleGun->GetParticleMomentum())*(particle->GetPDGCharge())/GeV;
	    }
	    pAnalysis->thein[count] = particleGun->GetParticleMomentumDirection().theta();
	    pAnalysis->phiin[count] = particleGun->GetParticleMomentumDirection().phi();
	  } // if (count < (int)pAnalysis->ngenmx) {
	} // if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {
      } // for (int count=0; count<iC_npart; count++) {
    } // if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
    // Corsika Event By Event
  } else if (InputFlag==4) {
    if (initialiseCor==0) {
      FileFLUX->cd();
      muFlux = (TH3F*)FileFLUX->Get("muFlux");
      mupFlux = (TH3F*)FileFLUX->Get("mupFlux");
      munFlux = (TH3F*)FileFLUX->Get("munFlux");
      initialiseCor=1;
    } // if (initialiseCor==0) {
    if(1) {
      pAnalysis->ievt		= g_nevt;
      pAnalysis->ievt_wt	= 1.;
      pAnalysis->ngent = particleGun->GetNumberOfParticles();	// number of particles
      // cout << " npart " <<  pAnalysis->ngent << endl;
      for (unsigned int count=0; count<pAnalysis->ngent; count++) {
	if(1) {
	  
	  double Pxx,Pyy,Pzz;
	  muFlux->GetRandom3(Pxx,Pyy,Pzz);
	  int gBin = muFlux->FindBin(Pxx,Pyy,Pzz);
	  double mupCnt = mupFlux->GetBinContent(gBin);
	  double munCnt = munFlux->GetBinContent(gBin);

	  int partID;
	  if(G4UniformRand()*(mupCnt+munCnt)>munCnt) {
	    partID = -13;
	  } else {
	    partID = 13;
	  }
	  // cout << " partID " << partID << endl;
	  
	  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	  G4ParticleDefinition* particle = particleTable->FindParticle(partID);
	  particleGun->SetParticleDefinition(particle);
	  particleGun->SetParticlePosition(G4ThreeVector(0,0,0));

	  /*** Check this properly ---> MeV or GeV ***/
	  // G4ThreeVector tmp3v(iC_cpx[count]*MeV,iC_cpy[count]*MeV,iC_cpz[count]*MeV);
	  G4ThreeVector tmp3v(Pxx*GeV,Pyy*GeV,Pzz*GeV);
	  
	  double Point2[3];
	  // double energy;
	  double vertexX;
	  double vertexY;
	  double vertexZ;
	  double Ini_Theta = 0;
	  double Ini_Phi = 0;
	  G4ThreeVector ini_Dir(tmp3v.unit());
	  double theta, phi;
	  int brkpt = 1;
	  int brkcnt = 0;

	  pAnalysis->ngenerated = 0;
	  pAnalysis->naperture = 0;

	  while(brkpt) {
	    vx = 1.2*pargas[0]*(2*G4UniformRand()-1.0);
	    vy = 1.2*pargas[1]*(2*G4UniformRand()-1.0);
	    vz = RPCLayerPosZ[toptrgly];
	    
	    if(abs(vx)<pargas[0] && abs(vy)<pargas[1]) {
	      pAnalysis->ngenerated++;
	    }
	    phi = ini_Dir.phi();
	    theta = ini_Dir.theta();
	    double Line1[6];
	    double Plane1[6];
	    double Point1[3] = {-100000.,-100000.,-100000.};
	    Line1[0] = vx;
	    Line1[1] = vy;
	    Line1[2] = vz;
	    Line1[3] = -sin(theta)*cos(phi);
	    Line1[4] = -sin(theta)*sin(phi);
	    Line1[5] = -cos(theta);
	    Plane1[0] =  0;
	    Plane1[1] =  0;
	    Plane1[2] =  RPCLayerPosZ[bottomtrgly];;
	    Plane1[3] = 0;
	    Plane1[4] = 0;
	    Plane1[5] = 1;
	    int trgCheck = LinePlaneInt(Line1,Plane1,Point1);
	    if(trgCheck == 1) {
	      if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
		pAnalysis->naperture++;
		double Line2[6];
		double Plane2[6];
		for(int lmn=0;lmn<3;lmn++) {Point2[lmn] = Point1[lmn];}
		Line2[0] = StackPosInWorld[0] + vx;
		Line2[1] = StackPosInWorld[1] + vy;
		Line2[2] = StackPosInWorld[2] + vz;
		Line2[3] = -sin(theta)*cos(phi);
		Line2[4] = -sin(theta)*sin(phi);
		Line2[5] = -cos(theta);
		
		Plane2[0] =  0;
		Plane2[1] =  0;
		Plane2[2] =  WorldZDim - 1*mm;
		Plane2[3] = 0;
		Plane2[4] = 0;
		Plane2[5] = 1;
		int TopPlane = LinePlaneInt(Line2,Plane2,Point2);
		if(TopPlane ==1) {
		  if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
		    vertexX = Point2[0];
		    vertexY = Point2[1];
		    vertexZ = Point2[2];
		    Ini_Theta = theta;
		    Ini_Phi = phi;
		    if(Ini_Phi < -pivalGA) {
		      Ini_Phi = Ini_Phi + 2*pivalGA;
		    } else if(Ini_Phi > pivalGA) {
		      Ini_Phi = Ini_Phi - 2*pivalGA;
		    }
		    brkpt = 0;

		  } // if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
		} // if(TopPlane ==1) {
	      } // if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	    }// if(trgCheck == 1) {
	    brkcnt++;
	    if(brkcnt>1000) {
	      cout<<"brkcnt "<<brkcnt<<endl;
	      vertexX = 0;
	      vertexY = 0; //Point2[1];
	      vertexZ = WorldZDim - 1*mm; //Point2[2];
	      Ini_Theta = theta;
	      Ini_Phi = phi;
	      if(Ini_Phi < -pivalGA) {
		Ini_Phi = Ini_Phi + 2*pivalGA;
	      } else if(Ini_Phi > pivalGA) {
		Ini_Phi = Ini_Phi - 2*pivalGA;
	      }
	      brkpt = 0;
	      cout<<"brkt "<<brkcnt<<" "<<Ini_Theta*180/pivalGA <<" "<<Ini_Phi*180/pivalGA<<endl;
	    }
	  } // while(brkpt) {	
	  // cout << " count " << count << " vertexX " << vertexX << " vertexY " << vertexY << " vertexZ " << vertexZ << endl;
	  // cout << " mom " << tmp3v.mag() << " " << Ini_Theta << " " << Ini_Phi << endl;
	  
	  // Ini_Enrgy = enrgy;
	  // double Ini_Enrgy = (G4UniformRand()*EUpLim+ELowLim)*MeV;
	  // cout << " Ini_Enrgy " << Ini_Enrgy << endl;
	  // cout << " ELowLim " << ELowLim << " EUpLim " << EUpLim << endl;	  
	  
	  ini_Dir.setTheta(Ini_Theta);
	  ini_Dir.setPhi(Ini_Phi);
	  particleGun->SetParticleMomentumDirection(ini_Dir);
	  particleGun->SetParticleMomentum(tmp3v.mag());
	  particleGun->SetParticlePosition(G4ThreeVector(vertexX, vertexY, vertexZ));
	  particleGun->GeneratePrimaryVertex(anEvent);
	  if (count < (int)pAnalysis->ngenmx) {
	    pAnalysis->pidin [count]= particleGun->GetParticleDefinition()->GetPDGEncoding();
	    pAnalysis->posxin[count]= particleGun->GetParticlePosition().x();
	    pAnalysis->posyin[count]= particleGun->GetParticlePosition().y();
	    pAnalysis->poszin[count]= particleGun->GetParticlePosition().z();
	    if(particle->GetPDGCharge()==0){
	      pAnalysis->momin[count] = particleGun->GetParticleMomentum()/GeV;
	    } else {
	      pAnalysis->momin[count] = (particleGun->GetParticleMomentum())*(particle->GetPDGCharge())/GeV;
	    }
	    pAnalysis->thein[count] = particleGun->GetParticleMomentumDirection().theta();
	    pAnalysis->phiin[count] = particleGun->GetParticleMomentumDirection().phi();
	  } // if (count < (int)pAnalysis->ngenmx) {
	} // if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {
      } // for (int count=0; count<iC_npart; count++) {
    } // if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
  } else if (InputFlag==5) {
    if (initialiseCor==0) {
      cout << " InputFlag==5 " << endl;
      FileFLUX->cd();
      muFlux = (TH3F*)FileFLUX->Get("muFlux");
      mupFlux = (TH3F*)FileFLUX->Get("mupFlux");
      munFlux = (TH3F*)FileFLUX->Get("munFlux");
      initialiseCor=1;
      cout << " InputFlag==5 " << endl;
    } // if (initialiseCor==0) {
    if(1) {
      pAnalysis->ievt		= g_nevt;
      pAnalysis->ievt_wt	= 1.;
      pAnalysis->ngent = particleGun->GetNumberOfParticles();	// number of particles
      // cout << " npart " <<  pAnalysis->ngent << endl;
      for (unsigned int count=0; count<pAnalysis->ngent; count++) {
	if(1) {
	  
	  double Pxx,Pyy,Pzz;
	  muFlux->GetRandom3(Pxx,Pyy,Pzz);
	  int gBin = muFlux->FindBin(Pxx,Pyy,Pzz);
	  double mupCnt = mupFlux->GetBinContent(gBin);
	  double munCnt = munFlux->GetBinContent(gBin);

	  int partID;
	  if(G4UniformRand()*(mupCnt+munCnt)>munCnt) {
	    partID = -13;
	  } else {
	    partID = 13;
	  }
	  // cout << " partID " << partID << endl;
	  
	  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	  G4ParticleDefinition* particle = particleTable->FindParticle(partID);
	  particleGun->SetParticleDefinition(particle);
	  particleGun->SetParticlePosition(G4ThreeVector(0,0,0));

	  /*** Check this properly ---> MeV or GeV ***/
	  // G4ThreeVector tmp3v(iC_cpx[count]*MeV,iC_cpy[count]*MeV,iC_cpz[count]*MeV);
	  G4ThreeVector tmp3v(Pxx*GeV,Pyy*GeV,Pzz*GeV);
	  
	  double Point2[3];
	  // double energy;
	  double vertexX;
	  double vertexY;
	  double vertexZ;
	  double Ini_Theta = 0;
	  double Ini_Phi = 0;
	  G4ThreeVector ini_Dir(tmp3v.unit());
	  double theta, phi;
	  int brkpt = 1;
	  int brkcnt = 0;

	  pAnalysis->ngenerated = 0;
	  pAnalysis->naperture = 0;

	  while(brkpt) {
	    vx = 1.2*pargas[0]*(2*G4UniformRand()-1.0);
	    vy = 1.2*pargas[1]*(2*G4UniformRand()-1.0);
	    vz = RPCLayerPosZ[toptrgly];
	    
	    if(abs(vx)<pargas[0] && abs(vy)<pargas[1]) {
	      pAnalysis->ngenerated++;
	    }
	    phi = ini_Dir.phi();
	    theta = ini_Dir.theta();
	    double Line1[6];
	    double Plane1[6];
	    double Point1[3] = {-100000.,-100000.,-100000.};
	    Line1[0] = vx;
	    Line1[1] = vy;
	    Line1[2] = vz;
	    Line1[3] = -sin(theta)*cos(phi);
	    Line1[4] = -sin(theta)*sin(phi);
	    Line1[5] = -cos(theta);
	    Plane1[0] =  0;
	    Plane1[1] =  0;
	    Plane1[2] =  RPCLayerPosZ[bottomtrgly];;
	    Plane1[3] = 0;
	    Plane1[4] = 0;
	    Plane1[5] = 1;
	    int trgCheck = LinePlaneInt(Line1,Plane1,Point1);
	    if(trgCheck == 1) {
	      if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
		pAnalysis->naperture++;
		double Line2[6];
		double Plane2[6];
		for(int lmn=0;lmn<3;lmn++) {Point2[lmn] = Point1[lmn];}
		Line2[0] = StackPosInWorld[0] + vx;
		Line2[1] = StackPosInWorld[1] + vy;
		Line2[2] = StackPosInWorld[2] + vz;
		Line2[3] = -sin(theta)*cos(phi);
		Line2[4] = -sin(theta)*sin(phi);
		Line2[5] = -cos(theta);
		
		Plane2[0] =  0;
		Plane2[1] =  0;
		Plane2[2] =  WorldZDim - 1*mm;
		Plane2[3] = 0;
		Plane2[4] = 0;
		Plane2[5] = 1;
		int TopPlane = LinePlaneInt(Line2,Plane2,Point2);
		if(TopPlane ==1) {
		  if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
		    vertexX = Point2[0];
		    vertexY = Point2[1];
		    vertexZ = Point2[2];
		    Ini_Theta = theta;
		    Ini_Phi = phi;
		    if(Ini_Phi < -pivalGA) {
		      Ini_Phi = Ini_Phi + 2*pivalGA;
		    } else if(Ini_Phi > pivalGA) {
		      Ini_Phi = Ini_Phi - 2*pivalGA;
		    }
		    brkpt = 0;

		  } // if(abs(Point2[0])<WorldXDim && abs(Point2[1])<WorldYDim) {
		} // if(TopPlane ==1) {
	      } // if(abs(Point1[0])<pargas[0] && abs(Point1[1])<pargas[1]) {
	    }// if(trgCheck == 1) {
	    brkcnt++;
	    if(brkcnt>1000) {
	      cout<<"brkcnt "<<brkcnt<<endl;
	      vertexX = 0;
	      vertexY = 0; //Point2[1];
	      vertexZ = WorldZDim - 1*mm; //Point2[2];
	      Ini_Theta = theta;
	      Ini_Phi = phi;
	      if(Ini_Phi < -pivalGA) {
		Ini_Phi = Ini_Phi + 2*pivalGA;
	      } else if(Ini_Phi > pivalGA) {
		Ini_Phi = Ini_Phi - 2*pivalGA;
	      }
	      brkpt = 0;
	      cout<<"brkt "<<brkcnt<<" "<<Ini_Theta*180/pivalGA <<" "<<Ini_Phi*180/pivalGA<<endl;
	    }
	  } // while(brkpt) {	
	  // cout << " count " << count << " vertexX " << vertexX << " vertexY " << vertexY << " vertexZ " << vertexZ << endl;
	  // cout << " mom " << tmp3v.mag() << " " << Ini_Theta << " " << Ini_Phi << endl;
	  
	  // Ini_Enrgy = enrgy;
	  double Ini_Enrgy = (G4UniformRand()*EUpLim+ELowLim)*MeV;
	  // cout << " Ini_Enrgy " << Ini_Enrgy << endl;
	  // cout << " ELowLim " << ELowLim << " EUpLim " << EUpLim << endl;	  
	  
	  // Ini_Phi = pivalGA*(2*G4UniformRand()-1)*rad;
	  
	  ini_Dir.setTheta(Ini_Theta);
	  ini_Dir.setPhi(Ini_Phi);
	  particleGun->SetParticleMomentumDirection(ini_Dir);
	  // particleGun->SetParticleMomentum(tmp3v.mag());
	  particleGun->SetParticleMomentum(Ini_Enrgy);
	  particleGun->SetParticlePosition(G4ThreeVector(vertexX, vertexY, vertexZ));
	  particleGun->GeneratePrimaryVertex(anEvent);
	  if (count < (int)pAnalysis->ngenmx) {
	    pAnalysis->pidin [count]= particleGun->GetParticleDefinition()->GetPDGEncoding();
	    pAnalysis->posxin[count]= particleGun->GetParticlePosition().x();
	    pAnalysis->posyin[count]= particleGun->GetParticlePosition().y();
	    pAnalysis->poszin[count]= particleGun->GetParticlePosition().z();
	    if(particle->GetPDGCharge()==0){
	      pAnalysis->momin[count] = particleGun->GetParticleMomentum()/GeV;
	    } else {
	      pAnalysis->momin[count] = (particleGun->GetParticleMomentum())*(particle->GetPDGCharge())/GeV;
	    }
	    pAnalysis->thein[count] = particleGun->GetParticleMomentumDirection().theta();
	    pAnalysis->phiin[count] = particleGun->GetParticleMomentumDirection().phi();
	  } // if (count < (int)pAnalysis->ngenmx) {
	} // if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {
      } // for (int count=0; count<iC_npart; count++) {
    } // if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
  } else if (InputFlag==6) {
    if (initialiseCor==0) {
      TreeCORSIKA = (TTree*)FileCORSIKA->Get("corsikaTreeAll");
      TreeCORSIKA->SetBranchAddress("iC_nevt",&iC_nevt);
      TreeCORSIKA->SetBranchAddress("iC_eventweight",&iC_eventweight);
      TreeCORSIKA->SetBranchAddress("iC_npart",&iC_npart);
      TreeCORSIKA->SetBranchAddress("iC_cpid",iC_cpid);
      TreeCORSIKA->SetBranchAddress("iC_cvx",iC_cvx);
      TreeCORSIKA->SetBranchAddress("iC_cvy",iC_cvy);
      TreeCORSIKA->SetBranchAddress("iC_cpx",iC_cpx);
      TreeCORSIKA->SetBranchAddress("iC_cpy",iC_cpy);
      TreeCORSIKA->SetBranchAddress("iC_cpz",iC_cpz);
      initialiseCor=1;
    } // if (initialiseCor==0) {
    if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
      TreeCORSIKA->GetEntry(FirstEvt+g_nevt);
      pAnalysis->ievt		= iC_nevt;
      pAnalysis->ievt_wt	= iC_eventweight;
      pAnalysis->ngent = (iC_npart < (int)pAnalysis->ngenmx) ? iC_npart : pAnalysis->ngenmx;
      for (unsigned int count=0; count<pAnalysis->ngent; count++) {
	if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {
	  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	  G4ParticleDefinition* particle = particleTable->FindParticle(iC_cpid[count]);
	  particleGun->SetParticleDefinition(particle);
	  particleGun->SetParticlePosition(G4ThreeVector(0,0,0));
	  G4ThreeVector tmp3v(iC_cpx[count]*GeV,iC_cpy[count]*GeV,iC_cpz[count]*GeV);
	  // cout<<"Momin "<<tmp3v.mag()<<endl;
	  double Point2[3];
	  // double energy;
	  double vertexX;
	  double vertexY;
	  double vertexZ;
	  double Ini_Theta = 0;
	  double Ini_Phi = 0;
	  G4ThreeVector ini_Dir(tmp3v.unit());
	  double theta, phi;
	  int brkpt = 1;
	  int brkcnt = 0;

	  pAnalysis->ngenerated = 0;
	  pAnalysis->naperture = 0;

	  vx = 0;//500*(2*G4UniformRand()-1.0);
	  vy = 0;//500*(2*G4UniformRand()-1.0);
	  vz = RPCLayerPosZ[toptrgly] + 101*mm;
	  
	  vertexX = vx + StackPosInWorld[0];
	  vertexY = vy + StackPosInWorld[1];
	  vertexZ = vz + StackPosInWorld[2];
	  
	  phi = ini_Dir.phi();
	  theta = ini_Dir.theta();

	  Ini_Theta = theta;
	  Ini_Phi = phi;
	  if(Ini_Phi < -pivalGA) {
	    Ini_Phi = Ini_Phi + 2*pivalGA;
	  } else if(Ini_Phi > pivalGA) {
	    Ini_Phi = Ini_Phi - 2*pivalGA;
	  }
	  	  
	  // Ini_Enrgy = enrgy;
	    
	  ini_Dir.setTheta(Ini_Theta);
	  ini_Dir.setPhi(Ini_Phi);
	  particleGun->SetParticleMomentumDirection(ini_Dir);
	  particleGun->SetParticleMomentum(tmp3v.mag());
	  particleGun->SetParticlePosition(G4ThreeVector(vertexX, vertexY, vertexZ));
	  particleGun->GeneratePrimaryVertex(anEvent);
	  // cout<<"particleGun->PartPos "<<particleGun->GetParticlePosition()<<endl;
	  // cout<<"particleGun->Mom "<<particleGun->GetParticleMomentum()/GeV<<endl;
	  if (count < (int)pAnalysis->ngenmx) {
	    pAnalysis->pidin [count]= particleGun->GetParticleDefinition()->GetPDGEncoding();//_cpid[count];
	    pAnalysis->posxin[count]= particleGun->GetParticlePosition().x();
	    pAnalysis->posyin[count]= particleGun->GetParticlePosition().y();
	    pAnalysis->poszin[count]= particleGun->GetParticlePosition().z();
	    if(particle->GetPDGCharge()==0){
	      pAnalysis->momin[count] = particleGun->GetParticleMomentum()/GeV;
	    } else {
	      pAnalysis->momin[count] = (particleGun->GetParticleMomentum())*(particle->GetPDGCharge())/GeV;
	    }
	    pAnalysis->thein[count] = particleGun->GetParticleMomentumDirection().theta();
	    pAnalysis->phiin[count] = particleGun->GetParticleMomentumDirection().phi();
	  } // if (count < (int)pAnalysis->ngenmx) {
	} // if(iC_cpid[count]!=0 && abs(iC_cpid[count])<1000000) {
      } // for (int count=0; count<iC_npart; count++) {
    } // if (FirstEvt+g_nevt<TreeCORSIKA->GetEntries()) {
    // Corsika Event By Event
  } else {
    //Future NuGenerators
  }
  
}

double micalPrimaryGeneratorAction::GetCosmicEnergy(double ELimLow, double ELimUp) {

  double selectedEnergy = 0.001;
  bool energyBool = true;
  double Pmax =0.00364709;//17.7098;// 0.00364709(Alkofer);//168.23;// 0.019372(alkofer);//168.23;
  
  double E1 = ELimLow*0.001;
  double E2 = ELimUp*0.001;

  while(energyBool) {
    double x1 = E1 + (G4UniformRand()*(E2-E1));
    double x2 = G4UniformRand();
    double f1 = energy_func(x1);
    double f2 = x2*Pmax;
    // cout<<"E "<<E1<<" "<<E2<<endl;
    // cout<<"f2 "<<f2<<" "<<x2<<" "<<f1<<" "<<x1<<endl;
    if(f2 < f1) {
      selectedEnergy = x1;
      energyBool = false;
    }
  }
  
  // return 1.5*1000;
  return selectedEnergy*1000.;
  
}

int micalPrimaryGeneratorAction::LinePlaneInt(double* Line, double* Plane, double* Point) { //, double &Dist) {
  double a, b, Dist;
  int ok = 0;
  b = Line[3]*Plane[3]  + Line[4]*Plane[4] + Line[5]*Plane[5];
  ok= (fabs(b) > 1e-10) ? 1 : 0;
  
  if(ok==1){
    a=(Plane[0]-Line[0])*Plane[3] + (Plane[1]-Line[1])*Plane[4] + (Plane[2]-Line[2])*Plane[5];
    Dist = a/b;
    Point[0] = Line[0] + Line[3]*Dist;
    Point[1] = Line[1] + Line[4]*Dist;
    Point[2] = Line[2] + Line[5]*Dist;
  }else{Point[0]=0; Point[1]=0; Point[2]=0;}
  
  return ok;
}

double micalPrimaryGeneratorAction::energy_func(double xx) {
  //double paren[7]={1.44021e-02,4.09983e-01,2.96593e-01,7.46643e-03,-2.20450e+00,1.15441e+00,-4.19923e+00};
  double paren[7]={9.20933e-03,4.57797e-01,1.00738e+00 ,4.82396e-03,-1.54754e+00,1.60996e-01 ,-2.90163e+00}; //Pmax = 0.00364709;
  if (xx<1.5) {
    return paren[0]*(TMath::Gaus(xx, paren[1], paren[2], kTRUE));
  } else if (xx<14.0) {
    return paren[3]*pow(xx, paren[4]);
  } else {
    return paren[5]*pow(xx, paren[6]);
  }
} 

void micalPrimaryGeneratorAction::OpenFileCORSIKA() {
  G4String infile;
  infile = CorsikaFileDir;
  infile.append(CorsikaFileName);
  cout<<"CorsikaFileName "<<CorsikaFileName<<endl;
  FileCORSIKA = new TFile(infile,"READ","input file");
}

void micalPrimaryGeneratorAction::CloseFileCORSIKA() {
  
  FileCORSIKA->Close();
  delete FileCORSIKA;
}

void micalPrimaryGeneratorAction::OpenFileFLUX() {
  G4String infile;
  infile = CorsikaFileDir;
  infile.append(FluxFileName);
  cout<<"FluxFileName "<<FluxFileName<<endl;
  FileFLUX = new TFile(infile,"READ","flux file");
}

void micalPrimaryGeneratorAction::CloseFileFLUX() {
  FileFLUX->Close();
  delete FileFLUX;
}

void micalPrimaryGeneratorAction::SetInputFlag(G4int p) {
  InputFlag = p;
  // cout<<"Setting Input Flag..."<<InputFlag<<endl;
  
}

void micalPrimaryGeneratorAction::SetIncPosition(G4ThreeVector p) {
  incPosition = p;
  // cout<<"Setting new inc pos "<<incPosition<<endl;
}
