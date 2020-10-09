#include <cassert>
#include <cmath>
#include <iomanip>

#include "TMinuit.h"
#include "TRandom.h"

#include "InoTrackFinder.h"

#include "TObjArray.h"
#include "TObjArray.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/DiagMatrix.h"
#include "CLHEP/Matrix/Vector.h"

InoTrackFinder::InoTrackFinder() :
  NumModules(1), ModuleType(1), PECut(0.0001), PECut2(0.001) //GMA Need proper vlaue, particularly PECut(2)
			      //  NumModules(1), PlanesInModule(150), ModuleType(1), PECut(0.05), PECut2(0.1) //GMA Need proper vlaue, particularly PECut(2)
{
  inoHit_pointer = new InoHit_Manager();
  inoRPC_pointer = InoRPCStrip_Manager::APointer; //new InoRPCStrip_Manager();

  pAnalysis = MultiSimAnalysis::AnPointer;        //asm 
  inoCluster_pointer = new InoCluster_Manager();
  
  inoTrack_pointer = new InoTrack_Manager();
  
  paradef = micalDetectorParameterDef::AnPointer;
  NumModules= 1;
  
  StripXWidth = (1/m)*paradef->GetXStrwd();
  StripYWidth = (1/m)*paradef->GetYStrwd();
  DigiToTimeConv = pAnalysis->GetTimeToDigiConvVal();
  SignalSpeed = pAnalysis->GetSignalSpeedVal();
  // cout<<"StripXWidth "<<StripXWidth<<", StripYWidth "<<StripYWidth<<endl;
  nXStrips = paradef->GetnXStrip();
  nYStrips = paradef->GetnYStrip();
  GasChmX = 2*paradef->GetPargas(0)/1000;
  GasChmY = 2*paradef->GetPargas(1)/1000;
  LayerThickness = (1/m)*2*(paradef->GetParlay(2)+paradef->GetParirlay(2)); //(1/m)*2*paradef->GetParlay(2);
  PlanesInModule = paradef->GetnLayer(); //GMAAAA

  
}

InoTrackFinder::~InoTrackFinder() {
  //GMA need to clear after evergy events;
  //  inoHit_pointer->InoHit_list.clear();
  
  for (unsigned int ij=0; ij<inoCluster_pointer->InoCluster_list.size(); ij++) {
    if (inoCluster_pointer->InoCluster_list[ij]) {
      delete inoCluster_pointer->InoCluster_list[ij];
      inoCluster_pointer->InoCluster_list[ij]=0;
    }
  }
  inoCluster_pointer->InoCluster_list.clear(); delete inoCluster_pointer;
  
  for (unsigned int ij=0; ij<inoTrack_pointer->InoTrack_list.size(); ij++)	{
    if (inoTrack_pointer->InoTrack_list[ij]) {
      delete inoTrack_pointer->InoTrack_list[ij];
      inoTrack_pointer->InoTrack_list[ij]=0;
    }
  }
  inoTrack_pointer->InoTrack_list.clear(); delete inoTrack_pointer;
  ClearUp();
  // cout<<"Hello:"<<endl;
}

void InoTrackFinder::RunTheFinder() {
  
  // Configure algorithm for the relevant detector
  // For ND, initial algorithm is applied only to calorimeter,
  // spectrometer is treated separately
  // cout<< "FINDING TRACKS IN SLICE" << endl;
  
  // Run the methods
  FormTheHits();// slice);
  // cout<< "FINDING TRACKS IN SLICE1" << endl;
  FormTheClusters();
  // cout<< "FINDING TRACKS IN SLICE2" << endl;
  IDTrkAndShwClusters();
  // cout<< "FINDING TRACKS IN SLICE3" << endl;
  FormTriplets();
  // cout<< "FINDING TRACKS IN SLICE4" << endl;
  FindAllAssociations();
  // cout<< "FINDING TRACKS IN SLICE5" << endl;
  FindPreferredJoins();
  // cout<< "FINDING TRACKS IN SLICE6" << endl;
  FindMatchedJoins();
  // cout<< "FINDING TRACKS IN SLICE7" << endl;
  FirstComparison();
  //  if(ModuleType==2) {NearDetectorTriplets();}
  FormTracks();
  // cout<< "FINDING TRACKS IN SLICE8" << endl;
  JoinTracks();
  // cout<< "FINDING TRACKS IN SLICE9" << endl;
  FormFinalTracks(); //(slice);
  // cout<< "FINDING TRACKS IN SLICE10" << endl;
  JoinCurvedTrack();  
  // cout<< "FINDING TRACKS IN SLICE11" << endl;
  CleanAndFilled();
  // cout<< "FINDING TRACKS IN SLICE12" << endl;
  return;
}
//===================================================================================================================================
void InoTrackFinder::FormTheHits() {
  // cout<<"void InoTrackFinder::FormTheHits() {..."<<endl;
  //MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;        //asm
  // Create and (store pointers to) InoHit objects for all strips in slice
  // meeting PH and XTalk requirements. InoHit class is essentially a wrapper
  // for CandStripHandles
  micalDetectorParameterDef* paradef = micalDetectorParameterDef::AnPointer;
  InoStripX_Manager *pstripx = InoStripX_Manager::APointer;
  InoStripY_Manager *pstripy = InoStripY_Manager::APointer;
  
  //  vector <int> iXfill;
  //  vector <int> iYFill;
  //GMA need to sort out the definition of vector<int> it out
  int iXFill[5000]; //GMA14 need to return back to vecttor, same is true for some other cases aslo (in other classes)
  int iYFill[5000];
  //-----------------------------------------------------------------------------------------------------------

  int RPCXStrip[2*3*150*8*8]={0};
  int RPCYStrip[2*3*150*8*8]={0};

  if (pstripx) {
    for (unsigned ix=0; ix<pstripx->InoStripX_list.size(); ix++) {
      InoStrip* XStrip = pstripx->InoStripX_list[ix];

      int xStrpId = XStrip->GetId();
      xStrpId>>=8;
      int nInX = xStrpId & 0x7F;
      //cout<<"InoTrackFinder nInX = "<<nInX<<endl;
      xStrpId>>=7;
      int nInCH = xStrpId & 0x03;
      xStrpId>>=3;
      int nInMO = xStrpId & 0x03;
      xStrpId>>=3;
      int nInLA = xStrpId & 0x4F;
      xStrpId>>=8;
      int nInDT = xStrpId & 0x03;
      
      // if(1) {
      // 	//cout<<"Noisy X-Strip, Noisy Strip nInDT, nInLA, nInMO, nInCH, nInX = "<<nInDT<<", "<<nInLA<<", "<<nInMO<<", "<<nInCH<<", "<<nInX<<endl;
      // 	continue;
      // }
      RPCXStrip[XStrip->GetRPCmod()]++;
      pAnalysis->NoisyStripX->Fill(nInX);
      // cout<<"pAnalysis->NoisyStripX->Fill(nInX);="<<nInX<<endl;
      if (XStrip->GetPulse()<PECut) continue;
      //GMA Use random number to use efficiency 
      //Also use Poission distribution of avearage number of 
      // primary electron-ion pair
      
      //noisy strip
      
      if (pstripy) {
	for (unsigned jy=0; jy<pstripy->InoStripY_list.size(); jy++) {
	  InoStrip* YStrip = pstripy->InoStripY_list[jy];
	  
	  int yStrpId = YStrip->GetId();
	  yStrpId>>=8;
	  int nInY = yStrpId & 0x7F;
	  //cout<<"InoTrackFinder nInY = "<<nInY<<endl;
	  yStrpId>>=7;
	  int nninch = yStrpId & 0x03;
	  yStrpId>>=3;
	  int nninmo = yStrpId & 0x03;
	  yStrpId>>=3;
	  int nninla = yStrpId & 0x4F;
	  yStrpId>>=8;
	  int nnindt = yStrpId & 0x03;
	  
	  // if(1) {
	  //   //  cout<<"Noisy Y-Strip, Noisy Strip nnindt, nninla, nninmo, nninch, nInY = "<<nnindt<<", "<<nninla<<", "<<nninmo<<", "<<nninch<<", "<<nInY<<endl;
	  //   continue;
	  // }

	  if (ix==0) { 
	    RPCYStrip[YStrip->GetRPCmod()]++; 
	  }

	  pAnalysis->NoisyStripY->Fill(nInY);
	  // cout<<"pAnalysis->NoisyStripY->Fill(nInY);="<<nInY<<endl;
      
	  if (YStrip->GetPulse()<PECut) continue;
	  if (XStrip->GetRPCmod() != YStrip->GetRPCmod()) continue;
	  //GMA Check it
	  double Xtime = DigiToTimeConv*XStrip->GetSmrTime() - (nInY+0.5)*SignalSpeed; //(sigXspeed);
	  double Ytime = DigiToTimeConv*YStrip->GetSmrTime() - (nInX+0.5)*SignalSpeed; //(sigYspeed);

	  // double Xtime = XStrip->GetTimeX() - (nInY*GasChmY*5.0*ns/nYStrips);
	  // double Ytime = YStrip->GetTimeY() - (nInX*GasChmX*5.0*ns/nXStrips);
	  double diffTime = Xtime - Ytime;

	  pAnalysis->DiffTime->Fill(diffTime);
	  //	  cout <<"ixxxx "<< ix<<" "<<jy<<" "<<Xtime<<" "<< Ytime<<endl;
	  if (abs(diffTime)>5.0) continue; //Already converted this to ns (anyhow ns=1.0)
	  //	  if (fabs(XStrip->GetSmrTime() - YStrip->GetSmrTime()) > 5.0 * ns) continue;
	  
	  // Mark for editting
	  
	  InoHit* tmphit = new InoHit(XStrip, YStrip);
	  inoHit_pointer->InoHit_list.push_back(tmphit);


	 

	  HitBank[XStrip->GetPlane()].push_back(tmphit);
	  
	  iXFill[ix] = 1;
	  iYFill[jy] = 1;
	}
      }
    }
  }
  
  typedef pair<int,int> xystrp;
  for (int ij=0; ij<2*3*150*8*8; ij++) {
    if (RPCXStrip[ij]>0 || RPCYStrip[ij]>0) {
      xystrp abc(ij, 100*RPCXStrip[ij]+ RPCYStrip[ij]);
      inoRPC_pointer->InoRPCStrip.push_back(abc); // pair(ij, 100*RPCXStrip[ij]+ RPCYStrip[ij]));
    }
  }

  //-------------------------------------------------------------------------------------------------------------
  int pi_event=1;
  
  if(pi_event) {
    //    int vtxPlane=-1;
    
    pAnalysis->inohits_old	=0;
    pAnalysis->orighits_old	=0;
    pAnalysis->x_hits_old	=0;
    pAnalysis->y_hits_old	=0;
    
    int ixstripno = -1;//HitBank[i][j]->GetXStripNum();
    int iystripno = -1;//HitBank[i][j]->GetYStripNum();
    
    int Inohits=0;
    int xxnhits=0;
    int yynhits=0;
    int Orighits=0;
    
    vector<int> xshwstrip;
    vector<int> yshwstrip;

    for(int ijx=0; ijx<250; ijx++) {
      if(HitBank[ijx].size()) {
	for(unsigned int jky=0; jky<HitBank[ijx].size(); jky++) {
	  double hxTrueTime = HitBank[ijx][jky]->GetXTrueTime();
	  double hxSmrTime = HitBank[ijx][jky]->GetXTime();
	  double hxSmrTimeCorr = HitBank[ijx][jky]->GetXTimeCorr();
	  double hyTrueTime = HitBank[ijx][jky]->GetYTrueTime();
	  double hySmrTime = HitBank[ijx][jky]->GetYTime();
	  double hySmrTimeCorr = HitBank[ijx][jky]->GetYTimeCorr();
	  double hitTime = HitBank[ijx][jky]->GetTime();
	  
	  double DiffTimeX1 = hxSmrTime - hxTrueTime;
	  double DiffTimeY1 = hySmrTime - hyTrueTime;
	  double DiffTimeX2 = hitTime - hxTrueTime;
	  double DiffTimeY2 = hitTime - hyTrueTime;
	  double DiffTimeX3 = hxSmrTimeCorr - hxTrueTime;
	  double DiffTimeY3 = hySmrTimeCorr - hyTrueTime;
	  pAnalysis->strpXtime->Fill(DiffTimeX1);
	  pAnalysis->strpYtime->Fill(DiffTimeY1);

	  pAnalysis->hitXtime->Fill(DiffTimeX2);
	  pAnalysis->hitYtime->Fill(DiffTimeY2);
	  
	  pAnalysis->strpXtimeCorr->Fill(DiffTimeX3);
	  pAnalysis->strpYtimeCorr->Fill(DiffTimeY3);
 
	} // for(int jky=0; jky<HitBank[ij].size(); jky++) {
      }
    } // for(int ijx=0; ijx<250; ijx++) {

    for(int ij=0; ij<250; ++ij) {
      Inohits = HitBank[ij].size();
      xxnhits=0; yynhits=0;  Orighits=0;
      if(HitBank[ij].size()) {
	for(unsigned int jk=0; jk<HitBank[ij].size(); ++jk) {
	  ixstripno = HitBank[ij][jk]->GetXStripNum();
	  iystripno = HitBank[ij][jk]->GetYStripNum();
	  
	  for (unsigned nn=0; nn< xshwstrip.size(); nn++) {
	    if (ixstripno == xshwstrip[nn]) {ixstripno = -1; break;}
	  }
	  for (unsigned nn=0; nn< yshwstrip.size(); nn++) {
	    if (iystripno == yshwstrip[nn]) {iystripno = -1; break;}
	  }
	  if ((ixstripno>=0 || iystripno>=0)) {
	    if(ixstripno>=0) {
	      xshwstrip.push_back(ixstripno);
	    }
	    if(iystripno>=0) {
	      yshwstrip.push_back(iystripno);
	    }
	  }
	}
	
	xxnhits+=xshwstrip.size();
	yynhits+=yshwstrip.size();
	
	Orighits=max(xxnhits,yynhits);
	
	pAnalysis->x_hits_old	+= xxnhits;
	pAnalysis->y_hits_old	+= yynhits;
	pAnalysis->orighits_old	+=Orighits;
	pAnalysis->inohits_old	+= Inohits;
      }
      xshwstrip.clear();
      yshwstrip.clear();
    }
  }
  
  //--------------------------------------------------------------------------------------------------------------
  
  if (pAnalysis->ihist < pAnalysis->nhistmx-1 && pAnalysis->isVisOut>=2) {
    for(int ij=0; ij<500; ++ij) { //GMA14 change this static array to vector all 500 to ..
      for(unsigned int jk=0; jk<HitBank[ij].size(); ++jk) {
	pAnalysis->gens_list[1][pAnalysis->ihist]->Fill(HitBank[ij][jk]->GetXPos(),HitBank[ij][jk]->GetYPos(),HitBank[ij][jk]->GetZPos()+0.01);
	cout<<"pAnalysis->gens_list[1][pAnalysis->ihst]->Fill();"<<endl;
	vectGr  tmpgr;
	tmpgr.x = HitBank[ij][jk]->GetXPos();
	tmpgr.y = HitBank[ij][jk]->GetYPos();
	tmpgr.z = HitBank[ij][jk]->GetZPos()+0.01;               //asm ?+0.01
	
	tmpgr.dx = 0;
	tmpgr.dy = 0;
	tmpgr.dz = 0;
	
	if (pAnalysis->isVisOut==3) pAnalysis->gens_vect[1].push_back(tmpgr);
      }
    }
  }
  
  //----------------------------------------------------------------------------------------------cout

  if(pAnalysis->isXtermOut==1) {
  // if(1) {
    cout <<"TrkFinder HitBank "<<endl;
    
    for(int Module=0; Module<NumModules; ++Module) {
      for(int Plane=0; Plane<PlanesInModule+1; ++Plane) {
	int ij=Plane + Module*(PlanesInModule+1);
	cout<<"HitBank["<<ij<<"].size()"<<HitBank[ij].size()<<endl;
	for(unsigned int jk=0; jk<HitBank[ij].size(); ++jk) {
	  HitBank[ij][jk]->Print();
	}
      } // for(int Plane=0; Plane<PlanesInModule+1; ++Plane) {
    } // for(int Module=0; Module<NumModules; ++Module) {
    cout<<"......................................................................"<<endl;

  }
  
  //---------------------------------------------------------------------------------------------ascii_output
  int nlayer=0;
  for(int ij=0; ij<500; ++ij) {
    if(HitBank[ij].size()>0){nlayer++;}
  }
  pAnalysis->nLayer=nlayer;
  
  if (pAnalysis->isVisOut==1) {
    //  cout << " pAnalysis->H->NHits "<< pAnalysis->H->NHits<<endl;
    int ll=0;
    for(int ij=0; ij<500; ++ij) {
      if(HitBank[ij].size()>0) {
	pAnalysis->H->NHits= pAnalysis->H->NHits + HitBank[ij].size();
      }
      for(unsigned int jk=0; jk<HitBank[ij].size(); ++jk) {
	pAnalysis->Hp= pAnalysis->H->AddHits(0,0); // add a track object //VALGRIND
	pAnalysis->Hp->TrackType=-1;// Track Type: -1: hits, -2: clulster, -3: triplet, -4: track
	pAnalysis->Hp->HitNum=ll;// Hit Number
	pAnalysis->Hp->ZZ=HitBank[ij][jk]->GetZPlane();
	pAnalysis->Hp->XX=HitBank[ij][jk]->GetXPos();
	pAnalysis->Hp->YY=HitBank[ij][jk]->GetYPos();
	ll++;
      }
    }
  }
  // cout<<"void InoTrackFinder::FormTheHits() completed."<<endl;  
  //--------------------------------------------------------------------------------------------
  // return;
}
//===================================================================================================================================
void InoTrackFinder::FormTheClusters() {
  // cout<<"void InoTrackFinder::FormTheClusters() {..."<<endl;
  // MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
  // For each plane where we stored hits, we look to form 1D clusters
  // of these hits. This is controlled by the IsHitAssoc method in InoCluster.
  // Pointers to the InoCluster objects are stored in plane order in a
  // clusterbank
  
  bool AddingHits;
  //  int i; //GMA14
  
  for(int Module=0; Module<NumModules; ++Module) {
    for(int Plane=0; Plane<PlanesInModule+1; ++Plane) {
      int ij=Plane + Module*(PlanesInModule+1); // GMA14
      for(unsigned int jk=0; jk<HitBank[ij].size(); ++jk) {
	//loop over hits in this plane
	if(HitBank[ij][jk]->GetUID()==0) {
	  //HitBank[ij][jk]::fUID = 0 b4 including the hit in a cluster
	  InoCluster* Clust = new InoCluster(HitBank[ij][jk]);
	  // Make a cluster for the first hit on the plane
	  
	  ClusterBank[HitBank[ij][jk]->GetZPlane()].push_back(Clust);
	  //ClusterBank[HitBank[ij][jk]->GetZPlane()].push_back(Clust);
	  
	  //Change UID from 0 to 1 to show that hit has been added
	  HitBank[ij][jk]->SetUID(1);
	  
	  AddingHits=true;
	  // Loop over other hits on plane to form clusters
	  while(AddingHits==true) {
	    AddingHits=false;
	    for(unsigned int kl=0; kl<HitBank[ij].size(); ++kl) {
	      //cout<<"Check "<<Clust->IsHitAssoc(HitBank[ij][kl])<<endl;
	      
	      //loop over the hits in that plane
	      if(HitBank[ij][kl]->GetUID()==0 && Clust->IsHitAssoc(HitBank[ij][kl])) {
		Clust->AddHit(HitBank[ij][kl]);
		//Add the hit if it satisfy IsHitAssoc()conditn for the cluster
		
		HitBank[ij][kl]->SetUID(1);
		//HitBank[ij][jk]::fUID = 1 when the hit is included as cluster
		AddingHits=true;
	      }
	    }
	  }
	  inoCluster_pointer->InoCluster_list.push_back(Clust);//All the clusters wil be stored in InoCluster_list
	}
      }
      // End loop over hits on plane
    }
    // End loop over planes in module
  }
  // End loop over modules
  //----------------------------------------------------------------------------------------------------
  
  if (pAnalysis->ihist < pAnalysis->nhistmx-1 && pAnalysis->isVisOut>=2) {
    for(int ij=0; ij<500; ++ij) {
      for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
	pAnalysis->gens_list[2][pAnalysis->ihist]->Fill(ClusterBank[ij][jk]->GetXPos(),ClusterBank[ij][jk]->GetYPos(),ClusterBank[ij][jk]->GetZPos()+0.02);
	cout<<"pAnalysis->gens_list[2][pAnalysis->ihst]->Fill();"<<endl;
	
	vectGr  tmpgr;
	tmpgr.x = ClusterBank[ij][jk]->GetXPos();
	tmpgr.y = ClusterBank[ij][jk]->GetYPos();
	tmpgr.z = ClusterBank[ij][jk]->GetZPos()+0.02;
	tmpgr.dx = 0.5*(ClusterBank[ij][jk]->GetEndXPos()-ClusterBank[ij][jk]->GetBegXPos());
	tmpgr.dy = 0.5*(ClusterBank[ij][jk]->GetEndYPos()-ClusterBank[ij][jk]->GetBegYPos());
	tmpgr.dz = 0;
	//	pAnalysis->clus_vect.push_back(tmpgr);
	if (pAnalysis->isVisOut==3) pAnalysis->gens_vect[2].push_back(tmpgr);
      }
    }
  }
  //-----------------------------------------------------------------------------------------------cout
  // Print out the list of cluster
  if(pAnalysis->isXtermOut==1) {
    cout <<"TrkFinder ClusterBank "<<endl;
    for(int ij=0; ij<500; ++ij) {
      for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
	cout<< "InoCluster "<<std::setw(3)<<jk
	    << " pln="    <<std::setw(3)<< ClusterBank[ij][jk]->GetZPlane()
	    << " beX=" <<std::setw(4)<< ClusterBank[ij][jk]->GetBegXStrip()
	    << "  -  " <<std::setw(4)<< ClusterBank[ij][jk]->GetEndXStrip()
	    << " beY=" <<std::setw(4)<< ClusterBank[ij][jk]->GetBegYStrip()
	    << "  -  " <<std::setw(4)<< ClusterBank[ij][jk]->GetEndYStrip()
	    << " Xpos="<<std::setw(8)<< ClusterBank[ij][jk]->GetXPos()
	    << " Ypos="<<std::setw(8)<< ClusterBank[ij][jk]->GetYPos()
	    << " chg=" <<std::setw(8)<< ClusterBank[ij][jk]->GetPulse()
	    << endl;
      }
    }	
  }
  //-------------------------------------------------------------------------------------------ascii_output
  if (pAnalysis->isVisOut==1) {
  }
  //------------------------------------------------------------------------------------------
  // cout<<"void InoTrackFinder::FormTheClusters() completed"<<endl;
  // return;
}

void InoTrackFinder::IDTrkAndShwClusters() {
  // Look at the 1D clusters formed and use simple techniques to get an idea
  // of how track-like or shower-like each cluster is.
  // Then look at all clusters on the plane and see how track-like or shower-like
  // the plane is.
  //  Detector::Detector_t Detector = vldc->GetDetector();
  
  //  int k0;
  int nclust0, nclust1, nhits0, nhits1, ShwAssocNum;
  
  vector<InoCluster*> TempClust0;
  vector<InoCluster*> TempClust1;
  
  //Set TrkFlag to 1, when the pulse value is above threshold (right now set to low value) or fdigit is >1 i.e more thatn two strips are hit
  //MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
  for(int Module=0; Module<NumModules; ++Module) {
    for(int Plane=0; Plane<PlanesInModule+1; ++Plane) {
      int ij=Plane + Module*(PlanesInModule+1);
      for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
	if(ClusterBank[ij][jk]->GetPulse()>PECut || ClusterBank[ij][jk]->GetDigits()>1 ) {
	  ClusterBank[ij][jk]->SetTrkFlag(1);
	}
      }
    }
  }
  
  // Identify the shower-like clusters
  for(int Module=0; Module<NumModules; ++Module) {
    for(int Plane=1; Plane<PlanesInModule+1; ++Plane) {
      int ij=Plane + Module*(PlanesInModule+1);
      
      for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
	InoCluster* Clust0 = ClusterBank[ij][jk];
	nhits0=0; nhits1=0; nclust0=0; nclust1=0;
	TempClust0.clear(); TempClust1.clear();
	
	// Look for shower associations in nearby planes
	// Loop over nearby planes
	for(int kl=-1; kl<2; ++kl) {
	  int kl0=ij+kl;
	  //ij+(2*kl);
	  // Check kl0 is a valid nearby plane
	  if(kl0>Module*(PlanesInModule+1) && kl0<(Module+1)*(PlanesInModule+1)) {
	    // Loop over clusters on plane kl0
	    for(unsigned int kl2=0; kl2<ClusterBank[kl0].size(); ++kl2) {
	      InoCluster* Clust1 = ClusterBank[kl0][kl2];
	      // Check shower-like associations
	      if(Clust0->GetPulse()>.1 && Clust1->GetPulse()>.1) {
		ShwAssocNum=Clust0->IsShwAssoc(Clust1);
		// Store the association data and the clusters
		if(ShwAssocNum>0) {
		  nhits0+=Clust1->GetHitEntries();
		  nclust0+=1;
		  TempClust0.push_back(Clust1);
		}
		if(ShwAssocNum>1) {
		  nhits1+=Clust1->GetHitEntries();
		  nclust1+=1;
		  TempClust1.push_back(Clust1);
		}
	      }
	    }
	    // End loop over clusters on nearby plane
	  }
	}
	// End loop over nearby planes
	// Use the above results to set the shower flags
	// Has important implications for finding tracks embedded in showers
	if(nclust1>2 && nhits1>4) {
	  //GMA optimise it
	  for(unsigned int kl=0; kl<TempClust1.size(); ++kl) {
	    TempClust1[kl]->SetShwFlag(1);
	    // TempClust1[kl]->SetTrkFlag(0);
	    // Clust0->SetTrkFlag(3); //asm: added new
	  }
	}
	
	if(nclust0>4 && nhits0>5) {
	  //GMA optimise it
	  for(unsigned int kl=0; kl<TempClust0.size(); ++kl) {
	    TempClust0[kl]->SetShwFlag(1);
	    // TempClust0[kl]->SetTrkFlag(0);
	    // Clust0->SetTrkFlag(3); //asm: added new
	  }
	}
      }
    }
  }
  
  // Identify track-like and shower-like planes
  double Charge;
  
  for(int Module=0; Module<NumModules; ++Module) {
    for(int Plane=0; Plane<PlanesInModule+1; ++Plane) {
      int ij=Plane + Module*(PlanesInModule+1);
      Charge=0.;
      
      // Get total charge in clusters on plane
      for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
	Charge+=ClusterBank[ij][jk]->GetPulse();
      }
      // asm: Conditions for shower flagging modified:
      // Set the plane flags for the clusters //GMA  should be <3, optimise it
      if(Charge>0.) {
	for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
	  InoCluster* Clust = ClusterBank[ij][jk];
	  if(Clust->GetHitEntries()<4 && (Clust->GetPulse()/Charge)>0.39) {Clust->SetTrkPlnFlag(1);}
	  if(Clust->GetHitEntries()<4 && (Clust->GetPulse()/Charge)<0.40) {Clust->SetTrkPlnFlag(1);}//Clust->SetTrkFlag(0);}
	  if(Clust->GetHitEntries()>3 && Clust->GetHitEntries()<7 && (Clust->GetPulse()/Charge)<0.6)
	    {Clust->SetShwPlnFlag(1);Clust->SetTrkPlnFlag(0);Clust->SetTrkFlag(0);}
	  if(Clust->GetHitEntries()>3 && Clust->GetHitEntries()<7 && (Clust->GetPulse()/Charge)>0.59)
	    {Clust->SetShwPlnFlag(1);Clust->SetTrkPlnFlag(1);Clust->SetTrkFlag(1);}
	  if(Clust->GetHitEntries()>6) {Clust->SetShwPlnFlag(1);Clust->SetTrkPlnFlag(0);Clust->SetTrkFlag(0);}
	  // if(Clust->GetHitEntries()<4) {Clust->SetTrkPlnFlag(1);}//Clust->SetTrkFlag(0); }
	  // if(Clust->GetHitEntries()>3 && Clust->GetHitEntries()<7) {Clust->SetShwPlnFlag(1);Clust->SetTrkPlnFlag(0);}//Clust->SetTrkFlag(0);}
	  // if(Clust->GetHitEntries()>3 && Clust->GetHitEntries()<7) {Clust->SetShwPlnFlag(1);Clust->SetTrkPlnFlag(0);}
	  // if(Clust->GetHitEntries()>6) {Clust->SetShwPlnFlag(1) ;Clust->SetTrkPlnFlag(0);}
	}
      }
    }
  }
  
  //----------------------------------------------------------------------------------------------------------------------cout
  // Print out list of hits and 1D clusters
  // cout <<" InoTrackFinder:ShwTrkFlagging" <<endl;
  
  pAnalysis->inoclust =0;
  for(int ij=0; ij<500; ++ij) {
    double Charge1=0;
    if(ClusterBank[ij].size()){pAnalysis->inoclust +=ClusterBank[ij].size();}
    for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
      Charge1+=ClusterBank[ij][jk]->GetPulse();
    }
    
    if (pAnalysis->isXtermOut==1) {
      for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
	cout<<"shwtrk"
	    << " pln=" << std::setw(3) << ClusterBank[ij][jk]->GetZPlane()
	    << "  beX="<< std::setw(8) << ClusterBank[ij][jk]->GetBegXPos()
	    << "<->"    << std::setw(8) << ClusterBank[ij][jk]->GetEndXPos()
	    << "  beY="<< std::setw(8) << ClusterBank[ij][jk]->GetBegYPos()
	    << "<->"    << std::setw(8) << ClusterBank[ij][jk]->GetEndYPos()
	    << "  nhits="    << std::setw(4) << ClusterBank[ij][jk]->GetHitEntries()
	    << "  puls  ="   << std::setw(10) << ClusterBank[ij][jk]->GetPulse()/Charge1
	    << "  trkf="     << std::setw(3) << ClusterBank[ij][jk]->GetTrkFlag()
	    << "  shwf="     << std::setw(3) << ClusterBank[ij][jk]->GetShwFlag()
	    << "  shwpl="    << std::setw(3) << ClusterBank[ij][jk]->GetShwPlnFlag()
	    << "  trkpl="    << std::setw(4) << ClusterBank[ij][jk]->GetTrkPlnFlag() << endl;
      }
    }// isXtermOut
  }
  //-----------------------------------------------------------------------------------------------------------------ascii_output
  if (pAnalysis->isVisOut==1) {
    int ll=0;
    pAnalysis->H->NClus=inoCluster_pointer->InoCluster_list.size();
    for(int ij=0; ij<500; ++ij) {
      for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
	for(unsigned int kl =0 ;kl<ClusterBank[ij][jk]->GetHitEntries(); ++kl) {
	  pAnalysis->Hp=  pAnalysis->H->AddHits(0,0); // add a track object  //VALGRIND
	  pAnalysis->Hp->TrackType=-2;// Track Type: -1: hits, -2: clulster, -3: trijlet, -4: track
	  pAnalysis->Hp->CluNum=ll;// Hit Number
	  pAnalysis->Hp->ZZ=ClusterBank[ij][jk]->GetZPlane();
	  pAnalysis->Hp->XX=ClusterBank[ij][jk]->GetHit(kl)->GetXPos();
	  pAnalysis->Hp->YY=ClusterBank[ij][jk]->GetHit(kl)->GetYPos();
	}
	ll++;
      }
    }
  }
  //---------------------------------------------------------------------------------------------------------------
  return;
}

void InoTrackFinder::FormTriplets() {
  
  //asm++++++++++++++++++++++++++++++++++++++++changed the conditions for forming triplets+++++++++++++++++++++++
  //asm[][][][][][][][][][][][][][][][][][][][]Refer to the changes in Cluster.cc file also[][][][][][][][][][][]
  
  // Treating U and V views separately, we consider associations between clusters
  // on nearby planes. We look for groups of three clusters, each on different 
  // planes, that can be joined together in a small track segment, or triplet. 
  
  // All possible triplets of the following form are created:
  
  // If we sit on plane 0, and there are no possible hits on a plane marked X, m 
  // means a lower plane number than our origin plane, and p means a higher plane 
  // number.
   
  // m 0 p         The simplest triplet
  // m X 0 p       One gap
  // m 0 X p
  // m X X 0 p     Two gaps
  // m X 0 X p
  // m 0 X X p
  //==================================================================================Look for triplets of the form m 0 p
  int TripAssocNum;
  bool AlternateTriplets;
  bool JoinFlag;
  //  int i; //GMA14
  //MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer; 
  // Begin loop to make the triplets
  for(int Module=0; Module<NumModules; ++Module) {
    for(int Plane=1; Plane<PlanesInModule-1; ++Plane) {
      int ij=Plane + Module*(PlanesInModule+1); //GMA14
      // Loop over clusters on plane. This is our plane 0, the origin.
      
      for(unsigned int k0=0; k0<ClusterBank[ij].size()&&ClusterBank[ij].size()<50; ++k0) {   //asm_310810
	if(ClusterBank[ij][k0]->GetTrkFlag()>0) {
	  JoinFlag=false;
	  
	  // Look for triplets of form m 0 p
	  if((ij-1)>=0 && ClusterBank[ij-1].size()>0 && ClusterBank[ij-1].size()<50 && (ij+1)<500 && ClusterBank[ij+1].size()>0 && ClusterBank[ij+1].size()<50) {
	    //asm_310810       
	    // Look at the clusters on plane ij-1...
	    for(unsigned int kl=0; kl<ClusterBank[ij-1].size(); ++kl) {
	      
	      // ...and the clusters on plane ij+1
	      for(unsigned int kp=0; kp<ClusterBank[ij+1].size(); ++kp) {
		
                if(ClusterBank[ij-1][kl]->GetTrkFlag()>0 && ClusterBank[ij+1][kp]->GetTrkFlag()>0) {
		  
                  // Check the track-like association of the three clusters
                  TripAssocNum=ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-1][kl],ClusterBank[ij+1][kp]);
		  
		  // cout <<"TripAssocNum "<<ij<<" "<<k0<<" "<<kl<<" "<<kp<<" "<<TripAssocNum<<endl;
                  // If the association is good, make the triplet
                  if(TripAssocNum> 2) {        //asm : was 0 
                    InoTrackSegment* seg0 = new InoTrackSegment(ClusterBank[ij-1][kl], ClusterBank[ij][k0], ClusterBank[ij+1][kp]);
		    
                    SegmentBank[ij].push_back(seg0);
                    JoinFlag=true;
		    
                    // Indicate that triplet is exceptionally track-like, by setting flag to 2.
                    if(TripAssocNum>2) 
		      {ClusterBank[ij][k0]->SetTrkFlag(2);}  //asm  : was 1
                  }
                }
	      }
	    }
	  }
	  
	  if((ij+1)<500 && ClusterBank[ij+1].size()>0&& ClusterBank[ij+1].size()<50) {         //asm_310810
	    //=================================================================================m X 0 p
	    
	    if( Plane>1 && (ij-2)>=0 && ClusterBank[ij-2].size()>0 &&ClusterBank[ij-2].size()<50) {
	      // Look at the clusters on plane ij+1...
	      for(unsigned int kp=0; kp<ClusterBank[ij+1].size()&&ClusterBank[ij+1].size()<50; ++kp) {
		
		// ...and the clusters on plane ij-2
		for(unsigned int kl=0; kl<ClusterBank[ij-2].size(); ++kl) {
		  
		  if(ClusterBank[ij-2][kl]->GetTrkFlag()>0 && ClusterBank[ij+1][kp]->GetTrkFlag()>0) {
		    
		    // Check the track-like association of the three clusters
		    TripAssocNum=ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-2][kl],ClusterBank[ij+1][kp]);
		    
		    // cout<<"ij m X 0 p "<< ij<<" "<<k0<<" "<<kl<<" "<<kp<<" "<<TripAssocNum<<endl;
		    
		    // If the association is good, check whether we can make the triplet
		    if(TripAssocNum>2) {                          //asm: was >0
		      AlternateTriplets=false;
		      
		      // Check to see if there is also an alternate triplet possibility with clusters on plane ij-2
		      // i.e. Want to make sure that the X in "m X 0 p" is correct.
		      
		      if(AlternateTriplets==false && ClusterBank[ij-1].size()>0 &&  ClusterBank[ij-1].size()<50) {
                        // Look at the clusters on plane ij-1
			for(unsigned int ktmp=0; ktmp<ClusterBank[ij-1].size(); ++ktmp) {
			  if(AlternateTriplets==false && ClusterBank[ij-1][ktmp]->GetTrkFlag()>0
			     && (ClusterBank[ij-1][ktmp]->IsTrkAssoc(ClusterBank[ij-2][kl],ClusterBank[ij][k0])>2            //ij-2,ij-1,ij
				 &&         ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-1][ktmp],ClusterBank[ij+1][kp])>2)              //ij-1,ij,ij+1 
			     ) {
			    AlternateTriplets=true;
			  }
			}
                      }
		      
		      // If everything is good, make the triplet
		      if(AlternateTriplets==false) {
			for(int jk=-1; jk<1; ++jk) {
			  InoTrackSegment* seg0 = new InoTrackSegment(ClusterBank[ij-2][kl],ClusterBank[ij][k0],ClusterBank[ij+1][kp]);
			  
			  SegmentBank[ij+jk].push_back(seg0);
			  JoinFlag=true;
			}
			
			// Indicate that triplet is exceptionally track-like, by setting flag to 2.
			if(TripAssocNum>2) 
			  {ClusterBank[ij][k0]->SetTrkFlag(2);}    //asm: was >1                                     
		      }
                    }
		  }
		}
	      }
	    }
	    //=================================================================== m X X 0 p
	    
	    //	    cout <<"trkfindertriplet ij21 "<<endl;
	    if( Plane>2 && (ij-3)>=0 && ClusterBank[ij-3].size()>0&& ClusterBank[ij-3].size()<50) {
	      
              // Look at the clusters on plane ij+1...
	      for(unsigned int kp=0; kp<ClusterBank[ij+1].size()&&ClusterBank[ij+1].size()<50; ++kp) {
		
                // ...and look at the clusters on plane ij-3
		for(unsigned int kl=0; kl<ClusterBank[ij-3].size()&& ClusterBank[ij-3].size()<50 ; ++kl) {
		  
                  if(ClusterBank[ij-3][kl]->GetTrkFlag()>0 && ClusterBank[ij+1][kp]->GetTrkFlag()>0) {
		    
                    // Check the track-like association of the three clusters
                    TripAssocNum=ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-3][kl],ClusterBank[ij+1][kp]);
		    
                    // If the association is good, check whether we can make the triplet
                    if(TripAssocNum>2) {
                      AlternateTriplets=false;
                      
                      // Check to see if there are alternate triplet possibilities with clusters on plane ij-2 or ij-4
                      // i.e. Want to make sure that the Xs in "m X X 0 p" are correct.
		      
                      // Check first X
                      if(AlternateTriplets==false && ClusterBank[ij-1].size()>0 && ClusterBank[ij-1].size()<50) {
                        // Look at any clusters on plane ij-1
                        for(unsigned int ktmp=0; ktmp<ClusterBank[ij-1].size(); ++ktmp) {
                          if(AlternateTriplets==false && ClusterBank[ij-1][ktmp]->GetTrkFlag()>0
                             &&( ClusterBank[ij-1][ktmp]->IsTrkAssoc(ClusterBank[ij-3][kl],ClusterBank[ij][k0])>2  //ij-3,ij-1,ij
				 && ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-1][ktmp],ClusterBank[ij+1][kp])>2) // ij-1,ij,ij+1
			     ) {
			    AlternateTriplets=true;
			  }
                        }
                      }
                      
                      // Check second X, plane ij-2
                      if(AlternateTriplets==false && ClusterBank[ij-2].size()>0&&ClusterBank[ij-2].size()<50) {
                        // Look at any clusters on plane ij-2
                        for(unsigned int ktmp=0; ktmp<ClusterBank[ij-2].size(); ++ktmp) {
                          if(AlternateTriplets==false && ClusterBank[ij-2][ktmp]->GetTrkFlag()>0
                             &&( ClusterBank[ij-2][ktmp]->IsTrkAssoc(ClusterBank[ij-3][kl],ClusterBank[ij][k0])>2 //ij-3,ij-2,ij
				 && ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-2][ktmp],ClusterBank[ij+1][kp])>2) //ij-2,ij,ij+1
			     ) {
			    AlternateTriplets=true;
			  }
                        }
                      }
                      
                      // Check both Xs together, planes ij-1 and ij-2
                      if(AlternateTriplets==false && ClusterBank[ij-1].size()>0 && ClusterBank[ij-1].size()<50&& ClusterBank[ij-2].size()>0 && ClusterBank[ij-2].size()<50) {
                        // Look again at any clusters on plane ij-2...
                        for(unsigned int ktmp=0; ktmp<ClusterBank[ij-2].size(); ++ktmp) {
			  
                          // ...and look again any clusters on plane ij-1
                          for(unsigned int ktmp1=0; ktmp1<ClusterBank[ij-1].size(); ++ktmp1) {
                            if(AlternateTriplets==false 
                               && ClusterBank[ij-2][ktmp]->GetTrkFlag()>0 && ClusterBank[ij-1][ktmp1]->GetTrkFlag()>0
                               &&( ClusterBank[ij-2][ktmp]->IsTrkAssoc(ClusterBank[ij-3][kl],ClusterBank[ij-1][ktmp1])>2 // ij-3,ij-2,ij-1
				   && ClusterBank[ij-1][ktmp1]->IsTrkAssoc(ClusterBank[ij-2][ktmp],ClusterBank[ij][k0])>2   // ij-2,ij-1,ij
				   && ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-1][ktmp1],ClusterBank[ij+1][kp])>2)    // ij-1,ij,ij+1
			       ) {
			      AlternateTriplets=true;
			    }
                          }
                        }
                      }
                      
                      // If everything is good, make the triplet
                      if(AlternateTriplets==false) {
                        // Store segment on empty planes too          
                        for(int jk=0; jk<1; ++jk){            //        for(int jk=-2; jk<1; ++jk) {
			  InoTrackSegment* seg0 = new InoTrackSegment(ClusterBank[ij-3][kl],ClusterBank[ij][k0],ClusterBank[ij+1][kp]);
                          SegmentBank[ij+jk].push_back(seg0);
                          JoinFlag=true;
                        }
			
                        // Indicate that triplet is exceptionally track-like, by setting flag to 2.
                        if(TripAssocNum>2) {ClusterBank[ij][k0]->SetTrkFlag(2);}
                      }
                    }
                  }
		}
	      }
	    }
	  } 
	  // Look for triplets of form m 0 <-> p
	  if((ij-1)>0 && ClusterBank[ij-1].size()>0&&ClusterBank[ij-1].size()<50) {
	    
	    //==================================================================================================m 0 X p        
	    
	    // Look for triplets of form m 0 X p  
	    if( Plane<PlanesInModule-1 && (ij+2)<500 && ClusterBank[ij+2].size()>0&&ClusterBank[ij+2].size()<50) {
	      
	      // Look at the clusters on plane ij-1...
	      for(unsigned int kl=0; kl<ClusterBank[ij-1].size(); ++kl) {
		
		// ... and the clusters on plane ij+2
		for(unsigned int kp=0; kp<ClusterBank[ij+2].size(); ++kp) {
		  
		  if(ClusterBank[ij-1][kl]->GetTrkFlag()>0 && ClusterBank[ij+2][kp]->GetTrkFlag()>0) {
		    
                    // Check the track-like association of the three clusters
		    TripAssocNum=ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-1][kl],ClusterBank[ij+2][kp]);
		    
		    // If the association is good, check whether we can make the triplet
		    if(TripAssocNum>2) {
		      AlternateTriplets=false;
		      
		      // Check to see if there is also an alternate triplet possibility with clusters on plane i+2
		      // i.e. Want to make sure that the X in "m 0 X p" is correct.
		      
		      if(AlternateTriplets==false && ClusterBank[ij+1].size()>0) {
			// Look at the clusters on plane ij+1
			for(unsigned int ktmp=0; ktmp<ClusterBank[ij+1].size(); ++ktmp) {
			  if(AlternateTriplets==false && ClusterBank[ij+1][ktmp]->GetTrkFlag()>0
			     &&( ClusterBank[ij+1][ktmp]->IsTrkAssoc(ClusterBank[ij][k0],ClusterBank[ij+2][kp])>2 
				 &&     ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-1][kl],ClusterBank[ij+1][ktmp])>2) 
			     ) {
			    AlternateTriplets=true;
			  }
			}
		      }
                      
		      
		      // If everything is good, make the triplet
		      if(AlternateTriplets==false) {
			// Store segment on empty plane too
			for(int jk=0; jk<2; ++jk) {  // for(int jk=0; jk<2; ++jk) {  //GMA why two layers, whereas previous case one layer
			  InoTrackSegment* seg0 = new InoTrackSegment(ClusterBank[ij-1][kl],ClusterBank[ij][k0],ClusterBank[ij+2][kp]);
			  
			  SegmentBank[ij+jk].push_back(seg0);
			  JoinFlag=true;
                        }
			
			// Indicate that triplet is exceptionally track-like, by setting flag to 2.
			if(TripAssocNum>2) {
			  ClusterBank[ij][k0]->SetTrkFlag(2);
			}
		      }
		    }
		  }
		}
	      }
	    }
	    //==========================================================================================m 0 X X p
	    // Look for triplets of form m 0 X X p  
	    if(Plane<PlanesInModule-2 && (ij+3)<500 && ClusterBank[ij+3].size()>0&& ClusterBank[ij+3].size()<50) {
	      
	      // Look at the clusters on plane ij-1...
	      for(unsigned int kl=0; kl<ClusterBank[ij-1].size() && ClusterBank[ij-1].size()<50; ++kl) {
		
		// ...and look at the clusters on plane ij+3
		for(unsigned int kp=0; kp<ClusterBank[ij+3].size(); ++kp) {
		  
		  if(ClusterBank[ij-1][kl]->GetTrkFlag()>0 && ClusterBank[ij+3][kp]->GetTrkFlag()>0) {
		    // Check the track-like association of the three clusters
		    TripAssocNum=ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-1][kl],ClusterBank[ij+3][kp]);
		    // If the association is good, check whether we can make the triplet
		    
		    if(TripAssocNum>2) {
		      AlternateTriplets=false;
		      
		      // Check to see if there are alternate triplet possibilities with clusters on plane ij+1 or ij+2
		      // i.e. Want to make sure that the Xs in "m 0 X X p" are correct.
		      
		      // Check first X, plane i+1
		      //		      cout <<"Alter "<<(int)AlternateTriplets<<" "<<ClusterBank[ij+1].size()<<endl;
		      if(AlternateTriplets==false && ClusterBank[ij+1].size()>0 && ClusterBank[ij+1].size()<50) {
			// Look at any clusters on plane ij+1
			for(unsigned int ktmp=0; ktmp<ClusterBank[ij+1].size(); ++ktmp) {
			  if(AlternateTriplets==false && ClusterBank[ij+1][ktmp]->GetTrkFlag()>0
			     &&( ClusterBank[ij+1][ktmp]->IsTrkAssoc(ClusterBank[ij][k0],ClusterBank[ij+3][kp])>2 
				 &&     ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-1][kl],ClusterBank[ij+1][ktmp])>2) 
			     ) {
			    AlternateTriplets=true;
			  }
			}
		      }
		      
		      // Check second X, plane ij+2
		      if(AlternateTriplets==false && ClusterBank[ij+2].size()>0 && ClusterBank[ij+2].size()<50) {
			// Look at any clusters on plane ij+2
			for(unsigned int ktmp=0; ktmp<ClusterBank[ij+2].size(); ++ktmp) {
			  if(AlternateTriplets==false && ClusterBank[ij+2][ktmp]->GetTrkFlag()>0
			     &&( ClusterBank[ij+2][ktmp]->IsTrkAssoc(ClusterBank[ij][k0],ClusterBank[ij+3][kp])>2 
				 &&     ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-1][kl],ClusterBank[ij+2][ktmp])>2) ) {
			    AlternateTriplets=true;
			  }
			}
		      }
		      
		      // Check both Xs together, planes ij+1 and ij+2
		      if(AlternateTriplets==false && ClusterBank[ij+1].size()>0  && ClusterBank[ij+1].size()<50 && ClusterBank[ij+2].size()>0&& ClusterBank[ij+2].size()<50) {
			// Look again at any clusters on plane ij+1...
			for(unsigned int ktmp=0; ktmp<ClusterBank[ij+1].size(); ++ktmp) {
			  
			  // ...and look again at any clusters on plane i+2
			  for(unsigned int ktmp1=0; ktmp1<ClusterBank[ij+2].size(); ++ktmp1) {
			    if(AlternateTriplets==false 
			       && ClusterBank[ij+1][ktmp]->GetTrkFlag()>0 &&ClusterBank[ij+2][ktmp1]->GetTrkFlag()>0
			       &&( ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-1][kl],ClusterBank[ij+1][ktmp])>2
				   &&     ClusterBank[ij+1][ktmp]->IsTrkAssoc(ClusterBank[ij][k0],ClusterBank[ij+2][ktmp1])>2
				   &&     ClusterBank[ij+2][ktmp1]->IsTrkAssoc(ClusterBank[ij+1][ktmp],ClusterBank[ij+3][kp])>2)
			       ) {
			      AlternateTriplets=true;
			    }
			  }
			}
		      }
		      
		      // If everything is good, make the triplet
		      if(AlternateTriplets==false) {
			// Store segment on empty planes too
			for(int jk=0; jk<3; ++jk) {
                          InoTrackSegment* seg0 = new InoTrackSegment(ClusterBank[ij-1][kl],ClusterBank[ij][k0],ClusterBank[ij+3][kp]);
                          
                          SegmentBank[ij+jk].push_back(seg0);
                          JoinFlag=true;
			}
			// Indicate that triplet is exceptionally track-like, by setting flag to 2.
			if(TripAssocNum>2) {
			  ClusterBank[ij][k0]->SetTrkFlag(2);
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  
	  //==================================================================================m X 0 X p
	  // Look for triplets of form m X 0 X p
	  if(Plane>1 && Plane<PlanesInModule-1 && (ij-2)>=0 && ClusterBank[ij-2].size()>0 &&ClusterBank[ij-2].size()<50 && (ij+2)<500 && ClusterBank[ij+2].size()>0&&ClusterBank[ij+2].size()<50) {
	    
	    // Look at clusters on planes ij-2 and ij+2
	    for(unsigned int kl=0; kl<ClusterBank[ij-2].size(); ++kl) {
	      for(unsigned int kp=0; kp<ClusterBank[ij+2].size(); ++kp) {
		
		if(ClusterBank[ij-2][kl]->GetTrkFlag()>0 && ClusterBank[ij+2][kp]->GetTrkFlag()>0) {
		  // Check the track-like association of the three clusters
		  TripAssocNum=ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-2][kl],ClusterBank[ij+2][kp]);
		  
		  // If the association is good, check whether we can make the triplet
		  if(TripAssocNum>2) {
		    AlternateTriplets=false;
		    
		    // Check to see if there are alternate triplet possibilities with clusters on plane ij-2 or ij+2
		    // i.e. Want to make sure that the Xs in "m X 0 X p" are correct.
		    
		    // Check first X, plane ij-1
		    if(AlternateTriplets==false && ClusterBank[ij-1].size()>0 && ClusterBank[ij-1].size()<50) {
		      // Look at any clusters on plane ij-1
		      for(unsigned int ktmp=0; ktmp<ClusterBank[ij-1].size(); ++ktmp) {
			if(AlternateTriplets==false && ClusterBank[ij-1][ktmp]->GetTrkFlag()>0
			   &&( ClusterBank[ij-1][ktmp]->IsTrkAssoc(ClusterBank[ij-2][kl],ClusterBank[ij][k0])>2 
			       &&     ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-1][ktmp],ClusterBank[ij+2][kp])>2)) {
			  AlternateTriplets=true;
			}
		      }
		    }
		    
		    // Check second X, plane ij+1
		    if(AlternateTriplets==false && ClusterBank[ij+1].size()>0 && ClusterBank[ij+1].size()<50) {
		      // Look at any clusters on plane ij+1
		      for(unsigned int ktmp=0; ktmp<ClusterBank[ij+1].size(); ++ktmp) {
			if(AlternateTriplets==false  && ClusterBank[ij+1][ktmp]->GetTrkFlag()>0
			   &&( ClusterBank[ij+1][ktmp]->IsTrkAssoc(ClusterBank[ij][k0],ClusterBank[ij+2][kp])>2 
			       &&    ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-2][kl],ClusterBank[ij+1][ktmp])>2)) {
			  AlternateTriplets=true;
			}
		      }
		    }
                    
		    // Check both Xs together, planes ij-1 and ij+1
		    if(AlternateTriplets==false && ClusterBank[ij-1].size()>0 && ClusterBank[ij-1].size()<50 && ClusterBank[ij+1].size()>0  && ClusterBank[ij+1].size()<50) {
		      
                      // Loop again at any clusters on plane ij-1
		      for(unsigned int ktmp=0; ktmp<ClusterBank[ij-1].size(); ++ktmp) {
			
			// Look again at any clusters on plane i+1
			for(unsigned int ktmp1=0; ktmp1<ClusterBank[ij+1].size(); ++ktmp1) {
			  if(AlternateTriplets==false 
			     && ClusterBank[ij-1][ktmp]->GetTrkFlag()>0 && ClusterBank[ij+1][ktmp1]->GetTrkFlag()>0
			     &&( ClusterBank[ij-1][ktmp]->IsTrkAssoc(ClusterBank[ij-2][kl],ClusterBank[ij][k0])>2 
				 &&     ClusterBank[ij][k0]->IsTrkAssoc(ClusterBank[ij-1][ktmp],ClusterBank[ij+1][ktmp1])>2 
				 &&     ClusterBank[ij+1][ktmp1]->IsTrkAssoc(ClusterBank[ij][k0],ClusterBank[ij+2][kp])>2)) {
			    AlternateTriplets=true;
			  }
			}
		      }
		    }
		    
		    // If everything is good, make the triplet
		    if(AlternateTriplets==false) {
		      // Store segment on empty planes too
		      for(int jk=0; jk<2; ++jk) { //for(int jk=-1; jk<2; ++jk) {
			InoTrackSegment* seg0 = new InoTrackSegment(ClusterBank[ij-2][kl],ClusterBank[ij][k0],ClusterBank[ij+2][kp]);
                        
			SegmentBank[ij+jk].push_back(seg0);
			JoinFlag=true;
		      }
                      // Indicate that triplet is exceptionally track-like, by setting flag to 2.
		      if(TripAssocNum>2) {
			ClusterBank[ij][k0]->SetTrkFlag(2);
		      }
		    }               
		  }
                }
	      }
	    }
	  }
	  
	  //	  cout <<"trkfindertriplet ij6 "<<endl;
	  // Set the exceptionally track-like flag for lone hits that form part of a triplet
	  //asmQ what it this part for?
	  if(JoinFlag==true && 
	     (ClusterBank[ij][k0]->GetEndXPos()==ClusterBank[ij][k0]->GetBegXPos()) &&
	     (ClusterBank[ij][k0]->GetEndYPos()==ClusterBank[ij][k0]->GetBegYPos())) {
            ClusterBank[ij][k0]->SetTrkFlag(2);
	  }
	}
	
      } // End loop over clusters on plane
    } // End loop over planes in module
  } // End loop over modules
  
  //---------------------------------------------------------------------------------------------------cout: cluster
  /*
  // Print out list of hits and 1D clusters
  //asmQ: is this cout required here? 
  if (pAnalysis->isXtermOut==1){
    cout << " InoTrackFinder : *** 2ND LIST OF CLUSTERS *** " << endl;
    for(int ij=0; ij<500; ++ij){
      for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
	cout 
	  << " plane="    << ClusterBank[ij][jk]->GetZPlane() 
	  << " beX"  << ClusterBank[ij][jk]->GetBegXPos() 
	  << " - "  << ClusterBank[ij][jk]->GetEndXPos() 
	  << " beY"  << ClusterBank[ij][jk]->GetBegYPos() 
	  << " - "  << ClusterBank[ij][jk]->GetEndYPos() 
	  << " trkflag="  << ClusterBank[ij][jk]->GetTrkFlag()
	  << " shwpln="   << ClusterBank[ij][jk]->GetShwPlnFlag() 
	  << " trkpln="   << ClusterBank[ij][jk]->GetTrkPlnFlag()
	  << " shwflag="  << ClusterBank[ij][jk]->GetShwFlag() << endl;
      }
    }
  }
*/
//----------------------------------------------------------------------------------------------------------cout
  // Print out list of triplets
  if (pAnalysis->isXtermOut==1){
    cout << " InoTrackFinder : *** LIST OF TRIPLETS *** " << endl;
    for(int ij=0; ij<500; ++ij) {
      for(unsigned int jk=0;jk<SegmentBank[ij].size(); ++jk) {
	cout << " SegmentBank " << jk << ", Pln " << ij <<endl;
	
	InoCluster* clr0 = SegmentBank[ij][jk]->GetCluster(0);
	InoCluster* clr1 = SegmentBank[ij][jk]->GetCluster(1);
	InoCluster* clr2 = SegmentBank[ij][jk]->GetCluster(2);
	
	cout
	  << "(" << clr0->GetZPlane() << "," << clr0->GetBegXPos() << "->" << clr0->GetEndXPos() << ") "
	  << "(" << clr1->GetZPlane() << "," << clr1->GetBegXPos() << "->" << clr1->GetEndXPos() << ") "
	  << "(" << clr2->GetZPlane() << "," << clr2->GetBegXPos() << "->" << clr2->GetEndXPos() << ") "
	  << " (" << clr0->GetZPlane() << "," << clr0->GetBegYPos() << "->" << clr0->GetEndYPos() << ") "
	  << "(" << clr1->GetZPlane() << "," << clr1->GetBegYPos() << "->" << clr1->GetEndYPos() << ") "
	  << "(" << clr2->GetZPlane() << "," << clr2->GetBegYPos() << "->" << clr2->GetEndYPos() << ") " << endl;
	//   cout<<endl; 
      }
    }
  }//(isXtermOut)
  
  //-----------------------------------------------------------------------------------------------------ascii_output
  if (pAnalysis->isVisOut==1){
    
    int n_clushit;               
    int ll=pAnalysis->H->NTrips; 
    for( int ij=0; ij<500; ++ij ) {
      for(unsigned int jk=0;jk<SegmentBank[ij].size(); ++jk) { 
	n_clushit =0;
	for (unsigned int kl=0; kl<SegmentBank[ij][jk]->ClustersInSegment.size(); ++kl){
	  n_clushit +=  SegmentBank[ij][jk]->GetCluster(kl)->HitsInCluster.size(); 
	}
	
	//           InoCluster* clr0 = SegmentBank[ij][jk]->GetCluster(k);
/*	pAnalysis->ascii_output << std::setw(12)<< "1"
				<< std::setw(12)<< "-3"
				<< std::setw(12)<< ll 
				<< std::setw(12)<< jk
				<< std::setw(14) << n_clushit
				<< std::setw(14)<< "-3"
				<< std::setw(14)<< "-3"
				<< std::setw(14)<< "-3"    
				<< endl;
*/
	int hh=0;     
	for (unsigned int kl=0; kl<SegmentBank[ij][jk]->ClustersInSegment.size(); ++kl) {
	  InoCluster* clr0 = SegmentBank[ij][jk]->GetCluster(kl);
	  /*   pAnalysis->ascii_output << std::setw(12)<< "1"
	      <<std::setw(12) << "-3" 
	      << std::setw(12)<< kl
	      << std::setw(12)
	      << std::setw(14)<< "-2"
	      << std::setw(14)<< "-2"
	      << std::setw(14)<< "-2"
	      << endl;  */   
          for (unsigned int lm=0;lm< clr0->HitsInCluster.size(); ++lm){ 
/*	    pAnalysis->ascii_output << std::setw(12)<< "1" 
				    << std::setw(12)<< "-3"
				    << std::setw(12)<< ll
				    << std::setw(12)<< hh
				    << std::setw(12)<< lm 
				    << std::setw(14)<<clr0->HitsInCluster[lm]->GetZPlane()
				    << std::setw(14)<<clr0->HitsInCluster[lm]->GetXPos()
				    << std::setw(14)<<clr0->HitsInCluster[lm]->GetYPos() 
				    << endl;
*/
	    pAnalysis->H->NTrips= pAnalysis->H->NTrips +1; //Number of Triplet events
	    pAnalysis->Hp= pAnalysis->H->AddHits(0,0); // add a track object //VALGRIND
	    pAnalysis->Hp->TrackType=-3;// Track Type: -1: hits, -2: clulster, -3: triplet, -4: track
	    pAnalysis->Hp->TriNum=ll;// Hit Number
	    pAnalysis->Hp->ZZ=clr0->HitsInCluster[lm]->GetZPlane();
	    pAnalysis->Hp->XX=clr0->HitsInCluster[lm]->GetXPos();
	    pAnalysis->Hp->YY=clr0->HitsInCluster[lm]->GetYPos();
	    
	    
	    hh++;
	  }
	}
	ll++;
      }
    }
  }
  //asm: here we display all hits in that cluster, we can move to clusters in the triplet.	  
  //--------------------------------------------------------------------------------------------
  return;
}

void InoTrackFinder::FindAllAssociations() {
  // For each triplet formed, we check nearby triplets and examine if
  // the triplets are associated. For each triplet, we simply find the
  // nearby other triplets that have a compatible beginning and/or end
  // position/direction. Later we look more that one triplet away, for
  // finding strings of "preferred associations" that might indicate a
  // track.
  
  //  int i = 0;
  //MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
  
  // Loop over the planes
  for(int Module = 0; Module < NumModules; ++Module) {
    for(int Plane = 1; Plane < PlanesInModule-1; ++Plane) {
      int ij = Plane + Module*(PlanesInModule + 1);
      if((ij + 1) < 500 && SegmentBank[ij].size() > 0 && SegmentBank[ij+1].size() > 0) {
	// Look at the segments on plane ij
	for(unsigned int k0 = 0; k0 < SegmentBank[ij].size(); ++k0) {
	  // For each of these, look at the segments on plane ij+1
	  for(int jk=1;jk<2;jk++) {
	    for(unsigned int kp=0; kp<SegmentBank[ij+jk].size(); ++kp) {
	      // Check associations.
	      if(SegmentBank[ij][k0]->IsAssoc(SegmentBank[ij+jk][kp])) {
		// If segments are associated, we add to list of ALL possible joins
		SegmentBank[ij][k0]->AddAssocSegToEnd(SegmentBank[ij+jk][kp]);
		SegmentBank[ij+jk][kp]->AddAssocSegToBeg(SegmentBank[ij][k0]);
	      }
	    }
	  }
	}
      }
    }// End loop over planes in module
  }// End loop over modules
  
  //asmS:---------------------------------------------------------------------------------------------------------
  // In the above code we are not inplementing the IsAssoc function with its complete strength,
  // let us see how, Every triplet that is formed is always kept in all othe layer except the last layer, now when
  // I look for triplets in layer i and i+1 it will always fall in the first condition of IsAssoc function where the
  // triplets overlap, but IsAssoc also look for association between triplets which do not have any layer in common.
  // which the above code we will never chck for such association.
  // May be that part of function IsAssoc is used in FormPreferedJoins function where the we check for association
  // when the segment is seperated by a gap.
  // infact over here if we select the loop variable properly we may not need to multiple copies of segment in empty
  // layers. After discussing with GM, it was decided that we will store the copy of the triplet in the empty layers
  // as well, but I fear if this in turm will make multiple copy of all the segments. I need to look into how this works
  //--------------------------------------------------------------------------------------------------------------
  
  // Print out list of all associations
  //  cout << " InoTrackFinder : *** LIST OF ALL SEG ASSOCIATIONS *** " << endl;
  if (pAnalysis->isXtermOut==1)	{
    for(int ij=0; ij<500; ++ij) {
      for(unsigned int jk=0;jk<SegmentBank[ij].size(); ++jk) {
	InoTrackSegment* segment = SegmentBank[ij][jk];
	cout << "--- Plane " << ij << ", Seg Number " << jk
	     << ", BegPlane " << segment->GetBegZPlane()
	     << ", EndPlane " << segment->GetEndZPlane()
	     << ", BegXPos " << segment->GetBegXPos()
	     << ", EndXPos " << segment->GetEndXPos()
	     << ", BegYPos " << segment->GetBegYPos()
	     << ", EndYPos " << segment->GetEndYPos()
	     << ", Entries " << segment->GetEntries()
	     <<endl;
	
	cout << " Beg: " <<endl;
	for(unsigned int kl=0; kl<segment->GetNAssocSegBeg(); ++kl) {
	  InoTrackSegment* segtmp = segment->GetAssocSegBeg(kl);
	  cout << "  Assoc number=" << kl
	       << "  begpln: " << segtmp->GetBegZPlane()
	       << "  begxpos: " << segtmp->GetBegXPos()
	       << "  begypos: " << segtmp->GetBegYPos()
	       << "  endpln: " << segtmp->GetEndZPlane()
	       << "  endxpos: " << segtmp->GetEndXPos()
	       << "  endypos: " << segtmp->GetEndYPos() 
	       << ", Entries " << segtmp->GetEntries()
	       << endl;
	}
	
	cout << " End: "<<endl;
	for(unsigned int kl=0; kl<segment->GetNAssocSegEnd(); ++kl) {
	  InoTrackSegment* segtmp = segment->GetAssocSegEnd(kl);
	  cout << "  Assoc number=" << kl
	       << "  begpln: " << segtmp->GetBegZPlane()
	       << "  begxpos: " << segtmp->GetBegXPos()
	       << "  begypos: " << segtmp->GetBegYPos()
	       << "  endpln: " << segtmp->GetEndZPlane()
	       << "  endxpos: " << segtmp->GetEndXPos()
	       << "  endypos: " << segtmp->GetEndYPos()<< endl;
	}
      }
    }
  }//if isXterm
  return;
}

void InoTrackFinder::FindPreferredJoins() {
  // Having made all possible associations between triplets above, let 
  // us select those preferred associations that are most track-like
  
  // For a given triplet, we know which triplets are associated with the
  // beginning and which are associated at the end. If there are associations
  // between these beginning and end triplets too, then we are quite likely 
  // to be considering track-like segments.
  
  //  int i;
  //MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
  vector<InoTrackSegment*> mTempSeg;
  vector<InoTrackSegment*> pTempSeg;
  
  bool JoinFlag; bool mJnFlag; bool pJnFlag;
  
  double GradientXTolerance=8.;
  double GradientYTolerance=8.;
  
  // Loop over the planes
  for(int Module = 0; Module < NumModules; ++Module) {
    for(int Plane = 1; Plane < PlanesInModule-1; ++Plane) {
      int ij = Plane + Module*(PlanesInModule + 1); //GMA14 "ij" should be local variable
      
      // SET THE TEMPORARY TRACK FLAGS (TmpTrkFlag)
      // Consider previous adjacent triplets   ( ij-1 <-> ) ij <-> ij+1
      //                                                  Seg
      // If there are associations at both ends of triplet, set triplet's 
      // temporary track flag to 1. 
      
      // If one of beginning associated triplets also starts at lower plane 
      // than current triplet, set the triplet's temporary track flag to 2
      
      mTempSeg.clear();
      //---------------------
      // Loop over the triplets centered on plane i
      //      cout <<"TrkFinder : preferjoin "<< ij <<" "<<SegmentBank[ij].size()<<endl;
      for(unsigned int jk=0; jk<SegmentBank[ij].size(); ++jk) {
	InoTrackSegment* Seg = SegmentBank[ij][jk];
	if(Seg->GetNAssocSegEnd()>0) {
	  //Check if any Associated triplets at end
	  mTempSeg.push_back(Seg);                  //if yes, store it in mTempSeg
	  
	  if(Seg->GetNAssocSegBeg()>0) {
	    //Check if any associated triplet at begin
	    Seg->SetTmpTrkFlag(1);                //When triplet has assocition on both the ends
	    
	    // Loop over the associations at beginning and set TmpTrackFlag to 2. if 
	    //there at least one SegBeg that starts from z plane lesser than the triplets Begplane
	    for(unsigned int kl=0; kl<Seg->GetNAssocSegBeg(); ++kl) {
	      if (Seg->GetTmpTrkFlag()<2 && (Seg->GetAssocSegBeg(kl)->GetBegZPlane() < Seg->GetBegZPlane())) {
		Seg->SetTmpTrkFlag(2);
		break;
	      }
	    }
	  }
	}
      }
      //-----------------------
      // Subsequent adjacent triplets   ij <-> ij+1 ( <-> ij+2 )
      //                                      Segp
      // If there are associations at both ends of triplet centered on 
      // plane ij+1, set triplet's temporary track flag to 1.
      
      // If one of the end associated triplets also ends at higher plane
      // than ij+1 triplet, set ij+1 triplet's temporary track flag to 2
      
      pTempSeg.clear();
      
      // Loop over the triplets centered on plane ij+1
      for(unsigned int jk = 0; jk < SegmentBank[ij+1].size(); ++jk) {
	InoTrackSegment* Segp = SegmentBank[ij+1][jk];
	
	// Check to see if there are associations at beginning
	if(Segp->GetNAssocSegBeg()>0) {
	  //Check if any Associated triplets at end
	  pTempSeg.push_back(Segp);        //if yes, store it in mTempSeg
	  
	  // Check to see if there are associations at end
	  if(Segp->GetNAssocSegEnd()>0) {
	    //Check if any associated triplet at begin
	    // We now know that triplet has possible associations at both ends
	    Segp->SetTmpTrkFlag(1);
	    
	    // Loop over associations at end and set temp track flags
	    // for segments on plane ij+1
	    for(unsigned int kl=0; kl<Segp->GetNAssocSegEnd(); ++kl) {
	      if(Segp->GetTmpTrkFlag()<2 && (Segp->GetAssocSegEnd(kl)->GetEndZPlane() > Segp->GetEndZPlane()) ) { 
		//asm: if BegAssoc seg is same as the triplet it self it will not do this flagging
		Segp->SetTmpTrkFlag(2); break;
	      }
	    }
	  }
	}
      }
      //----------------------------------------------------------
      // Look for preferred joins   ( ij-1 <-> ) ij <-> ij+1
      // Loop the over segments on plane ij+1, which we know have possible beginning associations
      for(unsigned int jk=0; jk<pTempSeg.size(); ++jk) {
	JoinFlag=false;
	
	// Loop over these beginning associations. If one beginning associated
	// triplet already has temp track flag 2, set JoinFlag to true.
	
	for(unsigned int kl=0; kl<pTempSeg[jk]->GetNAssocSegBeg(); ++kl) {
	  if(pTempSeg[jk]->GetAssocSegBeg(kl)->GetTmpTrkFlag()>1) {
	    JoinFlag=true; break;
	  }
	}
	
	// If JoinFlag is true, and a beginning assoc segment has temp track flag 1 (i.e. possible
	// beginning and end associations), set the temp track flag to 2 for this beg assoc segment
	if(JoinFlag==true) {
	  // Loop over beginning associations again and set temp track flags
	  // for segments on plane ij
	  for(unsigned int kl=0; kl<pTempSeg[jk]->GetNAssocSegBeg(); ++kl) {
	    if(pTempSeg[jk]->GetAssocSegBeg(kl)->GetTmpTrkFlag()==1) {
	      pTempSeg[jk]->GetAssocSegBeg(kl)->SetTmpTrkFlag(2); //260709 not found any difference with Asmita's code
	    }
	  }
	}
      }
      
      // Look for preferred joins    ij <-> ij+1 ( <-> ij+2 )
      // Loop over segments on plane ij, which we know have possible end associations
      for(unsigned int jk=0; jk<mTempSeg.size(); ++jk) {
	JoinFlag=false;
	
	// Loop over these end associations. If one end associated triplet already
	// has temp track flag 2, set JoinFlag to true.
	for(unsigned int kl=0; kl<mTempSeg[jk]->GetNAssocSegEnd(); ++kl) {
	  if(mTempSeg[jk]->GetAssocSegEnd(kl)->GetTmpTrkFlag()>1) {
	    JoinFlag=true;
	    break;
	  }
	}
        
	// If JoinFlag is true, and an end assoc segment has temp track flag 1 (i.e. possible 
	// beginning and end associations), set the temp track flag to 2 for this end assoc segment
	if(JoinFlag==true) {
	  // Loop over end associations again and set temp track flags
	  // for segments on plane i+1
	  for(unsigned int kl=0; kl<mTempSeg[jk]->GetNAssocSegEnd(); ++kl) {
	    if(mTempSeg[jk]->GetAssocSegEnd(kl)->GetTmpTrkFlag()==1) {
	      mTempSeg[jk]->GetAssocSegEnd(kl)->SetTmpTrkFlag(2);
	    }
	  }
	}
      }
      
      // Now have temporary track flags set to 2 for all segments on planes i and i+1
      // which look like very promising track segments
      
      // MAKE THE PREFERRED JOINS BETWEEN TRIPLETS ON PLANE i AND PLANE i+1
      if(mTempSeg.size()>0 && pTempSeg.size()>0) {
	// Look for doubly preferred joins    ( ij-1 <-> ) ij <-> ij+1 ( <-> ij+2 )
	// Loop over segments on plane ij, which we know have possible end associations
	for(unsigned int jk=0; jk<mTempSeg.size(); ++jk) {
	  InoTrackSegment* Seg = mTempSeg[jk];
	  // For these, loop over possible end associations
	  for(unsigned int kl=0; kl<Seg->GetNAssocSegEnd(); ++kl) {
	    InoTrackSegment* Segp = Seg->GetAssocSegEnd(kl);
	    // If both segments given temp track flag 2 above check assocations
	    // between triplets on either side
	    if(Seg->GetTmpTrkFlag()>1 && Segp->GetTmpTrkFlag()>1) {
	      // For each Segp, see if there is a Segm association
	      mJnFlag=false;
	      for(unsigned int lm=0; lm<Seg->GetNAssocSegBeg(); ++lm) {
		InoTrackSegment* Segm = Seg->GetAssocSegBeg(lm);
		if(Segm->IsAssoc(Segp) ) {
		  mJnFlag=true;
		  break;
		}
	      }
	      
	      // For each end assoc of Segp, see if it is assoc with Seg
	      pJnFlag=false;
	      for(unsigned int lm=0; lm<Segp->GetNAssocSegEnd(); ++lm) {
		InoTrackSegment* Segp2 = Segp->GetAssocSegEnd(lm);
		if(Seg->IsAssoc(Segp2) ) {
		  pJnFlag=true; break;
		}
	      }
	      
	      // Make the preferred assocations
	      if(mJnFlag==true && pJnFlag==true) {
		//GMA need to check this criteria and optimise it // asm:will fail which the triplets zlength is not the same
		if((fabs((Seg->GetBegXPos()-Seg->GetEndXPos())-(Segp->GetBegXPos()-Segp->GetEndXPos()))<GradientXTolerance*StripXWidth) && (fabs((Seg->GetBegYPos()-Seg->GetEndYPos())-(Segp->GetBegYPos()-Segp->GetEndYPos()))<GradientYTolerance*StripYWidth) && Seg->GetBegZPlane()<=Segp->GetBegZPlane() && Seg->GetEndZPlane()<=Segp->GetEndZPlane() ) {
		  Seg->AddPrefSegToEnd(Segp);
		  Segp->AddPrefSegToBeg(Seg);
		} else {
		  if (pAnalysis->isXtermOut==1) {
		    cout << "Segx " << Seg->GetBegZPlane() << " " 
			 << Seg->GetEndZPlane() << " "
			 << Seg->GetBegZPos() << " " 
			 << Seg->GetEndZPos() << " " 
			 << Seg->GetBegXPos() << " " 
			 << Seg->GetBegYPos() << " "
			 << Seg->GetEndXPos() << " " 
			 << Seg->GetEndYPos() <<" Segp " 
			 << Segp->GetBegZPlane() << " " 
			 << Segp->GetEndZPlane() << " "
			 << Segp->GetBegZPos() << " " 
			 << Segp->GetEndZPos() << " " 
			 << Segp->GetBegXPos() << " " 
			 << Segp->GetBegYPos() << " " 
			 << Segp->GetEndXPos() << " "
			 << Segp->GetEndYPos() << endl; 
		  }
		  //pAnalysis->isXtermOut)
		}
	      }
	    }
	  }
	}
	
	// Look for singly preferred joins (1)    ( ij-1 <-> ) ij <-> ij+1 
	
	// Loop over segments on plane ij, which we know have possible end associations
	for(unsigned int jk=0; jk<mTempSeg.size(); ++jk) {
	  InoTrackSegment* Seg = mTempSeg[jk];
	  // If segment given temp track flag 2 above
	  if(Seg->GetTmpTrkFlag()>1) {
	    JoinFlag=false;
	    // Check to see if we will already have made the preferred association above
	    for(unsigned int kl=0; kl<Seg->GetNAssocSegEnd(); ++kl) {
	      InoTrackSegment* Segp = Seg->GetAssocSegEnd(kl);
	      
	      if(Segp->GetTmpTrkFlag()>1) {
		JoinFlag=true;
		break;
	      }
	    }
	    
	    // If have not made the association, can check assocations
	    // between triplets on either side
	    
	    if(JoinFlag==false) {
	      for(unsigned int kl=0; kl<Seg->GetNAssocSegEnd(); ++kl) {
		InoTrackSegment* Segp = Seg->GetAssocSegEnd(kl);
		// For each Segp, see if there is a Segm associated
		mJnFlag=false;
		for(unsigned int lm=0; lm<Seg->GetNAssocSegBeg(); ++lm) {
		  InoTrackSegment* Segm = Seg->GetAssocSegBeg(lm);
		  if(mJnFlag==false && Segm->IsAssoc(Segp) ) {
		    mJnFlag=true;
		    break;
		  }
		}
		
		// Make the preferred assocations
		if(mJnFlag==true) {
		  //GMA need to check this criteria and optimise it
		  if((fabs((Seg->GetBegXPos()-Seg->GetEndXPos())-(Segp->GetBegXPos()-Segp->GetEndXPos()))<GradientXTolerance*StripXWidth) && (fabs((Seg->GetBegYPos()-Seg->GetEndYPos())-(Segp->GetBegYPos()-Segp->GetEndYPos()))<GradientYTolerance*StripYWidth) && Seg->GetBegZPlane()<=Segp->GetBegZPlane() && Seg->GetEndZPlane()<=Segp->GetEndZPlane() ) {
		    Seg->AddPrefSegToEnd(Segp);
		    Segp->AddPrefSegToBeg(Seg);
		  } else {
		    if (pAnalysis->isXtermOut==1) {
		      cout << "Segy " 
			   << Seg->GetBegZPlane() << " " 
			   << Seg->GetEndZPlane() << " "
			   << Seg->GetBegZPos() << " " 
			   << Seg->GetEndZPos() << " " 
			   << Seg->GetBegXPos() << " " 
			   << Seg->GetBegYPos() << " " 
			   << Seg->GetEndXPos() << " " 
			   << Seg->GetEndYPos() << " Segp " 
			   << Segp->GetBegZPlane() << " " 
			   << Segp->GetEndZPlane() << " "
			   << Segp->GetBegZPos() << " " 
			   << Segp->GetEndZPos() << " " 
			   << Segp->GetBegXPos() << " " 
			   << Segp->GetBegYPos() << " " 
			   << Segp->GetEndXPos() << " " 
			   << Segp->GetEndYPos() << endl;
		    }
		    //if isXterOut
		  }
		}
	      }
	    }
	  }
	}
	
	// Look for singly preferred joins (2)    ij <-> ij+1 ( <-> ij+2 )
	// Loop the over segments on plane ij+1, which we know have possible beginning associations
	for(unsigned int jk=0; jk<pTempSeg.size(); ++jk) {
	  InoTrackSegment* Segp = pTempSeg[jk];
	  // If segment given temp track flag 2 above
	  if(Segp->GetTmpTrkFlag()>1) {
	    JoinFlag=false;
	    // Check to see if we will already have made the preferred association above
	    for(unsigned int kl=0; kl<Segp->GetNAssocSegBeg(); ++kl) {
	      InoTrackSegment* Seg = Segp->GetAssocSegBeg(kl);
	      
	      if(Seg->GetTmpTrkFlag()>1) {
		JoinFlag=true;
		break;
	      }
	    }
	    
	    // If have not made the association, can check assocations
	    // between triplets on either side
	    
	    if(JoinFlag==false) {
	      for(unsigned int kl=0; kl<Segp->GetNAssocSegBeg(); ++kl) {
		InoTrackSegment* Seg = Segp->GetAssocSegBeg(kl);
		// For each Seg, see if there is a Segp2 associated
		pJnFlag=false;
		for(unsigned int lm=0; lm<Segp->GetNAssocSegEnd(); ++lm) {
		  InoTrackSegment* Segp2 = Segp->GetAssocSegEnd(lm);
		  if(pJnFlag==false && Seg->IsAssoc(Segp2) ) {
		    pJnFlag=true;
		    break;
		  }
		}
		
		// Make the preferred assocations
		if(pJnFlag==true) {
		  if((fabs((Seg->GetBegXPos()-Seg->GetEndXPos())-(Segp->GetBegXPos()-Segp->GetEndXPos()))<GradientXTolerance*StripXWidth) && (fabs((Seg->GetBegYPos()-Seg->GetEndYPos())-(Segp->GetBegYPos()-Segp->GetEndYPos()))<GradientYTolerance*StripYWidth) && Seg->GetBegZPlane()<=Segp->GetBegZPlane() && Seg->GetEndZPlane()<=Segp->GetEndZPlane() ) {
		    Seg->AddPrefSegToEnd(Segp);
		    Segp->AddPrefSegToBeg(Seg);
		  } else {
		    if (pAnalysis->isXtermOut==1) {
		      cout << "Segz " 
			   << Seg->GetBegZPlane() << " " 
			   << Seg->GetEndZPlane() << " "
			   << Seg->GetBegZPos() << " " 
			   << Seg->GetEndZPos() << " " 
			   << Seg->GetBegXPos() << " " 
			   << Seg->GetEndXPos() << " " 
			   << Seg->GetBegYPos() << " " 
			   << Seg->GetEndYPos() << " Segp " 
			   << Segp->GetBegZPlane() << " " 
			   << Segp->GetEndZPlane() << " "
			   << Segp->GetBegZPos() << " " 
			   << Segp->GetEndZPos() << " " 
			   << Segp->GetBegXPos() << " " 
			   << Segp->GetEndXPos() << " " 
			   << Segp->GetBegYPos() << " " 
			   << Segp->GetEndYPos() << endl;
		    }
		    // if (pAnalysis->isXtermOut){
		  }
		}
	      }
	    }
	  }
	}
	
	// Look for other joins we have missed    ij <-> ij+1
	// Loop over segments on plane ij, which we know have possible end associations
	for(unsigned int jk=0; jk<mTempSeg.size(); ++jk) {
	  InoTrackSegment* Seg = mTempSeg[jk];
	  // For these, loop over the end assocations
	  for(unsigned int kl=0; kl<Seg->GetNAssocSegEnd(); ++kl) {
	    InoTrackSegment* Segp = Seg->GetAssocSegEnd(kl);
	    // If segments and its end assoc segment were given temp track flag 2 above
	    if(Seg->GetTmpTrkFlag()<2 && Segp->GetTmpTrkFlag()<2) {
	      // Look at each possible end assoc for Seg to see if we've already
	      // made the assocation
	      mJnFlag=false;
	      for(unsigned int lm=0; lm<Seg->GetNAssocSegEnd(); ++lm) {
		if(Seg->GetAssocSegEnd(lm)->GetTmpTrkFlag()>1) {
		  mJnFlag=true;
		  break;
		}
	      }
	      
	      // Look at each possible beginning assoc for Segp to see if we've already
	      // made the assocation
	      
	      pJnFlag=false;
	      for(unsigned int lm=0; lm<Segp->GetNAssocSegBeg(); ++lm) {
		if(Segp->GetAssocSegBeg(lm)->GetTmpTrkFlag()>1) {
		  pJnFlag=true;
		  break;
		}
	      }
	      
	      // Make the preferred associations if it hasn't been done by one of the loops above
	      if(mJnFlag==false && pJnFlag==false) {
		if((fabs((Seg->GetBegXPos()-Seg->GetEndXPos())-(Segp->GetBegXPos()-Segp->GetEndXPos()))<GradientXTolerance*StripXWidth) && (fabs((Seg->GetBegYPos()-Seg->GetEndYPos())-(Segp->GetBegYPos()-Segp->GetEndYPos()))<GradientYTolerance*StripYWidth) && Seg->GetBegZPlane()<=Segp->GetBegZPlane() && Seg->GetEndZPlane()<=Segp->GetEndZPlane() ) {
		  Seg->AddPrefSegToEnd(Segp);
		  Segp->AddPrefSegToBeg(Seg);
		} else {
		  if (pAnalysis->isXtermOut==1) {
		    cout << "Seg " 
			 << Seg->GetBegZPlane() << " " 
			 << Seg->GetEndZPlane() << " "
			 << Seg->GetBegZPos() << " " 
			 << Seg->GetEndZPos() << " " 
			 << Seg->GetBegXPos() << " " 
			 << Seg->GetEndXPos() << " " 
			 << Seg->GetBegYPos() << " " 
			 << Seg->GetEndYPos() << " Segp " 
			 << Segp->GetBegZPlane() << " " 
			 << Segp->GetEndZPlane() << " "
			 << Segp->GetBegZPos() << " " 
			 << Segp->GetEndZPos() << " " 
			 << Segp->GetBegXPos() << " " 
			 << Segp->GetEndXPos() << " " 
			 << Segp->GetBegYPos() << " "
			 << Segp->GetEndYPos() << endl;
		  }
		  // if (pAnalysis->isXtermOut){
		}
	      }
	    }
	  }
	}
      }
      
      // Clear flags ready to consider next planes
      for(unsigned int jk=0; jk<mTempSeg.size(); ++jk) {
	mTempSeg[jk]->SetTmpTrkFlag(0);
      }
      
      for(unsigned int jk=0; jk<pTempSeg.size(); ++jk) {
	pTempSeg[jk]->SetTmpTrkFlag(0);
      }
    }
    // End loop over planes in module
  }
  // End loop over modules
  // Look for any other preferred associations we may have missed.
  // Find segments with no preferred end association and look a bit
  // further to find preferred associations.
  
  //--------------------------------------------------------------------------
  int NewPlane; int Increment;
  double SegXPos, NearbySegXPos, SegYPos, NearbySegYPos;
  bool AssocsStopHere;
  bool ConsiderSegment;
  
  for(int Module=0; Module<NumModules; ++Module) {
    for(int Plane=1; Plane<PlanesInModule-1; ++Plane) {
      int ij=Plane + Module*(PlanesInModule+1); //GMA14, "i", should be local variable
      
      // Loop over all segments
      for(unsigned int jk=0; jk<SegmentBank[ij].size(); ++jk) {
	InoTrackSegment* Seg0 = SegmentBank[ij][jk];
	
	// Loop over the segments with no preferred end association
	if(Seg0->GetNPrefSegEnd()==0 && Seg0->GetNPrefSegBeg()>0) {
	  // Look at nearby segments to see if these already continue the associations.
	  AssocsStopHere=true;
	  SegXPos=0.5*(Seg0->GetBegXPos()+Seg0->GetEndXPos());
	  SegYPos=0.5*(Seg0->GetBegYPos()+Seg0->GetEndYPos());
	  for(unsigned int kl=0; kl<SegmentBank[ij].size(); ++kl) {
	    InoTrackSegment* NearbySeg = SegmentBank[ij][kl];
	    
	    if(NearbySeg==Seg0) {continue;}
	    NearbySegXPos=0.5*(NearbySeg->GetBegXPos()+NearbySeg->GetEndXPos());
	    NearbySegYPos=0.5*(NearbySeg->GetBegYPos()+NearbySeg->GetEndYPos());
	    
	    if(fabs(NearbySegXPos-SegXPos)<1 && fabs(NearbySegYPos-SegYPos)<1 && NearbySeg->GetNPrefSegEnd()>0 && NearbySeg->GetNPrefSegBeg()>0) {
	      AssocsStopHere=false;
	      break;
	    }
	  }
	  if(AssocsStopHere==false) {
	    break;
	  }
	  
	  JoinFlag=false;
	  
	  // Consider the segment, plus all the segments, SegBeg, with which
	  // it has a preferred beginning association.
	  
	  for(int kl=-1; kl<int(Seg0->GetNPrefSegBeg()); ++kl) {
	    // First time, consider Seg0, rather than one of its preferred beginning associations.
	    InoTrackSegment* SegBeg = 0;
	    if(kl<0) {
	      SegBeg=Seg0;
	    } else {
	      SegBeg = Seg0->GetPrefSegBeg(kl);
	    }
	    
	    Increment = 2; // Original Increment=4;
	    // Look up to 6 planes ahead within the same module.
	    while(Increment<=7)	{
	      //GMA Optimise it
	      NewPlane=SegBeg->GetEndZPlane()+Increment;
	      if(NewPlane>=(Module+1)*(PlanesInModule+1)) {
		break;
	      }
	      
	      // See if SegBeg is associated with a segment on this plane, NextSeg.
	      // If so, also see if SegBeg is associated with one of NextSeg's
	      // preferred end associated segments.
	      
	      // If so, make the preferred associations between SegBeg and NextSeg.
	      for(unsigned int lm=0; lm<SegmentBank[NewPlane].size(); ++lm) {
		InoTrackSegment* NextSeg = SegmentBank[NewPlane][lm];
		
		// Check segments are not now linked in a chain of preferred associations.
		ConsiderSegment=true;
		for(unsigned int am=0; am<NextSeg->GetNPrefSegBeg(); ++am) {
		  InoTrackSegment* PrefSeg = NextSeg->GetPrefSegBeg(am);
		  if(PrefSeg->GetTmpTrkFlag()==1) {
		    NextSeg->SetTmpTrkFlag(1);
		    ConsiderSegment=false;
		    break;
		  }
		}
		if(ConsiderSegment==false) {
		  continue;
		}
		if( SegBeg->IsAssoc(NextSeg) ) {
		  for(unsigned int am=0; am<NextSeg->GetNPrefSegEnd(); ++am) {
		    InoTrackSegment* SegEnd = NextSeg->GetPrefSegEnd(am);
		    if( SegBeg->IsAssoc(SegEnd) ) {
		      JoinFlag=true;
		      SegBeg->AddPrefSegToEnd(NextSeg);
		      NextSeg->AddPrefSegToBeg(SegBeg);
		      SegBeg->SetTmpTrkFlag(1); NextSeg->SetTmpTrkFlag(1);
		      
		      // cout << "Made missing Pref assoc, SegBeg "
		      // << SegBeg->GetBegZPlane() << ", " << SegBeg->GetEndZPlane() << ", "
		      //			   << SegBeg->GetBegXPos() << ", " << SegBeg->GetEndXPos()<< ", " 
		      //			   << SegBeg->GetBegYPos() << ", " << SegBeg->GetEndYPos()
		      //			   << " NextSeg " << NextSeg->GetBegZPlane() << ", " 
		      //			   << NextSeg->GetEndZPlane() << ", " 
		      //			   << NextSeg->GetBegXPos() << ", " 
		      //			   << NextSeg->GetEndXPos() << ", " 
		      //			   << NextSeg->GetBegYPos() << ", " 
		      //			   << NextSeg->GetEndYPos() << endl;
		      break;
		    }
		  }
		}
	      }
	      Increment++;
	    }
	    // Reset TmpTrkFlags
	    for(int lm=ij; lm<(1+Module)*PlanesInModule; ++lm) {
	      for(unsigned int am=0; am<SegmentBank[lm].size(); ++am) {
		SegmentBank[lm][am]->SetTmpTrkFlag(0);
	      }
	    }
	    if(JoinFlag==true) {
	      break;
	    }
	  }
	}
	
	// Now loop over the segments with no preferred beginning association
	if(Seg0->GetNPrefSegBeg()==0 && Seg0->GetNPrefSegEnd()>0) {
	  // Look at nearby segments to see if these already continue the associations.
	  AssocsStopHere=true;
	  SegXPos=0.5*(Seg0->GetBegXPos()+Seg0->GetEndXPos());
	  SegYPos=0.5*(Seg0->GetBegYPos()+Seg0->GetEndYPos());
	  
	  for(unsigned int kl=0; kl<SegmentBank[ij].size(); ++kl) {
	    InoTrackSegment* NearbySeg = SegmentBank[ij][kl];
	    if(NearbySeg==Seg0) {continue;}
	    
	    NearbySegXPos=0.5*(NearbySeg->GetBegXPos()+NearbySeg->GetEndXPos());
	    NearbySegYPos=0.5*(NearbySeg->GetBegYPos()+NearbySeg->GetEndYPos());
	    
	    if(fabs(NearbySegXPos-SegXPos)<1 && fabs(NearbySegYPos-SegYPos)<1 && NearbySeg->GetNPrefSegEnd()>0 && NearbySeg->GetNPrefSegBeg()>0) {
	      AssocsStopHere=false;
	      break;
	    }
	  }
	  
	  if(AssocsStopHere==false) {
	    break;
	  }
	  
	  JoinFlag=false;
	  
	  // Consider the segment, plus all the segments, SegEnd, with which
	  // it has a preferred end association.
	  
	  for(int kl=-1; kl<int(Seg0->GetNPrefSegEnd()); ++kl) {
	    // First time, consider Seg0, rather than one of its preferred end associations.
	    InoTrackSegment* SegEnd = 0;
	    if(kl<0) {
	      SegEnd=Seg0;
	    } else {
	      SegEnd = Seg0->GetPrefSegEnd(kl);
	    }
	    Increment=2;
	    
	    // Look up to 6 planes behind within the same module
	    
	    while(Increment<=7) {
	      NewPlane=SegEnd->GetBegZPlane()-Increment;
	      if(NewPlane<=(Module)*(PlanesInModule+1)) {break;}
	      
	      // See if SegEnd is associated with a segment on this plane, PrevSeg.
	      // If so, also see if SegEnd is associated with one of PrevSeg's 
	      // preferred beginning associated segments.
	      // If so, make the preferred associations between SegEnd and PrevSeg.
	      
	      for(unsigned int lm=0; lm<SegmentBank[NewPlane].size(); ++lm) {
		InoTrackSegment* PrevSeg = SegmentBank[NewPlane][lm];
		
		// Check segments are not now linked in a chain of preferred associations.
		ConsiderSegment=true;
		for(unsigned int am=0; am<PrevSeg->GetNPrefSegEnd(); ++am) {
		  InoTrackSegment* PrefSeg = PrevSeg->GetPrefSegEnd(am);
		  if(PrefSeg->GetTmpTrkFlag()==1) {
		    PrevSeg->SetTmpTrkFlag(1);
		    ConsiderSegment=false;
		    break;
		  }
		}
		if(ConsiderSegment==false) {
		  continue;
		}
		
		if( PrevSeg->IsAssoc(SegEnd) ) {
		  for(unsigned int am=0; am<PrevSeg->GetNPrefSegBeg(); ++am) {
		    InoTrackSegment* SegBeg = PrevSeg->GetPrefSegBeg(am);
		    if( SegBeg->IsAssoc(SegEnd) ) {
		      JoinFlag=true;
		      SegEnd->AddPrefSegToBeg(PrevSeg);
		      PrevSeg->AddPrefSegToEnd(SegEnd);
		      
		      SegEnd->SetTmpTrkFlag(1); PrevSeg->SetTmpTrkFlag(1);
		      if (pAnalysis->isXtermOut==1) {
			cout << "Made missing Pref assoc, SegEnd "
			     << SegEnd->GetBegZPlane() << ", " << SegEnd->GetEndZPlane() << ", "
			     << SegEnd->GetBegXPos() << ", " << SegEnd->GetEndXPos()<<", "
			     << SegEnd->GetBegYPos() << ", " << SegEnd->GetEndYPos()
			     << " PrevSeg " << PrevSeg->GetBegZPlane() << ", " 
			     << PrevSeg->GetEndZPlane() << ", " 
			     << PrevSeg->GetBegXPos() <<", "
			     << PrevSeg->GetEndXPos() <<", "
			     << PrevSeg->GetBegYPos() <<", "
			     << PrevSeg->GetEndYPos() << endl;
		      }
		      // if (pAnalysis->isXtermOut){
		      break;
		    }
		  }
		}
	      }
	      Increment++;
	    }
	    // Reset TmpTrkFlags
	    for(int lm=ij; lm>Module*(PlanesInModule); --lm) {
	      for(unsigned int am=0; am<SegmentBank[lm].size(); ++am) {
		SegmentBank[lm][am]->SetTmpTrkFlag(0);
	      }
	    }
	    if(JoinFlag==true)
	      {
		break;
	      }
	  }
	}
      }
    }
  }
  
  if(pAnalysis->isXtermOut==1) {
    // Print out list of preferred associations
    cout << " InoTrackFinder : *** LIST OF PREF SEG ASSOCIATIONS *** " << endl;
    // Print out list of triplets
    for(int ij=0; ij<500; ++ij) {
      for(unsigned int jk=0;jk<SegmentBank[ij].size(); ++jk) {
	InoTrackSegment* segment = SegmentBank[ij][jk];
	cout << " ----Plane " << ij << ", Seg Number " << jk
	     << ", BegPlane " << segment->GetBegZPlane()
	     << ", EndPlane " << segment->GetEndZPlane()
	     << ", BegXPos " << segment->GetBegXPos()
	     << ", EndXPos " << segment->GetEndXPos()
	     << ", BegYPos " << segment->GetBegYPos()
	     << ", EndYPos " << segment->GetEndYPos()
	     <<endl;
	
	cout << " Beg: " << endl;
	for(unsigned int kl=0; kl<segment->GetNPrefSegBeg(); ++kl) {
	  InoTrackSegment* segtmp = segment->GetPrefSegBeg(kl);
	  cout << "  Pref number=" << kl
	       << "  begpln: " << segtmp->GetBegZPlane()
	       << "  begxpos: " << segtmp->GetBegXPos()
	       << "  begypos: " << segtmp->GetBegYPos()
	       << "  endpln: " << segtmp->GetEndZPlane()
	       << "  endxpos: " << segtmp->GetEndXPos()
	       << "  endypos: " << segtmp->GetEndYPos()
	       << endl;
	}
	
	cout << " End: " << endl;
	for(unsigned int kl=0; kl<segment->GetNPrefSegEnd(); ++kl) {
	  InoTrackSegment* segtmp = segment->GetPrefSegEnd(kl);
	  cout << "  Pref number=" << kl
	       << "  begpln: " << segtmp->GetBegZPlane()
	       << "  begxpos: " << segtmp->GetBegXPos()
	       << "  begypos: " << segtmp->GetBegYPos()
	       << "  endpln: " << segtmp->GetEndZPlane()
	       << "  endxpos: " << segtmp->GetEndXPos()
	       << "  endypos: " << segtmp->GetEndYPos()
	       << endl;
	}
      }
    }
  }
  return;
}

void InoTrackFinder::FindMatchedJoins() {
  // Having made all the preferred assocations, we actually match triplets
  // together to form longer segments of the track.
  // We use the idea of seed segments to identify the last segment in a chain
  // of segments.
  //MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
  
  //  int i;
  bool JoinFlag;
  int factX, factY;
  
  // Loop over the planes
  for(int Module = 0; Module < NumModules; ++Module) {
    for(int Plane = 1; Plane < PlanesInModule - 1; ++Plane) {
      int ij = Plane + Module * (PlanesInModule + 1); //GMA14 "i", need to bne defined as local
      
      // Loop over the triplets centered on plane i
      for(unsigned int jk = 0; jk < SegmentBank[ij].size(); ++jk) {
	//---------------------
	// Firstly, consider isolated segments. These must be their own seed segments.
	if(SegmentBank[ij][jk]->GetNPrefSegBeg()==0 && SegmentBank[ij][jk]->GetNPrefSegEnd()==0) {
	  ViewSegBank[0].push_back(SegmentBank[ij][jk]);
	  TempTrack[0].push_back(SegmentBank[ij][jk]);
	  SegmentBank[ij][jk]->SetSeedSegment(SegmentBank[ij][jk]);
	}
	
	//---------------------
	// Then, consider those segments which have some preferred associations
	if(SegmentBank[ij][jk]->GetNPrefSegBeg()>0 || SegmentBank[ij][jk]->GetNPrefSegEnd()>0) {
	  //asm: ijs the second condition required
	  JoinFlag=false;
	  
	  // If there is one pref assoc at beginning (always false for first segment in module)
	  if(SegmentBank[ij][jk]->GetNPrefSegBeg()==1) {
	    InoTrackSegment* Segm = SegmentBank[ij][jk]->GetPrefSegBeg(0);
	    
	    // Check that we are considering a simple chain of segments
	    if(Segm->GetNPrefSegEnd()==1) {
	      // Get the last segment in the chain of segments
	      InoTrackSegment* SeedSeg = Segm->GetSeedSegment();
	      SeedSeg->AddSegment(SegmentBank[ij][jk]);
	      SegmentBank[ij][jk]->SetSeedSegment(SeedSeg);
	      SegmentBank[ij][jk]->SetUID(-1);        		//fUID of the segment is set to -1 when it is 
	      
	      // cout << "trkfindmatchedjoin "<< ij<<" "<<SegmentBank[ij][jk]->GetUID()<<endl;
	      JoinFlag=true;
	    }
	  }
	  //----------------------------------
	  //Asm: not required
	  
	  if( SegmentBank[ij][jk]->GetNPrefSegBeg()==0 && SegmentBank[ij][jk]->GetNPrefSegEnd()==1) {
	    InoTrackSegment* Segp = SegmentBank[ij][jk]->GetPrefSegEnd(0);
	    
	    // Check that we are considering a simple chain of segments
	    if(Segp->GetNPrefSegBeg()==1&& Segp->GetNPrefSegEnd()==0) {
	      //(asm: 080311) 
	      // Get the last segment in the chain of segments
	      
	      SegmentBank[ij][jk]->SetSeedSegment(SegmentBank[ij][jk]);
	      InoTrackSegment* SeedSeg1 = SegmentBank[ij][jk]->GetSeedSegment();
	      
	      SeedSeg1->AddSegment(Segp);
	      Segp->SetSeedSegment(SeedSeg1);
	      Segp->SetUID(-1);                       	//fUID of the segment is set to -1 when it is
	      
	      JoinFlag=true;
	    }
	  }
	  //
	  //--------------------
	  //--------------------
	  //cout <<"trkfindmatchedjoin "<< ij<<" "<<SegmentBank[ij].size()<<endl;
	  // Now consider those triplets which aren't just one pref assoc at beg and one at end
	  
	  if(JoinFlag==false) {
	    ViewSegBank[0].push_back(SegmentBank[ij][jk]);
	    TempTrack[0].push_back(SegmentBank[ij][jk]);
	    SegmentBank[ij][jk]->SetSeedSegment(SegmentBank[ij][jk]);
	    
	    // Loop over the preferred associations at beginning
	    for(unsigned int kl=0; kl<SegmentBank[ij][jk]->GetNPrefSegBeg(); ++kl) {
	      InoTrackSegment* SeedSeg = SegmentBank[ij][jk]->GetPrefSegBeg(kl)->GetSeedSegment();
	      InoTrackSegment* MatchSeg = SegmentBank[ij][jk];

	      //asm:.....................condition in this part modified
	      if(SeedSeg->GetBegZPlane()<=MatchSeg->GetBegZPlane() && SeedSeg->GetEndZPlane()<=MatchSeg->GetEndZPlane()) {
		factX=((SeedSeg->GetBegXPos()-SeedSeg->GetEndXPos())/(SeedSeg->GetEndZPlane()-SeedSeg->GetBegZPlane()))<10*StripXWidth?4:10;
		factY=((SeedSeg->GetBegYPos()-SeedSeg->GetEndYPos())/(SeedSeg->GetEndZPlane()-SeedSeg->GetBegZPlane()))<10*StripYWidth?4:10;
		
		if( (SeedSeg->GetEntries()>3 && MatchSeg->GetEntries()>3) || (fabs( (SeedSeg->GetBegXPos()-SeedSeg->GetEndXPos())/(SeedSeg->GetEndZPlane()-SeedSeg->GetBegZPlane()) - (MatchSeg->GetBegXPos()-MatchSeg->GetEndXPos())/(MatchSeg->GetEndZPlane()-MatchSeg->GetBegZPlane()) ) < factX*StripXWidth && /*factor 10 should be replaced by no. this is function of the slope  <20 strip -4 * and for >20 strips - 10 */ fabs( (SeedSeg->GetBegYPos()-SeedSeg->GetEndYPos())/(SeedSeg->GetEndZPlane()-SeedSeg->GetBegZPlane()) - (MatchSeg->GetBegYPos()-MatchSeg->GetEndYPos())/(MatchSeg->GetEndZPlane()-MatchSeg->GetBegZPlane()) ) < factY*StripYWidth)) {
		  
		  //asm: ....................................Over.................................................................
		  
		  MatchSeg->AddMatchSegToBeg(SeedSeg);
		  SeedSeg->AddMatchSegToEnd(MatchSeg);
		}
	      }
	    }
	  }
	}
      }// End loop over triplets centered on plane ij
    }// End loop over planes in module
  }// End loop over modules
  
  //  cout << "---------------------------------------------------detgap------------------------------------- "<<endl;
  
  int detgap=2;                                    //this is just to set this part of the code  on and off while testing
  if (detgap==1)  {
    // GMA open this for coil and boundary etc
    // Form matched associations between non-adjacent segments
    vector<InoTrackSegment*> TempSeg1;
    vector<InoTrackSegment*> TempSeg2;
    
    const int npts0=4; int npts; double inv_npts=0; int Module; int Plane;
    bool CheckCoilHole=false;
    
    // Loop over the planeviews U and V
    for(int view=0; view<1; ++view) {
      TempSeg1.clear(); TempSeg2.clear();
      // Loop over entries stored in ViewSegBank, above
      cout<<" ------------------------------------ # of segments "<<ViewSegBank[view].size()<<endl;
      if(ViewSegBank[view].size()<2) continue;
      for(unsigned int ij=0; ij<ViewSegBank[view].size(); ++ij) {
	InoTrackSegment* Seg1 = ViewSegBank[view][ij];
	// If no matched entries at end of segment
	if(Seg1->GetNMatchSegEnd()==0) {
	  JoinFlag=false;
	  // For FD, attempt to match segments across the coil hole
	  // For ND, need to be careful making associations, so we don't join two separate tracks together
	  // Check which module we are in
	  Module=int(Seg1->GetEndZPlane()/(PlanesInModule+1));
	  npts=npts0;
	  //---------------------
	  // These are the tpos values required to locate the Magent coil hole
	  //GMA  put proper value
	  //asm: need to see if taking absolute value take care of the all possible case.
	  
	  double xx = abs(Seg1->GetEndXDir());
	  double yy = abs(Seg1->GetEndYDir());
	  bool rpcGap=false, MegModuleGap=false;
	  float  XXfrac =0; float YYfrac=0;
	  int    XXq=-11, YYq=-11;
	  //ASM //Using this condition first decide if the track is passing through module gap or detector gap.
	  //Set some flag and set the appropriate condition for both the falgs.
	  if((Seg1->GetEndXPos()>(-8.25-min(xx,1.52)*LayerThickness) && /*check this xx*/ Seg1->GetEndXPos()<(-7.95+min(xx,1.52)*LayerThickness)) || (Seg1->GetEndXPos()>( 7.95-min(xx,1.52)*LayerThickness) && Seg1->GetEndXPos()<( 8.25+min(xx,1.52)*LayerThickness))) {
	    MegModuleGap=true;
	    //cout<<"MegModuleGap " << Seg1->GetEndXPos()<<" "<< min(xx,1.52)*LayerThickness<<endl;
	  } else {
	    if(fabs(Seg1->GetEndXPos())<8.00)XXfrac= (int(fabs(Seg1->GetEndXPos()*1000))%2000)/1000.;
	    if(fabs(Seg1->GetEndXPos())>8.20)XXfrac= (int(fabs(Seg1->GetEndXPos()*1000)-200)%2000)/1000.;
	    //XXfrac= (int(fabs(Seg1->GetEndYPos()*1000))%2000)/1000.;
	    YYfrac= (int(fabs(Seg1->GetEndYPos()*1000))%2000)/1000.;
	    // if((fabs(XXfrac-1)>0.95-min(xx,1.52)*LayerThickness||fabs(YYfrac-1)>0.95-min(yy,1.52)*LayerThickness)){
	    if((fabs(XXfrac-1)>0.80||fabs(YYfrac-1)>0.80)) {
	      rpcGap=true;
	      //cout<<"rpcGap"<<endl;
	      if(fabs(XXfrac-1)>0.80) {
		if(fabs(Seg1->GetEndXPos())<8.00) XXq=int((fabs(Seg1->GetBegYPos()*1000))/2000);
		if(fabs(Seg1->GetEndXPos())>8.20) XXq=int((fabs(Seg1->GetEndXPos()*1000)-200)/2000);
	      }
	      if(fabs(YYfrac-1)>0.80)
		YYq=int((fabs(Seg1->GetEndYPos()*1000))/2000);
	    }
	  }
	  //---------------------
	  if(ModuleType==1 && ((MegModuleGap||rpcGap) || 1+(Seg1->GetEndZPlane()-Seg1->GetBegZPlane())/2>5) ) {
	    //asm: *2nd place  ||-> &&
	    if(MegModuleGap){
	      cout<<"pAnalysis->DGap->Fill(4);"<<endl;
	      pAnalysis->DGap->Fill(4);}
	    if(rpcGap) {
	      cout<<"pAnalysis->DGap->Fill(1);"<<endl;
	      pAnalysis->DGap->Fill(1);}
	    cout<<"pAnalysis->DGap->Fill(0);"<<endl;
	    pAnalysis->DGap->Fill(0);
	    //-----------------------------------------------------------------------------
            // Calculate range of possible association and configure
            // - worst case scenario - track crosses entire width of coil hole, tcoil = ~0.6m
            // - it will take a distance in z, zcoil = tcoil * dz/dt to do this
            // - zcoil is equivalent to nplanes = zcoil/0.096 = tcoil / ( 0.096 * dt/dz )
            // nplanes = zcoil/0.096
            // - npts = nplanes/2.0 = 1.0 / ( 0.198 * dt/dz)
            //.. nplanes/2.0 = 0.2 / ( 0.0852 * dt/dz * 2)
            // - where dt/dz = Seg1->GetEndD(X/Y)ir();
            //GMA put proper value

            //    double* xx = Seg1->GetEnd(X/Y)Dir();
            //    double xx = Seg1->GetEndXDir();
            //    double yy = Seg1->GetEndYDir();
	    
            //          inv_npts=0.198*(pow((xx*xx+yy*yy),0.5)); //GMA was inv_npts=0.5*Seg1->GetEnd(X/Y)Dir();
	    
	    //------------------------------------------------------------------------------
	    if(MegModuleGap) { 
	      inv_npts= 0.185*max(1/tan(xx),1/tan(1.55));          //0.852*xx;
	    } else if(rpcGap && XXq>=0) {
	      inv_npts= 0.3*max(1/tan(xx),1/tan(1.55));
	    } else if(rpcGap && YYq>=0) {
	      inv_npts= 0.3*max(1/tan(yy),1/tan(1.55));
	    }
	    
	    if(inv_npts<0) {inv_npts=-inv_npts;}
	    if(inv_npts<0.01) {inv_npts=0.01;}
	    npts=int(1./inv_npts);
	    
	    if(MegModuleGap&& npts<6) {
	      npts=6;
	    } else if (rpcGap&& npts<4) {   //replaced 8 with 13
	      npts=4;
	    }
	    //temp  if(npts>100) {npts=100;}  //replaced 18 with 30 - possibly this should be a function of the gradient of the track?
	  }
	  
	  // Loop over this range of possible association
	  //int inpts= npts-5<=Seg1->GetEndZPlane()?1:npts-5;
	  
	  for(int jk=1; jk<npts+4; ++jk) {
	    CheckCoilHole=false;
	    Plane=Seg1->GetEndZPlane()+jk;     //asm: (2*jk)---->jk;
	    if(JoinFlag==false && Plane>=0 &&Plane<((PlanesInModule+1)*(Module+1)) ) {
	      // Loop over segments on each possible plane
	      for(unsigned int kl=0; kl<SegmentBank[Plane].size(); ++kl) {
		InoTrackSegment* Seg2 = SegmentBank[Plane][kl];
		// Check association, first across coil-hole
		//GMA need proper values
		double xx2 = abs(Seg2->GetBegXDir());
		
		//--------------------------------------------------
		//asm==============condition changed====== even takes angle at which to track entres the gap region into account======
		
		if(MegModuleGap&&((Seg2->GetBegXPos()>(-8.25-min(xx2,1.52)*LayerThickness)&& Seg2->GetBegXPos()<(-7.95+min(xx2,1.52)*LayerThickness)&&Seg1->GetEndXPos()>(-8.25-min(xx2,1.52)*LayerThickness)&& Seg1->GetEndXPos()<(-7.95+min(xx2,1.52)*LayerThickness))||(Seg2->GetBegXPos()>(7.95-min(xx2,1.52)*LayerThickness) && Seg2->GetBegXPos()<(8.25+min(xx2,1.52)*LayerThickness)&&Seg1->GetEndXPos()>(7.95-min(xx2,1.52)*LayerThickness) && Seg1->GetEndXPos()<(8.25+min(xx2,1.52)*LayerThickness)))) {
		  CheckCoilHole=true;//cout<<"detcoilhole";
		  if(MegModuleGap) pAnalysis->DGap->Fill(5);
		} else if(rpcGap) {
		  float  XX1frac=0,  YY1frac=0;
		  int    XX1q=-20,  YY1q=-20;
		  if(fabs(Seg2->GetBegXPos())<8.00) {
		    XX1frac=(int(fabs(Seg2->GetBegXPos()*1000))%2000)/1000.;
		    XX1q=int((fabs(Seg2->GetBegXPos()*1000))/2000);
		  }
		  if(fabs(Seg2->GetBegXPos())>8.20) {
		    XX1frac=(int(fabs(Seg2->GetBegXPos()*1000)-200)%2000)/1000.;
		    XX1q=int((fabs(Seg2->GetBegXPos()*1000)-200)/2000);
		  }
		  YY1frac=(int(fabs(Seg2->GetBegYPos()*1000))%2000)/1000.;
		  YY1q=int((fabs(Seg2->GetBegYPos()*1000))/2000);
		  if(fabs(Seg2->GetBegXPos())<8.00)
		    cout<< " XXq-XX1q "<<fabs(Seg2->GetBegXPos()*1000)<<" "<<XX1q <<" "<<XXq<<" "<<
		      fabs(Seg2->GetBegYPos()*1000)<<" "<<YY1q <<" "<<YYq<<" "<< endl;
		  if((fabs(XX1frac-1)>0.80&&abs(XXq-XX1q)<2)||(fabs(YY1frac-1)>0.80&&abs(YYq-YY1q)<2)) {
		    CheckCoilHole=true;
		    rpcGap=true;
		    //cout<<"RPCGAP";
		    pAnalysis->DGap->Fill(2);
		  }
		}
		//----------------------------------------------------
		// Then main association checks
		//GMA need proper database
		//  cout<< "Is"; 
		if((Seg2->GetSeedSegment()==Seg2 &&  Seg2->GetNMatchSegBeg()==0) && (Seg2->GetBegZPlane()-Seg1->GetEndZPlane()>0 && (Seg1->GetEntries()>4 && Seg2->GetEntries()>4)) && (jk<npts0 || CheckCoilHole==true || (1+(Seg1->GetEndZPlane()-Seg1->GetBegZPlane())/2>5) ) ) {
		  if( Seg1->IsAssoc(Seg2) ) {
		    cout<< "Yes"<<endl;
		    // cout<<"1rpcGap"<<rpcGap<<endl;
		    // ND nearby preliminary track protection
		    //GMA not used here, but may need in future
		    
		    float atanX= atan(Seg1->GetEndXDir())-atan(Seg2->GetBegXDir());
		    float atanY= atan(Seg1->GetEndYDir())-atan(Seg2->GetBegYDir());
		    
		    //pAnalysis->ShwXw->Fill(atanX) ;
		    //pAnalysis->ShwYw->Fill(atanY);
		    if(CheckCoilHole){ 
		      cout<<"pAnalysis->ShwXw->Fill(atanX);\npAnalysis->ShwYw->Fill(atanY);"<<endl;
		      pAnalysis->ShwXw->Fill(atanX);
		      pAnalysis->ShwYw->Fill(atanY);
		      //cout<< "atan"<<atanX<< " "<<atanY<<endl;
		      //cout<<"rrpcGap"<<rpcGap<<endl;
		    }
		    
		    if(MegModuleGap) pAnalysis->DGap->Fill(6);
		    if(rpcGap)       pAnalysis->DGap->Fill(3);
		    
		    // asm" this condition may differ for both the gaps  can be set using the plot ShwXw
		    float atandiff = 0.8;
		    if(MegModuleGap) {
		      atandiff=0.5;
		    } else if(rpcGap) {
		      atandiff = 0.8;
		    }
		    
		    if(fabs(atan(Seg1->GetEndXDir())-atan(Seg2->GetBegXDir()))<atandiff && fabs(atan(Seg1->GetEndYDir())-atan(Seg2->GetBegYDir()))<atandiff && ((Seg1->GetEntries()<11 && Seg2->GetEntries()<11 && Seg2->GetBegZPlane()-Seg1->GetEndZPlane()<50) || Seg2->GetBegZPlane()-Seg1->GetEndZPlane()<100)) {
		      //correct the last two lines
		      TempSeg1.push_back(Seg1);
		      TempSeg2.push_back(Seg2);
		      JoinFlag=true;
		      cout<<" GAP JOINED******************************** " <<endl;
		      if(MegModuleGap) {pAnalysis->DGap->Fill(8);}
		      if(rpcGap) {pAnalysis->DGap->Fill(7);}
		    }//if( ModuleType!=2
		    
		  }  //  if( Seg1->IsAssoc(Seg2) )
		}//
		//  else{cout<<"NoAssoc"<<endl;}
		//-------------------------------------------------------------
	      }
	    }
	  }//jk
	}
      }// End loop over entries in ViewSegBank
      //cout <<"1=======trkfindmatchedjoin ======== "<< endl;
      
      // Make the associations
      for(unsigned int ij=0; ij<TempSeg1.size(); ++ij) {
	InoTrackSegment* Seg1 = TempSeg1[ij];
	InoTrackSegment* Seg2 = TempSeg2[ij];
	
	if(Seg1->GetBegZPlane()<=Seg2->GetBegZPlane() && Seg1->GetEndZPlane()<=Seg2->GetEndZPlane()) {
	  if( (Seg1->GetEntries()>3 && Seg2->GetEntries()>3) || (4.*fabs( (Seg1->GetBegXPos()-Seg1->GetEndXPos())/(Seg1->GetBegZPlane()-Seg1->GetEndZPlane()) - (Seg2->GetBegXPos()-Seg2->GetEndXPos())/(Seg2->GetBegZPlane()-Seg2->GetEndZPlane()) ) <  10*StripXWidth && 4.*fabs( (Seg1->GetBegYPos()-Seg1->GetEndYPos())/(Seg1->GetBegZPlane()-Seg1->GetEndZPlane()) - (Seg2->GetBegYPos()-Seg2->GetEndYPos())/(Seg2->GetBegZPlane()-Seg2->GetEndZPlane()) ) < 10*StripYWidth) ) {
	    Seg2->AddMatchSegToBeg(Seg1);
	    Seg1->AddMatchSegToEnd(Seg2);
	  }
	}
      }
    }// End loop over planeviews
    //  cout <<"2=======trkfindmatchedjoin ======== "<< endl;	
  } //if detgap
  if (detgap==2 ) {
    //----------------------
    // Form matched associations between non-adjacent segments
    vector<InoTrackSegment*> TempSeg1;
    vector<InoTrackSegment*> TempSeg2;
    const int npts0=4; int npts; double inv_npts=0; int Module; int Plane;
    bool CheckCoilHole=false;
    double xxmax = 1.52;
    double fmax =0.92; 
    
    //    int a0=0,a1=0;
    int a2=0,a3=0,b4=0,b5=0,b6=0,a7=0,b8=0;
    for(int view=0; view<1; ++view) {
      TempSeg1.clear(); TempSeg2.clear();
      if(ViewSegBank[view].size()<2) continue;
      
      //// please check this loop
      for(unsigned int ij=0; ij<ViewSegBank[view].size(); ++ij) {
	// Loop over entries stored in ViewSegBank, above
	InoTrackSegment* Seg1 = ViewSegBank[view][ij];
	if(Seg1->GetNMatchSegEnd()==0) {
	  // If no matched entries at end of segment
	  JoinFlag=false;
	  Module=int(Seg1->GetEndZPlane()/(PlanesInModule+1));
	  npts=npts0;
	  double xx = abs(Seg1->GetEndXDir());
	  double yy = abs(Seg1->GetEndYDir());
	  bool rpcGap=false, MegModuleGap=false;
	  //cout<<"  rpcGap=false, MegModuleGap=false " <<rpcGap << " "<<MegModuleGap<<endl;
	  float  XXfrac =0; float YYfrac=0;
	  int    XXq=-110, YYq=-110;
	  //-------------if begin
	  if((Seg1->GetEndXPos()>(-8.25-min(xx,xxmax)*LayerThickness) && /*check this xx*/ Seg1->GetEndXPos()<(-7.95+min(xx,xxmax)*LayerThickness))|| (Seg1->GetEndXPos()>( 7.95-min(xx,xxmax)*LayerThickness) && Seg1->GetEndXPos()<( 8.25+min(xx,xxmax)*LayerThickness) )) {
	    MegModuleGap=true;
	    //cout<<"MegModuleGap " << Seg1->GetEndXPos()<<" "<< min(xx,xxmax)*LayerThickness<<endl;
	  } else {
	    if(fabs(Seg1->GetEndXPos())<8.00) {
	      XXfrac= (int(fabs(Seg1->GetEndXPos()*1000))%2000)/1000.;
	      if(fabs(XXfrac-1)>fmax-fabs((int(fabs(min(xx,xxmax)*LayerThickness*1000))%2000)/1000.-1)) {
		XXq=int((Seg1->GetEndXPos()*1000)/2000)+24;
		rpcGap=true;
	      }
	    }
	    if(fabs(Seg1->GetEndXPos())>8.20) {
	      XXfrac= (int(fabs(Seg1->GetEndXPos()*1000)-200)%2000)/1000.;
	      if(fabs(XXfrac-1)>fmax-fabs((int(fabs(min(xx,xxmax)*LayerThickness*1000)-200)%2000)/1000.-1)) {
		if( Seg1->GetEndXPos()>0)XXq=int(((Seg1->GetEndXPos()*1000)+200)/2000)+24;
		else if( Seg1->GetEndXPos()<0)XXq=int(((Seg1->GetEndXPos()*1000)-200)/2000)+24;
		rpcGap=true;
	      }
	    }
	    
	    YYfrac= (int(fabs(Seg1->GetEndYPos()*1000))%2000)/1000.;
	    if(fabs(YYfrac-1)>fmax-fabs((int(fabs(min(yy,xxmax)*LayerThickness*1000))%2000)/1000.-1)) {
	      YYq=int((Seg1->GetEndYPos()*1000)/2000)+8;
	    }
	    //cout<<"XXq:YYq=" <<XXq<< " "<<YYq<<endl;
	  }
	  //if Ends
	  if(ModuleType==1 && ((MegModuleGap||rpcGap) || 1+(Seg1->GetEndZPlane()-Seg1->GetBegZPlane())/2>5) ) {   //asm: *2nd place  ||-> &&
	    if(MegModuleGap){
	      pAnalysis->DGap->Fill(4); b4++;
	    }
	    if(rpcGap) {
	      //pAnalysis->DGap->Fill(1);a1++;
	    }
	    // pAnalysis->DGap->Fill(0);a0++;
	    
	    if(MegModuleGap) {
	      inv_npts= 0.185*max(1/tan(xx),1/tan(1.55)) ;          //0.852*xx;
	    } else if(rpcGap && XXq>=0) {
	      inv_npts= 0.3*max(1/tan(xx),1/tan(1.55));
	    } else if(rpcGap && YYq>=0) {
	      inv_npts= 0.3*max(1/tan(yy),1/tan(1.55));
	    }
	    
	    if(inv_npts<0) {inv_npts=-inv_npts;}
	    if(inv_npts<0.01) {inv_npts=0.01;}
	    npts=int(1./inv_npts);
	    
	    if(MegModuleGap&& npts<6) {
	      npts=6;   //replaced 8 with 13
	    } else if(rpcGap&& npts<4) {
	      npts=4;    //replaced 8 with 13   //temp  if(npts>100) {npts=100;}  
	    }
	  }
	  
	  // Loop over this range of possible association
	  //int inpts= npts-5<=Seg1->GetEndZPlane()?1:npts-5;
	  for(int jk=1; jk<npts+4; ++jk) {
	    CheckCoilHole=false;
	    Plane=Seg1->GetEndZPlane()+jk;     //asm: (2*jk)---->jk;
	    
	    if( JoinFlag==false && Plane>=0 &&Plane<((PlanesInModule+1)*(Module+1)) ) {
	      for(unsigned int kl=0; kl<SegmentBank[Plane].size(); ++kl) {
		//for(unsigned int kl=0; kl<ViewSegBank[view].size(); ++kl) {
		InoTrackSegment* Seg2 = SegmentBank[Plane][kl];
		//InoTrackSegment* Seg2 = SegmentBank[Plane][kl]->GetSeedSegment();//ViewSegBank[view][k];
		double xx2 = abs(Seg2->GetBegXDir());
		double yy2 = abs(Seg2->GetBegYDir());
		// double xx1 = abs(Seg1->GetEndXDir());
		// double yy1 = abs(Seg2->GetEndYDir());
		//--------------------------------------------------
		//asm==============condition changed====== even takes angle at which to track entres the gap region into account======
		if(MegModuleGap&&((Seg2->GetBegXPos()>(-8.25-min(xx2,xxmax)*LayerThickness)&& Seg2->GetBegXPos()<(-7.95+min(xx2,xxmax)*LayerThickness)&&Seg1->GetEndXPos()>(-8.25-min(xx,xxmax)*LayerThickness)&& Seg1->GetEndXPos()<(-7.95+min(xx,xxmax)*LayerThickness))||(Seg2->GetBegXPos()>(7.95-min(xx2,xxmax)*LayerThickness) && Seg2->GetBegXPos()<(8.25+min(xx2,xxmax)*LayerThickness)&&Seg1->GetEndXPos()>(7.95-min(xx,xxmax)*LayerThickness) && Seg1->GetEndXPos()<(8.25+min(xx,xxmax)*LayerThickness)))) {
		  CheckCoilHole=true;
		  pAnalysis->DGap->Fill(5);b5++;
		} else if(rpcGap) {
		  float  XX1frac=0,  YY1frac=0;
		  int    XX1q=-100,  YY1q=-100;
		  
		  YY1frac=(int(fabs(Seg2->GetBegYPos()*1000))%2000)/1000.;
		  if(fabs(Seg2->GetBegXPos())<8.00) {
		    XX1frac=(int(fabs(Seg2->GetBegXPos()*1000))%2000)/1000.;
		    
		    if((fabs(XX1frac-1)>fmax-fabs((int(fabs(min(xx2,xxmax)*LayerThickness*1000))%2000)/1000.-1))||(fabs(YY1frac-1)>fmax-fabs((int(fabs(min(yy2,xxmax)*LayerThickness*1000))%2000)/1000.-1))) {
		      XX1q=int((Seg2->GetBegXPos()*1000)/2000)+24;
		      YY1q=int((Seg2->GetBegYPos()*1000)/2000)+8;
		    }
		  } else if(fabs(Seg2->GetBegXPos())>8.00) {
		    XX1frac=(int(fabs(Seg2->GetBegXPos()*1000)-200)%2000)/1000.;
		    if((fabs(XX1frac-1)>fmax-fabs((int(fabs(min(xx2,xxmax)*LayerThickness*1000-200))%2000)/1000.-1))||(fabs(YY1frac-1)>fmax-fabs((int(fabs(min(yy2,xxmax)*LayerThickness*1000))%2000)/1000.-1))) {
		      if(Seg2->GetBegXPos()>0) {
			XX1q=int(((Seg2->GetBegXPos()*1000)+200)/2000)+24;
		      } else if(Seg2->GetBegXPos()<0) {
			XX1q=int(((Seg2->GetBegXPos()*1000)-200)/2000)+24;
		      }
		      YY1q=int((Seg2->GetBegYPos()*1000)/2000)+8;
		    }
		  }
		  if((abs(XXq-XX1q)<2)||abs(YYq-YY1q)<2) {
		    CheckCoilHole=true;
		    //cout << "CheckCoilHole=true"<<endl;
		    pAnalysis->DGap->Fill(2);a2++ ;
		  }
		}//else if (rpcgap)
		//cout<< " IsAssoc Test starts"<<endl;
		//----------------------------------------------------
		// Then main association checks  //GMA need proper database
		
		if((Seg2->GetSeedSegment()==Seg2 &&  Seg2->GetNMatchSegBeg()==0) && (Seg2->GetBegZPlane()-Seg1->GetEndZPlane()>0 && (Seg1->GetSeedSegment()->GetEntries()>3 && Seg2->GetEntries()>3)) && (jk<npts0 || CheckCoilHole==true ||(1+(Seg1->GetEndZPlane()-Seg1->GetBegZPlane())/2>5)))  {
		  if( Seg1->IsAssoc(Seg2) ) {
		    //cout<< "Yes"<<endl;
		    float atanX= atan(Seg1->GetEndXDir())-atan(Seg2->GetBegXDir());
		    float atanY= atan(Seg1->GetEndYDir())-atan(Seg2->GetBegYDir());
		    
		    pAnalysis->ShwXw->Fill(atanX);
		    pAnalysis->ShwYw->Fill(atanY);
		    if(CheckCoilHole) {
		      //pAnalysis->ShwXw->Fill(atanX);
		      //pAnalysis->ShwYw->Fill(atanY);
		      //cout<< "atan"<<atanX<< " "<<atanY<<endl;
		    }
		    if(MegModuleGap) {
		      pAnalysis->DGap->Fill(6);b6++;
		    }
		    if(rpcGap) {
		      pAnalysis->DGap->Fill(3);a3++;
		    }
		    // asm" this condition may differ for both the gaps  can be set using the plot ShwXw
		    float atandiff = 0.8;
		    if(MegModuleGap) {
		      atandiff=0.5;
		    } else if(rpcGap) {
		      atandiff = 0.5;
		    }
		    if(fabs(atan(Seg1->GetEndXDir())-atan(Seg2->GetBegXDir()))<atandiff && fabs(atan(Seg1->GetEndYDir())-atan(Seg2->GetBegYDir()))<atandiff && (((Seg1->GetEntries()>11 || Seg2->GetEntries()>11) && Seg2->GetBegZPlane()-Seg1->GetEndZPlane()<60) || Seg2->GetBegZPlane()-Seg1->GetEndZPlane()<30 )) {
		      TempSeg1.push_back(Seg1);
		      TempSeg2.push_back(Seg2);
		      JoinFlag=true;              	
		      //cout<<" GAP JOINED******************************** " <<endl;
		      if(MegModuleGap) {
			pAnalysis->DGap->Fill(8);b8++;
		      }
		      //cout<< " rpcGap "<<rpcGap<<endl;
		      if(rpcGap) {
			pAnalysis->DGap->Fill(7);a7++;
		      }
		    }//if( ModuleType!=2
		  } else { //  if( Seg1->IsAssoc(Seg2) )
		    //cout<< "No "<<endl;
		  }
		  pAnalysis->DGap->Fill(3);a3++;
		}//
		//-------------------------------------------------------------
	      }//k
	    }//if
	  } //jk
	}// End loop over entries in ViewSegBank
      }
      //cout <<"1=======trkfindmatchedjoin ======== "<< endl;
      
      // Make the associations
      for(unsigned int ij=0; ij<TempSeg1.size(); ++ij) {
	InoTrackSegment* Seg1 = TempSeg1[ij];
	InoTrackSegment* Seg2 = TempSeg2[ij];
	if(Seg1->GetBegZPlane()<=Seg2->GetBegZPlane() && Seg1->GetEndZPlane()<=Seg2->GetEndZPlane()) {
	  if( (Seg1->GetEntries()>3 && Seg2->GetEntries()>3) || (4.*fabs( (Seg1->GetBegXPos()-Seg1->GetEndXPos())/(Seg1->GetBegZPlane()-Seg1->GetEndZPlane()) - (Seg2->GetBegXPos()-Seg2->GetEndXPos())/(Seg2->GetBegZPlane()-Seg2->GetEndZPlane()) ) < 10*StripXWidth && 4.*fabs( (Seg1->GetBegYPos()-Seg1->GetEndYPos())/(Seg1->GetBegZPlane()-Seg1->GetEndZPlane()) - (Seg2->GetBegYPos()-Seg2->GetEndYPos())/(Seg2->GetBegZPlane()-Seg2->GetEndZPlane()) ) < 10*StripYWidth) ) {
	    Seg2->AddMatchSegToBeg(Seg1);
	    Seg1->AddMatchSegToEnd(Seg2);
	  }
	}
      }
    }// End loop over planeviews
    //  cout <<"2=======trkfindmatchedjoin ======== "<< endl;
    //cout<< " aaa "<< a0<<" "<<a1<< " "<< a2<<" " << a3<< " "<<b4<<" "<<b5<<" " <<b6<<" " <<a7<<" "<<b8<<endl;;
    // cout << " ----------------------detgap end----------------------------- " <<endl;
  }//if detgap
  
  //========================================================== asm:Over=============================== 
  if(pAnalysis->isXtermOut==1) {
    // Print out list of segments
    //  cout<<"InoTrackFinder : *** LIST OF MATCHED SEG ASSOCIATIONS *** " << endl;
    
    for(int View=0; View<1; ++View) {
      cout<< "InoTrackFinder :View: " << View <<" "<<ViewSegBank[View].size()<< endl;
      for(unsigned int ij=0; ij<ViewSegBank[View].size(); ++ij) {
	InoTrackSegment* segment = (InoTrackSegment*)(ViewSegBank[View][ij]);
	cout<< "InoTrackFinder : ---  Seg0 number =" << ij
	    << "  begpln: "  << segment->GetBegZPlane()
	    << "  begxpos: " << segment->GetBegXPos()
	    << "  begypos: " << segment->GetBegYPos()
	    << "  endpln: "  << segment->GetEndZPlane()
	    << "  endxpos: " << segment->GetEndXPos()
	    << "  endypos: " << segment->GetEndYPos()<< endl;
	
	cout<< "InoTrackFinder :    Beg: " << endl;
	for(unsigned int jk=0; jk<segment->GetNMatchSegBeg(); ++jk) {
	  InoTrackSegment* segtmp = segment->GetMatchSegBeg(jk);
	  cout<<"InoTrackFinder :  Match number=" << jk
	      << "  begpln: "  << segtmp->GetBegZPlane()
	      << "  begxpos: " << segtmp->GetBegXPos()
	      << "  begypos: " << segtmp->GetBegYPos()
	      << "  endpln: "  << segtmp->GetEndZPlane()
	      << "  endxpos: "  << segtmp->GetEndXPos()
	      << "  endypos: "  << segtmp->GetEndYPos()<< endl;
	}
	
	cout << "InoTrackFinder :  End: " << endl;
	for(unsigned int jk=0; jk<segment->GetNMatchSegEnd(); ++jk) {
	  InoTrackSegment* segtmp = segment->GetMatchSegEnd(jk);
	  cout <<"InoTrackFinder :  Match number=" << jk
	       << "  begpln: "  << segtmp->GetBegZPlane()
	       << "  begxpos: " << segtmp->GetBegXPos()
	       << "  begypos: " << segtmp->GetBegYPos()
	       << "  endpln: "  << segtmp->GetEndZPlane()
	       << "  endxpos: " << segtmp->GetEndXPos()
	       << "  endypos: " << segtmp->GetEndYPos()<< endl;
	}
      }
    }
  }// if  isXtermOut
  return; 
}

void InoTrackFinder::FormTracks() {
  // Of the segments identified as good above, we identify the segments which
  // are good seeds for the track i.e. from which we can propagate back and forth
  // along matched segments to find a long track.
 
  // Then, for the first seed segment, we try to propagate backwards and forwards,
  // marking the segments we use with different TmpTrkFlag settings.

  // For a seed segment, Seg0, we firstly move in the direction of increasing plane
  // number. We mark all the segments we use with TmpTrkFlag 1. The seed segment is
  // labelled with TmpTrkFlag 3.

  // Some paths will lead further than others. We are interested in the longest paths, 
  // and move back towards Seg0 along the longest paths, marking the segments used 
  // with TmpTrkFlag 2.

  // Having done this, we can carry out a similar procedure, but initially moving 
  // backwards from Seg0 to lower plane numbers, before returning to Seg0 again.

  // So, this arrangement of segments...
  //
  //       Seg  Seg                                  Seg
  //                  Seg                      Seg
  //                        Seg   Seg0   Seg 
  //                  Seg                      Seg   Seg   Seg
  //                                                 
  //               Seg                               Seg
  //                                                       Seg   Seg

  // Is labelled as follows...
  //
  //       2     2                                    1 
  //                   2                        1 
  //                         2     3      2  
  //                   1                        2     2     2 
  //                                                 
  //                1                                 2 
  //                                                        2     2 
  //

  // It is possible that there are multiple paths along segments with TmpTrkFlag 2, 
  // as above. The best path is selected by obtaining a score for each, based on
  // straightness and length.

  // Define containers
  vector<InoTrackSegment*> BestSeedSegments;
  vector<InoTrackSegment*> SeedSegments;

  vector<InoTrackSegment*> Temp;
  vector<InoTrackSegment*> mTemp;
  vector<InoTrackSegment*> pTemp;
  vector<InoTrackSegment*> BegTemp;
  vector<InoTrackSegment*> EndTemp;
  vector<InoTrackSegment*> TempBank1;
  vector<InoTrackSegment*> TempBank2;

  vector <vector<InoTrackSegment*> > BegBank;
  vector <vector<InoTrackSegment*> > EndBank;

  vector <vector<InoTrackSegment*> > TempBegBank;
  vector <vector<InoTrackSegment*> > TempEndBank;


  bool Cont; bool tmpFlag; bool TrackFlag;
  int ntrks; int Counter; int nplane; int npts;
  int ObjCounter; int ObjVectorCounter; const int fObjVectorMax=100000;
  int id; double Score; double TopScore;

  unsigned int MostClusters;
  bool AlreadyAdded;


  // Main loop is over the views, first for U, then V

  if(ViewSegBank[0].size()>0) { // && ViewSegBank[1].size()>0) {
    for(int View=0; View<1; ++View) {
      
      Cont=true; ntrks=0;
      
      while(Cont==true) {
        Cont=false;
        BegBank.clear(); EndBank.clear();
	
        // First selection of SEED SEGMENTS - segment from which we can move
        // along matched associations to find majority of the track.
        //
        // First guesses at seed segments are stored in SeedSegments Container
        BestSeedSegments.clear(); SeedSegments.clear();
	
        // Loop over the entries in ViewSegBank
        for(unsigned int ij=0; ij<ViewSegBank[View].size(); ++ij) {
	  //	  cout <<"Trkfinderformtrk "<<ViewSegBank[0].size()<<" "<<ij<< endl;
          InoTrackSegment* Seg = ViewSegBank[View][ij];
	  
          if(Seg->GetUID()>0 && Seg->GetUID()<3) {
            Counter=0;
            
            // Loop over the constituent clusters and count the number that are track-like
            for(unsigned int jk=0; jk<Seg->GetEntries(); ++jk) {
              InoCluster* Clust = Seg->GetCluster(jk);
              if(Clust->GetTrkPlnFlag()>0 && Clust->GetTrkFlag()==2) {Counter++;}
            }
            
            // If we've already made a track and this segment doesn't contain 3 track-like
            // clusters, don't consider it any further
            if(ntrks>0 && Counter<3) {Seg->SetUID(0);}
          }

          // Temporarily store all the likely seed segments. Those with UID 2 are most track-like.
          if(Seg->GetUID()==2) {BestSeedSegments.push_back(Seg);}
          if(Seg->GetUID()==1) {SeedSegments.push_back(Seg);}

        }
	
	
	//	cout <<"Trkfinderformtrk 1x  "<<ViewSegBank[0].size()<< endl;
        // Hopefully we have segments with UID 2 (most track-like), else take those with UID 1
        // (just track-like) and move them from SeedSegments to BestSeedSegments. 
        // Can then empty SeedSegments.
        if(BestSeedSegments.size()==0) {BestSeedSegments=SeedSegments;}
        
        SeedSegments.clear();
	
        SeedSegments=BestSeedSegments;
	
        // Now, using our initial seed segments, try propagating back and forth
        // from each segment to refine our choice.
	
        // Final SEED SEGMENTS are stored BestSeedSegments Container
        BestSeedSegments.clear();


        // If we have some seed segments to work with, loop over them
	//	cout <<"Trkfinderformtrk : 2x "<<SeedSegments.size()<< endl;
        if(SeedSegments.size()>0) { 
          nplane=-1;
	  
          for(unsigned int ij=0; ij<SeedSegments.size(); ++ij) {
	    
            // See how far we can propagate backwards to lower plane numbers
            mTemp.clear(); BegTemp.clear(); 
	    
            // First put seed segment into BegTemp then begin while loop
            BegTemp.push_back(SeedSegments[ij]); TempBank1.clear(); 
	    //	    cout <<"Trkfinderformtrk : 3x "<<BegTemp.size()<< endl;

            while(BegTemp.size()>0) {
              Temp.clear();
              
              // Copy any segments into a temporary holder and loop over them,
              // gradually propagating back to lower plane numbers
              Temp=BegTemp; BegTemp.clear();
            
              for(unsigned int jk=0; jk<Temp.size(); ++jk) {
                tmpFlag=false;
                
                // For each segment, loop over their beginning matches
                for(unsigned int kl=0; kl<Temp[jk]->GetNMatchSegBeg(); ++kl) {
                  InoTrackSegment* Segb = Temp[jk]->GetMatchSegBeg(kl);
		  //		  cout <<"Trkfinderformtrk : 3x "<<BegTemp.size()<<" "<< jk<<" "<<kl<<" "<<endl;
                  // If matched segment is track-like and has zero TmpTrkFlag, 
                  // push it back into BegTemp, so the while loop continues.
                  if(Segb->GetUID()>-1 && Segb->GetUID()<3) {
                    tmpFlag=true;
                    
                    if(Segb->GetTmpTrkFlag()<1) {
                      BegTemp.push_back(Segb);
                      Segb->SetTmpTrkFlag(1); // Temporarily mark segments used
                      TempBank1.push_back(Segb);
		      
                    }
                  }
                }
                
                // Must have gone as far as we can go with this seed segment. Store segment at
                // lowest plane.
                if(tmpFlag==false) {mTemp.push_back(Temp[jk]);}
              }
            } // End while BegTemp.size()>0
	    
	    
            // Set the flags back to zero, for when we carry out actual propagation, below.
            for(unsigned int jk=0; jk<TempBank1.size(); ++jk) {TempBank1[jk]->SetTmpTrkFlag(0);}

            // Now see how far we can propagate forwards to higher plane numbers
            pTemp.clear(); EndTemp.clear();
	    
            // Put the seed segment into EndTemp then start while loop
            EndTemp.push_back(SeedSegments[ij]); TempBank1.clear();
	    //            cout <<"Trkfinderformtrk : 12x "<<EndTemp.size()<< endl;

            while(EndTemp.size()>0) {
              Temp.clear();
              
              // Copy any segments into a temporary holder and loop over them,
              // gradually propagating to higher plane numbers
              Temp=EndTemp; EndTemp.clear();
             
              for(unsigned int jk=0; jk<Temp.size(); ++jk) {
                tmpFlag=false;
                
                // For each segment, loop over their end matches
                for(unsigned int kl=0; kl<Temp[jk]->GetNMatchSegEnd(); ++kl) {
                  InoTrackSegment* Sege = Temp[jk]->GetMatchSegEnd(kl);
		  //		  cout <<"Trkfinderformtrk : 11x "<<EndTemp.size()<<" "<<jk<<" "<<kl<<" "<<Sege->GetUID()<< endl;


                  // If the matched segment is track-like and has zero TmpTrkFlag, 
                  // push back into EndTemp, so the while loop continues.
                  if(Sege->GetUID()>-1 && Sege->GetUID()<3) {
                    tmpFlag=true;
                    
                    if(Sege->GetTmpTrkFlag()<1) {
                      EndTemp.push_back(Sege);
                      Sege->SetTmpTrkFlag(1); // Temporarily mark segments used
                      TempBank1.push_back(Sege);
                    }
                  }
                }
                
                // Must have gone as far as we can go with this seed segment. Store segment at
                // highest plane.
                if(tmpFlag==false) {pTemp.push_back(Temp[jk]);}
              }
            } // End while EndTemp.size()>0
	    //            cout <<"Trkfinderformtrk : 13x "<< TempBank1.size()<<endl;

            // Set the flags back to zero, for when we carry out actual propagation, below.
            for(unsigned int jk=0; jk<TempBank1.size(); ++jk) {TempBank1[jk]->SetTmpTrkFlag(0);}
            
            // Find the maximum number of planes we can span by propagating the seed segment
            // back and forth.
            npts=-1;
            for(unsigned int jk=0; jk<mTemp.size(); ++jk) {
              for(unsigned int kl=0; kl<pTemp.size(); ++kl) {
                if((1+pTemp[kl]->GetEndZPlane()-mTemp[jk]->GetBegZPlane())>npts) {
                  npts=1+pTemp[kl]->GetEndZPlane()-mTemp[jk]->GetBegZPlane();
                }
              }
            }
            
            SeedSegments[ij]->SetNPlanes(npts);
            if(npts>nplane) {nplane=npts;}
	    
          } // End loop over seed segments stored above


          // Store the seed segments that lead to the biggest propagation in BestSeedSegments
          if(nplane>0) {
            for(unsigned int ij=0; ij<SeedSegments.size(); ++ij) {
              if(SeedSegments[ij]->GetNPlanes()>nplane-3) {BestSeedSegments.push_back(SeedSegments[ij]);}       
            }
          }
          
        } // End finding final seed segments  SeedSegments.size()>0

        // Order contents of BestSeedSegments
        SeedSegments=BestSeedSegments; BestSeedSegments.clear();
        
        for(unsigned int ij=0; ij<SeedSegments.size(); ++ij) {
          MostClusters=0;
          InoTrackSegment* LargestSeg = 0;        
          
          for(unsigned int jk=0; jk<SeedSegments.size(); ++jk) {
            AlreadyAdded=false;
            
            for(unsigned int kl=0; kl<BestSeedSegments.size(); ++kl) {
              if(BestSeedSegments[kl]==SeedSegments[jk]) {AlreadyAdded=true;}
            }
            if(AlreadyAdded) {continue;}

            if(SeedSegments[jk]->GetEntries()>MostClusters) {
              MostClusters=SeedSegments[jk]->GetEntries();
              LargestSeg=SeedSegments[jk];
            }
          }

          if(LargestSeg) {BestSeedSegments.push_back(LargestSeg);}
        }

        // Now use the final seed segments to actually PROPAGATE back and forth, 
        // setting the TmpTrkFlags for the segments used in the propagation.
        TempBank2.clear();

        // Loop over the seed segments, which are stored in BestSeedSegments
        for(unsigned int ij=0; ij<BestSeedSegments.size(); ++ij) {
          InoTrackSegment* Seg0 = BestSeedSegments[ij];
          // If we make a track, we will definitely include the seed segment, so set
          // it's TmpTrkFlag to 3.
          Seg0->SetTmpTrkFlag(3); 
          TrackFlag=false;
	  
          TempBank1.clear(); TempBank1.push_back(Seg0); 
	  
          // Track backwards (1)
          mTemp.clear(); BegTemp.clear();
	  
          // Push the seed segment into BegTemp and begin while loop
          BegTemp.push_back(Seg0);

          while(BegTemp.size()>0) {
            Temp.clear();

            // Copy any segments into a temporary holder and loop over them,
            // gradually propagating back to lower plane numbers
            Temp=BegTemp; BegTemp.clear();

            for(unsigned int jk=0; jk<Temp.size(); ++jk) {
              tmpFlag=false;
              InoTrackSegment* Segtmp = Temp[jk];
              
              // For each segment, loop over their beginning matches
              for(unsigned int kl=0; kl<Segtmp->GetNMatchSegBeg(); ++kl){
                InoTrackSegment* Segb = Segtmp->GetMatchSegBeg(kl);
                
                // If the matched segment is track-like and has zero TmpTrkFlag, 
                // push back into BegTemp, so the while loop continues.
                if(Segb->GetUID()>-1 && Segb->GetUID()<3) {
                  tmpFlag=true;
                  
                  if(Segb->GetTmpTrkFlag()<1) {
                    BegTemp.push_back(Segb);
                    Segb->SetTmpTrkFlag(1); // Mark the segments with TmpTrkFlag 1
                    TempBank1.push_back(Segb);
                  }
                }
              }
	      
              // Must have gone as far as we can go with this seed segment. Store segment at
              // lowest plane.
              if(tmpFlag==false) {mTemp.push_back(Segtmp);}
            }
          }
          
          // Note the first plane of the potential preliminary track
          nplane=999;
          for(unsigned int jk=0; jk<mTemp.size(); ++jk) {
            if(mTemp[jk]->GetBegZPlane()<nplane) {nplane=mTemp[jk]->GetBegZPlane();}
          }
	  
	  
          // For the longest possible paths, move back towards the seed segment,
          // setting TmpTrkFlags to 2 for the segments used.
	  
          // Find the segments reached by the longest propagations
          BegTemp.clear();
          for(unsigned int jk=0; jk<mTemp.size(); ++jk) {
            if(mTemp[jk]->GetBegZPlane()<nplane+3) {BegTemp.push_back(mTemp[jk]);}
          }

          // From these segments, propagate forwards and set the flags.
          while(BegTemp.size()>0) {
            Temp.clear();
            
            for(unsigned int jk=0; jk<BegTemp.size(); ++jk) {
              if(BegTemp[jk]->GetTmpTrkFlag()<2) {
                BegTemp[jk]->SetTmpTrkFlag(2); // Mark the segments with TmpTrkFlag 2
                Temp.push_back(BegTemp[jk]);
              }
              
              if(BegTemp[jk]->GetTrkFlag()<1) {
                BegTemp[jk]->SetTrkFlag(1);    // Also mark the segments with TrkFlag 1
                TempBank2.push_back(BegTemp[jk]);
                TrackFlag=true;
              }
            }
            BegTemp.clear();
            
            for(unsigned int jk=0; jk<Temp.size(); ++jk) {
              for(unsigned int kl=0; kl<Temp[jk]->GetNMatchSegEnd(); ++kl) {
                if(Temp[jk]->GetMatchSegEnd(kl)->GetTmpTrkFlag()>0) {
                  BegTemp.push_back(Temp[jk]->GetMatchSegEnd(kl));
                }
              }
            }
          }
          
          // Track forwards (1)
          pTemp.clear(); EndTemp.clear();
          
          // Push seed segment into EndTemp and begin while loop
          EndTemp.push_back(Seg0);
	  
          while(EndTemp.size()>0) {
            Temp.clear();
	    
            // Copy any segments into a temporary holder and loop over them,
            // gradually propagating to higher plane numbers
            Temp=EndTemp; EndTemp.clear();
	    
            for(unsigned int jk=0; jk<Temp.size(); ++jk) {
              tmpFlag=false;
              InoTrackSegment* Segtmp = Temp[jk];
              
              // For each segment, loop over their end matches
              for(unsigned int kl=0; kl<Segtmp->GetNMatchSegEnd(); ++kl){
                InoTrackSegment* Sege = Segtmp->GetMatchSegEnd(kl);
                
                // If the matched segment is track-like and has zero TmpTrkFlag, 
                // push back into BegTemp, so the while loop continues.
                if(Sege->GetUID()>-1 && Sege->GetUID()<3) {
                  tmpFlag=true;
                  
                  if(Sege->GetTmpTrkFlag()<1) {
                    EndTemp.push_back(Sege);
                    Sege->SetTmpTrkFlag(1); // Mark the segments with TmpTrkFlag 1
                    TempBank1.push_back(Sege);
                  }
                }
              }
	      
              // Must have gone as far as we can go with this seed segment. Store segment at
              // highest plane.
              if(tmpFlag==false) {pTemp.push_back(Segtmp);}
            }
          }
	  
          // Note the end plane of potential preliminary track     
          nplane=-999;
          for(unsigned int jk=0; jk<pTemp.size(); ++jk) {
            if(pTemp[jk]->GetEndZPlane()>nplane) {nplane=pTemp[jk]->GetEndZPlane();}
          }
	  
          // For the longest possible paths, move back towards the seed segment,
          // setting TmpTrkFlags to 2 for the segments used.
	  
          // Find the segments reached by the longest propagations
          EndTemp.clear();
          for(unsigned int jk=0; jk<pTemp.size(); ++jk) {
            if(pTemp[jk]->GetEndZPlane()>nplane-3) {EndTemp.push_back(pTemp[jk]);} //GMA original nplane-5 
          }
	  
          // From these segments, propagate back and set the flags.
          while(EndTemp.size()>0) {
            Temp.clear();
            
            for(unsigned int jk=0; jk<EndTemp.size(); ++jk) {
              if(EndTemp[jk]->GetTmpTrkFlag()<2) {
                EndTemp[jk]->SetTmpTrkFlag(2); // Mark the segments with TmpTrkFlag 2
                Temp.push_back(EndTemp[jk]);
              }
              
              if(EndTemp[jk]->GetTrkFlag()<1) {
                EndTemp[jk]->SetTrkFlag(1);    // Also mark the segments with TrkFlag 1
                TempBank2.push_back(EndTemp[jk]);
                TrackFlag=true;
              }
            }
            EndTemp.clear();
            
            for(unsigned int jk=0; jk<Temp.size(); ++jk) {
              for(unsigned int kl=0; kl<Temp[jk]->GetNMatchSegBeg(); ++kl) {
                if(Temp[jk]->GetMatchSegBeg(kl)->GetTmpTrkFlag()>0) {
                  EndTemp.push_back(Temp[jk]->GetMatchSegBeg(kl));
                }
              }
            }
          }

          // End setting TmpTrkFlags
	  
          // From paths of segments with TmpTrkFlag==2, find the best set of
          // segments at the beginning and the best set of segments at the end
	  
          // If propagation from the seed segment has given us a good potential preliminary track
	  //          cout<<"InoTrackFinder : TrackFlag " << TrackFlag << endl;
	  
          if(TrackFlag==true) { 
	    
            // Track backwards (2) and find best set of beginning segments
            TempBegBank.clear(); ObjVectorCounter=0;
            vector<InoTrackSegment*> TmpBegSegBank;
	    
            TmpBegSegBank.push_back(Seg0);
            TempBegBank.push_back(TmpBegSegBank);
            ObjCounter=1;
            
	    
            // Store the series of beginning segments that were flagged above during the propagation.
            // Each entry in TempBegbank will be a vectors of segments (a series of beginning segments 
            // flagged above).
	    
            // TempBegBank will contain vectors of segments representing every possible path back through
            // the segments marked with TmpTrkFlag 2 (possible beginning sections for the preliminary track)
            while(ObjCounter>0 && ObjVectorCounter<fObjVectorMax) {
              ObjCounter=0;
              npts=TempBegBank.size();
	      
              for(int jk=0; jk<npts; ++jk) {
                InoTrackSegment* Segtmp = TempBegBank[jk].back();
                mTemp.clear();
                
                for(unsigned int kl=0; kl<Segtmp->GetNMatchSegBeg(); ++kl) {
                  InoTrackSegment* Segbeg = Segtmp->GetMatchSegBeg(kl);
                  if(Segbeg->GetTmpTrkFlag()>1) {mTemp.push_back(Segbeg);};
                }

                if(mTemp.size()>0) {
                  for(unsigned int kl=1; kl<mTemp.size(); ++kl) {
                    if(ObjVectorCounter<fObjVectorMax) { 
                      vector<InoTrackSegment*> NewObj = TempBegBank[jk];
                      
                      NewObj.push_back(mTemp[kl]);
                      TempBegBank.push_back(NewObj);

                      ObjVectorCounter++;
                    } else {
		      cout<<"InoTrackFinder : *** CRAZY MEMORY GUZZLING *** " << endl; break;
		    }
                  }

                  if(ObjVectorCounter>=fObjVectorMax) {break;}

                  TempBegBank[jk].push_back(mTemp[0]); ObjCounter++;
                }
              }
            } // End while ObjCounter>0

            // Find beginning section with best score and store it
	    //            cout<< "InoTrackFinder : beginning score " << endl;

            id=-1; TopScore=-1.;
            for(unsigned int jk=0; jk<TempBegBank.size(); ++jk) {
              Score = Seg0->GetScore(&TempBegBank[jk],0);
              if(Score>TopScore) {TopScore=Score; id=jk;}
            }
	    
            if(id!=-1) {BegBank.push_back(TempBegBank[id]);}
            TempBegBank.clear();

            // Track forwards (2) and find best set of end segments
            TempEndBank.clear(); ObjVectorCounter=0;
            vector<InoTrackSegment*> TmpEndSegBank;
	    
            TmpEndSegBank.push_back(Seg0);
            TempEndBank.push_back(TmpEndSegBank);
            ObjCounter=1;
            
            // Store the series of end segments that were flagged above during the propagation.
            // Each entry in TempEndbank will be a vectors of segments (a series of end segments 
            // flagged above).

            // TempBegBank will contain vectors of segments representing every possible path forward through
            // the segments marked with TmpTrkFlag 2 (possible end sections for the preliminary track)
            while(ObjCounter>0 && ObjVectorCounter<fObjVectorMax) {
              ObjCounter=0;
              npts=TempEndBank.size();
              
              for(int jk=0; jk<npts; ++jk) {
                InoTrackSegment* Segtmp = TempEndBank[jk].back();
		
                pTemp.clear();
                
                for(unsigned int kl=0; kl<Segtmp->GetNMatchSegEnd(); ++kl) {
                  InoTrackSegment* Segend = Segtmp->GetMatchSegEnd(kl);
                  if(Segend->GetTmpTrkFlag()>1) {pTemp.push_back(Segend);};
                }
		
                if(pTemp.size()>0) {
                  for(unsigned int kl=1; kl<pTemp.size(); ++kl) {
                    if(ObjVectorCounter<fObjVectorMax) {
                      vector<InoTrackSegment*> NewObj = TempEndBank[jk];
                      
                      NewObj.push_back(pTemp[kl]);
                      TempEndBank.push_back(NewObj);

                      ObjVectorCounter++;
                    } else {
		      cout<<"InoTrackFinder : *** CRAZY MEMORY GUZZLING *** " << endl; break;
		    }
                  }

                  if(ObjVectorCounter>=fObjVectorMax) {break;}

                  TempEndBank[jk].push_back(pTemp[0]); ObjCounter++;
                }
              }
            } // End while ObjCounter>0

            // Find end section with best score and store it
	    //            cout << "end score " << endl;
	    
            id=-1; TopScore=-1.;
            for(unsigned int jk=0; jk<TempEndBank.size(); ++jk) {
              Score = Seg0->GetScore(0,&TempEndBank[jk]);
              if(Score>TopScore) {TopScore=Score; id=jk;} 
            }
	    
            if(id!=-1) {EndBank.push_back(TempEndBank[id]);}
            TempEndBank.clear(); // maybe clear individual vectors of segments
          }
	  
          
          // Clear up
          for(unsigned int jk=0; jk<TempBank1.size(); ++jk) {TempBank1[jk]->SetTmpTrkFlag(0);}
	  
        } // End loop over BestSeedSegments
        for(unsigned int jk=0; jk<TempBank2.size(); ++jk) {TempBank2[jk]->SetTrkFlag(0);}
	
        // Using the best beginning sections and the best end sections from all 
        // seed segments, find the BEST COMPLETE preliminary track
        if(BegBank.size()>0 && EndBank.size()>0) {
	  
          // Calculate a score, including the seed segment and best beginning/end paths,
          // to find the 'best' complete track
	  //          cout << "overall score " << endl;
	  
          id=-1; TopScore=-1.;
          for(unsigned int ij=0; ij<BegBank.size(); ++ij) {
            Score = BegBank[ij][0]->GetScore(&BegBank[ij],&EndBank[ij]);
            if(Score>TopScore) {TopScore=Score; id=ij;}
          }
	  
          // If we have found a best track, add together the segments to form
          // a segment that is the preliminary track, giving this UID 3. Mark segments  
          // used by setting their UID to -1.
          if(1+id>0) {
            // Get the seed segment
            InoTrackSegment* Seg0 = BegBank[id][0];
	    
            // Set segment UIDs to -1 if they are part of the track, so they
            // can't be used in another track.
            for(unsigned int ij=1; ij<BegBank[id].size(); ++ij) {
              InoTrackSegment* Segbeg = BegBank[id][ij];
              Seg0->AddSegment(Segbeg);
              Segbeg->SetUID(-1);
            }

            for(unsigned int ij=1; ij<EndBank[id].size(); ++ij) {
              InoTrackSegment* Segend = EndBank[id][ij];
              Seg0->AddSegment(Segend);
              Segend->SetUID(-1);
            }
            
            // Mark long segment as a preliminary track by setting its UID to 3
            Seg0->SetUID(3);

            // Mark all the clusters in preliminary track with TrkFlag 3
            for(unsigned int ij=0; ij<Seg0->GetEntries(); ++ij) {
              Seg0->GetCluster(ij)->SetTrkFlag(3);
            }

            // Store for use in Join tracks method
	    PossibleJoins[View].push_back(Seg0);
	    
            // Can now look for more preliminary tracks
            Cont=true; ntrks++;
          }
        }
	
        // Clean up
        BegBank.clear(); EndBank.clear();
	
        
      } // End while Cont==true loop
    } // End loop over planeviews
  } // End check to see if there is anything in ViewSegBank
  
  return;
}

void InoTrackFinder::JoinTracks() {
  //  cout << " *** JOIN preliminary TRACKS *** " << endl;
  
  //MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer; 
  int MxPlaneGap=15;  //GMA Maximum layer of allowed gap
  double ShortestDist; double LongestAllowedDist=3.; //GMA need optimisation   asm:LongestAllowedDist (in metres )
  double DeltaXPos, DeltaYPos, DeltaZPos, Dist;
  int DeltaPlane, Length, UID, NumAgree;
  
  InoTrackSegment* Seg1ToAdd=0;
  InoTrackSegment* Seg2ToAdd=0;
  bool PossibleAssocs, Cont; //, Overlap, DiffGrad;
  
  // First, make the best possible joins between the tracks
  for(int View=0; View<1; ++View) {
    if(PossibleJoins[View].size()<2) continue;
    PossibleAssocs=true;

     while(PossibleAssocs) {
      PossibleAssocs=false; 
      ShortestDist=LongestAllowedDist; Seg1ToAdd=0; Seg2ToAdd=0;

      // Consider all pairs of preliminary tracks and find the single best association.
      for(unsigned int ij=0; ij<PossibleJoins[View].size(); ++ij) {
	//if(PossibleJoins[View][ij].GetEntries()>10) LongestAllowedDist=4. ;   //asm:080311
	//ShortestDist=LongestAllowedDist; Seg1ToAdd=0; Seg2ToAdd=0; //asm:080311
        InoTrackSegment* Seg1 = PossibleJoins[View][ij];
	
        if(Seg1->GetUID()>2) {
	  
	  //          for(unsigned int jk=0; jk<PossibleJoins[View].size(); ++jk) {
          for(unsigned int jk=ij+1; jk<PossibleJoins[View].size(); ++jk) {
	    
            InoTrackSegment* Seg2 = PossibleJoins[View][jk];
	    
            if(Seg2->GetUID()>2 && Seg1!=Seg2) {         
	      
              // Calculate distance between segments
              DeltaXPos=fabs(Seg2->GetBegXPos()-Seg1->GetEndXPos());
	      DeltaYPos=fabs(Seg2->GetBegYPos()-Seg1->GetEndYPos());
              DeltaZPos=fabs(Seg2->GetBegZPos()-Seg1->GetEndZPos());
	      
              DeltaPlane=Seg2->GetBegZPlane()-Seg1->GetEndZPlane();
	      
              Dist=pow((DeltaXPos*DeltaXPos + DeltaYPos*DeltaYPos + DeltaZPos*DeltaZPos) ,0.5);

	      double seg1Xang =Seg1->GetBegXDir();
	      double seg2Xang =Seg2->GetEndXDir();
	      
	      double seg1Yang =Seg1->GetBegYDir();
	      double seg2Yang =Seg2->GetEndYDir();
	      
	      double Tr1EndX = 0, Tr1EndY = 0, Tr1EndZ = 0, Tr1EndTx = 0, Tr1EndTy = 0;
	      //		Tr1EndP = 0;
	      //	      double Tr2BegX = 0, Tr2BegY = 0, Tr2BegZ = 0, Tr2BegTx = 0, Tr2BegTy = 0, Tr2BegP = 0;
	      
	      Tr1EndX = Seg1->GetEndXPos();
	      Tr1EndY = Seg1->GetEndYPos();
	      Tr1EndZ = Seg1->GetEndZPos();
	      Tr1EndTx= Seg1->GetEndXDir();
	      Tr1EndTy= Seg1->GetEndYDir();
	      //Tr1EndP = fTrackCand->GetEndQP();
	      
	      double the1 = acos(1./pow(1+pow(seg1Xang,2.)+pow(seg1Yang,2.),0.5));
	      double phi1 = atan2(seg1Yang,seg1Xang);
	      
	      double the2 = acos(1./pow(1+pow(seg2Xang,2.)+pow(seg2Yang,2.),0.5));
	      double phi2 = atan2(seg2Yang,seg2Xang);
	      
	      double cosang = sin(the1)*sin(the2)*cos(phi1-phi2)+cos(the1)*cos(the2);
	      
	      // cout <<"InoTrackFinder::JoinTracks() "<<seg1Xang<<" "<<seg2Xang<<" "<<seg1Yang<<" "<<seg2Yang<<" "<<the1<<" "<<phi1<<" "<<the2<<" "<<phi2<<" "<<cosang<<endl;
	      
              if(DeltaPlane>=0 &&DeltaPlane <=MxPlaneGap && Dist<ShortestDist 
		 && cosang > 0.5 ) {  // 0.9 
		
                if(Seg1->IsAssoc(Seg2)) {Seg1ToAdd=Seg1; Seg2ToAdd=Seg2; ShortestDist=Dist;} //asm: see if IsAssoc takes care of all possible segment Joins
              } else if(DeltaPlane<=0 && DeltaPlane >=-MxPlaneGap && Dist<ShortestDist) {
                if(Seg2->IsAssoc(Seg1)) {
		  Seg1ToAdd=Seg1; Seg2ToAdd=Seg2; ShortestDist=Dist;
		}
              }
            }
          }
        }
      }
      
      if(ShortestDist<LongestAllowedDist && Seg2ToAdd!=0 && Seg1ToAdd!=0) { 
	if(pAnalysis->isXtermOut==1){
	  cout << "Missing preliminary track assocation: Seg1 " 
	       << Seg1ToAdd->GetBegZPlane() << " " << Seg1ToAdd->GetEndZPlane() 
	       << " (" << Seg1ToAdd->GetBegXPos() << "," << Seg1ToAdd->GetEndXPos() 
	       << " " << Seg1ToAdd->GetBegYPos() << "," << Seg1ToAdd->GetEndYPos() 
	       << "), Seg2 " << Seg2ToAdd->GetBegZPlane() 
	       << " " << Seg2ToAdd->GetEndZPlane() << " (" << Seg2ToAdd->GetBegXPos() 
	       << "," << Seg2ToAdd->GetEndXPos() 
	       << ", "<<Seg2ToAdd->GetBegYPos() 
	       << "," << Seg2ToAdd->GetEndYPos() << endl;
	} // if(isXtermOut)
	
	Seg1ToAdd->AddSegment(Seg2ToAdd); Seg2ToAdd->SetUID(-1); PossibleAssocs=true;
      }
     }
     
     // Now find the longest possible preliminary track and remove 2d tracks that overlap.
     Cont=true; 
     
     while(Cont) {
       Cont=false;
       
       InoTrackSegment* BestTrk = 0;
       Length=0;
       
       for(unsigned int ij=0; ij<PossibleJoins[View].size(); ++ij) {
	 InoTrackSegment* Seg = PossibleJoins[View][ij];
	 if(Seg->GetUID()>0 && Seg->GetUID()<4 && (Seg->GetEndZPlane()-Seg->GetBegZPlane())>Length) {
	   BestTrk=Seg; Length=Seg->GetEndZPlane()-Seg->GetBegZPlane();
	 }
       }
       
       if(BestTrk) {
	 BestTrk->SetUID(4);
	 
	 for(unsigned int ij=0; ij<BestTrk->GetEntries(); ++ij) {
	   InoCluster* BestTrkClust = BestTrk->GetCluster(ij);
	   
	   if(BestTrkClust->GetZPlane()==BestTrk->GetBegZPlane() || BestTrkClust->GetZPlane()==BestTrk->GetEndZPlane()) {continue;}
	   
	   for(unsigned int jk=0; jk<PossibleJoins[View].size(); ++jk) {
	     InoTrackSegment* Seg = PossibleJoins[View][jk];
	     NumAgree=0; 
	     
	     if(Seg->GetUID()<0 || Seg->GetUID()>3) {continue;}
	     
	     for(unsigned int kl=0; kl<Seg->GetEntries(); ++kl) {
	       InoCluster* SegClust = Seg->GetCluster(kl);   // asm: this condition needs some modification.... 
	       // should increment NumAgree if the clusters are just  seperated by one strips diference.  
	       
	       if(SegClust==BestTrkClust) {NumAgree++;}
	       
	       if(NumAgree>2) {Seg->SetUID(-1); break;}
	     }
	   }
	 }
       }
       
       // Conditions for continuing here
       for(unsigned int ij=0; ij<PossibleJoins[View].size(); ++ij) {
	 UID=PossibleJoins[View][ij]->GetUID();
	 
	 if(UID>0 && UID<4) {Cont=true;}
       }
       
     }
     
     
     // Reset UIDs
     for(unsigned int ij=0; ij<PossibleJoins[View].size(); ++ij) {
       if(PossibleJoins[View][ij]->GetUID()==4) {PossibleJoins[View][ij]->SetUID(3);}
     }
     
     PossibleJoins[View].clear();
  }
  
  
  // Print out the list of preliminary tracks
  for(int View=0; View<1; ++View) {
    //        cout << "InoTrackFinder ViewSegBank: " << View <<" "<<ViewSegBank[View].size()<< endl;
    for(unsigned ij=0; ij<ViewSegBank[View].size(); ++ij) {
      InoTrackSegment* Seg = ViewSegBank[View][ij];
      if(Seg->GetUID()>2) {
	if(pAnalysis->isXtermOut==1){
	  cout << "InoTrackFinder  ij=" << ij
	       << "  begpln: "  << Seg->GetBegZPlane()
	       << "  begxpos: " << Seg->GetBegXPos()
	       << "  begypos: " << Seg->GetBegYPos()
	       << "  endpln: "  << Seg->GetEndZPlane()
	       << "  endxpos: " << Seg->GetEndXPos() 
	       << "  endypos: " << Seg->GetEndYPos() 
	       << "  Entries: "<< Seg->GetEntries() << endl;
	}// if(isXtermOut){
	
      }
    }
    
    for(unsigned ij=0; ij<TempTrack[View].size(); ++ij) {
      InoTrackSegment* Seg = TempTrack[View][ij];
      if(Seg->GetUID()>2) {
	if(pAnalysis->isXtermOut==1){
	  cout << "InoTrackFinder TempTrack  ij=" << ij
	       << "  begpln: "  << Seg->GetBegZPlane()
	       << "  begxpos: " << Seg->GetBegXPos()
	       << "  begypos: " << Seg->GetBegYPos()
	       << "  endpln: "  << Seg->GetEndZPlane()
	       << "  endxpos: " << Seg->GetEndXPos() 
	       << "  endypos: " << Seg->GetEndYPos()
	       << endl;
	}//if(isXtermOut){   
      } else { //GMA 24/02/08
	TempTrack[View].erase(TempTrack[View].begin()+ij); ij--;
      }	
    }
  }
  return;
}

//==============================================================:FormFinalTracks=====================================================================

void InoTrackFinder::FormFinalTracks( )
{
  // The Final track segments are now stored in the TempTracks container, with
  // U and V segments stored separately. In this method we look at the consituent
  // clusters and hits in these segments, performing fits to choose the final
  // strips for the track.
  
//MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
//  cout << " *** Formation of Final tracks *** " << endl;

  int ktemp, View=0;
  double Xswx, Xswy, Xswx2, Xswxy, Xsw, ProjectedXPos, XPosWindow(0.);
  double Yswx, Yswy, Yswx2, Yswxy, Ysw, ProjectedYPos, YPosWindow(0.);

  double Xw, Xx, Zz, Xm, Xc;
  double Yw, Yy,     Ym, Yc;


  //  int UStrips, VStrips, TrackStrips, TrackPlanes;
  int TrackStrips, TrackPlanes;

  int begPlaneView[2], endPlaneView[2];
  int begplane1, begplane2, endplane1, endplane2, maxplane, Plane;
  double dScore, Score; int id;
  
  int ShwOrTrkList[500];
  int ViewList[500];
  vector<InoHit*> HitList[500];
  //  vector<InoCluster*> ClustList[500]; //GMA14
  //  vector<InoTrack*> FinalTrackTempHolder[2];
  
  for(int ij=0; ij<500; ++ij) {ShwOrTrkList[ij]=-1; ViewList[ij]=-1; ClustList[ij].clear();}
  for(int ij=0; ij<2; ++ij) {FinalTrackTempHolder[ij].clear();}

  // Main loop is over the pairs of track segments, stored above
  for(unsigned int ij=0; ij<TempTrack[0].size(); ++ij) {
    //    cout<< "InoTrackFinder making track " << ij << endl;
    
    // Find the beginning and end planes for each of this pair of segments.
    // Order them so that we have:
    // begplane2, begplane 1, endplane 1, endplane 2.
    begplane2=-1; endplane2=-1; begplane1=-1; endplane1=-1;
    
    for(View=0; View<1; ++View) {
      InoTrackSegment* seg = TempTrack[View][ij];
      
      if(begplane2<0 || seg->GetBegZXPlane()<begplane2) { 
	begplane1=begplane2; begplane2=seg->GetBegZXPlane(); 
      } else { 
	begplane1=seg->GetBegZPlane(); 
      }
      
      if(endplane2<0 || seg->GetEndZXPlane()>endplane2) { 
	endplane1=endplane2; endplane2=seg->GetEndZXPlane(); 
      } else { 
	endplane1=seg->GetEndZPlane(); 
      }
      //GMA for the time beig use only X plane, but need the combination of both X and Y
      if(begplane2<0 || seg->GetBegZYPlane()<begplane2) { 
	begplane1=begplane2; begplane2=seg->GetBegZYPlane(); 
      } else { 
	begplane1=seg->GetBegZPlane(); 
      }
      
      if(endplane2<0 || seg->GetEndZYPlane()>endplane2) { 
	endplane1=endplane2; endplane2=seg->GetEndZYPlane(); 
      } else { 
	endplane1=seg->GetEndZPlane(); 
      }
      
      begPlaneView[View]=-1; endPlaneView[View]=-1;
      
    }
    //    cout << "endplane1:" << endplane1 << "  endplane2:" << endplane2 <<endl;
    maxplane=1+endplane2;
    
 
    // Unpack the clusters from the segments.
    // For each track, should only have one cluster on a plane.
    for(View=0; View<1; ++View) {
      InoTrackSegment* Seg = TempTrack[View][ij];     // ij ? 
      
      for(unsigned int jk=0; jk<Seg->GetEntries(); ++jk) {
        InoCluster* Clust = Seg->GetCluster(jk);
	
        Plane=Clust->GetZPlane();
        if(Plane>=begplane2 && Plane<maxplane) {ClusterList[Plane].push_back(Clust);  //}
	  //	cout <<"final track "<< jk<<" "<< Plane<<" "<<Clust->GetEntries() <<" "<< Clust->GetBegXPos()<<" "<<Clust->GetBegYPos()<<endl;
        }
      }
    }
    
    // Trim the segments. Check we haven't tracked much further
    // in one view than the other, and adjust as necessary.
    
    // Maximum discrepancy allowed at each end is 3 planes. As only
    // clusters on planes greater than begplane1 are considered.
    int PrevPln, NextPln;
    bool LeavesAtEdge=false;
    
    // Check hasn't just left edge of detector in one view.
    //GMA  Modify it with proper geomtry file
    
    if(ModuleType==1 && (ClusterList[begplane1][0]->GetEndXPos()> 24.20 ||  // 33.9 ||
			 ClusterList[begplane1][0]->GetBegXPos()<-24.20 || //-33.9 ||
			 ClusterList[begplane1][0]->GetEndYPos()>  8 ||  // 33.9 || 33.9 ||
			 ClusterList[begplane1][0]->GetBegYPos()< -8 )  //-33.9 ||-33.9)
       ) {begplane1=begplane2+1; LeavesAtEdge=true;}
    
    
    if(ModuleType==1 && (ClusterList[endplane1][0]->GetEndXPos()> 24.20 ||  //33.9 ||
			 ClusterList[endplane1][0]->GetBegXPos()<-24.20 || //-33.9 ||
			 ClusterList[endplane1][0]->GetEndYPos()>  8 ||  // 33.9 ||
			 ClusterList[endplane1][0]->GetBegYPos()< -8 )   //-33.9  )
       ) {endplane1=endplane2-1; LeavesAtEdge=true;}
    
    //    cout <<"InoTrackFinder :"<< begplane2<< "  (" << begplane1 << ") -> (" << endplane1 << ") " << endplane2 << endl;
    //    cout <<"LeavesAtEdge "<<(int)LeavesAtEdge<<endl;
    
    if(!LeavesAtEdge) {
      // For the beginning
      // Jump back 2 planes before begplane1 (unless this takes us past begplane2)
      
      //-----------------------------------------------------------------------------------------------------    
      begplane1=begplane1-2; PrevPln=-1;
      Plane=begplane1; if(Plane<begplane2) {Plane=begplane2;}
      
      // Working forwards from this plane, find the first plane that has
      // a cluster.
      while(Plane<maxplane) {
        if(ClusterList[Plane].size()>0) {PrevPln=Plane; break;}
        else {Plane++;}
      }
      // Look at the next 2 planes in the same view as this cluster
      // and see if they have clusters too. Try and veto geometries such as
      //
      //                bp1
      // u       x   x   o   o
      // v         o   x   *   o
      //          bp2      * == shw-like, or gap
      //
      // In this case, begplane1 is set to begplane 2, so the lone hit at 
      // begplane 2 is later vetoed.
      if(PrevPln>-1 && PrevPln+2<maxplane && ClusterList[PrevPln+1].size()==0) {
        if(ClusterList[PrevPln+2].size()>0) {
          InoCluster* Clust = ClusterList[PrevPln+2][0];
	  
          if(Clust->GetTrkPlnFlag()>0 && Clust->GetShwPlnFlag()==0) {
            PrevPln=begplane1;
          }
        }
        begplane1=PrevPln;
      }
      
      //-----------------------------------------------------------------------------------------------------    
      // Similarly, for the end 
      // Jump forward 2 planes after endplane1 (unless this takes us past endplane2)
      endplane1=endplane1+2; NextPln=-1;
      Plane=endplane1; if(Plane>maxplane-1) {Plane=maxplane-1;} 
      
      //      cout <<"PLane "<<Plane<<" "<<begplane2<<" "<<NextPln<<" "<<maxplane<<endl;
      
      // Working backwards from this plane, find the first plane that has
      // a cluster.
      while(Plane>begplane2-1) {
        if(ClusterList[Plane].size()>0) {NextPln=Plane; break;}
        else {Plane--;}
      }
      
      // Look at the next 2 planes in the same view as this cluster
      // and see if they have clusters too, as above.
      if(NextPln>-1 && NextPln-2>(begplane2-1) && ClusterList[NextPln-1].size()==0) {
        if(ClusterList[NextPln-2].size()>0) {
          InoCluster* Clust = ClusterList[NextPln-2][0];
          if(Clust->GetTrkPlnFlag()>0 && Clust->GetShwPlnFlag()==0) {
            NextPln=endplane1;
          }
        }
        endplane1=NextPln; //if plane NextPln-2 has a trach\k like cluster then endPlane1 is not ignored
      }
    }
    
    //-----------------------------------------------------------------------------------------------------    
    
    // Count the strips and record planeviews, number of occupied planes, etc.
    //    UStrips=0; VStrips=0; 
    TrackStrips=0;
    
    //  cout <<"begplane2 "<<begplane2<<" "<<maxplane<<endl;
    // Look at the region over which the track is now defined
    for(int jk=begplane2; jk<maxplane; ++jk) {  
      
      if( //  (jk>=begplane1 && jk<=endplane1) &&
	 ClusterList[jk].size()>0 ) {
	
        InoCluster* Clust = ClusterList[jk][0];
        
	//        View=Clust->GetPlaneView();
	View = 0;      //note View is set to 0
	ViewList[jk]=View;              
        ShwOrTrkList[jk]=1;                  // to keep track of clusters
        
        // Record first and last planes in each view.
	if(begPlaneView[View]<0 || jk<begPlaneView[View]) {begPlaneView[View]=jk;}
	if(endPlaneView[View]<0 || jk>endPlaneView[View]) {endPlaneView[View]=jk;}
        
	//        if(View==0) {UStrips++;}
	//        if(View==1) {VStrips++;}
        
        if(Clust->GetTrkPlnFlag()>0) {TrackStrips++;}
      }
    }
  // print out clusterList at this point   asm  
    //-------------------------------------------------------------------------------------
    
    // If we have at least 3 strips in each view and at least 2 clusters
    // on track-like planes, proceed with the fits.
    //    if(UStrips>2 && VStrips>2 && (CambridgeAnalysis==false || TrackStrips>2) ) {
    if(TrackStrips>2 ) {
      
      // FLAG PLANES which we think are track-like
      for(View=0; View<1; ++View) {
	//	cout <<"InoTrackFinder : begplane2 "<<begplane2<<" "<<maxplane<<endl; 
        for(Plane=begplane2; Plane<maxplane; ++Plane) {
	  
          if(ShwOrTrkList[Plane]>0) {
            InoCluster* Clust0 = ClusterList[Plane][0];
	    
	    // Always set flag to 3 for the +/- 5 plane ND track segments
            if(Clust0->GetNDFlag()==2) {ShwOrTrkList[Plane]=3;}   // NDFlag =1 
	    
            // Find next planes, in the same view, where we have a cluster, both forwards
            // and backwards.
            ktemp=0; PrevPln=-1;
            while(Plane-ktemp>begplane2) {
              ktemp++;
	      if(ViewList[Plane-ktemp]==View) {PrevPln=ktemp; break;}
            }
	    
            ktemp=0; NextPln=-1;
            while(Plane+ktemp<maxplane-1) {
              ktemp++;
	      if(ViewList[Plane+ktemp]==View) {NextPln=ktemp; break;}
            }
	    
            // If there are planes forwards and backwards in the same view...
            if(PrevPln>0 && NextPln>0) {
              InoCluster* NextClust = ClusterList[Plane+NextPln][0];
              InoCluster* PrevClust = ClusterList[Plane-PrevPln][0];
              
              // If either both clusters are within range for this view, or they are the
              // next plane along and the current plane is a track plane.   
	      
	      //asm: with || (which is now replace by &&  ) this first condition will always be satisfied 
	      //asm: also the next condition will always be satisfied         
              if(Clust0->GetShwPlnFlag()>0)ShwOrTrkList[Plane]=2;            
              if( (PrevClust->GetZPlane()>begPlaneView[View] ||
		   (Clust0->GetZPlane()-PrevClust->GetZPlane()<3 && Clust0->GetTrkPlnFlag()>0)) &&     // asm 2 -> 3
		  (NextClust->GetZPlane()<endPlaneView[View] ||
		   (NextClust->GetZPlane()-Clust0->GetZPlane()<3 && Clust0->GetTrkPlnFlag()>0) )) {      //asm 2 -> 3
		
                // If the clusters have track-like association, identify plane as
                // being track like, by setting entry in ShwOrTrkList to 3.
                
                if(Clust0->IsTrkAssoc(PrevClust,NextClust)>2&&(PrevClust->GetTrkPlnFlag()&&NextClust->GetTrkPlnFlag()&&Clust0->GetTrkPlnFlag())){                   //asm 1 ->2
                  ShwOrTrkList[Plane]=3; 
		  //asm: new 2 lines

                  // If there is a good overlap of strip numbers between the three clusters         //asm: check this condition
                  if( (NextClust->GetBegXPos()-Clust0->GetEndXPos()>-1.1*StripXWidth &&
		       Clust0->GetBegXPos()-PrevClust->GetEndXPos()>-1.1*StripXWidth) ||
		      (NextClust->GetEndYPos()-Clust0->GetBegYPos()<1.1*StripYWidth &&
		       Clust0->GetEndYPos()-PrevClust->GetBegYPos()<1.1*StripYWidth) ) {
		    
                    // If forwards and backwards planes are adjacent in view,
                    // identify these planes as track-like too
		    if((Clust0->GetZPlane()-PrevClust->GetZPlane()<2 || Clust0->GetNDFlag()==2)&&PrevClust->GetTrkPlnFlag()>0){ShwOrTrkList[Plane-PrevPln]=3;}
		    if((NextClust->GetZPlane()-Clust0->GetZPlane()<2 || Clust0->GetNDFlag()==2)&&NextClust->GetTrkPlnFlag()>0){ShwOrTrkList[Plane+NextPln]=3;}
                  }
                } else {
		  // Otherwise, if cluster is small, can identify its plane as track-like.
		  // Don't set anything for forwards and backwards planes.
		  
                  if(Clust0->GetEndXPos()-Clust0->GetBegXPos()<2.1*StripXWidth &&
		     Clust0->GetEndYPos()-Clust0->GetBegYPos()<2.1*StripYWidth &&
		     Clust0->GetTrkPlnFlag()>0) {ShwOrTrkList[Plane]=3;}
                }
              }
            }
          }
        } // End loop over planes
      } // End loop over views

      //      cout<< "InoTrackFinder : *** SORT TRACK HITS *** " << endl;
      
      for(Plane=begplane2;Plane<maxplane; ++Plane) {
	if (pAnalysis->isXtermOut==1 && ViewList[Plane]>-1) {
	  InoCluster* Clust = ClusterList[Plane][0];
	  cout<< "InoTrackFinder : "<<ClusterList[Plane].size()<<" " << Clust->GetHitEntries()<<" "
	      << Clust->GetZPlane() <<" flg " << Clust->GetShwPlnFlag() << " " << Clust->GetTrkPlnFlag()
	      << " " << Clust->GetBegXPos() << " " << Clust->GetEndXPos()
	      << " " << Clust->GetBegYPos() << " " << Clust->GetEndYPos()
	      << " " << ShwOrTrkList[Plane] << endl;
	}
      }
      
      // Perform GLOBAL LINEAR FIT
      int ShowerStart, ShowerEnd, FitStart, FitEnd, PrevTrkPln, NextTrkPln;
      
      for(View=0; View<1; ++View) {
        for(Plane=begplane2; Plane<maxplane; ++Plane) {
	  
          // Currently ShwOrTrkList flag==1 means the plane was part of the Final track.
          // ShwOrTrkList flag==3 means that we also flagged the plane as being track-like, above.
          //
          // Here, we look at the planes we haven't flagged as track-like, find track-like planes
          // on either side and carry out a linear fit through this 'shower-like' region.
          // 
          // Using the fit, clusters in this shower-like region are identified, then hits.
	  if(ViewList[Plane]==View && ShwOrTrkList[Plane]<2) {
	    //	  if(ShwOrTrkList[Plane]<2) { 
	    //            cout<<"InoTrackFinder : *** GLOBAL FIT ***   Plane=" << Plane << endl;

            // Look for next planes in the view that we have flagged as track-like, 
            // forwards and backwards.
            ktemp=0; NextTrkPln=-1;
            while(Plane+ktemp<endplane1 && Plane+ktemp<maxplane-1) {
              ktemp++;
	      if(ViewList[Plane+ktemp]==View && ShwOrTrkList[Plane+ktemp]>2) {NextTrkPln=ktemp; break;}
	      //	      if(ShwOrTrkList[Plane+ktemp]>2) {NextTrkPln=ktemp; break;}
            }
            
            ktemp=0; PrevTrkPln=-1;
            while(Plane-ktemp>begplane1 && Plane-ktemp>begplane2) {
              ktemp++;
	      if(ViewList[Plane-ktemp]==View && ShwOrTrkList[Plane-ktemp]>2) {PrevTrkPln=ktemp; break;}
	      //              if(ShwOrTrkList[Plane-ktemp]>2) {PrevTrkPln=ktemp; break;}
            }
	    
            
            // Find the track-like planes which are just outside of the shower-like region ('ShowerStart' 
            // and 'ShowerEnd' planes away), as well as the planes which are up to 4 planes either side 
            // of the shower-like region ('FitStart' and 'FitEnd' planes away). These are the 
            // track used in the linear fit through the shower.
            if(PrevTrkPln>0) {
	      FitStart=-PrevTrkPln-4; ShowerStart=-PrevTrkPln;  //GMA need optimisation Orginial PrevTrkPln-7
	    } else {
	      FitStart=begplane1-Plane+1; ShowerStart=FitStart;
	    }
	    
            if(Plane+FitStart<begplane2) {FitStart=-Plane+begplane2;}
            if(Plane+ShowerStart<begplane2) {ShowerStart=-Plane+begplane2;}
	    
            if(NextTrkPln>0) {
	      FitEnd=1+NextTrkPln+4; ShowerEnd=1+NextTrkPln;
	    } else {
	      FitEnd=endplane1-Plane; ShowerEnd=FitEnd;
	    }
	    
            if(Plane+FitEnd>maxplane) {FitEnd=maxplane-Plane;}
            if(Plane+ShowerEnd>maxplane) {ShowerEnd=maxplane-Plane;}
	    
            InoCluster* ClustSeed = ClusterList[Plane][0];
	    
            Xm=0.; Xc=0.5*(ClustSeed->GetBegXPos()+ClustSeed->GetEndXPos());
            Xswx=0.; Xswy=0.; Xswx2=0.; Xswxy=0.; Xsw=0.;
	    
            Ym=0.; Yc=0.5*(ClustSeed->GetBegYPos()+ClustSeed->GetEndYPos());
            Yswx=0.; Yswy=0.; Yswx2=0.; Yswxy=0.; Ysw=0.;
	    
            // Carry out the linear fit through the shower, obtaining gradient and intercept.
            for(int kl=FitStart; kl<FitEnd; ++kl) {
	      if(ViewList[Plane+kl]==View) {
		InoCluster* ClustTemp = ClusterList[Plane+kl][0];
		
		for(unsigned int lm=0; lm<ClustTemp->GetHitEntries(); ++lm) {
		  InoHit* HitTemp = ClustTemp->GetHit(lm);
		  
		  Zz=HitTemp->GetZPos();
		  Xx=HitTemp->GetXPos();
		  Xw=HitTemp->GetXPulse()/ClustTemp->GetXPulse();
		  
		  if(ShwOrTrkList[Plane+kl]>2) {Xw=5.*Xw;}
		  else {Xw=0.5*Xw;}
		  Xswx+=Xw*Zz; Xswy+=Xw*Xx; Xswx2+=Xw*Zz*Zz; Xswxy+=Xw*Zz*Xx; Xsw+=Xw;
		  
		  Yy=HitTemp->GetYPos();
		  Yw=HitTemp->GetYPulse()/ClustTemp->GetYPulse();
		  
		  if(ShwOrTrkList[Plane+kl]>2) {Yw=5.*Yw;}
		  else {Yw=0.5*Yw;}
		  Yswx+=Yw*Zz; Yswy+=Yw*Yy; Yswx2+=Yw*Zz*Zz; Yswxy+=Yw*Zz*Yy; Ysw+=Yw;
		}
              }
            }
	    
            if(Xsw>0. && (Xsw*Xswx2-Xswx*Xswx)!=0.) {
              Xm=(Xsw*Xswxy-Xswx*Xswy)/(Xsw*Xswx2-Xswx*Xswx);
              Xc=(Xswy*Xswx2-Xswx*Xswxy)/(Xsw*Xswx2-Xswx*Xswx);
            }
	    
            if(Ysw>0. && (Ysw*Yswx2-Yswx*Yswx)!=0.) {
              Ym=(Ysw*Yswxy-Yswx*Yswy)/(Ysw*Yswx2-Yswx*Yswx);
              Yc=(Yswy*Yswx2-Yswx*Yswxy)/(Ysw*Yswx2-Yswx*Yswx);
            }	   
	    
	    //            cout<< " ShowerStart "<< ShowerStart+Plane <<" ShowerEnd "<<ShowerEnd+Plane-1<<endl;
            // Now loop over the planes in the shower-like region
            for(int kl=ShowerStart; kl<ShowerEnd; ++kl) {
	      if(ViewList[Plane+kl]==View && ShwOrTrkList[Plane+kl]<2) {
		//	      if(ShwOrTrkList[Plane+kl]<1) {	
                // ***First cluster in ClusterList is that from segment in Final track***
                InoCluster* ClustTemp = ClusterList[Plane+kl][0];
                
                // Projected strip number in shower
                ProjectedXPos = Xm*ClustTemp->GetZPos()+Xc;
		ProjectedYPos = Ym*ClustTemp->GetZPos()+Yc;
		
                // Calculate window around projected strip number
                if(Xm>=0.) {XPosWindow = (0.6*StripXWidth)+0.02*Xm;}
                if(Xm<=0.) {XPosWindow = (0.6*StripXWidth)-0.02*Xm;}
		
		if(Ym>=0.) {YPosWindow = (0.6*StripYWidth)+0.02*Ym;}
                if(Ym<=0.) {YPosWindow = (0.6*StripYWidth)-0.02*Ym;}

                // If  1. we have found track-like planes on either side of the shower-like region,
                // or  2. we haven't found any of these track-like planes (entire track in this view must be shower),
                // or  3. current cluster overlaps quite nicely with projected strip number from the fit,
                //
                // then proceed with trying to find the best hits from within the cluster.
                if( (PrevTrkPln>0 && NextTrkPln>0) || (PrevTrkPln<0 && NextTrkPln<0) ||
		    ((ClustTemp->GetEndXPos()>ProjectedXPos-(1.6*StripXWidth) &&
		      ClustTemp->GetBegXPos()<ProjectedXPos+(1.6*StripXWidth)) &&
		     (ClustTemp->GetEndYPos()>ProjectedYPos-(1.6*StripYWidth) &&
		      ClustTemp->GetBegYPos()<ProjectedYPos+(1.6*StripYWidth))) ) {
                  id=-1; dScore=0.; Score=999.9;
		  
                  // Loop over hits in cluster and get a 'score' for each. Score is based
                  // on distance from projected tpos, and we want to minimise this.
		  
		  // Create new cluster and store it. 
		  // ***Last cluster in ClusterList will be that found by the fit.***
		  //          InoCluster* ClustNew = new InoCluster(Hit);
		  ClusterList[Plane+kl].push_back(ClustTemp);
		  
		  ShwOrTrkList[Plane+kl]=2;
                }
              }
            }
          }
        }  // End loop over planes
      } // End loop over views
      
      //      cout <<"Now, Perform LOCAL LINEAR FITS"<<endl;
      
      // Now, Perform LOCAL LINEAR FITS
      for(View=0; View<1; ++View) {
        for(Plane=begplane2; Plane<maxplane; ++Plane) {
	  
          // If plane has been declared track-like
	  if(ViewList[Plane]==View && ShwOrTrkList[Plane]>1) {
	    //	    cout <<"InoTrackFinder : *** FIT (2) ***   Plane=" << Plane << endl;
	    
            // Get last cluster from ClusterList. This will have been flagged immediately as track-like,
            // or have been stored as a result of a linear fit through a shower-like region.
            // Using the cluster, calculate a projected strip number and window for use in the local fit.
            InoCluster* clr = ClusterList[Plane].back();
	    
            ProjectedXPos=0.5*(clr->GetBegXPos()+clr->GetEndXPos()); 
            XPosWindow=0.6*StripXWidth;
	    
	    ProjectedYPos=0.5*(clr->GetBegYPos()+clr->GetEndYPos()); 
            YPosWindow=0.6*StripYWidth;

	    //GMA 09/02/2009 Fill this with caution, here is just take all of them
	    

	    // 07/02/2009 If we want hit infor kep it, otherwise emove the following part
            // Loop over the hits in this cluster

	    const unsigned int nhits = clr->GetHitEntries();
	    //   if(nhits<2)ClustList[Plane].push_back(clr);
	    
            if(nhits){
              for(unsigned int kl=0; kl<nhits; ++kl) {clr->GetHit(kl)->SetTrkFlag(1);}
              
              Xswx=0.0; Xswy=0.0; Xswx2=0.0; Xswxy=0.0; Xsw=0.0;
	      Yswx=0.0; Yswy=0.0; Yswx2=0.0; Yswxy=0.0; Ysw=0.0;

              // For these, loop over the hits in nearby clusters (within 2 planes of the current plane)
              // and carry out the linear fit
              for(int kl=-2; kl<3; ++kl){
		//		// ViewList[Plane+kl]==View &&
                if( ViewList[Plane+kl]==View && (Plane+kl)>=begplane2 && (Plane+kl)<(maxplane) && ShwOrTrkList[Plane+kl]>1 ) {
		  
                  InoCluster* clrtmp = ClusterList[Plane+kl].back(); // Again get last entry 
                  for(unsigned int lm=0; lm<clrtmp->GetHitEntries(); ++lm){
                    InoHit* hittmp = clrtmp->GetHit(lm);
		    
		    Zz=hittmp->GetZPos();
		    Xx=hittmp->GetXPos();
		    Xw=hittmp->GetXPulse()/clrtmp->GetXPulse();
		
		    Xswx+=Xw*Zz; Xswy+=Xw*Xx; Xswx2+=Xw*Zz*Zz; Xswxy+=Xw*Zz*Xx; Xsw+=Xw;

		    Yy=hittmp->GetYPos();
		    Yw=hittmp->GetYPulse()/clrtmp->GetYPulse();
		
		    Yswx+=Yw*Zz; Yswy+=Yw*Yy; Yswx2+=Yw*Zz*Zz; Yswxy+=Yw*Zz*Yy; Ysw+=Yw;
		    
		    //  x=hittmp->GetZPos(); 
		    //  y=hittmp->GetTPos(); 
		    //  w=hittmp->GetPulse()/clrtmp->GetPulse();
		    //  swx+=w*x; swy+=w*y; swx2+=w*x*x; swxy+=w*x*y; sw+=w;
		    
                  }
                }
              }
	      
              if(Xsw>1.0 && (Xsw*Xswx2-Xswx*Xswx)!=0.){
                Xm=(Xsw*Xswxy-Xswx*Xswy)/(Xsw*Xswx2-Xswx*Xswx);
                Xc=(Xswy*Xswx2-Xswx*Xswxy)/(Xsw*Xswx2-Xswx*Xswx);
                ProjectedXPos=Xm*clr->GetZPos()+Xc; 
                if(Xm>=0.0) XPosWindow = (0.6*StripXWidth)+0.02*Xm; 
                if(Xm<=0) XPosWindow = (0.6*StripXWidth)-0.02*Xm; 
              }
	      
              if(Ysw>1.0 && (Ysw*Yswx2-Yswx*Yswx)!=0.){
                Ym=(Ysw*Yswxy-Yswx*Yswy)/(Ysw*Yswx2-Yswx*Yswx);
                Yc=(Yswy*Yswx2-Yswx*Yswxy)/(Ysw*Yswx2-Yswx*Yswx);
                ProjectedYPos=Ym*clr->GetZPos()+Yc; 
                if(Ym>=0.0) YPosWindow = (0.6*StripYWidth)+0.02*Ym; 
                if(Ym<=0) YPosWindow = (0.6*StripYWidth)-0.02*Ym; 
              }      
            }//nhits
            
            // Using the results of the fit, obtain a score for each hit. Again, this is
            // based on distance from the projected strip number. We want to pick the hit which
            // minimises this distance.
            id=-1; dScore=0.; Score=999.;
            
            for(unsigned int kl=0; kl<nhits; ++kl){
	      InoHit* hit = clr->GetHit(kl);  
     	      //              dScore=hit->GetTPos()-ProjectedTPos; if(dScore<0) dScore=-dScore;
	      double xx = hit->GetXPos()-ProjectedXPos;
     	      double yy = hit->GetYPos()-ProjectedYPos;
     	      dScore = pow((xx*xx+yy*yy),0.5);
	      if(dScore<Score){id=kl; Score=dScore;}
            }
            
            // If we have found a good hit, store it and flag the plane as part of the 
            // final track
            
            if(1+id>0){ 
	      InoHit* hit = clr->GetHit(id); 
	      //              HitList[Plane].push_back(hit);
	      //asmS:Oct: here identify the point closest to the predicted point and see if the clusters X and Y position is close to the 
	      //       the selected hit, if yes set ShowOrTrkList[Plane] to 4 
	      
	      //            InoCluster* ClustNew1 = new InoCluster(hit);
	      //            ClustList[Plane].push_back(ClustNew1); 
              // Include any other hits which also match projection within the window
	      //  for(unsigned int kl=0; kl<nhits; ++kl){
	      //    InoHit* hittmp = clr->GetHit(kl);
	      if( ( (clr->GetXPos()-hit->GetXPos()>=0 &&
		     clr->GetXPos()-hit->GetXPos()<XPosWindow) ||
		      (clr->GetXPos()-hit->GetXPos()<=0 && 
		       clr->GetXPos()-hit->GetXPos()>-XPosWindow) ) &&
		    ( (clr->GetYPos()-hit->GetYPos()>=0 &&
		       clr->GetYPos()-hit->GetYPos()<YPosWindow) ||
		      (clr->GetYPos()-hit->GetYPos()<=0 && 
		       clr->GetYPos()-hit->GetYPos()>-YPosWindow) )  ) {
		//             HitList[Plane].push_back(hittmp);
		//             ClustNew1->AddHit(hittmp); 
		ClustList[Plane].push_back(clr);
		ShwOrTrkList[Plane]=4;
	      }
	      //             }
            }//if
          } // End demand that plane has been flagged as track-like
	  
	} // End loop over planes
      } // End loop over views
      
      // *** FINALLY, FORM THE TRACKS ***
      // Count the number of planes flagged as being part of final track
      TrackPlanes=0;
      
      for(Plane=begplane2; Plane<maxplane; ++Plane) {
	//        cout<<"ShwOrTrkList[Plane]"<<Plane<<" "<< ClustList[Plane].size()<<" " <<ShwOrTrkList[Plane]<<" : ";
        if(ShwOrTrkList[Plane]>3) {TrackPlanes++;}           //asm:3->2
      }
      
      if(TrackPlanes>2) { //GMA original if(TrackPlanes>5) {
        for(View=0; View<1; ++View) {
	  
          InoTrack* Trk = new InoTrack(); //slice);
          FinalTrackTempHolder[View].push_back(Trk);
          
          for(Plane=begplane2; Plane<maxplane; ++Plane) {
	    if(ViewList[Plane]==View && ShwOrTrkList[Plane]>3) {
	      // if(ShwOrTrkList[Plane]>2)   //asm:3->2 
	      for(unsigned int kl=0; kl<ClustList[Plane].size(); ++kl) {
		ClustList[Plane][kl]->SetTrkFlag(2);             
                Trk->AddCluster(ClustList[Plane][kl]); //VALGRIND
              }	                                           
            }
          }
        }
	if (FinalTrackTempHolder[0].size()>0) {
	}
	
        InoTrack* Trku = FinalTrackTempHolder[0].back();
        
        for(unsigned int kl=0; kl<Trku->GetEntries(); ++kl) {ClustsInTracks[0].push_back(Trku->GetCluster(kl));}
      }
    } // End U/V/TrackStrips check
    
    for(int jk=0; jk<500; ++jk) {
      ShwOrTrkList[jk]=-1; 
      ViewList[jk]=-1;
      ClusterList[jk].clear();
      ClustList[jk].clear();
      HitList[jk].clear();
    }
  } // End loop over temporary tracks   
  
  // Check tracks don't contain same hits as a previous track
  bool GoodTrack; bool ThreePlanes[2];
  int DuplicateStrips[2]; int Planes[2]; int ThisPlane;
  //  vector<InoHit*> HitsInPrevTracks[2];
  vector<InoCluster*> ClustsInPrevTracks[2];
  //GMA first look for the largest track segments
  //Put on 24th Feb 2008
  for (unsigned int ij=0; ij<FinalTrackTempHolder[0].size(); ++ij) {
    unsigned int jk = 1+FinalTrackTempHolder[0].size();
    unsigned int mxentries = 0;
    for (unsigned int kl=0; kl<FinalTrackTempHolder[0].size(); ++kl) {
      
      if (FinalTrackTempHolder[0][kl]->GetUsed() <0 &&
	  FinalTrackTempHolder[0][kl]->GetEntries() > mxentries) {
	jk = kl; mxentries = FinalTrackTempHolder[0][kl]->GetEntries();
      }
    }

    if (jk >=FinalTrackTempHolder[0].size()) continue;
    InoTrack* Trku = FinalTrackTempHolder[0][jk];
    
    //  GMA All three are swithed off now, but open it later on
    //  if(ModuleType==1) {LookForHitsAcrossGap(Trku); LookForHitsAcrossGap(Trkv);}
    
    ExtendTrack(Trku);  //GMA  //asm commented out 
    GoodTrack=false; 
    
    if(Trku->GetEntries()>2) { // && Trkv->GetEntries()>2) {
      DuplicateStrips[0]=0; DuplicateStrips[1]=0;
      
      // For u
      Planes[0]=-999; Planes[1]=-999; ThreePlanes[0]=false;
      
      for(unsigned int kl=0; kl<Trku->GetEntries(); ++kl) {
	InoCluster* TrackClust = Trku->GetCluster(kl);
	
	for(unsigned int lm=0; lm<ClustsInPrevTracks[0].size(); ++lm) {
	  if(TrackClust->GetShwPlnFlag()==true)continue;
	  if(TrackClust==ClustsInPrevTracks[0][lm]) {DuplicateStrips[0]++; break;}
	}
	
	ThisPlane=TrackClust->GetZPlane();
	if(Planes[0]<0) {Planes[0]=ThisPlane;}
	ThreePlanes[0]=true;
	//        else if(Planes[1]<0 && Planes[0]!=ThisPlane) {Planes[1]=ThisPlane;}
	//        else if(Planes[0]!=ThisPlane && Planes[1]!=ThisPlane) {ThreePlanes[0]=true;}
	
	if(DuplicateStrips[0]>2) {break;}
      }
      
      ///GMA how many common hits are allowed ? Optimise this 
	if(ThreePlanes[0] && DuplicateStrips[0]<3 ) {GoodTrack=true;}    
    }
    // asm: for reordering of tracks...   //Hit to Cluster 260709 asm
    
    if(GoodTrack) {
      
      vector<InoCluster*>::iterator it;
      InoTrack* OrdTrk = new InoTrack();
      OrdTrk->AddCluster(Trku->GetCluster(0));
      for (unsigned int kl=1;kl<Trku->GetEntries();kl++){
	int CZplane=Trku->GetCluster(kl)->GetZPlane();
	if(CZplane >= OrdTrk->GetCluster(OrdTrk->GetEntries()-1)->GetZPlane()) {
	  it=OrdTrk->end(); OrdTrk->InsertCluster(it,Trku->GetCluster(kl));
	} else if(CZplane <=OrdTrk->GetCluster(0)->GetZPlane() ) {
	  it=OrdTrk->begin();OrdTrk->InsertCluster(it,Trku->GetCluster(kl));
	} else if (CZplane>OrdTrk->GetCluster(0)->GetZPlane() && 
		   CZplane < OrdTrk->GetCluster(OrdTrk->GetEntries()-1)->GetZPlane()) {
	  it=OrdTrk->begin();
	  for (int lm=0;OrdTrk->GetEntries();lm++){
	    if(CZplane<OrdTrk->GetCluster(lm)->GetZPlane()) it++;
	    else break; 
	  }
	  OrdTrk->InsertCluster(it,Trku->GetCluster(kl)) ;
	}
      }
      
      FinalTrackBank[0].push_back(OrdTrk);
      //      FinalTrackBank[0].push_back(Trku);
      FinalTrackTempHolder[0][jk]->SetUsed(1);
      
      for(unsigned int kl=0; kl<Trku->GetEntries(); ++kl) {
	//	cout <<"kl = "<<kl<<" "<< Trku->GetCluster(kl)->GetZPlane()<<endl;
	ClustsInPrevTracks[0].push_back(Trku->GetCluster(kl));
      }
    }
  }
  
  //  asm: for reordering of tracks...Over
  FinalTrackTempHolder[0].clear(); FinalTrackTempHolder[1].clear();
  ClustsInPrevTracks[0].clear(); ClustsInPrevTracks[1].clear();

  return;
}
//====================================================FormFinalTracks========================================================
       
void InoTrackFinder::JoinCurvedTrack() {            //asm:new function
  //  cout<<"InoTrackFinder::JoinCurvedTrack()"<<endl;
  
  CLHEP::HepMatrix MatA(4,4,0);
  CLHEP::HepMatrix MatB(4,1,0);
  CLHEP::HepMatrix MatP(4,1,0);
  CLHEP::HepMatrix MatP1(4,1,0);
  CLHEP::HepMatrix MatF(4,1,0);
  CLHEP::HepMatrix MatAI(4,4,0);
  /*
    HepMatrix A(4,4,0);
    HepMatrix B(4,1,0);
    HepMatrix P(4,1,0);
    HepMatrix P1(4,1,0);
    HepMatrix F(4,1,0);
    HepMatrix AI(4,4,0);
  */
  
  //InoTrack* CurvTrk1;
  //InoTrack* CurvTrk2;
  //bool JoinCurv=false;
  //int  Sign1, Sign2 ;
  //int iimax=0;
  InoTrack *CurvTrk;
  vector<InoCluster*>::iterator it;
  InoCluster* CurrentClust = 0;

  bool joinClust;
  int joinBeg;
  int strhit;
  int ierr;

  double minDiff; 
  //MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;        //asm 
  for(int View=0; View<1; ++View) {
    for(unsigned int ij=0; ij<FinalTrackBank[View].size(); ++ij) {
      CurvTrk = FinalTrackBank[View][ij];
      if(CurvTrk->GetEntries()<4) break;    //asm: make it 5 
      
      for (joinBeg=-1; joinBeg<2;joinBeg=joinBeg+2){
	joinClust=true;
	int rev=1;
	while(joinClust){ 
	  joinClust=false;
	  for(int jk=1; jk<5;jk++){MatF(jk,1)=0.;}  //GMA14 Dimension 4 or five ?
	  for(int jk=1; jk<5;jk++){
	    MatA(jk,4)=0.; MatA(jk,1)=0.; MatA(jk,2)=0.; MatA(jk,3)=0.;
	    MatB(jk,1)=0.;
	    MatP(jk,1)=0.;
	    MatAI(jk,4)=0.; MatAI(jk,1)=0.; MatAI(jk,2)=0.; MatAI(jk,3)=0.; 
	  } 
	  
          //Here should I give less weightage when the cluster is showerlike. 
	  //	  int kl=0;   //GMA14     
	  strhit= (joinBeg>0) ? CurvTrk->GetEntries()-1:0;
	  for (int jk=0; jk<4;jk++){
	    int kl=1;
	    MatA(jk+1,1)= -2*CurvTrk->GetCluster(strhit-joinBeg*jk*kl)->GetXPos();
	    MatA(jk+1,2)= -2*CurvTrk->GetCluster(strhit-joinBeg*jk*kl)->GetYPos();
	    MatA(jk+1,3)= -2*CurvTrk->GetCluster(strhit-joinBeg*jk*kl)->GetZPos();
	    MatA(jk+1,4)=  1. ;
	    MatB(jk+1,1)= pow(CurvTrk->GetCluster(strhit-joinBeg*jk*kl)->GetXPos(),2)+ 
	      pow(CurvTrk->GetCluster(strhit-joinBeg*jk*kl)->GetYPos(),2)+ 
	      pow(CurvTrk->GetCluster(strhit-joinBeg*jk*kl)->GetZPos(),2);
	    
	  } // for jk
          
	  double rTrk2;
	  for(int jk=0; jk<3;jk++){
	    MatAI=MatA.inverse(ierr);
	    MatP= MatAI*(MatF-MatB);
	    rTrk2 = pow(MatP[0][0],2)+pow(MatP[1][0],2)+pow(MatP[2][0],2)-MatP[3][0];   // will this term be negative
	    MatF= MatA*MatP+MatB; 
	  } //it
	  
          int PlaneView=0;
          int extPln=0;
          int nClust;  
          int mxplane;
	  
          do {
	    bool ConsiderClust;
	    int id=-1;
	    
	    //----------------------------------------------------------------------------------
	    
	    // These are the tpos values required to locate the Magent coil hole
	    //GMA  put proper value 
	    double xx, yy; double XPos, YPos;
	    bool rpcGap=false, MegModuleGap=false;
	    double xxmax=1.52;
	    double fmax =0.92;
	    
	    if(joinBeg>0) { xx = abs(CurvTrk->GetEndXDir());yy = abs(CurvTrk->GetEndYDir()); XPos= CurvTrk->GetEndXPos();YPos=CurvTrk->GetEndYPos();}
	    else          { xx = abs(CurvTrk->GetBegXDir());yy = abs(CurvTrk->GetBegYDir()); XPos= CurvTrk->GetBegXPos();YPos=CurvTrk->GetBegYPos();}
	    
	    float  XXfrac =0, YYfrac=0;
	    int    XXq=-100, YYq=-100;
	    //ASM //Using this condition first decide if the track is passing through module gap or detector gap.
	    //Set some flag and set the appropriate condition for both the falgs.
            
	    if((CurvTrk->GetCluster(strhit)->GetXPos()>(-8.25-min(xx,1.52)*LayerThickness) &&
		CurvTrk->GetCluster(strhit)->GetXPos()<(-7.95+min(xx,1.52)*LayerThickness))||
	       (CurvTrk->GetCluster(strhit)->GetXPos()>( 7.95-min(xx,1.52)*LayerThickness) &&
		CurvTrk->GetCluster(strhit)->GetXPos()<( 8.25+min(xx,1.52)*LayerThickness) )){
	      MegModuleGap=true;
	      //cout<<"JCMegModuleGap"<<endl;
            }else{
	      
	      if(fabs(CurvTrk->GetCluster(strhit)->GetXPos())<8.00) {
		XXfrac= (int(fabs(CurvTrk->GetCluster(strhit)->GetXPos()*1000))%2000)/1000.;
		if(fabs(XXfrac-1)>fmax-fabs((int(fabs(min(xx,xxmax)*LayerThickness*1000))%2000)/1000.-1)){
		  XXq=int(((CurvTrk->GetCluster(strhit)->GetXPos()*1000))/2000)+24;
		  rpcGap=true;
		}
	      }
	      if(fabs(CurvTrk->GetCluster(strhit)->GetXPos())>8.20){
		XXfrac= (int((CurvTrk->GetCluster(strhit)->GetXPos()*1000)-200)%2000)/1000.;
		if(fabs(XXfrac-1)>fmax-fabs((int(fabs(min(xx,xxmax)*LayerThickness*1000)-200)%2000)/1000.-1)){
		  if(CurvTrk->GetCluster(strhit)->GetXPos()>0)XXq=int(((CurvTrk->GetCluster(strhit)->GetXPos()*1000)+200)/2000)+24;
		  else if(CurvTrk->GetCluster(strhit)->GetXPos()<0)XXq=int(((CurvTrk->GetCluster(strhit)->GetXPos()*1000)-200)/2000)+24;
		  rpcGap=true;
		}
	      }
	      YYfrac= (int(fabs(CurvTrk->GetCluster(strhit)->GetYPos()*1000))%2000)/1000.;
	      if(fabs(YYfrac-1)>fmax-fabs((int(fabs(min(yy,xxmax)*LayerThickness*1000))%2000)/1000.-1)){
		YYq=int(((CurvTrk->GetCluster(strhit)->GetYPos()*1000))/2000)+8;
		rpcGap=true;
	      }
	    }//else
	    
	    if(MegModuleGap){
	      mxplane=40;minDiff=0.3;
	    } else if(rpcGap) {
	      mxplane=20;minDiff=0.3;
	    } else {
	      mxplane=3;minDiff=0.2;
	    }
	    //---------------------------------------------------------------------------------- 
	    if (rev==1) extPln++;
	    if (rev==-1)extPln--; 
	    int CurrentPlane = (CurvTrk->GetCluster(strhit)->GetZPlane()+joinBeg*extPln>0&&CurvTrk->GetCluster(strhit)->GetZPlane()+joinBeg*extPln<140)? CurvTrk->GetCluster(strhit)->GetZPlane()+joinBeg*extPln:-1;
	    if (CurrentPlane<0) break;     
	    
	    nClust=ClusterBank[CurrentPlane].size();                       //asm:8_7_09
            double ff;
            if(nClust<100){
	      // Loop over the hits on the plane
	      for( int jk=0; jk<nClust; ++jk) {
		int XX1q=-24; int YY1q=-24;
              // Don't add a hit that is already in another track in the slice
		ConsiderClust=true;
		
		if(ClusterBank[CurrentPlane][jk]->GetShwPlnFlag()){ConsiderClust=false;}
		for(unsigned int kl=0; kl<ClustsInTracks[PlaneView].size(); ++kl) {
		  //GMA 260709  AllHitBank
		  if(ClusterBank[CurrentPlane][jk]==ClustsInTracks[PlaneView][kl]) {ConsiderClust=false; break;}
		  if(rev==-1&&(fabs(ClusterBank[CurrentPlane][jk]->GetXPos()-ClustsInTracks[PlaneView][kl]->GetXPos())< 0.04 ||
			       fabs(ClusterBank[CurrentPlane][jk]->GetYPos()-ClustsInTracks[PlaneView][kl]->GetYPos())< 0.04 ||
			       fabs(ClusterBank[CurrentPlane][jk]->GetXPos()-ClustsInTracks[PlaneView][kl]->GetXPos())>10.0  ||
			       fabs(ClusterBank[CurrentPlane][jk]->GetYPos()-ClustsInTracks[PlaneView][kl]->GetYPos())>10.0  )) {
		    ConsiderClust=false;
		  }
		}
		
		if(MegModuleGap||rpcGap){
		  if(XXq!=-24){
		    if(fabs(ClusterBank[CurrentPlane][jk]->GetXPos())<8.00){
		      XX1q=int((ClusterBank[CurrentPlane][jk]->GetXPos()*1000))/2000+24;
		    }
		    if(fabs(ClusterBank[CurrentPlane][jk]->GetXPos())>8.20){
		      if(ClusterBank[CurrentPlane][jk]->GetXPos()>0)XX1q=int(((ClusterBank[CurrentPlane][jk]->GetXPos()*1000)+200)/2000)+24;
		      if(ClusterBank[CurrentPlane][jk]->GetXPos()<0)XX1q=int(((ClusterBank[CurrentPlane][jk]->GetXPos()*1000)-200)/2000)+24;
		    }
		  }
		  if(YYq!=-24){
		    YY1q=int((ClusterBank[CurrentPlane][jk]->GetYPos()*1000)/2000)+8;
		  }
		}
		if((MegModuleGap||rpcGap) && (abs(XXq-XX1q)>1 && abs(YYq-YY1q)>1)) { ConsiderClust=false;}
		if(ConsiderClust==false)continue;
		CurrentClust = ClusterBank[CurrentPlane][jk];
		if(CurvTrk->ContainsClust(CurrentClust)==true) continue;
		
		
		// Find the most likely hit to add
		double XXo= CurrentClust->GetXPos() - MatP[0][0]; 
		double YYo= CurrentClust->GetYPos() - MatP[1][0];
		double ZZo= CurrentClust->GetZPos() - MatP[2][0];
		ff = pow(XXo*XXo+YYo*YYo+ZZo*ZZo,0.5)- pow(rTrk2,0.5);
		
		// 260709 asm  
		//pAnalysis->RC->Fill(pow(fabs(rTrk2),0.5),ff);
		//pAnalysis->RC->Fill( extPln,ff);
		
		//asm:Oct09 Along with this condition also check if the cluster added is in the to the track is in the same direction as that of the track. 
		if(abs(ff)<minDiff){ 
		  id=jk;
		  minDiff=fabs(ff);//<minDiff?fabs(ff):minDiff;
		  joinClust=true;
		} 
		//else{ cout<< abs(ff) <<" < "<< minDiff << endl; } 
	      } //nClust
	    }
	    if(joinClust){ 
	      CurrentClust=NULL;
	      CurrentClust = ClusterBank[CurrentPlane][id];
	      
              it= joinBeg<0 ? CurvTrk->begin():CurvTrk->end();
	      
	      CurvTrk->InsertCluster(it,CurrentClust);
              //CurvTrk->AddCluster(it);
              //FinalTrackBank[View].erase(FinalTrackBank[View].begin()+l);            //asm
              //FinalTrackBank[View].insert(FinalTrackBank[View].begin()+l,CurvTrk);   //asm
              //CurrentClust->SetStraight(false); //11Nov2009
	      ClustsInTracks[PlaneView].push_back(CurrentClust);
              if(pAnalysis->isXtermOut==1){
		cout<<"   Track extended 1" << rev << " "  <<CurrentClust->GetXPos()<<" "
		    <<CurrentClust->GetYPos()<<" "
		    <<CurrentClust->GetZPlane() << "  extpln" <<mxplane<<" "<< extPln<< endl;
              }
	      //                cout<<  ClustsInTracks[PlaneView][ClustsInTracks[PlaneView].size()-1]->GetZPlane();
	      //  for (unsigned int jk=0;jk<CurvTrk->GetEntries();jk++){cout << "  " <<CurvTrk->GetCluster(jk)->GetZPlane();}
	    }   
	    if (extPln==mxplane) {extPln=1;rev=-1;}
	  } while(extPln<mxplane &&extPln>-1&&joinClust==false);
	} //joinClust 
      } //for ( joinBeg
    } //for l
  }  //for View
  
  int View=0;
  if(FinalTrackBank[View].size()>1){
    unsigned int TrkMax;
    unsigned int iid=0;                                           
    TrkMax=FinalTrackBank[View][0]->GetEntries();  
    for(unsigned int ij=0; ij<FinalTrackBank[View].size(); ++ij){
      if(FinalTrackBank[View][ij]->GetEntries()>TrkMax){
	TrkMax=FinalTrackBank[View][ij]->GetEntries();iid=ij;}
    }
    if(iid>0){
      FinalTrackBank[View].insert(FinalTrackBank[View].begin(),FinalTrackBank[View][iid]);
      FinalTrackBank[View].erase(FinalTrackBank[View].begin()+iid+1); 
    }
    
    //   for(unsigned int ij=0; ij<FinalTrackBank[View].size(); ++ij){cout<<" "<<FinalTrackBank[View][ij]->GetEntries()<< "  " ;}    
  }

  return;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//asm: ExtendedTrack
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void InoTrackFinder::ExtendTrack(InoTrack* Trk) {          //asm: modified
  
  //  cout<<"InoTrackFinder : Attempting to extend track at beginning and end " << endl;
  //MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
  int PlanesSinceAdded, CurrentPlane(0), Increment(0);
  unsigned int nClust;
  double TrkXDir(0.), TrkYDir(0.); 
  double TrkXPos(0.), TrkYPos(0.), TrkZPos(0.),  Xextrem1, Xextrem2, Yextrem1, Yextrem2;
  double StripXCentre, StripYCentre, CurrentZPos, MinDistanceToHit[2], StripXCharge, StripYCharge;
  double PredictedXPos, PredictedYPos;
  
  double XTolerance; double YTolerance; int MaxPlanes;
  int PlaneView=0; // Trk->GetPlaneView();
  bool ConsiderClust;
  
  if(ModuleType==1 || ModuleType==3) {
    MaxPlanes=2; XTolerance=10.*StripXWidth; YTolerance=10.*StripYWidth;
  } else {
    MaxPlanes=40; XTolerance=4.*StripXWidth; YTolerance=4.*StripYWidth;
  }
  
  // Extend beginning first, then the end.
  for(int direction=0; direction<2; ++direction) {
    PlanesSinceAdded=0;            //asm:
    if(direction==0) {CurrentPlane=Trk->GetBegZPlane()-1; Increment=-1;}
    else if(direction==1) {CurrentPlane=Trk->GetEndZPlane()+1; Increment=1;}
    
    while(PlanesSinceAdded<MaxPlanes && CurrentPlane>=0 && CurrentPlane<=PlanesInModule) { 
      // Don't do anything if there are no hits, or we aren't in the same view.
      nClust=ClusterBank[CurrentPlane].size();
      if(nClust>0) { // && ClusterBank[CurrentPlane][0]->GetPlaneView()==PlaneView) {
	
        if(direction==0) {
          TrkXDir=Trk->GetBegXDir(); 
	  TrkYDir=Trk->GetBegYDir(); 
          TrkXPos=Trk->GetBegXPos(); 
	  TrkYPos=Trk->GetBegYPos(); 
          TrkZPos=Trk->GetBegZPos(); 
        } else if(direction==1) {
          TrkXDir=Trk->GetEndXDir(); 
	  TrkYDir=Trk->GetEndYDir(); 
          TrkXPos=Trk->GetEndXPos(); 
	  TrkYPos=Trk->GetEndYPos(); 
          TrkZPos=Trk->GetEndZPos(); 
        }
	//==========================================================================================
	//asm: defined Maxplanes and XTolerance as a function of slope.  //Maxplane can be reduced.. 
	if(fabs(TrkXDir)>2||fabs(TrkYDir)>2) {
	  MaxPlanes=8;
	  if(fabs(TrkXDir)>3) {XTolerance=20.*StripXWidth;} else {XTolerance=10.*StripXWidth;}
	  if(fabs(TrkYDir)>3) {YTolerance=20.*StripYWidth;} else {YTolerance=10.*StripYWidth;}
	} else {
	  MaxPlanes=8; XTolerance=3.*StripXWidth; YTolerance=3.*StripYWidth;
	}
	
	//=========================================================================================
        InoCluster* CurrentClust = 0;
        CurrentZPos = ClusterBank[CurrentPlane][0]->GetZPos();
        MinDistanceToHit[0] = XTolerance; 
	MinDistanceToHit[1] = YTolerance;
        
        // Loop over the hits on the plane
        for(unsigned int ij=0; ij<nClust; ++ij) {
	  
          // Don't add a hit that is already in another track in the slice
          ConsiderClust=true;
          //unsigned int  temp;
	  CurrentClust = 0;
	  
	  if(ClusterBank[CurrentPlane][ij]->GetShwPlnFlag()){ConsiderClust=false;}
          for(unsigned int jk=0; jk<ClustsInTracks[PlaneView].size(); ++jk) {
	    if(ClusterBank[CurrentPlane][ij]==ClustsInTracks[PlaneView][jk]) {
	      ConsiderClust=false;  break;
	    }
          }
          if(ConsiderClust==false) {PlanesSinceAdded++;continue;}
	  
	  //asm: redefined Xextrem1/2 and Yextrem1/2 now  taking into account the planes since the last hit was added. 
	  
          // Find the most likely hit to add
          PredictedXPos = TrkXPos + TrkXDir*(CurrentZPos-TrkZPos);
          Xextrem1 = PredictedXPos + ((PlanesSinceAdded+1)*LayerThickness*TrkXDir)/2; //GMA 0.055     
          Xextrem2 = PredictedXPos - ((PlanesSinceAdded+1)*LayerThickness*TrkXDir)/2;
          
          StripXCentre=ClusterBank[CurrentPlane][ij]->GetXPos();
          StripXCharge=ClusterBank[CurrentPlane][ij]->GetXPulse();
	  
          PredictedYPos = TrkYPos + TrkYDir*(CurrentZPos-TrkZPos);
          Yextrem1 = PredictedYPos + ((PlanesSinceAdded+1)*LayerThickness*TrkYDir)/2;
          Yextrem2 = PredictedYPos - ((PlanesSinceAdded+1)*LayerThickness*TrkYDir)/2;
          
	  StripYCentre=ClusterBank[CurrentPlane][ij]->GetYPos();
	  StripYCharge=ClusterBank[CurrentPlane][ij]->GetYPulse();
	  
          if(fabs(StripXCentre-PredictedXPos)<MinDistanceToHit[0] &&
	     fabs(StripYCentre-PredictedYPos)<MinDistanceToHit[1]) {
            MinDistanceToHit[0] = fabs(StripXCentre-PredictedXPos);
	    MinDistanceToHit[1] = fabs(StripYCentre-PredictedYPos);
            CurrentClust = ClusterBank[CurrentPlane][ij];
          } else if(fabs(StripXCentre-Xextrem1)<MinDistanceToHit[0] &&
		    fabs(StripYCentre-Yextrem1)<MinDistanceToHit[1]) {
            MinDistanceToHit[0] = fabs(StripXCentre-Xextrem1);
	    MinDistanceToHit[1] = fabs(StripYCentre-Yextrem1);
            CurrentClust = ClusterBank[CurrentPlane][ij];
          } else if(fabs(StripXCentre-Xextrem2)<MinDistanceToHit[0] &&
		    fabs(StripYCentre-Yextrem2)<MinDistanceToHit[1]) {
            MinDistanceToHit[0] = fabs(StripXCentre-Xextrem2);
	    MinDistanceToHit[1] = fabs(StripYCentre-Yextrem2);
            CurrentClust = ClusterBank[CurrentPlane][ij];
          }
        } //for nClust
        // Add best hit we have found
        if(CurrentClust) {
          Trk->AddCluster(CurrentClust); //VALGRIND
          ClustsInTracks[PlaneView].push_back(CurrentClust);
          PlanesSinceAdded=0;
	  //        cout<<" InoTrackFinder : Track extended " << direction << " TPos " 
	  //	      << CurrentClust->GetXPos()<<" "<<CurrentClust->GetYPos() 
	  //	      << " ZPos " << CurrentClust->GetZPos() << endl;
        }
      } else { // if(nClust>0)
	PlanesSinceAdded++;
      }
      
      CurrentPlane+=Increment;
    }  //while
  }  //for plane
  return;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Extended Track
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


void InoTrackFinder::FillGapsInTrack(InoTrack* Trk) {
  //a  cout<<"InoTrackFinder : Attempting to fill gaps in track " << endl;
  
  // First store track hits in plane order
  vector<InoCluster*> TrackClusts[490];
  int PlaneView=0; //Trk->GetPlaneView();
  
  for(unsigned int ij=0; ij<Trk->GetEntries(); ++ij) {
    InoCluster* Clust = Trk->GetCluster(ij);
    TrackClusts[Clust->GetZPlane()].push_back(Clust);
  }
  
  // Find gaps
  int PrevPlaneWithClust=Trk->GetBegZPlane();
  int NextPlaneWithClust=Trk->GetEndZPlane();

  double PrevXPos, NextXPos, PrevYPos, NextYPos, PrevZPos, NextZPos, CurrentZPos;
  
  double PredictedGradient[2], PredictedXPos, MinDistanceToClust[2], StripXCharge;
  double PredictedYPos, StripYCharge;
  
  double XTolerance; double YTolerance; bool ConsiderClust;

  if(ModuleType==1 || ModuleType==3) {XTolerance=1.5*StripXWidth; YTolerance=1.5*StripYWidth;}
  else {XTolerance=3.*StripXWidth; YTolerance=3.*StripYWidth;}

  for(int ij=Trk->GetBegZPlane()+1; ij<Trk->GetEndZPlane(); ij++) {
    // Don't do anything if there are no clusters, or we aren't in the same view.
    if(ClusterBank[ij].size()>0) { // (&& ClusterBank[ij][0]->GetPlaneView()==PlaneView)
      
      // Find the cluster planes on either side of a gap.
      if(TrackClusts[ij].size()>0) {PrevPlaneWithClust=ij; continue;}
      
      for(int jk=ij+2; jk<=Trk->GetEndZPlane(); jk++) {
        if(TrackClusts[jk].size()>0) {NextPlaneWithClust=jk; break;}
      }
      
      // Now try and fill the gap. Limit to maximum gap size?
      if(PrevPlaneWithClust!=NextPlaneWithClust) {
        CurrentZPos=ClusterBank[ij][0]->GetZPos();
	
        PrevZPos=TrackClusts[PrevPlaneWithClust][0]->GetZPos();
        NextZPos=TrackClusts[NextPlaneWithClust][0]->GetZPos();

        PrevXPos=TrackClusts[PrevPlaneWithClust][0]->GetXPos();
        NextXPos=TrackClusts[NextPlaneWithClust][0]->GetXPos();

        PrevYPos=TrackClusts[PrevPlaneWithClust][0]->GetYPos();
        NextYPos=TrackClusts[NextPlaneWithClust][0]->GetYPos();

        PredictedGradient[0]=(NextXPos-PrevXPos)/(NextZPos-PrevZPos);
	PredictedGradient[1]=(NextYPos-PrevYPos)/(NextZPos-PrevZPos);
	
        PredictedXPos=PrevXPos+((CurrentZPos-PrevZPos)*PredictedGradient[0]);
	PredictedYPos=PrevYPos+((CurrentZPos-PrevZPos)*PredictedGradient[1]);

	InoCluster* CurrentClust = 0;

	MinDistanceToClust[0] = XTolerance; 
	MinDistanceToClust[1] = YTolerance;
	
        for(unsigned int kl=0; kl<ClusterBank[ij].size(); ++kl) {
          InoCluster* clust = ClusterBank[ij][kl];
	  
	  
          // Don't add a cluster that is already in another track in the slice
          ConsiderClust=true;
          for(unsigned int lm=0; lm<ClustsInTracks[PlaneView].size(); ++lm) {
            if(clust==ClustsInTracks[PlaneView][lm]) {ConsiderClust=false; break;}
          }
          if(ConsiderClust==false) {continue;}
          
	  
          // Find the most likely cluster to add
          StripXCharge=clust->GetXPulse();
	  StripYCharge=clust->GetYPulse();

          if(fabs(PredictedXPos-clust->GetXPos())<MinDistanceToClust[0] &&
	     fabs(PredictedYPos-clust->GetYPos())<MinDistanceToClust[1]) {
            MinDistanceToClust[0]=fabs(PredictedXPos-clust->GetXPos());
	    MinDistanceToClust[1]=fabs(PredictedYPos-clust->GetYPos());
            CurrentClust=clust; 
          }
        }
	
        if(CurrentClust) {
          Trk->AddCluster(CurrentClust); 
          TrackClusts[ij].push_back(CurrentClust);
          ClustsInTracks[PlaneView].push_back(CurrentClust);
        }
      }
    }
  }
  
  // Clean up
  for(int ij=0; ij<490; ++ij) {TrackClusts[ij].clear();}
  
  return;
}

void InoTrackFinder::LookForHitsAcrossGap(InoTrack* Trk) {
  // Try and add missing clusters that are just across the SM gap
  
  int PlanesSinceAdded, CurrentPlane(0), ClustsAdded(0);
  double CurrentZPos, TargetXPos, TargetYPos;
  double TrkBegZPos, TrkBegXPos, TrkBegYPos;
  double TrkEndZPos, TrkEndXPos, TrkEndYPos;

  double TrkBegDir[2]={0,0};
  double TrkEndDir[2]={0,0};
  
  double Tolerance[2]={0,0};
  
  double MinDistanceToClust[2], DistanceToCurrentClust[2];
  double Xcentre, Xextrem1, Xextrem2;
  double Ycentre, Yextrem1, Yextrem2;
  unsigned int nclusts; 

  double z, t, w, sw, swx, swx2, swy, swyx;
  bool ConsiderClust;
  
  int PlaneView= 0; // Trk->GetPlaneView();

  // If we want to look back into SM 1
  //GMA an arbitray optin, depending on detector configuration one need to change it
  if(Trk->GetBegZPlane()>1140 && Trk->GetBegZPlane()<=1156) {  //GMA14 fixenumber
    
    if(PlaneView==0) {CurrentPlane=248;} else if(PlaneView==1) {CurrentPlane=247;}
    
    PlanesSinceAdded=0;
    TrkBegZPos=Trk->GetBegZPos();
    TrkBegXPos=Trk->GetBegXPos();
    TrkBegYPos=Trk->GetBegYPos();
    
    // Get beginning direction
    z=0.; t=0.; sw=0.; swx=0.; swx2=0.; swy=0.; swyx=0.;
    
    // Include first few clusters in track
    for(unsigned int ij=0; ij<Trk->GetEntries(); ++ij) {
      InoCluster* clust = Trk->GetCluster(ij);
      
      if(clust->GetZPlane()<Trk->GetBegZPlane()+4) {
        z=clust->GetZPos(); t=clust->GetXPos(); w=clust->GetXPulse();
        sw+=w; swx+=w*z; swx2+=w*z*z; swy+=w*t; swyx+=w*t*z;  
      }
    }
    // Use result so far to define a tolerance
    if((swx*swx-sw*swx2)!=0) {
      TrkBegDir[0] = (swx*swy-sw*swyx)/(swx*swx-sw*swx2); 
      Tolerance[0] = 0.2 + fabs(TrkBegDir[0]);
    }
    
    z=0.; t=0.; sw=0.; swx=0.; swx2=0.; swy=0.; swyx=0.;
    
    // Include first few clusters in track
    for(unsigned int ij=0; ij<Trk->GetEntries(); ++ij) {
      InoCluster* clust = Trk->GetCluster(ij);
      
      if(clust->GetZPlane()<Trk->GetBegZPlane()+4) {
        z=clust->GetZPos(); t=clust->GetYPos(); w=clust->GetYPulse();
        sw+=w; swx+=w*z; swx2+=w*z*z; swy+=w*t; swyx+=w*t*z;  
      }
    }
    // Use result so far to define a tolerance
    if((swx*swx-sw*swx2)!=0) {
      TrkBegDir[1] = (swx*swy-sw*swyx)/(swx*swx-sw*swx2); 
      Tolerance[1] = 0.2 + fabs(TrkBegDir[1]);
    }
    
    // Look at some of the relevant clusters across the gap
    for(int ij=CurrentPlane; ij>CurrentPlane-4; ij-=1) {
      for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
        InoCluster* clust = ClusterBank[ij][jk];
	
        z=clust->GetZPos(); t=clust->GetXPos(); w=clust->GetXPulse();
	
        // If cluster is within tolerance region, include with weight 0.5
        if(fabs(t-(TrkBegXPos+TrkBegDir[0]*(z-TrkBegZPos)))<Tolerance[0]) {
          sw+=w; swx+=w*z; swx2+=w*z*z; swy+=w*t; swyx+=w*t*z;  
        }
      }
    }
    // Get results of total fit
    if((swx*swx-sw*swx2)!=0) {TrkBegDir[0] = (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
    
    // Look at some of the relevant clusters across the gap
    for(int ij=CurrentPlane; ij>CurrentPlane-4; ij-=1) {
      for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
        InoCluster* clust = ClusterBank[ij][jk];
	
        z=clust->GetZPos(); t=clust->GetYPos(); w=clust->GetYPulse();
	
        // If cluster is within tolerance region, include with weight 0.5
        if(fabs(t-(TrkBegYPos+TrkBegDir[1]*(z-TrkBegZPos)))<Tolerance[1]) {
          sw+=w; swx+=w*z; swx2+=w*z*z; swy+=w*t; swyx+=w*t*z;  
        }
      }
    }
    // Get results of total fit
    if((swx*swx-sw*swx2)!=0) {TrkBegDir[1] = (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
    
    // Look to add the clusters
    while(CurrentPlane>0 && PlanesSinceAdded<5 && ClustsAdded<2) {
      
      nclusts=ClusterBank[CurrentPlane].size();
      
      if(nclusts>0) {
        CurrentZPos=ClusterBank[CurrentPlane][0]->GetZPos();
	
        InoCluster* BestClust = 0;
        MinDistanceToClust[0] = 3.*StripXWidth; 
	MinDistanceToClust[1] = 3.*StripYWidth;

        for(unsigned int jk=0; jk<nclusts; ++jk) {
          // Don't add a cluster that is already in another track in the slice
          ConsiderClust=true;
          for(unsigned int kl=0; kl<ClustsInTracks[PlaneView].size(); ++kl) {
            if(ClusterBank[CurrentPlane][jk]==ClustsInTracks[PlaneView][kl]) {
	      ConsiderClust=false; break;
	    }
          }
          if(ConsiderClust==false) {continue;}
	  
	  
          // Compare transverse positions
          TargetXPos=TrkBegXPos+TrkBegDir[0]*(CurrentZPos-TrkBegZPos);
	  TargetYPos=TrkBegYPos+TrkBegDir[1]*(CurrentZPos-TrkBegZPos);

          Xcentre=ClusterBank[CurrentPlane][jk]->GetXPos();
          Xextrem1=Xcentre+(LayerThickness*TrkBegDir[0]); //GMA 0.0852 Exact value from final detector geometry
          Xextrem2=Xcentre-(LayerThickness*TrkBegDir[0]);
          
          Ycentre=ClusterBank[CurrentPlane][jk]->GetYPos();
          Yextrem1=Ycentre+(LayerThickness*TrkBegDir[1]); //GMA 0.0852 Exact value from final detector geometry
          Yextrem2=Ycentre-(LayerThickness*TrkBegDir[1]);

          DistanceToCurrentClust[0]=min(fabs(Xcentre-TargetXPos),
                                   min(fabs(Xextrem1-TargetXPos),fabs(Xextrem2-TargetXPos)));
	  
	  DistanceToCurrentClust[1]=min(fabs(Ycentre-TargetYPos),
                                   min(fabs(Yextrem1-TargetYPos),fabs(Yextrem2-TargetYPos)));


          if(DistanceToCurrentClust[0]<MinDistanceToClust[0] &&
	     DistanceToCurrentClust[1]<MinDistanceToClust[1]) {
            MinDistanceToClust[0]=DistanceToCurrentClust[0]; 
	    MinDistanceToClust[1]=DistanceToCurrentClust[1]; 
            BestClust=ClusterBank[CurrentPlane][jk]; 
            PlanesSinceAdded=0;
          }
        }
        if(BestClust) {Trk->AddCluster(BestClust);}
      } else {PlanesSinceAdded+=1;}
      CurrentPlane-=1;
    }
  }
  
  // If we want to look forward into SM 2
  if(Trk->GetEndZPlane()>=242 && Trk->GetEndZPlane()<250) { //GMA Just an arbitray number
    
    if(PlaneView==0) {CurrentPlane=251;} else if(PlaneView==1) {CurrentPlane=250;}
    PlanesSinceAdded=0;
    TrkEndZPos=Trk->GetEndZPos();
    TrkEndXPos=Trk->GetEndXPos();
    TrkEndYPos=Trk->GetEndYPos();    
    
    // Get end direction
    z=0.; t=0.; sw=0.; swx=0.; swx2=0.; swy=0.; swyx=0.;
    
    // Include last few clusters in track
    for(unsigned int ij=0; ij<Trk->GetEntries(); ++ij) {
      InoCluster* clust = Trk->GetCluster(ij);
      
      if(clust->GetZPlane()>Trk->GetEndZPlane()-4) {
        z=clust->GetZPos(); t=clust->GetXPos(); w=clust->GetXPulse();
        sw+=w; swx+=w*z; swx2+=w*z*z; swy+=w*t; swyx+=w*t*z; 
      }
    }
    // Use result so far to define a tolerance
    if((swx*swx-sw*swx2)!=0) {TrkEndDir[0] = (swx*swy-sw*swyx)/(swx*swx-sw*swx2); Tolerance[0] = 0.2 + fabs(TrkEndDir[0]);}
    
    // Get end direction
    z=0.; t=0.; sw=0.; swx=0.; swx2=0.; swy=0.; swyx=0.;
    
    // Include last few clusts in track
    for(unsigned int ij=0; ij<Trk->GetEntries(); ++ij) {
      InoCluster* clust = Trk->GetCluster(ij);
      
      if(clust->GetZPlane()>Trk->GetEndZPlane()-4) {
        z=clust->GetZPos(); t=clust->GetYPos(); w=clust->GetYPulse();
        sw+=w; swx+=w*z; swx2+=w*z*z; swy+=w*t; swyx+=w*t*z; 
      }
    }
    // Use result so far to define a tolerance
    if((swx*swx-sw*swx2)!=0) {TrkEndDir[1] = (swx*swy-sw*swyx)/(swx*swx-sw*swx2); Tolerance[1] = 0.2 + fabs(TrkEndDir[1]);}
    
    // Look at some of the relevant clusters across the gap
    for(int ij=CurrentPlane; ij<CurrentPlane+4; ij+=1) {
      for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
        InoCluster* clust = ClusterBank[ij][jk];
        
        z=clust->GetZPos(); t=clust->GetXPos(); w=clust->GetPulse();
	
        // If clust is within tolerance region, include with weight 0.5
        if(fabs(t-(TrkEndXPos+TrkEndDir[0]*(z-TrkEndZPos)))<Tolerance[0]) {
          sw+=w; swx+=w*z; swx2+=w*z*z; swy+=w*t; swyx+=w*t*z;  
        }
      }
    }
    
    // Get results of total fit
    if((swx*swx-sw*swx2)!=0) {TrkEndDir[0] = (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
    
    // Look at some of the relevant clusts across the gap
    for(int ij=CurrentPlane; ij<CurrentPlane+4; ij+=1) {
      for(unsigned int jk=0; jk<ClusterBank[ij].size(); ++jk) {
        InoCluster* clust = ClusterBank[ij][jk];
        
        z=clust->GetZPos(); t=clust->GetYPos(); w=clust->GetPulse();
	
        // If clust is within tolerance region, include with weight 0.5
        if(fabs(t-(TrkEndYPos+TrkEndDir[1]*(z-TrkEndZPos)))<Tolerance[1]) {
          sw+=w; swx+=w*z; swx2+=w*z*z; swy+=w*t; swyx+=w*t*z;  
        }
      }
    }
    
    // Get results of total fit
    if((swx*swx-sw*swx2)!=0) {TrkEndDir[1] = (swx*swy-sw*swyx)/(swx*swx-sw*swx2);}
    
    // Look at clusters in region
    while(CurrentPlane<500 && PlanesSinceAdded<5 && ClustsAdded<2) {
      
      nclusts=ClusterBank[CurrentPlane].size();
        
      if(nclusts>0) {
        CurrentZPos=ClusterBank[CurrentPlane][0]->GetZPos();

        InoCluster* BestClust = 0;
        MinDistanceToClust[0] = 3.*StripXWidth; 
	MinDistanceToClust[1] = 3.*StripYWidth;
                  
        for(unsigned int jk=0; jk<nclusts; ++jk) {
          // Don't add a cluster that is already in another track in the slice
          ConsiderClust=true;
          for(unsigned int kl=0; kl<ClustsInTracks[PlaneView].size(); ++kl) {
            if(ClusterBank[CurrentPlane][jk]==ClustsInTracks[PlaneView][kl]) {
	      ConsiderClust=false; break;
	    }
          }
          if(ConsiderClust==false) {continue;}
	  
          // Compare transverse positions
          TargetXPos=TrkEndXPos+TrkEndDir[0]*(CurrentZPos-TrkEndZPos);
	  
          Xcentre=ClusterBank[CurrentPlane][jk]->GetXPos();
          Xextrem1=Xcentre+(LayerThickness*TrkEndDir[0]);
          Xextrem2=Xcentre-(LayerThickness*TrkEndDir[0]);
          
          DistanceToCurrentClust[0]=min(fabs(Xcentre-TargetXPos),
					min(fabs(Xextrem1-TargetXPos),fabs(Xextrem2-TargetXPos)));
	  
	  
          TargetYPos=TrkEndYPos+TrkEndDir[0]*(CurrentZPos-TrkEndZPos);
	  
          Ycentre=ClusterBank[CurrentPlane][jk]->GetYPos();
          Yextrem1=Ycentre+(LayerThickness*TrkEndDir[1]);
          Yextrem2=Ycentre-(LayerThickness*TrkEndDir[1]);
          
          DistanceToCurrentClust[1]=min(fabs(Ycentre-TargetYPos),
					min(fabs(Yextrem1-TargetYPos),fabs(Yextrem2-TargetYPos)));
	  
          if(DistanceToCurrentClust[0]<MinDistanceToClust[0] &&
	     DistanceToCurrentClust[1]<MinDistanceToClust[1]) {
            MinDistanceToClust[0]=DistanceToCurrentClust[0]; 
	    MinDistanceToClust[1]=DistanceToCurrentClust[1]; 
            BestClust=ClusterBank[CurrentPlane][jk]; 
            PlanesSinceAdded=0;
          }
        }
        
        if(BestClust) {Trk->AddCluster(BestClust);}
      } else {PlanesSinceAdded+=1;}
      CurrentPlane+=1;
    }
  }
  return;
}

void InoTrackFinder::ClearUp() {
  
  // Reset containers and tidy memory usage
  //  unsigned int i, j;
  
  for (unsigned int ij=0; ij<500; ++ij) {
    for(unsigned int jk=0; jk<HitBank[ij].size(); ++jk) {
      if(HitBank[ij][jk]) {delete HitBank[ij][jk]; HitBank[ij][jk]=0;}
    }
    HitBank[ij].clear();

    ClusterBank[ij].clear();
    for(unsigned int jk=0; jk<ClusterList[ij].size(); ++jk) {
      if(ClusterList[ij][jk]){ delete ClusterList[ij][jk];ClusterList[ij][jk]=0;}
    }
    ClusterList[ij].clear();

    for(unsigned int jk=0; jk<SegmentBank[ij].size(); ++jk) { 
      if(SegmentBank[ij][jk] ) { delete SegmentBank[ij][jk]; SegmentBank[ij][jk]=0;}
    }
    SegmentBank[ij].clear();

    for(unsigned int jk=0; jk<ClustList[ij].size(); ++jk) {
      if(ClustList[ij][jk]){ delete ClustList[ij][jk];ClustList[ij][jk]=0;}
    }
    ClustList[ij].clear();   
  }
  
  for (unsigned int ij=0; ij<2; ++ij) {
    TempTrack[ij].clear();
    PossibleJoins[ij].clear();
    FinalTrackBank[ij].clear();
    ViewSegBank[ij].clear();
    ClustsInTracks[ij].clear();
    for(unsigned int jk=0; jk<FinalTrackTempHolder[ij].size(); ++jk) {
      if(FinalTrackTempHolder[ij][jk]){ delete FinalTrackTempHolder[ij][jk];FinalTrackTempHolder[ij][jk]=0;}
    }
    FinalTrackTempHolder[ij].clear(); 
  }
}


void InoTrackFinder::FirstComparison() {
  // Having made sections of the track, from 'matched' segments, we compare
  // the sections in the U view with those in the V view, looking for compatibity.

  // First all clusters are investigated to mark the segments as track-like or
  // shower-like. Then track-like segments which have 'overlapping' partners in the 
  // opposite view are marked by setting their UID values to 2

  //  cout <<"InoTrackFinder : *** Comparing views (1) *** " << endl;


  int shwctr, plnctr;
  unsigned int nsegments=0;
  //  unsigned int nsegmentsV=0;
  //  vector<InoTrackSegment*> mysegbnk[2];


  // Determine track-like segments
  for(int View=0; View<1; ++View) {

    // Loop over all segments in the current view
    nsegments = ViewSegBank[View].size();
    for(unsigned int ij=0; ij<nsegments; ++ij) {
      
      InoTrackSegment* Seg = ViewSegBank[View][ij];
      shwctr=0; plnctr=0;
      
      // Loop over all the clusters in the current segment
      const unsigned int ncluster = Seg->GetEntries();
      for(unsigned int jk=0;jk<ncluster;++jk) {
        InoCluster* cluster = Seg->GetCluster(jk);
        
        // Find how track-like or shower-like the cluster is
        if(cluster->GetTrkFlag()>1) { 
          if(cluster->GetShwFlag()<1) {shwctr++;}
        }
        if(cluster->GetTrkPlnFlag()>0) {plnctr++;}
      }
      
      // Set the UID accordingly. This depends greatly on the initial values 
      // of the shower flags used when identifying track and shower clusters.
      //      if(CambridgeAnalysis==false) {Seg->SetUID(1);}
      
      Seg->SetUID(1);
      
// asm://     cout <<" set "<< plnctr<<" "<<shwctr <<" "<<Seg->GetEntries()<<" "<<Seg->GetUID()<<endl; 

      if( plnctr>1 && Seg->GetEntries()>4 ) { // Just added
        if(shwctr>2) {Seg->SetUID(2);} 

//a	cout <<" set2 "<< plnctr<<" "<<shwctr <<" "<<Seg->GetEntries()<<" "<<Seg->GetUID()<<endl; 	
      }
    }
  }
  
  return;
}


void InoTrackFinder::Trace(const char * /* c */) const
{
}

void InoTrackFinder::CleanAndFilled() {

  //MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
  
  //11Nov2009 Moved from  void InoTrackFinder::FormFinalTracks( ) 
  /* */  // Print out list of final tracks           
  //  cout<< "InoTrackFinder : *** LIST OF TRACKS *** " <<FinalTrackBank[0].size()<< endl;
  for(int View=0; View<1; ++View) {                       //asm
    for(unsigned int ij=0; ij<FinalTrackBank[View].size(); ++ij) { 
      InoTrack* TrkTemp = FinalTrackBank[View][ij];
      
      TrkTemp->SetStraight();
      inoTrack_pointer->InoTrack_list.push_back(TrkTemp);
      /*  Reopen 11Nov2009 */
      //const int nrandom = 10;
      //float randvar[nrandom];
      for (unsigned int jk=0; jk<TrkTemp->GetEntries(); ++jk){
	if(jk>0) {
          //asm: commented out while releasing the code.
	  //if(TrkTemp->GetCluster(jk)->GetZPlane()==TrkTemp->GetCluster(jk-1)->GetZPlane())  cout<<"Double Z planes++++++++++++ "<<endl;
	}
	
	//GMA 090809 why ??????
	if (jk >0 && 
	    TrkTemp->GetCluster(jk)->GetZPlane() == TrkTemp->GetCluster(jk-1)->GetZPlane()) {
	  TrkTemp->ClustsInTrack.erase(TrkTemp->ClustsInTrack.begin()+jk);
	  jk--;
	}
      }
     // cout <<"tmp 3 "<< TrkTemp->ClustsInTrack.size()<<" "<<TrkTemp->GetEntries()<<endl;
      if(pAnalysis->isXtermOut==1){
	cout<< "\n InoTrackFinder_final tracks"
	    << "  " << View << " " << ij
	    << "  " << TrkTemp->GetBegZPlane() << "  " << TrkTemp->GetEndZPlane() << endl;
	
      } //if(isXtermOut){
      
    }
  }
}
