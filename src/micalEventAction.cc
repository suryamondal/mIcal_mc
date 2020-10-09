//
// $Id: micalEventAction.cc,v 1.24 2005/05/30 14:24:31 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
//GMAA Store begin and end hit position and muo energy at those points for the resolution of reconstruction algorithms, which is much better than true muon energy resolution

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "micalEventAction.hh"
#include "micalEventActionMessenger.hh"

#include "micalCal0Hit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"

#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

#include "Randomize.hh"
#include <iomanip>
#include <utility>
#include "G4SDManager.hh"
#include "G4Trajectory.hh"
//#include <iomanip.h>

#include "vect_manager.h"
//#include "InoPatternRecognition.h" //GMA14
#include "InoTrackFinder.h"
#include "InoTrackFitAlg.h"
#include "InoVertex.h"

#include "TStyle.h"
#include "TMatrixD.h"
#include "TMath.h"
//#include "TVectorD.h"
//#include "TMatrixTBase.h"
#include "TMatrixDEigen.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalEventAction::micalEventAction()
  :drawFlag("all"),printModulo(1000),eventMessenger(0)
				  //  :drawFlag("charged"),printModulo(1),eventMessenger(0)
{
  eventMessenger = new micalEventActionMessenger(this);
  cal0CollID = -1;
  nevent = 0;
  rang = 0;
  //VtxPlane[]={0};
  //EndPlane[]={0};
  //  cal1CollID = -1;  //  cal2CollID = -1;  //  c0 =0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

micalEventAction::~micalEventAction() {
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Orighit algo has been included as  a function so that we can use more than once
//only by calling this function
//added by SSE 09/15
int micalEventAction::orighit_calc(vector<InoHit*> tmphitcluster){
  vector<iXYZ> xShwstrip;
  vector<iXYZ> yShwstrip;

  xShwstrip.clear();yShwstrip.clear(); 
  // cout<< " xShwstrip.size "<< xShwstrip.size()<< " yShwstrip.size "<< yShwstrip.size() <<endl;
  for( unsigned int mn=0; mn<tmphitcluster.size(); mn++) {
    x_stripno = tmphitcluster[mn]->GetXStripNum();
    y_stripno = tmphitcluster[mn]->GetYStripNum();
    z_plane=tmphitcluster[mn]->GetZPlane();

    if (x_stripno && y_stripno) {
      //cout<<" 09/15 orig func chk: "<< " "<<z_plane << " "<<x_stripno<< " "<<y_stripno<< " UID "<<tmphitcluster[mn]->GetUID()<<endl;
      for (unsigned nn=0; nn< xShwstrip.size(); nn++) {
	if (x_stripno == xShwstrip[nn].second &&
	    z_plane == xShwstrip[nn].first) {
	  // cout<< "09/15: ix working"<<endl;
	  x_stripno = -1; break;
	}
      }
      for (unsigned nn=0; nn< yShwstrip.size(); nn++) {
	if (y_stripno == yShwstrip[nn].second &&
	    z_plane == yShwstrip[nn].first) {
	  //cout<< "09/15: iy working"<<endl;
	  y_stripno = -1; break;
	}
      }
      if ((x_stripno>=0 || y_stripno>=0) ) {
	if(x_stripno>=0){
	  //cout<< " 09/15: ixstripno "<< ixstripno<<endl;
	  iXYZ ZXstrp(z_plane,x_stripno);xShwstrip.push_back(ZXstrp);
	}
	if(y_stripno>=0){
	  //cout<< " 09/15: iystripno "<< iystripno<<endl;
	  iXYZ ZYstrp(z_plane,y_stripno);yShwstrip.push_back(ZYstrp);
	}
      }//end of if ((ixstripno>=0 || iystripno>=0) )
    }
  }//end of for( int mn=0; mn<tcluster.size(); mn++) loop on total hits
  // cout<< " xShwstrip.size "<< xShwstrip.size()<< " yShwstrip.size "<< yShwstrip.size() <<endl;
  int x_nhits =0;
  int y_nhits =0;
  int Orighits_all=0;
  for(int kl=0; kl<nLayer; ++kl){

    if(xShwstrip.size() && yShwstrip.size()){
      x_nhits =0;
      y_nhits =0;
      for(unsigned nn=0; nn< xShwstrip.size(); nn++) {
	if(xShwstrip[nn].first==int(kl)){
	  x_nhits++;
	}
      }
      for (unsigned nn=0; nn< yShwstrip.size(); nn++){
	if(yShwstrip[nn].first==int(kl)){
	  y_nhits++;
	}
      }
      Orighits_all +=max(x_nhits,y_nhits);
    }
  }
  return Orighits_all;
}//SSE 09/15

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void micalEventAction::BeginOfEventAction(const G4Event* evt) {  
  G4int evtNb = evt->GetEventID();
  if (evtNb%printModulo == 0) {
    
    G4cout << "\n---> Begin of event: " << evtNb << G4endl;
    //HepRandom::showEngineStatus();
  }
  
  rang = 0.0;
  //initialisation per event
  Energycal0 = EnergyGap = 0.;
  TrackLcal0 = TrackLGap = 0.;
  //Energycal1 = Energycal2 = 0.;
  //TrackLcal1 = TrackLcal2 = 0.;
  
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  //  if(cal0CollID<0||cal1CollID<0||cal2CollID<0)
  
  if(cal0CollID<0) {
    G4String colNam;
    cal0CollID = SDman->GetCollectionID(colNam="cal0Collect");
    //    cal1CollID = SDman->GetCollectionID(colNam="cal1Collect");
    //    cal2CollID = SDman->GetCollectionID(colNam="cal2Collect");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void micalEventAction::EndOfEventAction(const G4Event* evt) {
  paradef = micalDetectorParameterDef::AnPointer;
  //StripXWidth = (1/m)*paradef->GetXStrwd();
  //StripYWidth = (1/m)*paradef->GetYStrwd();
  nLayer      = paradef->GetnLayer();
  LayerThickness = (1/m)*2*(paradef->GetParlay(2)+paradef->GetParirlay(2));
  
  typedef pair<int,int> ixyz;
  nevent++;
  MultiSimAnalysis *pAnalysis = MultiSimAnalysis::AnPointer;
  InoRPCStrip_Manager* inoRPC_pointer = new InoRPCStrip_Manager(); //InoRPCStrip_Manager::APointer;

  pAnalysis->range = rang; //meghna
  pAnalysis->nhtcal0 = 0;
  
  G4int evtNb = evt->GetEventID();
  cout<<"# "<<evtNb<<endl;
  // pAnalysis->timeAsciiOutput<<"# "<<evtNb<<endl;
  if (evtNb%printModulo == 0) {
    if(pAnalysis->isXtermOut==1) {
      //    G4cout << "---> End of event: " << evtNb << G4endl;
      
      G4cout
	<< G4endl
	<< "    Absrober: total energy: " << std::setw(7)
	<< G4BestUnit(EnergyGap,"Energy")
	<< "        total track length: " << std::setw(7)
	<< G4BestUnit(TrackLGap,"Length")
	<< G4endl
	<< "RPC_gas: cal0 total energy: " << std::setw(7)
	<< G4BestUnit(Energycal0,"Energy")
	<< "        total track length: " << std::setw(7)
	<< G4BestUnit(TrackLcal0,"Length")
	<< G4endl;
    }
    pAnalysis->simtotabenr = EnergyGap;
    pAnalysis->simtotrpcenr= Energycal0;
    pAnalysis->simtotablen = TrackLGap;
    pAnalysis->simtotrpclen= TrackLcal0;
  }
  
  if (pAnalysis->InputOutput <=2) {
    G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
    if (pVisManager) {
      G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
      G4int n_trajectories = 0;
      if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
      
      for (G4int ij=0; ij<n_trajectories; ij++) {
	G4VTrajectory* trj = ((*(evt->GetTrajectoryContainer()))[ij]);
	if (drawFlag == "all") { pVisManager->Draw(*trj); // GMA14 ,100);
	} else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.)) {
	  pVisManager->Draw(*trj); // GMA14 ,100);
	} else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.)) {
	  pVisManager->Draw(*trj); //GMA14 ,100);
	}
      }
      
      if(cal0CollID<0) return;
      
      G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
      micalcal0HitsCollection* EHC0 = 0;
      //  micalcal1HitsCollection* EHC1 = 0;
      //  micalcal2HitsCollection* EHC2 = 0;
      if(HCE) {
	EHC0 = (micalcal0HitsCollection*)(HCE->GetHC(cal0CollID));
	//    EHC1 = (micalcal1HitsCollection*)(HCE->GetHC(cal1CollID));
	//    EHC2 = (micalcal2HitsCollection*)(HCE->GetHC(cal2CollID));
      }
      if(EHC0) {
	int n_hit = EHC0->entries();
	//    pAnalysis->ascii_output <<n_hit<<" ";
	
	G4double totE = 0;
	for (int ij=0; ij<n_hit; ij++) {
	  totE +=(*EHC0)[ij]->GetEdep();
	  //G4cout <<"Energy deposited in hit "<<i<<" "<<(*EHC0)[ij]->GetEdep()<< " at "<<(*EHC0)[ij]->GetPos()<<" at "<<(*EHC0)[ij]->GetTime()<<" "<<(*EHC0)[ij]->GetHitId() <<G4endl;
	  if (pAnalysis->nhtcal0 <pAnalysis->nmxhit) {
	    pAnalysis->calid0[pAnalysis->nhtcal0] = (*EHC0)[ij]->GetHitId();
	    pAnalysis->calen0[pAnalysis->nhtcal0] = (*EHC0)[ij]->GetEdep()/keV;
	    //G4ThreeVector MCpos = (*EHC0)[i]->GetPos();
	    //cout<<"x "<<1.e-3*MCpos.x()<<"     y "<<1.e-3*MCpos.x()<<"     z "<<1.e-3*MCpos.z()<<endl;
	    pAnalysis->nhtcal0++;
	    //pAnalysis->ascii_output <<std::setw(12)<<(*EHC0)[i]->GetHitId()<<" "<<std::setw(9)<<(*EHC0)[i]->GetEdep()/keV;
	  }
	}
	
	pAnalysis->caltot0 = totE/keV;
	//    pAnalysis->ascii_output <<"  "<<totE/keV<<G4endl;
	if(pAnalysis->isXtermOut==1) {
	  G4cout << "     " << n_hit
		 << " hits are stored in micalEcalHitsCollection with total Energy  "<<totE<< " KeV"<<endl; // G4BestUnit(totE,"Energy") << G4endl;
	}
      }
    } // if (pVisManager)
  } // if (pAnalysis->InputOutput <=2)
  
  if (pAnalysis->InputOutput==0 || pAnalysis->InputOutput==3 || pAnalysis->InputOutput==5) {
    //    InoPatternRecognition pattern; //GMA14
    //    pattern.RunPattRecog();
    
    InoTrackFinder trackfinder;
    // cout<<"Running the finder now ...."<<endl;
    trackfinder.RunTheFinder();	//VALGRIND
    // cout<<"Completed..."<<endl;
    //Now add other cluster close to vertex;
    vector <InoCluster*> totcluster = trackfinder.inoCluster_pointer->InoCluster_list;
    int totclustersize = totcluster.size();
    // cout<<"totalclustersize = "<<totclustersize<<endl;
    InoTrack_Manager *pinotrack = InoTrack_Manager::APointer;
    // pinotrack = 0;
    if (pinotrack) {
      // cout<<"pAnalysis->InoTrack_listsize->Fill(pinotrack->InoTrack_list.size()); "<<pinotrack->InoTrack_list.size()<<endl;
      pAnalysis->InoTrack_listsize->Fill(pinotrack->InoTrack_list.size());        //asm : from here
      // cout<<"Hola..."<<endl;
      //default version has isVisOut= 1, but we change that to 0 to obtain hit information in output without generating .inh file
      if (pAnalysis->isVisOut==1)  {
	pAnalysis->H->NFinders=pinotrack->InoTrack_list.size(); // number of finder sets
	for (unsigned ij=0; ij<pinotrack->InoTrack_list.size() ; ij++) {
	  for ( unsigned int jk =0; jk<pinotrack->InoTrack_list[ij]->ClustsInTrack.size();jk++) {
	    pAnalysis->Hp= pAnalysis->H->AddHits(0,0); // add a track object //VALGRIND
	    pAnalysis->Hp->TrackType=-4;// Track Type: -1: hits, -2: clulster, -3: triplet, -4: track
	    pAnalysis->Hp->FindNum=ij;// track Number
	    pAnalysis->Hp->ZZ=pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetZPlane();
	    pAnalysis->Hp->XX=pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetXPos();
	    pAnalysis->Hp->YY=pinotrack->InoTrack_list[ij]->ClustsInTrack[jk]->GetYPos();
	  }
	} //asm : upto here
	
	if (pinotrack->InoTrack_list.size()>0) {
	  for (unsigned int jk =0; jk<pinotrack->InoTrack_list[0]->ClustsInTrack.size();jk++) {
	    //GMA where is it ptinting ?
	    //cout<<"pAnalysis->TrkDist->Fill(pinotrack->InoTrack_list[0]->ClustsInTrack[jk]->GetZPlane());"<<endl;
	    pAnalysis->TrkDist->Fill(pinotrack->InoTrack_list[0]->ClustsInTrack[jk]->GetZPlane());
	    //cout<<"pAnalysis->EffDist->Fill(pinotrack->InoTrack_list[0]->ClustsInTrack[jk]->GetZPlane());"<<endl;
	    pAnalysis->EffDist->Fill(pinotrack->InoTrack_list[0]->ClustsInTrack[jk]->GetZPlane());
	  }
	}
      }
      // cout<<"Going to Track Fitter..."<<endl;
      InoTrackFitAlg trackfitter;
      // cout<<"Running Track Fit Alg ..."<<endl;
      trackfitter.RunAlg(); //VALGRIND
      // cout<<"Track Fit Alg Completed "<<endl;
      vector<ixyz> xtrkstrip;
      vector<ixyz> ytrkstrip;
      vector<ixyz> xshwstrip;
      vector<ixyz> yshwstrip;
      
      int tottrkXstrp = 0;
      int tottrkYstrp = 0;

      InoTrackCand_Manager *pfitTrack = InoTrackCand_Manager::APointer;
      
      if (pfitTrack) {
	//	G4cout <<"tmphitlist in eventaction "<< pinohit->InoHit_list.size()<<G4endl;
	// cout <<"tmptracklist in eventaction "<< nevent<<" "<<pfitTrack->InoTrackCand_list.size()<<" "<<pAnalysis->ihist<<G4endl;
	{
	  unsigned ij=0;
	
	  // cout<<"pAnalysis->ntrkmx "<<pfitTrack->InoTrackCand_list.size()<<endl;
	  for (unsigned jk=0; jk<pfitTrack->InoTrackCand_list.size() ; jk++) {
	    // cout <<"Identical "<< int(pfitTrack->InoTrackCand_list[jk]->ClustsInTrack[0]->isIdentical(trackfinder.inoCluster_pointer->InoCluster_list[0]));
	    G4HCofThisEvent * hce = evt->GetHCofThisEvent();
	    //	    micalcal0HitsCollection* ehco = 0;
	    G4ThreeVector MCpos;
	    if(hce) {
	      // SSE : Error while running code sepataterly (other than 0 option)
	      //  ehco = (micalcal0HitsCollection*)(hce->GetHC(cal0CollID));
	      //MCpos = (*ehco)[jk]->GetPos();
	    }
	  
	    if (ij <pAnalysis->ntrkmx) {
	      pAnalysis->itype[ij] =  pfitTrack->InoTrackCand_list[jk]->GetFitType();
	      pAnalysis->nhits[ij] = (pfitTrack->InoTrackCand_list[jk]->GetNDOF()+5)/2;
	      pAnalysis->chisq[ij] =  pfitTrack->InoTrackCand_list[jk]->GetChi2();
	      pAnalysis->cvalue[ij] = pfitTrack->InoTrackCand_list[jk]->Getcval();
	    
	      pAnalysis->fc_or_pc[ij] = pfitTrack->InoTrackCand_list[jk]->GetFCPC();
	      pAnalysis->trkmm[ij] = pfitTrack->InoTrackCand_list[jk]->GetMomentum();
	      pAnalysis->trkth[ij] = pfitTrack->InoTrackCand_list[jk]->GetTheta();
	      pAnalysis->trkph[ij] = pfitTrack->InoTrackCand_list[jk]->GetPhi();
	      // cout<<"-------------------------------------------------------------------------"<<endl;
	      // cout<<"Reconstructed P = "<<pAnalysis->trkmm[ij]<<"  |  theta "<<pAnalysis->trkth[ij]*180/3.1415<<"  |  phi "<<pAnalysis->trkph[ij]*180/3.1415<<endl;
	      // cout<<"-------------------------------------------------------------------------"<<endl;
	      pAnalysis->therr[ij] = pfitTrack->InoTrackCand_list[jk]->GetThErr();
	      pAnalysis->pherr[ij] = pfitTrack->InoTrackCand_list[jk]->GetPhErr();
	    
	      pAnalysis->momvx[ij] = pfitTrack->InoTrackCand_list[jk]->GetMomentumCurve();
	      pAnalysis->thevx[ij] = acos(pfitTrack->InoTrackCand_list[jk]->GetDirCosZ());
	      pAnalysis->phivx[ij] = atan2(pfitTrack->InoTrackCand_list[jk]->GetDirCosV(),
					   pfitTrack->InoTrackCand_list[jk]->GetDirCosU());
	      pAnalysis->posxvx[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxU();
	      pAnalysis->posyvx[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxV();
	      pAnalysis->poszvx[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxZ();
	    
	      pAnalysis->momend[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndMomentumCurve();
	      pAnalysis->theend[ij] = acos(pfitTrack->InoTrackCand_list[jk]->GetEndDirCosZ());
	      pAnalysis->phiend[ij] = atan2(pfitTrack->InoTrackCand_list[jk]->GetEndDirCosV(),
					    pfitTrack->InoTrackCand_list[jk]->GetEndDirCosU());
	      pAnalysis->tx_end[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndDirCosU();
	      pAnalysis->ty_end[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndDirCosV();
	    
	      pAnalysis->posxend[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndU();
	      pAnalysis->posyend[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndV();
	      pAnalysis->poszend[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndZ();
	      /*
		pAnalysis->momds[ij] = pfitTrack->InoTrackCand_list[jk]->GetMomentumdS()
		+ pfitTrack->InoTrackCand_list[ij]->GetdSExtra(); //GMA -ve value for lower momentum !!!!!!
		pAnalysis->momrg[jk] = pfitTrack->InoTrackCand_list[ij]->GetMomentumRange()
		+  pfitTrack->InoTrackCand_list[ij]->GetRangeExtra();
	      */
	      pAnalysis->momds[ij] = pfitTrack->InoTrackCand_list[jk]->GetMomentumdS();  //GMA14 was bug
	      pAnalysis->momrg[ij] = pfitTrack->InoTrackCand_list[jk]->GetMomentumRange();
	      
	      // cout<<"micalEventAction ij "<<ij<<endl;
	      // cout<<"pAnalysis->trkmm "<<pAnalysis->trkmm[ij]<<endl;
	      // cout<<"pAnalysis->momds "<<pAnalysis->momds[ij]<<endl;
	      // cout<<"pAnalysis->momrg "<<pAnalysis->momrg[ij]<<endl;
	      pAnalysis->vtxzplane[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxPlane();
	      pAnalysis->endzplane[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndPlane();
            
	      pAnalysis->xxerr[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxUError(); //GMA14 all these 15 tesrms
	      pAnalysis->yyerr[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxVError();
	      pAnalysis->txerr[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxdUError();
	      pAnalysis->tyerr[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxdVError();
	      pAnalysis->qperr[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxQPError();
	    
	      pAnalysis->xxenderr[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndUError();
	      pAnalysis->yyenderr[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndVError();
	      pAnalysis->txenderr[ij] = pfitTrack->InoTrackCand_list[jk]->GetEnddUError();
	      pAnalysis->tyenderr[ij] = pfitTrack->InoTrackCand_list[jk]->GetEnddVError();
	      pAnalysis->qpenderr[ij] = pfitTrack->InoTrackCand_list[jk]->GetEndQPError();
	    
	      pAnalysis->xxin[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxXX();
	      pAnalysis->yyin[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxYY();
	      pAnalysis->txin[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxTX();
	      pAnalysis->tyin[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxTY();
	      pAnalysis->tx[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxdU();
	      pAnalysis->ty[ij] = pfitTrack->InoTrackCand_list[jk]->GetVtxdV();
	    
	      pAnalysis->mcxgnvx[ij]=1.e-3*MCpos.x();
	      pAnalysis->mcygnvx[ij]=1.e-3*MCpos.y();
	      //Initialise with this mainly for noise track
	      
	      // cout<<"---------------------------------------------------------------"<<endl;
	      // cout<<"p = "<<pAnalysis->trkmm[ij]<<", momds = "<<pAnalysis->momds[ij]<<", E_mu = "<<pAnalysis->momin[0]<<endl;
	      // cout<<"---------------------------------------------------------------"<<endl;

	      pAnalysis->momgnvx[ij]  = pAnalysis->momgnend[ij] = 0.0;
	      pAnalysis->thegnvx[ij] = pAnalysis->thegnend[ij] = 10.;
	      pAnalysis->phignvx[ij] = pAnalysis->phignend[ij] = 10.;
	    
	      vector<InoCluster*> tmpclusts = pfitTrack->InoTrackCand_list[jk]->ClustsInTrack; // fTrack->ClustsInTrack;
	    
	      int plane = pfitTrack->InoTrackCand_list[jk]->GetVtxPlane();
	    
	      // cout <<"genplane "<<plane<<endl;
	      for (unsigned kl = 0; kl <tmpclusts.size(); kl++)	{
		if (tmpclusts[kl]->GetZPlane()==plane) {
		  pAnalysis->momgnvx[ij] = 0.001*tmpclusts[kl]->HitsInCluster[0]->GetMomentum();
		  pAnalysis->thegnvx[ij] = tmpclusts[kl]->HitsInCluster[0]->GetTheta();
		  pAnalysis->phignvx[ij] = tmpclusts[kl]->HitsInCluster[0]->GetPhi();
		  break;
		} // if (tmpclusts[kl]->GetZPlane()==plane)
	      } //  for (unsigned kl = 0; kl <tmpclusts.size(); kl++)
	      int trkclustersize = tmpclusts.size();
	      int tottrkhit=0;
	      int trkxstrp = 0;
	      int trkystrp = 0;
	      for (int kl=0; kl<trkclustersize; kl++) {
		vector<InoHit*>  tmphit = tmpclusts[kl]->HitsInCluster;
		int nhits = tmphit.size();
		int zplane =tmphit[0]->GetZPlane();
	      
		for (int lm=0; lm<nhits; lm++) {
		  int ixstripno = tmphit[lm]->GetXStripNum();
		  int iystripno = tmphit[lm]->GetYStripNum();
		  if (tmphit[lm]->GetXStrip()) {
		    for (unsigned nn=0; nn< xtrkstrip.size(); nn++) {
		      //if (istripno == xtrkstrip[nn]) {istripno = -1; break;}
		      if (ixstripno == xtrkstrip[nn].second &&
			  zplane == xtrkstrip[nn].first) {ixstripno = -1; break;}
		    }
		  }
		  if (tmphit[lm]->GetYStrip()) {
		    for (unsigned nn=0; nn< ytrkstrip.size(); nn++) {
		      //if (istripno == ytrkstrip[nn]) {istripno = -1; break;}
		      if (iystripno == xtrkstrip[nn].second &&
			  zplane == xtrkstrip[nn].first) {iystripno = -1; break;}
		    }
		  }
		  if (ixstripno>=0 && iystripno>=0) {
		    trkxstrp++;  trkystrp++;
		    if (pAnalysis->cvalue[ij]>0){ixyz Zxstrip(zplane,ixstripno) ; xtrkstrip.push_back(Zxstrip);}
		    if (pAnalysis->cvalue[ij]>0){ixyz Zystrip(zplane,iystripno) ; ytrkstrip.push_back(Zystrip);}
		  }
		} // for (int lm=0; lm<nhits; lm++)
		//tmphit.clear();
		tottrkhit +=nhits;
	      } // for (int kl=0; kl<trkclustersize; kl++)
	    
	      pAnalysis->ntrkcl[ij] = 1000*trkclustersize + tottrkhit;
	      pAnalysis->ntrkst[ij] = 1000*trkxstrp + trkystrp;
	      tottrkXstrp +=trkxstrp;
	      tottrkYstrp +=trkystrp;
	    
	      plane = pfitTrack->InoTrackCand_list[jk]->GetEndPlane();
	    
	      for (unsigned kl = 0; kl <tmpclusts.size(); kl++) {
		if (tmpclusts[kl]->GetZPlane()==plane) {
		  pAnalysis->momgnend[ij] = 0.001*tmpclusts[kl]->HitsInCluster[0]->GetMomentum();
		  pAnalysis->thegnend[ij] = tmpclusts[kl]->HitsInCluster[0]->GetTheta();
		  pAnalysis->phignend[ij] = tmpclusts[kl]->HitsInCluster[0]->GetPhi();
		  break;
		}
	      }
	      //tmpclusts.clear();
	    } // if (ij <pAnalysis->ntrkmx)
	    ij++;
	  } //for (unsigned jk=0; jk<pfitTrack->InoTrackCand_list.size() ; jk++)
      
	  //	  cout <<"ntrack "<< ij<<endl;
	  pAnalysis->ntrkt = ij; //(pfitTrack->InoTrackCand_list.size() <=pAnalysis->ntrkmx) ? pfitTrack->InoTrackCand_list.size() : pAnalysis->ntrkmx;
	
	  //Now add other cluster close to vertex;
	  //	  vector <InoCluster*> totcluster = trackfinder.inoCluster_pointer->InoCluster_list;
	  //	  int totclustersize = totcluster.size();
	  //cout<<"totclustersize "<<totclustersize<<endl;
	  //int trackhits=0; //asm  temp variable
	  //Tag all clusters which belongs to any leading track

	  for (int jk=0; jk<totclustersize; jk++) {
	  
	    totcluster[jk]->SetInTrack(0);
	    totcluster[jk]->SetInShower(0);
	  
	    for (unsigned kl=0; kl<pfitTrack->InoTrackCand_list.size() ; kl++) {
	      if(((pfitTrack->InoTrackCand_list[0]->GetChi2()/pfitTrack->InoTrackCand_list[0]->GetNDOF())<30 ||
		  (abs(pfitTrack->InoTrackCand_list[0]->GetMomentum())>0. &&
		   abs(pfitTrack->InoTrackCand_list[0]->GetMomentum())<1000 ))) { // ||
		//		(pfitTrack->InoTrackCand_list[0]->Getcval()>0))) {

		for (unsigned lm=0; lm< pfitTrack->InoTrackCand_list[kl]->GetClusterEntries(); lm++) {
		  if (totcluster[jk]->isIdentical(pfitTrack->InoTrackCand_list[kl]->ClustsInTrack[lm])) {
		  
		    if((kl>0&&abs(pfitTrack->InoTrackCand_list[kl]->GetVtxZ()-pfitTrack->InoTrackCand_list[0]->GetVtxZ())<.1)) {
		      totcluster[jk]->SetInShower(1);
		    
		      vector<InoHit*> tmphit = totcluster[jk]->HitsInCluster;
		      int nhits = tmphit.size();
		      int zplane =tmphit[0]->GetZPlane();
		    
		      for (int mn=0; mn<nhits; mn++) {
			int ixstripno = tmphit[mn]->GetXStripNum();
			int iystripno = tmphit[mn]->GetYStripNum();
			//		      if ((ixstripno==0 && iystripno==0) ||( ixstripno>=0 && iystripno>=0) ) {
			if (ixstripno>=0 && iystripno>=0) {
			  ixyz ZXstrp(zplane,ixstripno); ixyz ZYstrp(zplane,iystripno);
			  xshwstrip.push_back(ZXstrp); yshwstrip.push_back(ZYstrp);
			}
			if (tmphit[mn]->GetXStrip() && tmphit[mn]->GetYStrip() ) {
			  //tmphit[mn]->SetUID(-442);// SSE 08/15
			  if (pAnalysis->isVisOut) {
			    pAnalysis->H->NShowerHits++;
			    pAnalysis->Hp= pAnalysis->H->AddHits(0,0); // add a track object  //VALGRIND
			    pAnalysis->Hp->ShowerHitNum++;
			    pAnalysis->Hp->ZZ=  zplane;//pfitTrack->InoTrackCand_list[i]->ClustsInTrack[j]->GetZPlane(); //os();
			    pAnalysis->Hp->XX=tmphit[mn]->GetXPos();
			    pAnalysis->Hp->YY=tmphit[mn]->GetYPos();
			    pAnalysis->Hp->TrackType=-442; // track type up -44 + 0 FIT UP
			    pAnalysis->Hp->Fup=2;// Shower hits from an event
			  } //isVis
			}
		      }
		    }  else {
		      totcluster[jk]->SetInTrack(1);
		    }
		  }
		}
	      }
	    }
	  } // for (int jk=0; jk<totclustersize; jk++)
	
	  {int nc_cl=0; int nc_ht=0; int ccc_cl=0; int ccc_ht=0 ; int t_ht=0;
	    for (int jk=0; jk<totclustersize; ++jk) {
	      if(totcluster[jk]->GetInTrack()==1&&totcluster[jk]->GetInShower()!=1) {
		nc_cl++; nc_ht += totcluster[jk]->GetHitEntries();
	      } else {
		ccc_cl++; ccc_ht+=totcluster[jk]->GetHitEntries();
	      }
	      t_ht+=totcluster[jk]->GetHitEntries();
	    }
	    //	    cout <<"before implemeting : " <<"hits not in track:"<< ccc_ht
	    //	 << " + track hits: "<<  nc_ht <<" total " << t_ht <<endl;//SSE 29 Oct 2015
	    if(pAnalysis->isVisOut==1){
	      pAnalysis->ascii_output <<endl;
	      pAnalysis->ascii_output <<"before implemeting : " <<"hits not in track:"<< ccc_ht << " + track hits: "<<  nc_ht <<" total " << t_ht <<endl;
	    }
	  }

	  //	  cout<< "totclustersize  " <<totclustersize<<  endl;
	  //====================================================================
	  //  calculation of hadron energy //SSE mofied and added here 291015
	  //Calculation of hadron momentum and direction  from simulation input Oct8, 15
	  //====================================================================
	
	  double E_had;
	  float Px_nu=0.0, Py_nu=0.0, Pz_nu=0.0; //SSE 081015
	  float Px_muon=0.0, Py_muon=0.0, Pz_muon=0.0; //SSE 081015
	  float  thetain_had=0.0, phiin_had=0.0;//,pmagin_had=0.0; //SSE 081015
	  //float  costhetain_had=0.0; //SSE 081015
	  //float costheta_had=0.0;
	  float dotAngle=0.0;
	  G4ThreeVector tmp3pin_nu(0.,0.,0.);
	  G4ThreeVector tmp3pin_muon(0.,0.,0.);
	  G4ThreeVector tmp3pin_had(0.,0.,0.);
	  if (abs(pAnalysis->pidin[1])==13||abs(pAnalysis->pidin[1])==14){
	    Px_muon=(pAnalysis->momin[1])*(TMath::Sin(pAnalysis->thein[1]))*(TMath::Cos(pAnalysis->phiin[1]));
	    Py_muon=(pAnalysis->momin[1])*(TMath::Sin(pAnalysis->thein[1]))*(TMath::Sin(pAnalysis->phiin[1]));
	    Pz_muon=(pAnalysis->momin[1])*(TMath::Cos(pAnalysis->thein[1]));
	    E_had=abs(pAnalysis->momin[0])- abs(pAnalysis->momin[1]);
	  } else {
	    E_had=pAnalysis->momin[0];
	  }

	  //=========================cal. of  hadron momentum and direction=================================
	
	  Px_nu=(pAnalysis->momin[0])*(TMath::Sin(pAnalysis->thein[0]))*(TMath::Cos(pAnalysis->phiin[0]));
	  Py_nu=(pAnalysis->momin[0])*(TMath::Sin(pAnalysis->thein[0]))*(TMath::Sin(pAnalysis->phiin[0]));
	  Pz_nu=(pAnalysis->momin[0])*(TMath::Cos(pAnalysis->thein[0]));
	
	  //G4ThreeVector 
	  tmp3pin_nu=G4ThreeVector (Px_nu,Py_nu, Pz_nu);
	  //G4ThreeVector 
	  tmp3pin_muon= G4ThreeVector (Px_muon,Py_muon, Pz_muon);
	  //G4ThreeVector 
	  tmp3pin_had=tmp3pin_nu-tmp3pin_muon;
	
	  //pmagin_had=tmp3pin_had.mag();
	  thetain_had=tmp3pin_had.theta();
	  //costhetain_had=tmp3pin_had.cosTheta();
	  phiin_had=tmp3pin_had.phi();
	
	  pAnalysis->theta_hadron_in=0;//SSE 081015
	  //pAnalysis->costheta_hadron_in=0;//SSE 081015
	  pAnalysis->phi_hadron_in=0;//SSE 081015
	  pAnalysis->theta_hadron_in=thetain_had ;//SSE 081015
	  //pAnalysis->costheta_hadron_in=costhetain_had ;//SSE 081015
	  pAnalysis->phi_hadron_in=phiin_had ;//SSE 081015
	
	  // cout<< " nu p in: x, y, z "<< tmp3pin_nu[0]<< " "<<tmp3pin_nu[1]<<" "<< tmp3pin_nu[2]<<endl;
	  // cout<< " mu p in: x, y, z "<< tmp3pin_muon[0]<< " "<<tmp3pin_muon[1]<<" "<< tmp3pin_muon[2]<<endl;
	  // cout<< " had p in: x, y, z "<< tmp3pin_had[0]<< " "<<tmp3pin_had[1]<<" "<< tmp3pin_had[2]<<endl;
	  // cout<< " had theta phi  input: x, y, z "<< thetain_had<< " "<<phiin_had<<endl;
	  //============= end of cal. of  hadron momentum and direction ==================================
	
	  // cout<<"E_had"<<" "<<E_had<<" "<<pAnalysis->momin[0]<<" "<<pAnalysis->momin[1]<<endl;//SSE 291015
	
	  //======================================================================
	  //Add flag through SetUID (Aug 2015; SSE) added here 291015
	  //
	  //
	  //=======================================================================
	  // cout<<"pinotrack->InoTrack_list.size() = "<<pinotrack->InoTrack_list.size()<<endl;
	  for (unsigned ixj=0; ixj<pinotrack->InoTrack_list.size() ; ixj++) {
	    // cout<<"Inside Loop 1, pinotrack->InoTrack_list[ixj]->ClustsInTrack.size() = "<<pinotrack->InoTrack_list[ixj]->ClustsInTrack.size()<<endl;
	    for ( unsigned int jxk =0; jxk<pinotrack->InoTrack_list[ixj]->ClustsInTrack.size();jxk++) {	
	      // cout<<"Inside Loop 2, pinotrack->InoTrack_list[ixj]->ClustsInTrack[jxk]->GetHitEntries() = "<<pinotrack->InoTrack_list[ixj]->ClustsInTrack[jxk]->GetHitEntries()<<endl;
	      for(unsigned int mxn=0;mxn<pinotrack->InoTrack_list[ixj]->ClustsInTrack[jxk]->GetHitEntries();mxn++){
		// cout<<"Inside Loop 3"<<endl;
		pinotrack->InoTrack_list[ixj]->ClustsInTrack[jxk]->HitsInCluster[mxn]->SetUID(-4);//SSE 08/15
	      }
	    }
	  }
	
	  // cout<<"pfitTrack->InoTrackCand_list.size() = "<<pfitTrack->InoTrackCand_list.size()<<endl;
	  for (unsigned kxl=0; kxl< pfitTrack->InoTrackCand_list.size() ; kxl++) {
	    // cout<<"Inside Loop 1, pfitTrack->InoTrackCand_list[kxl]->GetClusterEntries() = "<<pfitTrack->InoTrackCand_list[kxl]->GetClusterEntries()<<endl;
	    for ( unsigned int lxm =0; lxm<pfitTrack->InoTrackCand_list[kxl]->GetClusterEntries(); lxm++){
	      // cout<<"Inside Loop 2, pfitTrack->InoTrackCand_list[kxl]->ClustsInTrack[lxm]->GetHitEntries() = "<<pfitTrack->InoTrackCand_list[kxl]->ClustsInTrack[lxm]->GetHitEntries()<<endl;
	      for(unsigned int mxn=0; mxn<pfitTrack->InoTrackCand_list[kxl]->ClustsInTrack[lxm]->GetHitEntries();mxn++){
		// cout<<"Inside Loop 3"<<endl;
		if(pfitTrack->InoTrackCand_list[kxl]->GetFitType()==1){
		  // cout<<"UID == 10"<<endl;
		  pfitTrack->InoTrackCand_list[kxl]->ClustsInTrack[lxm]->HitsInCluster[mxn]->SetUID(10);	//SSE 08/15
		} else {
		  // cout<<"UID == 11"<<endl;
		  pfitTrack->InoTrackCand_list[kxl]->ClustsInTrack[lxm]->HitsInCluster[mxn]->SetUID(11);
		}
	      }
	    }
	  }
	  // cout<<"Outside Loop"<<endl;
	  //=========================================================


	  int Orighits_cluster=0;//SSE Hits in cluster + strip finding algo
	  vector<InoHit*> tmphit_trape;////SSE 301015  
	  int Hit_wo_ghst=0;
	    
	  int Orighits_trape=0;//hits in trapezoid + strip finding algo 
	  int wogh_Orighits=0;
	  int num_hadron=0;//SSE 250915
	  
	  float nhits_large_clust_selected=0;// SSE 091015
	  
	  double matrix_XX=0.0, matrix_YY=0.0, matrix_ZZ=0.0;//SSE 031115
	  double Matrix_xx=0.0, Matrix_yy=0.0, Matrix_zz=0.0;//SSE 031115
	  double Matrix_xy=0.0, Matrix_yz=0.0, Matrix_zx=0.0;//SSE 031115
	  
	  
	  int minplane_cluster=0;//SSE 091015
	  
	  
    
	  int allremhit = 0;
	  int allremcls = 0;
	  double x0=0,y0=0, z0=0;
	  if(pfitTrack->InoTrackCand_list.size() >0) {
	    x0 = pfitTrack->InoTrackCand_list[0]->GetVtxU();  // asm:  vertex of longest track
	    y0 = pfitTrack->InoTrackCand_list[0]->GetVtxV();   //asm
	    z0 = pfitTrack->InoTrackCand_list[0]->GetVtxZ();  
	    // cout<<" chi2 "<< (pfitTrack->InoTrackCand_list[0]->GetChi2()/pfitTrack->InoTrackCand_list[0]->GetNDOF())<<endl;
	    // cout<<	" abs(pfitTrack->InoTrackCand_list[0]->GetMomentum())  "<< abs(pfitTrack->InoTrackCand_list[0]->GetMomentum())<<endl;	
	    // cout<<	" pfitTrack->InoTrackCand_list[0]->Getcval()  "<< pfitTrack->InoTrackCand_list[0]->Getcval()<<endl;	
	  }
	  // cout<<" x0 "<< x0<< " y0 "<<y0<<" "<<pfitTrack->InoTrackCand_list.size()<<endl;//SSE 281015
	  float mXX=0; float  mYY=0; float mZZ=0;
	  float Mxx =0; float Myy=0; float Mzz =0;
	  float Mxy =0; float Mxz=0; float Myz =0;
	  //	commented by SSE because cval= 0 or -0

	  //   &&(pfitTrack->InoTrackCand_list[0]->Getcval()>0)){
	  if( (x0!=0 || y0!=0 ) && 
	      (pfitTrack->InoTrackCand_list.size()>0) &&
	      (pfitTrack->InoTrackCand_list[0]->GetChi2()/pfitTrack->InoTrackCand_list[0]->GetNDOF())<30 &&  //add a function for nhits
	      abs(pfitTrack->InoTrackCand_list[0]->GetMomentum())>0. && 
	      abs(pfitTrack->InoTrackCand_list[0]->GetMomentum())<1000000) {
	    
	    // cout<<"in loop:  x0 "<< x0<< " y0 "<<y0<<endl;//SSE 281015
	    for (unsigned kl=0; kl<1 ; kl++) {
	      int plane = pfitTrack->InoTrackCand_list[kl]->fVertex->GetPlane();
	      //TVector3 direction = pfitTrack->InoTrackCand_list[kl]->fVertex->GetDirCosine();
	    
	      double theta = pfitTrack->InoTrackCand_list[0]->GetTheta(); if (abs(theta<1.e-4)) {theta=0;};
	      double phi = pfitTrack->InoTrackCand_list[0]->GetPhi();   if (abs(phi<1.e-4)) {phi=0;};
	      double dxdz = tan(theta)*cos(phi);     if( abs(dxdz)<1.e-4) {dxdz=0;}
	      double dydz = tan(theta)*sin(phi);     if( abs(dydz)<1.e-4) {dydz=0;}
	      // cout<< " theta "<< theta<<" phi  "<< phi<< " dxdz "<< dxdz<<" dydz "<<dydz<<endl;//SSE 28Oct 
	    
	      int dir = (cos(theta)>0) ? 1 : -1;  //asm: check this
	      int plusplane=5;
	      int minusplane=-1;
	      double basewindow1 = 0.05, basewindow=0.05;// 10 cm
	    
	      double slopex = 1.0  ; //dz/dx(dy) = 1
	      double slopey = 1.0  ;
	    
	      if (abs(dxdz)>1) {slopex = 1.0;} else {slopex = 2.0;}
	      if (abs(dydz)>1) {slopey = 1.0;} else {slopey = 2.0;}
	      if(pAnalysis->isXtermOut==2) {
		pAnalysis->ascii_output <<"-------------------------------------------------------------"<<endl;
		pAnalysis->ascii_output <<" x0:y0 "<< x0 <<":"<<y0<< "  dxdz:dydz "<<  dxdz <<":" <<dydz<<endl;
	      }
	      int collected_hits=0;
	      do {
		collected_hits=0;//Trk->SetStraight();
		//inoTrack_pointer->InoTrack_list.push_back(Trk);
		basewindow1 +=0.05;
		if (basewindow1<0.2){
		  plusplane=5; minusplane=-2;
		} else if (basewindow1>=0.2){
		  plusplane+=1; minusplane = -3;
		}
	      
		for (int jxk=minusplane; jxk<plusplane; jxk++) {
		  int newplane = plane + jxk*dir;
		  if (newplane<0 && newplane>=nLayer) continue;
		
		  if (jxk<0) {
		    basewindow =basewindow1+0.05;
		  } else if(jxk>=0) {
		    basewindow=basewindow1;
		  }
		
		  for (int lm=0; lm<totclustersize; lm++) {
		    // cout<<"Loop lm at 739 = "<<lm<<endl;
		    if (totcluster[lm]->GetZPlane() !=newplane) continue;
		    if (totcluster[lm]->GetInTrack()) continue; //remove clusters in tracks
		  
		    double xmn = x0 +  LayerThickness*(jxk+1)*dir*(dxdz - dir*slopex) - basewindow;
		    double xmx = x0 +  LayerThickness*(jxk+1)*dir*(dxdz + dir*slopex) + basewindow;
		  
		    double ymn = y0 + LayerThickness*(jxk+1)*dir*(dydz - dir*slopey) - basewindow;
		    double ymx = y0 + LayerThickness*(jxk+1)*dir*(dydz + dir*slopey) + basewindow;
		  
		    if (totcluster[lm]->GetBegXPos() > xmn &&
			totcluster[lm]->GetEndXPos() < xmx &&
			totcluster[lm]->GetBegYPos() > ymn &&
			totcluster[lm]->GetEndYPos() < ymx &&
			!totcluster[lm]->GetInShower() && 
			!totcluster[lm]->GetInTrack()) {
		    
		      allremhit +=totcluster[lm]->GetHitEntries();
		      allremcls +=1;   totcluster[lm]->SetInShower(1);
		    
		      vector<InoHit*> tmphit = totcluster[lm]->HitsInCluster;
		      int nhits = tmphit.size();
		      int zplane =tmphit[0]->GetZPlane();
		    
		    
		      for (int mn=0; mn<nhits; mn++) {
			// cout<<" UID " <<tmphit[mn]->GetUID()<<endl;
			int ixstripno = tmphit[mn]->GetXStripNum();
			int iystripno = tmphit[mn]->GetYStripNum();
		      
			if (tmphit[mn]->GetXStrip() && tmphit[mn]->GetYStrip() ) {
			  if (tmphit[mn]->GetXPos() < xmn || tmphit[mn]->GetXPos() > xmx) continue;
			  if (tmphit[mn]->GetYPos() < ymn || tmphit[mn]->GetYPos() > ymx) continue;

			  for (unsigned nn=0; nn< xtrkstrip.size(); nn++) {
			    if (ixstripno == xtrkstrip[nn].second &&
				zplane == xtrkstrip[nn].first) {ixstripno = -1; break;}
			  }
			
			  for (unsigned nn=0; nn< ytrkstrip.size(); nn++) {
			    if (iystripno == ytrkstrip[nn].second &&
				zplane == ytrkstrip[nn].first) {
			      iystripno = -1; break;
			    }
			  }
			  for (unsigned nn=0; nn< xshwstrip.size(); nn++) {
			    if (ixstripno == xshwstrip[nn].second &&
				zplane == xshwstrip[nn].first) {
			      ixstripno = -1; break;
			    }
			  }
			 
			  for (unsigned nn=0; nn< yshwstrip.size(); nn++) {
			    if (iystripno == yshwstrip[nn].second &&
				zplane == yshwstrip[nn].first) {
			      iystripno = -1; break;
			    }
			  }
			
			  if((ixstripno>=0 || iystripno>=0) ) {
			    tmphit_trape.push_back(tmphit[mn]);//SSE 291015

			    if(ixstripno>=0) {
			      ixyz ZXstrp(zplane,ixstripno);xshwstrip.push_back(ZXstrp);
			    }
			    if(iystripno>=0){
			      ixyz ZYstrp(zplane,iystripno);yshwstrip.push_back(ZYstrp);
			    }
			    collected_hits++;
			    
			  }
		       
			  //tmphit[mn]->SetUID(-442); //SSE 08/15	
			  if (pAnalysis->isVisOut) {
			    pAnalysis->H->NShowerHits++;
			    pAnalysis->Hp= pAnalysis->H->AddHits(0,0); // add a track object  //VALGRIND
			    pAnalysis->Hp->ShowerHitNum++;
			    pAnalysis->Hp->ZZ=  zplane;//pfitTrack->InoTrackCand_list[i]->ClustsInTrack[j]->GetZPlane(); //os();
			    pAnalysis->Hp->XX=tmphit[mn]->GetXPos();
			    pAnalysis->Hp->YY=tmphit[mn]->GetYPos();
			    pAnalysis->Hp->TrackType=-442; // track type up -44 + 0 FIT UP
			    pAnalysis->Hp->Fup=2;// Shower hits from an event
			    mXX=tmphit[mn]->GetXPos();
			    mYY=tmphit[mn]->GetYPos();
			    mZZ=tmphit[mn]->GetZPos();
			    Mxx += mXX*mXX;
			    Mxy += mXX*mYY;
			    Mxz += mXX*mZZ;
			    Myy += mYY*mYY;
			    Myz += mYY*mZZ;
			    Mzz += mZZ*mZZ;
			  
			  } //isVis
			} //if (tmphit[mn]->GetXStrip() && tmphit[mn]->GetYStrip() )
		      } //for (int mn=0; mn<nhits; mn++)
		      //	      tmphit.clear();
		    } //if (totcluster[lm]->GetBegXPos() > xmn .....)
		    // cout<< " trapezoidal cluster size Oct 27, 2015  size x & y"<< xshwstrip.size ()<< " "<<yshwstrip.size () <<endl;
		  } //for (int lm=0; lm<totclustersize; lm++)//asm: Dec7 : need to cleat x/ytrkstrip  for the next plane
		} //for (int jxk=minusplane; jxk<plusplane; jxk++)
		
	      } while(collected_hits>4||basewindow<=0.35);
	    
	    } //for (unsigned kl=0; kl<1 ; kl++) c: //sse 
	    // WE need to call orighit_calc function here
	    
	     Orighits_trape = orighit_calc(tmphit_trape);//SSE 30102015
	    
	    //  histogram filling for Orighits
	   // pAnalysis->hist_orighits_mod_E->Fill(Orighits_mod, E_had);
	    
	    
	    
	   
	    //=================================================
	    //Added on 031115 
	    //to Find hadron shower direction
	    //=================================================
	    
	    for( unsigned int ijk=0; ijk<tmphit_trape.size(); ijk++) {
	      //need min two planes including muon vetrex plane
	      if (ijk>0 && pfitTrack->InoTrackCand_list[0]->fVertex->GetPlane()!=tmphit_trape[ijk]->GetZPlane()){
		minplane_cluster++;
		break;
	      }
	    } 
	    
	    if (minplane_cluster>0){// lagest cluster must contains min two planes.
	      for( unsigned int mn=0; mn<tmphit_trape.size(); mn++) {
		
		matrix_XX=(tmphit_trape[mn]->GetXPos() - x0);
		matrix_YY=(tmphit_trape[mn]->GetYPos() - y0);
		matrix_ZZ=(tmphit_trape[mn]->GetZPos() - z0);

			// cout<< " xx  "<< tmphit_trape[mn]->GetXPos()
			// << " yy  "<< tmphit_trape[mn]->GetYPos()
			// << " zz  "<< " "<<tmphit_trape[mn]->GetZPos()<<endl<<endl;
		
			// cout<< " matrix_XX "<< matrix_XX<< " matrix_YY "<< matrix_YY<< "  matrix_ZZ "<<  matrix_ZZ<<endl;			 
        
		Matrix_xx +=matrix_XX*matrix_XX;
		Matrix_xy +=matrix_XX*matrix_YY;
		Matrix_yy +=matrix_YY*matrix_YY;
		Matrix_yz +=matrix_YY*matrix_ZZ;
		Matrix_zz +=matrix_ZZ*matrix_ZZ;
		Matrix_zx +=matrix_ZZ*matrix_XX;
		
		//cout<< " Matrix_xx  " << Matrix_xx<< " Matrix_xy "<<Matrix_xy <<
		//	" Matrix_yy " << Matrix_yy<< " Matrix_yz "<<Matrix_yz <<
		//	" Matrix_zz " << Matrix_zz<< " Matrix_zx "<<Matrix_zx << endl;
		
	      }


	      // cout<< " Matrix_xx  " << Matrix_xx<< " Matrix_xy "<<Matrix_xy <<
	      // 	" Matrix_yy " << Matrix_yy<< " Matrix_yz "<<Matrix_yz <<
	      // 	" Matrix_zz " << Matrix_zz<< " Matrix_zx "<<Matrix_zx << endl;
	

		
	      }// end of if (minplane_cluster>0
	    
	  } else { // of (x0!=0||y0!=0)&&(pfitTrack->InoTrackCand_list.size()!=0) // finding a good recosntructed muon 
	    
	    //========================================================================
	    //  To remove event if a layer has more than 100 hits //SSE 09/15
	    //========================================================================
	    
	  

	    InoHit_Manager* tmphitlist = InoHit_Manager::APointer;
	    // cout<<"2E_had"<<" "<<E_had<<endl;//SSE
	    int totalhits = 0; // tmphitlist->InoHit_list.size();
	    vector<InoHit*> tmphadcluster1;
	    
	    for (int ixj =0; ixj<nLayer;ixj++){
	      HitBank_All[ixj].clear();
	    }
	    // cout <<"xx1 "<<endl;
	    // cout<<"tmphitlist->InoHit_list.size() = "<<tmphitlist->InoHit_list.size()<<endl;
	    for( unsigned int iyj=0; iyj<tmphitlist->InoHit_list.size(); iyj++) {
	      // cout<<"Inside Loop "<<endl;
	      // cout <<"iyj "<< iyj<<endl;
	      int nxstrip=0, nystrip=0;
	      
	      int irpcmod = tmphitlist->InoHit_list[iyj]->GetRPCmod();
	      // cout <<"ixxx "<<irpcmod <<endl;
	      // cout<<"<inoRPC_pointer- "<<(inoRPC_pointer)<<endl;
	      // cout<<inoRPC_pointer->InoRPCStrip.size()<<endl;

	      for (unsigned int ix=0; ix<inoRPC_pointer->InoRPCStrip.size(); ix++) {
		// cout <<"ix "<< ix<<" "<<irpcmod <<endl;
		if (inoRPC_pointer->InoRPCStrip[ix].first==irpcmod) {
		  nxstrip = inoRPC_pointer->InoRPCStrip[ix].second/100;
		  nystrip = inoRPC_pointer->InoRPCStrip[ix].second%100;
		  break;
		}
	      }
	      // cout <<"x2x2 "<<endl;
	      if (nxstrip*nystrip <100) {
		int plane = tmphitlist->InoHit_list[iyj]->GetZPlane();
		HitBank_All[plane].push_back(tmphitlist->InoHit_list[iyj]);
		tmphadcluster1.push_back(tmphitlist->InoHit_list[iyj]);
	      }
	    }
	    // cout <<"xx2 "<<endl;
	  //  pAnalysis->total_inohits      =0;
	   // pAnalysis->orighits_new	=0;//SSE 250915
	    //pAnalysis->orighits_mod	=0;//SSE 250915
	    //pAnalysis->hit_wo_ghst	=0;//SSE 250915
	    //pAnalysis->hit_wogh_orighits       =0;//SSE 250915
	    //pAnalysis->e_hadron=0;//SSE 250915
	    //pAnalysis->nhits_largest_cluster=0;//SSE 250915
	    //pAnalysis->theta_hadron_shw=0;//SSE 071015
	   // pAnalysis->costheta_hadron_shw=0;//SSE 081015
	    //pAnalysis->phi_hadron_shw=0;//SSE 071015
	    //pAnalysis-> dot_angle_had_shw=0.0;
	   // pAnalysis-> nhits_largest_cluster_selected=0; //SSE 091015	
	    vector<InoHit*> tmphadcluster;
	    vector < vector <InoHit* > > allhadcluster;
	    
	    //pAnalysis->total_inohits=tmphadcluster1.size();
	    // cout<< " before removal :  "<<tmphitlist->InoHit_list.size()<<endl;
	    

	    // cout<< " tmphadcluster1.size() "<<tmphadcluster1.size()<<endl;
	    for (unsigned int ixx=0; ixx<tmphadcluster1.size(); ixx++) {
	      // cout<<" tmphadcluster1[ixx]->GetUID()  "<< tmphadcluster1[ixx]->GetUID()<<endl;
	      if (tmphadcluster1[ixx]->GetUID()==10 || tmphadcluster1[ixx]->GetUID()==11) {
		tmphadcluster1.erase(tmphadcluster1.begin()+ixx);
		ixx--;
	      }
	    }
	    // cout<< " tmphadcluster1.size() "<<tmphadcluster1.size()<<endl;
	    //-------------------------------------------------
	    //orighit calculation 
	    //
	    //-------------------------------------------------
	    //Orighits_new = orighit_calc(tmphadcluster1);
	    // cout<<" orighit function: new orighit "<<Orighits_new<<endl;

	    //cout<< " after removal :  "<<tmphitlist->InoHit_list.size()<<" totalhits before removal "<< totalhits<<endl;
	    //totalhits = tmphitlist->InoHit_list.size();	
	    totalhits=tmphadcluster1.size();
	    // cout<<"tmphadcluster1.size() "<<tmphadcluster1.size()<<endl;
	    //----------------------------------------------
	    //Begin cluster algo
	    //
	    //---------------------------------------------
	    
	    int tmphitcnt = 0;
	    // cout <<"tmphitcnt "<< tmphitcnt<<" "<<totalhits<<endl;
	    while (totalhits>tmphitcnt) {
	      //cout<<"total hits before for (int ix=0; ix<totalhits; ix++) loop "<<totalhits<<endl;
	      for (int ix=0; ix<totalhits; ix++) {
		tmphadcluster.clear();
		InoHit* tmphit = tmphadcluster1[ix];//SSE
		//		cout <<"ixxxxx "<<ix<<" "<< tmphit->GetUID()<<endl;
		
		if (tmphit->GetUID()>1) continue;
		tmphadcluster.push_back(tmphit);
		tmphit->SetUID(2);
		tmphitcnt++;
		// cout<<"tmphitcint "<<tmphitcnt<<endl;
		// cout <<"ix "<<ix<<"  tmphit->GetZPlane() "<< tmphit->GetZPlane()<<" "<<tmphit->GetUID()<<" "<<tmphit->GetXStripNum()<<" "<<tmphit->GetYStripNum()<<endl;
	    
		int initialclustersize = 0;
		int finalclustersize = tmphadcluster.size();
		while (initialclustersize != finalclustersize) {
		  initialclustersize = finalclustersize;
		  // cout<< " initial cluster size "<<initialclustersize<<endl;//SSE
		  for (int iy=0; iy<initialclustersize; iy++) {
		    InoHit* tmphit2 = tmphadcluster[iy];
		    //cout <<"iy "<< iy<<" "<<tmphit2->GetZPlane()<<" "<<tmphit2->GetUID()<<" X strp "<<tmphit2->GetXStripNum()<<" YStrp "<<tmphit2->GetYStripNum()<<" XPOS " <<tmphit2->GetXPos()<<" YPOS " <<tmphit2->GetYPos()<<endl;
		    if (tmphit2->GetUID() >2) continue;
		    tmphit2->SetUID(3);
		    
		    for (int iz=0; iz<totalhits; iz++) {
		      //InoHit* tmphit3 = tmphitlist->InoHit_list[iz];////this line is replaced by next line
		      InoHit* tmphit3 = tmphadcluster1[iz];//SSE
		      // cout <<"iz "<< iz<<" "<<tmphit3->GetUID()<<endl;
		      if (tmphit3->GetUID() >1) continue;
		      //cout <<"iz "<< iz<<" "<<tmphit3->GetUID()<<" "<<tmphit2->GetUID()<<" "<<tmphit3->GetXStripNum()<<" "<<tmphit3->GetYStripNum()<<endl;
		      //  cout<<"dif "<< tmphit2->GetZPlane()<<" "<<tmphit3->GetZPlane()<<" "
		      //   <<tmphit2->GetXPos()<<" "<<tmphit3->GetXPos()<<" "
		      //  <<tmphit2->GetYPos()<<" "<<tmphit3->GetYPos()<<endl;
		      if (abs(tmphit2->GetZPlane()- tmphit3->GetZPlane()) < 5 &&
			  abs(tmphit2->GetXPos()- tmphit3->GetXPos())<0.15 &&
			  abs(tmphit2->GetYPos()- tmphit3->GetYPos())<0.15) {// strip distance 15 cm and 5 layers
			tmphadcluster.push_back(tmphit3);
			tmphit3->SetUID(2);
			tmphitcnt++;
			finalclustersize++;
			// cout <<"2: iz "<< iz<<" "<<tmphit3->GetUID()<<" "<<tmphit2->GetUID()<<" "<<tmphit3->GetXStripNum()<<" "<<tmphit3->GetYStripNum()<<endl;
			
			//   cout <<" tmphadcluster size tmphitcnt final ini size "<< tmphadcluster.size()<<" "<< tmphitcnt<<" "<<finalclustersize<<" "<< initialclustersize<<endl;
		      }
		    }
		  }
		  // cout<< finalclustersize <<" "<<initialclustersize<<endl;
		}
		allhadcluster.push_back(tmphadcluster);
		// cout<<"allhadcluster "<<allhadcluster.size()<<endl;
	      }
	      // cout<<"for (int ix=0; ix<totalhits; ix++) out of loop"<<endl;
	    }
	    // cout <<"tmphitcnt2 "<< tmphitcnt<<" "<<totalhits<<endl;
	    //According to the number of hit, sort them out
	    int nclhit=allhadcluster.size();//no. of clusters 
	    // cout<<" total no. of cluster"<<nclhit<<endl;
	    if (nclhit>1000) nclhit=1000;
	    int in_r[1000];
	    int clhits[1000];
	    for (int ix=0; ix<1000;ix++){
	      clhits[ix]=0;//added by SSE, initialisation needed, otherwise gives problem during storing num_hadron		
	    }
	    
	    for (int ix=0; ix<nclhit; ix++) {
	      in_r[ix] = ix;//store cluster number
	      clhits[ix] = allhadcluster[ix].size();//no of hits in "ix th" cluster
	      // cout<< "ix :"<< ix<< " hit in ix th cluster  "<< clhits[ix]<<endl;//SSE
	    }

	    for (int ix=0; ix<nclhit; ix++){
	      for (int iy=ix+1; iy<nclhit; iy++) {
		if (clhits[iy]>clhits[ix]) {
		  int tmpid = in_r[ix];
		  in_r[ix] = in_r[iy];// 0 th element of in_r stores the largest cluster's serial no
		  in_r[iy] = tmpid;
		  
		  tmpid =clhits[ix];
		  clhits[ix] = clhits[iy];
		  clhits[iy] = tmpid;
		  //	cout<< "clhits[ix] "<<clhits[ix]<<" clhits[iy] "<<clhits[iy]<<endl;	
		}
	      }
	    }

	    num_hadron=clhits[0]; 
	    vector<InoHit*> tmphadcluster2;//vector to store hits from largest cluster
	    float cluster_plane[1000];
	    int ii=0;
	    for (unsigned int ixx=0; ixx<allhadcluster.size(); ixx++) {
	      int ix= in_r[ixx];
	      //cout<<"ixxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx "<< ix<<" size  "<<allhadcluster[ix].size()<<endl;
	      for (unsigned int iy=0; iy<allhadcluster[ix].size(); iy++) {
		
		if (ixx==0) {
		  allhadcluster[ix][iy]->SetUID(3333);
		  ii++;
		  tmphadcluster2.push_back(allhadcluster[ix][iy]);
		  cluster_plane[ii]=allhadcluster[ix][iy]->GetZPlane(); 
		}
		//cout<<"ix "<<ix<<" "<< iy<<" "<<allhadcluster[ix][iy]->GetZPlane()<<" "
		//  <<allhadcluster[ix][iy]->GetXPos()<<" "
		// <<allhadcluster[ix][iy]->GetYPos()<<" "
		
		
		// <<allhadcluster[ix][iy]->GetXStripNum()<< " "
		//<<allhadcluster[ix][iy]->GetYStripNum()<<endl;
	      }
	    }

	    //pAnalysis->hist_nhits_LargestCluster_E->Fill(num_hadron, E_had);//SSE

	    //==================================================================
	    // orighit calculation taking the output of cluster algorithm as 
	    // input of orighit algo  
	    //SSE 09/15
	    //==================================================================
	    Orighits_cluster = orighit_calc(tmphadcluster2);
	    
	    //histogram filling for Orighits
	    //pAnalysis->hist_orighits_mod_E->Fill(Orighits_cluster, E_had);
	    //
	    //==================================================================
	    //choosing the largest cluster contains atleast two plane
	    //=================================================================
	    
	    for( unsigned int ijk=0; ijk<tmphadcluster2.size(); ijk++) {
	      if (ijk>0 && cluster_plane[0]!=cluster_plane[ijk]) {
		minplane_cluster++;
		break;
	      }
	    } 
	    //cout<< " mincluster_plane "<<minplane_cluster<<endl;
	    
	    //===================================================================
	    //Algorithm to find direction of hadronic shower
	    //=================================================================== 
	    
	    if (minplane_cluster>0){// lagest cluster must contains min two planes.
	      nhits_large_clust_selected=num_hadron;
	      
	      //Calculation of center of largsest cluster	
	      
	      float center_x=0.0, center_y=0.0, center_z=0.0;
	      
	      for( unsigned int mn=0; mn<tmphadcluster2.size(); mn++) {
		center_x +=tmphadcluster2[mn]->GetXPos();
		center_y +=tmphadcluster2[mn]->GetYPos();
		center_z +=tmphadcluster2[mn]->GetZPos();
	      }
	      center_x= center_x/(tmphadcluster2.size()); 
	      center_y= center_y/(tmphadcluster2.size()); 
	      center_z= center_z/(tmphadcluster2.size()); 
	      //cout<< " center_x "<< center_x<< " center_y "<< center_y<<" center_z "<< center_z<<endl;
	      // and Calculating matrix element to find hadron direction.
	      
	      for( unsigned int mn=0; mn<tmphadcluster2.size(); mn++) {
		
		matrix_XX=(tmphadcluster2[mn]->GetXPos() - center_x);
		matrix_YY=(tmphadcluster2[mn]->GetYPos() - center_y);
		matrix_ZZ=(tmphadcluster2[mn]->GetZPos() - center_z);
		
		//	cout<< " xx  "<< tmphadcluster2[mn]->GetXPos()
		//	<< " yy  "<< tmphadcluster2[mn]->GetYPos()
		//	<< " zz  "<< " "<<tmphadcluster2[mn]->GetZPos()<<endl<<endl;
		
		//	cout<< " matrix_XX "<< matrix_XX<< " matrix_YY "<< matrix_YY<< "  matrix_ZZ "<<  matrix_ZZ<<endl;			 
        
		Matrix_xx +=matrix_XX*matrix_XX;
		Matrix_xy +=matrix_XX*matrix_YY;
		Matrix_yy +=matrix_YY*matrix_YY;
		Matrix_yz +=matrix_YY*matrix_ZZ;
		Matrix_zz +=matrix_ZZ*matrix_ZZ;
		Matrix_zx +=matrix_ZZ*matrix_XX;
		
		//cout<< " Matrix_xx  " << Matrix_xx<< " Matrix_xy "<<Matrix_xy <<
		//	" Matrix_yy " << Matrix_yy<< " Matrix_yz "<<Matrix_yz <<
		//	" Matrix_zz " << Matrix_zz<< " Matrix_zx "<<Matrix_zx << endl;
		
	      }
	      
	      // cout<< " Matrix_xx  " << Matrix_xx<< " Matrix_xy "<<Matrix_xy <<
	      // 	" Matrix_yy " << Matrix_yy<< " Matrix_yz "<<Matrix_yz <<
	      // 	" Matrix_zz " << Matrix_zz<< " Matrix_zx "<<Matrix_zx << endl;
	      
	    } //if (minplane_cluster>0)// SSE 03 Nov 2015
	    
	    //==================================================================================
	    //				Algorithm of ghost hits removal
	    //=================================================================================
	    //create hitbank using hits from largest cluster from cluster algo.
	    //initialise the HitBank_largestClust otherwise some info stored init
	    for (int ixj =0; ixj<nLayer;ixj++){
	      HitBank_largestClust[ixj].clear();
	    }
	    
	    // Filling HIt Bank (needed in GHR algo)
	    for( unsigned int mn=0; mn<tmphadcluster2.size(); mn++) {
	      HitBank_largestClust[tmphadcluster2[mn]->GetZPlane()].push_back(tmphadcluster2[mn]); //filling HitBank
	    }
	    
	    int iFlag_had[500][300];
	    //Initialisation
	    int iXHFill[500][500];// 1st array is for layer, 2nd on n_had (n_had no. of layer have one hit) 
	    int iYHFill[500][500];
	    
	    int diff_old=1000;
	    int diff=1000;
	    int n_had=0,n_had1=0;
	    int layer_prev=-10000;
	    int layer_post=-10000;
	    
	    int count_hit_layer[250];
	    
	    //initilisation of iFlag_had
	    for(int ixj=0; ixj<500; ++ixj) {
	      for(unsigned int jxk=0; jxk<HitBank_largestClust[ixj].size(); ++jxk) {
		iFlag_had[ixj][jxk]=0;				
	      }
	    }
	    //initialisation of iXHFill[ixj][jxk] and iYHFill[ixj][jxk], array of layer and n_had		
	    for(unsigned int ixj=0; ixj<500; ++ixj) {
	      for(int jxk=0; jxk<500; ++jxk) {
		iXHFill[ixj][jxk]=-100; 	 iYHFill[ixj][jxk]=-100; 				
	      }
	    }
	    for (int iter=0;iter<1;iter++){
	      layer_prev=-10000;//initialise previous layer
	      layer_post=-10000;//initialise post layer
	      
	      //To find which layer has only one hit
	      for(int ixj=0; ixj<nLayer; ++ixj) {//loop on layer number 
		count_hit_layer[ixj]=0;
		for(unsigned int jxk=0;jxk<HitBank_largestClust[ixj].size(); ++jxk){
		  if (iFlag_had[ixj][jxk]==0 && HitBank_largestClust[ixj][jxk]->GetUID()==3333){
		    count_hit_layer[ixj]=count_hit_layer[ixj]+1;//count number of hit in a layer
		  }
		}
		if (count_hit_layer[ixj]==1) {
		  // cout<<ij<< " layer has  count "<< count_hit_layer[ij]<<endl;
		}
	      }
	      
	      //Start finding ghost hit 
	      //	Start loop in upward direction 
	      
	      for(int iyj=0; iyj<nLayer; ++iyj) {//loop on layer number
		if((count_hit_layer[iyj]==1 && iyj>=0) || (count_hit_layer[iyj]>1 && iyj==layer_post)) {
		  //if a layer having one hit or connected to the hit of post layer
		    //start the loop on hits of that layer
		    //start connecting to hit in next layer and also identify ghost hits
		    //
		    for(unsigned int jyk=0; jyk<HitBank_largestClust[iyj].size(); ++jyk) {
		      if (count_hit_layer[iyj]==1 && iFlag_had[iyj][jyk]==0) {
			n_had=n_had+1; 
			iFlag_had[iyj][jyk]=(1+10*iter);//Flag changing
			iXHFill[iyj][n_had]= HitBank_largestClust[iyj][jyk]->GetXStripNum() ;
			iYHFill[iyj][n_had]= HitBank_largestClust[iyj][jyk]->GetYStripNum() ;
			//cout<< "layer "<<iyj << " has 1 hit " <<iXHFill[iyj][n_had] << " "<<iYHFill[iyj][n_had] << " n_had  "<< n_had<<endl; 
		      }

		      diff_old=1000;		
		      //cout	<< " start connecting hit (of next layer with)  iXHFill[iyj][n_had]" 
		      //	<<iXHFill[iyj][n_had]<<" "<<iYHFill[iyj][n_had]<<endl;	
	     		
		      //Start the if on connected hit
		      if(HitBank_largestClust[iyj][jyk]->GetXStripNum()==iXHFill[iyj][n_had] &&
			 HitBank_largestClust[iyj][jyk]->GetYStripNum()==iYHFill[iyj][n_had] &&
			 (iFlag_had[iyj][jyk]==0 ||iFlag_had[iyj][jyk]==(1+10*iter))) {

			if (iFlag_had[iyj][jyk]==0) {iFlag_had[iyj][jyk]=(2+10*iter);} //Flag changing
			for(unsigned int kl=0; iyj<nLayer-1 && kl<HitBank_largestClust[iyj+1].size(); ++kl) {
			  //start loop on hit points of consecutive post layer
			  //    cout	<< "next layer: "<< (iyj+1)
			  //		<< " hit no. "<< kl
			  //		<<" X: "<<HitBank_largestClust[iyj+1][kl]->GetXStripNum() << "Y: "<< HitBank_largestClust[iyj+1][kl]->GetYStripNum()<<endl;
			  //if next layer's hit is non-connected then find nearest hit
			  if(iFlag_had[iyj+1][kl]==0) {           
			    diff=abs(HitBank_largestClust[iyj][jyk]->GetXStripNum()-HitBank_largestClust[iyj+1][kl]->GetXStripNum())
			      +abs(HitBank_largestClust[iyj][jyk]->GetYStripNum()-HitBank_largestClust[iyj+1][kl]->GetYStripNum());
			    // 		cout<<" post: diff "<< diff<< " diff _old "<< diff_old<<endl;
			    if(diff<diff_old){
			      diff_old=diff;
			      layer_post=iyj+1;		
			      iXHFill[iyj][n_had]= HitBank_largestClust[iyj][jyk]->GetXStripNum();
			      iYHFill[iyj][n_had]= HitBank_largestClust[iyj][jyk]->GetYStripNum(); 
			      iXHFill[iyj+1][n_had]= HitBank_largestClust[iyj+1][kl]->GetXStripNum() ;
			      iYHFill[iyj+1][n_had]= HitBank_largestClust[iyj+1][kl]->GetYStripNum();
			    }
			  }//end of if(iFlag_had[iyj-1][kl]==0)
			}//end loop on hit points of  post layer
	
			//   cout	<<" *layer "<<iyj<<" n_had "<< n_had<<" "<<iXHFill[iyj][n_had]<<" "<<iYHFill[iyj][n_had]
			//		<<" post layer "<< iyj+1<< " "<<iXHFill[iyj+1][n_had]<< " "<<iYHFill[iyj+1][n_had]<<endl;
		      } else { //end if on connected hit
			int irpcmod =HitBank_largestClust[iyj][jyk]->GetRPCmod();
			int nxstrip=0;
			int nystrip=0;
			for (unsigned int ix=0; ix<inoRPC_pointer->InoRPCStrip.size(); ix++) {
			  if (inoRPC_pointer->InoRPCStrip[ix].first==irpcmod) {
			    nxstrip = inoRPC_pointer->InoRPCStrip[ix].second/100;
			    nystrip = inoRPC_pointer->InoRPCStrip[ix].second%100;
			    break;
			  }
			}
			
			if (nxstrip>1 && nystrip>1 &&
			    ((HitBank_largestClust[iyj][jyk]->GetXStripNum()==iXHFill[iyj][n_had]&&
			      HitBank_largestClust[iyj][jyk]->GetYStripNum()!=iYHFill[iyj][n_had])||
			     (HitBank_largestClust[iyj][jyk]->GetXStripNum()!=iXHFill[iyj][n_had]&& 
			      HitBank_largestClust[iyj][jyk]->GetYStripNum()==iYHFill[iyj][n_had]))&& 
			    iFlag_had[iyj][jyk]==0){//if on unconnected hit
			  iFlag_had[iyj][jyk]=-1-10*iter;//flag realted to ghost hits
			  // cout<<"layer : "<<iyj<< " "<< HitBank_largestClust[iyj][jyk]->GetXStripNum() <<" "<<HitBank_largestClust[iyj][jyk]->GetYStripNum()<< "is ghost hit"<<endl;
			}//endi of else if
			//		}//if(HitBank_largestClust[iyj][jyk]->GetUID()==3333)
		      }
		    }//end of loop on hit bank on the layer which has one hit or connected hit
		}//end of if(HitBank_largestClust[iyj].size()==1 && iyj>0)
	      }//end ....loop on layer no.
	      //cout<<"End of loop on post-layer, no. of layer having one hit " <<n_had<<endl<<endl;


	      //-------------------------------------------------
	      //start the loop on the downwards direction
	      
	      for(int iyj=nLayer-1; iyj>=0; --iyj) {//loop on layer number
		//if((HitBank[iyj].size()==1 && iyj<250)||(HitBank[iyj].size()>1 && iyj==layer_prev)) {//if on layer having one hit
		if((count_hit_layer[iyj]==1 && iyj<nLayer)||(count_hit_layer[iyj]>1 && iyj==layer_prev)) {//if on layer having one hit
		  // cout<<" Enter SEcond loop for dowards layer"<<endl;
		  for(unsigned int jyk=0; jyk<HitBank_largestClust[iyj].size(); ++jyk) {
		    //	if(HitBank_largestClust[iyj][jyk]->GetUID()==3333){		
		    //if (HitBank_largestClust[iyj].size()==1){
		    if (count_hit_layer[iyj]==1 && iFlag_had[iyj][jyk]==1+10*iter) {
		      n_had1=n_had1+1; iFlag_had[iyj][jyk]=1+10*iter;
		      iXHFill[iyj][n_had1]= HitBank_largestClust[iyj][jyk]->GetXStripNum() ;	
		      iYHFill[iyj][n_had1]= HitBank_largestClust[iyj][jyk]->GetYStripNum() ;
		      //  cout<< "layer "<<iyj << " have one hit "<<iXHFill[iyj][n_had1] << " "<<iYHFill[iyj][n_had1] <<endl; 
		    }	
		    //checking previous layer 
		    diff_old=1000;
		    // cout<< "connect hit of prev layer with  iXHFill[iyj][n_had]" <<iXHFill[iyj][n_had1]<<" "<<iYHFill[iyj][n_had1]<<endl;	
		    if (HitBank_largestClust[iyj][jyk]->GetXStripNum()==iXHFill[iyj][n_had1]&& HitBank_largestClust[iyj][jyk]->GetYStripNum()==iYHFill[iyj][n_had1] 
			&& (iFlag_had[iyj][jyk]==0|| iFlag_had[iyj][jyk]==(1+10*iter)||iFlag_had[iyj][jyk]==(2+10*iter))){//if on connected hit
		      if (iFlag_had[iyj][jyk]==0)iFlag_had[iyj][jyk]=3+10*iter;
		      for(unsigned int kl=0; iyj>0 && kl<HitBank_largestClust[iyj-1].size(); ++kl) {
			//  cout<< "prev. layer: "<< iyj-1<< " hit no. "<< kl
			//	  <<" "<<HitBank_largestClust[iyj-1][kl]->GetXStripNum() 
			//	  << " "<< HitBank_largestClust[iyj-1][kl]->GetYStripNum()<<"iFlag"<< iFlag_had[iyj-1][kl]<<endl;
			
			if(iFlag_had[iyj-1][kl]==0){    
			  diff=abs(HitBank_largestClust[iyj][jyk]->GetXStripNum()-HitBank_largestClust[iyj-1][kl]->GetXStripNum())
			    +abs(HitBank_largestClust[iyj][jyk]->GetYStripNum()-HitBank_largestClust[iyj-1][kl]->GetYStripNum());
			  //cout<<" calculate difference "<<endl;
			  //cout<<"pre : diff "<< diff<< "diff _old "<< diff_old<<endl;
			  if(diff<diff_old){
			    diff_old=diff;
			    layer_prev=iyj-1;		
			    iXHFill[iyj][n_had1]=HitBank_largestClust[iyj][jyk]->GetXStripNum() ; 
			    iYHFill[iyj][n_had1]= HitBank_largestClust[iyj][jyk]->GetYStripNum(); 
			    iXHFill[iyj-1][n_had1]=HitBank_largestClust[iyj-1][kl]->GetXStripNum();
			    iYHFill[iyj-1][n_had1]=HitBank_largestClust[iyj-1][kl]->GetYStripNum();
			  }
			} else if (iFlag_had[iyj-1][kl]==1+10*iter || iFlag_had[iyj-1][kl]==2+10*iter) {
			  //if two hits are already connected by upwards case
			  //then take those are nearest
			  //which help to avoid the case same hit of a particular layer
			  //connect two diiferent hits in previous(or succesive layer)
			  //during upwards & downwards direction checking
			  
			  iXHFill[iyj][n_had1]=HitBank_largestClust[iyj][jyk]->GetXStripNum() ; 
			  iYHFill[iyj][n_had1]= HitBank_largestClust[iyj][jyk]->GetYStripNum(); 
			  iXHFill[iyj-1][n_had1]=HitBank_largestClust[iyj-1][kl]->GetXStripNum();
			  iYHFill[iyj-1][n_had1]=HitBank_largestClust[iyj-1][kl]->GetYStripNum();
			  diff_old=-0.1;//arbitrary -ve number so that diff from prev. if never less than that.
			  layer_prev=iyj-1;
			  //cout<<"forced to particular hit"<<endl;
			  //continue;
			  //break;
			}
			
		      }//loop on hit points of pre layer
		      //cout	<<" ***layer "<<iyj<<" n_had1 "<<n_had1 <<" "<<iXHFill[iyj][n_had1]<<" "<<iYHFill[iyj][n_had1]
		      //		<<" prev layer "<< iyj-1<< " "<<iXHFill[iyj-1][n_had1]<< " "<<iYHFill[iyj-1][n_had1]<<endl;
		    } else { //end if on connected hit
		      int irpcmod =HitBank_largestClust[iyj][jyk]->GetRPCmod();
		      int nxstrip=0;
		      int nystrip=0;
		      for (unsigned int ix=0; ix<inoRPC_pointer->InoRPCStrip.size(); ix++) {
			if (inoRPC_pointer->InoRPCStrip[ix].first==irpcmod) {
			  nxstrip = inoRPC_pointer->InoRPCStrip[ix].second/100;
			  nystrip = inoRPC_pointer->InoRPCStrip[ix].second%100;
			  break;
			}
		      }
		      if (nxstrip>1 && nystrip>1 &&
			  ((HitBank_largestClust[iyj][jyk]->GetXStripNum()==iXHFill[iyj][n_had1]&& 
			    HitBank_largestClust[iyj][jyk]->GetYStripNum()!=iYHFill[iyj][n_had1])||
			   (HitBank_largestClust[iyj][jyk]->GetXStripNum()!=iXHFill[iyj][n_had1]&& 
			    HitBank_largestClust[iyj][jyk]->GetYStripNum()==iYHFill[iyj][n_had1])) && 
			  iFlag_had[iyj][jyk]==0){//if on unconnected hit
			iFlag_had[iyj][jyk]=-1-10*iter;
			// cout<< "layer: "<< iyj<< " "<<HitBank_largestClust[iyj][jyk]->GetXStripNum() 
			//     <<" "<<HitBank_largestClust[iyj][jyk]->GetYStripNum()<< "is ghost hit"<<endl;
		      }//end of else if
		    }
		    //	}//if(HitBank_largestClust[iyj][jyk]->GetUID()==3333)
		  }//end of loop on hit bank	
		}//end of if(HitBank[iyj].size()==1 && iyj>0)
		//cout<<"layer "<<iyj<<" "<<iXHFill[iyj][n_had]<<" "<<iYHFill[iyj][n_had]<<" prev "<< iXHFill[iyj-1][n_had]<< " "<<iYHFill[iyj-1][n_had] <<endl;
	      }//end ....loop on layer no.
	      // cout<<" one  in a layer " <<n_had1<<endl;
	    } //for loop on iter //End of ghost removal algorithm
	    
	    int count_h=0;//counter to calculate remaining hit after ghost hit removal algo
	    vector<InoHit*> tmphadcluster3;//Create vector of InoHit with hits remained after removing ghost hits
	    for(int iyj=0; iyj<nLayer; ++iyj) {//loop on layer number 
	      count_hit_layer[iyj]=0;
	      for(unsigned int jyk=0;jyk<HitBank_largestClust[iyj].size(); ++jyk){
		// cout<< "layer "<<iyj<< "   x: "<<HitBank_largestClust[iyj][jyk]->GetXStripNum() <<" y:"<<HitBank_largestClust[iyj][jyk]->GetYStripNum()
		//	  << "    iFlag     "<< iFlag_had[iyj][jyk] <<" GetUID() "<<HitBank_largestClust[iyj][jyk]->GetUID() <<endl;
		
		if (iFlag_had[iyj][jyk]>=0) {
		  HitBank_largestClust[iyj][jyk]->SetUID(4444);
		  tmphadcluster3.push_back(HitBank_largestClust[iyj][jyk]);
		  count_h=count_h+1;
		}
		//cout<< " GetUID() "<<HitBank_largestClust[iyj][jyk]->GetUID() <<endl;
	      }
	    }

	    Hit_wo_ghst=count_h;
	    //cout<< " Hit_wo_ghst "<< count_h<< " "<<tmphadcluster3.size()<<endl;
	    //==========================================================================
	    //			APPLY orighit algorithm after removing ghost hits
	    //==========================================================================
	    wogh_Orighits = orighit_calc(tmphadcluster3);
	    // cout<<" chk: orighit function after ghost hit removal"<<wogh_Orighits<<endl;
	    
	    //cout<< " hit_wo_ghst  "<< " "<<Hit_wo_ghst<<" pAnalysis->momin[0] "<<pAnalysis->momin[0]<< " pAnalysis->momin[1] "<<pAnalysis->momin[1] <<" E_had"<<E_had<<endl;
	    //pAnalysis->hh_woghst_E->Fill(Hit_wo_ghst, E_had);//SSE
	    
	    //histogram filling for Orighits 
	   // pAnalysis->hist_wogh_orighits_E->Fill(wogh_Orighits , E_had);
	     cout<<" Orighits_cluster "<<Orighits_cluster<<"  "<<Hit_wo_ghst<<" wogh_Orighits  "<<wogh_Orighits<<endl;
	    
	  } //else x0= // loop jyk  // else of finding a good recosntructed muon 
	  
	  
	  //=========================================================================
	  //hadron shower direction
	  //calculation of eigen value and eigen vector
	  //==============================================================================
	  double Rx_had=0.0,  Ry_had=0.0,  Rz_had=0.0;
	  float R_had=0.0, theta_had=0.0, phi_had=0.0; //SSE 071015
	  
	  double matrix_elemnt[3][3]={{Matrix_xx,Matrix_xy,Matrix_zx}, {Matrix_xy, Matrix_yy,  Matrix_yz}, {Matrix_zx, Matrix_yz,  Matrix_zz}};
	  //double matrix_elemnt[3][3]={{1,0,0}, {0,0,1}, {0,1,0}}; //for testing
	  TMatrixD mat_dir(3,3);
	  for (int row=0; row<3; row++) { 
	    for (int column=0;column<3;column++){ 
	      mat_dir[row][column] = matrix_elemnt[row][column];
	    } 
	  }
	  
	  // mat_dir.Print();
	  const TMatrixDEigen eigen(mat_dir);
	  //const TVectorD eigenVal = eigen.GetEigenValues();
	  const TMatrixD eigenVal = eigen.GetEigenValues();
	  const TMatrixD eigenVect=eigen.GetEigenVectors();
	  // eigenVal.Print();
	  
	  //cout<< " eigen values " <<eigenVal[0]<< " "<<eigenVal[1] << " " << eigenVal[2] << endl;
	  // eigenVect.Print();
	  
	  //have to check whether sorting required
	  for (int row=0; row<3; row++) {
	    for (int column=0;column<1;column++){
	      if (row==0)Rx_had=eigenVal[0][0]*eigenVect[row][column];
	      if (row==1)Ry_had=eigenVal[0][0]*eigenVect[row][column];
	      if (row==2)Rz_had=eigenVal[0][0]*eigenVect[row][column];
	      //cout<< eigenVect[row][column]<<endl;
	    }
	  }
	  G4ThreeVector tmp3v(Rx_had,Ry_had,Rz_had);
	  R_had = tmp3v.mag();
	  theta_had = tmp3v.theta();
	  phi_had = tmp3v.phi();
	  //costheta_had=tmp3v.cosTheta();
	  dotAngle= tmp3v.cosTheta(tmp3pin_had);
	  //cout<< " Rx_had "<< Rx_had << " Ry_had "<< Ry_had<< " Rz_had "<<Rz_had<<endl;
	  //cout<< " R "<<R_had  <<" had shower theta "<< theta_had<< " had shower phi  "<<phi_had <<endl;
	  //cout<< " had shower costheta "<< costheta_had<< " dotAngle "<< dotAngle<<endl;
	  
	  pAnalysis->e_hadron=E_had;//SSE 250915
	  pAnalysis->orighits_trape=Orighits_trape; //SSE 250915
	  pAnalysis->nhits_largest_cluster=num_hadron;//SSE 250915
	  pAnalysis->orighits_cluster=Orighits_cluster;////SSE 250915
	  pAnalysis->hit_wo_ghst=Hit_wo_ghst;//SSE 250915
	  pAnalysis->hit_wogh_orighits=wogh_Orighits ;//SSE 250915
	  if(theta_had!=0)pAnalysis->theta_hadron_shw=theta_had ;//SSE 071015
	  //pAnalysis->costheta_hadron_shw=costheta_had ;//SSE 071015
	  if(phi_had!=0)pAnalysis->phi_hadron_shw=phi_had ;//SSE 071015
	  if( dotAngle!=0)pAnalysis-> dot_angle_had_shw=dotAngle ;//SSE 081015
	  pAnalysis->nhits_largest_cluster_selected=nhits_large_clust_selected;//SSE 
	  if(theta_had!=0) {
	    for(int eval=0; eval<3; eval++) {
	      pAnalysis->had_eigen_val[eval] = eigenVal[eval][eval];
	    }
	  }
	  pAnalysis->ntotcl = 1000*allremcls + allremhit;
	  pAnalysis->ntotst = 1000*(xtrkstrip.size()-tottrkXstrp) + (ytrkstrip.size()-tottrkYstrp);
	  //pAnalysis->inohits =0;
	  pAnalysis->orighits =0;
	  pAnalysis->x_hits =0;
	  pAnalysis->y_hits =0;
	
	  pAnalysis->inohits =0;
	  int nc_cl=0; int nc_ht=0; int ccc_ht=0 ; int t_ht=0;
	  for (int kl=0; kl<totclustersize; ++kl){
	    
	    if(totcluster[kl]->GetInShower()==0&&totcluster[kl]->GetInTrack()==0){
	      nc_cl++ ; nc_ht += totcluster[kl]->GetHitEntries();
	    } else if(totcluster[kl]->GetInShower()==1){
	      ccc_ht+=totcluster[kl]->GetHitEntries();
	    } else if(totcluster[kl]->GetInTrack()==1) {
	      t_ht+=totcluster[kl]->GetHitEntries();
	    }
	  }
	
	  pAnalysis->inohits =ccc_ht;
	  if(pAnalysis->isXtermOut==2){
	    pAnalysis->ascii_output <<"hadron shower collected: " <<"shower:"<< ccc_ht << " + : track "<<  t_ht <<" + not selected " <<nc_ht <<endl;
	  }
	}

//===========================================
//previous orighit calc
//=====================================

	int Orighits=0;
	int xxnhits =0;
	int yynhits =0;
	// cout<< " xshwstrip.size "<< xshwstrip.size()<< " yshwstrip.size "<< yshwstrip.size() <<endl;	  
	for(int kl=0; kl<nLayer; ++kl){
	  
	  if(xshwstrip.size() && yshwstrip.size()){
	    xxnhits =0;
	    yynhits =0;
	    for(unsigned nn=0; nn< xshwstrip.size(); nn++) {
	      if(xshwstrip[nn].first==int(kl)){
		xxnhits++;
	      }
	    }
	    for (unsigned nn=0; nn< yshwstrip.size(); nn++){
	      if(yshwstrip[nn].first==int(kl)){
		yynhits++;
	      }
	    }
	    //if(showerXHit[kl].size()||showerYHit[kl].size()){
	    //xxnhits += showerXHit[kl].size();
	    //yynhits += showerYHit[kl].size();
	    //pAnalysis->x_hits += xxnhits;
	    //pAnalysis->y_hits += yynhits;
	    Orighits +=max(xxnhits,yynhits);
	  }
	}
	
	// cout<<" orighits "<<Orighits<<endl;
	//histogram filling for Orighits //SSE
	//pAnalysis->hh_E->Fill(Orighits, E_had);//SSE
	//	TProfile* profx =pAnalysis->hh_E->ProfileX("Profile_orighits",-10, 90, "o");//firstybin, lastybin, "o" : original range//SSE
	//		profx->Draw();//SSE
	
	if( abs(pAnalysis->x_hits) <100000) pAnalysis->x_hits=xshwstrip.size();else{pAnalysis->x_hits=-1;}//asm//temporarily commented out
	if( abs(pAnalysis->y_hits) <100000) pAnalysis->y_hits=yshwstrip.size();else{pAnalysis->y_hits=-1;}
	if( abs(pAnalysis->x_hits)<100000 && abs(pAnalysis->y_hits) <100000)pAnalysis->orighits=Orighits;else{pAnalysis->orighits=-1;}

	if(pAnalysis->isXtermOut==2){
	  pAnalysis->ascii_output <<  "(x_hits y_hit orighits inohits)  (" << pAnalysis->x_hits << " " <<pAnalysis->y_hits << " " << pAnalysis->orighits<< " " << pAnalysis->inohits<< ")"<<endl;
	}
      } else { //if (pfitTrack)
	cout <<"XXXXXXXXXXXXXXXXXXXXXXX InoTrackCand_Manager::APointer is not found "<<endl; 
	
      } //loop jk
      pAnalysis->pEventTree->Fill(); //VALGRIND
      int nUp=0;
      int nDown=0;
      
      if (pAnalysis->isVisOut) {
	for (unsigned kl=0; kl< pfitTrack->InoTrackCand_list.size() ; kl++) {
	  pAnalysis->H->NRecTracks++;
	  pAnalysis->Hp= pAnalysis->H->AddHits(0,0); // add a track object/
	  pAnalysis->Hp->TrackType=-140;// Track Type: -1: hits, -2: clulster, -3: triplet, -4: track -14: particle info
	  pAnalysis->Hp->ZZ =((((0.5*(paradef->GetnIRLayer()+1))*LayerThickness)+pAnalysis->poszvx[kl]*cm)/(LayerThickness*m));// vertex z incase of particle ifo
	  
	  pAnalysis->Hp->XX=pAnalysis->posxvx[kl]/100; // vertex x incase of particle info
	  pAnalysis->Hp->YY=pAnalysis->posyvx[kl]/100; // vertex y incase of particle info
	  pAnalysis->Hp->pmag=pAnalysis->trkmm[kl]; // vertex y incase of particle info
	  pAnalysis->Hp->pt=pAnalysis->trkth[kl]; // vertex y incase of particle info
	  pAnalysis->Hp->pp=pAnalysis->trkph[kl]; // vertex y incase of particle info
	  
	  for ( unsigned int lm =0; lm<pfitTrack->InoTrackCand_list[kl]->GetClusterEntries(); lm++){
	    
	    pAnalysis->Hp= pAnalysis->H->AddHits(0,0); // add a track object  //VALGRIND
	    for(unsigned int mn=0;mn<pfitTrack->InoTrackCand_list[kl]->ClustsInTrack[lm]->GetHitEntries();mn++){
	      
	      int zplane = pfitTrack->InoTrackCand_list[kl]->ClustsInTrack[lm]->HitsInCluster[mn]->GetZPlane();
	      if (zplane >=250) {zplane -=250;}
	      pAnalysis->Hp->ZZ=  zplane;//pfitTrack->InoTrackCand_list[kl]->ClustsInTrack[lm]->GetZPlane(); //os();
	      pAnalysis->Hp->XX=pfitTrack->InoTrackCand_list[kl]->ClustsInTrack[lm]->HitsInCluster[mn]->GetXPos();
	      pAnalysis->Hp->YY=pfitTrack->InoTrackCand_list[kl]->ClustsInTrack[lm]->HitsInCluster[mn]->GetYPos();
	      if(pfitTrack->InoTrackCand_list[kl]->GetFitType()==1){
		
		pAnalysis->Hp->TrackType=-440; // track type up -44 + 0 FIT UP
		pAnalysis->Hp->Fup=0;// Finder set up
		pAnalysis->Hp->FitUpNum = nUp;
		pAnalysis->H->NFitUp=nUp+1;
	      } else {
		pAnalysis->Hp->TrackType=-441; // track type down -44 + 1 FIT DOWN
		pAnalysis->Hp->Fup=1;//Finder set down
		pAnalysis->Hp->FitDownNum = nDown;
		pAnalysis->H->NFitDown=nDown+1;
	      }
	    }
	    //--------------------------
	  }
	  if (pfitTrack->InoTrackCand_list[kl]->GetFitType()==1) {
	    nUp++;
	  } else {
	    nDown++;
	  }
	}
      } // if (pAnalysis->isVisOut
    } //if (pinotrack)
  } //if (pAnalysis->InputOutput==0 || pAnalysis->InputOutput==3 || pAnalysis->InputOutput==5) 
  
  if (inoRPC_pointer) { delete inoRPC_pointer; inoRPC_pointer=0;}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
