
//////////////////////////////////////////////////////////////////////////////
//
//
//
//
////////////////////////////////////////////////////////////////////////////////

#include "MultiSimAnalysis.hh"
#include "micalPrimaryGeneratorAction.hh"

MultiSimAnalysis *MultiSimAnalysis::AnPointer;

MultiSimAnalysis::MultiSimAnalysis(G4String sprefix) {
  AnPointer = this;
  
  text_inputFile = sprefix.data();  // file name without .root
  paradef = micalDetectorParameterDef::AnPointer;

  h2dDimX = paradef->GetParino(0);
  h2dDimY = paradef->GetParino(1);
  numberInLA = paradef->GetnLayer();
  cout<<"numberInLA "<<numberInLA<<endl;
  nbinxMag2d = (int)(h2dDimX/10);
  nbinyMag2d = (int)(h2dDimY/10);

  magZmax = 2*(paradef->GetParirlay(2) + paradef->GetParlay(2));
  nbinxMagZ = (int)(magZmax);
  //  isVisOut = true;
  isVisOut	= 0; //1:ascii file, 2:histogram, 3:histogram+Tgraph
  isXtermOut	= 0;
  nmxhit		= 10000;
  nhtcal0		= 0;
  ievt		= 0;
  ihist		= 0;
  ievent		= 0;
  
  ievt_wt=0;
  intxn_id=0;
  
  H	=0;
  Hp	=0;
  
  pVisFile	=0;
  pRootFile	=0;
  inputRootFile=0;
  
  pEventTree	=0;
  inputEventTree=0;
  
  pDepthDose	=0;
  pRadialMomentum=0;
  pEnergyDeposit =0;
  
  rtree=0;

  // double pPosLimXY = 52.5;
  // double pPosLimZ = 2.5;

  double pPosLimXY = 105;
  double pPosLimZ = 20.5;

  pPosX = new TH1F("deltaX", "#Delta X (cm)", 105, -pPosLimXY, pPosLimXY);
  pPosY = new TH1F("deltaY", "#Delta Y (cm)", 105, -pPosLimXY, pPosLimXY);
  pPosZ = new TH1F("deltaZ", "#Delta Z (cm)", 100, -pPosLimZ, pPosLimZ);
  
  char title[100];
  for (int ij=0; ij<20; ij++) {
    sprintf(title, "dedz_%i", ij);
    pdedz[ij] = new TH1F(title, title, 4000, -600., 600.);
  }
  
  pPosXX = new TH2F("deltaXX", "#Delta XX (cm)", 100, -2500, 2500, 100, -pPosLimXY, pPosLimXY);
  pPosYY = new TH2F("deltaYY", "#Delta YY (cm)", 100, -500, 500, 100, -pPosLimXY, pPosLimXY);
  pPosZZ = new TH2F("deltaZZ", "#Delta ZZ (cm)", 100, -1000, 1000, 100, -pPosLimZ, pPosLimZ);
  
  
  hitDist  = new TH1F("HitDist","This is the total rawhit distribution across layers",182,0,182);     //asm
  TrkDist  = new TH1F("TrkDist","This is the total track-hit distribution across layers",182,0,182);     //asm
  EffDist  = new TH1F("EffDist","This is the efficiency distribution across layers",182,0,182);     //asm
  
  InoTrack_listsize =new TH1F("TrksizeDist","Distribution of tracks per event ",10,0,10); //asm
  
  ShwXw	= new TH1F("ShwXw ","This is a distribution for X position resolution in shower region " ,80 , -3.14 , 3.14);  //asm
  ShwYw	= new TH1F("ShwYw ","This is a distribution for Y position resolution in shower region " ,80 , -3.14 , 3.14);  //asm
  RC		= new TH2F("RCvsff","This is a distribution for raduis of curvature at the edges " ,1000 , 0 , 1000,100,-50,50);  //asm
  DGap	= new TH1F("DetGap"," Plot to see the performance of code at detgap ",10,0,10);  //asm
  //  nhtcal1 = 0;
  //  nhtcal2 = 0;
  DeadStripX = new TH1F("DeadStripX","Plot to see dead strips on X-plane",64,-0.5,63.5);
  NoisyStripX = new TH1F("NoisyStripX","Plot to see noisy strips on X-Plane",64,-0.5,63.5);
  DeadStripY = new TH1F("DeadStripY","Plot to see dead strips on Y-Plane",64,-0.5,63.5);
  NoisyStripY = new TH1F("NoisyStripY","Plot to see noisy strips on Y-Plane",64,-0.5,63.5);
  DiffTime = new TH1F("DiffTime","DiffTime",100,-100.,100.);
  
  strpXtime = new TH1F("strpXtime","xsmrtime - xtruetime",100,-40,40.);
  strpYtime = new TH1F("strpYtime","ysmrtime - ytruetime",100,-40,40.);
  strpXtimeCorr = new TH1F("strpXtimeCorr","xsmrtimeCorr - xtruetime",100,-20,20.);
  strpYtimeCorr = new TH1F("strpYtimeCorr","ysmrtimeCorr - ytruetime",100,-20,20.);
  hitXtime = new TH1F("hitXtime","hittime - xtruetime",100,-20,20.);
  hitYtime = new TH1F("hitYtime","hittime - ytruetime",100,-20,20.);

  smagFieldX = new TH1D("smagfieldinx","magfieldinx simu",nbinxMagZ,0,magZmax);
  smagFieldY = new TH1D("smagfieldiny","magfieldiny simu",nbinxMagZ,0,magZmax);
  smag2dX = new TH2D("smag2dinx","magfieldinx 2D simu",nbinxMagZ,0,magZmax,180,-2.0,2.0);
  smag2dY = new TH2D("smag2diny","magfieldiny 2D simu",nbinxMagZ,0,magZmax,180,-2.0,2.0);

  rmagFieldX = new TH1D("rmagfieldinx","magfieldinx in reco",nbinxMagZ,0,magZmax);
  rmagFieldY = new TH1D("rmagfieldiny","magfieldiny in reco",nbinxMagZ,0,magZmax);
  rmag2dX = new TH2D("rmag2dinx","magfieldinx 2D in reco",nbinxMagZ,0,magZmax,180,-2.0,2.0);
  rmag2dY = new TH2D("rmag2diny","magfieldiny 2D in reco",nbinxMagZ,0,magZmax,180,-2.0,2.0);
  // int nbin2d = 1600;
  smag2dXYpixel_iron = new TH2D("smag2dxypixeliron","magfield in xy pixel in sim in iron",nbinxMag2d,-h2dDimX,h2dDimX,nbinyMag2d,-h2dDimY,h2dDimY);
  smag2dXYpixel_air = new TH2D("smag2dxypixelair","magfield in xy pixel in sim in air",nbinxMag2d,-h2dDimX,h2dDimX,nbinyMag2d,-h2dDimY,h2dDimY);
  rmag2dXYpixel_iron = new TH2D("rmag2dxypixeliron","magfield in xy pixel in reco in iron",nbinxMag2d,-h2dDimX,h2dDimX,nbinyMag2d,-h2dDimY,h2dDimY);
  rmag2dXYpixel_air = new TH2D("rmag2dxypixelair","magfield in xy pixel in reco in air",nbinxMag2d,-h2dDimX,h2dDimX,nbinyMag2d,-h2dDimY,h2dDimY);

  xyvsbxin = new TH2D("xyvsbxin", "xyvsbxin", 321, -8025, 8025, 321, -8025, 8025);
  xyvsbyin = new TH2D("xyvsbyin", "xyvsbyin", 321, -8025, 8025, 321, -8025, 8025);
  xyvsbxdiff = new TH2D("xyvsbxdiff", "xyvsbxdiff", 321, -8025, 8025, 321, -8025, 8025);
  xyvsbydiff = new TH2D("xyvsbydiff", "xyvsbydiff", 321, -8025, 8025, 321, -8025, 8025);
  // xyvsbxindiff = new TH2D("xyvsbxindiff", "xyvsbxindiff", 319, -7975, 7975, 319, -7975, 7975);
  // xyvsbyindiff = new TH2D("xyvsbyindiff", "xyvsbyindiff", 319, -7975, 7975, 319, -7975, 7975);
  xyvsbxout = new TH2D("xyvsbxout", "xyvsbxout", 81, -2025, 2025, 81, -2025, 2025);
  xyvsbyout = new TH2D("xyvsbyout", "xyvsbyout", 81, -2025, 2025, 81, -2025, 2025);
  // xyvsbxout = new TH2D("xyvsbxout", "xyvsbxout", 969, -24225, 24225, 321, -8025, 8025);
  // xyvsbyout =  new TH2D("xyvsbyout", "xyvsbyout", 969, -24225, 24225, 321, -8025, 8025);

  //G4int hit_wo_ghst;//SSE
  //hprof  = new TProfile("hprof","Profile",2,0,100,0.,1000.);//SSE
 // hh_E 	= new TH2D("hh_E","E_had vs orighits;orighits;E_had (GeV)",25, 0. ,50.,100, 0., 100.);  //SSE
  //hh_woghst_E 	= new TH2D("hh_woghst_E","E_had vs hits;hits;E_had (GeV)",50, 0. ,500.,100, 0., 100.);  //SSE
  //hist_nhits_LargestCluster_E 	= new TH2D("hist_nhits_LargestCluster_E","E_had vs no of hits in largest cluster;hits;E_had (GeV)",50, 0. ,500.,100, 0., 100.);  //SSE
  //hist_orighits_new_E= new TH2D("hist_orighits_new_E","E_had vs orighits (new) ;orighits_new;E_had (GeV)",50, 0. ,100.,100, 0., 100.);//SSE 09/15
 //hist_orighits_mod_E= new TH2D("hist_orighits_mod_E","E_had vs orighits (modified) ;orighits_mod;E_had (GeV)",25, 0. ,50.,100, 0., 100.);//SSE 09/15
 //hist_wogh_orighits_E= new TH2D("hist_wogh_orighits_E","E_had vs wogh_orighits ;Ghosthit_removed_orighit;E_had (GeV)",25, 0. ,50.,100, 0., 100.);//SSE 09/15


}

void MultiSimAnalysis::OpenRootfiles(G4String infile,  G4String outfile, G4String collFile) {
  // cout<<"void MultiSimAnalysis::OpenRootfiles(G4String infile,  G4String outfile) {..."<<endl;
  cout<<"my infile is "<<infile<<endl;
  cout<<"my collated is "<<collFile<<endl;
  //  inputRootFile = new TFile(infile, "read");
  
  //  G4String fname(sprefix.data());  // file name without .root
  ievent=0;
  G4String titleinh= outfile;//text_inputFile;
  G4String titledat= outfile;
  G4String titlefld= outfile;
  G4String titleop= outfile;
  //NOT defined yet
  // G4String titleTimdat = outfile;

  // titleTimdat.append("time.dat");
  // timeAsciiOutput.open(titleTimdat);

  if (isXtermOut ==2 && (InputOutput==1 || InputOutput==2|| InputOutput==3)) {
    titlefld.append("B_field.dat");
    B_ascii_output.open(titlefld);
  }
  if ((InputOutput==0 ||InputOutput==3 || InputOutput==5)) {
    if (isXtermOut==2) {
      titledat.append("_ascii.dat");
      ascii_output.open(titledat);
    }
    if(isVisOut==1) {
      titleinh.append(".inh");
      pVisFile = new TFile(titleinh, "RECREATE"); //VALGRIND
      if (!pVisFile) {
	G4cout << "Error opening .inh file !" << G4endl;
	exit(-1);
      }
    }
  }

  if(collatedIn) {
    collatedRootFile = new TFile(collFile, "read");
    if(!collatedRootFile) {
      cout << "Error opening collated file !" << endl;
      exit(-1);
    }

    /* Please check this size of the arrays in .hh file */
    char namex[200];
    // for(int iki=0; iki<numberInLA; iki++) {
    //   sprintf(namex,"inefficiency_cor_m0_xr0_yr0_l%i",iki);
    //   inefficiency_corx[iki] = (TH2D*)collatedRootFile->Get(namex);
    //   sprintf(namex,"inefficiency_unc_m0_xr0_yr0_x%i",iki);
    //   inefficiency_uncx[iki] = (TH2D*)collatedRootFile->Get(namex);
    //   sprintf(namex,"inefficiency_unc_m0_xr0_yr0_y%i",iki);
    //   inefficiency_uncy[iki] = (TH2D*)collatedRootFile->Get(namex);
    //   sprintf(namex,"triggereffi_evt_m0_xr0_yr0_x%i",iki);
    //   triggereffi_xevt[iki] = (TH2D*)collatedRootFile->Get(namex);
    //   sprintf(namex,"triggereffi_evt_m0_xr0_yr0_y%i",iki);
    //   triggereffi_yevt[iki] = (TH2D*)collatedRootFile->Get(namex);
    //   // sprintf(namex,"strp_xmulsim_cor_l%i",iki);
    //   // strp_xmulsim_cor[iki] = (TH2D*)collatedRootFile->Get(namex);
    //   // sprintf(namex,"strp_ymulsim_cor_l%i",iki);
    //   // strp_ymulsim_cor[iki] = (TH2D*)collatedRootFile->Get(namex);
    // }
    
    // for(int iki=0; iki<numberInLA; iki++) {
    //   cout<<"rading mul hist "<<iki<<endl;
    //   for(int ikj=0; ikj<64; ikj++) {
    // 	for(int ikk=0; ikk<64; ikk++) {
    // 	  sprintf(namex,"strp_xmullaysim_m0_xr0_yr0_l%i_%i_%i",iki,ikj,ikk);
    // 	  // block_xmulsim[iki][ikj][ikk] = (TH2D*)((TH2D*)collatedRootFile->Get(namex))->Clone(namex);
    // 	  block_xmulsim[iki][ikj][ikk] = (TH2D*)collatedRootFile->Get(namex);
    // 	  sprintf(namex,"strp_ymullaysim_m0_xr0_yr0_l%i_%i_%i",iki,ikj,ikk);
    // 	  // block_ymulsim[iki][ikj][ikk] = (TH2D*)((TH2D*)collatedRootFile->Get(namex))->Clone(namex);
    // 	  block_ymulsim[iki][ikj][ikk] = (TH2D*)collatedRootFile->Get(namex);
    // 	}
    //   }
    // } // for(int iki=0; iki<numberInLA; iki++) {
    
    // for(int iki=0; iki<numberInLA; iki++) {
    //   for(int ikj=0; ikj<16; ikj++) {
    // 	for(int ikk=0; ikk<16; ikk++) {
    // 	  sprintf(namex,"blk_xmullaysim_m0_xr0_yr0_l%i_%i_%i",iki,ikj,ikk);
    // 	  block_xmulsim[iki][ikj][ikk] = (TH2D*)collatedRootFile->Get(namex);
    // 	  sprintf(namex,"blk_ymullaysim_m0_xr0_yr0_l%i_%i_%i",iki,ikj,ikk);
    // 	  block_ymulsim[iki][ikj][ikk] = (TH2D*)collatedRootFile->Get(namex);
    // 	}
    //   }
    // } // for(int iki=0; iki<numberInLA; iki++) {

    if(InputOutput==1 || InputOutput==3 || InputOutput==4) {
      cout<<" reading strp sensitivity "<<endl;
      system("free");
      // for(int iki=0; iki<numberInLA; iki++) {
      //   for(int nj=0;nj<2;nj++) {
      // 	sprintf(namex,"laymul_2D_m0_xr0_yr0_%s%i",nj==0?"x":"y",iki);
      // 	// laymul_2D[nj][iki] = (TH3D*)((TH3D*)collatedRootFile->Get(namex))->Clone();
      // 	laymul_2D[nj][iki] = (TH3S*)collatedRootFile->Get(namex);
      // 	system("free");
      //   }} // for(int iki=0; iki<numberInLA; iki++) {
      for(int iki=0; iki<numberInLA; iki++) {
	for(int nj=0;nj<2;nj++) {
	  // sprintf(namex,"laymul_2Dpos_m0_xr0_yr0_%s%i",nj==0?"x":"y",iki);
	  // laymul_2D[nj][iki] = (TH3S*)collatedRootFile->Get(namex);
	  for(int ns=0;ns<64;ns++) {
	    sprintf(namex,"laymul_2Dstrp_m0_xr0_yr0_%s%i_s%i",nj==0?"x":"y",iki,ns);
	    // laymul_2Dstrp_m0_xr0_yr0_y9_s63
	    laymul_2Dstrp[nj][iki][ns] = (TH2S*)collatedRootFile->Get(namex);
	    // system("free");
	  }
	}} // for(int iki=0; iki<numberInLA; iki++) {
      cout<<" reading strp sensitivity : end "<<endl;
      system("free");
    } // if(InputOutput==1 || InputOutput==3 || InputOutput==4) {
    
    // fnMul[0] = new TF1("fn2",
    // 		       "[0]+exp([1]+[2]*pow(x,2))",
    // 		       -0.5,0.5);
    // fnMul[1] = new TF1("fn2inv",
    // 		       "[0]+exp([1]+[2]*pow(0.5-fabs(x),2))",
    // 		       -0.5,0.5);    
    // for(int iki=0;iki<nparMul;iki++) {
    //   sprintf(namex,"proportion_vs_mul_par%i",iki);
    //   profoutMul[iki] = (TH2D*)collatedRootFile->Get(namex);}
    // for(int iki=0; iki<numberInLA; iki++) {
    //   sprintf(namex,"layblocktotmul_m0_xr0_yr0_x%i",iki);
    //   layblocktotmul[iki][0] = (TH3D*)collatedRootFile->Get(namex);
    //   sprintf(namex,"layblocktotmul_m0_xr0_yr0_y%i",iki);
    //   layblocktotmul[iki][1] = (TH3D*)collatedRootFile->Get(namex);}
    
    
    // for(int iki=0; iki<numberInLA; iki++) {
    //   sprintf(namex,"inefficiency_corx_l%i",iki);
    //   inefficiency_corx[iki] = (TH2D*)collatedRootFile->Get(namex);
    //   sprintf(namex,"inefficiency_uncx_l%i",iki);
    //   inefficiency_uncx[iki] = (TH2D*)collatedRootFile->Get(namex);
    //   cout<<"inefficiency_uncx["<<iki<<"] "<<inefficiency_uncx[iki]<<endl;
    //   sprintf(namex,"inefficiency_uncy_l%i",iki);
    //   inefficiency_uncy[iki] = (TH2D*)collatedRootFile->Get(namex);
    //   cout<<"inefficiency_uncy["<<iki<<"] "<<inefficiency_uncy[iki]<<endl;
    //   sprintf(namex,"triggereffi_xevt_l%i",iki);
    //   triggereffi_xevt[iki] = (TH2D*)collatedRootFile->Get(namex);
    //   sprintf(namex,"triggereffi_yevt_l%i",iki);
    //   triggereffi_yevt[iki] = (TH2D*)collatedRootFile->Get(namex);
    //   sprintf(namex,"strp_xmulsim_cor_l%i",iki);
    //   strp_xmulsim_cor[iki] = (TH2D*)collatedRootFile->Get(namex);
    //   sprintf(namex,"strp_ymulsim_cor_l%i",iki);
    //   strp_ymulsim_cor[iki] = (TH2D*)collatedRootFile->Get(namex);
    // }
    // for(int iki=0; iki<numberInLA; iki++) {
    //   for(int ikj=0; ikj<16; ikj++) {
    // 	for(int ikk=0; ikk<16; ikk++) {
    // 	  sprintf(namex,"strp_xmullaysim_l%i_%i_%i",iki,ikj,ikk);
    // 	  block_xmulsim[iki][ikj][ikk] = (TH2D*)collatedRootFile->Get(namex);
    // 	  sprintf(namex,"strp_ymullaysim_l%i_%i_%i",iki,ikj,ikk);
    // 	  block_ymulsim[iki][ikj][ikk] = (TH2D*)collatedRootFile->Get(namex);
    // 	}
    //   }
    // } // for(int iki=0; iki<numberInLA; iki++) {
    
  } 
  
  micalPrimaryGeneratorAction *pgPointer = micalPrimaryGeneratorAction::AnPointer;
  
  if(pgPointer->InputFlag>1) {
    if((InputOutput==0 || InputOutput==1|| InputOutput==2) &&
       (pgPointer->InputFlag<=3)) {
      pgPointer->OpenFileCORSIKA();
    }
  }
  if(pgPointer->InputFlag==4 || pgPointer->InputFlag==5) {
    pgPointer->OpenFileFLUX();
  }

  if(pgPointer->InputFlag==1) {
    if( InputOutput==0 || InputOutput==1|| InputOutput==2) {
      // pgPointer->OpenNuanceFile(infile);
    }
  } else if(InputOutput==3 || InputOutput==4) {
    //Input is simulation file
    inputRootFile = new TFile(infile, "read");
    cout<< "Data is read from simulation output rootfile :"<< infile <<endl;
    
    // inputRootFile = new TFile("simulation.root", "read");
    inputEventTree = (TTree*)inputRootFile->Get("T3");
    
    inputEventTree->SetBranchAddress("nsimht", &nsimht);
    inputEventTree->SetBranchAddress("detid", detid);
    inputEventTree->SetBranchAddress("simpdgid", simpdgid);
    inputEventTree->SetBranchAddress("simtime", simtime);
    inputEventTree->SetBranchAddress("simenr", simenr);
    inputEventTree->SetBranchAddress("simvx", simvx);
    inputEventTree->SetBranchAddress("simvy", simvy);
    inputEventTree->SetBranchAddress("simvz", simvz);
    inputEventTree->SetBranchAddress("simpx", simpx);
    inputEventTree->SetBranchAddress("simpy", simpy);
    inputEventTree->SetBranchAddress("simpz", simpz);
    inputEventTree->SetBranchAddress("simlocvx", simlocvx);
    inputEventTree->SetBranchAddress("simlocvy", simlocvy);
  } else if (InputOutput==5) {
    // Input is digitisation file
    // inputRootFile = new TFile("digitisation.root", "read");
    inputRootFile = new TFile(infile, "read");
    cout<< "Data is read from digitization output file"<< infile <<endl;
    
    inputEventTree = (TTree*)inputRootFile->Get("T2");
    
    inputEventTree->SetBranchAddress("ndigiht", &ndigiht);
    inputEventTree->SetBranchAddress("stripid", stripid);
    inputEventTree->SetBranchAddress("digipdgid", digipdgid);
    inputEventTree->SetBranchAddress("digitime", digitime);
    inputEventTree->SetBranchAddress("digitruetime", digitruetime);
    inputEventTree->SetBranchAddress("digienr", digienr);
    inputEventTree->SetBranchAddress("digivx", digivx);
    inputEventTree->SetBranchAddress("digivy", digivy);
    inputEventTree->SetBranchAddress("digivz", digivz);
    inputEventTree->SetBranchAddress("digipx", digipx);
    inputEventTree->SetBranchAddress("digipy", digipy);
    inputEventTree->SetBranchAddress("digipz", digipz);
  } if (InputOutput==3 || InputOutput==4 || InputOutput==5) {
    inputEventTree->SetBranchAddress("irun",&irun);
    inputEventTree->SetBranchAddress("ievt",&ievt);
    
    inputEventTree->SetBranchAddress("ngent",&ngent);
    inputEventTree->SetBranchAddress("pidin",pidin);
    inputEventTree->SetBranchAddress("momin",momin);
    inputEventTree->SetBranchAddress("thein",thein);
    inputEventTree->SetBranchAddress("phiin",phiin);
    inputEventTree->SetBranchAddress("posxin",posxin);
    inputEventTree->SetBranchAddress("posyin",posyin);
    inputEventTree->SetBranchAddress("poszin",poszin); 
  }

  titleop.append(".root");
  //cout <<"title top=========================="<<titleop<<endl;
  pRootFile = new TFile(titleop, "RECREATE"); //VALGRIND
  
  if (!pRootFile) {
    G4cout << "Error opening histogram file !" << G4endl;
    exit(-1);
  } else {
    cout<< "Output stored in root file:"<< titleop <<endl;
  }
  
  if(!pPosX)pPosX = new TH1F("deltaX", "#Delta X (cm)", 100, -52.5, 52.5);
  //else  pPosX->Clear();
  if(!pPosY)pPosY = new TH1F("deltaY", "#Delta Y (cm)", 100, -52.5, 52.5);
  //else pPosY->Clear();
  if(!pPosZ)pPosZ = new TH1F("deltaZ", "#Delta Z (cm)", 100, -2.5, 2.5);
  //else pPosZ->Clear();
  char title[100];
  for (int i=0; i<20; i++) {
    sprintf(title, "dedz_%i", i);
    if(!pdedz[i]) {
      pdedz[i] = new TH1F(title, title, 4000, -600., 600.);
    }
    //else pdedz[i]->Clear();
  }
  if(!pPosXX)pPosXX = new TH2F("deltaXX", "#Delta XX (cm)", 100, -2500, 2500, 100, -52.5, 52.5);
  //else pPosXX->Clear();
  if(!pPosYY)pPosYY = new TH2F("deltaYY", "#Delta YY (cm)", 100, -500, 500, 100, -52.5, 52.5);
  //else pPosYY->Clear();
  if(!pPosZZ)pPosZZ = new TH2F("deltaZZ", "#Delta ZZ (cm)", 100, -1000, 1000, 100, -2.5, 2.5);
  //else pPosZZ->Clear();
  
  if(!hitDist)  hitDist  = new TH1F("HitDist","This is the total rawhit distribution across layers",182,0,182);     //asm
  // else hitDist->Clear();    //asm
  if(!TrkDist)TrkDist  = new TH1F("TrkDist","This is the total track-hit distribution across layers",182,0,182);     //asm
  // else   TrkDist->Clear();     //asm
  if(!EffDist)EffDist  = new TH1F("EffDist","This is the efficiency distribution across layers",182,0,182);     //asm
  //else  EffDist->Clear();     //asm
  if(!InoTrack_listsize)InoTrack_listsize =new TH1F("TrksizeDist","Distribution of tracks per event ",10,0,10); //asm
  //else  InoTrack_listsize->Clear(); //asm
  
  if(!ShwXw)ShwXw    = new TH1F("ShwXw ","This is a distribution for X position resolution in shower region " ,80 , -3.14 , 3.14);  //asm
  //else ShwXw->Clear(); //asm
  if(!ShwYw)ShwYw    = new TH1F("ShwYw ","This is a distribution for Y position resolution in shower region " ,80 , -3.14 , 3.14);  //asm
  // else ShwYw->Clear();  //asm
  if(!RC)RC       = new TH2F("RCvsff","This is a distribution for raduis of curvature at the edges " ,1000 , 0 , 1000,100,-50,50);  //asm
  //else RC->Clear();  //asm
  if(!DGap) DGap      = new TH1F("DetGap"," Plot to see the performance of code at detgap ",10,0,10);  //asm
  //else DGap->Clear();  //asm

//if (!hh_E)hh_E   = new TH2D("hh_E","E_had vs orighits;orighits;E_had (GeV)",25 , 0. , 50.,100, 0., 100.);  //SSE
 //if(!hh_woghst_E)  hh_woghst_E= new TH2D("hh_woghst_E","E_had vs hits;hits;E_had (GeV)",50, 0. ,500.,100, -10., 90.);  //SSE
 //if(!hist_nhits_LargestCluster_E) hist_nhits_LargestCluster_E= new TH2D("hist_nhits_LargestCluster_E","E_had vs no of hits in largest cluster;hits;E_had (GeV)",50, 0. ,500.,100, -10., 90.);  //SSE
//if(!hist_orighits_new_E)hist_orighits_new_E= new TH2D("hist_orighits_new_E","E_had vs orighits (new) ;orighits_new;E_had (GeV)",50, 0. ,100.,100, 0., 100.);//SSE 09/15
  //if(!hist_orighits_mod_E)hist_orighits_mod_E= new TH2D("hist_orighits_mod_E","E_had vs orighits (modified) ;orighits_mod;E_had (GeV)",25, 0. ,50.,100, 0., 100.);//SSE
 //if(!hist_wogh_orighits_E)hist_wogh_orighits_E= new TH2D("hist_wogh_orighits_E","E_had vs wogh_orighits;Ghosthit_removed_orighit;E_had (GeV)",25, 0. ,50.,100, 0., 100.);//SSE 09/15
//if(!hprof)hprof =new TProfile("hprof","Profile",2,0,100,0.,100.);//SSE 

  



  if (InputOutput==0 || InputOutput==3 || InputOutput==5) {
    //Reco output
    pEventTree = new TTree("T1", "INORECO");
  }  else if (InputOutput==1 || InputOutput==4) {
    //DigiOutput
    pEventTree = new TTree("T2", "INODIGI");
  } else {
    //SImulation output
    pEventTree = new TTree("T3", "INOSIM");
  }
  
  pEventTree->Branch("irun",&irun,"irun/i"); //VALGRIND
  pEventTree->Branch("ievt",&ievt,"ievt/i");
  
  pEventTree->Branch("ngent",&ngent,"ngent/i");
  pEventTree->Branch("pidin",pidin,"pidin[ngent]/I");
  pEventTree->Branch("ievt_wt",&ievt_wt,"ievt_wt/F");		// GMa
  pEventTree->Branch("intxn_id",&intxn_id,"intxn_id/I");		// GMa
  pEventTree->Branch("momin",momin,"momin[ngent]/F"); //VALGRIND
  pEventTree->Branch("thein",thein,"thein[ngent]/F");
  pEventTree->Branch("phiin",phiin,"phiin[ngent]/F");
  pEventTree->Branch("posxin",posxin,"posxin[ngent]/F");
  pEventTree->Branch("posyin",posyin,"posyin[ngent]/F");
  pEventTree->Branch("poszin",poszin,"poszin[ngent]/F");
  pEventTree->Branch("ngenerated",&ngenerated,"ngenerated/i");
  pEventTree->Branch("naperture",&naperture,"naperture/i");
  
  
  if (InputOutput==0 || InputOutput==3 || InputOutput==5) {
    //Reco output
    if (InputOutput==0)	{
      pEventTree->Branch("caltot0",&caltot0,"caltot0/F");
      // pEventTree->Branch("caltot1",&caltot1,"caltot1/F");
      // pEventTree->Branch("caltot2",&caltot2,"caltot2/F");
      
      pEventTree->Branch("nhtcal0",&nhtcal0,"nhtcal0/I");
      pEventTree->Branch("calid0",calid0,"calid0[nhtcal0]/i");
      pEventTree->Branch("calen0",calen0,"calen0[nhtcal0]/F");
    }
    pEventTree->Branch("ntotcl",&ntotcl,"ntotcl/I");
    pEventTree->Branch("ntotst",&ntotst,"ntotst/I");  
    pEventTree->Branch("inohits",&inohits,"inohits/I");
    pEventTree->Branch("orighits",&orighits,"orighits/I");
    pEventTree->Branch("inoclust",&inoclust,"inoclust/I");
    pEventTree->Branch("origclust",&origclust,"origclust/I");
    
    pEventTree->Branch("hPathlength",&hPathlength,"hPathlength/F");
    pEventTree->Branch("x_hits",&x_hits,"x_hits/I");
    pEventTree->Branch("y_hits",&y_hits,"y_hits/I");
    // pEventTree->Branch("nhtcal1",&nhtcal1,"nhtcal1/I");
    pEventTree->Branch("inohits_old",&inohits_old,"inohits_old/I");
    pEventTree->Branch("orighits_old",&orighits_old,"orighits_old/I");
    pEventTree->Branch("x_hits_old",&x_hits_old,"x_hits_old/I");
    pEventTree->Branch("y_hits_old",&y_hits_old,"y_hits_old/I");
    
    //pEventTree->Branch("total_inohits",&total_inohits,"total_inohits/I");//SSE
    pEventTree->Branch("hit_wo_ghst",&hit_wo_ghst,"hit_wo_ghst/I");//SSE
    pEventTree->Branch("e_hadron",&e_hadron,"e_hadron/F");//SSE
    pEventTree->Branch("nhits_largest_cluster",&nhits_largest_cluster,"nhits_largest_cluster/I");//SSE
    pEventTree->Branch("orighits_trape",&orighits_trape,"orighits_trape/I");//SSE 09/15
    pEventTree->Branch("orighits_cluster",&orighits_cluster,"orighits_cluster/I");//SSE 09/15
    pEventTree->Branch("hit_wogh_orighits",&hit_wogh_orighits,"hit_wogh_orighits/I");//SSE 09/15
    pEventTree->Branch("theta_hadron_shw",&theta_hadron_shw,"theta_hadron_shw/F");// 07Oct15 SSE
    pEventTree->Branch("had_eigen_val",had_eigen_val,"had_eigen_val[3]/F");
    //pEventTree->Branch("costheta_hadron_shw",&costheta_hadron_shw,"costheta_hadron_shw/F");// 08Oct15 SSE
    pEventTree->Branch("phi_hadron_shw",&phi_hadron_shw,"phi_hadron_shw/F");// 07Oct15 SSE
    pEventTree->Branch("theta_hadron_in",&theta_hadron_in,"theta_hadron_in/F");// 08Oct15 SSE
    //pEventTree->Branch("costheta_hadron_in",&costheta_hadron_in,"costheta_hadron_in/F");// 08Oct15 SSE
    pEventTree->Branch("phi_hadron_in",&phi_hadron_in,"phi_hadron_in/F");// 08Oct15 SSE
    pEventTree->Branch("dot_angle_had_shw",&dot_angle_had_shw,"dot_angle_had_shw/F");// 08Oct15 SSE
    pEventTree->Branch("nhits_largest_cluster_selected",&nhits_largest_cluster_selected,"nhits_largest_cluster_selected/I");//SSE


    //  pEventTree->Branch("nhtcal1",&nhtcal1,"nhtcal1/I");
    //  pEventTree->Branch("calid1",calid1,"calid2[nhtcal1]/I");
    //  pEventTree->Branch("calen1",calen1,"calen2[nhtcal1]/F");
    
    //  pEventTree->Branch("nhtcal2",&nhtcal2,"nhtcal2/I");
    //  pEventTree->Branch("calid2",calid2,"calid2[nhtcal2]/I");
    //  pEventTree->Branch("calen2",calen2,"calen2[nhtcal2]/F");
    
    pEventTree->Branch("ntrkt",&ntrkt,"ntrkt/i");
    pEventTree->Branch("itype",itype,"itype[ntrkt]/I");
    pEventTree->Branch("nLayer",&nLayer,"nLayer/I");
    pEventTree->Branch("nhits",nhits,"nhits[ntrkt]/I");
    pEventTree->Branch("nhits_finder",nhits_finder,"nhits_finder[ntrkt]/I");
    pEventTree->Branch("chisq",chisq,"chisq[ntrkt]/F");
    pEventTree->Branch("cvalue",cvalue,"cvalue[ntrkt]/F");
    pEventTree->Branch("fc_or_pc",fc_or_pc,"fc_or_pc[ntrkt]/I");
    pEventTree->Branch("trkmm",trkmm,"trkmm[ntrkt]/F");
    pEventTree->Branch("trkth",trkth,"trkth[ntrkt]/F");
    pEventTree->Branch("trkph",trkph,"trkph[ntrkt]/F");
    
    pEventTree->Branch("momvx",momvx,"momvx[ntrkt]/F");
    pEventTree->Branch("thevx",thevx,"thevx[ntrkt]/F");
    pEventTree->Branch("phivx",phivx,"phivx[ntrkt]/F");
    pEventTree->Branch("posxvx",posxvx,"posxvx[ntrkt]/F");
    pEventTree->Branch("posyvx",posyvx,"posyvx[ntrkt]/F");
    pEventTree->Branch("poszvx",poszvx,"poszvx[ntrkt]/F");
    
    pEventTree->Branch("momend",momend,"momend[ntrkt]/F");
    pEventTree->Branch("theend",theend,"theend[ntrkt]/F");
    pEventTree->Branch("phiend",phiend,"phiend[ntrkt]/F");
    pEventTree->Branch("posxend",posxend,"posxend[ntrkt]/F");
    pEventTree->Branch("posyend",posyend,"posyend[ntrkt]/F");
    pEventTree->Branch("poszend",poszend,"poszend[ntrkt]/F");
    pEventTree->Branch("tx_end",tx_end,"tx_end[ntrkt]/F");
    pEventTree->Branch("ty_end",ty_end,"ty_end[ntrkt]/F");
    
    pEventTree->Branch("momds"   ,momds   ,"momds[ntrkt]/F");
    pEventTree->Branch("momrg"   ,momrg   ,"momrg[ntrkt]/F");
    
    pEventTree->Branch("mcxgnvx" ,mcxgnvx ,"mcxgnvx[ntrkt]/F");
    pEventTree->Branch("mcygnvx" ,mcygnvx ,"mcygnvx[ntrkt]/F");
    pEventTree->Branch("momgnvx" ,momgnvx ,"momgnvx[ntrkt]/F");
    pEventTree->Branch("thegnvx" ,thegnvx ,"thegnvx[ntrkt]/F");
    pEventTree->Branch("phignvx" ,phignvx ,"phignvx[ntrkt]/F");
    pEventTree->Branch("momgnend",momgnend,"momgnend[ntrkt]/F");
    pEventTree->Branch("thegnend",thegnend,"thegnend[ntrkt]/F");
    pEventTree->Branch("phignend",phignend,"phignend[ntrkt]/F");
    
    pEventTree->Branch("vtxzplane",vtxzplane,"vtxzplane[ntrkt]/I");
    pEventTree->Branch("endzplane",endzplane,"endzplane[ntrkt]/I");
    pEventTree->Branch("ntrkcl",ntrkcl,"ntrkcl[ntrkt]/I");
    pEventTree->Branch("ntrkst",ntrkst,"ntrkst[ntrkt]/I");
    pEventTree->Branch("range", &range, "range/F");
    
    pEventTree->Branch("tx", &tx, "tx/F");
    pEventTree->Branch("ty", &ty, "ty/F");
    pEventTree->Branch("xxin", &xxin, "xxin/F");
    pEventTree->Branch("yyin", &yyin, "yyin/F");
    pEventTree->Branch("txin", &txin, "txin/F");
    pEventTree->Branch("tyin", &tyin, "tyin/F");
    
    pEventTree->Branch("therr", &therr, "therr/F");
    pEventTree->Branch("pherr", &pherr, "pherr/F");
    pEventTree->Branch("xxerr", &xxerr, "xxerr/F");
    pEventTree->Branch("yyerr", &yyerr, "yyerr/F");
    pEventTree->Branch("txerr", &txerr, "txerr/F");
    pEventTree->Branch("tyerr", &tyerr, "tyerr/F");
    pEventTree->Branch("qperr", &qperr, "qperr/F");
    
    pEventTree->Branch("xxenderr", &xxenderr, "xxenderr/F");
    pEventTree->Branch("yyenderr", &yyenderr, "yyenderr/F");
    pEventTree->Branch("txenderr", &txenderr, "txenderr/F");
    pEventTree->Branch("tyenderr", &tyenderr, "tyenderr/F");
    pEventTree->Branch("qpenderr", &qpenderr, "qpenderr/F");
    
    //Arnab
    /*
      if(!H) H= new Hits();
      if(!Hp) Hp = new HitPos();
      pEventTree->Branch("Hits_Branch","Hits",&H,1600000,2);
      pEventTree->Branch("HitPos_Branch","HitPos",&Hp,1600000,2);
      //H->ENum = 0;
      */
  } else if (InputOutput==1 || InputOutput==4) {
    //DigiOutput
    pEventTree->Branch("ndigiht", &ndigiht, "ndigiht/i");
    pEventTree->Branch("trigx", &trigx, "trigx/I");
    pEventTree->Branch("trigy", &trigy, "trigy/I");
    pEventTree->Branch("stripid", stripid, "stripid[ndigiht]/i");
    pEventTree->Branch("digipdgid", digipdgid, "digipdgid[ndigiht]/I");
    pEventTree->Branch("digitime", digitime, "digitime[ndigiht]/I");
    pEventTree->Branch("digitruetime", digitruetime, "digitruetime[ndigiht]/I");
    pEventTree->Branch("digienr", digienr, "digienr[ndigiht]/F");
    pEventTree->Branch("digivx", digivx, "digivx[ndigiht]/F"); 
    pEventTree->Branch("digivy", digivy, "digivy[ndigiht]/F");
    pEventTree->Branch("digivz", digivz, "digivz[ndigiht]/F");
    pEventTree->Branch("digipx", digipx, "digipx[ndigiht]/F");
    pEventTree->Branch("digipy", digipy, "digipy[ndigiht]/F");
    pEventTree->Branch("digipz", digipz, "digipz[ndigiht]/F");
  } else if (InputOutput==2)  { //Sim output
    //SImulation output
    pEventTree->Branch("nsimht", &nsimht, "nsimht/i");
    pEventTree->Branch("detid", detid, "detid[nsimht]/i");
    pEventTree->Branch("simpdgid", simpdgid, "simpdgid[nsimht]/I");
    pEventTree->Branch("simtime", simtime, "simtime[nsimht]/F");
    pEventTree->Branch("simenr", simenr, "simenr[nsimht]/F");
    pEventTree->Branch("simvx", simvx, "simvx[nsimht]/F");
    pEventTree->Branch("simvy", simvy, "simvy[nsimht]/F");
    pEventTree->Branch("simvz", simvz, "simvz[nsimht]/F");
    pEventTree->Branch("simpx", simpx, "simpx[nsimht]/F");
    pEventTree->Branch("simpy", simpy, "simpy[nsimht]/F");
    pEventTree->Branch("simpz", simpz, "simpz[nsimht]/F");
    pEventTree->Branch("simlocvx", simlocvx, "simlocvx[nsimht]/F");
    pEventTree->Branch("simlocvy", simlocvy, "simlocvy[nsimht]/F");

    //pEventTree->Branch("range", &range, "range/F");
    pEventTree->Branch("simtotabenr",&simtotabenr, "simtotabenr/F");
    pEventTree->Branch("simtotrpcenr",&simtotrpcenr, "simtotrpcenr/F");
    pEventTree->Branch("simtotablen",&simtotablen, "simtotablen/F");
    pEventTree->Branch("simtotrpclen",&simtotrpclen, "simtotrpclen/F");
  }
  if (!pEventTree) {
    G4cout << "Error allocating Tree !" << G4endl;
    exit(-1);
  } 
  
  if(isVisOut==1 && (InputOutput==0 ||InputOutput==3 || InputOutput==5)) {
    pVisFile->cd();
    EveCnt=0;//Event Counter
    if(!H) H = new Hits(); //VALGRIND
    if(!Hp) Hp= new HitPos();
    nloops=0;// number of tree fills
    if(!rtree){rtree = new TTree("Hitstree","Geant Hits File");}
    rtree->Branch("Hits_Branch","Hits",&H,1600000,2);
    //H->Clear();
    //Hp->Clear();
    //H->ClearTracks();
    H->ENum=0;
  }
  // cout<<"Closing . . .void MultiSimAnalysis::OpenRootfiles(G4String infile,  G4String outfile) {"<<endl;
}

/*
void MultiSimAnalysis::CloseRootfiles(){
cout<< "MultiSimAnalysis::CloseRootfiles("<<endl;
 if (pRootFile)
  {pRootFile->cd();
    EffDist->Divide(hitDist);   //asm
    G4cout << "Histogram file writen !" << G4endl;
    pRootFile->Write(); //VALGRIND
    pRootFile->Close();

    if (pEventTree) { delete pEventTree; pEventTree=0;}
    delete pRootFile;
      }
  else {
    G4cout << "Histogram file not made !" << G4endl;
  }

if(isVisOut==1&&(InputOutput==0 ||InputOutput==3 || InputOutput==5)){
if (pVisFile) {
    pVisFile->cd();
    G4cout << "Hit Display Tree File writen !" << G4endl;
    pVisFile->Write(); //VALGRIND
    pVisFile->Close();

 if (rtree) {G4cout << "rtree !" << G4endl;delete rtree; rtree=0;}
    if (Hp) {delete Hp; Hp=0;}
    if (H) {delete H; H=0;}
  } else {
    G4cout << "No output Hit Display Tree !" << G4endl;
  }
}
  if (inputRootFile)
  {
    inputRootFile->cd();
    delete inputEventTree; inputEventTree=0;
    delete inputRootFile; inputRootFile=0;
  }  else {
    G4cout << "No inputRootFile !" << G4endl;
  }
  if (InputOutput==0 ||InputOutput==3 || InputOutput==5) {
    ascii_output.close();
  }
  if (InputOutput==0 ||InputOutput==1 || InputOutput==2) {
  micalPrimaryGeneratorAction *pgPointer = micalPrimaryGeneratorAction::AnPointer;
  pgPointer->CloseNuanceFile();

}

}
*/
void MultiSimAnalysis::CloseRootfiles() {
  // cout<< "MultiSimAnalysis::CloseRootfiles("<<endl;
  if (pRootFile) {
    pRootFile->cd();
    EffDist->Divide(hitDist);   //asm
    
    cout << "Histogram file writen !" << endl;
    pRootFile->Write(); //VALGRIND
    cout << "Histogram file writen !" << endl;

    if (pEventTree) { delete pEventTree; pEventTree=0;}
    
    
    if (pPosX) {pPosX->Write(); delete pPosX; pPosX=0;}
    if (pPosY) {pPosY->Write(); delete pPosY; pPosY=0;}
    if (pPosZ) {pPosZ->Write(); delete pPosZ; pPosZ=0;}
    
    if (pPosXX) {pPosXX->Write(); delete pPosXX; pPosXX=0;}
    if (pPosYY) {pPosYY->Write(); delete pPosYY; pPosYY=0;}
    if (pPosZZ) {pPosZZ->Write(); delete pPosZZ; pPosZZ=0;}
    
    for (int i=0; i<20; i++) {if (pdedz[i]) {pdedz[i]->Write();  delete pdedz[i]; pdedz[i]=0;}};
    
    if (hitDist) {hitDist->Write(); delete hitDist; hitDist=0;}   //asm
    if (TrkDist) {TrkDist->Write(); delete TrkDist; TrkDist=0;}   //asm
    if (EffDist) {EffDist->Write(); delete EffDist; EffDist=0;}   //asm
    if (InoTrack_listsize) {InoTrack_listsize->Write(); delete InoTrack_listsize;  InoTrack_listsize=0;} //asm
    
    if (ShwXw) {ShwXw->Write(); delete ShwXw; ShwXw=0;}
    if (ShwYw) {ShwYw->Write(); delete ShwYw; ShwYw=0;}
    if (RC) {RC->Write(); delete RC; RC=0;}
    if (DGap) {DGap->Write(); delete DGap; DGap=0;}
    
    if (DeadStripX) {DeadStripX->Write(); delete DeadStripX; DeadStripX=0;}
    if (DeadStripY) {DeadStripY->Write(); delete DeadStripY; DeadStripY=0;}
    if (NoisyStripX) {NoisyStripX->Write(); delete NoisyStripX; NoisyStripX=0;}
    if (NoisyStripY) {NoisyStripY->Write(); delete NoisyStripY; NoisyStripY=0;}
    if (DiffTime) {DiffTime->Write(); delete DiffTime; DiffTime=0;}
    if (strpXtime) {strpXtime->Write(); delete strpXtime; strpXtime=0;}    
    if (strpYtime) {strpYtime->Write(); delete strpYtime; strpYtime=0;}    
    if (strpXtimeCorr) {strpXtimeCorr->Write(); delete strpXtimeCorr; strpXtimeCorr=0;}    
    if (strpYtimeCorr) {strpYtimeCorr->Write(); delete strpYtimeCorr; strpYtimeCorr=0;}    
    if (hitXtime) {hitXtime->Write(); delete hitXtime; hitXtime=0;}    
    if (hitYtime) {hitYtime->Write(); delete hitYtime; hitYtime=0;}    
    if (smagFieldX) {smagFieldX->Write(); delete smagFieldX; smagFieldX=0;} 
    if (smagFieldY) {smagFieldY->Write(); delete smagFieldY; smagFieldY=0;} 
    if (smag2dX) {smag2dX->Write(); delete smag2dX; smag2dX=0;} 
    if (smag2dY) {smag2dY->Write(); delete smag2dY; smag2dY=0;} 
    if (rmagFieldX) {rmagFieldX->Write(); delete rmagFieldX; rmagFieldX=0;} 
    if (rmagFieldY) {rmagFieldY->Write(); delete rmagFieldY; rmagFieldY=0;} 
    if (rmag2dX) {rmag2dX->Write(); delete rmag2dX; rmag2dX=0;} 
    if (rmag2dY) {rmag2dY->Write(); delete rmag2dY; rmag2dY=0;} 
    if (smag2dXYpixel_iron) {smag2dXYpixel_iron->Write(); delete smag2dXYpixel_iron; smag2dXYpixel_iron=0;}
    if (smag2dXYpixel_air) {smag2dXYpixel_air->Write(); delete smag2dXYpixel_air; smag2dXYpixel_air=0;}
    if (rmag2dXYpixel_iron) {rmag2dXYpixel_iron->Write(); delete rmag2dXYpixel_iron; rmag2dXYpixel_iron=0;}
    if (rmag2dXYpixel_air) {rmag2dXYpixel_air->Write(); delete rmag2dXYpixel_air; rmag2dXYpixel_air=0;}
    if (xyvsbxin) {xyvsbxin->Write(); delete xyvsbxin; xyvsbxin=0;}
    if (xyvsbyin) {xyvsbyin->Write(); delete xyvsbyin; xyvsbyin=0;}
    if (xyvsbxdiff) {xyvsbxdiff->Write(); delete xyvsbxdiff; xyvsbxdiff=0;}
    if (xyvsbydiff) {xyvsbydiff->Write(); delete xyvsbydiff; xyvsbydiff=0;}
    // if (xyvsbxindiff) {xyvsbxindiff->Write(); delete xyvsbxindiff; xyvsbxindiff=0;}
    // if (xyvsbyindiff) {xyvsbyindiff->Write(); delete xyvsbyindiff; xyvsbyindiff=0;}
    if (xyvsbxout) {xyvsbxout->Write(); delete xyvsbxout; xyvsbxout=0;}
    if (xyvsbyout) {xyvsbyout->Write(); delete xyvsbyout; xyvsbyout=0;}
    //if (hh_E) {hh_E->Write(); delete hh_E; hh_E=0;}//SSE 
    //if (hh_woghst_E) {hh_woghst_E->Write(); delete hh_woghst_E; hh_woghst_E=0;}//SSE
    //if (hist_nhits_LargestCluster_E) {hist_nhits_LargestCluster_E->Write(); delete hist_nhits_LargestCluster_E; hist_nhits_LargestCluster_E=0;}//SSE
    //if (hist_orighits_new_E) {hist_orighits_new_E->Write(); delete hist_orighits_new_E; hist_orighits_new_E=0;}//SSE 09/15
    //if (hist_orighits_mod_E) {hist_orighits_mod_E->Write(); delete hist_orighits_mod_E; hist_orighits_mod_E=0;}//SSE 09/15
   // if (hist_wogh_orighits_E) {hist_wogh_orighits_E->Write(); delete hist_wogh_orighits_E; hist_wogh_orighits_E=0;}//SSE 09/15
    
    //  if (hprof) {hprof->Write(); delete hprof; hprof=0;}//SSE

    // pRootFile->ls();
    pRootFile->Close();
    delete pRootFile; pRootFile=0;
    
    cout << "Histogram file  made !" << endl;
  } else {
    cout << "Histogram file not made !" << endl;
  }
  
  if (pVisFile) {
    pVisFile->cd();
    rtree->Fill(); // fill tree
    H->Clear();
    Hp->Clear();
    H->ClearTracks();
    pVisFile->Write(); //VALGRIND
    cout<< "Write to .inh file  done"<<endl;
    if (rtree) {delete rtree; rtree=0;}
    pVisFile->Close();
    if (Hp) {delete Hp; Hp=0;}
    if (H) {delete H; H=0;}
    delete pVisFile; pVisFile=0;
  } else {
    cout << "No output Hit Display Tree !" << endl;
  }
  
  if (inputRootFile) {
    inputRootFile->cd();
    delete inputEventTree; inputEventTree=0;
    delete inputRootFile; inputRootFile=0;
    cout<<"Input root file closed."<<endl;
  } else {
    cout << "No input histograms and Tree !" << endl;
  }

  if (collatedRootFile) {
    collatedRootFile->cd();
    for(int iki=0; iki<numberInLA; iki++) {
      delete inefficiency_corx[iki];
      delete inefficiency_uncx[iki];
      delete inefficiency_uncy[iki];
      delete triggereffi_xevt[iki];
      delete triggereffi_yevt[iki];
      delete strp_xmulsim_cor[iki];
      delete strp_ymulsim_cor[iki];
    }
    collatedRootFile->Close();
    delete collatedRootFile; collatedRootFile=0;
    cout<<"Collated root file closed."<<endl;
  } else {
    cout << "No collated histograms !" << endl;
  }

  
  if ( isXtermOut==2&&(InputOutput==0 ||InputOutput==3 || InputOutput==5)) {
    ascii_output.close();
  }
  if( isXtermOut==2 &&(InputOutput==0 ||InputOutput==1 || InputOutput==2)) {
    B_ascii_output.close();
  }
  
  // timeAsciiOutput.close();
  
  micalPrimaryGeneratorAction *pgPointer = micalPrimaryGeneratorAction::AnPointer;
  // if(pgPointer->InputFlag==2)  {
  //   if(InputOutput==0) {
  //     pgPointer->CloseFileCORSIKA();
  //   }
  // }
  if((InputOutput==0 || InputOutput==1|| InputOutput==2) &&
     (pgPointer->InputFlag<=3)) {
    pgPointer->CloseFileCORSIKA();
  }
  // if(pgPointer->InputFlag==1 && (InputOutput==0 || InputOutput==1|| InputOutput==2)) {
  //   pgPointer->CloseNuanceFile();
  // }
  // pgPointer->initialise=0;
  //cout <<"end of Multisimanalysis detructor "<<endl;
}

MultiSimAnalysis::~MultiSimAnalysis() {
  // cout <<"end of Multisimanalysis detructor "<<endl;

}
void MultiSimAnalysis::SetCorrTimeError(G4double val) {
  cout<<"void MultiSimAnalysis::SetCorrTimeError(G4double "<<val<<")"<<endl;
  CorrTimeError = val;
}

void MultiSimAnalysis::SetUnCorrTimeError(G4double val) {
  cout<<"void MultiSimAnalysis::SetUnCorrTimeError(G4double "<<val<<")"<<endl;
  UnCorrTimeError = val;
}

void MultiSimAnalysis::SetTimeToDigiConvVal(G4double val) {
  cout<<"void MultiSimAnalysis::SetTimeToDigiConvVal(G4double "<<val<<")"<<endl;
  TimeToDigiConv = val;
}

void MultiSimAnalysis::SetSignalSpeedVal(G4double val) {
  cout<<"void MultiSimAnalysis::SetSignalSpeedVal(G4double "<<val<<")"<<endl;
  SignalSpeed = val;
}
