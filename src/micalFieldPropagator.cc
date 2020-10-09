#include "micalFieldPropagator.hh"

#define interpolate 1
micalFieldPropagator *micalFieldPropagator::FdPointer;
micalFieldPropagator::micalFieldPropagator() {   
  FdPointer=this;
  paradef = micalDetectorParameterDef::AnPointer;
  pAnalysis= MultiSimAnalysis::AnPointer;
  geometryIn = (InoGeometry_Manager::APointer)->icalGeometry;
  slotsize = 50;
  irlayZ = paradef->GetParirlay(2);
  rpclayZ = paradef->GetParlay(2);
  ironrpcZ = 2*(irlayZ+rpclayZ);
  nLayer = paradef->GetnLayer();
  nIRLayer = paradef->GetnIRLayer();
  gapino = 0.0;//paradef->GetGapino();
  //  nCoil =  paradef->GetnCoil(); 
  for(int ij=0;ij<3;ij++) {
    parino[ij] = paradef->GetParino(ij);
    parirlay[ij] = paradef->GetParirlay(ij);
    StackShift[ij] = paradef->GetStackPosInRoom(ij) + paradef->GetINOroomPos(ij);
    // cout<<"StackShift "<<ij<<" "<<StackShift[ij]<<" "<<paradef->GetStackPosInRoom(ij)<<" "<<paradef->GetINOroomPos(ij)<<endl;
  }
  for(int ij=0; ij<nIRLayer; ij++) {
    IRONLayerPosZ[ij] = paradef->GetIRONLayerPosZ(ij);
  }
  for(int ij=0; ij<nLayer; ij++) {
    RPCLayerPosZ[ij] = paradef->GetRPCLayerPosZ(ij);
  }
  // cout<<"--------------------------------"<<endl;
  // cout<<"irlayZ "<<irlayZ<<endl;
  // cout<<"rpclayZ "<<rpclayZ<<endl;
  // cout<<"ironrpcZ "<<ironrpcZ<<endl;
  // cout<<"nLayer "<<nLayer<<endl;
  // cout<<"nIRLayer "<<nIRLayer<<endl;
  // cout<<"gapino "<<gapino<<endl;
  // for(int ij=0;ij<3;ij++) {
  //   cout<<"parino["<<ij<<"] "<<parino[ij]<<endl;
  //   cout<<"parirlay["<<ij<<"] "<<parirlay[ij]<<endl;
  //   cout<<"StackShift["<<ij<<"] "<<StackShift[ij]<<endl;
  // }
  // cout<<"--------------------------------"<<endl;

  pMagFile = new TFile("B_mical_hist.root","read");
  
  if (!pMagFile) {
    cout << "Error Field map root opening file !" << endl;
    exit(-1);
  } else {
    cout<< " Field Map file being read opened" <<endl;
  }

  fieldxin = (TH2D*)pMagFile->Get("xyvsbxin");
  fieldyin = (TH2D*)pMagFile->Get("xyvsbyin");

  PrintFieldMap();
}

micalFieldPropagator::~micalFieldPropagator() {
  if(pMagFile) {
    pMagFile->Close();
    delete pMagFile;
    pMagFile=0;
  }
}

void micalFieldPropagator::PrintFieldMap() {
  int ixmax, iymax;
  ixmax = 80; iymax=80;
  double tmpxyz[3];
  double posx1[3];
  tmpxyz[2] = IRONLayerPosZ[5];
  double tmpbx1[2];
  cout<<"Checking field map file..."<<endl;
  for (int ix11=0; ix11<=ixmax; ix11++) {
    for (int iy11=0; iy11<=iymax; iy11++) {
      tmpxyz[0] = -parirlay[0] + ix11*50.0 + 25;// + StackShift[0];
      tmpxyz[1] = -parirlay[1] + iy11*50.0 + 25;// + StackShift[1];
      for(int ss1=0; ss1<3; ss1++) {
	posx1[ss1] = tmpxyz[ss1] + StackShift[ss1];
      }
      if((abs(tmpxyz[0])<parirlay[0]-10) && (abs(tmpxyz[1])<parirlay[1]-10)) {
	ElectroMagneticField(posx1,tmpbx1[0],tmpbx1[1],0);
	// cout<<"ij "<<ix11<<" "<<iy11<<" "<<tmpxyz[0]<<" "<<tmpxyz[1]<<" "<<posx1[0]<<" "<<posx1[1]<<" "<<tmpbx1[0]/tesla<<" "<<tmpbx1[1]/tesla<<endl;
	pAnalysis->xyvsbxout->Fill(tmpxyz[0],tmpxyz[1],tmpbx1[0]/tesla);
	pAnalysis->xyvsbyout->Fill(tmpxyz[0],tmpxyz[1],tmpbx1[1]/tesla);
      }
    }
  }
  cout<<"Checking complete field map file..."<<endl;
}

void micalFieldPropagator::ElectroMagneticField(const double xyzc[3], double &Bx, double &By, int ftype) {
  //translating the xyzc[] position to miniICAL 0-4m module. //apoorva
  int igrid[3]={0};
  double Bxf=0;
  double Byf=0;
  double tmppos[3];
  for(int iw=0; iw<3; iw++) {tmppos[iw] = xyzc[iw] - StackShift[iw];}
  if(ftype==0) {
    // cout<<"xyzc "<<xyzc[0]<<" "<<xyzc[1]<<" "<<xyzc[2]<<endl;
    // cout<<"tmppos "<<tmppos[0]<<" "<<tmppos[1]<<" "<<tmppos[2]<<endl;
  }
  igrid[0]= (int)((tmppos[0]+2000)*mm);
  igrid[1]= (int)((tmppos[1]+2000)*mm);
  igrid[2]= (int)(tmppos[2]*mm);
  
  double local_pos1[3], local_dir1[3];
  local_pos1[0] = xyzc[0]/10;
  local_pos1[1] = xyzc[1]/10;
  local_pos1[2] = xyzc[2]/10;
  
  local_dir1[0] = 0.0;
  local_dir1[1] = 0.0;
  local_dir1[2] = 1.0;
  
  geometryIn->InitTrack(local_pos1, local_dir1);
  // cout<<"xyzc[0] "<<xyzc[0]<<" "<<xyzc[1]<<" "<<xyzc[2]<<" "<<geometryIn->GetCurrentVolume()->GetName()<<endl;
  if(strstr(geometryIn->GetCurrentVolume()->GetName(),"IRLAYElog")) {
    // cout<<"xyzc[0] 111 "<<xyzc[0]<<" "<<xyzc[1]<<" "<<xyzc[2]<<" "<<geometryIn->GetCurrentVolume()->GetName()<<endl;
    F2int(igrid,Bxf,Byf);
    Bx = Bxf*tesla;
    By = Byf*tesla;
  } else {
    Bx = 0*tesla;
    By = 0*tesla;
  }
}

void micalFieldPropagator::F2int(int* ag_f2i, double& bx_f2i, double& by_f2i) {
  bx_f2i = 0.0;
  by_f2i = 0.0;
  double fgridx_loc[4],fgridy_loc[4],pgrid_loc[4];
  int binx = ag_f2i[0]/slotsize;
  int biny = ag_f2i[1]/slotsize;
  // cout<<"bins "<<ag_f2i[0]<<" "<<ag_f2i[1]<<" "<<binx<<" "<<biny<<endl;
  pgrid_loc[0] = fieldxin->GetXaxis()->GetBinLowEdge(binx+1)/1000;
  pgrid_loc[1] = fieldxin->GetXaxis()->GetBinLowEdge(binx+2)/1000;
  pgrid_loc[2] = fieldxin->GetYaxis()->GetBinLowEdge(biny+1)/1000;
  pgrid_loc[3] = fieldxin->GetYaxis()->GetBinLowEdge(biny+2)/1000;
  fgridx_loc[0] = fieldxin->GetBinContent(binx+2,biny+2);
  fgridx_loc[1] = fieldxin->GetBinContent(binx+1,biny+2);
  fgridx_loc[2] = fieldxin->GetBinContent(binx+2,biny+1);
  fgridx_loc[3] = fieldxin->GetBinContent(binx+1,biny+1);
  fgridy_loc[0] = fieldyin->GetBinContent(binx+2,biny+2);
  fgridy_loc[1] = fieldyin->GetBinContent(binx+1,biny+2);
  fgridy_loc[2] = fieldyin->GetBinContent(binx+2,biny+1);
  fgridy_loc[3] = fieldyin->GetBinContent(binx+1,biny+1);
  bx_f2i = bilinearInterpolation(fgridx_loc,pgrid_loc,ag_f2i);
  by_f2i = bilinearInterpolation(fgridy_loc,pgrid_loc,ag_f2i);
  if (fabs(bx_f2i)< 2 && fabs(bx_f2i) >0) {
    // return bx_f2i;
  } else {
    cout<<"nan  : "<< ag_f2i[0]<<" "<<ag_f2i[1] <<" "<< bx_f2i<<endl;
    bx_f2i = 0;
    // return 0;
  }
  if (fabs(by_f2i)< 2 && fabs(by_f2i) >0) {
    // return by_f2i;
  } else {
    cout<<"nan  : "<< ag_f2i[0]<<" "<<ag_f2i[1] <<" "<< by_f2i<<endl;
    by_f2i = 0;
    // return 0;
  }
}

double micalFieldPropagator::bilinearInterpolation(double* f_bli,double* arg_bli,int* ag_bli) {
  double x_loc1= arg_bli[0];
  double x_loc2= arg_bli[1];
  double y_loc1= arg_bli[2];
  double y_loc2= arg_bli[3];
  double f11_loc = f_bli[0];
  double f12_loc = f_bli[2];
  double f21_loc = f_bli[1];
  double f22_loc = f_bli[3];
  double x_loc= (double)ag_bli[0]/1000;
  double y_loc= (double)ag_bli[1]/1000;
  double term1_loc = f11_loc*(x_loc2-x_loc)*(y_loc2-y_loc)/((x_loc2-x_loc1)*(y_loc2-y_loc1));  
  double term2_loc = f21_loc*(x_loc1-x_loc)*(y_loc2-y_loc)/((x_loc1-x_loc2)*(y_loc2-y_loc1));
  double term3_loc = f12_loc*(x_loc2-x_loc)*(y_loc1-y_loc)/((x_loc2-x_loc1)*(y_loc1-y_loc2));
  double term4_loc = f22_loc*(x_loc1-x_loc)*(y_loc1-y_loc)/((x_loc1-x_loc2)*(y_loc1-y_loc2));
  // cout << "x_loc " << x_loc << " " << x_loc1 << " " << x_loc2 << endl;
  // cout << "y_loc " << y_loc << " " << y_loc1 << " " << y_loc2 << endl;
  // cout << "f11_loc " << f11_loc << " f12_loc " << f12_loc << " f21_loc " << f21_loc << " f22_loc " << f22_loc << endl;
  // cout << "term1_loc " << term1_loc << " term2_loc " << term2_loc << " term3_loc " << term3_loc << " term4_loc " << term4_loc << endl;
  // cout << "bli " << term1_loc+term2_loc+term3_loc+term4_loc << endl;
  return term1_loc+term2_loc+term3_loc+term4_loc;
}
