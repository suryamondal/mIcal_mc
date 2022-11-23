
void plotFluxThetaPhi() {

  auto FileFLUX = new TFile("flux_Pethu_FLUKA_SIBYLL.root","read");

  auto muFlux  = (TH3F*)FileFLUX->Get("muFlux");
  auto mupFlux = (TH3F*)FileFLUX->Get("mupFlux");
  auto munFlux = (TH3F*)FileFLUX->Get("munFlux");

  TH1D *pidin = new TH1D("pidin", "pidin", 41, -20.5, 20.5);
  TH1D *theinTo   = new TH1D("theinTo", "#theta : Going to",
			     45, 0., TMath::Pi());
  TH1D *phiinTo   = new TH1D("phiinTo", "#phi : Going to",
			     45, -TMath::Pi(), TMath::Pi());
  TH1D *theinFrom = new TH1D("theinFrom", "#theta : Coming from",
			     45, 0., TMath::Pi());
  TH1D *phiinFrom = new TH1D("phiinFrom", "#phi : Coming from",
			     45, -TMath::Pi(), TMath::Pi());

  {
    // SetBinLabel
    int gBinN = phiinTo->FindBin(0);
    int gBinW = phiinTo->FindBin(TMath::Pi()/2.);
    int gBinE = phiinTo->FindBin(-TMath::Pi()/2.);
    phiinTo->GetXaxis()->SetBinLabel(gBinN, "N");
    phiinTo->GetXaxis()->SetBinLabel(gBinW, "W");
    phiinTo->GetXaxis()->SetBinLabel(gBinE, "E");
    phiinFrom->GetXaxis()->SetBinLabel(gBinN, "N");
    phiinFrom->GetXaxis()->SetBinLabel(gBinW, "W");
    phiinFrom->GetXaxis()->SetBinLabel(gBinE, "E");

    phiinTo->GetXaxis()->SetLabelSize(0.075);
    phiinFrom->GetXaxis()->SetLabelSize(0.075);
  }
  
  for (int count=0; count<1000000; count++) {

    double Pxx,Pyy,Pzz;
    muFlux->GetRandom3(Pxx,Pyy,Pzz);
    int gBin = muFlux->FindBin(Pxx,Pyy,Pzz);
    double mupCnt = mupFlux->GetBinContent(gBin);
    double munCnt = munFlux->GetBinContent(gBin);

    int partID;
    if(gRandom->Uniform()*(mupCnt+munCnt)>munCnt) {
      partID = -13;
    } else {
      partID = 13;
    }
    pidin->Fill(partID);
    // cout << " partID " << partID << endl;

    TVector3 tmp3v(Pxx,Pyy,Pzz);
    theinTo->Fill(tmp3v.Theta());
    phiinTo->Fill(tmp3v.Phi());

    TVector3 tmp3vFrom = -1.*tmp3v;
    theinFrom->Fill(tmp3vFrom.Theta());
    phiinFrom->Fill(tmp3vFrom.Phi());
  }

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","c1", 1200, 400);
  c1->Divide(3,2);

  c1->cd(1);
  pidin->DrawNormalized();
  c1->cd(2);
  theinTo->DrawNormalized();
  c1->cd(3);
  phiinTo->DrawNormalized();

  c1->cd(5);
  theinFrom->DrawNormalized();
  c1->cd(6);
  phiinFrom->DrawNormalized();
}

  
