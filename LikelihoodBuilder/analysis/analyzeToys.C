#include <iostream>
#include <fstream>

#include "TString.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TFile.h"

using namespace std;

double tZeroLeptonTopWJetsYieldTotal=0;
double tZeroLeptonQCDYieldTotal=0;
double tZeroLeptonZtoNuNuYieldTotal=0;

int trueN = 0; 
TString trueXsec_string = "0.0";
double trueXsec=0;


TTree* treeToys = new TTree("treeToys","treeToys");

bool drawTrue = true;

void makeTree(int trueN_in, TString inFile, TString fitter = "LB") {

  
  TString branchDescriptor; 
  if(fitter == "OAK") {
    branchDescriptor = "zeroLeptonSignalYieldTotal:zeroLeptonTopWJetsYieldTotal:zeroLeptonQCDYieldTotal:zeroLeptonZtoNuNuYieldTotal";
  }
  else {
    branchDescriptor = "signalCrossSection:signalCrossSectionError";
    branchDescriptor += ":lowerLimit68:upperLimit68:lowerLimit95:upperLimit95";
    branchDescriptor += ":signalUncertainty";
    branchDescriptor += ":zeroLeptonSignalYieldTotal:zeroLeptonSignalYieldTotalError";
    branchDescriptor += ":zeroLeptonTopWJetsYieldTotal:zeroLeptonQCDYieldTotal:zeroLeptonZtoNuNuYieldTotal";
  }

  treeToys->ReadFile(inFile,branchDescriptor);

  //TString foutName = ""; foutName += "tree_"; foutName += trueN; foutName += "_"; foutName += fitter; foutName += ".root";
  //TFile fout(foutName, "RECREATE");
  //treeToys->Write();
  //fout.Close();

}


void drawTest() {
  treeToys->Draw("zeroLeptonTopWJetsYieldTotal");
}


void makeHists(TString fitter = "LB") {
  
  
  /*
  TH1D* hSignalCrossSection = new TH1D("hSignalCrossSection", "Signal cross section", 50, 0, 200);
  TH1D* hPull = new TH1D("hPull", "Pull using fit errors", 50, -5, 5);
  TH1D* hPullPL = new TH1D("hPullPL", "Pull using PL errors", 50, -5, 5);
  TH1D* hZeroLeptonSignalYieldTotal = new TH1D("hZeroLeptonSignalYieldTotal", "Zero lepton signal yield", 50, 0, 200);
  TH1D* hZeroLeptonTopWJetsYieldTotal = new TH1D("hZeroLeptonTopWJetsYieldTotal", "Zero lepton top W-jets yield", 50, 1500, 2200);
  TH1D* hZeroLeptonQCDYieldTotal = new TH1D("hZeroLeptonQCDYieldTotal", "Zero lepton QCD yield", 50, 0, 500);
  TH1D* hZeroLeptonZtoNuNuYieldTotal = new TH1D("hZeroLeptonZtoNuNuYieldTotal", "Zero lepton Z-invisible yield", 50, 150, 450);
  */
  TH1D* hSignalCrossSection = new TH1D("hSignalCrossSection", "Signal cross section", 50, 0, 100);
  TH1D* hPull = new TH1D("hPull", "Pull using fit errors", 50, -5, 5);
  TH1D* hPullPL = new TH1D("hPullPL", "Pull using PL errors", 50, -5, 5);
  TH1D* hZeroLeptonSignalYieldTotal = new TH1D("hZeroLeptonSignalYieldTotal", "Zero lepton signal yield", 50, 0, 500);
  TH1D* hZeroLeptonTopWJetsYieldTotal = new TH1D("hZeroLeptonTopWJetsYieldTotal", "Zero lepton top W-jets yield", 50, 15000, 22000);
  TH1D* hZeroLeptonQCDYieldTotal = new TH1D("hZeroLeptonQCDYieldTotal", "Zero lepton QCD yield", 50, 15000, 22000);
  TH1D* hZeroLeptonZtoNuNuYieldTotal = new TH1D("hZeroLeptonZtoNuNuYieldTotal", "Zero lepton Z-invisible yield", 50, 1000, 2000);
  
  hSignalCrossSection->SetFillColor(6);
  hZeroLeptonSignalYieldTotal->SetFillColor(6);  
  hZeroLeptonTopWJetsYieldTotal->SetFillColor(kBlue-9);
  hZeroLeptonQCDYieldTotal->SetFillColor(2);
  hZeroLeptonZtoNuNuYieldTotal->SetFillColor(kGreen-3);

  treeToys->Project("hSignalCrossSection","signalCrossSection");
  treeToys->Project("hZeroLeptonSignalYieldTotal", "zeroLeptonSignalYieldTotal");
  treeToys->Project("hZeroLeptonTopWJetsYieldTotal", "zeroLeptonTopWJetsYieldTotal");
  treeToys->Project("hZeroLeptonQCDYieldTotal", "zeroLeptonQCDYieldTotal");
  treeToys->Project("hZeroLeptonZtoNuNuYieldTotal", "zeroLeptonZtoNuNuYieldTotal");

  treeToys->Project("hPullPL", "(signalCrossSection-"+trueXsec_string+")/( (upperLimit68-lowerLimit68)/2.0 )");

  //outfile
  TString foutName = ""; foutName += "output_"; foutName += trueN; foutName += "_"; foutName += fitter; foutName += ".dat";
  ofstream fout;
  fout.open(foutName.Data(), ios::out | ios::trunc);
  assert(fout.is_open());
  fout << trueN << " " << hZeroLeptonSignalYieldTotal->GetMean() << " " << hZeroLeptonSignalYieldTotal->GetRMS() << " ";
  fout << hZeroLeptonTopWJetsYieldTotal->GetMean() << " " << hZeroLeptonTopWJetsYieldTotal->GetRMS() << " " ;
  fout << hZeroLeptonQCDYieldTotal->GetMean() << " " << hZeroLeptonQCDYieldTotal->GetRMS() << " " ;
  fout << hZeroLeptonZtoNuNuYieldTotal->GetMean() << " " << hZeroLeptonZtoNuNuYieldTotal->GetRMS() << endl;
  fout.close();


  //Plots
  TString trueN_string = ""; trueN_string+=trueN;
  TCanvas* cBackground = new TCanvas("cBackground", "Backgrounds", 1200, 330);
  cBackground->Divide(3,1);
  cBackground->cd(1);
  hZeroLeptonTopWJetsYieldTotal->Draw();
  TLine * lttwj = new TLine(tZeroLeptonTopWJetsYieldTotal, 0, tZeroLeptonTopWJetsYieldTotal, hZeroLeptonTopWJetsYieldTotal->GetMaximum());
  lttwj->SetLineWidth(3);
  lttwj->SetLineColor(kOrange-3);
  lttwj->SetLineStyle(2);
  if(drawTrue) lttwj->Draw();

  cBackground->cd(2);
  hZeroLeptonQCDYieldTotal->Draw();
  TLine * lqcd = new TLine(tZeroLeptonQCDYieldTotal, 0, tZeroLeptonQCDYieldTotal, hZeroLeptonQCDYieldTotal->GetMaximum());
  lqcd->SetLineWidth(3);
  lqcd->SetLineColor(kOrange-3);
  lqcd->SetLineStyle(2);
  if(drawTrue) lqcd->Draw();

  cBackground->cd(3);
  hZeroLeptonZtoNuNuYieldTotal->Draw();
  TLine * lznn = new TLine(tZeroLeptonZtoNuNuYieldTotal, 0, tZeroLeptonZtoNuNuYieldTotal, hZeroLeptonZtoNuNuYieldTotal->GetMaximum());
  lznn->SetLineWidth(3);
  lznn->SetLineColor(kOrange-3);
  lznn->SetLineStyle(2);
  if(drawTrue) lznn->Draw();

  cBackground->Print("cBackground_"+trueN_string+"_"+fitter+".pdf");

  TCanvas* cSignalCrossSection = new TCanvas("cSignalCrossSection", "Signal cross section", 800, 330);
  cSignalCrossSection->Divide(2,1);
  cSignalCrossSection->cd(1);
  hSignalCrossSection->Draw();
  TLine * lsig = new TLine(trueXsec, 0, trueXsec, hSignalCrossSection->GetMaximum());
  lsig->SetLineWidth(3);
  lsig->SetLineColor(kOrange-3);
  lsig->SetLineStyle(2);
  if(drawTrue) lsig->Draw();

  cSignalCrossSection->cd(2);
  hPullPL->Draw();
  cSignalCrossSection->Print("cSignalCrossSection_"+trueN_string+"_"+fitter+".pdf");

  TCanvas * cZeroLeptonSignalYieldTotal = new TCanvas("cZeroLeptonSignalYieldTotal", "Zero lepton signal yield", 400, 330);
  cZeroLeptonSignalYieldTotal->cd();
  hZeroLeptonSignalYieldTotal->Draw();
  TLine * lsigy = new TLine(trueN, 0, trueN, hZeroLeptonSignalYieldTotal->GetMaximum());
  lsigy->SetLineWidth(3);
  lsigy->SetLineColor(kOrange-3);
  lsigy->SetLineStyle(2);
  if(drawTrue) lsigy->Draw();

  cZeroLeptonSignalYieldTotal->Print("cZeroLeptonSignalYieldTotal_"+trueN_string+"_"+fitter+".pdf");
  
  return;
}


void initialize(int trueN_in, TString inFile, TString fitter = "LB") {

  trueN = trueN_in;
  makeTree(trueN, inFile, fitter);

  double trueN_double = trueN;

  //signal point 850 600, lumi = 15/fb
  trueXsec = 1000.0 * trueN_double / 1.5 / 2119.0 ; 
  //trueXsec = 52;//3x3 with 5/fb
  trueXsec_string = "";
  trueXsec_string += trueXsec;

  //true ttbar, qcd, znn, lumi = 15/fb, useExpected0lep
  tZeroLeptonTopWJetsYieldTotal = 18825.470703;
  tZeroLeptonQCDYieldTotal = 18169.785156;
  tZeroLeptonZtoNuNuYieldTotal = 1451.139893;

}

