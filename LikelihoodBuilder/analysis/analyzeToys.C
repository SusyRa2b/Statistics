#include <iostream>
#include <fstream>

#include "TStyle.h"
#include "TROOT.h"
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
  else if(fitter == "test") {
    //branchDescriptor = "zeroLeptonCountTotal";
    branchDescriptor = "XXsignalCrossSection:signalCrossSectionError:signalCrossSection";
  }
  else if(fitter == "sensitivity") {
    branchDescriptor = "signalCrossSection/D:signalCrossSectionError";
    branchDescriptor += ":lowerLimit95:upperLimit95";
    branchDescriptor += ":susyFloatingNLL:susyFixedToZeroNLL";
  }
  else {
    branchDescriptor = "signalCrossSection/D:signalCrossSectionError";
    branchDescriptor += ":lowerLimit68:upperLimit68:lowerLimit95:upperLimit95";
    branchDescriptor += ":signalUncertainty";
    branchDescriptor += ":zeroLeptonSignalYieldTotal:zeroLeptonSignalYieldTotalError";
    branchDescriptor += ":zeroLeptonTopWJetsYieldTotal:zeroLeptonQCDYieldTotal:zeroLeptonZtoNuNuYieldTotal:zeroLeptonDibosonYieldTotal:zeroLeptonCountTotal";
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
  
  gStyle->SetStatY(0.98);
  gStyle->SetStatX(0.98);
  gStyle->SetStatW(0.3);
  gStyle->SetStatH(0.3);

  /*
  TH1D* hSignalCrossSection = new TH1D("hSignalCrossSection", "Signal cross section", 50, 0, 200);
  TH1D* hPull = new TH1D("hPull", "Pull using fit errors", 50, -5, 5);
  TH1D* hPullPL = new TH1D("hPullPL", "Pull using PL errors", 50, -5, 5);
  TH1D* hZeroLeptonSignalYieldTotal = new TH1D("hZeroLeptonSignalYieldTotal", "Zero lepton signal yield", 50, 0, 200);
  TH1D* hZeroLeptonTopWJetsYieldTotal = new TH1D("hZeroLeptonTopWJetsYieldTotal", "Zero lepton top W-jets yield", 50, 1500, 2200);
  TH1D* hZeroLeptonQCDYieldTotal = new TH1D("hZeroLeptonQCDYieldTotal", "Zero lepton QCD yield", 50, 0, 500);
  TH1D* hZeroLeptonZtoNuNuYieldTotal = new TH1D("hZeroLeptonZtoNuNuYieldTotal", "Zero lepton Z-invisible yield", 50, 150, 450);
  */
  TH1D* hSignalCrossSection = new TH1D("hSignalCrossSection", "Signal cross section", 50, 0, 200);
  TH1D* hPull = new TH1D("hPull", "Pull using fit errors", 50, -5, 5);
  TH1D* hPullPL = new TH1D("hPullPL", "Pull using PL errors", 50, -5, 5);
  TH1D* hZeroLeptonSignalYieldTotal = new TH1D("hZeroLeptonSignalYieldTotal", "Zero lepton signal yield", 50, 0, 500);
  TH1D* hZeroLeptonTopWJetsYieldTotal = new TH1D("hZeroLeptonTopWJetsYieldTotal", "Zero lepton top W-jets yield", 50, 9000, 16000);
  TH1D* hZeroLeptonQCDYieldTotal = new TH1D("hZeroLeptonQCDYieldTotal", "Zero lepton QCD yield", 50, 15000, 22000);
  TH1D* hZeroLeptonZtoNuNuYieldTotal = new TH1D("hZeroLeptonZtoNuNuYieldTotal", "Zero lepton Z-invisible yield", 50, 1000, 2000);
  TH1D* hZeroLeptonYieldTotal = new TH1D("hZeroLeptonYieldTotal", "Zero lepton yield", 50, 25000, 28000);
  TH1D* hZeroLeptonCountTotal = new TH1D("hZeroLeptonCountTotal", "Zero lepton count", 50, 25000, 28000);
  TH1D* hUpperLimit95 = new TH1D("hUpperLimit95", "95% Confidence Level Upper Limit", 50, 0, 200);
  TH1D* hSignificance = new TH1D("hSignificance", "Significance", 50, 0, 5);


  hSignalCrossSection->SetFillColor(6);
  hZeroLeptonSignalYieldTotal->SetFillColor(6);  
  hZeroLeptonTopWJetsYieldTotal->SetFillColor(kBlue-9);
  hZeroLeptonQCDYieldTotal->SetFillColor(2);
  hZeroLeptonZtoNuNuYieldTotal->SetFillColor(kGreen-3);
  hZeroLeptonYieldTotal->SetFillColor(kGray+1);
  hZeroLeptonCountTotal->SetFillColor(kGray+1);
  hUpperLimit95->SetFillColor(kGray+1);
  hSignificance->SetFillColor(kGray+1);

  treeToys->Project("hSignalCrossSection","signalCrossSection");
  treeToys->Project("hZeroLeptonSignalYieldTotal", "zeroLeptonSignalYieldTotal");
  treeToys->Project("hZeroLeptonTopWJetsYieldTotal", "zeroLeptonTopWJetsYieldTotal");
  treeToys->Project("hZeroLeptonQCDYieldTotal", "zeroLeptonQCDYieldTotal");
  treeToys->Project("hZeroLeptonZtoNuNuYieldTotal", "zeroLeptonZtoNuNuYieldTotal");
  treeToys->Project("hZeroLeptonYieldTotal", "zeroLeptonSignalYieldTotal+zeroLeptonTopWJetsYieldTotal+zeroLeptonQCDYieldTotal+zeroLeptonZtoNuNuYieldTotal+zeroLeptonDibosonYieldTotal");
  treeToys->Project("hZeroLeptonCountTotal", "zeroLeptonCountTotal");

  treeToys->Project("hPullPL", "(signalCrossSection-"+trueXsec_string+")/( (upperLimit68-lowerLimit68)/2.0 )");
  treeToys->Project("hPull", "(signalCrossSection-"+trueXsec_string+")/( signalCrossSectionError )");

  treeToys->Project("hUpperLimit95", "upperLimit95");
  treeToys->Project("hSignificance", "sqrt(2.0*(susyFixedToZeroNLL - susyFloatingNLL))");

  cout << "CrossSection = " << hSignalCrossSection->GetMean() << " +- " << hSignalCrossSection->GetRMS() << endl;
  cout << "Yield Total = " << hZeroLeptonYieldTotal->GetMean() << " +- " << hZeroLeptonYieldTotal->GetRMS() << endl;
  cout << "Count Total = " << hZeroLeptonCountTotal->GetMean() << " +- " << hZeroLeptonCountTotal->GetRMS() << endl;

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
  //hPullPL->Draw();
  hPull->Draw();
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
  
  TCanvas* cZeroLeptonYieldTotal = new TCanvas("cZeroLeptonYieldTotal", "Zero lepton yield", 2*400, 330);
  cZeroLeptonYieldTotal->Divide(2,1);
  cZeroLeptonYieldTotal->cd(1);
  hZeroLeptonYieldTotal->Draw();
  cZeroLeptonYieldTotal->cd(2);
  hZeroLeptonCountTotal->Draw();
  cZeroLeptonYieldTotal->Print("cZeroLeptonYield_"+trueN_string+"_"+fitter+".pdf");

  TCanvas* cSensitivity = new TCanvas("cSensitivity", "Sensitivity", 800, 330);
  cSensitivity->Divide(2,1);
  cSensitivity->cd(1);
  hUpperLimit95->Draw();
  cSensitivity->cd(2);
  hSignificance->Draw();
  cSensitivity->Print("cSensitivity_"+trueN_string+"_"+fitter+".pdf");


  return;
}


void initialize(int trueN_in, TString inFile, TString fitter = "LB") {

  trueN = trueN_in;
  makeTree(trueN, inFile, fitter);

  double trueN_double = trueN;

  //last number should not included skipped bins

  //signal point t1bbbb 850 600, lumi = 19.399/fb
  trueXsec = 1000.0 * trueN_double / 1.9399 / 2651.34 ; 
  
  //signal point t1tttt 1175 400, lumi = 19.399/fb
  //trueXsec = 1000.0 * trueN_double / 1.9399 / 13532.1 ;

  cout << "true xsec = " << trueXsec << endl;

  trueXsec_string = "";
  trueXsec_string += trueXsec;

  //true ttbar, qcd, znn, lumi = 15/fb, useExpected0lep
  tZeroLeptonTopWJetsYieldTotal = 18825.470703;
  tZeroLeptonQCDYieldTotal = 18169.785156;
  tZeroLeptonZtoNuNuYieldTotal = 1451.139893;

}

