#include <iostream>
#include <fstream>

#include "TString.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLine.h"

using namespace std;

void drawSummary() {

  double tZeroLeptonTopWJetsYieldTotal = 18825.470703;
  double tZeroLeptonQCDYieldTotal = 18169.785156;
  double tZeroLeptonZtoNuNuYieldTotal = 1451.139893;


  std::vector<double> vSignalInjected;
  std::vector<double> vSignalInjectedError;
  std::vector<double> vSignal;
  std::vector<double> vSignalRMS;

  std::vector<double> vTopWJets;
  std::vector<double> vTopWJetsRMS;
  std::vector<double> vQCD;
  std::vector<double> vQCDRMS;
  std::vector<double> vZtoNuNu;
  std::vector<double> vZtoNuNuRMS;

  fstream f;
  f.open("output_LB_cat.dat");
  assert(f.is_open());

  int n=0;

  double nSignalInjected, nSignal, nSignalRMS, nTopWJets, nTopWJetsRMS, nQCD, nQCDRMS, nZtoNuNu, nZtoNuNuRMS;
  while(f>>nSignalInjected>>nSignal>>nSignalRMS>>nTopWJets>>nTopWJetsRMS>>nQCD>>nQCDRMS>>nZtoNuNu>>nZtoNuNuRMS) {
    n++;
    vSignalInjectedError.push_back(0);

    cout << nSignalInjected << " " << nSignal << " " << nSignalRMS << endl;
    vSignalInjected.push_back(nSignalInjected);
    vSignal.push_back(nSignal);
    vSignalRMS.push_back(nSignalRMS);
    vTopWJets.push_back(nTopWJets);
    vTopWJetsRMS.push_back(nTopWJetsRMS);
    vQCD.push_back(nQCD);
    vQCDRMS.push_back(nQCDRMS);
    vZtoNuNu.push_back(nZtoNuNu);
    vZtoNuNuRMS.push_back(nZtoNuNu);

  }

  TCanvas * cSignal = new TCanvas("cSignal", "Signal");
  cSignal->cd();
  TGraphErrors * gSignal = new TGraphErrors(n,&vSignalInjected[0],&vSignal[0],&vSignalInjectedError[0],&vSignalRMS[0]);
  gSignal->SetTitle("");
  gSignal->SetMarkerStyle(20);
  gSignal->SetMarkerColor(6);
  gSignal->Draw("AP");
  gSignal->SetMinimum(0);
  gSignal->SetMaximum(250);

  TF1 * f1 = new TF1("f1", "x", 0,2e2);
  f1->SetLineColor(kOrange-3);
  f1->SetLineStyle(9);
  f1->Draw("SAME");
 
  cSignal->Print("cSignal.pdf");

  
  TCanvas * cBackgrounds = new TCanvas("cBackgrounds", "Backgrounds",1200,330);
  cBackgrounds->Divide(3,1);
  cBackgrounds->cd(1);
  TGraphErrors * gTopWJets = new TGraphErrors(n,&vSignalInjected[0],&vTopWJets[0],&vSignalInjectedError[0],&vTopWJetsRMS[0]);
  gTopWJets->SetTitle("");
  gTopWJets->SetMarkerStyle(20);
  gTopWJets->SetMarkerColor(kBlue-9);
  gTopWJets->Draw("AP");
  gTopWJets->SetMinimum(17600);
  gTopWJets->SetMaximum(19300);
  TLine * lttwj = new TLine(0, tZeroLeptonTopWJetsYieldTotal, 110, tZeroLeptonTopWJetsYieldTotal);
  lttwj->SetLineWidth(3);
  lttwj->SetLineColor(kOrange-3);
  lttwj->SetLineStyle(2);
  lttwj->Draw();

  cBackgrounds->cd(2);
  TGraphErrors * gQCD = new TGraphErrors(n,&vSignalInjected[0],&vQCD[0],&vSignalInjectedError[0],&vQCDRMS[0]);
  gQCD->SetTitle("");
  gQCD->SetMarkerStyle(20);
  gQCD->SetMarkerColor(2);
  gQCD->Draw("AP");
  gQCD->SetMinimum(17700);
  gQCD->SetMaximum(19200);
  TLine * lqcd = new TLine(0, tZeroLeptonQCDYieldTotal, 110, tZeroLeptonQCDYieldTotal);
  lqcd->SetLineWidth(3);
  lqcd->SetLineColor(kOrange-3);
  lqcd->SetLineStyle(2);
  lqcd->Draw();


  cBackgrounds->cd(3);
  TGraphErrors * gZtoNuNu = new TGraphErrors(n,&vSignalInjected[0],&vZtoNuNu[0],&vSignalInjectedError[0],&vZtoNuNuRMS[0]);
  gZtoNuNu->SetTitle("");
  gZtoNuNu->SetMarkerStyle(20);
  gZtoNuNu->SetMarkerColor(kGreen-3);
  gZtoNuNu->Draw("AP");
  gZtoNuNu->SetMinimum(0);
  gZtoNuNu->SetMaximum(3200);
  TLine * lznn = new TLine(0, tZeroLeptonZtoNuNuYieldTotal, 110, tZeroLeptonZtoNuNuYieldTotal);
  lznn->SetLineWidth(3);
  lznn->SetLineColor(kOrange-3);
  lznn->SetLineStyle(2);
  lznn->Draw();
  

cBackgrounds->Print("cBackgrounds.pdf");

  

}
