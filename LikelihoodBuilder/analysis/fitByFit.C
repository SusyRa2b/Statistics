#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"

using namespace std;

void fitByFit(){

  TChain * tLB = new TChain("treeToys");
  tLB->Add("tree_100_LB.root");
  cout << "Number of LB entries: " << tLB->GetEntries() << endl;

  TChain * tOAK = new TChain("treeToys");
  tOAK->Add("tree_100_OAK.root");
  cout << "Number of OAK entries: " << tOAK->GetEntries() << endl;
  
  unsigned int ntot =  tOAK->GetEntries();


  //Variables from tree
  float zeroLeptonSignalYieldTotal_LB, zeroLeptonSignalYieldTotal_OAK;
  float zeroLeptonTopWJetsYieldTotal_LB, zeroLeptonTopWJetsYieldTotal_OAK;
  float zeroLeptonQCDYieldTotal_LB, zeroLeptonQCDYieldTotal_OAK;
  float zeroLeptonZtoNuNuYieldTotal_LB, zeroLeptonZtoNuNuYieldTotal_OAK;
  tLB->SetBranchAddress("zeroLeptonSignalYieldTotal", &zeroLeptonSignalYieldTotal_LB);
  tOAK->SetBranchAddress("zeroLeptonSignalYieldTotal", &zeroLeptonSignalYieldTotal_OAK);
  tLB->SetBranchAddress("zeroLeptonTopWJetsYieldTotal", &zeroLeptonTopWJetsYieldTotal_LB);
  tOAK->SetBranchAddress("zeroLeptonTopWJetsYieldTotal", &zeroLeptonTopWJetsYieldTotal_OAK);
  tLB->SetBranchAddress("zeroLeptonQCDYieldTotal", &zeroLeptonQCDYieldTotal_LB);
  tOAK->SetBranchAddress("zeroLeptonQCDYieldTotal", &zeroLeptonQCDYieldTotal_OAK);
  tLB->SetBranchAddress("zeroLeptonZtoNuNuYieldTotal", &zeroLeptonZtoNuNuYieldTotal_LB);
  tOAK->SetBranchAddress("zeroLeptonZtoNuNuYieldTotal", &zeroLeptonZtoNuNuYieldTotal_OAK);


  //Histograms
  TH1F* hZeroLeptonSignalYieldTotalDifference = new TH1F("hZeroLeptonSignalYieldTotalDifference", "Signal Absolute Difference", 50, -10, 10);
  TH1F* hZeroLeptonSignalYieldTotalPercentDifference = new TH1F("hZeroLeptonSignalYieldTotalPercentDifference", "Signal Percent Difference", 50, -250, 250);
  
  TH1F* hZeroLeptonTopWJetsYieldTotalDifference = new TH1F("hZeroLeptonTopWJetsYieldTotalDifference", "TopWJets Absolute Difference", 50, -100, 100);
  TH1F* hZeroLeptonQCDYieldTotalDifference = new TH1F("hZeroLeptonQCDYieldTotalDifference", "QCD Absolute Difference", 50, -100, 100);
  TH1F* hZeroLeptonZtoNuNuYieldTotalDifference = new TH1F("hZeroLeptonZtoNuNuYieldTotalDifference", "ZtoNuNu Absolute Difference", 50, -20, 20);
  
  TH1F* hZeroLeptonTopWJetsYieldTotalPercentDifference = new TH1F("hZeroLeptonTopWJetsYieldTotalPercentDifference", "TopWJets Percent Difference", 50, -1, 1);
  TH1F* hZeroLeptonQCDYieldTotalPercentDifference = new TH1F("hZeroLeptonQCDYieldTotalPercentDifference", "QCD Percent Difference", 50, -1, 1);
  TH1F* hZeroLeptonZtoNuNuYieldTotalPercentDifference = new TH1F("hZeroLeptonZtoNuNuYieldTotalPercentDifference", "ZtoNuNu Percent Difference", 50, -1, 1);
  
  TH1F* hZeroLeptonSignalYieldTotal_LB = new TH1F("hZeroLeptonSignalYieldTotal_LB", "ZeroLeptonSignalYieldTotal", 50, 0, 200);
  TH1F* hZeroLeptonSignalYieldTotal_OAK = new TH1F("hZeroLeptonSignalYieldTotal_OAK", "ZeroLeptonSignalYieldTotal", 50, 0, 200);

  //Style
  hZeroLeptonSignalYieldTotalDifference->SetFillColor(6);
  hZeroLeptonSignalYieldTotalPercentDifference->SetFillColor(6);
 
  hZeroLeptonTopWJetsYieldTotalDifference->SetFillColor(kBlue-9);
  hZeroLeptonTopWJetsYieldTotalPercentDifference->SetFillColor(kBlue-9);
  hZeroLeptonQCDYieldTotalDifference->SetFillColor(2);
  hZeroLeptonQCDYieldTotalPercentDifference->SetFillColor(2);
  hZeroLeptonZtoNuNuYieldTotalDifference->SetFillColor(kGreen-3);
  hZeroLeptonZtoNuNuYieldTotalPercentDifference->SetFillColor(kGreen-3);
  
  hZeroLeptonSignalYieldTotal_LB->SetLineWidth(2);
  hZeroLeptonSignalYieldTotal_OAK->SetLineWidth(2);
  hZeroLeptonSignalYieldTotal_LB->SetLineColor(6);
  hZeroLeptonSignalYieldTotal_OAK->SetLineColor(kGray+1);


  //Loop
  for(unsigned int i = 0; i<ntot; i++) {
    tLB->GetEvent(i);
    tOAK->GetEvent(i);
    
    hZeroLeptonSignalYieldTotal_LB->Fill(zeroLeptonSignalYieldTotal_LB);
    hZeroLeptonSignalYieldTotal_OAK->Fill(zeroLeptonSignalYieldTotal_OAK);
    
    float differenceSignal = zeroLeptonSignalYieldTotal_LB - zeroLeptonSignalYieldTotal_OAK;
    float percentDifferenceSignal =  100*(zeroLeptonSignalYieldTotal_LB - zeroLeptonSignalYieldTotal_OAK)/( 0.5*(zeroLeptonSignalYieldTotal_LB + zeroLeptonSignalYieldTotal_OAK) );
    float differenceTopWJets = zeroLeptonTopWJetsYieldTotal_LB - zeroLeptonTopWJetsYieldTotal_OAK;
    float percentDifferenceTopWJets =  100*(zeroLeptonTopWJetsYieldTotal_LB - zeroLeptonTopWJetsYieldTotal_OAK)/( 0.5*(zeroLeptonTopWJetsYieldTotal_LB + zeroLeptonTopWJetsYieldTotal_OAK) );
    float differenceQCD = zeroLeptonQCDYieldTotal_LB - zeroLeptonQCDYieldTotal_OAK;
    float percentDifferenceQCD =  100*(zeroLeptonQCDYieldTotal_LB - zeroLeptonQCDYieldTotal_OAK)/( 0.5*(zeroLeptonQCDYieldTotal_LB + zeroLeptonQCDYieldTotal_OAK) );
    float differenceZtoNuNu = zeroLeptonZtoNuNuYieldTotal_LB - zeroLeptonZtoNuNuYieldTotal_OAK;
    float percentDifferenceZtoNuNu =  100*(zeroLeptonZtoNuNuYieldTotal_LB - zeroLeptonZtoNuNuYieldTotal_OAK)/( 0.5*(zeroLeptonZtoNuNuYieldTotal_LB + zeroLeptonZtoNuNuYieldTotal_OAK) );


    //cout << zeroLeptonSignalYieldTotal_LB << " " << zeroLeptonSignalYieldTotal_OAK << " " << difference << " " << percentDifference << endl;
    cout << zeroLeptonTopWJetsYieldTotal_LB << " " << zeroLeptonTopWJetsYieldTotal_OAK << endl;
    
    hZeroLeptonSignalYieldTotalDifference->Fill( differenceSignal );
    hZeroLeptonSignalYieldTotalPercentDifference->Fill( percentDifferenceSignal );
    
    hZeroLeptonTopWJetsYieldTotalDifference->Fill( differenceTopWJets );
    hZeroLeptonTopWJetsYieldTotalPercentDifference->Fill( percentDifferenceTopWJets );
    
    hZeroLeptonQCDYieldTotalDifference->Fill( differenceQCD );
    hZeroLeptonQCDYieldTotalPercentDifference->Fill( percentDifferenceQCD );
    
    hZeroLeptonZtoNuNuYieldTotalDifference->Fill( differenceZtoNuNu );
    hZeroLeptonZtoNuNuYieldTotalPercentDifference->Fill( percentDifferenceZtoNuNu );
    
  }
  
  TCanvas * cDifference = new TCanvas("cDifference", "Difference", 2*640, 480);
  cDifference->Divide(2,1);
  cDifference->cd(1);
  hZeroLeptonSignalYieldTotalDifference->Draw();
  cDifference->cd(2);
  hZeroLeptonSignalYieldTotalPercentDifference->Draw();
  cDifference->Print("cDifference.pdf");

  TCanvas * cDifferenceBackgrounds = new TCanvas("cDifferenceBackgrounds", "cDifferenceBackgrounds", 3*640, 480);
  cDifferenceBackgrounds->Divide(3,1);
  cDifferenceBackgrounds->cd(1);
  hZeroLeptonTopWJetsYieldTotalDifference->Draw();
  cDifferenceBackgrounds->cd(2);
  hZeroLeptonQCDYieldTotalDifference->Draw();
  cDifferenceBackgrounds->cd(3);
  hZeroLeptonZtoNuNuYieldTotalDifference->Draw();
  cDifferenceBackgrounds->Print("cDifferenceBackgrounds.pdf");

  TCanvas * cPercentDifferenceBackgrounds = new TCanvas("cPercentDifferenceBackgrounds", "cPercentDifferenceBackgrounds", 3*640, 480);
  cPercentDifferenceBackgrounds->Divide(3,1);
  cPercentDifferenceBackgrounds->cd(1);
  hZeroLeptonTopWJetsYieldTotalPercentDifference->Draw();
  cPercentDifferenceBackgrounds->cd(2);
  hZeroLeptonQCDYieldTotalPercentDifference->Draw();
  cPercentDifferenceBackgrounds->cd(3);
  hZeroLeptonZtoNuNuYieldTotalPercentDifference->Draw();
  cPercentDifferenceBackgrounds->Print("cPercentDifferenceBackgrounds.pdf");

  TCanvas * cDists = new TCanvas("cDists", "Distributions", 640, 480);
  cDists->cd();
  hZeroLeptonSignalYieldTotal_OAK->SetMaximum(50);
  hZeroLeptonSignalYieldTotal_OAK->Draw();
  hZeroLeptonSignalYieldTotal_LB->Draw("SAMES");
  cDists->Print("cDists.pdf");

  


  return;
}
