#include <cassert>
#include <iostream>

#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "TPRegexp.h"
#include "RooAddition.h"
#include "RooProduct.h"

#include "RooStats/ModelConfig.h"


using namespace RooFit ;
using namespace RooStats ; 
using namespace std;


void integratedTotals(TString workspaceFile = "test.root", TString name =  "", TString binFilesFile = "binFilesFile.dat", TString datFile = "")
{
  
  TFile* wstf = new TFile ( workspaceFile );
  
  RooWorkspace* ws = (RooWorkspace*)wstf->Get("workspace");
  
  RooRealVar * signalCrossSection = ws->var("signalCrossSection");
    

  ifstream binStream;
  binStream.open(binFilesFile.Data(),fstream::in);


  TString index;
  string fileLine;
  std::vector<TString> binNames;
  while(!binStream.eof()) {
    getline(binStream,fileLine);
    TString thisLine(fileLine.c_str());

    TStringToken nameAndNumber(thisLine," ");
    nameAndNumber.NextToken();
    index = nameAndNumber;
    if(index=="") continue;
    
    binNames.push_back(index);
  }
  binStream.close();

  cout << "nbins: " << binNames.size() << endl;
  
  //Now get ttwj qcd znn -- do it by simply incrementing doubles since we don't care about the error for now
  double ttwjtot=0, qcdtot=0, znntot=0;
  for(unsigned int i = 0; i<binNames.size(); i++) {
    
    RooAbsReal * ttwj = ws->function("zeroLepton_"+binNames.at(i)+"_TopWJetsYield");
    RooAbsReal * qcd = ws->function("zeroLepton_"+binNames.at(i)+"_QCDYield");
    RooAbsReal * znn = ws->function("zeroLepton_"+binNames.at(i)+"_ZtoNuNuYield");
    
    assert (  (ttwj != NULL) && (qcd != NULL) && (znn !=NULL) );

    ttwjtot += ttwj->getVal();
    qcdtot += qcd->getVal();
    znntot += znn->getVal();
    
  }


  //for signal, we care about the error, so let's do it the Roo way
  RooArgSet* signalFractionSet = new RooArgSet("signalFractions");
  for(unsigned int i =0; i<binNames.size(); i++) {
    signalFractionSet->add( *(ws->arg("zeroLeptonSignalYieldFraction_"+binNames.at(i)+"_BetaPrimeInverseCDF")) );
  }
  RooAddition* signalFractionAddition = new RooAddition("signalFractionAddition", "signalFractionAddition", *signalFractionSet);
  RooProduct*  signalYield = new RooProduct("signalYield","signalYield", RooArgSet(*signalFractionAddition,*signalCrossSection));;;

  RooFitResult *fitResult = (RooFitResult*)wstf->Get("fitresult_likelihood_data");
  
  double sig = signalYield->getVal();
  double sigerr = signalYield->getPropagatedError(*fitResult);
  cout << "analyzeFitOutput: " << name << " " << sig << "+-" << sigerr << " " << ttwjtot << " " << qcdtot << " " << znntot << endl; 
  cout << "DEBUG: datFile: " << datFile << endl;
  
  if(datFile != "") {
    cout << "DEBUG: printing to file" << endl;
    ofstream myfile;
    myfile.open(datFile.Data(), ios::out | ios::app);
    assert(myfile.is_open());
    myfile << sig << " " << sigerr << " " << ttwjtot << " " << qcdtot << " " << znntot << " ";
    myfile.close();
  }


  wstf->Close();

  return;
}


/*
  //this version takes the yield fractions  from the input file (hardcoded  in scale object)
void integratedTotals(TString workspaceFile = "test.root", TString name = "")
{
  
  TFile* wstf = new TFile ( workspaceFile );
  
  RooWorkspace* ws = (RooWorkspace*)wstf->Get("workspace");
  
  RooRealVar * signalCrossSection = ws->var("signalCrossSection");
    
  //Now get ttwj qcd znn
  double ttwjtot=0, qcdtot=0, znntot=0;
  for(int i = 1; i<=27; i++) {
    
    TString num = "";
    num+=i;
    
    RooAbsReal * ttwj = ws->function("zeroLepton_bin"+num+"_TopWJetsYield");
    RooAbsReal * qcd = ws->function("zeroLepton_bin"+num+"_QCDYield");
    RooAbsReal * znn = ws->function("zeroLepton_bin"+num+"_ZtoNuNuYield");
    
    assert (  (ttwj != NULL) && (qcd != NULL) && (znn !=NULL) );

    ttwjtot += ttwj->getVal();
    qcdtot += qcd->getVal();
    znntot += znn->getVal();
    
  }

  double scale = (5.0*1934.0)/10000.0 ;
  
  double sig = scale*signalCrossSection->getVal();
  double sigerr = scale*signalCrossSection->getError();

  cout  << "Susy cross section scaled by " <<  scale << endl;
  cout << "analyzeFitOutput: " << name << " " << sig << "+-" << sigerr << " " << ttwjtot << " " << qcdtot << " " << znntot << endl; 


  wstf->Close();

  return;
}
*/

struct s3d {
  int nM;
  int nH;
  int nB;
};



void owenPlots(TString workspaceFile = "test.root", TString name = "",  TString binFilesFile = "binFilesFile.dat") {
  //This function assumes bin dat files are named like bin_Mx_Hy_zb

  gStyle->SetOptStat(0);

  int maxNM=1, maxNH=1, maxNB=1; 
  
  ifstream binStream;
  binStream.open(binFilesFile.Data(),fstream::in);
  
  map<TString,s3d> binMap;
  TString index, fileName;
  string fileLine;
  while(!binStream.eof()) {
    getline(binStream,fileLine);
    TString thisLine(fileLine.c_str());

    TStringToken nameAndNumber(thisLine," ");
    nameAndNumber.NextToken();
    index = nameAndNumber;
    if(index=="") continue;
    nameAndNumber.NextToken();
    fileName = nameAndNumber;
    fileName.Remove( fileName.First('.'), fileName.Length() - fileName.First('.') );
    
    TObjArray *subStrL = TPRegexp("(bin_M)([1-99])(_H)([1-99])(_)([1-99])(b)").MatchS(fileName);
    s3d binInfo; 
    binInfo.nM = ( ((TObjString *)subStrL->At(2))->GetString() ).Atoi();
    binInfo.nH = ( ((TObjString *)subStrL->At(4))->GetString() ).Atoi();
    binInfo.nB = ( ((TObjString *)subStrL->At(6))->GetString() ).Atoi();
    
    if(binInfo.nM > maxNM) maxNM = binInfo.nM;
    if(binInfo.nH > maxNH) maxNH = binInfo.nH;
    if(binInfo.nB > maxNB) maxNB = binInfo.nB;

    binMap[index] = binInfo;
  }
  binStream.close();
  cout << "Total number of bins: " << (int)binMap.size() << endl;
  cout << "Total number of bins in MET=" << maxNM << ", in HT=" << maxNH << ", in BTAGS=" << maxNB << endl; 

  TFile* wstf = new TFile ( workspaceFile );
  
  RooWorkspace* ws = (RooWorkspace*)wstf->Get("workspace");
  
  //each entry will be btag multiplicity
  std::vector<TH1D*> hZL_sig;  std::vector<TH1D*> hZL_ttwj;  std::vector<TH1D*> hZL_qcd;  std::vector<TH1D*> hZL_znn; std::vector<TH1D*> hZL_data;
  std::vector<TH1D*> hSL_sig;  std::vector<TH1D*> hSL_ttwj;  std::vector<TH1D*> hSL_qcd;  std::vector<TH1D*> hSL_znn;  std::vector<TH1D*> hSL_data;
  std::vector<TH1D*> hLDP_sig;  std::vector<TH1D*> hLDP_ttwj;  std::vector<TH1D*> hLDP_qcd;  std::vector<TH1D*> hLDP_znn; std::vector<TH1D*> hLDP_data;

  std::vector<THStack*> hsZL;  std::vector<THStack*> hsSL;   std::vector<THStack*> hsLDP; 


  int nbins = maxNM*(maxNH+1)+1;

  TH1D* hZee_z = new TH1D("Zee", "Zee", nbins, 0.5, nbins+0.5); //assume no btag dependence
  TH1D* hZmm_z = new TH1D("Zmm", "Zmm", nbins, 0.5, nbins+0.5); //assume no btag dependence
  hZee_z->SetFillColor(kGreen);
  hZmm_z->SetFillColor(kGreen);
  TH1D* hZee_d = new TH1D("Zee_data", "Zee_data", nbins, 0.5, nbins+0.5); //assume no btag dependence
  TH1D* hZmm_d = new TH1D("Zmm_data", "Zmm_data", nbins, 0.5, nbins+0.5); //assume no btag dependence
  hZee_d->SetMarkerStyle(20);
  hZmm_d->SetMarkerStyle(20);

  for(int i = 1; i<=maxNB; i++) {
    TString b = "";
    b+=i;

    hZL_data.push_back(new TH1D("ZL_data_"+b+"b","ZL_data_"+b+"b",nbins,0.5,nbins+0.5));
    hSL_data.push_back(new TH1D("SL_data_"+b+"b","SL_data_"+b+"b",nbins,0.5,nbins+0.5));
    hLDP_data.push_back(new TH1D("LDP_data_"+b+"b","LDP_data_"+b+"b",nbins,0.5,nbins+0.5));
    hZL_data[i-1]->SetMarkerStyle(20);
    hSL_data[i-1]->SetMarkerStyle(20);
    hLDP_data[i-1]->SetMarkerStyle(20);

    hZL_sig.push_back(new TH1D("ZL_sig_"+b+"b","ZL_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hSL_sig.push_back(new TH1D("SL_sig_"+b+"b","SL_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hLDP_sig.push_back(new TH1D("LDP_sig_"+b+"b","LDP_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hZL_sig[i-1]->SetFillColor(6);
    hSL_sig[i-1]->SetFillColor(6);
    hLDP_sig[i-1]->SetFillColor(6);
    
    hZL_ttwj.push_back(new TH1D("ZL_ttwj_"+b+"b","ZL_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hSL_ttwj.push_back(new TH1D("SL_ttwj_"+b+"b","SL_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hLDP_ttwj.push_back(new TH1D("LDP_ttwj_"+b+"b","LDP_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hZL_ttwj[i-1]->SetFillColor(kBlue-9);
    hSL_ttwj[i-1]->SetFillColor(kBlue-9);
    hLDP_ttwj[i-1]->SetFillColor(kBlue-9);    

    hZL_qcd.push_back(new TH1D("ZL_qcd_"+b+"b","ZL_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hSL_qcd.push_back(new TH1D("SL_qcd_"+b+"b","SL_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hLDP_qcd.push_back(new TH1D("LDP_qcd_"+b+"b","LDP_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hZL_qcd[i-1]->SetFillColor(2);
    hSL_qcd[i-1]->SetFillColor(2);
    hLDP_qcd[i-1]->SetFillColor(2);    

    hZL_znn.push_back(new TH1D("ZL_znn_"+b+"b","ZL_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hSL_znn.push_back(new TH1D("SL_znn_"+b+"b","SL_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hLDP_znn.push_back(new TH1D("LDP_znn_"+b+"b","LDP_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hZL_znn[i-1]->SetFillColor(kGreen-3);
    hSL_znn[i-1]->SetFillColor(kGreen-3);
    hLDP_znn[i-1]->SetFillColor(kGreen-3);
 

    hsZL.push_back(new THStack("ZL_"+b+"b","ZL_"+b+"b"));
    hsZL[i-1]->Add(hZL_znn[i-1]);   
    hsZL[i-1]->Add(hZL_qcd[i-1]);
    hsZL[i-1]->Add(hZL_ttwj[i-1]);
    hsZL[i-1]->Add(hZL_sig[i-1]);
    
    hsSL.push_back(new THStack("SL_"+b+"b","SL_"+b+"b"));
    hsSL[i-1]->Add(hSL_znn[i-1]);
    hsSL[i-1]->Add(hSL_qcd[i-1]);    
    hsSL[i-1]->Add(hSL_ttwj[i-1]);
    hsSL[i-1]->Add(hSL_sig[i-1]);
    
    hsLDP.push_back(new THStack("LDP_"+b+"b","LDP_"+b+"b"));
    hsLDP[i-1]->Add(hLDP_znn[i-1]);
    hsLDP[i-1]->Add(hLDP_qcd[i-1]);
    hsLDP[i-1]->Add(hLDP_ttwj[i-1]);
    hsLDP[i-1]->Add(hLDP_sig[i-1]);
    
  }


  //take data out of workspace and fill histograms
  for(map<TString,s3d>::iterator it = binMap.begin(); it != binMap.end(); it++) {
    TString binName = it->first;
    s3d binInfo = it->second;
    
    cout << binName << " M" << binInfo.nM << " H" << binInfo.nH << " " << binInfo.nB << "b" << endl; 
    int binIndex = 1 + (maxNH+1)*(binInfo.nM-1) + binInfo.nH;
    int bBin = binInfo.nB;
    
    hZL_data[bBin-1]->SetBinContent( binIndex,(ws->function( "zeroLepton_"+binName+"_Count" ))->getVal() );
    hSL_data[bBin-1]->SetBinContent( binIndex,(ws->function( "oneMuon_"+binName+"_Count" ))->getVal() );//BEN FIXME -- just muon now!
    hLDP_data[bBin-1]->SetBinContent( binIndex,(ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_Count" ))->getVal() );

    hZL_sig[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLepton_"+binName+"_SignalYield" ))->getVal() );
    hSL_sig[bBin-1]->SetBinContent( binIndex, (ws->function( "oneLepton_"+binName+"_SignalYield" ))->getVal() );
    hLDP_sig[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_SignalYield" ))->getVal() );
    
    hZL_ttwj[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLepton_"+binName+"_TopWJetsYield" ))->getVal() );
    hSL_ttwj[bBin-1]->SetBinContent( binIndex, (ws->function( "oneLepton_"+binName+"_TopWJetsYield" ))->getVal() );
    hLDP_ttwj[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_TopWJetsYield" ))->getVal() );
    
    hZL_qcd[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLepton_"+binName+"_QCDYield" ))->getVal() );
    hSL_qcd[bBin-1]->SetBinContent( binIndex, 0);
    hLDP_qcd[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_QCDYield" ))->getVal() );
    
    hZL_znn[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLepton_"+binName+"_ZtoNuNuYield" ))->getVal() );
    hSL_znn[bBin-1]->SetBinContent( binIndex, 0);
    hLDP_znn[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_ZtoNuNuYield" ))->getVal() );
    
  }
   

  for(int i=1; i<=maxNM; i++) {
    for(int j=1; j<=maxNH; j++) {
      TString binLabel = "M"; binLabel+=i; binLabel+="_H"; binLabel+=j;
      int binIndex = 1 + (maxNH+1)*(i-1) + j;
      
      hZee_z->SetBinContent( binIndex, (ws->function( "diElectron_"+binLabel+"_Yield" ))->getVal() );
      hZee_z->GetXaxis()->SetBinLabel( binIndex, binLabel);
      hZmm_z->SetBinContent( binIndex, (ws->function( "diMuon_"+binLabel+"_Yield" ))->getVal() );
      hZmm_z->GetXaxis()->SetBinLabel( binIndex, binLabel);

      hZee_d->SetBinContent( binIndex, (ws->function( "diElectron_"+binLabel+"_Count" ))->getVal() );
      hZmm_d->SetBinContent( binIndex, (ws->function( "diMuon_"+binLabel+"_Count" ))->getVal() );
      
    }
  }


  double maxS = 1.25;
  hZee_z->SetMaximum(maxS*hZee_d->GetMaximum());
  hZmm_z->SetMaximum(maxS*hZmm_d->GetMaximum());

  for(int i = 1; i<=maxNB; i++) {
    hsZL[i-1]->SetMaximum(maxS*hZL_data[i-1]->GetMaximum());
    hsSL[i-1]->SetMaximum(maxS*hSL_data[i-1]->GetMaximum());
    hsLDP[i-1]->SetMaximum(maxS*hLDP_data[i-1]->GetMaximum());
  }


  TCanvas * c = new TCanvas("c", "c", 1000, 1000);
  c->Divide(3,4);
  c->cd(1);
  hsZL[0]->Draw();
  hZL_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(2);
  hsZL[1]->Draw();
  hZL_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);  
  c->cd(3);
  hsZL[2]->Draw();
  hZL_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(4);
  hsSL[0]->Draw();
  hSL_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);  
  c->cd(5);
  hsSL[1]->Draw();
  hSL_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(6);
  hsSL[2]->Draw();
  hSL_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(7);
  hsLDP[0]->Draw();
  hLDP_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(8);
  hsLDP[1]->Draw();
  hLDP_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(9);
  hsLDP[2]->Draw();
  hLDP_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(10);
  hZee_z->Draw();
  hZee_d->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(11);
  hZmm_z->Draw();
  hZmm_d->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);


  //label axes
  for(map<TString,s3d>::iterator it = binMap.begin(); it != binMap.end(); it++) {
    s3d binInfo = it->second;
    TString binLabel = "M"; binLabel+=binInfo.nM; binLabel+="_H"; binLabel+=binInfo.nH; binLabel+="_"; binLabel+=binInfo.nB; binLabel+="b";
    int binIndex = 1 + (maxNH+1)*(binInfo.nM-1) + binInfo.nH;
    int bBin = binInfo.nB;
    
    hsZL[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
    hsSL[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
    hsLDP[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
      
  }

  for(int i=1; i<=maxNB; i++) {
    hsZL[i-1]->GetXaxis()->LabelsOption("v");
    hsSL[i-1]->GetXaxis()->LabelsOption("v");
    hsLDP[i-1]->GetXaxis()->LabelsOption("v");
  }
  hZee_z->GetXaxis()->LabelsOption("v");
  hZmm_z->GetXaxis()->LabelsOption("v");



  gPad->Modified();


  c->SaveAs("test.pdf");
  
  
  wstf->Close();

}


void analyzeFit(TString workspaceFile = "test.root", TString name = "", TString binFilesFile = "", TString datFile= "") {

  integratedTotals(workspaceFile, name, binFilesFile, datFile);
  //integratedTotals(workspaceFile, name);
  //owenPlots();

}
