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

bool skipBinInAnalysis(TString binName) {
  
  return false;

  //FIXME hardcoded
  if(binName=="bin37" || binName=="bin38" || binName=="bin39") return true;
  //skipBin(binName);//from likelihoodBuilder

  //bool partial = true;
  bool partial = false;
  if(partial) {
    //if(!(binName).Contains("3b")) return true;;
    if(!(binName).Contains("2b")) return true;
  }
  
  return false;
}

void extractFromWorkspace(TString workspaceFile = "test.root", TString datFile = "", bool skip = true)
{

  if(skip) {
    if(datFile != "") {
      ofstream myfile;
      myfile.open(datFile.Data(), ios::out | ios::app);
      assert(myfile.is_open());
      myfile << "-99" << " ";
      myfile.close();
    }
    return;
  }
  
  cout <<  "Starting extractFromWorkspace" << endl;
  TFile* wstf = new TFile ( workspaceFile );
  
  RooWorkspace* ws = (RooWorkspace*)wstf->Get("workspace");
  
  RooRealVar* signalUncertainty = ws->var("signalUncertainty"); 
  double sigU = signalUncertainty->getVal();
  
  if(datFile != "") {
    ofstream myfile;
    myfile.open(datFile.Data(), ios::out | ios::app);
    assert(myfile.is_open());
    myfile << sigU << " ";
    myfile.close();
  }

  wstf->Close();

  return;
  
}

void integratedTotals(TString workspaceFile = "test.root", TString binFilesFile = "binFilesFile.dat", TString datFile = "", bool includeTriggerEfficiency = false)
{
  cout << "Starting integratedTotals" << endl;
  
  includeTriggerEfficiency = false;
  TString addName = "";
  if(includeTriggerEfficiency) addName = "Data";

  TFile* wstf = new TFile ( workspaceFile );
  RooWorkspace* ws = (RooWorkspace*)wstf->Get("workspace");
  
  RooRealVar * signalCrossSection = ws->var("signalCrossSection");
  RooRealVar * luminosity = ws->var("luminosity");
  
  ifstream binStream;
  binStream.open(binFilesFile.Data(),fstream::in);
  
  TString index;
  string fileLine;
  std::vector<TString> binNames;
  std::vector<TString> binFiles;  
  while(!binStream.eof()) {
    getline(binStream,fileLine);
    TString thisLine(fileLine.c_str());
    
    TStringToken nameAndNumber(thisLine," ");
    nameAndNumber.NextToken();
    index = nameAndNumber;
    if(index=="") continue;
    
    binNames.push_back(index);
    
    nameAndNumber.NextToken();
    index = nameAndNumber;
    binFiles.push_back(index);
  }
  binStream.close();

  cout << "nbins: " << binNames.size() << endl;
  
  //Now get ttwj qcd znn -- do it by simply incrementing doubles since we don't care about the error for now
  double ttwjtot=0, qcdtot=0, znntot=0, vvtot=0;
  double dattot=0;
  for(unsigned int i = 0; i<binNames.size(); i++) {
    
    if( skipBinInAnalysis(binNames.at(i)) ) continue;
    cout << "binName: " << binNames.at(i) << endl;
    cout << "binFile: " << binFiles.at(i) << endl;
    

    RooAbsReal * ttwj = ws->function("zeroLepton_"+binNames.at(i)+"_TopWJets"+addName+"Yield");
    RooAbsReal * qcd = ws->function("zeroLepton_"+binNames.at(i)+"_QCD"+addName+"Yield");
    RooAbsReal * znn = ws->function("zeroLepton_"+binNames.at(i)+"_ZtoNuNu"+addName+"Yield");
    RooAbsReal * vv = ws->function("zeroLepton_"+binNames.at(i)+"_Diboson"+addName+"Yield");
    
    RooRealVar * dat = ws->var("zeroLepton_"+binNames.at(i)+"_Count");
    
    assert (  (ttwj != NULL) && (qcd != NULL) && (znn !=NULL) && (vv !=NULL) && (dat !=NULL) );

    ttwjtot += ttwj->getVal();
    qcdtot += qcd->getVal();
    znntot += znn->getVal();
    vvtot += vv->getVal();
    dattot += dat->getVal();
    
  }
  
  
  //for signal, we care about the error, so let's do it the Roo way
  RooArgSet* signalYieldSet = new RooArgSet("signalYields");
  for(unsigned int i =0; i<binNames.size(); i++) {
    if( skipBinInAnalysis(binNames.at(i)) ) continue;
    signalYieldSet->add( *(ws->arg("zeroLepton_"+binNames.at(i)+"_Signal"+addName+"Yield")) );
    cout << (ws->function("zeroLepton_"+binNames.at(i)+"_Signal"+addName+"Yield"))->getVal() << endl;
    cout << binNames.at(i) << " signal yield pointer = " << (ws->arg("zeroLepton_"+binNames.at(i)+"_Signal"+addName+"Yield"))  << endl;
  }
  RooAddition* signalYield = new RooAddition("signalYield", "signalYield", *signalYieldSet);
  
  RooFitResult *fitResult = (RooFitResult*)wstf->Get("fitresult_likelihood_data");
  
  cout << "fit result pointer = " << fitResult << endl;

  double sig = signalYield->getVal();
  double sigerr = signalYield->getPropagatedError(*fitResult);
  cout << "analyzeFitOutput: " << sig << "+-" << sigerr << " " << ttwjtot << " " << qcdtot << " " << znntot << " " << vvtot << " " << dattot << endl; 
  
  if(datFile != "") {
    cout << "DEBUG: printing to file" << endl;
    ofstream myfile;
    myfile.open(datFile.Data(), ios::out | ios::app);
    assert(myfile.is_open());
    myfile << sig << " " << sigerr << " " << ttwjtot << " " << qcdtot << " " << znntot << " " << vvtot << " " << dattot << " ";
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



void owenPlots(TString workspaceFile = "test.root", TString binFilesFile = "binFilesFile.dat", bool logy=false, bool ABCD = false) {
  //This function assumes bin dat files are named like bin_Mx_Hy_zb
  bool debug = true;
  if(debug) cout << "Starting owenPlots" << endl;

  gStyle->SetOptStat(0);

  int maxNM=1, maxNH=1, maxNB=1; 
  
  if(debug) cout <<  binFilesFile << endl;
  ifstream binStream;
  binStream.open(binFilesFile.Data(),fstream::in);
  assert(binStream.is_open());

  map<TString,s3d> binMap;
  TString index, fileName;
  string fileLine;
  while(!binStream.eof()) {
    getline(binStream,fileLine);
    TString thisLine(fileLine.c_str());

    TStringToken nameAndNumber(thisLine," ");
    nameAndNumber.NextToken();
    index = nameAndNumber;
    if(debug) cout << "index: " << index << endl;
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
  std::vector<TH1D*> hZL_sig;  std::vector<TH1D*> hZL_ttwj;  std::vector<TH1D*> hZL_qcd;  std::vector<TH1D*> hZL_znn; std::vector<TH1D*> hZL_vv; std::vector<TH1D*> hZL_data;
  std::vector<TH1D*> hSL_sig;  std::vector<TH1D*> hSL_ttwj;  std::vector<TH1D*> hSL_qcd;  std::vector<TH1D*> hSL_znn; std::vector<TH1D*> hSL_vv; std::vector<TH1D*> hSL_data;
  std::vector<TH1D*> hLDP_sig;  std::vector<TH1D*> hLDP_ttwj;  std::vector<TH1D*> hLDP_qcd;  std::vector<TH1D*> hLDP_znn; std::vector<TH1D*> hLDP_vv; std::vector<TH1D*> hLDP_data;

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

    hZL_vv.push_back(new TH1D("ZL_vv_"+b+"b","ZL_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hSL_vv.push_back(new TH1D("SL_vv_"+b+"b","SL_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hLDP_vv.push_back(new TH1D("LDP_vv_"+b+"b","LDP_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hZL_vv[i-1]->SetFillColor(kOrange-3);
    hSL_vv[i-1]->SetFillColor(kOrange-3);
    hLDP_vv[i-1]->SetFillColor(kOrange-3);
 

    hsZL.push_back(new THStack("ZL_"+b+"b","ZL_"+b+"b"));
    hsZL[i-1]->Add(hZL_vv[i-1]);   
    hsZL[i-1]->Add(hZL_znn[i-1]);   
    hsZL[i-1]->Add(hZL_qcd[i-1]);   
    hsZL[i-1]->Add(hZL_ttwj[i-1]);
    hsZL[i-1]->Add(hZL_sig[i-1]);
    
    hsSL.push_back(new THStack("SL_"+b+"b","SL_"+b+"b"));
    hsSL[i-1]->Add(hSL_vv[i-1]);
    hsSL[i-1]->Add(hSL_znn[i-1]);
    hsSL[i-1]->Add(hSL_qcd[i-1]);    
    hsSL[i-1]->Add(hSL_ttwj[i-1]);
    hsSL[i-1]->Add(hSL_sig[i-1]);
    
    hsLDP.push_back(new THStack("LDP_"+b+"b","LDP_"+b+"b"));
    hsLDP[i-1]->Add(hLDP_vv[i-1]);
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
    if(ABCD) hSL_data[bBin-1]->SetBinContent( binIndex,(ws->function( "oneLepton_"+binName+"_Count" ))->getVal() );
    hLDP_data[bBin-1]->SetBinContent( binIndex,(ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_Count" ))->getVal() );

    hZL_sig[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLepton_"+binName+"_SignalDataYield" ))->getVal() );
    if(ABCD) hSL_sig[bBin-1]->SetBinContent( binIndex, (ws->function( "oneLepton_"+binName+"_SignalDataYield" ))->getVal() );
    hLDP_sig[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_SignalDataYield" ))->getVal() );
    cout << (ws->function( "zeroLepton_"+binName+"_SignalDataYield" ))->getVal() << " " << (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_SignalDataYield" ))->getVal() << endl;

    hZL_ttwj[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLepton_"+binName+"_TopWJetsDataYield" ))->getVal() );
    if(ABCD) hSL_ttwj[bBin-1]->SetBinContent( binIndex, (ws->function( "oneLepton_"+binName+"_TopWJetsDataYield" ))->getVal() );
    hLDP_ttwj[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_TopWJetsDataYield" ))->getVal() );
    
    
    hZL_qcd[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLepton_"+binName+"_QCDDataYield" ))->getVal() );
    if(ABCD) hSL_qcd[bBin-1]->SetBinContent( binIndex, 0);
    hLDP_qcd[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_QCDDataYield" ))->getVal() );
    
    hZL_znn[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLepton_"+binName+"_ZtoNuNuDataYield" ))->getVal() );
    if(ABCD) hSL_znn[bBin-1]->SetBinContent( binIndex, 0);
    hLDP_znn[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_ZtoNuNuDataYield" ))->getVal() );

    hZL_vv[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLepton_"+binName+"_DibosonDataYield" ))->getVal() );
    if(ABCD) hSL_vv[bBin-1]->SetBinContent( binIndex, (ws->function( "oneLepton_"+binName+"_DibosonDataYield" ))->getVal() );
    hLDP_vv[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_DibosonDataYield" ))->getVal() );
    
  }
   
  
  for(int i=1; i<=maxNM; i++) {//hardcoded
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

  if(logy) {
    hZee_z->SetMinimum(.1);
    hZmm_z->SetMinimum(.1);

    for(int i = 1; i<=maxNB; i++) {
      hsZL[i-1]->SetMinimum(.1);
      hsSL[i-1]->SetMinimum(.1);
      hsLDP[i-1]->SetMinimum(.1);
    }
  }


  TCanvas * c = new TCanvas("c", "c", 1000, 1000);
  c->Divide(3,4);
  c->cd(1);
  if(hsZL[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsZL[0]->Draw();				
  hZL_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(2);
  if(hsZL[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsZL[1]->Draw();
  hZL_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);  
  c->cd(3);
  if(hsZL[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsZL[2]->Draw();
  hZL_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(4);
  if(hsSL[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsSL[0]->Draw();
  hSL_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);  
  c->cd(5);
  if(hsSL[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsSL[1]->Draw();
  hSL_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(6);
  if(hsSL[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsSL[2]->Draw();
  hSL_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(7);
  if(hsLDP[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLDP[0]->Draw();
  hLDP_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(8);
  if(hsLDP[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLDP[1]->Draw();
  hLDP_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(9);
  if(hsLDP[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLDP[2]->Draw();
  hLDP_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(10);
  if(hZee_z->GetMaximum()>0) gPad->SetLogy(logy);
  hZee_z->Draw();
  hZee_d->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  c->cd(11);
  if(hZmm_z->GetMaximum()>0) gPad->SetLogy(logy);
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

  TString logstring="";
  if(logy) logstring="_log";
  c->SaveAs("stack"+logstring+".pdf");
  
  
  wstf->Close();

}


void analyzeFit(TString workspaceFile, TString binFilesFile, TString datFile) {
  cout << "Starting analyzeFit" << endl;

  
  extractFromWorkspace(workspaceFile, datFile);

  integratedTotals(workspaceFile, binFilesFile, datFile);
  //integratedTotals(workspaceFile);
  //owenPlots();

}
