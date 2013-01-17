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
#include "RooAbsReal.h"
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

bool skipMHBin(int met, int ht)
{
  if (met==4 && ht==1) return true;
  return false;
}

void extractFromWorkspace(TString workspaceFile = "test.root", TString datFile = "", bool skip = false)
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

  //For MR
  std::vector<TH1D*> hLL1_sig;  std::vector<TH1D*> hLL1_ttwj;  std::vector<TH1D*> hLL1_qcd;  std::vector<TH1D*> hLL1_znn; std::vector<TH1D*> hLL1_vv; std::vector<TH1D*> hLL1_data;
  std::vector<TH1D*> hLL2_sig;  std::vector<TH1D*> hLL2_ttwj;  std::vector<TH1D*> hLL2_qcd;  std::vector<TH1D*> hLL2_znn; std::vector<TH1D*> hLL2_vv; std::vector<TH1D*> hLL2_data;
  std::vector<TH1D*> hLL3_sig;  std::vector<TH1D*> hLL3_ttwj;  std::vector<TH1D*> hLL3_qcd;  std::vector<TH1D*> hLL3_znn; std::vector<TH1D*> hLL3_vv; std::vector<TH1D*> hLL3_data;
  std::vector<TH1D*> hLL4_sig;  std::vector<TH1D*> hLL4_ttwj;  std::vector<TH1D*> hLL4_qcd;  std::vector<TH1D*> hLL4_znn; std::vector<TH1D*> hLL4_vv; std::vector<TH1D*> hLL4_data;
  std::vector<TH1D*> hLL5_sig;  std::vector<TH1D*> hLL5_ttwj;  std::vector<TH1D*> hLL5_qcd;  std::vector<TH1D*> hLL5_znn; std::vector<TH1D*> hLL5_vv; std::vector<TH1D*> hLL5_data;
  std::vector<TH1D*> hTM1_sig;  std::vector<TH1D*> hTM1_ttwj;  std::vector<TH1D*> hTM1_qcd;  std::vector<TH1D*> hTM1_znn; std::vector<TH1D*> hTM1_vv; std::vector<TH1D*> hTM1_data;
  std::vector<TH1D*> hTM2_sig;  std::vector<TH1D*> hTM2_ttwj;  std::vector<TH1D*> hTM2_qcd;  std::vector<TH1D*> hTM2_znn; std::vector<TH1D*> hTM2_vv; std::vector<TH1D*> hTM2_data;
  std::vector<TH1D*> hTM3_sig;  std::vector<TH1D*> hTM3_ttwj;  std::vector<TH1D*> hTM3_qcd;  std::vector<TH1D*> hTM3_znn; std::vector<TH1D*> hTM3_vv; std::vector<TH1D*> hTM3_data;
  std::vector<TH1D*> hTM4_sig;  std::vector<TH1D*> hTM4_ttwj;  std::vector<TH1D*> hTM4_qcd;  std::vector<TH1D*> hTM4_znn; std::vector<TH1D*> hTM4_vv; std::vector<TH1D*> hTM4_data;
  std::vector<TH1D*> hTM5_sig;  std::vector<TH1D*> hTM5_ttwj;  std::vector<TH1D*> hTM5_qcd;  std::vector<TH1D*> hTM5_znn; std::vector<TH1D*> hTM5_vv; std::vector<TH1D*> hTM5_data;
  std::vector<TH1D*> h2LL_sig;  std::vector<TH1D*> h2LL_ttwj;  std::vector<TH1D*> h2LL_qcd;  std::vector<TH1D*> h2LL_znn; std::vector<TH1D*> h2LL_vv; std::vector<TH1D*> h2LL_data;
  std::vector<TH1D*> h2TM_sig;  std::vector<TH1D*> h2TM_ttwj;  std::vector<TH1D*> h2TM_qcd;  std::vector<TH1D*> h2TM_znn; std::vector<TH1D*> h2TM_vv; std::vector<TH1D*> h2TM_data;

  std::vector<THStack*> hsLL1;  std::vector<THStack*> hsLL2;   std::vector<THStack*> hsLL3;  std::vector<THStack*> hsLL4;  std::vector<THStack*> hsLL5; 
  std::vector<THStack*> hsTM1;  std::vector<THStack*> hsTM2;   std::vector<THStack*> hsTM3;  std::vector<THStack*> hsTM4;  std::vector<THStack*> hsTM5; 
  std::vector<THStack*> hs2LL;  std::vector<THStack*> hs2TM;   

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

    hLL1_data.push_back(new TH1D("LL1_data_"+b+"b","LL1_data_"+b+"b",nbins,0.5,nbins+0.5));
    hLL2_data.push_back(new TH1D("LL2_data_"+b+"b","LL2_data_"+b+"b",nbins,0.5,nbins+0.5));
    hLL3_data.push_back(new TH1D("LL3_data_"+b+"b","LL3_data_"+b+"b",nbins,0.5,nbins+0.5));
    hLL4_data.push_back(new TH1D("LL4_data_"+b+"b","LL4_data_"+b+"b",nbins,0.5,nbins+0.5));
    hLL5_data.push_back(new TH1D("LL5_data_"+b+"b","LL5_data_"+b+"b",nbins,0.5,nbins+0.5));
    hTM1_data.push_back(new TH1D("TM1_data_"+b+"b","TM1_data_"+b+"b",nbins,0.5,nbins+0.5));
    hTM2_data.push_back(new TH1D("TM2_data_"+b+"b","TM2_data_"+b+"b",nbins,0.5,nbins+0.5));
    hTM3_data.push_back(new TH1D("TM3_data_"+b+"b","TM3_data_"+b+"b",nbins,0.5,nbins+0.5));
    hTM4_data.push_back(new TH1D("TM4_data_"+b+"b","TM4_data_"+b+"b",nbins,0.5,nbins+0.5));
    hTM5_data.push_back(new TH1D("TM5_data_"+b+"b","TM5_data_"+b+"b",nbins,0.5,nbins+0.5));
    h2LL_data.push_back(new TH1D("2LL_data_"+b+"b","2LL_data_"+b+"b",nbins,0.5,nbins+0.5));
    h2TM_data.push_back(new TH1D("2TM_data_"+b+"b","2TM_data_"+b+"b",nbins,0.5,nbins+0.5));
    hLL1_data[i-1]->SetMarkerStyle(20);
    hLL2_data[i-1]->SetMarkerStyle(20);
    hLL3_data[i-1]->SetMarkerStyle(20);
    hLL4_data[i-1]->SetMarkerStyle(20);
    hLL5_data[i-1]->SetMarkerStyle(20);
    hTM1_data[i-1]->SetMarkerStyle(20);
    hTM2_data[i-1]->SetMarkerStyle(20);
    hTM3_data[i-1]->SetMarkerStyle(20);
    hTM4_data[i-1]->SetMarkerStyle(20);
    hTM5_data[i-1]->SetMarkerStyle(20);
    h2LL_data[i-1]->SetMarkerStyle(20);
    h2TM_data[i-1]->SetMarkerStyle(20);

    hZL_sig.push_back(new TH1D("ZL_sig_"+b+"b","ZL_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hSL_sig.push_back(new TH1D("SL_sig_"+b+"b","SL_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hLDP_sig.push_back(new TH1D("LDP_sig_"+b+"b","LDP_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hZL_sig[i-1]->SetFillColor(6);
    hSL_sig[i-1]->SetFillColor(6);
    hLDP_sig[i-1]->SetFillColor(6);
    
    hLL1_sig.push_back(new TH1D("LL1_sig_"+b+"b","LL1_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hLL2_sig.push_back(new TH1D("LL2_sig_"+b+"b","LL2_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hLL3_sig.push_back(new TH1D("LL3_sig_"+b+"b","LL3_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hLL4_sig.push_back(new TH1D("LL4_sig_"+b+"b","LL4_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hLL5_sig.push_back(new TH1D("LL5_sig_"+b+"b","LL5_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hTM1_sig.push_back(new TH1D("TM1_sig_"+b+"b","TM1_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hTM2_sig.push_back(new TH1D("TM2_sig_"+b+"b","TM2_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hTM3_sig.push_back(new TH1D("TM3_sig_"+b+"b","TM3_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hTM4_sig.push_back(new TH1D("TM4_sig_"+b+"b","TM4_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hTM5_sig.push_back(new TH1D("TM5_sig_"+b+"b","TM5_sig_"+b+"b",nbins,0.5,nbins+0.5));
    h2LL_sig.push_back(new TH1D("2LL_sig_"+b+"b","2LL_sig_"+b+"b",nbins,0.5,nbins+0.5));
    h2TM_sig.push_back(new TH1D("2TM_sig_"+b+"b","2TM_sig_"+b+"b",nbins,0.5,nbins+0.5));
    hLL1_sig[i-1]->SetFillColor(6);
    hLL2_sig[i-1]->SetFillColor(6);
    hLL3_sig[i-1]->SetFillColor(6);
    hLL4_sig[i-1]->SetFillColor(6);
    hLL5_sig[i-1]->SetFillColor(6);
    hTM1_sig[i-1]->SetFillColor(6);
    hTM2_sig[i-1]->SetFillColor(6);
    hTM3_sig[i-1]->SetFillColor(6);
    hTM4_sig[i-1]->SetFillColor(6);
    hTM5_sig[i-1]->SetFillColor(6);
    h2LL_sig[i-1]->SetFillColor(6);
    h2TM_sig[i-1]->SetFillColor(6);

    hZL_ttwj.push_back(new TH1D("ZL_ttwj_"+b+"b","ZL_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hSL_ttwj.push_back(new TH1D("SL_ttwj_"+b+"b","SL_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hLDP_ttwj.push_back(new TH1D("LDP_ttwj_"+b+"b","LDP_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hZL_ttwj[i-1]->SetFillColor(kBlue-9);
    hSL_ttwj[i-1]->SetFillColor(kBlue-9);
    hLDP_ttwj[i-1]->SetFillColor(kBlue-9);    

    hLL1_ttwj.push_back(new TH1D("LL1_ttwj_"+b+"b","LL1_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hLL2_ttwj.push_back(new TH1D("LL2_ttwj_"+b+"b","LL2_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hLL3_ttwj.push_back(new TH1D("LL3_ttwj_"+b+"b","LL3_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hLL4_ttwj.push_back(new TH1D("LL4_ttwj_"+b+"b","LL4_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hLL5_ttwj.push_back(new TH1D("LL5_ttwj_"+b+"b","LL5_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hTM1_ttwj.push_back(new TH1D("TM1_ttwj_"+b+"b","TM1_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hTM2_ttwj.push_back(new TH1D("TM2_ttwj_"+b+"b","TM2_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hTM3_ttwj.push_back(new TH1D("TM3_ttwj_"+b+"b","TM3_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hTM4_ttwj.push_back(new TH1D("TM4_ttwj_"+b+"b","TM4_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hTM5_ttwj.push_back(new TH1D("TM5_ttwj_"+b+"b","TM5_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    h2LL_ttwj.push_back(new TH1D("2LL_ttwj_"+b+"b","2LL_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    h2TM_ttwj.push_back(new TH1D("2TM_ttwj_"+b+"b","2TM_ttwj_"+b+"b",nbins,0.5,nbins+0.5));
    hLL1_ttwj[i-1]->SetFillColor(kBlue-9);
    hLL2_ttwj[i-1]->SetFillColor(kBlue-9);
    hLL3_ttwj[i-1]->SetFillColor(kBlue-9);
    hLL4_ttwj[i-1]->SetFillColor(kBlue-9);
    hLL5_ttwj[i-1]->SetFillColor(kBlue-9);
    hTM1_ttwj[i-1]->SetFillColor(kBlue-9);
    hTM2_ttwj[i-1]->SetFillColor(kBlue-9);
    hTM3_ttwj[i-1]->SetFillColor(kBlue-9);
    hTM4_ttwj[i-1]->SetFillColor(kBlue-9);
    hTM5_ttwj[i-1]->SetFillColor(kBlue-9);
    h2LL_ttwj[i-1]->SetFillColor(kBlue-9);
    h2TM_ttwj[i-1]->SetFillColor(kBlue-9);

    hZL_qcd.push_back(new TH1D("ZL_qcd_"+b+"b","ZL_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hSL_qcd.push_back(new TH1D("SL_qcd_"+b+"b","SL_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hLDP_qcd.push_back(new TH1D("LDP_qcd_"+b+"b","LDP_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hZL_qcd[i-1]->SetFillColor(2);
    hSL_qcd[i-1]->SetFillColor(2);
    hLDP_qcd[i-1]->SetFillColor(2);    

    hLL1_qcd.push_back(new TH1D("LL1_qcd_"+b+"b","LL1_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hLL2_qcd.push_back(new TH1D("LL2_qcd_"+b+"b","LL2_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hLL3_qcd.push_back(new TH1D("LL3_qcd_"+b+"b","LL3_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hLL4_qcd.push_back(new TH1D("LL4_qcd_"+b+"b","LL4_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hLL5_qcd.push_back(new TH1D("LL5_qcd_"+b+"b","LL5_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hTM1_qcd.push_back(new TH1D("TM1_qcd_"+b+"b","TM1_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hTM2_qcd.push_back(new TH1D("TM2_qcd_"+b+"b","TM2_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hTM3_qcd.push_back(new TH1D("TM3_qcd_"+b+"b","TM3_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hTM4_qcd.push_back(new TH1D("TM4_qcd_"+b+"b","TM4_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hTM5_qcd.push_back(new TH1D("TM5_qcd_"+b+"b","TM5_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    h2LL_qcd.push_back(new TH1D("2LL_qcd_"+b+"b","2LL_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    h2TM_qcd.push_back(new TH1D("2TM_qcd_"+b+"b","2TM_qcd_"+b+"b",nbins,0.5,nbins+0.5));
    hLL1_qcd[i-1]->SetFillColor(2);
    hLL2_qcd[i-1]->SetFillColor(2);
    hLL3_qcd[i-1]->SetFillColor(2);
    hLL4_qcd[i-1]->SetFillColor(2);
    hLL5_qcd[i-1]->SetFillColor(2);
    hTM1_qcd[i-1]->SetFillColor(2);
    hTM2_qcd[i-1]->SetFillColor(2);
    hTM3_qcd[i-1]->SetFillColor(2);
    hTM4_qcd[i-1]->SetFillColor(2);
    hTM5_qcd[i-1]->SetFillColor(2);
    h2LL_qcd[i-1]->SetFillColor(2);
    h2TM_qcd[i-1]->SetFillColor(2);

    hZL_znn.push_back(new TH1D("ZL_znn_"+b+"b","ZL_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hSL_znn.push_back(new TH1D("SL_znn_"+b+"b","SL_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hLDP_znn.push_back(new TH1D("LDP_znn_"+b+"b","LDP_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hZL_znn[i-1]->SetFillColor(kGreen-3);
    hSL_znn[i-1]->SetFillColor(kGreen-3);
    hLDP_znn[i-1]->SetFillColor(kGreen-3);

    hLL1_znn.push_back(new TH1D("LL1_znn_"+b+"b","LL1_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hLL2_znn.push_back(new TH1D("LL2_znn_"+b+"b","LL2_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hLL3_znn.push_back(new TH1D("LL3_znn_"+b+"b","LL3_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hLL4_znn.push_back(new TH1D("LL4_znn_"+b+"b","LL4_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hLL5_znn.push_back(new TH1D("LL5_znn_"+b+"b","LL5_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hTM1_znn.push_back(new TH1D("TM1_znn_"+b+"b","TM1_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hTM2_znn.push_back(new TH1D("TM2_znn_"+b+"b","TM2_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hTM3_znn.push_back(new TH1D("TM3_znn_"+b+"b","TM3_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hTM4_znn.push_back(new TH1D("TM4_znn_"+b+"b","TM4_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hTM5_znn.push_back(new TH1D("TM5_znn_"+b+"b","TM5_znn_"+b+"b",nbins,0.5,nbins+0.5));
    h2LL_znn.push_back(new TH1D("2LL_znn_"+b+"b","2LL_znn_"+b+"b",nbins,0.5,nbins+0.5));
    h2TM_znn.push_back(new TH1D("2TM_znn_"+b+"b","2TM_znn_"+b+"b",nbins,0.5,nbins+0.5));
    hLL1_znn[i-1]->SetFillColor(kGreen-3);
    hLL2_znn[i-1]->SetFillColor(kGreen-3);
    hLL3_znn[i-1]->SetFillColor(kGreen-3);
    hLL4_znn[i-1]->SetFillColor(kGreen-3);
    hLL5_znn[i-1]->SetFillColor(kGreen-3);
    hTM1_znn[i-1]->SetFillColor(kGreen-3);
    hTM2_znn[i-1]->SetFillColor(kGreen-3);
    hTM3_znn[i-1]->SetFillColor(kGreen-3);
    hTM4_znn[i-1]->SetFillColor(kGreen-3);
    hTM5_znn[i-1]->SetFillColor(kGreen-3);
    h2LL_znn[i-1]->SetFillColor(kGreen-3);
    h2TM_znn[i-1]->SetFillColor(kGreen-3);

    hZL_vv.push_back(new TH1D("ZL_vv_"+b+"b","ZL_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hSL_vv.push_back(new TH1D("SL_vv_"+b+"b","SL_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hLDP_vv.push_back(new TH1D("LDP_vv_"+b+"b","LDP_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hZL_vv[i-1]->SetFillColor(kOrange-3);
    hSL_vv[i-1]->SetFillColor(kOrange-3);
    hLDP_vv[i-1]->SetFillColor(kOrange-3);
 
    hLL1_vv.push_back(new TH1D("LL1_vv_"+b+"b","LL1_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hLL2_vv.push_back(new TH1D("LL2_vv_"+b+"b","LL2_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hLL3_vv.push_back(new TH1D("LL3_vv_"+b+"b","LL3_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hLL4_vv.push_back(new TH1D("LL4_vv_"+b+"b","LL4_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hLL5_vv.push_back(new TH1D("LL5_vv_"+b+"b","LL5_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hTM1_vv.push_back(new TH1D("TM1_vv_"+b+"b","TM1_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hTM2_vv.push_back(new TH1D("TM2_vv_"+b+"b","TM2_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hTM3_vv.push_back(new TH1D("TM3_vv_"+b+"b","TM3_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hTM4_vv.push_back(new TH1D("TM4_vv_"+b+"b","TM4_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hTM5_vv.push_back(new TH1D("TM5_vv_"+b+"b","TM5_vv_"+b+"b",nbins,0.5,nbins+0.5));
    h2LL_vv.push_back(new TH1D("2LL_vv_"+b+"b","2LL_vv_"+b+"b",nbins,0.5,nbins+0.5));
    h2TM_vv.push_back(new TH1D("2TM_vv_"+b+"b","2TM_vv_"+b+"b",nbins,0.5,nbins+0.5));
    hLL1_vv[i-1]->SetFillColor(kOrange-3);
    hLL2_vv[i-1]->SetFillColor(kOrange-3);
    hLL3_vv[i-1]->SetFillColor(kOrange-3);
    hLL4_vv[i-1]->SetFillColor(kOrange-3);
    hLL5_vv[i-1]->SetFillColor(kOrange-3);
    hTM1_vv[i-1]->SetFillColor(kOrange-3);
    hTM2_vv[i-1]->SetFillColor(kOrange-3);
    hTM3_vv[i-1]->SetFillColor(kOrange-3);
    hTM4_vv[i-1]->SetFillColor(kOrange-3);
    hTM5_vv[i-1]->SetFillColor(kOrange-3);
    h2LL_vv[i-1]->SetFillColor(kOrange-3);
    h2TM_vv[i-1]->SetFillColor(kOrange-3);

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

    //MR
    hsLL1.push_back(new THStack("LL1_"+b+"b","LL1_"+b+"b"));
    hsLL1[i-1]->Add(hLL1_vv[i-1]);
    hsLL1[i-1]->Add(hLL1_znn[i-1]);
    hsLL1[i-1]->Add(hLL1_qcd[i-1]);  
    hsLL1[i-1]->Add(hLL1_ttwj[i-1]);
    hsLL1[i-1]->Add(hLL1_sig[i-1]);
    hsLL2.push_back(new THStack("LL2_"+b+"b","LL2_"+b+"b"));
    hsLL2[i-1]->Add(hLL2_vv[i-1]);
    hsLL2[i-1]->Add(hLL2_znn[i-1]);
    hsLL2[i-1]->Add(hLL2_qcd[i-1]);  
    hsLL2[i-1]->Add(hLL2_ttwj[i-1]);
    hsLL2[i-1]->Add(hLL2_sig[i-1]);
    hsLL3.push_back(new THStack("LL3_"+b+"b","LL3_"+b+"b"));
    hsLL3[i-1]->Add(hLL3_vv[i-1]);
    hsLL3[i-1]->Add(hLL3_znn[i-1]);
    hsLL3[i-1]->Add(hLL3_qcd[i-1]);  
    hsLL3[i-1]->Add(hLL3_ttwj[i-1]);
    hsLL3[i-1]->Add(hLL3_sig[i-1]);
    hsLL4.push_back(new THStack("LL4_"+b+"b","LL4_"+b+"b"));
    hsLL4[i-1]->Add(hLL4_vv[i-1]);
    hsLL4[i-1]->Add(hLL4_znn[i-1]);
    hsLL4[i-1]->Add(hLL4_qcd[i-1]);  
    hsLL4[i-1]->Add(hLL4_ttwj[i-1]);
    hsLL4[i-1]->Add(hLL4_sig[i-1]);
    hsLL5.push_back(new THStack("LL5_"+b+"b","LL5_"+b+"b"));
    hsLL5[i-1]->Add(hLL5_vv[i-1]);
    hsLL5[i-1]->Add(hLL5_znn[i-1]);
    hsLL5[i-1]->Add(hLL5_qcd[i-1]);  
    hsLL5[i-1]->Add(hLL5_ttwj[i-1]);
    hsLL5[i-1]->Add(hLL5_sig[i-1]);
    hsTM1.push_back(new THStack("TM1_"+b+"b","TM1_"+b+"b"));
    hsTM1[i-1]->Add(hTM1_vv[i-1]);
    hsTM1[i-1]->Add(hTM1_znn[i-1]);
    hsTM1[i-1]->Add(hTM1_qcd[i-1]);  
    hsTM1[i-1]->Add(hTM1_ttwj[i-1]);
    hsTM1[i-1]->Add(hTM1_sig[i-1]);
    hsTM2.push_back(new THStack("TM2_"+b+"b","TM2_"+b+"b"));
    hsTM2[i-1]->Add(hTM2_vv[i-1]);
    hsTM2[i-1]->Add(hTM2_znn[i-1]);
    hsTM2[i-1]->Add(hTM2_qcd[i-1]);  
    hsTM2[i-1]->Add(hTM2_ttwj[i-1]);
    hsTM2[i-1]->Add(hTM2_sig[i-1]);
    hsTM3.push_back(new THStack("TM3_"+b+"b","TM3_"+b+"b"));
    hsTM3[i-1]->Add(hTM3_vv[i-1]);
    hsTM3[i-1]->Add(hTM3_znn[i-1]);
    hsTM3[i-1]->Add(hTM3_qcd[i-1]);  
    hsTM3[i-1]->Add(hTM3_ttwj[i-1]);
    hsTM3[i-1]->Add(hTM3_sig[i-1]);
    hsTM4.push_back(new THStack("TM4_"+b+"b","TM4_"+b+"b"));
    hsTM4[i-1]->Add(hTM4_vv[i-1]);
    hsTM4[i-1]->Add(hTM4_znn[i-1]);
    hsTM4[i-1]->Add(hTM4_qcd[i-1]);  
    hsTM4[i-1]->Add(hTM4_ttwj[i-1]);
    hsTM4[i-1]->Add(hTM4_sig[i-1]);
    hsTM5.push_back(new THStack("TM5_"+b+"b","TM5_"+b+"b"));
    hsTM5[i-1]->Add(hTM5_vv[i-1]);
    hsTM5[i-1]->Add(hTM5_znn[i-1]);
    hsTM5[i-1]->Add(hTM5_qcd[i-1]);  
    hsTM5[i-1]->Add(hTM5_ttwj[i-1]);
    hsTM5[i-1]->Add(hTM5_sig[i-1]);
    hs2LL.push_back(new THStack("2LL_"+b+"b","2LL_"+b+"b"));
    hs2LL[i-1]->Add(h2LL_vv[i-1]);
    hs2LL[i-1]->Add(h2LL_znn[i-1]);
    hs2LL[i-1]->Add(h2LL_qcd[i-1]);  
    hs2LL[i-1]->Add(h2LL_ttwj[i-1]);
    hs2LL[i-1]->Add(h2LL_sig[i-1]);
    hs2TM.push_back(new THStack("2TM_"+b+"b","2TM_"+b+"b"));
    hs2TM[i-1]->Add(h2TM_vv[i-1]);
    hs2TM[i-1]->Add(h2TM_znn[i-1]);
    hs2TM[i-1]->Add(h2TM_qcd[i-1]);  
    hs2TM[i-1]->Add(h2TM_ttwj[i-1]);
    hs2TM[i-1]->Add(h2TM_sig[i-1]);
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

    //MR - data
    if(!ABCD)
      {
	hLL1_data[bBin-1]->SetBinContent(binIndex,(ws->function( "oneLooseLep_"+binName+"_Theta1_Count" ))->getVal());
	hLL2_data[bBin-1]->SetBinContent(binIndex,(ws->function( "oneLooseLep_"+binName+"_Theta2_Count" ))->getVal());
	hLL3_data[bBin-1]->SetBinContent(binIndex,(ws->function( "oneLooseLep_"+binName+"_Theta3_Count" ))->getVal());
	hLL4_data[bBin-1]->SetBinContent(binIndex,(ws->function( "oneLooseLep_"+binName+"_Theta4_Count" ))->getVal());
	hLL5_data[bBin-1]->SetBinContent(binIndex,(ws->function( "oneLooseLep_"+binName+"_Theta5_Count" ))->getVal());
	hTM1_data[bBin-1]->SetBinContent(binIndex,(ws->function( "oneTightMu_"+binName+"_Theta1_Count" ))->getVal());
	hTM2_data[bBin-1]->SetBinContent(binIndex,(ws->function( "oneTightMu_"+binName+"_Theta2_Count" ))->getVal());
	hTM3_data[bBin-1]->SetBinContent(binIndex,(ws->function( "oneTightMu_"+binName+"_Theta3_Count" ))->getVal());
	hTM4_data[bBin-1]->SetBinContent(binIndex,(ws->function( "oneTightMu_"+binName+"_Theta4_Count" ))->getVal());
	hTM5_data[bBin-1]->SetBinContent(binIndex,(ws->function( "oneTightMu_"+binName+"_Theta5_Count" ))->getVal());
	h2LL_data[bBin-1]->SetBinContent(binIndex,(ws->function( "twoLooseLep_"+binName+"_Count" ))->getVal());
	h2TM_data[bBin-1]->SetBinContent(binIndex,(ws->function( "twoTightMu_"+binName+"_Count" ))->getVal());
      }

    hZL_sig[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLepton_"+binName+"_SignalDataYield" ))->getVal() );
    if(ABCD) hSL_sig[bBin-1]->SetBinContent( binIndex, (ws->function( "oneLepton_"+binName+"_SignalDataYield" ))->getVal() );
    hLDP_sig[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_SignalDataYield" ))->getVal() );
    cout << (ws->function( "zeroLepton_"+binName+"_SignalDataYield" ))->getVal() << " " << (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_SignalDataYield" ))->getVal() << endl;

    hZL_ttwj[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLepton_"+binName+"_TopWJetsDataYield" ))->getVal() );
    if(ABCD) hSL_ttwj[bBin-1]->SetBinContent( binIndex, (ws->function( "oneLepton_"+binName+"_TopWJetsDataYield" ))->getVal() );
    hLDP_ttwj[bBin-1]->SetBinContent( binIndex, (ws->function( "zeroLeptonLowDeltaPhiN_"+binName+"_TopWJetsDataYield" ))->getVal() );
    
    //MR - ttwj
    if(!ABCD)
      {
	hLL1_ttwj[bBin-1]->SetBinContent(binIndex,(ws->function( "oneLooseLep_"+binName+"_Theta1_TopWJetsDataYield" ))->getVal());
	hLL2_ttwj[bBin-1]->SetBinContent(binIndex,(ws->function( "oneLooseLep_"+binName+"_Theta2_TopWJetsDataYield" ))->getVal());
	hLL3_ttwj[bBin-1]->SetBinContent(binIndex,(ws->function( "oneLooseLep_"+binName+"_Theta3_TopWJetsDataYield" ))->getVal());
	hLL4_ttwj[bBin-1]->SetBinContent(binIndex,(ws->function( "oneLooseLep_"+binName+"_Theta4_TopWJetsDataYield" ))->getVal());
	hLL5_ttwj[bBin-1]->SetBinContent(binIndex,(ws->function( "oneLooseLep_"+binName+"_Theta5_TopWJetsDataYield" ))->getVal());
	hTM1_ttwj[bBin-1]->SetBinContent(binIndex,(ws->function( "oneTightMu_"+binName+"_Theta1_TopWJetsDataYield" ))->getVal());
	hTM2_ttwj[bBin-1]->SetBinContent(binIndex,(ws->function( "oneTightMu_"+binName+"_Theta2_TopWJetsDataYield" ))->getVal());
	hTM3_ttwj[bBin-1]->SetBinContent(binIndex,(ws->function( "oneTightMu_"+binName+"_Theta3_TopWJetsDataYield" ))->getVal());
	hTM4_ttwj[bBin-1]->SetBinContent(binIndex,(ws->function( "oneTightMu_"+binName+"_Theta4_TopWJetsDataYield" ))->getVal());
	hTM5_ttwj[bBin-1]->SetBinContent(binIndex,(ws->function( "oneTightMu_"+binName+"_Theta5_TopWJetsDataYield" ))->getVal());
	h2LL_ttwj[bBin-1]->SetBinContent(binIndex,(ws->function( "twoLooseLep_"+binName+"_TopWJetsDataYield" ))->getVal());
	h2TM_ttwj[bBin-1]->SetBinContent(binIndex,(ws->function( "twoTightMu_"+binName+"_TopWJetsDataYield" ))->getVal());
      }

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
   
  //Dielectron and dimuon
  for(int i=1; i<=maxNM; i++) {//hardcoded
    for(int j=1; j<=maxNH; j++) {
      TString binLabel = "M"; binLabel+=i; binLabel+="_H"; binLabel+=j;
      int binIndex = 1 + (maxNH+1)*(i-1) + j;
      
      if(skipMHBin(i, j)) continue;
      
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

    //MR
    if(!ABCD)
      {
	hsLL1[i-1]->SetMaximum(maxS*hLL1_data[i-1]->GetMaximum());
	hsLL2[i-1]->SetMaximum(maxS*hLL2_data[i-1]->GetMaximum());
	hsLL3[i-1]->SetMaximum(maxS*hLL3_data[i-1]->GetMaximum());
	hsLL4[i-1]->SetMaximum(maxS*hLL4_data[i-1]->GetMaximum());
	hsLL5[i-1]->SetMaximum(maxS*hLL5_data[i-1]->GetMaximum());
	hsTM1[i-1]->SetMaximum(maxS*hTM1_data[i-1]->GetMaximum());
	hsTM2[i-1]->SetMaximum(maxS*hTM2_data[i-1]->GetMaximum());
	hsTM3[i-1]->SetMaximum(maxS*hTM3_data[i-1]->GetMaximum());
	hsTM4[i-1]->SetMaximum(maxS*hTM4_data[i-1]->GetMaximum());
	hsTM5[i-1]->SetMaximum(maxS*hTM5_data[i-1]->GetMaximum());
	hs2LL[i-1]->SetMaximum(maxS*h2LL_data[i-1]->GetMaximum());
	hs2TM[i-1]->SetMaximum(maxS*h2TM_data[i-1]->GetMaximum());
      }
  }

  if(logy) {
    hZee_z->SetMinimum(.1);
    hZmm_z->SetMinimum(.1);

    for(int i = 1; i<=maxNB; i++) {
      hsZL[i-1]->SetMinimum(.1);
      hsSL[i-1]->SetMinimum(.1);
      hsLDP[i-1]->SetMinimum(.1);

      //MR
      if(!ABCD)
	{
	  hsLL1[i-1]->SetMinimum(.1);
	  hsLL2[i-1]->SetMinimum(.1);
	  hsLL3[i-1]->SetMinimum(.1);
	  hsLL4[i-1]->SetMinimum(.1);
	  hsLL5[i-1]->SetMinimum(.1);
	  hsTM1[i-1]->SetMinimum(.1);
	  hsTM2[i-1]->SetMinimum(.1);
	  hsTM3[i-1]->SetMinimum(.1);
	  hsTM4[i-1]->SetMinimum(.1);
	  hsTM5[i-1]->SetMinimum(.1);
	  hs2LL[i-1]->SetMinimum(.1);
	  hs2TM[i-1]->SetMinimum(.1);
	}
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


  TCanvas * cMRLL = new TCanvas("cMRLL", "cMRLL", 1000, 1000);
  cMRLL->Divide(3,5);
  cMRLL->cd(1);
  if(hsLL1[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL1[0]->Draw();				
  hLL1_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRLL->cd(2);
  if(hsLL1[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL1[1]->Draw();				
  hLL1_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRLL->cd(3);
  if(hsLL1[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL1[2]->Draw();				
  hLL1_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  //
  cMRLL->cd(4);
  if(hsLL2[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL2[0]->Draw();				
  hLL2_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRLL->cd(5);
  if(hsLL2[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL2[1]->Draw();				
  hLL2_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRLL->cd(6);
  if(hsLL2[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL2[2]->Draw();				
  hLL2_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  //
  cMRLL->cd(7);
  if(hsLL3[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL3[0]->Draw();				
  hLL3_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRLL->cd(8);
  if(hsLL3[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL3[1]->Draw();				
  hLL3_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRLL->cd(9);
  if(hsLL3[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL3[2]->Draw();				
  hLL3_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  //
  cMRLL->cd(10);
  if(hsLL4[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL4[0]->Draw();				
  hLL4_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRLL->cd(11);
  if(hsLL4[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL4[1]->Draw();				
  hLL4_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRLL->cd(12);
  if(hsLL4[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL4[2]->Draw();				
  hLL4_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  //
  cMRLL->cd(13);
  if(hsLL5[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL5[0]->Draw();				
  hLL5_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRLL->cd(14);
  if(hsLL5[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL5[1]->Draw();				
  hLL5_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRLL->cd(15);
  if(hsLL5[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsLL5[2]->Draw();				
  hLL5_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);


  TCanvas * cMRTM = new TCanvas("cMRTM", "cMRTM", 1000, 1000);
  cMRTM->Divide(3,5);
  cMRTM->cd(1);
  if(hsTM1[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM1[0]->Draw();				
  hTM1_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRTM->cd(2);
  if(hsTM1[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM1[1]->Draw();				
  hTM1_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRTM->cd(3);
  if(hsTM1[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM1[2]->Draw();				
  hTM1_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  //
  cMRTM->cd(4);
  if(hsTM2[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM2[0]->Draw();				
  hTM2_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRTM->cd(5);
  if(hsTM2[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM2[1]->Draw();				
  hTM2_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRTM->cd(6);
  if(hsTM2[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM2[2]->Draw();				
  hTM2_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  //
  cMRTM->cd(7);
  if(hsTM3[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM3[0]->Draw();				
  hTM3_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRTM->cd(8);
  if(hsTM3[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM3[1]->Draw();				
  hTM3_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRTM->cd(9);
  if(hsTM3[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM3[2]->Draw();				
  hTM3_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  //
  cMRTM->cd(10);
  if(hsTM4[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM4[0]->Draw();				
  hTM4_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRTM->cd(11);
  if(hsTM4[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM4[1]->Draw();				
  hTM4_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRTM->cd(12);
  if(hsTM4[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM4[2]->Draw();				
  hTM4_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  //
  cMRTM->cd(13);
  if(hsTM5[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM5[0]->Draw();				
  hTM5_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRTM->cd(14);
  if(hsTM5[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM5[1]->Draw();				
  hTM5_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMRTM->cd(15);
  if(hsTM5[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hsTM5[2]->Draw();				
  hTM5_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);


  TCanvas * cMR2 = new TCanvas("cMR2", "cMR2", 1000, 400);
  cMR2->Divide(3,2);
  cMR2->cd(1);
  if(hs2LL[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hs2LL[0]->Draw();				
  h2LL_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMR2->cd(2);
  if(hs2LL[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hs2LL[1]->Draw();				
  h2LL_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMR2->cd(3);
  if(hs2LL[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hs2LL[2]->Draw();				
  h2LL_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  //
  cMR2->cd(4);
  if(hs2TM[0]->GetMaximum()>0) gPad->SetLogy(logy);
  hs2TM[0]->Draw();				
  h2TM_data[0]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMR2->cd(5);
  if(hs2TM[1]->GetMaximum()>0) gPad->SetLogy(logy);
  hs2TM[1]->Draw();				
  h2TM_data[1]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  cMR2->cd(6);
  if(hs2TM[2]->GetMaximum()>0) gPad->SetLogy(logy);
  hs2TM[2]->Draw();				
  h2TM_data[2]->Draw("HISTPE SAME");
  gPad->SetBottomMargin(0.15);
  



  //label axes
  for(map<TString,s3d>::iterator it = binMap.begin(); it != binMap.end(); it++) {
    s3d binInfo = it->second;
    TString binLabel = "M"; binLabel+=binInfo.nM; binLabel+="_H"; binLabel+=binInfo.nH; binLabel+="_"; binLabel+=binInfo.nB; binLabel+="b";
    int binIndex = 1 + (maxNH+1)*(binInfo.nM-1) + binInfo.nH;
    int bBin = binInfo.nB;

    c->cd();
    hsZL[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
    hsSL[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
    hsLDP[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );

    //MR
    if(!ABCD)
      {
	cMRLL->cd();
	hsLL1[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
	hsLL2[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
	hsLL3[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
	hsLL4[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
	hsLL5[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
	cMRTM->cd();
	hsTM1[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
	hsTM2[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
	hsTM3[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
	hsTM4[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
	hsTM5[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
	cMR2->cd();
	hs2LL[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
	hs2TM[bBin-1]->GetXaxis()->SetBinLabel( binIndex, binLabel );
      }
  }

  for(int i=1; i<=maxNB; i++) {
    hsZL[i-1]->GetXaxis()->LabelsOption("v");
    hsSL[i-1]->GetXaxis()->LabelsOption("v");
    hsLDP[i-1]->GetXaxis()->LabelsOption("v");

    //MR
    if(!ABCD)
      {
	hsLL1[i-1]->GetXaxis()->LabelsOption("v");
	hsLL2[i-1]->GetXaxis()->LabelsOption("v");
	hsLL3[i-1]->GetXaxis()->LabelsOption("v");
	hsLL4[i-1]->GetXaxis()->LabelsOption("v");
	hsLL5[i-1]->GetXaxis()->LabelsOption("v");
	hsTM1[i-1]->GetXaxis()->LabelsOption("v");
	hsTM2[i-1]->GetXaxis()->LabelsOption("v");
	hsTM3[i-1]->GetXaxis()->LabelsOption("v");
	hsTM4[i-1]->GetXaxis()->LabelsOption("v");
	hsTM5[i-1]->GetXaxis()->LabelsOption("v");
	hs2LL[i-1]->GetXaxis()->LabelsOption("v");
	hs2TM[i-1]->GetXaxis()->LabelsOption("v");
      }
  }
  hZee_z->GetXaxis()->LabelsOption("v");
  hZmm_z->GetXaxis()->LabelsOption("v");



  gPad->Modified();

  TString logstring="";
  if(logy) logstring="_log";
  c->SaveAs("stack"+logstring+".pdf");
  cMRLL->SaveAs("stack_MRLL"+logstring+".pdf");
  cMRTM->SaveAs("stack_MRTM"+logstring+".pdf");
  cMR2->SaveAs("stack_MR2"+logstring+".pdf");
  
  
  
  wstf->Close();

}


void analyzeFit(TString workspaceFile, TString binFilesFile, TString datFile) {
  cout << "Starting analyzeFit" << endl;

  
  extractFromWorkspace(workspaceFile, datFile, true);

  integratedTotals(workspaceFile, binFilesFile, datFile);
  //integratedTotals(workspaceFile);
  //owenPlots();

}
