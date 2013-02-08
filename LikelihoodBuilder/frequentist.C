#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooProfileLL.h"
#include "RooFitResult.h"

#include "RooStats/HypoTestResult.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/ProfileLikelihoodTestStat.h"


#include <iostream>
#include <fstream>
#include <time.h>

using namespace RooFit;
using namespace RooStats;
using namespace std;


//void RunToyScan5(TString fileName, double startVal, double stopVal, TString outFile) {
void frequentist(TString fileName) {
  cout << "Starting frequentist " << time(NULL) << endl;
  double startVal = 0;
  double stopVal = 200;
  TString outFile = "";

  int nToys = 1 ;
  int nscanpoints = 2 ;

  /*
  gROOT->LoadMacro("RooBetaPdf.cxx+") ;
  gROOT->LoadMacro("RooRatio.cxx+") ;
  gROOT->LoadMacro("RooPosDefCorrGauss.cxx+") ;
  */

  // get relevant objects out of the "ws" file

  TFile *file = TFile::Open(fileName);
  if(!file){
    cout <<"file not found" << endl;
    return;
  } 

  RooWorkspace* w = (RooWorkspace*) file->Get("workspace");
  if(!w){
    cout <<"workspace not found" << endl;
    return;
  }

  ModelConfig* mc = (ModelConfig*) w->obj("S+B_model");
  RooAbsData* data = w->data("data");

  if( !data || !mc ){
    w->Print();
    cout << "data or ModelConfig was not found" <<endl;
    return;
  }

  RooRealVar* myPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();
  myPOI->setRange(0, 1000.);

  ModelConfig* bModel = (ModelConfig*) w->obj("B_model");
  ModelConfig* sbModel = (ModelConfig*) w->obj("S+B_model");

  ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
  profll.SetPrintLevel(2);
  profll.SetOneSided(1);
  TestStatistic * testStat = &profll;

  HypoTestCalculatorGeneric *  hc = 0;
  hc = new FrequentistCalculator(*data, *bModel, *sbModel);
  
  ToyMCSampler *toymcs = (ToyMCSampler*)hc->GetTestStatSampler();
  toymcs->SetMaxToys(10000);
  toymcs->SetNEventsPerToy(1);
  toymcs->SetTestStatistic(testStat);


  ((FrequentistCalculator *)hc)->SetToys(nToys,nToys);
  
  HypoTestInverter calc(*hc);
  calc.SetConfidenceLevel(0.95);
  calc.UseCLs(true);
  //calc.SetVerbose(true);
  calc.SetVerbose(2);

  cout << "About to set fixed scan " << time(NULL) << endl;
  calc.SetFixedScan(nscanpoints,startVal,stopVal);
  cout << "About to do inverter " << time(NULL) << endl;
  HypoTestInverterResult * res_toysCLs_calculator = calc.GetInterval();

  cout << "CLs = " << res_toysCLs_calculator->UpperLimit() 
	    << "   CLs_exp = " << res_toysCLs_calculator->GetExpectedUpperLimit(0) 
	    << "   CLs_exp(-1s) = " << res_toysCLs_calculator->GetExpectedUpperLimit(-1) 
	    << "   CLs_exp(+1s) = " << res_toysCLs_calculator->GetExpectedUpperLimit(1) << endl ;

  /*
  // dump results string to output file
  ofstream outStream ;
  outStream.open(outFile,ios::app) ;
  
  outStream << "CLs = " << res_toysCLs_calculator->UpperLimit() 
	    << "   CLs_exp = " << res_toysCLs_calculator->GetExpectedUpperLimit(0) 
	    << "   CLs_exp(-1s) = " << res_toysCLs_calculator->GetExpectedUpperLimit(-1) 
	    << "   CLs_exp(+1s) = " << res_toysCLs_calculator->GetExpectedUpperLimit(1) << endl ;
  
  outStream.close() ;
  */


  cout << "End of frequentist " << time(NULL) << endl;
  return ;

}
