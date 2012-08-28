#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"
#include "RooProfileLL.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooCurve.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "THStack.h"
#include "TPRegexp.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/ToyMCSampler.h"

#include <iostream>
#include <map>
#include <string>
#include <fstream>

using namespace RooFit;
using namespace RooStats;
using namespace std;


void profileLikelihoodLimit(const char * fileName = "test.root",
                            const char * wsName = "workspace",
                            const char * modelSBName = "S+B_model",
                            const char * modelBName = "B_model",
                            const char * dataName = "data",
                            const char * modelName = "",
			    double percentInterval = 0.95, 
			    double signalCrossSectionLow = 0.0,
			    double signalCrossSectionHigh = 500.0,
                            bool doPlots = false)

{
  cout <<"Starting Profile Likelihood Calculator"<<endl;

  TFile *file = TFile::Open(fileName);
  // if input file was specified byt not found, quit
  if(!file){
    cout <<"file not found" << endl;
    return;
  }
  
  // get the workspace out of the file
  RooWorkspace* w = (RooWorkspace*) file->Get(wsName);
  if(!w){
    cout <<"workspace not found" << endl;
    return;
  }
  
  // get the modelConfig out of the file
  ModelConfig* mc = (ModelConfig*) w->obj(modelSBName);
  
  // get the modelConfig out of the file
  RooAbsData* data = w->data(dataName);
  RooAbsPdf* pdf = mc->GetPdf();
  
  // make sure ingredients are found
  if(!data || !mc){
    w->Print();
    cout << "data or ModelConfig was not found" <<endl;
    return;
  }
  

  RooRealVar* firstPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();
  TString firstPOIName = firstPOI->GetName();
  firstPOI->setRange(signalCrossSectionLow, signalCrossSectionHigh);
  
  ProfileLikelihoodCalculator pl(*data,*pdf, RooArgSet(*firstPOI));
  pl.SetTestSize( 1.0 - percentInterval );
  //pl.SetTestSize( 0.05 ) ; // 95%
  //pl.SetTestSize( 1. - .682 ) ; //68% 
  //pl.SetConfidenceLevel(0.975); // 95% one sided limit
  
  LikelihoodInterval* interval = pl.GetInterval();
    
  cout << endl;
  cout << percentInterval << " interval on " <<firstPOI->GetName()<<" is : ["<<
    interval->LowerLimit(*firstPOI) << ", "<<
    interval->UpperLimit(*firstPOI) <<"] "<<endl;
  return;

  
  if(doPlots){
    cout << "Making a plot of the profile likelihood function ....(if it is taking a lot of time use less points or the TF1 drawing option)\n";

    double asymptoticUpperLimit = interval->UpperLimit(*firstPOI);
    double asymptoticLowerLimit = interval->LowerLimit(*firstPOI);   
    firstPOI->setRange(0.8*asymptoticLowerLimit, 1.2*asymptoticUpperLimit);

    TCanvas* c = new TCanvas("c");;    
    LikelihoodIntervalPlot plot(interval);
    plot.SetNPoints(100);  // do not use too many points, it could become very slow for some models
    plot.SetMaximum(5);
    plot.SetRange(0.8*asymptoticLowerLimit,1.2*asymptoticUpperLimit);
    plot.Draw("");  // use option TF1 if too slow (plot.Draw("tf1")
    c->SaveAs(TString(modelName)+"_profileLikelihoodScan.pdf");
    delete c;
  }
  
  
  return;
}



void minimalProfileLikelihood() 
{

  profileLikelihoodLimit("test.root", "workspace", "S+B_model", "B_model", "data", "modelName", 0.95, 0.0, 500.0, true);
  
  return;
}
