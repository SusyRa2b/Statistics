#include <iostream>

#include "TFile.h"
#include "TString.h"

#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooFitResult.h"

#include "RooStats/ModelConfig.h"


using namespace RooFit ;
using namespace RooStats ;
using namespace std;


void minimalFit(TString workspaceFile = "test.root", double signalCrossSectionGuess = 50.0, double signalCrossSectionLow = 0.0, double signalCrossSectionHigh = 1000.0, bool updateWS = false) 
{
  
  TFile* wstf = 0;  
  if( updateWS == true ) wstf = new TFile ( workspaceFile, "UPDATE", "", 1);
  else wstf = new TFile ( workspaceFile, "READ", "", 1);
    
  
  RooWorkspace* ws = (RooWorkspace*)wstf->Get("workspace");

  RooDataSet* rds = (RooDataSet*) ws->obj( "data" );
  rds->Print("v");

  RooAbsPdf* likelihood = ws->pdf( "likelihood" );

  RooRealVar* signalCrossSection = ws->var( "signalCrossSection" );
  signalCrossSection->setVal( signalCrossSectionGuess );
  signalCrossSection->setRange( signalCrossSectionLow, signalCrossSectionHigh );

  RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(0) );
  fitResult->Print();
  cout << "RooFitResult status = " << fitResult->status() << endl;
  cout << "RooFitResult minNll = " << fitResult->minNll() << endl;

  if(updateWS) ws->Write();
  wstf->Close();

  return;
}
