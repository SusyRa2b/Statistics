#include <iostream>
#include <fstream>

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


int minimalFit(TString workspaceFile = "test.root", double signalCrossSectionGuess = 50.0, double signalCrossSectionLow = 0.0, double signalCrossSectionHigh = 1000.0, bool updateWS = false, TString datFile = "") 
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

  double sig = signalCrossSection->getVal();
  double sigerr = signalCrossSection->getError();
  cout << "signalCrossSection = " << sig << " +- " <<  sigerr << endl;
  cout << "RooFitResult status = " << fitResult->status() << endl;
  cout << "RooFitResult minNll = " << fitResult->minNll() << endl;

  //if(fitResult->status() != 0) return fitResult->status();
  if(fabs((sig-signalCrossSectionGuess)/signalCrossSectionGuess)<1e-5) return 1;

  if(datFile != "") {
    ofstream myfile;
    myfile.open(datFile.Data(), ios::out | ios::app);
    assert(myfile.is_open());
    myfile << sig << " " << sigerr << " ";
    myfile.close();
  }

  if(updateWS) {
    //ws->import(fitResult->GetName());//if i do this, how do i get it out of the workspace later??
    fitResult->Write();
  }

  if(updateWS) ws->Write();
  wstf->Close();

  //return fitResult->status();
  return 0;
}
