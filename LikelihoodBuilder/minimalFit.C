#include <iostream>
#include <string.h>
#include <complex>
#include <map>
#include <cassert>

#include "TCanvas.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TText.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TH2F.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TStringLong.h"
#include "TString.h"
#include "TPRegexp.h"

#include "RooArgSet.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooTrace.h"
#include "RooUniform.h"
#include "RooAddition.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooRatio.h"
#include "RooAddition.h"
#include "RooPoisson.h"
#include "RooFitResult.h"

#include "RooStats/ModelConfig.h"

#include "TMath.h"

using namespace RooFit ;
using namespace RooStats ;
using namespace std;

void minimalFit() 
{

  TFile* wstf = new TFile ( "test.root" );

  RooWorkspace* ws = (RooWorkspace*)wstf->Get("workspace");

  //ModelConfig* modelConfig = (ModelConfig*) ws->obj("S+B_model");

  RooDataSet* rds = (RooDataSet*) ws->obj( "data" );

  rds->Print("v");

  //RooAbsPdf* likelihood = modelConfig->GetPdf() ;
  RooAbsPdf* likelihood = ws->pdf( "likelihood" );

  RooRealVar* signalCrossSection = ws->var( "signalCrossSection" );

  //RooAbsReal* pll = likelihood->createProfile(*signalCrossSection);

  //pll->getVal();

  signalCrossSection->setConstant();
  
  RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(0) );

  fitResult->Print();

  delete fitResult;

  signalCrossSection->setConstant(kFALSE);

  fitResult = likelihood->fitTo( *rds, Save(true) , PrintLevel(0));
  
  cout << "RooFitResult status = " << fitResult->status() << endl;
  cout << "RooFitResult minNll = " << fitResult->minNll() << endl;

  return;
}
