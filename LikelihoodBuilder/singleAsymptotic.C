#include "TFile.h"
#include "TString.h"
#include "RooWorkspace.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"
//#include "RooStats/AsymptoticCalculator.h"
#include "AsymptoticCalculatorNestedSimPdf.h"
#include "RooStats/ProfileLikelihoodTestStat.h"

#include <iostream>
#include <time.h>

using namespace RooFit;
using namespace RooStats;
using namespace std;

void singleAsymptotic(TString fileName = "test.root")
{

  TString dataName = "data";
  //TString fileName = "test.root";
  TString wsName = "workspace";
  TString modelBName = "B_model";
  TString modelSBName = "S+B_model";

  cout << "Starting Single Asymptotic Calculator " << time(NULL) << endl;

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
  
  //Get the B and S+B model configs 
  ModelConfig* bModel = (ModelConfig*) w->obj(modelBName);
  ModelConfig* sbModel = (ModelConfig*) w->obj(modelSBName);

  // get the data and pdf
  RooAbsData* data = w->data(dataName);
  RooAbsPdf* pdf = sbModel->GetPdf();

  // make sure ingredients are found
  if(!bModel || !sbModel || !pdf || !data){
    w->Print();
    cout << "ModelConfig, pdf, or data was not found" <<endl;
    return;
  }

  
  cout << "Got what you needed.  Now set up asymptotic CLs" << endl;

  ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
  profll.SetOneSided(1);

  //AsymptoticCalculator asymptotic_calculator(*data, *bModel, *sbModel);
  AsymptoticCalculatorNestedSimPdf asymptotic_calculator(*data, *bModel, *sbModel);
  HypoTestInverter calc_asymptotic_calculator(asymptotic_calculator);
  
  calc_asymptotic_calculator.SetConfidenceLevel(0.95);
  calc_asymptotic_calculator.SetVerbose(true);
  calc_asymptotic_calculator.UseCLs(true);

  int nScanPoints = 1; double poimin = 75; double poimax = 75;
  calc_asymptotic_calculator.SetFixedScan( nScanPoints , poimin , poimax) ;


  cout << "Asymptotic CLs set up.  Now get result." << endl;

  HypoTestInverterResult* res_asymptotic_calculator ;
  res_asymptotic_calculator = calc_asymptotic_calculator.GetInterval();

  cout << "Done with the asymptotic interval calculation " << time(NULL) << endl;
  
  for(int j = 0 ; j < res_asymptotic_calculator->ArraySize() ; j++)
    {
      cout << "CLb(" << j << ")          : " <<res_asymptotic_calculator->CLb(j)           << endl;
      cout << "CLbError(" << j << ")     : " <<res_asymptotic_calculator->CLbError(j)      << endl;
      cout << "CLs(" << j << ")          : " <<res_asymptotic_calculator->CLs(j)            << endl;
      cout << "CLsError(" << j << ")     : " <<res_asymptotic_calculator->CLsError(j)      << endl;
      cout << "CLsplusb(" << j << ")     : " <<res_asymptotic_calculator->CLsplusb(j)      << endl;
      cout << "CLsplusbError(" << j << "): " <<res_asymptotic_calculator->CLsplusbError(j) << endl;
    }

  
  


  
}
