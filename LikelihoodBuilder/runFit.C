#include <iostream>
#include <fstream>

#include "TString.h"
#include "RooMsgService.h"

#include "likelihoodBuilder.C"
#include "minimalFit.C"
#include "minimalProfileLikelihood.C"
#include "singleAsymptotic.C"
#include "frequentist.C"
#include "analyzeFit.C"

using namespace std;

void runFit(TString inpath, TString outpath, TString outname = "", TString option = "")
{

  //RooFit::PrintLevel(0); 
  //ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(3) ;
  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("OldMinuit","BenKreis"); //put trash in for algorithm to make sure it doesn't matter
  //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Migrad"); 
  //RooMsgService::instance().addStream(INFO,OutputStream(cout));


  TString datFile = outpath; datFile+="dat_"; datFile+=outname; datFile+=".dat";
  if(datFile != "") {
    ofstream myfile;
    myfile.open(datFile.Data(), ios::out | ios::trunc);
    assert(myfile.is_open());
    myfile.close();
  }
  
  likelihoodBuilder(inpath+"setupFile.dat", inpath+"binFilesFile.dat", inpath, inpath+"sig1/", "workspace", outpath+"likelihood_"+outname+".root", inpath+"binFilesFileMR.dat", inpath+"countsMR/", option);

  return;

  int status=0;
  //status = minimalFit(outpath+"likelihood_"+outname+".root", 0, 0, 1000, true, outpath+"dat_"+outname+".dat", true); //fix susy to zero
  //if(status != 0) return;
  //status = minimalFit(outpath+"likelihood_"+outname+".root", 5, 0, 1000, true, outpath+"dat_"+outname+".dat"); //susy floating 
  //if(status != 0) return;
  //int status = minimalFit(outpath+"likelihood_"+outname+".root", 56.628, 0, 1000, true, outpath+"dat_"+outname+".dat", true); //fix susy to 56.628


  //singleAsymptotic(outpath+"likelihood_"+outname+".root");
  frequentist(outpath+"likelihood_"+outname+".root");
  return;


  //profileLikelihoodLimit(outpath+"likelihood_"+outname+".root", "workspace", "S+B_model", "B_model", "data", "modelName", 0.682, 0.0, 1000.0, false, outpath+"dat_"+outname+".dat");
  profileLikelihoodLimit(outpath+"likelihood_"+outname+".root", "workspace", "S+B_model", "B_model", "data", "modelName", 0.95, 0.0, 1000.0, false, outpath+"dat_"+outname+".dat");
  return;
  
  analyzeFit(outpath+"likelihood_"+outname+".root", inpath+"binFilesFile.dat", outpath+"dat_"+outname+".dat"); 
  
  if(datFile != "") {
    ofstream myfile;
    myfile.open(datFile.Data(), ios::out | ios::app);
    assert(myfile.is_open());
    myfile << endl;
    myfile.close();
  }

}

