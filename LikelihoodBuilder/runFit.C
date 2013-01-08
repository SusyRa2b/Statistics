#include <iostream>
#include <fstream>

#include "TString.h"

#include "likelihoodBuilder.C"
#include "minimalFit.C"
#include "minimalProfileLikelihood.C"
#include "analyzeFit.C"

using namespace std;

void runFit(TString inpath, TString outpath, TString outname = "")
{

  TString datFile = outpath; datFile+="dat_"; datFile+=outname; datFile+=".dat";
  if(datFile != "") {
    ofstream myfile;
    myfile.open(datFile.Data(), ios::out | ios::trunc);
    assert(myfile.is_open());
    myfile.close();
  }

  likelihoodBuilder(inpath+"setupFile.dat", inpath+"binFilesFile.dat", inpath, inpath+"sig1/", "workspace", outpath+"likelihood_"+outname+".root", inpath+"binFilesFileMR.dat", inpath+"countsMR/");
   
  //int status=0;
  int status = minimalFit(outpath+"likelihood_"+outname+".root", 5, 0, 1000, true, outpath+"dat_"+outname+".dat"); //susy floating
  //int status = minimalFit(outpath+"likelihood_"+outname+".root", 0, 0, 1000, true, outpath+"dat_"+outname+".dat", true); //fix susy to zero
  //int status = minimalFit(outpath+"likelihood_"+outname+".root", 56.628, 0, 1000, true, outpath+"dat_"+outname+".dat", true); //fix susy to 56.628
  if(status != 0) return;


  return;

  profileLikelihoodLimit(outpath+"likelihood_"+outname+".root", "workspace", "S+B_model", "B_model", "data", "modelName", 0.682, 0.0, 1000.0, false, outpath+"dat_"+outname+".dat");
  profileLikelihoodLimit(outpath+"likelihood_"+outname+".root", "workspace", "S+B_model", "B_model", "data", "modelName", 0.95, 0.0, 1000.0, false, outpath+"dat_"+outname+".dat");

  analyzeFit(outpath+"likelihood_"+outname+".root", inpath+"binFilesFile.dat", outpath+"dat_"+outname+".dat"); 

  
  if(datFile != "") {
    ofstream myfile;
    myfile.open(datFile.Data(), ios::out | ios::app);
    assert(myfile.is_open());
    myfile << endl;
    myfile.close();
  }

}

