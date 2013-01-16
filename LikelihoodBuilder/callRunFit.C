#include "setup.C"
#include "setupAnalysis.C"

void callRunFit(TString inpath, TString outpath, TString outname = "", TString option="allWidths"){

  setup();
  setupAnalysis();

  runFit(inpath, outpath, outname, option);

}
