#include "setup.C"
#include "setupAnalysis.C"

void callRunFit(TString inpath, TString outpath, TString outname = ""){

  setup();
  setupAnalysis();

  runFit(inpath, outpath, outname);

}
