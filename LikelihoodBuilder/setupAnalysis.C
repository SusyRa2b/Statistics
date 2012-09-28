#include "TSystem.h"

void setupAnalysis(){
  gSystem->CompileMacro("minimalFit.C"                      ,"kO") ;
  gSystem->CompileMacro("minimalProfileLikelihood.C"        ,"kO") ;
  gSystem->CompileMacro("analyzeFit.C"                      ,"kO") ;
  gSystem->CompileMacro("runFit.C"                         ,"kO") ;
  return;
}
