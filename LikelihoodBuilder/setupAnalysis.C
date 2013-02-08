#include "TSystem.h"

void setupAnalysis(){
  gSystem->CompileMacro("AsymptoticCalculatorNestedSimPdf.cxx", "kO") ;
  gSystem->CompileMacro("minimalFit.C"                        , "kO") ;
  gSystem->CompileMacro("minimalProfileLikelihood.C"          , "kO") ;
  gSystem->CompileMacro("singleAsymptotic.C"                  , "kO") ;
  gSystem->CompileMacro("frequentist.C"                       , "kO") ;
  gSystem->CompileMacro("analyzeFit.C"                        , "kO") ;
  gSystem->CompileMacro("runFit.C"                            , "kO") ;
  return;
}
