#include "TSystem.h"

void setup(){
  gSystem->CompileMacro("RooProdPdfLogSum.cxx"                        ,"kO") ;
  gSystem->CompileMacro("RooPoissonLogEval.cxx"                       ,"kO") ;
  gSystem->CompileMacro("RooRatio.cxx"                                ,"kO") ;
  gSystem->CompileMacro("RooBetaPdf.cxx"                              ,"kO") ;
  gSystem->CompileMacro("RooBetaPrimePdf.cxx"                         ,"kO") ;
  gSystem->CompileMacro("betaHelperFunctions.h"                       ,"kO") ;
  gSystem->CompileMacro("RooNormalFromFlatPdf.cxx"                    ,"kO") ;
  gSystem->CompileMacro("RooBetaInverseCDF.cxx"                       ,"kO") ;
  gSystem->CompileMacro("RooBetaPrimeInverseCDF.cxx"                  ,"kO") ;
  gSystem->CompileMacro("RooCorrelatedBetaGeneratorHelper.cxx"        ,"kO") ;
  gSystem->CompileMacro("RooCorrelatedBetaPrimeGeneratorHelper.cxx"   ,"kO") ;
  gSystem->CompileMacro("rooFitBetaHelperFunctions.h"                 ,"kO") ;
  gSystem->CompileMacro("../3Dcode/RooPosDefCorrGauss.cxx"            ,"kO") ;
  gSystem->CompileMacro("rooFitGaussianHelperFunctions.h"             ,"kO") ;
  gSystem->CompileMacro("rooFitLogNormalHelperFunctions.h"            ,"kO") ;
  gSystem->CompileMacro("fixParameters.C"                             ,"kO") ; 
  gSystem->CompileMacro("metReweightingBuilderSIMPLETAU.C"            ,"kO") ;
  gSystem->CompileMacro("likelihoodBuilder.C"                         ,"kO") ;

  return;
}
