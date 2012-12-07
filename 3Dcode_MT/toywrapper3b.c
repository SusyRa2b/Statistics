
void toywrapper3b( const char* input_datfile = "datfiles/Input-met4-ht4-v15-standard.dat",
		   const char* input_susyfile = "datfiles/T1tttt-met4-ht4-v15-njets3.dat",
		   double input_mgl=1025, double input_mlsp=500,
		   const char* input_deffdbtagfile = "datfiles/dummy_DeffDbtag-met4-ht4.dat",
		   double input_nSusy0lep = 90.,
		   const char* input_outputDir = "FTest_exp0lep_P3_FLOAT_nSusy90",
		   int nToy = 100,
		   const char* input_mcvals_rootfile = "rootfiles/gi-plots-met4-ht4-v15.root",
		   bool useExpected0lep = true,
		   int qcdModelIndex = 4,
		   bool input_doSignif = false,
		   bool input_doUL = false,
		   const char* blindBinsList = "null",
		   bool input_inputObservablesArePostTrigger = true,
		   bool constrainBjetShape = false,
		   bool floatSLSigRatios = true,
		   const char* systfilename = "datfiles/DummySyst.txt"
		   ) {
  gROOT->LoadMacro("RooRatio.cxx+") ;
  gROOT->LoadMacro("RooBetaPdf.cxx+") ;
  gROOT->LoadMacro("RooPosDefCorrGauss.cxx+") ;
  gROOT->LoadMacro("RooProdPdfLogSum.cxx+") ;
  gROOT->LoadMacro("toymc3b.c+") ;

  toymc3b( input_datfile, input_susyfile, input_mgl, input_mlsp, input_deffdbtagfile, input_nSusy0lep, input_outputDir, nToy, input_mcvals_rootfile, useExpected0lep, qcdModelIndex, input_doSignif, input_doUL, blindBinsList, input_inputObservablesArePostTrigger, constrainBjetShape, floatSLSigRatios, systfilename ) ;
  
}

