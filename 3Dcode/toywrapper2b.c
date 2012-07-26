
  void toywrapper2b( const char* input_datfile = "datfiles/Input-met4-ht4-newMC-wsyst.dat",
                const char* input_susyfile = "datfiles/Susy-range1-met4-ht4.dat",
                double input_mgl=850, double input_mlsp=600.,
                const char* input_deffdbtagfile = "datfiles/dummy_DeffDbtag-met4-ht4.dat",
                double input_nSusy0lep = 100.,
                const char* input_outputDir = "output-toymc2b-mgl850-mlsp600-100evts-newMC-wsyst-test2d",
                int nToy = 50,
                const char* input_mcvals_rootfile = "rootfiles/gi-plots-met4-ht4-newMC.root"
      ) {
     gROOT->LoadMacro("RooRatio.cxx+") ;
     gROOT->LoadMacro("toymc2b.c+") ;
     toymc2b( input_datfile, input_susyfile, input_mgl, input_mlsp, input_deffdbtagfile, input_nSusy0lep, input_outputDir, nToy, input_mcvals_rootfile ) ;
  }

