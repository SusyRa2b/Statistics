
  void toywrapper2b( const char* input_datfile = "datfiles/Input-met3-ht3-v1-wsyst1.dat",
                const char* input_susyfile = "datfiles/T1bbbb-met3-ht3-v1.dat",
                double input_mgl=800, double input_mlsp=300.,
                const char* input_deffdbtagfile = "datfiles/dummy_DeffDbtag-met3-ht3.dat",
                double input_nSusy0lep = 100.,
<<<<<<< toywrapper2b.c
                const char* input_outputDir = "output-toymc2b-mgl800-mlsp300-100evts-wsyst-test6e-useEZLtrue",
=======
                const char* input_outputDir = "output-toymc2b-mgl800-mlsp300-100evts-wsyst-test6d-useEZLtrue",
>>>>>>> 1.7
                int nToy = 100,
                const char* input_mcvals_rootfile = "rootfiles/gi-plots-met3-ht3.root"
      ) {
     gROOT->LoadMacro("RooRatio.cxx+") ;
     gROOT->LoadMacro("toymc2b.c+") ;
     toymc2b( input_datfile, input_susyfile, input_mgl, input_mlsp, input_deffdbtagfile, input_nSusy0lep, input_outputDir, nToy, input_mcvals_rootfile ) ;
  }

