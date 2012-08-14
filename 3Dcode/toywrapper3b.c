
  void toywrapper3b( const char* input_datfile = "datfiles/Input-met3-ht3-v1-wclosuresyst-wttwjcorr-wqcdcorr.dat",
                const char* input_susyfile = "datfiles/T1bbbb-met3-ht3-v1.dat",
                double input_mgl=800, double input_mlsp=700.,
                const char* input_deffdbtagfile = "datfiles/dummy_DeffDbtag-met3-ht3.dat",
                double input_nSusy0lep = 100.,
                const char* input_outputDir = "output-toymc3b-mgl800-mlsp700-100evts-wsyst-useEZLtrue-test2g",
                int nToy = 400,
                const char* input_mcvals_rootfile = "rootfiles/gi-plots-met3-ht3-v1.root",
                bool useExpected0lep = true
      ) {
     gROOT->LoadMacro("RooRatio.cxx+") ;
     gROOT->LoadMacro("toymc3b.c+") ;
     toymc3b( input_datfile, input_susyfile, input_mgl, input_mlsp, input_deffdbtagfile, input_nSusy0lep, input_outputDir, nToy, input_mcvals_rootfile, useExpected0lep ) ;
  }

