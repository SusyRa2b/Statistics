
  void toywrapper3b( const char* input_datfile = "datfiles/Input-met3-ht3-v1-wclosuresyst-wttwjcorr-wqcdcorr.dat",
                const char* input_susyfile = "datfiles/T1bbbb-met3-ht3-v1-flat10perr.dat",
                double input_mgl=800, double input_mlsp=600.,
                const char* input_deffdbtagfile = "datfiles/dummy_DeffDbtag-met3-ht3.dat",
                double input_nSusy0lep = 37.,
                const char* input_outputDir = "output-toymc3b-800-600-37evts-flat10perr_posdef",
                int nToy = 400,
                const char* input_mcvals_rootfile = "rootfiles/gi-plots-met3-ht3-v1.root",
                bool useExpected0lep = true
      ) {
     gROOT->LoadMacro("RooRatio.cxx+") ;
     gROOT->LoadMacro("RooPosDefCorrGauss.cxx+") ;
     gROOT->LoadMacro("toymc3b.c+") ;
     toymc3b( input_datfile, input_susyfile, input_mgl, input_mlsp, input_deffdbtagfile, input_nSusy0lep, input_outputDir, nToy, input_mcvals_rootfile, useExpected0lep ) ;
  }

