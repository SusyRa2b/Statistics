
  void toywrapper2b( const char* input_datfile = "datfiles/Input-met4-ht4-wsyst1.dat",
                const char* input_susyfile = "datfiles/Susy-range1-met4-ht4.dat",
                double input_mgl=900, double input_mlsp=300.,
                const char* input_deffdbtagfile = "datfiles/dummy_DeffDbtag-met4-ht4.dat",
                double input_nSusy0lep = 80.,
                const char* input_outputDir = "output-toymc2b-mgl900-mlsp300-80evts-wsyst1-exp0lep",
                int nToy = 400
      ) {
     gROOT->LoadMacro("RooRatio.cxx+") ;
     gROOT->LoadMacro("toymc2b.c+") ;
     toymc2b( input_datfile, input_susyfile, input_mgl, input_mlsp, input_deffdbtagfile, input_nSusy0lep, input_outputDir, nToy ) ;
  }

