
  void toywrapper1( const char* input_datfile = "datfiles/Input-met4-ht4-wsyst1.dat",
                const char* input_susyfile = "datfiles/Susy-mgl900-mlsp300-met4-ht4.dat",
                double input_mgl=900, double input_mlsp=300.,
                const char* input_deffdbtagfile = "datfiles/dummy_DeffDbtag-met4-ht4.dat",
                double input_nSusy0lep = 60.,
                const char* input_outputDir = "output-toymc1",
                int nToy = 10
      ) {
     gROOT->LoadMacro("RooRatio.cxx+") ;
     gROOT->LoadMacro("toymc1.c+") ;
     toymc1( input_datfile, input_susyfile, input_mgl, input_mlsp, input_deffdbtagfile, input_nSusy0lep, input_outputDir, nToy ) ;
  }

