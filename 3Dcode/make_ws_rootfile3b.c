


   void make_ws_rootfile3b( const char* input_datfile = "output-toymc3b-debug1/Input-met3-ht3-v5-wclosuresyst-wttwjcorr-wqcdcorr-toy0000.dat",
                           const char* input_susyfile = "datfiles/T1bbbb-met3-ht3-v5.dat",
                           double mgl=1100, double mlsp=700.,
                           const char* input_deffdbtagfile = "datfiles/dummy_DeffDbtag-met3-ht3.dat",
                           const char* input_deffdbtagfile_lightflavor = "datfiles/dummy_DeffDbtag-met3-ht3.dat",
                           int qcdModelIndex = 3,
                           const char* wsrootfilename = "ws3b.root",
                           const char* blindBinsList = "null",
                           const char* systfilename = "datfiles/sigsystematics.T1bbbb.JES.txt",
                           const char* pdf_syst_file = "datfiles/sigsystematics.T1bbbb.PDF.txt",
                           const char* wjets_xsec_shapesyst_file = "blah.txt",
                           const char* singletop_xsec_shapesyst_file = "blah.txt"
                         ) {

       gROOT->LoadMacro("RooRatio.cxx+") ;
       gROOT->LoadMacro("RooBetaPdf.cxx+") ;
       gROOT->LoadMacro("RooPosDefCorrGauss.cxx+") ;
       gROOT->LoadMacro("ra2bRoostatsClass3D_3b.c+") ;

       ra2bRoostatsClass3D_3b ra2b ;

       ra2b.initialize( input_datfile, input_susyfile, mgl, mlsp, false, 0., input_deffdbtagfile, input_deffdbtagfile_lightflavor, qcdModelIndex, wsrootfilename, blindBinsList, systfilename, pdf_syst_file, wjets_xsec_shapesyst_file, singletop_xsec_shapesyst_file ) ;


   }

