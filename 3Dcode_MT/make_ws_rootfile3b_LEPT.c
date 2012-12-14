


   void make_ws_rootfile3b_LEPT( const char* input_datfile = "output-toymc3b-debug1/Input-met3-ht3-v5-wclosuresyst-wttwjcorr-wqcdcorr-toy0000.dat",
				 const char* input_susyfile = "datfiles/T1bbbb-met3-ht3-v5.dat",
				 double mgl=1100, double mlsp=700.,
				 const char* input_deffdbtagfile = "datfiles/dummy_DeffDbtag-met3-ht3.dat",
				 const char* input_lightmistagfile = "datfiles/dummy_DeffDbtag-met3-ht3.dat",
				 int qcdModelIndex = 3,
				 const char* wsrootfilename = "ws3b.root",
				 const char* blindBinsList = "null",
				 bool constrainBjetShape = false,
				 bool floatSLSigRatios = false,
				 const char* systfilename = "datfiles/sigsystematics.T1bbbb.JES.txt",
				 const char* systfilename2 = "datfiles/sigsystematics.T1bbbb.PDF.txt"
				 ) {

     gROOT->LoadMacro("RooRatio.cxx+") ;
     gROOT->LoadMacro("RooBetaPdf.cxx+") ;
     gROOT->LoadMacro("RooPosDefCorrGauss.cxx+") ;
     gROOT->LoadMacro("RooProdPdfLogSum.cxx+") ;
     gROOT->LoadMacro("ra2bRoostatsClass3D_3b_LEPT.c+") ;

     ra2bRoostatsClass3D_3b_LEPT ra2b ;

     ra2b.initialize( input_datfile, input_susyfile, mgl, mlsp, false, 0., input_deffdbtagfile, input_lightmistagfile, qcdModelIndex, 
		      wsrootfilename, blindBinsList, constrainBjetShape, floatSLSigRatios, systfilename, systfilename2 ) ;


   }

