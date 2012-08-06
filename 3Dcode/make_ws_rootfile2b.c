


   void make_ws_rootfile2b( const char* input_datfile = "datfiles/Input-met3-ht3-v1-wclosuresyst.dat",
                           const char* input_susyfile = "datfiles/T1bbbb-met3-ht3-v1.dat",
                           double mgl=800, double mlsp=700.,
                           const char* input_deffdbtagfile = "datfiles/dummy_DeffDbtag-met3-ht3.dat"
                         ) {

       gROOT->LoadMacro("RooRatio.cxx+") ;
       gROOT->LoadMacro("ra2bRoostatsClass3D_2b.c+") ;

       ra2bRoostatsClass3D_2b ra2b ;

       int qcdModelIndex = 2 ;
       ra2b.initialize( input_datfile, input_susyfile, mgl, mlsp, false, 0., input_deffdbtagfile, qcdModelIndex ) ;

       printf("\n\n Renaming ws.root as ws2.root\n\n") ;
       gSystem->Exec("mv ws.root ws2.root") ;


   }

