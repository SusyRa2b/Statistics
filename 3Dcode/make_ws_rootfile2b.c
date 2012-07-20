


   void make_ws_rootfile2b( const char* input_datfile = "datfiles/Input-met4-ht4-wsyst1.dat",
                           const char* input_susyfile = "datfiles/Susy-mgl900-mlsp300-met4-ht4.dat",
                           double mgl=900, double mlsp=300.,
                           const char* input_deffdbtagfile = "datfiles/dummy_DeffDbtag-met4-ht4.dat"
                         ) {

       gROOT->LoadMacro("RooRatio.cxx+") ;
       gROOT->LoadMacro("ra2bRoostatsClass3D_2b.c+") ;

       ra2bRoostatsClass3D_2b ra2b ;

       int qcdModelIndex = 2 ;
       ra2b.initialize( input_datfile, input_susyfile, mgl, mlsp, false, 0., input_deffdbtagfile, qcdModelIndex ) ;

       printf("\n\n Renaming ws.root as ws2.root\n\n") ;
       gSystem->Exec("mv ws.root ws2.root") ;


   }

