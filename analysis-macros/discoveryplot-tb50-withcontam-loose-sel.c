

   void discoveryplot_tb50_withcontam_loose_sel() {

      gROOT->LoadMacro("ra2bRoostatsClass3.c+") ;

      // ra2bRoostatsClass3( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass3 rfit(1,0) ;
      rfit.initialize("input-files/mc-inputs-SMonly-loose-sel.txt") ;
      rfit.discoveryScanWithContam("input-files/susyScanInput4-tanb50-SSVHPT-nj3-nb1-ht350-sb150_200-sig200.txt") ;


   }
