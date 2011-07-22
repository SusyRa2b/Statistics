

   void discoveryplot_tb10_withcontam_tight_sel() {

      gROOT->LoadMacro("ra2bRoostatsClass3.c+") ;

      // ra2bRoostatsClass3( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass3 rfit(1,0,2) ;
      rfit.initialize("input-files/mc-inputs-SMonly-tight-sel.txt") ;
      rfit.discoveryScanWithContam("input-files/susyScanInput4-tanb10-SSVHPT-nj3-nb1-ht500-sb150_200-sig300.txt") ;


   }
