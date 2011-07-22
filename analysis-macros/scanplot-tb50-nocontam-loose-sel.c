

   void scanplot_tb50_nocontam_loose_sel() {

      gROOT->LoadMacro("ra2bRoostatsClass3.c+") ;

      float low,high ;

      // ra2bRoostatsClass3( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass3 rfit(1,0) ;
      rfit.initialize("input-files/mc-inputs-SMonly-loose-sel.txt") ;
      rfit.doFit() ;
      rfit.parameterSnapshot() ;
      rfit.susyScanNoContam("input-files/susyScanInput4-tanb50-SSVHPT-nj3-nb1-ht350-sb150_200-sig200.txt") ;


   }
