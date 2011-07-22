

   void scanplot_unblind_tb10_nocontam_tight_sel() {

      gROOT->LoadMacro("ra2bRoostatsClass3.c+") ;

      float low,high ;

      // ra2bRoostatsClass3( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass3 rfit(1,0,2) ;
      rfit.initialize("input-files/rsinput-v4-SSVHPT-nj3-nb1-ht500-sb150_200-sig300.txt") ;
      rfit.doFit() ;
      rfit.parameterSnapshot() ;
      rfit.susyScanNoContam("input-files/susyScanInput4-tanb10-SSVHPT-nj3-nb1-ht500-sb150_200-sig300.txt","-scanplot-nocontam-unblind.png" ) ;


   }
