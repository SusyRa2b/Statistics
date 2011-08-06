

   void an_scanplot_unblind_tb50_withcontam_ge1b_tight() {

      gROOT->LoadMacro("ra2bRoostatsClass3.c+") ;

      float low,high ;

      // ra2bRoostatsClass3( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass3 rfit(1,0,2) ;
      rfit.initialize("an-11-257-v2-files/input-files/byhand-data-ge1b-tight.txt") ;
      rfit.doFit() ;
      rfit.parameterSnapshot() ;
      rfit.susyScanWithContam("an-11-257-v2-files/input-files/susyScanInput4-tanb50-SSVHPT-nj3-nb1-ht500-sb150_200-sig300.txt",
                              "an-11-257-v2-files/output-files/an-scanplot-unblind-tb50-withcontam-ge1b-tight" ) ;


   }
