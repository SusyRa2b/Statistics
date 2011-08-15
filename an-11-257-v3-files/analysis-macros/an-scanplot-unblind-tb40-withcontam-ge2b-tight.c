

   void an_scanplot_unblind_tb40_withcontam_ge2b_tight() {

      gROOT->LoadMacro("ra2bRoostatsClass4.c+") ;

      float low,high ;

      // ra2bRoostatsClass4( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass4 rfit(0,1,1) ;
      rfit.initialize("an-11-257-v3-files/input-files/byhand-data-ge2b-tight.txt") ;
      rfit.doFit() ;
      rfit.parameterSnapshot() ;
      rfit.susyScanWithContam("an-11-257-v3-files/input-files/signalSyst.mSUGRAtanb40.ge2btight.dat",
                              "an-11-257-v3-files/output-files/an-scanplot-unblind-tb40-withcontam-ge2b-tight" ) ;


   }
