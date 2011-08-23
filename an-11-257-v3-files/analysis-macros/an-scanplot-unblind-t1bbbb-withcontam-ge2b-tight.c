

   void an_scanplot_unblind_t1bbbb_withcontam_ge2b_tight() {

      gROOT->LoadMacro("ra2bRoostatsClass4ln.c+") ;

      float low,high ;

      // ra2bRoostatsClass4ln( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass4ln rfit(0,1,2) ; //-- important!  3rd argument must be 2 for "tight" Znn model.
      rfit.initialize("an-11-257-v3-files/input-files/byhand-data-ge2b-tight.txt") ;
      rfit.doFit() ;
      rfit.parameterSnapshot() ;
      rfit.susyScanWithContam("an-11-257-v3-files/input-files/signalSyst.T1bbbb-preliminary.ge2btight.dat",
                              "an-11-257-v3-files/output-files/an-scanplot-unblind-t1bbbb-withcontam-ge2b-tight", 1 ) ;


   }
