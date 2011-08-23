

   void an_scanplot_unblind_t1bbbb_withcontam_ge2b_loose() {

      gROOT->LoadMacro("ra2bRoostatsClass4ln.c+") ;

      float low,high ;

      ra2bRoostatsClass4ln rfit(0,1,2) ; //-- important!  3rd argument must be 2 for "tight" Znn model.
      rfit.initialize("an-11-257-v3-files/input-files/byhand-data-ge2b-loose.txt") ;
      rfit.doFit() ;
      rfit.parameterSnapshot() ;
      rfit.susyScanWithContam("an-11-257-v3-files/input-files/signalSyst.T1bbbb-preliminary.ge2bLoose.dat",
                              "an-11-257-v3-files/output-files/an-scanplot-unblind-t1bbbb-withcontam-ge2b-loose", 1 ) ;


   }
