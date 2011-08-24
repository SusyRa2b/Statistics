

   void an_scanplot_unblind_t2bb_withcontam_ge1b_loose() {

      gROOT->LoadMacro("ra2bRoostatsClass4ln.c+") ;

      float low,high ;

      ra2bRoostatsClass4ln rfit(0,1,1) ;
      rfit.initialize("an-11-257-v3-files/input-files/byhand-data-ge1b-loose.txt") ;
      rfit.doFit() ;
      rfit.parameterSnapshot() ;
      rfit.susyScanWithContam("an-11-257-v3-files/input-files/signalSyst.T2bb.ge1bLoose.dat",
                              "an-11-257-v3-files/output-files/an-scanplot-unblind-t2bb-withcontam-ge1b-loose", 1 ) ;


   }
