

   void an_scanplot_unblind_t1bbbb_withcontam_ge1b_loose() {

      gROOT->LoadMacro("ra2bRoostatsClass4.c+") ;

      float low,high ;


      ra2bRoostatsClass4 rfit(0,1,1) ;
      rfit.initialize("an-11-257-v3-files/input-files/byhand-data-ge1b-loose.txt") ;
      rfit.doFit() ;
      rfit.parameterSnapshot() ;
      rfit.susyScanWithContam("an-11-257-v3-files/input-files/signalSyst.T1bbbb-preliminary.ge1bLoose.dat",
                              "an-11-257-v3-files/output-files/an-scanplot-unblind-t1bbbb-withcontam-ge1b-loose", 1 ) ;


   }
