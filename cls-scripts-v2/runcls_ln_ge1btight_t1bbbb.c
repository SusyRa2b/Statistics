

   void runcls_ln_ge1btight_t1bbbb( int m0 = 560, int m12 = 300, double t1bbbbXsec = 1.0 ) {

      gROOT->LoadMacro("ra2bRoostatsClass4ln.c+") ;



      ra2bRoostatsClass4ln rfit(1,0,2) ; // --- Important!   3rd argument must be 2 (for alternative Znn model)

      rfit.initialize("an-11-257-v3-files/input-files/byhand-data-ge1b-tight.txt") ;
      rfit.setSusyScanPoint("an-11-257-v3-files/input-files/signalSyst.T1bbbb.ge1bTight.dat", m0, m12, true, t1bbbbXsec ) ;






      rfit.doFit() ;
      printf("\n\n ======== best fit, susy floating ==============\n\n") ;
      rfit.parameterSnapshot() ;

      double maxLogL = rfit.getLogLikelihoodValue() ;







      rfit.reinitialize() ;
      rfit.setAndFixSusySigToPredictedValue() ;

      rfit.doFit() ;
      printf("\n\n ======== best fit, susy fixed to predicted value ==============\n\n") ;
      rfit.parameterSnapshot() ;

      double splusbLogL = rfit.getLogLikelihoodValue() ;





      double data_q = 2.0*( splusbLogL - maxLogL ) ;
      printf("\n\n ++++++ Data value of test statistic : %8.2f\n", data_q ) ;





 ///  char outfname[10000] ;
 ///  sprintf( outfname, "output-files-ge1btight-t1bbbb/clstest_ge1btight_tb40_%d_%d.root", m0, m12 ) ;
 ///  TFile clsroot( outfname,"recreate") ;









      rfit.reinitialize() ;
      rfit.setAndFixSusySig(0.) ;

      rfit.doFit() ;
      printf("\n\n ======== best fit, no susy included ==============\n\n") ;
      rfit.parameterSnapshot() ;

      rfit.saveToymeanSnapshot() ;

      double bgonlyPvalue = rfit.doToyStudy(1000,1, data_q ) ;

      printf(" +++++++++++ p-value of BG-only hypothesis : %8.3f\n", bgonlyPvalue ) ;

      //------------------





      rfit.reinitialize() ;
      rfit.setAndFixSusySigToPredictedValue() ;

      rfit.doFit() ;
      printf("\n\n ======== best fit, susy fixed to predicted value ==============\n\n") ;
      rfit.parameterSnapshot() ;



      rfit.saveToymeanSnapshot() ;

      double splusbPvalue = rfit.doToyStudy(1000,0, data_q ) ;

      printf(" +++++++++++ p-value of SIG+BG hypothesis : %8.3f\n", splusbPvalue ) ;

      double cls = 0. ;
      if ( bgonlyPvalue > 0. ) { cls = splusbPvalue / bgonlyPvalue ; }

      printf("\n\n +++++++ CLs value = %8.3f\n\n\n", cls ) ;

 ///  clsroot.Close() ;

      printf("\n\n  final result  %3d  %3d  %8.3f  %8.3f  %8.3f  %9.2f\n",
               m0, m12, cls, splusbPvalue, bgonlyPvalue, t1bbbbXsec ) ;


   }



