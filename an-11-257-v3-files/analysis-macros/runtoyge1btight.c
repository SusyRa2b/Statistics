

   void runtoyge1btight() {

      gROOT->LoadMacro("ra2bRoostatsClass4ln.c+") ;


      // ra2bRoostatsClass4ln( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass4ln rfit(1,0,2) ;  //+++ VERY IMPORTANT!  3rd argument must be 2 for tight sel!
      rfit.initialize("an-11-257-v3-files/input-files/mc-inputs-SMonly-tight-sel.txt") ;
      rfit.setSusyScanPoint("an-11-257-v3-files/input-files/signalSyst.LM9.ge1bTight.1143invpb.dat",0,0) ;






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





      TFile clsroot("an-11-257-v3-files/output-files/toyge1btight-smonly-mctest.root","recreate") ;











      rfit.reinitialize() ;
      rfit.setAndFixSusySig(0.) ;

      rfit.doFit() ;
      printf("\n\n ======== best fit, no susy included ==============\n\n") ;
      rfit.parameterSnapshot() ;

      rfit.saveToymeanSnapshot() ;

      rfit.doToyStudyNoSusyInFit(1000, "an-11-257-v3-files/input-files/mc-trueval-SMonly-tight-sel.txt" ) ;

      //------------------








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

      clsroot.Close() ;



   }



