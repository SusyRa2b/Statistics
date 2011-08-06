

   void an_data_ge2b_tight() {

      gROOT->LoadMacro("ra2bRoostatsClass3.c+") ;

      float low,high ;
      float susylow, susyhigh ;

      // ra2bRoostatsClass3( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass3 rfitNoSusy(1,0,2) ;
      rfitNoSusy.initialize("an-11-257-v2-files/input-files/byhand-data-ge2b-tight.txt") ;
      rfitNoSusy.setAndFixSusySig(0.) ;

      rfitNoSusy.doFit() ;
      printf("\n\n ======== best fit, no susy included ==============\n\n") ;
      rfitNoSusy.parameterSnapshot() ;
      rfitNoSusy.fitQualityPlot(0,"an-11-257-v2-files/output-files/fitQual-unblind-BestFit-NoSusy-ge2b-tight-sel.png") ;
      rfitNoSusy.fitQualityPlot(1,"an-11-257-v2-files/output-files/fitQualNorm-unblind-BestFit-NoSusy-ge2b-tight-sel.png") ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profilettwjSig(low,high,1,"an-11-257-v2-files/output-files/ttwjSigProf-unblind-NoSusy-ge2b-tight-sel.png", 20. ) ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profileqcdSig(low,high,1,"an-11-257-v2-files/output-files/qcdSigProf-unblind-NoSusy-ge2b-tight-sel.png", 5. ) ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profileZnnSig(low,high,1,"an-11-257-v2-files/output-files/znnSigProf-unblind-NoSusy-ge2b-tight-sel.png", 10. ) ;



   //-------------------------

      ra2bRoostatsClass3 rfitSusyFloat(1,0,2) ;
      rfitSusyFloat.initialize("an-11-257-v2-files/input-files/byhand-data-ge2b-tight.txt") ;
      rfitSusyFloat.setSusyScanPoint("an-11-257-v2-files/input-files/susy-lm9-nb2-ht500-sb150_200-sig300.txt",0,0) ;

      rfitSusyFloat.doFit() ;
      printf("\n\n ======== best fit, susy floating ==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v2-files/output-files/fitQual-unblind-BestFit-susyFloat-lm9-ge2b-tight-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v2-files/output-files/fitQualNorm-unblind-BestFit-susyFloat-lm9-ge2b-tight-sel.png") ;
      rfitSusyFloat.profileSusySig(susylow,susyhigh,1,"an-11-257-v2-files/output-files/susySigProf-unblind-susyFloat-lm9-ge2b-tight-sel.png", 20. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profilettwjSig(low,high,1,"an-11-257-v2-files/output-files/ttwjSigProf-unblind-susyFloat-lm9-ge2b-tight-sel.png", 20. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profileqcdSig(low,high,1,"an-11-257-v2-files/output-files/qcdSigProf-unblind-susyFloat-lm9-ge2b-tight-sel.png", 5. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profileZnnSig(low,high,1,"an-11-257-v2-files/output-files/znnSigProf-unblind-susyFloat-lm9-ge2b-tight-sel.png", 10. ) ;


      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.setAndFixSusySig(susylow) ;
      rfitSusyFloat.doFit() ;
      printf("\n\n ======== lower limit , susy floating==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v2-files/output-files/fitQual-unblind-susyFloat-LowerLimit-lm9-ge2b-tight-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v2-files/output-files/fitQualNorm-unblind-susyFloat-LowerLimit-lm9-ge2b-tight-sel.png") ;

      rfitSusyFloat.setAndFixSusySig(susyhigh) ;
      rfitSusyFloat.doFit() ;
      printf("\n\n ======== upper limit , susy floating==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v2-files/output-files/fitQual-unblind-susyFloat-UpperLimit-lm9-ge2b-tight-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v2-files/output-files/fitQualNorm-unblind-susyFloat-UpperLimit-lm9-ge2b-tight-sel.png") ;


   //-------------------------

      ra2bRoostatsClass3 rfitSusyFix(1,0,2) ;
      rfitSusyFix.initialize("an-11-257-v2-files/input-files/byhand-data-ge2b-tight.txt") ;
      rfitSusyFix.setSusyScanPoint("an-11-257-v2-files/input-files/susy-lm9-nb2-ht500-sb150_200-sig300.txt",0,0) ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;

      rfitSusyFix.doFit() ;
      printf("\n\n ======== best fit, susy fixed ==============\n\n") ;
      rfitSusyFix.parameterSnapshot() ;
      rfitSusyFix.fitQualityPlot(0,"an-11-257-v2-files/output-files/fitQual-unblind-BestFit-susyFixed-lm9-ge2b-tight-sel.png") ;
      rfitSusyFix.fitQualityPlot(1,"an-11-257-v2-files/output-files/fitQualNorm-unblind-BestFit-susyFixed-lm9-ge2b-tight-sel.png") ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profilettwjSig(low,high,1,"an-11-257-v2-files/output-files/ttwjSigProf-unblind-susyFixed-lm9-ge2b-tight-sel.png", 20. ) ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profileqcdSig(low,high,1,"an-11-257-v2-files/output-files/qcdSigProf-unblind-susyFixed-lm9-ge2b-tight-sel.png", 5. ) ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profileZnnSig(low,high,1,"an-11-257-v2-files/output-files/znnSigProf-unblind-susyFixed-lm9-ge2b-tight-sel.png", 10. ) ;









   }
