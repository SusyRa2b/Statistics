

   void MCSMsusy_lm9_tight_sel() {

      gROOT->LoadMacro("ra2bRoostatsClass3.c+") ;

      float low,high ;
      float susylow, susyhigh ;

      // ra2bRoostatsClass3( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass3 rfitNoSusy(1,0,2) ;
      rfitNoSusy.initialize("an-11-257-v2-files/input-files/mc-inputs-SM+lm9-tight-sel.txt") ;
      rfitNoSusy.setAndFixSusySig(0.) ;

      rfitNoSusy.doFit() ;
      printf("\n\n ======== best fit, no susy included ==============\n\n") ;
      rfitNoSusy.parameterSnapshot() ;
      rfitNoSusy.fitQualityPlot(0,"an-11-257-v2-files/output-files/fitQual-MCSMsusy-BestFit-NoSusy-tight-sel.png") ;
      rfitNoSusy.fitQualityPlot(1,"an-11-257-v2-files/output-files/fitQualNorm-MCSMsusy-BestFit-NoSusy-tight-sel.png") ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profilettwjSig(low,high,1,"an-11-257-v2-files/output-files/ttwjSigProf-MCSMsusy-NoSusy-tight-sel.png", 45. ) ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profileqcdSig(low,high,1,"an-11-257-v2-files/output-files/qcdSigProf-MCSMsusy-NoSusy-tight-sel.png", 5. ) ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profileZnnSig(low,high,1,"an-11-257-v2-files/output-files/znnSigProf-MCSMsusy-NoSusy-tight-sel.png", 10. ) ;









   //-------------------------

      ra2bRoostatsClass3 rfitSusyFloat(1,0,2) ;
      rfitSusyFloat.initialize("an-11-257-v2-files/input-files/mc-inputs-SM+lm9-tight-sel.txt") ;
      rfitSusyFloat.setSusyScanPoint("an-11-257-v2-files/input-files/susy-lm9-nb1-ht500-sb150_200-sig300.txt",0,0) ;

      rfitSusyFloat.doFit() ;
      printf("\n\n ======== best fit, susy floating ==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v2-files/output-files/fitQual-MCSMsusy-BestFit-susyFloat-lm9-tight-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v2-files/output-files/fitQualNorm-MCSMsusy-BestFit-susyFloat-lm9-tight-sel.png") ;
      rfitSusyFloat.profileSusySig(susylow,susyhigh,1,"an-11-257-v2-files/output-files/susySigProf-MCSMsusy-susyFloat-lm9-tight-sel.png", 45. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profilettwjSig(low,high,1,"an-11-257-v2-files/output-files/ttwjSigProf-MCSMsusy-susyFloat-lm9-tight-sel.png", 45. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profileqcdSig(low,high,1,"an-11-257-v2-files/output-files/qcdSigProf-MCSMsusy-susyFloat-lm9-tight-sel.png", 5. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profileZnnSig(low,high,1,"an-11-257-v2-files/output-files/znnSigProf-MCSMsusy-susyFloat-lm9-tight-sel.png", 10. ) ;


      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.setAndFixSusySig(susylow) ;
      rfitSusyFloat.doFit() ;
      printf("\n\n ======== lower limit , susy floating==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v2-files/output-files/fitQual-MCSMsusy-susyFloat-LowerLimit-lm9-tight-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v2-files/output-files/fitQualNorm-MCSMsusy-susyFloat-LowerLimit-lm9-tight-sel.png") ;

      rfitSusyFloat.setAndFixSusySig(susyhigh) ;
      rfitSusyFloat.doFit() ;
      printf("\n\n ======== upper limit , susy floating==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v2-files/output-files/fitQual-MCSMsusy-susyFloat-UpperLimit-lm9-tight-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v2-files/output-files/fitQualNorm-MCSMsusy-susyFloat-UpperLimit-lm9-tight-sel.png") ;


   //-------------------------

      ra2bRoostatsClass3 rfitSusyFix(1,0,2) ;
      rfitSusyFix.initialize("an-11-257-v2-files/input-files/mc-inputs-SM+lm9-tight-sel.txt") ;
      rfitSusyFix.setSusyScanPoint("an-11-257-v2-files/input-files/susy-lm9-nb1-ht500-sb150_200-sig300.txt",0,0) ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;

      rfitSusyFix.doFit() ;
      printf("\n\n ======== best fit, susy fixed ==============\n\n") ;
      rfitSusyFix.parameterSnapshot() ;
      rfitSusyFix.fitQualityPlot(0,"an-11-257-v2-files/output-files/fitQual-MCSMsusy-BestFit-susyFixed-lm9-tight-sel.png") ;
      rfitSusyFix.fitQualityPlot(1,"an-11-257-v2-files/output-files/fitQualNorm-MCSMsusy-BestFit-susyFixed-lm9-tight-sel.png") ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profilettwjSig(low,high,1,"an-11-257-v2-files/output-files/ttwjSigProf-MCSMsusy-susyFixed-lm9-tight-sel.png", 45. ) ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profileqcdSig(low,high,1,"an-11-257-v2-files/output-files/qcdSigProf-MCSMsusy-susyFixed-lm9-tight-sel.png", 5. ) ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profileZnnSig(low,high,1,"an-11-257-v2-files/output-files/znnSigProf-MCSMsusy-susyFixed-lm9-tight-sel.png", 10. ) ;




   }
