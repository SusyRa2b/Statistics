

   void MCSMonly_lm9_loose_sel() {

      gROOT->LoadMacro("ra2bRoostatsClass4ln.c+") ;

      float low,high ;
      float susylow, susyhigh ;

      // ra2bRoostatsClass4ln( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass4ln rfitNoSusy(1,0,1) ;
      rfitNoSusy.initialize("an-11-257-v3-files/input-files/mc-inputs-SMonly-loose-sel.txt") ;
      rfitNoSusy.setAndFixSusySig(0.) ;

      rfitNoSusy.doFit() ;
      printf("\n\n ======== best fit, no susy included ==============\n\n") ;
      rfitNoSusy.parameterSnapshot() ;
      rfitNoSusy.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-MCSMonly-BestFit-NoSusy-loose-sel.png") ;
      rfitNoSusy.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-MCSMonly-BestFit-NoSusy-loose-sel.png") ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profilettwjSig(low,high,1,"an-11-257-v3-files/output-files/ttwjSigProf-MCSMonly-NoSusy-loose-sel.png", 200. ) ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profileqcdSig(low,high,1,"an-11-257-v3-files/output-files/qcdSigProf-MCSMonly-NoSusy-loose-sel.png", 100. ) ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profileZnnSig(low,high,1,"an-11-257-v3-files/output-files/znnSigProf-MCSMonly-NoSusy-loose-sel.png", 80. ) ;









   //-------------------------

      ra2bRoostatsClass4ln rfitSusyFloat(1,0,1) ;
      rfitSusyFloat.initialize("an-11-257-v3-files/input-files/mc-inputs-SMonly-loose-sel.txt") ;
      rfitSusyFloat.setSusyScanPoint("an-11-257-v3-files/input-files/signalSyst.LM9.ge1bLoose.1143invpb.dat",0,0) ;

      rfitSusyFloat.doFit() ;
      printf("\n\n ======== best fit, susy floating ==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-MCSMonly-BestFit-susyFloat-lm9-loose-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-MCSMonly-BestFit-susyFloat-lm9-loose-sel.png") ;
      rfitSusyFloat.profileSusySig(susylow,susyhigh,1,"an-11-257-v3-files/output-files/susySigProf-MCSMonly-susyFloat-lm9-loose-sel.png", 200. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profilettwjSig(low,high,1,"an-11-257-v3-files/output-files/ttwjSigProf-MCSMonly-susyFloat-lm9-loose-sel.png", 200. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profileqcdSig(low,high,1,"an-11-257-v3-files/output-files/qcdSigProf-MCSMonly-susyFloat-lm9-loose-sel.png", 100. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profileZnnSig(low,high,1,"an-11-257-v3-files/output-files/znnSigProf-MCSMonly-susyFloat-lm9-loose-sel.png", 80. ) ;


      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.setAndFixSusySig(susylow) ;
      rfitSusyFloat.doFit() ;
      printf("\n\n ======== lower limit , susy floating==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-MCSMonly-susyFloat-LowerLimit-lm9-loose-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-MCSMonly-susyFloat-LowerLimit-lm9-loose-sel.png") ;

      rfitSusyFloat.setAndFixSusySig(susyhigh) ;
      rfitSusyFloat.doFit() ;
      printf("\n\n ======== upper limit , susy floating==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-MCSMonly-susyFloat-UpperLimit-lm9-loose-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-MCSMonly-susyFloat-UpperLimit-lm9-loose-sel.png") ;


   //-------------------------

      ra2bRoostatsClass4ln rfitSusyFix(1,0,1) ;
      rfitSusyFix.initialize("an-11-257-v3-files/input-files/mc-inputs-SMonly-loose-sel.txt") ;
      rfitSusyFix.setSusyScanPoint("an-11-257-v3-files/input-files/signalSyst.LM9.ge1bLoose.1143invpb.dat",0,0) ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;

      rfitSusyFix.doFit() ;
      printf("\n\n ======== best fit, susy fixed ==============\n\n") ;
      rfitSusyFix.parameterSnapshot() ;
      rfitSusyFix.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-MCSMonly-BestFit-susyFixed-lm9-loose-sel.png") ;
      rfitSusyFix.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-MCSMonly-BestFit-susyFixed-lm9-loose-sel.png") ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profilettwjSig(low,high,1,"an-11-257-v3-files/output-files/ttwjSigProf-MCSMonly-susyFixed-lm9-loose-sel.png", 200. ) ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profileqcdSig(low,high,1,"an-11-257-v3-files/output-files/qcdSigProf-MCSMonly-susyFixed-lm9-loose-sel.png", 100. ) ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profileZnnSig(low,high,1,"an-11-257-v3-files/output-files/znnSigProf-MCSMonly-susyFixed-lm9-loose-sel.png", 80. ) ;




   }
