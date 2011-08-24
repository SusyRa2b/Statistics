

   void MCSMsusy_lm9_loose_sel() {

      gROOT->LoadMacro("ra2bRoostatsClass4ln.c+") ;

      float low,high ;
      float susylow, susyhigh ;

      // ra2bRoostatsClass4ln( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass4ln rfitNoSusy(1,0,1) ;
      rfitNoSusy.initialize("an-11-257-v3-files/input-files/mc-inputs-SM+lm9-loose-sel.txt") ;
      rfitNoSusy.setAndFixSusySig(0.) ;

      rfitNoSusy.doFit() ;
      printf("\n\n ======== best fit, no susy included ==============\n\n") ;
      rfitNoSusy.parameterSnapshot() ;
      rfitNoSusy.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-MCSMsusy-BestFit-NoSusy-loose-sel.png") ;
      rfitNoSusy.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-MCSMsusy-BestFit-NoSusy-loose-sel.png") ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profilettwjSig(low,high,1,"an-11-257-v3-files/output-files/ttwjSigProf-MCSMsusy-NoSusy-loose-sel.png", 450. ) ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profileqcdSig(low,high,1,"an-11-257-v3-files/output-files/qcdSigProf-MCSMsusy-NoSusy-loose-sel.png", 100. ) ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profileZnnSig(low,high,1,"an-11-257-v3-files/output-files/znnSigProf-MCSMsusy-NoSusy-loose-sel.png", 100. ) ;









   //-------------------------

      ra2bRoostatsClass4ln rfitSusyFloat(1,0,1) ;
      rfitSusyFloat.initialize("an-11-257-v3-files/input-files/mc-inputs-SM+lm9-loose-sel.txt") ;
      rfitSusyFloat.setSusyScanPoint("an-11-257-v3-files/input-files/signalSyst.LM9.ge1bLoose.1143invpb.dat",0,0) ;

      rfitSusyFloat.doFit() ;
      printf("\n\n ======== best fit, susy floating ==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-MCSMsusy-BestFit-susyFloat-lm9-loose-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-MCSMsusy-BestFit-susyFloat-lm9-loose-sel.png") ;
      rfitSusyFloat.profileSusySig(susylow,susyhigh,1,"an-11-257-v3-files/output-files/susySigProf-MCSMsusy-susyFloat-lm9-loose-sel.png", 450. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profilettwjSig(low,high,1,"an-11-257-v3-files/output-files/ttwjSigProf-MCSMsusy-susyFloat-lm9-loose-sel.png", 450. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profileqcdSig(low,high,1,"an-11-257-v3-files/output-files/qcdSigProf-MCSMsusy-susyFloat-lm9-loose-sel.png", 100. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profileZnnSig(low,high,1,"an-11-257-v3-files/output-files/znnSigProf-MCSMsusy-susyFloat-lm9-loose-sel.png", 100. ) ;


      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.setAndFixSusySig(susylow) ;
      rfitSusyFloat.doFit() ;
      printf("\n\n ======== lower limit , susy floating==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-MCSMsusy-susyFloat-LowerLimit-lm9-loose-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-MCSMsusy-susyFloat-LowerLimit-lm9-loose-sel.png") ;

      rfitSusyFloat.setAndFixSusySig(susyhigh) ;
      rfitSusyFloat.doFit() ;
      printf("\n\n ======== upper limit , susy floating==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-MCSMsusy-susyFloat-UpperLimit-lm9-loose-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-MCSMsusy-susyFloat-UpperLimit-lm9-loose-sel.png") ;


   //-------------------------

      ra2bRoostatsClass4ln rfitSusyFix(1,0,1) ;
      rfitSusyFix.initialize("an-11-257-v3-files/input-files/mc-inputs-SM+lm9-loose-sel.txt") ;
      rfitSusyFix.setSusyScanPoint("an-11-257-v3-files/input-files/signalSyst.LM9.ge1bLoose.1143invpb.dat",0,0) ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;

      rfitSusyFix.doFit() ;
      printf("\n\n ======== best fit, susy fixed ==============\n\n") ;
      rfitSusyFix.parameterSnapshot() ;
      rfitSusyFix.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-MCSMsusy-BestFit-susyFixed-lm9-loose-sel.png") ;
      rfitSusyFix.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-MCSMsusy-BestFit-susyFixed-lm9-loose-sel.png") ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profilettwjSig(low,high,1,"an-11-257-v3-files/output-files/ttwjSigProf-MCSMsusy-susyFixed-lm9-loose-sel.png", 450. ) ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profileqcdSig(low,high,1,"an-11-257-v3-files/output-files/qcdSigProf-MCSMsusy-susyFixed-lm9-loose-sel.png", 100. ) ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profileZnnSig(low,high,1,"an-11-257-v3-files/output-files/znnSigProf-MCSMsusy-susyFixed-lm9-loose-sel.png", 100. ) ;




   }
