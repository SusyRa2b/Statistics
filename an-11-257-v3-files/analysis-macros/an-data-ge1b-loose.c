

   void an_data_ge1b_loose() {

      gROOT->LoadMacro("ra2bRoostatsClass4.c+") ;

      float low,high ;
      float susylow, susyhigh ;

      // ra2bRoostatsClass4( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass4 rfitNoSusy(1,0,1) ;
      rfitNoSusy.initialize("an-11-257-v3-files/input-files/byhand-data-ge1b-loose.txt") ;
      rfitNoSusy.setAndFixSusySig(0.) ;

      rfitNoSusy.doFit() ;
      printf("\n\n ======== best fit, no susy included ==============\n\n") ;
      rfitNoSusy.parameterSnapshot() ;
      rfitNoSusy.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-unblind-BestFit-NoSusy-ge1b-loose-sel.png") ;
      rfitNoSusy.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-unblind-BestFit-NoSusy-ge1b-loose-sel.png") ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profilettwjSig(low,high,1,"an-11-257-v3-files/output-files/ttwjSigProf-unblind-NoSusy-ge1b-loose-sel.png", 200. ) ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profileqcdSig(low,high,1,"an-11-257-v3-files/output-files/qcdSigProf-unblind-NoSusy-ge1b-loose-sel.png", 40. ) ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profileZnnSig(low,high,1,"an-11-257-v3-files/output-files/znnSigProf-unblind-NoSusy-ge1b-loose-sel.png", 80. ) ;



   //-------------------------

      ra2bRoostatsClass4 rfitSusyFloat(1,0,1) ;
      rfitSusyFloat.initialize("an-11-257-v3-files/input-files/byhand-data-ge1b-loose.txt") ;
      rfitSusyFloat.setSusyScanPoint("an-11-257-v3-files/input-files/signalSyst.LM9.ge1bLoose.1143invpb.dat",0,0) ;

      rfitSusyFloat.doFit() ;
      printf("\n\n ======== best fit, susy floating ==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-unblind-BestFit-susyFloat-lm9-ge1b-loose-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-unblind-BestFit-susyFloat-lm9-ge1b-loose-sel.png") ;
      rfitSusyFloat.profileSusySig(susylow,susyhigh,1,"an-11-257-v3-files/output-files/susySigProf-unblind-susyFloat-lm9-ge1b-loose-sel.png", 200. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profilettwjSig(low,high,1,"an-11-257-v3-files/output-files/ttwjSigProf-unblind-susyFloat-lm9-ge1b-loose-sel.png", 200. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profileqcdSig(low,high,1,"an-11-257-v3-files/output-files/qcdSigProf-unblind-susyFloat-lm9-ge1b-loose-sel.png", 40. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profileZnnSig(low,high,1,"an-11-257-v3-files/output-files/znnSigProf-unblind-susyFloat-lm9-ge1b-loose-sel.png", 80. ) ;


      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.setAndFixSusySig(susylow) ;
      rfitSusyFloat.doFit() ;
      printf("\n\n ======== lower limit , susy floating==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-unblind-susyFloat-LowerLimit-lm9-ge1b-loose-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-unblind-susyFloat-LowerLimit-lm9-ge1b-loose-sel.png") ;

      rfitSusyFloat.setAndFixSusySig(susyhigh) ;
      rfitSusyFloat.doFit() ;
      printf("\n\n ======== upper limit , susy floating==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-unblind-susyFloat-UpperLimit-lm9-ge1b-loose-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-unblind-susyFloat-UpperLimit-lm9-ge1b-loose-sel.png") ;


   //-------------------------

      ra2bRoostatsClass4 rfitSusyFix(1,0,1) ;
      rfitSusyFix.initialize("an-11-257-v3-files/input-files/byhand-data-ge1b-loose.txt") ;
      rfitSusyFix.setSusyScanPoint("an-11-257-v3-files/input-files/signalSyst.LM9.ge1bLoose.1143invpb.dat",0,0) ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;

      rfitSusyFix.doFit() ;
      printf("\n\n ======== best fit, susy fixed ==============\n\n") ;
      rfitSusyFix.parameterSnapshot() ;
      rfitSusyFix.fitQualityPlot(0,"an-11-257-v3-files/output-files/fitQual-unblind-BestFit-susyFixed-lm9-ge1b-loose-sel.png") ;
      rfitSusyFix.fitQualityPlot(1,"an-11-257-v3-files/output-files/fitQualNorm-unblind-BestFit-susyFixed-lm9-ge1b-loose-sel.png") ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profilettwjSig(low,high,1,"an-11-257-v3-files/output-files/ttwjSigProf-unblind-susyFixed-lm9-ge1b-loose-sel.png", 200. ) ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profileqcdSig(low,high,1,"an-11-257-v3-files/output-files/qcdSigProf-unblind-susyFixed-lm9-ge1b-loose-sel.png", 40. ) ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySigToPredictedValue() ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profileZnnSig(low,high,1,"an-11-257-v3-files/output-files/znnSigProf-unblind-susyFixed-lm9-ge1b-loose-sel.png", 80. ) ;









   }
