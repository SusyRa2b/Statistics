

   void unblind_tb50_m0_500_m12_250_loose_sel() {

      gROOT->LoadMacro("ra2bRoostatsClass3.c+") ;

      float low,high ;
      float susylow, susyhigh ;

      // ra2bRoostatsClass3( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

      ra2bRoostatsClass3 rfitNoSusy(1,0,1) ;
      rfitNoSusy.initialize("input-files/rsinput-v4-SSVHPT-nj3-nb1-ht350-sb150_200-sig200.txt") ;
      rfitNoSusy.setAndFixSusySig(0.) ;

      rfitNoSusy.doFit() ;
      printf("\n\n ======== best fit, no susy included ==============\n\n") ;
      rfitNoSusy.parameterSnapshot() ;
      rfitNoSusy.fitQualityPlot(0,"output-files/fitQual-unblind-BestFit-NoSusy-loose-sel.png") ;
      rfitNoSusy.fitQualityPlot(1,"output-files/fitQualNorm-unblind-BestFit-NoSusy-loose-sel.png") ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profilettwjSig(low,high,1,"output-files/ttwjSigProf-unblind-NoSusy-loose-sel.png", 120. ) ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profileqcdSig(low,high,1,"output-files/qcdSigProf-unblind-NoSusy-loose-sel.png", 120. ) ;
      rfitNoSusy.reinitialize() ;
      rfitNoSusy.setAndFixSusySig(0.) ;
      rfitNoSusy.doFit() ;
      rfitNoSusy.profileZnnSig(low,high,1,"output-files/znnSigProf-unblind-NoSusy-loose-sel.png", 120. ) ;









   //-------------------------

      ra2bRoostatsClass3 rfitSusyFloat(1,0,1) ;
      rfitSusyFloat.initialize("input-files/rsinput-v4-SSVHPT-nj3-nb1-ht350-sb150_200-sig200.txt") ;
      rfitSusyFloat.setSusyScanPoint("input-files/susyScanInput4-tanb50-SSVHPT-nj3-nb1-ht350-sb150_200-sig200.txt",500,250) ;

      rfitSusyFloat.doFit() ;
      printf("\n\n ======== best fit, susy floating ==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"output-files/fitQual-unblind-BestFit-susyFloat-tb50-m0-500-m12-250-loose-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"output-files/fitQualNorm-unblind-BestFit-susyFloat-tb50-m0-500-m12-250-loose-sel.png") ;
      rfitSusyFloat.profileSusySig(susylow,susyhigh,1,"output-files/susySigProf-unblind-susyFloat-tb50-m0-500-m12-250-loose-sel.png", 120. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profilettwjSig(low,high,1,"output-files/ttwjSigProf-unblind-susyFloat-tb50-m0-500-m12-250-loose-sel.png", 120. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profileqcdSig(low,high,1,"output-files/qcdSigProf-unblind-susyFloat-tb50-m0-500-m12-250-loose-sel.png", 120. ) ;
      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.profileZnnSig(low,high,1,"output-files/znnSigProf-unblind-susyFloat-tb50-m0-500-m12-250-loose-sel.png", 120. ) ;


      rfitSusyFloat.reinitialize() ;
      rfitSusyFloat.doFit() ;
      rfitSusyFloat.setAndFixSusySig(susylow) ;
      rfitSusyFloat.doFit() ;
      printf("\n\n ======== lower limit , susy floating==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"output-files/fitQual-unblind-susyFloat-LowerLimit-tb50-m0-500-m12-250-loose-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"output-files/fitQualNorm-unblind-susyFloat-LowerLimit-tb50-m0-500-m12-250-loose-sel.png") ;

      rfitSusyFloat.setAndFixSusySig(susyhigh) ;
      rfitSusyFloat.doFit() ;
      printf("\n\n ======== upper limit , susy floating==============\n\n") ;
      rfitSusyFloat.parameterSnapshot() ;
      rfitSusyFloat.fitQualityPlot(0,"output-files/fitQual-unblind-susyFloat-UpperLimit-tb50-m0-500-m12-250-loose-sel.png") ;
      rfitSusyFloat.fitQualityPlot(1,"output-files/fitQualNorm-unblind-susyFloat-UpperLimit-tb50-m0-500-m12-250-loose-sel.png") ;


   //-------------------------

      ra2bRoostatsClass3 rfitSusyFix(1,0,1) ;
      rfitSusyFix.initialize("input-files/rsinput-v4-SSVHPT-nj3-nb1-ht350-sb150_200-sig200.txt") ;
      rfitSusyFix.setSusyScanPoint("input-files/susyScanInput4-tanb50-SSVHPT-nj3-nb1-ht350-sb150_200-sig200.txt",500,250) ;
      rfitSusyFix.setAndFixSusySig(117.) ;

      rfitSusyFix.doFit() ;
      printf("\n\n ======== best fit, susy fixed ==============\n\n") ;
      rfitSusyFix.parameterSnapshot() ;
      rfitSusyFix.fitQualityPlot(0,"output-files/fitQual-unblind-BestFit-susyFixed-tb50-m0-500-m12-250-loose-sel.png") ;
      rfitSusyFix.fitQualityPlot(1,"output-files/fitQualNorm-unblind-BestFit-susyFixed-tb50-m0-500-m12-250-loose-sel.png") ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySig(117.) ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profilettwjSig(low,high,1,"output-files/ttwjSigProf-unblind-susyFixed-tb50-m0-500-m12-250-loose-sel.png", 120. ) ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySig(117.) ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profileqcdSig(low,high,1,"output-files/qcdSigProf-unblind-susyFixed-tb50-m0-500-m12-250-loose-sel.png", 120. ) ;
      rfitSusyFix.reinitialize() ;
      rfitSusyFix.setAndFixSusySig(117.) ;
      rfitSusyFix.doFit() ;
      rfitSusyFix.profileZnnSig(low,high,1,"output-files/znnSigProf-unblind-susyFixed-tb50-m0-500-m12-250-loose-sel.png", 120. ) ;




   }
