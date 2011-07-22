

   void MCSMonly_limitcomp_loose() {

      gROOT->LoadMacro("ra2bRoostatsClass3.c+") ;

      float low,high ;
      float susylow, susyhigh ;

      // ra2bRoostatsClass3( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;







   //-------------------------

      ra2bRoostatsClass3 rfitSusyFloatNC(1,0,1) ;
      rfitSusyFloatNC.initialize("input-files/mc-inputs-SMonly-loose-sel.txt") ;
      rfitSusyFloatNC.doFit() ;
      rfitSusyFloatNC.parameterSnapshot() ;
      printf("\n\n ======== best fit , susy floating, no contam ==============\n\n") ;
      rfitSusyFloatNC.fitQualityPlot(0,"output-files/fitQual-MCSMonly-BestFit-susyFloat-tb50-m0-500-m12-250-loose-sel-nocontam.png") ;
      rfitSusyFloatNC.fitQualityPlot(1,"output-files/fitQualNorm-MCSMonly-BestFit-susyFloat-tb50-m0-500-m12-250-loose-sel-nocontam.png") ;
      rfitSusyFloatNC.profileSusySig(susylow,susyhigh,1,"output-files/susySigProf-MCSMonly-susyFloat-tb50-m0-500-m12-250-loose-sel-nocontam.png", 120. ) ;

      rfitSusyFloatNC.setAndFixSusySig(susyhigh) ;
      rfitSusyFloatNC.doFit() ;
      printf("\n\n ======== susy upper limit , no contam ==============\n\n") ;
      rfitSusyFloatNC.parameterSnapshot() ;
      rfitSusyFloatNC.fitQualityPlot(0,"output-files/fitQual-MCSMonly-UpperLimit-tb50-m0-500-m12-250-loose-sel-nocontam.png") ;
      rfitSusyFloatNC.fitQualityPlot(1,"output-files/fitQualNorm-MCSMonly-UpperLimit-tb50-m0-500-m12-250-loose-sel-nocontam.png") ;





      ra2bRoostatsClass3 rfitSusyFloatWC(1,0,1) ;
      rfitSusyFloatWC.initialize("input-files/mc-inputs-SMonly-loose-sel.txt") ;
      rfitSusyFloatWC.setSusyScanPoint("input-files/susyScanInput4-tanb50-SSVHPT-nj3-nb1-ht350-sb150_200-sig200.txt",500,250) ;
      rfitSusyFloatWC.doFit() ;
      rfitSusyFloatWC.parameterSnapshot() ;
      printf("\n\n ======== best fit , susy floating, no contam ==============\n\n") ;
      rfitSusyFloatWC.fitQualityPlot(0,"output-files/fitQual-MCSMonly-BestFit-susyFloat-tb50-m0-500-m12-250-loose-sel-withcontam.png") ;
      rfitSusyFloatWC.fitQualityPlot(1,"output-files/fitQualNorm-MCSMonly-BestFit-susyFloat-tb50-m0-500-m12-250-loose-sel-withcontam.png") ;
      rfitSusyFloatWC.profileSusySig(susylow,susyhigh,1,"output-files/susySigProf-MCSMonly-susyFloat-tb50-m0-500-m12-250-loose-sel-withcontam.png", 120. ) ;

      rfitSusyFloatWC.setAndFixSusySig(susyhigh) ;
      rfitSusyFloatWC.doFit() ;
      printf("\n\n ======== susy upper limit , no contam ==============\n\n") ;
      rfitSusyFloatWC.parameterSnapshot() ;
      rfitSusyFloatWC.fitQualityPlot(0,"output-files/fitQual-MCSMonly-UpperLimit-tb50-m0-500-m12-250-loose-sel-withcontam.png") ;
      rfitSusyFloatWC.fitQualityPlot(1,"output-files/fitQualNorm-MCSMonly-UpperLimit-tb50-m0-500-m12-250-loose-sel-withcontam.png") ;


   }
