


   void run_ws_cls_hybrid2n() {

      gROOT->LoadMacro("workspace-macros-v1/ws_cls_hybrid2n.c+") ;

      TFile* outf = new TFile("cls-newfit-test-1000.root","recreate") ;

      gDirectory->pwd() ;

      bool isBgonlyStudy ;
      double poiVal ;
      int nToys ;
      bool makeTtree ;
      int verbLevel = 1 ;
      bool oneSidedTestStat = false ;

      makeTtree = true ;

      nToys = 1000 ;

    //-------

      for ( int poii=0; poii<1; poii++ ) {

         //poiVal = 15. + 15*poii ;
         poiVal = 0 ;

         isBgonlyStudy = false ;
         ws_cls_hybrid2n( "output-files/ws-newfit-lm9-1BL.root", isBgonlyStudy, poiVal, nToys, makeTtree, verbLevel, oneSidedTestStat ) ;

         isBgonlyStudy = true ;
         ws_cls_hybrid2n( "output-files/ws-newfit-lm9-1BL.root", isBgonlyStudy, poiVal, nToys, makeTtree, verbLevel, oneSidedTestStat ) ;

      }



      outf->ls() ;


      outf->Close() ;

   }




