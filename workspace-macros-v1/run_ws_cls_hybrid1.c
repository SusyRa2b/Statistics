


   void run_ws_cls_hybrid1() {

      gROOT->LoadMacro("workspace-macros-v1/ws_cls_hybrid1.c+") ;

      TFile* outf = new TFile("cls-expected-2BL-sbvars.root","recreate") ;

      gDirectory->pwd() ;

      bool isBgonlyStudy ;
      double poiVal ;
      int nToys ;
      bool makeTtree ;
      int verbLevel = 0 ;

      makeTtree = true ;

      nToys = 1000 ;

    //-------

      for ( int poii=0; poii<8; poii++ ) {

         poiVal = 25. + 25*poii ;

         isBgonlyStudy = false ;
         ws_cls_hybrid1( "output-files/expected-ws-lm9-2BL-sbvars.root", isBgonlyStudy, poiVal, nToys, makeTtree, verbLevel ) ;

         isBgonlyStudy = true ;
         ws_cls_hybrid1( "output-files/expected-ws-lm9-2BL-sbvars.root", isBgonlyStudy, poiVal, nToys, makeTtree, verbLevel ) ;

      }



      outf->ls() ;


      outf->Close() ;

   }




