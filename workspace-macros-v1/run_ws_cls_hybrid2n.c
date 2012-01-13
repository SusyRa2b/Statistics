


   void run_ws_cls_hybrid2n() {

      gROOT->LoadMacro("workspace-macros-v1/ws_cls_hybrid2n.c+") ;

      TFile* outf = new TFile("cls-newfit-Loose-1000toys-per-point.root","recreate") ;

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

         poiVal = 4. + 4*poii ;

         isBgonlyStudy = false ;
         ws_cls_hybrid2n( "output-files/ws-newfit-lm9-1BL.root", isBgonlyStudy, poiVal, nToys, makeTtree, verbLevel ) ;

         isBgonlyStudy = true ;
         ws_cls_hybrid2n( "output-files/ws-newfit-lm9-1BL.root", isBgonlyStudy, poiVal, nToys, makeTtree, verbLevel ) ;

      }



      outf->ls() ;


      outf->Close() ;

   }




