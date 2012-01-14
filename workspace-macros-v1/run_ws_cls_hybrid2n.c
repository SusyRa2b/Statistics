


   void run_ws_cls_hybrid2n() {

      gROOT->LoadMacro("workspace-macros-v1/ws_cls_hybrid2n.c+") ;

      TFile* outf = new TFile("cls-newfit-Loose-withLM9-xSec4.0-400toys-per-point.root","recreate") ;

      gDirectory->pwd() ;

      bool isBgonlyStudy ;
      double poiVal ;
      int nToys ;
      bool makeTtree ;
      int verbLevel = 0 ;

      makeTtree = true ;

      nToys = 400 ;

    //-------

      for ( int poii=0; poii<8; poii++ ) {

         poiVal = 15. + 15*poii ;

         isBgonlyStudy = false ;
         ws_cls_hybrid2n( "output-files/ws-newfit-lm9-1BL-withLM9-xSec4.0.root", isBgonlyStudy, poiVal, nToys, makeTtree, verbLevel ) ;

         isBgonlyStudy = true ;
         ws_cls_hybrid2n( "output-files/ws-newfit-lm9-1BL-withLM9-xSec4.0.root", isBgonlyStudy, poiVal, nToys, makeTtree, verbLevel ) ;

      }



      outf->ls() ;


      outf->Close() ;

   }




