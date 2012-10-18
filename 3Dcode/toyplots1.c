
#include "TFile.h"
#include "TChain.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TLine.h"


   const int  ncomps(4) ;
   char compname[ncomps][100] = { "susy", "ttwj", "qcd", "znn" } ;
   int  compcolor[ncomps] = { 6, kBlue-9, 2, kGreen-3 } ;

   const int nBinsMET(4) ;
   const int nBinsHT(4) ;

// const int nBinsMET(3) ;
// const int nBinsHT(3) ;


   const int nBinsBjets(3) ;


   const int nCanvas(9) ;

   char plotname[nCanvas][100] = { 
       "nbjet",
       "met",
       "ht",
       "met_nbjet1",
       "met_nbjet2",
       "met_nbjet3",
       "ht_nbjet1",
       "ht_nbjet2",
       "ht_nbjet3"
       } ;

   char daformat[nCanvas][10000] = {
        "Sum$(%s_%s_0lep_3da[][][%d]):%d>>+h%s_%s_%s"
       ,"Sum$(%s_%s_0lep_3da[%d][][]):%d>>+h%s_%s_%s"
       ,"Sum$(%s_%s_0lep_3da[][%d][]):%d>>+h%s_%s_%s"
       ,"Sum$(%s_%s_0lep_3da[%d][][0]):%d>>+h%s_%s_%s"
       ,"Sum$(%s_%s_0lep_3da[%d][][1]):%d>>+h%s_%s_%s"
       ,"Sum$(%s_%s_0lep_3da[%d][][2]):%d>>+h%s_%s_%s"
       ,"Sum$(%s_%s_0lep_3da[][%d][0]):%d>>+h%s_%s_%s"
       ,"Sum$(%s_%s_0lep_3da[][%d][1]):%d>>+h%s_%s_%s"
       ,"Sum$(%s_%s_0lep_3da[][%d][2]):%d>>+h%s_%s_%s"
        } ;

   int  svmin[nCanvas] = { 1, 1, 1, 1, 1, 1, 1, 1, 1 } ;
   int  svmax[nCanvas] = { nBinsBjets, nBinsMET, nBinsHT, nBinsMET, nBinsMET, nBinsMET, nBinsHT, nBinsHT, nBinsHT } ;

  //----------------------------

   void toyplots1( const char* infilename = "toy1200300-wSUSY.root", bool doDiff=false, bool doFrac=false, float sethmax=0.0 ) {

      gStyle->SetOptStat(0) ;
      gStyle->SetLabelSize(0.08,"y") ;
      gStyle->SetNdivisions(504,"y") ;
      gStyle->SetPadLeftMargin(0.15) ;
      gStyle->SetEndErrorSize(5) ;

      gDirectory->Delete("h*") ;

      TChain* toytt = new TChain("toytt") ;
      toytt->Add( infilename ) ;

      toytt->Print("toponly") ;

      TLine* line = new TLine() ;
      line -> SetLineStyle(2) ;



      for ( int cani=0; cani<nCanvas; cani++ ) {

         printf("\n\n Canvas %s\n\n", plotname[cani] ) ;

         char cname[100] ;
         sprintf(cname, "c_%s", plotname[cani] ) ;
         TCanvas* can = (TCanvas*) gDirectory->FindObject( cname ) ;
         if ( can == 0x0 ) {
            can = new TCanvas( cname, cname, 370, 370 ) ;
            can -> SetWindowPosition(1+30*cani,1+25*cani) ;
         }
         can->Clear() ;
         can->Divide(2,2) ;

         for ( int ci=0; ci<ncomps; ci++ ) {


            can -> cd(ci+1) ;

            char hname[1000] ;
            char htitle[1000] ;
            sprintf( htitle, "Fit %s vs %s", compname[ci], plotname[cani] ) ;
            sprintf( hname, "hfit_%s_%s", compname[ci], plotname[cani] ) ;
            TProfile* hfit   = new TProfile( hname, htitle, (svmax[cani]-svmin[cani]+1), svmin[cani]-0.5, svmax[cani]+0.5, "s" ) ;
            sprintf( hname, "hfitem_%s_%s", compname[ci], plotname[cani] ) ;
            TProfile* hfitem = new TProfile( hname, htitle, (svmax[cani]-svmin[cani]+1), svmin[cani]-0.5, svmax[cani]+0.5 ) ;
            sprintf( hname, "hmcval_%s_%s", compname[ci], plotname[cani] ) ;
            sprintf( htitle, "MCval %s vs %s", compname[ci], plotname[cani] ) ;
            TProfile* hmcval   = new TProfile( hname, htitle, (svmax[cani]-svmin[cani]+1), svmin[cani]-0.5, svmax[cani]+0.5 ) ;

            for ( int svi=svmin[cani]; svi<=svmax[cani]; svi++ ) {

               char drawarg1[10000] ;

               sprintf( drawarg1, daformat[cani], "fit", compname[ci], svi-1, svi, "fit", compname[ci], plotname[cani] ) ;
               printf( "%s\n", drawarg1 ) ;
               toytt -> Draw( drawarg1, "fit_covqual_susyfloat==3" ) ;

               sprintf( drawarg1, daformat[cani], "fit", compname[ci], svi-1, svi, "fitem", compname[ci], plotname[cani] ) ;
               printf( "%s\n", drawarg1 ) ;
               toytt -> Draw( drawarg1, "fit_covqual_susyfloat==3" ) ;

               sprintf( drawarg1, daformat[cani], "mcval", compname[ci], svi-1, svi, "mcval", compname[ci], plotname[cani] ) ;
               printf( "%s\n", drawarg1 ) ;
               toytt -> Draw( drawarg1, "fit_covqual_susyfloat==3" ) ;

            } // svi.

            if ( !doDiff ) {

               double hfitmax = 1.1 * (hfit->GetBinContent( hfit->GetMaximumBin() ) + hfit->GetBinError( hfit->GetMaximumBin() ) ) ;
               double hmcvalmax = 1.1 * (hmcval->GetBinContent( hmcval->GetMaximumBin() ) ) ;
               if ( hfitmax > hmcvalmax ) {
                  hfit->SetMaximum( hfitmax ) ;
               } else {
                  hfit->SetMaximum( hmcvalmax ) ;
               }
               hmcval -> SetLineColor( compcolor[ci] ) ;
               hmcval -> SetLineWidth(2) ;
               hfit   -> SetLineWidth(2) ;
               hfit   -> SetLineColor(4) ;
               hfitem -> SetLineWidth(2) ;

               hfit -> SetMarkerStyle(20) ;
               hfit -> Draw("e1") ;
               hmcval -> Draw("histsame") ;
               hfit -> Draw("same") ;
               hfitem->Draw("samee1") ;

            } else {

               sprintf( hname, "hdiff_%s_%s", compname[ci], plotname[cani] ) ;
               sprintf( htitle, "Ave fit-true diff %s vs %s", compname[ci], plotname[cani] ) ;
               TH1F* hdiff = new TH1F( hname, htitle, (svmax[cani]-svmin[cani]+1), svmin[cani]-0.5, svmax[cani]+0.5 ) ;

               sprintf( hname, "hdiffem_%s_%s", compname[ci], plotname[cani] ) ;
               sprintf( htitle, "Ave fit-true diff %s vs %s", compname[ci], plotname[cani] ) ;
               TH1F* hdiffem = new TH1F( hname, htitle, (svmax[cani]-svmin[cani]+1), svmin[cani]-0.5, svmax[cani]+0.5 ) ;


               double drawmax(0.) ;
               double drawmin(0.) ;
               for ( int bi=1; bi<=(svmax[cani]-svmin[cani]+1); bi++ ) {
                  double diff = hfit->GetBinContent( bi ) - hmcval->GetBinContent( bi ) ;
                  double err  = hfit->GetBinError( bi ) ;
                  double errm = hfitem->GetBinError( bi ) ;
                  if ( doFrac ) {
                     if ( hmcval->GetBinContent( bi ) > 0. ) {
                        diff = diff / ( hmcval->GetBinContent( bi ) ) ;
                        err  = err  / ( hmcval->GetBinContent( bi ) ) ;
                        errm = errm / ( hmcval->GetBinContent( bi ) ) ;
                     }
                     printf(" fit=%7.1f, MC=%7.1f, frac diff = %5.3f +/- %5.3f\n", hfit->GetBinContent( bi ), hmcval->GetBinContent( bi ), diff, err ) ;
                  }
                  hdiff -> SetBinContent( bi, diff ) ;
                  hdiff -> SetBinError( bi, err ) ;
                  hdiffem -> SetBinContent( bi, diff ) ;
                  hdiffem -> SetBinError( bi, errm ) ;
                  if ( 1.1*(diff+err) > drawmax ) { drawmax = 1.1*(diff+err) ; }
                  if ( diff<0 && 1.1*(diff-err) < drawmin ) { drawmin = 1.1*(diff-err) ; }
               } // bi.

               if ( sethmax <= 0 ) {
                  if ( fabs(drawmax) > fabs(drawmin) ) {
                     hdiff->SetMaximum(drawmax) ;
                     hdiff->SetMinimum(-drawmax) ;
                  } else {
                     hdiff->SetMaximum(fabs(drawmin)) ;
                     hdiff->SetMinimum(-fabs(drawmin)) ;
                  }
               } else {
                  hdiff->SetMaximum(sethmax) ;
                  hdiff->SetMinimum(-sethmax) ;
               }

               hdiff -> SetMarkerStyle(20) ;
               hdiff -> SetLineWidth(2) ;
               hdiff -> SetLineColor( compcolor[ci] ) ;
               hdiff -> Draw("e1") ;
               hdiffem -> Draw("samee1") ;
               line->DrawLine( svmin[cani]-0.5, 0.,    svmax[cani]+0.5, 0. ) ;

            }

            printf("\n") ;

         } // ci.


      } // cani.

   } // toyplots1

  //=========================================================================================================

   void positionWindows1() {

      Int_t x,y ;
      UInt_t w,h ;
      gVirtualX->GetGeometry(-1,x,y,w,h) ;
      printf("\n\n X11 screen size: w=%d, h=%d\n\n", w,h ) ;

      int wmargin(40), hmargin(100) ;

      TCanvas* c ;
      char cname[1000] ;
      int width, height ;

      sprintf( cname, "c_nbjet" ) ;
      c = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( c == 0x0 ) { printf("\n\n Can't find canvas %s\n\n", cname ) ; return ; }
      width  = c -> GetWindowWidth() ;
      height = c -> GetWindowHeight() ;
      c->SetWindowPosition( w-width-wmargin, hmargin ) ;

      sprintf( cname, "c_met" ) ;
      c = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( c == 0x0 ) { printf("\n\n Can't find canvas %s\n\n", cname ) ; return ; }
      width  = c -> GetWindowWidth() ;
      height = c -> GetWindowHeight() ;
      c->SetWindowPosition( w-2*width-wmargin, hmargin ) ;

      sprintf( cname, "c_ht" ) ;
      c = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( c == 0x0 ) { printf("\n\n Can't find canvas %s\n\n", cname ) ; return ; }
      width  = c -> GetWindowWidth() ;
      height = c -> GetWindowHeight() ;
      c->SetWindowPosition( w-3*width-wmargin, hmargin ) ;


      sprintf( cname, "c_met_nbjet3" ) ;
      c = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( c == 0x0 ) { printf("\n\n Can't find canvas %s\n\n", cname ) ; return ; }
      width  = c -> GetWindowWidth() ;
      height = c -> GetWindowHeight() ;
      c->SetWindowPosition( w-width-wmargin, h-1.05*height-hmargin ) ;

      sprintf( cname, "c_met_nbjet2" ) ;
      c = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( c == 0x0 ) { printf("\n\n Can't find canvas %s\n\n", cname ) ; return ; }
      width  = c -> GetWindowWidth() ;
      height = c -> GetWindowHeight() ;
      c->SetWindowPosition( w-2*width-wmargin, h-1.05*height-hmargin ) ;

      sprintf( cname, "c_met_nbjet1" ) ;
      c = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( c == 0x0 ) { printf("\n\n Can't find canvas %s\n\n", cname ) ; return ; }
      width  = c -> GetWindowWidth() ;
      height = c -> GetWindowHeight() ;
      c->SetWindowPosition( w-3*width-wmargin, h-1.05*height-hmargin ) ;



      sprintf( cname, "c_ht_nbjet3" ) ;
      c = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( c == 0x0 ) { printf("\n\n Can't find canvas %s\n\n", cname ) ; return ; }
      width  = c -> GetWindowWidth() ;
      height = c -> GetWindowHeight() ;
      c->SetWindowPosition( w-width-wmargin, h-2.1*height-hmargin ) ;

      sprintf( cname, "c_ht_nbjet2" ) ;
      c = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( c == 0x0 ) { printf("\n\n Can't find canvas %s\n\n", cname ) ; return ; }
      width  = c -> GetWindowWidth() ;
      height = c -> GetWindowHeight() ;
      c->SetWindowPosition( w-2*width-wmargin, h-2.1*height-hmargin ) ;

      sprintf( cname, "c_ht_nbjet1" ) ;
      c = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( c == 0x0 ) { printf("\n\n Can't find canvas %s\n\n", cname ) ; return ; }
      width  = c -> GetWindowWidth() ;
      height = c -> GetWindowHeight() ;
      c->SetWindowPosition( w-3*width-wmargin, h-2.1*height-hmargin ) ;







   } // positionWindows1.

  //=========================================================================================================



