
#include "TFile.h"
#include "TChain.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TStyle.h"


   const int  ncomps(4) ;
   char compname[ncomps][100] = { "susy", "ttwj", "qcd", "znn" } ;
   int  compcolor[ncomps] = { 6, kBlue-9, 2, kGreen-3 } ;

// const int nBinsMET(4) ;
// const int nBinsHT(4) ;

   const int nBinsMET(3) ;
   const int nBinsHT(3) ;


   const int nBinsBjets(3) ;


   const int nCanvas(9) ;

   char plotname[nCanvas][100] = { "nbjet", "met", "ht", "met_nbjet1", "met_nbjet2", "met_nbjet3", "ht_nbjet1", "ht_nbjet2", "ht_nbjet3" } ;

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

   void toyplots1( const char* infilename = "output-toymc2b-mgl850-mlsp600-100evts-newMC-wsyst-test5a-useE0Ltrue/toy-results.root" ) {

      gStyle->SetOptStat(0) ;

      gDirectory->Delete("h*") ;

      TChain* toytt = new TChain("toytt") ;
      toytt->Add( infilename ) ;

      toytt->Print("toponly") ;


      for ( int cani=0; cani<nCanvas; cani++ ) {

         printf("\n\n Canvas %s\n\n", plotname[cani] ) ;

         char cname[100] ;
         sprintf(cname, "c_%s", plotname[cani] ) ;
         TCanvas* can = (TCanvas*) gDirectory->FindObject( cname ) ;
         if ( can == 0x0 ) {
            can = new TCanvas( cname, cname, 500, 500 ) ;
            can -> SetWindowPosition(1+30*cani,1+25*cani) ;
         }
         can->Clear() ;
         can->Divide(2,2) ;

         for ( int ci=0; ci<ncomps; ci++ ) {


            can -> cd(ci+1) ;

            char hname[1000] ;
            char htitle[1000] ;
            sprintf( hname, "hfit_%s_%s", compname[ci], plotname[cani] ) ;
            sprintf( htitle, "Fit %s vs %s", compname[ci], plotname[cani] ) ;
            TProfile* hfit   = new TProfile( hname, htitle, (svmax[cani]-svmin[cani]+1), svmin[cani]-0.5, svmax[cani]+0.5, "s" ) ;
            sprintf( hname, "hmcval_%s_%s", compname[ci], plotname[cani] ) ;
            sprintf( htitle, "MCval %s vs %s", compname[ci], plotname[cani] ) ;
            TProfile* hmcval   = new TProfile( hname, htitle, (svmax[cani]-svmin[cani]+1), svmin[cani]-0.5, svmax[cani]+0.5 ) ;

            for ( int svi=svmin[cani]; svi<=svmax[cani]; svi++ ) {

               char drawarg1[10000] ;

               sprintf( drawarg1, daformat[cani], "fit", compname[ci], svi-1, svi, "fit", compname[ci], plotname[cani] ) ;
               printf( "%s\n", drawarg1 ) ;
               toytt -> Draw( drawarg1 ) ;

               sprintf( drawarg1, daformat[cani], "mcval", compname[ci], svi-1, svi, "mcval", compname[ci], plotname[cani] ) ;
               printf( "%s\n", drawarg1 ) ;
               toytt -> Draw( drawarg1 ) ;

            } // svi.

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

            hfit -> SetMarkerStyle(20) ;
            hfit -> Draw() ;
            hmcval -> Draw("histsame") ;
            hfit -> Draw("same") ;

            printf("\n") ;

         } // ci.


      } // cani.




   }

