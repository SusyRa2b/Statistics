
#include "TChain.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPaveStats.h"


   const int  ncomps(4) ;
   char compname[ncomps][100] = { "susy", "ttwj", "qcd", "znn" } ;
   int  compcolor[ncomps] = { 6, kBlue-9, 2, kGreen-3 } ;

    //--- For the bnb, bmet, and bht variables, 0 means integrate over that dimension
    //    while >0 means set the bin to that index, where 1 is the lowest bin.

   void histtoy1(
                  int bmet = 0, int bht = 0, int bnb = 0,
                  const char* file1="batch-3x3v1-ntoy400/toy-results-realsyst-wttwjcorr-useEZL-mgl800-mlsp700-nsig100.root",
                  const char* file2="batch-3x3v1-ntoy400/toy-results-tinysyst-wttwjcorr-useEZL-mgl800-mlsp700-nsig100.root"
   ) {

      gStyle->SetOptStat("mr") ;
      gStyle->SetNdivisions(506,"x") ;
      gStyle->SetNdivisions(506,"y") ;
      gStyle->SetTitleH(0.08) ;


      TChain ch1("toytt") ;
      int nadded1 = ch1.Add( file1 ) ;
      printf(" n added for file 1 %s: %d\n", file1, nadded1 ) ;
      if ( nadded1 <= 0 ) { printf("\n\n *** No TTrees in %s.\n\n", file1 ) ; return ; }

      TChain ch2("toytt") ;
      int nadded2 = ch2.Add( file2, 0 ) ;
      bool ch2good(true) ;
      printf(" n added for file 2 %s: %d\n", file2, nadded2 ) ;
      if ( nadded2 <= 0 ) { printf("\n\n *** No TTrees in %s.\n\n", file2 ) ; ch2good = false ; }

      char cname[1000] ;
      sprintf( cname, "chist_met%d_ht%d_mb%d", bmet, bht, bnb ) ;

      char ctitle[1000] ;
      char metstr[100] ;
      char htstr[100] ;
      char nbstr[100] ;
      if ( bmet>0 ) { sprintf( metstr, "MET bin %d", bmet ) ; } else { sprintf( metstr, "Sum over MET bins" ) ; }
      if ( bht >0 ) { sprintf( htstr , "HT  bin %d", bht )  ; } else { sprintf( htstr , "Sum over HT bins" ) ; }
      if ( bnb >0 ) { sprintf( nbstr , "nb  bin %d", bnb )  ; } else { sprintf( nbstr , "Sum over nb bins" ) ; }
      sprintf( ctitle, "%s, %s, %s", metstr, htstr, nbstr ) ;
      TCanvas* can = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( can == 0x0 ) {
         can = new TCanvas( cname, ctitle, 900, 900 ) ;
      }

      can -> Clear() ;
      can -> Divide(2,2) ;

      char bmetstr[10] ;
      char bhtstr[10] ;
      char bnbstr[10] ;
      if ( bmet > 0 ) { sprintf( bmetstr, "%d", bmet-1 ) ; } else { sprintf( bmetstr, "" ) ; }
      if ( bht  > 0 ) { sprintf( bhtstr , "%d", bht -1 ) ; } else { sprintf( bhtstr , "" ) ; }
      if ( bnb  > 0 ) { sprintf( bnbstr , "%d", bnb -1 ) ; } else { sprintf( bnbstr , "" ) ; }


      int nhbins(30) ;

      for ( int ci=0; ci<ncomps; ci++ ) {

         printf("\n") ;

         can -> cd(ci+1) ;

         char drawarg[1000] ;

         sprintf( drawarg, "Sum$(fit_%s_0lep_3da[%s][%s][%s])>>h1%s(%d)", compname[ci], bmetstr, bhtstr, bnbstr, compname[ci], nhbins ) ;
         printf( "%s\n", drawarg ) ;
         ch1.Draw( drawarg, "" ) ;

         sprintf( drawarg, "Sum$(fit_%s_0lep_3da[%s][%s][%s])>>h2%s(%d)", compname[ci], bmetstr, bhtstr, bnbstr, compname[ci], nhbins ) ;
         printf( "%s\n", drawarg ) ;
         if ( ch2good ) { ch2.Draw( drawarg, "", "same" ) ; }

         sprintf( drawarg, "Sum$(mcval_%s_0lep_3da[%s][%s][%s])>>hmc%s(%d)", compname[ci], bmetstr, bhtstr, bnbstr, compname[ci], 10*nhbins ) ;
         printf( "%s\n", drawarg ) ;
         ch1.Draw( drawarg, "" ) ;

         char hname[1000] ;

         sprintf( hname, "h1%s", compname[ci] ) ;
         TH1F* h1 = (TH1F*) gDirectory->FindObject( hname ) ;
         if ( h1 == 0x0 ) { printf("\n\n *** can't find histogram %s\n\n", hname ) ; return ; }

         TH1F* h2(0x0) ;
         if ( ch2good ) {
            sprintf( hname, "h2%s", compname[ci] ) ;
            h2 = (TH1F*) gDirectory->FindObject( hname ) ;
            if ( h2 == 0x0 ) { printf("\n\n *** can't find histogram %s\n\n", hname ) ; return ; }
         }

         sprintf( hname, "hmc%s", compname[ci] ) ;
         TH1F* hmc = (TH1F*) gDirectory->FindObject( hname ) ;
         if ( hmc == 0x0 ) { printf("\n\n *** can't find histogram %s\n\n", hname ) ; return ; }

         float h1max = 1.1*(h1->GetMaximum()) ;
         float h2max = 0. ;
         if ( ch2good ) { h2max = 1.1*(h2->GetMaximum()) ; }
         float hmax ;
         if ( h1max > h2max ) { hmax = h1max ; } else { hmax = h2max ; }

         h1->SetLineColor(1) ;
         hmc->SetLineColor(4) ;
         h1->SetLineWidth(2) ;
         hmc->SetLineWidth(2) ;
         if ( ch2good ) {
            h2->SetLineColor(2) ;
            h2->SetLineWidth(2) ;
         }

         char htitle[1000] ;
         if ( bmet>0 ) { sprintf( metstr, "MET bin %d", bmet ) ; } else { sprintf( metstr, " " ) ; }
         if ( bht >0 ) { sprintf( htstr , "HT  bin %d", bht )  ; } else { sprintf( htstr , " " ) ; }
         if ( bnb >0 ) { sprintf( nbstr , "nb  bin %d", bnb )  ; } else { sprintf( nbstr , " " ) ; }
         sprintf( ctitle, "%s %s %s", metstr, htstr, nbstr ) ;
         sprintf( htitle, "%s : %s", compname[ci], ctitle ) ;
         h1->SetTitle( htitle ) ;
         if ( ch2good ) { h2->SetTitle( htitle ) ; }

         h1->SetMaximum( hmax ) ;

         h1->Draw() ;
         if ( ch2good ) { h2->Draw("sames") ; }
         hmc->Draw("same") ;

         gPad->Update() ;
         TPaveStats* st1 = (TPaveStats*) h1->FindObject("stats") ;
         TPaveStats* st2(0x0) ;
         if ( ch2good ) { st2 = (TPaveStats*) h2->FindObject("stats") ; }
         st1->SetX1NDC(0.75) ;
         st1->SetY1NDC(0.8) ;
         st1->SetX2NDC(0.99) ;
         st1->SetY2NDC(0.95) ;
         st1->SetTextColor(1) ;
         if ( ch2good ) {
            st2->SetX1NDC(0.75) ;
            st2->SetY1NDC(0.6) ;
            st2->SetX2NDC(0.99) ;
            st2->SetY2NDC(0.75) ;
            st2->SetTextColor(2) ;
         }

         h1->DrawCopy() ;
         if ( ch2good ) { h2->DrawCopy("sames") ; }
         hmc->DrawCopy("same") ;



      } // ci.



   } // histtoy1.
