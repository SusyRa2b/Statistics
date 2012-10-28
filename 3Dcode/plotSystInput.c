

#include "TSystem.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>



    bool plotSystInput( const char* infile,
                        double target_mgl = 0., double target_mlsp = 0.,
                        double scaleMax = 0.2 ) {


      int nBinsMET(4), nBinsHT(4), nBinsBtag(3) ;

      gStyle->SetOptStat(0) ;
      gStyle->SetPadRightMargin(0.2) ;

      gDirectory->Delete("h_syst*") ;

      char command[1000] ;
      sprintf( command, "ls %s >& /dev/null", infile ) ;
      int returnstat = gSystem->Exec( command ) ;
      if ( returnstat != 0 ) {
         printf("\n\n ***  Input file doesn't exist: %s\n\n", infile ) ;
         return false ;
      }


      sprintf( command, "head -1 %s | awk '{print NF}'", infile ) ;
      const char* nfields_str = gSystem->GetFromPipe( command ) ;
      int nfields ;
      sscanf( nfields_str, "%d", &nfields ) ;
      printf(" Nfields is %d\n", nfields ) ;

      bool hasSL(false) ;

      int ArraySize ;
      if ( nfields == 2+2*(nBinsMET*nBinsHT*nBinsBtag) ) {
         hasSL = false ;
         printf("\n\n Format is consistent with no SL observables.\n") ;
         ArraySize = nfields ;
      } else if ( nfields == 2+3*(nBinsMET*nBinsHT*nBinsBtag) ) {
         hasSL = true ;
         printf("\n\n Format is consistent with including SL observables.\n") ;
         ArraySize = nfields ;
      } else {
         printf("\n\n I don't know what to do with nfields = %d\n\n", nfields ) ;
         return false ;
      }


      ifstream infq ;
      infq.open(infile) ;
      if ( !infq.good() ) {
         printf("\n\n *** setupShapeSyst: Problem opening input file: %s.\n\n", infile ) ;
         return false ;
      }

      bool found = false ;

      double syst_zl[10][10][10] ;
      double syst_sl[10][10][10] ;
      double syst_ldp[10][10][10] ;

      double minSyst(0.) ;
      double maxSyst(0.) ;

      while ( infq.good() ) {

         int nBins = nBinsMET*nBinsHT*nBinsBtag ;

         double ArrayContent[ArraySize] ;
         for ( int i=0; infq && i<ArraySize; ++ i) {
            infq >> ArrayContent[i] ;
         }

         double mgl  = ArrayContent[0] ;
         double mlsp = ArrayContent[1] ;

         if ( !(fabs( mgl-target_mgl ) < 10. && fabs( mlsp - target_mlsp ) < 10 ) ) continue ;

         found = true ;

         printf("\n\n Found mgl=%.0f, mlsp=%.0f\n\n", mgl, mlsp ) ;

         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
               for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {
                  if ( hasSL ) {
                     syst_zl [mbi][hbi][bbi]  = ArrayContent[2 + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
                     syst_sl [mbi][hbi][bbi]  = ArrayContent[2 + nBins + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
                     syst_ldp[mbi][hbi][bbi]  = ArrayContent[2 + 2*nBins + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
                  } else {
                     syst_zl [mbi][hbi][bbi]  = ArrayContent[2 + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
                     syst_sl [mbi][hbi][bbi]  = 0. ;
                     syst_ldp[mbi][hbi][bbi]  = ArrayContent[2 + nBins + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
                  }
                  if ( syst_zl [mbi][hbi][bbi] > maxSyst ) maxSyst = syst_zl [mbi][hbi][bbi] ;
                  if ( syst_sl [mbi][hbi][bbi] > maxSyst ) maxSyst = syst_sl [mbi][hbi][bbi] ;
                  if ( syst_ldp[mbi][hbi][bbi] > maxSyst ) maxSyst = syst_ldp[mbi][hbi][bbi] ;
                  if ( syst_zl [mbi][hbi][bbi] < minSyst ) minSyst = syst_zl [mbi][hbi][bbi] ;
                  if ( syst_sl [mbi][hbi][bbi] < minSyst ) minSyst = syst_sl [mbi][hbi][bbi] ;
                  if ( syst_ldp[mbi][hbi][bbi] < minSyst ) minSyst = syst_ldp[mbi][hbi][bbi] ;
               } // bbi.
            } // hbi.
         } // mbi.

         break ;

      } // reading file?

      if ( !found ) {
         printf("\n\n *** Did not find target point: mgl=%.0f, mlsp=%.0f\n\n", target_mgl, target_mlsp ) ;
         return false ;
      }

      printf("\n\n Min syst = %6.2f, Max syst = %6.2f\n\n", minSyst, maxSyst ) ;

      TH2F* h_syst_zl[10] ;
      TH2F* h_syst_sl[10] ;
      TH2F* h_syst_ldp[10] ;

      for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {

         char hname[100] ;
         char htitle[100] ;
         char binlabel[100] ;

         sprintf( hname, "h_syst_zl_%db", bbi+1 ) ;
         sprintf( htitle, "ZL syst, MET vs HT, nB=%d", bbi+1 ) ;
         h_syst_zl[bbi] = new TH2F( hname, htitle, nBinsHT, 0.5, 0.5+nBinsHT,   nBinsMET, 0.5, nBinsMET+0.5 ) ;
         sprintf( hname, "h_syst_sl_%db", bbi+1 ) ;
         sprintf( htitle, "SL syst, MET vs HT, nB=%d", bbi+1 ) ;
         h_syst_sl[bbi] = new TH2F( hname, htitle, nBinsHT, 0.5, 0.5+nBinsHT,   nBinsMET, 0.5, nBinsMET+0.5 ) ;
         sprintf( hname, "h_syst_ldp_%db", bbi+1 ) ;
         sprintf( htitle, "LDP syst, MET vs HT, nB=%d", bbi+1 ) ;
         h_syst_ldp[bbi] = new TH2F( hname, htitle, nBinsHT, 0.5, 0.5+nBinsHT,   nBinsMET, 0.5, nBinsMET+0.5 ) ;

         if ( scaleMax > 0 ) {
            h_syst_zl [bbi] -> SetMinimum( -1.*scaleMax ) ;
            h_syst_sl [bbi] -> SetMinimum( -1.*scaleMax ) ;
            h_syst_ldp[bbi] -> SetMinimum( -1.*scaleMax ) ;
            h_syst_zl [bbi] -> SetMaximum(     scaleMax ) ;
            h_syst_sl [bbi] -> SetMaximum(     scaleMax ) ;
            h_syst_ldp[bbi] -> SetMaximum(     scaleMax ) ;
         }

         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            sprintf( binlabel, "HT%d", hbi+1 ) ;
            h_syst_zl [bbi] -> GetXaxis() -> SetBinLabel( hbi+1, binlabel ) ;
            h_syst_sl [bbi] -> GetXaxis() -> SetBinLabel( hbi+1, binlabel ) ;
            h_syst_ldp[bbi] -> GetXaxis() -> SetBinLabel( hbi+1, binlabel ) ;
         } // hbi
         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            sprintf( binlabel, "MET%d", mbi+1 ) ;
            h_syst_zl [bbi] -> GetYaxis() -> SetBinLabel( nBinsMET-mbi, binlabel ) ;
            h_syst_sl [bbi] -> GetYaxis() -> SetBinLabel( nBinsMET-mbi, binlabel ) ;
            h_syst_ldp[bbi] -> GetYaxis() -> SetBinLabel( nBinsMET-mbi, binlabel ) ;
         } // hbi

         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
               h_syst_zl [bbi] -> SetBinContent( hbi+1, nBinsMET-mbi, syst_zl [mbi][hbi][bbi] ) ;
               h_syst_sl [bbi] -> SetBinContent( hbi+1, nBinsMET-mbi, syst_sl [mbi][hbi][bbi] ) ;
               h_syst_ldp[bbi] -> SetBinContent( hbi+1, nBinsMET-mbi, syst_ldp[mbi][hbi][bbi] ) ;
            } // hbi.
         } // mbi.

      } // bbi.


      TCanvas* csyst = (TCanvas*) gDirectory->FindObject("csyst") ;
      if ( csyst==0x0 ) {
         csyst = new TCanvas( "csyst", "Systematics", 900, 900 ) ;
      }

      csyst->Clear() ;
      csyst->Divide(3,3) ;

      csyst->cd(1) ;
      h_syst_zl[0] -> Draw("colz") ;

      csyst->cd(2) ;
      h_syst_zl[1] -> Draw("colz") ;

      csyst->cd(3) ;
      h_syst_zl[2] -> Draw("colz") ;


      csyst->cd(4) ;
      h_syst_sl[0] -> Draw("colz") ;

      csyst->cd(5) ;
      h_syst_sl[1] -> Draw("colz") ;

      csyst->cd(6) ;
      h_syst_sl[2] -> Draw("colz") ;


      csyst->cd(7) ;
      h_syst_ldp[0] -> Draw("colz") ;

      csyst->cd(8) ;
      h_syst_ldp[1] -> Draw("colz") ;

      csyst->cd(9) ;
      h_syst_ldp[2] -> Draw("colz") ;



      return true ;


   } // plotSystInput




