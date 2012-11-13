

#include "TSystem.h"
#include "TFile.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>



    bool plotSigInput( const char* infile,
                        double target_mgl = 0., double target_mlsp = 0.
                        ) {


      gStyle->SetLabelSize(0.06,"x") ;
      gStyle->SetLabelSize(0.06,"y") ;
      gStyle->SetPadLeftMargin(0.11) ;

      double dataLumi(12.0) ; // 1/fb.

      TFile prospino("referenceXSecs.root");
      TH1F *gluinoxsec8TeV = (TH1F*) prospino.Get("gluino8TeV_NLONLL") ;
      if ( gluinoxsec8TeV == 0 ) {
         printf("\n\n *** Can't find histogram with name gluino8TeV_NLONLL in referenceXSecs.root.\n\n") ;
         return false ;
      }

      int theBin8TeV = gluinoxsec8TeV->FindBin( target_mgl ) ;
      if ( theBin8TeV <=0 || theBin8TeV > gluinoxsec8TeV->GetNbinsX() ) {
         printf("\n\n *** can't find bin for mgl=%.0f.  Returned %d\n\n", target_mgl, theBin8TeV ) ;
         return false ;
      }
      double xsec8TeV = gluinoxsec8TeV->GetBinContent( theBin8TeV ) ;

      prospino.Close() ;

         // each point has 10k events generated. The sigma is in pb and I want to normalized to 15 fb-1. 
         // so multiple cross section by 1.5 to get events in 15 fb-1

      double scaleFactor = (dataLumi/10.)*xsec8TeV ;
      printf("\n Will rescale SUSY MC counts by %g\n\n", scaleFactor ) ;

      gStyle->SetPaintTextFormat(".2f") ;

      int nBinsMET(4), nBinsHT(4), nBinsBtag(3) ;

      gStyle->SetOptStat(0) ;
      gStyle->SetPadRightMargin(0.2) ;

      gDirectory->Delete("h_*") ;

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


      int ArraySize = 3 + 4*(nBinsMET*nBinsHT*nBinsBtag) ;
      if ( nfields == ArraySize ) {
         printf("\n\n Found expected array size: %d\n\n", nfields ) ;
      } else {
         printf("\n\n I don't know what to do with nfields = %d.  Expected %d.\n\n", nfields, ArraySize ) ;
         return false ;
      }


      ifstream infq ;
      infq.open(infile) ;
      if ( !infq.good() ) {
         printf("\n\n *** setupShapeSyst: Problem opening input file: %s.\n\n", infile ) ;
         return false ;
      }

      bool found = false ;

      double n_0l_raw[10][10][10] ;
      double n_1l_raw[10][10][10] ;
      double n_ldp_raw[10][10][10] ;

      double n_0l_err[10][10][10] ;
      double n_1l_err[10][10][10] ;
      double n_ldp_err[10][10][10] ;


      double binMax(0.) ;

      while ( infq.good() ) {

         //int nBins = nBinsMET*nBinsHT*nBinsBtag ;

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

                  n_0l_raw[mbi][hbi][bbi]  = ArrayContent[3 + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
                  n_1l_raw[mbi][hbi][bbi]  = 0. ;
                  n_ldp_raw[mbi][hbi][bbi] = ArrayContent[3 + (nBinsMET*nBinsHT*nBinsBtag) + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;

                  if ( scaleFactor * n_0l_raw[mbi][hbi][bbi] > binMax ) { binMax = scaleFactor * n_0l_raw[mbi][hbi][bbi] ; }
                  if ( scaleFactor * n_ldp_raw[mbi][hbi][bbi] > binMax ) { binMax = scaleFactor * n_ldp_raw[mbi][hbi][bbi] ; }

                  n_0l_err[mbi][hbi][bbi]  = ArrayContent[3 + 2*(nBinsMET*nBinsHT*nBinsBtag) + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
                  n_1l_err[mbi][hbi][bbi]  = 0.1 ;
                  n_ldp_err[mbi][hbi][bbi] = ArrayContent[3 + 3*(nBinsMET*nBinsHT*nBinsBtag) + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;

               } // bbi.
            } // hbi.
         } // mbi.

         break ;

      } // reading file?

      if ( !found ) {
         printf("\n\n *** Did not find target point: mgl=%.0f, mlsp=%.0f\n\n", target_mgl, target_mlsp ) ;
         return false ;
      }


      TH2F* h_sig_zl[10] ;
      TH2F* h_sig_sl[10] ;
      TH2F* h_sig_ldp[10] ;

      TH2F* h_staterr_zl[10] ;
      TH2F* h_staterr_sl[10] ;
      TH2F* h_staterr_ldp[10] ;


      for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {

         printf("\n\n\n ======= nB = %d,  Event counts :\n\n", bbi+1 ) ;

         char hname[100] ;
         char htitle[100] ;
         char binlabel[100] ;

         sprintf( hname, "h_sig_zl_%db", bbi+1 ) ;
         sprintf( htitle, "ZL signal, MET vs HT, nB=%d", bbi+1 ) ;
         h_sig_zl[bbi] = new TH2F( hname, htitle, nBinsHT, 0.5, 0.5+nBinsHT,   nBinsMET, 0.5, nBinsMET+0.5 ) ;
         sprintf( hname, "h_sig_sl_%db", bbi+1 ) ;
         sprintf( htitle, "SL signal, MET vs HT, nB=%d", bbi+1 ) ;
         h_sig_sl[bbi] = new TH2F( hname, htitle, nBinsHT, 0.5, 0.5+nBinsHT,   nBinsMET, 0.5, nBinsMET+0.5 ) ;
         sprintf( hname, "h_sig_ldp_%db", bbi+1 ) ;
         sprintf( htitle, "LDP signal, MET vs HT, nB=%d", bbi+1 ) ;
         h_sig_ldp[bbi] = new TH2F( hname, htitle, nBinsHT, 0.5, 0.5+nBinsHT,   nBinsMET, 0.5, nBinsMET+0.5 ) ;


         h_sig_zl[bbi]  -> SetFillColor(6) ;
         h_sig_sl[bbi]  -> SetFillColor(6) ;
         h_sig_ldp[bbi] -> SetFillColor(6) ;

         h_sig_zl[bbi]  -> SetMarkerSize(2.5) ;
         h_sig_sl[bbi]  -> SetMarkerSize(2.5) ;
         h_sig_ldp[bbi] -> SetMarkerSize(2.5) ;

         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            sprintf( binlabel, "HT%d", hbi+1 ) ;
            h_sig_zl [bbi] -> GetXaxis() -> SetBinLabel( hbi+1, binlabel ) ;
            h_sig_sl [bbi] -> GetXaxis() -> SetBinLabel( hbi+1, binlabel ) ;
            h_sig_ldp[bbi] -> GetXaxis() -> SetBinLabel( hbi+1, binlabel ) ;
         } // hbi
         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            sprintf( binlabel, "MET%d", mbi+1 ) ;
            h_sig_zl [bbi] -> GetYaxis() -> SetBinLabel( nBinsMET-mbi, binlabel ) ;
            h_sig_sl [bbi] -> GetYaxis() -> SetBinLabel( nBinsMET-mbi, binlabel ) ;
            h_sig_ldp[bbi] -> GetYaxis() -> SetBinLabel( nBinsMET-mbi, binlabel ) ;
         } // hbi

         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
               printf("  m,h,b (%d,%d,%d) : zl_raw =%5.1f, zl_weighted =%6.2f;   sl_raw =%5.1f, sl_weighted =%6.2f;   ldp_raw = %5.1f, ldp_weighted = %6.2f\n",
                    mbi+1, hbi+1, bbi+1,
                    n_0l_raw [mbi][hbi][bbi], scaleFactor * n_0l_raw [mbi][hbi][bbi],
                    n_1l_raw [mbi][hbi][bbi], scaleFactor * n_1l_raw [mbi][hbi][bbi],
                    n_ldp_raw [mbi][hbi][bbi], scaleFactor * n_ldp_raw [mbi][hbi][bbi]
                    ) ;
               h_sig_zl [bbi] -> SetBinContent( hbi+1, nBinsMET-mbi, scaleFactor * n_0l_raw [mbi][hbi][bbi] ) ;
               h_sig_sl [bbi] -> SetBinContent( hbi+1, nBinsMET-mbi, scaleFactor * n_1l_raw [mbi][hbi][bbi] ) ;
               h_sig_ldp[bbi] -> SetBinContent( hbi+1, nBinsMET-mbi, scaleFactor * n_ldp_raw[mbi][hbi][bbi] ) ;
            } // hbi.
         } // mbi.

         h_sig_zl[bbi]  -> SetMaximum(binMax) ;
         h_sig_sl[bbi]  -> SetMaximum(binMax) ;
         h_sig_ldp[bbi] -> SetMaximum(binMax) ;


      //---------

         printf("\n\n\n ======= nB = %d,  Stat errors :\n\n", bbi+1 ) ;

         sprintf( hname, "h_staterr_zl_%db", bbi+1 ) ;
         sprintf( htitle, "ZL signal stat err, MET vs HT, nB=%d", bbi+1 ) ;
         h_staterr_zl[bbi] = new TH2F( hname, htitle, nBinsHT, 0.5, 0.5+nBinsHT,   nBinsMET, 0.5, nBinsMET+0.5 ) ;
         sprintf( hname, "h_staterr_sl_%db", bbi+1 ) ;
         sprintf( htitle, "SL signal stat err, MET vs HT, nB=%d", bbi+1 ) ;
         h_staterr_sl[bbi] = new TH2F( hname, htitle, nBinsHT, 0.5, 0.5+nBinsHT,   nBinsMET, 0.5, nBinsMET+0.5 ) ;
         sprintf( hname, "h_staterr_ldp_%db", bbi+1 ) ;
         sprintf( htitle, "LDP signal stat err, MET vs HT, nB=%d", bbi+1 ) ;
         h_staterr_ldp[bbi] = new TH2F( hname, htitle, nBinsHT, 0.5, 0.5+nBinsHT,   nBinsMET, 0.5, nBinsMET+0.5 ) ;

         h_staterr_zl[bbi]  -> SetFillColor(18) ;
         h_staterr_sl[bbi]  -> SetFillColor(18) ;
         h_staterr_ldp[bbi] -> SetFillColor(18) ;

         h_staterr_zl[bbi]  -> SetMarkerSize(2.5) ;
         h_staterr_sl[bbi]  -> SetMarkerSize(2.5) ;
         h_staterr_ldp[bbi] -> SetMarkerSize(2.5) ;

         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            sprintf( binlabel, "HT%d", hbi+1 ) ;
            h_staterr_zl [bbi] -> GetXaxis() -> SetBinLabel( hbi+1, binlabel ) ;
            h_staterr_sl [bbi] -> GetXaxis() -> SetBinLabel( hbi+1, binlabel ) ;
            h_staterr_ldp[bbi] -> GetXaxis() -> SetBinLabel( hbi+1, binlabel ) ;
         } // hbi
         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            sprintf( binlabel, "MET%d", mbi+1 ) ;
            h_staterr_zl [bbi] -> GetYaxis() -> SetBinLabel( nBinsMET-mbi, binlabel ) ;
            h_staterr_sl [bbi] -> GetYaxis() -> SetBinLabel( nBinsMET-mbi, binlabel ) ;
            h_staterr_ldp[bbi] -> GetYaxis() -> SetBinLabel( nBinsMET-mbi, binlabel ) ;
         } // hbi

         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            for ( int hbi=0; hbi<nBinsHT; hbi++ ) {

               double err_0l(0.), err_1l(0.), err_ldp(0.) ;

               if ( n_0l_raw [mbi][hbi][bbi]  > 0. ) { err_0l  = n_0l_err[mbi][hbi][bbi] / n_0l_raw [mbi][hbi][bbi]  ; }
               if ( n_1l_raw [mbi][hbi][bbi]  > 0. ) { err_1l  = n_1l_err[mbi][hbi][bbi] / n_1l_raw [mbi][hbi][bbi]  ; }
               if ( n_ldp_raw [mbi][hbi][bbi] > 0. ) { err_ldp = n_ldp_err[mbi][hbi][bbi] / n_ldp_raw [mbi][hbi][bbi]  ; }

               printf("  m,h,b (%d,%d,%d) : zl_raw = %5.1f, zl=%.2f;   sl_raw = %5.1f, sl=%.2f;    ldp_raw = %5.1f  ldp=%.2f\n",
                   mbi+1, hbi+1, bbi+1,
                   n_0l_err[mbi][hbi][bbi], err_0l,
                   n_1l_err[mbi][hbi][bbi], err_1l,
                   n_ldp_err[mbi][hbi][bbi], err_ldp
                   ) ;

               h_staterr_zl [bbi] -> SetBinContent( hbi+1, nBinsMET-mbi,  err_0l  ) ;
               h_staterr_sl [bbi] -> SetBinContent( hbi+1, nBinsMET-mbi,  err_1l  ) ;
               h_staterr_ldp[bbi] -> SetBinContent( hbi+1, nBinsMET-mbi,  err_ldp ) ;
            } // hbi.
         } // mbi.

         h_staterr_zl[bbi]  -> SetMaximum(1.0) ;
         h_staterr_sl[bbi]  -> SetMaximum(1.0) ;
         h_staterr_ldp[bbi] -> SetMaximum(1.0) ;


         printf("\n\n\n") ;


      } // bbi.






      TCanvas* csig = (TCanvas*) gDirectory->FindObject("csig") ;
      if ( csig==0x0 ) {
         csig = new TCanvas( "csig", "signal counts", 900, 900 ) ;
      }

      TCanvas* cstaterr = (TCanvas*) gDirectory->FindObject("cstaterr") ;
      if ( cstaterr==0x0 ) {
         cstaterr = new TCanvas( "cstaterr", "signal stat err", 900, 900 ) ;
      }


      char drawoptions[100] ;
      sprintf( drawoptions, "boxtext" ) ;
  //  sprintf( drawoptions, "colz" ) ;



    //------

      csig->Clear() ;
      csig->Divide(3,3) ;

      csig->cd(1) ;
      h_sig_zl[0] -> DrawCopy( drawoptions ) ;

      csig->cd(2) ;
      h_sig_zl[1] -> DrawCopy( drawoptions ) ;

      csig->cd(3) ;
      h_sig_zl[2] -> DrawCopy( drawoptions ) ;


      csig->cd(4) ;
      h_sig_sl[0] -> DrawCopy( drawoptions ) ;

      csig->cd(5) ;
      h_sig_sl[1] -> DrawCopy( drawoptions ) ;

      csig->cd(6) ;
      h_sig_sl[2] -> DrawCopy( drawoptions ) ;


      csig->cd(7) ;
      h_sig_ldp[0] -> DrawCopy( drawoptions ) ;

      csig->cd(8) ;
      h_sig_ldp[1] -> DrawCopy( drawoptions ) ;

      csig->cd(9) ;
      h_sig_ldp[2] -> DrawCopy( drawoptions ) ;


      csig->Update() ; csig->Draw() ;




    //------

      cstaterr->Clear() ;
      cstaterr->Divide(3,3) ;

      cstaterr->cd(1) ;
      h_staterr_zl[0] -> DrawCopy( drawoptions ) ;

      cstaterr->cd(2) ;
      h_staterr_zl[1] -> DrawCopy( drawoptions ) ;

      cstaterr->cd(3) ;
      h_staterr_zl[2] -> DrawCopy( drawoptions ) ;


      cstaterr->cd(4) ;
      h_staterr_sl[0] -> DrawCopy( drawoptions ) ;

      cstaterr->cd(5) ;
      h_staterr_sl[1] -> DrawCopy( drawoptions ) ;

      cstaterr->cd(6) ;
      h_staterr_sl[2] -> DrawCopy( drawoptions ) ;


      cstaterr->cd(7) ;
      h_staterr_ldp[0] -> DrawCopy( drawoptions ) ;

      cstaterr->cd(8) ;
      h_staterr_ldp[1] -> DrawCopy( drawoptions ) ;

      cstaterr->cd(9) ;
      h_staterr_ldp[2] -> DrawCopy( drawoptions ) ;


      cstaterr->Update() ; cstaterr->Draw() ;





      return true ;


   } // plotstaterrInput




