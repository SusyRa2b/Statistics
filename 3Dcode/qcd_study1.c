
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TFile.h"
#include "TLine.h"
#include "TLegend.h"

#include "updateFileValue.c"


#include <iostream>

   void saveHist(const char* filename, const char* pat) ;
   void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;
   TH1F* bookHist(const char* hname, const char* htitle, int nBinsMET, int nBinsHT, int sampleIndex ) ;

   //-------

   void qcd_study1( bool fillHists=false, bool savePlots=false, bool interactive=true, const char* datfile = "null" ) {

      TLine* line = new TLine() ;

      gStyle->SetPalette(1) ;
      gStyle->SetPadBottomMargin(0.20) ;

      gDirectory->Delete("h*") ;

      int nBinsBjets = 3 ;

      const int nBinsMET = 3 ;
      const int nBinsHT  = 3 ;
      float Mbins[nBinsMET+1] = { 125, 200,  350, 99999. } ;
      float Hbins[nBinsHT+1]  = { 400, 600, 1000, 99999. } ;

      TH2F* hdummy = new TH2F("hdummy","",2, Mbins[0], 1500., 2, Hbins[0], 1500. ) ;

      const int nQcdSamples(9) ;

      char inputfile[9][1000] = {
        "files15fb_8TeV/nn-qcd-0120-to-0170.root"
       ,"files15fb_8TeV/nn-qcd-0170-to-0300.root"
       ,"files15fb_8TeV/nn-qcd-0300-to-0470.root"
       ,"files15fb_8TeV/nn-qcd-0470-to-0600.root"
       ,"files15fb_8TeV/nn-qcd-0600-to-0800.root"
       ,"files15fb_8TeV/nn-qcd-0800-to-1000.root"
       ,"files15fb_8TeV/nn-qcd-1000-to-1400.root"
       ,"files15fb_8TeV/nn-qcd-1400-to-1800.root"
       ,"files15fb_8TeV/nn-qcd-1800.root"
      } ;

      char samplename[9][100] = {
        "qcd_0120_to_0170"
       ,"qcd_0170_to_0300"
       ,"qcd_0300_to_0470"
       ,"qcd_0470_to_0600"
       ,"qcd_0600_to_0800"
       ,"qcd_0800_to_1000"
       ,"qcd_1000_to_1400"
       ,"qcd_1400_to_1800"
       ,"qcd_1800_to_9999"
      } ;

 //   int samplecolor[9] = {
 //      kOrange+8,
 //      kRed,
 //      KMagenta+2,
 //      KBlue+1,
 //      KAzure+2,
 //      KGreen+2,
 //      kPink-2,
 //      kViolet+2,
 //      kOrange+2
 //   } ;

      int samplecolor[9] = {
         800+8,
         632,
         616+2,
         600+1,
         860+2,
         416+2,
         900-2,
         880+2,
         800+2
      } ;

      TH2F*   h0lep[nQcdSamples][nBinsBjets] ;
      TH2F*   hldp [nQcdSamples][nBinsBjets] ;

      TCanvas* cqcd = (TCanvas*) gDirectory->FindObject("cqcd") ;
      if ( cqcd == 0x0 ) {
         cqcd = new TCanvas("cqcd", "qcd study", 700, 900 ) ;
      }

      char hname[1000] ;
      char htitle[1000] ;

      if ( fillHists ) {

         TChain* qcdch[nQcdSamples] ;

         printf("\n\n") ;
         for ( int si=0; si<nQcdSamples; si++ ) {

            qcdch[si] = new TChain("tree") ;
            printf(" %2d : connecting to %s\n", si, inputfile[si] ) ;
            qcdch[si] -> Add( inputfile[si] ) ;

            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               sprintf( hname, "h_0lep_%db_%s", bbi+1, samplename[si] ) ;
               sprintf( htitle, "QCD 0lep yield, nb=%d, %s", bbi+1, samplename[si] ) ;
               printf("         booking hist %s : %s\n", hname, htitle ) ;
               h0lep[si][bbi] = new TH2F( hname, htitle, nBinsMET, Mbins, nBinsHT, Hbins ) ;
               h0lep[si][bbi] -> Sumw2() ;
               sprintf( hname, "h_ldp_%db_%s", bbi+1, samplename[si] ) ;
               sprintf( htitle, "QCD  LDP yield, nb=%d, %s", bbi+1, samplename[si] ) ;
               printf("         booking hist %s  : %s\n", hname, htitle ) ;
               hldp [si][bbi] = new TH2F( hname, htitle, nBinsMET, Mbins, nBinsHT, Hbins ) ;
               hldp [si][bbi] -> Sumw2() ;

            } // bbi.

         } // si.
         printf("\n\n") ;





         char basecuts[10000] ;
         sprintf( basecuts, "pt_1st_leadJet>50&&pt_2nd_leadJet>50&&pt_3rd_leadJet>50&&HT>400&&MET>125&&nMu==0&&nEl==0&&nB>=1&&nJets>=3" ) ;

         char bcut[3][100] = { "nB==1", "nB==2", "nB>=3" } ;


         for ( int si=0; si<nQcdSamples; si++ ) {

            printf(" %2d : %s : 0lep\n", si, samplename[si] ) ; cout << flush ;
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               char arg1[1000] ;

               char cuts0lep[10000] ;
               sprintf( cuts0lep, "(%s)&&(minDelPhiN>4&&%s)", basecuts, bcut[bbi] ) ;
               printf("     %db, 0lep cuts : %s\n", bbi+1, cuts0lep ) ;
               sprintf( arg1, "HT:MET>>h_0lep_%db_%s", bbi+1, samplename[si] ) ;
               qcdch[si] -> Draw( arg1, cuts0lep ) ;
               hdummy->Draw() ;
               h0lep[si][bbi]->Draw("samecolz") ;
               cqcd->Update() ; cqcd->Draw() ;


               char cutsldp[10000] ;
               sprintf( cutsldp, "(%s)&&(minDelPhiN<=4&&%s)", basecuts, bcut[bbi] ) ;
               printf("     %db, ldp  cuts : %s\n", bbi+1, cutsldp  ) ;
               sprintf( arg1, "HT:MET>>h_ldp_%db_%s", bbi+1, samplename[si] ) ;
               qcdch[si] -> Draw( arg1, cutsldp, "colz" ) ;
               hdummy->Draw() ;
               hldp[si][bbi]->Draw("samecolz") ;
               cqcd->Update() ; cqcd->Draw() ;

            } // bbi.

         } // si.

         saveHist( "rootfiles/qcd-study1.root", "h*" ) ;

      }

      //--------

      gStyle->SetOptStat(0) ;

      gDirectory->Delete("h*") ;
      loadHist( "rootfiles/qcd-study1.root" ) ;

      for ( int si=0; si<nQcdSamples; si++ ) {
         for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

            sprintf( hname, "h_0lep_%db_%s", bbi+1, samplename[si] ) ;
            printf("  loading %s\n", hname ) ;
            h0lep[si][bbi] = (TH2F*) gDirectory->FindObject( hname ) ;
            if ( h0lep[si][bbi] == 0x0 ) { printf("\n\n *** %s missing.\n\n", hname ) ; return ; }
            sprintf( hname, "h_ldp_%db_%s", bbi+1, samplename[si] ) ;
            printf("  loading %s\n", hname ) ;
            hldp [si][bbi] = (TH2F*) gDirectory->FindObject( hname ) ;
            if ( hldp [si][bbi] == 0x0 ) { printf("\n\n *** %s missing.\n\n", hname ) ; return ; }

         } // bbi.
      } // si.






     //--- Flatten 2D histos and compute 0lep/LDP ratios for each sample.

      int nbins = nBinsMET*(nBinsHT+1) + 1 ;

      TH1F* hflat_0lep[nQcdSamples][nBinsBjets] ;
      TH1F* hflat_ldp [nQcdSamples][nBinsBjets] ;
      TH1F* hflat_0lepldp_ratio[nQcdSamples][nBinsBjets] ;

      sprintf( hname, "hflat_0lep_all" ) ;
      sprintf( htitle, "All QCD 0lep events" ) ;
      TH1F* hflat_0lep_all = bookHist( hname, htitle, nBinsMET, nBinsHT, 5 ) ;
      sprintf( hname, "hflat_ldp_all" ) ;
      sprintf( htitle, "All QCD ldp events" ) ;
      TH1F* hflat_ldp_all = bookHist( hname, htitle, nBinsMET, nBinsHT, 5 ) ;


      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
         for ( int si=0; si<nQcdSamples; si++ ) {

            sprintf( hname, "hflat_0lep_%db_%s", bbi+1, samplename[si] ) ;
            sprintf( htitle, "QCD 0lep events, nb=%d, %s", bbi+1, samplename[si] ) ;
            hflat_0lep[si][bbi] = bookHist( hname, htitle, nBinsMET, nBinsHT, si ) ;
            hflat_0lep[si][bbi]->SetFillColor(11) ;
            hldp[si][bbi]->Print("all") ;

            sprintf( hname, "hflat_ldp_%db_%s", bbi+1, samplename[si] ) ;
            sprintf( htitle, "QCD LDP events, nb=%d, %s", bbi+1, samplename[si] ) ;
            hflat_ldp[si][bbi] = bookHist( hname, htitle, nBinsMET, nBinsHT, si ) ;
            hflat_ldp[si][bbi]->SetFillColor(11) ;

            sprintf( hname, "hflat_0lepldp_ratio_%db_%s", bbi+1, samplename[si] ) ;
            sprintf( htitle, "QCD 0lep/LDP ratio, nb=%d, %s", bbi+1, samplename[si] ) ;
            hflat_0lepldp_ratio[si][bbi] = bookHist( hname, htitle, nBinsMET, nBinsHT, si ) ;

            for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
               for ( int hbi=0; hbi<nBinsHT; hbi++ ) {

                  int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

                  float n0lep     = h0lep[si][bbi] -> GetBinContent( mbi+1, hbi+1 ) ;
                  float n0lep_err = h0lep[si][bbi] -> GetBinError( mbi+1, hbi+1 ) ;

                  float nldp      = hldp [si][bbi] -> GetBinContent( mbi+1, hbi+1 ) ;
                  float nldp_err  = hldp [si][bbi] -> GetBinError( mbi+1, hbi+1 ) ;

                  hflat_0lep[si][bbi] -> SetBinContent( histbin, n0lep ) ;
                  hflat_0lep[si][bbi] -> SetBinError( histbin, n0lep_err ) ;

                  hflat_ldp[si][bbi]  -> SetBinContent( histbin, nldp ) ;
                  hflat_ldp[si][bbi]  -> SetBinError( histbin, nldp_err ) ;

                  double ratio(0.) ;
                  if ( nldp > 0 ) { ratio = n0lep / nldp ; }
                  double ratio_err(0.) ;
                  if ( n0lep_err > 0. && nldp_err > 0. ) {
                     ratio_err = ratio * sqrt( pow(n0lep_err/n0lep,2) + pow(nldp_err/nldp,2) ) ;
                  }

                  hflat_0lepldp_ratio[si][bbi] -> SetBinContent( histbin, ratio ) ;
                  hflat_0lepldp_ratio[si][bbi] -> SetBinError( histbin, ratio_err ) ;


               } // hbi.
            } // mbi.

            hflat_0lep_all -> Add( hflat_0lep[si][bbi] ) ;
            hflat_ldp_all  -> Add( hflat_ldp [si][bbi] ) ;

            cqcd->Clear() ;
            cqcd->Divide(1,3) ;

            cqcd->cd(1) ;
            hflat_0lep[si][bbi] -> DrawCopy("histe") ;
            hflat_0lep[si][bbi] -> DrawCopy("samee") ;

            cqcd->cd(2) ;
            hflat_ldp[si][bbi] -> DrawCopy("histe") ;
            hflat_ldp[si][bbi] -> DrawCopy("samee") ;

            cqcd->cd(3) ;
            hflat_0lepldp_ratio[si][bbi]->SetMinimum(-0.1) ;
            hflat_0lepldp_ratio[si][bbi]->SetMaximum(0.6) ;
            hflat_0lepldp_ratio[si][bbi]->SetMarkerStyle(20) ;
            hflat_0lepldp_ratio[si][bbi]->SetLineWidth(2) ;
            hflat_0lepldp_ratio[si][bbi]->DrawCopy() ;
            gPad->SetGridx(1) ;
            gPad->SetGridy(1) ;
            line->DrawLine(0.5,0,nbins+0.5,0) ;


            cqcd->Update() ; cqcd->Draw() ;

         } // si.
      } // bbi.

      TH1F* hflat_0lep_all_sqrtNerrs = (TH1F*) hflat_0lep_all->Clone("hflat_0lep_all_sqrtNerrs") ;
      TH1F* hflat_ldp_all_sqrtNerrs  = (TH1F*) hflat_ldp_all ->Clone("hflat_ldp_all_sqrtNerrs") ;

      for ( int bi=1; bi<=nbins; bi++ ) {
         double val ;
         val = hflat_0lep_all_sqrtNerrs->GetBinContent( bi ) ;
         hflat_0lep_all_sqrtNerrs -> SetBinError( bi, sqrt(val) ) ;
         val = hflat_ldp_all_sqrtNerrs->GetBinContent( bi ) ;
         hflat_ldp_all_sqrtNerrs -> SetBinError( bi, sqrt(val) ) ;
      }












     //--- Compute sample average of 0lep/LDP ratio


      TH1F* hflat_0lepldp_ratio_ave[nBinsBjets] ;

      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

         sprintf( hname, "hflat_0lepldp_ratio_ave_%db", bbi+1 ) ;
         sprintf( htitle, "QCD 0lep/LDP average ratio, nb=%d", bbi+1 ) ;
         hflat_0lepldp_ratio_ave[bbi] = bookHist( hname, htitle, nBinsMET, nBinsHT, 5 ) ;

         hflat_0lepldp_ratio_ave[bbi] -> SetMarkerStyle(24) ;
         hflat_0lepldp_ratio_ave[bbi] -> SetMarkerSize(2.0) ;
         hflat_0lepldp_ratio_ave[bbi] -> SetLineWidth(3) ;

         for ( int binind=1; binind<=nbins; binind++ ) {

            double wvsum(0.) ;
            double wsum(0.) ;

            for ( int si=0; si<nQcdSamples; si++ ) {

               double val = hflat_0lepldp_ratio[si][bbi] -> GetBinContent( binind ) ;
               double err = hflat_0lepldp_ratio[si][bbi] -> GetBinError( binind ) ;
               if ( err > 0 ) {
                  double weight = 1./(err*err) ;
                  wvsum += val*weight ;
                  wsum  += weight ;
               }
            } // si.

            if ( wsum > 0. ) {
               hflat_0lepldp_ratio_ave[bbi] -> SetBinContent( binind, wvsum/wsum ) ;
               hflat_0lepldp_ratio_ave[bbi] -> SetBinError( binind, 1./sqrt(wsum) ) ;
            }

         } // binind

      } // bbi.








     //--- Compute RMS spread of samples and total uncertainty


      TH1F* hflat_0lepldp_ratio_ave_withRMSerror[nBinsBjets] ;

      printf("\n\n") ;
      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

         sprintf( hname, "hflat_0lepldp_ratio_ave_%db_withRMSerror", bbi+1 ) ;
         sprintf( htitle, "QCD 0lep/LDP average ratio, RMS included in error, nb=%d", bbi+1 ) ;
         hflat_0lepldp_ratio_ave_withRMSerror[bbi] = bookHist( hname, htitle, nBinsMET, nBinsHT, 5 ) ;

         hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetMarkerStyle(24) ;
         hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetMarkerSize(2.0) ;
         hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetLineWidth(3) ;

         for ( int binind=1; binind<=nbins; binind++ ) {

            double sumsquare(0.) ;
            int    nsum(0) ;

            double ave     = hflat_0lepldp_ratio_ave[bbi] -> GetBinContent( binind ) ;
            double ave_err = hflat_0lepldp_ratio_ave[bbi] -> GetBinError( binind ) ;

            for ( int si=0; si<nQcdSamples; si++ ) {

               double val = hflat_0lepldp_ratio[si][bbi] -> GetBinContent( binind ) ;
               double err = hflat_0lepldp_ratio[si][bbi] -> GetBinError( binind ) ;

               if ( err <= 0 ) { continue ; }
               if ( val <= 0 ) { continue ; }
               if ( err/val > 0.5 ) { continue ; } //-- do not include points with very large errors.


               double diff = val - ave ;

               sumsquare += diff*diff ;

               nsum ++ ;

            } // si.

            if ( nsum > 1 ) {
               double rms = sqrt(sumsquare/nsum) ;
               double totalerr = sqrt( ave_err*ave_err + rms*rms ) ;
               const char* binlabel = hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> GetXaxis() -> GetBinLabel( binind ) ;
               printf(" %s : %s : ratio = %5.3f, stat err = %5.3f, RMS = %5.3f,  total err = %5.3f\n",
                   hname, binlabel, ave, ave_err, rms, totalerr ) ;
               hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetBinContent( binind, ave ) ;
               hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetBinError( binind, totalerr ) ;
            } else {
               hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetBinContent( binind, ave ) ;
               hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetBinError( binind, ave_err ) ;
            }

         } // binind
         printf("\n") ;

      } // bbi.
      printf("\n\n") ;






     //--- compute global average 0lep/LDP ratio.

      double all0lep(0.) ;
      double allldp(0.) ;

      printf("\n\n") ;
      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
         for ( int si=0; si<nQcdSamples; si++ ) {

            double sample_all0lep = h0lep[si][bbi] -> Integral() ;
            double sample_allldp  = hldp [si][bbi] -> Integral() ;
            double sample_ratio(0.) ;
            if ( sample_allldp > 0 ) { sample_ratio = sample_all0lep / sample_allldp ; }

            printf( "%s, nb=%d : 0lep = %8.1f,  LDP = %8.1f, ratio = %5.3f\n", samplename[si], bbi+1, sample_all0lep, sample_allldp, sample_ratio ) ;

            all0lep += sample_all0lep ;
            allldp += sample_allldp ;

         } // si.
         printf("\n") ;
      } // bbi.
      printf("\n\n") ;

      double global_0lepldp_ratio(0.) ;
      if ( allldp > 0.) {
         global_0lepldp_ratio = all0lep / allldp ;
      } else {
         printf("\n\n *** You screwed up.\n\n") ; return ;
      }
      printf("  overall : 0lep = %9.1f,  LDP = %9.1f, ratio = %5.3f\n", all0lep, allldp, global_0lepldp_ratio ) ;










     //--- Compute scale factors.  SF_i = (0lep/LDP)_i / global_(0lep/LDP)

      TH1F* hflat_scale_factor[nBinsBjets] ;
      TH1F* hflat_scale_factor_withRMSerror[nBinsBjets] ;

      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

         sprintf( hname, "hflat_scale_factor_%db", bbi+1 ) ;
         sprintf( htitle, "QCD scale factor, nb=%d", bbi+1 ) ;
         hflat_scale_factor[bbi] = (TH1F*) hflat_0lepldp_ratio_ave[bbi]->Clone( hname ) ;
         hflat_scale_factor[bbi] -> SetTitle( htitle ) ;
         hflat_scale_factor[bbi] -> Scale( 1./ global_0lepldp_ratio ) ;

         sprintf( hname, "hflat_scale_factor_%db_withRMSerror", bbi+1 ) ;
         sprintf( htitle, "QCD scale factor, nb=%d, RMS error included ", bbi+1 ) ;
         hflat_scale_factor_withRMSerror[bbi] = (TH1F*) hflat_0lepldp_ratio_ave_withRMSerror[bbi]->Clone( hname ) ;
         hflat_scale_factor_withRMSerror[bbi] -> SetTitle( htitle ) ;
         hflat_scale_factor_withRMSerror[bbi] -> Scale( 1./ global_0lepldp_ratio ) ;

      } // bbi.









     //--- set Scale Factor values in dat file if one is provided.

      if ( strcmp( datfile, "null" ) != 0 ) {

         printf("\n\n\n Setting QCD scale factors in %s\n\n", datfile ) ;

         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            for ( int hbi=0; hbi<nBinsMET; hbi++ ) {
               for ( int bbi=0; bbi<nBinsMET; bbi++ ) {
                  char parname[1000] ;
                  int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
                  double val = hflat_scale_factor_withRMSerror[bbi]->GetBinContent( histbin ) ;
                  double err = hflat_scale_factor_withRMSerror[bbi]->GetBinError( histbin ) ;
                  if ( err <= 0 ) {
                     val = 1.0 ;
                     err = 3.0 ;
                  }
                  sprintf( parname, "sf_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  updateFileValue( datfile, parname, val ) ;
                  sprintf( parname, "sf_qcd_M%d_H%d_%db_err", mbi+1, hbi+1, bbi+1 ) ;
                  updateFileValue( datfile, parname, err ) ;
               } // bbi.
            } // hbi.
         } // mbi.

      }







     //===========  End calculations.  Plots and other output below here.  ===================================


     //--- Plot all samples together

      TCanvas* cqcd2 = (TCanvas*) gDirectory->FindObject("cqcd2") ;
      if ( cqcd2 == 0x0 ) {
         cqcd2 = new TCanvas("cqcd2", "qcd study", 1700, 600 ) ;
      }

      TLegend* l2 = new TLegend( 0.79, 0.70,  0.99, 0.99 ) ;
      l2->AddEntry( hflat_0lepldp_ratio[0][0], "120 to 170" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[1][0], "170 to 300" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[2][0], "300 to 470" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[3][0], "470 to 600" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[4][0], "600 to 800" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[5][0], "800 to 1000" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[6][0], "1000 to 1400" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[7][0], "1400 to 1800" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[8][0], "> 1800" ) ;
      l2->AddEntry( hflat_0lepldp_ratio_ave[0], "Average ratio" ) ;


      cqcd2 -> Clear() ;
      cqcd2 -> Divide(3,1) ;

      cqcd2 -> cd(1) ;
      hflat_0lepldp_ratio[0][0] -> SetTitle("QCD 0lep/LDP ratio, nb=1") ;
      for ( int si=0; si<nQcdSamples; si++ ) {
         hflat_0lepldp_ratio[si][0] -> SetMarkerStyle(20+si) ;
         hflat_0lepldp_ratio[si][0] -> SetLineColor(samplecolor[si]) ;
         hflat_0lepldp_ratio[si][0] -> SetMarkerColor(samplecolor[si]) ;
         if ( si == 0 ) { hflat_0lepldp_ratio[si][0] -> DrawCopy() ; } else { hflat_0lepldp_ratio[si][0] -> DrawCopy("same") ; }
      }
      hflat_0lepldp_ratio_ave[0]->DrawCopy("same") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      l2->Draw() ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd2 -> cd(2) ;
      hflat_0lepldp_ratio[0][1] -> SetTitle("QCD 0lep/LDP ratio, nb=2") ;
      for ( int si=0; si<nQcdSamples; si++ ) {
         hflat_0lepldp_ratio[si][1] -> SetMarkerStyle(20+si) ;
         hflat_0lepldp_ratio[si][1] -> SetLineColor(samplecolor[si]) ;
         hflat_0lepldp_ratio[si][1] -> SetMarkerColor(samplecolor[si]) ;
         if ( si == 0 ) { hflat_0lepldp_ratio[si][1] -> DrawCopy() ; } else { hflat_0lepldp_ratio[si][1] -> DrawCopy("same") ; }
      }
      hflat_0lepldp_ratio_ave[1]->DrawCopy("same") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      l2->Draw() ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd2 -> cd(3) ;
      hflat_0lepldp_ratio[0][2] -> SetTitle("QCD 0lep/LDP ratio, nb>=3") ;
      for ( int si=0; si<nQcdSamples; si++ ) {
         hflat_0lepldp_ratio[si][2] -> SetMarkerStyle(20+si) ;
         hflat_0lepldp_ratio[si][2] -> SetLineColor(samplecolor[si]) ;
         hflat_0lepldp_ratio[si][2] -> SetMarkerColor(samplecolor[si]) ;
         if ( si == 0 ) { hflat_0lepldp_ratio[si][2] -> DrawCopy() ; } else { hflat_0lepldp_ratio[si][2] -> DrawCopy("same") ; }
      }
      hflat_0lepldp_ratio_ave[2]->DrawCopy("same") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      l2->Draw() ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;


      cqcd2->Update() ; cqcd2->Draw() ;

      if ( savePlots ) {
         cqcd2 -> SaveAs("outputfiles/qcd-study-allonone.png") ;
         cqcd2 -> SaveAs("outputfiles/qcd-study-allonone.pdf") ;
      }








     //--- Go through the samples, 1 by 1, and draw numerator, denominator, and ratio

      TCanvas* cqcd3 = (TCanvas*) gDirectory->FindObject("cqcd3") ;
      if ( cqcd3 == 0x0 ) {
         cqcd3 = new TCanvas("cqcd3", "qcd study", 1200, 900 ) ;
      }

      char oldtitle[1000] ;
      sprintf( oldtitle, "QCD 0lep/LDP ratio, nb=1, %s", samplename[0] ) ;
      hflat_0lepldp_ratio[0][0] -> SetTitle(oldtitle) ;
      sprintf( oldtitle, "QCD 0lep/LDP ratio, nb=2, %s", samplename[0] ) ;
      hflat_0lepldp_ratio[0][1] -> SetTitle(oldtitle) ;
      sprintf( oldtitle, "QCD 0lep/LDP ratio, nb>=3, %s", samplename[0] ) ;
      hflat_0lepldp_ratio[0][2] -> SetTitle(oldtitle) ;

      cqcd3->Clear() ;
      cqcd3->Divide(3,3) ;

      for ( int si=0; si<nQcdSamples; si++ ) {
         TLegend* l3 = new TLegend( 0.65, 0.78, 0.99, 0.92) ;
         l3->AddEntry( hflat_0lepldp_ratio[si][0], samplename[si] ) ;
         l3->AddEntry( hflat_0lepldp_ratio_ave[0], "Average") ;
         for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

            cqcd3->cd(1+bbi) ;
            hflat_0lep[si][bbi] -> DrawCopy("histe") ;
            hflat_0lep[si][bbi] -> DrawCopy("samee") ;

            cqcd3->cd(4+bbi) ;
            hflat_ldp[si][bbi] -> DrawCopy("histe") ;
            hflat_ldp[si][bbi] -> DrawCopy("samee") ;


            cqcd3->cd(7+bbi) ;
            gPad->SetGridx(1) ;
            gPad->SetGridy(1) ;
            hflat_0lepldp_ratio[si][bbi]->DrawCopy() ;
            hflat_0lepldp_ratio_ave[bbi]->DrawCopy("same") ;
            line->DrawLine(0.5,0,nbins+0.5,0) ;
            l3->Draw() ;

         } // bbi.

         cqcd3->Update() ; cqcd3->Draw() ;

         if ( savePlots ) {
            char filename[1000] ;
            sprintf( filename, "outputfiles/qcd-study-%s.pdf", samplename[si] ) ;
            cqcd3 -> SaveAs( filename ) ;
            sprintf( filename, "outputfiles/qcd-study-%s.png", samplename[si] ) ;
            cqcd3 -> SaveAs( filename ) ;
         }
         if ( interactive ) {
            char a = getchar() ;
            if ( a == 'q') { return ; }
         }

      } // si.







     //--- Draw the average 0lep/LDP ratio with MC stat error and total error.

      TCanvas* cqcd4 = (TCanvas*) gDirectory->FindObject("cqcd4") ;
      if ( cqcd4 == 0x0 ) {
         cqcd4 = new TCanvas("cqcd4", "qcd study", 1700, 600 ) ;
      }

      gStyle->SetEndErrorSize(5) ;

      cqcd4->Clear() ;
      cqcd4->Divide(3,1) ;

      cqcd4 -> cd(1) ;
      hflat_0lepldp_ratio_ave[0] -> SetTitle("Average QCD 0lep/LDP ratio, nb=1") ;
      hflat_0lepldp_ratio_ave[0] -> SetMarkerSize(1.2) ;
      hflat_0lepldp_ratio_ave[0] -> SetMarkerStyle(20) ;
      hflat_0lepldp_ratio_ave[0] -> SetMaximum(0.6) ;
      hflat_0lepldp_ratio_ave[0] -> SetMinimum(-0.1) ;
      hflat_0lepldp_ratio_ave[0]->DrawCopy("e1") ;
      hflat_0lepldp_ratio_ave_withRMSerror[0]->SetMarkerStyle() ;
      hflat_0lepldp_ratio_ave_withRMSerror[0]->SetLineColor(4) ;
      hflat_0lepldp_ratio_ave_withRMSerror[0]->DrawCopy("samee1") ;
      hflat_0lepldp_ratio_ave[0]->DrawCopy("samee1") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd4 -> cd(2) ;
      hflat_0lepldp_ratio_ave[1] -> SetTitle("Average QCD 0lep/LDP ratio, nb=2") ;
      hflat_0lepldp_ratio_ave[1] -> SetMarkerSize(1.2) ;
      hflat_0lepldp_ratio_ave[1] -> SetMarkerStyle(20) ;
      hflat_0lepldp_ratio_ave[1] -> SetMaximum(0.6) ;
      hflat_0lepldp_ratio_ave[1] -> SetMinimum(-0.1) ;
      hflat_0lepldp_ratio_ave[1]->DrawCopy("e1") ;
      hflat_0lepldp_ratio_ave_withRMSerror[1]->SetMarkerStyle() ;
      hflat_0lepldp_ratio_ave_withRMSerror[1]->SetLineColor(4) ;
      hflat_0lepldp_ratio_ave_withRMSerror[1]->DrawCopy("samee1") ;
      hflat_0lepldp_ratio_ave[1]->DrawCopy("samee1") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd4 -> cd(3) ;
      hflat_0lepldp_ratio_ave[2] -> SetTitle("Average QCD 0lep/LDP ratio, nb>=3") ;
      hflat_0lepldp_ratio_ave[2] -> SetMarkerSize(1.2) ;
      hflat_0lepldp_ratio_ave[2] -> SetMarkerStyle(20) ;
      hflat_0lepldp_ratio_ave[2] -> SetMaximum(0.6) ;
      hflat_0lepldp_ratio_ave[2] -> SetMinimum(-0.1) ;
      hflat_0lepldp_ratio_ave[2]->DrawCopy("e1") ;
      hflat_0lepldp_ratio_ave_withRMSerror[2]->SetMarkerStyle() ;
      hflat_0lepldp_ratio_ave_withRMSerror[2]->SetLineColor(4) ;
      hflat_0lepldp_ratio_ave_withRMSerror[2]->DrawCopy("samee1") ;
      hflat_0lepldp_ratio_ave[2]->DrawCopy("samee1") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;


      cqcd4->Update() ; cqcd4->Draw() ;

      if ( savePlots ) {
         cqcd4 -> SaveAs("outputfiles/qcd-study-averatio.png") ;
         cqcd4 -> SaveAs("outputfiles/qcd-study-averatio.pdf") ;
      }












     //--- Draw the QCD scale factor with MC stat error and total error.

      TCanvas* cqcd5 = (TCanvas*) gDirectory->FindObject("cqcd5") ;
      if ( cqcd5 == 0x0 ) {
         cqcd5 = new TCanvas("cqcd5", "qcd study", 1700, 600 ) ;
      }

      gStyle->SetEndErrorSize(5) ;

      cqcd5->Clear() ;
      cqcd5->Divide(3,1) ;

      cqcd5 -> cd(1) ;
      hflat_scale_factor[0] -> SetTitle("QCD Scale Factor, nb=1") ;
      hflat_scale_factor[0] -> SetMarkerSize(1.2) ;
      hflat_scale_factor[0] -> SetMarkerStyle(20) ;
      hflat_scale_factor[0] -> SetMaximum(2.6) ;
      hflat_scale_factor[0] -> SetMinimum(-0.1) ;
      hflat_scale_factor[0]->DrawCopy("e1") ;
      hflat_scale_factor_withRMSerror[0]->SetMarkerStyle() ;
      hflat_scale_factor_withRMSerror[0]->SetLineColor(4) ;
      hflat_scale_factor_withRMSerror[0]->DrawCopy("samee1") ;
      hflat_scale_factor[0]->DrawCopy("samee1") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd5 -> cd(2) ;
      hflat_scale_factor[1] -> SetTitle("QCD Scale Factor, nb=2") ;
      hflat_scale_factor[1] -> SetMarkerSize(1.2) ;
      hflat_scale_factor[1] -> SetMarkerStyle(20) ;
      hflat_scale_factor[1] -> SetMaximum(2.6) ;
      hflat_scale_factor[1] -> SetMinimum(-0.1) ;
      hflat_scale_factor[1]->DrawCopy("e1") ;
      hflat_scale_factor_withRMSerror[1]->SetMarkerStyle() ;
      hflat_scale_factor_withRMSerror[1]->SetLineColor(4) ;
      hflat_scale_factor_withRMSerror[1]->DrawCopy("samee1") ;
      hflat_scale_factor[1]->DrawCopy("samee1") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd5 -> cd(3) ;
      hflat_scale_factor[2] -> SetTitle("QCD Scale Factor, nb>=3") ;
      hflat_scale_factor[2] -> SetMarkerSize(1.2) ;
      hflat_scale_factor[2] -> SetMarkerStyle(20) ;
      hflat_scale_factor[2] -> SetMaximum(2.6) ;
      hflat_scale_factor[2] -> SetMinimum(-0.1) ;
      hflat_scale_factor[2]->DrawCopy("e1") ;
      hflat_scale_factor_withRMSerror[2]->SetMarkerStyle() ;
      hflat_scale_factor_withRMSerror[2]->SetLineColor(4) ;
      hflat_scale_factor_withRMSerror[2]->DrawCopy("samee1") ;
      hflat_scale_factor[2]->DrawCopy("samee1") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;


      cqcd5->Update() ; cqcd5->Draw() ;

      if ( savePlots ) {
         cqcd5 -> SaveAs("outputfiles/qcd-study-scalefactor.png") ;
         cqcd5 -> SaveAs("outputfiles/qcd-study-scalefactor.pdf") ;
      }







   } // qcd_study.

  //==========================================================================================
  
  
    TH1F* bookHist(const char* hname, const char* htitle, int nBinsMET, int nBinsHT, int sampleIndex ) {
  
       int nbins = nBinsMET*(nBinsHT+1) + 1 ;
  
       TH1F* retVal = new TH1F( hname, htitle, nbins, 0.5+0.05*(sampleIndex-5), nbins+0.5+0.05*(sampleIndex-5) ) ;
       TAxis* xaxis = retVal->GetXaxis() ;
  
       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
             int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
             char binlabel[1000] ;
             sprintf( binlabel, "M%d_H%d", mbi+1, hbi+1 ) ;
             xaxis->SetBinLabel( histbin, binlabel ) ;
          } // hbi.
       } // mbi.
  
       retVal->SetLabelSize(0.055,"x") ;
       xaxis->LabelsOption("v") ;
  
       return retVal ;
  
    }
  
  //==========================================================================================


void saveHist(const char* filename, const char* pat)
{

  cout << "\n\n Saving histograms matching " << pat << " in file " << filename << "\n\n" << flush ;

  TList* list = gDirectory->GetList() ;
  TIterator* iter = list->MakeIterator();

  TRegexp re(pat,kTRUE) ;

  TFile outf(filename,"RECREATE") ;
  TObject* obj ;
  while((obj=iter->Next())) {
    if (TString(obj->GetName()).Index(re)>=0) {
      obj->Write() ;
      std::cout << "." ;
    }
  }
  std::cout << std::endl ;
  outf.Close() ;

  delete iter ;
}

//==========================================================================================


void loadHist(const char* filename, const char* pfx, const char* pat, Bool_t doAdd, Double_t scaleFactor)
{
  cout << " Reading histograms from file: " << filename << endl << flush ;
  TFile inf(filename) ;
  //inf.ReadAll() ;
  TList* list = inf.GetListOfKeys() ;
  TIterator* iter = list->MakeIterator();

  TRegexp re(pat,kTRUE) ;
  std::cout << "pat = " << pat << std::endl ;

  gDirectory->cd("Rint:") ;

  TObject* obj ;
  TKey* key ;
  std::cout << "doAdd = " << (doAdd?"T":"F") << std::endl ;
  std::cout << "loadHist: reading." ;
  while((key=(TKey*)iter->Next())) {
   
    Int_t ridx = TString(key->GetName()).Index(re) ;    
    if (ridx==-1) {
      continue ;
    }

    obj = inf.Get(key->GetName()) ;
    TObject* clone ;
    if (pfx) {

      // Find existing TH1-derived objects
      TObject* oldObj = 0 ;
      if (doAdd){
	oldObj = gDirectory->Get(Form("%s_%s",pfx,obj->GetName())) ;
	if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	  oldObj = 0 ;
	}
      }
      if (oldObj) {
	clone = oldObj ;
        if ( scaleFactor > 0 ) {
           ((TH1*)clone)->Sumw2() ;
           ((TH1*)clone)->Add((TH1*)obj, scaleFactor) ;
        } else {
           ((TH1*)clone)->Add((TH1*)obj) ;
        }
      } else {
	clone = obj->Clone(Form("%s_%s",pfx,obj->GetName())) ;
      }


    } else {

      // Find existing TH1-derived objects
      TObject* oldObj = 0 ;
      if (doAdd){
	oldObj = gDirectory->Get(key->GetName()) ;
	if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	  oldObj = 0 ;
	}
      }

      if (oldObj) {
	clone = oldObj ;
        if ( scaleFactor > 0 ) {
           ((TH1*)clone)->Sumw2() ;
           ((TH1*)clone)->Add((TH1*)obj, scaleFactor) ;
        } else {
           ((TH1*)clone)->Add((TH1*)obj) ;
        }
      } else {
	clone = obj->Clone() ;
      }
    }
    if ( scaleFactor > 0 && !doAdd ) {
       ((TH1*) clone)->Sumw2() ;
       ((TH1*) clone)->Scale(scaleFactor) ;
    }
    if (!gDirectory->GetList()->FindObject(clone)) {
      gDirectory->Append(clone) ;
    }
    std::cout << "." ;
    std::cout.flush() ;
  }
  std::cout << std::endl;
  inf.Close() ;
  delete iter ;
}

//==========================================================================================


