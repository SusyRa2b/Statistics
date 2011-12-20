
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TDirectory.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"

    void plot_ws_cls_hybrid1( const char* infile="cls-expected-2BT.root", int nStep=10, double minPoi=5., double deltaPoi=5. ) {

        gStyle->SetOptStat(0) ;
        gStyle->SetOptTitle(0) ;
        gStyle->SetLabelSize( 0.08, "x") ;

        TCanvas* canv_tstat = new TCanvas("canv_tstat", "Test Statistic distributions", 800,1000) ;
        TCanvas* canv_cls   = new TCanvas("canv_cls"  , "CLs vs mu_susy_sig", 900, 600 ) ;

        canv_tstat->Clear() ;
        canv_tstat->Divide(2,nStep/2) ;

        TText* label = new TText() ;
        label->SetTextSize(0.09) ;

        TText* pvallabel = new TText() ;
        pvallabel->SetTextSize(0.06) ;
        pvallabel->SetTextAlign(31) ;

        double bgopvals[100] ;
        double spbpvals[100] ;
        double clsvals[100] ;
        double medianclsvals[100] ;
        double poivals[100] ;
        double medianclserrlow1[100] ;
        double medianclserrlow2[100] ;
        double medianclserrhigh1[100] ;
        double medianclserrhigh2[100] ;
        double xerr[100] ;


        for ( int si=0; si<nStep; si++ ) {

           canv_tstat->cd(si+1) ;

           double poiVal = minPoi + si * deltaPoi ;
           poivals[si] = poiVal ;

           char bgottname[1000] ;
           char spbttname[1000] ;

           sprintf( bgottname, "toytt_%.0f_bgo", poiVal ) ;
           sprintf( spbttname, "toytt_%.0f_spb", poiVal ) ;

           TChain bgochain( bgottname ) ;
           TChain spbchain( spbttname ) ;

           bgochain.Add( infile ) ;
           spbchain.Add( infile ) ;

           bgochain.SetLineColor(4) ;
           bgochain.SetLineWidth(2) ;
           spbchain.SetLineColor(1) ;
           spbchain.SetFillColor(2) ;

           double nbgo_ok = bgochain.Draw("testStat","testStat>=0") ;
           double nbgo_worse = bgochain.Draw("testStat","testStat>=dataTestStat") ;
           double bgo_pval = nbgo_worse / nbgo_ok ;

           double nspb_ok = spbchain.Draw("testStat","testStat>=0") ;
           double nspb_worse = spbchain.Draw("testStat","testStat>=dataTestStat") ;
           double spb_pval = nspb_worse / nspb_ok ;

           printf("\n") ;
           printf(" %3.0f : BG-only  p-value : %4.0f / %4.0f = %5.3f\n", poiVal, nbgo_worse, nbgo_ok, bgo_pval ) ;
           printf(" %3.0f : S-plus-B p-value : %4.0f / %4.0f = %5.3f\n", poiVal, nspb_worse, nspb_ok, spb_pval ) ;
           printf(" %3.0f : CLs value        : %4.0f / %4.0f = %5.3f\n", poiVal, spb_pval, bgo_pval, spb_pval / bgo_pval  ) ;


           bgopvals[si] = bgo_pval ;
           spbpvals[si] = spb_pval ;
           clsvals[si] = spb_pval / bgo_pval ;

        //--- draw the pretty plot.

           char drawcmd[1000] ;

           sprintf( drawcmd, "testStat>>hbgo%d", si ) ;
           bgochain.Draw( drawcmd,"testStat>=0") ;
           sprintf( drawcmd, "testStat>>hspb%d", si ) ;
           spbchain.Draw( drawcmd,"","same") ;
           bgochain.Draw("testStat","testStat>=0","same") ;
           bgochain.SetLineColor(1) ;
           bgochain.SetFillColor(1) ;
           bgochain.Draw("dataTestStat","","same") ;

           char histname[1000] ;
           sprintf( histname, "hbgo%d", si ) ;
           TH1D* hbgo = (TH1D*) gDirectory->FindObject( histname ) ;
           sprintf( histname, "hspb%d", si ) ;
           TH1D* hspb = (TH1D*) gDirectory->FindObject( histname ) ;

           TLegend* legend = new TLegend(0.65,0.65, 0.9, 0.9) ;
           legend->SetFillColor(kWhite) ;
           legend->AddEntry( hbgo, "BG-only") ;
           legend->AddEntry( hspb, "S-plus-B") ;
           legend->Draw() ;

           char labeltext[1000] ;
           sprintf( labeltext, "mu_susy_sig = %.0f", poiVal ) ;
           label->DrawTextNDC(0.1,0.93, labeltext ) ;

           char pvaltext[1000] ;
           sprintf( pvaltext, "p-val: %4.0f / %4.0f = %5.3f", nbgo_worse, nbgo_ok, bgo_pval ) ;
           pvallabel->DrawTextNDC( 0.63, 0.80, pvaltext ) ;
           sprintf( pvaltext, "p-val: %4.0f / %4.0f = %5.3f", nspb_worse, nspb_ok, spb_pval ) ;
           pvallabel->DrawTextNDC( 0.63, 0.70, pvaltext ) ;

           canv_tstat->Update() ;



         //--- Make the expected CLs histogram for the background hypothesis.

           char clshname[1000] ;
           sprintf( clshname, "hcls_expected_%.0f", poiVal ) ;
           TH1F* hcls = new TH1F( clshname, clshname, 10000, 0., 2. ) ;
           int nBins = hbgo->GetNbinsX() ;
           for ( int bi=1; bi<=(nBins+1); bi++) {
              double bgo_pv = ( hbgo->Integral(bi,nBins+1) ) / nbgo_ok ;
              double spb_pv = ( hspb->Integral(bi,nBins+1) ) / nspb_ok ;
              double cls(99.) ;
              if ( bgo_pv > 0. ) { cls = spb_pv / bgo_pv ; }
              hcls->Fill( cls, hbgo->GetBinContent(bi) ) ;
           } // bi.

           double low2sig_cls(-10.) ;
           double low1sig_cls(-10.) ;
           double median_cls(-10.) ;
           double high1sig_cls(-10.) ;
           double high2sig_cls(-10.) ;

         //--- get the numbers for the expected CLs

           for ( int bi=1; bi<=(hcls->GetNbinsX()+1); bi++) {
              double frac = (hcls->Integral(1,bi)) / nbgo_ok ;
              double cls = hcls->GetBinCenter(bi) ;
              if ( frac >= 0.0455 && low2sig_cls < 0 ) { low2sig_cls = cls ; }
              if ( frac >= 0.3173 && low1sig_cls < 0 ) { low1sig_cls = cls ; }
              if ( frac >= 0.5000 && median_cls  < 0 ) { median_cls  = cls ; }
              if ( frac >= (1-0.3173) && high1sig_cls < 0 ) { high1sig_cls = cls ; }
              if ( frac >= (1-0.0455) && high2sig_cls < 0 ) { high2sig_cls = cls ; }
           } // bi.


           medianclsvals[si] = median_cls ;

           medianclserrlow1[si] = median_cls - low1sig_cls ;
           medianclserrlow2[si] = median_cls - low2sig_cls ;
           medianclserrhigh1[si] = high1sig_cls - median_cls ;
           medianclserrhigh2[si] = high2sig_cls - median_cls ;


           xerr[si] = deltaPoi/2. ;


           printf(" low2sig  (cls) :  %5.3f\n",  low2sig_cls ) ;
           printf(" low1sig  (cls) :  %5.3f\n",  low1sig_cls ) ;
           printf(" median   (cls) :  %5.3f\n",  median_cls ) ;
           printf(" high1sig (cls) :  %5.3f\n",  high1sig_cls ) ;
           printf(" high2sig (cls) :  %5.3f\n",  high2sig_cls ) ;

        } // si.

        canv_tstat->SaveAs("teststat-distributions.png") ;
        canv_tstat->SaveAs("teststat-distributions.pdf") ;


       //--- cls vs mu_susy_sig plot

        gStyle->SetLabelSize( 0.04, "x") ;

        TH2F* hdummy = new TH2F("hdummy","",2, minPoi-deltaPoi, minPoi+nStep*deltaPoi, 2, 0., 1. ) ;
        hdummy->SetXTitle("mu_susy_sig") ;
        hdummy->SetYTitle("p-value") ;

        TGraph* bgo_gr = new TGraph( nStep, poivals, bgopvals ) ;
        TGraph* spb_gr = new TGraph( nStep, poivals, spbpvals ) ;
        TGraph* cls_gr = new TGraph( nStep, poivals, clsvals ) ;
        TGraph* mediancls_gr = new TGraph( nStep, poivals, medianclsvals ) ;
        TGraphAsymmErrors* onesigband_gr = new TGraphAsymmErrors( nStep, poivals, medianclsvals, xerr, xerr, medianclserrlow1, medianclserrhigh1 ) ;
        TGraphAsymmErrors* twosigband_gr = new TGraphAsymmErrors( nStep, poivals, medianclsvals, xerr, xerr, medianclserrlow2, medianclserrhigh2 ) ;

        bgo_gr->SetMarkerStyle(20) ;
        spb_gr->SetMarkerStyle(20) ;
        cls_gr->SetMarkerStyle(20) ;
        spb_gr->SetMarkerColor(4) ;
        cls_gr->SetMarkerColor(2) ;
        bgo_gr->SetLineWidth(2) ;
        spb_gr->SetLineWidth(2) ;
        cls_gr->SetLineWidth(2) ;
        cls_gr->SetLineStyle(1) ;
        spb_gr->SetLineStyle(3) ;
        mediancls_gr->SetLineStyle(2) ;
        mediancls_gr->SetLineWidth(2) ;
        onesigband_gr->SetFillColor(3) ;
        twosigband_gr->SetFillColor(5) ;

        canv_cls->cd() ;

        hdummy->Draw() ;
        twosigband_gr->Draw("3") ;
        onesigband_gr->Draw("3") ;
        bgo_gr->Draw("LP") ;
        spb_gr->Draw("LP") ;
        cls_gr->Draw("LP") ;
        mediancls_gr->Draw("L") ;

        TLine* line = new TLine() ;
        line->SetLineColor(2) ;

        line->DrawLine(minPoi-deltaPoi, 0.05, minPoi+nStep*deltaPoi, 0.05) ;


        canv_cls->SaveAs("cls-vs-mu_susy_sig.png") ;
        canv_cls->SaveAs("cls-vs-mu_susy_sig.pdf") ;


    }





