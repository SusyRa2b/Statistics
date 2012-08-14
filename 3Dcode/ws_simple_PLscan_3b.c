
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooStats/ModelConfig.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "TAxis.h"
#include "TLine.h"
#include "TText.h"

#include <iostream>
#include <sstream>

  using namespace RooFit;
  using namespace RooStats;


// simple minded PL scan (integrating SUSY signal over all the bins)

void ws_simple_PLscan_3b( TString wsfile = "ws.root", double scanHigh = -1., double scanLow = 0. ) {


    TFile* wstf = new TFile( wsfile ) ;
    RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
    ws->Print() ;

    RooDataSet* rds = (RooDataSet*) ws->obj( "ra2b_observed_rds" ) ;
    cout << "\n\n\n  ===== RooDataSet ====================\n\n" << endl ;
    rds->Print() ;
    rds->printMultiline(cout, 1, kTRUE, "") ;



    TString parName = "mu_susy_all0lep" ;

    cout << "\n\n\n  ===== Grabbing " << parName << " rrv ====================\n\n" ;
    RooRealVar* rrv_susy_poi = ws->var( parName ) ;
    if ( rrv_susy_poi == 0x0 ) {
      cout << "\n\n\n *** can't find " << parName << " in workspace.  Quitting.\n\n\n" ; return ;
    } else {
      cout << " current value is : " << rrv_susy_poi->getVal() << endl ; cout << flush ;
    }


    cout << "\n\n\n  ===== Grabbing likelihood pdf ====================\n\n" ;
    RooAbsPdf* likelihood = ws->pdf("likelihood") ;
    if ( likelihood == 0x0 ) {
      cout << "\n\n\n *** can't find likelihood pdf in workspace.  Quitting.\n\n\n" ; return ;
    } else {
      cout << "\n\n likelihood pdf: \n\n" ;
      likelihood->Print() ;
    }


    cout << "\n\n\n  ===== Doing a fit ====================\n\n" ;

    rrv_susy_poi->setConstant(kFALSE);

    RooArgSet ras ;
    ras.add( *rrv_susy_poi ) ;
    RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(0), Hesse(true), Minos(ras), Strategy(1) ) ;
    double susy_poi_atMinNll = rrv_susy_poi->getVal() ;
    double susy_poi_plusErr = rrv_susy_poi->getErrorHi() ;
    double susy_poi_minusErr = rrv_susy_poi->getErrorLo() ;
    double susy_poi_err = rrv_susy_poi->getError() ;

    double minNllSusyFloat = fitResult->minNll() ;

    printf("\n\n  Best fit value of %s : %7.1f +/- %5.1f (+%5.1f, %5.1f)\n\n",
       parName.Data(), susy_poi_atMinNll, susy_poi_err, susy_poi_plusErr, susy_poi_minusErr ) ;


   if ( scanHigh < 0. ) {

    printf("\n\n ===== Finding the UL =====================\n\n") ;

      //-- Search for the 1-sided 95% CL UL in two steps of linear interpolation.
      //   Initially, use testStat at susy = +1.5 sigma and +2.0 sigma.

      float susy_poi_a = susy_poi_atMinNll + 1.5 * susy_poi_plusErr ;
      float susy_poi_b = susy_poi_atMinNll + 2.0 * susy_poi_plusErr ;

      rrv_susy_poi->setVal( susy_poi_a ) ;
      rrv_susy_poi->setConstant( kTRUE ) ;
      RooFitResult* fitResult_a = likelihood -> fitTo( *rds, Save(true), Hesse(false), Strategy(1), PrintLevel(0) ) ;
      double minNll_a = fitResult_a->minNll() ;
      double testStat_a = 2.*( minNll_a - minNllSusyFloat ) ;
      delete fitResult_a ;

      rrv_susy_poi->setVal( susy_poi_b ) ;
      rrv_susy_poi->setConstant( kTRUE ) ;
      RooFitResult* fitResult_b = likelihood -> fitTo( *rds, Save(true), Hesse(false), Strategy(1), PrintLevel(0) ) ;
      double minNll_b = fitResult_b->minNll() ;
      double testStat_b = 2.*( minNll_b - minNllSusyFloat ) ;
      delete fitResult_b ;




      double susy_poi_g1 = susy_poi_a + (2.70 - testStat_a) * (susy_poi_b-susy_poi_a)/(testStat_b-testStat_a) ;
      printf("  Test stat at +1.5 sigma = %5.2f, +2.0 sigma = %5.2f.  Will try %5.2f sigma.\n",
         testStat_a, testStat_b, (susy_poi_g1-susy_poi_atMinNll)/susy_poi_plusErr ) ;


      rrv_susy_poi->setVal( susy_poi_g1 ) ;
      rrv_susy_poi->setConstant( kTRUE ) ;
      RooFitResult* fitResult_g1 = likelihood -> fitTo( *rds, Save(true), Hesse(false), Strategy(1), PrintLevel(0) ) ;
      double minNll_g1 = fitResult_g1->minNll() ;
      double testStat_g1 = 2.*( minNll_g1 - minNllSusyFloat ) ;
      delete fitResult_g1 ;
      printf("  Test stat at first guess is %5.2f\n", testStat_g1 ) ;





      double susy_poi_g2 = susy_poi_a + (2.70 - testStat_a) * (susy_poi_g1-susy_poi_a)/(testStat_g1-testStat_a) ;
      printf("  Test stat at +1.5 sigma = %5.2f, +%5.2f sigma = %5.2f.  Will try %5.2f sigma.\n",
         testStat_a,
         (susy_poi_g1-susy_poi_atMinNll)/susy_poi_plusErr,
         testStat_g1, (susy_poi_g2-susy_poi_atMinNll)/susy_poi_plusErr ) ;

      rrv_susy_poi->setVal( susy_poi_g2 ) ;
      rrv_susy_poi->setConstant( kTRUE ) ;
      RooFitResult* fitResult_g2 = likelihood -> fitTo( *rds, Save(true), Hesse(false), Strategy(1), PrintLevel(0) ) ;
      double minNll_g2 = fitResult_g2->minNll() ;
      double testStat_g2 = 2.*( minNll_g2 - minNllSusyFloat ) ;
      delete fitResult_g2 ;
      printf("  Test stat at second guess is %5.2f\n", testStat_g2 ) ;

      double fit_susy_ul = susy_poi_g2 ;

      printf("  SUSY upper limit is : %5.2f 0lep events.\n", fit_susy_ul ) ;

      scanHigh = 1.2* fit_susy_ul ;

   } // no preset scanHigh?



      printf("\n\n ========== Doing the scan ==================== \n\n") ;

      int nScanPoints(20) ;

      double poiVals[20] ;
      double testStatVals[20] ;

      for ( int spi=0; spi<nScanPoints; spi++ ) {

         poiVals[spi] = scanLow + spi*(scanHigh-scanLow)/(nScanPoints-1.) ;

         rrv_susy_poi->setVal( poiVals[spi] ) ;
         rrv_susy_poi->setConstant( kTRUE ) ;


         RooFitResult* fitResult_sp = likelihood -> fitTo( *rds, Save(true), Hesse(false), Strategy(1), PrintLevel(0) ) ;
         double minNll_sp = fitResult_sp->minNll() ;
         testStatVals[spi] = 2.*( minNll_sp - minNllSusyFloat ) ;
         delete fitResult_sp ;

         printf("  %3d : %s = %6.1f,  test stat = %7.3f\n", spi, parName.Data(), poiVals[spi], testStatVals[spi] ) ;


      } // spi.


    TGraph *ScanPlot = new TGraph( nScanPoints, poiVals, testStatVals ) ;

    //--- Use interpolation to find +/- 1 sigma and UL points.

      double tsAtPoiMin(99.) ;
      double poiMin(-1) ;
      double poiMinusOneSigma(-1.) ;
      double poiPlusOneSigma(-1.) ;
      double poiUL(-1.) ;
      int nsearch(1000) ;
      printf("\n\n") ;
      for ( int i=0; i<nsearch; i++ ) {
         double poi = scanLow + i*(scanHigh-scanLow)/(nsearch-1.) ;
         double tsVal = ScanPlot->Eval( poi, 0, "S" ) ;
         if ( poi < susy_poi_atMinNll ) {
            if ( tsVal <= 1.0 && poiMinusOneSigma < 0. ) {
               poiMinusOneSigma = poi ;
               printf(" -1 sigma   : %s = %7.2f\n", parName.Data(), poi) ;
            }
         } else {
            if ( tsVal >= 1.0 && poiPlusOneSigma < 0. ) {
               poiPlusOneSigma = poi ;
               printf(" +1 sigma   : %s = %7.2f\n", parName.Data(), poi) ;
            }
            if ( tsVal >= 2.71 && poiUL < 0. ) {
               poiUL = poi ;
               printf(" Upper lim. : %s = %7.2f\n", parName.Data(), poi) ;
            }
         }
         if ( tsVal < tsAtPoiMin ) {
            tsAtPoiMin = tsVal ;
            poiMin = poi ;
         }
         /// printf("  %4d : %s = %6.1f,  test stat = %7.3f\n", i, parName.Data(), poi, tsVal ) ;
      }
      printf(" Minimum    : %s = %7.2f\n", parName.Data(), poiMin ) ;
      double poiSignif(-1.) ;
      if ( scanLow <= 0 ) {
         double ts = ScanPlot->Eval( 0., 0, "S" ) ;
         if ( ts >=0. ) { poiSignif = sqrt( ts ) ; }
         printf(" Signif     : %5.2f\n", poiSignif ) ;
      }
      printf("\n\n") ;



    ScanPlot->GetXaxis()->SetTitle("N SUSY 0lep");
    ScanPlot->GetYaxis()->SetTitle("test statistic");

    ScanPlot->SetMarkerStyle(20);
    ScanPlot->SetMarkerColor(2);
    ScanPlot->SetLineColor(2);
    ScanPlot->SetLineWidth(3);

    ScanPlot->SetTitle("Profile likelihood scan of SUSY 0lep yield") ;

    TLine* line = new TLine() ;
    line->SetLineColor(4) ;

    TCanvas *c0 = new TCanvas("c0","simple PL scan",700,600);

    ScanPlot->Draw("AC");

    line->DrawLine(scanLow, 1.0, poiPlusOneSigma, 1.0) ;
    line->DrawLine(scanLow, 2.71, poiUL, 2.71) ;
    line->DrawLine(poiMinusOneSigma,0.,poiMinusOneSigma,1.0) ;
    line->DrawLine(poiPlusOneSigma,0.,poiPlusOneSigma,1.0) ;
    line->DrawLine(poiUL,0.,poiUL,2.71) ;

    TText* text = new TText() ;
    text->SetTextSize(0.04) ;
    char tstring[1000] ;

    sprintf( tstring, "SUSY 0lep yield = %5.1f +%5.1f, -%5.1f", poiMin, poiPlusOneSigma-poiMin, poiMin-poiMinusOneSigma ) ;
    text->DrawTextNDC( 0.15, 0.85, tstring ) ;
    sprintf( tstring, "SUSY 0lep 95%% CL 1-sided UL = %5.1f ", poiUL ) ;
    text->DrawTextNDC( 0.15, 0.80, tstring ) ;
    if ( scanLow == 0. && poiSignif>-1. ) {
       sprintf( tstring, "Significance = %5.2f", poiSignif ) ;
       text->DrawTextNDC( 0.15, 0.75, tstring ) ;
    }

    gPad->SetGridx(1) ;
    gPad->SetGridy(1) ;


    TString savestr( wsfile ) ;
    savestr.ReplaceAll(".root","") ;

    char savename[1000] ;


    sprintf( savename, "%s-susy-PL-scan.pdf", savestr.Data() ) ;
    printf(" Saving as %s\n", savename ) ;
    c0->SaveAs( savename );

    sprintf( savename, "%s-susy-PL-scan.gif", savestr.Data() ) ;
    printf(" Saving as %s\n", savename ) ;
    c0->SaveAs( savename );



      printf("\n\n") ;
      for ( int spi=0; spi<nScanPoints; spi++ ) {
         printf("  %3d : %s = %6.1f,  test stat = %7.3f\n", spi, parName.Data(), poiVals[spi], testStatVals[spi] ) ;
      } // spi.
      printf("\n\n") ;

    printf("\n\n  Best fit value of %s : %7.1f +/- %5.1f (+%5.1f, %5.1f)\n\n",
       parName.Data(), susy_poi_atMinNll, susy_poi_err, susy_poi_plusErr, susy_poi_minusErr ) ;


    return ;

}
