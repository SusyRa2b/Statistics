
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TAxis.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TGAxis.h"
#include "TText.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooStats/ModelConfig.h"

#include <string.h>

  using namespace RooFit;
  using namespace RooStats;

  //------
  //
  // Note: If the 2nd argument (mu_susy_sig_val) is negative, mu_susy_sig will be floated
  //        in the fit.  If it's zero or positive, mu_susy_sig will be fixed to that value
  //        in the fit.
  //
  //------

   void ws_fitqual_plots2n_nuisancepars( const char* wsfile = "output-files/ws-newfit-lm9-1BL.root",
                            double mu_susy_sig_1b_val = 75., double yMax = 3.5 ) {


       char sel[100] ;
       if ( strstr( wsfile, "1BL" ) != 0 ) {
          sprintf( sel, "Loose HT,MET" ) ;
       } else if ( strstr( wsfile, "1BT" ) != 0 ) {
          sprintf( sel, "1BT" ) ;
       } else if ( strstr( wsfile, "2BL" ) != 0 ) {
          sprintf( sel, "2BL" ) ;
       } else if ( strstr( wsfile, "2BT" ) != 0 ) {
          sprintf( sel, "2BT" ) ;
       } else if ( strstr( wsfile, "3B" ) != 0 ) {
          sprintf( sel, "3B" ) ;
       } else {
          printf("\n\n\n *** can't figure out which selection this is.  I quit.\n\n" ) ;
          return ;
       }
       printf("\n\n selection is %s\n\n", sel ) ;

       if ( strstr( wsfile, "newfit") == 0 ) {
          printf("\n\n\n *** Input ws file is not new simultaneous fit of bjet bins.\n\n\n") ;
          return ;
       }



       TFile* wstf = new TFile( wsfile ) ;

       RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
       ws->Print() ;






       ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;

       printf("\n\n\n  ===== SbModel ====================\n\n") ;
       modelConfig->Print() ;







       RooDataSet* rds = (RooDataSet*) ws->obj( "ra2b_observed_rds" ) ;
       printf("\n\n\n  ===== RooDataSet ====================\n\n") ;

       rds->Print() ;
       rds->printMultiline(cout, 1, kTRUE, "") ;





       RooAbsPdf* likelihood = modelConfig->GetPdf() ;

       RooRealVar* rrv_mu_susy_sig_1b = ws->var("mu_susy_sig_1b") ;
       if ( rrv_mu_susy_sig_1b == 0x0 ) {
          printf("\n\n\n *** can't find mu_susy_sig_1b in workspace.  Quitting.\n\n\n") ;
          return ;
       } else {
          printf(" current value is : %8.3f\n", rrv_mu_susy_sig_1b->getVal() ) ; cout << flush ;
          rrv_mu_susy_sig_1b->setConstant(kFALSE) ;
       }


       printf("\n\n\n  ===== Doing a fit with SUSY component floating ====================\n\n") ;

       RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true) ) ;
       double logLikelihoodSusyFloat = fitResult->minNll() ;

       double logLikelihoodSusyFixed(0.) ;
       double testStatVal(-1.) ;
       if ( mu_susy_sig_1b_val >= 0. ) {
          printf("\n\n\n  ===== Doing a fit with SUSY fixed ====================\n\n") ;
          printf(" fixing mu_susy_sig_1b to %8.2f.\n", mu_susy_sig_1b_val ) ;
          rrv_mu_susy_sig_1b->setVal( mu_susy_sig_1b_val ) ;
          rrv_mu_susy_sig_1b->setConstant(kTRUE) ;

          fitResult = likelihood->fitTo( *rds, Save(true) ) ;
          logLikelihoodSusyFixed = fitResult->minNll() ;
          testStatVal = 2.*(logLikelihoodSusyFixed - logLikelihoodSusyFloat) ;
          printf("\n\n\n ======= test statistic : -2 * ln (L_fixed / L_max) = %8.3f\n\n\n", testStatVal ) ;
       }


       printf("\n ==== Final floating parameter values\n\n") ;
       const RooArgList fitFloatVals = fitResult->floatParsFinal() ;
       {
         TIterator* parIter = fitFloatVals.createIterator() ;
         while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {
            printf(" %20s : %8.2f\n", par->GetName(), par->getVal() ) ;
         }
       }


       printf("\n ==== Constant parameter values\n\n") ;
       const RooArgList fitConstVals = fitResult->constPars() ;
       {
         TIterator* parIter = fitConstVals.createIterator() ;
         while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {
            printf(" %20s : %8.2f\n", par->GetName(), par->getVal() ) ;
         }
       }

       printf("\n ==== Function values\n\n") ;
       RooArgSet funcs = ws->allFunctions() ;
       TIterator* funcIter = funcs.createIterator() ;
       while ( RooFormulaVar* func = (RooFormulaVar*) funcIter->Next() ) {
          printf(" %20s : %8.2f\n", func->GetName(), func->getVal() ) ;
       }

       printf("\n\n") ;





      gStyle->SetPadBottomMargin(0.35) ;
      gStyle->SetPadRightMargin(0.20) ;


       int nbins(36) ;

       TH1F* hfitqual_susy = new TH1F("hfitqual_susy", "RA2b likelihood fit results, susy", nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_ttwj = new TH1F("hfitqual_ttwj", "RA2b likelihood fit results, ttwj", nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_qcd  = new TH1F("hfitqual_qcd" , "RA2b likelihood fit results, qcd" , nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_znn  = new TH1F("hfitqual_znn" , "RA2b likelihood fit results, znn" , nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_mc   = new TH1F("hfitqual_mc"  , "RA2b likelihood fit results, mc" , nbins, 0.5, nbins+0.5 ) ;

       hfitqual_ttwj  -> SetFillColor(kBlue-9) ;
       hfitqual_qcd   -> SetFillColor(2) ;
       hfitqual_znn   -> SetFillColor(kGreen-3) ;
       hfitqual_susy  -> SetFillColor(6) ;
       hfitqual_mc  -> SetFillColor(kOrange+1) ;



       TAxis* xaxis = hfitqual_ttwj->GetXaxis() ;


       char binLabel[1000] ;
       int  binIndex ;

       double value ;

       binIndex = 1 ;


     //-- SUSY ----------------------------------------------------------------------

       //-- Eff sf

       sprintf( binLabel, "SUSY,eff" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("eff_sf_prim")) -> getVal() ;

       hfitqual_susy -> SetBinContent( binIndex, value ) ;


       //-- btag Eff sf

       sprintf( binLabel, "SUSY,btag eff" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("btageff_sf_prim")) -> getVal() ;

       hfitqual_susy -> SetBinContent( binIndex, value ) ;





       binIndex++ ;

     //-- ttwj ----------------------------------------------------------------------

       //-- ttwj, SIG, 1b

       sprintf( binLabel, "ttwj,SIG - 1b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_ttwj_sig_1b_prim")) -> getVal() ;

       hfitqual_ttwj -> SetBinContent( binIndex, value ) ;


       //-- ttwj, SIG, 2b

       sprintf( binLabel, "ttwj,SIG - 2b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_ttwj_sig_2b_prim")) -> getVal() ;

       hfitqual_ttwj -> SetBinContent( binIndex, value ) ;


       //-- ttwj, SIG, 3b

       sprintf( binLabel, "ttwj,SIG - 3b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_ttwj_sig_3b_prim")) -> getVal() ;

       hfitqual_ttwj -> SetBinContent( binIndex, value ) ;





       //-- ttwj, SB, 2b

       sprintf( binLabel, "ttwj,SB - 2b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_ttwj_sb_2b_prim")) -> getVal() ;

       hfitqual_ttwj -> SetBinContent( binIndex, value ) ;


       //-- ttwj, SB, 3b

       sprintf( binLabel, "ttwj,SB - 3b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_ttwj_sb_3b_prim")) -> getVal() ;

       hfitqual_ttwj -> SetBinContent( binIndex, value ) ;



       binIndex++ ;

     //-- Znn ----------------------------------------------------------------------

       //-- Znn, ll

       sprintf( binLabel, "Znn, ll" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_ll_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;


       //-- Znn, knn, SIG, 1b

       sprintf( binLabel, "Znn, F,SIG - 1b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("knn_sig_1b_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;



       //-- Znn, knn, SIG, 2b

       sprintf( binLabel, "Znn, F,SIG - 2b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("knn_sig_2b_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;



       //-- Znn, knn, SIG, 3b

       sprintf( binLabel, "Znn, F,SIG - 3b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("knn_sig_3b_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;



       //-- Znn, knn, SB, 1b

       sprintf( binLabel, "Znn, F,SB - 1b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("knn_sb_1b_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;



       //-- Znn, knn, SB, 2b

       sprintf( binLabel, "Znn, F,SB - 2b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("knn_sb_2b_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;



       //-- Znn, knn, SB, 3b

       sprintf( binLabel, "Znn, F,SB - 3b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("knn_sb_3b_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;



       //-- Znn, Acc, SIG, ee

       sprintf( binLabel, "Znn, Acc,SIG,ee" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("acc_ee_sig_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;



       //-- Znn, Acc, SIG, mm

       sprintf( binLabel, "Znn, Acc,SIG,mm" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("acc_mm_sig_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;




       //-- Znn, Acc, SB, ee

       sprintf( binLabel, "Znn, Acc,SB,ee" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("acc_ee_sb_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;



       //-- Znn, Acc, SB, mm

       sprintf( binLabel, "Znn, Acc,SB,mm" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("acc_mm_sb_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;



       //-- Znn, Eff, ee

       sprintf( binLabel, "Znn, Eff,ee" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("eff_ee_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;



       //-- Znn, Eff, mm

       sprintf( binLabel, "Znn, Eff,mm" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("eff_mm_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;



       //-- Znn, Purity,ee

       sprintf( binLabel, "Znn, Purity,ee" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("fsig_ee_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;



       //-- Znn, Purity,mm

       sprintf( binLabel, "Znn, Purity,mm" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("fsig_mm_prim")) -> getVal() ;

       hfitqual_znn -> SetBinContent( binIndex, value ) ;






       binIndex++ ;

     //-- QCD ----------------------------------------------------------------------

       //-- QCD, Rlsb

       sprintf( binLabel, "QCD, Rlsb" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("Rlsb_passfail_prim")) -> getVal() ;

       hfitqual_qcd -> SetBinContent( binIndex, value ) ;





       //-- QCD, SIG,1b

       sprintf( binLabel, "QCD, SIG - 1b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_qcd_sig_1b_prim")) -> getVal() ;

       hfitqual_qcd -> SetBinContent( binIndex, value ) ;



       //-- QCD, SIG,2b

       sprintf( binLabel, "QCD, SIG - 2b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_qcd_sig_2b_prim")) -> getVal() ;

       hfitqual_qcd -> SetBinContent( binIndex, value ) ;



       //-- QCD, SIG,3b

       sprintf( binLabel, "QCD, SIG - 3b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_qcd_sig_3b_prim")) -> getVal() ;

       hfitqual_qcd -> SetBinContent( binIndex, value ) ;






       //-- QCD, SB,1b

       sprintf( binLabel, "QCD, SB - 1b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_qcd_sb_1b_prim")) -> getVal() ;

       hfitqual_qcd -> SetBinContent( binIndex, value ) ;



       //-- QCD, SB,2b

       sprintf( binLabel, "QCD, SB - 2b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_qcd_sb_2b_prim")) -> getVal() ;

       hfitqual_qcd -> SetBinContent( binIndex, value ) ;



       //-- QCD, SB,3b

       sprintf( binLabel, "QCD, SB - 3b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_qcd_sb_3b_prim")) -> getVal() ;

       hfitqual_qcd -> SetBinContent( binIndex, value ) ;






       binIndex++ ;

     //-- MC ----------------------------------------------------------------------

       //-- MC sf

       sprintf( binLabel, "MC, sf" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;

       value = ((RooRealVar*) ws->obj("sf_mc_prim")) -> getVal() ;

       hfitqual_mc -> SetBinContent( binIndex, value ) ;


    //-------------------------------------------------------------------------------




      printf("\n\n\n") ;


     //--- final formatting and drawing.

      TCanvas* cfitqual = new TCanvas("cfitqualnp","RA2b fit quality, nuisance pars", 1000, 550 ) ;

      gStyle->SetOptStat(0) ;
      gStyle->SetOptTitle(0) ;

      hfitqual_ttwj->SetLabelSize(0.055,"x") ;
      hfitqual_ttwj->GetXaxis()->LabelsOption("v") ;

      hfitqual_ttwj->SetMinimum( -yMax ) ;
      hfitqual_ttwj->SetMaximum(  yMax ) ;

      hfitqual_ttwj->Draw() ;
      hfitqual_susy->Draw("same") ;
      hfitqual_znn->Draw("same") ;
      hfitqual_qcd->Draw("same") ;
      hfitqual_mc->Draw("same") ;

      TLegend* legend = new TLegend(0.85,0.2,0.97,0.9) ;

      legend->AddEntry( hfitqual_susy,  "SUSY" ) ;
      legend->AddEntry( hfitqual_ttwj, "ttwj" ) ;
      legend->AddEntry( hfitqual_znn,   "Znunu" ) ;
      legend->AddEntry( hfitqual_qcd,   "QCD" ) ;
      legend->AddEntry( hfitqual_mc,   "MC" ) ;




      gPad->SetTicks(1,0) ;


      gPad->SetGridy(1) ;
      legend->Draw() ;

      cfitqual->Update() ;





      TText* fittext = new TText() ;
      fittext->SetTextSize(0.045) ;




      TPad* bigPad = (TPad*) gPad ;




   //--- title

      bigPad->cd() ;

      char titletext[10000] ;
      if ( mu_susy_sig_1b_val < 0. ) {
         sprintf( titletext, "%s : SUSY floating", sel ) ;
      } else {
         sprintf( titletext, "%s : SUSY fixed to %.1f", sel, mu_susy_sig_1b_val ) ;
      }
      fittext->SetTextSize(0.050) ;
      fittext->DrawTextNDC( 0.02, 0.93, titletext ) ;

      if ( testStatVal >=0 ) {
         sprintf( titletext, "-2 * ln (L_fixed / L_max) = %8.3f\n", testStatVal ) ;
         fittext->SetTextSize(0.035) ;
         fittext->DrawTextNDC( 0.12, 0.85, titletext ) ;
      }

      cfitqual->Update() ;

      cfitqual->SaveAs("fitqual-nuisancepars.png") ;
      cfitqual->SaveAs("fitqual-nuisancepars.pdf") ;


   }






