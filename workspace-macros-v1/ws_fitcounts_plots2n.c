
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
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

   void ws_fitcounts_plots2n( const char* wsfile = "output-files/ws-newfit-lm9-1BL.root",
                            double mu_susy_sig_1b_val = -1 ) {

      gStyle->SetPaintTextFormat("6.1f") ;
      gStyle->SetOptStat(0) ;


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







      //--- unpack observables.

       int dataNsig_1b(0) ;
       int dataNsb_1b(0) ;
       int dataNsig_sl_1b(0) ;
       int dataNsb_sl_1b(0) ;
       int dataNsig_ldp_1b(0) ;
       int dataNsb_ldp_1b(0) ;

       int dataNsig_2b(0) ;
       int dataNsb_2b(0) ;
       int dataNsig_sl_2b(0) ;
       int dataNsb_sl_2b(0) ;
       int dataNsig_ldp_2b(0) ;
       int dataNsb_ldp_2b(0) ;

       int dataNsig_3b(0) ;
       int dataNsb_3b(0) ;
       int dataNsig_sl_3b(0) ;
       int dataNsb_sl_3b(0) ;
       int dataNsig_ldp_3b(0) ;
       int dataNsb_ldp_3b(0) ;

       int dataNsig_ee(0) ;
       int dataNsb_ee(0) ;
       int dataNsig_mm(0) ;
       int dataNsb_mm(0) ;



       const RooArgSet* dsras = rds->get() ;
       TIterator* obsIter = dsras->createIterator() ;
       while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {

          if ( strcmp( obs->GetName(), "Nsig_1b"     ) == 0 ) { dataNsig_1b     = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_1b"      ) == 0 ) { dataNsb_1b      = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsig_sl_1b"  ) == 0 ) { dataNsig_sl_1b  = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_sl_1b"   ) == 0 ) { dataNsb_sl_1b   = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsig_ldp_1b" ) == 0 ) { dataNsig_ldp_1b = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_ldp_1b"  ) == 0 ) { dataNsb_ldp_1b  = obs->getVal() ; }

          if ( strcmp( obs->GetName(), "Nsig_2b"     ) == 0 ) { dataNsig_2b     = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_2b"      ) == 0 ) { dataNsb_2b      = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsig_sl_2b"  ) == 0 ) { dataNsig_sl_2b  = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_sl_2b"   ) == 0 ) { dataNsb_sl_2b   = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsig_ldp_2b" ) == 0 ) { dataNsig_ldp_2b = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_ldp_2b"  ) == 0 ) { dataNsb_ldp_2b  = obs->getVal() ; }

          if ( strcmp( obs->GetName(), "Nsig_3b"     ) == 0 ) { dataNsig_3b     = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_3b"      ) == 0 ) { dataNsb_3b      = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsig_sl_3b"  ) == 0 ) { dataNsig_sl_3b  = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_sl_3b"   ) == 0 ) { dataNsb_sl_3b   = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsig_ldp_3b" ) == 0 ) { dataNsig_ldp_3b = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_ldp_3b"  ) == 0 ) { dataNsb_ldp_3b  = obs->getVal() ; }

          if ( strcmp( obs->GetName(), "Nsig_ee"  ) == 0 ) { dataNsig_ee  = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_ee"  ) == 0 )  { dataNsb_ee  = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsig_mm"  ) == 0 ) { dataNsig_mm  = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_mm"  ) == 0 )  { dataNsb_mm  = obs->getVal() ; }

       }





       TH2F* hfitcounts_data = new TH2F("hfitqual_data", "RA2b likelihood fit results, data", 6, 0.5, 6.5,  3, 0.5, 3.5 ) ;

       TAxis* xaxis = hfitcounts_data->GetXaxis() ;
       xaxis->SetBinLabel( 1, "SIG" ) ;
       xaxis->SetBinLabel( 2, "SB" ) ;
       xaxis->SetBinLabel( 3, "SL,SIG" ) ;
       xaxis->SetBinLabel( 4, "SL,SB" ) ;
       xaxis->SetBinLabel( 5, "LDP,SIG" ) ;
       xaxis->SetBinLabel( 6, "LDP,SB" ) ;
       TAxis* yaxis = hfitcounts_data->GetYaxis() ;
       yaxis->SetBinLabel( 1, "= 1 b") ;
       yaxis->SetBinLabel( 2, "= 2 b") ;
       yaxis->SetBinLabel( 3, ">= 3 b") ;

       TH2F* hfitcounts_ttwj = (TH2F*) hfitcounts_data->Clone("hfitcounts_ttwj") ;
       TH2F* hfitcounts_qcd  = (TH2F*) hfitcounts_data->Clone("hfitcounts_qcd") ;
       TH2F* hfitcounts_znn  = (TH2F*) hfitcounts_data->Clone("hfitcounts_znn") ;
       TH2F* hfitcounts_susy = (TH2F*) hfitcounts_data->Clone("hfitcounts_susy") ;
       TH2F* hfitcounts_fit  = (TH2F*) hfitcounts_data->Clone("hfitcounts_fit") ;
       TH2F* hfitcounts_allsm  = (TH2F*) hfitcounts_data->Clone("hfitcounts_allsm") ;





       double n_sig_1b     = ((RooRealVar*) ws->obj("n_sig_1b"    )) -> getVal() ;
       double mu_ttwj_sig_1b     = ((RooRealVar*) ws->obj("mu_ttwj_sig_1b"    )) -> getVal() ;
       double mu_znn_sig_1b      = ((RooRealVar*) ws->obj("mu_znn_sig_1b"     )) -> getVal() ;
       double mu_qcd_sig_1b      = ((RooRealVar*) ws->obj("mu_qcd_sig_1b"     )) -> getVal() ;
       double mu_susy_sig_1b     = ((RooRealVar*) ws->obj("mu_susy_sig_1b"    )) -> getVal() ;
       double eff_sf_sig_1b      = ((RooRealVar*) ws->obj("eff_sf_sig_1b"     )) -> getVal() ;
       double btageff_sf_sig_1b      = ((RooRealVar*) ws->obj("btageff_sf_sig_1b"     )) -> getVal() ;

       double n_sig_2b     = ((RooRealVar*) ws->obj("n_sig_2b"    )) -> getVal() ;
       double mu_ttwj_sig_2b     = ((RooRealVar*) ws->obj("mu_ttwj_sig_2b"    )) -> getVal() ;
       double mu_znn_sig_2b      = ((RooRealVar*) ws->obj("mu_znn_sig_2b"     )) -> getVal() ;
       double mu_qcd_sig_2b      = ((RooRealVar*) ws->obj("mu_qcd_sig_2b"     )) -> getVal() ;
       double mu_susy_sig_2b     = ((RooRealVar*) ws->obj("mu_susy_sig_2b"    )) -> getVal() ;
       double eff_sf_sig_2b      = ((RooRealVar*) ws->obj("eff_sf_sig_2b"     )) -> getVal() ;
       double btageff_sf_sig_2b      = ((RooRealVar*) ws->obj("btageff_sf_sig_2b"     )) -> getVal() ;

       double n_sig_3b     = ((RooRealVar*) ws->obj("n_sig_3b"    )) -> getVal() ;
       double mu_ttwj_sig_3b     = ((RooRealVar*) ws->obj("mu_ttwj_sig_3b"    )) -> getVal() ;
       double mu_znn_sig_3b      = ((RooRealVar*) ws->obj("mu_znn_sig_3b"     )) -> getVal() ;
       double mu_qcd_sig_3b      = ((RooRealVar*) ws->obj("mu_qcd_sig_3b"     )) -> getVal() ;
       double mu_susy_sig_3b     = ((RooRealVar*) ws->obj("mu_susy_sig_3b"    )) -> getVal() ;
       double eff_sf_sig_3b      = ((RooRealVar*) ws->obj("eff_sf_sig_3b"     )) -> getVal() ;
       double btageff_sf_sig_3b      = ((RooRealVar*) ws->obj("btageff_sf_sig_3b"     )) -> getVal() ;



       double n_sb_1b     = ((RooRealVar*) ws->obj("n_sb_1b"    )) -> getVal() ;
       double mu_ttwj_sb_1b     = ((RooRealVar*) ws->obj("mu_ttwj_sb_1b"    )) -> getVal() ;
       double mu_znn_sb_1b      = ((RooRealVar*) ws->obj("mu_znn_sb_1b"     )) -> getVal() ;
       double mu_qcd_sb_1b      = ((RooRealVar*) ws->obj("mu_qcd_sb_1b"     )) -> getVal() ;
       double mu_susy_sb_1b     = ((RooRealVar*) ws->obj("mu_susy_sb_1b"    )) -> getVal() ;
       double eff_sf_sb_1b      = ((RooRealVar*) ws->obj("eff_sf_sb_1b"     )) -> getVal() ;
       double btageff_sf_sb_1b      = ((RooRealVar*) ws->obj("btageff_sf_sb_1b"     )) -> getVal() ;

       double n_sb_2b     = ((RooRealVar*) ws->obj("n_sb_2b"    )) -> getVal() ;
       double mu_ttwj_sb_2b     = ((RooRealVar*) ws->obj("mu_ttwj_sb_2b"    )) -> getVal() ;
       double mu_znn_sb_2b      = ((RooRealVar*) ws->obj("mu_znn_sb_2b"     )) -> getVal() ;
       double mu_qcd_sb_2b      = ((RooRealVar*) ws->obj("mu_qcd_sb_2b"     )) -> getVal() ;
       double mu_susy_sb_2b     = ((RooRealVar*) ws->obj("mu_susy_sb_2b"    )) -> getVal() ;
       double eff_sf_sb_2b      = ((RooRealVar*) ws->obj("eff_sf_sb_2b"     )) -> getVal() ;
       double btageff_sf_sb_2b      = ((RooRealVar*) ws->obj("btageff_sf_sb_2b"     )) -> getVal() ;

       double n_sb_3b     = ((RooRealVar*) ws->obj("n_sb_3b"    )) -> getVal() ;
       double mu_ttwj_sb_3b     = ((RooRealVar*) ws->obj("mu_ttwj_sb_3b"    )) -> getVal() ;
       double mu_znn_sb_3b      = ((RooRealVar*) ws->obj("mu_znn_sb_3b"     )) -> getVal() ;
       double mu_qcd_sb_3b      = ((RooRealVar*) ws->obj("mu_qcd_sb_3b"     )) -> getVal() ;
       double mu_susy_sb_3b     = ((RooRealVar*) ws->obj("mu_susy_sb_3b"    )) -> getVal() ;
       double eff_sf_sb_3b      = ((RooRealVar*) ws->obj("eff_sf_sb_3b"     )) -> getVal() ;
       double btageff_sf_sb_3b      = ((RooRealVar*) ws->obj("btageff_sf_sb_3b"     )) -> getVal() ;



       double n_sig_sl_1b     = ((RooRealVar*) ws->obj("n_sig_sl_1b"    )) -> getVal() ;
       double mu_ttwj_sig_sl_1b     = ((RooRealVar*) ws->obj("mu_ttwj_sig_sl_1b"    )) -> getVal() ;
       double mu_susy_sig_sl_1b     = ((RooRealVar*) ws->obj("mu_susy_sig_sl_1b"    )) -> getVal() ;
       double eff_sf_sig_sl_1b      = ((RooRealVar*) ws->obj("eff_sf_sig_sl_1b"     )) -> getVal() ;
       double btageff_sf_sig_sl_1b      = ((RooRealVar*) ws->obj("btageff_sf_sig_sl_1b"     )) -> getVal() ;

       double n_sig_sl_2b     = ((RooRealVar*) ws->obj("n_sig_sl_2b"    )) -> getVal() ;
       double mu_ttwj_sig_sl_2b     = ((RooRealVar*) ws->obj("mu_ttwj_sig_sl_2b"    )) -> getVal() ;
       double mu_susy_sig_sl_2b     = ((RooRealVar*) ws->obj("mu_susy_sig_sl_2b"    )) -> getVal() ;
       double eff_sf_sig_sl_2b      = ((RooRealVar*) ws->obj("eff_sf_sig_sl_2b"     )) -> getVal() ;
       double btageff_sf_sig_sl_2b      = ((RooRealVar*) ws->obj("btageff_sf_sig_sl_2b"     )) -> getVal() ;

       double n_sig_sl_3b     = ((RooRealVar*) ws->obj("n_sig_sl_3b"    )) -> getVal() ;
       double mu_ttwj_sig_sl_3b     = ((RooRealVar*) ws->obj("mu_ttwj_sig_sl_3b"    )) -> getVal() ;
       double mu_susy_sig_sl_3b     = ((RooRealVar*) ws->obj("mu_susy_sig_sl_3b"    )) -> getVal() ;
       double eff_sf_sig_sl_3b      = ((RooRealVar*) ws->obj("eff_sf_sig_sl_3b"     )) -> getVal() ;
       double btageff_sf_sig_sl_3b      = ((RooRealVar*) ws->obj("btageff_sf_sig_sl_3b"     )) -> getVal() ;





       double n_sb_sl_1b     = ((RooRealVar*) ws->obj("n_sb_sl_1b"    )) -> getVal() ;
       double mu_ttwj_sb_sl_1b     = ((RooRealVar*) ws->obj("mu_ttwj_sb_sl_1b"    )) -> getVal() ;
       double mu_susy_sb_sl_1b     = ((RooRealVar*) ws->obj("mu_susy_sb_sl_1b"    )) -> getVal() ;
       double eff_sf_sb_sl_1b      = ((RooRealVar*) ws->obj("eff_sf_sb_sl_1b"     )) -> getVal() ;
       double btageff_sf_sb_sl_1b      = ((RooRealVar*) ws->obj("btageff_sf_sb_sl_1b"     )) -> getVal() ;

       double n_sb_sl_2b     = ((RooRealVar*) ws->obj("n_sb_sl_2b"    )) -> getVal() ;
       double mu_ttwj_sb_sl_2b     = ((RooRealVar*) ws->obj("mu_ttwj_sb_sl_2b"    )) -> getVal() ;
       double mu_susy_sb_sl_2b     = ((RooRealVar*) ws->obj("mu_susy_sb_sl_2b"    )) -> getVal() ;
       double eff_sf_sb_sl_2b      = ((RooRealVar*) ws->obj("eff_sf_sb_sl_2b"     )) -> getVal() ;
       double btageff_sf_sb_sl_2b      = ((RooRealVar*) ws->obj("btageff_sf_sb_sl_2b"     )) -> getVal() ;

       double n_sb_sl_3b     = ((RooRealVar*) ws->obj("n_sb_sl_3b"    )) -> getVal() ;
       double mu_ttwj_sb_sl_3b     = ((RooRealVar*) ws->obj("mu_ttwj_sb_sl_3b"    )) -> getVal() ;
       double mu_susy_sb_sl_3b     = ((RooRealVar*) ws->obj("mu_susy_sb_sl_3b"    )) -> getVal() ;
       double eff_sf_sb_sl_3b      = ((RooRealVar*) ws->obj("eff_sf_sb_sl_3b"     )) -> getVal() ;
       double btageff_sf_sb_sl_3b      = ((RooRealVar*) ws->obj("btageff_sf_sb_sl_3b"     )) -> getVal() ;





       double sf_mc      = ((RooRealVar*) ws->obj("sf_mc"     )) -> getVal() ;

       double n_sig_ldp_1b     = ((RooRealVar*) ws->obj("n_sig_ldp_1b"    )) -> getVal() ;
       double mu_ttwj_sig_ldp_1b     = ((RooRealVar*) ws->obj("mu_ttwj_sig_ldp_1b"    )) -> getVal() ;
       double mu_znn_sig_ldp_1b      = ((RooRealVar*) ws->obj("mu_znn_sig_ldp_1b"     )) -> getVal() ;
       double mu_qcd_sig_ldp_1b      = ((RooRealVar*) ws->obj("mu_qcd_sig_ldp_1b"     )) -> getVal() ;
       double mu_susy_sig_ldp_1b     = ((RooRealVar*) ws->obj("mu_susy_sig_ldp_1b"    )) -> getVal() ;
       double eff_sf_sig_ldp_1b      = ((RooRealVar*) ws->obj("eff_sf_sig_ldp_1b"     )) -> getVal() ;
       double btageff_sf_sig_ldp_1b      = ((RooRealVar*) ws->obj("btageff_sf_sig_ldp_1b"     )) -> getVal() ;

       double n_sig_ldp_2b     = ((RooRealVar*) ws->obj("n_sig_ldp_2b"    )) -> getVal() ;
       double mu_ttwj_sig_ldp_2b     = ((RooRealVar*) ws->obj("mu_ttwj_sig_ldp_2b"    )) -> getVal() ;
       double mu_znn_sig_ldp_2b      = ((RooRealVar*) ws->obj("mu_znn_sig_ldp_2b"     )) -> getVal() ;
       double mu_qcd_sig_ldp_2b      = ((RooRealVar*) ws->obj("mu_qcd_sig_ldp_2b"     )) -> getVal() ;
       double mu_susy_sig_ldp_2b     = ((RooRealVar*) ws->obj("mu_susy_sig_ldp_2b"    )) -> getVal() ;
       double eff_sf_sig_ldp_2b      = ((RooRealVar*) ws->obj("eff_sf_sig_ldp_2b"     )) -> getVal() ;
       double btageff_sf_sig_ldp_2b      = ((RooRealVar*) ws->obj("btageff_sf_sig_ldp_2b"     )) -> getVal() ;

       double n_sig_ldp_3b     = ((RooRealVar*) ws->obj("n_sig_ldp_3b"    )) -> getVal() ;
       double mu_ttwj_sig_ldp_3b     = ((RooRealVar*) ws->obj("mu_ttwj_sig_ldp_3b"    )) -> getVal() ;
       double mu_znn_sig_ldp_3b      = ((RooRealVar*) ws->obj("mu_znn_sig_ldp_3b"     )) -> getVal() ;
       double mu_qcd_sig_ldp_3b      = ((RooRealVar*) ws->obj("mu_qcd_sig_ldp_3b"     )) -> getVal() ;
       double mu_susy_sig_ldp_3b     = ((RooRealVar*) ws->obj("mu_susy_sig_ldp_3b"    )) -> getVal() ;
       double eff_sf_sig_ldp_3b      = ((RooRealVar*) ws->obj("eff_sf_sig_ldp_3b"     )) -> getVal() ;
       double btageff_sf_sig_ldp_3b      = ((RooRealVar*) ws->obj("btageff_sf_sig_ldp_3b"     )) -> getVal() ;




       double n_sb_ldp_1b     = ((RooRealVar*) ws->obj("n_sb_ldp_1b"    )) -> getVal() ;
       double mu_ttwj_sb_ldp_1b     = ((RooRealVar*) ws->obj("mu_ttwj_sb_ldp_1b"    )) -> getVal() ;
       double mu_znn_sb_ldp_1b      = ((RooRealVar*) ws->obj("mu_znn_sb_ldp_1b"     )) -> getVal() ;
       double mu_qcd_sb_ldp_1b      = ((RooRealVar*) ws->obj("mu_qcd_sb_ldp_1b"     )) -> getVal() ;
       double mu_susy_sb_ldp_1b     = ((RooRealVar*) ws->obj("mu_susy_sb_ldp_1b"    )) -> getVal() ;
       double eff_sf_sb_ldp_1b      = ((RooRealVar*) ws->obj("eff_sf_sb_ldp_1b"     )) -> getVal() ;
       double btageff_sf_sb_ldp_1b      = ((RooRealVar*) ws->obj("btageff_sf_sb_ldp_1b"     )) -> getVal() ;

       double n_sb_ldp_2b     = ((RooRealVar*) ws->obj("n_sb_ldp_2b"    )) -> getVal() ;
       double mu_ttwj_sb_ldp_2b     = ((RooRealVar*) ws->obj("mu_ttwj_sb_ldp_2b"    )) -> getVal() ;
       double mu_znn_sb_ldp_2b      = ((RooRealVar*) ws->obj("mu_znn_sb_ldp_2b"     )) -> getVal() ;
       double mu_qcd_sb_ldp_2b      = ((RooRealVar*) ws->obj("mu_qcd_sb_ldp_2b"     )) -> getVal() ;
       double mu_susy_sb_ldp_2b     = ((RooRealVar*) ws->obj("mu_susy_sb_ldp_2b"    )) -> getVal() ;
       double eff_sf_sb_ldp_2b      = ((RooRealVar*) ws->obj("eff_sf_sb_ldp_2b"     )) -> getVal() ;
       double btageff_sf_sb_ldp_2b      = ((RooRealVar*) ws->obj("btageff_sf_sb_ldp_2b"     )) -> getVal() ;

       double n_sb_ldp_3b     = ((RooRealVar*) ws->obj("n_sb_ldp_3b"    )) -> getVal() ;
       double mu_ttwj_sb_ldp_3b     = ((RooRealVar*) ws->obj("mu_ttwj_sb_ldp_3b"    )) -> getVal() ;
       double mu_znn_sb_ldp_3b      = ((RooRealVar*) ws->obj("mu_znn_sb_ldp_3b"     )) -> getVal() ;
       double mu_qcd_sb_ldp_3b      = ((RooRealVar*) ws->obj("mu_qcd_sb_ldp_3b"     )) -> getVal() ;
       double mu_susy_sb_ldp_3b     = ((RooRealVar*) ws->obj("mu_susy_sb_ldp_3b"    )) -> getVal() ;
       double eff_sf_sb_ldp_3b      = ((RooRealVar*) ws->obj("eff_sf_sb_ldp_3b"     )) -> getVal() ;
       double btageff_sf_sb_ldp_3b      = ((RooRealVar*) ws->obj("btageff_sf_sb_ldp_3b"     )) -> getVal() ;





       hfitcounts_ttwj->SetBinContent(1,1, mu_ttwj_sig_1b ) ;
       hfitcounts_ttwj->SetBinContent(1,2, mu_ttwj_sig_2b ) ;
       hfitcounts_ttwj->SetBinContent(1,3, mu_ttwj_sig_3b ) ;

       hfitcounts_ttwj->SetBinContent(2,1, mu_ttwj_sb_1b ) ;
       hfitcounts_ttwj->SetBinContent(2,2, mu_ttwj_sb_2b ) ;
       hfitcounts_ttwj->SetBinContent(2,3, mu_ttwj_sb_3b ) ;

       hfitcounts_ttwj->SetBinContent(3,1, mu_ttwj_sig_sl_1b ) ;
       hfitcounts_ttwj->SetBinContent(3,2, mu_ttwj_sig_sl_2b ) ;
       hfitcounts_ttwj->SetBinContent(3,3, mu_ttwj_sig_sl_3b ) ;

       hfitcounts_ttwj->SetBinContent(4,1, mu_ttwj_sb_sl_1b ) ;
       hfitcounts_ttwj->SetBinContent(4,2, mu_ttwj_sb_sl_2b ) ;
       hfitcounts_ttwj->SetBinContent(4,3, mu_ttwj_sb_sl_3b ) ;

       hfitcounts_ttwj->SetBinContent(5,1, btageff_sf_sig_ldp_1b * eff_sf_sig_ldp_1b * sf_mc * mu_ttwj_sig_ldp_1b ) ;
       hfitcounts_ttwj->SetBinContent(5,2, btageff_sf_sig_ldp_2b * eff_sf_sig_ldp_2b * sf_mc * mu_ttwj_sig_ldp_2b ) ;
       hfitcounts_ttwj->SetBinContent(5,3, btageff_sf_sig_ldp_3b * eff_sf_sig_ldp_3b * sf_mc * mu_ttwj_sig_ldp_3b ) ;

       hfitcounts_ttwj->SetBinContent(6,1, btageff_sf_sb_ldp_1b * eff_sf_sb_ldp_1b * sf_mc * mu_ttwj_sb_ldp_1b ) ;
       hfitcounts_ttwj->SetBinContent(6,2, btageff_sf_sb_ldp_2b * eff_sf_sb_ldp_2b * sf_mc * mu_ttwj_sb_ldp_2b ) ;
       hfitcounts_ttwj->SetBinContent(6,3, btageff_sf_sb_ldp_3b * eff_sf_sb_ldp_3b * sf_mc * mu_ttwj_sb_ldp_3b ) ;






       hfitcounts_znn->SetBinContent(1,1, mu_znn_sig_1b ) ;
       hfitcounts_znn->SetBinContent(1,2, mu_znn_sig_2b ) ;
       hfitcounts_znn->SetBinContent(1,3, mu_znn_sig_3b ) ;

       hfitcounts_znn->SetBinContent(2,1, mu_znn_sb_1b ) ;
       hfitcounts_znn->SetBinContent(2,2, mu_znn_sb_2b ) ;
       hfitcounts_znn->SetBinContent(2,3, mu_znn_sb_3b ) ;

       hfitcounts_znn->SetBinContent(3,1, 0 ) ;
       hfitcounts_znn->SetBinContent(3,2, 0 ) ;
       hfitcounts_znn->SetBinContent(3,3, 0 ) ;

       hfitcounts_znn->SetBinContent(4,1, 0 ) ;
       hfitcounts_znn->SetBinContent(4,2, 0 ) ;
       hfitcounts_znn->SetBinContent(4,3, 0 ) ;

       hfitcounts_znn->SetBinContent(5,1, btageff_sf_sig_ldp_1b * eff_sf_sig_ldp_1b * sf_mc * mu_znn_sig_ldp_1b ) ;
       hfitcounts_znn->SetBinContent(5,2, btageff_sf_sig_ldp_2b * eff_sf_sig_ldp_2b * sf_mc * mu_znn_sig_ldp_2b ) ;
       hfitcounts_znn->SetBinContent(5,3, btageff_sf_sig_ldp_3b * eff_sf_sig_ldp_3b * sf_mc * mu_znn_sig_ldp_3b ) ;

       hfitcounts_znn->SetBinContent(6,1, btageff_sf_sb_ldp_1b * eff_sf_sb_ldp_1b * sf_mc * mu_znn_sb_ldp_1b ) ;
       hfitcounts_znn->SetBinContent(6,2, btageff_sf_sb_ldp_2b * eff_sf_sb_ldp_2b * sf_mc * mu_znn_sb_ldp_2b ) ;
       hfitcounts_znn->SetBinContent(6,3, btageff_sf_sb_ldp_3b * eff_sf_sb_ldp_3b * sf_mc * mu_znn_sb_ldp_3b ) ;






       hfitcounts_qcd->SetBinContent(1,1, mu_qcd_sig_1b ) ;
       hfitcounts_qcd->SetBinContent(1,2, mu_qcd_sig_2b ) ;
       hfitcounts_qcd->SetBinContent(1,3, mu_qcd_sig_3b ) ;

       hfitcounts_qcd->SetBinContent(2,1, mu_qcd_sb_1b ) ;
       hfitcounts_qcd->SetBinContent(2,2, mu_qcd_sb_2b ) ;
       hfitcounts_qcd->SetBinContent(2,3, mu_qcd_sb_3b ) ;

       hfitcounts_qcd->SetBinContent(3,1, 0 ) ;
       hfitcounts_qcd->SetBinContent(3,2, 0 ) ;
       hfitcounts_qcd->SetBinContent(3,3, 0 ) ;

       hfitcounts_qcd->SetBinContent(4,1, 0 ) ;
       hfitcounts_qcd->SetBinContent(4,2, 0 ) ;
       hfitcounts_qcd->SetBinContent(4,3, 0 ) ;

       hfitcounts_qcd->SetBinContent(5,1, mu_qcd_sig_ldp_1b ) ;
       hfitcounts_qcd->SetBinContent(5,2, mu_qcd_sig_ldp_2b ) ;
       hfitcounts_qcd->SetBinContent(5,3, mu_qcd_sig_ldp_3b ) ;

       hfitcounts_qcd->SetBinContent(6,1, mu_qcd_sb_ldp_1b ) ;
       hfitcounts_qcd->SetBinContent(6,2, mu_qcd_sb_ldp_2b ) ;
       hfitcounts_qcd->SetBinContent(6,3, mu_qcd_sb_ldp_3b ) ;





       hfitcounts_susy->SetBinContent(1,1, btageff_sf_sig_1b * eff_sf_sig_1b * mu_susy_sig_1b ) ;
       hfitcounts_susy->SetBinContent(1,2, btageff_sf_sig_2b * eff_sf_sig_2b * mu_susy_sig_2b ) ;
       hfitcounts_susy->SetBinContent(1,3, btageff_sf_sig_3b * eff_sf_sig_3b * mu_susy_sig_3b ) ;

       hfitcounts_susy->SetBinContent(2,1, btageff_sf_sb_1b * eff_sf_sb_1b * mu_susy_sb_1b ) ;
       hfitcounts_susy->SetBinContent(2,2, btageff_sf_sb_2b * eff_sf_sb_2b * mu_susy_sb_2b ) ;
       hfitcounts_susy->SetBinContent(2,3, btageff_sf_sb_3b * eff_sf_sb_3b * mu_susy_sb_3b ) ;

       hfitcounts_susy->SetBinContent(3,1, btageff_sf_sig_sl_1b * eff_sf_sig_sl_1b * mu_susy_sig_sl_1b ) ;
       hfitcounts_susy->SetBinContent(3,2, btageff_sf_sig_sl_2b * eff_sf_sig_sl_2b * mu_susy_sig_sl_2b ) ;
       hfitcounts_susy->SetBinContent(3,3, btageff_sf_sig_sl_3b * eff_sf_sig_sl_3b * mu_susy_sig_sl_3b ) ;

       hfitcounts_susy->SetBinContent(4,1, btageff_sf_sb_sl_1b * eff_sf_sb_sl_1b * mu_susy_sb_sl_1b ) ;
       hfitcounts_susy->SetBinContent(4,2, btageff_sf_sb_sl_2b * eff_sf_sb_sl_2b * mu_susy_sb_sl_2b ) ;
       hfitcounts_susy->SetBinContent(4,3, btageff_sf_sb_sl_3b * eff_sf_sb_sl_3b * mu_susy_sb_sl_3b ) ;

       hfitcounts_susy->SetBinContent(5,1, btageff_sf_sig_ldp_1b * eff_sf_sig_ldp_1b * mu_susy_sig_ldp_1b ) ;
       hfitcounts_susy->SetBinContent(5,2, btageff_sf_sig_ldp_2b * eff_sf_sig_ldp_2b * mu_susy_sig_ldp_2b ) ;
       hfitcounts_susy->SetBinContent(5,3, btageff_sf_sig_ldp_3b * eff_sf_sig_ldp_3b * mu_susy_sig_ldp_3b ) ;

       hfitcounts_susy->SetBinContent(6,1, btageff_sf_sb_ldp_1b * eff_sf_sb_ldp_1b * mu_susy_sb_ldp_1b ) ;
       hfitcounts_susy->SetBinContent(6,2, btageff_sf_sb_ldp_2b * eff_sf_sb_ldp_2b * mu_susy_sb_ldp_2b ) ;
       hfitcounts_susy->SetBinContent(6,3, btageff_sf_sb_ldp_3b * eff_sf_sb_ldp_3b * mu_susy_sb_ldp_3b ) ;



       hfitcounts_allsm->SetBinContent(1,1,  hfitcounts_ttwj->GetBinContent(1,1) + hfitcounts_znn->GetBinContent(1,1) + hfitcounts_qcd->GetBinContent(1,1) ) ;
       hfitcounts_allsm->SetBinContent(2,1,  hfitcounts_ttwj->GetBinContent(2,1) + hfitcounts_znn->GetBinContent(2,1) + hfitcounts_qcd->GetBinContent(2,1) ) ;
       hfitcounts_allsm->SetBinContent(3,1,  hfitcounts_ttwj->GetBinContent(3,1) + hfitcounts_znn->GetBinContent(3,1) + hfitcounts_qcd->GetBinContent(3,1) ) ;
       hfitcounts_allsm->SetBinContent(4,1,  hfitcounts_ttwj->GetBinContent(4,1) + hfitcounts_znn->GetBinContent(4,1) + hfitcounts_qcd->GetBinContent(4,1) ) ;
       hfitcounts_allsm->SetBinContent(5,1,  hfitcounts_ttwj->GetBinContent(5,1) + hfitcounts_znn->GetBinContent(5,1) + hfitcounts_qcd->GetBinContent(5,1) ) ;
       hfitcounts_allsm->SetBinContent(6,1,  hfitcounts_ttwj->GetBinContent(6,1) + hfitcounts_znn->GetBinContent(6,1) + hfitcounts_qcd->GetBinContent(6,1) ) ;

       hfitcounts_allsm->SetBinContent(1,2,  hfitcounts_ttwj->GetBinContent(1,2) + hfitcounts_znn->GetBinContent(1,2) + hfitcounts_qcd->GetBinContent(1,2) ) ;
       hfitcounts_allsm->SetBinContent(2,2,  hfitcounts_ttwj->GetBinContent(2,2) + hfitcounts_znn->GetBinContent(2,2) + hfitcounts_qcd->GetBinContent(2,2) ) ;
       hfitcounts_allsm->SetBinContent(3,2,  hfitcounts_ttwj->GetBinContent(3,2) + hfitcounts_znn->GetBinContent(3,2) + hfitcounts_qcd->GetBinContent(3,2) ) ;
       hfitcounts_allsm->SetBinContent(4,2,  hfitcounts_ttwj->GetBinContent(4,2) + hfitcounts_znn->GetBinContent(4,2) + hfitcounts_qcd->GetBinContent(4,2) ) ;
       hfitcounts_allsm->SetBinContent(5,2,  hfitcounts_ttwj->GetBinContent(5,2) + hfitcounts_znn->GetBinContent(5,2) + hfitcounts_qcd->GetBinContent(5,2) ) ;
       hfitcounts_allsm->SetBinContent(6,2,  hfitcounts_ttwj->GetBinContent(6,2) + hfitcounts_znn->GetBinContent(6,2) + hfitcounts_qcd->GetBinContent(6,2) ) ;

       hfitcounts_allsm->SetBinContent(1,3,  hfitcounts_ttwj->GetBinContent(1,3) + hfitcounts_znn->GetBinContent(1,3) + hfitcounts_qcd->GetBinContent(1,3) ) ;
       hfitcounts_allsm->SetBinContent(2,3,  hfitcounts_ttwj->GetBinContent(2,3) + hfitcounts_znn->GetBinContent(2,3) + hfitcounts_qcd->GetBinContent(2,3) ) ;
       hfitcounts_allsm->SetBinContent(3,3,  hfitcounts_ttwj->GetBinContent(3,3) + hfitcounts_znn->GetBinContent(3,3) + hfitcounts_qcd->GetBinContent(3,3) ) ;
       hfitcounts_allsm->SetBinContent(4,3,  hfitcounts_ttwj->GetBinContent(4,3) + hfitcounts_znn->GetBinContent(4,3) + hfitcounts_qcd->GetBinContent(4,3) ) ;
       hfitcounts_allsm->SetBinContent(5,3,  hfitcounts_ttwj->GetBinContent(5,3) + hfitcounts_znn->GetBinContent(5,3) + hfitcounts_qcd->GetBinContent(5,3) ) ;
       hfitcounts_allsm->SetBinContent(6,3,  hfitcounts_ttwj->GetBinContent(6,3) + hfitcounts_znn->GetBinContent(6,3) + hfitcounts_qcd->GetBinContent(6,3) ) ;




       hfitcounts_fit->SetBinContent(1,1, n_sig_1b ) ;
       hfitcounts_fit->SetBinContent(1,2, n_sig_2b ) ;
       hfitcounts_fit->SetBinContent(1,3, n_sig_3b ) ;

       hfitcounts_fit->SetBinContent(2,1, n_sb_1b ) ;
       hfitcounts_fit->SetBinContent(2,2, n_sb_2b ) ;
       hfitcounts_fit->SetBinContent(2,3, n_sb_3b ) ;

       hfitcounts_fit->SetBinContent(3,1, n_sig_sl_1b ) ;
       hfitcounts_fit->SetBinContent(3,2, n_sig_sl_2b ) ;
       hfitcounts_fit->SetBinContent(3,3, n_sig_sl_3b ) ;

       hfitcounts_fit->SetBinContent(4,1, n_sb_sl_1b ) ;
       hfitcounts_fit->SetBinContent(4,2, n_sb_sl_2b ) ;
       hfitcounts_fit->SetBinContent(4,3, n_sb_sl_3b ) ;

       hfitcounts_fit->SetBinContent(5,1, n_sig_ldp_1b ) ;
       hfitcounts_fit->SetBinContent(5,2, n_sig_ldp_2b ) ;
       hfitcounts_fit->SetBinContent(5,3, n_sig_ldp_3b ) ;

       hfitcounts_fit->SetBinContent(6,1, n_sb_ldp_1b ) ;
       hfitcounts_fit->SetBinContent(6,2, n_sb_ldp_2b ) ;
       hfitcounts_fit->SetBinContent(6,3, n_sb_ldp_3b ) ;




       hfitcounts_data->SetBinContent(1,1, dataNsig_1b ) ;
       hfitcounts_data->SetBinContent(1,2, dataNsig_2b ) ;
       hfitcounts_data->SetBinContent(1,3, dataNsig_3b ) ;

       hfitcounts_data->SetBinContent(2,1, dataNsb_1b ) ;
       hfitcounts_data->SetBinContent(2,2, dataNsb_2b ) ;
       hfitcounts_data->SetBinContent(2,3, dataNsb_3b ) ;

       hfitcounts_data->SetBinContent(3,1, dataNsig_sl_1b ) ;
       hfitcounts_data->SetBinContent(3,2, dataNsig_sl_2b ) ;
       hfitcounts_data->SetBinContent(3,3, dataNsig_sl_3b ) ;

       hfitcounts_data->SetBinContent(4,1, dataNsb_sl_1b ) ;
       hfitcounts_data->SetBinContent(4,2, dataNsb_sl_2b ) ;
       hfitcounts_data->SetBinContent(4,3, dataNsb_sl_3b ) ;

       hfitcounts_data->SetBinContent(5,1, dataNsig_ldp_1b ) ;
       hfitcounts_data->SetBinContent(5,2, dataNsig_ldp_2b ) ;
       hfitcounts_data->SetBinContent(5,3, dataNsig_ldp_3b ) ;

       hfitcounts_data->SetBinContent(6,1, dataNsb_ldp_1b ) ;
       hfitcounts_data->SetBinContent(6,2, dataNsb_ldp_2b ) ;
       hfitcounts_data->SetBinContent(6,3, dataNsb_ldp_3b ) ;

      double msize(2.5) ;
      double xlabelsize(0.06) ;
      double ylabelsize(0.08) ;

      hfitcounts_znn->SetLabelSize( xlabelsize, "x" ) ;
      hfitcounts_qcd->SetLabelSize( xlabelsize, "x" ) ;
      hfitcounts_susy->SetLabelSize( xlabelsize, "x" ) ;
      hfitcounts_allsm->SetLabelSize( xlabelsize, "x" ) ;
      hfitcounts_ttwj->SetLabelSize( xlabelsize, "x" ) ;
      hfitcounts_data->SetLabelSize( xlabelsize, "x" ) ;
      hfitcounts_fit->SetLabelSize( xlabelsize, "x" ) ;

      hfitcounts_znn->SetLabelSize( ylabelsize, "y" ) ;
      hfitcounts_qcd->SetLabelSize( ylabelsize, "y" ) ;
      hfitcounts_susy->SetLabelSize( ylabelsize, "y" ) ;
      hfitcounts_allsm->SetLabelSize( ylabelsize, "y" ) ;
      hfitcounts_ttwj->SetLabelSize( ylabelsize, "y" ) ;
      hfitcounts_data->SetLabelSize( ylabelsize, "y" ) ;
      hfitcounts_fit->SetLabelSize( ylabelsize, "y" ) ;

      hfitcounts_znn->SetMarkerSize( msize ) ;
      hfitcounts_qcd->SetMarkerSize( msize ) ;
      hfitcounts_susy->SetMarkerSize( msize ) ;
      hfitcounts_allsm->SetMarkerSize( msize ) ;
      hfitcounts_ttwj->SetMarkerSize( msize ) ;
      hfitcounts_data->SetMarkerSize( msize ) ;
      hfitcounts_fit->SetMarkerSize( msize ) ;

      hfitcounts_znn->SetFillColor(kGreen-3) ;
      hfitcounts_qcd->SetFillColor(2) ;
      hfitcounts_susy->SetFillColor(6) ;
      hfitcounts_allsm->SetFillColor(kMagenta+1) ;
      hfitcounts_ttwj->SetFillColor(kBlue-9) ;
      hfitcounts_data->SetFillColor(18) ;
      hfitcounts_fit->SetFillColor(11) ;

      hfitcounts_allsm->SetTitle("Fit, ttwj+QCD+Znn") ;
      hfitcounts_ttwj->SetTitle("Fit ttbar + Wjets") ;
      hfitcounts_qcd->SetTitle("Fit QCD") ;
      hfitcounts_znn->SetTitle("Fit Zinvis") ;
      hfitcounts_susy->SetTitle("Fit SUSY") ;
      hfitcounts_data->SetTitle("Nobs") ;
      hfitcounts_fit->SetTitle("Fit, all components") ;

      gStyle->SetPaintTextFormat("6.1f") ;


      TCanvas* canvas = new TCanvas("canvas","N bjet vs met bin, fit", 1700, 700 ) ;

      canvas->Divide(4,2) ;

      canvas->cd(1) ;
      hfitcounts_susy->Draw("box") ;
      hfitcounts_susy->Draw("textsame") ;

      canvas->cd(5) ;
      hfitcounts_data->Draw("box") ;
      hfitcounts_data->Draw("textsame") ;

      canvas->cd(2) ;
      hfitcounts_allsm->Draw("box") ;
      hfitcounts_allsm->Draw("textsame") ;

      canvas->cd(6) ;
      hfitcounts_ttwj->Draw("box") ;
      hfitcounts_ttwj->Draw("textsame") ;

      canvas->cd(4) ;

      canvas->cd(8) ;
      hfitcounts_fit->Draw("box") ;
      hfitcounts_fit->Draw("textsame") ;

      canvas->cd(3) ;
      hfitcounts_znn->Draw("box") ;
      hfitcounts_znn->Draw("textsame") ;

      canvas->cd(7) ;
      hfitcounts_qcd->Draw("box") ;
      hfitcounts_qcd->Draw("textsame") ;


      canvas->SaveAs("nbjetvsmetbin-fit.png") ;
      canvas->SaveAs("nbjetvsmetbin-fit.pdf") ;





   }






