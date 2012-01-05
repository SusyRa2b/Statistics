
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

   void ws_fitqual_plots2n( const char* wsfile = "output-files/ws-newfit-lm9-1BL.root",
                            double mu_susy_sig_1b_val = 0., bool doNorm = false ) {

      double hmax = 1.5 ;

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








       double ttwjSig1bVal, ttwjSb1bVal ;

       double qcdSig1bVal,  qcdSb1bVal  ;

       double ttwjSig1bErr ;

       double qcdSig1bErr ;


     //-- Note: Depending on the first 2 arguments used in ra2bRoostatsClass*.c, you either
     //          have the SIG or SB parameter (but not both) for ttwj and qcd.
     //
     //          To look at the sig vars, use ra2bRoostatsClass7( true, false )
     //          To look at the sb  vars, use ra2bRoostatsClass7( false, true )
     //
     //         These are set in the make_*_ws_*.c macros
     //         (e.g. make_lm9_ws_ge1btight.c).
     //

     //--- sort out which configuration was used (SIG or SB vars for ttwj and qcd).

       TObject* ttwj_sig_1b_obj = ws->obj("mu_ttwj_sig_1b") ;
       TObject* ttwj_sb_1b_obj = ws->obj("mu_ttwj_sb_1b") ;

       if ( ttwj_sig_1b_obj->IsA()->InheritsFrom(RooRealVar::Class()) )    {
          ttwjSig1bVal = ((RooRealVar*) ttwj_sig_1b_obj) -> getVal() ;
          ttwjSig1bErr = ((RooRealVar*) ttwj_sig_1b_obj) -> getError() ;
          ttwjSb1bVal = ((RooFormulaVar*) ttwj_sb_1b_obj) -> getVal() ;
          printf(" mu_ttwj_sig_1b is a RooRealVar.\n" ) ; 
       } else if ( ttwj_sig_1b_obj->IsA()->InheritsFrom(RooFormulaVar::Class()) ) {
          ttwjSig1bVal = ((RooFormulaVar*) ttwj_sig_1b_obj) -> getVal() ;
          ttwjSig1bErr = -1. ;
          ttwjSb1bVal = ((RooRealVar*) ttwj_sb_1b_obj) -> getVal() ;
          printf(" mu_ttwj_sig_1b is a RooFormulaVar\n" ) ;
       } else {
          printf("\n\n\n *** what kind of class is mu_ttwj_sig_1b ???\n\n\n\n") ;
          return ;
       }
       printf(" mu_ttwj_sig_1b = %8.2f\n", ttwjSig1bVal ) ;
       printf(" mu_ttwj_sb_1b  = %8.2f\n", ttwjSb1bVal ) ;



       TObject* qcd_sig_1b_obj = ws->obj("mu_qcd_sig_1b") ;
       TObject* qcd_sb_1b_obj = ws->obj("mu_qcd_sb_1b") ;

       if ( qcd_sig_1b_obj->IsA()->InheritsFrom(RooRealVar::Class()) )    {
          qcdSig1bVal = ((RooRealVar*) qcd_sig_1b_obj) -> getVal() ;
          qcdSig1bErr = ((RooRealVar*) qcd_sig_1b_obj) -> getError() ;
          qcdSb1bVal = ((RooFormulaVar*) qcd_sb_1b_obj) -> getVal() ;
          printf(" mu_qcd_sig is a RooRealVar.\n" ) ; 
       } else if ( qcd_sig_1b_obj->IsA()->InheritsFrom(RooFormulaVar::Class()) ) {
          qcdSig1bVal = ((RooFormulaVar*) qcd_sig_1b_obj) -> getVal() ;
          qcdSig1bErr = -1. ;
          qcdSb1bVal = ((RooRealVar*) qcd_sb_1b_obj) -> getVal() ;
          printf(" mu_qcd_sig_1b is a RooFormulaVar\n" ) ;
       } else {
          printf("\n\n\n *** what kind of class is mu_qcd_sig_1b???\n\n\n\n") ;
          return ;
       }
       printf(" mu_qcd_sig_1b  = %8.2f\n", qcdSig1bVal ) ;
       printf(" mu_qcd_sb_1b   = %8.2f\n", qcdSb1bVal ) ;





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



      gStyle->SetPadBottomMargin(0.25) ;
      gStyle->SetPadRightMargin(0.20) ;


       int nbins(30) ;

       TH1F* hfitqual_data = new TH1F("hfitqual_data", "RA2b likelihood fit results, data", nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_susy = new TH1F("hfitqual_susy", "RA2b likelihood fit results, susy", nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_ttwj = new TH1F("hfitqual_ttwj", "RA2b likelihood fit results, ttwj", nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_qcd  = new TH1F("hfitqual_qcd" , "RA2b likelihood fit results, qcd" , nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_znn  = new TH1F("hfitqual_znn" , "RA2b likelihood fit results, znn" , nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_np   = new TH1F("hfitqual_np"  , "RA2b likelihood fit results, np"  , 1, 0., 1. ) ;

       hfitqual_ttwj  -> SetFillColor(kBlue-9) ;
       hfitqual_qcd   -> SetFillColor(2) ;
       hfitqual_znn   -> SetFillColor(kGreen-3) ;
       hfitqual_susy  -> SetFillColor(6) ;
       hfitqual_np    -> SetFillColor(kOrange+1) ;

       hfitqual_data->SetMarkerStyle(20) ;
       hfitqual_data->SetLineWidth(2) ;

       THStack* hfitqual_fit = new THStack( "hfitqual_fit", "RA2b likelihood fit results, fit" ) ;

       TAxis* xaxis = hfitqual_data->GetXaxis() ;


       char binLabel[1000] ;
       int  binIndex ;

       double dataVal(0.) ;
       double dataErr(0.) ;
       double ttwjVal(0.) ;
       double qcdVal(0.) ;
       double znnVal(0.) ;
       double susyVal(0.) ;
       double lhtotalVal(0.) ;

       double eff_sf_sig_1b(0.) ;
       double eff_sf_sb_1b(0.) ;
       double eff_sf_sig_sl_1b(0.) ;
       double eff_sf_sb_sl_1b(0.) ;
       double eff_sf_sig_ldp_1b(0.) ;
       double eff_sf_sb_ldp_1b(0.) ;

       double eff_sf_sig_2b(0.) ;
       double eff_sf_sb_2b(0.) ;
       double eff_sf_sig_sl_2b(0.) ;
       double eff_sf_sb_sl_2b(0.) ;
       double eff_sf_sig_ldp_2b(0.) ;
       double eff_sf_sb_ldp_2b(0.) ;

       double eff_sf_sig_3b(0.) ;
       double eff_sf_sb_3b(0.) ;
       double eff_sf_sig_sl_3b(0.) ;
       double eff_sf_sb_sl_3b(0.) ;
       double eff_sf_sig_ldp_3b(0.) ;
       double eff_sf_sb_ldp_3b(0.) ;

       double sf_mc(0.) ;

       double znnSigVal ;
       double znnSigErr ;
       double susySigVal ;
       double susySigErr ;


     //-- SIG ---------------------------------------------------------------------------

       //-- 1b

       sprintf( binLabel, "SIG - 1b" ) ;
       binIndex = 2 ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_1b ;
       eff_sf_sig_1b = ((RooFormulaVar*) ws->obj("eff_sf_sig_1b")) -> getVal() ;
       susyVal = eff_sf_sig_1b * ( ((RooRealVar*) ws->obj("mu_susy_sig_1b")) -> getVal() ) ;
       ttwjVal = ttwjSig1bVal ;
       qcdVal  = qcdSig1bVal ;
       znnVal  = ((RooRealVar*) ws->obj("mu_znn_sig_1b"))  -> getVal() ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;

    //-- get values for later.
       susySigVal = ((RooRealVar*) ws->obj("mu_susy_sig_1b")) -> getVal() ;
       susySigErr = ((RooRealVar*) ws->obj("mu_susy_sig_1b")) -> getError() ;
       znnSigVal = ((RooRealVar*) ws->obj("mu_znn_sig_1b"))  -> getVal() ;
       znnSigErr = ((RooRealVar*) ws->obj("mu_znn_sig_1b"))  -> getError() ;

       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;



     //-- SIG ---------------------------------------------------------------------------

       //-- 2b

       sprintf( binLabel, "SIG - 2b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_2b ;
       eff_sf_sig_2b = ((RooFormulaVar*) ws->obj("eff_sf_sig_2b")) -> getVal() ;
       susyVal = eff_sf_sig_2b * ( ((RooRealVar*) ws->obj("mu_susy_sig_2b")) -> getVal() ) ;
       ttwjVal = ((RooFormulaVar*) ws->obj("mu_ttwj_sig_2b")) -> getVal()  ;
       qcdVal  = ((RooFormulaVar*) ws->obj("mu_qcd_sig_2b")) -> getVal()  ;
       znnVal  = ((RooRealVar*) ws->obj("mu_znn_sig_2b"))  -> getVal() ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;

       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;



     //-- SIG ---------------------------------------------------------------------------

       //-- 3b

       sprintf( binLabel, "SIG - 3b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_3b ;
       eff_sf_sig_3b = ((RooFormulaVar*) ws->obj("eff_sf_sig_3b")) -> getVal() ;
       susyVal = eff_sf_sig_3b * ( ((RooRealVar*) ws->obj("mu_susy_sig_3b")) -> getVal() ) ;
       ttwjVal = ((RooFormulaVar*) ws->obj("mu_ttwj_sig_3b")) -> getVal()  ;
       qcdVal  = ((RooFormulaVar*) ws->obj("mu_qcd_sig_3b")) -> getVal()  ;
       znnVal  = ((RooRealVar*) ws->obj("mu_znn_sig_3b"))  -> getVal() ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;

       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;



       binIndex++ ;

     //-- SB ---------------------------------------------------------------------------

       //-- 1b

       sprintf( binLabel, "SB - 1b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_1b ;
       eff_sf_sb_1b = ((RooFormulaVar*) ws->obj("eff_sf_sb_1b")) -> getVal() ;
       susyVal = eff_sf_sb_1b * ( ((RooRealVar*) ws->obj("mu_susy_sb_1b")) -> getVal() ) ;
       ttwjVal = ttwjSb1bVal ;
       qcdVal  = qcdSb1bVal ;
       znnVal  = ((RooRealVar*) ws->obj("mu_znn_sb_1b"))  -> getVal() ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;



     //-- SB ---------------------------------------------------------------------------

       //-- 2b

       sprintf( binLabel, "SB - 2b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_2b ;
       eff_sf_sb_2b = ((RooFormulaVar*) ws->obj("eff_sf_sb_2b")) -> getVal() ;
       susyVal = eff_sf_sb_2b * ( ((RooRealVar*) ws->obj("mu_susy_sb_2b")) -> getVal() ) ;
       ttwjVal = ((RooRealVar*) ws->obj("mu_ttwj_sb_2b")) -> getVal() ;
       qcdVal  = ((RooRealVar*) ws->obj("mu_qcd_sb_2b")) -> getVal() ;
       znnVal  = ((RooRealVar*) ws->obj("mu_znn_sb_2b"))  -> getVal() ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;



     //-- SB ---------------------------------------------------------------------------

       //-- 3b

       sprintf( binLabel, "SB - 3b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_3b ;
       eff_sf_sb_3b = ((RooFormulaVar*) ws->obj("eff_sf_sb_3b")) -> getVal() ;
       susyVal = eff_sf_sb_3b * ( ((RooRealVar*) ws->obj("mu_susy_sb_3b")) -> getVal() ) ;
       ttwjVal = ((RooRealVar*) ws->obj("mu_ttwj_sb_3b")) -> getVal() ;
       qcdVal  = ((RooRealVar*) ws->obj("mu_qcd_sb_3b")) -> getVal() ;
       znnVal  = ((RooRealVar*) ws->obj("mu_znn_sb_3b"))  -> getVal() ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;


       binIndex++ ;

     //-- SIG, SL ---------------------------------------------------------------------------

       //-- 1b
       //
       sprintf( binLabel, "SIG,SL - 1b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_sl_1b ;
       eff_sf_sig_sl_1b = ((RooFormulaVar*) ws->obj("eff_sf_sig_sl_1b")) -> getVal() ;
       susyVal = eff_sf_sig_sl_1b * ( ((RooRealVar*) ws->obj("mu_susy_sig_sl_1b")) -> getVal() ) ;
       ttwjVal = ((RooRealVar*) ws->obj("mu_ttwj_sig_sl_1b")) -> getVal() ;
       qcdVal  = 0. ;
       znnVal  = 0. ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;



     //-- SIG, SL ---------------------------------------------------------------------------

       //-- 2b
       //
       sprintf( binLabel, "SIG,SL - 2b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_sl_2b ;
       eff_sf_sig_sl_2b = ((RooFormulaVar*) ws->obj("eff_sf_sig_sl_2b")) -> getVal() ;
       susyVal = eff_sf_sig_sl_2b * ( ((RooRealVar*) ws->obj("mu_susy_sig_sl_2b")) -> getVal() ) ;
       ttwjVal = ((RooRealVar*) ws->obj("mu_ttwj_sig_sl_2b")) -> getVal() ;
       qcdVal  = 0. ;
       znnVal  = 0. ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;



     //-- SIG, SL ---------------------------------------------------------------------------

       //-- 3b
       //
       sprintf( binLabel, "SIG,SL - 3b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_sl_3b ;
       eff_sf_sig_sl_3b = ((RooFormulaVar*) ws->obj("eff_sf_sig_sl_3b")) -> getVal() ;
       susyVal = eff_sf_sig_sl_3b * ( ((RooRealVar*) ws->obj("mu_susy_sig_sl_3b")) -> getVal() ) ;
       ttwjVal = ((RooRealVar*) ws->obj("mu_ttwj_sig_sl_3b")) -> getVal() ;
       qcdVal  = 0. ;
       znnVal  = 0. ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;



       binIndex++ ;

     //-- SB, SL ---------------------------------------------------------------------------

       //--- 1b
       //
       sprintf( binLabel, "SB,SL - 1b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_sl_1b ;
       eff_sf_sb_sl_1b = ((RooFormulaVar*) ws->obj("eff_sf_sb_sl_1b")) -> getVal() ;
       susyVal = eff_sf_sb_sl_1b * ( ((RooRealVar*) ws->obj("mu_susy_sb_sl_1b")) -> getVal() ) ;
       ttwjVal = ((RooRealVar*) ws->obj("mu_ttwj_sb_sl_1b")) -> getVal() ;
       qcdVal  = 0. ;
       znnVal  = 0. ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;




     //-- SB, SL ---------------------------------------------------------------------------

       //--- 2b
       //
       sprintf( binLabel, "SB,SL - 2b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_sl_2b ;
       eff_sf_sb_sl_2b = ((RooFormulaVar*) ws->obj("eff_sf_sb_sl_2b")) -> getVal() ;
       susyVal = eff_sf_sb_sl_2b * ( ((RooRealVar*) ws->obj("mu_susy_sb_sl_2b")) -> getVal() ) ;
       ttwjVal = ((RooRealVar*) ws->obj("mu_ttwj_sb_sl_2b")) -> getVal() ;
       qcdVal  = 0. ;
       znnVal  = 0. ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;




     //-- SB, SL ---------------------------------------------------------------------------

       //--- 3b
       //
       sprintf( binLabel, "SB,SL - 3b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_sl_3b ;
       eff_sf_sb_sl_3b = ((RooFormulaVar*) ws->obj("eff_sf_sb_sl_3b")) -> getVal() ;
       susyVal = eff_sf_sb_sl_3b * ( ((RooRealVar*) ws->obj("mu_susy_sb_sl_3b")) -> getVal() ) ;
       ttwjVal = ((RooRealVar*) ws->obj("mu_ttwj_sb_sl_3b")) -> getVal() ;
       qcdVal  = 0. ;
       znnVal  = 0. ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;





       binIndex++ ;

     //-- SIG, LDP ---------------------------------------------------------------------------

       //-- 1b
       //
       sprintf( binLabel, "SIG,LDP - 1b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_ldp_1b ;
       eff_sf_sig_ldp_1b = ((RooFormulaVar*) ws->obj("eff_sf_sig_ldp_1b")) -> getVal() ;
       susyVal = eff_sf_sig_ldp_1b * ( ((RooRealVar*) ws->obj("mu_susy_sig_ldp_1b")) -> getVal() ) ;
       sf_mc = ((RooFormulaVar*) ws->obj("sf_mc")) -> getVal() ;
       ttwjVal = eff_sf_sig_ldp_1b * sf_mc * ( ( ((RooRealVar*) ws->obj("mu_ttbarsingletopzjetsmc_sig_ldp_1b")) -> getVal() )  + ( ((RooRealVar*) ws->obj("mu_WJmc_sig_ldp_1b")) -> getVal() )  ) ;
       qcdVal  = ((RooRealVar*) ws->obj("mu_qcd_sig_ldp_1b")) -> getVal() ;
       znnVal  = eff_sf_sig_ldp_1b * sf_mc * ( ((RooRealVar*) ws->obj("mu_Znnmc_sig_ldp_1b")) -> getVal()   ) ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;



     //-- SIG, LDP ---------------------------------------------------------------------------

       //-- 2b
       //
       sprintf( binLabel, "SIG,LDP - 2b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_ldp_2b ;
       eff_sf_sig_ldp_2b = ((RooFormulaVar*) ws->obj("eff_sf_sig_ldp_2b")) -> getVal() ;
       susyVal = eff_sf_sig_ldp_2b * ( ((RooRealVar*) ws->obj("mu_susy_sig_ldp_2b")) -> getVal() ) ;
       sf_mc = ((RooFormulaVar*) ws->obj("sf_mc")) -> getVal() ;
       ttwjVal = eff_sf_sig_ldp_2b * sf_mc * ( ( ((RooRealVar*) ws->obj("mu_ttbarsingletopzjetsmc_sig_ldp_2b")) -> getVal() )  + ( ((RooRealVar*) ws->obj("mu_WJmc_sig_ldp_2b")) -> getVal() )  ) ;
       qcdVal  = ((RooRealVar*) ws->obj("mu_qcd_sig_ldp_2b")) -> getVal() ;
       znnVal  = eff_sf_sig_ldp_2b * sf_mc * ( ((RooRealVar*) ws->obj("mu_Znnmc_sig_ldp_2b")) -> getVal()   ) ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;



     //-- SIG, LDP ---------------------------------------------------------------------------

       //-- 3b
       //
       sprintf( binLabel, "SIG,LDP - 3b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_ldp_3b ;
       eff_sf_sig_ldp_3b = ((RooFormulaVar*) ws->obj("eff_sf_sig_ldp_3b")) -> getVal() ;
       susyVal = eff_sf_sig_ldp_3b * ( ((RooRealVar*) ws->obj("mu_susy_sig_ldp_3b")) -> getVal() ) ;
       sf_mc = ((RooFormulaVar*) ws->obj("sf_mc")) -> getVal() ;
       ttwjVal = eff_sf_sig_ldp_3b * sf_mc * ( ( ((RooRealVar*) ws->obj("mu_ttbarsingletopzjetsmc_sig_ldp_3b")) -> getVal() )  + ( ((RooRealVar*) ws->obj("mu_WJmc_sig_ldp_3b")) -> getVal() )  ) ;
       qcdVal  = ((RooRealVar*) ws->obj("mu_qcd_sig_ldp_3b")) -> getVal() ;
       znnVal  = eff_sf_sig_ldp_3b * sf_mc * ( ((RooRealVar*) ws->obj("mu_Znnmc_sig_ldp_3b")) -> getVal()   ) ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;




       binIndex++ ;

     //-- SB, LDP ---------------------------------------------------------------------------

       //-- 1b
       sprintf( binLabel, "SB,LDP - 1b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_ldp_1b ;
       eff_sf_sb_ldp_1b = ((RooFormulaVar*) ws->obj("eff_sf_sb_ldp_1b")) -> getVal() ;
       susyVal = eff_sf_sb_ldp_1b * ( ((RooRealVar*) ws->obj("mu_susy_sb_ldp_1b")) -> getVal() ) ;
       sf_mc = ((RooFormulaVar*) ws->obj("sf_mc")) -> getVal() ;
       ttwjVal = eff_sf_sb_ldp_1b * sf_mc * ( ( ((RooRealVar*) ws->obj("mu_ttbarsingletopzjetsmc_sb_ldp_1b")) -> getVal() )  + ( ((RooRealVar*) ws->obj("mu_WJmc_sb_ldp_1b")) -> getVal() )  ) ;
       qcdVal  = ((RooRealVar*) ws->obj("mu_qcd_sb_ldp_1b")) -> getVal() ;
       znnVal  = eff_sf_sb_ldp_1b * sf_mc * ( ((RooRealVar*) ws->obj("mu_Znnmc_sb_ldp_1b")) -> getVal()   ) ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;





     //-- SB, LDP ---------------------------------------------------------------------------

       //-- 2b
       sprintf( binLabel, "SB,LDP - 2b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_ldp_2b ;
       eff_sf_sb_ldp_2b = ((RooFormulaVar*) ws->obj("eff_sf_sb_ldp_2b")) -> getVal() ;
       susyVal = eff_sf_sb_ldp_2b * ( ((RooRealVar*) ws->obj("mu_susy_sb_ldp_2b")) -> getVal() ) ;
       sf_mc = ((RooFormulaVar*) ws->obj("sf_mc")) -> getVal() ;
       ttwjVal = eff_sf_sb_ldp_2b * sf_mc * ( ( ((RooRealVar*) ws->obj("mu_ttbarsingletopzjetsmc_sb_ldp_2b")) -> getVal() )  + ( ((RooRealVar*) ws->obj("mu_WJmc_sb_ldp_2b")) -> getVal() )  ) ;
       qcdVal  = ((RooRealVar*) ws->obj("mu_qcd_sb_ldp_2b")) -> getVal() ;
       znnVal  = eff_sf_sb_ldp_2b * sf_mc * ( ((RooRealVar*) ws->obj("mu_Znnmc_sb_ldp_2b")) -> getVal()   ) ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;





     //-- SB, LDP ---------------------------------------------------------------------------

       //-- 3b
       sprintf( binLabel, "SB,LDP - 3b" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_ldp_3b ;
       eff_sf_sb_ldp_3b = ((RooFormulaVar*) ws->obj("eff_sf_sb_ldp_3b")) -> getVal() ;
       susyVal = eff_sf_sb_ldp_3b * ( ((RooRealVar*) ws->obj("mu_susy_sb_ldp_3b")) -> getVal() ) ;
       sf_mc = ((RooFormulaVar*) ws->obj("sf_mc")) -> getVal() ;
       ttwjVal = eff_sf_sb_ldp_3b * sf_mc * ( ( ((RooRealVar*) ws->obj("mu_ttbarsingletopzjetsmc_sb_ldp_3b")) -> getVal() )  + ( ((RooRealVar*) ws->obj("mu_WJmc_sb_ldp_3b")) -> getVal() )  ) ;
       qcdVal  = ((RooRealVar*) ws->obj("mu_qcd_sb_ldp_3b")) -> getVal() ;
       znnVal  = eff_sf_sb_ldp_3b * sf_mc * ( ((RooRealVar*) ws->obj("mu_Znnmc_sb_ldp_3b")) -> getVal()   ) ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      susy : %8.2f\n", binLabel, susyVal ) ;
       printf(" %8s      ttwj : %8.2f\n", binLabel, ttwjVal ) ;
       printf(" %8s      qcd  : %8.2f\n", binLabel, qcdVal ) ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          susyVal = susyVal / dataVal ;
          ttwjVal = ttwjVal / dataVal ;
          qcdVal  = qcdVal / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_susy -> SetBinContent( binIndex, susyVal ) ;
       hfitqual_ttwj -> SetBinContent( binIndex, ttwjVal ) ;
       hfitqual_qcd  -> SetBinContent( binIndex, qcdVal ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;










       binIndex++ ;

     //-- SIG, ee ---------------------------------------------------------------------------

       sprintf( binLabel, "SIG,ee" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_ee ;
       znnVal  = ( ((RooRealVar*) ws->obj("mu_zee_sig")) -> getVal() ) / ( ((RooRealVar*) ws->obj("fsig_ee"))->getVal() ) ;
       lhtotalVal = znnVal  ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;



     //-- SB, ee ---------------------------------------------------------------------------

       sprintf( binLabel, "SB,ee" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_ee ;
       znnVal  = ( ((RooRealVar*) ws->obj("mu_zee_sb")) -> getVal() ) / ( ((RooRealVar*) ws->obj("fsig_ee"))->getVal() ) ;
       lhtotalVal = znnVal  ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;




     //-- SIG, mm ---------------------------------------------------------------------------

       sprintf( binLabel, "SIG,mm" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_mm ;
       znnVal  = ( ((RooRealVar*) ws->obj("mu_zmm_sig")) -> getVal() ) / ( ((RooRealVar*) ws->obj("fsig_mm"))->getVal() ) ;
       lhtotalVal = znnVal  ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;



     //-- SB, mm ---------------------------------------------------------------------------

       sprintf( binLabel, "SB,mm" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_mm ;
       znnVal  = ( ((RooRealVar*) ws->obj("mu_zmm_sb")) -> getVal() ) / ( ((RooRealVar*) ws->obj("fsig_mm"))->getVal() ) ;
       lhtotalVal = znnVal  ;

       dataErr = sqrt(dataVal) ;


       printf("\n\n") ;
       printf(" %8s      znn  : %8.2f\n", binLabel, znnVal ) ;
       printf(" %8s  LH total : %8.2f\n", binLabel, lhtotalVal ) ;
       printf(" %8s      data : %5.0f\n"  , binLabel, dataVal ) ;


       if ( doNorm && dataVal > 0. ) {
          dataErr = dataErr / dataVal ;
          znnVal  = znnVal / dataVal ;
          dataVal = 1. ;
       }

       hfitqual_data -> SetBinContent( binIndex, dataVal ) ;
       hfitqual_data -> SetBinError( binIndex, dataErr ) ;
       hfitqual_znn  -> SetBinContent( binIndex, znnVal ) ;











//   //-- Eff sf ---------

//    float eff_sf_prim   =  ((RooRealVar*) ws->obj("eff_sf_prim")) -> getVal() ;

//    hfitqual_np -> SetBinContent( 1, eff_sf_prim ) ;










      printf("\n\n\n") ;


     //--- final formatting and drawing.


      gStyle->SetOptStat(0) ;
      gStyle->SetOptTitle(0) ;

      hfitqual_fit->Add( hfitqual_znn ) ;
      hfitqual_fit->Add( hfitqual_qcd ) ;
      hfitqual_fit->Add( hfitqual_ttwj ) ;
      hfitqual_fit->Add( hfitqual_susy ) ;

      TLegend* legend = new TLegend(0.85,0.2,0.97,0.9) ;

      legend->AddEntry( hfitqual_data,  "data" ) ;
      legend->AddEntry( hfitqual_susy,  "SUSY" ) ;
      legend->AddEntry( hfitqual_ttwj, "ttwj" ) ;
      legend->AddEntry( hfitqual_qcd,   "QCD" ) ;
      legend->AddEntry( hfitqual_znn,   "Znunu" ) ;
    //legend->AddEntry( hfitqual_np,  "Eff PG" ) ;


      if ( doNorm ) {
         hfitqual_data->SetMaximum( hmax ) ;
      } else {
         hfitqual_data->SetMaximum( hmax*(hfitqual_data->GetMaximum()) ) ;
      }

      TCanvas* cfitqual = new TCanvas("cfitqual","RA2b fit quality", 1000, 550 ) ;

      gPad->SetTicks(1,0) ;


      hfitqual_data->SetLabelSize(0.055,"x") ;
      hfitqual_data->GetXaxis()->LabelsOption("v") ;

      hfitqual_data->Draw("histpe") ;
      hfitqual_fit->Draw("same") ;
      hfitqual_data->Draw("same") ;
      gPad->SetGridy(1) ;
      legend->Draw() ;

      cfitqual->Update() ;




// //-- Give SIG box values and uncertainties.

      TText* fittext = new TText() ;
      fittext->SetTextSize(0.045) ;
//    double fitval ;
//    double fiterr ;
//    double tx = 0.80 ;
//    double ty = 0.71 ;
//    double dy = 0.115 ;
//    char fitvalchars[1000] ;



//    if ( mu_susy_sig_1b_val < 0. ) {
//       fitval = susySigVal ;
//       fiterr = susySigErr ;
//       if ( fitval<10.) {
//          sprintf( fitvalchars, "%3.1f +/- %3.1f", fitval, fiterr ) ;
//       } else {
//          sprintf( fitvalchars, "%3.0f +/- %3.0f", fitval, fiterr ) ;
//       }
//       fittext->DrawTextNDC( tx, ty, fitvalchars ) ;
//    } else {
//       fitval = susySigVal ;
//       if ( fitval<10.) {
//          sprintf( fitvalchars, "%3.1f", fitval ) ;
//       } else {
//          sprintf( fitvalchars, "%3.0f", fitval ) ;
//       }
//       fittext->DrawTextNDC( tx, ty, fitvalchars ) ;
//    }



//    if ( ttwjSigErr > 0. ) {
//       fitval = ttwjSigVal ;
//       fiterr = ttwjSigErr ;
//       if ( fitval<10.) {
//          sprintf( fitvalchars, "%3.1f +/- %3.1f", fitval, fiterr ) ;
//       } else {
//          sprintf( fitvalchars, "%3.0f +/- %3.0f", fitval, fiterr ) ;
//       }
//    } else {
//       fitval = ttwjSigVal ;
//       sprintf( fitvalchars, "%3.0f", fitval ) ;
//       if ( fitval<10.) {
//          sprintf( fitvalchars, "%3.1f", fitval ) ;
//       } else {
//          sprintf( fitvalchars, "%3.0f", fitval ) ;
//       }
//    }
//    ty = ty - dy ;
//    fittext->DrawTextNDC( tx, ty, fitvalchars ) ;




//    if ( qcdSigErr > 0. ) {
//       fitval = qcdSigVal ;
//       fiterr = qcdSigErr ;
//       if ( fitval<10.) {
//          sprintf( fitvalchars, "%3.1f +/- %3.1f", fitval, fiterr ) ;
//       } else {
//          sprintf( fitvalchars, "%3.0f +/- %3.0f", fitval, fiterr ) ;
//       }
//    } else {
//       fitval = qcdSigVal ;
//       if ( fitval<10.) {
//          sprintf( fitvalchars, "%3.1f", fitval ) ;
//       } else {
//          sprintf( fitvalchars, "%3.0f", fitval ) ;
//       }
//    }
//    ty = ty - dy ;
//    fittext->DrawTextNDC( tx, ty, fitvalchars ) ;




//    fitval = znnSigVal ;
//    fiterr = znnSigErr ;
//    if ( fitval<10.) {
//       sprintf( fitvalchars, "%3.1f +/- %3.1f", fitval, fiterr ) ;
//    } else {
//       sprintf( fitvalchars, "%3.0f +/- %3.0f", fitval, fiterr ) ;
//    }
//    ty = ty - dy ;
//    fittext->DrawTextNDC( tx, ty, fitvalchars ) ;




//    fitval =  eff_sf_prim ;
//    sprintf( fitvalchars, "%4.2f", fitval ) ;
//    ty = ty - dy ;
//    fittext->DrawTextNDC( tx, ty, fitvalchars ) ;




      TPad* bigPad = (TPad*) gPad ;


// //--- Efficiency scale factor primary Gaussian value.

//    hfitqual_np->SetMinimum(-5.) ;
//    hfitqual_np->SetMaximum( 5.) ;

//    hfitqual_np->SetNdivisions(101,"x") ;
//    hfitqual_np->SetNdivisions(101,"y") ;
//    hfitqual_np->SetLabelOffset(99,"y") ;



//    TPad* tp = new TPad("tp","tp",0.52,0.,0.61,1.0) ;

//    tp->SetRightMargin(0.4) ;



//    tp->Draw() ;
//    tp->cd() ;
//    hfitqual_np->SetLabelSize(0.5,"x") ;
//    xaxis = hfitqual_np->GetXaxis() ;
//    xaxis->SetBinLabel(1,"Eff PG") ;
//    hfitqual_np->GetXaxis()->LabelsOption("v") ;
//    hfitqual_np->Draw() ;

//    cfitqual->Update() ;


//    TGaxis* axis = new TGaxis() ;
//    axis->SetLabelOffset(0.1) ;
//    axis->SetLabelSize(0.30) ;
//    axis->SetTickSize(0.2) ;
//    axis->DrawAxis( 1.0, -5., 1.0, 5., -5., 5., 510, "+LS") ;

//    cfitqual->Update() ;



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

      cfitqual->SaveAs("fitqual.png") ;


   }






