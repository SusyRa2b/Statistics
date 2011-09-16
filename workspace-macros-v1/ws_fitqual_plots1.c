
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

   void ws_fitqual_plots1( const char* wsfile = "ws-lm9-ge1btight.root", double mu_susy_sig_val = 0., bool doNorm = false ) {

      double hmax = 1.5 ;

       char sel[100] ;
       if ( strstr( wsfile, "ge1bloose" ) != 0 ) {
          sprintf( sel, "ge1bloose" ) ;
       } else if ( strstr( wsfile, "ge1btight" ) != 0 ) {
          sprintf( sel, "ge1btight" ) ;
       } else if ( strstr( wsfile, "ge2bloose" ) != 0 ) {
          sprintf( sel, "ge2bloose" ) ;
       } else if ( strstr( wsfile, "ge2btight" ) != 0 ) {
          sprintf( sel, "ge2btight" ) ;
       } else {
          printf("\n\n\n *** can't figure out which selection this is.  I quit.\n\n" ) ;
          return ;
       }
       printf("\n\n selection is %s\n\n", sel ) ;




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

       RooRealVar* rrv_mu_susy_sig = ws->var("mu_susy_sig") ;
       if ( rrv_mu_susy_sig == 0x0 ) {
          printf("\n\n\n *** can't find mu_susy_sig in workspace.  Quitting.\n\n\n") ;
          return ;
       } else {
          printf(" current value is : %8.3f\n", rrv_mu_susy_sig->getVal() ) ; cout << flush ;
          rrv_mu_susy_sig->setConstant(kFALSE) ;
       }


       printf("\n\n\n  ===== Doing a fit with SUSY component floating ====================\n\n") ;

       RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true) ) ;
       double logLikelihoodSusyFloat = fitResult->minNll() ;

       double logLikelihoodSusyFixed(0.) ;
       double testStatVal(-1.) ;
       if ( mu_susy_sig_val >= 0. ) {
          printf("\n\n\n  ===== Doing a fit with SUSY fixed ====================\n\n") ;
          printf(" fixing mu_susy_sig to %8.2f.\n", mu_susy_sig_val ) ;
          rrv_mu_susy_sig->setVal( mu_susy_sig_val ) ;
          rrv_mu_susy_sig->setConstant(kTRUE) ;

          fitResult = likelihood->fitTo( *rds, Save(true) ) ;
          logLikelihoodSusyFixed = fitResult->minNll() ;
          testStatVal = 2.*(logLikelihoodSusyFixed - logLikelihoodSusyFloat) ;
          printf("\n\n\n ======= test statistic : -2 * ln (L_fixed / ln L_max) = %8.3f\n\n\n", testStatVal ) ;
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








       double ttwjSigVal, ttwjSbVal ;
       double qcdSigVal,  qcdSbVal  ;

       double ttwjSigErr ;
       double qcdSigErr ;


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

       TObject* ttwj_sig_obj = ws->obj("mu_ttwj_sig") ;
       TObject* ttwj_sb_obj = ws->obj("mu_ttwj_sb") ;

       if ( ttwj_sig_obj->IsA()->InheritsFrom(RooRealVar::Class()) )    {
          ttwjSigVal = ((RooRealVar*) ttwj_sig_obj) -> getVal() ;
          ttwjSigErr = ((RooRealVar*) ttwj_sig_obj) -> getError() ;
          ttwjSbVal = ((RooFormulaVar*) ttwj_sb_obj) -> getVal() ;
          printf(" mu_ttwj_sig is a RooRealVar.\n" ) ; 
       } else if ( ttwj_sig_obj->IsA()->InheritsFrom(RooFormulaVar::Class()) ) {
          ttwjSigVal = ((RooFormulaVar*) ttwj_sig_obj) -> getVal() ;
          ttwjSigErr = -1. ;
          ttwjSbVal = ((RooRealVar*) ttwj_sb_obj) -> getVal() ;
          printf(" mu_ttwj_sig is a RooFormulaVar\n" ) ;
       } else {
          printf("\n\n\n *** what kind of class is mu_ttwj_sig???\n\n\n\n") ;
          return ;
       }
       printf(" mu_ttwj_sig = %8.2f\n", ttwjSigVal ) ;
       printf(" mu_ttwj_sb  = %8.2f\n", ttwjSbVal ) ;



       TObject* qcd_sig_obj = ws->obj("mu_qcd_sig") ;
       TObject* qcd_sb_obj = ws->obj("mu_qcd_sb") ;

       if ( qcd_sig_obj->IsA()->InheritsFrom(RooRealVar::Class()) )    {
          qcdSigVal = ((RooRealVar*) qcd_sig_obj) -> getVal() ;
          qcdSigErr = ((RooRealVar*) qcd_sig_obj) -> getError() ;
          qcdSbVal = ((RooFormulaVar*) qcd_sb_obj) -> getVal() ;
          printf(" mu_qcd_sig is a RooRealVar.\n" ) ; 
       } else if ( qcd_sig_obj->IsA()->InheritsFrom(RooFormulaVar::Class()) ) {
          qcdSigVal = ((RooFormulaVar*) qcd_sig_obj) -> getVal() ;
          qcdSigErr = -1. ;
          qcdSbVal = ((RooRealVar*) qcd_sb_obj) -> getVal() ;
          printf(" mu_qcd_sig is a RooFormulaVar\n" ) ;
       } else {
          printf("\n\n\n *** what kind of class is mu_qcd_sig???\n\n\n\n") ;
          return ;
       }
       printf(" mu_qcd_sig  = %8.2f\n", qcdSigVal ) ;
       printf(" mu_qcd_sb   = %8.2f\n", qcdSbVal ) ;





      //--- unpack observables.

       int dataNsig(0) ;
       int dataNsb(0) ;
       int dataNsig_sl(0) ;
       int dataNsb_sl(0) ;
       int dataNsig_ldp(0) ;
       int dataNsb_ldp(0) ;

       const RooArgSet* dsras = rds->get() ;
       TIterator* obsIter = dsras->createIterator() ;
       while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {
          if ( strcmp( obs->GetName(), "Nsig"     ) == 0 ) { dataNsig     = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb"      ) == 0 ) { dataNsb      = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsig_sl"  ) == 0 ) { dataNsig_sl  = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_sl"   ) == 0 ) { dataNsb_sl   = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsig_ldp" ) == 0 ) { dataNsig_ldp = obs->getVal() ; }
          if ( strcmp( obs->GetName(), "Nsb_ldp"  ) == 0 ) { dataNsb_ldp  = obs->getVal() ; }
       }



      gStyle->SetPadBottomMargin(0.20) ;
      gStyle->SetPadRightMargin(0.50) ;


       int nbins(8) ;

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

       double eff_sf_sig(0.) ;
       double eff_sf_sb(0.) ;
       double eff_sf_sig_sl(0.) ;
       double eff_sf_sb_sl(0.) ;
       double eff_sf_sig_ldp(0.) ;
       double eff_sf_sb_ldp(0.) ;

       double sf_mc(0.) ;


     //-- SIG ---------------------------------------------------------------------------

       double znnSigVal ;
       double znnSigErr ;
       double susySigVal ;
       double susySigErr ;

       sprintf( binLabel, "SIG" ) ;
       binIndex = 2 ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig ;
       eff_sf_sig = ((RooFormulaVar*) ws->obj("eff_sf_sig")) -> getVal() ;
       susyVal = eff_sf_sig * ( ((RooRealVar*) ws->obj("mu_susy_sig")) -> getVal() ) ;
       ttwjVal = ttwjSigVal ;
       qcdVal  = qcdSigVal ;
       znnVal  = ((RooRealVar*) ws->obj("mu_znn_sig"))  -> getVal() ;
       lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

       dataErr = sqrt(dataVal) ;

    //-- get values for later.
       susySigVal = ((RooRealVar*) ws->obj("mu_susy_sig")) -> getVal() ;
       susySigErr = ((RooRealVar*) ws->obj("mu_susy_sig")) -> getError() ;
       znnSigVal = ((RooRealVar*) ws->obj("mu_znn_sig"))  -> getVal() ;
       znnSigErr = ((RooRealVar*) ws->obj("mu_znn_sig"))  -> getError() ;

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

       sprintf( binLabel, "SB" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb ;
       eff_sf_sb = ((RooFormulaVar*) ws->obj("eff_sf_sb")) -> getVal() ;
       susyVal = eff_sf_sb * ( ((RooRealVar*) ws->obj("mu_susy_sb")) -> getVal() ) ;
       ttwjVal = ttwjSbVal ;
       qcdVal  = qcdSbVal ;
       znnVal  = ((RooRealVar*) ws->obj("mu_znn_sb"))  -> getVal() ;
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

       sprintf( binLabel, "SIG,SL" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_sl ;
       eff_sf_sig_sl = ((RooFormulaVar*) ws->obj("eff_sf_sig_sl")) -> getVal() ;
       susyVal = eff_sf_sig_sl * ( ((RooRealVar*) ws->obj("mu_susy_sig_sl")) -> getVal() ) ;
       ttwjVal = ((RooRealVar*) ws->obj("mu_ttwj_sig_sl")) -> getVal() ;
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

       sprintf( binLabel, "SB,SL" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_sl ;
       eff_sf_sb_sl = ((RooFormulaVar*) ws->obj("eff_sf_sb_sl")) -> getVal() ;
       susyVal = eff_sf_sb_sl * ( ((RooRealVar*) ws->obj("mu_susy_sb_sl")) -> getVal() ) ;
       ttwjVal = ((RooRealVar*) ws->obj("mu_ttwj_sb_sl")) -> getVal() ;
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






     //-- SIG, LDP ---------------------------------------------------------------------------

       sprintf( binLabel, "SIG,LDP" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsig_ldp ;
       eff_sf_sig_ldp = ((RooFormulaVar*) ws->obj("eff_sf_sig_ldp")) -> getVal() ;
       susyVal = eff_sf_sig_ldp * ( ((RooRealVar*) ws->obj("mu_susy_sig_ldp")) -> getVal() ) ;
       sf_mc = ((RooFormulaVar*) ws->obj("sf_mc")) -> getVal() ;
       ttwjVal = eff_sf_sig_ldp * sf_mc * ( ( ((RooRealVar*) ws->obj("mu_ttbarmc_sig_ldp")) -> getVal() )  + ( ((RooRealVar*) ws->obj("mu_WJmc_sig_ldp")) -> getVal() )  ) ;
       qcdVal  = ((RooRealVar*) ws->obj("mu_qcd_sig_ldp")) -> getVal() ;
       znnVal  = eff_sf_sig_ldp * sf_mc * ( ((RooRealVar*) ws->obj("mu_Znnmc_sig_ldp")) -> getVal()   ) ;
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

       sprintf( binLabel, "SB,LDP" ) ;
       binIndex++ ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;


       dataVal = dataNsb_ldp ;
       eff_sf_sb_ldp = ((RooFormulaVar*) ws->obj("eff_sf_sb_ldp")) -> getVal() ;
       susyVal = eff_sf_sb_ldp * ( ((RooRealVar*) ws->obj("mu_susy_sb_ldp")) -> getVal() ) ;
       sf_mc = ((RooFormulaVar*) ws->obj("sf_mc")) -> getVal() ;
       ttwjVal = eff_sf_sb_ldp * sf_mc * ( ( ((RooRealVar*) ws->obj("mu_ttbarmc_sb_ldp")) -> getVal() )  + ( ((RooRealVar*) ws->obj("mu_WJmc_sb_ldp")) -> getVal() )  ) ;
       qcdVal  = ((RooRealVar*) ws->obj("mu_qcd_sb_ldp")) -> getVal() ;
       znnVal  = eff_sf_sb_ldp * sf_mc * ( ((RooRealVar*) ws->obj("mu_Znnmc_sb_ldp")) -> getVal()   ) ;
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









     //-- Eff sf ---------

      float eff_sf_prim   =  ((RooRealVar*) ws->obj("eff_sf_prim")) -> getVal() ;

      hfitqual_np -> SetBinContent( 1, eff_sf_prim ) ;










      printf("\n\n\n") ;


     //--- final formatting and drawing.


      gStyle->SetOptStat(0) ;
      gStyle->SetOptTitle(0) ;

      hfitqual_fit->Add( hfitqual_znn ) ;
      hfitqual_fit->Add( hfitqual_qcd ) ;
      hfitqual_fit->Add( hfitqual_ttwj ) ;
      hfitqual_fit->Add( hfitqual_susy ) ;

      TLegend* legend = new TLegend(0.65,0.2,0.77,0.9) ;

      legend->AddEntry( hfitqual_data,  "data" ) ;
      legend->AddEntry( hfitqual_susy,  "SUSY" ) ;
      legend->AddEntry( hfitqual_ttwj, "ttwj" ) ;
      legend->AddEntry( hfitqual_qcd,   "QCD" ) ;
      legend->AddEntry( hfitqual_znn,   "Znunu" ) ;
      legend->AddEntry( hfitqual_np,  "Eff PG" ) ;


      if ( doNorm ) {
         hfitqual_data->SetMaximum( hmax ) ;
      } else {
         hfitqual_data->SetMaximum( hmax*(hfitqual_data->GetMaximum()) ) ;
      }

      TCanvas* cfitqual = new TCanvas("cfitqual","RA2b fit quality", 780, 550 ) ;

      gPad->SetTicks(1,0) ;


      hfitqual_data->SetLabelSize(0.055,"x") ;
      hfitqual_data->GetXaxis()->LabelsOption("v") ;

      hfitqual_data->Draw("histpe") ;
      hfitqual_fit->Draw("same") ;
      hfitqual_data->Draw("same") ;
      gPad->SetGridy(1) ;
      legend->Draw() ;

      cfitqual->Update() ;




   //-- Give SIG box values and uncertainties.

      TText* fittext = new TText() ;
      fittext->SetTextSize(0.045) ;
      double fitval ;
      double fiterr ;
      double tx = 0.80 ;
      double ty = 0.71 ;
      double dy = 0.115 ;
      char fitvalchars[1000] ;



      if ( mu_susy_sig_val < 0. ) {
         fitval = susySigVal ;
         fiterr = susySigErr ;
         if ( fitval<10.) {
            sprintf( fitvalchars, "%3.1f +/- %3.1f", fitval, fiterr ) ;
         } else {
            sprintf( fitvalchars, "%3.0f +/- %3.0f", fitval, fiterr ) ;
         }
         fittext->DrawTextNDC( tx, ty, fitvalchars ) ;
      } else {
         fitval = susySigVal ;
         if ( fitval<10.) {
            sprintf( fitvalchars, "%3.1f", fitval ) ;
         } else {
            sprintf( fitvalchars, "%3.0f", fitval ) ;
         }
         fittext->DrawTextNDC( tx, ty, fitvalchars ) ;
      }



      if ( ttwjSigErr > 0. ) {
         fitval = ttwjSigVal ;
         fiterr = ttwjSigErr ;
         if ( fitval<10.) {
            sprintf( fitvalchars, "%3.1f +/- %3.1f", fitval, fiterr ) ;
         } else {
            sprintf( fitvalchars, "%3.0f +/- %3.0f", fitval, fiterr ) ;
         }
      } else {
         fitval = ttwjSigVal ;
         sprintf( fitvalchars, "%3.0f", fitval ) ;
         if ( fitval<10.) {
            sprintf( fitvalchars, "%3.1f", fitval ) ;
         } else {
            sprintf( fitvalchars, "%3.0f", fitval ) ;
         }
      }
      ty = ty - dy ;
      fittext->DrawTextNDC( tx, ty, fitvalchars ) ;




      if ( qcdSigErr > 0. ) {
         fitval = qcdSigVal ;
         fiterr = qcdSigErr ;
         if ( fitval<10.) {
            sprintf( fitvalchars, "%3.1f +/- %3.1f", fitval, fiterr ) ;
         } else {
            sprintf( fitvalchars, "%3.0f +/- %3.0f", fitval, fiterr ) ;
         }
      } else {
         fitval = qcdSigVal ;
         if ( fitval<10.) {
            sprintf( fitvalchars, "%3.1f", fitval ) ;
         } else {
            sprintf( fitvalchars, "%3.0f", fitval ) ;
         }
      }
      ty = ty - dy ;
      fittext->DrawTextNDC( tx, ty, fitvalchars ) ;




      fitval = znnSigVal ;
      fiterr = znnSigErr ;
      if ( fitval<10.) {
         sprintf( fitvalchars, "%3.1f +/- %3.1f", fitval, fiterr ) ;
      } else {
         sprintf( fitvalchars, "%3.0f +/- %3.0f", fitval, fiterr ) ;
      }
      ty = ty - dy ;
      fittext->DrawTextNDC( tx, ty, fitvalchars ) ;




      fitval =  eff_sf_prim ;
      sprintf( fitvalchars, "%4.2f", fitval ) ;
      ty = ty - dy ;
      fittext->DrawTextNDC( tx, ty, fitvalchars ) ;




      TPad* bigPad = (TPad*) gPad ;


   //--- Efficiency scale factor primary Gaussian value.

      hfitqual_np->SetMinimum(-5.) ;
      hfitqual_np->SetMaximum( 5.) ;

      hfitqual_np->SetNdivisions(101,"x") ;
      hfitqual_np->SetNdivisions(101,"y") ;
      hfitqual_np->SetLabelOffset(99,"y") ;



      TPad* tp = new TPad("tp","tp",0.52,0.,0.61,1.0) ;

      tp->SetRightMargin(0.4) ;



      tp->Draw() ;
      tp->cd() ;
      hfitqual_np->SetLabelSize(0.5,"x") ;
      xaxis = hfitqual_np->GetXaxis() ;
      xaxis->SetBinLabel(1,"Eff PG") ;
      hfitqual_np->GetXaxis()->LabelsOption("v") ;
      hfitqual_np->Draw() ;

      cfitqual->Update() ;


      TGaxis* axis = new TGaxis() ;
      axis->SetLabelOffset(0.1) ;
      axis->SetLabelSize(0.30) ;
      axis->SetTickSize(0.2) ;
      axis->DrawAxis( 1.0, -5., 1.0, 5., -5., 5., 510, "+LS") ;

      cfitqual->Update() ;



   //--- title

      bigPad->cd() ;

      char titletext[10000] ;
      if ( mu_susy_sig_val < 0. ) {
         sprintf( titletext, "%s : SUSY floating", sel ) ;
      } else {
         sprintf( titletext, "%s : SUSY fixed to %.1f", sel, mu_susy_sig_val ) ;
      }
      fittext->SetTextSize(0.050) ;
      fittext->DrawTextNDC( 0.02, 0.93, titletext ) ;

      if ( testStatVal >=0 ) {
         sprintf( titletext, "-2 * ln (L_fixed / ln L_max) = %8.3f\n", testStatVal ) ;
         fittext->SetTextSize(0.035) ;
         fittext->DrawTextNDC( 0.12, 0.85, titletext ) ;
      }

      cfitqual->Update() ;


   }






