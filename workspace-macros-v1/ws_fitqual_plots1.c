
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
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

   void ws_fitqual_plots1( const char* wsfile = "ws-lm9-ge1btight.root", double mu_susy_sig_val = 0. ) {


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

       const RooArgSet* observables = modelConfig->GetObservables() ;

       const RooArgSet* nuisanceParameters = modelConfig->GetNuisanceParameters() ;

       RooRealVar* rrv_mu_susy_sig = ws->var("mu_susy_sig") ;
       if ( rrv_mu_susy_sig == 0x0 ) {
          printf("\n\n\n *** can't find mu_susy_sig in workspace.  Quitting.\n\n\n") ;
          return ;
       } else {
          printf(" current value is : %8.3f\n", rrv_mu_susy_sig->getVal() ) ; cout << flush ;
          if ( mu_susy_sig_val < 0. ) {
             printf(" allowing mu_susy_sig to float in the fit.\n") ;
             rrv_mu_susy_sig->setConstant(kFALSE) ;
          } else {
             printf(" fixing mu_susy_sig to %8.2f.\n", mu_susy_sig_val ) ;
             rrv_mu_susy_sig->setVal( mu_susy_sig_val ) ;
             rrv_mu_susy_sig->setConstant(kTRUE) ;
          }
       }

       printf("\n\n\n  ===== Doing a fit ====================\n\n") ;

       RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true) ) ;

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

     //-- Note: Depending on the first 2 arguments used in ra2bRoostatsClass*.c, you either
     //          have the SIG or SB parameter (but not both) for ttwj and qcd.
     //
     //          To look at the sig vars, use ra2bRoostatsClass7( true, false )
     //          To look at the sb  vars, use ra2bRoostatsClass7( false, true )
     //
     //         These are set in the make_*_ws_*.c macros
     //         (e.g. make_lm9_ws_ge1btight.c).
     //


       TObject* ttwj_sig_obj = ws->obj("mu_ttwj_sig") ;
       TObject* ttwj_sb_obj  = ws->obj("mu_ttwj_sb") ;

       if ( ttwj_sig_obj->IsA()->InheritsFrom(RooRealVar::Class()) )    { printf(" mu_ttwj_sig is a RooRealVar\n" ) ; }
       if ( ttwj_sig_obj->IsA()->InheritsFrom(RooFormulaVar::Class()) ) { printf(" mu_ttwj_sig is a RooFormulaVar\n" ) ; }
       if ( ttwj_sb_obj->IsA()->InheritsFrom(RooRealVar::Class()) )    { printf(" mu_ttwj_sb is a RooRealVar\n" ) ; }
       if ( ttwj_sb_obj->IsA()->InheritsFrom(RooFormulaVar::Class()) ) { printf(" mu_ttwj_sb is a RooFormulaVar\n" ) ; }





       int nbins(10) ;

       TH1F* hfitqual_data = new TH1F("hfitqual_data", "RA2b likelihood fit results, data", nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_susy = new TH1F("hfitqual_susy", "RA2b likelihood fit results, susy", nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_ttwj = new TH1F("hfitqual_ttwj", "RA2b likelihood fit results, ttwj", nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_qcd  = new TH1F("hfitqual_qcd" , "RA2b likelihood fit results, qcd" , nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_znn  = new TH1F("hfitqual_znn" , "RA2b likelihood fit results, znn" , nbins, 0.5, nbins+0.5 ) ;
       TH1F* hfitqual_np   = new TH1F("hfitqual_np"  , "RA2b likelihood fit results, np"  , nbins, 0.5, nbins+0.5 ) ;

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

       double dataVal ;
       double ttwjVal ;
       double qcdVal ;
       double znnVal ;
       double susyVal ;

       double eff_sf_sig ;
       double eff_sf_sb ;
       double eff_sf_sig_sl ;
       double eff_sf_sb_sl ;
       double eff_sf_sig_ldp ;
       double eff_sf_sb_ldp ;

       double sf_mc ;



     //-- SIG ---------------------------------------------------------------------------

       sprintf( binLabel, "SIG" ) ;
       binIndex = 2 ;
       xaxis->SetBinLabel(binIndex, binLabel ) ;




//     char parName[1000] ;

//     RooRealVar* rrv_ttwj(0x0) ;
//     RooRealVar* rrv_qcd(0x0) ;
//     RooRealVar* rrv_znn(0x0) ;

//     sprintf( parName, "mu_ttwj_sig" ) ;
//     rrv_ttwj = ws->var( parName ) ;
//     if ( rrv_ttwj == 0x0 ) {
//        printf("\n\n\n *** can't find %s in workspace.  Checking other var.\n\n\n", parName ) ;
//        sprintf( parName, "mu_ttwj_sb" ) ;
//        rrv_ttwj = ws->var( parName ) ;
//        if ( rrv_ttwj == 0x0 ) {
//           printf("\n\n\n *** can't find %s either.  I quit.\n\n\n", parName ) ;
//           return ;
//        }
//     }

//     sprintf( parName, "mu_qcd_sig" ) ;
//     rrv_qcd = ws->var( parName ) ;
//     if ( rrv_qcd == 0x0 ) {
//        printf("\n\n\n *** can't find %s in workspace.  Checking other var.\n\n\n", parName ) ;
//        sprintf( parName, "mu_qcd_sb" ) ;
//        rrv_qcd = ws->var( parName ) ;
//        if ( rrv_qcd == 0x0 ) {
//           printf("\n\n\n *** can't find %s either.  I quit.\n\n\n", parName ) ;
//           return ;
//        }
//     }

//     sprintf( parName, "mu_znn_sig" ) ;
//     rrv_znn = ws->var( parName ) ;
//     if ( rrv_znn == 0x0 ) {
//        printf("\n\n\n *** can't find %s in workspace.  I quit.\n\n\n", parName ) ;
//        return ;
//     }




//     printf("\n\n\n  ===== Begin Toy study ====================\n\n") ;

//     RooMCStudy mcs( *likelihood,
//                     *observables,
//                     Constrain( *nuisanceParameters ),
//                     Silence(),
//                     FitOptions(PrintLevel(-1),PrintEvalErrors(-1))
//                     ) ;



//    mcs.generate( nToys, 1, true, "" ) ;

//    int    gen_nsig ;
//    double fit_tru_ttwj ;
//    double fit_tru_qcd ;
//    double fit_tru_znn ;
//    double fit_val_ttwj ;
//    double fit_val_qcd ;
//    double fit_val_znn ;
//    double fit_err_ttwj ;
//    double fit_err_qcd ;
//    double fit_err_znn ;
//    int    fit_cov_qual ;

//    fit_tru_ttwj = rrv_ttwj->getVal() ;
//    fit_tru_qcd  = rrv_qcd->getVal() ;
//    fit_tru_znn  = rrv_znn->getVal() ;



//    TTree* tt(0x0) ;

//    tt = new TTree("toytt", "Toy TTree" ) ;

//    tt->Branch( "gen_nsig",  &gen_nsig,  "gen_nsig/I"  ) ;
//    tt->Branch( "fit_val_ttwj",  &fit_val_ttwj,  "fit_val_ttwj/D"  ) ;
//    tt->Branch( "fit_val_qcd",  &fit_val_qcd,  "fit_val_qcd/D"  ) ;
//    tt->Branch( "fit_val_znn",  &fit_val_znn,  "fit_val_znn/D"  ) ;
//    tt->Branch( "fit_err_ttwj",  &fit_err_ttwj,  "fit_err_ttwj/D"  ) ;
//    tt->Branch( "fit_err_qcd",  &fit_err_qcd,  "fit_err_qcd/D"  ) ;
//    tt->Branch( "fit_err_znn",  &fit_err_znn,  "fit_err_znn/D"  ) ;
//    tt->Branch( "fit_tru_ttwj",  &fit_tru_ttwj,  "fit_tru_ttwj/D"  ) ;
//    tt->Branch( "fit_tru_qcd",  &fit_tru_qcd,  "fit_tru_qcd/D"  ) ;
//    tt->Branch( "fit_tru_znn",  &fit_tru_znn,  "fit_tru_znn/D"  ) ;
//    tt->Branch( "fit_cov_qual",  &fit_cov_qual,  "fit_cov_qual/I"  ) ;



//    for ( int ti=0; ti<nToys; ti++ ) {


//       //-- initialize all floating parameters to the values from the pre fit.
//       TIterator* parIter = preFitFloatVals.createIterator() ;
//       while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {
//          RooRealVar* rrv = ws->var( par->GetName() ) ;
//          rrv->setVal( par->getVal() ) ;
//          //printf(" %20s : %8.2f\n", rrv->GetName(), rrv->getVal() ) ;
//       }
//       //printf("\n\n") ;

//       RooAbsData* toyrds = (RooAbsData*) mcs.genData(ti) ;

//       const RooArgSet* dsras = toyrds->get() ;
//       TIterator* obsIter = dsras->createIterator() ;
//       while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {
//          if ( strcmp( obs->GetName(), "Nsig" ) == 0 ) { 
//             gen_nsig = obs->getVal() ;
//             cout << "   found Nsig : " << gen_nsig << endl << flush ;
//             break ;
//          }
//       }

//       RooFitResult* rfr = likelihood->fitTo( *toyrds, Save(true), PrintLevel(-1), PrintEvalErrors(-1) ) ;


//       fit_val_ttwj = rrv_ttwj->getVal() ;
//       fit_val_qcd = rrv_qcd->getVal() ;
//       fit_val_znn = rrv_znn->getVal() ;
//       fit_err_ttwj = rrv_ttwj->getError() ;
//       fit_err_qcd = rrv_qcd->getError() ;
//       fit_err_znn = rrv_znn->getError() ;
//       fit_cov_qual = rfr->covQual() ;


//       tt->Fill() ;

//       printf("\n\n\n  ===== RooDataSet for toy %d ====================\n\n", ti) ;

//       toyrds->Print() ;
//       toyrds->printMultiline(cout, 1, kTRUE, "") ;
//       printf("\n\n\n Toy %4d Fit covariance quality : %d\n", ti, fit_cov_qual ) ;

//       printf("  generated Nsig : %3d\n", gen_nsig ) ;
//       printf("    %20s : %8.2f +/- %8.2f\n", rrv_ttwj->GetName(), rrv_ttwj->getVal(), rrv_ttwj->getError() ) ;
//       printf("    %20s : %8.2f +/- %8.2f\n", rrv_qcd->GetName(), rrv_qcd->getVal(), rrv_qcd->getError() ) ;
//       printf("    %20s : %8.2f +/- %8.2f\n", rrv_znn->GetName(), rrv_znn->getVal(), rrv_znn->getError() ) ;

//       printf("\n\n") ;

//       delete rfr ;


//   } // ti.

//   char toyRootFile[10000] ;
//   gSystem->Exec("mkdir -p output-files") ;
//   sprintf( toyRootFile, "output-files/toy-data-%s.root", sel ) ;
//   TFile f( toyRootFile, "recreate") ;
//   tt->Write() ;



   }






