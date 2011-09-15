
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

   void ws_toystudy1( const char* wsfile = "ws-lm9-ge1btight.root", int nToys=100 ) {


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
          printf(" fixing to zero.\n") ;
          rrv_mu_susy_sig->setVal(0.) ;
          rrv_mu_susy_sig->setConstant(kTRUE) ;
       }

       printf("\n\n\n  ===== Doing a fit ====================\n\n") ;

       RooFitResult* preFitResult = likelihood->fitTo( *rds, Save(true) ) ;
       const RooArgList preFitFloatVals = preFitResult->floatParsFinal() ;
       {
         TIterator* parIter = preFitFloatVals.createIterator() ;
         while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {
            printf(" %20s : %8.2f\n", par->GetName(), par->getVal() ) ;
         }
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


       char parName[1000] ;

       RooRealVar* rrv_ttwj(0x0) ;
       RooRealVar* rrv_qcd(0x0) ;
       RooRealVar* rrv_znn(0x0) ;

       sprintf( parName, "mu_ttwj_sig" ) ;
       rrv_ttwj = ws->var( parName ) ;
       if ( rrv_ttwj == 0x0 ) {
          printf("\n\n\n *** can't find %s in workspace.  Checking other var.\n\n\n", parName ) ;
          sprintf( parName, "mu_ttwj_sb" ) ;
          rrv_ttwj = ws->var( parName ) ;
          if ( rrv_ttwj == 0x0 ) {
             printf("\n\n\n *** can't find %s either.  I quit.\n\n\n", parName ) ;
             return ;
          }
       }

       sprintf( parName, "mu_qcd_sig" ) ;
       rrv_qcd = ws->var( parName ) ;
       if ( rrv_qcd == 0x0 ) {
          printf("\n\n\n *** can't find %s in workspace.  Checking other var.\n\n\n", parName ) ;
          sprintf( parName, "mu_qcd_sb" ) ;
          rrv_qcd = ws->var( parName ) ;
          if ( rrv_qcd == 0x0 ) {
             printf("\n\n\n *** can't find %s either.  I quit.\n\n\n", parName ) ;
             return ;
          }
       }

       sprintf( parName, "mu_znn_sig" ) ;
       rrv_znn = ws->var( parName ) ;
       if ( rrv_znn == 0x0 ) {
          printf("\n\n\n *** can't find %s in workspace.  I quit.\n\n\n", parName ) ;
          return ;
       }




       printf("\n\n\n  ===== Begin Toy study ====================\n\n") ;

       RooMCStudy mcs( *likelihood,
                       *observables,
                       Constrain( *nuisanceParameters ),
                       Silence(),
                       FitOptions(PrintLevel(-1),PrintEvalErrors(-1))
                       ) ;



      mcs.generate( nToys, 1, true, "" ) ;

      int    gen_nsig ;
      double fit_tru_ttwj ;
      double fit_tru_qcd ;
      double fit_tru_znn ;
      double fit_val_ttwj ;
      double fit_val_qcd ;
      double fit_val_znn ;
      double fit_err_ttwj ;
      double fit_err_qcd ;
      double fit_err_znn ;
      int    fit_cov_qual ;

      fit_tru_ttwj = rrv_ttwj->getVal() ;
      fit_tru_qcd  = rrv_qcd->getVal() ;
      fit_tru_znn  = rrv_znn->getVal() ;



      TTree* tt(0x0) ;

      tt = new TTree("toytt", "Toy TTree" ) ;

      tt->Branch( "gen_nsig",  &gen_nsig,  "gen_nsig/I"  ) ;
      tt->Branch( "fit_val_ttwj",  &fit_val_ttwj,  "fit_val_ttwj/D"  ) ;
      tt->Branch( "fit_val_qcd",  &fit_val_qcd,  "fit_val_qcd/D"  ) ;
      tt->Branch( "fit_val_znn",  &fit_val_znn,  "fit_val_znn/D"  ) ;
      tt->Branch( "fit_err_ttwj",  &fit_err_ttwj,  "fit_err_ttwj/D"  ) ;
      tt->Branch( "fit_err_qcd",  &fit_err_qcd,  "fit_err_qcd/D"  ) ;
      tt->Branch( "fit_err_znn",  &fit_err_znn,  "fit_err_znn/D"  ) ;
      tt->Branch( "fit_tru_ttwj",  &fit_tru_ttwj,  "fit_tru_ttwj/D"  ) ;
      tt->Branch( "fit_tru_qcd",  &fit_tru_qcd,  "fit_tru_qcd/D"  ) ;
      tt->Branch( "fit_tru_znn",  &fit_tru_znn,  "fit_tru_znn/D"  ) ;
      tt->Branch( "fit_cov_qual",  &fit_cov_qual,  "fit_cov_qual/I"  ) ;



      for ( int ti=0; ti<nToys; ti++ ) {


         //-- initialize all floating parameters to the values from the pre fit.
         TIterator* parIter = preFitFloatVals.createIterator() ;
         while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {
            RooRealVar* rrv = ws->var( par->GetName() ) ;
            rrv->setVal( par->getVal() ) ;
            //printf(" %20s : %8.2f\n", rrv->GetName(), rrv->getVal() ) ;
         }
         //printf("\n\n") ;

         RooAbsData* toyrds = (RooAbsData*) mcs.genData(ti) ;

         const RooArgSet* dsras = toyrds->get() ;
         TIterator* obsIter = dsras->createIterator() ;
         while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {
            if ( strcmp( obs->GetName(), "Nsig" ) == 0 ) { 
               gen_nsig = obs->getVal() ;
               cout << "   found Nsig : " << gen_nsig << endl << flush ;
               break ;
            }
         }

         RooFitResult* rfr = likelihood->fitTo( *toyrds, Save(true), PrintLevel(-1), PrintEvalErrors(-1) ) ;


         fit_val_ttwj = rrv_ttwj->getVal() ;
         fit_val_qcd = rrv_qcd->getVal() ;
         fit_val_znn = rrv_znn->getVal() ;
         fit_err_ttwj = rrv_ttwj->getError() ;
         fit_err_qcd = rrv_qcd->getError() ;
         fit_err_znn = rrv_znn->getError() ;
         fit_cov_qual = rfr->covQual() ;


         tt->Fill() ;

         printf("\n\n\n  ===== RooDataSet for toy %d ====================\n\n", ti) ;

         toyrds->Print() ;
         toyrds->printMultiline(cout, 1, kTRUE, "") ;
         printf("\n\n\n Toy %4d Fit covariance quality : %d\n", ti, fit_cov_qual ) ;

         printf("  generated Nsig : %3d\n", gen_nsig ) ;
         printf("    %20s : %8.2f +/- %8.2f\n", rrv_ttwj->GetName(), rrv_ttwj->getVal(), rrv_ttwj->getError() ) ;
         printf("    %20s : %8.2f +/- %8.2f\n", rrv_qcd->GetName(), rrv_qcd->getVal(), rrv_qcd->getError() ) ;
         printf("    %20s : %8.2f +/- %8.2f\n", rrv_znn->GetName(), rrv_znn->getVal(), rrv_znn->getError() ) ;

         printf("\n\n") ;

         delete rfr ;


     } // ti.

     char toyRootFile[10000] ;
     gSystem->Exec("mkdir -p output-files") ;
     sprintf( toyRootFile, "output-files/toy-data-%s.root", sel ) ;
     TFile f( toyRootFile, "recreate") ;
     tt->Write() ;



   }






