
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"

  using namespace RooFit;
  using namespace RooStats;

   void ws_profile_interval1( const char* wsfile = "ws-test1.root", const char* parName = "mu_susy_sig", double alpha = 0.10, double mu_susy_sig_val = 0., double xmax = -1. ) {




       TFile* wstf = new TFile( wsfile ) ;

       RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
       ws->Print() ;






 ////  ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;

 ////  printf("\n\n\n  ===== SbModel ====================\n\n") ;
 ////  modelConfig->Print() ;







       RooDataSet* rds = (RooDataSet*) ws->obj( "ra2b_observed_rds" ) ;
       printf("\n\n\n  ===== RooDataSet ====================\n\n") ;

       rds->Print() ;
       rds->printMultiline(cout, 1, kTRUE, "") ;









       printf("\n\n\n  ===== Grabbing %s rrv ====================\n\n", parName ) ;
       RooRealVar* rrv_par = ws->var( parName ) ;
       if ( rrv_par == 0x0 ) {
          printf("\n\n\n *** can't find %s in workspace.  Quitting.\n\n\n", parName ) ; return ;
       } else {
          printf(" current value is : %8.3f\n", rrv_par->getVal() ) ; cout << flush ;
       }
       if ( xmax > 0 ) { rrv_par->setMax( xmax ) ; }

       printf("\n\n\n  ===== Grabbing mu_susy_sig rrv ====================\n\n") ;
       RooRealVar* rrv_mu_susy_sig = ws->var("mu_susy_sig") ;
       if ( rrv_mu_susy_sig == 0x0 ) {
          printf("\n\n\n *** can't find mu_susy_sig in workspace.  Quitting.\n\n\n") ; return ;
       }
       if ( strcmp( parName, "mu_susy_sig" ) != 0 ) {
          if ( mu_susy_sig_val >= 0. ) {
             printf(" current value is : %8.3f\n", rrv_mu_susy_sig->getVal() ) ; cout << flush ;
             printf(" fixing to %8.2f.\n", mu_susy_sig_val ) ;
             rrv_mu_susy_sig->setVal( mu_susy_sig_val ) ;
             rrv_mu_susy_sig->setConstant(kTRUE) ;
          } else {
             printf(" current value is : %8.3f\n", rrv_mu_susy_sig->getVal() ) ; cout << flush ;
             printf(" allowing mu_susy_sig to float.\n") ;
             rrv_mu_susy_sig->setConstant(kFALSE) ;
          }
       } else {
          printf("\n\n profile plot parameter is mu_susy_sig.\n") ;
          rrv_mu_susy_sig->setConstant(kFALSE) ;
       }

       printf("\n\n\n  ===== Grabbing likelihood pdf ====================\n\n") ;
       RooAbsPdf* likelihood = ws->pdf("likelihood") ;
       if ( likelihood == 0x0 ) {
          printf("\n\n\n *** can't find likelihood pdf in workspace.  Quitting.\n\n\n") ; return ;
       } else {
          printf("\n\n likelihood pdf: \n\n") ;
          likelihood->Print() ;
       }












       printf("\n\n\n  ===== Doing a fit ====================\n\n") ;

       likelihood->fitTo( *rds ) ;

       double mlValue = rrv_par->getVal() ;
       printf("  Maximum likelihood value of %s : %8.3f +/- %8.3f\n",
            parName, rrv_par->getVal(), rrv_par->getError() ) ;






       printf("\n\n ========== Creating ProfileLikelihoodCalculator\n\n" ) ; cout << flush ;

    // ProfileLikelihoodCalculator plc( *rds, *modelConfig ) ;
       ProfileLikelihoodCalculator plc( *rds, *likelihood, RooArgSet( *rrv_par ) ) ;

       plc.SetTestSize( alpha ) ;
       ConfInterval* plinterval = plc.GetInterval() ;
       double low  = ((LikelihoodInterval*) plinterval)->LowerLimit(*rrv_par) ;
       double high = ((LikelihoodInterval*) plinterval)->UpperLimit(*rrv_par) ;

       printf("\n\n  Limits: %8.3f,  %8.3f\n\n", low, high ) ;








       printf("\n\n ========= Making profile likelihood plot\n\n") ; cout << flush ;

       LikelihoodIntervalPlot* profPlot = new LikelihoodIntervalPlot((LikelihoodInterval*)plinterval) ;

       TCanvas* cplplot = new TCanvas("cplplot","cplplot", 500, 400) ;
       profPlot->Draw() ;
       gPad->SetGridy(1) ;
       char plotname[10000] ;
       sprintf( plotname, "plplot-%s.png", parName ) ;
       cplplot->SaveAs( plotname ) ;



       if ( alpha > 0.3 ) {
          printf("\n\n\n 1 standard-deviation errors for %s : %8.2f   + %8.2f  - %8.2f\n\n\n",
              parName, mlValue, high-mlValue, mlValue-low ) ;
       }



   }






