
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooPlot.h"

//----------
#include "RooStats/LikelihoodIntervalPlot.h"
//#include "LikelihoodIntervalPlot.cxx"
//----------

  using namespace RooFit;
  using namespace RooStats;

   void ws_profile_interval2n( const char* wsfile = "ws-test1.root",
                              const char* parName = "mu_susy_sig_1b",
                              double alpha = 0.10,
                              double mu_susy_sig_1b_val = 0.,
                              double xmax = -1.,
                              double efficiency = -1.,
                              double integratedLumi = 4684.
                              ) {

       //double efficiency = (3039. / 413270.) * 0.93 = 0.006839 ; // 1BL



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

       printf("\n\n\n  ===== Grabbing mu_susy_sig_1b rrv ====================\n\n") ;
       RooRealVar* rrv_mu_susy_sig_1b = ws->var("mu_susy_sig_1b") ;
       if ( rrv_mu_susy_sig_1b == 0x0 ) {
          printf("\n\n\n *** can't find mu_susy_sig_1b in workspace.  Quitting.\n\n\n") ; return ;
       }
       if ( strcmp( parName, "mu_susy_sig_1b" ) != 0 ) {
          if ( mu_susy_sig_1b_val >= 0. ) {
             printf(" current value is : %8.3f\n", rrv_mu_susy_sig_1b->getVal() ) ; cout << flush ;
             printf(" fixing to %8.2f.\n", mu_susy_sig_val_1b ) ;
             rrv_mu_susy_sig_1b->setVal( mu_susy_sig_1b_val ) ;
             rrv_mu_susy_sig_1b->setConstant(kTRUE) ;
          } else {
             printf(" current value is : %8.3f\n", rrv_mu_susy_sig_1b->getVal() ) ; cout << flush ;
             printf(" allowing mu_susy_sig_1b to float.\n") ;
             rrv_mu_susy_sig_1b->setConstant(kFALSE) ;
          }
       } else {
          printf("\n\n profile plot parameter is mu_susy_sig_1b.\n") ;
          rrv_mu_susy_sig_1b->setConstant(kFALSE) ;
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
  ///  gDirectory->ls() ;

  ///  TFile fout("fout.root","recreate") ;
  ///  profPlot->Write("theProfilePlot") ;
  ///  fout.Close() ;

       gPad->SetGridy(1) ;
       char plotname[10000] ;
       sprintf( plotname, "plplot-%s.png", parName ) ;
       cplplot->SaveAs( plotname ) ;



       if ( alpha > 0.3 ) {
          printf("\n\n\n 1 standard-deviation errors for %s : %8.2f   + %8.2f  - %8.2f\n\n\n",
              parName, mlValue, high-mlValue, mlValue-low ) ;
       }

       RooPlot* rp = (RooPlot*) gDirectory->FindObject("LIP_profile") ;
       if ( rp != 0x0 ) {

          if ( efficiency > 0. ) {

             rp->Dump() ;
             cout << flush ;
             RooCurve* rc = (RooCurve*) rp->findObject("nll_likelihood_ra2b_observed_rds_with_constr_Profile[mu_susy_sig_1b]_Norm[mu_susy_sig_1b]") ;

             if ( rc != 0x0 ) {

                printf("\n\n Found RooCurve.\n\n") ; cout << flush ;
             ///TCanvas* rccanv = new TCanvas("rccanv","rccanv") ;
             ///rc->Draw("AC*") ;
             ///rccanv->SaveAs("rccanv.png") ;

                Double_t xsec[10000] ;
                Double_t* mu_susy_sig = rc->GetX() ;
                Double_t* delta_log_L = rc->GetY() ;
                int nPoints = rc->GetN() ;


                printf("\n\n Number of points in TGraph: %d\n\n", nPoints ) ; cout << flush ;
                double xsec_for_vline(-1.) ;
                double dll_for_vline(-1.) ;
                for ( int pi=0; pi<nPoints; pi++ ) {
                   xsec[pi] = mu_susy_sig[pi] / ( integratedLumi * efficiency ) ;
                   printf(" %3d : dll=%8.3f   mu=%8.1f   xSec=%8.3f\n", pi, delta_log_L[pi],
                           mu_susy_sig[pi], xsec[pi] ) ;
                   if ( delta_log_L[pi] > 1.355 && xsec_for_vline<0 ) {
                      xsec_for_vline = xsec[pi] ;
                      dll_for_vline = delta_log_L[pi] ;
                   }
                } // pi.

                TLine* line = new TLine() ;
                line->SetLineColor(kGreen) ;
                TGraph* newgraph = new TGraph( nPoints-1, xsec, delta_log_L ) ; //-- truncate last point.
                TCanvas* xseccanv = new TCanvas("xseccanv","xseccanv") ;
                newgraph->SetLineColor(4) ;
                newgraph->SetLineWidth(3) ;
                newgraph->Draw("AC") ;
                line->DrawLine( xsec_for_vline, 0., xsec_for_vline, dll_for_vline ) ;
                newgraph->GetXaxis()->SetRangeUser( 0., xsec[nPoints-2] ) ;
                newgraph->GetYaxis()->SetRangeUser( 0., 2.0 ) ;
                line->DrawLine( 0., 1.355, xsec[nPoints-2], 1.355 ) ;
                xseccanv->Update() ;
                gPad->SetGridy(1) ;
                xseccanv->SaveAs("deltall-vs-xsec.png") ;


             } else {
                printf("\n\n *** Can't find RooCurve with name : nll_likelihood_ra2b_observed_rds_with_constr_Profile[mu_susy_sig_1b]_Norm[mu_susy_sig_1b]\n\n") ; cout << flush ;
             }
             printf("\n\n Deleteing RooPlot object.\n\n") ; cout << flush ;

          }

          delete rp ;
          printf("\n\n After deleteing RooPlot object.\n\n") ; cout << flush ;

       } else {
          printf("\n\n\n *** Can't find LIP_profile object in memory.\n\n\n") ; cout << flush ;
       }

       printf("\n\n Exiting now.\n\n") ; cout << flush ;


   }






