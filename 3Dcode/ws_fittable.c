
#include "TMath.h"
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
#include "TGaxis.h"
#include "TText.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooConstVar.h"
#include "RooStats/ModelConfig.h"
#include "TRegexp.h"

#include <string.h>
#include <sstream>

  using namespace RooFit;
  using namespace RooStats;


  //------
  //
  // Note: If the 2nd argument (mu_susy_sig_val) is negative, mu_susy_sig will be floated
  //        in the fit.  If it's zero or positive, mu_susy_sig will be fixed to that value
  //        in the fit.
  //
  //------

   void ws_fittable( const char* wsfile = "rootfiles/ws-data-halfblind.root",
                            double mu_susy_sig_val = 0. ) {





     TFile* wstf = new TFile( wsfile ) ;

     RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
     ws->Print() ;


     // figure out the binning from the variables in the likelihood
     int nBinsMET(1) ;
     for ( int bi=2; bi<20; bi++ ) {
        char vname[1000] ;
        sprintf( vname, "mu_ttwj_sl_M%d_H2_1b", bi ) ;
        RooRealVar* rrv = ws->var( vname ) ;
        if ( rrv == 0x0 ) break ;
        nBinsMET = bi ;
     }

     int nBinsHT(1) ;
     for ( int bi=2; bi<20; bi++ ) {
        char vname[1000] ;
        sprintf( vname, "mu_ttwj_sl_M1_H%d_1b", bi ) ;
        RooRealVar* rrv = ws->var( vname ) ;
        if ( rrv == 0x0 ) break ;
        nBinsHT = bi ;
     }

      //-- Hardwire in to ignore the highest MET bin in the lowest HT bin.
      bool ignoreBin[4][4] ;
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            ignoreBin[mbi][hbi] = false ;
         } // hbi
      } // mbi
      ignoreBin[nBinsMET-1][0] = true ;

      printf("\n\n *** Ignoring these bins in the analysis.\n\n") ;
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            if ( ignoreBin[mbi][hbi] ) {
               printf("  MET %d, HT %d\n", mbi+1, hbi+1 ) ;
            }
         } // hbi
      } // mbi
      printf("\n\n\n") ;

     printf("\n\n Binning: nBinsMET = %d, nBinsHT = %d\n\n", nBinsMET, nBinsHT ) ;




     ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;

     printf("\n\n\n  ===== SbModel ====================\n\n") ;
     modelConfig->Print() ;



     RooDataSet* rds = (RooDataSet*) ws->obj( "ra2b_observed_rds" ) ;
     printf("\n\n\n  ===== RooDataSet ====================\n\n") ;

     rds->Print() ;
     rds->printMultiline(cout, 1, kTRUE, "") ;


     RooAbsPdf* likelihood = modelConfig->GetPdf() ;


     RooRealVar* rrv_mu_susy_sig = ws->var("mu_susy_all0lep") ;
     if ( rrv_mu_susy_sig == 0x0 ) {
       printf("\n\n\n *** can't find mu_susy_all0lep in workspace.  Quitting.\n\n\n") ;
       return ;
     } else {
       printf(" current value is : %8.3f\n", rrv_mu_susy_sig->getVal() ) ; cout << flush ;
       rrv_mu_susy_sig->setConstant(kFALSE) ;
     }


     printf("\n\n\n  ===== Doing a fit with SUSY component floating ====================\n\n") ;

     RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(0) ) ;
     double logLikelihoodSusyFloat = fitResult->minNll() ;
     
     double logLikelihoodSusyFixed(0.) ;
     double testStatVal(-1.) ;
     if ( mu_susy_sig_val >= 0. ) {
       printf("\n\n\n  ===== Doing a fit with SUSY fixed ====================\n\n") ;
       printf(" fixing mu_susy_sig to %8.2f.\n", mu_susy_sig_val ) ;
       rrv_mu_susy_sig->setVal( mu_susy_sig_val ) ;
       rrv_mu_susy_sig->setConstant(kTRUE) ;
       
       fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(0) ) ;
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

     printf("\n\n") ; cout << flush ;







     //--- Extract model and observed values for backgrounds in key bins


     double n_0lep_2b_model[4] ;
     double n_0lep_2b_obs[4] ;

     for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
        char pname[1000] ;
        sprintf( pname, "n_M4_H%d_2b", hbi+1 ) ;
        RooAbsReal* par = (RooAbsReal*) ws->obj( pname ) ;
        if ( par == 0x0 ) {
           printf(" %s missing.\n", pname) ;
           n_0lep_2b_model[hbi] = 0. ;
        } else {
           printf(" %s found\n", pname ) ;
           printf(" %s = %g\n", pname, par->getVal() ) ;
           n_0lep_2b_model[hbi] = par->getVal() ;
        }
        sprintf( pname, "N_0lep_M4_H%d_2b", hbi+1 ) ;
        par = (RooAbsReal*) ws->obj( pname ) ;
        if ( par == 0x0 ) {
           printf(" %s missing.\n", pname) ;
           n_0lep_2b_obs[hbi] = 0. ;
        } else {
           printf(" %s found\n", pname ) ;
           printf(" %s = %g\n", pname, par->getVal() ) ;
           n_0lep_2b_obs[hbi] = par->getVal() ;
        }
     } // hbi.


     double n_0lep_3b_model[4][4] ;
     double n_0lep_3b_obs[4][4] ;

     for ( int mbi=1; mbi<nBinsMET; mbi++ ) {
        for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
           char pname[1000] ;
           sprintf( pname, "n_M%d_H%d_3b", mbi+1, hbi+1 ) ;
           RooAbsReal* par = (RooAbsReal*) ws->obj( pname ) ;
           if ( par == 0x0 ) {
              printf(" %s missing.\n", pname) ;
              n_0lep_3b_model[mbi][hbi] = 0. ;
           } else {
              printf(" %s found\n", pname ) ;
              printf(" %s = %g\n", pname, par->getVal() ) ;
              n_0lep_3b_model[mbi][hbi] = par->getVal() ;
           }
           sprintf( pname, "N_0lep_M%d_H%d_3b", mbi+1, hbi+1 ) ;
           par = (RooAbsReal*) ws->obj( pname ) ;
           if ( par == 0x0 ) {
              printf(" %s missing.\n", pname) ;
              n_0lep_2b_obs[hbi] = 0. ;
           } else {
              printf(" %s found\n", pname ) ;
              printf(" %s = %g\n", pname, par->getVal() ) ;
              n_0lep_3b_obs[mbi][hbi] = par->getVal() ;
           }
        } // hbi.
     } // mbi





   //--- print pretty tables

    {

     printf("\n\n ===== nB == 2 =====================\n\n") ;

     printf("\n\n\n") ;
     printf("   Model values\n") ;
     printf("             H1      H2      H3      H4       H1-4\n") ;
     printf("---------------------------------------------------------\n") ;
     printf("  M4   | " ) ;
     double model_hsum(0.) ;
     for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
        printf(" %6.1f ", n_0lep_2b_model[hbi] ) ;
        model_hsum += n_0lep_2b_model[hbi] ;
     } // hbi.
     printf(" | %6.1f\n", model_hsum ) ;


     printf("\n\n\n") ;
     printf("   Observed values\n") ;
     printf("             H1      H2      H3      H4       H1-4\n") ;
     printf("---------------------------------------------------------\n") ;
     printf("  M4   | " ) ;
     double obs_hsum(0.) ;
     for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
        printf(" %6.1f ", n_0lep_2b_obs[hbi] ) ;
        obs_hsum += n_0lep_2b_obs[hbi] ;
     } // hbi.
     printf(" | %6.1f\n", obs_hsum ) ;



    }



    {

     printf("\n\n ===== nB >= 3 =====================\n\n") ;

     double model_metsum[4] = {0., 0., 0., 0.} ;
     printf("\n\n\n") ;
     printf("   Model values\n") ;
     printf("             H1      H2      H3      H4       H1-4\n") ;
     printf("---------------------------------------------------------\n") ;
     for ( int mbi=1; mbi<nBinsMET; mbi++ ) {
        printf("  M%d   | ", mbi+1 ) ;
        double hsum(0.) ;
        for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
           printf(" %6.1f ", n_0lep_3b_model[mbi][hbi] ) ;
           hsum += n_0lep_3b_model[mbi][hbi] ;
           model_metsum[hbi] += n_0lep_3b_model[mbi][hbi] ;
        } // hbi.
        printf(" | %6.1f\n", hsum ) ;
     } // mbi
     printf("---------------------------------------------------------\n") ;
     printf("  M2-4 | ") ;
     double model_allsum(0.) ;
     for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
        printf(" %6.1f ", model_metsum[hbi] ) ;
        model_allsum += model_metsum[hbi] ;
     } // hbi.
     printf(" | %6.1f\n", model_allsum ) ;


     double obs_metsum[4] = {0., 0., 0., 0.} ;
     printf("\n\n\n") ;
     printf("   Observed values\n") ;
     printf("             H1      H2      H3      H4       H1-4\n") ;
     printf("---------------------------------------------------------\n") ;
     for ( int mbi=1; mbi<nBinsMET; mbi++ ) {
        printf("  M%d   | ", mbi+1 ) ;
        double hsum(0.) ;
        for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
           printf(" %4.0f   ", n_0lep_3b_obs[mbi][hbi] ) ;
           hsum += n_0lep_3b_obs[mbi][hbi] ;
           obs_metsum[hbi] += n_0lep_3b_obs[mbi][hbi] ;
        } // hbi.
        printf(" | %4.0f\n", hsum ) ;
     } // mbi
     printf("---------------------------------------------------------\n") ;
     printf("  M2-4 | ") ;
     double obs_allsum(0.) ;
     for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
        printf(" %4.0f   ", obs_metsum[hbi] ) ;
        obs_allsum += obs_metsum[hbi] ;
     } // hbi.
     printf(" | %4.0f\n", obs_allsum ) ;

    }











   }
