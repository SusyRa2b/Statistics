
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
#include "RooMinuit.h"
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

   void ws_fittable_halfblind( const char* wsfile = "rootfiles/ws-data-unblind.root",
                            double mu_susy_sig_val = 0.
                             ) {





     TFile* wstf = new TFile( wsfile ) ;

     RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
     ws->Print() ;

     int nBinsBtag(3) ;
     int nBinsMET(4) ;
     int nBinsHT(4) ;

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





    //----------------------------------------------------------------------------------------------

       //RooMinuit* rminuit( 0x0 ) ;

       //RooFormulaVar* plot_var( 0x0 ) ;

       char ignore_observable_name[50][100] ;

       sprintf( ignore_observable_name[0], "N_0lep_M4_H2_2b" ) ;
       sprintf( ignore_observable_name[1], "N_0lep_M4_H3_2b" ) ;
       sprintf( ignore_observable_name[2], "N_0lep_M4_H4_2b" ) ;

       sprintf( ignore_observable_name[3], "N_0lep_M2_H1_3b" ) ;
       sprintf( ignore_observable_name[4], "N_0lep_M2_H2_3b" ) ;
       sprintf( ignore_observable_name[5], "N_0lep_M2_H3_3b" ) ;
       sprintf( ignore_observable_name[6], "N_0lep_M2_H4_3b" ) ;

       sprintf( ignore_observable_name[7], "N_0lep_M3_H1_3b" ) ;
       sprintf( ignore_observable_name[8], "N_0lep_M3_H2_3b" ) ;
       sprintf( ignore_observable_name[9], "N_0lep_M3_H3_3b" ) ;
       sprintf( ignore_observable_name[10],"N_0lep_M3_H4_3b" ) ;

       sprintf( ignore_observable_name[11],"N_0lep_M4_H2_3b" ) ;
       sprintf( ignore_observable_name[12],"N_0lep_M4_H3_3b" ) ;
       sprintf( ignore_observable_name[13],"N_0lep_M4_H4_3b" ) ;

       int n_ignore_observable(14) ;

       RooAbsReal* rar_ignore_pdf[100] ;
       RooAbsReal* rar_ignore_obs[100] ;

       char ignoreTermFormula[10000] ;
       RooArgSet  ignorePdfList ;

       for ( int ii=0; ii < n_ignore_observable; ii++ ) {

          char name[100] ;

          sprintf( name, "pdf_%s", ignore_observable_name[ii] ) ;
          rar_ignore_pdf[ii] = ws -> pdf( name ) ;
          if ( rar_ignore_pdf[ii] == 0x0 ) {
             printf("\n\n\n *** Told to ignore %s but can't find %s\n\n", ignore_observable_name[ii], name ) ;
             return ;
          }

          rar_ignore_obs[ii] = ws -> var( ignore_observable_name[ii] ) ;
          ( (RooRealVar*) rar_ignore_obs[ii] ) -> setConstant(kTRUE) ; // probably not necessary, but can't hurt.

          if ( ii==0 ) {
             sprintf( ignoreTermFormula, "log(@%d)", ii ) ;
          } else {
             char buffer[10000] ;
             sprintf( buffer, "%s+log(@%d)", ignoreTermFormula, ii ) ;
             sprintf( ignoreTermFormula, "%s", buffer ) ;
          }

          ignorePdfList.add( *rar_ignore_pdf[ii] ) ;

       } // ii

       printf("\n\n Creating ignore formula var with : %s\n", ignoreTermFormula ) ;

       RooFormulaVar* rfv_ignorePdfTerm = new RooFormulaVar("ignorePdfTerm", ignoreTermFormula, ignorePdfList ) ;
       rfv_ignorePdfTerm -> Print() ;


       RooAbsReal* nll = likelihood -> createNLL( *rds, Verbose(true) ) ;

       char minuit_formula_unbiased_unconstrained[100000] ;
       sprintf( minuit_formula_unbiased_unconstrained, "%s+%s", nll->GetName(), rfv_ignorePdfTerm->GetName() ) ;
       RooFormulaVar* rfv_minuitvar_unbiased_unconstrained = new RooFormulaVar( "minuitvar_unbiased_unconstrained",
         minuit_formula_unbiased_unconstrained,
         RooArgList( *nll, *rfv_ignorePdfTerm ) ) ;

       RooMinuit* rminuit_ub_uc = new RooMinuit( *rfv_minuitvar_unbiased_unconstrained  ) ;

       int verbLevel = 0 ;

       rminuit_ub_uc->setPrintLevel(verbLevel-1) ;
       rminuit_ub_uc->setNoWarn() ;

       // rminuit_ub_uc->migrad() ;
       // rminuit_ub_uc->hesse() ;

       RooFitResult* rfr_ub_uc = rminuit_ub_uc->fit("mr") ;

       double floatParInitVal[10000] ;
       char   floatParName[10000][100] ;
       int nFloatParInitVal(0) ;
       RooArgList ral_floats = rfr_ub_uc->floatParsFinal() ;
       TIterator* floatParIter = ral_floats.createIterator() ;
       {
          RooRealVar* par ;
          while ( (par = (RooRealVar*) floatParIter->Next()) ) {
             sprintf( floatParName[nFloatParInitVal], "%s", par->GetName() ) ;
             floatParInitVal[nFloatParInitVal] = par->getVal() ;
             nFloatParInitVal++ ;
          }
       }


       //printf("\n\n Unbiased best value for new POI %s is : %7.1f\n\n", new_poi_rar->GetName(), new_poi_rar->getVal() ) ;

       //double best_unbiased_poi_val = new_poi_rar->getVal() ;

    //----------------------------------------------------------------------------------------------







 //  printf("\n ==== Final floating parameter values\n\n") ;
 //  const RooArgList fitFloatVals = fitResult->floatParsFinal() ;
 //  {
 //    TIterator* parIter = fitFloatVals.createIterator() ;
 //    while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {
 //        printf(" %20s : %8.2f\n", par->GetName(), par->getVal() ) ;
 //    }
 //  }


 //  printf("\n ==== Constant parameter values\n\n") ;
 //  const RooArgList fitConstVals = fitResult->constPars() ;
 //  {
 //    TIterator* parIter = fitConstVals.createIterator() ;
 //    while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {
 //        printf(" %20s : %8.2f\n", par->GetName(), par->getVal() ) ;
 //    }
 //  }

 //  printf("\n ==== Function values\n\n") ;
 //  RooArgSet funcs = ws->allFunctions() ;
 //  TIterator* funcIter = funcs.createIterator() ;
 //  while ( RooFormulaVar* func = (RooFormulaVar*) funcIter->Next() ) {
 //    printf(" %20s : %8.2f\n", func->GetName(), func->getVal() ) ;
 //  }

 //  printf("\n\n") ; cout << flush ;










   //--- collect everything from workspace.

     printf("\n\n\n\n") ; cout << flush ;


     int ncomp(4) ;
     char comp    [4][100] = {       "ttwj",     "qcd",        "znn",          "vv" } ;
     char trigtype[4][100] = { "trigeff_sl", "trigeff", "trigeff_sl",  "trigeff_sl" } ;

     int nsel(3) ;
     char ws_selname    [3][100] = {     "",  "_sl", "_ldp" } ;
     char output_selname[3][100] = {   "zl",   "sl",  "ldp" } ;
     char obs_selname   [3][100] = { "0lep", "1lep",  "ldp" } ;

     double mu_raw  [4][4][3][3][4] ; // [mbi][hbi][bbi][si][ci] : met, ht, nb, sel, comp
     double mu_wtrig[4][4][3][3][4] ; // [mbi][hbi][bbi][si][ci] : met, ht, nb, sel, comp
     int    N_obs   [4][4][3][3]    ; // [mbi][hbi][bbi][si]     : met, ht, nb, sel

     bool comp_included[3][4] ; // [si][ci] : sel, comp
     comp_included[0][0] = true ; // zl, ttwj
     comp_included[0][1] = true ; // zl, qcd
     comp_included[0][2] = true ; // zl, znn
     comp_included[0][3] = true ; // zl, vv
     comp_included[1][0] = true  ; // sl, ttwj
     comp_included[1][1] = false ; // sl, qcd
     comp_included[1][2] = false ; // sl, znn
     comp_included[1][3] = false ; // sl, vv
     comp_included[2][0] = true ; // ldp, ttwj
     comp_included[2][1] = true ; // ldp, qcd
     comp_included[2][2] = true ; // ldp, znn
     comp_included[2][3] = false ; // ldp, vv

     for ( int si=0; si<nsel; si++ ) {
        for ( int ci=0; ci<ncomp; ci++ ) {

           if ( !comp_included[si][ci] ) continue ;

           for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {
              for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
                 for ( int hbi=0; hbi<nBinsHT; hbi++ ) {

                    if ( ignoreBin[mbi][hbi] ) continue ;

                    char pname[1000] ;
                    RooAbsReal* rar(0x0) ;

                    mu_raw  [mbi][hbi][bbi][si][ci] = 0. ;
                    mu_wtrig[mbi][hbi][bbi][si][ci] = 0. ;

                    sprintf( pname, "N_%s_M%d_H%d_%db", obs_selname[si], mbi+1, hbi+1, bbi+1 ) ;
                    rar = (RooAbsReal*) ws->obj(pname) ;
                    //if ( rar == 0x0 ) { printf("\n\n *** %s missing from ws.\n\n", pname ) ; return ; }
                    if ( rar == 0x0 ) { printf("\n\n *** %s missing from ws.\n\n", pname ) ; continue ; }
                    int observed_events = TMath::Nint( rar->getVal() ) ;

                    sprintf( pname, "%s_M%d_H%d", trigtype[ci], mbi+1, hbi+1 ) ;
                    rar = (RooAbsReal*) ws->obj(pname) ;
                    //if ( rar == 0x0 ) { printf("\n\n *** %s missing from ws.\n\n", pname ) ; return ; }
                    if ( rar == 0x0 ) { printf("\n\n *** %s missing from ws.\n\n", pname ) ; continue ; }
                    double trig_eff = rar->getVal() ;

                    sprintf( pname, "mu_%s%s_M%d_H%d_%db", comp[ci], ws_selname[si], mbi+1, hbi+1, bbi+1 ) ;
                    rar = (RooAbsReal*) ws->obj(pname) ;
                    //if ( rar == 0x0 ) { printf("\n\n *** %s missing from ws.\n\n", pname ) ; return ; }
                    if ( rar == 0x0 ) { printf("\n\n *** %s missing from ws.\n\n", pname ) ; continue ; }
                    double mu_val = rar->getVal() ;

                    N_obs   [mbi][hbi][bbi][si]     = observed_events ;
                    mu_raw  [mbi][hbi][bbi][si][ci] = mu_val ;
                    mu_wtrig[mbi][hbi][bbi][si][ci] = trig_eff * mu_val ;

                    printf(" %s :  mu = %7.1f,  trig_eff = %6.3f,  N_obs = %6d\n", pname, mu_val, trig_eff, observed_events ) ;


                 } // hbi.
                 printf("\n") ;
              } // mbi.
              printf("\n\n") ;
           } // bbi.
           printf("\n ---------------- \n") ;
        } // ci.
        printf("\n ============================== \n") ;
     } // si.




   //--- print out all observables and fit totals.

     printf("\n\n\n\n") ; cout << flush ;

     for ( int si=0; si<nsel; si++ ) {

        for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {
           for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
              for ( int hbi=0; hbi<nBinsHT; hbi++ ) {

                 if ( ignoreBin[mbi][hbi] ) continue ;

                 double fitsum(0.) ;

                 for ( int ci=0; ci<ncomp; ci++ ) {

                    if ( !comp_included[si][ci] ) continue ;

                    fitsum += mu_wtrig[mbi][hbi][bbi][si][ci] ;

                 } // ci.

                 printf( "%3s M%d, H%d, %db :  obs=%6d,  fit=%8.1f  : ",
                         output_selname[si], mbi+1, hbi+1, bbi+1,
                         N_obs[mbi][hbi][bbi][si],
                         fitsum ) ;

                 for ( int ci=0; ci<ncomp; ci++ ) {
                    if ( !comp_included[si][ci] ) continue ;
                    printf( "%4s = %7.1f  ", comp[ci], mu_wtrig[mbi][hbi][bbi][si][ci] ) ;
                 } // ci.
                 printf("\n") ;

              } // hbi.
              printf("\n") ;
           } // mbi.
           printf("\n\n") ;
        } // bbi.
        printf("\n ---------------- \n") ;
     } // si.







   //--- compute all relevant integrals for components and save in a file.

     printf("\n\n\n\n") ; cout << flush ;

  // TString ws_fname( wsfile ) ;
  // ws_fname.ReplaceAll("rootfiles/","") ;
  // ws_fname.ReplaceAll(".root","") ;

     char command[1000] ;
     sprintf( command, "basename %s", wsfile ) ;
     TString ws_fname = gSystem->GetFromPipe( command ) ;
     ws_fname.ReplaceAll(".root","") ;


     char outfilename[10000] ;
     if ( mu_susy_sig_val < 0. ) {
        sprintf( outfilename, "outputfiles/fitresults-halfblind-%s-susyFloat.txt", ws_fname.Data() ) ;
     } else {
        sprintf( outfilename, "outputfiles/fitresults-halfblind-%s-susyFixed%.1f.txt", ws_fname.Data(), mu_susy_sig_val ) ;
     }
     printf(" ws_fname : %s\n", ws_fname.Data() ) ; cout << flush ;
     printf("\n\n Opening output file : %s\n\n", outfilename ) ; cout << flush ;
     FILE* outfile = fopen( outfilename, "w" ) ;

     for ( int si=0; si<nsel; si++ ) {
        for ( int ci=0; ci<ncomp; ci++ ) {

           if ( !comp_included[si][ci] ) continue ;

           double all(0.) ;

           for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {

              double this_nb(0.) ;

              for ( int mbi=0; mbi<nBinsMET; mbi++ ) {

                 double this_nb_met(0.) ;

                 for ( int hbi=0; hbi<nBinsHT; hbi++ ) {

                    if ( ignoreBin[mbi][hbi] ) continue ;

                    this_nb_met += mu_wtrig[mbi][hbi][bbi][si][ci] ;
                    this_nb     += mu_wtrig[mbi][hbi][bbi][si][ci] ;
                    all         += mu_wtrig[mbi][hbi][bbi][si][ci] ;

                 } // hbi.

                 fprintf( outfile, "%s_%s_wt_M%d_%db %.2f\n", comp[ci], output_selname[si], mbi+1, bbi+1, this_nb_met ) ;

              } // mbi.

              fprintf( outfile, "%s_%s_wt_%db %.2f\n", comp[ci], output_selname[si], bbi+1, this_nb ) ;

              for ( int hbi=0; hbi<nBinsHT; hbi++ ) {

                 double this_nb_ht(0.) ;

                 for ( int mbi=0; mbi<nBinsMET; mbi++ ) {

                    if ( ignoreBin[mbi][hbi] ) continue ;

                    this_nb_ht  += mu_wtrig[mbi][hbi][bbi][si][ci] ;

                 } // mbi.

                 fprintf( outfile, "%s_%s_wt_H%d_%db %.2f\n", comp[ci], output_selname[si], hbi+1, bbi+1, this_nb_ht ) ;

              } // hbi.

           } // bbi.

           fprintf( outfile, "%s_%s_wt %.2f\n", comp[ci], output_selname[si], all ) ;

        } // ci.
     } // si.

     fclose( outfile ) ;


















 //  //--- Extract model and observed values for backgrounds in key bins


 //  double n_0lep_2b_model[4] ;
 //  double n_0lep_2b_obs[4] ;

 //  for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
 //     char pname[1000] ;
 //     sprintf( pname, "n_M4_H%d_2b", hbi+1 ) ;
 //     RooAbsReal* par = (RooAbsReal*) ws->obj( pname ) ;
 //     if ( par == 0x0 ) {
 //        printf(" %s missing.\n", pname) ;
 //        n_0lep_2b_model[hbi] = 0. ;
 //     } else {
 //        printf(" %s found\n", pname ) ;
 //        printf(" %s = %g\n", pname, par->getVal() ) ;
 //        n_0lep_2b_model[hbi] = par->getVal() ;
 //     }
 //     sprintf( pname, "N_0lep_M4_H%d_2b", hbi+1 ) ;
 //     par = (RooAbsReal*) ws->obj( pname ) ;
 //     if ( par == 0x0 ) {
 //        printf(" %s missing.\n", pname) ;
 //        n_0lep_2b_obs[hbi] = 0. ;
 //     } else {
 //        printf(" %s found\n", pname ) ;
 //        printf(" %s = %g\n", pname, par->getVal() ) ;
 //        n_0lep_2b_obs[hbi] = par->getVal() ;
 //     }
 //  } // hbi.


 //  double n_0lep_3b_model[4][4] ;
 //  double n_0lep_3b_obs[4][4] ;

 //  for ( int mbi=1; mbi<nBinsMET; mbi++ ) {
 //     for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
 //        char pname[1000] ;
 //        sprintf( pname, "n_M%d_H%d_3b", mbi+1, hbi+1 ) ;
 //        RooAbsReal* par = (RooAbsReal*) ws->obj( pname ) ;
 //        if ( par == 0x0 ) {
 //           printf(" %s missing.\n", pname) ;
 //           n_0lep_3b_model[mbi][hbi] = 0. ;
 //        } else {
 //           printf(" %s found\n", pname ) ;
 //           printf(" %s = %g\n", pname, par->getVal() ) ;
 //           n_0lep_3b_model[mbi][hbi] = par->getVal() ;
 //        }
 //        sprintf( pname, "N_0lep_M%d_H%d_3b", mbi+1, hbi+1 ) ;
 //        par = (RooAbsReal*) ws->obj( pname ) ;
 //        if ( par == 0x0 ) {
 //           printf(" %s missing.\n", pname) ;
 //           n_0lep_2b_obs[hbi] = 0. ;
 //        } else {
 //           printf(" %s found\n", pname ) ;
 //           printf(" %s = %g\n", pname, par->getVal() ) ;
 //           n_0lep_3b_obs[mbi][hbi] = par->getVal() ;
 //        }
 //     } // hbi.
 //  } // mbi





 ////--- print pretty tables

 // {

 //  printf("\n\n ===== nB == 2 =====================\n\n") ;

 //  printf("\n\n\n") ;
 //  printf("   Model values\n") ;
 //  printf("             H1      H2      H3      H4       H1-4\n") ;
 //  printf("---------------------------------------------------------\n") ;
 //  printf("  M4   | " ) ;
 //  double model_hsum(0.) ;
 //  for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
 //     printf(" %6.1f ", n_0lep_2b_model[hbi] ) ;
 //     model_hsum += n_0lep_2b_model[hbi] ;
 //  } // hbi.
 //  printf(" | %6.1f\n", model_hsum ) ;


 //  printf("\n\n\n") ;
 //  printf("   Observed values\n") ;
 //  printf("             H1      H2      H3      H4       H1-4\n") ;
 //  printf("---------------------------------------------------------\n") ;
 //  printf("  M4   | " ) ;
 //  double obs_hsum(0.) ;
 //  for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
 //     printf(" %6.1f ", n_0lep_2b_obs[hbi] ) ;
 //     obs_hsum += n_0lep_2b_obs[hbi] ;
 //  } // hbi.
 //  printf(" | %6.1f\n", obs_hsum ) ;



 // }



 // {

 //  printf("\n\n ===== nB >= 3 =====================\n\n") ;

 //  double model_metsum[4] = {0., 0., 0., 0.} ;
 //  printf("\n\n\n") ;
 //  printf("   Model values\n") ;
 //  printf("             H1      H2      H3      H4       H1-4\n") ;
 //  printf("---------------------------------------------------------\n") ;
 //  for ( int mbi=1; mbi<nBinsMET; mbi++ ) {
 //     printf("  M%d   | ", mbi+1 ) ;
 //     double hsum(0.) ;
 //     for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
 //        printf(" %6.1f ", n_0lep_3b_model[mbi][hbi] ) ;
 //        hsum += n_0lep_3b_model[mbi][hbi] ;
 //        model_metsum[hbi] += n_0lep_3b_model[mbi][hbi] ;
 //     } // hbi.
 //     printf(" | %6.1f\n", hsum ) ;
 //  } // mbi
 //  printf("---------------------------------------------------------\n") ;
 //  printf("  M2-4 | ") ;
 //  double model_allsum(0.) ;
 //  for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
 //     printf(" %6.1f ", model_metsum[hbi] ) ;
 //     model_allsum += model_metsum[hbi] ;
 //  } // hbi.
 //  printf(" | %6.1f\n", model_allsum ) ;


 //  double obs_metsum[4] = {0., 0., 0., 0.} ;
 //  printf("\n\n\n") ;
 //  printf("   Observed values\n") ;
 //  printf("             H1      H2      H3      H4       H1-4\n") ;
 //  printf("---------------------------------------------------------\n") ;
 //  for ( int mbi=1; mbi<nBinsMET; mbi++ ) {
 //     printf("  M%d   | ", mbi+1 ) ;
 //     double hsum(0.) ;
 //     for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
 //        printf(" %4.0f   ", n_0lep_3b_obs[mbi][hbi] ) ;
 //        hsum += n_0lep_3b_obs[mbi][hbi] ;
 //        obs_metsum[hbi] += n_0lep_3b_obs[mbi][hbi] ;
 //     } // hbi.
 //     printf(" | %4.0f\n", hsum ) ;
 //  } // mbi
 //  printf("---------------------------------------------------------\n") ;
 //  printf("  M2-4 | ") ;
 //  double obs_allsum(0.) ;
 //  for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
 //     printf(" %4.0f   ", obs_metsum[hbi] ) ;
 //     obs_allsum += obs_metsum[hbi] ;
 //  } // hbi.
 //  printf(" | %4.0f\n", obs_allsum ) ;

 // }








     printf("\n\n\n\n") ; cout << flush ;



   } // ws_fittable









