
#include "TRandom2.h"
#include "TLegend.h"
#include "TFile.h"
#include "TLine.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TString.h"
#include "TStyle.h"
#include "TText.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooStats/ModelConfig.h"
#include "RooMinuit.h"

#include <string.h>

  using namespace RooFit;
  using namespace RooStats;

  //-- Note: 2nd argument is list of bins to sum, separated by semicolons.
  //         Here's an example: "n_M2_H4_3b:n_M3_H4_3b:n_M4_H4_3b"

  //==============================================================================================

   void ws_halfblind_binsum_profile3D( const char* wsfile = "rootfiles/ws-data-unblind.root",
                                   const char* new_poi_list = "n_M2_H4_3b:n_M3_H4_3b:n_M4_H4_3b",
                                   const char* new_poi_name = "n_M234_H4_3b",
                                   int npoiPoints = 20,
                                   double poiMinVal = 0.,
                                   double poiMaxVal = 20.,
                                   double constraintWidth = 1.5,
                                   double ymax = 10.,
                                   int verbLevel=0 ) {


     gStyle->SetOptStat(0) ;

     //--- make output directory.

     //// TString outputdir( wsfile ) ;
     //// outputdir.ReplaceAll("rootfiles/","outputfiles/scans-") ;
     //// outputdir.ReplaceAll(".root","") ;

     char command[10000] ;
     sprintf( command, "basename %s", wsfile ) ;
     TString wsfilenopath = gSystem->GetFromPipe( command ) ;
     wsfilenopath.ReplaceAll(".root","") ;
     char outputdirstr[1000] ;
     sprintf( outputdirstr, "outputfiles/scans-%s", wsfilenopath.Data() ) ;
     TString outputdir( outputdirstr ) ;


     printf("\n\n Creating output directory: %s\n\n", outputdir.Data() ) ;
     sprintf(command, "mkdir -p %s", outputdir.Data() ) ;
     gSystem->Exec( command ) ;


     //--- Tell RooFit to shut up about anything less important than an ERROR.
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;



       if ( verbLevel > 0 ) { printf("\n\n Verbose level : %d\n\n", verbLevel) ; }


       TFile* wstf = new TFile( wsfile ) ;

       RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );

       if ( verbLevel > 0 ) { ws->Print() ; }






       RooDataSet* rds = (RooDataSet*) ws->obj( "ra2b_observed_rds" ) ;

       if ( verbLevel > 0 ) {
          printf("\n\n\n  ===== RooDataSet ====================\n\n") ;
          rds->Print() ;
          rds->printMultiline(cout, 1, kTRUE, "") ;
       }





       ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;
       RooAbsPdf* likelihood = modelConfig->GetPdf() ;

       RooRealVar* rrv_mu_susy_all0lep = ws->var("mu_susy_all0lep") ;
       if ( rrv_mu_susy_all0lep == 0x0 ) {
          printf("\n\n\n *** can't find mu_susy_all0lep in workspace.  Quitting.\n\n\n") ;
          return ;
       }





       //-- do BG only for now.
       rrv_mu_susy_all0lep->setVal(0.) ;
       rrv_mu_susy_all0lep->setConstant( kTRUE ) ;


  //// //-- Include floating susy
  //// rrv_mu_susy_all0lep->setVal(0.) ;
  //// rrv_mu_susy_all0lep->setConstant( kFALSE ) ;








       //-- do a prefit.

       printf("\n\n\n ====== Pre fit with unmodified nll var.\n\n") ;

       RooFitResult* dataFitResultSusyFixed = likelihood->fitTo(*rds, Save(true),Hesse(false),Minos(false),Strategy(1),PrintLevel(verbLevel));
       int dataSusyFixedFitCovQual = dataFitResultSusyFixed->covQual() ;
       if ( dataSusyFixedFitCovQual < 2 ) { printf("\n\n\n *** Failed fit!  Cov qual %d.  Quitting.\n\n", dataSusyFixedFitCovQual ) ; return ; }
       double dataFitSusyFixedNll = dataFitResultSusyFixed->minNll() ;

       if ( verbLevel > 0 ) {
          dataFitResultSusyFixed->Print("v") ;
       }

       printf("\n\n Nll value, from fit result : %.3f\n\n", dataFitSusyFixedNll ) ;

       delete dataFitResultSusyFixed ;






       //-- Construct the new POI parameter.
       RooAbsReal* new_poi_rar(0x0) ;

  //   new_poi_rar = ws->var( new_poi_name ) ;
  //   if ( new_poi_rar == 0x0 ) {
  //      printf("\n\n New POI %s is not a variable.  Trying function.\n\n", new_poi_name ) ;
  //      new_poi_rar = ws->function( new_poi_name ) ;
  //      if ( new_poi_rar == 0x0 ) {
  //         printf("\n\n New POI %s is not a function.  I quit.\n\n", new_poi_name ) ;
  //         return ;
  //      }
  //   } else {
  //      printf("\n\n     New POI %s is a variable with current value %.1f.\n\n", new_poi_name, new_poi_rar->getVal() ) ;
  //   }



       TString new_poi_list_ts( new_poi_list ) ;
       TObjArray* poi_obs_name_list = new_poi_list_ts.Tokenize(":") ;

       int n_poi_observables = poi_obs_name_list -> GetEntries() ;

       printf(" New poi list : %d : %s\n", n_poi_observables, new_poi_list ) ;


       RooAbsReal* rar_new_poi_obs[100] ;

       char new_poi_formula_string[10000] ;
       sprintf( new_poi_formula_string, "@0" ) ;

       RooArgSet ras_poi_list ;

       for ( int pi=0; pi<n_poi_observables; pi++ ) {
          TObjString* os = (TObjString*) (*poi_obs_name_list)[pi] ;
          TString obs_str = os->GetString() ;
          rar_new_poi_obs[pi] = ws -> function( obs_str.Data() ) ;
          if ( rar_new_poi_obs[pi] == 0x0 ) {
             rar_new_poi_obs[pi] = ws -> var( obs_str.Data() ) ;
             if ( rar_new_poi_obs[pi] == 0x0 ) {
                printf("\n\n *** %s missing from workspace.\n\n", obs_str.Data() ) ;
                return ;
             }
          }
          if ( pi>0 ) {
             char buffer[10000] ;
             sprintf( buffer, "%s+@%d", new_poi_formula_string, pi ) ;
             sprintf( new_poi_formula_string, "%s", buffer ) ;
          }
          ras_poi_list.add( *rar_new_poi_obs[pi] ) ;
          printf(" %2d : %s  %.1f\n", pi, obs_str.Data(), rar_new_poi_obs[pi]->getVal() ) ;
       } // pi

       new_poi_rar = new RooFormulaVar( new_poi_name, new_poi_formula_string, ras_poi_list ) ;
       new_poi_rar -> Print() ;





       if ( npoiPoints <=0 ) {
          printf("\n\n Quitting now.\n\n" ) ;
          return ;
       }


       double startPoiVal = new_poi_rar->getVal() ;



      //--- The RooNLLVar is NOT equivalent to what minuit uses.
  //   RooNLLVar* nll = new RooNLLVar("nll","nll", *likelihood, *rds ) ;
  //   printf("\n\n Nll value, from construction : %.3f\n\n", nll->getVal() ) ;

      //--- output of createNLL IS what minuit uses, so use that.
       RooAbsReal* nll = likelihood -> createNLL( *rds, Verbose(true) ) ;

       RooRealVar* rrv_poiValue = new RooRealVar( "poiValue", "poiValue", 0., -10000., 10000. ) ;
   /// rrv_poiValue->setVal( poiMinVal ) ;
   /// rrv_poiValue->setConstant(kTRUE) ;

       RooRealVar* rrv_constraintWidth = new RooRealVar("constraintWidth","constraintWidth", 0.1, 0.1, 1000. ) ;
       rrv_constraintWidth -> setVal( constraintWidth ) ;
       rrv_constraintWidth -> setConstant(kTRUE) ;




       if ( verbLevel > 0 ) {
          printf("\n\n ======= debug likelihood print\n\n") ;
          likelihood->Print("v") ;
          printf("\n\n ======= debug nll print\n\n") ;
          nll->Print("v") ;
       }






    //----------------------------------------------------------------------------------------------

       RooMinuit* rminuit( 0x0 ) ;

       RooFormulaVar* plot_var( 0x0 ) ;

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


       char minuit_formula_unbiased_unconstrained[100000] ;
       sprintf( minuit_formula_unbiased_unconstrained, "%s+%s", nll->GetName(), rfv_ignorePdfTerm->GetName() ) ;
       RooFormulaVar* rfv_minuitvar_unbiased_unconstrained = new RooFormulaVar( "minuitvar_unbiased_unconstrained",
         minuit_formula_unbiased_unconstrained,
         RooArgList( *nll, *rfv_ignorePdfTerm ) ) ;

       RooMinuit* rminuit_ub_uc = new RooMinuit( *rfv_minuitvar_unbiased_unconstrained  ) ;

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
       while ( RooRealVar* par = (RooRealVar*) floatParIter->Next() ) {
          sprintf( floatParName[nFloatParInitVal], "%s", par->GetName() ) ;
          floatParInitVal[nFloatParInitVal] = par->getVal() ;
          nFloatParInitVal++ ;
       }


       printf("\n\n Unbiased best value for new POI %s is : %7.1f\n\n", new_poi_rar->GetName(), new_poi_rar->getVal() ) ;

       double best_unbiased_poi_val = new_poi_rar->getVal() ;

     //-------


       char minuit_formula[10000] ;
       sprintf( minuit_formula, "%s+%s*(%s-%s)*(%s-%s)+%s",
         nll->GetName(),
         rrv_constraintWidth->GetName(),
         new_poi_rar->GetName(), rrv_poiValue->GetName(),
         new_poi_rar->GetName(), rrv_poiValue->GetName(),
         rfv_ignorePdfTerm->GetName() ) ;

       printf("\n\n Creating new minuit variable with formula: %s\n\n", minuit_formula ) ;
       RooFormulaVar* new_minuit_var = new RooFormulaVar("new_minuit_var", minuit_formula,
           RooArgList( *nll,
                       *rrv_constraintWidth,
                       *new_poi_rar, *rrv_poiValue,
                       *new_poi_rar, *rrv_poiValue,
                       *rfv_ignorePdfTerm
                       ) ) ;

       printf("\n\n Current value is %.2f\n\n",
            new_minuit_var->getVal() ) ;

       rminuit = new RooMinuit( *new_minuit_var ) ;

       char plotvar_formula[10000] ;
       sprintf( plotvar_formula, "%s+%s",
              nll->GetName(),
              rfv_ignorePdfTerm->GetName()
              ) ;

       printf("\n\n Creating new plot variable with formula: %s\n\n", plotvar_formula ) ;
       plot_var = new RooFormulaVar("plot_var", plotvar_formula,
           RooArgList( *nll,
                       *rfv_ignorePdfTerm
                       ) ) ;

       printf("\n\n Current value is %.2f\n\n",
            plot_var->getVal() ) ;




       rminuit->setPrintLevel(verbLevel-1) ;
       if ( verbLevel <=0 ) { rminuit->setNoWarn() ; }

    //----------------------------------------------------------------------------------------------

       //-- If POI range is -1 to -1, automatically determine the range using the set value.

       if ( poiMinVal < 0. && poiMaxVal < 0. ) {

          printf("\n\n Automatic determination of scan range.\n\n") ;

          if ( startPoiVal <= 0. ) {
             printf("\n\n *** POI starting value zero or negative %g.  Quit.\n\n\n", startPoiVal ) ;
             return ;
          }

          poiMinVal = startPoiVal - 3.5 * sqrt(startPoiVal) ;
          poiMaxVal = startPoiVal + 6.0 * sqrt(startPoiVal) ;

          if ( poiMinVal < 0. ) { poiMinVal = 0. ; }

          printf("    Start val = %g.   Scan range:   %g  to  %g\n\n", startPoiVal, poiMinVal, poiMaxVal ) ;


       }



    //----------------------------------------------------------------------------------------------


       double poiVals_scanDown[1000] ;
       double nllVals_scanDown[1000] ;

       //-- Do scan down from best value.

       printf("\n\n +++++ Starting scan down from best value.\n\n") ;

       double minNllVal(1.e9) ;

       for ( int poivi=0; poivi < npoiPoints/2 ; poivi++ ) {

          ////double poiValue = poiMinVal + poivi*(poiMaxVal-poiMinVal)/(1.*(npoiPoints-1)) ;
          double poiValue = best_unbiased_poi_val - poivi*(best_unbiased_poi_val-poiMinVal)/(1.*(npoiPoints/2-1)) ;

          rrv_poiValue -> setVal( poiValue ) ;
          rrv_poiValue -> setConstant( kTRUE ) ;


       //+++++++++++++++++++++++++++++++++++

          rminuit->migrad() ;
          rminuit->hesse() ;
          RooFitResult* rfr = rminuit->save() ;

       //+++++++++++++++++++++++++++++++++++


          if ( verbLevel > 0 ) { rfr->Print("v") ; }


          float fit_minuit_var_val = rfr->minNll() ;

          printf(" %02d : poi constraint = %.2f : allvars : MinuitVar, createNLL, PV, POI :    %.5f   %.5f   %.5f   %.5f\n",
                poivi, rrv_poiValue->getVal(), fit_minuit_var_val, nll->getVal(), plot_var->getVal(), new_poi_rar->getVal() ) ;
          cout << flush ;



       // TString binstring( new_poi_name ) ;
       // binstring.ReplaceAll("n_","") ;
       // TString trigbinstring( binstring ) ;
       // trigbinstring.Resize(5) ;

          //char pname[1000] ;
          //RooAbsReal* par(0x0) ;


          poiVals_scanDown[poivi] = new_poi_rar->getVal() ;
          nllVals_scanDown[poivi] = plot_var->getVal() ;

          if ( nllVals_scanDown[poivi] < minNllVal ) { minNllVal = nllVals_scanDown[poivi] ; }

          delete rfr ;


       } // poivi


      //-- Refit for best unbiased value.


       printf("\n\n +++++ Resetting floats to best unbiased fit values.\n\n") ;

       for ( int pi=0; pi<nFloatParInitVal; pi++ ) {
          RooRealVar* par = ws->var( floatParName[pi] ) ;
          par->setVal( floatParInitVal[pi] ) ;
       } // pi.

       printf("\n\n +++++ Starting scan up from best value.\n\n") ;

      //-- Now do scan up.

       double poiVals_scanUp[1000] ;
       double nllVals_scanUp[1000] ;

       for ( int poivi=0; poivi < npoiPoints/2 ; poivi++ ) {

          double poiValue = best_unbiased_poi_val + poivi*(poiMaxVal-best_unbiased_poi_val)/(1.*(npoiPoints/2-1)) ;
      //  printf(" best_unbiased_poi_val=%g, poivi=%d, poiMaxVal=%g, npoiPoints=%d, poiValue=%g\n",
      //      best_unbiased_poi_val, poivi, poiMaxVal, npoiPoints, poiValue ) ;

          rrv_poiValue -> setVal( poiValue ) ;
          rrv_poiValue -> setConstant( kTRUE ) ;

      //  printf("  debug rrv_poiValue = %g\n", rrv_poiValue->getVal() ) ;


       //+++++++++++++++++++++++++++++++++++

          rminuit->migrad() ;
          rminuit->hesse() ;
          RooFitResult* rfr = rminuit->save() ;

       //+++++++++++++++++++++++++++++++++++


          if ( verbLevel > 0 ) { rfr->Print("v") ; }


          float fit_minuit_var_val = rfr->minNll() ;

          printf(" %02d : poi constraint = %.2f : allvars : MinuitVar, createNLL, PV, POI :    %.5f   %.5f   %.5f   %.5f\n",
                poivi, rrv_poiValue->getVal(), fit_minuit_var_val, nll->getVal(), plot_var->getVal(), new_poi_rar->getVal() ) ;
          cout << flush ;



       // TString binstring( new_poi_name ) ;
       // binstring.ReplaceAll("n_","") ;
       // TString trigbinstring( binstring ) ;
       // trigbinstring.Resize(5) ;

          //char pname[1000] ;
          //RooAbsReal* par(0x0) ;


          poiVals_scanUp[poivi] = new_poi_rar->getVal() ;
          nllVals_scanUp[poivi] = plot_var->getVal() ;

          if ( nllVals_scanUp[poivi] < minNllVal ) { minNllVal = nllVals_scanUp[poivi] ; }

          delete rfr ;


       } // poivi





       double poiVals[1000] ;
       double nllVals[1000] ;

       int pointCount(0) ;
       for ( int pi=0; pi<npoiPoints/2; pi++ ) {
          poiVals[pi] = poiVals_scanDown[(npoiPoints/2-1)-pi] ;
          nllVals[pi] = nllVals_scanDown[(npoiPoints/2-1)-pi] ;
          pointCount++ ;
       }
       for ( int pi=1; pi<npoiPoints/2; pi++ ) {
          poiVals[pointCount] = poiVals_scanUp[pi] ;
          nllVals[pointCount] = nllVals_scanUp[pi] ;
          pointCount++ ;
       }
       npoiPoints = pointCount ;

       printf("\n\n --- TGraph arrays:\n") ;
       for ( int i=0; i<npoiPoints; i++ ) {
          printf("  %2d : poi = %6.1f, nll = %g\n", i, poiVals[i], nllVals[i] ) ;
       }
       printf("\n\n") ;

       double nllDiffVals[1000] ;

       double poiAtMinlnL(-1.) ;
       double poiAtMinusDelta2(-1.) ;
       double poiAtPlusDelta2(-1.) ;
       for ( int poivi=0; poivi < npoiPoints ; poivi++ ) {
          nllDiffVals[poivi] = 2.*(nllVals[poivi] - minNllVal) ;
          double poiValue = poiMinVal + poivi*(poiMaxVal-poiMinVal)/(1.*npoiPoints) ;
          if ( nllDiffVals[poivi] < 0.01 ) { poiAtMinlnL = poiValue ; }
          if ( poiAtMinusDelta2 < 0. && nllDiffVals[poivi] < 2.5 ) { poiAtMinusDelta2 = poiValue ; }
          if ( poiAtMinlnL > 0. && poiAtPlusDelta2 < 0. && nllDiffVals[poivi] > 2.0 ) { poiAtPlusDelta2 = poiValue ; }
       } // poivi

       printf("\n\n Estimates for poi at delta ln L = -2, 0, +2:  %g ,   %g ,   %g\n\n", poiAtMinusDelta2, poiAtMinlnL, poiAtPlusDelta2 ) ;




      //--- Main canvas

       TCanvas* cscan = (TCanvas*) gDirectory->FindObject("cscan") ;
       if ( cscan == 0x0 ) {
          printf("\n Creating canvas.\n\n") ;
          cscan = new TCanvas("cscan","Delta nll") ;
       }


       char gname[1000] ;

       TGraph* graph = new TGraph( npoiPoints, poiVals, nllDiffVals ) ;
       sprintf( gname, "scan_%s", new_poi_name ) ;
       graph->SetName( gname ) ;

       double poiBest(-1.) ;
       double poiMinus1stdv(-1.) ;
       double poiPlus1stdv(-1.) ;
       double poiMinus2stdv(-1.) ;
       double poiPlus2stdv(-1.) ;
       double twoDeltalnLMin(1e9) ;

       int nscan(1000) ;
       for ( int xi=0; xi<nscan; xi++ ) {

          double x = poiVals[0] + xi*(poiVals[npoiPoints-1]-poiVals[0])/(nscan-1) ;

          double twoDeltalnL = graph -> Eval( x, 0, "S" ) ;

          if ( poiMinus1stdv < 0. && twoDeltalnL < 1.0 ) { poiMinus1stdv = x ; printf(" set m1 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}
          if ( poiMinus2stdv < 0. && twoDeltalnL < 4.0 ) { poiMinus2stdv = x ; printf(" set m2 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}
          if ( twoDeltalnL < twoDeltalnLMin ) { poiBest = x ; twoDeltalnLMin = twoDeltalnL ; }
          if ( twoDeltalnLMin < 0.3 && poiPlus1stdv < 0. && twoDeltalnL > 1.0 ) { poiPlus1stdv = x ; printf(" set p1 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}
          if ( twoDeltalnLMin < 0.3 && poiPlus2stdv < 0. && twoDeltalnL > 4.0 ) { poiPlus2stdv = x ; printf(" set p2 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}

          if ( xi%100 == 0 ) { printf( " %4d : poi=%6.2f,  2DeltalnL = %6.2f\n", xi, x, twoDeltalnL ) ; }

       }
       printf("\n\n POI estimate :  %g  +%g  -%g    [%g,%g],   two sigma errors: +%g  -%g   [%g,%g]\n\n",
               poiBest,
               (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv), poiMinus1stdv, poiPlus1stdv,
               (poiPlus2stdv-poiBest), (poiBest-poiMinus2stdv), poiMinus2stdv, poiPlus2stdv
               ) ;

       printf(" %s val,pm1sig,pm2sig: %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n",
          new_poi_name, poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv), (poiPlus2stdv-poiBest), (poiBest-poiMinus2stdv) ) ;

       char htitle[1000] ;
       sprintf(htitle, "%s profile likelihood scan: -2ln(L/Lm)", new_poi_name ) ;
       TH1F* hscan = new TH1F("hscan", htitle, 10, poiMinVal, poiMaxVal ) ;
       hscan->SetMinimum(0.) ;
       hscan->SetMaximum(ymax) ;


       hscan->DrawCopy() ;
       graph->SetLineColor(4) ;
       graph->SetLineWidth(3) ;
       graph->Draw("CP") ;
       gPad->SetGridx(1) ;
       gPad->SetGridy(1) ;
       cscan->Update() ;

       TLine* line = new TLine() ;
       line->SetLineColor(2) ;
       line->DrawLine(poiMinVal, 1., poiPlus1stdv, 1.) ;
       line->DrawLine(poiMinus1stdv,0., poiMinus1stdv, 1.) ;
       line->DrawLine(poiPlus1stdv ,0., poiPlus1stdv , 1.) ;

       TText* text = new TText() ;
       text->SetTextSize(0.04) ;
       char tstring[1000] ;

       sprintf( tstring, "%s = %.1f +%.1f -%.1f", new_poi_name, poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv) ) ;
       text -> DrawTextNDC( 0.15, 0.85, tstring ) ;

       sprintf( tstring, "68%% interval [%.1f,  %.1f]", poiMinus1stdv, poiPlus1stdv ) ;
       text -> DrawTextNDC( 0.15, 0.78, tstring ) ;


       char hname[1000] ;
       sprintf( hname, "hscanout_%s", new_poi_name ) ;
       TH1F* hsout = new TH1F( hname,"scan results",4,0.,4.) ;
       double obsVal(-1.) ;
       hsout->SetBinContent(1, obsVal ) ;
       hsout->SetBinContent(2, poiPlus1stdv ) ;
       hsout->SetBinContent(3, poiBest ) ;
       hsout->SetBinContent(4, poiMinus1stdv ) ;
       TAxis* xaxis = hsout->GetXaxis() ;
       xaxis->SetBinLabel(1,"Observed val.") ;
       xaxis->SetBinLabel(2,"Model+1sd") ;
       xaxis->SetBinLabel(3,"Model") ;
       xaxis->SetBinLabel(4,"Model-1sd") ;

       char outrootfile[10000] ;
       sprintf( outrootfile, "%s/scan-hb-bs-%s.root", outputdir.Data(), new_poi_name ) ;

       char outpdffile[10000] ;
       sprintf( outpdffile, "%s/scan-hb-bs-%s.pdf", outputdir.Data(), new_poi_name ) ;

       cscan->Update() ; cscan->Draw() ;

       printf("\n Saving %s\n", outpdffile ) ;
       cscan->SaveAs( outpdffile ) ;



     //--- save in root file

       printf("\n Saving %s\n", outrootfile ) ;
       TFile fout(outrootfile,"recreate") ;
       graph->Write() ;
       hsout->Write() ;
       fout.Close() ;

       delete ws ;
       wstf->Close() ;

   }

