
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

  //==============================================================================================

   void ws_n_susy_profile3D( const char* wsfile = "rootfiles/ws-data-unblind-1200-600.root",
                                   const char* bintag = "M4_H4_3b",
                                   const char* ignore_observable_name = "none",
                                   int npoiPoints = 20,
                                   double poiMinVal = 0.,
                                   double poiMaxVal = 8.,
                                   double constraintWidth = 40.,
                                   double ymax = 5.,
                                   int verbLevel=0
                                   ) {


     int metbin(-1), htbin(-1), nbbin(-1) ;
     sscanf( bintag, "M%d_H%d_%db", &metbin, &htbin, &nbbin ) ;
     if ( metbin < 0 || htbin < 0 || nbbin < 0 ) {
        printf("\n\n\n bad bintag: %s\n\n\n", bintag ) ; return ;
     }

     gStyle->SetOptStat(0) ;

     //--- make output directory.

     TString outputdir( wsfile ) ;
     outputdir.ReplaceAll("rootfiles/","outputfiles/scans-") ;
     outputdir.ReplaceAll(".root","") ;
     printf("\n\n Creating output directory: %s\n\n", outputdir.Data() ) ;
     char command[10000] ;
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





       rrv_mu_susy_all0lep->setVal(0.) ;
       rrv_mu_susy_all0lep->setConstant( kFALSE ) ;










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






       //-- construct the new POI parameter, which is the number of susy events in the given bin .


       char parname[1000] ;

       sprintf( parname, "mu_susy_%s", bintag ) ;
       RooAbsReal* rar_mu_susy = (RooAbsReal*) ws -> obj( parname ) ;
       if ( rar_mu_susy == 0x0 ) { printf("\n\n\n *** Can't find %s in workspace.\n\n", parname ) ; return ; }

       sprintf( parname, "eff_sf_%s", bintag ) ;
       RooAbsReal* rar_eff_sf = (RooAbsReal*) ws -> obj( parname ) ;
       if ( rar_eff_sf == 0x0 ) { printf("\n\n\n *** Can't find %s in workspace.\n\n", parname ) ; return ; }

       sprintf( parname, "all_gu" ) ;
       RooAbsReal* rar_all_gu = (RooAbsReal*) ws -> obj( parname ) ;
       if ( rar_all_gu == 0x0 ) { printf("\n\n\n *** Can't find %s in workspace.\n\n", parname ) ; return ; }

       sprintf( parname, "shapesyst_prod_%s", bintag ) ;
       RooAbsReal* rar_shapesyst_prod = (RooAbsReal*) ws -> obj( parname ) ;
       if ( rar_shapesyst_prod == 0x0 ) { printf("\n\n\n *** Can't find %s in workspace.\n\n", parname ) ; return ; }

       sprintf( parname, "trigeff_sl_M%d_H%d", metbin, htbin ) ;
       RooAbsReal* rar_trigeff_sl = (RooAbsReal*) ws -> obj( parname ) ;
       if ( rar_trigeff_sl == 0x0 ) { printf("\n\n\n *** Can't find %s in workspace.\n\n", parname ) ; return ; }

       sprintf( parname, "n_susy_%s", bintag ) ;
       RooFormulaVar* rfv_n_susy = new RooFormulaVar( parname, "@0 * @1 * @2 * @3 * @4",
            RooArgList( *rar_mu_susy, *rar_eff_sf, *rar_all_gu, *rar_shapesyst_prod, *rar_trigeff_sl )) ;

       printf("\n\n new POI is %s\n", parname ) ;
       rfv_n_susy -> Print() ;



       if ( npoiPoints <=0 ) {
          printf("\n\n Quitting now.\n\n" ) ;
          return ;
       }


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


       RooRealVar* rrv_obs(0x0) ;

       if ( strcmp( ignore_observable_name, "none" ) != 0 ) {

          char pdfName[1000] ;
          sprintf( pdfName, "pdf_%s", ignore_observable_name ) ;
          RooAbsPdf* pdf_ignore = ws -> pdf( pdfName ) ;

          if ( pdf_ignore == 0x0 ) {
             printf("\n\n *** Told to ignore %s but can't find pdf %s.\n\n", ignore_observable_name, pdfName ) ;
             return ;
          }

          rrv_obs = ws -> var( ignore_observable_name ) ;
          if ( rrv_obs == 0x0 ) {
             printf("\n\n *** Can't find RooRealVar for observable %s\n\n", ignore_observable_name ) ;
             return ;
          }
          rrv_obs->setConstant(kTRUE) ;


          printf("\n\n Ignoring observable %s by removing pdf %s\n\n", ignore_observable_name, pdfName ) ;

          char minuit_formula[10000] ;
          sprintf( minuit_formula, "%s+%s*(%s-%s)*(%s-%s)+log(%s)",
                 nll->GetName(),
                 rrv_constraintWidth->GetName(),
                 rfv_n_susy->GetName(), rrv_poiValue->GetName(),
                 rfv_n_susy->GetName(), rrv_poiValue->GetName(),
                 pdf_ignore->GetName()
                 ) ;

          printf("\n\n Creating new minuit variable with formula: %s\n\n", minuit_formula ) ;
          RooFormulaVar* new_minuit_var = new RooFormulaVar("new_minuit_var", minuit_formula,
              RooArgList( *nll,
                          *rrv_constraintWidth,
                          *rfv_n_susy, *rrv_poiValue,
                          *rfv_n_susy, *rrv_poiValue,
                          *pdf_ignore
                          ) ) ;

          printf("\n\n Current value is %.2f\n\n",
               new_minuit_var->getVal() ) ;

          rminuit = new RooMinuit( *new_minuit_var ) ;




          char plotvar_formula[10000] ;
          sprintf( plotvar_formula, "%s+log(%s)",
                 nll->GetName(),
                 pdf_ignore->GetName()
                 ) ;

          printf("\n\n Creating new plot variable with formula: %s\n\n", plotvar_formula ) ;
          plot_var = new RooFormulaVar("plot_var", plotvar_formula,
              RooArgList( *nll,
                          *pdf_ignore
                          ) ) ;

          printf("\n\n Current value is %.2f\n\n",
               plot_var->getVal() ) ;


       } else {

          RooFormulaVar* new_minuit_var(0x0) ;

          char minuit_formula[10000] ;
          sprintf( minuit_formula, "%s+%s*(%s-%s)*(%s-%s)",
              nll->GetName(),
              rrv_constraintWidth->GetName(),
              rfv_n_susy->GetName(), rrv_poiValue->GetName(),
              rfv_n_susy->GetName(), rrv_poiValue->GetName()
              ) ;
          printf("\n\n Creating new minuit variable with formula: %s\n\n", minuit_formula ) ;
          new_minuit_var = new RooFormulaVar("new_minuit_var", minuit_formula,
              RooArgList( *nll,
                       *rrv_constraintWidth,
                       *rfv_n_susy, *rrv_poiValue,
                       *rfv_n_susy, *rrv_poiValue
                       ) ) ;


          printf("\n\n Current value is %.2f\n\n",
               new_minuit_var->getVal() ) ;

          rminuit = new RooMinuit( *new_minuit_var ) ;




          char plotvar_formula[10000] ;
          sprintf( plotvar_formula, "%s",
                 nll->GetName()
                 ) ;

          printf("\n\n Creating new plot variable with formula: %s\n\n", plotvar_formula ) ;
          plot_var = new RooFormulaVar("plot_var", plotvar_formula,
              RooArgList( *nll
                          ) ) ;

          printf("\n\n Current value is %.2f\n\n",
               plot_var->getVal() ) ;


       }


       rminuit->setPrintLevel(verbLevel-1) ;
       if ( verbLevel <=0 ) { rminuit->setNoWarn() ; }

 // //----------------------------------------------------------------------------------------------

 //    //-- If POI range is -1 to -1, automatically determine the range using the set value.

 //    if ( poiMinVal < 0. && poiMaxVal < 0. ) {

 //       printf("\n\n Automatic determination of scan range.\n\n") ;

 //       if ( startPoiVal <= 0. ) {
 //          printf("\n\n *** POI starting value zero or negative %g.  Quit.\n\n\n", startPoiVal ) ;
 //          return ;
 //       }

 //       poiMinVal = startPoiVal - 3.5 * sqrt(startPoiVal) ;
 //       poiMaxVal = startPoiVal + 6.0 * sqrt(startPoiVal) ;

 //       if ( poiMinVal < 0. ) { poiMinVal = 0. ; }

 //       printf("    Start val = %g.   Scan range:   %g  to  %g\n\n", startPoiVal, poiMinVal, poiMaxVal ) ;


 //    }



 // //----------------------------------------------------------------------------------------------


       double poiVals[1000] ;
       double nllVals[1000] ;
       double minNllVal(1.e9) ;



       for ( int poivi=0; poivi < npoiPoints ; poivi++ ) {

          double poiValue = poiMinVal + poivi*(poiMaxVal-poiMinVal)/(1.*(npoiPoints-1)) ;

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
                poivi, rrv_poiValue->getVal(), fit_minuit_var_val, nll->getVal(), plot_var->getVal(), rfv_n_susy->getVal() ) ;
          cout << flush ;

          poiVals[poivi] = rfv_n_susy->getVal() ;
          nllVals[poivi] = plot_var->getVal() ;

          if ( nllVals[poivi] < minNllVal ) { minNllVal = nllVals[poivi] ; }

          delete rfr ;


       } // poivi

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
       sprintf( gname, "scan_%s", rfv_n_susy->GetName() ) ;
       graph->SetName( gname ) ;


       double poiBest(-1.) ;
       double poiMinus1stdv(-1.) ;
       double poiPlus1stdv(-1.) ;
       double twoDeltalnLMin(1e9) ;

       int nscan(1000) ;
       for ( int xi=0; xi<nscan; xi++ ) {

          double x = poiVals[0] + xi*(poiVals[npoiPoints-1]-poiVals[0])/(nscan-1) ;

          double twoDeltalnL = graph -> Eval( x ) ;

          if ( poiMinus1stdv < 0. && twoDeltalnL < 1.0 ) { poiMinus1stdv = x ; printf(" set m1 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}
          if ( twoDeltalnL < twoDeltalnLMin ) { poiBest = x ; twoDeltalnLMin = twoDeltalnL ; }
          if ( twoDeltalnLMin < 0.3 && poiPlus1stdv < 0. && twoDeltalnL > 1.0 ) { poiPlus1stdv = x ; printf(" set p1 : %d, x=%g, 2dnll=%g\n", xi, x, twoDeltalnL) ;}

          if ( xi%10 == 0 ) { printf( " %4d : poi=%6.2f,  2DeltalnL = %6.2f\n", xi, x, twoDeltalnL ) ; }

       }
       printf("\n\n POI estimate :  %g  +%g  -%g    [%g,%g]\n\n",
               poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv), poiMinus1stdv, poiPlus1stdv ) ;

       char htitle[1000] ;
       sprintf(htitle, "%s profile likelihood scan: -2ln(L/Lm)", rfv_n_susy->GetName() ) ;
       TH1F* hscan = new TH1F("hscan", htitle, 10, poiMinVal, poiMaxVal ) ;
       hscan->SetMinimum(0.) ;
       hscan->SetMaximum(ymax) ;


       hscan->Draw() ;
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

       sprintf( tstring, "N susy in %s = %.1f +%.1f -%.1f", bintag, poiBest, (poiPlus1stdv-poiBest), (poiBest-poiMinus1stdv) ) ;
       text -> DrawTextNDC( 0.15, 0.85, tstring ) ;

       sprintf( tstring, "68%% interval [%.1f,  %.1f]", poiMinus1stdv, poiPlus1stdv ) ;
       text -> DrawTextNDC( 0.15, 0.78, tstring ) ;





       char hname[1000] ;
       sprintf( hname, "hscanout_%s", rfv_n_susy->GetName() ) ;
       TH1F* hsout = new TH1F( hname,"scan results",4,0.,4.) ;
       double obsVal(-1.) ;
       if ( rrv_obs != 0 ) { obsVal = rrv_obs->getVal() ; }
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
       sprintf( outrootfile, "%s/scan-%s.root", outputdir.Data(), rfv_n_susy->GetName() ) ;

       char outpdffile[10000] ;
       sprintf( outpdffile, "%s/scan-%s.pdf", outputdir.Data(), rfv_n_susy->GetName() ) ;

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

