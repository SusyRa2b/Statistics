
#include "TRandom2.h"
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

  //==============================================================================================

    //--- prototypes

    void printObservables() ;

    void saveObservables() ;

    void resetObservables() ;

    void generateObservables() ;


  //==============================================================================================

    //--- Global variables.


       //--- pointers to the observables.

       RooRealVar* rrv_Nsig(0x0) ;
       RooRealVar* rrv_Nsb(0x0) ;
       RooRealVar* rrv_Nsig_sl(0x0) ;
       RooRealVar* rrv_Nsb_sl(0x0) ;
       RooRealVar* rrv_Nsig_ldp(0x0) ;
       RooRealVar* rrv_Nsb_ldp(0x0) ;
       RooRealVar* rrv_Nsig_ee(0x0) ;
       RooRealVar* rrv_Nsb_ee(0x0) ;
       RooRealVar* rrv_Nsig_mm(0x0) ;
       RooRealVar* rrv_Nsb_mm(0x0) ;


       //--- pointers to likelihood model predictions for the observables.

       RooAbsReal* rfv_n_sig(0x0)       ;
       RooAbsReal* rfv_n_sb(0x0)        ;
       RooAbsReal* rfv_n_sig_sl(0x0)    ;
       RooAbsReal* rfv_n_sb_sl(0x0)     ;
       RooAbsReal* rfv_n_sig_ldp(0x0)   ;
       RooAbsReal* rfv_n_sb_ldp(0x0)    ;
       RooAbsReal* rfv_n_sig_ee(0x0)    ;
       RooAbsReal* rfv_n_sb_ee(0x0)     ;
       RooAbsReal* rfv_n_sig_mm(0x0)    ;
       RooAbsReal* rfv_n_sb_mm(0x0)     ;


       //--- values of the observables.

       int actual_Nsig(0) ;
       int actual_Nsb(0) ;
       int actual_Nsig_sl(0) ;
       int actual_Nsb_sl(0) ;
       int actual_Nsig_ldp(0) ;
       int actual_Nsb_ldp(0) ;
       int actual_Nsig_ee(0) ;
       int actual_Nsb_ee(0) ;
       int actual_Nsig_mm(0) ;
       int actual_Nsb_mm(0) ;



       //--- nuisance parameter initial values.

       int np_count(0) ;
       double np_initial_val[1000] ;




       TRandom2* random_ng ;




  //==============================================================================================

   void ws_cls_hybrid1( const char* wsfile = "output-files/expected-ws-lm9-2BL.root", bool isBgonlyStudy=false, double poiVal = 150.0, int nToys=100, bool makeTtree=true, int verbLevel=0 ) {



       TTree* toytt(0x0) ;
       TFile* ttfile(0x0) ;

       int    tt_gen_Nsig ;
       double tt_testStat ;
       double tt_dataTestStat ;
       double tt_hypo_mu_susy_sig ;
       char ttname[1000] ;
       char tttitle[1000] ;

       if ( makeTtree ) {

          ttfile = gDirectory->GetFile() ;
          if ( ttfile == 0x0 ) { printf("\n\n\n *** asked for a ttree but no open file???\n\n") ; return ; }


          if ( isBgonlyStudy ) {
             sprintf( ttname, "toytt_%.0f_bgo", poiVal ) ;
             sprintf( tttitle, "Toy study for background only, mu_susy_sig = %.0f", poiVal ) ;
          } else {
             sprintf( ttname, "toytt_%.0f_spb", poiVal ) ;
             sprintf( tttitle, "Toy study for signal+background, mu_susy_sig = %.0f", poiVal ) ;
          }

          printf("\n\n Creating TTree : %s : %s\n\n", ttname, tttitle ) ;

          gDirectory->pwd() ;
          gDirectory->ls() ;

          toytt = new TTree( ttname, tttitle ) ;

          gDirectory->ls() ;

          toytt -> Branch(  "gen_Nsig"         ,       &tt_gen_Nsig         ,      "gen_Nsig/I"         ) ;
          toytt -> Branch(  "testStat"         ,       &tt_testStat         ,      "testStat/D"         ) ;
          toytt -> Branch(  "dataTestStat"     ,       &tt_dataTestStat     ,      "dataTestStat/D"     ) ;
          toytt -> Branch(  "hypo_mu_susy_sig" ,       &tt_hypo_mu_susy_sig ,      "hypo_mu_susy_sig/D" ) ;

       }


     //--- Tell RooFit to shut up about anything less important than an ERROR.
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;


       random_ng = new TRandom2(12345) ;

   /// char sel[100] ;
   /// if ( strstr( wsfile, "1BL" ) != 0 ) {
   ///    sprintf( sel, "1BL" ) ;
   /// } else if ( strstr( wsfile, "2BL" ) != 0 ) {
   ///    sprintf( sel, "2BL" ) ;
   /// } else if ( strstr( wsfile, "3B" ) != 0 ) {
   ///    sprintf( sel, "3B" ) ;
   /// } else if ( strstr( wsfile, "1BT" ) != 0 ) {
   ///    sprintf( sel, "1BT" ) ;
   /// } else if ( strstr( wsfile, "2BT" ) != 0 ) {
   ///    sprintf( sel, "2BT" ) ;
   /// } else {
   ///    printf("\n\n\n *** can't figure out which selection this is.  I quit.\n\n" ) ;
   ///    return ;
   /// }
   /// printf("\n\n selection is %s\n\n", sel ) ;




       TFile* wstf = new TFile( wsfile ) ;

       RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );

       ws->Print() ;






       RooDataSet* rds = (RooDataSet*) ws->obj( "ra2b_observed_rds" ) ;
       printf("\n\n\n  ===== RooDataSet ====================\n\n") ;

       rds->Print() ;
       rds->printMultiline(cout, 1, kTRUE, "") ;





       ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;
       RooAbsPdf* likelihood = modelConfig->GetPdf() ;

       const RooArgSet* nuisanceParameters = modelConfig->GetNuisanceParameters() ;

       RooRealVar* rrv_mu_susy_sig = ws->var("mu_susy_sig") ;
       if ( rrv_mu_susy_sig == 0x0 ) {
          printf("\n\n\n *** can't find mu_susy_sig in workspace.  Quitting.\n\n\n") ;
          return ;
       }






 ////  printf("\n\n\n  ===== Doing a fit ====================\n\n") ;

 ////  RooFitResult* preFitResult = likelihood->fitTo( *rds, Save(true) ) ;
 ////  const RooArgList preFitFloatVals = preFitResult->floatParsFinal() ;
 ////  {
 ////    TIterator* parIter = preFitFloatVals.createIterator() ;
 ////    while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {
 ////       printf(" %20s : %8.2f\n", par->GetName(), par->getVal() ) ;
 ////    }
 ////  }







       //--- Get pointers to the model predictions of the observables.

       rfv_n_sig       = ws->function("n_sig") ;
       rfv_n_sb        = ws->function("n_sb") ;
       rfv_n_sig_sl    = ws->function("n_sig_sl") ;
       rfv_n_sb_sl     = ws->function("n_sb_sl") ;
       rfv_n_sig_ldp   = ws->function("n_sig_ldp") ;
       rfv_n_sb_ldp    = ws->function("n_sb_ldp") ;
       rfv_n_sig_ee    = ws->function("n_sig_ee") ;
       rfv_n_sb_ee     = ws->function("n_sb_ee") ;
       rfv_n_sig_mm    = ws->function("n_sig_mm") ;
       rfv_n_sb_mm     = ws->function("n_sb_mm") ;

       if ( rfv_n_sig         == 0x0 ) { printf("\n\n\n *** can't find n_sig       in workspace.  Quitting.\n\n\n") ; return ; }
       if ( rfv_n_sb          == 0x0 ) { printf("\n\n\n *** can't find n_sb        in workspace.  Quitting.\n\n\n") ; return ; }
       if ( rfv_n_sig_sl      == 0x0 ) { printf("\n\n\n *** can't find n_sig_sl    in workspace.  Quitting.\n\n\n") ; return ; }
       if ( rfv_n_sb_sl       == 0x0 ) { printf("\n\n\n *** can't find n_sb_sl     in workspace.  Quitting.\n\n\n") ; return ; }
       if ( rfv_n_sig_ldp     == 0x0 ) { printf("\n\n\n *** can't find n_sig_ldp   in workspace.  Quitting.\n\n\n") ; return ; }
       if ( rfv_n_sb_ldp      == 0x0 ) { printf("\n\n\n *** can't find n_sb_ldp    in workspace.  Quitting.\n\n\n") ; return ; }
       if ( rfv_n_sig_ee      == 0x0 ) { printf("\n\n\n *** can't find n_sig_ee    in workspace.  Quitting.\n\n\n") ; return ; }
       if ( rfv_n_sb_ee       == 0x0 ) { printf("\n\n\n *** can't find n_sb_ee     in workspace.  Quitting.\n\n\n") ; return ; }
       if ( rfv_n_sig_mm      == 0x0 ) { printf("\n\n\n *** can't find n_sig_mm    in workspace.  Quitting.\n\n\n") ; return ; }
       if ( rfv_n_sb_mm       == 0x0 ) { printf("\n\n\n *** can't find n_sb_mm     in workspace.  Quitting.\n\n\n") ; return ; }







       //--- Get pointers to the observables.

       const RooArgSet* dsras = rds->get() ;
       TIterator* obsIter = dsras->createIterator() ;
       while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {
          if ( strcmp( obs->GetName(), "Nsig"     ) == 0 ) { rrv_Nsig      = obs ; }
          if ( strcmp( obs->GetName(), "Nsb"      ) == 0 ) { rrv_Nsb       = obs ; }
          if ( strcmp( obs->GetName(), "Nsig_sl"  ) == 0 ) { rrv_Nsig_sl   = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_sl"   ) == 0 ) { rrv_Nsb_sl    = obs ; }
          if ( strcmp( obs->GetName(), "Nsig_ldp" ) == 0 ) { rrv_Nsig_ldp  = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_ldp"  ) == 0 ) { rrv_Nsb_ldp   = obs ; }
          if ( strcmp( obs->GetName(), "Nsig_ee"  ) == 0 ) { rrv_Nsig_ee   = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_ee"   ) == 0 ) { rrv_Nsb_ee    = obs ; }
          if ( strcmp( obs->GetName(), "Nsig_mm"  ) == 0 ) { rrv_Nsig_mm   = obs ; }
          if ( strcmp( obs->GetName(), "Nsb_mm"   ) == 0 ) { rrv_Nsb_mm    = obs ; }
       }

       if ( rrv_Nsig       == 0x0 ) { printf("\n\n\n *** can't find Nsig       in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb        == 0x0 ) { printf("\n\n\n *** can't find Nsb        in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsig_sl    == 0x0 ) { printf("\n\n\n *** can't find Nsig_sl    in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_sl     == 0x0 ) { printf("\n\n\n *** can't find Nsb_sl     in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsig_ldp   == 0x0 ) { printf("\n\n\n *** can't find Nsig_ldp   in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_ldp    == 0x0 ) { printf("\n\n\n *** can't find Nsb_ldp    in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsig_ee    == 0x0 ) { printf("\n\n\n *** can't find Nsig_ee    in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_ee     == 0x0 ) { printf("\n\n\n *** can't find Nsb_ee     in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsig_mm    == 0x0 ) { printf("\n\n\n *** can't find Nsig_mm    in dataset.  Quitting.\n\n\n") ; return ; }
       if ( rrv_Nsb_mm     == 0x0 ) { printf("\n\n\n *** can't find Nsb_mm     in dataset.  Quitting.\n\n\n") ; return ; }






       printf("\n\n\n === Model values for observables\n\n") ;

       printObservables() ;



      //--- save the actual values of the observables.

       saveObservables() ;











       //--- evaluate the test stat on the data: fit with susy floating.

       rrv_mu_susy_sig->setVal( poiVal ) ;
       rrv_mu_susy_sig->setConstant( kTRUE ) ;

       printf("\n\n\n ====== Fitting the data with susy fixed.\n\n") ;

       RooFitResult* dataFitResultSusyFixed = likelihood->fitTo(*rds, Save(true));
       int dataSusyFixedFitCovQual = dataFitResultSusyFixed->covQual() ;
       if ( dataSusyFixedFitCovQual != 3 ) { printf("\n\n\n *** Failed fit!  Cov qual %d.  Quitting.\n\n", dataSusyFixedFitCovQual ) ; return ; }
       double dataFitSusyFixedNll = dataFitResultSusyFixed->minNll() ;


       rrv_mu_susy_sig->setVal( 0.0 ) ;
       rrv_mu_susy_sig->setConstant( kFALSE ) ;

       printf("\n\n\n ====== Fitting the data with susy floating.\n\n") ;

       RooFitResult* dataFitResultSusyFloat = likelihood->fitTo(*rds, Save(true));
       int dataSusyFloatFitCovQual = dataFitResultSusyFloat->covQual() ;
       if ( dataSusyFloatFitCovQual != 3 ) { printf("\n\n\n *** Failed fit!  Cov qual %d.  Quitting.\n\n", dataSusyFloatFitCovQual ) ; return ; }
       double dataFitSusyFloatNll = dataFitResultSusyFloat->minNll() ;

       double dataTestStat = 2.*( dataFitSusyFixedNll - dataFitSusyFloatNll) ;

       printf("\n\n\n Data value of test stat : %8.2f\n", dataTestStat ) ;












       printf("\n\n\n === Nuisance parameters\n\n") ;

       {
          int npi(0) ;
          TIterator* npIter = nuisanceParameters->createIterator() ;
          while ( RooRealVar* np_rrv = (RooRealVar*) npIter->Next() ) {

             np_initial_val[npi] = np_rrv->getVal() ; //--- I am assuming that the order of the NPs in the iterator does not change.

             TString npname( np_rrv->GetName() ) ;
             npname.ReplaceAll("_prim","") ;
             RooAbsReal* np_rfv = ws->function( npname ) ;

             TString pdfname( np_rrv->GetName() ) ;
             pdfname.ReplaceAll("_prim","") ;
             pdfname.Prepend("pdf_") ;
             RooAbsPdf* np_pdf = ws->pdf( pdfname ) ;
             if ( np_pdf == 0x0 ) { printf("\n\n *** Can't find nuisance parameter pdf with name %s.\n\n", pdfname.Data() ) ; }

             if ( np_rfv != 0x0 ) {
                printf(" %20s : %8.2f , %20s, %8.2f\n", np_rrv->GetName(), np_rrv->getVal(), np_rfv->GetName(), np_rfv->getVal() ) ;
             } else {
                printf(" %20s : %8.2f\n", np_rrv->GetName(), np_rrv->getVal() ) ;
             }

             npi++ ;
          } // np_rrv iterator.

          np_count = npi ;

       }







       tt_dataTestStat = dataTestStat ;
       tt_hypo_mu_susy_sig = poiVal ;














       printf("\n\n\n === Doing the toys\n\n") ;

       int nToyOK(0) ;
       int nToyWorseThanData(0) ;

       for ( int ti=0; ti<nToys; ti++ ) {

          printf("\n\n\n ======= Toy %4d\n\n\n", ti ) ;





          //--- 1) pick values for the nuisance parameters from the PDFs and fix them.

          {
             TIterator* npIter = nuisanceParameters->createIterator() ;
             while ( RooRealVar* np_rrv = (RooRealVar*) npIter->Next() ) {

                TString pdfname( np_rrv->GetName() ) ;
                pdfname.ReplaceAll("_prim","") ;
                pdfname.Prepend("pdf_") ;
                RooAbsPdf* np_pdf = ws->pdf( pdfname ) ;
                if ( np_pdf == 0x0 ) { printf("\n\n *** Can't find nuisance parameter pdf with name %s.\n\n", pdfname.Data() ) ; return ; }

                RooDataSet* nprds = np_pdf->generate( RooArgSet(*np_rrv) ,1) ;
                const RooArgSet* npdsras = nprds->get() ;
                TIterator* valIter = npdsras->createIterator() ;
                RooRealVar* val = (RooRealVar*) valIter->Next() ;

                //--- reset the value of the nuisance parameter and fix it for the toy model definition fit.
                np_rrv->setVal( val->getVal() ) ;
                np_rrv->setConstant( kTRUE ) ;


                TString npname( np_rrv->GetName() ) ;
                npname.ReplaceAll("_prim","") ;
                RooAbsReal* np_rfv = ws->function( npname ) ;

                if ( verbLevel > 0 ) {
                   if ( np_rfv != 0x0 ) {
                      printf(" %20s : %8.2f , %15s, %8.3f\n", val->GetName(), val->getVal(), np_rfv->GetName(), np_rfv->getVal() ) ;
                   } else if ( strstr( npname.Data(), "eff_sf" ) != 0 ) {
                      np_rfv = ws->function( "eff_sf_sig" ) ;
                      RooAbsReal* np_rfv2 = ws->function( "eff_sf_sb" ) ;
                      printf(" %20s : %8.2f , %15s, %8.3f , %15s, %8.3f\n", val->GetName(), val->getVal(), np_rfv->GetName(), np_rfv->getVal(), np_rfv2->GetName(), np_rfv2->getVal() ) ;
                   } else if ( strstr( npname.Data(), "sf_ll" ) != 0 ) {
                      np_rfv = ws->function( "sf_ee" ) ;
                      RooAbsReal* np_rfv2 = ws->function( "sf_mm" ) ;
                      printf(" %20s : %8.2f , %15s, %8.3f , %15s, %8.3f\n", val->GetName(), val->getVal(), np_rfv->GetName(), np_rfv->getVal(), np_rfv2->GetName(), np_rfv2->getVal() ) ;
                   } else {
                      printf(" %20s : %8.2f\n", val->GetName(), val->getVal() ) ;
                   }
                }

                delete nprds ;

             } // np_rrv iterator
          }






          //--- 2) Fit the dataset with these values for the nuisance parameters.

          if ( isBgonlyStudy ) {
            //-- fit with susy yield fixed to zero.
             rrv_mu_susy_sig -> setVal( 0. ) ;
             if ( verbLevel > 0 ) { printf("\n Setting mu_susy_sig to zero.\n\n") ; }
          } else {
            //-- fit with susy yield fixed to predicted value.
             rrv_mu_susy_sig -> setVal( poiVal ) ;
             if ( verbLevel > 0 ) { printf("\n Setting mu_susy_sig to %8.1f.\n\n", poiVal) ; }
          }
          rrv_mu_susy_sig->setConstant( kTRUE ) ;

          if ( verbLevel > 0 ) {
             printf("\n\n") ;
             printf("  Fitting with these values for the observables to define the model for toy generation.\n") ;
             rds->printMultiline(cout, 1, kTRUE, "") ;
             printf("\n\n") ;
          }

          RooFitResult* toyModelDefinitionFitResult(0x0) ;
          if ( verbLevel < 2 ) {
             toyModelDefinitionFitResult = likelihood->fitTo(*rds, Save(true), PrintLevel(-1));
          } else {
             toyModelDefinitionFitResult = likelihood->fitTo(*rds, Save(true));
          }

          int toyModelDefFitCovQual = toyModelDefinitionFitResult->covQual() ;
          if ( verbLevel > 0 ) { printf("\n fit covariance matrix quality: %d\n\n", toyModelDefFitCovQual ) ; }
          if ( toyModelDefFitCovQual != 3 ) {
             printf("\n\n\n *** Bad toy model definition fit.  Cov qual %d.  Aborting this toy.\n\n\n", toyModelDefFitCovQual ) ;
             continue ;
          }

          delete toyModelDefinitionFitResult ;

          if ( verbLevel > 0 ) {
             printf("\n\n\n === Model values for observables.  These will be used to generate the toy dataset.\n\n") ;
             printObservables() ;
          }









          //--- 3) Generate a new set of observables based on this model.

          generateObservables() ;

          printf("\n\n\n   Generated dataset\n") ;
          rds->Print() ;
          rds->printMultiline(cout, 1, kTRUE, "") ;

          //--- Apparently, I need to make a new RooDataSet...  Resetting the
          //    values in the old one doesn't stick.  If you do likelihood->fitTo(*rds), it
          //    uses the original values, not the reset ones, in the fit.

          RooArgSet toyFitobservedParametersList ;
          toyFitobservedParametersList.add( *rrv_Nsig        ) ;
          toyFitobservedParametersList.add( *rrv_Nsb         ) ;
          toyFitobservedParametersList.add( *rrv_Nsig_sl     ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_sl      ) ;
          toyFitobservedParametersList.add( *rrv_Nsig_ldp    ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_ldp     ) ;
          toyFitobservedParametersList.add( *rrv_Nsig_ee     ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_ee      ) ;
          toyFitobservedParametersList.add( *rrv_Nsig_mm     ) ;
          toyFitobservedParametersList.add( *rrv_Nsb_mm      ) ;


          RooDataSet* toyFitdsObserved = new RooDataSet("toyfit_ra2b_observed_rds", "RA2b toy observed data values",
                                         toyFitobservedParametersList ) ;
          toyFitdsObserved->add( toyFitobservedParametersList ) ;





          //--- 4) Reset and free the nuisance parameters.

          {
             if ( verbLevel > 0 ) { printf("\n\n") ; }
             int npi(0) ;
             TIterator* npIter = nuisanceParameters->createIterator() ;
             while ( RooRealVar* np_rrv = (RooRealVar*) npIter->Next() ) {
                np_rrv -> setVal( np_initial_val[npi] ) ; // assuming that the order in the iterator does not change.
                np_rrv -> setConstant( kFALSE ) ;
                npi++ ;
                if ( verbLevel > 0 ) { printf("    reset %20s to %8.2f and freed it.\n", np_rrv->GetName() , np_rrv->getVal() ) ; }
             } // np_rrv iterator.
             if ( verbLevel > 0 ) { printf("\n\n") ; }
          }





          //--- 5a) Evaluate the test statistic: Fit with susy yield floating to get the absolute maximum log likelihood.

          if ( verbLevel > 0 ) { printf("\n\n  Evaluating the test statistic for this toy.  Fitting with susy floating.\n\n") ; }

          rrv_mu_susy_sig->setVal( 0.0 ) ;
          rrv_mu_susy_sig->setConstant( kFALSE ) ;

          if ( verbLevel > 0 ) {
             printf("\n toy dataset\n\n") ;
             toyFitdsObserved->printMultiline(cout, 1, kTRUE, "") ;
          }

     /////---- nfg.  Need to create a new dataset  ----------
     /////RooFitResult* maxLikelihoodFitResult = likelihood->fitTo(*rds, Save(true), PrintLevel(-1));
     /////RooFitResult* maxLikelihoodFitResult = likelihood->fitTo(*rds, Save(true));
     /////--------------

          RooFitResult* maxLikelihoodFitResult(0x0) ;
          if ( verbLevel < 2 ) {
             maxLikelihoodFitResult = likelihood->fitTo(*toyFitdsObserved, Save(true), PrintLevel(-1));
          } else {
             maxLikelihoodFitResult = likelihood->fitTo(*toyFitdsObserved, Save(true));
          }

          if ( verbLevel > 0 ) { printObservables() ; }

          int mlFitCovQual = maxLikelihoodFitResult->covQual() ;
          if ( verbLevel > 0 ) { printf("\n fit covariance matrix quality: %d , -log likelihood %f\n\n", mlFitCovQual, maxLikelihoodFitResult->minNll() ) ; }
          if ( mlFitCovQual != 3 ) {
             printf("\n\n\n *** Bad maximum likelihood fit (susy floating).  Cov qual %d.  Aborting this toy.\n\n\n", mlFitCovQual ) ;
             continue ;
          }
          double maxL_susyFloat = maxLikelihoodFitResult->minNll() ;
          double maxL_mu_susy_sig = rrv_mu_susy_sig->getVal() ;

          delete maxLikelihoodFitResult ;






          //--- 5b) Evaluate the test statistic: Fit with susy yield fixed to hypothesis value.
          //        This is only necessary if the maximum likelihood fit value of the susy yield
          //        is less than the hypothesis value (to get a one-sided limit).


          double testStat(0.0) ;
          double maxL_susyFixed(0.0) ;

          if ( maxL_mu_susy_sig < poiVal ) {

             if ( verbLevel > 0 ) { printf("\n\n  Evaluating the test statistic for this toy.  Fitting with susy fixed to %8.2f.\n\n", poiVal ) ; }

             rrv_mu_susy_sig->setVal( poiVal ) ;
             rrv_mu_susy_sig->setConstant( kTRUE ) ;

             if ( verbLevel > 0 ) {
                printf("\n toy dataset\n\n") ;
                rds->printMultiline(cout, 1, kTRUE, "") ;
             }

         ////--------- nfg.  need to make a new dataset  ---------------
         ////RooFitResult* susyFixedFitResult = likelihood->fitTo(*rds, Save(true), PrintLevel(-1));
         ////RooFitResult* susyFixedFitResult = likelihood->fitTo(*rds, Save(true));
         ////-----------------------------

             RooFitResult* susyFixedFitResult(0x0) ;
             if ( verbLevel < 2 ) {
                susyFixedFitResult = likelihood->fitTo(*toyFitdsObserved, Save(true), PrintLevel(-1));
             } else {
                susyFixedFitResult = likelihood->fitTo(*toyFitdsObserved, Save(true));
             }

             if ( verbLevel > 0 ) { printObservables() ; }

             int susyFixedFitCovQual = susyFixedFitResult->covQual() ;
             if ( verbLevel > 0 ) { printf("\n fit covariance matrix quality: %d , -log likelihood %f\n\n", susyFixedFitCovQual, susyFixedFitResult->minNll()  ) ; }
             if ( susyFixedFitCovQual != 3 ) {
                printf("\n\n\n *** Bad maximum likelihood fit (susy fixed).  Cov qual %d.  Aborting this toy.\n\n\n", susyFixedFitCovQual ) ;
                continue ;
             }
             maxL_susyFixed = susyFixedFitResult->minNll() ;
             testStat = 2. * (maxL_susyFixed - maxL_susyFloat) ;


             delete susyFixedFitResult ;


          } else {

             if ( verbLevel > 0 ) { printf("\n\n  Floating value of susy yield greater than hypo value (%8.2f > %8.2f).  Setting test stat to zero.\n\n", maxL_mu_susy_sig, poiVal ) ; }

             testStat = 0.0 ;

          }

          printf("   --- test stat for toy %4d : %8.2f\n", ti, testStat ) ;





          nToyOK++ ;

          if ( testStat > dataTestStat ) { nToyWorseThanData++ ; }


          if ( makeTtree ) {

             tt_testStat = testStat ;
             tt_gen_Nsig = rrv_Nsig->getVal() ;

             toytt->Fill() ;

          }





          //--- *) reset things for the next toy.

          resetObservables() ;

          delete toyFitdsObserved ;




       } // ti.

       wstf->Close() ;

       printf("\n\n\n") ;

       if ( nToyOK == 0 ) { printf("\n\n\n *** All toys bad !?!?!\n\n\n") ; return ; }

       double pValue = (1.0*nToyWorseThanData) / (1.0*nToyOK) ;

       if ( isBgonlyStudy ) {
          printf("\n\n\n p-value result, BG-only , poi=%3.0f : %4d / %4d = %6.3f\n\n\n\n", poiVal, nToyWorseThanData, nToyOK, pValue ) ;
       } else {
          printf("\n\n\n p-value result, S-plus-B, poi=%3.0f : %4d / %4d = %6.3f\n\n\n\n", poiVal, nToyWorseThanData, nToyOK, pValue ) ;
       }


       if ( makeTtree ) {
          printf("\n\n Writing TTree : %s : %s\n\n", ttname, tttitle ) ;
          ttfile->cd() ;
          toytt->Write() ;
       }


   } // ws_cls_hybrid1


  //==============================================================================================


   void printObservables() {

       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig     ->GetName(), rfv_n_sig     ->getVal(),  rrv_Nsig     ->GetName(), rrv_Nsig     ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb      ->GetName(), rfv_n_sb      ->getVal(),  rrv_Nsb      ->GetName(), rrv_Nsb      ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_sl  ->GetName(), rfv_n_sig_sl  ->getVal(),  rrv_Nsig_sl  ->GetName(), rrv_Nsig_sl  ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_sl   ->GetName(), rfv_n_sb_sl   ->getVal(),  rrv_Nsb_sl   ->GetName(), rrv_Nsb_sl   ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_ldp ->GetName(), rfv_n_sig_ldp ->getVal(),  rrv_Nsig_ldp ->GetName(), rrv_Nsig_ldp ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_ldp  ->GetName(), rfv_n_sb_ldp  ->getVal(),  rrv_Nsb_ldp  ->GetName(), rrv_Nsb_ldp  ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_ee  ->GetName(), rfv_n_sig_ee  ->getVal(),  rrv_Nsig_ee  ->GetName(), rrv_Nsig_ee  ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_ee   ->GetName(), rfv_n_sb_ee   ->getVal(),  rrv_Nsb_ee   ->GetName(), rrv_Nsb_ee   ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sig_mm  ->GetName(), rfv_n_sig_mm  ->getVal(),  rrv_Nsig_mm  ->GetName(), rrv_Nsig_mm  ->getVal() ) ;
       printf(" %20s : %8.2f ,  %10s : %8.0f\n", rfv_n_sb_mm   ->GetName(), rfv_n_sb_mm   ->getVal(),  rrv_Nsb_mm   ->GetName(), rrv_Nsb_mm   ->getVal() ) ;

   }


  //==============================================================================================



   void saveObservables() {

       actual_Nsig     =  rrv_Nsig     ->getVal()  ;
       actual_Nsb      =  rrv_Nsb      ->getVal()  ;
       actual_Nsig_sl  =  rrv_Nsig_sl  ->getVal()  ;
       actual_Nsb_sl   =  rrv_Nsb_sl   ->getVal()  ;
       actual_Nsig_ldp =  rrv_Nsig_ldp ->getVal()  ;
       actual_Nsb_ldp  =  rrv_Nsb_ldp  ->getVal()  ;
       actual_Nsig_ee  =  rrv_Nsig_ee  ->getVal()  ;
       actual_Nsb_ee   =  rrv_Nsb_ee   ->getVal()  ;
       actual_Nsig_mm  =  rrv_Nsig_mm  ->getVal()  ;
       actual_Nsb_mm   =  rrv_Nsb_mm   ->getVal()  ;

   }


  //==============================================================================================



   void resetObservables() {

       rrv_Nsig     ->setVal(actual_Nsig    )  ;
       rrv_Nsb      ->setVal(actual_Nsb     )  ;
       rrv_Nsig_sl  ->setVal(actual_Nsig_sl )  ;
       rrv_Nsb_sl   ->setVal(actual_Nsb_sl  )  ;
       rrv_Nsig_ldp ->setVal(actual_Nsig_ldp)  ;
       rrv_Nsb_ldp  ->setVal(actual_Nsb_ldp )  ;
       rrv_Nsig_ee  ->setVal(actual_Nsig_ee )  ;
       rrv_Nsb_ee   ->setVal(actual_Nsb_ee  )  ;
       rrv_Nsig_mm  ->setVal(actual_Nsig_mm )  ;
       rrv_Nsb_mm   ->setVal(actual_Nsb_mm  )  ;

   }


  //==============================================================================================


   void generateObservables() {

       rrv_Nsig     ->setVal( random_ng->Poisson( rfv_n_sig    ->getVal() ) )  ;
       rrv_Nsb      ->setVal( random_ng->Poisson( rfv_n_sb     ->getVal() ) )  ;
       rrv_Nsig_sl  ->setVal( random_ng->Poisson( rfv_n_sig_sl ->getVal() ) )  ;
       rrv_Nsb_sl   ->setVal( random_ng->Poisson( rfv_n_sb_sl  ->getVal() ) )  ;
       rrv_Nsig_ldp ->setVal( random_ng->Poisson( rfv_n_sig_ldp->getVal() ) )  ;
       rrv_Nsb_ldp  ->setVal( random_ng->Poisson( rfv_n_sb_ldp ->getVal() ) )  ;
       rrv_Nsig_ee  ->setVal( random_ng->Poisson( rfv_n_sig_ee ->getVal() ) )  ;
       rrv_Nsb_ee   ->setVal( random_ng->Poisson( rfv_n_sb_ee  ->getVal() ) )  ;
       rrv_Nsig_mm  ->setVal( random_ng->Poisson( rfv_n_sig_mm ->getVal() ) )  ;
       rrv_Nsb_mm   ->setVal( random_ng->Poisson( rfv_n_sb_mm  ->getVal() ) )  ;

   }

  //==============================================================================================





