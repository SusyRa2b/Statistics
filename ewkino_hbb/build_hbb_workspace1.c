
#include "TMath.h"
#include "TSystem.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooUniform.h"
#include "RooStats/ModelConfig.h"


#include "getFileValue.c"

   using namespace RooFit ;
   using namespace RooStats ;

   const int bins_of_nb(3) ;
   const int max_bins_of_met(50) ;




   void build_hbb_workspace1( const char* infile = "outputfiles/input-file.txt", const char* outfile = "outputfiles/ws.root" ) {


    //-------------------------------------------------------------------------

     //-- Create workspace and other RooStats things.

      printf("\n\n Creating workspace.\n\n") ;

      RooWorkspace workspace("ws") ;
      workspace.autoImportClassCode(true) ;

   // RooArgSet* globalObservables      = new RooArgSet("globalObservables");
   // RooArgSet* allNuisances           = new RooArgSet("allNuisances");
   // RooArgSet* allNuisancePdfs        = new RooArgSet("allNuisancePdfs");
      RooArgSet* observedParametersList = new RooArgSet("observables") ;




    //-------------------------------------------------------------------------

      printf("\n\n Reading input file: %s\n\n", infile ) ;

      float fileVal ;
      char pname[1000] ;
      char formula[1000] ;


      sprintf( pname, "bins_of_met" ) ;
      if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
      int bins_of_met = TMath::Nint( fileVal ) ;

      //-- save bins_of_met in the workspace for convenience.
      RooRealVar bom( "bins_of_met", "bins_of_met", bins_of_met, 0., 1000. ) ;
      bom.setConstant(kTRUE) ;
      workspace.import(bom) ;


      //-- save bins_of_nb in the workspace for convenience.
      RooRealVar bonb( "bins_of_nb", "bins_of_nb", bins_of_nb, 0., 1000. ) ;
      bonb.setConstant(kTRUE) ;
      workspace.import(bonb) ;


      RooRealVar* rv_N_msig[bins_of_nb][max_bins_of_met] ; // first index is number of btags, second is met bin.
      RooRealVar* rv_N_msb[bins_of_nb][max_bins_of_met]  ; // first index is number of btags, second is met bin.

      RooRealVar* rv_smc_msig[bins_of_nb][max_bins_of_met] ; // first index is number of btags, second is met bin.
      RooRealVar* rv_smc_msb[bins_of_nb][max_bins_of_met]  ; // first index is number of btags, second is met bin.


      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            sprintf( pname, "N_%db_msig_met%d", nbi+2, mbi+1 ) ;
            if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
            rv_N_msig[nbi][mbi] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
            rv_N_msig[nbi][mbi] -> setVal( TMath::Nint(fileVal) ) ;
            rv_N_msig[nbi][mbi] -> setConstant( kTRUE ) ;
            observedParametersList -> add( *rv_N_msig[nbi][mbi] ) ;

            sprintf( pname, "N_%db_msb_met%d", nbi+2, mbi+1 ) ;
            if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
            rv_N_msb[nbi][mbi] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
            rv_N_msb[nbi][mbi] -> setVal( TMath::Nint(fileVal) ) ;
            rv_N_msb[nbi][mbi] -> setConstant( kTRUE ) ;
            observedParametersList -> add( *rv_N_msb[nbi][mbi] ) ;

            sprintf( pname, "smc_%db_msig_met%d", nbi+2, mbi+1 ) ;
            if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
            rv_smc_msig[nbi][mbi] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
            rv_smc_msig[nbi][mbi] -> setVal( TMath::Nint(fileVal) ) ;
            rv_smc_msig[nbi][mbi] -> setConstant( kTRUE ) ;

            sprintf( pname, "smc_%db_msb_met%d", nbi+2, mbi+1 ) ;
            if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
            rv_smc_msb[nbi][mbi] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
            rv_smc_msb[nbi][mbi] -> setVal( TMath::Nint(fileVal) ) ;
            rv_smc_msb[nbi][mbi] -> setConstant( kTRUE ) ;

         } // mbi.

      } // nbi.

     //-- Finished reading input from file.

    //-------------------------------------------------------------------------

      printf("\n\n Creating and importing dataset into workspace.\n\n") ;

      RooDataSet* dsObserved = new RooDataSet("hbb_observed_rds", "hbb observed data values", *observedParametersList ) ;
      dsObserved -> add( *observedParametersList ) ;
      workspace.import( *dsObserved ) ;

    //-------------------------------------------------------------------------

     //-- Define all floats.

      printf("\n\n Defining all unconstrained floats (Ratios, signal strength).\n\n") ;

      double R_msigmsb_initialval(0.15) ;

      RooRealVar* rv_R_msigmsb[50] ;

      for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

         sprintf( pname, "R_msigmsb_met%d", mbi+1 ) ;
         printf( "  %s\n", pname ) ;
         rv_R_msigmsb[mbi] = new RooRealVar( pname, pname, R_msigmsb_initialval, 0., 3. ) ;
         rv_R_msigmsb[mbi] -> setConstant( kFALSE ) ;
         rv_R_msigmsb[mbi] -> Print() ;

      } // mbi.

      printf("\n") ;

      sprintf( pname, "sig_strength" ) ;
      RooRealVar* rv_sig_strength = new RooRealVar( pname, pname, 1.0, 0., 10. ) ;
      rv_sig_strength -> setConstant(kFALSE) ;
      rv_sig_strength -> Print() ;
      printf("  %s\n\n", pname ) ;

    //-------------------------------------------------------------------------

     //-- Define all mu parameters.

      printf("\n\n Defining mu parameters.\n\n") ;

      RooAbsReal* rv_mu_bg_msig[bins_of_nb][max_bins_of_met] ;  // first index is number of btags, second is met bin.
      RooAbsReal* rv_mu_bg_msb[bins_of_nb][max_bins_of_met]  ;  // first index is number of btags, second is met bin.

      RooAbsReal* rv_mu_sig_msig[bins_of_nb][max_bins_of_met] ; // first index is number of btags, second is met bin.
      RooAbsReal* rv_mu_sig_msb[bins_of_nb][max_bins_of_met]  ; // first index is number of btags, second is met bin.

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            sprintf( pname, "mu_bg_%db_msb_met%d", nbi+2, mbi+1 ) ;
            printf( "  %s\n", pname ) ;
            rv_mu_bg_msb[nbi][mbi] = new RooRealVar( pname, pname, rv_N_msb[nbi][mbi] -> getVal(), 0., 1.e6 ) ;
            rv_mu_bg_msb[nbi][mbi] -> Print() ;


            sprintf( formula, "@0 * @1" ) ;

            sprintf( pname, "mu_bg_%db_msig_met%d", nbi+2, mbi+1 ) ;
            printf( "  %s\n", pname ) ;
            rv_mu_bg_msig[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_R_msigmsb[mbi], *rv_mu_bg_msb[nbi][mbi] ) ) ;
            rv_mu_bg_msig[nbi][mbi] -> Print() ;

            sprintf( pname, "mu_sig_%db_msig_met%d", nbi+2, mbi+1 ) ;
            printf( "  %s\n", pname ) ;
            rv_mu_sig_msig[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_sig_strength, *rv_smc_msig[nbi][mbi] ) ) ;
            rv_mu_sig_msig[nbi][mbi] -> Print() ;

            sprintf( pname, "mu_sig_%db_msb_met%d", nbi+2, mbi+1 ) ;
            printf( "  %s\n", pname ) ;
            rv_mu_sig_msb[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_sig_strength, *rv_smc_msb[nbi][mbi] ) ) ;
            rv_mu_sig_msb[nbi][mbi] -> Print() ;


         } // mbi.

      } // nbi.

     //-- Finished defining mu parameters.

    //-------------------------------------------------------------------------

     //-- Defining small n's

     printf("\n\n Defining small n's.\n\n") ;

     RooAbsReal* rv_n_msig[bins_of_nb][max_bins_of_met] ;  // first index is number of btags, second is met bin.
     RooAbsReal* rv_n_msb[bins_of_nb][max_bins_of_met]  ;  // first index is number of btags, second is met bin.

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            sprintf( formula, "@0 + @1" ) ;

            sprintf( pname, "n_%db_msig_met%d", nbi+2, mbi+1 ) ;
            printf( "  %s\n", pname ) ;
            rv_n_msig[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_mu_sig_msig[nbi][mbi], *rv_mu_bg_msig[nbi][mbi] ) ) ;
            rv_n_msig[nbi][mbi] -> Print() ;
            workspace.import( *rv_n_msig[nbi][mbi] ) ;

            sprintf( pname, "n_%db_msb_met%d", nbi+2, mbi+1 ) ;
            printf( "  %s\n", pname ) ;
            rv_n_msb[nbi][mbi] = new RooFormulaVar( pname, formula, RooArgSet( *rv_mu_sig_msb[nbi][mbi], *rv_mu_bg_msb[nbi][mbi] ) ) ;
            rv_n_msb[nbi][mbi] -> Print() ;
            workspace.import( *rv_n_msb[nbi][mbi] ) ;

         } // mbi.

      } // nbi.

    //-------------------------------------------------------------------------

     //-- Define the Poisson pdfs for the observables.

      printf("\n\n Defining Poisson pdfs for the observables.\n\n") ;

      RooAbsReal* rv_pdf_msig[bins_of_nb][max_bins_of_met] ;  // first index is number of btags, second is met bin.
      RooAbsReal* rv_pdf_msb[bins_of_nb][max_bins_of_met]  ;  // first index is number of btags, second is met bin.

      RooArgSet pdflist ;

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {

            sprintf( pname, "pdf_%db_msig_met%d", nbi+2, mbi+1 ) ;
            printf( "  %s\n", pname ) ;
            rv_pdf_msig[nbi][mbi] = new RooPoisson( pname, pname, *rv_N_msig[nbi][mbi], *rv_n_msig[nbi][mbi] ) ;
            rv_pdf_msig[nbi][mbi] -> Print() ;

            pdflist.add( *rv_pdf_msig[nbi][mbi] ) ;

            sprintf( pname, "pdf_%db_msb_met%d", nbi+2, mbi+1 ) ;
            printf( "  %s\n", pname ) ;
            rv_pdf_msb[nbi][mbi] = new RooPoisson( pname, pname, *rv_N_msb[nbi][mbi], *rv_n_msb[nbi][mbi] ) ;
            rv_pdf_msb[nbi][mbi] -> Print() ;

            pdflist.add( *rv_pdf_msb[nbi][mbi] ) ;

         } // mbi.

      } // nbi.

    //-------------------------------------------------------------------------

     //-- Build the likelihood.

      printf("\n\n Building the likelihood.\n\n") ;

      pdflist.Print() ;
      printf("\n") ;

      RooProdPdf* likelihood = new RooProdPdf( "likelihood", "hbb likelihood", pdflist ) ;
      likelihood->Print() ;


    //-------------------------------------------------------------------------


  //  printf("\n\n Running a test fit.\n\n") ;


  //  dsObserved -> Print() ;
  //  dsObserved -> printMultiline(cout, 1, kTRUE, "") ;


  //  printf("\n\n =============================================\n\n") ;
  //  likelihood -> fitTo( *dsObserved, PrintLevel(3), Hesse(0), Minos(0) ) ;
  //  printf("\n\n =============================================\n\n") ;







     //-- Set up RooStats models.

      printf("\n\n Setting up S+B model.\n\n") ;

      RooArgSet poi( *rv_sig_strength, "poi" ) ;
      RooUniform signal_prior( "signal_prior", "signal_prior", *rv_sig_strength ) ;

      ModelConfig sbModel ("SbModel");
      sbModel.SetWorkspace( workspace ) ;
      sbModel.SetPdf( *likelihood ) ;
      sbModel.SetParametersOfInterest( poi );
      sbModel.SetPriorPdf(signal_prior);
      sbModel.SetObservables( *observedParametersList );
   // sbModel.SetNuisanceParameters( *allNuisances );
   // sbModel.SetGlobalObservables( *globalObservables );

      workspace.Print() ;

      printf("\n\n Doing fit for S+B model.\n" ) ; fflush(stdout) ;

      RooAbsReal* pNll = sbModel.GetPdf()->createNLL(*dsObserved);
      RooAbsReal* pProfile = pNll->createProfile(RooArgSet());
      pProfile->getVal();
      RooArgSet* pPoiAndNuisance = new RooArgSet();
      pPoiAndNuisance->add(*sbModel.GetParametersOfInterest());
      if(sbModel.GetNuisanceParameters()) pPoiAndNuisance->add(*sbModel.GetNuisanceParameters());
      printf("\n\n Will save these parameter points that correspond to the fit to data.\n\n") ; fflush(stdout) ;
      pPoiAndNuisance->Print("v");
      sbModel.SetSnapshot(*pPoiAndNuisance);
      workspace.import (sbModel);

      delete pProfile ;
      delete pNll ;
      delete pPoiAndNuisance ;

      printf("\n\n Setting up BG-only model.\n\n") ;

      ModelConfig bModel (*(RooStats::ModelConfig *)workspace.obj("SbModel"));
      bModel.SetName("BModel");
      bModel.SetWorkspace(workspace);

      printf("\n\n Doing fit for BG-only model.\n" ) ; fflush(stdout) ;
      pNll = bModel.GetPdf()->createNLL(*dsObserved);
      pProfile = pNll->createProfile(*bModel.GetParametersOfInterest());
      ((RooRealVar *)(bModel.GetParametersOfInterest()->first()))->setVal(0.);
      pProfile->getVal();
      pPoiAndNuisance = new RooArgSet();
      pPoiAndNuisance->add(*bModel.GetParametersOfInterest());
      if(bModel.GetNuisanceParameters()) pPoiAndNuisance->add(*bModel.GetNuisanceParameters());
      printf("\n\n Should use these parameter points to generate pseudo data for bkg only.\n\n") ; fflush(stdout) ;
      pPoiAndNuisance->Print("v");
      bModel.SetSnapshot(*pPoiAndNuisance);
      workspace.import (bModel);

      delete pProfile ;
      delete pNll ;
      delete pPoiAndNuisance ;

      workspace.Print() ;

      printf("\n\n Saving workspace in : %s\n\n", outfile ) ;

      gSystem->Exec(" mkdir -p outputfiles " ) ;

      workspace.writeToFile( outfile ) ;




   } // build_hbb_workspace1.










