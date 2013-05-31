
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooStats/ModelConfig.h"

#include "TMath.h"

#include "getFileValue.c"


   void build_hbb_workspace1( const char* infile = "outputfiles/input-file.txt" ) {

      printf("\n\n Reading input file: %s\n\n", infile ) ;

      float fileVal ;
      char pname[1000] ;


      sprintf( pname, "bins_of_met" ) ;
      if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
      int bins_of_met = TMath::Nint( fileVal ) ;

      RooRealVar* rrv_N_4b_msig_met[50] ;
      RooRealVar* rrv_N_3b_msig_met[50] ;
      RooRealVar* rrv_N_2b_msig_met[50] ;
      RooRealVar* rrv_N_4b_msb_met[50] ;
      RooRealVar* rrv_N_3b_msb_met[50] ;
      RooRealVar* rrv_N_2b_msb_met[50] ;

      RooRealVar* rrv_smc_4b_msig_met[50] ;
      RooRealVar* rrv_smc_3b_msig_met[50] ;
      RooRealVar* rrv_smc_2b_msig_met[50] ;
      RooRealVar* rrv_smc_4b_msb_met[50] ;
      RooRealVar* rrv_smc_3b_msb_met[50] ;
      RooRealVar* rrv_smc_2b_msb_met[50] ;

      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {

         sprintf( pname, "N_4b_msig_met%d", mbi ) ;
         if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
         rrv_N_4b_msig_met[mbi-1] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rrv_N_4b_msig_met[mbi-1] -> setVal( fileVal ) ;
         rrv_N_4b_msig_met[mbi-1] -> setConstant( kTRUE ) ;

         sprintf( pname, "N_3b_msig_met%d", mbi ) ;
         if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
         rrv_N_3b_msig_met[mbi-1] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rrv_N_3b_msig_met[mbi-1] -> setVal( fileVal ) ;
         rrv_N_3b_msig_met[mbi-1] -> setConstant( kTRUE ) ;

         sprintf( pname, "N_2b_msig_met%d", mbi ) ;
         if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
         rrv_N_2b_msig_met[mbi-1] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rrv_N_2b_msig_met[mbi-1] -> setVal( fileVal ) ;
         rrv_N_2b_msig_met[mbi-1] -> setConstant( kTRUE ) ;


         sprintf( pname, "N_4b_msb_met%d", mbi ) ;
         if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
         rrv_N_4b_msb_met[mbi-1] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rrv_N_4b_msb_met[mbi-1] -> setVal( fileVal ) ;
         rrv_N_4b_msb_met[mbi-1] -> setConstant( kTRUE ) ;

         sprintf( pname, "N_3b_msb_met%d", mbi ) ;
         if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
         rrv_N_3b_msb_met[mbi-1] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rrv_N_3b_msb_met[mbi-1] -> setVal( fileVal ) ;
         rrv_N_3b_msb_met[mbi-1] -> setConstant( kTRUE ) ;

         sprintf( pname, "N_2b_msb_met%d", mbi ) ;
         if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
         rrv_N_2b_msb_met[mbi-1] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rrv_N_2b_msb_met[mbi-1] -> setVal( fileVal ) ;
         rrv_N_2b_msb_met[mbi-1] -> setConstant( kTRUE ) ;


         sprintf( pname, "smc_4b_msig_met%d", mbi ) ;
         if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
         rrv_smc_4b_msig_met[mbi-1] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rrv_smc_4b_msig_met[mbi-1] -> setVal( fileVal ) ;
         rrv_smc_4b_msig_met[mbi-1] -> setConstant( kTRUE ) ;

         sprintf( pname, "smc_3b_msig_met%d", mbi ) ;
         if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
         rrv_smc_3b_msig_met[mbi-1] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rrv_smc_3b_msig_met[mbi-1] -> setVal( fileVal ) ;
         rrv_smc_3b_msig_met[mbi-1] -> setConstant( kTRUE ) ;

         sprintf( pname, "smc_2b_msig_met%d", mbi ) ;
         if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
         rrv_smc_2b_msig_met[mbi-1] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rrv_smc_2b_msig_met[mbi-1] -> setVal( fileVal ) ;
         rrv_smc_2b_msig_met[mbi-1] -> setConstant( kTRUE ) ;


         sprintf( pname, "smc_4b_msb_met%d", mbi ) ;
         if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
         rrv_smc_4b_msb_met[mbi-1] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rrv_smc_4b_msb_met[mbi-1] -> setVal( fileVal ) ;
         rrv_smc_4b_msb_met[mbi-1] -> setConstant( kTRUE ) ;

         sprintf( pname, "smc_3b_msb_met%d", mbi ) ;
         if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
         rrv_smc_3b_msb_met[mbi-1] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rrv_smc_3b_msb_met[mbi-1] -> setVal( fileVal ) ;
         rrv_smc_3b_msb_met[mbi-1] -> setConstant( kTRUE ) ;

         sprintf( pname, "smc_2b_msb_met%d", mbi ) ;
         if ( !getFileValue( infile, pname, fileVal ) ) { printf("\n\n *** Error.  Can't find %s\n\n", pname ) ; return ; }
         rrv_smc_2b_msb_met[mbi-1] = new RooRealVar( pname, pname, 0., 1.e6 ) ;
         rrv_smc_2b_msb_met[mbi-1] -> setVal( fileVal ) ;
         rrv_smc_2b_msb_met[mbi-1] -> setConstant( kTRUE ) ;


      } // mbi.



   } // build_hbb_workspace1.

