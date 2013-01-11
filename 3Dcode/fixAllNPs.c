
#include "RooArgList.h"
#include "TIterator.h"
#include "TSystem.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

  //==============================================================
  //
  //  The definition of a Nuisance Parameter (NP) for this
  //  code is any floating parameter that is NOT in the
  //  list of unconstrained parameters: unconstrained-pars.txt
  //
  //  Get the allFloats list from a RooFitResult using floatParsFinal()
  //
  //==============================================================


  bool fixAllNPs( const RooArgList& allFloats,
                  RooWorkspace* workspace,
                  const char* uncon_pars_file = "unconstrained-pars.txt",
                  bool resetValsToMean = true
                ) {

     printf("\n\n === Fixing all NPs.  Excluding pars listed in %s\n", uncon_pars_file ) ;
     printf("       Length of RooArgList of floats : %d\n\n", allFloats.getSize() ) ;


     int pcount(0) ;
     TIterator* parIter = allFloats.createIterator() ;
     while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {

        pcount++ ;
        char command[1000] ;
        sprintf( command, "grep -q %s %s", par->GetName(), uncon_pars_file ) ;
        int status = gSystem->Exec( command ) ;
        if ( status == 0 ) {
           printf(" %3d: %30s is unconstrained.  NOT fixing it.\n", pcount, par->GetName() ) ;
        } else {
           printf(" %3d: %30s is constrained.  Fixing it.\n", pcount, par->GetName() ) ;
           ////--- this doesn't work.  Need to dig the parameter out of workspace. ---
           ////  par->setConstant(kTRUE) ;
           ////-----------------------------------------------------------------------
           RooRealVar* rrv = workspace->var( par->GetName() ) ;
           if ( rrv == 0x0 ) {
              printf("\n\n *** fixAllNPs : can't find %s in workspace.\n\n", par->GetName() ) ;
              return false ;
           }
           if ( resetValsToMean ) {
              //-- find parameter that corresponds to the mean of the constraint pdf.
              char meanparname[1000] ;
              TString pname = par->GetName() ;
              if ( pname.Contains("prim") ) {
                 TString tsmpn = pname ;
                 tsmpn.ReplaceAll("prim_","prim_mean_") ;
                 sprintf( meanparname, "%s", tsmpn.Data() ) ;
              } else {
                 sprintf( meanparname, "mean_%s", par->GetName() ) ;
              }
              RooAbsReal* rar_mean = (RooAbsReal*) workspace->obj( meanparname ) ;
              if ( rar_mean == 0x0 ) {
                 printf("\n\n *** fixAllNPs : can't find mean parameter : %s\n\n", meanparname ) ;
                 return false ;
              }
              printf("     resetting %s to %g\n", rrv->GetName(), rar_mean->getVal() ) ;
              rrv->setVal( rar_mean->getVal() ) ;
           }
           rrv->setConstant(kTRUE) ;
        }

     } // par.

     printf("\n\n") ;

     return true ;

  } // fixAllNPs.

