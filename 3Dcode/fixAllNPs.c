
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
                  const char* uncon_pars_file = "datfiles/unconstrained-pars.txt",
                  bool resetValsToMean = true
                ) {

     printf("\n\n === Fixing all NPs.  Excluding pars listed in %s\n\n", uncon_pars_file ) ;

     TIterator* parIter = allFloats.createIterator() ;
     while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {

        char command[1000] ;
        sprintf( command, "grep -q %s %s", par->GetName(), uncon_pars_file ) ;
        int status = gSystem->Exec( command ) ;
        if ( status == 0 ) {
           printf(" %30s is unconstrained.  NOT fixing it.\n", par->GetName() ) ;
        } else {
           printf(" %30s is constrained.  Fixing it.\n", par->GetName() ) ;
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
              TString pname = par->GetName() ;
           }
           rrv->setConstant(kTRUE) ;
        }

     } // par.

     printf("\n\n") ;

     return true ;

  } // fixAllNPs.

