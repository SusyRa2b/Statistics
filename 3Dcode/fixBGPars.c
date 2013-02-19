
#include "RooArgList.h"
#include "TIterator.h"
#include "TSystem.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"

  //==============================================================
  //
  //  Get the allFloats list from a RooFitResult using floatParsFinal()
  //
  //==============================================================


  bool fixBGPars( const RooArgList& allFloats,
                  RooWorkspace* workspace
                ) {

     printf("\n\n === Fixing all Background Parameters.\n\n") ;


     TIterator* parIter = allFloats.createIterator() ;
     while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {

        TString parname( par->GetName() ) ;

        if ( parname.Contains("mu_ttwj_sl_M") ||
             parname.Contains("mu_qcd_ldp_M") ||
             parname.Contains("mu_znn_M") ||
             parname.Contains("qcd_0lepLDP_ratio_H") ||
             parname.Contains("SFqcd_met2") ||
             parname.Contains("SFqcd_nb2") ||
             parname.Contains("ttwj_0lep1lep_ratio")
           ) {
           printf(" %30s is a BG parameter.  Fixing it.  Current value is %g\n", par->GetName(), par->getVal() ) ;
           ////--- this doesn't work.  Need to dig the parameter out of workspace. ---
           ////  par->setConstant(kTRUE) ;
           ////-----------------------------------------------------------------------
           RooRealVar* rrv = workspace->var( par->GetName() ) ;
           if ( rrv == 0x0 ) {
              printf("\n\n *** fixBGPars : can't find %s in workspace.\n\n", par->GetName() ) ;
              return false ;
           }
           rrv->setConstant(kTRUE) ;
        }

     } // par.

     printf("\n\n") ;

     return true ;

  } // fixBGPars.

