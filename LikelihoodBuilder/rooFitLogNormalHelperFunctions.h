#include <iostream>

#include "TString.h"

#include "RooAbsArg.h"
#include "RooWorkspace.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooGaussian.h"
#include "RooArgSet.h"

using namespace RooFit;

RooAbsArg* getLognormalConstraint( RooWorkspace& ws, TString NP_name, const TString binName, 
				    double NP_val, double NP_err,
				    const TString observables, const TString nuisances) {
  
  NP_name = NP_name + binName;

  RooAbsArg* constrained_check = ws.arg(NP_name);
  if(constrained_check) return constrained_check;
  
  if ( NP_err <= 0. ) {
    RooRealVar* rrc = new RooRealVar( NP_name, NP_name, NP_val ) ;
    rrc->setConstant();
    ws.import(*rrc);
    return ws.var(rrc->GetName());
  }
  
  RooRealVar* np_prim_rrv = new RooRealVar( "prim_"+NP_name, "prim_"+NP_name, 0., -6., 6. ) ;
  np_prim_rrv -> setVal( 0. ) ;
  np_prim_rrv -> setConstant( kFALSE ) ;

  RooRealVar* np_prim_mean = new RooRealVar( "prim_mean_"+NP_name, "prim_mean_"+NP_name, 0., -6., 6. ) ;
  np_prim_mean->setConstant(kTRUE) ;
  
  RooConstVar* np_prim_sigma = new RooConstVar( "prim_sigma_"+NP_name, "prim_sigma_"+NP_name, 1. ) ;
  
  RooGaussian* np_prim_pdf = new RooGaussian( "pdf_prim_"+NP_name, "pdf_prim_"+NP_name, *np_prim_rrv, *np_prim_mean, *np_prim_sigma ) ;

  ws.import( *np_prim_rrv );
  ws.import( *np_prim_pdf, RecycleConflictNodes() );
  //allNuisances -> add( *np_prim_rrv ) ;
  //allNuisancePdfs -> add( *np_prim_pdf ) ;
  //globalObservables -> add( *np_prim_mean ) ;

  //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

  RooConstVar* g_mean  = new RooConstVar( "mean_"+NP_name, "mean_"+NP_name, NP_val ) ;
  RooConstVar* g_sigma = new RooConstVar( "sigma_"+NP_name, "sigma_"+NP_name, NP_err ) ;

  //-- compute the log-normal-distributed parameter from the primary parameter.

  RooFormulaVar* np_rfv = new RooFormulaVar( NP_name, "@0 * pow( exp( @1/@0 ), @2)",
					     RooArgSet( *g_mean, *g_sigma, *np_prim_rrv ) ) ;

  
  ws.import( *np_rfv , RecycleConflictNodes() );
  return np_rfv ;


} // makeLognormalConstraint.




RooAbsReal* getCorrelatedLogNormalConstraint( RooWorkspace& ws, TString NP_name, TString binName, 
					       double NP_val, double NP_err, 
					       const TString observables, const TString nuisances,
					       const TString NP_base_name, 
					       bool changeSign ) {
  
  NP_name = NP_name + binName;
  
  RooAbsReal* constrained_check = ws.var(NP_name);
  if(constrained_check) return constrained_check;
  
  if ( NP_err <= 0. ) {
    RooRealVar* rrc = new RooRealVar( NP_name, NP_name, NP_val ) ;
    rrc->setConstant();
    ws.import(*rrc);
    return ws.var(rrc->GetName());
  }
  
  RooRealVar* rrv_np_base_par = ws.var( NP_base_name ) ;
  
  if ( rrv_np_base_par == 0x0 ) {
    
    rrv_np_base_par = new RooRealVar( NP_base_name, NP_base_name, -6.0, 6.0 ) ;
    rrv_np_base_par -> setVal( 0. ) ;
    rrv_np_base_par -> setConstant( kFALSE ) ;

    ws.import( *rrv_np_base_par );
    //allNuisances -> add( *rrv_np_base_par ) ;

    RooRealVar* g_mean = new RooRealVar( "mean_"+NP_base_name, "mean_"+NP_base_name, 0.0,-10.,10. ) ;
    g_mean->setConstant(kTRUE);
    RooConstVar* g_sigma = new RooConstVar( "sigma_"+NP_base_name, "sigma_"+NP_base_name, 1.0 ) ;
    
    RooGaussian* base_np_pdf = new RooGaussian( "pdf_"+NP_base_name, "pdf_"+NP_base_name, *rrv_np_base_par, *g_mean, *g_sigma ) ;
    
    ws.import( *base_np_pdf, RecycleConflictNodes() );
    //allNuisancePdfs -> add( *base_np_pdf ) ;
    //globalObservables -> add( *g_mean ) ;
    
  }
  
  //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.
  
  RooConstVar* ln_mean  = new RooConstVar( "mean_"+NP_name, "mean_"+NP_name, NP_val ) ;
  
  RooConstVar* ln_sigma = new RooConstVar( "sigma_"+NP_name, "sigma_"+NP_name, NP_err ) ;
  
  RooAbsReal* rar(0x0) ;
  
  char formula[1000] ;
  
  if ( !changeSign ) {
    sprintf( formula, "@0 * pow( exp( @1/@0 ), @2 )" ) ;
  } else {
    sprintf( formula, "@0 * pow( exp( @1/@0 ), -1.0 * @2 )" ) ;
  }
  
  rar = new RooFormulaVar( NP_name, formula, RooArgSet( *ln_mean, *ln_sigma, *rrv_np_base_par ) ) ;
  
  ws.import ( *rar, RecycleConflictNodes() );

  return rar ;
  
} // makeCorrelatedLognormalConstraint.

