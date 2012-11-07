#include <iostream>
#include <string.h>
#include <complex>
#include <map>

#include "TString.h"
#include "TPRegexp.h"

#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"

#include "RooGaussian.h"

#include "RooFormulaVar.h"
#include "RooRatio.h"
#include "../3Dcode/RooPosDefCorrGauss.h"

using namespace RooFit ;


RooAbsArg* getGaussianConstraint( RooWorkspace& ws, TString NP_name, const TString binName,
				    const double NP_val, const double NP_err,
				    const TString observables, const TString nuisances,
				    bool allowNegative = false ) {
  
  NP_name = NP_name + binName;

  RooAbsArg* constrained_check = ws.arg(NP_name);
  if(constrained_check) return constrained_check;

  if ( NP_err <= 0. ) {
    //RooRealVar* rrc = new RooRealVar( NP_name+"_zeroErrorGaussian", NP_name+"_zeroErrorGaussian", NP_val ) ;
    RooRealVar* rrc = new RooRealVar( NP_name, NP_name, NP_val ) ;
    rrc->setConstant();
    ws.import(*rrc);
    return ws.var(rrc->GetName());
  }
  
  double max = NP_val + 6.*NP_err ;
  double min = NP_val - 6.*NP_err ;

  if ( min < 0. && !allowNegative ) { min = 1e-5 ; }

  RooRealVar* np_rrv = new RooRealVar(NP_name, NP_name, min, max ) ;
  np_rrv -> setVal( NP_val ) ;
  np_rrv -> setConstant( kFALSE ) ;

  RooConstVar* g_mean = new RooConstVar( "mean_"+NP_name, "mean_"+NP_name, NP_val ) ;
  RooConstVar* g_sigma = new RooConstVar( "sigma_"+NP_name, "sigma_"+NP_name, NP_err ) ;
  RooGaussian* np_pdf = new RooGaussian( "pdf_"+NP_name, "pdf_"+NP_name, *np_rrv, *g_mean, *g_sigma ) ;

  //allNuisances -> add( *np_rrv ) ;//BEN FIXME -- add to nuisances?
  //allNuisancePdfs -> add( *np_pdf ) ;//BEN FIXME -- add to nuisances?
  ws.import( *np_rrv );
  ws.import( *np_pdf );
  
  return np_rrv ;


}



RooAbsReal* getCorrelatedGaussianConstraint( RooWorkspace& ws, TString NP_name, const TString binName,
					      double NP_val, double NP_err, 
					      const TString observables, const TString nuisances,
					      const TString NP_base_name, 
					      bool changeSign = false, bool allowNegative = false) {
  
  NP_name = NP_name + binName;
  
  RooAbsReal* constrained_check = ws.var(NP_name);
  if(constrained_check) return constrained_check;

  if ( NP_err <= 0. ) {
    //RooRealVar* rrc = new RooRealVar( NP_name+"_zeroErrorCorrelatedGaussian", NP_name+"_zeroErrorCorrelatedGaussian", NP_val ) ;
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
    
    //allNuisances -> add( *rrv_np_base_par ) ;//BEN FIXME -- add to nuisances
    ws.import( *rrv_np_base_par ) ;

    RooConstVar* g_mean = new RooConstVar( "mean_"+NP_base_name, "mean_"+NP_base_name, 0.0 ) ;
    RooConstVar* g_sigma = new RooConstVar( "sigma_"+NP_base_name, "sigma_"+NP_base_name, 1.0 ) ;

    RooGaussian* base_np_pdf = new RooGaussian( "pdf_"+NP_base_name, "pdf_"+NP_base_name, *rrv_np_base_par, *g_mean, *g_sigma ) ;
    
    //allNuisancePdfs -> add( *base_np_pdf ) ;//BEN FIXME -- add to nuisances 
    ws.import( *base_np_pdf );

  }

  //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.
  RooConstVar* g_mean = new RooConstVar( "mean_"+NP_name, "mean_"+NP_name, NP_val ) ;
  RooConstVar* g_sigma = new RooConstVar( "sigma_"+NP_name, "sigma_"+NP_name, NP_err ) ;

  RooAbsReal* rar(0x0) ;

  if ( allowNegative ) {

    char formula[1000] ;

    if ( !changeSign ) {
      sprintf( formula, "@0+@1*@2" ) ;
    } else {
      sprintf( formula, "@0-@1*@2" ) ;
    }

    rar = new RooFormulaVar( NP_name, formula, RooArgSet( *g_mean, *g_sigma, *rrv_np_base_par ) ) ;

  } else {
    const char* NP_name_char = NP_name.Data();
    rar = new RooPosDefCorrGauss( NP_name_char, NP_name_char, *g_mean, *g_sigma, *rrv_np_base_par, changeSign ) ;
  }

  ws.import(*rar, RecycleConflictNodes());
  cout << " returning corr gauss " << endl;
  return rar ;
  
}




