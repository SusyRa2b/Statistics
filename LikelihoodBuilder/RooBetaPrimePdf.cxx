/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooBetaPrimePdf.cxx,v 1.4 2012/10/12 14:12:49 kreis Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//
// BEGIN_HTML
// Plain BetaPdf p.d.f
// END_HTML
//

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "RooBetaPrimePdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "Math/DistFunc.h"
#include "Math/SpecFuncMathCore.h"

#include "betaHelperFunctions.h"

ClassImp(RooBetaPrimePdf)


//_____________________________________________________________________________
RooBetaPrimePdf::RooBetaPrimePdf(const char *name, const char *title,
			 RooAbsReal& _x, RooAbsReal& _alpha,
			 RooAbsReal& _beta) :
  RooAbsPdf(name,title),
  x("x","Observable",this,_x),
  alpha("alpha","Alpha",this,_alpha),
  beta("beta","Width",this,_beta)
{
  _alpha0 = alpha;
  _beta0 = beta;
  _x0 = -1;
  _sigma0 = sqrt( alpha * (alpha + beta - 1 ) / ( pow(beta - 1 , 2 ) * (beta - 2) ) );
}



//_____________________________________________________________________________
RooBetaPrimePdf::RooBetaPrimePdf(const RooBetaPrimePdf& other, const char* name) : 
  RooAbsPdf(other,name), x("x",this,other.x), alpha("alpha",this,other.alpha),
  beta("beta",this,other.beta)
{
  _alpha0 = alpha;
  _beta0 = beta;
  _x0 = other._x0;
  _sigma0 = other._sigma0;
}



//_____________________________________________________________________________
Double_t RooBetaPrimePdf::evaluate() const
{

  if (x<=0)
    {
      _logValue = log(0.0);
      return 0;
    }
  else {
    Double_t lnval = (alpha-1)*log(x) - (alpha + beta)*log(1+x) - (ROOT::Math::lgamma(alpha)+ROOT::Math::lgamma(beta)-ROOT::Math::lgamma(alpha+beta));
    //cout << "the log value is " << lnval << endl;
    _logValue = lnval;
    return exp(lnval);
  }

  //return pow(x,alpha-1)*pow(1+x,-alpha-beta)/ROOT::Math::beta(alpha,beta); // this version has problems with large alpha/beta
}



Double_t RooBetaPrimePdf::getValV(const RooArgSet* nset) const
{
  // Return current value, normalizated by integrating over
  // the observables in 'nset'. If 'nset' is 0, the unnormalized value. 
  // is returned. All elements of 'nset' must be lvalues
  //
  // Unnormalized values are not cached
  // Doing so would be complicated as _norm->getVal() could
  // spoil the cache and interfere with returning the cached
  // return value. Since unnormalized calls are typically
  // done in integration calls, there is no performance hit.

  // Fast-track processing of clean-cache objects
  //   if (_operMode==AClean) {
  //     cout << "RooAbsPdf::getValV(" << this << "," << GetName() << ") CLEAN  value = " << _value << endl ;
  //     return _value ;
  //   }

  // Special handling of case without normalization set (used in numeric integration of pdfs)
  if (!nset) {
    RooArgSet* tmp = _normSet ;
    _normSet = 0 ;
    Double_t val = evaluate() ;
    _normSet = tmp ;
    Bool_t error = traceEvalPdf(val) ;

    if (error) {
//       raiseEvalError() ;
      return 0 ;
    }
    return val ;
  }


  // Process change in last data set used
  Bool_t nsetChanged(kFALSE) ;
  if (nset!=_normSet || _norm==0) {
    nsetChanged = syncNormalization(nset) ;
  }

  // Return value of object. Calculated if dirty, otherwise cached value is returned.
  if (isValueDirty() || nsetChanged || _norm->isValueDirty()) {

    // Evaluate numerator
    Double_t rawVal = evaluate() ;
    Bool_t error = traceEvalPdf(rawVal) ; // Error checking and printing

    // Evaluate denominator
    Double_t normVal(_norm->getVal()) ;
    
    if (normVal<=0.) {
      error=kTRUE ;
      logEvalError("p.d.f normalization integral is zero or negative") ;  
    }

    // Raise global error flag if problems occur
    if (error) {
//       raiseEvalError() ;
      _value = 0 ;
      _logValue = log(_value);
    } else {
      _value = rawVal / normVal ;
      _logValue = _logValue - log(normVal);
//       cout << "RooAbsPdf::getValV(" << GetName() << ") writing _value = " << _value << endl ;
    }

    clearValueAndShapeDirty() ; //setValueDirty(kFALSE) ;
  } 

  return _value ;
}



Double_t RooBetaPrimePdf::getLogVal(const RooArgSet* nset) const 
{
  Double_t prob = getVal(nset) ;
  //evaluate();
  //cout << "returning log value is " << _logValue << endl;
  return _logValue;
  //if(prob < 0) {
  //
  //  logEvalError("getLogVal() top-level p.d.f evaluates to a negative number") ;
  //
  //  return 0;
  //}
  //if(prob == 0) {
  //
  //  logEvalError("getLogVal() top-level p.d.f evaluates to zero") ;
  //
  //  return log((double)0);
  //}
  //
  //if (TMath::IsNaN(prob)) {
  //  logEvalError("getLogVal() top-level p.d.f evaluates to NaN") ;
  //
  //  return log((double)0);
  //  
  //}
  //return log(prob);
}



//_____________________________________________________________________________
Int_t RooBetaPrimePdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  if (matchArgs(allVars,analVars,x)) return 1 ;
  return 0 ;
}



//_____________________________________________________________________________
Double_t RooBetaPrimePdf::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1) ;
  double retVal = 1.0;//ROOT::Math::beta(alpha,beta);
  if( x.min(rangeName) <=  0. && x.max(rangeName) == RooNumber::infinity()) return retVal;
  if( x.min(rangeName) <=  0.) return retVal*ROOT::Math::beta_cdf( x.max(rangeName) / (1+x.max(rangeName)), alpha , beta ) ;
  return ( retVal*(ROOT::Math::beta_cdf( x.max(rangeName) / (1+x.max(rangeName)), alpha , beta ) - ROOT::Math::beta_cdf( x.min(rangeName) / (1+x.min(rangeName)), alpha , beta ) ) ) ;
}

//_____________________________________________________________________________
Int_t RooBetaPrimePdf::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  if (matchArgs(directVars,generateVars,x)) return 1 ;  
  if (matchArgs(directVars,generateVars,alpha,beta)) return 2 ;
  return 0 ;
}

//_____________________________________________________________________________
void RooBetaPrimePdf::generateEvent(Int_t code)
{
  assert(code==1||code==2) ;
  if(code==1){
    double xgen,x1,x2;
    while(true)
      {
	x1 = generateGamma(alpha);
	x2 = generateGamma(beta);
	xgen = x1/x2;
	if (xgen<x.max() && xgen>x.min()) {
	  x = xgen ;
	  break;
	}
      }
    return;
  }
  if(code==2){
    if(_x0 != x) {
      //cout << "_x0 = " << _x0 << ", x = " << x << ", _sigma0 = " << _sigma0 << ", _alpha0 = " << _alpha0 << "_beta0 = " << _beta0 << endl;
      betaPrimeModeTransform(x , _sigma0 , _alpha0 , _beta0 );
      _x0 = x;
    }
    while(true)
      {
	double alphagen = RooRandom::randomGenerator()->Poisson(_alpha0);
	double betagen = RooRandom::randomGenerator()->Poisson(_beta0);
	
	if( alphagen<alpha.max() && alphagen>alpha.min() && betagen<beta.max() && betagen>beta.min() ) {
	  alpha = alphagen;
	  beta = betagen;
	  break;
	}
      }
    return;
  }
}

//_____________________________________________________________________________
double RooBetaPrimePdf::generateGamma(const double& a)
{
  //algorithm adapted from code example in:
  //Marsaglia, G. and Tsang, W. W.
  //A Simple Method for Generating Gamma Variables
  //ACM Transactions on Mathematical Software, Vol. 26, No. 3, September 2000
  //
  //The speed of this algorithm depends on the speed of generating normal variates.
  //The algorithm is limited to $\gamma \geq 0$ !

   /* assume a > 0 */

  //stolen from GSL

  if(a < 1) return generateGamma(a+1)*pow(RooRandom::randomGenerator()->Uniform(),1/a);

  double b = 1.0;
  double x_gen, v, u;
  double d = a - 1.0 / 3.0;
  double c = (1.0 / 3.0) / sqrt (d);
 
  while (1)
    {
      do
	{
	  x_gen = RooRandom::randomGenerator()->Gaus();
	  v = 1.0 + c * x_gen;
	}
      while (v <= 0);
      
      v = v * v * v;
      u = RooRandom::randomGenerator()->Uniform();
      
      if (u < 1 - 0.0331 * x_gen * x_gen * x_gen * x_gen) 
	break;
      
      if (log (u) < 0.5 * x_gen * x_gen + d * (1 - v + log (v)))
	break;
    }
  
  return b * d * v;
}
  
