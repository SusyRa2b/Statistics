/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooBetaPrimePdf.cxx,v 1.1 2012/07/20 13:22:48 kreis Exp $
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
}



//_____________________________________________________________________________
RooBetaPrimePdf::RooBetaPrimePdf(const RooBetaPrimePdf& other, const char* name) : 
  RooAbsPdf(other,name), x("x",this,other.x), alpha("alpha",this,other.alpha),
  beta("beta",this,other.beta)
{
}



//_____________________________________________________________________________
Double_t RooBetaPrimePdf::evaluate() const
{
  if(x<=0) return 0;
  return pow(x,alpha-1)*pow(1+x,-alpha-beta)/ROOT::Math::beta(alpha,beta);
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
  return 0 ;
}

//_____________________________________________________________________________
void RooBetaPrimePdf::generateEvent(Int_t code)
{
  assert(code==1) ;
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
  
