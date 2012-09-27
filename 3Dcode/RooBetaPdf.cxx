/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooBetaPdf.cxx,v 1.1 2012/07/20 13:22:18 kreis Exp $
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

#include "RooBetaPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "Math/DistFunc.h"

ClassImp(RooBetaPdf)


//_____________________________________________________________________________
RooBetaPdf::RooBetaPdf(const char *name, const char *title,
			 RooAbsReal& _x, RooAbsReal& _alpha,
			 RooAbsReal& _beta) :
  RooAbsPdf(name,title),
  x("x","Observable",this,_x),
  alpha("alpha","Alpha",this,_alpha),
  beta("beta","Width",this,_beta)
{
}



//_____________________________________________________________________________
RooBetaPdf::RooBetaPdf(const RooBetaPdf& other, const char* name) : 
  RooAbsPdf(other,name), x("x",this,other.x), alpha("alpha",this,other.alpha),
  beta("beta",this,other.beta)
{
}



//_____________________________________________________________________________
Double_t RooBetaPdf::evaluate() const
{
  if(x<=0 || x>=1) return 0;
  return ROOT::Math::beta_pdf(x,alpha,beta);
}



//_____________________________________________________________________________
Int_t RooBetaPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  if (matchArgs(allVars,analVars,x)) return 1 ;
  return 0 ;
}



//_____________________________________________________________________________
Double_t RooBetaPdf::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  assert(code==1) ;
  if( x.min(rangeName) <=  0. && x.max(rangeName) >= 1.) return 1.;
  if( x.min(rangeName) <=  0.) return ROOT::Math::beta_cdf( x.max(rangeName) , alpha , beta ) ;
  if( x.max(rangeName) >=  1.) return (1. - ROOT::Math::beta_cdf( x.min(rangeName) , alpha , beta ) ) ;
  return ( ROOT::Math::beta_cdf( x.max(rangeName) , alpha , beta ) - ROOT::Math::beta_cdf( x.min(rangeName) , alpha , beta ) ) ;
}

//_____________________________________________________________________________
Int_t RooBetaPdf::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  if (matchArgs(directVars,generateVars,x)) return 1 ;  
  return 0 ;
}

//_____________________________________________________________________________
void RooBetaPdf::generateEvent(Int_t code)
{
  assert(code==1) ;
  double xgen,x1,x2;
  while(true)
    {
      x1 = generateGamma(alpha);
      x2 = generateGamma(beta);
      xgen = x1/(x1+x2);
      if (xgen<x.max() && xgen>x.min()) {
	x = xgen ;
	break;
      }
    }
  return;
}

//_____________________________________________________________________________
double RooBetaPdf::generateGamma(const double& a)
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
  
