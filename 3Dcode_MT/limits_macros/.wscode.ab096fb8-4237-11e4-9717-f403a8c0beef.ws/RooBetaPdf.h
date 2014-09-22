/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooBetaPdf.h,v 1.1 2012/09/27 19:29:31 owen Exp $
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
#ifndef ROO_BETAPDF
#define ROO_BETAPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooBetaPdf : public RooAbsPdf {
public:
  RooBetaPdf() {} ;
  RooBetaPdf(const char *name, const char *title,
	      RooAbsReal& _x, RooAbsReal& _alpha, RooAbsReal& _beta);
  RooBetaPdf(const RooBetaPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooBetaPdf(*this,newname); }
  inline virtual ~RooBetaPdf() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  void generateEvent(Int_t code);

protected:

  RooRealProxy x ;
  RooRealProxy alpha ;
  RooRealProxy beta ;
  
  Double_t evaluate() const ;

  double generateGamma(const double& a);

private:

  ClassDef(RooBetaPdf,1) // BetaPdf PDF
};

#endif
