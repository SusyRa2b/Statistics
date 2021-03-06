/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooBetaPdf.h,v 1.3 2012/11/02 13:31:53 kreis Exp $
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

  virtual Double_t getValV(const RooArgSet* set=0) const ;
  virtual Double_t getLogVal(const RooArgSet* set=0) const ;

  void setSigma0(Double_t newValue) { _sigma0 = newValue; _x0 = -1; }

protected:

  RooRealProxy x ;
  RooRealProxy alpha ;
  RooRealProxy beta ;
  
  Double_t evaluate() const ;

  double generateGamma(const double& a);

  mutable Double_t _logValue ;
  
  Double_t _sigma0 ;
  Double_t _alpha0 ;
  Double_t _beta0 ;
  Double_t _x0 ;

private:

  ClassDef(RooBetaPdf,1) // BetaPdf PDF
};

#endif
