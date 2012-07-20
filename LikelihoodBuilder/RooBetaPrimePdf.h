/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id: RooBetaPrimePdf.h,v 1.16 2007/07/12 20:30:49 wouter Exp $
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
#ifndef ROO_BETAPRIMEPDF
#define ROO_BETAPRIMEPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooBetaPrimePdf : public RooAbsPdf {
public:
  RooBetaPrimePdf() {} ;
  RooBetaPrimePdf(const char *name, const char *title,
	      RooAbsReal& _x, RooAbsReal& _alpha, RooAbsReal& _beta);
  RooBetaPrimePdf(const RooBetaPrimePdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooBetaPrimePdf(*this,newname); }
  inline virtual ~RooBetaPrimePdf() { }

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

  ClassDef(RooBetaPrimePdf,1) // BetaPdf PDF
};

#endif
