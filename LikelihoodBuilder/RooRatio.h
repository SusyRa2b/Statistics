/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROORATIO
#define ROORATIO

#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooRatio : public RooAbsReal {
public:
  RooRatio() {} ; 
  RooRatio(const char *name, const char *title,
	      RooAbsReal& _numerator,
	      RooAbsReal& _denominator);
  RooRatio(const RooRatio& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooRatio(*this,newname); }
  inline virtual ~RooRatio() { }

protected:

  RooRealProxy numerator ;
  RooRealProxy denominator ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooRatio,1) // Your description goes here...
};
 
#endif
