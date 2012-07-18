#include "betaHelperFunctions.h"

#include <iostream>
#include <string.h>
#include <complex>
#include <map>

#include "TString.h"
#include "TPRegexp.h"

#include "RooArgSet.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooAddition.h"
#include "RooWorkspace.h"
#include "RooUniform.h"
#include "RooPoisson.h"

#include "RooRatio.h"
#include "RooNormalFromFlatPdf.h"
#include "RooBetaInverseCDF.h"
#include "RooBetaPrimeInverseCDF.h"
#include "RooCorrelatedBetaGeneratorHelper.h"
#include "RooCorrelatedBetaPrimeGeneratorHelper.h"

using namespace RooFit ;

RooAbsArg* getCorrelatedBetaPrimeConstraint(RooWorkspace& ws,const TString varName,const TString binName,
					    const double value ,const double error,
					    const TString observables, const TString nuisances,
					    const TString correlatedName,
                                            Bool_t useComplement=false )
{
  //handle case with 0 error
  if(error == 0.)
    {
      RooRealVar valueVar ( varName+binName+"_value" , varName+binName+"_value" ,value );
      valueVar.setConstant();
      ws.import(valueVar);
      //NB this should not be on the 'observables' or 'nuisances' list; it should just be a constant
      return ws.arg(varName+binName+"_value");
    }
  //
  double alpha,beta;
  betaPrimeModeTransform(value , error , alpha , beta );

  RooRealVar alphaVar ( varName+binName+"_alpha" , varName+binName+"_alpha" ,alpha );
  RooRealVar betaVar ( varName+binName+"_beta" , varName+binName+"_beta" , beta );

  alphaVar.setConstant();
  betaVar.setConstant();

  ws.import(alphaVar);
  ws.extendSet(observables,alphaVar.GetName());

  ws.import(betaVar);
  ws.extendSet(observables,betaVar.GetName());

  RooAbsReal* correlationParameter = ws.function(correlatedName);
  if(correlationParameter == NULL)
    {
      RooRealVar parameter(correlatedName,correlatedName,0.5,0.,1.);
      ws.import(parameter);
      ws.extendSet(nuisances,parameter.GetName());
      RooNormalFromFlatPdf constraint(correlatedName+"_Constraint",correlatedName+"_Constraint",parameter);
      ws.import(constraint,RecycleConflictNodes());
      correlationParameter = ws.function(correlatedName);
    }
  assert(correlationParameter != NULL);

  RooCorrelatedBetaPrimeGeneratorHelper helper(varName+binName+"_GeneratorHelper",varName+binName+"_GeneratorHelper",*correlationParameter,alphaVar,betaVar);
  ws.import(helper,RecycleConflictNodes());

  RooBetaPrimeInverseCDF inverseCDF(varName+binName+"_BetaPrimeInverseCDF",varName+binName+"_BetaPrimeInverseCDF",*correlationParameter,alphaVar,betaVar, useComplement);
  ws.import(inverseCDF,RecycleConflictNodes());

  return ws.arg(varName+binName+"_BetaPrimeInverseCDF");
}


RooAbsArg* getCorrelatedBetaConstraint(RooWorkspace& ws,const TString varName,const TString binName,
				       const double value ,const double error,
				       const TString observables, const TString nuisances,
				       const TString correlatedName,
                                       Bool_t useComplement=false )
{
  //handle case with 0 error
  if(error == 0.)
    {
      RooRealVar valueVar ( varName+binName+"_value" , varName+binName+"_value" ,value );
      valueVar.setConstant();
      ws.import(valueVar);
      //NB this should not be on the 'observables' or 'nuisances' list; it should just be a constant
      return ws.arg(varName+binName+"_value");
    }
  //
  double alpha,beta;
  betaModeTransform(value , error , alpha , beta );

  RooRealVar alphaVar ( varName+binName+"_alpha" , varName+binName+"_alpha" ,alpha );
  RooRealVar betaVar ( varName+binName+"_beta" , varName+binName+"_beta" , beta );

  alphaVar.setConstant();
  betaVar.setConstant();

  ws.import(alphaVar);
  ws.extendSet(observables,alphaVar.GetName());

  ws.import(betaVar);
  ws.extendSet(observables,betaVar.GetName());

  RooAbsReal* correlationParameter = ws.function(correlatedName);
  if(correlationParameter == NULL)
    {
      RooRealVar parameter(correlatedName,correlatedName,0.5,0.,1.);
      ws.import(parameter);
      ws.extendSet(nuisances,parameter.GetName());
      RooNormalFromFlatPdf constraint(correlatedName+"_Constraint",correlatedName+"_Constraint",parameter);
      ws.import(constraint,RecycleConflictNodes());
      correlationParameter = ws.function(correlatedName);
    }
  assert(correlationParameter != NULL);

  RooCorrelatedBetaGeneratorHelper helper(varName+binName+"_GeneratorHelper",varName+binName+"_GeneratorHelper",*correlationParameter,alphaVar,betaVar);
  ws.import(helper,RecycleConflictNodes());

  RooBetaInverseCDF inverseCDF(varName+binName+"_BetaInverseCDF",varName+binName+"_BetaInverseCDF",*correlationParameter,alphaVar,betaVar, useComplement);
  ws.import(inverseCDF,RecycleConflictNodes());

  return ws.arg(varName+binName+"_BetaInverseCDF");
}

RooAbsArg* getBetaPrimeConstraint(RooWorkspace& ws,const TString varName,const TString binName,
				  const double value ,const double error,
				  const TString observables, const TString nuisances,
				  TString* passObs = NULL , TString* failObs = NULL,
				  TString* passPar = NULL , TString* failPar = NULL)
{
  //handle case with 0 error
  if(error == 0.)
    {
      RooRealVar valueVar ( varName+binName+"_value" , varName+binName+"_value" ,value );
      valueVar.setConstant();
      ws.import(valueVar);
      //NB this should not be on the 'observables' or 'nuisances' list; it should just be a constant
      return ws.arg(varName+binName+"_value");
    }
  //
  double alpha,beta;
  betaPrimeModeTransform(value , error , alpha , beta );

  RooRealVar passObservable (varName+binName+"_PassObs", varName+binName+"_PassObs", alpha-1,1e-5,1e5);
  passObservable.setConstant();
  RooRealVar failObservable (varName+binName+"_FailObs", varName+binName+"_FailObs", beta -1,1e-5,1e5);
  failObservable.setConstant();

  ws.import(passObservable);
  ws.import(failObservable);
  ws.extendSet(observables,passObservable.GetName());
  ws.extendSet(observables,failObservable.GetName());

  RooRealVar passParameter (varName+binName+"_PassPar", varName+binName+"_PassPar", alpha-1,1e-5,1e5);
  RooRealVar failParameter (varName+binName+"_FailPar", varName+binName+"_FailPar", beta +1,1e-5,1e5);

  ws.import(passParameter);
  ws.import(failParameter);
  ws.extendSet(nuisances,passParameter.GetName());
  ws.extendSet(nuisances,failParameter.GetName());

  if(passObs) *passObs = passObservable.GetName();
  if(failObs) *failObs = failObservable.GetName();
  if(passPar) *passPar = passParameter.GetName();
  if(failPar) *failPar = failParameter.GetName();

  RooPoisson passConstraint (varName+binName+"_PassConstraint" , varName+binName+"_PassConstraint", passObservable, passParameter);
  RooPoisson failConstraint (varName+binName+"_FailConstraint" , varName+binName+"_FailConstraint", failObservable, failParameter);

  ws.import(passConstraint, RecycleConflictNodes());
  ws.import(failConstraint, RecycleConflictNodes());

  RooRatio ratio(varName+binName+"_Ratio",varName+binName+"_Ratio",passParameter,failParameter);
  ws.import(ratio, RecycleConflictNodes());
  return ws.arg(varName+binName+"_Ratio");

}

RooAbsArg* getBetaConstraint(RooWorkspace& ws,const TString varName,const TString binName,
			     const double value ,const double error,
			     const TString observables, const TString nuisances,
			     TString* passObs = NULL , TString* failObs = NULL,
			     TString* passPar = NULL , TString* failPar = NULL)
{
  //handle case with 0 error
  if(error == 0.)
    {
      RooRealVar valueVar ( varName+binName+"_value" , varName+binName+"_value" ,value );
      valueVar.setConstant();
      ws.import(valueVar);
      //NB this should not be on the 'observables' or 'nuisances' list; it should just be a constant
      return ws.arg(varName+binName+"_value");
    }
  //
  double alpha,beta;
  betaModeTransform(value , error , alpha , beta ); 

  RooRealVar passObservable (varName+binName+"_PassObs", varName+binName+"_PassObs", alpha-1,1e-5,1e5);
  passObservable.setConstant();
  RooRealVar failObservable (varName+binName+"_FailObs", varName+binName+"_FailObs", beta -1,1e-5,1e5);
  failObservable.setConstant();

  ws.import(passObservable);
  ws.import(failObservable);
  ws.extendSet(observables,passObservable.GetName());
  ws.extendSet(observables,failObservable.GetName());

  RooRealVar passParameter (varName+binName+"_PassPar", varName+binName+"_PassPar", alpha-1,1e-5,1e5);
  RooRealVar failParameter (varName+binName+"_FailPar", varName+binName+"_FailPar", beta -1,1e-5,1e5);

  ws.import(passParameter);
  ws.import(failParameter);
  ws.extendSet(nuisances,passParameter.GetName());
  ws.extendSet(nuisances,failParameter.GetName());

  if(passObs) *passObs = passObservable.GetName();
  if(failObs) *failObs = failObservable.GetName();
  if(passPar) *passPar = passParameter.GetName();
  if(failPar) *failPar = failParameter.GetName();

  RooPoisson passConstraint (varName+binName+"_PassConstraint" , varName+binName+"_PassConstraint", passObservable, passParameter);
  RooPoisson failConstraint (varName+binName+"_FailConstraint" , varName+binName+"_FailConstraint", failObservable, failParameter);

  ws.import(passConstraint, RecycleConflictNodes());
  ws.import(failConstraint, RecycleConflictNodes());

  RooAddition passPlusFail(varName+binName+"_PassPlusFail",varName+binName+"_PassPlusFail",RooArgSet(passParameter,failParameter));

  RooRatio ratio(varName+binName+"_Ratio",varName+binName+"_Ratio",passParameter,passPlusFail);
  ws.import(ratio, RecycleConflictNodes());
  return ws.arg(varName+binName+"_Ratio");

}

RooAbsArg* getInverseBetaConstraint(RooWorkspace& ws,const TString varName,const TString binName,
				    const double value ,const double error,
				    const TString observables, const TString nuisances,
				    TString* passObs = NULL , TString* failObs = NULL,
				    TString* passPar = NULL , TString* failPar = NULL)
{
  //handle case with 0 error
  if(error == 0.)
    {
      RooRealVar valueVar ( varName+binName+"_invValue" , varName+binName+"_invValue" , 1.0/value );
      valueVar.setConstant();
      ws.import(valueVar);
      //NB this should not be on the 'observables' or 'nuisances' list; it should just be a constant
      return ws.arg(varName+binName+"_value");
    }
  //

  double alpha,beta;
  betaModeTransform(value , error , alpha , beta );

  RooRealVar passObservable (varName+binName+"_PassObs", varName+binName+"_PassObs", alpha-1,1e-5,1e5);
  passObservable.setConstant();
  RooRealVar failObservable (varName+binName+"_FailObs", varName+binName+"_FailObs", beta -1,1e-5,1e5);
  failObservable.setConstant();

  ws.import(passObservable);
  ws.import(failObservable);
  ws.extendSet(observables,passObservable.GetName());
  ws.extendSet(observables,failObservable.GetName());

  RooRealVar passParameter (varName+binName+"_PassPar", varName+binName+"_PassPar", alpha-1,1e-5,1e5);
  RooRealVar failParameter (varName+binName+"_FailPar", varName+binName+"_FailPar", beta -1,1e-5,1e5);

  ws.import(passParameter);
  ws.import(failParameter);
  ws.extendSet(nuisances,passParameter.GetName());
  ws.extendSet(nuisances,failParameter.GetName());

  if(passObs) *passObs = passObservable.GetName();
  if(failObs) *failObs = failObservable.GetName();
  if(passPar) *passPar = passParameter.GetName();
  if(failPar) *failPar = failParameter.GetName();

  RooPoisson passConstraint (varName+binName+"_PassConstraint" , varName+binName+"_PassConstraint", passObservable, passParameter);
  RooPoisson failConstraint (varName+binName+"_FailConstraint" , varName+binName+"_FailConstraint", failObservable, failParameter);

  ws.import(passConstraint, RecycleConflictNodes());
  ws.import(failConstraint, RecycleConflictNodes());

  RooAddition passPlusFail(varName+binName+"_PassPlusFail",varName+binName+"_PassPlusFail",RooArgSet(passParameter,failParameter));

  RooRatio ratio(varName+binName+"_InverseRatio",varName+binName+"_InverseRatio",passPlusFail,passParameter);
  ws.import(ratio, RecycleConflictNodes());
  return ws.arg(varName+binName+"_InverseRatio");

}


