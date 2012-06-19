#include <iostream>
#include <string.h>
#include <complex>
#include <map>

#include "TCanvas.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TText.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TH2F.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TStringLong.h"
#include "TString.h"
#include "TPRegexp.h"

#include "RooArgSet.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooTrace.h"
#include "RooUniform.h"
#include "RooAddition.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooRatio.h"
#include "RooAddition.h"
#include "RooPoisson.h"

#include "RooStats/ModelConfig.h"

#include "TMath.h"

using namespace RooFit ;
using namespace RooStats ;

void betaPrimeModeTransform(const double& mu_measured , const double& sigma_measured ,
			    double& alpha , double& beta ) 
 {
  complex<double> mu ( mu_measured, 0.0 ) ;
  complex<double> sigma2 ( sigma_measured * sigma_measured , 0.0 ) ;
  //complex<double> complex_beta = sigma2*2.;
  //complex<double> cubeRoot = pow( sqrt(pow(sigma2,3)+pow(sigma2,2)*105.+ sigma2*3.-1.) / (sigma2*sqrt(3)) + (pow(sigma2,3)+pow(sigma2,2)*168.+sigma2*93.+8.)/(pow(sigma2,3)*27.) , 1.0/3.0);
  //complex<double> complex_beta = cubeRoot + (pow(sigma2,2) + sigma2*31. + 4.)/(pow(sigma2,2)*9.*cubeRoot) + (sigma2*4.+2.)/(sigma2*3.);
  complex<double> complex_beta = (-0.3333333333333333*(-1.*mu - 1.*pow(mu,2) - 4.*sigma2))/sigma2 - (0.41997368329829105*
      (-1.*pow(-1.*mu - 1.*pow(mu,2) - 4.*sigma2,2) + 3.*sigma2*(-1. - 2.*mu - 2.*pow(mu,2) + 5.*sigma2)))/
    (sigma2*pow(2.*pow(mu,3) + 6.*pow(mu,4) + 6.*pow(mu,5) + 2.*pow(mu,6) + 9.*mu*sigma2 + 51.*pow(mu,2)*sigma2 + 84.*pow(mu,3)*sigma2 + 42.*pow(mu,4)*sigma2 + 
        36.*pow(sigma2,2) + 150.*mu*pow(sigma2,2) + 150.*pow(mu,2)*pow(sigma2,2) + 2.*pow(sigma2,3) + 
        sqrt(pow(2.*pow(mu,3) + 6.*pow(mu,4) + 6.*pow(mu,5) + 2.*pow(mu,6) + 9.*mu*sigma2 + 51.*pow(mu,2)*sigma2 + 84.*pow(mu,3)*sigma2 + 42.*pow(mu,4)*sigma2 + 
            36.*pow(sigma2,2) + 150.*mu*pow(sigma2,2) + 150.*pow(mu,2)*pow(sigma2,2) + 2.*pow(sigma2,3),2) + 
          4.*pow(-1.*pow(-1.*mu - 1.*pow(mu,2) - 4.*sigma2,2) + 3.*sigma2*(-1. - 2.*mu - 2.*pow(mu,2) + 5.*sigma2),3)),0.3333333333333333)) + 
   (0.26456684199469993*pow(2.*pow(mu,3) + 6.*pow(mu,4) + 6.*pow(mu,5) + 2.*pow(mu,6) + 9.*mu*sigma2 + 51.*pow(mu,2)*sigma2 + 84.*pow(mu,3)*sigma2 + 
        42.*pow(mu,4)*sigma2 + 36.*pow(sigma2,2) + 150.*mu*pow(sigma2,2) + 150.*pow(mu,2)*pow(sigma2,2) + 2.*pow(sigma2,3) + 
        sqrt(pow(2.*pow(mu,3) + 6.*pow(mu,4) + 6.*pow(mu,5) + 2.*pow(mu,6) + 9.*mu*sigma2 + 51.*pow(mu,2)*sigma2 + 84.*pow(mu,3)*sigma2 + 42.*pow(mu,4)*sigma2 + 
            36.*pow(sigma2,2) + 150.*mu*pow(sigma2,2) + 150.*pow(mu,2)*pow(sigma2,2) + 2.*pow(sigma2,3),2) + 
	     4.*pow(-1.*pow(-1.*mu - 1.*pow(mu,2) - 4.*sigma2,2) + 3.*sigma2*(-1. - 2.*mu - 2.*pow(mu,2) + 5.*sigma2),3)),0.3333333333333333))/sigma2;
  beta = complex_beta.real();
  alpha = 1. + mu_measured + mu_measured*beta;
  cout << "Beta Prime Mode Transform " << endl;
  cout << "The calculated value of alpha is : " << alpha << endl;
  cout << "The calculated value of beta is  : " << beta << endl;
  cout << "The input value of mu is         : " << mu_measured << endl;
  cout << "The input value of sigma is      : " << sigma_measured << endl;
  cout << "The calculated value of mu is    : " << (alpha - 1) / ( beta + 1 ) << endl;
  cout << "The calculated value of sigma is : " << sqrt( alpha * (alpha + beta - 1 ) / ( pow(beta - 1 , 2 ) * (beta - 2) ) ) << endl;
}

void betaModeTransform(const double& mu_measured , const double& sigma_measured ,
		   double& alpha , double& beta ) 
{
  //double sigma2 = sigma_measured*sigma_measured;
  //double mu2 = mu_measured*mu_measured;
  //double mu3 = mu_measured*mu2;
  //alpha = (mu2 - mu3 - mu_measured * sigma2) / sigma2 ;
  //beta = (alpha - alpha * mu_measured) / mu_measured ;
  complex<double> mu ( mu_measured, 0.0 ) ;
  complex<double> mu2 = mu * mu ;
  complex<double> mu3 = mu2 * mu ;
  complex<double> muMinusOne ( mu_measured - 1 , 0.0 ) ;
  complex<double> muMinusOne2 = muMinusOne * muMinusOne ;
  complex<double> muMinusOne3 = muMinusOne2 * muMinusOne ;
  complex<double> muMinusOne6 = muMinusOne3 * muMinusOne3 ;
  complex<double> sigma2 ( sigma_measured * sigma_measured , 0.0 ) ;
  complex<double> sigma4 = sigma2 * sigma2 ;
  complex<double> sigma6 = sigma4 * sigma2 ;
  complex<double> common_expression = pow( muMinusOne3 * ( muMinusOne3 * mu3 * -2. - 
							   muMinusOne * mu * 3. * ( muMinusOne * mu * 14. + 3. ) * sigma2 - 
							   ( mu * 5. - 3. ) * ( mu * 5. - 2. ) * sigma4 * 6. - sigma6 * 2. ) + 
					   sqrt( muMinusOne6 * pow( sigma2 - 2. * mu * sigma2 , 2 ) *
						 ( -1. * muMinusOne2 * mu2 + 
						   4. * ( -1. + ( muMinusOne ) * mu * ( -4. + 3. * ( muMinusOne ) * mu ) ) * sigma2 + 
						   4. * ( 11. + 47. * ( muMinusOne ) * mu ) * sigma4 + 4. * sigma6 ) ) * 5.196152422706632 ,
					   0.3333333333333333 );

  complex<double> complex_beta = ( ( muMinusOne2 * mu * 4. + ( mu * 7. - 4. ) * sigma2 * 4. + 
				     ( complex<double> ( 2.519842099789747 , 4.364494543886885 ) * muMinusOne2 *
				       ( muMinusOne2 * mu2 + ( muMinusOne * mu * 14. + 3. ) * sigma2 + sigma4 ) ) /
				     common_expression + complex<double> ( 1.5874010519681996 , -2.7494592739972052 ) * 
				     common_expression ) * 0.08333333333333333 ) / sigma2;
  beta = complex_beta.real();
  alpha = (2 * mu_measured - beta * mu_measured - 1) / ( mu_measured - 1 ) ;

  if( !(alpha >= 1) && !(beta >= 1) )
    {
      cout << "mode and variance impossible with beta pdf, setting uniform distribution" << endl;
      alpha = 1;
      beta = 1;
    }
  cout << "Beta Mode Transform " << endl;
  cout << "The calculated value of alpha is : " << alpha << endl;
  cout << "The calculated value of beta is  : " << beta << endl;
  cout << "The input value of mu is         : " << mu_measured << endl;
  cout << "The input value of sigma is      : " << sigma_measured << endl;
  cout << "The calculated value of mu is    : " << (alpha - 1) / ( alpha + beta - 2 ) << endl;
  cout << "The calculated value of sigma is : " << sqrt( alpha * beta / ( pow(alpha + beta,2) * (alpha + beta + 1) ) ) << endl;
}

struct channels
{
  double zeroLepton;
  double zeroLeptonLowDeltaPhiN;
  double oneMuon;
  double oneElectron;
} ;

struct yields
{
  double zeroLepton;
  double zeroLeptonLowDeltaPhiN;
  double oneLepton;
} ;

struct abcdBinParameters
{
  double zeroLeptonTriggerEfficiency;
  double zeroLeptonTriggerEfficiencyError;
  double zeroLeptonLowDeltaPhiNTriggerEfficiency;
  double zeroLeptonLowDeltaPhiNTriggerEfficiencyError;
  double oneElectronTriggerEfficiency;
  double oneElectronTriggerEfficiencyError;
  double oneMuonTriggerEfficiency;
  double oneMuonTriggerEfficiencyError;
  double zeroLeptonLowDeltaPhiNMC;
  double ZtoNuNubTagScaling;
  double ZtoNuNubTagScalingError;
  double ZtollAcceptance;
  double ZtollAcceptanceError;
  double qcdClosure;
  double qcdClosureError;
  double topWJetsClosure;
  double topWJetsClosureError;
  double ZtoNuNuClosure;
  double ZtoNuNuClosureError;
} ;

struct allBinNames
{
  TString ZtoNuNu;
  TString deltaPhiNRatio;
  TString singleLeptonScaling;
  TString MCUncertainty;
  TString signalCrossSection;
  TString signalCrossSectionPassObs;
  TString signalCrossSectionFailObs;
  TString signalCrossSectionPassPar;
  TString signalCrossSectionFailPar;
  TString observables;
  TString nuisances;
  TString qcdClosurePassObs;
  TString qcdClosureFailObs;
  TString qcdClosurePassPar;
  TString qcdClosureFailPar;
  TString topWJetsClosurePassObs;
  TString topWJetsClosureFailObs;
  TString topWJetsClosurePassPar;
  TString topWJetsClosureFailPar;
  TString ZtoNuNuClosurePassObs;
  TString ZtoNuNuClosureFailObs;
  TString ZtoNuNuClosurePassPar;
  TString ZtoNuNuClosureFailPar;
} ;

struct allBins
{
  double diElectronCount;
  double diMuonCount;
  double ZtollOverZtoNuNuRatio;
  double ZtoeePurity;
  double ZtoeePurityError;
  double ZtomumuPurity;
  double ZtomumuPurityError;
  double ZtoeeEfficiency;
  double ZtoeeEfficiencyError;
  double ZtomumuEfficiency;
  double ZtomumuEfficiencyError;
  double deltaPhiNRatio;
  double deltaPhiNRatioError;
  double qcdClosure;
  double qcdClosureError;
  double topWJetsClosure;
  double topWJetsClosureError;
  double ZtoNuNuClosure;
  double ZtoNuNuClosureError;
  double MCUncertainty;
  double signal;
  double signalError;
} ;

RooAbsArg* getCorrelatedBetaPrimeConstraint(RooWorkspace& ws,const TString varName,const TString binName,
					    const double value ,const double error,
					    const TString correlatedPassObs,const TString correlatedFailObs,
					    const TString correlatedPassPar,const TString correlatedFailPar)
{
  double alpha,beta;
  betaPrimeModeTransform(value , error , alpha , beta );

  RooAbsReal* passObs = ws.function(correlatedPassObs);
  RooAbsReal* failObs = ws.function(correlatedFailObs);
  
  RooRealVar passScale ( varName+binName+"_PassScale" , varName+binName+"_PassScale" , (alpha-1) / passObs->getVal() );
  RooRealVar failScale ( varName+binName+"_FailScale" , varName+binName+"_FailScale" , (beta -1) / failObs->getVal() );

  passScale.setConstant();
  failScale.setConstant();

  RooAbsReal* passPar = ws.function(correlatedPassPar);
  RooAbsReal* failPar = ws.function(correlatedFailPar);

  RooProduct pass(varName+binName+"_Pass" , varName+binName+"_Pass" , RooArgSet(*passPar,passScale));
  RooProduct fail(varName+binName+"_Fail" , varName+binName+"_Fail" , RooArgSet(*failPar,failScale));

  RooRatio ratio(varName+binName+"_Ratio",varName+binName+"_Ratio",pass,fail);
  ws.import(ratio, RecycleConflictNodes());
  return ws.arg(varName+binName+"_Ratio");
}

RooAbsArg* getBetaPrimeConstraint(RooWorkspace& ws,const TString varName,const TString binName,
				  const double value ,const double error,
				  const TString observables, const TString nuisances,
				  TString* passObs = NULL , TString* failObs = NULL,
				  TString* passPar = NULL , TString* failPar = NULL)
{

  double alpha,beta;
  betaPrimeModeTransform(value , error , alpha , beta );

  RooRealVar passObservable (varName+binName+"_PassObs", varName+binName+"_PassObs", alpha-1,1e-5,1e5);
  passObservable.setConstant();
  RooRealVar failObservable (varName+binName+"_FailObs", varName+binName+"_FailObs", alpha-1,1e-5,1e5);
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

RooAbsArg* getCorrelatedBetaConstraint(RooWorkspace& ws,const TString varName,const TString binName,
				       const double value ,const double error,
				       const TString correlatedPassObs,const TString correlatedFailObs,
				       const TString correlatedPassPar,const TString correlatedFailPar)
{
  double alpha,beta;
  betaPrimeModeTransform(value , error , alpha , beta );

  RooAbsReal* passObs = ws.function(correlatedPassObs);
  RooAbsReal* failObs = ws.function(correlatedFailObs);
  
  RooRealVar passScale ( varName+binName+"_PassScale" , varName+binName+"_PassScale" , (alpha-1) / passObs->getVal() );
  RooRealVar failScale ( varName+binName+"_FailScale" , varName+binName+"_FailScale" , (beta -1) / failObs->getVal() );

  passScale.setConstant();
  failScale.setConstant();

  RooAbsReal* passPar = ws.function(correlatedPassPar);
  RooAbsReal* failPar = ws.function(correlatedFailPar);

  RooProduct pass(varName+binName+"_Pass" , varName+binName+"_Pass" , RooArgSet(*passPar,passScale));
  RooProduct fail(varName+binName+"_Fail" , varName+binName+"_Fail" , RooArgSet(*failPar,failScale));

  RooAddition passPlusFail(varName+binName+"_PassPlusFail",varName+binName+"_PassPlusFail",RooArgSet(pass,fail));

  RooRatio ratio(varName+binName+"_Ratio",varName+binName+"_Ratio",pass,passPlusFail);
  ws.import(ratio, RecycleConflictNodes());
  return ws.arg(varName+binName+"_Ratio");
}

RooAbsArg* getBetaConstraint(RooWorkspace& ws,const TString varName,const TString binName,
			     const double value ,const double error,
			     const TString observables, const TString nuisances,
			     TString* passObs = NULL , TString* failObs = NULL,
			     TString* passPar = NULL , TString* failPar = NULL)
{
  double alpha,beta;
  betaPrimeModeTransform(value , error , alpha , beta );

  RooRealVar passObservable (varName+binName+"_PassObs", varName+binName+"_PassObs", alpha-1,1e-5,1e5);
  passObservable.setConstant();
  RooRealVar failObservable (varName+binName+"_FailObs", varName+binName+"_FailObs", alpha-1,1e-5,1e5);
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

  double alpha,beta;
  betaPrimeModeTransform(value , error , alpha , beta );

  RooRealVar passObservable (varName+binName+"_PassObs", varName+binName+"_PassObs", alpha-1,1e-5,1e5);
  passObservable.setConstant();
  RooRealVar failObservable (varName+binName+"_FailObs", varName+binName+"_FailObs", alpha-1,1e-5,1e5);
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

bool makeOneBin(RooWorkspace& ws , TString& binName , allBinNames& names , channels& observed , abcdBinParameters& abcd , yields& signal , yields& signalError , double& oneLeptonTotal, double& zeroLeptonTopWJetsGuess)
{

  //Make sure that names are unique to this bin

  TString zeroLeptonName("zeroLepton_");
  zeroLeptonName+=binName;
  TString zeroLeptonLowDeltaPhiNName("zeroLeptonLowDeltaPhiN_");
  zeroLeptonLowDeltaPhiNName+=binName;
  TString oneLeptonName("oneLepton_");
  oneLeptonName+=binName;
  TString oneMuonName("oneMuon_");
  oneMuonName+=binName;
  TString oneElectronName("oneElectron_");
  oneElectronName+=binName;

  //Define counts and add them to the workspace

  RooRealVar zeroLeptonCount(zeroLeptonName+"_Count",zeroLeptonName+"_Count",observed.zeroLepton);
  RooRealVar zeroLeptonLowDeltaPhiNCount(zeroLeptonLowDeltaPhiNName+"_Count",zeroLeptonLowDeltaPhiNName+"_Count",observed.zeroLeptonLowDeltaPhiN);
  RooRealVar oneMuonCount(oneMuonName+"_Count",oneMuonName+"_Count",observed.oneMuon);
  RooRealVar oneElectronCount(oneElectronName+"_Count",oneElectronName+"_Count",observed.oneElectron);

  zeroLeptonCount.setConstant();
  zeroLeptonLowDeltaPhiNCount.setConstant();
  oneMuonCount.setConstant();
  oneElectronCount.setConstant();

  ws.import(zeroLeptonCount);
  ws.import(zeroLeptonLowDeltaPhiNCount);
  ws.import(oneMuonCount);
  ws.import(oneElectronCount);
  ws.extendSet(names.observables,zeroLeptonCount.GetName());
  ws.extendSet(names.observables,zeroLeptonLowDeltaPhiNCount.GetName());
  ws.extendSet(names.observables,oneMuonCount.GetName());
  ws.extendSet(names.observables,oneElectronCount.GetName());

  //Define QCD component

  RooAbsArg* deltaPhiNRatio = ws.arg(names.deltaPhiNRatio);
  double qcdGuess = observed.zeroLeptonLowDeltaPhiN - abcd.zeroLeptonLowDeltaPhiNMC;
  if(qcdGuess < 0) qcdGuess = 0;
  RooRealVar zeroLeptonLowDeltaPhiNQCDYield(zeroLeptonLowDeltaPhiNName+"_QCDYield",zeroLeptonLowDeltaPhiNName+"_QCDYield",qcdGuess,1e-5,10000);
  ws.import(zeroLeptonLowDeltaPhiNQCDYield);
  ws.extendSet(names.nuisances,zeroLeptonLowDeltaPhiNQCDYield.GetName());
  RooAbsArg* zeroLeptonQCDClosure = 
    getCorrelatedBetaPrimeConstraint(ws,"zeroLeptonQCDClosure_",binName,
				     abcd.qcdClosure,abcd.qcdClosureError,
				     names.qcdClosurePassObs,names.qcdClosureFailObs,
				     names.qcdClosurePassPar,names.qcdClosureFailPar);
  RooProduct zeroLeptonQCDYield(zeroLeptonName+"_QCDYield",zeroLeptonName+"_QCDYield",RooArgSet(*deltaPhiNRatio,*zeroLeptonQCDClosure,zeroLeptonLowDeltaPhiNQCDYield));

  //Define top and W+jets component:

  RooRealVar* singleLeptonScaling = ws.var(names.singleLeptonScaling);
  RooRealVar  oneLeptonTopWJetsYield(oneLeptonName+"_TopWJetsYield",oneLeptonName+"_TopWJetsYield",0.5*(observed.oneMuon+observed.oneElectron),1e-5,100000);
  ws.import(oneLeptonTopWJetsYield);
  ws.extendSet(names.nuisances,oneLeptonTopWJetsYield.GetName());
  oneLeptonTotal += oneLeptonTopWJetsYield.getVal();
  RooAbsArg*  zeroLeptonTopWJetsClosure = 
    getCorrelatedBetaPrimeConstraint(ws,"zeroLeptonTopWJetsClosure_",binName,
				     abcd.topWJetsClosure,abcd.topWJetsClosureError,
				     names.topWJetsClosurePassObs,names.topWJetsClosureFailObs,
				     names.topWJetsClosurePassPar,names.topWJetsClosureFailPar);
  RooProduct  zeroLeptonTopWJetsYield(zeroLeptonName+"_TopWJetsYield",zeroLeptonName+"_TopWJetsYield",RooArgSet(*singleLeptonScaling,*zeroLeptonTopWJetsClosure,oneLeptonTopWJetsYield));

  //Define Z to invisible component:

  RooAbsArg* unscaledZtoNuNu = ws.arg(names.ZtoNuNu);
  RooAbsArg* zeroLeptonZtoNuNubTagScaling = 
    getBetaConstraint(ws,"zeroLeptonZtoNuNubTagScaling_",binName,
		      abcd.ZtoNuNubTagScaling,abcd.ZtoNuNubTagScalingError,
		      names.observables,names.nuisances);
  RooAbsArg* zeroLeptonZtollInvAcceptance = 
    getInverseBetaConstraint(ws,"zeroLeptonZtollInvAcceptance_",binName,
			     abcd.ZtollAcceptance,abcd.ZtollAcceptanceError,
			     names.observables,names.nuisances);
  RooAbsArg*  zeroLeptonZtoNuNuClosure = 
    getCorrelatedBetaPrimeConstraint(ws,"zeroLeptonZtoNuNuClosure_",binName,
				     abcd.ZtoNuNuClosure,abcd.ZtoNuNuClosureError,
				     names.ZtoNuNuClosurePassObs,names.ZtoNuNuClosureFailObs,
				     names.ZtoNuNuClosurePassPar,names.ZtoNuNuClosureFailPar);
  RooProduct zeroLeptonZtoNuNuYield(zeroLeptonName+"_ZtoNuNuYield",zeroLeptonName+"_ZtoNuNuYield",RooArgSet(*zeroLeptonZtoNuNuClosure,*zeroLeptonZtoNuNubTagScaling,*zeroLeptonZtollInvAcceptance,*unscaledZtoNuNu));

  // Monte Carlo yield (summed over all channels) in zero lepton low delta phi_N region

  RooAbsArg* MCUncertainty = ws.arg(names.MCUncertainty);
  RooRealVar zeroLeptonLowDeltaPhiNMCCount(zeroLeptonLowDeltaPhiNName+"_MCCount",zeroLeptonLowDeltaPhiNName+"_MCCount",abcd.zeroLeptonLowDeltaPhiNMC);
  zeroLeptonLowDeltaPhiNMCCount.setConstant();
  RooProduct  zeroLeptonLowDeltaPhiNMCYield(zeroLeptonLowDeltaPhiNName+"_MCYield",zeroLeptonLowDeltaPhiNName+"_MCYield",RooArgSet(*MCUncertainty,zeroLeptonLowDeltaPhiNMCCount));

  // Setup signal yields

  RooRealVar* signalCrossSection = ws.var(names.signalCrossSection);

  RooAbsArg*  zeroLeptonSignalYieldFraction = 
    getCorrelatedBetaPrimeConstraint(ws,"zeroLeptonSignalYieldFraction_",binName,
				     signal.zeroLepton,signalError.zeroLepton,
				     names.signalCrossSectionPassObs,names.signalCrossSectionFailObs,
				     names.signalCrossSectionPassPar,names.signalCrossSectionFailPar);

  RooAbsArg*  zeroLeptonLowDeltaPhiNSignalYieldFraction = 
    getCorrelatedBetaPrimeConstraint(ws,"zeroLeptonLowDeltaPhiNSignalYieldFraction_",binName,
				     signal.zeroLeptonLowDeltaPhiN,signalError.zeroLeptonLowDeltaPhiN,
				     names.signalCrossSectionPassObs,names.signalCrossSectionFailObs,
				     names.signalCrossSectionPassPar,names.signalCrossSectionFailPar);

  RooAbsArg*  oneLeptonSignalYieldFraction = 
    getCorrelatedBetaPrimeConstraint(ws,"oneLeptonSignalYieldFraction_",binName,
				     signal.oneLepton,signalError.oneLepton,
				     names.signalCrossSectionPassObs,names.signalCrossSectionFailObs,
				     names.signalCrossSectionPassPar,names.signalCrossSectionFailPar);
  
  RooProduct zeroLeptonSignalYield(zeroLeptonName+"_SignalYield",zeroLeptonName+"_SignalYield",RooArgSet(*signalCrossSection,*zeroLeptonSignalYieldFraction));
  RooProduct zeroLeptonLowDeltaPhiNSignalYield(zeroLeptonLowDeltaPhiNName+"_SignalYield",zeroLeptonLowDeltaPhiNName+"_SignalYield",RooArgSet(*signalCrossSection,*zeroLeptonLowDeltaPhiNSignalYieldFraction));
  RooProduct oneLeptonSignalYield(oneLeptonName+"_SignalYield",oneLeptonName+"_SignalYield",RooArgSet(*signalCrossSection,*oneLeptonSignalYieldFraction));
  
  // Setup yields in all bins

  RooAddition zeroLeptonYieldSum(zeroLeptonName+"_YieldSum",zeroLeptonName+"_YieldSum",RooArgSet(zeroLeptonSignalYield,zeroLeptonZtoNuNuYield,zeroLeptonTopWJetsYield,zeroLeptonQCDYield));
  double topGuess = observed.zeroLepton - zeroLeptonZtoNuNuYield.getVal() - zeroLeptonQCDYield.getVal();
  if(topGuess > 0 ) zeroLeptonTopWJetsGuess += topGuess;
  RooAddition zeroLeptonLowDeltaPhiNYieldSum(zeroLeptonLowDeltaPhiNName+"_YieldSum",zeroLeptonLowDeltaPhiNName+"_YieldSum",RooArgSet(zeroLeptonLowDeltaPhiNSignalYield,zeroLeptonLowDeltaPhiNQCDYield,zeroLeptonLowDeltaPhiNMCYield));
  RooAddition oneLeptonYieldSum(oneLeptonName+"_YieldSum",oneLeptonName+"_YieldSum",RooArgSet(oneLeptonSignalYield,oneLeptonTopWJetsYield));

  // Setup trigger efficiencies

  RooAbsArg* zeroLeptonTriggerEfficiency = 
    getBetaConstraint(ws,"zeroLeptonTriggerEfficiency_",binName,
		      abcd.zeroLeptonTriggerEfficiency,abcd.zeroLeptonTriggerEfficiencyError,
		      names.observables,names.nuisances);

  RooAbsArg* zeroLeptonLowDeltaPhiNTriggerEfficiency = 
    getBetaConstraint(ws,"zeroLeptonLowDeltaPhiNTriggerEfficiency_",binName,
		      abcd.zeroLeptonLowDeltaPhiNTriggerEfficiency,abcd.zeroLeptonLowDeltaPhiNTriggerEfficiencyError,
		      names.observables,names.nuisances);

  RooAbsArg* oneMuonTriggerEfficiency = 
    getBetaConstraint(ws,"oneMuonTriggerEfficiency_",binName,
		      abcd.oneMuonTriggerEfficiency,abcd.oneMuonTriggerEfficiencyError,
		      names.observables,names.nuisances);

  RooAbsArg* oneElectronTriggerEfficiency = 
    getBetaConstraint(ws,"oneElectronTriggerEfficiency_",binName,
		      abcd.oneElectronTriggerEfficiency,abcd.oneElectronTriggerEfficiencyError,
		      names.observables,names.nuisances);

  // Total Yields in bins

  RooProduct zeroLeptonYield(zeroLeptonName+"_Yield",zeroLeptonName+"_Yield",RooArgSet(*zeroLeptonTriggerEfficiency,zeroLeptonYieldSum));
  RooProduct zeroLeptonLowDeltaPhiNYield(zeroLeptonLowDeltaPhiNName+"_Yield",zeroLeptonLowDeltaPhiNName+"_Yield",RooArgSet(*zeroLeptonLowDeltaPhiNTriggerEfficiency,zeroLeptonLowDeltaPhiNYieldSum));
  RooProduct oneMuonYield(oneMuonName+"_Yield",oneMuonName+"_Yield",RooArgSet(*oneMuonTriggerEfficiency,oneLeptonYieldSum));
  RooProduct oneElectronYield(oneElectronName+"_Yield",oneElectronName+"_Yield",RooArgSet(*oneElectronTriggerEfficiency,oneLeptonYieldSum));

  // Define poisson constraints

  RooPoisson zeroLeptonConstraint(zeroLeptonName+"_Constraint",zeroLeptonName+"_Constraint",zeroLeptonYield,zeroLeptonCount);
  RooPoisson zeroLeptonLowDeltaPhiNConstraint(zeroLeptonLowDeltaPhiNName+"_Constraint",zeroLeptonLowDeltaPhiNName+"_Constraint",zeroLeptonLowDeltaPhiNYield,zeroLeptonLowDeltaPhiNCount);
  RooPoisson oneMuonConstraint(oneMuonName+"_Constraint",oneMuonName+"_Constraint",oneMuonYield,oneMuonCount);
  RooPoisson oneElectronConstraint(oneElectronName+"_Constraint",oneElectronName+"_Constraint",oneElectronYield,oneElectronCount);

  // Cleanup

  ws.import(zeroLeptonConstraint, RecycleConflictNodes());
  ws.import(zeroLeptonLowDeltaPhiNConstraint, RecycleConflictNodes());
  ws.import(oneMuonConstraint, RecycleConflictNodes());
  ws.import(oneElectronConstraint, RecycleConflictNodes());

  return true;

}

void setupUnderlyingLikelihood(RooWorkspace& ws ,allBinNames& names, allBins& numbers)
{
  ws.defineSet("observables","");
  names.observables = "observables";

  ws.defineSet("nuisances","");
  names.nuisances = "nuisances";

  //Universal parameters

  RooAbsArg* deltaPhiNRatio = 
    getBetaConstraint(ws,"deltaPhiNRatio","",
		      numbers.deltaPhiNRatio,numbers.deltaPhiNRatioError,
		      names.observables,names.nuisances);
  names.deltaPhiNRatio = deltaPhiNRatio->GetName();

  RooRealVar MCUncertainty("MCUncertainty","MCUncertainty",numbers.MCUncertainty);
  MCUncertainty.setConstant();
  ws.import(MCUncertainty);
  names.MCUncertainty = MCUncertainty.GetName();

  RooRealVar signalCrossSection("signalCrossSection","signalCrossSection",0.0,0.0,1e5);
  ws.import(signalCrossSection);
  names.signalCrossSection = signalCrossSection.GetName();

  RooAbsArg* signalUncertainty = 
    getBetaPrimeConstraint(ws,"signalUncertainty","",
			   numbers.signal,numbers.signalError,
			   names.observables, names.nuisances,
			   &names.signalCrossSectionPassObs,&names.signalCrossSectionFailObs,
			   &names.signalCrossSectionPassPar,&names.signalCrossSectionFailPar);

  RooRealVar singleLeptonScaling("singleLeptonScaling","singleLeptonScaling",0.0,1e5);
  ws.import(singleLeptonScaling);
  names.singleLeptonScaling = singleLeptonScaling.GetName();

  //Closure Tests:

  RooAbsArg* qcdClosure = 
    getBetaPrimeConstraint(ws,"qcdClosure","",
			   numbers.qcdClosure,numbers.qcdClosureError,
			   names.observables, names.nuisances,
			   &names.qcdClosurePassObs,&names.qcdClosureFailObs,
			   &names.qcdClosurePassPar,&names.qcdClosureFailPar);

  RooAbsArg* topWJetsClosure = 
    getBetaPrimeConstraint(ws,"topWJetsClosure","",
			   numbers.topWJetsClosure,numbers.topWJetsClosureError,
			   names.observables, names.nuisances,
			   &names.topWJetsClosurePassObs,&names.topWJetsClosureFailObs,
			   &names.topWJetsClosurePassPar,&names.topWJetsClosureFailPar);

  RooAbsArg* ZtoNuNuClosure = 
    getBetaPrimeConstraint(ws,"ZtoNuNuClosure","",
			   numbers.ZtoNuNuClosure,numbers.ZtoNuNuClosureError,
			   names.observables, names.nuisances,
			   &names.ZtoNuNuClosurePassObs,&names.ZtoNuNuClosureFailObs,
			   &names.ZtoNuNuClosurePassPar,&names.ZtoNuNuClosureFailPar);

  //dileptons for Z->invisible background:

  //Define Z->ll observables

  RooRealVar diMuonCount("diMuon_Count","diMuon_Count",numbers.diMuonCount);
  RooRealVar diElectronCount("diElectron_Count","diElectron_Count",numbers.diElectronCount);

  ws.import(diMuonCount);
  ws.import(diElectronCount);
  ws.extendSet(names.observables,diMuonCount.GetName());
  ws.extendSet(names.observables,diElectronCount.GetName());

  //Define Z->ll yields
  
  RooRealVar ZtollOverZtoNuNuRatio("ZtollOverZtoNuNuRatio","ZtollOverZtoNuNuRatio",numbers.ZtollOverZtoNuNuRatio);
  ZtollOverZtoNuNuRatio.setConstant();
  RooRealVar ZtoNuNu("ZtoNuNu","ZtoNuNu",numbers.ZtollOverZtoNuNuRatio*0.5*(numbers.diMuonCount*numbers.ZtomumuEfficiency/numbers.ZtomumuPurity + numbers.diElectronCount*numbers.ZtoeeEfficiency/numbers.ZtoeePurity),0.0,1e5);
  ws.import(ZtoNuNu);
  ws.extendSet(names.nuisances,ZtoNuNu.GetName());
  names.ZtoNuNu = ZtoNuNu.GetName();

  RooProduct Ztoll("Ztoll","Ztoll",RooArgSet(ZtoNuNu,ZtollOverZtoNuNuRatio));

  RooAbsArg* ZtomumuInvPurity = 
    getInverseBetaConstraint(ws,"ZtomumuInvPurity","",
			     numbers.ZtomumuPurity,numbers.ZtomumuPurityError,
			     names.observables,names.nuisances);

  RooAbsArg* ZtomumuEfficiency = 
    getBetaConstraint(ws,"Efficiency","",
		      numbers.ZtomumuEfficiency,numbers.ZtomumuEfficiencyError,
		      names.observables,names.nuisances);

  RooAbsArg* ZtoeeInvPurity = 
    getInverseBetaConstraint(ws,"ZtoeeInvPurity","",
			     numbers.ZtoeePurity,numbers.ZtoeePurityError,
			     names.observables,names.nuisances);

  RooAbsArg* ZtoeeEfficiency = 
    getBetaConstraint(ws,"Efficiency","",
		      numbers.ZtoeeEfficiency,numbers.ZtoeeEfficiencyError,
		      names.observables,names.nuisances);

  RooProduct diMuonYield("diMuon_Yield","diMuon_Yield",RooArgSet(Ztoll,*ZtomumuEfficiency,*ZtomumuInvPurity));
  RooProduct diElectronYield("diElectron_Yield","diElectron_Yield",RooArgSet(Ztoll,*ZtoeeEfficiency,*ZtoeeInvPurity));

  // Define dilepton poisson constraints

  RooPoisson diMuonConstraint("diMuon_Constraint","diMuon_Constraint",diMuonYield,diMuonCount);
  RooPoisson diElectronConstraint("diElectron_Constraint","diElectron_Constraint",diElectronYield,diElectronCount);

  ws.import(diMuonConstraint, RecycleConflictNodes());
  ws.import(diElectronConstraint, RecycleConflictNodes());

}

void setupSignalModel(vector<TString> binNames, TString signalModelFileName, int signalModelLineNumber, map<TString,yields>& signal,map<TString,yields>& signalError, const double& luminosity)
{
  cout << "getting the file: " << signalModelFileName << endl;

  ifstream signalFile;  
  
  signalFile.open(signalModelFileName.Data(),fstream::in);
  
  int thisLineNumber = 0;
  string fileLine;

  while( thisLineNumber < signalModelLineNumber && !signalFile.eof() )
    {
      thisLineNumber++;
      getline(signalFile,fileLine);
    }
  
  double m0,m12,susyGenerated;

  signalFile>>m0;
  cout << "m0 : " << m0 << endl;
  signalFile>>m12;
  cout << "m12 : " << m12 << endl;
  signalFile>>susyGenerated;
  cout << "SUSY Generated : " << susyGenerated << endl;

  // Luminosity*sigma = number of events 

  double crossSectionScaling = luminosity/susyGenerated;

  for(vector<TString>::iterator thisBin = binNames.begin() ; thisBin != binNames.end() ; thisBin++)
    {
      yields thisSignal,thisSignalError;
      double valueHolder;
      signalFile>>valueHolder ; thisSignal.zeroLepton = valueHolder*crossSectionScaling ;
      signalFile>>valueHolder ; thisSignalError.zeroLepton = valueHolder*crossSectionScaling ;
      signalFile>>valueHolder ; thisSignal.zeroLeptonLowDeltaPhiN = valueHolder*crossSectionScaling ;
      signalFile>>valueHolder ; thisSignalError.zeroLeptonLowDeltaPhiN = valueHolder*crossSectionScaling ;
      signalFile>>valueHolder ; thisSignal.oneLepton = valueHolder*crossSectionScaling ;
      signalFile>>valueHolder ; thisSignalError.oneLepton = valueHolder*crossSectionScaling ;
      signal[*thisBin]=thisSignal;
      signalError[*thisBin]=thisSignalError;
      cout << "For selection " << *thisBin << " signal is " << endl;
      cout << "zero lepton bin               : " << thisSignal.zeroLepton             << endl ;
      cout << "zero lepton low delta phi bin : " << thisSignal.zeroLeptonLowDeltaPhiN << endl ;
      cout << "one lepton bin                : " << thisSignal.oneLepton              << endl ;
      cout << "For selection " << *thisBin << " signal error is " << endl;
      cout << "zero lepton bin               : " << thisSignalError.zeroLepton             << endl ;
      cout << "zero lepton low delta phi bin : " << thisSignalError.zeroLeptonLowDeltaPhiN << endl ;
      cout << "one lepton bin                : " << thisSignalError.oneLepton              << endl ;
    }

  signalFile.close();

}

void setupObservations(TString binName, TString binFileName, map<TString,abcdBinParameters>& bins, map<TString,channels>& observations)
{
  channels counts;
  abcdBinParameters abcd;
  
  ifstream binFile;
       
  cout << "getting the file: " << binFileName << " to fill bin " << binName << endl;
       
  binFile.open(binFileName.Data(),fstream::in);
  
  string fileLine;
  
  TString index;
  double value;
       
  while(!binFile.eof())
    {
      getline(binFile,fileLine);
      TString thisLine(fileLine.c_str());

      TStringToken nameAndNumber(thisLine," ");
      nameAndNumber.NextToken();
      index = nameAndNumber;
      nameAndNumber.NextToken();
      value = nameAndNumber.Atof();
      cout << index << " : " << nameAndNumber << endl;
      if(index == "zeroLeptonCount"                                   ) counts.zeroLepton = value;		 
      else if(index == "zeroLeptonLowDeltaPhiNCount"                  ) counts.zeroLeptonLowDeltaPhiN = value;
      else if(index == "oneMuonCount"                                 ) counts.oneMuon = value;		 
      else if(index == "oneElectronCount"                             ) counts.oneElectron = value;           
      else if(index == "zeroLeptonTriggerEfficiency"		      ) abcd.zeroLeptonTriggerEfficiency = value;			
      else if(index == "zeroLeptonTriggerEfficiencyError"	      ) abcd.zeroLeptonTriggerEfficiencyError = value;		
      else if(index == "zeroLeptonLowDeltaPhiNTriggerEfficiency"      ) abcd.zeroLeptonLowDeltaPhiNTriggerEfficiency = value;	
      else if(index == "zeroLeptonLowDeltaPhiNTriggerEfficiencyError" ) abcd.zeroLeptonLowDeltaPhiNTriggerEfficiencyError = value;	
      else if(index == "oneElectronTriggerEfficiency"		      ) abcd.oneElectronTriggerEfficiency = value;			
      else if(index == "oneElectronTriggerEfficiencyError"    	      ) abcd.oneElectronTriggerEfficiencyError = value;		
      else if(index == "oneMuonTriggerEfficiency"	              ) abcd.oneMuonTriggerEfficiency = value;			
      else if(index == "oneMuonTriggerEfficiencyError"		      ) abcd.oneMuonTriggerEfficiencyError = value;			
      else if(index == "zeroLeptonLowDeltaPhiNMC"		      ) abcd.zeroLeptonLowDeltaPhiNMC = value;			
      else if(index == "ZtoNuNubTagScaling"			      ) abcd.ZtoNuNubTagScaling = value;				
      else if(index == "ZtoNuNubTagScalingError"		      ) abcd.ZtoNuNubTagScalingError = value;			
      else if(index == "ZtollAcceptance"			      ) abcd.ZtollAcceptance = value;				
      else if(index == "ZtollAcceptanceError"		              ) abcd.ZtollAcceptanceError = value;				
      else if(index == "qcdClosure"		        	      ) abcd.qcdClosure = value;					
      else if(index == "qcdClosureError"			      ) abcd.qcdClosureError = value;				
      else if(index == "topWJetsClosure"		              ) abcd.topWJetsClosure = value;				
      else if(index == "topWJetsClosureError"			      ) abcd.topWJetsClosureError = value;				
      else if(index == "ZtoNuNuClosure"				      ) abcd.ZtoNuNuClosure = value;					
      else if(index == "ZtoNuNuClosureError"                          ) abcd.ZtoNuNuClosureError = value;                            
    }

  binFile.close();

  bins[binName] = abcd;
  observations[binName] = counts;
}

void setupUnderlyingModel(map<TString,TString>& binFileNames, vector<TString>& binNames, TString& modelFileName , TString& binFilesFileName , allBins& numbers, double& luminosity)
{
  ifstream setupFile;
       
  cout << "getting the file: " << modelFileName << " to setup likelihood " << endl;
  
  string fileLine;
  
  TString index;
  double value;
       
  setupFile.open(modelFileName.Data(),fstream::in);

  while(!setupFile.eof())
    {
      getline(setupFile,fileLine);
      TString thisLine(fileLine.c_str());

      TStringToken nameAndNumber(thisLine," ");
      nameAndNumber.NextToken();
      index = nameAndNumber;
      nameAndNumber.NextToken();
      value = nameAndNumber.Atof();
      cout << index << " : " << value << endl;
      if(index == "diElectronCount"             ) numbers.diElectronCount = value;	     
      else if(index == "diMuonCount"	        ) numbers.diMuonCount = value;	      
      else if(index == "ZtollOverZtoNuNuRatio"  ) numbers.ZtollOverZtoNuNuRatio = value;
      else if(index == "ZtoeePurity"	        ) numbers.ZtoeePurity = value;	      	
      else if(index == "ZtoeePurityError"       ) numbers.ZtoeePurityError = value;     
      else if(index == "ZtomumuPurity"	        ) numbers.ZtomumuPurity = value;	      
      else if(index == "ZtomumuPurityError"     ) numbers.ZtomumuPurityError = value;   
      else if(index == "ZtoeeEfficiency"        ) numbers.ZtoeeEfficiency = value;	     
      else if(index == "ZtoeeEfficiencyError"   ) numbers.ZtoeeEfficiencyError = value; 
      else if(index == "ZtomumuEfficiency"      ) numbers.ZtomumuEfficiency = value;    
      else if(index == "ZtomumuEfficiencyError" ) numbers.ZtomumuEfficiencyError = value;	
      else if(index == "deltaPhiNRatio"	        ) numbers.deltaPhiNRatio = value;	     
      else if(index == "deltaPhiNRatioError"    ) numbers.deltaPhiNRatioError = value;  	  
      else if(index == "MCUncertainty" 	        ) numbers.MCUncertainty = value;	      	   
      else if(index == "Luminosity"   	        ) luminosity = value;	      	   
    }

  setupFile.close();

  setupFile.open(binFilesFileName.Data(),fstream::in);

  TString fileName;

  while(!setupFile.eof())
    {
      getline(setupFile,fileLine);
      TString thisLine(fileLine.c_str());

      TStringToken nameAndNumber(thisLine," ");
      nameAndNumber.NextToken();
      index = nameAndNumber;
      nameAndNumber.NextToken();
      fileName = nameAndNumber;
      cout << index << " : " << fileName << endl;
      binNames.push_back(index);
      binFileNames[index] = fileName;
    }
  setupFile.close();
}

void chooseUnderlyingUncertainties(allBins& numbers , map<TString,abcdBinParameters> bins, map<TString,yields> signal, map<TString,yields> signalError)
{
  double qcdAlpha(0.);
  double qcdBeta(0.);
  double topWJetsAlpha(0.);
  double topWJetsBeta(0.);
  double ZtoNuNuAlpha(0.);
  double ZtoNuNuBeta(0.);
  for(map<TString,abcdBinParameters>::iterator thisBin = bins.begin(); thisBin!= bins.end(); thisBin++)
    {
      double alpha,beta;
      betaPrimeModeTransform(thisBin->second.qcdClosure,thisBin->second.qcdClosureError,alpha,beta);
      if(alpha > qcdAlpha) qcdAlpha = alpha;
      if(beta  > qcdBeta ) qcdBeta  = beta ;
      betaPrimeModeTransform(thisBin->second.topWJetsClosure,thisBin->second.topWJetsClosureError,alpha,beta);
      if(alpha > topWJetsAlpha) topWJetsAlpha = alpha;
      if(beta  > topWJetsBeta ) topWJetsBeta  = beta ;
      betaPrimeModeTransform(thisBin->second.ZtoNuNuClosure,thisBin->second.ZtoNuNuClosureError,alpha,beta);
      if(alpha > ZtoNuNuAlpha) ZtoNuNuAlpha = alpha;
      if(beta  > ZtoNuNuBeta ) ZtoNuNuBeta  = beta ;
    }
  numbers.qcdClosure = (qcdAlpha - 1) / ( qcdBeta + 1 );
  numbers.qcdClosureError = sqrt( qcdAlpha * (qcdAlpha + qcdBeta - 1 ) / ( pow(qcdBeta - 1 , 2 ) * (qcdBeta - 2) ) );
  numbers.topWJetsClosure = (topWJetsAlpha - 1) / ( topWJetsBeta + 1 );
  numbers.topWJetsClosureError = sqrt( topWJetsAlpha * (topWJetsAlpha + topWJetsBeta - 1 ) / ( pow(topWJetsBeta - 1 , 2 ) * (topWJetsBeta - 2) ) );
  numbers.ZtoNuNuClosure = (ZtoNuNuAlpha - 1) / ( ZtoNuNuBeta + 1 );
  numbers.ZtoNuNuClosureError = sqrt( ZtoNuNuAlpha * (ZtoNuNuAlpha + ZtoNuNuBeta - 1 ) / ( pow(ZtoNuNuBeta - 1 , 2 ) * (ZtoNuNuBeta - 2) ) );

  map<TString,yields>::iterator thisSignal      = signal     .begin();
  map<TString,yields>::iterator thisSignalError = signalError.begin();

  double signalAlpha(0.);
  double signalBeta(0.);

  for(; thisSignal != signal.end() && thisSignalError != signalError.end(); thisSignal++, thisSignalError++)
    {
      double alpha,beta;
      betaPrimeModeTransform(thisSignal->second.zeroLepton,thisSignalError->second.zeroLepton,alpha,beta);
      if(alpha > signalAlpha) signalAlpha = alpha;
      if(beta  > signalBeta ) signalBeta  = beta ;
      betaPrimeModeTransform(thisSignal->second.zeroLeptonLowDeltaPhiN,thisSignalError->second.zeroLeptonLowDeltaPhiN,alpha,beta);
      if(alpha > signalAlpha) signalAlpha = alpha;
      if(beta  > signalBeta ) signalBeta  = beta ;
      betaPrimeModeTransform(thisSignal->second.oneLepton,thisSignalError->second.oneLepton,alpha,beta);
      if(alpha > signalAlpha) signalAlpha = alpha;
      if(beta  > signalBeta ) signalBeta  = beta ;
    }
  numbers.signal = (signalAlpha - 1) / ( signalBeta + 1 );
  numbers.signalError = sqrt( signalAlpha * (signalAlpha + signalBeta - 1 ) / ( pow(signalBeta - 1 , 2 ) * (signalBeta - 2) ) );
}



void buildLikelihood( TString setupFileName, TString binFilesFileName, TString signalModelFileName, int signalModelFileLine, TString workspaceName, TString outputFileName ) 
{

  double luminosity(1.);
  RooWorkspace ws (workspaceName) ;
  ws.autoImportClassCode(true);
  vector<TString> binNames;
  map<TString,TString> binFileNames;
  allBinNames names;
  allBins numbers;
  map<TString,abcdBinParameters> bins;
  map<TString,channels> observations;
  map<TString,yields> signal;
  map<TString,yields> signalError;
  double oneLeptonTotal(0.);
  double zeroLeptonTopWJetsGuess(0.);

  setupUnderlyingModel(binFileNames, binNames, setupFileName , binFilesFileName , numbers , luminosity);
  for(map<TString,TString>::iterator thisBin = binFileNames.begin(); thisBin != binFileNames.end() ; thisBin++)
    {
      setupObservations(thisBin->first , thisBin->second , bins, observations);
    }

  setupSignalModel(binNames , signalModelFileName , signalModelFileLine , signal , signalError , luminosity);

  chooseUnderlyingUncertainties(numbers , bins, signal, signalError);

  setupUnderlyingLikelihood(ws , names, numbers);

  for(vector<TString>::iterator thisBin = binNames.begin(); thisBin != binNames.end() ; thisBin++)
    {
      makeOneBin(ws , *thisBin , names , observations[*thisBin] , bins[*thisBin] , signal[*thisBin] , signalError[*thisBin] , oneLeptonTotal, zeroLeptonTopWJetsGuess );
    }

  ws.var(names.singleLeptonScaling)->setVal(zeroLeptonTopWJetsGuess/oneLeptonTotal);

  RooArgSet allpdfs = ws.allPdfs();
  
  RooProdPdf likelihood("likelihood","likelihood",allpdfs);
  
  ws.import(likelihood, RecycleConflictNodes());

  ws.defineSet("poi",names.signalCrossSection);  

  RooDataSet data("data","data",*ws.set(names.observables));

  ws.import(data);

  cout << "setting up models" << endl;

  ModelConfig sbModel("S+B_model",&ws);
  sbModel.SetPdf(*ws.pdf("likelihood"));
  sbModel.SetObservables(*ws.set(names.observables));
  sbModel.SetNuisanceParameters(*ws.set(names.nuisances));
  sbModel.SetParametersOfInterest(*ws.set("poi"));
  sbModel.SetProtoData(*ws.data("data"));

  ModelConfig bModel("B_model",&ws);
  bModel.SetPdf(*ws.pdf("likelihood"));
  bModel.SetObservables(*ws.set(names.observables));
  bModel.SetNuisanceParameters(*ws.set(names.nuisances));
  bModel.SetParametersOfInterest(*ws.set("poi"));
  bModel.SetProtoData(*ws.data("data"));
  ws.var(names.signalCrossSection)->setVal(0.0);
  bModel.SetSnapshot(*ws.set("poi"));

  ws.import (sbModel);
  ws.import (bModel);

  ws.Print() ;
  ws.writeToFile(outputFileName) ;

}


