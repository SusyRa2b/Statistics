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
  TString diMuonName;
  double diMuon;
  TString diElectronName;
  double diElectron;
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
  TString ZtoNuNubTagScalingName;
  double ZtoNuNubTagScaling;
  double ZtoNuNubTagScalingError;
  TString ZtomumuAcceptanceName;
  double ZtomumuAcceptance;
  double ZtomumuAcceptanceError;
  TString ZtoeeAcceptanceName;
  double ZtoeeAcceptance;
  double ZtoeeAcceptanceError;
  TString deltaPhiNRatioName;
  double deltaPhiNRatio;
  double deltaPhiNRatioError;
  double qcdClosure;
  double qcdClosureError;
  double topWJetsClosure;
  double topWJetsClosureError;
  TString ZtoeeSystematicName;
  double ZtoeeSystematic;
  double ZtoeeSystematicError;
  TString ZtomumuSystematicName;
  double ZtomumuSystematic;
  double ZtomumuSystematicError;
} ;

struct allBinNames
{
  TString ZtoeeInvPurity;
  TString ZtomumuInvPurity;
  TString ZtoeeEfficiency;
  TString ZtomumuEfficiency;
  TString ZtollOverZtoNuNuRatio;
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
} ;

struct allBins
{
  double ZtollOverZtoNuNuRatio;
  double ZtoeePurity;
  double ZtoeePurityError;
  double ZtomumuPurity;
  double ZtomumuPurityError;
  double ZtoeeEfficiency;
  double ZtoeeEfficiencyError;
  double ZtomumuEfficiency;
  double ZtomumuEfficiencyError;
  double qcdClosure;
  double qcdClosureError;
  double topWJetsClosure;
  double topWJetsClosureError;
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

RooAbsArg* getCorrelatedBetaConstraint(RooWorkspace& ws,const TString varName,const TString binName,
				       const double value ,const double error,
				       const TString correlatedPassObs,const TString correlatedFailObs,
				       const TString correlatedPassPar,const TString correlatedFailPar)
{
  double alpha,beta;
  betaModeTransform(value , error , alpha , beta );

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

bool makeOneBin(RooWorkspace& ws , TString& binName , allBinNames& names , const allBins& numbers, channels& observed , abcdBinParameters& abcd , yields& signal , yields& signalError , double& oneLeptonTotal, double& zeroLeptonTopWJetsGuess)
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

  //Define unique counts and add them to the workspace

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

  RooAbsArg* deltaPhiNRatio = ws.arg("deltaPhiNRatio_"+abcd.deltaPhiNRatioName+"_Ratio");
  if(deltaPhiNRatio == NULL) 
    {
      deltaPhiNRatio = 
	getBetaConstraint(ws,"deltaPhiNRatio_",abcd.deltaPhiNRatioName,
			  abcd.deltaPhiNRatio,abcd.deltaPhiNRatioError,
			  names.observables,names.nuisances);
    }
  double qcdGuess = observed.zeroLeptonLowDeltaPhiN - abcd.zeroLeptonLowDeltaPhiNMC;
  if(qcdGuess < 0) qcdGuess = 0;
  RooRealVar zeroLeptonLowDeltaPhiNQCDYield(zeroLeptonLowDeltaPhiNName+"_QCDYield",zeroLeptonLowDeltaPhiNName+"_QCDYield",qcdGuess,1e-5,10000);
  ws.import(zeroLeptonLowDeltaPhiNQCDYield);
  ws.extendSet(names.nuisances,zeroLeptonLowDeltaPhiNQCDYield.GetName());
  RooAbsArg* zeroLeptonQCDClosure = 
    getBetaPrimeConstraint(ws,"zeroLeptonQCDClosure_", binName,
			   abcd.qcdClosure,abcd.qcdClosureError,
			   names.observables,names.nuisances);
  //RooAbsArg* zeroLeptonQCDClosure = getCorrelatedBetaPrimeConstraint(ws,"zeroLeptonQCDClosure_",binName, abcd.qcdClosure,abcd.qcdClosureError, names.qcdClosurePassObs, names.qcdClosureFailObs, names.qcdClosurePassPar,names.qcdClosureFailPar);//Correlated
  RooProduct zeroLeptonQCDYield(zeroLeptonName+"_QCDYield",zeroLeptonName+"_QCDYield",RooArgSet(*deltaPhiNRatio,*zeroLeptonQCDClosure,zeroLeptonLowDeltaPhiNQCDYield));
  
  //Define top and W+jets component:
  
  RooRealVar* singleLeptonScaling = ws.var(names.singleLeptonScaling);
  RooRealVar  oneLeptonTopWJetsYield(oneLeptonName+"_TopWJetsYield",oneLeptonName+"_TopWJetsYield",0.5*(observed.oneMuon+observed.oneElectron),1e-5,100000);
  ws.import(oneLeptonTopWJetsYield);
  ws.extendSet(names.nuisances,oneLeptonTopWJetsYield.GetName());
  oneLeptonTotal += oneLeptonTopWJetsYield.getVal();
  RooAbsArg*  zeroLeptonTopWJetsClosure = 
    getBetaPrimeConstraint(ws,"zeroLeptonTopWJetsClosure_", binName,
			   abcd.topWJetsClosure,abcd.topWJetsClosureError,
			   names.observables,names.nuisances);
  //RooAbsArg*  zeroLeptonTopWJetsClosure = getCorrelatedBetaPrimeConstraint(ws,"zeroLeptonTopWJetsClosure_",binName, abcd.topWJetsClosure,abcd.topWJetsClosureError, names.topWJetsClosurePassObs,names.topWJetsClosureFailObs, names.topWJetsClosurePassPar,names.topWJetsClosureFailPar);//Correlated
  RooProduct  zeroLeptonTopWJetsYield(zeroLeptonName+"_TopWJetsYield",zeroLeptonName+"_TopWJetsYield",RooArgSet(*singleLeptonScaling,*zeroLeptonTopWJetsClosure,oneLeptonTopWJetsYield));
  
  //Define Z to invisible component:
  
  //-----Define Z->ll observables
  
  RooRealVar* diElectronCount = (RooRealVar*)ws.arg("diElectron_"+observed.diElectronName+"_Count");
  if(diElectronCount == NULL) 
    {
      diElectronCount = new RooRealVar("diElectron_"+observed.diElectronName+"_Count","diElectron_"+observed.diElectronName+"_Count",observed.diElectron); 
      diElectronCount->setConstant();
      ws.import(*diElectronCount);
      ws.extendSet(names.observables,diElectronCount->GetName());
    }
  RooRealVar* diMuonCount = (RooRealVar*)ws.arg("diMuon_"+observed.diMuonName+"_Count");
  if(diMuonCount == NULL) 
    {
      diMuonCount = new RooRealVar("diMuon_"+observed.diMuonName+"_Count","diMuon_"+observed.diMuonName+"_Count",observed.diMuon);
      diMuonCount->setConstant();
      ws.import(*diMuonCount);
      ws.extendSet(names.observables,diMuonCount->GetName());
    }
  
  //-----Define Z->ll yields
  
  RooAbsArg* ZtollOverZtoNuNuRatio = ws.arg(names.ZtollOverZtoNuNuRatio);
  
  RooRealVar* ZtoNuNu = (RooRealVar*)ws.arg("ZtoNuNu_"+observed.diMuonName);//use MuonName, which should be the same as ElectronName
  if(ZtoNuNu == NULL) 
    {
      ZtoNuNu = new RooRealVar("ZtoNuNu_"+observed.diMuonName,"ZtoNuNu_"+observed.diMuonName,numbers.ZtollOverZtoNuNuRatio*0.5*(observed.diMuon*numbers.ZtomumuEfficiency*abcd.ZtomumuAcceptance/numbers.ZtomumuPurity + observed.diElectron*numbers.ZtoeeEfficiency*abcd.ZtoeeAcceptance/numbers.ZtoeePurity),0.0,1e5);
      ws.import(*ZtoNuNu);
      ws.extendSet(names.nuisances,ZtoNuNu->GetName());
  }
  
  RooProduct* Ztoll = (RooProduct*)ws.arg("Ztoll_"+observed.diMuonName);//use MuonName, which should be the same as ElectronName
  if(Ztoll == NULL) 
    {
      Ztoll = new RooProduct("Ztoll_"+observed.diMuonName,"Ztoll_"+observed.diMuonName,RooArgSet(*ZtoNuNu,*ZtollOverZtoNuNuRatio));
      ws.import(*Ztoll, RecycleConflictNodes());
    }
  
  RooAbsArg* ZtomumuEfficiency = ws.arg(names.ZtomumuEfficiency);
  RooAbsArg* ZtomumuInvPurity = ws.arg(names.ZtomumuInvPurity);
  RooAbsArg* ZtoeeEfficiency = ws.arg(names.ZtoeeEfficiency);
  RooAbsArg* ZtoeeInvPurity = ws.arg(names.ZtoeeInvPurity);
  
  RooAbsArg* ZtoeeAcceptance = ws.arg("ZtoeeAcceptance_"+abcd.ZtoeeAcceptanceName+"_Ratio");
  if(ZtoeeAcceptance == NULL) 
    {
      ZtoeeAcceptance = 
	getBetaConstraint(ws,"ZtoeeAcceptance_",abcd.ZtoeeAcceptanceName,
			  abcd.ZtoeeAcceptance,abcd.ZtoeeAcceptanceError,
			  names.observables,names.nuisances);
    }
  RooAbsArg* ZtomumuAcceptance = ws.arg("ZtomumuAcceptance_"+abcd.ZtomumuAcceptanceName+"_Ratio");
  if(ZtomumuAcceptance == NULL) 
    {
      ZtomumuAcceptance = 
	getBetaConstraint(ws,"ZtomumuAcceptance_",abcd.ZtomumuAcceptanceName,
			  abcd.ZtomumuAcceptance,abcd.ZtomumuAcceptanceError,
			  names.observables,names.nuisances);
    }
  
  /*
  //BEN FIXME - systematics not implemented yet
  RooAbsArg*  ZtoeeSystematic = ws.arg("ZtoeeSystematic_"+abcd.ZtoeeSystematicName+"_Ratio");
  if(ZtoeeSystematic == NULL) 
    {
      ZtoeeSystematic = 
	getBetaPrimeConstraint(ws,"ZtoeeSystematic_",abcd.ZtoeeSystematicName,
			       abcd.ZtoeeSystematic,abcd.ZtoeeSystematicError,
			       names.observables,names.nuisances);
    }
  RooAbsArg*  ZtomumuSystematic = ws.arg("ZtomumuSystematic_"+abcd.ZtomumuSystematicName+"_Ratio");
  if(ZtomumuSystematic == NULL) 
    {
      ZtomumuSystematic = 
	getBetaPrimeConstraint(ws,"ZtomumuSystematic_",abcd.ZtomumuSystematicName,
			       abcd.ZtomumuSystematic,abcd.ZtomumuSystematicError,
			       names.observables,names.nuisances);
    }
  */

  RooProduct* diMuonYield = (RooProduct*)ws.arg("diMuon_"+observed.diMuonName+"_Yield");//Assumes acceptance is only binned in zero or more dimensions of the count.
  if(diMuonYield == NULL) 
    {
      diMuonYield = new RooProduct("diMuon_"+observed.diMuonName+"_Yield","diMuon_"+observed.diMuonName+"_Yield",RooArgSet(*Ztoll,*ZtomumuAcceptance,*ZtomumuEfficiency,*ZtomumuInvPurity));
      ws.import(*diMuonYield, RecycleConflictNodes());
    }
  
  RooProduct* diElectronYield = (RooProduct*)ws.arg("diElectron_"+observed.diElectronName+"_Yield");//Assumes acceptance is only binned in zero or more dimensions of the count.
  if(diElectronYield == NULL) 
    {
      diElectronYield = new RooProduct("diElectron_"+observed.diElectronName+"_Yield","diElectron_"+observed.diElectronName+"_Yield",RooArgSet(*Ztoll,*ZtoeeAcceptance,*ZtoeeEfficiency,*ZtoeeInvPurity));
      ws.import(*diElectronYield, RecycleConflictNodes());
    }
  
  //-----Define Z->ll Poisson constraints
  
  RooPoisson* diMuonConstraint = (RooPoisson*)ws.arg("diMuon_"+observed.diMuonName+"_Constraint");
  if(diMuonConstraint == NULL) 
    {
      diMuonConstraint = new RooPoisson("diMuon_"+observed.diMuonName+"_Constraint","diMuon_"+observed.diMuonName+"_Constraint",*diMuonYield,*diMuonCount);
      ws.import(*diMuonConstraint, RecycleConflictNodes());
    }
  
  RooPoisson* diElectronConstraint = (RooPoisson*)ws.arg("diElectron_"+observed.diElectronName+"_Constraint");
  if(diElectronConstraint == NULL) 
    {
      diElectronConstraint = new RooPoisson("diElectron_"+observed.diElectronName+"_Constraint","diElectron_"+observed.diElectronName+"_Constraint",*diElectronYield,*diElectronCount);
      ws.import(*diElectronConstraint, RecycleConflictNodes());
    }
  
  //-----Define Z->nunu yield
  
  RooAbsArg* zeroLeptonZtoNuNubTagScaling = ws.arg("zeroLeptonZtoNuNubTagScaling_"+abcd.ZtoNuNubTagScalingName+"_Ratio");
  if(zeroLeptonZtoNuNubTagScaling == NULL) 
    {
      zeroLeptonZtoNuNubTagScaling = getBetaConstraint(ws,"zeroLeptonZtoNuNubTagScaling_",abcd.ZtoNuNubTagScalingName,
						       abcd.ZtoNuNubTagScaling,abcd.ZtoNuNubTagScalingError,
						       names.observables,names.nuisances);
    }
  
  RooProduct zeroLeptonZtoNuNuYield(zeroLeptonName+"_ZtoNuNuYield",zeroLeptonName+"_ZtoNuNuYield",RooArgSet(*zeroLeptonZtoNuNubTagScaling,*ZtoNuNu));
  
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

  diElectronCount = 0;
  diMuonCount = 0;
  ZtoNuNu = 0;
  Ztoll = 0;
  diElectronYield = 0;
  diMuonYield = 0;
  diElectronConstraint = 0;
  diMuonYield = 0;
  delete diElectronCount;
  delete diMuonCount;
  delete ZtoNuNu;
  delete Ztoll;
  delete diElectronYield;
  delete diMuonYield;
  delete diElectronConstraint;
  delete diMuonYield;

  return true;

}

void setupUnderlyingLikelihood(RooWorkspace& ws ,allBinNames& names, allBins& numbers)
{
  ws.defineSet("observables","");
  names.observables = "observables";

  ws.defineSet("nuisances","");
  names.nuisances = "nuisances";

  //Universal parameters

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
  /*
  //Commented out for now to assume they are uncorrelated
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
  */

  //Objects for Z->invisible background:

  RooRealVar ZtollOverZtoNuNuRatio("ZtollOverZtoNuNuRatio","ZtollOverZtoNuNuRatio",numbers.ZtollOverZtoNuNuRatio);
  ZtollOverZtoNuNuRatio.setConstant();
  ws.import(ZtollOverZtoNuNuRatio);
  names.ZtollOverZtoNuNuRatio = ZtollOverZtoNuNuRatio.GetName();

  RooAbsArg* ZtomumuInvPurity = 
    getInverseBetaConstraint(ws,"ZtomumuInvPurity","",
			     numbers.ZtomumuPurity,numbers.ZtomumuPurityError,
			     names.observables,names.nuisances);
  names.ZtomumuInvPurity = ZtomumuInvPurity->GetName();
  
  RooAbsArg* ZtomumuEfficiency = 
    getBetaConstraint(ws,"ZtomumuEfficiency","",
		      numbers.ZtomumuEfficiency,numbers.ZtomumuEfficiencyError,
		      names.observables,names.nuisances);
  names.ZtomumuEfficiency = ZtomumuEfficiency->GetName();
  
  RooAbsArg* ZtoeeInvPurity = 
    getInverseBetaConstraint(ws,"ZtoeeInvPurity","",
			     numbers.ZtoeePurity,numbers.ZtoeePurityError,
			     names.observables,names.nuisances);
  names.ZtoeeInvPurity = ZtoeeInvPurity->GetName();

  RooAbsArg* ZtoeeEfficiency = 
    getBetaConstraint(ws,"ZtoeeEfficiency","",
		      numbers.ZtoeeEfficiency,numbers.ZtoeeEfficiencyError,
		      names.observables,names.nuisances);
  names.ZtoeeEfficiency = ZtoeeEfficiency->GetName();

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
  TString value;
  
  while(!binFile.eof())
    {
      getline(binFile,fileLine);
      TString thisLine(fileLine.c_str());

      TStringToken nameAndNumber(thisLine," ");
      nameAndNumber.NextToken();
      index = nameAndNumber;
      nameAndNumber.NextToken();
      value = nameAndNumber;
      if(value.IsDigit() == false) 
	{
	  //hack to get rid of potential whitespace in strings
	  string str = (string)value;
	  for(unsigned int i=0; i<str.length(); i++) 
	    {
	    if(str[i] == '\t') str.erase(i,1);
	    if(str[i] == ' ') str.erase(i,1);
	    }
	  value = (TString)str;
	}
      cout << index << " : " << value << endl;
      if(index == "zeroLeptonCount"                                   ) counts.zeroLepton = value.Atof();		 
      else if(index == "zeroLeptonLowDeltaPhiNCount"                  ) counts.zeroLeptonLowDeltaPhiN = value.Atof();
      else if(index == "oneMuonCount"                                 ) counts.oneMuon = value.Atof();		 
      else if(index == "oneElectronCount"                             ) counts.oneElectron = value.Atof();           
      else if(index == "diElectronCountName"                          ) counts.diElectronName = value;	     
      else if(index == "diElectronCount"                              ) counts.diElectron = value.Atof();	     
      else if(index == "diMuonCountName"	                      ) counts.diMuonName = value;	      
      else if(index == "diMuonCount"	                              ) counts.diMuon = value.Atof();	      
      else if(index == "zeroLeptonTriggerEfficiency"		      ) abcd.zeroLeptonTriggerEfficiency = value.Atof();			
      else if(index == "zeroLeptonTriggerEfficiencyError"	      ) abcd.zeroLeptonTriggerEfficiencyError = value.Atof();		
      else if(index == "zeroLeptonLowDeltaPhiNTriggerEfficiency"      ) abcd.zeroLeptonLowDeltaPhiNTriggerEfficiency = value.Atof();	
      else if(index == "zeroLeptonLowDeltaPhiNTriggerEfficiencyError" ) abcd.zeroLeptonLowDeltaPhiNTriggerEfficiencyError = value.Atof();	
      else if(index == "oneElectronTriggerEfficiency"		      ) abcd.oneElectronTriggerEfficiency = value.Atof();			
      else if(index == "oneElectronTriggerEfficiencyError"    	      ) abcd.oneElectronTriggerEfficiencyError = value.Atof();		
      else if(index == "oneMuonTriggerEfficiency"	              ) abcd.oneMuonTriggerEfficiency = value.Atof();			
      else if(index == "oneMuonTriggerEfficiencyError"		      ) abcd.oneMuonTriggerEfficiencyError = value.Atof();			
      else if(index == "zeroLeptonLowDeltaPhiNMC"		      ) abcd.zeroLeptonLowDeltaPhiNMC = value.Atof();			
      else if(index == "ZtoNuNubTagScalingName"			      ) abcd.ZtoNuNubTagScalingName = value;				
      else if(index == "ZtoNuNubTagScaling"			      ) abcd.ZtoNuNubTagScaling = value.Atof();				
      else if(index == "ZtoNuNubTagScalingError"		      ) abcd.ZtoNuNubTagScalingError = value.Atof();			
      else if(index == "ZtomumuAcceptanceName"			      ) abcd.ZtomumuAcceptanceName = value;				
      else if(index == "ZtomumuAcceptance"			      ) abcd.ZtomumuAcceptance = value.Atof();				
      else if(index == "ZtomumuAcceptanceError"		              ) abcd.ZtomumuAcceptanceError = value.Atof();				
      else if(index == "ZtoeeAcceptanceName"			      ) abcd.ZtoeeAcceptanceName = value;				
      else if(index == "ZtoeeAcceptance"			      ) abcd.ZtoeeAcceptance = value.Atof();				
      else if(index == "ZtoeeAcceptanceError"		              ) abcd.ZtoeeAcceptanceError = value.Atof();				
      else if(index == "deltaPhiNRatioName"                           ) abcd.deltaPhiNRatioName = value;	     
      else if(index == "deltaPhiNRatio"	                              ) abcd.deltaPhiNRatio = value.Atof();	     
      else if(index == "deltaPhiNRatioError"                          ) abcd.deltaPhiNRatioError = value.Atof();  	  
      else if(index == "qcdClosure"		        	      ) abcd.qcdClosure = value.Atof();					
      else if(index == "qcdClosureError"			      ) abcd.qcdClosureError = value.Atof();				
      else if(index == "topWJetsClosure"		              ) abcd.topWJetsClosure = value.Atof();				
      else if(index == "topWJetsClosureError"			      ) abcd.topWJetsClosureError = value.Atof();				
      else if(index == "ZtoeeSystematicName"		       	      ) abcd.ZtoeeSystematicName = value;					
      else if(index == "ZtoeeSystematic"			      ) abcd.ZtoeeSystematic = value.Atof();					
      else if(index == "ZtoeeSystematicError"                         ) abcd.ZtoeeSystematicError = value.Atof();                            
      else if(index == "ZtomumuSystematicName"		       	      ) abcd.ZtomumuSystematicName = value;					
      else if(index == "ZtomumuSystematic"			      ) abcd.ZtomumuSystematic = value.Atof();					
      else if(index == "ZtomumuSystematicError"                       ) abcd.ZtomumuSystematicError = value.Atof();                            
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
      if(index == "ZtollOverZtoNuNuRatio"  ) numbers.ZtollOverZtoNuNuRatio = value;
      else if(index == "ZtoeePurity"	        ) numbers.ZtoeePurity = value;	      	
      else if(index == "ZtoeePurityError"       ) numbers.ZtoeePurityError = value;     
      else if(index == "ZtomumuPurity"	        ) numbers.ZtomumuPurity = value;	      
      else if(index == "ZtomumuPurityError"     ) numbers.ZtomumuPurityError = value;   
      else if(index == "ZtoeeEfficiency"        ) numbers.ZtoeeEfficiency = value;	     
      else if(index == "ZtoeeEfficiencyError"   ) numbers.ZtoeeEfficiencyError = value; 
      else if(index == "ZtomumuEfficiency"      ) numbers.ZtomumuEfficiency = value;    
      else if(index == "ZtomumuEfficiencyError" ) numbers.ZtomumuEfficiencyError = value;	
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
  //Stuff for correlated closure systematics 
  //--Currently not used because we assume they are uncorrelated!
  double qcdAlpha(0.);
  double qcdBeta(0.);
  double topWJetsAlpha(0.);
  double topWJetsBeta(0.);
   for(map<TString,abcdBinParameters>::iterator thisBin = bins.begin(); thisBin!= bins.end(); thisBin++)
    {
      double alpha,beta;
      betaPrimeModeTransform(thisBin->second.qcdClosure,thisBin->second.qcdClosureError,alpha,beta);
      if(alpha > qcdAlpha) qcdAlpha = alpha;
      if(beta  > qcdBeta ) qcdBeta  = beta ;
      betaPrimeModeTransform(thisBin->second.topWJetsClosure,thisBin->second.topWJetsClosureError,alpha,beta);
      if(alpha > topWJetsAlpha) topWJetsAlpha = alpha;
      if(beta  > topWJetsBeta ) topWJetsBeta  = beta ;
    }
  numbers.qcdClosure = (qcdAlpha - 1) / ( qcdBeta + 1 );
  numbers.qcdClosureError = sqrt( qcdAlpha * (qcdAlpha + qcdBeta - 1 ) / ( pow(qcdBeta - 1 , 2 ) * (qcdBeta - 2) ) );
  numbers.topWJetsClosure = (topWJetsAlpha - 1) / ( topWJetsBeta + 1 );
  numbers.topWJetsClosureError = sqrt( topWJetsAlpha * (topWJetsAlpha + topWJetsBeta - 1 ) / ( pow(topWJetsBeta - 1 , 2 ) * (topWJetsBeta - 2) ) );
  
  //Signal
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
      makeOneBin(ws , *thisBin , names , numbers, observations[*thisBin] , bins[*thisBin] , signal[*thisBin] , signalError[*thisBin] , oneLeptonTotal, zeroLeptonTopWJetsGuess );
    }

  ws.var(names.singleLeptonScaling)->setVal(zeroLeptonTopWJetsGuess/oneLeptonTotal);

  RooArgSet allpdfs = ws.allPdfs();
  
  cout << "*#*#*#*#*#* allpdfs *#*#*#*#*#*" << endl;
  cout << "size: " << allpdfs.getSize() << endl;
  allpdfs.Print("v");
  
  RooProdPdf likelihood("likelihood","likelihood",allpdfs);
  
  ws.import(likelihood, RecycleConflictNodes());

  ws.defineSet("poi",names.signalCrossSection);  

  cout << "*#*#*#*#*#* poi *#*#*#*#*#*" << endl;
  cout << "size: " << (*ws.set("poi")).getSize() << endl;
  (*ws.set("poi")).Print("v");

  cout << "*#*#*#*#*#* data *#*#*#*#*#*" << endl;
  cout << "size: " << (*ws.set(names.observables)).getSize() << endl;
  (*ws.set(names.observables)).Print("v");
  //same output as allpdfs, except Constraint->Obs or Count

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


