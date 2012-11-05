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


//===================================================================================================================================

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

void betaPrimeMeanTransform(const double& mu_measured , const double& sigma_measured ,
			    double& alpha , double& beta ) 
 {
   double factor = (mu_measured*(mu_measured+1))/(sigma_measured*sigma_measured); // common factor of mu(mu+1)/sigma**2
   beta = factor + 2;
   alpha =  mu_measured*(factor + 1);
  cout << "Beta Prime Mean Transform " << endl;
  cout << "The calculated value of alpha is : " << alpha << endl;
  cout << "The calculated value of beta is  : " << beta << endl;
  cout << "The input value of mu is         : " << mu_measured << endl;
  cout << "The input value of sigma is      : " << sigma_measured << endl;
  cout << "The calculated value of mu is    : " << (alpha ) / ( beta - 1 ) << endl;
  cout << "The calculated value of sigma is : " << sqrt( alpha * (alpha + beta - 1 ) / ( pow(beta - 1 , 2 ) * (beta - 2) ) ) << endl;
 }

void gammaMeanTransform(const double& mu_measured , const double& sigma_measured ,
                        double& k , double& theta )
 {
  double mu2 = mu_measured*mu_measured; 					         
  double sigma2 = sigma_measured*sigma_measured;				         
  // from Wikipedia's http://en.wikipedia.org/wiki/Gamma_distribution parameterization.  
  theta = sigma2/mu_measured;							         
  k = mu2/sigma2;								         

  cout << "For mu    : " << mu_measured << endl;				         
  cout << "For sigma : " << sigma_measured << endl;				         
  cout << "k is      : " << k << endl;						         
  cout << "theta is  : " << theta << endl;	
  
  cout << "Gamma Mean Transform " << endl;
  cout << "The calculated value of alpha is : " << k << endl;
  cout << "The calculated value of beta is  : " << theta << endl;
  cout << "The input value of mu is         : " << mu_measured << endl;
  cout << "The input value of sigma is      : " << sigma_measured << endl;
  cout << "The calculated value of mu is    : " << (k - 1) / ( theta + 1 ) << endl;
  cout << "The calculated value of sigma is : " << sqrt( k * (k + theta - 1 ) / ( pow(theta - 1 , 2 ) * (theta - 2) ) ) << endl;

  				         
}

//===================================================================================================================================

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


void betaMeanTransform(const double& mu_measured , const double& sigma_measured ,
                   double& alpha , double& beta )
{
  double mu2 = mu_measured*mu_measured;
  double mu3 = mu_measured*mu2;
  double sigma2 = sigma_measured*sigma_measured;

  alpha = (mu2-mu3-mu_measured*sigma2)/sigma2;
  beta = (mu_measured-1.)*(mu2-mu_measured+sigma2)/sigma2;

  if( !(alpha >= 1) && !(beta >= 1) )
    {
      cout << "mean and variance un-physical with beta pdf, try setting mode to central value instead" << endl;
      betaModeTransform(mu_measured, sigma_measured, alpha, beta );
    }
  cout << "Beta Mean Transform " << endl;
  cout << "The calculated value of alpha is : " << alpha << endl;
  cout << "The calculated value of beta is  : " << beta << endl;
  cout << "The input value of mu is         : " << mu_measured << endl;
  cout << "The input value of sigma is      : " << sigma_measured << endl;
  cout << "The calculated value of mu from mean is    : " << alpha/(alpha+beta) << endl;
  cout << "The calculated value of mu from mode is    : " << (alpha - 1) / ( alpha + beta - 2 ) << endl;
  cout << "The calculated value of sigma is : " << sqrt( alpha * beta / ( pow(alpha + beta,2) * (alpha + beta + 1) ) ) << endl;



}


//===================================================================================================================================

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
  RooRealVar failObservable (varName+binName+"_FailObs", varName+binName+"_FailObs", beta +1,1e-5,1e5);
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

//===================================================================================================================================


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
  RooRealVar failObservable (varName+binName+"_FailObs", varName+binName+"_FailObs", alpha+beta-2,1e-5,1e5);
  failObservable.setConstant();

  ws.import(passObservable);
  ws.import(failObservable);
  ws.extendSet(observables,passObservable.GetName());
  ws.extendSet(observables,failObservable.GetName());

  RooRealVar passParameter (varName+binName+"_PassPar", varName+binName+"_PassPar", alpha-1,1e-5,1e5);
  RooRealVar failParameter (varName+binName+"_FailPar", varName+binName+"_FailPar", alpha+beta -2,1e-5,1e5);

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

//===================================================================================================================================


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

//===================================================================================================================================

