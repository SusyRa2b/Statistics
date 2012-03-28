//
//   Owen Long, UCR
//   Harrison Prosper, FSU
//   Sezen Sekmen, FSU
//
//   a.g. : implementing trigger efficiency corrections 
//

#include "ra2bRoostatsClass7_noNSig_withKristen.h"

#include <iostream>
#include <string.h>
#include <complex>

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

#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooTrace.h"
#include "RooUniform.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooGammaPdf.h"
#include "RooBetaPdf.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"

#include "RooCFunction4Binding.h" 
#include "RooCFunction3Binding.h"
#include "RooTFnBinding.h" 
#include "TMath.h"

#include "TF1.h" 
#include "Math/WrappedTF1.h"
#include "Math/RootFinderAlgorithms.h"


//////////#include "LikelihoodIntervalPlot.cxx"

using namespace RooFit ;
using namespace RooStats ;
using namespace ROOT::Math ;


class gammaScalingFunctor
{
private:
  
  double currentSigma,desiredSigma;
  
public:

  gammaScalingFunctor(const double& currentSigma_,const double& desiredSigma_):currentSigma(currentSigma_),desiredSigma(desiredSigma_){};
  double sigmaModeCalc(double alpha)
  {
    long double scale = pow(-1 + sqrt(1 + 4*pow(currentSigma,2)),2*alpha)*pow(4,-alpha);
    long double linGamma = tgammal((long double) ((2 + 4*alpha + (1 + sqrt(1 + 4*pow(currentSigma,2)))/pow(currentSigma,2))/2.))/
      tgammal((long double) ((1 + 2*pow(currentSigma,2) + sqrt(1 + 4*pow(currentSigma,2)))/(2.*pow(currentSigma,2))));
    long double quadGamma = tgammal((long double) (alpha + (1 + 2*pow(currentSigma,2) + sqrt(1 + 4*pow(currentSigma,2)))/(2.*pow(currentSigma,2))))/
      tgammal((long double)((1 + 2*pow(currentSigma,2) + sqrt(1 + 4*pow(currentSigma,2)))/(2.*pow(currentSigma,2))));
    return sqrt(scale*(linGamma-quadGamma*quadGamma));
  };
  double sigmaModeIntersect(double* alpha, double*)
  {
    return sigmaModeCalc(*alpha) - desiredSigma;
  };

  double sigmaMeanCalc(double alpha)
  {
    long double scale = pow(pow(currentSigma,2),2*alpha) ;
    long double linGamma = tgammal(2*alpha + pow(currentSigma,-2))/tgammal(pow(currentSigma,-2));
    long double quadGamma = tgammal(alpha + pow(currentSigma,-2))/tgammal(pow(currentSigma,-2));
    return sqrt(scale*(linGamma-quadGamma*quadGamma));
  };
  double sigmaMeanIntersect(double* alpha, double*)
  {
    return sigmaMeanCalc(*alpha) - desiredSigma;
  };

  void setCurrentSigma(double currentSigma_){currentSigma=currentSigma_;};
  void setDesiredSigma(double desiredSigma_){desiredSigma=desiredSigma_;};
};

void gammaModeTransform(const double& mu_measured , const double& sigma_measured ,
		    double& k , double& theta ) 
{
  double mu2 = mu_measured*mu_measured;
  double sigma2 = sigma_measured*sigma_measured;
  // from Wikipedia's http://en.wikipedia.org/wiki/Gamma_distribution parameterization.
  theta = 0.5*(sqrt(mu2+4*sigma2)-mu_measured);
  k = sigma2/(theta*theta);

  cout << "For mu    : " << mu_measured << endl;
  cout << "For sigma : " << sigma_measured << endl;
  cout << "k is      : " << k << endl;
  cout << "theta is  : " << theta << endl;
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
}

void gammaCousinsTransform(const double& mu_measured , const double& sigma_measured ,
			double& k , double& theta ) 
{
  double mu2 = mu_measured*mu_measured;
  double sigma2 = sigma_measured*sigma_measured;
  // from Wikipedia's http://en.wikipedia.org/wiki/Gamma_distribution parameterization.
  theta = sigma2/mu_measured;
  k = mu2/sigma2+1;

  cout << "For mu    : " << mu_measured << endl;
  cout << "For sigma : " << sigma_measured << endl;
  cout << "k is      : " << k << endl;
  cout << "theta is  : " << theta << endl;
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

  //cout << "The calculated value of alpha is : " << alpha << endl;
  //cout << "The calculated value of beta is  : " << beta << endl;
  //cout << "The input value of mu is         : " << mu_measured << endl;
  //cout << "The input value of sigma is      : " << sigma_measured << endl;
  //cout << "The calculated value of mu is    : " << (alpha - 1) / ( alpha + beta - 2 ) << endl;
  //cout << "The calculated value of sigma is : " << sqrt( alpha * beta / ( pow(alpha + beta,2) * (alpha + beta + 1) ) ) << endl;

  //beta = exp(sigma_measured);
  //alpha = (mu_measured);
  //
  //cout << "The input value of mu is         : " << mu_measured << endl;
  //cout << "The input value of sigma is      : " << sigma_measured << endl;
  //cout << "The input value of sigma/mu is   : " << sigma_measured/mu_measured << endl;
  //
  //cout << "The calculated value of alpha is : " << log(alpha) << endl;
  //cout << "The calculated value of beta is  : " << log(beta) << endl;

}

void betaMeanTransform(const double& mu_measured , const double& sigma_measured ,
		   double& alpha , double& beta ) 
{
  double mu2 = mu_measured*mu_measured;
  double mu3 = mu_measured*mu2;
  double sigma2 = sigma_measured*sigma_measured;

  alpha = (mu2-mu3-mu_measured*sigma2)/sigma2;
  beta = (mu_measured-1.)*(mu2 - mu_measured + sigma2)/sigma2;

  if( alpha<1 && beta<1 )
    {
      cout << "mean and variance un-physical with beta pdf, setting uniform distribution" << endl;
      alpha = 1;
      beta = 1;
    }

}

// to be used for binding

//double twoSidedGaussian(double x, double mu, double sigmaL, double sigmaR)
//{
//  if( x < 0 )
//    {
//      return x * sigmaL + mu;
//    }
//  else if( x > 0 )
//    {
//      return x * sigmaR + mu;
//    }
//  return mu;
//}
//
//double brokenLogNormal(double x, double mu, double sigmaL, double sigmaR)
//{
//  cout <<"using broken log-normal function"<<endl;
//  if( x < 0 ) {
//    double thisSigma = log(mu) - log(mu-sigmaL);
//    if( thisSigma < 0 ) x*=-1;
//    cout << "returning " << exp( thisSigma*x + log(mu) ) << endl;
//    return exp( thisSigma*x + log(mu) );
//  }
//  else {
//    double thisSigma = log(mu+sigmaR) - log(mu);
//    if( thisSigma < 0 ) x*=-1;
//    cout << "returning " << exp( thisSigma*x + log(mu) ) << endl;
//    return exp( thisSigma*x + log(mu) );
//  }
//}

  //=====================================================================================================


   ra2bRoostatsClass7::ra2bRoostatsClass7( bool ArgUseSigTtwjVar, bool ArgUseLdpVars ) {

      gStyle->SetOptStat(0) ;

      useSigTtwjVar = ArgUseSigTtwjVar ;
      useLdpVars = ArgUseLdpVars ;


     //--- Tell RooFit to shut up about anything less important than an ERROR.
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
  /// RooMsgService::instance().addStream(DEBUG,Topic(Tracing),OutputFile("debug.log")) ;
      printf("\n\n ==== RooFit output configuration ===============\n") ;
      RooMsgService::instance().Print("v") ;
      printf("\n\n ================================================\n") ;

      varsAtFitVals = false ;
      initialized = false ;




   }





  //===================================================================

   ra2bRoostatsClass7::~ra2bRoostatsClass7() { }



  //===================================================================

    bool ra2bRoostatsClass7::initialize( const char* infile ,
                                         const char* inputScanFile, double m0, double m12, bool isT1bbbb, double t1bbbbXsec,
					 const char* outfile, bool effInitFixed, bool constantNonPoisson , bool constantSigEff ) {


       printf( "\n\n Opening input file : %s\n\n", infile ) ;

       FILE* infp ;
       if ( (infp=fopen( infile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
          return false ;
       }

       sprintf( initializeFile, "%s", infile ) ;


       int    Nsig                  ;
       int    Nsb                   ;
       int    Nsig_sl               ;
       int    Nsb_sl                ;
       int    Nsig_ldp              ;
       int    Nsb_ldp               ;
   /// int    Nhtonlytrig_lsb_0b      ;
   /// int    Nhtonlytrig_lsb_0b_ldp  ;
       int    Nsb_ee                  ;
       int    Nsig_ee                 ;
       int    Nsb_mm                  ;
       int    Nsig_mm                 ;
       float  Nttbarmc_sig_ldp      ;
       float  Nttbarmc_sb_ldp       ;
       float  NWJmc_sig_ldp         ;
       float  NWJmc_sb_ldp          ;
       float  NZnnmc_sig_ldp        ;
       float  NZnnmc_sb_ldp         ;
       float  NEwomc_sig_ldp        ;
       float  NEwomc_sb_ldp         ;
       float  NEwomc_sig  ;
       float  NEwomc_sb   ;


       float  sf_mc            ;
       float  sf_mc_err        ;
       float  sf_qcd_sb        ;
       float  sf_qcd_sb_err    ;
       float  sf_qcd_sig       ;
       float  sf_qcd_sig_err   ;
       float  sf_ttwj_sig      ;
       float  sf_ttwj_sig_err  ;
       float  sf_ee            ;
       float  sf_ee_err        ;
       float  sf_mm            ;
       float  sf_mm_err        ;

       float  Rlsb_passfail     ;
       float  Rlsb_passfail_err ;

       // trigger efficiency corrections
       float eps_sb_mean ;
       float eps_sig_mean ;
       float eps_sb_ldp_mean ;
       float eps_sb_sl_mean ;
       float eps_sig_sl_mean ;

       // trigger efficiency scale factors (with asymmetric errors)
       float epsSF_sb ;
       float epsSF_sb_errp ;
       float epsSF_sb_errm ;
       float epsSF_sig ;
       float epsSF_sig_errp ;
       float epsSF_sig_errm ;
       float epsSF_sb_ldp ;
       float epsSF_sb_ldp_errp ;
       float epsSF_sb_ldp_errm ;
       float epsSF_sb_sl ;
       float epsSF_sb_sl_errp ;
       float epsSF_sb_sl_errm ;
       float epsSF_sig_sl ;
       float epsSF_sig_sl_errp ;
       float epsSF_sig_sl_errm ;

       //Kristen's Measurement for T1bbbb
       float kristen_ttwj_sig_mean;
       float kristen_ttwj_sig_stat;
       float kristen_ttwj_sig_syst;

       //--- read in description line.
       printf("\n\n") ;
       char c(0) ;
       while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
       printf("\n\n") ;


       char label[1000] ;

      //-----    The order here matters!

       fscanf( infp, "%s %d", label, &Nsig                  ) ;   printf( "%s %d\n", label, Nsig                  ) ;
       fscanf( infp, "%s %d", label, &Nsb                   ) ;   printf( "%s %d\n", label, Nsb                   ) ;
       fscanf( infp, "%s %d", label, &Nsig_sl               ) ;   printf( "%s %d\n", label, Nsig_sl               ) ;
       fscanf( infp, "%s %d", label, &Nsb_sl                ) ;   printf( "%s %d\n", label, Nsb_sl                ) ;
       fscanf( infp, "%s %d", label, &Nsig_ldp              ) ;   printf( "%s %d\n", label, Nsig_ldp              ) ;
       fscanf( infp, "%s %d", label, &Nsb_ldp               ) ;   printf( "%s %d\n", label, Nsb_ldp               ) ;
       fscanf( infp, "%s %d", label, &Nsb_ee                    ) ;   printf( "%s %d\n", label, Nsb_ee                    ) ;
       fscanf( infp, "%s %d", label, &Nsig_ee                   ) ;   printf( "%s %d\n", label, Nsig_ee                   ) ;
       fscanf( infp, "%s %d", label, &Nsb_mm                    ) ;   printf( "%s %d\n", label, Nsb_mm                    ) ;
       fscanf( infp, "%s %d", label, &Nsig_mm                   ) ;   printf( "%s %d\n", label, Nsig_mm                   ) ;
       fscanf( infp, "%s %g", label, &Rlsb_passfail         ) ;   printf( "%s %g\n", label, Rlsb_passfail         ) ;
       fscanf( infp, "%s %g", label, &Rlsb_passfail_err     ) ;   printf( "%s %g\n", label, Rlsb_passfail_err     ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sig_ldp      ) ;   printf( "%s %g\n", label, Nttbarmc_sig_ldp      ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sb_ldp       ) ;   printf( "%s %g\n", label, Nttbarmc_sb_ldp       ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig_ldp         ) ;   printf( "%s %g\n", label, NWJmc_sig_ldp         ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb_ldp          ) ;   printf( "%s %g\n", label, NWJmc_sb_ldp          ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig_ldp        ) ;   printf( "%s %g\n", label, NZnnmc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb_ldp         ) ;   printf( "%s %g\n", label, NZnnmc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig_ldp        ) ;   printf( "%s %g\n", label, NEwomc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb_ldp         ) ;   printf( "%s %g\n", label, NEwomc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig        ) ;   printf( "%s %g\n", label, NEwomc_sig        ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb         ) ;   printf( "%s %g\n", label, NEwomc_sb         ) ;
       fscanf( infp, "%s %g", label, &DataLumi                  ) ;   printf( "%s %g\n", label, DataLumi                  ) ;
       fscanf( infp, "%s %g", label, &acc_ee_sig_mean              ) ;   printf( "%s %g\n", label, acc_ee_sig_mean               ) ;
       fscanf( infp, "%s %g", label, &acc_ee_sig_err               ) ;   printf( "%s %g\n", label, acc_ee_sig_err                ) ;
       fscanf( infp, "%s %g", label, &acc_ee_sb_mean               ) ;   printf( "%s %g\n", label, acc_ee_sb_mean                ) ;
       fscanf( infp, "%s %g", label, &acc_ee_sb_err                ) ;   printf( "%s %g\n", label, acc_ee_sb_err                 ) ;
       fscanf( infp, "%s %g", label, &acc_mm_sig_mean              ) ;   printf( "%s %g\n", label, acc_mm_sig_mean               ) ;
       fscanf( infp, "%s %g", label, &acc_mm_sig_err               ) ;   printf( "%s %g\n", label, acc_mm_sig_err                ) ;
       fscanf( infp, "%s %g", label, &acc_mm_sb_mean               ) ;   printf( "%s %g\n", label, acc_mm_sb_mean                ) ;
       fscanf( infp, "%s %g", label, &acc_mm_sb_err                ) ;   printf( "%s %g\n", label, acc_mm_sb_err                 ) ;
       fscanf( infp, "%s %g", label, &eff_ee_mean               ) ;   printf( "%s %g\n", label, eff_ee_mean               ) ;
       fscanf( infp, "%s %g", label, &eff_ee_err                ) ;   printf( "%s %g\n", label, eff_ee_err                ) ;
       fscanf( infp, "%s %g", label, &eff_mm_mean               ) ;   printf( "%s %g\n", label, eff_mm_mean               ) ;
       fscanf( infp, "%s %g", label, &eff_mm_err                ) ;   printf( "%s %g\n", label, eff_mm_err                ) ;
       fscanf( infp, "%s %g", label, &Ztoll_lumi                ) ;   printf( "%s %g\n", label, Ztoll_lumi                ) ;
       fscanf( infp, "%s %g", label, &knn_ee_sig_mean              ) ;   printf( "%s %g\n", label, knn_ee_sig_mean              ) ;
       fscanf( infp, "%s %g", label, &knn_ee_sig_err               ) ;   printf( "%s %g\n", label, knn_ee_sig_err               ) ;
       fscanf( infp, "%s %g", label, &knn_ee_sb_mean               ) ;   printf( "%s %g\n", label, knn_ee_sb_mean               ) ;
       fscanf( infp, "%s %g", label, &knn_ee_sb_err                ) ;   printf( "%s %g\n", label, knn_ee_sb_err                ) ;
       fscanf( infp, "%s %g", label, &knn_mm_sig_mean              ) ;   printf( "%s %g\n", label, knn_mm_sig_mean              ) ;
       fscanf( infp, "%s %g", label, &knn_mm_sig_err               ) ;   printf( "%s %g\n", label, knn_mm_sig_err               ) ;
       fscanf( infp, "%s %g", label, &knn_mm_sb_mean               ) ;   printf( "%s %g\n", label, knn_mm_sb_mean               ) ;
       fscanf( infp, "%s %g", label, &knn_mm_sb_err                ) ;   printf( "%s %g\n", label, knn_mm_sb_err                ) ;
       fscanf( infp, "%s %g", label, &fsig_ee_mean              ) ;   printf( "fsig_ee_mean    %s %g\n", label, fsig_ee_mean              ) ;
       fscanf( infp, "%s %g", label, &fsig_ee_err               ) ;   printf( "fsig_ee_err     %s %g\n", label, fsig_ee_err               ) ;
       fscanf( infp, "%s %g", label, &fsig_mm_mean              ) ;   printf( "fsig_mm_mean    %s %g\n", label, fsig_mm_mean              ) ;
       fscanf( infp, "%s %g", label, &fsig_mm_err               ) ;   printf( "fsig_mm_err     %s %g\n", label, fsig_mm_err               ) ;
       fscanf( infp, "%s %g", label, &sf_mc                     ) ;   printf( "sf_mc           %s %g\n", label, sf_mc                     ) ;
       fscanf( infp, "%s %g", label, &sf_mc_err                 ) ;   printf( "sf_mc_err       %s %g\n", label, sf_mc_err                 ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb                 ) ;   printf( "sf_qcd_sb       %s %g\n", label, sf_qcd_sb                 ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb_err             ) ;   printf( "sf_qcd_sb_err   %s %g\n", label, sf_qcd_sb_err             ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig                ) ;   printf( "sf_qcd_sig      %s %g\n", label, sf_qcd_sig                ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig_err            ) ;   printf( "sf_qcd_sig_err  %s %g\n", label, sf_qcd_sig_err            ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig               ) ;   printf( "sf_ttwj_sig     %s %g\n", label, sf_ttwj_sig               ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig_err           ) ;   printf( "sf_ttwj_sig_err %s %g\n", label, sf_ttwj_sig_err           ) ;
       fscanf( infp, "%s %g", label, &sf_ee                     ) ;   printf( "sf_ee           %s %g\n", label, sf_ee                     ) ;
       fscanf( infp, "%s %g", label, &sf_ee_err                 ) ;   printf( "sf_ee_err       %s %g\n", label, sf_ee_err                 ) ;
       fscanf( infp, "%s %g", label, &sf_mm                     ) ;   printf( "sf_mm           %s %g\n", label, sf_mm                     ) ;
       fscanf( infp, "%s %g", label, &sf_mm_err                 ) ;   printf( "sf_mm_err       %s %g\n", label, sf_mm_err                 ) ;
       fscanf( infp, "%s %g", label, &eps_sb_mean               ) ;   printf( "eps_sb_mean     %s %g\n", label, eps_sb_mean               ) ;
       fscanf( infp, "%s %g", label, &eps_sig_mean              ) ;   printf( "eps_sig_mean    %s %g\n", label, eps_sig_mean              ) ;
       fscanf( infp, "%s %g", label, &eps_sb_ldp_mean           ) ;   printf( "eps_sb_ldp_mean %s %g\n", label, eps_sb_ldp_mean           ) ;
       fscanf( infp, "%s %g", label, &eps_sb_sl_mean            ) ;   printf( "eps_sb_sl_mean  %s %g\n", label, eps_sb_sl_mean            ) ;
       fscanf( infp, "%s %g", label, &eps_sig_sl_mean           ) ;   printf( "eps_sig_sl_mean %s %g\n", label, eps_sig_sl_mean           ) ;
       fscanf( infp, "%s %g", label, &epsSF_sb                  ) ;   printf( "epsSF_sb          %s %g\n", label, epsSF_sb                ) ;
       fscanf( infp, "%s %g", label, &epsSF_sb_errp             ) ;   printf( "epsSF_sb_errp     %s %g\n", label, epsSF_sb_errp           ) ;
       fscanf( infp, "%s %g", label, &epsSF_sb_errm             ) ;   printf( "epsSF_sb_errm     %s %g\n", label, epsSF_sb_errm           ) ;
       fscanf( infp, "%s %g", label, &epsSF_sig                 ) ;   printf( "epsSF_sig         %s %g\n", label, epsSF_sig               ) ;
       fscanf( infp, "%s %g", label, &epsSF_sig_errp            ) ;   printf( "epsSF_sig_errp    %s %g\n", label, epsSF_sig_errp          ) ;
       fscanf( infp, "%s %g", label, &epsSF_sig_errm            ) ;   printf( "epsSF_sig_errm    %s %g\n", label, epsSF_sig_errm          ) ;
       fscanf( infp, "%s %g", label, &epsSF_sb_ldp              ) ;   printf( "epsSF_sb_ldp      %s %g\n", label, epsSF_sb_ldp            ) ;
       fscanf( infp, "%s %g", label, &epsSF_sb_ldp_errp         ) ;   printf( "epsSF_sb_ldp_errp %s %g\n", label, epsSF_sb_ldp_errp       ) ;
       fscanf( infp, "%s %g", label, &epsSF_sb_ldp_errm         ) ;   printf( "epsSF_sb_ldp_errm %s %g\n", label, epsSF_sb_ldp_errm       ) ;
       fscanf( infp, "%s %g", label, &epsSF_sb_sl               ) ;   printf( "epsSF_sb_sl       %s %g\n", label, epsSF_sb_sl             ) ;
       fscanf( infp, "%s %g", label, &epsSF_sb_sl_errp          ) ;   printf( "epsSF_sb_sl_errp  %s %g\n", label, epsSF_sb_sl_errp        ) ;
       fscanf( infp, "%s %g", label, &epsSF_sb_sl_errm          ) ;   printf( "epsSF_sb_sl_errm  %s %g\n", label, epsSF_sb_sl_errm        ) ;
       fscanf( infp, "%s %g", label, &epsSF_sig_sl              ) ;   printf( "epsSF_sig_sl      %s %g\n", label, epsSF_sig_sl            ) ;
       fscanf( infp, "%s %g", label, &epsSF_sig_sl_errp         ) ;   printf( "epsSF_sig_sl_errp %s %g\n", label, epsSF_sig_sl_errp       ) ;
       fscanf( infp, "%s %g", label, &epsSF_sig_sl_errm         ) ;   printf( "epsSF_sig_sl_errm %s %g\n", label, epsSF_sig_sl_errm       ) ;
       fscanf( infp, "%s %g", label, &kristen_ttwj_sig_mean    ) ;   printf( "kristen_ttwj_sig_mean %s %g\n", label, kristen_ttwj_sig_mean       ) ;
       fscanf( infp, "%s %g", label, &kristen_ttwj_sig_stat     ) ;   printf( "kristen_ttwj_sig_stat %s %g\n", label, kristen_ttwj_sig_stat       ) ;
       fscanf( infp, "%s %g", label, &kristen_ttwj_sig_syst     ) ;   printf( "kristen_ttwj_sig_syst %s %g\n", label, kristen_ttwj_sig_syst       ) ;

       printf("\n Done reading in %s\n\n", infile ) ;
       fclose( infp ) ;

    //---- calculations for determining initial values for floating parameters.

     //-- Znunu stuff

       float initialval_znn_sig_ee(2.) ;
       float initialval_znn_sig_mm(2.) ;
       float initialval_znn_sig(2.) ;
       float initialval_znn_sb_ee(2.) ;
       float initialval_znn_sb_mm(2.) ;
       float initialval_znn_sb(2.) ;


       initialval_znn_sig_ee = (Nsig_ee) * ( 5.95 * DataLumi * fsig_ee_mean * knn_ee_sig_mean ) / ( acc_ee_sig_mean * eff_ee_mean * Ztoll_lumi ) ;
       initialval_znn_sb_ee  = (Nsb_ee ) * ( 5.95 * DataLumi * fsig_ee_mean * knn_ee_sb_mean  ) / ( acc_ee_sb_mean  * eff_ee_mean * Ztoll_lumi ) ;

       initialval_znn_sig_mm = (Nsig_mm) * ( 5.95 * DataLumi * fsig_mm_mean * knn_mm_sig_mean ) / ( acc_mm_sig_mean * eff_mm_mean * Ztoll_lumi ) ;
       initialval_znn_sb_mm  = (Nsb_mm ) * ( 5.95 * DataLumi * fsig_mm_mean * knn_mm_sb_mean  ) / ( acc_mm_sb_mean  * eff_mm_mean * Ztoll_lumi ) ;


       //-- really dumb ave.
       initialval_znn_sig = 0.5 * ( initialval_znn_sig_ee + initialval_znn_sig_mm ) ;
       initialval_znn_sb  = 0.5 * ( initialval_znn_sb_ee  + initialval_znn_sb_mm ) ;


  /////--- QCD lsb stuff

  ///  double initialval_qcd_lsb_0b     = Nhtonlytrig_lsb_0b     ;
  ///  double initialval_qcd_lsb_0b_ldp = Nhtonlytrig_lsb_0b_ldp ;


     //--- QCD ldp stuff

       double initialval_qcd_sig_ldp = ( Nsig_ldp / eps_sig_mean ) - ( Nttbarmc_sig_ldp + NWJmc_sig_ldp + NZnnmc_sig_ldp + NEwomc_sig_ldp ) ;
       double initialval_qcd_sb_ldp  = ( Nsb_ldp  / eps_sb_ldp_mean ) - ( Nttbarmc_sb_ldp  + NWJmc_sb_ldp  + NZnnmc_sb_ldp  + NEwomc_sb_ldp  ) ;


     //--- QCD sig and sb stuff

   /// if ( Nhtonlytrig_lsb_0b_ldp <= 0 ) {
   ///    printf("\n\n\n *** Nhtonlytrig_lsb_0b_ldp has a crazy value (%d).  I quit.\n\n\n", Nhtonlytrig_lsb_0b_ldp ) ;
   ///    return false ;
   /// }
   /// double Rldp = initialval_qcd_lsb_0b / initialval_qcd_lsb_0b_ldp ;
       double initialval_qcd_sig = Rlsb_passfail * initialval_qcd_sig_ldp ;
       double initialval_qcd_sb  = Rlsb_passfail * initialval_qcd_sb_ldp  ;


     //--- ttwj SL

       double initialval_ttwj_sig_sl = Nsig_sl / eps_sig_sl_mean ;
       double initialval_ttwj_sb_sl  = Nsb_sl / eps_sb_sl_mean ;


     //--- ttwj sig and sb

       if ( initialval_ttwj_sb_sl <= 0 ) {
          printf("\n\n\n *** initialval_ttwj_sb_sl has a crazy value (%.1f).  I quit.\n\n\n", initialval_ttwj_sb_sl ) ;
          return false ;
       }
       double initialval_ttwj_sb     = ( Nsb / eps_sb_mean ) - ( initialval_qcd_sb + initialval_znn_sb ) ;
       double initialval_ttwj_sig    = kristen_ttwj_sig_mean ;



       printf("\n\n\n --------- Observables and floating parameter initial values. ------------\n\n") ;

       printf("          |  Nobs   ||  ttwj  |  QCD  |  Znn  |\n") ;
       printf(" SIG      | %5d   || %5.1f | %5.1f | %5.1f |\n", Nsig, initialval_ttwj_sig, initialval_qcd_sig, initialval_znn_sig ) ;
       printf(" SB       | %5d   || %5.1f | %5.1f | %5.1f |\n", Nsb , initialval_ttwj_sb , initialval_qcd_sb , initialval_znn_sb  ) ;
       printf(" SIG,LDP  | %5d   ||*%5.1f | %5.1f |*%5.1f |\n", Nsig_ldp, (Nttbarmc_sig_ldp+NWJmc_sig_ldp+NEwomc_sig_ldp), initialval_qcd_sig_ldp, NZnnmc_sig_ldp ) ;
       printf(" SB ,LDP  | %5d   ||*%5.1f | %5.1f |*%5.1f |\n", Nsb_ldp , (Nttbarmc_sb_ldp +NWJmc_sb_ldp +NEwomc_sb_ldp) , initialval_qcd_sb_ldp , NZnnmc_sb_ldp  ) ;
       printf("\n * means fixed MC value.\n\n\n") ;




     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


       printf(" --- Defining observables.\n" ) ;


      rv_Nsig        = new RooRealVar( "Nsig"        , "Nsig"        , 0.0, 1000000. ) ;
      rv_Nsb         = new RooRealVar( "Nsb"         , "Nsb"         , 0.0, 1000000. ) ;

 ///  rv_Nsig_sl     = new RooRealVar( "Nsig_sl"     , "Nsig_sl"     , 0.0, 1000000. ) ;
 ///  rv_Nsb_sl      = new RooRealVar( "Nsb_sl"      , "Nsb_sl"      , 0.0, 1000000. ) ;

      rv_Nsig_ldp    = new RooRealVar( "Nsig_ldp"    , "Nsig_ldp"    , 0.0, 1000000. ) ;
      rv_Nsb_ldp     = new RooRealVar( "Nsb_ldp"     , "Nsb_ldp"     , 0.0, 1000000. ) ;

 ///  rv_Nlsb_0b     = new RooRealVar( "Nlsb_0b"     , "Nlsb_0b"     , 0.0, 1000000. ) ;
 ///  rv_Nlsb_0b_ldp = new RooRealVar( "Nlsb_0b_ldp" , "Nlsb_0b_ldp" , 0.0, 1000000. ) ;


      rv_Nsb_ee      = new RooRealVar( "Nsb_ee"      ,"Nsb_ee"       , 0., 10000. ) ;
      rv_Nsig_ee     = new RooRealVar( "Nsig_ee"     ,"Nsig_ee"      , 0., 10000. ) ;

      rv_Nsb_mm      = new RooRealVar( "Nsb_mm"      ,"Nsb_mm"       , 0., 10000. ) ;
      rv_Nsig_mm     = new RooRealVar( "Nsig_mm"     ,"Nsig_mm"      , 0., 10000. ) ;






      rv_Nsig        -> setVal( Nsig ) ;
      rv_Nsb         -> setVal( Nsb ) ;

 ///  rv_Nsig_sl     -> setVal( Nsig_sl ) ;
 ///  rv_Nsb_sl      -> setVal( Nsb_sl ) ;

      rv_Nsig_ldp    -> setVal( Nsig_ldp ) ;
      rv_Nsb_ldp     -> setVal( Nsb_ldp ) ;

  /// rv_Nlsb_0b     -> setVal( Nhtonlytrig_lsb_0b ) ;
  /// rv_Nlsb_0b_ldp -> setVal( Nhtonlytrig_lsb_0b_ldp ) ;


      rv_Nsb_ee      -> setVal( Nsb_ee ) ;
      rv_Nsig_ee     -> setVal( Nsig_ee ) ;

      rv_Nsb_mm      -> setVal( Nsb_mm ) ;
      rv_Nsig_mm     -> setVal( Nsig_mm ) ;








    //++++++++ Parameters of the likelihood +++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining parameters.\n" ) ;

    //--- Systematics and other nuisance parameters
      RooArgSet globalObservables ("globalObservables");
      RooArgSet allNuisances ("allNuisances");
      RooArgSet allPoissonNuisances ("allPoissonNuisances");
      RooArgSet allNonPoissonNuisances ("allNonPoissonNuisances");
      RooArgSet allNuisancePdfs ("allNuisancePdfs");



    //____ Counts in SIG ______________________

      //----------------Kristen's Constraints---------------------//

      double k,theta;

      //Want to take the statistical error that is not sqrt(N) and put it in the systematic
      //So, want to have new_syst = sqrt(old_syst^2 + old_stat^2 - mean)
      //Then mean +/- sqrt(mean) +/- new_syst is ~right.  Nothing else is correlated with this so should be ok.
      //double kristen_ttwj_sig_new_syst = kristen_ttwj_sig_syst*kristen_ttwj_sig_syst + kristen_ttwj_sig_stat*kristen_ttwj_sig_stat - kristen_ttwj_sig_mean;
      double kristen_ttwj_sig_new_syst = kristen_ttwj_sig_syst*kristen_ttwj_sig_syst + kristen_ttwj_sig_stat*kristen_ttwj_sig_stat;
      kristen_ttwj_sig_new_syst = (kristen_ttwj_sig_new_syst > 0) ? sqrt(kristen_ttwj_sig_new_syst) : 0;

      gammaModeTransform(kristen_ttwj_sig_mean,kristen_ttwj_sig_new_syst,k,theta);
      RooRealVar kristen_mu_ttwj_sig ("mu_ttwj_sig","mu_ttwj_sig",kristen_ttwj_sig_mean,0,1e5);
      RooRealVar mu_ttwj_sig_k ("mu_ttwj_sig_k", "mu_ttwj_sig_k", k);
      RooRealVar mu_ttwj_sig_theta ("mu_ttwj_sig_theta", "mu_ttwj_sig_theta", theta);
      RooGammaPdf pdf_mu_ttwj_sig ("pdf_mu_ttwj_sig" , "pdf_mu_ttwj_sig", kristen_mu_ttwj_sig, mu_ttwj_sig_k, mu_ttwj_sig_theta, RooConst(0));
      mu_ttwj_sig_k.setConstant();
      mu_ttwj_sig_theta.setConstant();
      allNonPoissonNuisances.add (kristen_mu_ttwj_sig);
      if(constantNonPoisson) kristen_mu_ttwj_sig.setConstant();
      allNuisancePdfs.add (pdf_mu_ttwj_sig);

      //cout << "The mean is " << kristen_ttwj_sig_mean << endl;
      //cout << "The stat is " << kristen_ttwj_sig_stat << endl;
      //cout << "The syst is " << kristen_ttwj_sig << endl;
      //cout << "total error^2 is " << kristen_ttwj_sig*kristen_ttwj_sig + kristen_ttwj_sig_stat*kristen_ttwj_sig_stat << endl;
      //cout << "the new syst is: " << kristen_ttwj_sig_new << endl;
      //if(kristen_ttwj_sig_new_syst > 0)
      //	  {
      //	    allNuisances.add (kristen_mu_ttwj_sig);
      //	    allNuisancePdfs.add (pdf_mu_ttwj_sig);
      //	  }

      if ( useSigTtwjVar ) {
	

	if(kristen_ttwj_sig_new_syst > 0)
	  {
	    rv_mu_ttwj_sig = &kristen_mu_ttwj_sig;
	  }
	else{
	  //cout << "using RooRealVar" << endl;
	  rv_mu_ttwj_sig = &kristen_mu_ttwj_sig;
	  ((RooRealVar*)rv_mu_ttwj_sig)->setConstant();
	}
	
      }
      if ( !useLdpVars ) {
	cout <<"This option is invalid, please set useLdpVars to true)"<<endl;
	return false;
      }


      //-- Note: Ewo is rfv

      rv_mu_znn_sig      = new RooRealVar( "mu_znn_sig"    , "mu_znn_sig"    , 0.0, 300. ) ;
      allPoissonNuisances.add(*rv_mu_znn_sig);

      float maxSusySig = 4.0*Nsig ;
      rv_mu_susy_sig     = new RooRealVar( "mu_susy_sig"   , "mu_susy_sig"   , 0.0, maxSusySig ) ;


      rv_mu_znn_sig   -> setVal( initialval_znn_sig ) ;  //-- this is a starting value only.
      rv_mu_susy_sig    -> setVal( 0. ) ;  //-- this is a starting value only.




    //____ Counts in SB  ______________________


      if ( !useSigTtwjVar ) {
	cout << "This option is not valid with this likelihood.  Please set useSigTtwjVar to true." << endl;
	return false;
      }
      if ( !useLdpVars ) {
	cout <<"This option is invalid, please set useLdpVars to true"<<endl;
	return false;
      }

      //-- Note: QCD is rfv
      //-- Note: Ewo is rfv
      //-- Note: SUSY is rfv


      rrv_mu_znn_sb       = new RooRealVar( "mu_znn_sb"     , "mu_znn_sb"     , 0.0, 350. ) ;

      rrv_mu_znn_sb   -> setVal( initialval_znn_sb ) ;  //-- this is a starting value only.

      rv_mu_znn_sb = rrv_mu_znn_sb ;
      allPoissonNuisances.add(*rv_mu_znn_sb);
      //-- Note: Znn is rfv in Znn model 2.






    //____ Counts in SIG, SL  ______________________

      // rv_mu_ttwj_sig_sl  = new RooRealVar( "mu_ttwj_sig_sl"    , "mu_ttwj_sig_sl"    , 0.0, 2500. ) ;

      //-- Note: QCD, Ewo, and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      // rv_mu_ttwj_sig_sl  -> setVal( initialval_ttwj_sig_sl ) ;  //-- this is a starting value only.







    //____ Counts in SB, SL  ______________________

      // rv_mu_ttwj_sb_sl  = new RooRealVar( "mu_ttwj_sb_sl"    , "mu_ttwj_sb_sl"    , 0.0, 3000. ) ;

      //-- Note: QCD, Ewo, and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      // rv_mu_ttwj_sb_sl  -> setVal( initialval_ttwj_sb_sl ) ;  //-- this is a starting value only.







    //____ Counts in SIG, LDP  ______________________

      cout << " --- SIG, LDP" << endl;

      if ( useLdpVars ) {
         rrv_mu_qcd_sig_ldp  = new RooRealVar( "mu_qcd_sig_ldp"    , "mu_qcd_sig_ldp"    , 0.0, 3500. ) ;
         rv_mu_qcd_sig_ldp = rrv_mu_qcd_sig_ldp ;
	 allPoissonNuisances.add(*rv_mu_qcd_sig_ldp);
         rrv_mu_qcd_sig_ldp  -> setVal( initialval_qcd_sig_ldp ) ; //-- this is a starting value only.
      }

      //-- Note: Ewo is assumed to be negligible and is not explicitly included.
      //-- Note: Znn is rfv (MC)
      //-- Note: ttwj is rfv (MC)
      //-- Note: SUSY is rfv




    //____ Counts in SB, LDP  ______________________

      cout << " --- SB, LDP" << endl;

      if ( useLdpVars ) {
         rrv_mu_qcd_sb_ldp  = new RooRealVar( "mu_qcd_sb_ldp"    , "mu_qcd_sb_ldp"    , 0.0, 3000. ) ;
         rv_mu_qcd_sb_ldp = rrv_mu_qcd_sb_ldp ;
	 allPoissonNuisances.add(*rv_mu_qcd_sb_ldp);
         rrv_mu_qcd_sb_ldp  -> setVal( initialval_qcd_sb_ldp ) ; //-- this is a starting value only.
      }

      //-- Note: Ewo is assumed to be negligible and is not explicitly included.
      //-- Note: Znn is rfv (MC)
      //-- Note: ttwj is rfv (MC)
      //-- Note: SUSY is rfv








 // //____ Counts in LSB, 0b  ______________________

 //   rv_mu_qcd_lsb_0b  = new RooRealVar( "mu_qcd_lsb_0b"    ,  "mu_qcd_lsb_0b" ,  0.    ,  10000. ) ;

 //   //-- Note: The 0btag LSB is assumed to be 100% QCD.

 //   rv_mu_qcd_lsb_0b  -> setVal( initialval_qcd_lsb_0b ) ;  //-- this is a starting value only.






 // //____ Counts in LSB, 0b, LDP  ______________________

 //   rv_mu_qcd_lsb_0b_ldp  = new RooRealVar( "mu_qcd_lsb_0b_ldp"    ,  "mu_qcd_lsb_0b_ldp",   0. ,  10000. ) ;

 //   //-- Note: The 0btag LSB, LDP is assumed to be 100% QCD.

 //   rv_mu_qcd_lsb_0b_ldp  -> setVal( initialval_qcd_lsb_0b_ldp ) ;  //-- this is a starting value only.






    //____ Counts in SB, ee  ______________________


      //-- Note: The Z to ee sample is assumed to be 100% Z to ee.
      //-- Note: zee is rfv






    //____ Counts in SIG, ee  ______________________


      //-- Note: The Z to ee sample is assumed to be 100% Z to ee.
      //-- Note: zee is rfv








    //____ Counts in SB, mm  ______________________


      //-- Note: The Z to mm sample is assumed to be 100% Z to mm.
      //-- Note: zmm is rfv






    //____ Counts in SIG, mm  ______________________


      //-- Note: The Z to mm sample is assumed to be 100% Z to mm.
      //-- Note: zmm is rfv






    //____ MC inputs _______________________________

      printf(" --- Defining MC parameters.\n" ) ;

     //-- SUSY

      rv_mu_susymc_sig      = new RooRealVar( "mu_susymc_sig"     , "mu_susymc_sig"     , 0.0, 100000. ) ;
      ///rv_mu_susymc_sb       = new RooRealVar( "mu_susymc_sb"      , "mu_susymc_sb"      , 0.0, 100000. ) ;
      ///rv_mu_susymc_sig_sl   = new RooRealVar( "mu_susymc_sig_sl"  , "mu_susymc_sig_sl"  , 0.0, 100000. ) ;
      ///rv_mu_susymc_sb_sl    = new RooRealVar( "mu_susymc_sb_sl"   , "mu_susymc_sb_sl"   , 0.0, 100000. ) ;
      rv_mu_susymc_sig_ldp  = new RooRealVar( "mu_susymc_sig_ldp" , "mu_susymc_sig_ldp" , 0.0, 100000. ) ;
      rv_mu_susymc_sb_ldp   = new RooRealVar( "mu_susymc_sb_ldp"  , "mu_susymc_sb_ldp"  , 0.0, 100000. ) ;

      rv_mu_susymc_sig     -> setVal( 0.1 ) ;
      ///rv_mu_susymc_sb      -> setVal( 0. ) ;
      ///rv_mu_susymc_sig_sl  -> setVal( 0. ) ;
      ///rv_mu_susymc_sb_sl   -> setVal( 0. ) ;
      rv_mu_susymc_sig_ldp -> setVal( 0. ) ;
      rv_mu_susymc_sb_ldp  -> setVal( 0. ) ;

      rv_mu_susymc_sig     -> setConstant(kTRUE) ;
      ///rv_mu_susymc_sb      -> setConstant(kTRUE) ;
      ///rv_mu_susymc_sig_sl  -> setConstant(kTRUE) ;
      ///rv_mu_susymc_sb_sl   -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_ldp -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_ldp  -> setConstant(kTRUE) ;



     //-- SIG, LDP

      rv_mu_ttbarmc_sig_ldp   = new RooRealVar( "mu_ttbarmc_sig_ldp" ,"mu_ttbarmc_sig_ldp" , 0., 1000. ) ;
      rv_mu_WJmc_sig_ldp      = new RooRealVar( "mu_WJmc_sig_ldp"    ,"mu_WJmc_sig_ldp"    , 0., 1000. ) ;
      rv_mu_Znnmc_sig_ldp     = new RooRealVar( "mu_Znnmc_sig_ldp"   ,"mu_Znnmc_sig_ldp"   , 0., 1000. ) ;
      rv_mu_Ewomc_sig_ldp     = new RooRealVar( "mu_Ewomc_sig_ldp"   ,"mu_Ewomc_sig_ldp"   , 0., 1000. ) ;

      rv_mu_ttbarmc_sig_ldp  -> setVal( Nttbarmc_sig_ldp ) ;
      rv_mu_WJmc_sig_ldp     -> setVal( NWJmc_sig_ldp ) ;
      rv_mu_Znnmc_sig_ldp    -> setVal( NZnnmc_sig_ldp ) ;
      rv_mu_Ewomc_sig_ldp    -> setVal( NEwomc_sig_ldp ) ;

      rv_mu_ttbarmc_sig_ldp  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sig_ldp     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sig_ldp    -> setConstant( kTRUE ) ;
      rv_mu_Ewomc_sig_ldp    -> setConstant( kTRUE ) ;


     //-- SB, LDP

      rv_mu_ttbarmc_sb_ldp   = new RooRealVar( "mu_ttbarmc_sb_ldp" ,"mu_ttbarmc_sb_ldp" , 0., 1000. ) ;
      rv_mu_WJmc_sb_ldp      = new RooRealVar( "mu_WJmc_sb_ldp"    ,"mu_WJmc_sb_ldp"    , 0., 1000. ) ;
      rv_mu_Znnmc_sb_ldp     = new RooRealVar( "mu_Znnmc_sb_ldp"   ,"mu_Znnmc_sb_ldp"   , 0., 1000. ) ;
      rv_mu_Ewomc_sb_ldp     = new RooRealVar( "mu_Ewomc_sb_ldp"   ,"mu_Ewomc_sb_ldp"   , 0., 1000. ) ;

      rv_mu_ttbarmc_sb_ldp  -> setVal( Nttbarmc_sb_ldp ) ;
      rv_mu_WJmc_sb_ldp     -> setVal( NWJmc_sb_ldp ) ;
      rv_mu_Znnmc_sb_ldp    -> setVal( NZnnmc_sb_ldp ) ;
      rv_mu_Ewomc_sb_ldp    -> setVal( NEwomc_sb_ldp ) ;

      rv_mu_ttbarmc_sb_ldp  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sb_ldp     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sb_ldp    -> setConstant( kTRUE ) ;
      rv_mu_Ewomc_sb_ldp    -> setConstant( kTRUE ) ;


     //-- SIG

      rv_mu_Ewomc_sig     = new RooRealVar( "mu_Ewomc_sig"   ,"mu_Ewomc_sig"   , 0., 1000. ) ;
      rv_mu_Ewomc_sig    -> setVal( NEwomc_sig ) ;
      rv_mu_Ewomc_sig    -> setConstant( kTRUE ) ;

     //-- SB

      ///rv_mu_Ewomc_sb     = new RooRealVar( "mu_Ewomc_sb"   ,"mu_Ewomc_sb"   , 0., 1000. ) ;
      ///rv_mu_Ewomc_sb    -> setVal( NEwomc_sb ) ;
      ///rv_mu_Ewomc_sb    -> setConstant( kTRUE ) ;





    //+++++++ Gaussian constraints ++++++++++++++++++++++++++++++++

      printf(" --- Defining Gaussian constraint and constant parameters.\n" ) ;

    //_______ Efficiency scale factor.  Applied to SUSY and all MC inputs _______________

    //   August 23, 2011:  switching to correlated log-normal PDFs
    //
    //
    //  double pmin, pmax ;


    //--- mean parameters.
      rv_mean_eff_sf_sig     = new RooRealVar( "mean_eff_sf_sig"    , "mean_eff_sf_sig", 0., 10. ) ;
      ///rv_mean_eff_sf_sb      = new RooRealVar( "mean_eff_sf_sb"     , "mean_eff_sf_sb", 0., 10. ) ;
      ///rv_mean_eff_sf_sig_sl  = new RooRealVar( "mean_eff_sf_sig_sl" , "mean_eff_sf_sig_sl", 0., 10. ) ;
      ///rv_mean_eff_sf_sb_sl   = new RooRealVar( "mean_eff_sf_sb_sl"  , "mean_eff_sf_sb_sl", 0., 10. ) ;
      rv_mean_eff_sf_sig_ldp = new RooRealVar( "mean_eff_sf_sig_ldp", "mean_eff_sf_sig_ldp", 0., 10. ) ;
      rv_mean_eff_sf_sb_ldp  = new RooRealVar( "mean_eff_sf_sb_ldp" , "mean_eff_sf_sb_ldp", 0., 10. ) ;

    //--- width parameters.
      rv_width_eff_sf_sig     = new RooRealVar( "width_eff_sf_sig"    , "width_eff_sf_sig", 0., 10. ) ;
      ///rv_width_eff_sf_sb      = new RooRealVar( "width_eff_sf_sb"     , "width_eff_sf_sb", 0., 10. ) ;
      ///rv_width_eff_sf_sig_sl  = new RooRealVar( "width_eff_sf_sig_sl" , "width_eff_sf_sig_sl", 0., 10. ) ;
      ///rv_width_eff_sf_sb_sl   = new RooRealVar( "width_eff_sf_sb_sl"  , "width_eff_sf_sb_sl", 0., 10. ) ;
      rv_width_eff_sf_sig_ldp = new RooRealVar( "width_eff_sf_sig_ldp", "width_eff_sf_sig_ldp", 0., 10. ) ;
      rv_width_eff_sf_sb_ldp  = new RooRealVar( "width_eff_sf_sb_ldp" , "width_eff_sf_sb_ldp", 0., 10. ) ;


      rv_mean_eff_sf_sig     -> setVal( 1.00 ) ;
      ///rv_mean_eff_sf_sb      -> setVal( 1.00 ) ;
      ///rv_mean_eff_sf_sig_sl  -> setVal( 1.00 ) ;
      ///rv_mean_eff_sf_sb_sl   -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_ldp -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_ldp  -> setVal( 1.00 ) ;

      rv_mean_eff_sf_sig     -> setConstant( kTRUE ) ;
      ///rv_mean_eff_sf_sb      -> setConstant( kTRUE ) ;
      ///rv_mean_eff_sf_sig_sl  -> setConstant( kTRUE ) ;
      ///rv_mean_eff_sf_sb_sl   -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_ldp -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_ldp  -> setConstant( kTRUE ) ;


      //--- Initialize all width parameters to 15%.  They will be reset in the susy scan.
      rv_width_eff_sf_sig     -> setVal( 0.15 ) ;
      ///rv_width_eff_sf_sb      -> setVal( 0.15 ) ;
      ///rv_width_eff_sf_sig_sl  -> setVal( 0.15 ) ;
      ///rv_width_eff_sf_sb_sl   -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_ldp -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_ldp  -> setVal( 0.15 ) ;

      rv_width_eff_sf_sig     -> setConstant( kTRUE ) ;
      ///rv_width_eff_sf_sb      -> setConstant( kTRUE ) ;
      ///rv_width_eff_sf_sig_sl  -> setConstant( kTRUE ) ;
      ///rv_width_eff_sf_sb_sl   -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_ldp -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_ldp  -> setConstant( kTRUE ) ;


   //
   // Owen : Sept 7, 2011: need to set these here so that it gets into the workspace.
   //
   //++++++++++
      setSusyScanPoint( inputScanFile,  m0,  m12,  isT1bbbb,  t1bbbbXsec ) ;
   //++++++++++

      gammaModeTransform(Rlsb_passfail,Rlsb_passfail_err,k,theta);
      RooRealVar rrv_Rlsb_passfail ("Rlsb_passfail","Rlsb_passfail",Rlsb_passfail,0,1e5);
      RooRealVar Rlsb_passfail_k ("Rlsb_passfail_k", "Rlsb_passfail_k", k);
      RooRealVar Rlsb_passfail_theta ("Rlsb_passfail_theta", "Rlsb_passfail_theta", theta);
      RooGammaPdf pdf_Rlsb_passfail ("pdf_Rlsb_passfail" , "pdf_Rlsb_passfail", rrv_Rlsb_passfail, Rlsb_passfail_k, Rlsb_passfail_theta, RooConst(0));
      Rlsb_passfail_k.setConstant();
      Rlsb_passfail_theta.setConstant();
      allNonPoissonNuisances.add (rrv_Rlsb_passfail);
      if(constantNonPoisson) rrv_Rlsb_passfail.setConstant();
      allNuisancePdfs.add (pdf_Rlsb_passfail);

    //--- Nov 24, 2011:  sf_ee and sf_mm are now derived from a common function.

      double sigma_avg = 0.5*(sf_ee_err + sf_mm_err);
      double mu_avg = 0.5*(sf_ee + sf_mm);

      gammaModeTransform(mu_avg,sigma_avg,k,theta);

      RooRealVar rrv_sf_ll ("sf_ll","sf_ll",mu_avg,0,1e5);
      RooRealVar sf_ll_k ("sf_ll_k", "sf_ll_k", k);
      RooRealVar sf_ll_theta ("sf_ll_theta", "sf_ll_theta", theta);
      RooGammaPdf pdf_sf_ll ("pdf_sf_ll" , "pdf_sf_ll", rrv_sf_ll, sf_ll_k, sf_ll_theta, RooConst(0));
      sf_ll_k.setConstant();
      sf_ll_theta.setConstant();
      allNonPoissonNuisances.add (rrv_sf_ll);
      if(constantNonPoisson) rrv_sf_ll.setConstant();
      allNuisancePdfs.add (pdf_sf_ll);

      gammaScalingFunctor alphaFinding(sigma_avg,1.);
  
      TF1 alphaFindingFunction("alphaFindingFunction",&alphaFinding,&gammaScalingFunctor::sigmaModeIntersect,0,1,0,"alphaFindingFunction","alphaFindingFunction");
      WrappedTF1 WrappedAlphaFindingFunction(alphaFindingFunction);
      ROOT::Math::Roots::Brent alphaRoot;
      alphaRoot.SetFunction(WrappedAlphaFindingFunction, 1e-9, 5*max(sf_mm_err,sf_ee_err)/sigma_avg);
      alphaFinding.setDesiredSigma(sf_ee_err);
      
      alphaRoot.Solve();

      RooRealVar sf_ee_alpha("sf_ee_alpha","sf_ee_alpha",alphaRoot.Root());
      sf_ee_alpha.setConstant();

      alphaFinding.setDesiredSigma(sf_mm_err);
      alphaRoot.SetFunction(WrappedAlphaFindingFunction, 1e-9, 5*max(sf_mm_err,sf_ee_err)/sigma_avg);
      alphaRoot.Solve();
      RooRealVar sf_mm_alpha("sf_mm_alpha","sf_mm_alpha",alphaRoot.Root());
      sf_mm_alpha.setConstant();

      cout << "sf_ll section -- calculation of power law" << endl;
      cout << "sf_ll section -- sf_ee has alpha of: " << sf_ee_alpha.getVal() << endl;
      cout << "sf_ll section -- sf_mm has alpha of: " << sf_mm_alpha.getVal() << endl;

      RooFormulaVar fv_sf_ee ("sf_ee", "pow(@0,@1)" , RooArgList(rrv_sf_ll,sf_ee_alpha));
      RooFormulaVar fv_sf_mm ("sf_mm", "pow(@0,@1)" , RooArgList(rrv_sf_ll,sf_mm_alpha));

      gammaModeTransform(sf_mc,sf_mc_err,k,theta);
      RooRealVar rrv_sf_mc ("sf_mc","sf_mc",sf_mc,0,1e5);
      RooRealVar sf_mc_k ("sf_mc_k", "sf_mc_k", k);
      RooRealVar sf_mc_theta ("sf_mc_theta", "sf_mc_theta", theta);
      RooGammaPdf pdf_sf_mc ("pdf_sf_mc" , "pdf_sf_mc", rrv_sf_mc, sf_mc_k, sf_mc_theta, RooConst(0));
      sf_mc_k.setConstant();
      sf_mc_theta.setConstant();
      allNonPoissonNuisances.add (rrv_sf_mc);
      if(constantNonPoisson) rrv_sf_mc.setConstant();
      allNuisancePdfs.add (pdf_sf_mc);

      gammaModeTransform(sf_qcd_sig,sf_qcd_sig_err,k,theta);
      RooRealVar rrv_sf_qcd_sig ("sf_qcd_sig","sf_qcd_sig",sf_qcd_sig,0,1e5);
      RooRealVar sf_qcd_sig_k ("sf_qcd_sig_k", "sf_qcd_sig_k", k);
      RooRealVar sf_qcd_sig_theta ("sf_qcd_sig_theta", "sf_qcd_sig_theta", theta);
      RooGammaPdf pdf_sf_qcd_sig ("pdf_sf_qcd_sig" , "pdf_sf_qcd_sig", rrv_sf_qcd_sig, sf_qcd_sig_k, sf_qcd_sig_theta, RooConst(0));
      sf_qcd_sig_k.setConstant();
      sf_qcd_sig_theta.setConstant();
      allNonPoissonNuisances.add (rrv_sf_qcd_sig);
      if(constantNonPoisson) rrv_sf_qcd_sig.setConstant();
      allNuisancePdfs.add (pdf_sf_qcd_sig);

      //gammaModeTransform(sf_ttwj_sig,sf_ttwj_sig_err,k,theta);
      //RooRealVar rrv_sf_ttwj_sig ("sf_ttwj_sig","sf_ttwj_sig",sf_ttwj_sig,0,1e5);
      //RooRealVar sf_ttwj_sig_k ("sf_ttwj_sig_k", "sf_ttwj_sig_k", k);
      //RooRealVar sf_ttwj_sig_theta ("sf_ttwj_sig_theta", "sf_ttwj_sig_theta", theta);
      //RooGammaPdf pdf_sf_ttwj_sig ("pdf_sf_ttwj_sig" , "pdf_sf_ttwj_sig", rrv_sf_ttwj_sig, sf_ttwj_sig_k, sf_ttwj_sig_theta, RooConst(0));
      //sf_ttwj_sig_k.setConstant();
      //sf_ttwj_sig_theta.setConstant();
      //allNuisances.add (rrv_sf_ttwj_sig);
      //allNuisancePdfs.add (pdf_sf_ttwj_sig);

      //Beta Pdf Distributed Variables (efficiencies and acceptances) 

      double alpha,beta;

      betaModeTransform(acc_ee_sig_mean,acc_ee_sig_err,alpha,beta);
      RooRealVar rrv_acc_ee_sig ("acc_ee_sig","acc_ee_sig",acc_ee_sig_mean,0,1);
      RooRealVar acc_ee_sig_alpha ("acc_ee_sig_alpha", "acc_ee_sig_alpha", alpha);
      RooRealVar acc_ee_sig_beta ("acc_ee_sig_beta", "acc_ee_sig_beta", beta);
      RooBetaPdf pdf_acc_ee_sig ("pdf_acc_ee_sig" , "pdf_acc_ee_sig", rrv_acc_ee_sig, acc_ee_sig_alpha, acc_ee_sig_beta);
      acc_ee_sig_alpha.setConstant();
      acc_ee_sig_beta.setConstant();
      allNonPoissonNuisances.add (rrv_acc_ee_sig);
      if(constantNonPoisson) rrv_acc_ee_sig.setConstant();
      allNuisancePdfs.add (pdf_acc_ee_sig);

      betaModeTransform(acc_mm_sig_mean,acc_mm_sig_err,alpha,beta);
      RooRealVar rrv_acc_mm_sig ("acc_mm_sig","acc_mm_sig",acc_mm_sig_mean,0,1);
      RooRealVar acc_mm_sig_alpha ("acc_mm_sig_alpha", "acc_mm_sig_alpha", alpha);
      RooRealVar acc_mm_sig_beta ("acc_mm_sig_beta", "acc_mm_sig_beta", beta);
      RooBetaPdf pdf_acc_mm_sig ("pdf_acc_mm_sig" , "pdf_acc_mm_sig", rrv_acc_mm_sig, acc_mm_sig_alpha, acc_mm_sig_beta);
      acc_mm_sig_alpha.setConstant();
      acc_mm_sig_beta.setConstant();
      allNonPoissonNuisances.add (rrv_acc_mm_sig);
      if(constantNonPoisson) rrv_acc_mm_sig.setConstant();
      allNuisancePdfs.add (pdf_acc_mm_sig);

      betaModeTransform(eff_ee_mean,eff_ee_err,alpha,beta);
      RooRealVar rrv_eff_ee ("eff_ee","eff_ee",eff_ee_mean,0,1);
      RooRealVar eff_ee_alpha ("eff_ee_alpha", "eff_ee_alpha", alpha);
      RooRealVar eff_ee_beta ("eff_ee_beta", "eff_ee_beta", beta);
      RooBetaPdf pdf_eff_ee ("pdf_eff_ee" , "pdf_eff_ee", rrv_eff_ee, eff_ee_alpha, eff_ee_beta);
      eff_ee_alpha.setConstant();
      eff_ee_beta.setConstant();
      allNonPoissonNuisances.add (rrv_eff_ee);
      if(constantNonPoisson) rrv_eff_ee.setConstant();
      allNuisancePdfs.add (pdf_eff_ee);

      betaModeTransform(eff_mm_mean,eff_mm_err,alpha,beta);
      RooRealVar rrv_eff_mm ("eff_mm","eff_mm",eff_mm_mean,0,1);
      RooRealVar eff_mm_alpha ("eff_mm_alpha", "eff_mm_alpha", alpha);
      RooRealVar eff_mm_beta ("eff_mm_beta", "eff_mm_beta", beta);
      RooBetaPdf pdf_eff_mm ("pdf_eff_mm" , "pdf_eff_mm", rrv_eff_mm, eff_mm_alpha, eff_mm_beta);
      eff_mm_alpha.setConstant();
      eff_mm_beta.setConstant();
      allNonPoissonNuisances.add (rrv_eff_mm);
      if(constantNonPoisson) rrv_eff_mm.setConstant();
      allNuisancePdfs.add (pdf_eff_mm);

      betaModeTransform(fsig_ee_mean,fsig_ee_err,alpha,beta);
      RooRealVar rrv_fsig_ee ("fsig_ee","fsig_ee",fsig_ee_mean,0,1);
      RooRealVar fsig_ee_alpha ("fsig_ee_alpha", "fsig_ee_alpha", alpha);
      RooRealVar fsig_ee_beta ("fsig_ee_beta", "fsig_ee_beta", beta);
      RooBetaPdf pdf_fsig_ee ("pdf_fsig_ee" , "pdf_fsig_ee", rrv_fsig_ee, fsig_ee_alpha, fsig_ee_beta);
      fsig_ee_alpha.setConstant();
      fsig_ee_beta.setConstant();
      allNonPoissonNuisances.add (rrv_fsig_ee);
      if(constantNonPoisson) rrv_fsig_ee.setConstant();
      allNuisancePdfs.add (pdf_fsig_ee);

      betaModeTransform(fsig_mm_mean,fsig_mm_err,alpha,beta);
      RooRealVar rrv_fsig_mm ("fsig_mm","fsig_mm",fsig_mm_mean,0,1);
      RooRealVar fsig_mm_alpha ("fsig_mm_alpha", "fsig_mm_alpha", alpha);
      RooRealVar fsig_mm_beta ("fsig_mm_beta", "fsig_mm_beta", beta);
      RooBetaPdf pdf_fsig_mm ("pdf_fsig_mm" , "pdf_fsig_mm", rrv_fsig_mm, fsig_mm_alpha, fsig_mm_beta);
      fsig_mm_alpha.setConstant();
      fsig_mm_beta.setConstant();
      allNonPoissonNuisances.add (rrv_fsig_mm);
      if(constantNonPoisson) rrv_fsig_mm.setConstant();
      allNuisancePdfs.add (pdf_fsig_mm);

      betaModeTransform(knn_ee_sig_mean,knn_ee_sig_err,alpha,beta);
      RooRealVar rrv_knn_ee_sig ("knn_ee_sig","knn_ee_sig",knn_ee_sig_mean,0,1);
      RooRealVar knn_ee_sig_alpha ("knn_ee_sig_alpha", "knn_ee_sig_alpha", alpha);
      RooRealVar knn_ee_sig_beta ("knn_ee_sig_beta", "knn_ee_sig_beta", beta);
      RooBetaPdf pdf_knn_ee_sig ("pdf_knn_ee_sig" , "pdf_knn_ee_sig", rrv_knn_ee_sig, knn_ee_sig_alpha, knn_ee_sig_beta);
      knn_ee_sig_alpha.setConstant();
      knn_ee_sig_beta.setConstant();
      allNonPoissonNuisances.add (rrv_knn_ee_sig);
      if(constantNonPoisson) rrv_knn_ee_sig.setConstant();
      allNuisancePdfs.add (pdf_knn_ee_sig);

      betaModeTransform(knn_mm_sig_mean,knn_mm_sig_err,alpha,beta);
      RooRealVar rrv_knn_mm_sig ("knn_mm_sig","knn_mm_sig",knn_mm_sig_mean,0,1);
      RooRealVar knn_mm_sig_alpha ("knn_mm_sig_alpha", "knn_mm_sig_alpha", alpha);
      RooRealVar knn_mm_sig_beta ("knn_mm_sig_beta", "knn_mm_sig_beta", beta);
      RooBetaPdf pdf_knn_mm_sig ("pdf_knn_mm_sig" , "pdf_knn_mm_sig", rrv_knn_mm_sig, knn_mm_sig_alpha, knn_mm_sig_beta);
      knn_mm_sig_alpha.setConstant();
      knn_mm_sig_beta.setConstant();
      allNonPoissonNuisances.add (rrv_knn_mm_sig);
      if(constantNonPoisson) rrv_knn_mm_sig.setConstant();
      allNuisancePdfs.add (pdf_knn_mm_sig);

    //-------  add trigger efficiencies as nuisance parameters with Beta Funtions

      betaModeTransform(eps_sig_mean-epsSF_sig_errm+0.5*(epsSF_sig_errm+epsSF_sig_errp),epsSF_sig_errm+epsSF_sig_errp,alpha,beta);
      RooRealVar rrv_epsSF_sig ("epsSF_sig","epsSF_sig",eps_sig_mean-epsSF_sig_errm+0.5*(epsSF_sig_errm+epsSF_sig_errp),0,1);
      RooRealVar epsSF_sig_alpha ("epsSF_sig_alpha", "epsSF_sig_alpha", alpha);
      RooRealVar epsSF_sig_beta ("epsSF_sig_beta", "epsSF_sig_beta", beta);
      RooBetaPdf pdf_epsSF_sig ("pdf_epsSF_sig" , "pdf_epsSF_sig", rrv_epsSF_sig, epsSF_sig_alpha, epsSF_sig_beta);
      epsSF_sig_alpha.setConstant();
      epsSF_sig_beta.setConstant();
      allNonPoissonNuisances.add (rrv_epsSF_sig);
      if(constantNonPoisson) rrv_epsSF_sig.setConstant();
      allNuisancePdfs.add (pdf_epsSF_sig);

      betaModeTransform(eps_sb_ldp_mean-epsSF_sb_ldp_errm+0.5*(epsSF_sb_ldp_errm+epsSF_sb_ldp_errp),epsSF_sb_ldp_errm+epsSF_sb_ldp_errp,alpha,beta);
      RooRealVar rrv_epsSF_sb_ldp ("epsSF_sb_ldp","epsSF_sb_ldp",eps_sb_ldp_mean-epsSF_sb_ldp_errm+0.5*(epsSF_sb_ldp_errm+epsSF_sb_ldp_errp),0,1);
      RooRealVar epsSF_sb_ldp_alpha ("epsSF_sb_ldp_alpha", "epsSF_sb_ldp_alpha", alpha);
      RooRealVar epsSF_sb_ldp_beta ("epsSF_sb_ldp_beta", "epsSF_sb_ldp_beta", beta);
      RooBetaPdf pdf_epsSF_sb_ldp ("pdf_epsSF_sb_ldp" , "pdf_epsSF_sb_ldp", rrv_epsSF_sb_ldp, epsSF_sb_ldp_alpha, epsSF_sb_ldp_beta);
      epsSF_sb_ldp_alpha.setConstant();
      epsSF_sb_ldp_beta.setConstant();
      allNonPoissonNuisances.add (rrv_epsSF_sb_ldp);
      if(constantNonPoisson) rrv_epsSF_sb_ldp.setConstant();
      allNuisancePdfs.add (pdf_epsSF_sb_ldp);

    //--------------------------------------


      // Underlying Gamma Distribution for signal systematics -- taken from average of all errors.

      double sigma_eff_sf_sig      = rv_width_eff_sf_sig    ->getVal();
      //double sigma_eff_sf_sb       = rv_width_eff_sf_sb	    ->getVal();
      //double sigma_eff_sf_sig_sl   = rv_width_eff_sf_sig_sl ->getVal();
      //double sigma_eff_sf_sb_sl    = rv_width_eff_sf_sb_sl  ->getVal();
      double sigma_eff_sf_sig_ldp  = rv_width_eff_sf_sig_ldp->getVal();
      double sigma_eff_sf_sb_ldp   = rv_width_eff_sf_sb_ldp ->getVal();

      double sigma_max = max(sigma_eff_sf_sig,
			     max(sigma_eff_sf_sig_ldp,
				 sigma_eff_sf_sb_ldp));

      sigma_avg =(sigma_eff_sf_sig + sigma_eff_sf_sig_ldp + sigma_eff_sf_sb_ldp)/3.;
      //sigma_avg = sigma_eff_sf_sig;
      mu_avg = 1.;

      gammaModeTransform(mu_avg,sigma_avg,k,theta);

      //RooRealVar rrv_eff_sf ("eff_sf","eff_sf",mu_avg,0,1e5);
      //RooRealVar eff_sf_k ("eff_sf_k", "eff_sf_k", k);
      //RooRealVar eff_sf_theta ("eff_sf_theta", "eff_sf_theta", theta);
      //RooFormulaVar inv_mu_susy_sig_expected("inv_mu_susy_sig_expected","(@0>0)*1/@0",RooArgSet(*rv_mu_susy_sig));
      //cout << "value of rv_mu_susy_sig is " << rv_mu_susy_sig->getVal() << endl;
      //cout << "value of inv_mu_susy_sig_expected is " << inv_mu_susy_sig_expected.getVal() << endl;
      ////inv_mu_susy_sig_expected.setConstant();
      //RooProduct signalEfficiency_eff_sf("eff_sf","eff_sf",RooArgSet(*rv_mu_susy_sig,inv_mu_susy_sig_expected));
      //RooAbsReal* eff_sf_pointer;
      //if(/*sigPOI*/true) eff_sf_pointer = &rrv_eff_sf;
      //else eff_sf_pointer = &signalEfficiency_eff_sf ;
      //RooGammaPdf pdf_eff_sf ("pdf_eff_sf" , "pdf_eff_sf", *eff_sf_pointer, eff_sf_k, eff_sf_theta, RooConst(0));
      //eff_sf_k.setConstant();
      //eff_sf_theta.setConstant();
      //if(!effInitFixed) allNuisances.add (rrv_eff_sf); 
      //if(constantSigEff) rrv_eff_sf.setConstant();
      //allNuisancePdfs.add (pdf_eff_sf);

      RooRealVar rrv_eff_sf ("eff_sf","eff_sf",mu_avg,0,1e5);
      RooRealVar eff_sf_k ("eff_sf_k", "eff_sf_k", k);
      RooRealVar eff_sf_theta ("eff_sf_theta", "eff_sf_theta", theta);
      RooGammaPdf pdf_eff_sf ("pdf_eff_sf" , "pdf_eff_sf", rrv_eff_sf, eff_sf_k, eff_sf_theta, RooConst(0));
      eff_sf_k.setConstant();
      eff_sf_theta.setConstant();
      if(!effInitFixed) allNonPoissonNuisances.add (rrv_eff_sf); 
      if(constantSigEff) rrv_eff_sf.setConstant();
      allNuisancePdfs.add (pdf_eff_sf);

      alphaFinding.setCurrentSigma(sigma_avg);
  
      if(sigma_eff_sf_sig > 0){
	alphaFinding.setDesiredSigma(sigma_eff_sf_sig);
	alphaRoot.SetFunction(WrappedAlphaFindingFunction, 1e-6, 5*sigma_max/sigma_avg);
	alphaRoot.Solve();
      }
      RooRealVar eff_sf_sig_alpha("eff_sf_sig_alpha","eff_sf_sig_alpha", (sigma_eff_sf_sig > 0) ? alphaRoot.Root() : 0 );
      eff_sf_sig_alpha.setConstant();

      if(sigma_eff_sf_sig_ldp > 0){
	alphaFinding.setDesiredSigma(sigma_eff_sf_sig_ldp);
	alphaRoot.SetFunction(WrappedAlphaFindingFunction, 1e-6, 5*sigma_max/sigma_avg);
	alphaRoot.Solve();
      }
      RooRealVar eff_sf_sig_ldp_alpha("eff_sf_sig_ldp_alpha","eff_sf_sig_ldp_alpha", (sigma_eff_sf_sig_ldp > 0) ? alphaRoot.Root() : 0 );
      eff_sf_sig_ldp_alpha.setConstant();

      if(sigma_eff_sf_sb_ldp > 0){
	alphaFinding.setDesiredSigma(sigma_eff_sf_sb_ldp);
	alphaRoot.SetFunction(WrappedAlphaFindingFunction, 1e-6, 5*sigma_max/sigma_avg);
	alphaRoot.Solve();
      }
      RooRealVar eff_sf_sb_ldp_alpha("eff_sf_sb_ldp_alpha","eff_sf_sb_ldp_alpha", (sigma_eff_sf_sb_ldp > 0) ? alphaRoot.Root() : 0 );
      eff_sf_sb_ldp_alpha.setConstant();

      cout << "eff_sf section -- calculation of power law" << endl;
      cout << "eff_sf section -- eff_sf_sig     has alpha of: " << eff_sf_sig_alpha    .getVal() << endl;
      cout << "eff_sf section -- eff_sf_sig_ldp has alpha of: " << eff_sf_sig_ldp_alpha.getVal() << endl;
      cout << "eff_sf section -- eff_sf_sb_ldp  has alpha of: " << eff_sf_sb_ldp_alpha .getVal() << endl;

      //-- Parametric relations between correlated signal efficiency scale factors.

      RooFormulaVar fv_eff_sf_sig     ("eff_sf_sig"     , "pow(@0,@1)" , RooArgList( rrv_eff_sf , eff_sf_sig_alpha      ));
      RooFormulaVar fv_eff_sf_sig_ldp ("eff_sf_sig_ldp" , "pow(@0,@1)" , RooArgList( rrv_eff_sf , eff_sf_sig_ldp_alpha  ));
      RooFormulaVar fv_eff_sf_sb_ldp  ("eff_sf_sb_ldp"  , "pow(@0,@1)" , RooArgList( rrv_eff_sf , eff_sf_sb_ldp_alpha   ));

    //-- Z to nunu stuff



      rv_znnoverll_bfratio = new RooRealVar( "znnoverll_bfratio", "znnoverll_bfratio", 0.1, 10. ) ;
      rv_znnoverll_bfratio -> setVal( 5.95 ) ;
      rv_znnoverll_bfratio -> setConstant( kTRUE ) ;

      rv_dataoverll_lumiratio = new RooRealVar( "dataoverll_lumiratio", "dataoverll_lumiratio", 0.1, 10.0 ) ;
      rv_dataoverll_lumiratio  -> setVal( DataLumi / Ztoll_lumi ) ;
      rv_dataoverll_lumiratio  -> setConstant( kTRUE ) ;




     //+++++++++++++++++ Relationships between parameters ++++++++++++++++++++++++++++++++++++++++++

       printf(" --- Defining relationships between parameters.\n" ) ;




    //-- ttwj

      if ( useSigTtwjVar ) {

      } else {
	cout << "This option is not valid with this likelihood.  Please set useSigTtwjVar to true." << endl;
	return false;

      }




      rv_mu_ttwj_sig_ldp = new RooFormulaVar("mu_ttwj_sig_ldp",
                              "mu_ttbarmc_sig_ldp + mu_WJmc_sig_ldp",
                              RooArgSet( *rv_mu_ttbarmc_sig_ldp,
                                         *rv_mu_WJmc_sig_ldp ) ) ;


      rv_mu_ttwj_sb_ldp = new RooFormulaVar("mu_ttwj_sb_ldp",
                              "mu_ttbarmc_sb_ldp + mu_WJmc_sb_ldp",
                              RooArgSet( *rv_mu_ttbarmc_sb_ldp,
                                         *rv_mu_WJmc_sb_ldp ) ) ;






    //-- QCD

   //-------------------------------
   ///if ( useLdpVars ) {

   ///   rfv_mu_qcd_sig = new RooFormulaVar( "mu_qcd_sig",
   ///                               "mu_qcd_sig_ldp * sf_qcd_sig * ( mu_qcd_lsb_0b / mu_qcd_lsb_0b_ldp )",
   ///                               RooArgSet( *rv_mu_qcd_sig_ldp, fv_sf_qcd_sig, *rv_mu_qcd_lsb_0b, *rv_mu_qcd_lsb_0b_ldp ) ) ;
   ///   rv_mu_qcd_sig = rfv_mu_qcd_sig ;

   ///   rfv_mu_qcd_sb = new RooFormulaVar( "mu_qcd_sb",
   ///                               "mu_qcd_sb_ldp * sf_qcd_sb * ( mu_qcd_lsb_0b / mu_qcd_lsb_0b_ldp )",
   ///                               RooArgSet( *rv_mu_qcd_sb_ldp, fv_sf_qcd_sb, *rv_mu_qcd_lsb_0b, *rv_mu_qcd_lsb_0b_ldp ) ) ;
   ///   rv_mu_qcd_sb = rfv_mu_qcd_sb ;

   ///} else {

   ///   rfv_mu_qcd_sig_ldp = new RooFormulaVar( "mu_qcd_sig_ldp",
   ///                               "mu_qcd_sig * (1.0/sf_qcd_sig) * ( mu_qcd_lsb_0b_ldp / mu_qcd_lsb_0b )",
   ///                               RooArgSet( *rv_mu_qcd_sig, fv_sf_qcd_sig, *rv_mu_qcd_lsb_0b_ldp, *rv_mu_qcd_lsb_0b ) ) ;
   ///   rv_mu_qcd_sig_ldp = rfv_mu_qcd_sig_ldp ;

   ///   rfv_mu_qcd_sb_ldp = new RooFormulaVar( "mu_qcd_sb_ldp",
   ///                               "mu_qcd_sb * (1.0/sf_qcd_sb) * ( mu_qcd_lsb_0b_ldp / mu_qcd_lsb_0b )",
   ///                               RooArgSet( *rv_mu_qcd_sb, fv_sf_qcd_sb, *rv_mu_qcd_lsb_0b_ldp, *rv_mu_qcd_lsb_0b ) ) ;
   ///   rv_mu_qcd_sb_ldp = rfv_mu_qcd_sb_ldp ;

   ///}
   //-------------------------------

      if ( useLdpVars ) {

         rfv_mu_qcd_sig = new RooFormulaVar("mu_qcd_sig",
					     "mu_qcd_sig_ldp * sf_qcd_sig * Rlsb_passfail * (1.0/epsSF_sig)",
					     RooArgSet( *rv_mu_qcd_sig_ldp, rrv_sf_qcd_sig, rrv_Rlsb_passfail, rrv_epsSF_sig ) ) ;
         rv_mu_qcd_sig = rfv_mu_qcd_sig ;


      } else {
	cout <<"This option is invalid, please set useLdpVars to true)"<<endl;
	return false;

      }

   //-------------------------------



    //-- SUSY


      rv_mu_susy_sig_ldp = new RooFormulaVar("mu_susy_sig_ldp",
					      "mu_susymc_sig_ldp * (mu_susy_sig/mu_susymc_sig)",
					      RooArgSet( *rv_mu_susymc_sig_ldp, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sb_ldp = new RooFormulaVar("mu_susy_sb_ldp",
					     "mu_susymc_sb_ldp * (mu_susy_sig/mu_susymc_sig) * (1.0/epsSF_sb_ldp)",
					     RooArgSet( *rv_mu_susymc_sb_ldp, *rv_mu_susy_sig, *rv_mu_susymc_sig, rrv_epsSF_sb_ldp ) ) ;






    //-- Z to nu nu

      rv_mu_znn_sig_ldp = new RooFormulaVar( "mu_znn_sig_ldp",
                              "mu_Znnmc_sig_ldp",
                              RooArgSet( *rv_mu_Znnmc_sig_ldp ) ) ;

      rv_mu_znn_sb_ldp = new RooFormulaVar( "mu_znn_sb_ldp",
                              "mu_Znnmc_sb_ldp",
                              RooArgSet( *rv_mu_Znnmc_sb_ldp ) ) ;


   //---------------------------------
   // if ( znnModel == 1 ) {

   //    rv_mu_zee_sb_ee = new RooFormulaVar( "mu_zee_sb_ee",
   //                                    "mu_znn_sb * sf_ee * ( acc_ee * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
   //                                    RooArgSet( *rv_mu_znn_sb, fv_sf_ee, fv_acc_ee, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

   //    rv_mu_zee_sig_ee = new RooFormulaVar( "mu_zee_sig_ee",
   //                                    "mu_znn_sig * sf_ee * ( acc_ee * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
   //                                    RooArgSet( *rv_mu_znn_sig, fv_sf_ee, fv_acc_ee, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

   //    rv_mu_zmm_sb_mm = new RooFormulaVar( "mu_zmm_sb_mm",
   //                                    "mu_znn_sb * sf_mm * ( acc_mm * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
   //                                    RooArgSet( *rv_mu_znn_sb, fv_sf_mm, fv_acc_mm, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

   //    rv_mu_zmm_sig_mm = new RooFormulaVar( "mu_zmm_sig_mm",
   //                                    "mu_znn_sig * sf_mm * ( acc_mm * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
   //                                    RooArgSet( *rv_mu_znn_sig, fv_sf_mm, fv_acc_mm, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

   // } else if ( znnModel == 2 ) {

   //    rfv_mu_znn_sb = new RooFormulaVar( "mu_znn_sb",
   //                                     "mu_znn_sig * ( knn_sb / knn_sig )",
   //                                     RooArgSet( *rv_mu_znn_sig, fv_knn_sb, fv_knn_sig ) ) ;

   //    rv_mu_znn_sb = rfv_mu_znn_sb ;

   //    rv_mu_zee_sigsb_ee = new RooFormulaVar( "mu_zee_sigsb_ee",
   //                                  "( mu_znn_sig / knn_sig ) * sf_ee * ( (acc_ee * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
   //                                    RooArgSet( *rv_mu_znn_sig, fv_knn_sig, fv_sf_ee, fv_acc_ee, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

   //    rv_mu_zmm_sigsb_mm = new RooFormulaVar( "mu_zmm_sigsb_mm",
   //                                  "( mu_znn_sig / knn_sig ) * sf_mm * ( (acc_mm * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
   //                                    RooArgSet( *rv_mu_znn_sig, fv_knn_sig, fv_sf_mm, fv_acc_mm, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

   // }
   //---------------------------------

      rv_mu_zee_sig_ee = new RooFormulaVar("mu_zee_sig_ee",
                                    "( mu_znn_sig / knn_ee_sig ) * sf_ee * ( (acc_ee_sig * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                      RooArgSet( *rv_mu_znn_sig, rrv_knn_ee_sig, fv_sf_ee, rrv_acc_ee_sig, rrv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;


      rv_mu_zmm_sig_mm = new RooFormulaVar("mu_zmm_sig_mm",
                                    "( mu_znn_sig / knn_mm_sig ) * sf_mm * ( (acc_mm_sig * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                      RooArgSet( *rv_mu_znn_sig, rrv_knn_mm_sig, fv_sf_mm, rrv_acc_mm_sig, rrv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;


    //-- EWO

      rv_mu_ewo_sig     = new RooFormulaVar( "mu_ewo_sig"     , "mu_Ewomc_sig"     , RooArgSet( *rv_mu_Ewomc_sig     ) ) ;
      rv_mu_ewo_sig_ldp = new RooFormulaVar( "mu_ewo_sig_ldp" , "mu_Ewomc_sig_ldp" , RooArgSet( *rv_mu_Ewomc_sig_ldp ) ) ;
      rv_mu_ewo_sb_ldp  = new RooFormulaVar( "mu_ewo_sb_ldp"  , "mu_Ewomc_sb_ldp"  , RooArgSet( *rv_mu_Ewomc_sb_ldp  ) ) ;



    //+++++++++++++ Expected counts for observables in terms of parameters ++++++++++++++++++

       printf(" --- Defining expected counts in terms of parameters.\n" ) ;

      rv_n_sig         = new RooFormulaVar("n_sig",
					    "mu_ttwj_sig + mu_qcd_sig + mu_znn_sig + sf_mc * ( mu_ewo_sig ) + eff_sf_sig * ( mu_susy_sig )",
					    RooArgSet( *rv_mu_ttwj_sig, *rv_mu_qcd_sig, *rv_mu_znn_sig, rrv_sf_mc, *rv_mu_ewo_sig, fv_eff_sf_sig, *rv_mu_susy_sig ) ) ;
      

      rv_n_sig_ldp     = new RooFormulaVar("n_sig_ldp",
					    "mu_qcd_sig_ldp + sf_mc * (mu_ttwj_sig_ldp + mu_znn_sig_ldp + mu_ewo_sig_ldp) + eff_sf_sig_ldp*(mu_susy_sig_ldp)",
					    RooArgSet( *rv_mu_qcd_sig_ldp, fv_eff_sf_sig_ldp, rrv_sf_mc, *rv_mu_ttwj_sig_ldp, *rv_mu_znn_sig_ldp, *rv_mu_ewo_sig_ldp, *rv_mu_susy_sig_ldp ) ) ;

      rv_n_sb_ldp      = new RooFormulaVar("n_sb_ldp",
					    "mu_qcd_sb_ldp + sf_mc * (mu_ttwj_sb_ldp + mu_znn_sb_ldp + mu_ewo_sb_ldp) + eff_sf_sb_ldp*(mu_susy_sb_ldp)",
					    RooArgSet( *rv_mu_qcd_sb_ldp, fv_eff_sf_sb_ldp, rrv_sf_mc, *rv_mu_ttwj_sb_ldp, *rv_mu_znn_sb_ldp, *rv_mu_ewo_sb_ldp, *rv_mu_susy_sb_ldp ) ) ;
      
      rv_n_sig_ee      = new RooFormulaVar("n_sig_ee",
                                     "mu_zee_sig_ee / fsig_ee",
                                     RooArgSet( *rv_mu_zee_sig_ee, rrv_fsig_ee ) ) ;

      rv_n_sig_mm      = new RooFormulaVar("n_sig_mm",
                                     "mu_zmm_sig_mm / fsig_mm",
                                     RooArgSet( *rv_mu_zmm_sig_mm, rrv_fsig_mm ) ) ;

   //++++++++++++ PDFs for the likelihood +++++++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining PDFs of the likelihood.\n" ) ;

      pdf_Nsig        = new RooPoisson( "pdf_Nsig"        , "Nsig Poisson PDF"        , *rv_Nsig        , *rv_n_sig ) ;
      pdf_Nsig_ldp    = new RooPoisson( "pdf_Nsig_ldp"    , "Nsig_ldp Poisson PDF"    , *rv_Nsig_ldp    , *rv_n_sig_ldp ) ;
      pdf_Nsb_ldp     = new RooPoisson( "pdf_Nsb_ldp"     , "Nsb_ldp Poisson PDF"     , *rv_Nsb_ldp     , *rv_n_sb_ldp ) ;

      pdf_Nsig_ee     = new RooPoisson( "pdf_Nsig_ee"     , "Nsig_ee Poisson PDF"     , *rv_Nsig_ee     , *rv_n_sig_ee ) ;
      pdf_Nsig_mm     = new RooPoisson( "pdf_Nsig_mm"     , "Nsig_mm Poisson PDF"     , *rv_Nsig_mm     , *rv_n_sig_mm ) ;

      {
         RooArgSet pdflist ;
         pdflist.add( *pdf_Nsig        ) ;
         pdflist.add( *pdf_Nsig_ldp    ) ;
         pdflist.add( *pdf_Nsb_ldp     ) ;
         pdflist.add( *pdf_Nsig_ee     ) ;
         pdflist.add( *pdf_Nsig_mm     ) ;

	 pdflist.add(allNuisancePdfs);

         likelihood = new RooProdPdf( "likelihood", "ra2b likelihood", pdflist ) ;

         RooArgSet pdflist_noNSig ;
         pdflist_noNSig.add( *pdf_Nsig_ldp    ) ;
         pdflist_noNSig.add( *pdf_Nsb_ldp     ) ;
         pdflist_noNSig.add( *pdf_Nsig_ee     ) ;
         pdflist_noNSig.add( *pdf_Nsig_mm     ) ;

	 pdflist_noNSig.add(allNuisancePdfs);

         likelihood_noNSig = new RooProdPdf( "likelihood_noNSig", "ra2b likelihood w/out NSig", pdflist_noNSig ) ;

      }


     //---- Define the list of observables.

       observedParametersList.add( *rv_Nsig        ) ;
       observedParametersList.add( *rv_Nsig_ldp    ) ;
       observedParametersList.add( *rv_Nsb_ldp     ) ;
       observedParametersList.add( *rv_Nsig_ee     ) ;
       observedParametersList.add( *rv_Nsig_mm     ) ;



       dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
                                      observedParametersList ) ;
       dsObserved->add( observedParametersList ) ;


       RooWorkspace workspace ("ws") ;
       workspace.autoImportClassCode(true);

       workspace.import(*dsObserved);

       // parameters of interest
       RooArgSet poi(*rv_mu_susy_sig, "poi");
       // flat prior for POI
       RooUniform signal_prior ("signal_prior","signal_prior",*rv_mu_susy_sig);

       allNuisances.add(allPoissonNuisances);
       allNuisances.add(allNonPoissonNuisances);

       // signal+background model
       ModelConfig sbModel ("SbModel");
       sbModel.SetWorkspace(workspace);
       sbModel.SetPdf(*likelihood);
       sbModel.SetParametersOfInterest(poi);
       sbModel.SetPriorPdf(signal_prior);
       sbModel.SetNuisanceParameters(allNuisances);
       sbModel.SetObservables(observedParametersList);
       sbModel.SetGlobalObservables(globalObservables);

       // find global maximum with the signal+background model
       // with conditional MLEs for nuisance parameters
       // and save the parameter point snapshot in the Workspace
       //  - safer to keep a default name because some RooStats calculators
       //    will anticipate it
       RooAbsReal * pNll = sbModel.GetPdf()->createNLL(*dsObserved);
       RooAbsReal * pProfile = pNll->createProfile(RooArgSet());
       pProfile->getVal(); // this will do fit and set POI and nuisance parameters to fitted values
       RooArgSet * pPoiAndNuisance = new RooArgSet();
       pPoiAndNuisance->add(*sbModel.GetParametersOfInterest());
       if(sbModel.GetNuisanceParameters()) pPoiAndNuisance->add(*sbModel.GetNuisanceParameters());
       cout << "\nWill save these parameter points that correspond to the fit to data" << endl;
       pPoiAndNuisance->Print("v");
       sbModel.SetSnapshot(*pPoiAndNuisance);
       workspace.import (sbModel);

       delete pProfile;
       delete pNll;
       delete pPoiAndNuisance;


       // background-only model
       // use the same PDF as s+b, with xsec=0
       // POI value under the background hypothesis
       ModelConfig bModel (*(RooStats::ModelConfig *)workspace.obj("SbModel"));
       bModel.SetName("BModel");
       bModel.SetWorkspace(workspace);

       // Find a parameter point for generating pseudo-data
       // with the background-only data.
       // Save the parameter point snapshot in the Workspace
       pNll = bModel.GetPdf()->createNLL(*dsObserved);
       // bug discovered by Fedor on Sep 6th 2011:
       //pProfile = pNll->createProfile(poi);
       //((RooRealVar *)poi.first())->setVal(0.); // set signal = 0
       pProfile = pNll->createProfile(*bModel.GetParametersOfInterest());
       ((RooRealVar *)(bModel.GetParametersOfInterest()->first()))->setVal(0.); // set signal = 0
       pProfile->getVal(); // this will do fit and set nuisance parameters to profiled values
       pPoiAndNuisance = new RooArgSet();
       pPoiAndNuisance->add(*bModel.GetParametersOfInterest());
       if(bModel.GetNuisanceParameters()) pPoiAndNuisance->add(*bModel.GetNuisanceParameters());
       cout << "\nShould use these parameter points to generate pseudo data for bkg only" << endl;
       pPoiAndNuisance->Print("v");
       bModel.SetSnapshot(*pPoiAndNuisance);
       workspace.import (bModel);

       workspace.import(*likelihood_noNSig, RecycleConflictNodes());

       RooArgSet nuiscancePriors;
       //nuiscancePriors.add(gammaNuisancePdfs);
       //nuiscancePriors.add(betaNuisancePdfs);
       nuiscancePriors.add(allNuisancePdfs);

       RooProdPdf nuisancePrior("nuisancePrior", "nuisances we have priors for", nuiscancePriors ) ;
       workspace.import(nuisancePrior, RecycleConflictNodes());

       delete pProfile;
       delete pNll;
       delete pPoiAndNuisance;

       workspace.defineSet("allNonPoissonNuisances",allNonPoissonNuisances);

       workspace.Print() ;
       workspace.writeToFile(outfile);

       return true ;


    } // initialize.


   //===================================================================================================================================


    bool ra2bRoostatsClass7::setSusyScanPoint( const char* inputScanFile, double m0, double m12, bool isT1bbbb, double t1bbbbXsec ) {


       //--- Aug 15, 2011: updated to new format for AN, v3.


       printf("\n\n Opening SUSY scan input file : %s\n", inputScanFile ) ;

       FILE* infp ;
       if ( (infp=fopen( inputScanFile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", inputScanFile ) ;
          return false ;
       }

       double deltaM0(0.) ;
       double deltaM12(0.) ;

       if ( !isT1bbbb ) {
          deltaM0 = 20 ;
          deltaM12 = 20 ;
       } else {
          deltaM0 = 25 ;
          deltaM12 = 25 ;
       }

       bool found(false) ;

       //--- Loop over the scan points.
       while ( !feof( infp ) ) {

          float pointM0 ;
          float pointM12 ;



   //+++ ORL: Aug 14, 2011 ++++++++++++++++++++++++++
          float n_sig_raw ;
          float n_sb_raw ;
          float n_sig_sl_raw ;
          float n_sb_sl_raw ;
          float n_sig_ldp_raw ;
          float n_sb_ldp_raw ;

          float n_sig_correction ;
          float n_sb_correction ;
          float n_sig_sl_correction ;
          float n_sb_sl_correction ;
          float n_sig_ldp_correction ;
          float n_sb_ldp_correction ;

          float n_sig_error ;
          float n_sb_error ;
          float n_sig_sl_error ;
          float n_sb_sl_error ;
          float n_sig_ldp_error ;
          float n_sb_ldp_error ;

          int nGen ;

          fscanf( infp, "%f %f %d  %f %f %f %f %f %f  %f %f %f %f %f %f %f %f %f %f %f %f",
            &pointM0, &pointM12, &nGen,
            &n_sig_raw, &n_sb_raw, &n_sig_sl_raw, &n_sb_sl_raw, &n_sig_ldp_raw, &n_sb_ldp_raw,
            &n_sig_correction, &n_sb_correction, &n_sig_sl_correction, &n_sb_sl_correction, &n_sig_ldp_correction, &n_sb_ldp_correction,
            &n_sig_error, &n_sb_error, &n_sig_sl_error, &n_sb_sl_error, &n_sig_ldp_error, &n_sb_ldp_error ) ;

          if ( feof(infp) ) break ;
          if ( n_sig_raw < 0.00001 ) continue ;
   //--If you are asking for it, I'll assume it's good.  Josh is using 0 for ngen dummy in LM9.
   /////  if ( nGen != 10000 ) continue ; // get rid of bad scan points.

          if (    fabs( pointM0 - m0 ) <= deltaM0/2.
               && fabs( pointM12 - m12 ) <= deltaM12/2. ) {

              double nGenPerPoint = 10000 ; // for t1bbbb

             printf("\n\n") ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_error     = %6.1f %%\n", n_sig_error     ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_error      = %6.1f %%\n", n_sb_error      ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_sl_error  = %6.1f %%\n", n_sig_sl_error  ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_sl_error   = %6.1f %%\n", n_sb_sl_error   ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_ldp_error = %6.1f %%\n", n_sig_ldp_error ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_ldp_error  = %6.1f %%\n", n_sb_ldp_error  ) ;
             printf("\n\n") ;

      //--Include the stat error on the efficiency for t1bbbb.
             if ( isT1bbbb ) {

              //-- absolute raw eff
                 float n_sig_raw_eff     = n_sig_raw     / nGenPerPoint ;
                 float n_sb_raw_eff      = n_sb_raw      / nGenPerPoint ;
                 float n_sig_sl_raw_eff  = n_sig_sl_raw  / nGenPerPoint ;
                 float n_sb_sl_raw_eff   = n_sb_sl_raw   / nGenPerPoint ;
                 float n_sig_ldp_raw_eff = n_sig_ldp_raw / nGenPerPoint ;
                 float n_sb_ldp_raw_eff  = n_sb_ldp_raw  / nGenPerPoint ;


               //-- absolute stat err.
                 float n_sig_stat_error     =  sqrt(  n_sig_raw_eff    * ( 1.0 - n_sig_raw_eff     ) / nGenPerPoint ) ;
                 float n_sb_stat_error      =  sqrt(  n_sb_raw_eff     * ( 1.0 - n_sb_raw_eff      ) / nGenPerPoint ) ;
                 float n_sig_sl_stat_error  =  sqrt(  n_sig_sl_raw_eff * ( 1.0 - n_sig_sl_raw_eff  ) / nGenPerPoint ) ;
                 float n_sb_sl_stat_error   =  sqrt(  n_sb_sl_raw_eff  * ( 1.0 - n_sb_sl_raw_eff   ) / nGenPerPoint ) ;
                 float n_sig_ldp_stat_error =  sqrt(  n_sig_ldp_raw_eff* ( 1.0 - n_sig_ldp_raw_eff ) / nGenPerPoint ) ;
                 float n_sb_ldp_stat_error  =  sqrt(  n_sb_ldp_raw_eff * ( 1.0 - n_sb_ldp_raw_eff  ) / nGenPerPoint ) ;

               //-- relative stat err in percent.
                 if ( n_sig_raw_eff     > 0 ) { n_sig_stat_error     = 100.* n_sig_stat_error     / n_sig_raw_eff     ; } else { n_sig_stat_error     = 0. ; }
                 if ( n_sb_raw_eff      > 0 ) { n_sb_stat_error      = 100.* n_sb_stat_error      / n_sb_raw_eff      ; } else { n_sb_stat_error      = 0. ; }
                 if ( n_sig_sl_raw_eff  > 0 ) { n_sig_sl_stat_error  = 100.* n_sig_sl_stat_error  / n_sig_sl_raw_eff  ; } else { n_sig_sl_stat_error  = 0. ; }
                 if ( n_sb_sl_raw_eff   > 0 ) { n_sb_sl_stat_error   = 100.* n_sb_sl_stat_error   / n_sb_sl_raw_eff   ; } else { n_sb_sl_stat_error   = 0. ; }
                 if ( n_sig_ldp_raw_eff > 0 ) { n_sig_ldp_stat_error = 100.* n_sig_ldp_stat_error / n_sig_ldp_raw_eff ; } else { n_sig_ldp_stat_error = 0. ; }
                 if ( n_sb_ldp_raw_eff  > 0 ) { n_sb_ldp_stat_error  = 100.* n_sb_ldp_stat_error  / n_sb_ldp_raw_eff  ; } else { n_sb_ldp_stat_error  = 0. ; }

                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_stat_error     = %6.1f %%\n", n_sig_stat_error     ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_stat_error      = %6.1f %%\n", n_sb_stat_error      ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_sl_stat_error  = %6.1f %%\n", n_sig_sl_stat_error  ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_sl_stat_error   = %6.1f %%\n", n_sb_sl_stat_error   ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_ldp_stat_error = %6.1f %%\n", n_sig_ldp_stat_error ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_ldp_stat_error  = %6.1f %%\n", n_sb_ldp_stat_error  ) ;

               //-- total err in percent.
                 n_sig_error     = sqrt( pow( n_sig_error    , 2) + pow( n_sig_stat_error    , 2) ) ;
                 n_sb_error      = sqrt( pow( n_sb_error     , 2) + pow( n_sb_stat_error     , 2) ) ;
                 n_sig_sl_error  = sqrt( pow( n_sig_sl_error , 2) + pow( n_sig_sl_stat_error , 2) ) ;
                 n_sb_sl_error   = sqrt( pow( n_sb_sl_error  , 2) + pow( n_sb_sl_stat_error  , 2) ) ;
                 n_sig_ldp_error = sqrt( pow( n_sig_ldp_error, 2) + pow( n_sig_ldp_stat_error, 2) ) ;
                 n_sb_ldp_error  = sqrt( pow( n_sb_ldp_error , 2) + pow( n_sb_ldp_stat_error , 2) ) ;

                 printf("\n\n") ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_error     = %6.1f %%\n", n_sig_error     ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_error      = %6.1f %%\n", n_sb_error      ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_sl_error  = %6.1f %%\n", n_sig_sl_error  ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_sl_error   = %6.1f %%\n", n_sb_sl_error   ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_ldp_error = %6.1f %%\n", n_sig_ldp_error ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_ldp_error  = %6.1f %%\n", n_sb_ldp_error  ) ;
                 printf("\n\n") ;

             }

       //--- Not needed with log-normal
        ///  //-- enforce a maximum efficiency uncertainty (to avoid negative scale factors).
        ///  if ( n_sig_error     > 35. ) { n_sig_error     = 35. ; }
        ///  if ( n_sb_error      > 35. ) { n_sb_error      = 35. ; }
        ///  if ( n_sig_sl_error  > 35. ) { n_sig_sl_error  = 35. ; }
        ///  if ( n_sb_sl_error   > 35. ) { n_sb_sl_error   = 35. ; }
        ///  if ( n_sig_ldp_error > 35. ) { n_sig_ldp_error = 35. ; }
        ///  if ( n_sb_ldp_error  > 35. ) { n_sb_ldp_error  = 35. ; }

      //++++++++++++++++++++++++++++++++++++++++++++++++



             double setVal_n_sig(0.) ;
             double setVal_n_sb(0.) ;
             double setVal_n_sig_sl(0.) ;
             double setVal_n_sb_sl(0.) ;
             double setVal_n_sig_ldp(0.) ;
             double setVal_n_sb_ldp(0.) ;


             if ( !isT1bbbb ) {
                //-- tanb40
                setVal_n_sig     = n_sig_raw     * n_sig_correction     ;
                setVal_n_sb      = n_sb_raw      * n_sb_correction      ;
                setVal_n_sig_sl  = n_sig_sl_raw  * n_sig_sl_correction  ;
                setVal_n_sb_sl   = n_sb_sl_raw   * n_sb_sl_correction   ;
                setVal_n_sig_ldp = n_sig_ldp_raw * n_sig_ldp_correction ;
                setVal_n_sb_ldp  = n_sb_ldp_raw  * n_sb_ldp_correction  ;
             } else {
                //-- t1bbbb
                setVal_n_sig     = DataLumi * t1bbbbXsec * (( n_sig_raw      * n_sig_correction     )/ nGenPerPoint ) ;
                setVal_n_sb      = DataLumi * t1bbbbXsec * (( n_sb_raw       * n_sb_correction      )/ nGenPerPoint ) ;
                setVal_n_sig_sl  = DataLumi * t1bbbbXsec * (( n_sig_sl_raw   * n_sig_sl_correction  )/ nGenPerPoint ) ;
                setVal_n_sb_sl   = DataLumi * t1bbbbXsec * (( n_sb_sl_raw    * n_sb_sl_correction   )/ nGenPerPoint ) ;
                setVal_n_sig_ldp = DataLumi * t1bbbbXsec * (( n_sig_ldp_raw  * n_sig_ldp_correction )/ nGenPerPoint ) ;
                setVal_n_sb_ldp  = DataLumi * t1bbbbXsec * (( n_sb_ldp_raw   * n_sb_ldp_correction  )/ nGenPerPoint ) ;
             }

             rv_mu_susymc_sig       -> setVal( setVal_n_sig      ) ;
             //rv_mu_susymc_sb        -> setVal( setVal_n_sb       ) ;
             //rv_mu_susymc_sig_sl    -> setVal( setVal_n_sig_sl   ) ;
             //rv_mu_susymc_sb_sl     -> setVal( setVal_n_sb_sl    ) ;
             rv_mu_susymc_sig_ldp   -> setVal( setVal_n_sig_ldp  ) ;
             rv_mu_susymc_sb_ldp    -> setVal( setVal_n_sb_ldp   ) ;

             rv_width_eff_sf_sig     -> setVal( n_sig_error     / 100. ) ;
             //rv_width_eff_sf_sb      -> setVal( n_sb_error      / 100. ) ;
             //rv_width_eff_sf_sig_sl  -> setVal( n_sig_sl_error  / 100. ) ;
             //rv_width_eff_sf_sb_sl   -> setVal( n_sb_sl_error   / 100. ) ;
             rv_width_eff_sf_sig_ldp -> setVal( n_sig_ldp_error / 100. ) ;
             rv_width_eff_sf_sb_ldp  -> setVal( n_sb_ldp_error  / 100. ) ;

	     eff_susy_sig     =  new RooRealVar( "eff_susy_sig"     , "eff_susy_sig"     , n_sig_correction    );
	     eff_susy_sb      =  new RooRealVar( "eff_susy_sb"      , "eff_susy_sb"      , n_sb_correction     );
	     eff_susy_sig_sl  =  new RooRealVar( "eff_susy_sig_sl"  , "eff_susy_sig_sl"  , n_sig_sl_correction );
	     eff_susy_sb_sl   =  new RooRealVar( "eff_susy_sb_sl"   , "eff_susy_sb_sl"   , n_sb_sl_correction  );
	     eff_susy_sig_ldp =  new RooRealVar( "eff_susy_sig_ldp" , "eff_susy_sig_ldp" , n_sig_ldp_correction);
	     eff_susy_sb_ldp  =  new RooRealVar( "eff_susy_sb_ldp"  , "eff_susy_sb_ldp"  , n_sb_ldp_correction );

	     eff_susy_sig     ->setConstant();
	     eff_susy_sb      ->setConstant();
	     eff_susy_sig_sl  ->setConstant();
	     eff_susy_sb_sl   ->setConstant();
	     eff_susy_sig_ldp ->setConstant();
	     eff_susy_sb_ldp  ->setConstant();


             if ( !isT1bbbb ) {
                printf("\n\n Found point m0 = %4.0f,  m1/2 = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, setVal_n_sig ) ;
             } else {
                printf("\n\n Found point mGluino = %4.0f,  mLSP = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, setVal_n_sig ) ;
             }


             printf("\n\n") ;
             printf(" Setting susy N_sig     to  %7.1f\n", setVal_n_sig       ) ;
             //printf(" Setting susy N_sb      to  %7.1f\n", setVal_n_sb        ) ;
             //printf(" Setting susy N_sig_sl  to  %7.1f\n", setVal_n_sig_sl    ) ;
             //printf(" Setting susy N_sb_sl   to  %7.1f\n", setVal_n_sb_sl     ) ;
             printf(" Setting susy N_sig_ldp to  %7.1f\n", setVal_n_sig_ldp   ) ;
             printf(" Setting susy N_sb_ldp  to  %7.1f\n", setVal_n_sb_ldp    ) ;
             printf("\n\n") ;

             found = true ;

             break ;

          } // point match?

       } // not eof ?

       fclose( infp ) ;

       if ( found ) {
          return true ;
       } else {
          printf("\n\n *** Point not found in scan.\n\n" ) ;
          return false ;
       }

    } // setSusyScanPoint.

