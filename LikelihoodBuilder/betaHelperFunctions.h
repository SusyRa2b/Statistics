#include <iostream>
#include <complex>

#ifndef BETAHELPERFUNCTIONS
#define BETAHELPERFUNCTIONS



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


#endif
