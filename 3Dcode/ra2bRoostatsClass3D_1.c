#include "ra2bRoostatsClass3D_1.h"

#include <iostream>
#include <sstream>
#include <string.h>

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

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"

  using namespace RooFit ;
  using namespace RooStats ;


  //=====================================================================================================


   ra2bRoostatsClass3D_1::ra2bRoostatsClass3D_1() {

      gStyle->SetOptStat(0) ;

     //--- Tell RooFit to shut up about anything less important than an ERROR.
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
      printf("\n\n ==== RooFit output configuration ===============\n") ;
      RooMsgService::instance().Print("v") ;
      printf("\n\n ================================================\n") ;

      varsAtFitVals = false ;
      initialized = false ;

   }

  //===================================================================

   ra2bRoostatsClass3D_1::~ra2bRoostatsClass3D_1() { }



  //===================================================================

    bool ra2bRoostatsClass3D_1::initialize( const char* infile ,
					    const char* inputScanFile,
					    double m0, double m12, bool isT1bbbb, double t1bbbbXsec,
					    const char* inputSusy_deff_dbtageff_file
					    ) {

      printf( "\n\n Opening input file : %s\n\n", infile ) ;

      FILE* infp ;
      if ( (infp=fopen( infile,"r"))==NULL ) {
	printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
	return false ;
      }
      
      sprintf( initializeFile, "%s", infile ) ;


      // define inputs using arrays:
      //
      // 1st index: MET bin
      // 2nd index: HT bin
      // 3rd index: b-multiplicity bin


      // observables:

      int N_0lep[nBinsMET][nBinsHT][nBinsBtag] ;
      int N_1lep[nBinsMET][nBinsHT][nBinsBtag] ;
      int N_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

      int N_Zee[nBinsMET][nBinsHT] ;
      int N_Zmm[nBinsMET][nBinsHT] ;


      // MC inputs:

      float Nttbarsingletopzjetsmc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
      float NWJmc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
      float NZnnmc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

      
      // other parameters:

      float R_lsb[nBinsHT][nBinsBtag] ;
      float R_lsb_err[nBinsHT][nBinsBtag] ;

      float acc_Zee[nBinsMET] ;
      float acc_Zee_err[nBinsMET] ;
      float acc_Zmm[nBinsMET] ;
      float acc_Zmm_err[nBinsMET] ;

      float eff_Zee ;
      float eff_Zee_err ;
      float eff_Zmm ;
      float eff_Zmm_err ;
      
      float knn[nBinsBtag] ;
      float knn_err[nBinsBtag] ;

      float pur_Zee;
      float pur_Zee_err;
      float pur_Zmm;
      float pur_Zmm_err;

      float sf_mc            ;
      float sf_mc_err        ;
      float sf_qcd[nBinsMET][nBinsHT][nBinsBtag] ;
      float sf_qcd_err[nBinsMET][nBinsHT][nBinsBtag] ;
      float sf_ttwj[nBinsMET][nBinsHT][nBinsBtag] ;
      float sf_ttwj_err[nBinsMET][nBinsHT][nBinsBtag] ;
      float sf_ee[nBinsBtag];
      float sf_ee_err[nBinsBtag];       
      float sf_mm[nBinsBtag];
      float sf_mm_err[nBinsBtag];       
      
      float btageff_err ;



      //--- read in description lines.
      printf("\n\n") ;
      char c(0) ;
      while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
      printf("\n\n") ;


      char label[1000] ;
   // char equal;
      
      // read all observables and other parameters
      // first create the TString stubs needed to read the text file:
     
      TString sMbins[nBinsMET];
      TString sHbins[nBinsHT];
      TString sBbins[3] = {"_1b","_2b","_3b"};
      
      for (int i = 0 ; i < nBinsMET ; i++) {
	TString base = "_M";
	stringstream sbin;
	sbin << i+1;
	base += sbin.str();
	sMbins[i] = base;
      }
      
      for (int j = 0 ; j < nBinsHT ; j++) {
	TString base = "_H";
	stringstream sbin;
	sbin << j+1;
	base += sbin.str();
	sHbins[j] = base;
      }


      // 0lep observables

      for ( int i = 0 ; i < nBinsMET ; i++ ) {
	for ( int j = 0 ; j < nBinsHT ; j++ ) {
	  for ( int k = 0 ; k < nBinsBtag ; k++ ) {

	    TString inPar = "N_0lep"+sMbins[i]+sHbins[j]+sBbins[k] ;
            float value ;
	    fscanf( infp, "%s %g", label, &value ) ;
            N_0lep[i][j][k] = TMath::Nint( value ) ;

	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << N_0lep[i][j][k] << endl ;

	  }
	}
      }


      // 1lep observables

      for ( int i = 0 ; i < nBinsMET ; i++ ) {
	for ( int j = 0 ; j < nBinsHT ; j++ ) {
	  for ( int k = 0 ; k < nBinsBtag ; k++ ) {

	    TString inPar = "N_1lep"+sMbins[i]+sHbins[j]+sBbins[k] ;
            float value ;
	    fscanf( infp, "%s %g", label, &value ) ;
            N_1lep[i][j][k] = TMath::Nint( value ) ;

	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << N_1lep[i][j][k] << endl ;

	  }
	}
      }


      // ldp observables

      for ( int i = 0 ; i < nBinsMET ; i++ ) {
	for ( int j = 0 ; j < nBinsHT ; j++ ) {
	  for ( int k = 0 ; k < nBinsBtag ; k++ ) {

	    TString inPar = "N_ldp"+sMbins[i]+sHbins[j]+sBbins[k] ;
            float value ;
	    fscanf( infp, "%s %g", label, &value ) ;
            N_ldp[i][j][k] = TMath::Nint( value ) ;

	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << N_ldp[i][j][k] << endl ;

	  }
	}
      }


      // R_lsb (and its error)

      for ( int j = 0 ; j < nBinsHT ; j++ ) {
	for ( int k = 0 ; k < nBinsBtag ; k++ ) {

	  TString inPar = "R_lsb"+sHbins[j]+sBbins[k] ;
	  fscanf( infp, "%s %g", label, &R_lsb[j][k] ) ;

	  if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	  cout << inPar << " = " << R_lsb[j][k] << endl ;

	  inPar += "_err" ;
	  fscanf( infp, "%s %g", label, &R_lsb_err[j][k] ) ;

	  cout << inPar << " = " << R_lsb_err[j][k] << endl ;

	}
      }


      // Z -> ll VL events:

      for ( int i = 0 ; i < nBinsMET ; i++ ) {
	for ( int j = 0 ; j < nBinsHT ; j++ ) {

	  TString inPar = "N_Zee"+sMbins[i]+sHbins[j] ;
          float value ;
	  fscanf( infp, "%s %g", label, &value ) ;
          N_Zee[i][j] = TMath::Nint( value ) ;

	  if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	  cout << inPar << " = " << N_Zee[i][j] << endl ;	  

	}
      }

      for ( int i = 0 ; i < nBinsMET ; i++ ) {
	for ( int j = 0 ; j < nBinsHT ; j++ ) {

	  TString inPar = "N_Zmm"+sMbins[i]+sHbins[j] ;
          float value ;
	  fscanf( infp, "%s %g", label, &value ) ;
          N_Zmm[i][j] = TMath::Nint( value ) ;

	  if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	  cout << inPar << " = " << N_Zmm[i][j] << endl ;	  

	}
      }


      // Nttbarsingletopzjetsmc_ldp

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {
	    
	    TString inPar = "N_ttbarsingletopzjetsmc_ldp"+sMbins[i]+sHbins[j]+sBbins[k] ;
	    fscanf( infp, "%s %g", label, &Nttbarsingletopzjetsmc_ldp[i][j][k] ) ;
	    
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << Nttbarsingletopzjetsmc_ldp[i][j][k] << endl ;	  

	  }
	}
      }


      // NWJmc_ldp

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {
	    
	    TString inPar = "N_WJmc_ldp"+sMbins[i]+sHbins[j]+sBbins[k] ;
	    fscanf( infp, "%s %g", label, &NWJmc_ldp[i][j][k] ) ;
	    
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << NWJmc_ldp[i][j][k] << endl ;	  

	  }
	}
      }


      // NZnnmc_ldp

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {
	    
	    TString inPar = "N_Znnmc_ldp"+sMbins[i]+sHbins[j]+sBbins[k] ;
	    fscanf( infp, "%s %g", label, &NZnnmc_ldp[i][j][k] ) ;
	    
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << NZnnmc_ldp[i][j][k] << endl ;	  

	  }
	}
      }

      
      // Z -> ll acceptance (+error)
      
      for (int i = 0 ; i < nBinsMET ; i++) {
	    
	TString inPar = "acc_Zee"+sMbins[i] ;
	fscanf( infp, "%s %g", label, &acc_Zee[i] ) ;
	
	if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	cout << inPar << " = " << acc_Zee[i] << endl ;	  
	
	inPar += "_err" ;
	fscanf( infp, "%s %g", label, &acc_Zee_err[i] ) ;
	cout << inPar << " = " << acc_Zee_err[i] << endl ;	  

      }
      
      for (int i = 0 ; i < nBinsMET ; i++) {
	    
	TString inPar = "acc_Zmm"+sMbins[i] ;
	fscanf( infp, "%s %g", label, &acc_Zmm[i] ) ;
	
	if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	cout << inPar << " = " << acc_Zmm[i] << endl ;	  
	
	inPar += "_err" ;
	fscanf( infp, "%s %g", label, &acc_Zmm_err[i] ) ;
	cout << inPar << " = " << acc_Zmm_err[i] << endl ;	  

      }


      // Z -> ll efficiencies

      fscanf( infp, "%s %g", label, &eff_Zee ) ; cout << "Z_ee_eff" << " = " << eff_Zee << endl ;
      fscanf( infp, "%s %g", label, &eff_Zee_err ) ; cout << "Z_ee_eff_err" << " = " << eff_Zee_err << endl ;
      fscanf( infp, "%s %g", label, &eff_Zmm ) ; cout << "Z_mm_eff" << " = " << eff_Zmm << endl ;
      fscanf( infp, "%s %g", label, &eff_Zmm_err ) ; cout << "Z_mm_eff_err" << " = " << eff_Zmm_err << endl ;


      // VL -> nominal scale factors
      
      for (int k = 0 ; k < nBinsBtag ; k++) {
	    
	TString inPar = "knn"+sBbins[k] ;
	fscanf( infp, "%s %g", label, &knn[k] ) ;
	
	if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	cout << inPar << " = " << knn[k] << endl ;	  
	
	inPar += "_err" ;
	fscanf( infp, "%s %g", label, &knn_err[k] ) ;
	cout << inPar << " = " << knn_err[k] << endl ;	  

      }
      

      // Z -> ll purities

      fscanf( infp, "%s %g", label, &pur_Zee ) ; cout << "Z_ee_pur" << " = " << pur_Zee << endl ;
      fscanf( infp, "%s %g", label, &pur_Zee_err ) ; cout << "Z_ee_pur_err" << " = " << pur_Zee_err << endl ;
      fscanf( infp, "%s %g", label, &pur_Zmm ) ; cout << "Z_mm_pur" << " = " << pur_Zmm << endl ;
      fscanf( infp, "%s %g", label, &pur_Zmm_err ) ; cout << "Z_mm_pur_err" << " = " << pur_Zmm_err << endl ;


      // sf MC

      fscanf( infp, "%s %g", label, &sf_mc ) ; cout << "sf_mc" << " = " << sf_mc << endl ;
      fscanf( infp, "%s %g", label, &sf_mc_err ) ; cout << "sf_mc_err" << " = " << sf_mc_err << endl ;


      // sf QCD

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {

	    TString inPar = "sf_qcd"+sMbins[i]+sHbins[j]+sBbins[k] ;
	    fscanf( infp, "%s %g", label, &sf_qcd[i][j][k] ) ;
	
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << sf_qcd[i][j][k] << endl ;	  
	    
	    inPar += "_err" ;
	    fscanf( infp, "%s %g", label, &sf_qcd_err[i][j][k] ) ;
	    cout << inPar << " = " << sf_qcd_err[i][j][k] << endl ;	  

	  }
	}
      }      


      // sf TTWJ

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {

	    TString inPar = "sf_ttwj"+sMbins[i]+sHbins[j]+sBbins[k] ;
	    fscanf( infp, "%s %g", label, &sf_ttwj[i][j][k] ) ;
	
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << sf_ttwj[i][j][k] << endl ;	  
	    
	    inPar += "_err" ;
	    fscanf( infp, "%s %g", label, &sf_ttwj_err[i][j][k] ) ;
	    cout << inPar << " = " << sf_ttwj_err[i][j][k] << endl ;	  

	  }
	}
      }      

      
      // sf Z -> ll
      
      for (int k = 0 ; k < nBinsBtag ; k++) {
	    
	TString inPar = "sf_ee"+sBbins[k] ;
	fscanf( infp, "%s %g", label, &sf_ee[k] ) ;
	
	if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	cout << inPar << " = " << sf_ee[k] << endl ;	  
	
	inPar += "_err" ;
	fscanf( infp, "%s %g", label, &sf_ee_err[k] ) ;
	cout << inPar << " = " << sf_ee_err[k] << endl ;	  

      }
      
      for (int k = 0 ; k < nBinsBtag ; k++) {
	    
	TString inPar = "sf_mm"+sBbins[k] ;
	fscanf( infp, "%s %g", label, &sf_mm[k] ) ;
	
	if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	cout << inPar << " = " << sf_mm[k] << endl ;	  
	
	inPar += "_err" ;
	fscanf( infp, "%s %g", label, &sf_mm_err[k] ) ;
	cout << inPar << " = " << sf_mm_err[k] << endl ;	  

      }

      
      // btag efficiency error

      fscanf( infp, "%s %g", label, &btageff_err ) ; cout << "btageff_err" << " = " << btageff_err << endl ;


      printf("\n Done reading in %s\n\n", infile ) ;
      fclose( infp ) ;
      
      
      //---- calculations for determining initial values for floating parameters.
      
      double initialval_znn_ee[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_znn_mm[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_znn[nBinsMET][nBinsHT][nBinsBtag] ;

      double initialval_qcd_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_qcd[nBinsMET][nBinsHT][nBinsBtag] ;

      double initialval_ttwj_sl[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_ttwj[nBinsMET][nBinsHT][nBinsBtag] ;
      double ratio_ttwj[nBinsMET][nBinsHT][nBinsBtag] ;


      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {     

	    // Z -> invis stuff :
	    
	    initialval_znn_ee[i][j][k] = (N_Zee[i][j]) * ( 5.95 * pur_Zee * knn[k] ) / ( acc_Zee[i] * eff_Zee ) ;
	    initialval_znn_mm[i][j][k] = (N_Zmm[i][j]) * ( 5.95 * pur_Zmm * knn[k] ) / ( acc_Zmm[i] * eff_Zmm ) ;
	    
	    // simple average
	    initialval_znn[i][j][k] = 0.5 * ( initialval_znn_ee[i][j][k] + initialval_znn_mm[i][j][k] ) ;


	    // QCD stuff:

	    initialval_qcd_ldp[i][j][k] = (N_ldp[i][j][k]) - ( Nttbarsingletopzjetsmc_ldp[i][j][k] + NWJmc_ldp[i][j][k] + NZnnmc_ldp[i][j][k] ) ;
	    initialval_qcd[i][j][k] = R_lsb[j][k] * initialval_qcd_ldp[i][j][k] ;

            if ( initialval_qcd_ldp[i][j][k] < 0.) { initialval_qcd_ldp[i][j][k] = 0. ; }
            if ( initialval_qcd[i][j][k]     < 0.) { initialval_qcd[i][j][k]     = 0. ; }
	    
	    // TTWJ stuff. Reference bin is (1,1,1):

	    ratio_ttwj[i][j][k] = (double)N_1lep[i][j][k] / N_1lep[0][0][0] ;
	    
	    initialval_ttwj_sl[i][j][k] = N_1lep[i][j][k] ;
	    initialval_ttwj[i][j][k] = ( N_0lep[0][0][0] - initialval_znn[0][0][0] - initialval_qcd[0][0][0] ) * ratio_ttwj[i][j][k] ;
	    //initialval_ttwj[i][j][k] = N_0lep[i][j][k] - initialval_znn[i][j][k] - initialval_qcd[i][j][k] ; // this is just for debugging !!!

	  }
	}
      }


      // print out initial values:

      printf("\n * means fixed MC value.\n\n\n") ;

      printf("\n\n\n --------- Observables and floating parameter initial values. ------------\n\n") ;
      
      printf("           |  Nobs   |  Model |    PDF        ||  ttwj  |  QCD  |  Znn  |\n") ;
      printf("-----------+---------+--------+---------------++--------+-------+-------+\n") ;

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {     

            double model_0lep = initialval_ttwj[i][j][k] + initialval_qcd[i][j][k] + initialval_znn[i][j][k] ;
            double model_ldp  = Nttbarsingletopzjetsmc_ldp[i][j][k] + initialval_qcd_ldp[i][j][k] + NZnnmc_ldp[i][j][k] ;
            double pdf_0lep = TMath::PoissonI( N_0lep[i][j][k], model_0lep ) ;
            double pdf_ldp  = TMath::PoissonI( N_ldp[i][j][k] , model_ldp  ) ;
            char warning0lep[4] ;
            char warningldp[4] ;
            sprintf( warning0lep,"   ") ;
            sprintf( warningldp,"   ") ;
            if ( pdf_0lep < 0.0001   ) { sprintf( warning0lep, "*  ") ; }
            if ( pdf_0lep < 0.00001  ) { sprintf( warning0lep, "** ") ; }
            if ( pdf_0lep < 0.000001 ) { sprintf( warning0lep, "***") ; }
            if ( pdf_ldp < 0.0001   ) { sprintf( warningldp, "*  ") ; }
            if ( pdf_ldp < 0.00001  ) { sprintf( warningldp, "** ") ; }
            if ( pdf_ldp < 0.000001 ) { sprintf( warningldp, "***") ; }
	    cout << " MET bin " << i+1 << ", HT bin " << j+1 << ", Btag bin " << k+1 << endl;
	    printf(" 0-lep     | %5d   | %6.1f | %8.6f %4s ||  %5.1f | %5.1f | %5.1f |\n", N_0lep[i][j][k], model_0lep, pdf_0lep, warning0lep,
                                                                                       initialval_ttwj[i][j][k], initialval_qcd[i][j][k], initialval_znn[i][j][k] ) ;
	    printf(" ldp       | %5d   | %6.1f | %8.6f %4s || *%5.1f | %5.1f |*%5.1f |\n", N_ldp[i][j][k], model_ldp, pdf_ldp, warningldp,
                                                                                       Nttbarsingletopzjetsmc_ldp[i][j][k], initialval_qcd_ldp[i][j][k], NZnnmc_ldp[i][j][k] ) ;
            printf("-----------+---------+--------+---------------++--------+-------+-------+\n") ;

	  }
	}
      }

      printf("\n * means fixed MC value.\n\n\n") ;


      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining observables and parameters of the likelihood.\n" ) ;


      rv_mu_susy_M1_H1_1b = new RooRealVar( "mu_susy_M1_H1_1b", "mu_susy_M1_H1_1b", 0., 100000. ) ;
      rv_mu_susy_M1_H1_1b->setVal( 1. ) ;
      rv_mu_susy_M1_H1_1b->setConstant( kTRUE ) ;

            
      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {     

	    TString zlString  = "N_0lep";
	    TString slString  = "N_1lep";
	    TString ldpString = "N_ldp";

	    zlString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    slString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    ldpString += sMbins[i]+sHbins[j]+sBbins[k] ;

	    rv_0lep[i][j][k] = new RooRealVar( zlString,  zlString,  0., 1000000.);
	    rv_1lep[i][j][k] = new RooRealVar( slString,  slString,  0., 1000000.);
	    rv_ldp[i][j][k]  = new RooRealVar( ldpString, ldpString, 0., 1000000.);

	    // set values:
	    rv_0lep[i][j][k]->setVal( N_0lep[i][j][k] );
	    rv_1lep[i][j][k]->setVal( N_1lep[i][j][k] );
	    rv_ldp[i][j][k]->setVal( N_ldp[i][j][k] );

	    
	    // parameters of the likelihood:
	    TString muTtString     = "mu_ttwj";
	    TString muTtSlString   = "mu_ttwj_sl";
	    TString muQcdString    = "mu_qcd";
	    TString muQcdLdpString = "mu_qcd_ldp";
	    TString muZnnString    = "mu_znn";

	    TString muSusMcString    = "mu_susymc";
	    TString muSusSlMcString  = "mu_susymc_sl";
	    TString muSusLdpMcString = "mu_susymc_ldp";

	    TString muTtLdpString  = "mu_ttbarsingletopzjetsmc";
	    TString muWjLdpString  = "mu_WJmc";
	    TString muZnnLdpString = "mu_Znnmc";

	    TString MEffSfString     = "mean_eff_sf";
	    TString WEffSfString     = "width_eff_sf";
	    TString dEffdBtString    = "deff_dbtageff";

	    TString MEffSfSlString   = "mean_eff_sf_sl";
	    TString WEffSfSlString   = "width_eff_sf_sl";
	    TString dEffdBtSlString  = "deff_dbtageff_sl";

	    TString MEffSfLdpString  = "mean_eff_sf_ldp";
	    TString WEffSfLdpString  = "width_eff_sf_ldp";
	    TString dEffdBtLdpString = "deff_dbtageff_ldp";


	    muTtString     += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muTtSlString   += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muQcdString    += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muQcdLdpString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muZnnString    += sMbins[i]+sHbins[j]+sBbins[k] ;

	    muSusMcString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muSusSlMcString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muSusLdpMcString  += sMbins[i]+sHbins[j]+sBbins[k] ;

	    muTtLdpString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muWjLdpString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muZnnLdpString += sMbins[i]+sHbins[j]+sBbins[k] ;

	    MEffSfString   += sMbins[i]+sHbins[j]+sBbins[k] ;
	    WEffSfString   += sMbins[i]+sHbins[j]+sBbins[k] ;
	    dEffdBtString  += sMbins[i]+sHbins[j]+sBbins[k] ;


	    rrv_mu_ttwj[i][j][k] = new RooRealVar( muTtString, muTtString, 0., 100000. ) ;
	    rv_mu_ttwj[i][j][k] = rrv_mu_ttwj[i][j][k];
	    rrv_mu_ttwj[i][j][k]->setVal( initialval_ttwj[i][j][k] ) ;         // this is a starting value only

	    rv_mu_ttwj_sl[i][j][k] = new RooRealVar( muTtSlString, muTtSlString, 0., 100000. ) ;
	    rv_mu_ttwj_sl[i][j][k]->setVal( initialval_ttwj_sl[i][j][k] ) ;    // this is a starting value only

	    rrv_mu_qcd[i][j][k] = new RooRealVar( muQcdString, muQcdString, 0., 100000. ) ;
	    rv_mu_qcd[i][j][k] = rrv_mu_qcd[i][j][k];
	    rrv_mu_qcd[i][j][k]->setVal( initialval_qcd[i][j][k] ) ;           // this is a starting value only

	    rrv_mu_qcd_ldp[i][j][k] = new RooRealVar( muQcdLdpString, muQcdLdpString, 0., 100000. ) ;
	    rv_mu_qcd_ldp[i][j][k] = rrv_mu_qcd_ldp[i][j][k];
	    rrv_mu_qcd_ldp[i][j][k]->setVal( initialval_qcd_ldp[i][j][k] ) ;   // this is a starting value only

	    rrv_mu_znn[i][j][k] = new RooRealVar( muZnnString, muZnnString, 0., 100000. ) ;
	    rv_mu_znn[i][j][k] = rrv_mu_znn[i][j][k];
	    rrv_mu_znn[i][j][k]->setVal( initialval_znn[i][j][k] ) ;           // this is a starting value only


	    // MC inputs

	    rv_mu_susymc[i][j][k] = new RooRealVar( muSusMcString, muSusMcString, 0., 100000. ) ;
	    rv_mu_susymc[i][j][k]->setVal( 1. ) ;
	    rv_mu_susymc[i][j][k]->setConstant( kTRUE ) ;

	    rv_mu_susymc_sl[i][j][k] = new RooRealVar( muSusSlMcString, muSusSlMcString, 0., 100000. ) ;
	    rv_mu_susymc_sl[i][j][k]->setVal( 1. ) ;
	    rv_mu_susymc_sl[i][j][k]->setConstant( kTRUE ) ;

	    rv_mu_susymc_ldp[i][j][k] = new RooRealVar( muSusLdpMcString, muSusLdpMcString, 0., 100000. ) ;
	    rv_mu_susymc_ldp[i][j][k]->setVal( 1. ) ;
	    rv_mu_susymc_ldp[i][j][k]->setConstant( kTRUE ) ;

	    rv_mu_ttbarsingletopzjetsmc_ldp[i][j][k] = new RooRealVar( muTtLdpString, muTtLdpString, 0., 100000. ) ;
	    rv_mu_ttbarsingletopzjetsmc_ldp[i][j][k]->setVal( Nttbarsingletopzjetsmc_ldp[i][j][k] ) ;
	    rv_mu_ttbarsingletopzjetsmc_ldp[i][j][k]->setConstant( kTRUE ) ;

	    rv_mu_WJmc_ldp[i][j][k] = new RooRealVar( muWjLdpString, muWjLdpString, 0., 100000. ) ;
	    rv_mu_WJmc_ldp[i][j][k]->setVal( NWJmc_ldp[i][j][k] ) ;
	    rv_mu_WJmc_ldp[i][j][k]->setConstant( kTRUE ) ;

	    rv_mu_Znnmc_ldp[i][j][k] = new RooRealVar( muZnnLdpString, muZnnLdpString, 0., 100000. ) ;
	    rv_mu_Znnmc_ldp[i][j][k]->setVal( NZnnmc_ldp[i][j][k] ) ;
	    rv_mu_Znnmc_ldp[i][j][k]->setConstant( kTRUE ) ;

	    
	    // gaussian constraints

	    rv_mean_eff_sf[i][j][k] = new RooRealVar( MEffSfString, MEffSfString, 0., 10. );
	    rv_mean_eff_sf[i][j][k]->setVal( 1.0 ) ;
	    rv_mean_eff_sf[i][j][k]->setConstant( kTRUE ) ;

	    rv_width_eff_sf[i][j][k] = new RooRealVar( WEffSfString, WEffSfString, 0., 10. );
	    rv_width_eff_sf[i][j][k]->setVal( 0.15 ) ;            // arbitrarily set the uncertainty to 15%
	    rv_width_eff_sf[i][j][k]->setConstant( kTRUE ) ;   


	    rv_mean_eff_sf_sl[i][j][k] = new RooRealVar( MEffSfSlString, MEffSfSlString, 0., 10. );
	    rv_mean_eff_sf_sl[i][j][k]->setVal( 1.0 ) ;
	    rv_mean_eff_sf_sl[i][j][k]->setConstant( kTRUE ) ;

	    rv_width_eff_sf_sl[i][j][k] = new RooRealVar( WEffSfSlString, WEffSfSlString, 0., 10. );
	    rv_width_eff_sf_sl[i][j][k]->setVal( 0.15 ) ;            // arbitrarily set the uncertainty to 15%
	    rv_width_eff_sf_sl[i][j][k]->setConstant( kTRUE ) ;   


	    rv_mean_eff_sf_ldp[i][j][k] = new RooRealVar( MEffSfLdpString, MEffSfLdpString, 0., 10. );
	    rv_mean_eff_sf_ldp[i][j][k]->setVal( 1.0 ) ;
	    rv_mean_eff_sf_ldp[i][j][k]->setConstant( kTRUE ) ;

	    rv_width_eff_sf_ldp[i][j][k] = new RooRealVar( WEffSfLdpString, WEffSfLdpString, 0., 10. );
	    rv_width_eff_sf_ldp[i][j][k]->setVal( 0.15 ) ;            // arbitrarily set the uncertainty to 15%
	    rv_width_eff_sf_ldp[i][j][k]->setConstant( kTRUE ) ;   

	    
	    // btag efficiency derivatives

	    rv_deff_dbtageff[i][j][k] = new RooRealVar( dEffdBtString,  dEffdBtString, -10., 10. );
	    rv_deff_dbtageff[i][j][k]->setVal( 0.1 ) ;
	    rv_deff_dbtageff[i][j][k]->setConstant( kTRUE ) ;

	    rv_deff_dbtageff_sl[i][j][k] = new RooRealVar( dEffdBtSlString,  dEffdBtSlString, -10., 10. );
	    rv_deff_dbtageff_sl[i][j][k]->setVal( 0.1 ) ;
	    rv_deff_dbtageff_sl[i][j][k]->setConstant( kTRUE ) ;

	    rv_deff_dbtageff_ldp[i][j][k] = new RooRealVar( dEffdBtLdpString,  dEffdBtLdpString, -10., 10. );
	    rv_deff_dbtageff_ldp[i][j][k]->setVal( 0.1 ) ;
	    rv_deff_dbtageff_ldp[i][j][k]->setConstant( kTRUE ) ;

	  }
	  
	  TString zeeString = "N_Zee";
	  TString zmmString = "N_Zmm";
	  zeeString += sMbins[i]+sHbins[j];
	  zmmString += sMbins[i]+sHbins[j];

	  rv_Zee[i][j] = new RooRealVar( zeeString, zeeString, 0., 100000.);
	  rv_Zmm[i][j] = new RooRealVar( zmmString, zmmString, 0., 100000.);

	  // set values:
	  rv_Zee[i][j]->setVal( N_Zee[i][j] );
	  rv_Zmm[i][j]->setVal( N_Zmm[i][j] );

	}
      }


      RooRealVar* rv_btageff_err = new RooRealVar( "btageff_err", "btageff_err", 0., 10. ) ;
      rv_btageff_err->setVal( btageff_err ) ;
      rv_btageff_err->setConstant( kTRUE ) ;


      //++++++++++

      bool ssspOk = setSusyScanPoint( inputScanFile,  m0,  m12,  isT1bbbb,  t1bbbbXsec, inputSusy_deff_dbtageff_file ) ;

      if ( !ssspOk ) {
         printf("\n\n\n *** setSusyScanPoint failed.  I quit.\n\n\n") ;
         return false ;
      }
      
      //++++++++++

      cout << "\n\n Back from setSusyScanPoint.  Now defining parameters.\n\n" << flush ;

      //--- Systematics and other nuisance parameters
      // THIS WILL NEED TO BE CAREFULLY REVISED !!!!
      // STILL USING THE LOG-NORMALS TO GET STARTED

      char formula[1024];
      RooArgSet globalObservables ("globalObservables");
      RooArgSet allNuisances ("allNuisances");
      RooArgSet allNuisancePdfs ("allNuisancePdfs");


      // Rlsb_passfail

      RooRealVar* Rlsb_passfail_prim[nBinsHT][nBinsBtag];
      RooRealVar* Rlsb_passfail_nom[nBinsHT][nBinsBtag];
      RooGaussian* pdf_Rlsb_passfail[nBinsHT][nBinsBtag];
      RooFormulaVar* fv_Rlsb_passfail[nBinsHT][nBinsBtag];

      for (int j = 0 ; j < nBinsHT ; j++) {
	for (int k = 0 ; k < nBinsBtag ; k++) {     

	  TString RlsbPfString     = "Rlsb_passfail";
	  TString RlsbPfPrimString = "Rlsb_passfail_prim";
	  TString RlsbPfNomString  = "Rlsb_passfail_nom";
	  TString RlsbPfPdfString  = "pdf_Rlsb_passfail";
	  
	  RlsbPfString     += sHbins[j]+sBbins[k] ;
	  RlsbPfPrimString += sHbins[j]+sBbins[k] ;
	  RlsbPfNomString  += sHbins[j]+sBbins[k] ;
	  RlsbPfPdfString  += sHbins[j]+sBbins[k] ;

	  Rlsb_passfail_prim[j][k] = new RooRealVar( RlsbPfPrimString, RlsbPfPrimString, 0., -5., 5.);
	  Rlsb_passfail_nom[j][k] = new RooRealVar( RlsbPfNomString, RlsbPfNomString, 0., -5., 5.);
	  pdf_Rlsb_passfail[j][k] = new RooGaussian( RlsbPfPdfString, RlsbPfPdfString, *Rlsb_passfail_prim[j][k], *Rlsb_passfail_nom[j][k], RooConst(1) );
	  sprintf (formula, "%f*pow(%f,@0)", R_lsb[j][k], exp(R_lsb_err[j][k]/R_lsb[j][k]));
	  fv_Rlsb_passfail[j][k] = new RooFormulaVar(RlsbPfString, formula, RooArgList(*Rlsb_passfail_prim[j][k]));
	  Rlsb_passfail_nom[j][k]->setConstant();
	  globalObservables.add (*Rlsb_passfail_nom[j][k]);
	  allNuisances.add (*Rlsb_passfail_prim[j][k]);
	  allNuisancePdfs.add (*pdf_Rlsb_passfail[j][k]);

	}
      }


      // Z -> ll acceptance

      RooRealVar* acc_Zee_prim[nBinsMET];
      RooRealVar* acc_Zee_nom[nBinsMET];
      RooGaussian* pdf_acc_Zee[nBinsMET];
      RooFormulaVar* fv_acc_Zee[nBinsMET];

      RooRealVar* acc_Zmm_prim[nBinsMET];
      RooRealVar* acc_Zmm_nom[nBinsMET];
      RooGaussian* pdf_acc_Zmm[nBinsMET];
      RooFormulaVar* fv_acc_Zmm[nBinsMET];

      for (int i = 0 ; i < nBinsMET ; i++) {

	TString acc_ZeeString     = "acc_Zee";
	TString acc_ZeePrimString = "acc_Zee_prim";
	TString acc_ZeeNomString  = "acc_Zee_nom";
	TString acc_ZeePdfString  = "pdf_acc_Zee";
	
	acc_ZeeString     += sMbins[i] ;
	acc_ZeePrimString += sMbins[i] ;
	acc_ZeeNomString  += sMbins[i] ;
	acc_ZeePdfString  += sMbins[i] ;
	
	acc_Zee_prim[i] = new RooRealVar( acc_ZeePrimString, acc_ZeePrimString, 0., -5., 5.);
	acc_Zee_nom[i] = new RooRealVar( acc_ZeeNomString, acc_ZeeNomString, 0., -5., 5.);
	pdf_acc_Zee[i] = new RooGaussian( acc_ZeePdfString, acc_ZeePdfString, *acc_Zee_prim[i], *acc_Zee_nom[i], RooConst(1) );
	sprintf (formula, "%f*pow(%f,@0)", acc_Zee[i], exp(acc_Zee_err[i]/acc_Zee[i]));
	fv_acc_Zee[i] = new RooFormulaVar(acc_ZeeString, formula, RooArgList(*acc_Zee_prim[i]));
	acc_Zee_nom[i]->setConstant();
	globalObservables.add (*acc_Zee_nom[i]);
	allNuisances.add (*acc_Zee_prim[i]);
	allNuisancePdfs.add (*pdf_acc_Zee[i]);
	
	TString acc_ZmmString     = "acc_Zmm";
	TString acc_ZmmPrimString = "acc_Zmm_prim";
	TString acc_ZmmNomString  = "acc_Zmm_nom";
	TString acc_ZmmPdfString  = "pdf_acc_Zmm";
	
	acc_ZmmString     += sMbins[i] ;
	acc_ZmmPrimString += sMbins[i] ;
	acc_ZmmNomString  += sMbins[i] ;
	acc_ZmmPdfString  += sMbins[i] ;

	acc_Zmm_prim[i] = new RooRealVar( acc_ZmmPrimString, acc_ZmmPrimString, 0., -5., 5.);
	acc_Zmm_nom[i] = new RooRealVar( acc_ZmmNomString, acc_ZmmNomString, 0., -5., 5.);
	pdf_acc_Zmm[i] = new RooGaussian( acc_ZmmPdfString, acc_ZmmPdfString, *acc_Zmm_prim[i], *acc_Zmm_nom[i], RooConst(1) );
	sprintf (formula, "%f*pow(%f,@0)", acc_Zmm[i], exp(acc_Zmm_err[i]/acc_Zmm[i]));
	fv_acc_Zmm[i] = new RooFormulaVar(acc_ZmmString, formula, RooArgList(*acc_Zmm_prim[i]));
	acc_Zmm_nom[i]->setConstant();
	globalObservables.add (*acc_Zmm_nom[i]);
	allNuisances.add (*acc_Zmm_prim[i]);
	allNuisancePdfs.add (*pdf_acc_Zmm[i]);
	
      }


      // Z -> ll efficiencies
	
      RooRealVar *eff_Zee_prim = new RooRealVar( "eff_Zee_prim", "eff_Zee_prim", 0., -5., 5.);
      RooRealVar *eff_Zee_nom = new RooRealVar( "eff_Zee_nom", "eff_Zee_nom", 0., -5., 5.);
      RooGaussian *pdf_eff_Zee = new RooGaussian( "pdf_eff_Zee", "pdf_eff_Zee", *eff_Zee_prim, *eff_Zee_nom, RooConst(1) );
      sprintf (formula, "%f*pow(%f,@0)", eff_Zee, exp(eff_Zee_err/eff_Zee));
      RooFormulaVar* fv_eff_Zee = new RooFormulaVar("fv_eff_Zee", formula, RooArgList(*eff_Zee_prim));
      eff_Zee_nom->setConstant();
      globalObservables.add (*eff_Zee_nom);
      allNuisances.add (*eff_Zee_prim);
      allNuisancePdfs.add (*pdf_eff_Zee);
	
      RooRealVar *eff_Zmm_prim = new RooRealVar( "eff_Zmm_prim", "eff_Zmm_prim", 0., -5., 5.);
      RooRealVar *eff_Zmm_nom = new RooRealVar( "eff_Zmm_nom", "eff_Zmm_nom", 0., -5., 5.);
      RooGaussian *pdf_eff_Zmm = new RooGaussian( "pdf_eff_Zmm", "pdf_eff_Zmm", *eff_Zmm_prim, *eff_Zmm_nom, RooConst(1) );
      sprintf (formula, "%f*pow(%f,@0)", eff_Zmm, exp(eff_Zmm_err/eff_Zmm));
      RooFormulaVar* fv_eff_Zmm = new RooFormulaVar("fv_eff_Zmm", formula, RooArgList(*eff_Zmm_prim));
      eff_Zmm_nom->setConstant();
      globalObservables.add (*eff_Zmm_nom);
      allNuisances.add (*eff_Zmm_prim);
      allNuisancePdfs.add (*pdf_eff_Zmm);


      // Z -> ll purities
      
      RooRealVar *pur_Zee_prim = new RooRealVar( "pur_Zee_prim", "pur_Zee_prim", 0., -5., 5.);
      RooRealVar *pur_Zee_nom = new RooRealVar( "pur_Zee_nom", "pur_Zee_nom", 0., -5., 5.);
      RooGaussian *pdf_pur_Zee = new RooGaussian( "pdf_pur_Zee", "pdf_pur_Zee", *pur_Zee_prim, *pur_Zee_nom, RooConst(1) );
      sprintf (formula, "%f*pow(%f,@0)", pur_Zee, exp(pur_Zee_err/pur_Zee));
      RooFormulaVar* fv_pur_Zee = new RooFormulaVar("fv_pur_Zee", formula, RooArgList(*pur_Zee_prim));
      pur_Zee_nom->setConstant();
      globalObservables.add (*pur_Zee_nom);
      allNuisances.add (*pur_Zee_prim);
      allNuisancePdfs.add (*pdf_pur_Zee);
	
      RooRealVar *pur_Zmm_prim = new RooRealVar( "pur_Zmm_prim", "pur_Zmm_prim", 0., -5., 5.);
      RooRealVar *pur_Zmm_nom = new RooRealVar( "pur_Zmm_nom", "pur_Zmm_nom", 0., -5., 5.);
      RooGaussian *pdf_pur_Zmm = new RooGaussian( "pdf_pur_Zmm", "pdf_pur_Zmm", *pur_Zmm_prim, *pur_Zmm_nom, RooConst(1) );
      sprintf (formula, "%f*pow(%f,@0)", pur_Zmm, exp(pur_Zmm_err/pur_Zmm));
      RooFormulaVar* fv_pur_Zmm = new RooFormulaVar("fv_pur_Zmm", formula, RooArgList(*pur_Zmm_prim));
      pur_Zmm_nom->setConstant();
      globalObservables.add (*pur_Zmm_nom);
      allNuisances.add (*pur_Zmm_prim);
      allNuisancePdfs.add (*pdf_pur_Zmm);


      // VL -> nominal scale factors

      RooRealVar* knn_prim[nBinsBtag];
      RooRealVar* knn_nom[nBinsBtag];
      RooGaussian* pdf_knn[nBinsBtag];
      RooFormulaVar* fv_knn[nBinsBtag];

      for (int k = 0 ; k < nBinsBtag ; k++) {

	TString knnString     = "knn";
	TString knnPrimString = "knn_prim";
	TString knnNomString  = "knn_nom";
	TString knnPdfString  = "pdf_knn";
	
	knnString     += sBbins[k] ;
	knnPrimString += sBbins[k] ;
	knnNomString  += sBbins[k] ;
	knnPdfString  += sBbins[k] ;
	
	knn_prim[k] = new RooRealVar( knnPrimString, knnPrimString, 0., -5., 5.);
	knn_nom[k] = new RooRealVar( knnNomString, knnNomString, 0., -5., 5.);
	pdf_knn[k] = new RooGaussian( knnPdfString, knnPdfString, *knn_prim[k], *knn_nom[k], RooConst(1) );
	sprintf (formula, "%f*pow(%f,@0)", knn[k], exp(knn_err[k]/knn[k]));
	fv_knn[k] = new RooFormulaVar(knnString, formula, RooArgList(*knn_prim[k]));
	knn_nom[k]->setConstant();
	globalObservables.add (*knn_nom[k]);
	allNuisances.add (*knn_prim[k]);
	allNuisancePdfs.add (*pdf_knn[k]);

      }


      // sf_ee and sf_mm derived from a common underlying gaussian

      RooRealVar* sf_ll_prim = new RooRealVar( "sf_ll_prim", "sf_ll_prim", 0, -5, 5);
      RooRealVar* sf_ll_nom = new RooRealVar( "sf_ll_nom", "sf_ll_nom", 0, -5, 5);
      RooGaussian* pdf_sf_ll = new RooGaussian("pdf_sf_ll" , "pdf_sf_ll", *sf_ll_prim, *sf_ll_nom, RooConst(1));
      sf_ll_nom->setConstant();
      globalObservables.add (*sf_ll_nom);
      allNuisances.add (*sf_ll_prim);
      allNuisancePdfs.add (*pdf_sf_ll);

      RooFormulaVar* fv_sf_ee[nBinsBtag];
      RooFormulaVar* fv_sf_mm[nBinsBtag];

      for (int k = 0 ; k < nBinsBtag ; k++) {

	TString sf_eeString = "sf_ee";
	TString sf_mmString = "sf_mm";
	
	sf_eeString += sBbins[k] ;
	sf_mmString += sBbins[k] ;

	sprintf (formula, "%f*pow(%f,@0)", sf_ee[k], exp(sf_ee_err[k]/sf_ee[k]));
	fv_sf_ee[k] = new RooFormulaVar( sf_eeString, formula, RooArgList(*sf_ll_prim));

	sprintf (formula, "%f*pow(%f,@0)", sf_mm[k], exp(sf_mm_err[k]/sf_mm[k]));
	fv_sf_mm[k] = new RooFormulaVar( sf_mmString, formula, RooArgList(*sf_ll_prim));

      }


      // MC scale factor

      RooRealVar* sf_mc_prim = new RooRealVar( "sf_mc_prim", "sf_mc_prim", 0, -5, 5);
      RooRealVar* sf_mc_nom = new RooRealVar( "sf_mc_nom", "sf_mc_nom", 0, -5, 5);
      RooGaussian* pdf_sf_mc = new RooGaussian("pdf_sf_mc" , "pdf_sf_mc", *sf_mc_prim, *sf_mc_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_mc, exp(sf_mc_err/sf_mc));
      RooFormulaVar* fv_sf_mc = new RooFormulaVar("sf_mc", formula, RooArgList(*sf_mc_prim));
      sf_mc_nom->setConstant();
      globalObservables.add (*sf_mc_nom);
      allNuisances.add (*sf_mc_prim);
      allNuisancePdfs.add (*pdf_sf_mc);


      // QCD and TTWJ scale factors

      RooRealVar* sf_qcd_prim[nBinsMET][nBinsHT][nBinsBtag];
      RooRealVar* sf_qcd_nom[nBinsMET][nBinsHT][nBinsBtag];
      RooGaussian* pdf_sf_qcd[nBinsMET][nBinsHT][nBinsBtag];
      RooFormulaVar* fv_sf_qcd[nBinsMET][nBinsHT][nBinsBtag];

      RooRealVar* sf_ttwj_prim[nBinsMET][nBinsHT][nBinsBtag];
      RooRealVar* sf_ttwj_nom[nBinsMET][nBinsHT][nBinsBtag];
      RooGaussian* pdf_sf_ttwj[nBinsMET][nBinsHT][nBinsBtag];
      RooFormulaVar* fv_sf_ttwj[nBinsMET][nBinsHT][nBinsBtag];

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {     

	    TString sfQcdString     = "sf_qcd";
	    TString sfQcdPrimString = "sf_qcd_prim";
	    TString sfQcdNomString  = "sf_qcd_nom";
	    TString sfQcdPdfString  = "pdf_sf_qcd";
	    
	    sfQcdString     += sMbins[i]+sHbins[j]+sBbins[k] ;
	    sfQcdPrimString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    sfQcdNomString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    sfQcdPdfString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    
	    sf_qcd_prim[i][j][k] = new RooRealVar( sfQcdPrimString, sfQcdPrimString, 0., -5., 5.);
	    sf_qcd_nom[i][j][k] = new RooRealVar( sfQcdNomString, sfQcdNomString, 0., -5., 5.);
	    pdf_sf_qcd[i][j][k] = new RooGaussian( sfQcdPdfString, sfQcdPdfString, *sf_qcd_prim[i][j][k], *sf_qcd_nom[i][j][k], RooConst(1) );
	    sprintf (formula, "%f*pow(%f,@0)", sf_qcd[i][j][k], exp(sf_qcd_err[i][j][k]/sf_qcd[i][j][k]));
	    fv_sf_qcd[i][j][k] = new RooFormulaVar(sfQcdString, formula, RooArgList(*sf_qcd_prim[i][j][k]));
	    sf_qcd_nom[i][j][k]->setConstant();
	    globalObservables.add (*sf_qcd_nom[i][j][k]);
	    allNuisances.add (*sf_qcd_prim[i][j][k]);
	    allNuisancePdfs.add (*pdf_sf_qcd[i][j][k]);


	    TString sfTtwjString     = "sf_ttwj";
	    TString sfTtwjPrimString = "sf_ttwj_prim";
	    TString sfTtwjNomString  = "sf_ttwj_nom";
	    TString sfTtwjPdfString  = "pdf_sf_ttwj";
	    
	    sfTtwjString     += sMbins[i]+sHbins[j]+sBbins[k] ;
	    sfTtwjPrimString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    sfTtwjNomString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    sfTtwjPdfString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    
	    sf_ttwj_prim[i][j][k] = new RooRealVar( sfTtwjPrimString, sfTtwjPrimString, 0., -5., 5.);
	    sf_ttwj_nom[i][j][k] = new RooRealVar( sfTtwjNomString, sfTtwjNomString, 0., -5., 5.);
	    pdf_sf_ttwj[i][j][k] = new RooGaussian( sfTtwjPdfString, sfTtwjPdfString, *sf_ttwj_prim[i][j][k], *sf_ttwj_nom[i][j][k], RooConst(1) );
	    sprintf (formula, "%f*pow(%f,@0)", sf_ttwj[i][j][k], exp(sf_ttwj_err[i][j][k]/sf_ttwj[i][j][k]));
	    fv_sf_ttwj[i][j][k] = new RooFormulaVar(sfTtwjString, formula, RooArgList(*sf_ttwj_prim[i][j][k]));
	    sf_ttwj_nom[i][j][k]->setConstant();
	    globalObservables.add (*sf_ttwj_nom[i][j][k]);
	    allNuisances.add (*sf_ttwj_prim[i][j][k]);
	    allNuisancePdfs.add (*pdf_sf_ttwj[i][j][k]);

	  }
	}
      }


      RooRealVar* eff_sf_prim = new RooRealVar( "eff_sf_prim", "eff_sf_prim", 0, -5, 5);
      RooRealVar* eff_sf_nom = new RooRealVar( "eff_sf_nom", "eff_sf_nom", 0, -5, 5);
      RooGaussian* pdf_eff_sf = new RooGaussian( "pdf_eff_sf" , "pdf_eff_sf", *eff_sf_prim, *eff_sf_nom, RooConst(1));
      eff_sf_nom->setConstant();
      globalObservables.add (*eff_sf_nom);
      allNuisances.add (*eff_sf_prim);
      allNuisancePdfs.add (*pdf_eff_sf);


      RooRealVar* btageff_sf_prim = new RooRealVar( "btageff_sf_prim", "btageff_sf_prim", 0, -5, 5);
      RooRealVar* btageff_sf_nom = new RooRealVar( "btageff_sf_nom", "btageff_sf_nom", 0, -5, 5);
      RooGaussian* pdf_btageff_sf = new RooGaussian( "pdf_btageff_sf", "pdf_btageff_sf", *btageff_sf_prim, *btageff_sf_nom, RooConst(1));
      btageff_sf_nom->setConstant();
      globalObservables.add (*btageff_sf_nom);
      allNuisances.add (*btageff_sf_prim);
      allNuisancePdfs.add (*pdf_btageff_sf);


      //+++++++++++++++++ Relationships between parameters ++++++++++++++++++++++++++++++++++++++++++++
      
      printf(" --- Defining relationships between parameters.\n" ) ;
      cout << flush ;

      RooArgSet pdflist ;

      RooRealVar* rv_znnoverll_bfratio = new RooRealVar( "znnoverll_bfratio", "znnoverll_bfratio", 0.1, 10. ) ;
      rv_znnoverll_bfratio -> setVal( 5.95 ) ;
      rv_znnoverll_bfratio -> setConstant( kTRUE ) ;

      RooRealVar* rv_dataoverll_lumiratio = new RooRealVar( "dataoverll_lumiratio", "dataoverll_lumiratio", 0.1, 10.0 ) ;
      rv_dataoverll_lumiratio  -> setVal( 1.0 ) ;         // might have to change it later on
      rv_dataoverll_lumiratio  -> setConstant( kTRUE ) ;


      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {     

	    //---- TTWJ

	    TString ttwjString   = "mu_ttwj";
	    ttwjString   += sMbins[i]+sHbins[j]+sBbins[k] ;

	    // for the normalization, float bin (1,1,1) and define the others using the SL ratios
	    if ( !(i == 0 && j == 0 && k == 0) ) {
	       
	      TString rfvString =  " @0 * @1 * ( @2 / @3 )" ;
	      
	      rfv_mu_ttwj[i][j][k] = new RooFormulaVar( ttwjString, rfvString, 
							RooArgSet( *rv_mu_ttwj[0][0][0], *fv_sf_ttwj[i][j][k], *rv_mu_ttwj_sl[i][j][k], *rv_mu_ttwj_sl[0][0][0] )) ;
	      
	      rv_mu_ttwj[i][j][k] = rfv_mu_ttwj[i][j][k] ;

	    }
	    
	    TString ttwjLdpString  = "mu_ttwj_ldp" ;
	    ttwjLdpString  += sMbins[i]+sHbins[j]+sBbins[k] ;

	    TString rfvString = "@0 + @1" ;

	    rv_mu_ttwj_ldp[i][j][k] = new RooFormulaVar( ttwjLdpString, rfvString, 
							 RooArgSet( *rv_mu_ttbarsingletopzjetsmc_ldp[i][j][k], *rv_mu_WJmc_ldp[i][j][k] )) ;
	    

	    //---- QCD

	    TString qcdString    = "mu_qcd" ;
	    qcdString    += sMbins[i]+sHbins[j]+sBbins[k] ;

	    TString rfvQcdString = "@0 * @1 * @2" ;

	    rfv_mu_qcd[i][j][k] = new RooFormulaVar( qcdString, rfvQcdString, 
						     RooArgSet( *rv_mu_qcd_ldp[i][j][k], *fv_sf_qcd[i][j][k], *fv_Rlsb_passfail[j][k] ) ) ;

	    rv_mu_qcd[i][j][k] = rfv_mu_qcd[i][j][k] ;


	    //---- SUSY

	    TString susyString    = "mu_susy" ;
	    TString susySlString  = "mu_susy_sl" ;
	    TString susyLdpString = "mu_susy_ldp" ;
	    
	    susyString    += sMbins[i]+sHbins[j]+sBbins[k] ;
	    susySlString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    susyLdpString += sMbins[i]+sHbins[j]+sBbins[k] ;

	    rv_mu_susy[0][0][0] = rv_mu_susy_M1_H1_1b ;

	    // for the normalization, float bin (1,1,1) and define the others using the MC ratios
	    if ( !(i == 0 && j == 0 && k == 0) ) {

	      TString rfvSusyString = "@0 * ( @1 / @2)" ;

	      rv_mu_susy[i][j][k] = new RooFormulaVar( susyString, rfvSusyString,
						       RooArgSet( *rv_mu_susymc[i][j][k], *rv_mu_susy_M1_H1_1b, *rv_mu_susymc[0][0][0] ) ) ;
	    }

	    TString rfvSusySlString  = "@0 * ( @1 / @2)" ;

	    rv_mu_susy_sl[i][j][k] = new RooFormulaVar( susySlString, rfvSusySlString,
							RooArgSet( *rv_mu_susymc_sl[i][j][k], *rv_mu_susy_M1_H1_1b, *rv_mu_susymc[0][0][0] ) ) ;

	    TString rfvSusyLdpString = "@0 * ( @1 / @2 )" ;

	    rv_mu_susy_ldp[i][j][k] = new RooFormulaVar( susyLdpString, rfvSusyLdpString,
							 RooArgSet( *rv_mu_susymc_ldp[i][j][k], *rv_mu_susy_M1_H1_1b, *rv_mu_susymc[0][0][0] ) ) ;


	    //---- Z -> nunu
	    
	    TString znnLdpString   = "mu_znn_ldp" ;
	    znnLdpString   += sMbins[i]+sHbins[j]+sBbins[k] ;

	    TString znnMcLdpString = "@0" ;

	    rv_mu_znn_ldp[i][j][k] = new RooFormulaVar( znnLdpString, znnMcLdpString, RooArgSet( *rv_mu_Znnmc_ldp[i][j][k] ) ) ;


	    //-- Float the Znn 1b vars and derive 2b and 3b vars using knn ratios


	    TString muZnnString = "mu_znn" ;
	    TString muZeeString = "mu_zee" ;
	    TString muZmmString = "mu_zmm" ;

	    muZnnString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muZeeString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muZmmString += sMbins[i]+sHbins[j]+sBbins[k] ;


	    if ( k == 0 ) {

	      TString rfvZeeString = "( @0 / @1 ) * @2 * ( ( @3 * @4 ) / ( @5 * @6 ) )" ;

	      rv_mu_zee[i][j] = new RooFormulaVar( muZeeString, rfvZeeString,
						   RooArgSet( *rv_mu_znn[i][j][k], *fv_knn[k], *fv_sf_ee[k], *fv_acc_Zee[i], 
							      *fv_eff_Zee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;


	      TString rfvZmmString = "( @0 / @1 ) * @2 * ( ( @3 * @4 ) / ( @5 * @6 ) )" ;

	      rv_mu_zmm[i][j] = new RooFormulaVar( muZmmString, rfvZmmString,
						   RooArgSet( *rv_mu_znn[i][j][k], *fv_knn[k], *fv_sf_mm[k], *fv_acc_Zmm[i], 
							      *fv_eff_Zmm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

	    }
	    else {
	      
	      TString rfvZnnString = "@0 * ( @1 / @2 )" ;

	      rv_mu_znn[i][j][k] = new RooFormulaVar( muZnnString, rfvZnnString,
						      RooArgSet( *rv_mu_znn[i][j][0], *fv_knn[k], *fv_knn[0] ) ) ;
	      
	    }

      
	    //-- Parametric relations between correlated efficiency scale factors.

	    TString effSfString      = "eff_sf" ;
	    TString effSfSlString    = "eff_sf_sl" ;
	    TString effSfLdpString   = "eff_sf_ldp" ;

	    effSfString      += sMbins[i]+sHbins[j]+sBbins[k] ;
	    effSfSlString    += sMbins[i]+sHbins[j]+sBbins[k] ;
	    effSfLdpString   += sMbins[i]+sHbins[j]+sBbins[k] ;

	    
	    TString rfvSfString = "@0 * pow( exp( @1 / @0 ), @2 )" ;

	    rv_eff_sf[i][j][k] = new RooFormulaVar( effSfString, rfvSfString,
						    RooArgSet( *rv_mean_eff_sf[i][j][k], *rv_width_eff_sf[i][j][k], 
							       *rv_mean_eff_sf[i][j][k], *eff_sf_prim ) ) ;


	    TString rfvSfSlString = "@0 * pow( exp( @1 / @0 ), @2 )" ;

	    rv_eff_sf_sl[i][j][k] = new RooFormulaVar( effSfSlString, rfvSfSlString,
						       RooArgSet( *rv_mean_eff_sf_sl[i][j][k], *rv_width_eff_sf_sl[i][j][k], 
								  *rv_mean_eff_sf_sl[i][j][k], *eff_sf_prim ) ) ;


	    TString rfvSfLdpString = "@0 * pow( exp( @1 / @0 ), @2 )" ;

	    rv_eff_sf_ldp[i][j][k] = new RooFormulaVar( effSfLdpString, rfvSfLdpString,
							RooArgSet( *rv_mean_eff_sf_ldp[i][j][k], *rv_width_eff_sf_ldp[i][j][k], 
								   *rv_mean_eff_sf_ldp[i][j][k], *eff_sf_prim ) ) ;


	    // b-tag efficiency scale factors

	    TString bteffString    = "btageff_sf" ;
	    TString bteffSlString  = "btageff_sf_sl" ;
	    TString bteffLdpString = "btageff_sf_ldp" ;
	    
	    bteffString    += sMbins[i]+sHbins[j]+sBbins[k] ;
	    bteffSlString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    bteffLdpString += sMbins[i]+sHbins[j]+sBbins[k] ;


	    // 0-lep

	    if ( rv_deff_dbtageff[i][j][k]->getVal() > 0. ) {
	      TString rfvDeffDbtString = "pow( exp( @0 * @1 ), @2 )" ;
	      rv_btageff_sf[i][j][k] = new RooFormulaVar( bteffString, rfvDeffDbtString,
							  RooArgSet( *rv_btageff_err, *rv_deff_dbtageff[i][j][k], *btageff_sf_prim ) ) ;
	    } else {
	      rv_deff_dbtageff[i][j][k]->setVal( -1.0 * rv_deff_dbtageff[i][j][k]->getVal() ) ;
	      TString rfvDeffDbtString = "pow( exp( @0 * @1 ), -1.0*@2 )" ;
	      rv_btageff_sf[i][j][k] = new RooFormulaVar( bteffString, rfvDeffDbtString,
							  RooArgSet( *rv_btageff_err, *rv_deff_dbtageff[i][j][k], *btageff_sf_prim ) ) ;
	    }


	    // 1-lep

	    if ( rv_deff_dbtageff_sl[i][j][k]->getVal() > 0. ) {
	      TString rfvDeffDbtSlString = "pow( exp( @0 * @1 ), @2 )" ;
	      rv_btageff_sf_sl[i][j][k] = new RooFormulaVar( bteffSlString, rfvDeffDbtSlString,
							  RooArgSet( *rv_btageff_err, *rv_deff_dbtageff_sl[i][j][k], *btageff_sf_prim ) ) ;
	    } else {
	      rv_deff_dbtageff_sl[i][j][k]->setVal( -1.0 * rv_deff_dbtageff_sl[i][j][k]->getVal() ) ;
	      TString rfvDeffDbtSlString = "pow( exp( @0 * @1 ), -1.0*@2 )" ;
	      rv_btageff_sf_sl[i][j][k] = new RooFormulaVar( bteffSlString, rfvDeffDbtSlString,
							  RooArgSet( *rv_btageff_err, *rv_deff_dbtageff_sl[i][j][k], *btageff_sf_prim ) ) ;
	    }


	    // ldp

	    if ( rv_deff_dbtageff_ldp[i][j][k]->getVal() > 0. ) {
	      TString rfvDeffDbtLdpString = "pow( exp( @0 * @1 ), @2 )" ;
	      rv_btageff_sf_ldp[i][j][k] = new RooFormulaVar( bteffLdpString, rfvDeffDbtLdpString,
							  RooArgSet( *rv_btageff_err, *rv_deff_dbtageff_ldp[i][j][k], *btageff_sf_prim ) ) ;
	    } else {
	      rv_deff_dbtageff_ldp[i][j][k]->setVal( -1.0 * rv_deff_dbtageff_ldp[i][j][k]->getVal() ) ;
	      TString rfvDeffDbtLdpString = "pow( exp( @0 * @1 ), -1.0*@2 )" ;
	      rv_btageff_sf_ldp[i][j][k] = new RooFormulaVar( bteffLdpString, rfvDeffDbtLdpString,
							  RooArgSet( *rv_btageff_err, *rv_deff_dbtageff_ldp[i][j][k], *btageff_sf_prim ) ) ;
	    }


	    //+++++++++++++ Expected counts for observables in terms of parameters ++++++++++++++++++

	    TString nString    = "n" ;
	    TString nSlString  = "n_sl" ;
	    TString nLdpString = "n_ldp" ;

	    nString    += sMbins[i]+sHbins[j]+sBbins[k] ;
	    nSlString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    nLdpString += sMbins[i]+sHbins[j]+sBbins[k] ;

	    TString rfvNString =  "@0 + @1 + @2 + (@3 * @4 * @5)" ;

	    rv_n[i][j][k] = new RooFormulaVar( nString, rfvNString,
	    				       RooArgSet( *rv_mu_ttwj[i][j][k], *rv_mu_qcd[i][j][k], *rv_mu_znn[i][j][k],  
							  *rv_btageff_sf[i][j][k], *rv_eff_sf[i][j][k], *rv_mu_susy[i][j][k] ) ) ;


	    TString rfvNSlString = "@0 + (@1 * @2 * @3)" ;

	    rv_n_sl[i][j][k] = new RooFormulaVar( nSlString, rfvNSlString,
						  RooArgSet( *rv_mu_ttwj_sl[i][j][k], *rv_btageff_sf_sl[i][j][k], 
							     *rv_eff_sf_sl[i][j][k], *rv_mu_susy_sl[i][j][k] ) ) ;

	    
	    TString rfvNLdpString = "@0 + @1 * @2 * ( @3 * ( @4 + @5 ) + @6 )" ;
	      
	    rv_n_ldp[i][j][k] = new RooFormulaVar( nLdpString, rfvNLdpString,
						   RooArgSet( *rv_mu_qcd_ldp[i][j][k], *rv_btageff_sf_ldp[i][j][k], *rv_eff_sf_ldp[i][j][k], 
							      *fv_sf_mc, *rv_mu_ttwj_ldp[i][j][k], *rv_mu_znn_ldp[i][j][k], *rv_mu_susy_ldp[i][j][k] ) ) ;


	    // pdf's

	    TString pdfN0lepString = "pdf_N_0lep";
	    TString pdfN1lepString = "pdf_N_1lep";
	    TString pdfNldpString  = "pdf_N_ldp";

	    pdfN0lepString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    pdfN1lepString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    pdfNldpString  += sMbins[i]+sHbins[j]+sBbins[k] ;

	    pdf_N_0lep[i][j][k] = new RooPoisson( pdfN0lepString, pdfN0lepString, *rv_0lep[i][j][k], *rv_n[i][j][k] ) ;
	    pdf_N_1lep[i][j][k] = new RooPoisson( pdfN1lepString, pdfN1lepString, *rv_1lep[i][j][k], *rv_n_sl[i][j][k] ) ;
	    pdf_N_ldp[i][j][k]  = new RooPoisson( pdfNldpString,  pdfNldpString,  *rv_ldp[i][j][k],  *rv_n_ldp[i][j][k] ) ;

	    pdflist.add( *pdf_N_0lep[i][j][k] ) ;
	    pdflist.add( *pdf_N_1lep[i][j][k] ) ;
	    pdflist.add( *pdf_N_ldp[i][j][k] ) ;

	    observedParametersList.add( *rv_0lep[i][j][k] ) ;
	    observedParametersList.add( *rv_1lep[i][j][k] ) ;
	    observedParametersList.add( *rv_ldp[i][j][k] ) ;


	  }

	    
	  TString nEeString  = "n_ee" ;
	  TString nMmString  = "n_mm" ;

          nEeString += sMbins[i]+sHbins[j] ;
          nMmString += sMbins[i]+sHbins[j] ;

	  TString rfvNeeString = "@0 / @1" ;

	  rv_n_ee[i][j] = new RooFormulaVar( nEeString, rfvNeeString,
					     RooArgSet( *rv_mu_zee[i][j], *fv_pur_Zee ) ) ;

	  
	  TString rfvNmmString = "@0 / @1" ;

	  rv_n_mm[i][j] = new RooFormulaVar( nMmString, rfvNmmString,
					     RooArgSet( *rv_mu_zmm[i][j], *fv_pur_Zmm ) ) ;


	  // pdf's

	  TString pdfZeeString = "pdf_N_Zee";
	  TString pdfZmmString = "pdf_N_Zmm";
	  
	  pdfZeeString += sMbins[i]+sHbins[j] ;
	  pdfZmmString += sMbins[i]+sHbins[j] ;

	  pdf_N_Zee[i][j] = new RooPoisson( pdfZeeString, pdfZeeString, *rv_Zee[i][j], *rv_n_ee[i][j] ) ;
	  pdf_N_Zmm[i][j] = new RooPoisson( pdfZmmString, pdfZmmString, *rv_Zmm[i][j], *rv_n_mm[i][j] ) ;

	  pdflist.add( *pdf_N_Zee[i][j] ) ;
	  pdflist.add( *pdf_N_Zmm[i][j] ) ;

	  observedParametersList.add( *rv_Zee[i][j] ) ;
	  observedParametersList.add( *rv_Zmm[i][j] ) ;

	}
      }      

      printf(" --- Constructing likelihood.\n" ) ; cout << flush ;

      pdflist.add(allNuisancePdfs);
	    
      likelihood = new RooProdPdf( "likelihood", "ra2b likelihood", pdflist ) ;

      dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
				  observedParametersList ) ;
      dsObserved->add( observedParametersList ) ;


      printf(" --- Creating workspace.\n" ) ; cout << flush ;

      RooWorkspace workspace ("ws") ;

      printf(" --- Importing dataset.\n" ) ; cout << flush ;

      workspace.import(*dsObserved);


      //--- Do a simple pre-fit to tune initial values of key parameters.
      printf(" --- Performing simple pre-fit of 0-lep ttwj normalization.\n") ; cout << flush ;

      { // begin scoping bracket.
         double initialGuess = ((RooRealVar*) rv_mu_ttwj[0][0][0]) -> getVal() ;

         // printf("\n\n Susy set to: %6.1f\n\n", ((RooRealVar*)rv_mu_susy[0][0][0]) -> getVal() ) ;

         ((RooRealVar*)rv_mu_susy[0][0][0]) -> setVal(0.) ;

           {
            double logL(0.) ;
            for (int i = 0 ; i < nBinsMET ; i++) {
               for (int j = 0 ; j < nBinsHT ; j++) {
                  for (int k = 0 ; k < nBinsBtag ; k++) {     
                     double pdfVal = pdf_N_0lep[i][j][k] -> getVal() ;
                     if ( pdfVal > 0. ) { logL += log( pdfVal ) ; } else { printf(" *** PDF %s evaluates to %g\n", pdf_N_0lep[i][j][k] -> GetName(), pdfVal ) ; }
                     printf("       PDF %20s, val=%8.6f,  N=%5.0f,  n=%6.1f\n", pdf_N_0lep[i][j][k] -> GetName(), pdfVal, rv_0lep[i][j][k]->getVal(), rv_n[i][j][k]->getVal() ) ;
                  }
               }
            }
            printf( "\n val = %6.1f, ln L = %g\n\n\n", initialGuess, logL ) ;
           }

         printf("\n\n Initial guess for 0-lep ttwj first bin: %6.1f\n\n", initialGuess ) ;

         double scanLow  = initialGuess - 6 * sqrt( initialGuess ) ;
         double scanHigh = initialGuess + 6 * sqrt( initialGuess ) ;
         double bestVal = initialGuess ;
         double bestlnL( -1.e99 ) ;
         int nScanPoints(20) ;
         if ( scanLow < 0. ) { scanLow = 1. ; }
         for ( int spi=0 ; spi<nScanPoints; spi++ ) {
            double logL(0.) ;
            double scanVal = scanLow + (scanHigh-scanLow)/(nScanPoints-1.)*spi ;
            ((RooRealVar*)rv_mu_ttwj[0][0][0]) -> setVal( scanVal ) ;
            for (int i = 0 ; i < nBinsMET ; i++) {
               for (int j = 0 ; j < nBinsHT ; j++) {
                  for (int k = 0 ; k < nBinsBtag ; k++) {     
                     double pdfVal = pdf_N_0lep[i][j][k] -> getVal() ;
                     if ( pdfVal > 0. ) { logL += log( pdfVal ) ; } else { printf(" *** PDF %s evaluates to %g\n", pdf_N_0lep[i][j][k] -> GetName(), pdfVal ) ; }
                ///  printf("       PDF %20s, val=%8.6f,  N=%5.0f,  n=%6.1f\n", pdf_N_0lep[i][j][k] -> GetName(), pdfVal, rv_0lep[i][j][k]->getVal(), rv_n[i][j][k]->getVal() ) ;
                  }
               }
            }
            if ( logL > bestlnL ) {
               bestlnL = logL ;
               bestVal = scanVal ;
            }
            printf( " val = %6.1f, ln L = %g\n", scanVal, logL ) ;
         } // spi.

         printf("\n\n  Best val = %6.1f, ln L = %g\n\n\n", bestVal, bestlnL ) ;

         ((RooRealVar*) rv_mu_ttwj[0][0][0]) -> setVal(bestVal) ;

      } // end scoping bracket.

      // parameters of interest
      RooArgSet poi(*rv_mu_susy[0][0][0], "poi");
      // flat prior for POI
      RooUniform signal_prior ("signal_prior","signal_prior",*rv_mu_susy[0][0][0]);

      printf(" --- Setting up S+B model.\n" ) ; cout << flush ;

      // signal+background model
      ModelConfig sbModel ("SbModel");
      sbModel.SetWorkspace(workspace);
      sbModel.SetPdf(*likelihood);
      sbModel.SetParametersOfInterest(poi);
      sbModel.SetPriorPdf(signal_prior);
      sbModel.SetNuisanceParameters(allNuisances);
      sbModel.SetObservables(observedParametersList);
      sbModel.SetGlobalObservables(globalObservables);


      printf(" --- Doing fit for S+B model.\n" ) ; cout << flush ;
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
      cout << "\n\n  Will save these parameter points that correspond to the fit to data" << endl << flush ;
      pPoiAndNuisance->Print("v");
      sbModel.SetSnapshot(*pPoiAndNuisance);
      workspace.import (sbModel);
      
      delete pProfile;
      delete pNll;
      delete pPoiAndNuisance;


      printf(" --- Setting up BG-only model.\n" ) ; cout << flush ;
      // background-only model
      // use the same PDF as s+b, with xsec=0
      // POI value under the background hypothesis
      ModelConfig bModel (*(RooStats::ModelConfig *)workspace.obj("SbModel"));
      bModel.SetName("BModel");
      bModel.SetWorkspace(workspace);

      printf(" --- Doing fit for BG-only model.\n" ) ; cout << flush ;
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
      cout << "\n\n  Should use these parameter points to generate pseudo data for bkg only" << endl << flush ;
      pPoiAndNuisance->Print("v");
      bModel.SetSnapshot(*pPoiAndNuisance);
      workspace.import (bModel);
      
      delete pProfile;
      delete pNll;
      delete pPoiAndNuisance;
      

      workspace.Print() ;
      printf("\n\n Creating output root file: ws.root\n\n\n") ;
      workspace.writeToFile("ws.root");

      return true;

    }  // end of initialize


   //=====================================================================================================================

    bool ra2bRoostatsClass3D_1::setSusyScanPoint( const char* inputScanFile,
						  double m0, double m12, bool isT1bbbb, double t1bbbbXsec,
						  const char* inputSusy_deff_dbtageff_file
						  ) {

      printf("\n\n Opening SUSY scan input file : %s\n", inputScanFile ) ;

      ifstream infp ;
      infp.open(inputScanFile) ;
      if ( !infp.good() ) {
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

     
      TString sMbins[nBinsMET];
      TString sHbins[nBinsHT];
      TString sBbins[3] = {"_1b","_2b","_3b"};
      
      for (int i = 0 ; i < nBinsMET ; i++) {
	TString base = "M";
	stringstream sbin;
	sbin << i+1;
	base += sbin.str();
	sMbins[i] = base;
      }
      
      for (int j = 0 ; j < nBinsHT ; j++) {
	TString base = "_H";
	stringstream sbin;
	sbin << j+1;
	base += sbin.str();
	sHbins[j] = base;
      }


      //--- Loop over the scan points.
      while ( infp.good() ) {

	float pointM0 ;
	float pointM12 ;
	
	float n_0l_raw[nBinsMET][nBinsHT][nBinsBtag] ;
	float n_1l_raw[nBinsMET][nBinsHT][nBinsBtag] ;
	float n_ldp_raw[nBinsMET][nBinsHT][nBinsBtag] ;
	
	float n_0l_correction[nBinsMET][nBinsHT][nBinsBtag] ;
	float n_1l_correction[nBinsMET][nBinsHT][nBinsBtag] ;
	float n_ldp_correction[nBinsMET][nBinsHT][nBinsBtag] ;
	
	float n_0l_error[nBinsMET][nBinsHT][nBinsBtag] ;
	float n_1l_error[nBinsMET][nBinsHT][nBinsBtag] ;
	float n_ldp_error[nBinsMET][nBinsHT][nBinsBtag] ;

	int nGen ;

	int ArraySize = 3 + 9*(nBinsMET*nBinsHT*nBinsBtag) ;
	double ArrayContent[ArraySize] ;

	for (int i = 0; infp && i < ArraySize; ++ i) {
	  infp >> ArrayContent[i];
	}

	pointM0  = ArrayContent[0] ;
	pointM12 = ArrayContent[1] ;

	nGen = (int)ArrayContent[2] ;

	int nBins = nBinsMET*nBinsHT*nBinsBtag*3 ;

	for (int i = 0 ; i < nBinsMET ; i++) {
	  for (int j = 0 ; j < nBinsHT ; j++) {
	    for (int k = 0 ; k < nBinsBtag ; k++) {     

	      n_0l_raw[i][j][k]  = ArrayContent[3 + i*(nBinsHT*nBinsBtag*3) + j*(nBinsBtag*3) + k*3] ;
	      n_1l_raw[i][j][k]  = ArrayContent[4 + i*(nBinsHT*nBinsBtag*3) + j*(nBinsBtag*3) + k*3] ;
	      n_ldp_raw[i][j][k] = ArrayContent[5 + i*(nBinsHT*nBinsBtag*3) + j*(nBinsBtag*3) + k*3] ;

	      n_0l_correction[i][j][k]  = ArrayContent[3 + nBins + i*(nBinsHT*nBinsBtag*3) + j*(nBinsBtag*3) + k*3] ;
	      n_1l_correction[i][j][k]  = ArrayContent[4 + nBins + i*(nBinsHT*nBinsBtag*3) + j*(nBinsBtag*3) + k*3] ;
	      n_ldp_correction[i][j][k] = ArrayContent[5 + nBins + i*(nBinsHT*nBinsBtag*3) + j*(nBinsBtag*3) + k*3] ;

	      n_0l_error[i][j][k]  = ArrayContent[3 + 2*nBins + i*(nBinsHT*nBinsBtag*3) + j*(nBinsBtag*3) + k*3] ;
	      n_1l_error[i][j][k]  = ArrayContent[4 + 2*nBins + i*(nBinsHT*nBinsBtag*3) + j*(nBinsBtag*3) + k*3] ;
	      n_ldp_error[i][j][k] = ArrayContent[5 + 2*nBins + i*(nBinsHT*nBinsBtag*3) + j*(nBinsBtag*3) + k*3] ;

	    }
	  }
	}


	printf("  pointM0 = %g , pointM12 = %g\n", pointM0, pointM12 ) ;

	if (    fabs( pointM0 - m0 ) <= deltaM0/2 && fabs( pointM12 - m12 ) <= deltaM12/2 ) {
	  
	  double nGenPerPoint = 10000 ; // for SMS's
	  
	  cout << "Systematic uncertainties (in %):" << endl;

	  for (int i = 0 ; i < nBinsMET ; i++) {
	    for (int j = 0 ; j < nBinsHT ; j++) {
	      for (int k = 0 ; k < nBinsBtag ; k++) {     

		TString binString = "";
		binString += sMbins[i]+sHbins[j]+sBbins[k] ;
		
		cout << binString + " - 0 lep  = " << n_0l_error[i][j][k] << endl ; 
		cout << binString + " - 1 lep  = " << n_1l_error[i][j][k] << endl ; 
		cout << binString + " - ldp    = " << n_ldp_error[i][j][k] << endl ; 
		
	      }
	    }
	  }
	  

	  //--Include the stat error on the efficiency for SMS's

	  if ( isT1bbbb ) {

	    float n_0l_raw_eff[nBinsMET][nBinsHT][nBinsBtag] ;
	    float n_1l_raw_eff[nBinsMET][nBinsHT][nBinsBtag] ;
	    float n_ldp_raw_eff[nBinsMET][nBinsHT][nBinsBtag] ;

	    float n_0l_stat_error[nBinsMET][nBinsHT][nBinsBtag] ;
	    float n_1l_stat_error[nBinsMET][nBinsHT][nBinsBtag] ;
	    float n_ldp_stat_error[nBinsMET][nBinsHT][nBinsBtag] ;


	    for (int i = 0 ; i < nBinsMET ; i++) {
	      for (int j = 0 ; j < nBinsHT ; j++) {
		for (int k = 0 ; k < nBinsBtag ; k++) {     

		  TString binString = "";
		  binString += sMbins[i]+sHbins[j]+sBbins[k] ;

		  // absolute raw efficiency
		  n_0l_raw_eff[i][j][k]  = n_0l_raw[i][j][k] / nGenPerPoint ;
		  n_1l_raw_eff[i][j][k]  = n_1l_raw[i][j][k] / nGenPerPoint ;
		  n_ldp_raw_eff[i][j][k] = n_ldp_raw[i][j][k] / nGenPerPoint ;

		  // absolute stat error
		  n_0l_stat_error[i][j][k]  = sqrt( n_0l_raw_eff[i][j][k]  * ( 1.0 - n_0l_raw_eff[i][j][k]  ) / nGenPerPoint ) ;
		  n_1l_stat_error[i][j][k]  = sqrt( n_1l_raw_eff[i][j][k]  * ( 1.0 - n_1l_raw_eff[i][j][k]  ) / nGenPerPoint ) ;
		  n_ldp_stat_error[i][j][k] = sqrt( n_ldp_raw_eff[i][j][k] * ( 1.0 - n_ldp_raw_eff[i][j][k] ) / nGenPerPoint ) ;

		  // relative stat err in percent.
		  if ( n_0l_raw_eff[i][j][k] > 0 ) { n_0l_stat_error[i][j][k] = 100.* n_0l_stat_error[i][j][k] / n_0l_raw_eff[i][j][k] ; } 
		  else { n_0l_stat_error[i][j][k] = 0. ; }

		  if ( n_1l_raw_eff[i][j][k] > 0 ) { n_1l_stat_error[i][j][k] = 100.* n_1l_stat_error[i][j][k] / n_1l_raw_eff[i][j][k] ; } 
		  else { n_1l_stat_error[i][j][k] = 0. ; }

		  if ( n_ldp_raw_eff[i][j][k] > 0 ) { n_ldp_stat_error[i][j][k] = 100.* n_ldp_stat_error[i][j][k] / n_ldp_raw_eff[i][j][k] ; } 
		  else { n_ldp_stat_error[i][j][k] = 0. ; }
		  
		  cout << binString + " - 0 lep - SUSY statistical uncertainty (%) = " << n_0l_stat_error[i][j][k] << endl ; 
		  cout << binString + " - 1 lep - SUSY statistical uncertainty (%) = " << n_1l_stat_error[i][j][k] << endl ; 
		  cout << binString + " - ldp   - SUSY statistical uncertainty (%) = " << n_ldp_stat_error[i][j][k] << endl ; 

		  // total error
		  n_0l_error[i][j][k]  = sqrt( pow( n_0l_error[i][j][k],  2) + pow( n_0l_stat_error[i][j][k], 2) ) ;
		  n_1l_error[i][j][k]  = sqrt( pow( n_1l_error[i][j][k],  2) + pow( n_1l_stat_error[i][j][k], 2) ) ;
		  n_ldp_error[i][j][k] = sqrt( pow( n_ldp_error[i][j][k], 2) + pow( n_ldp_stat_error[i][j][k], 2) ) ;
		  
		  cout << binString + " - 0 lep - SUSY total uncertainty (%) = " << n_0l_error[i][j][k] << endl ; 
		  cout << binString + " - 1 lep - SUSY total uncertainty (%) = " << n_1l_error[i][j][k] << endl ; 
		  cout << binString + " - ldp   - SUSY total uncertainty (%) = " << n_ldp_error[i][j][k] << endl ; 

		}
	      }
	    }
	  } // end of if (isT1bbbb)


	  double setVal_n_0l[nBinsMET][nBinsHT][nBinsBtag] ;
	  double setVal_n_1l[nBinsMET][nBinsHT][nBinsBtag] ;
	  double setVal_n_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

	  for (int i = 0 ; i < nBinsMET ; i++) {
	    for (int j = 0 ; j < nBinsHT ; j++) {
	      for (int k = 0 ; k < nBinsBtag ; k++) {     
		
		if (!isT1bbbb) {
		  setVal_n_0l[i][j][k]  = n_0l_raw[i][j][k]  * n_0l_correction[i][j][k] ;
		  setVal_n_1l[i][j][k]  = n_1l_raw[i][j][k]  * n_1l_correction[i][j][k] ;
		  setVal_n_ldp[i][j][k] = n_ldp_raw[i][j][k] * n_ldp_correction[i][j][k] ;
		}
		else {
		  setVal_n_0l[i][j][k]  = DataLumi * t1bbbbXsec * (( n_0l_raw[i][j][k]  * n_0l_correction[i][j][k] ) / nGenPerPoint ) ;
		  setVal_n_1l[i][j][k]  = DataLumi * t1bbbbXsec * (( n_1l_raw[i][j][k]  * n_1l_correction[i][j][k] ) / nGenPerPoint ) ;
		  setVal_n_ldp[i][j][k] = DataLumi * t1bbbbXsec * (( n_ldp_raw[i][j][k] * n_ldp_correction[i][j][k] ) / nGenPerPoint ) ;
		}


		rv_mu_susymc[i][j][k]     -> setVal( setVal_n_0l[i][j][k] ) ;
		rv_mu_susymc_sl[i][j][k]  -> setVal( setVal_n_1l[i][j][k] ) ;
		rv_mu_susymc_ldp[i][j][k] -> setVal( setVal_n_ldp[i][j][k] ) ;

		rv_width_eff_sf[i][j][k]     -> setVal( n_0l_error[i][j][k] / 100. ) ;
		rv_width_eff_sf_sl[i][j][k]  -> setVal( n_1l_error[i][j][k] / 100. ) ;
		rv_width_eff_sf_ldp[i][j][k] -> setVal( n_ldp_error[i][j][k] / 100. ) ;

	      }
	    }
	  }


	  if ( !isT1bbbb ) {
	    printf("\n\n Found point m0 = %4.0f,  m1/2 = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, setVal_n_0l[0][0][0] ) ;
	  } else {
	    printf("\n\n Found point mGluino = %4.0f,  mLSP = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, setVal_n_0l[0][0][0]) ;
	  }
	  

	  // print out signal values

	  cout << "\n\nSetting susy signal: \n" << endl ;

          double total0lep(0.) ;
	  
	  for (int i = 0 ; i < nBinsMET ; i++) {
	    for (int j = 0 ; j < nBinsHT ; j++) {
	      for (int k = 0 ; k < nBinsBtag ; k++) {     

		TString binString = "";
		binString += sMbins[i]+sHbins[j]+sBbins[k] ;	  

		cout << binString + " - 0 lep - setting susy signal to " << setVal_n_0l[i][j][k] << endl ;
		cout << binString + " - 1 lep - setting susy signal to " << setVal_n_1l[i][j][k] << endl ;
		cout << binString + " - ldp   - setting susy signal to " << setVal_n_ldp[i][j][k] << endl ;

                total0lep += setVal_n_0l[i][j][k] ;

	      }
	    }
	  }

          printf("\n\n SUSY 0lep total: %7.2f\n\n", total0lep ) ;

	  found = true ;

	  break ;

	}

      }


      infp.close() ;

      if ( !found ) {
	printf("\n\n *** Point not found in scan.  Check this file %s\n  to see if this point is there: mgl=%4.0f, mlsp=%4.0f\n\n", inputScanFile, m0, m12 ) ;
	return false ;
      }


      //----- Now, read in the deff_dbtageff numbers.

      printf("\n\n Opening SUSY deff_dbtageff input file : %s\n", inputSusy_deff_dbtageff_file ) ;

      ifstream infb ;
      infb.open(inputSusy_deff_dbtageff_file) ;
      if ( !infb.good() ) {
	printf("\n\n *** Problem opening input file: %s.\n\n", inputSusy_deff_dbtageff_file ) ;
	return false ;
      } 

      found = false ;

      while ( infb.good() ) {

	float pointM0 ;
	float pointM12 ;

	float deff_dbtageff_0l[nBinsMET][nBinsHT][nBinsBtag] ;
	float deff_dbtageff_1l[nBinsMET][nBinsHT][nBinsBtag] ;
	float deff_dbtageff_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
	
	int ArraySize = 2 + 3*(nBinsMET*nBinsHT*nBinsBtag) ;
	double ArrayContent[ArraySize] ;

	for (int i = 0; infb && i < ArraySize; ++ i) {
	  infb >> ArrayContent[i];
	}

	pointM0  = ArrayContent[0] ;
	pointM12 = ArrayContent[1] ;

	int nBins = nBinsMET*nBinsHT*nBinsBtag ;

	for (int i = 0 ; i < nBinsMET ; i++) {
	  for (int j = 0 ; j < nBinsHT ; j++) {
	    for (int k = 0 ; k < nBinsBtag ; k++) {     

	      deff_dbtageff_0l[i][j][k]  = ArrayContent[2 + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;
	      deff_dbtageff_1l[i][j][k]  = ArrayContent[2 + nBins + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;
	      deff_dbtageff_ldp[i][j][k] = ArrayContent[2 + 2*nBins + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;

	    }
	  }
	}

	printf("  pointM0 = %g , pointM12 = %g\n", pointM0, pointM12 ) ;

	if (    fabs( pointM0 - m0 ) <= deltaM0/2. && fabs( pointM12 - m12 ) <= deltaM12/2. ) {
	 
	  cout << "Setting susy deff_dbtag derivatives: " << endl ;

	  for (int i = 0 ; i < nBinsMET ; i++) {
	    for (int j = 0 ; j < nBinsHT ; j++) {
	      for (int k = 0 ; k < nBinsBtag ; k++) {     
		
		TString binString = "";
		binString += sMbins[i]+sHbins[j]+sBbins[k] ;	  	  

		cout << binString + " - 0 lep - setting susy deff_dbtag to " << deff_dbtageff_0l[i][j][k] << endl ;
		cout << binString + " - 1 lep - setting susy deff_dbtag to " << deff_dbtageff_1l[i][j][k] << endl ;
		cout << binString + " - ldp   - setting susy deff_dbtag to " << deff_dbtageff_ldp[i][j][k] << endl ;

		rv_deff_dbtageff[i][j][k]     -> setVal ( deff_dbtageff_0l[i][j][k] ) ;
		rv_deff_dbtageff_sl[i][j][k]  -> setVal ( deff_dbtageff_1l[i][j][k] ) ;
		rv_deff_dbtageff_ldp[i][j][k] -> setVal ( deff_dbtageff_ldp[i][j][k] ) ;

	      }
	    }
	  }

	  found = true ;

	  break ;
 
	}

      }


      if ( found ) {
	return true ;
      } else {
	printf("\n\n *** Point not found in scan.  Check this file %s\n  to see if this point is there: mgl=%4.0f, mlsp=%4.0f\n\n", inputSusy_deff_dbtageff_file, m0, m12 ) ;
	return false ;
      }

    }  // end of setSusyScanPoint

   //==============================================================================================================

    void ra2bRoostatsClass3D_1::mismatchErr( char* label, TString inPar ) {

      cout << "Mismatch in input file:" << endl;
      cout << "Expecting: " << inPar << endl;
      cout << "Reading:   " << label << endl;
      cout << "\nCheck binning!\n" << endl;

      return;

    }

   //==============================================================================================================








