#include "ra2bRoostatsClass3D_2.h"

//
//  Owen: The difference between ra2bRoostatsClass3D_1 and ra2bRoostatsClass3D_2
//        is in the implementation of the nuisance parameters.
//
//          ra2bRoostatsClass3D_1 : uses log normals
//
//          ra2bRoostatsClass3D_2 : uses beta and beta prime distributions.
//
//


#include <iostream>
#include <sstream>
#include <string.h>

#include "TRoot.h"
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

#include "RooRatio.h"
#include "betaPrimeConstraint.c"

  using namespace RooFit ;
  using namespace RooStats ;


  //=====================================================================================================


   ra2bRoostatsClass3D_2::ra2bRoostatsClass3D_2() {

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

   ra2bRoostatsClass3D_2::~ra2bRoostatsClass3D_2() { }



  //===================================================================

    bool ra2bRoostatsClass3D_2::initialize( const char* infile ,
					    const char* inputScanFile,
					    double m0, double m12, bool isT1bbbb, double t1bbbbXsec,
					    const char* inputSusy_deff_dbtageff_file,
                                            int   qcdModelIndex
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










      printf(" --- Creating workspace.\n" ) ; cout << flush ;

      RooWorkspace workspace ("ws") ;
      workspace.autoImportClassCode(true);


      globalObservables      = new RooArgSet("globalObservables");
      allNuisances           = new RooArgSet("allNuisances");
      allNuisancePdfs        = new RooArgSet("allNuisancePdfs");
      observedParametersList = new RooArgSet("observables") ;









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

	    observedParametersList -> add( *rv_0lep[i][j][k] ) ;
	    observedParametersList -> add( *rv_1lep[i][j][k] ) ;
	    observedParametersList -> add( *rv_ldp[i][j][k] ) ;
	    
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

	  observedParametersList -> add( *rv_Zee[i][j] ) ;
	  observedParametersList -> add( *rv_Zmm[i][j] ) ;

	}
      }

      RooDataSet* dsObserved ;
      dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
				  *observedParametersList ) ;
      dsObserved->add( *observedParametersList ) ;

      printf(" --- Importing dataset.\n" ) ; cout << flush ;

      workspace.import(*dsObserved);



      char NP_name[1000] ;
      char NP_base_name[1000] ;




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



      //--- QCD model 1:  Use LSB pass/fail ratio for QCD 0lep/LDP ratio.
      RooRealVar* Rlsb_passfail_prim[nBinsHT][nBinsBtag];
      RooRealVar* Rlsb_passfail_nom[nBinsHT][nBinsBtag];
      RooGaussian* pdf_Rlsb_passfail[nBinsHT][nBinsBtag];
      RooFormulaVar* fv_Rlsb_passfail[nBinsHT][nBinsBtag];

      //--- QCD model 2:  Fit for QCD 0lep/LDP ratio in each HT bin.
      RooRealVar* rv_qcd_0lepLDP_ratio[20] ;

      if ( qcdModelIndex == 1 ) {

         printf("\n\n  QCD Model 1 : Use LSB pass/fail ratio for QCD 0lep/LDP ratio.\n\n") ;

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
             globalObservables -> add (*Rlsb_passfail_nom[j][k]);
             allNuisances -> add (*Rlsb_passfail_prim[j][k]);
             allNuisancePdfs -> add (*pdf_Rlsb_passfail[j][k]);

           }
         }

      } else if ( qcdModelIndex == 2 ) {

         printf("\n\n  QCD Model 2 : One floating 0lep/LDP ratio for each HT bin..\n\n") ;

         for ( int htbi=0; htbi < nBinsHT ; htbi++ ) {
            char vname[1000] ;
            sprintf( vname, "qcd_0lepLDP_ratio_H%d", htbi+1 ) ;
            rv_qcd_0lepLDP_ratio[htbi] = new RooRealVar( vname, vname, 0.05, 0., 10. ) ;
         }

      } else {

         printf("\n\n *** Unknown QCD model index : %d\n\n", qcdModelIndex ) ;
         return false ;

      }






      // Z -> ll acceptance

      RooAbsReal* rar_acc_Zee[nBinsMET] ;
      RooAbsReal* rar_acc_Zmm[nBinsMET] ;

      for (int i = 0 ; i < nBinsMET ; i++) {

         sprintf( NP_name, "acc_Zee_M%d", i+1 ) ;
         rar_acc_Zee[i] = makeBetaConstraint( NP_name, acc_Zee[i], acc_Zee_err[i] ) ;

         sprintf( NP_name, "acc_Zmm_M%d", i+1 ) ;
         rar_acc_Zmm[i] = makeBetaConstraint( NP_name, acc_Zmm[i], acc_Zmm_err[i] ) ;

      }












      // Z -> ll efficiencies

       sprintf( NP_name, "eff_Zee" ) ;
       RooAbsReal* rar_eff_Zee = makeBetaConstraint( NP_name, eff_Zee, eff_Zee_err ) ;

       sprintf( NP_name, "eff_Zmm" ) ;
       RooAbsReal* rar_eff_Zmm = makeBetaConstraint( NP_name, eff_Zmm, eff_Zmm_err ) ;











      // Z -> ll purities

       sprintf( NP_name, "pur_Zee" ) ;
       RooAbsReal* rar_pur_Zee = makeBetaConstraint( NP_name, pur_Zee, pur_Zee_err ) ;

       sprintf( NP_name, "pur_Zmm" ) ;
       RooAbsReal* rar_pur_Zmm = makeBetaConstraint( NP_name, pur_Zmm, pur_Zmm_err ) ;










      // VL -> nominal scale factors

      RooAbsReal* rar_knn[nBinsBtag] ;

      for (int k = 0 ; k < nBinsBtag ; k++) {

       sprintf( NP_name, "knn_%db", k+1 ) ;
       rar_knn[k] = makeBetaPrimeConstraint( NP_name, knn[k], knn_err[k] ) ;

      }









      // sf_ee and sf_mm derived from a common underlying parameter

      sprintf( NP_name, "sf_ll" ) ;
      makeBetaPrimeConstraint( NP_name, 1.0, 0.1 ) ;

      RooAbsReal* rar_sf_ee[nBinsBtag] ;
      RooAbsReal* rar_sf_mm[nBinsBtag] ;

      sprintf( NP_base_name, "sf_ll" ) ;

      for (int k = 0 ; k < nBinsBtag ; k++) {

         sprintf( NP_name, "sf_ee_%db", k+1 ) ;
         rar_sf_ee[k] = makeCorrelatedBetaPrimeConstraint( NP_name, sf_ee[k], sf_ee_err[k], NP_base_name ) ;

         sprintf( NP_name, "sf_mm_%db", k+1 ) ;
         rar_sf_mm[k] = makeCorrelatedBetaPrimeConstraint( NP_name, sf_mm[k], sf_mm_err[k], NP_base_name ) ;

      }







      // MC scale factor

      sprintf( NP_name, "sf_mc" ) ;
      RooAbsReal* rar_sf_mc = makeBetaPrimeConstraint( NP_name, sf_mc, sf_mc_err ) ;








      // QCD and TTWJ scale factors

      RooAbsReal* rar_sf_qcd [nBinsMET][nBinsHT][nBinsBtag];
      RooAbsReal* rar_sf_ttwj[nBinsMET][nBinsHT][nBinsBtag];

      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBtag ; k++) {


             //--- QCD

             sprintf( NP_name, "sf_qcd_M%d_H%d_%db", i+1, j+1, k+1 ) ;
             rar_sf_qcd[i][j][k] = makeBetaPrimeConstraint( NP_name, sf_qcd[i][j][k], sf_qcd_err[i][j][k] ) ;


             //--- ttwj

             if ( ! (i==0 && j==0 && k==0 ) ) {

                sprintf( NP_name, "sf_ttwj_M%d_H%d_%db", i+1, j+1, k+1 ) ;
                rar_sf_ttwj[i][j][k] = makeBetaPrimeConstraint( NP_name, sf_ttwj[i][j][k], sf_ttwj_err[i][j][k] ) ;

             }

             printf("--------\n") ;

          } // i (met)
             printf("++++++++++++++++\n") ;
        } // j (ht)
             printf("=========================\n") ;
      } // k (nbtag)











      //zzzzz needs updating to beta zzzzzzzzzz

      RooRealVar* eff_sf_prim = new RooRealVar( "eff_sf_prim", "eff_sf_prim", 0, -5, 5);
      RooRealVar* eff_sf_nom = new RooRealVar( "eff_sf_nom", "eff_sf_nom", 0, -5, 5);
      RooGaussian* pdf_eff_sf = new RooGaussian( "pdf_eff_sf" , "pdf_eff_sf", *eff_sf_prim, *eff_sf_nom, RooConst(1));
      eff_sf_nom->setConstant();
      globalObservables -> add (*eff_sf_nom);
      allNuisances -> add (*eff_sf_prim);
      allNuisancePdfs -> add (*pdf_eff_sf);


      RooRealVar* btageff_sf_prim = new RooRealVar( "btageff_sf_prim", "btageff_sf_prim", 0, -5, 5);
      RooRealVar* btageff_sf_nom = new RooRealVar( "btageff_sf_nom", "btageff_sf_nom", 0, -5, 5);
      RooGaussian* pdf_btageff_sf = new RooGaussian( "pdf_btageff_sf", "pdf_btageff_sf", *btageff_sf_prim, *btageff_sf_nom, RooConst(1));
      btageff_sf_nom->setConstant();
      globalObservables -> add (*btageff_sf_nom);
      allNuisances -> add (*btageff_sf_prim);
      allNuisancePdfs -> add (*pdf_btageff_sf);

















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
	                      RooArgSet( *rv_mu_ttwj[0][0][0], *rar_sf_ttwj[i][j][k], *rv_mu_ttwj_sl[i][j][k], *rv_mu_ttwj_sl[0][0][0] )) ;
	      
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

            if ( qcdModelIndex == 1 ) {

               rfv_mu_qcd[i][j][k] = new RooFormulaVar( qcdString, rfvQcdString, 
                                                        RooArgSet( *rv_mu_qcd_ldp[i][j][k], *rar_sf_qcd[i][j][k], *fv_Rlsb_passfail[j][k] ) ) ;


            } else if ( qcdModelIndex == 2 ) {

               rfv_mu_qcd[i][j][k] = new RooFormulaVar( qcdString, rfvQcdString, 
                                                        RooArgSet( *rv_mu_qcd_ldp[i][j][k], *rar_sf_qcd[i][j][k], *rv_qcd_0lepLDP_ratio[j] ) ) ;


            }

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
						   RooArgSet( *rv_mu_znn[i][j][k], *rar_knn[k], *rar_sf_ee[k], *rar_acc_Zee[i], 
							      *rar_eff_Zee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;


	      TString rfvZmmString = "( @0 / @1 ) * @2 * ( ( @3 * @4 ) / ( @5 * @6 ) )" ;

	      rv_mu_zmm[i][j] = new RooFormulaVar( muZmmString, rfvZmmString,
						   RooArgSet( *rv_mu_znn[i][j][k], *rar_knn[k], *rar_sf_mm[k], *rar_acc_Zmm[i], 
							      *rar_eff_Zmm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

	    }
	    else {
	      
	      TString rfvZnnString = "@0 * ( @1 / @2 )" ;

	      rv_mu_znn[i][j][k] = new RooFormulaVar( muZnnString, rfvZnnString,
						      RooArgSet( *rv_mu_znn[i][j][0], *rar_knn[k], *rar_knn[0] ) ) ;
	      
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
							      *rar_sf_mc, *rv_mu_ttwj_ldp[i][j][k], *rv_mu_znn_ldp[i][j][k], *rv_mu_susy_ldp[i][j][k] ) ) ;


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



	  }

	    
	  TString nEeString  = "n_ee" ;
	  TString nMmString  = "n_mm" ;

          nEeString += sMbins[i]+sHbins[j] ;
          nMmString += sMbins[i]+sHbins[j] ;

	  TString rfvNeeString = "@0 / @1" ;

	  rv_n_ee[i][j] = new RooFormulaVar( nEeString, rfvNeeString,
					     RooArgSet( *rv_mu_zee[i][j], *rar_pur_Zee ) ) ;

	  
	  TString rfvNmmString = "@0 / @1" ;

	  rv_n_mm[i][j] = new RooFormulaVar( nMmString, rfvNmmString,
					     RooArgSet( *rv_mu_zmm[i][j], *rar_pur_Zmm ) ) ;


	  // pdf's

	  TString pdfZeeString = "pdf_N_Zee";
	  TString pdfZmmString = "pdf_N_Zmm";
	  
	  pdfZeeString += sMbins[i]+sHbins[j] ;
	  pdfZmmString += sMbins[i]+sHbins[j] ;

	  pdf_N_Zee[i][j] = new RooPoisson( pdfZeeString, pdfZeeString, *rv_Zee[i][j], *rv_n_ee[i][j] ) ;
	  pdf_N_Zmm[i][j] = new RooPoisson( pdfZmmString, pdfZmmString, *rv_Zmm[i][j], *rv_n_mm[i][j] ) ;

	  pdflist.add( *pdf_N_Zee[i][j] ) ;
	  pdflist.add( *pdf_N_Zmm[i][j] ) ;


	}
      }      

      printf(" --- Constructing likelihood.\n" ) ; cout << flush ;

      pdflist.add( *allNuisancePdfs );
	    
      pdflist.Print() ;

      likelihood = new RooProdPdf( "likelihood", "ra2b likelihood", pdflist ) ;



      /// likelihood->printMultiline( cout, 1, kTRUE, "" ) ;


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
      sbModel.SetNuisanceParameters( *allNuisances );
      sbModel.SetObservables( *observedParametersList );
      sbModel.SetGlobalObservables( *globalObservables );

      workspace.Print() ;

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

    bool ra2bRoostatsClass3D_2::setSusyScanPoint( const char* inputScanFile,
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

    void ra2bRoostatsClass3D_2::mismatchErr( char* label, TString inPar ) {

      cout << "Mismatch in input file:" << endl;
      cout << "Expecting: " << inPar << endl;
      cout << "Reading:   " << label << endl;
      cout << "\nCheck binning!\n" << endl;

      return;

    }

   //==============================================================================================================


    RooAbsReal* ra2bRoostatsClass3D_2::makeBetaPrimeConstraint( const char* NP_name, double NP_val, double NP_err ) {

       double alpha, beta ;
       char varname[1000] ;

       RooRealVar *rrv_passObs, *rrv_failObs, *rrv_passPar, *rrv_failPar ;
       RooPoisson *passConstraint, *failConstraint ;

       double parVal, parErr, upperLimit, lowerLimit ;

       betaPrimeModeTransform( NP_val, NP_err, alpha, beta ) ;

       sprintf( varname, "passObs_%s", NP_name ) ;
       rrv_passObs = new RooRealVar( varname, varname, alpha-1., 1e-5, 1e5 ) ;
       rrv_passObs -> setConstant( kTRUE ) ;

       sprintf( varname, "failObs_%s", NP_name ) ;
       rrv_failObs = new RooRealVar( varname, varname, beta -1., 1e-5, 1e5 ) ;
       rrv_failObs -> setConstant( kTRUE ) ;

       globalObservables -> add( *rrv_passObs ) ;
       globalObservables -> add( *rrv_failObs ) ;


       parVal = alpha-1. ;
       lowerLimit = parVal - 6*sqrt(parVal) ;
       if ( lowerLimit <= 0. ) { lowerLimit = 1e-5 ; }
       upperLimit = parVal + 6*sqrt(parVal) ;
       if ( parVal > 0. ) { parErr = sqrt( parVal ) ; } else { parErr = 0. ; }
       sprintf( varname, "passPar_%s", NP_name ) ;
       rrv_passPar = new RooRealVar( varname, varname, parVal, lowerLimit, upperLimit ) ;
       rrv_passPar -> setError( parErr ) ;
       rrv_passPar->setConstant( kFALSE ) ;
       printf(" floating nuisance parameter: %s = %g +/- %g, [%g, %g]\n",
           varname, parVal, parErr, lowerLimit, upperLimit ) ;

       parVal = beta-1. ;
       lowerLimit = parVal - 6*sqrt(parVal) ;
       if ( lowerLimit <= 0. ) { lowerLimit = 1e-5 ; }
       upperLimit = parVal + 6*sqrt(parVal) ;
       if ( parVal > 0. ) { parErr = sqrt( parVal ) ; } else { parErr = 0. ; }
       sprintf( varname, "failPar_%s", NP_name ) ;
       rrv_failPar = new RooRealVar( varname, varname, parVal, lowerLimit, upperLimit ) ;
       rrv_failPar -> setError( parErr ) ;
       rrv_failPar->setConstant( kFALSE ) ;
       printf(" floating nuisance parameter: %s = %g +/- %g, [%g, %g]\n",
           varname, parVal, parErr, lowerLimit, upperLimit ) ;

       allNuisances -> add( *rrv_passPar ) ;
       allNuisances -> add( *rrv_failPar ) ;


       sprintf( varname, "passConstraint_%s", NP_name ) ;
       passConstraint = new RooPoisson( varname, varname, *rrv_passObs, *rrv_passPar ) ;
       printf("  Created constraint PDF : %s, val = %g, logval = %g\n",
         varname, passConstraint->getVal(), passConstraint->getLogVal() ) ;

       sprintf( varname, "failConstraint_%s", NP_name ) ;
       failConstraint = new RooPoisson( varname, varname, *rrv_failObs, *rrv_failPar ) ;
       printf("  Created constraint PDF : %s, val = %g, logval = %g\n",
         varname, failConstraint->getVal(), failConstraint->getLogVal() ) ;

       allNuisancePdfs -> add( *passConstraint ) ;
       allNuisancePdfs -> add( *failConstraint ) ;

       sprintf( varname, "%s", NP_name ) ;
       RooAbsReal* rar = new RooRatio( varname, varname, *rrv_passPar, *rrv_failPar ) ;
       printf("  Created nuisance parameter %s : val = %g\n", varname, rar -> getVal() ) ;

       return rar ;

    } // makeBetaPrimeConstraint.


   //==============================================================================================================



    RooAbsReal* ra2bRoostatsClass3D_2::makeBetaConstraint( const char* NP_name, double NP_val, double NP_err ) {

       if ( NP_val >=1 ) {
          printf("\n\n *** ra2bRoostatsClass3D_2::makeBetaConstraint:  warning.  Input efficiency value is %g.  Resetting to 0.999.\n\n", NP_val) ;
          NP_val = 0.999 ;
       }

       double alpha, beta ;
       char varname[1000] ;

       RooRealVar *rrv_passObs, *rrv_failObs, *rrv_passPar, *rrv_failPar ;
       RooPoisson *passConstraint, *failConstraint ;

       double parVal, parErr, upperLimit, lowerLimit ;

       betaModeTransform( NP_val, NP_err, alpha, beta ) ;

       sprintf( varname, "passObs_%s", NP_name ) ;
       rrv_passObs = new RooRealVar( varname, varname, alpha-1., 1e-5, 1e5 ) ;
       rrv_passObs -> setConstant( kTRUE ) ;

       sprintf( varname, "failObs_%s", NP_name ) ;
       ///// rrv_failObs = new RooRealVar( varname, varname, alpha+beta -2., 1e-5, 1e5 ) ;
       ///// rrv_failObs = new RooRealVar( varname, varname, beta +1., 1e-5, 1e5 ) ;
       rrv_failObs = new RooRealVar( varname, varname, beta -1., 1e-5, 1e5 ) ;
       rrv_failObs -> setConstant( kTRUE ) ;

       globalObservables -> add( *rrv_passObs ) ;
       globalObservables -> add( *rrv_failObs ) ;


       parVal = alpha-1. ;
       lowerLimit = parVal - 6*sqrt(parVal) ;
       if ( lowerLimit <= 0. ) { lowerLimit = 1e-5 ; }
       upperLimit = parVal + 6*sqrt(parVal) ;
       if ( parVal > 0. ) { parErr = sqrt( parVal ) ; } else { parErr = 0. ; }
       sprintf( varname, "passPar_%s", NP_name ) ;
       rrv_passPar = new RooRealVar( varname, varname, parVal, lowerLimit, upperLimit ) ;
       rrv_passPar -> setError( parErr ) ;
       rrv_passPar->setConstant( kFALSE ) ;
       printf(" floating nuisance parameter: %s = %g +/- %g, [%g, %g]\n",
           varname, parVal, parErr, lowerLimit, upperLimit ) ;

       ////// parVal = alpha+beta-2. ;
       ////// parVal = beta+1. ;
       parVal = beta-1. ;
       lowerLimit = parVal - 6*sqrt(parVal) ;
       if ( lowerLimit <= 0. ) { lowerLimit = 1e-5 ; }
       upperLimit = parVal + 6*sqrt(parVal) ;
       if ( parVal > 0. ) { parErr = sqrt( parVal ) ; } else { parErr = 0. ; }
       sprintf( varname, "failPar_%s", NP_name ) ;
       rrv_failPar = new RooRealVar( varname, varname, parVal, lowerLimit, upperLimit ) ;
       rrv_failPar -> setError( parErr ) ;
       rrv_failPar->setConstant( kFALSE ) ;
       printf(" floating nuisance parameter: %s = %g +/- %g, [%g, %g]\n",
           varname, parVal, parErr, lowerLimit, upperLimit ) ;

       allNuisances -> add( *rrv_passPar ) ;
       allNuisances -> add( *rrv_failPar ) ;


       sprintf( varname, "passConstraint_%s", NP_name ) ;
       passConstraint = new RooPoisson( varname, varname, *rrv_passObs, *rrv_passPar ) ;
       printf("  Created constraint PDF : %s, val = %g, logval = %g\n",
         varname, passConstraint->getVal(), passConstraint->getLogVal() ) ;

       sprintf( varname, "failConstraint_%s", NP_name ) ;
       failConstraint = new RooPoisson( varname, varname, *rrv_failObs, *rrv_failPar ) ;
       printf("  Created constraint PDF : %s, val = %g, logval = %g\n",
         varname, failConstraint->getVal(), failConstraint->getLogVal() ) ;

       allNuisancePdfs -> add( *passConstraint ) ;
       allNuisancePdfs -> add( *failConstraint ) ;

       sprintf( varname, "%s_passPlusFail", NP_name ) ;
       RooAddition* passPlusFail = new RooAddition( varname, varname, RooArgSet(*rrv_passPar, *rrv_failPar) ) ;


       sprintf( varname, "%s", NP_name ) ;
       RooAbsReal* rar = new RooRatio( varname, varname, *rrv_passPar, *passPlusFail ) ;
       printf("  Created nuisance parameter %s : val = %g\n", varname, rar -> getVal() ) ;

       return rar ;

    } // makeBetaConstraint.


   //==============================================================================================================




    RooAbsReal* ra2bRoostatsClass3D_2::makeCorrelatedBetaPrimeConstraint( const char* NP_name, double NP_val, double NP_err, const char* NP_base_name ) {

       double alpha, beta ;
       char varname[1000] ;

       betaPrimeModeTransform( NP_val, NP_err, alpha, beta ) ;




       sprintf( varname, "failPar_%s", NP_base_name ) ;
       RooAbsReal* baseFailPar = (RooAbsReal*) allNuisances -> find( varname ) ;
       if ( baseFailPar != 0x0 ) {
          printf("  base fail parameter : %s = %g\n", baseFailPar->GetName(), baseFailPar->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_2::makeCorrelatedBetaPrimeConstraint: can't find base fail parameter.\n\n") ;
          return 0x0 ;
       }

       sprintf( varname, "passPar_%s", NP_base_name ) ;
       RooAbsReal* basePassPar = (RooAbsReal*) allNuisances -> find( varname ) ;
       if ( basePassPar != 0x0 ) {
          printf("  base pass parameter : %s = %g\n", basePassPar->GetName(), basePassPar->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_2::makeCorrelatedBetaPrimeConstraint: can't find base pass parameter.\n\n") ;
          return 0x0 ;
       }





       sprintf( varname, "failObs_%s", NP_base_name ) ;
       RooAbsReal* baseFailObs = (RooAbsReal*) globalObservables -> find( varname ) ;
       if ( baseFailObs != 0x0 ) {
          printf("  base fail obs : %s = %g\n", baseFailObs->GetName(), baseFailObs->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_2::makeCorrelatedBetaPrimeConstraint: can't find base fail obs.\n\n") ;
          return 0x0 ;
       }

       sprintf( varname, "passObs_%s", NP_base_name ) ;
       RooAbsReal* basePassObs = (RooAbsReal*) globalObservables -> find( varname ) ;
       if ( basePassObs != 0x0 ) {
          printf("  base pass obs : %s = %g\n", basePassObs->GetName(), basePassObs->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_2::makeCorrelatedBetaPrimeConstraint: can't find base pass obs.\n\n") ;
          return 0x0 ;
       }





       sprintf( varname, "passScale_%s", NP_name ) ;
       RooRealVar* passScale = new RooRealVar( varname, varname, (alpha-1)/(basePassObs->getVal()) ) ;
       passScale->setConstant(kTRUE) ;

       sprintf( varname, "failScale_%s", NP_name ) ;
       RooRealVar* failScale = new RooRealVar( varname, varname, (beta-1)/(baseFailObs->getVal()) ) ;
       failScale->setConstant(kTRUE) ;




       sprintf( varname, "passPar_%s", NP_name ) ;
       RooProduct* passPar = new RooProduct( varname, varname, RooArgSet( *basePassPar, *passScale ) ) ;

       sprintf( varname, "failPar_%s", NP_name ) ;
       RooProduct* failPar = new RooProduct( varname, varname, RooArgSet( *baseFailPar, *failScale ) ) ;






       sprintf( varname, "%s", NP_name ) ;
       RooAbsReal* rar = new RooRatio( varname, varname, *passPar, *failPar ) ;
       printf("  Created nuisance parameter %s : val = %g\n", varname, rar -> getVal() ) ;

       return rar ;

    } // makeCorrelatedBetaPrimeConstraint.


   //==============================================================================================================


    RooAbsReal* ra2bRoostatsClass3D_2::makeCorrelatedBetaConstraint( const char* NP_name, double NP_val, double NP_err, const char* NP_base_name ) {

       double alpha, beta ;
       char varname[1000] ;

       betaModeTransform( NP_val, NP_err, alpha, beta ) ;




       sprintf( varname, "failPar_%s", NP_base_name ) ;
       RooAbsReal* baseFailPar = (RooAbsReal*) allNuisances -> find( varname ) ;
       if ( baseFailPar != 0x0 ) {
          printf("  base fail parameter : %s = %g\n", baseFailPar->GetName(), baseFailPar->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_2::makeCorrelatedBetaPrimeConstraint: can't find base fail parameter.\n\n") ;
          return 0x0 ;
       }

       sprintf( varname, "passPar_%s", NP_base_name ) ;
       RooAbsReal* basePassPar = (RooAbsReal*) allNuisances -> find( varname ) ;
       if ( basePassPar != 0x0 ) {
          printf("  base pass parameter : %s = %g\n", basePassPar->GetName(), basePassPar->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_2::makeCorrelatedBetaPrimeConstraint: can't find base pass parameter.\n\n") ;
          return 0x0 ;
       }





       sprintf( varname, "failObs_%s", NP_base_name ) ;
       RooAbsReal* baseFailObs = (RooAbsReal*) globalObservables -> find( varname ) ;
       if ( baseFailObs != 0x0 ) {
          printf("  base fail obs : %s = %g\n", baseFailObs->GetName(), baseFailObs->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_2::makeCorrelatedBetaPrimeConstraint: can't find base fail obs.\n\n") ;
          return 0x0 ;
       }

       sprintf( varname, "passObs_%s", NP_base_name ) ;
       RooAbsReal* basePassObs = (RooAbsReal*) globalObservables -> find( varname ) ;
       if ( basePassObs != 0x0 ) {
          printf("  base pass obs : %s = %g\n", basePassObs->GetName(), basePassObs->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_2::makeCorrelatedBetaPrimeConstraint: can't find base pass obs.\n\n") ;
          return 0x0 ;
       }





       sprintf( varname, "passScale_%s", NP_name ) ;
       RooRealVar* passScale = new RooRealVar( varname, varname, (alpha-1)/(basePassObs->getVal()) ) ;
       passScale->setConstant(kTRUE) ;

       sprintf( varname, "failScale_%s", NP_name ) ;
       RooRealVar* failScale = new RooRealVar( varname, varname, (beta-1)/(baseFailObs->getVal()) ) ;
       failScale->setConstant(kTRUE) ;




       sprintf( varname, "passPar_%s", NP_name ) ;
       RooProduct* passPar = new RooProduct( varname, varname, RooArgSet( *basePassPar, *passScale ) ) ;

       sprintf( varname, "failPar_%s", NP_name ) ;
       RooProduct* failPar = new RooProduct( varname, varname, RooArgSet( *baseFailPar, *failScale ) ) ;






       sprintf( varname, "%s_passPlusFail", NP_name ) ;
       RooAddition* passPlusFail = new RooAddition( varname, varname, RooArgSet(*passPar, *failPar) ) ;


       sprintf( varname, "%s", NP_name ) ;
       RooAbsReal* rar = new RooRatio( varname, varname, *passPar, *passPlusFail ) ;
       printf("  Created nuisance parameter %s : val = %g\n", varname, rar -> getVal() ) ;

       return rar ;

    } // makeCorrelatedBetaPrimeConstraint.


   //==============================================================================================================


