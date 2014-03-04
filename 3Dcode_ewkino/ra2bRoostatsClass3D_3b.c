#include "ra2bRoostatsClass3D_3b.h"

//
//  Version for the ewkino inclusive analysis.
//  
//  Four components are used:
//
//  1) signal
//  2) ttwj
//  3) QCD (taken from MC)
//  4) Z -> inv. (shapes taken from MC), but it's possible to 
//     constrain normalization with VL counts
//


#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>

#include "TROOT.h"
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
#include "TSystem.h"

#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooTrace.h"
#include "RooUniform.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"

#include "RooRatio.h"
#include "RooPosDefCorrGauss.h"
#include "RooBetaPdf.h"
#include "betaPrimeConstraint.c"
#include "RooProdPdfLogSum.h"

  using namespace RooFit ;
  using namespace RooStats ;


  //=====================================================================================================


   ra2bRoostatsClass3D_3b::ra2bRoostatsClass3D_3b() {

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

   ra2bRoostatsClass3D_3b::~ra2bRoostatsClass3D_3b() { }



  //===================================================================

    bool ra2bRoostatsClass3D_3b::initialize( const char* infile ,
					     const char* inputScanFile,
					     double m0, double m12, bool isT1bbbb, double t1bbbbXsec,
					     const char* inputSusy_deff_dbtageff_file,
                                             const char* inputSusy_deff_dbtageff_lightflavor_file,
					     const char* wsrootfilename,
					     const char* blindBinsList,
					     bool constrainBjetShape,
					     bool floatSLSigRatios,
					     int qcdModelIndex,
					     bool floatQcdShape,
					     int znnModelIndex,
					     const char* systFile1,
					     const char* pdf_syst_file,
					     const char* isr_syst_file,
					     const char* wjets_xsec_shapesyst_file,
                                             const char* singletop_xsec_shapesyst_file
					     ) {




      bool useBeta(false) ;
      bool useLognormal(true) ;
      
      
      // check that the qcd and znnModelIndex make sense

      if ( qcdModelIndex < 1 || qcdModelIndex > 2 ) {
	cout << "qcdModelIndex = " << qcdModelIndex << " is not implemented. Exiting ... " << endl ;
	return false ;
      }

      if ( znnModelIndex < 1 || znnModelIndex > 2 ) {
	cout << "znnModelIndex = " << znnModelIndex << " is not implemented. Exiting ... " << endl ;
	return false ;
      }

      // for now ... 
      if ( nBinsBtag > 3 ) {
	cout << "As of now nBinsBtag > 3 is not implemented. Exiting ... " << endl ;
      }


      //-- add the possibility to ignore particular bins in the analysis
      //   (for now don't ignore anything)

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            ignoreBin[mbi][hbi] = false ;
         } // hbi
      } // mbi
      //ignoreBin[nBinsMET-1][0] = true ;

      printf("\n\n *** Ignoring these bins in the analysis.\n\n") ;
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            if ( ignoreBin[mbi][hbi] ) {
               printf("  MET %d, HT %d\n", mbi+1, hbi+1 ) ;
            }
         } // hbi
      } // mbi
      printf("\n\n\n") ;




      //--- read in blind bins, if file given.

      bool blind0lepBin[nBinsMET][nBinsHT][nBinsBtag] ;
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
	for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
	  for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {
	    blind0lepBin[mbi][hbi][bbi] = false ;
	  } // bbi.
	} // hbi.
      } // mbi.
      
      bool blindStudy(false) ;
      
      if ( strcmp( blindBinsList, "null" ) != 0 ) {
	char command[1000] ;
	sprintf( command, "ls %s >& /dev/null", blindBinsList ) ;
	int returnstat = gSystem->Exec( command ) ;
	if ( returnstat != 0 ) {
	  printf("\n\n *** Given non-null blindBinsList file, but it doesn't exist: %s\n\n", blindBinsList ) ;
	  return false ;
	}
	FILE* bbl_file ;
	if ( (bbl_file = fopen( blindBinsList, "r"))==NULL ) {
	  printf("\n\n *** Problem opening blindBinsList file: %s\n\n", blindBinsList ) ;
	  return false ;
	}
	blindStudy = true ;
	printf("\n\n Reading blindBinsList file: %s\n\n", blindBinsList ) ;
	while ( ! feof(bbl_file) ) {
	  char binlabel[1000] ;
	  fscanf( bbl_file, "%s", binlabel ) ;
	  if ( feof(bbl_file) ) break ;
	  int bmb(-1), bhb(-1), bbb(-1) ;
	  int rv = sscanf( binlabel, "N_0lep_M%d_H%d_%db", &bmb, &bhb, &bbb ) ;
	  if ( rv == 0 ) {
	    printf("\n\n *** bad bin label format: %s\n\n", binlabel ) ;
	    return false ;
	  }
	  if ( bmb < 1 || bhb < 1 || bbb < 1 || bmb > nBinsMET || bhb > nBinsHT || bbb > nBinsBtag ) {
	    printf("\n\n *** bin index out of range: met=%d, ht=%d, nb=%d\n\n", bmb, bhb, bbb ) ;
	    return false ;
	  }
	  blind0lepBin[bmb-1][bhb-1][bbb-1] = true ;
	  printf("   bin label %s : will blind met=%d, ht=%d, nb=%d 0lep bin.\n", binlabel, bmb, bhb, bbb ) ;
	} // still reading?
      } // process blindBinsList file?





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
      int N_1lepSig[nBinsMET][nBinsHT][nBinsBtag] ;
      int N_1lep[nBinsMET][nBinsHT][nBinsBtag] ;

      int N_ZeeVL[nBinsMET][nBinsHT] ;
      int N_ZmmVL[nBinsMET][nBinsHT] ;


      // MC inputs:

      float NQCDmc_0lep[nBinsMET][nBinsHT][nBinsBtag] ;
      float NZnnmc_0lep[nBinsMET][nBinsHT][nBinsBtag] ;

      
      // other parameters:

      float sf_ttwj[nBinsMET][nBinsHT][nBinsBtag] ;
      float sf_ttwj_err[nBinsMET][nBinsHT][nBinsBtag] ;
      float sf_ttwj_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
      float sf_ttwj_slSig_err[nBinsMET][nBinsHT][nBinsBtag] ;

      float trigeff_0L[nBinsMET][nBinsHT];
      float trigefferr_0L[nBinsMET][nBinsHT];
      float trigeff_1L[nBinsMET][nBinsHT];
      float trigefferr_1L[nBinsMET][nBinsHT];

      int   n_global_uncertainties(0) ;
      float global_uncertainty_val[100] ;
      char  global_uncertainty_name[100][100] ;


      //--- read in description lines.
      printf("\n\n") ;
      char c(0) ;
      while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
      printf("\n\n") ;


      char label[1000] ;
      
      // read all observables and other parameters
      // first create the TString stubs needed to read the text file:
     
      TString sMbins[nBinsMET];
      TString sHbins[nBinsHT];
      TString sBbins[4] = {"_1b","_2b","_3b","_4b"};
      
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

            if ( ignoreBin[i][j] ) {
               N_0lep[i][j][k] = 0. ;
               continue ;
            }

	    cout << inPar << " = " << N_0lep[i][j][k] << endl ;

	  }
	}
      }


      // 1lepSig observables

      for ( int i = 0 ; i < nBinsMET ; i++ ) {
	for ( int j = 0 ; j < nBinsHT ; j++ ) {
	  for ( int k = 0 ; k < nBinsBtag ; k++ ) {

	    TString inPar = "N_1lepSig"+sMbins[i]+sHbins[j]+sBbins[k] ;
            float value ;
	    fscanf( infp, "%s %g", label, &value ) ;
            N_1lepSig[i][j][k] = TMath::Nint( value ) ;

	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }

            if ( ignoreBin[i][j] ) {
               N_1lepSig[i][j][k] = 0. ;
               continue ;
            }

	    cout << inPar << " = " << N_1lepSig[i][j][k] << endl ;

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

            if ( ignoreBin[i][j] ) {
               N_1lep[i][j][k] = 0. ;
               continue ;
            }

	    cout << inPar << " = " << N_1lep[i][j][k] << endl ;

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


      // sf TTWJ slSig

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {

	    TString inPar = "sf_ttwj_slSig"+sMbins[i]+sHbins[j]+sBbins[k] ;
	    fscanf( infp, "%s %g", label, &sf_ttwj_slSig[i][j][k] ) ;
	
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << sf_ttwj_slSig[i][j][k] << endl ;	  
	    
	    inPar += "_err" ;
	    fscanf( infp, "%s %g", label, &sf_ttwj_slSig_err[i][j][k] ) ;
	    cout << inPar << " = " << sf_ttwj_slSig_err[i][j][k] << endl ;	  

	  }
	}
      }      


      // Next get 3b/2b and 4b/2b ratios for ttwj

      float      ttwj_3b2b_ratio[nBinsMET][nBinsHT] ;
      float ttwjSlSig_3b2b_ratio[nBinsMET][nBinsHT] ;

      float      ttwj_3b2b_ratio_err[nBinsMET][nBinsHT] ;
      float ttwjSlSig_3b2b_ratio_err[nBinsMET][nBinsHT] ;

      float      ttwj_4b2b_ratio[nBinsMET][nBinsHT] ;
      float ttwjSlSig_4b2b_ratio[nBinsMET][nBinsHT] ;

      float      ttwj_4b2b_ratio_err[nBinsMET][nBinsHT] ;
      float ttwjSlSig_4b2b_ratio_err[nBinsMET][nBinsHT] ;


      // 0lep sample

      if ( nBinsBtag > 2 ) {

	printf("\n\n") ;
	for (int i = 0 ; i < nBinsMET ; i++) {
	  for (int j = 0 ; j < nBinsHT ; j++) {
	    char expectedlabel[1000] ;                                                                                                          
	    sprintf( expectedlabel, "ttwj_mc_3bover2b_ratio_M%d_H%d", i+1, j+1 ) ;                                                      
	    char vallabel[1000] ;                                                                                                               
	    char errlabel[1000] ;                                                                                                               
	    fscanf( infp, "%s %g", vallabel, &ttwj_3b2b_ratio[i][j] ) ;                                                                 
	    if ( strcmp( vallabel, expectedlabel ) != 0 ) { mismatchErr( label, expectedlabel ) ; return false ; }                              
	    fscanf( infp, "%s %g", errlabel, &ttwj_3b2b_ratio_err[i][j] ) ;                                                             
	    printf(" ttwj 3b/2b ratio : %s  %6.3f +/- %5.3f\n", vallabel, ttwj_3b2b_ratio[i][j], ttwj_3b2b_ratio_err[i][j] ) ;   
	  } // j
	} // i
	printf("\n\n") ;
	
	if ( nBinsBtag > 3 ) {

	  printf("\n\n") ;
	  for (int i = 0 ; i < nBinsMET ; i++) {
	    for (int j = 0 ; j < nBinsHT ; j++) {
	      char expectedlabel[1000] ;                                                                                                          
	      sprintf( expectedlabel, "ttwj_mc_4bover2b_ratio_M%d_H%d", i+1, j+1 ) ;                                                      
	      char vallabel[1000] ;                                                                                                               
	      char errlabel[1000] ;                                                                                                               
	      fscanf( infp, "%s %g", vallabel, &ttwj_4b2b_ratio[i][j] ) ;                                                                 
	      if ( strcmp( vallabel, expectedlabel ) != 0 ) { mismatchErr( label, expectedlabel ) ; return false ; }                              
	      fscanf( infp, "%s %g", errlabel, &ttwj_4b2b_ratio_err[i][j] ) ;                                                             
	      printf(" ttwj 4b/2b ratio : %s  %6.3f +/- %5.3f\n", vallabel, ttwj_4b2b_ratio[i][j], ttwj_4b2b_ratio_err[i][j] ) ;   
	    } // j
	  } // i
	  printf("\n\n") ;
	  
	}
      }


      // SlSig sample

      if ( nBinsBtag > 2 ) {

	printf("\n\n") ;
	for (int i = 0 ; i < nBinsMET ; i++) {
	  for (int j = 0 ; j < nBinsHT ; j++) {
	    char expectedlabel[1000] ;                                                                                                          
	    sprintf( expectedlabel, "ttwjSlSig_mc_3bover2b_ratio_M%d_H%d", i+1, j+1 ) ;                                                      
	    char vallabel[1000] ;                                                                                                               
	    char errlabel[1000] ;                                                                                                               
	    fscanf( infp, "%s %g", vallabel, &ttwjSlSig_3b2b_ratio[i][j] ) ;                                                                 
	    if ( strcmp( vallabel, expectedlabel ) != 0 ) { mismatchErr( label, expectedlabel ) ; return false ; }                              
	    fscanf( infp, "%s %g", errlabel, &ttwjSlSig_3b2b_ratio_err[i][j] ) ;                                                             
	    printf(" ttwjSlSig 3b/2b ratio : %s  %6.3f +/- %5.3f\n", vallabel, ttwjSlSig_3b2b_ratio[i][j], ttwjSlSig_3b2b_ratio_err[i][j] ) ;   
	  } // j
	} // i
	printf("\n\n") ;
	
	if ( nBinsBtag > 3 ) {

	  printf("\n\n") ;
	  for (int i = 0 ; i < nBinsMET ; i++) {
	    for (int j = 0 ; j < nBinsHT ; j++) {
	      char expectedlabel[1000] ;                                                                                                          
	      sprintf( expectedlabel, "ttwjSlSig_mc_4bover2b_ratio_M%d_H%d", i+1, j+1 ) ;                                                      
	      char vallabel[1000] ;                                                                                                               
	      char errlabel[1000] ;                                                                                                               
	      fscanf( infp, "%s %g", vallabel, &ttwjSlSig_4b2b_ratio[i][j] ) ;                                                                 
	      if ( strcmp( vallabel, expectedlabel ) != 0 ) { mismatchErr( label, expectedlabel ) ; return false ; }                              
	      fscanf( infp, "%s %g", errlabel, &ttwjSlSig_4b2b_ratio_err[i][j] ) ;                                                             
	      printf(" ttwjSlSig 4b/2b ratio : %s  %6.3f +/- %5.3f\n", vallabel, ttwjSlSig_4b2b_ratio[i][j], ttwjSlSig_4b2b_ratio_err[i][j] ) ;   
	    } // j
	  } // i
	  printf("\n\n") ;
	  
	}
      }


      // NQCDmc_0lep

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {
	    
	    TString inPar = "N_QCDmc_0lep"+sMbins[i]+sHbins[j]+sBbins[k] ;
	    fscanf( infp, "%s %g", label, &NQCDmc_0lep[i][j][k] ) ;
	    
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << NQCDmc_0lep[i][j][k] << endl ;	  

	  }
	}
      }
      

      // Next get 1b/2b and 3b/2b ratios for qcd

      float qcd_2b1b_ratio ;
      float qcd_2b1b_ratio_err ;
      float qcd_3b1b_ratio ;
      float qcd_3b1b_ratio_err ;
      float QCD_SF ;
      float QCD_SF_err ;

      fscanf( infp, "%s %g", label, &qcd_2b1b_ratio ) ;     cout << "qcd_2b1b_ratio    " << " = " << qcd_2b1b_ratio << endl ;
      fscanf( infp, "%s %g", label, &qcd_2b1b_ratio_err ) ; cout << "qcd_2b1b_ratio_err" << " = " << qcd_2b1b_ratio_err << endl ;
      fscanf( infp, "%s %g", label, &qcd_3b1b_ratio ) ;     cout << "qcd_3b1b_ratio    " << " = " << qcd_3b1b_ratio << endl ;
      fscanf( infp, "%s %g", label, &qcd_3b1b_ratio_err ) ; cout << "qcd_3b1b_ratio_err" << " = " << qcd_3b1b_ratio_err << endl ;
      fscanf( infp, "%s %g", label, &QCD_SF ) ;     cout << "QCD_SF    " << " = " << QCD_SF << endl ;
      fscanf( infp, "%s %g", label, &QCD_SF_err ) ; cout << "QCD_SF_err" << " = " << QCD_SF_err << endl ;


      // NZnnmc_0lep

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {
	    
	    TString inPar = "N_ZNNmc_0lep"+sMbins[i]+sHbins[j]+sBbins[k] ;
	    fscanf( infp, "%s %g", label, &NZnnmc_0lep[i][j][k] ) ;
	    
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << NZnnmc_0lep[i][j][k] << endl ;	  

	  }
	}
      }

      
      // NZeeVL_0lep

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	    
	  TString inPar = "N_ZeeVL_0lep"+sMbins[i]+sHbins[j] ;
	  fscanf( infp, "%s %d", label, &N_ZeeVL[i][j] ) ;
	  
	  if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	  cout << inPar << " = " << N_ZeeVL[i][j] << endl ;	  
	  
	}
      }

      
      // NZmmVL_0lep

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	    
	  TString inPar = "N_ZmmVL_0lep"+sMbins[i]+sHbins[j] ;
	  fscanf( infp, "%s %d", label, &N_ZmmVL[i][j] ) ;
	  
	  if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	  cout << inPar << " = " << N_ZmmVL[i][j] << endl ;	  
	  
	}
      }


      // Next get 1b/2b and 3b/2b ratios for Znn

      float Znn_2b1b_ratio ;
      float Znn_2b1b_ratio_err ;
      float Znn_3b1b_ratio ;
      float Znn_3b1b_ratio_err ;
      float Znn_SF ;
      float Znn_SF_err ;

      fscanf( infp, "%s %g", label, &Znn_2b1b_ratio ) ;     cout << "Znn_2b1b_ratio    " << " = " << Znn_2b1b_ratio << endl ;
      fscanf( infp, "%s %g", label, &Znn_2b1b_ratio_err ) ; cout << "Znn_2b1b_ratio_err" << " = " << Znn_2b1b_ratio_err << endl ;
      fscanf( infp, "%s %g", label, &Znn_3b1b_ratio ) ;     cout << "Znn_3b1b_ratio    " << " = " << Znn_3b1b_ratio << endl ;
      fscanf( infp, "%s %g", label, &Znn_3b1b_ratio_err ) ; cout << "Znn_3b1b_ratio_err" << " = " << Znn_3b1b_ratio_err << endl ;
      fscanf( infp, "%s %g", label, &Znn_SF ) ;     cout << "Znn_SF    " << " = " << Znn_SF << endl ;
      fscanf( infp, "%s %g", label, &Znn_SF_err ) ; cout << "Znn_SF_err" << " = " << Znn_SF_err << endl ;



      // next get trigger efficiency values
      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	    
	    TString inPar = "trigeff_val_0L"+sMbins[i]+sHbins[j];
	    fscanf( infp, "%s %g", label, &trigeff_0L[i][j] ) ;
	    
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << trigeff_0L[i][j] << endl ;	  

	    inPar = "trigeff_err_0L"+sMbins[i]+sHbins[j];
	    fscanf( infp, "%s %g", label, &trigefferr_0L[i][j] ) ;
	    
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << trigefferr_0L[i][j] << endl ;	  

        }
      }

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	    
	    TString inPar = "trigeff_val_1L"+sMbins[i]+sHbins[j];
	    fscanf( infp, "%s %g", label, &trigeff_1L[i][j] ) ;
	    
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << trigeff_1L[i][j] << endl ;	  

	    inPar = "trigeff_err_1L"+sMbins[i]+sHbins[j];
	    fscanf( infp, "%s %g", label, &trigefferr_1L[i][j] ) ;
	    
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << trigefferr_1L[i][j] << endl ;	  

        }
      }





      //--- Oct 31, 2012: read in global uncertainties.

      char target_label[100] ;

      sprintf( target_label, "GU_luminosity" ) ;
      fscanf( infp, "%s %g", label, &global_uncertainty_val[n_global_uncertainties] ) ;
      if ( strcmp( label, target_label ) != 0 ) { printf("\n\n *** expecting %s, found %s\n\n", target_label, label ) ; return false ; }
      sprintf( global_uncertainty_name[n_global_uncertainties], "%s", label ) ;
      n_global_uncertainties++ ;

      sprintf( target_label, "GU_metcleaning" ) ;
      fscanf( infp, "%s %g", label, &global_uncertainty_val[n_global_uncertainties] ) ;
      if ( strcmp( label, target_label ) != 0 ) { printf("\n\n *** expecting %s, found %s\n\n", target_label, label ) ; return false ; }
      sprintf( global_uncertainty_name[n_global_uncertainties], "%s", label ) ;
      n_global_uncertainties++ ;

      sprintf( target_label, "GU_JER" ) ;
      fscanf( infp, "%s %g", label, &global_uncertainty_val[n_global_uncertainties] ) ;
      if ( strcmp( label, target_label ) != 0 ) { printf("\n\n *** expecting %s, found %s\n\n", target_label, label ) ; return false ; }
      sprintf( global_uncertainty_name[n_global_uncertainties], "%s", label ) ;
      n_global_uncertainties++ ;

      sprintf( target_label, "GU_unclMET" ) ;
      fscanf( infp, "%s %g", label, &global_uncertainty_val[n_global_uncertainties] ) ;
      if ( strcmp( label, target_label ) != 0 ) { printf("\n\n *** expecting %s, found %s\n\n", target_label, label ) ; return false ; }
      sprintf( global_uncertainty_name[n_global_uncertainties], "%s", label ) ;
      n_global_uncertainties++ ;
 



      //---- calculations for determining initial values for floating parameters.

      double initialval_qcd[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_znn[nBinsMET][nBinsHT][nBinsBtag] ;

      double initialval_ttwj_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_ttwj_sl[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_ttwj[nBinsMET][nBinsHT][nBinsBtag] ;


      double initialguess_ttwj_0lep1lep_ratio = 1.40 ;


      // set initial values of 1lepSig just by rescaling 1lep
      // take 1lepSig/1lep ratio from bin (1,1,1)

      float ratio_slSigsl = (double)N_1lepSig[0][0][0] / N_1lep[0][0][0] ;


      printf("\n\n initial guess for ttwj 0lep/1lep ratio : %5.2f\n\n", initialguess_ttwj_0lep1lep_ratio ) ;

      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {

          if ( ignoreBin[i][j] ) continue ;

          for (int k = 0 ; k < nBinsBtag ; k++) {     

            // TTWJ stuff.

            initialval_ttwj_sl[i][j][k] = N_1lep[i][j][k] / trigeff_1L[i][j] ;
            initialval_ttwj[i][j][k] = sf_ttwj[i][j][k] * initialguess_ttwj_0lep1lep_ratio * initialval_ttwj_sl[i][j][k] ;
	    initialval_ttwj_slSig[i][j][k] = sf_ttwj_slSig[i][j][k] * initialval_ttwj_sl[i][j][k] * ratio_slSigsl ;

            // QCD and ZNN stuff:

            initialval_qcd[i][j][k] = NQCDmc_0lep[i][j][k] ;
            initialval_znn[i][j][k] = NZnnmc_0lep[i][j][k] ;

          }
        }
      }



      // print out initial values:


      printf("\n\n\n --------- Observables and floating parameter initial values. ------------\n\n") ;

      printf("           |   Nobs   |  Model  |    PDF        ||   ttwj   |   QCD   |   Znn  |\n") ;
      printf("-----------+----------+---------+---------------++----------+---------+--------+\n") ;

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {

          if ( ignoreBin[i][j] ) continue ;

	  for (int k = 0 ; k < nBinsBtag ; k++) {     

            double model_0lep = trigeff_0L[i][j] * initialval_qcd[i][j][k] + trigeff_1L[i][j] * (initialval_ttwj[i][j][k]+initialval_znn[i][j][k]) ;
            double pdf_0lep = TMath::PoissonI( N_0lep[i][j][k], model_0lep ) ;
            char warning0lep[4] ;
            sprintf( warning0lep,"   ") ;
            if ( pdf_0lep < 0.0001   ) { sprintf( warning0lep, "*  ") ; }
            if ( pdf_0lep < 0.00001  ) { sprintf( warning0lep, "** ") ; }
            if ( pdf_0lep < 0.000001 ) { sprintf( warning0lep, "***") ; }
	    cout << " MET bin " << i+1 << ", HT bin " << j+1 << ", Btag bin " << k+1 << endl;
	    printf(" 0-lep     | %6d   | %7.1f | %8.6f %4s ||  %7.1f | %7.1f | %7.1f | \n", N_0lep[i][j][k], model_0lep, pdf_0lep, warning0lep,
		   trigeff_1L[i][j] * initialval_ttwj[i][j][k], trigeff_0L[i][j] * initialval_qcd[i][j][k], trigeff_1L[i][j] * initialval_znn[i][j][k] ) ;
            printf("-----------+----------+---------+---------------++----------+---------+---------+\n") ;

	  }
	}
      }






      printf(" --- Creating workspace.\n" ) ; cout << flush ;

      RooWorkspace workspace ("ws") ;
      workspace.autoImportClassCode(true);


      globalObservables      = new RooArgSet("globalObservables");
      allNuisances           = new RooArgSet("allNuisances");
      allNuisancePdfs        = new RooArgSet("allNuisancePdfs");
      observedParametersList = new RooArgSet("observables") ;






      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining observables and parameters of the likelihood.\n" ) ;


      rv_mu_susy_all0lep = new RooRealVar( "mu_susy_all0lep", "mu_susy_all0lep", 0., 100000. ) ;
      rv_mu_susy_all0lep->setVal( 1. ) ;
      rv_mu_susy_all0lep->setConstant( kTRUE ) ;


      rv_mean_qcd_sf = new RooRealVar( "mean_qcd_sf", "mean_qcd_sf", 0., 10. );			  
      rv_mean_qcd_sf->setVal( 1.0 ) ;		   // for now hardcode 1.0, this can be set on the input file						  
      rv_mean_qcd_sf->setConstant( kTRUE ) ;							  
      rv_width_qcd_sf = new RooRealVar( "width_qcd_sf", "width_qcd_sf", 0., 10. );  		  
      rv_width_qcd_sf->setVal( 1.00 ) ; 	   // for now hardcode 1.0, this can be set on the input file						  
      rv_width_qcd_sf->setConstant( kTRUE ) ;							  


      rv_mean_znn_sf = new RooRealVar( "mean_znn_sf", "mean_znn_sf", 0., 10. );			  
      rv_mean_znn_sf->setVal( 1.0 ) ;		   // for now hardcode 1.0, this can be set on the input file						  
      rv_mean_znn_sf->setConstant( kTRUE ) ;							  
      rv_width_znn_sf = new RooRealVar( "width_znn_sf", "width_znn_sf", 0., 10. );  		  
      rv_width_znn_sf->setVal( 1.00 ) ; 	   // for now hardcode 1.0, this can be set on the input file						  
      rv_width_znn_sf->setConstant( kTRUE ) ;							  


            
      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {

          if ( ignoreBin[i][j] ) continue ;

	  for (int k = 0 ; k < nBinsBtag ; k++) {     

	    TString zlString    = "N_0lep";
	    TString slSigString = "N_1lepSig";
	    TString slString    = "N_1lep";

	    zlString    += sMbins[i]+sHbins[j]+sBbins[k] ;
	    slSigString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    slString    += sMbins[i]+sHbins[j]+sBbins[k] ;

	    rv_0lep[i][j][k]    = new RooRealVar( zlString,    zlString,    0., 1000000.);
	    rv_1lepSig[i][j][k] = new RooRealVar( slSigString, slSigString, 0., 1000000.);
	    rv_1lep[i][j][k]    = new RooRealVar( slString,    slString,    0., 1000000.);

	    // set values:
	    rv_0lep[i][j][k]->setVal( N_0lep[i][j][k] );
	    rv_1lepSig[i][j][k]->setVal( N_1lepSig[i][j][k] );
	    rv_1lep[i][j][k]->setVal( N_1lep[i][j][k] );

	    observedParametersList -> add( *rv_0lep[i][j][k] ) ;
	    //observedParametersList -> add( *rv_1lepSig[i][j][k] ) ;
	    observedParametersList -> add( *rv_1lep[i][j][k] ) ;
	    
	    // parameters of the likelihood:
	    TString muTtString      = "mu_ttwj";
	    TString muTtSlSigString = "mu_ttwj_slSig";
	    TString muTtSlString    = "mu_ttwj_sl";
	    TString muQcdString     = "mu_qcd";
	    TString muZnnString     = "mu_znn";

	    TString muSusMcString      = "mu_susymc";
	    TString muSusSlSigMcString = "mu_susymc_slSig";
	    TString muSusSlMcString    = "mu_susymc_sl";

	    TString MEffSfString     = "mean_eff_sf";
	    TString WEffSfString     = "width_eff_sf";

	    TString MEffSfSlSigString = "mean_eff_sf_slSig";
	    TString WEffSfSlSigString = "width_eff_sf_slSig";

	    TString MEffSfSlString   = "mean_eff_sf_sl";
	    TString WEffSfSlString   = "width_eff_sf_sl";


	    muTtString      += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muTtSlSigString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muTtSlString    += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muQcdString     += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muZnnString     += sMbins[i]+sHbins[j]+sBbins[k] ;

	    muSusMcString      += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muSusSlSigMcString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    muSusSlMcString    += sMbins[i]+sHbins[j]+sBbins[k] ;

	    MEffSfString   += sMbins[i]+sHbins[j]+sBbins[k] ;
	    WEffSfString   += sMbins[i]+sHbins[j]+sBbins[k] ;

	    MEffSfSlSigString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    WEffSfSlSigString += sMbins[i]+sHbins[j]+sBbins[k] ;

	    MEffSfSlString   += sMbins[i]+sHbins[j]+sBbins[k] ;
	    WEffSfSlString   += sMbins[i]+sHbins[j]+sBbins[k] ;


            rrv_mu_ttwj_slSig[i][j][k] = new RooRealVar( muTtSlSigString, muTtSlSigString, 0., 100000. ) ;
            rv_mu_ttwj_slSig[i][j][k] = rrv_mu_ttwj_slSig[i][j][k] ;
            rrv_mu_ttwj_slSig[i][j][k]->setVal( initialval_ttwj_slSig[i][j][k] ) ;    // this is a starting value only
            //if ( !(ignoreBin[i][j]) ) { allNuisances -> add( *rv_mu_ttwj_slSig[i][j][k] ) ; } //-- is this correct?

            rrv_mu_ttwj_sl[i][j][k] = new RooRealVar( muTtSlString, muTtSlString, 0., 100000. ) ;
            rv_mu_ttwj_sl[i][j][k] = rrv_mu_ttwj_sl[i][j][k] ;
            rrv_mu_ttwj_sl[i][j][k]->setVal( initialval_ttwj_sl[i][j][k] ) ;    // this is a starting value only
            if ( !(ignoreBin[i][j]) ) { allNuisances -> add( *rv_mu_ttwj_sl[i][j][k] ) ; } //-- is this correct?

	    rrv_mu_qcd[i][j][k] = new RooRealVar( muQcdString, muQcdString, 0., 100000. ) ;
	    rv_mu_qcd[i][j][k] = rrv_mu_qcd[i][j][k];
	    rrv_mu_qcd[i][j][k]->setVal( initialval_qcd[i][j][k] ) ;          
	    rrv_mu_qcd[i][j][k]->setConstant(kTRUE);

	    rrv_mu_znn[i][j][k] = new RooRealVar( muZnnString, muZnnString, 0., 100000. ) ;
	    rv_mu_znn[i][j][k] = rrv_mu_znn[i][j][k];
	    rrv_mu_znn[i][j][k]->setVal( initialval_znn[i][j][k] ) ;          
	    rrv_mu_znn[i][j][k]->setConstant(kTRUE);


	    // MC inputs

	    rv_mu_susymc[i][j][k] = new RooRealVar( muSusMcString, muSusMcString, 0., 100000. ) ;
	    rv_mu_susymc[i][j][k]->setVal( 1. ) ;
	    rv_mu_susymc[i][j][k]->setConstant( kTRUE ) ;

	    rv_mu_susymc_slSig[i][j][k] = new RooRealVar( muSusSlSigMcString, muSusSlSigMcString, 0., 100000. ) ;
	    rv_mu_susymc_slSig[i][j][k]->setVal( 1. ) ;
	    rv_mu_susymc_slSig[i][j][k]->setConstant( kTRUE ) ;

	    rv_mu_susymc_sl[i][j][k] = new RooRealVar( muSusSlMcString, muSusSlMcString, 0., 100000. ) ;
	    rv_mu_susymc_sl[i][j][k]->setVal( 1. ) ;
	    rv_mu_susymc_sl[i][j][k]->setConstant( kTRUE ) ;

	    
	    // gaussian constraints

	    rv_mean_eff_sf[i][j][k] = new RooRealVar( MEffSfString, MEffSfString, 0., 10. );
	    rv_mean_eff_sf[i][j][k]->setVal( 1.0 ) ;
	    rv_mean_eff_sf[i][j][k]->setConstant( kTRUE ) ;

	    rv_width_eff_sf[i][j][k] = new RooRealVar( WEffSfString, WEffSfString, 0., 10. );
	    rv_width_eff_sf[i][j][k]->setVal( 0.15 ) ;            // arbitrarily set the uncertainty to 15%
	    rv_width_eff_sf[i][j][k]->setConstant( kTRUE ) ;   


	    rv_mean_eff_sf_slSig[i][j][k] = new RooRealVar( MEffSfSlSigString, MEffSfSlSigString, 0., 10. );
	    rv_mean_eff_sf_slSig[i][j][k]->setVal( 1.0 ) ;
	    rv_mean_eff_sf_slSig[i][j][k]->setConstant( kTRUE ) ;

	    rv_width_eff_sf_slSig[i][j][k] = new RooRealVar( WEffSfSlSigString, WEffSfSlSigString, 0., 10. );
	    rv_width_eff_sf_slSig[i][j][k]->setVal( 0.15 ) ;            // arbitrarily set the uncertainty to 15%
	    rv_width_eff_sf_slSig[i][j][k]->setConstant( kTRUE ) ;   


	    rv_mean_eff_sf_sl[i][j][k] = new RooRealVar( MEffSfSlString, MEffSfSlString, 0., 10. );
	    rv_mean_eff_sf_sl[i][j][k]->setVal( 1.0 ) ;
	    rv_mean_eff_sf_sl[i][j][k]->setConstant( kTRUE ) ;

	    rv_width_eff_sf_sl[i][j][k] = new RooRealVar( WEffSfSlString, WEffSfSlString, 0., 10. );
	    rv_width_eff_sf_sl[i][j][k]->setVal( 0.15 ) ;            // arbitrarily set the uncertainty to 15%
	    rv_width_eff_sf_sl[i][j][k]->setConstant( kTRUE ) ;   

	  }

	}
      }


      //--- New in 3 : floating ttwj_0lep1lep_ratio

      RooRealVar* rv_ttwj_0lep1lep_ratio = new RooRealVar( "ttwj_0lep1lep_ratio", "ttwj_0lep1lep_ratio", 0., 5. ) ;
      rv_ttwj_0lep1lep_ratio -> setVal( initialguess_ttwj_0lep1lep_ratio ) ; // initial guess
      rv_ttwj_0lep1lep_ratio -> setConstant( kFALSE ) ;

      allNuisances -> add( *rv_ttwj_0lep1lep_ratio ) ; //-- is this correct?

      
      //--- New for T1tttt : slSig/sl ratio

      RooRealVar* rv_ttwj_slSigsl_ratio = new RooRealVar( "ttwj_slSigsl_ratio", "ttwj_slSigsl_ratio", 0., 5. ) ;
      rv_ttwj_slSigsl_ratio -> setVal( ratio_slSigsl ) ; // initial guess
      rv_ttwj_slSigsl_ratio -> setConstant( kFALSE ) ;

      //allNuisances -> add( *rv_ttwj_slSigsl_ratio ) ;


      RooRealVar* rv_ttwj_slSigDD_ratio[nBinsMET][nBinsHT] ;

      if ( floatSLSigRatios ) {
	for (int i = 0 ; i < nBinsMET ; i++) {
	  for (int j = 0 ; j < nBinsHT ; j++) {

	    TString SlSigDD_ratioS = "slSigDD_ratio";
	    SlSigDD_ratioS += sMbins[i]+sHbins[j];

	    rv_ttwj_slSigDD_ratio[i][j] = new RooRealVar( SlSigDD_ratioS, SlSigDD_ratioS, 0., 5. ) ;
	    rv_ttwj_slSigDD_ratio[i][j] -> setVal( 1. ) ; // initial guess
	    rv_ttwj_slSigDD_ratio[i][j] -> setConstant( kFALSE ) ;

	    //allNuisances -> add( *rv_ttwj_slSigDD_ratio[i][j] ) ;

	  }
	}
      }


      // float QCD shape (if specified)

      RooRealVar* rv_qcd_shape_ratio[nBinsMET][nBinsHT] ;

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  
	  TString qcdShape_ratioS = "qcdShape_ratio";
	  qcdShape_ratioS += sMbins[i]+sHbins[j];
	  
	  rv_qcd_shape_ratio[i][j] = new RooRealVar( qcdShape_ratioS, qcdShape_ratioS, 0., 5. ) ;
	  rv_qcd_shape_ratio[i][j] -> setVal( 1. ) ; 
	  rv_qcd_shape_ratio[i][j] -> setConstant( kTRUE ) ;

	  if ( floatQcdShape ) {

	    if ( i == 0 && j == 0 ) continue ;    // set ratios wrt lowest (Var1,Var2) bin

	    rv_qcd_shape_ratio[i][j] -> setConstant( kFALSE ) ;
	    allNuisances -> add( *rv_qcd_shape_ratio[i][j] ) ;
	    
	  }

	}
      }

      

      RooDataSet* dsObserved ;
      dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
				  *observedParametersList ) ;
      dsObserved->add( *observedParametersList ) ;

      printf(" --- Importing dataset.\n" ) ; cout << flush ;

      workspace.import(*dsObserved);



      char NP_name[1000] ;





      //++++++++++

      bool ssspOk = setSusyScanPoint( inputScanFile,  m0,  m12 ) ;

      if ( !ssspOk ) {
         printf("\n\n\n *** setSusyScanPoint failed.  I quit.\n\n\n") ;
         return false ;
      }

      RooRealVar* rv_mu_susymc_all0lep = new RooRealVar( "mu_susymc_all0lep", "mu_susymc_all0lep", 0., 100000. ) ;

      double susymc_all0lep(0.) ;
      for (int i = 0 ; i < nBinsMET ; i++) {
         for (int j = 0 ; j < nBinsHT ; j++) {
            if ( ignoreBin[i][j] ) continue ;
            for (int k = 0 ; k < nBinsBtag ; k++) {
               susymc_all0lep += rv_mu_susymc[i][j][k] -> getVal() ;
            }
         }
      }

      rv_mu_susymc_all0lep -> setVal( susymc_all0lep ) ;
      rv_mu_susymc_all0lep -> setConstant( kTRUE ) ;

      //++++++++++

      cout << "\n\n Back from setSusyScanPoint.  Now defining parameters.\n\n" << flush ;

      //--- Systematics and other nuisance parameters

      char formula[1024];


      printf("\n\n Ttwj scale factors.\n\n") ; cout << flush ;

      // TTWJ scale factors

      RooAbsReal* rar_sf_ttwj[nBinsMET][nBinsHT][nBinsBtag];
      RooAbsReal* rar_sf_ttwj_slSig[nBinsMET][nBinsHT][nBinsBtag];


      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          if ( ignoreBin[i][j] ) continue ;
          for (int k = 0 ; k < nBinsBtag ; k++) {


             //--- ttwj


             sprintf( NP_name, "sf_ttwj_M%d_H%d_%db", i+1, j+1, k+1 ) ;
	     if ( ! (blindStudy && blind0lepBin[i][j][k]) ) {
	       if ( useLognormal ) {
		 rar_sf_ttwj[i][j][k] = makeLognormalConstraint( NP_name, sf_ttwj[i][j][k], sf_ttwj_err[i][j][k] ) ;
	       } else {
		 rar_sf_ttwj[i][j][k] = makeGaussianConstraint( NP_name, sf_ttwj[i][j][k], sf_ttwj_err[i][j][k] ) ;
	       }
             } else {
	       rar_sf_ttwj[i][j][k] = makeGaussianConstraint( NP_name, sf_ttwj[i][j][k], 0. ) ; //-- this just creates a constant.
             }



             //--- ttwj slSig

             sprintf( NP_name, "sf_ttwj_slSig_M%d_H%d_%db", i+1, j+1, k+1 ) ;
	     if ( ! (blindStudy && blind0lepBin[i][j][k]) ) {
	       if ( useLognormal ) {
		 rar_sf_ttwj_slSig[i][j][k] = makeLognormalConstraint( NP_name, sf_ttwj_slSig[i][j][k], sf_ttwj_slSig_err[i][j][k] ) ;
	       } else {
		 rar_sf_ttwj_slSig[i][j][k] = makeGaussianConstraint( NP_name, sf_ttwj_slSig[i][j][k], sf_ttwj_slSig_err[i][j][k] ) ;
	       }
             } else {
	       rar_sf_ttwj_slSig[i][j][k] = makeGaussianConstraint( NP_name, sf_ttwj_slSig[i][j][k], 0. ) ; //-- this just creates a constant.
             }

             printf("--------\n") ;

          } // i (met)
             printf("++++++++++++++++\n") ;
        } // j (ht)
             printf("=========================\n") ;
      } // k (nbtag)



      printf("\n\n Ttwj nB ratio factors.\n\n") ; cout << flush ;

      // TTWJ nB ratio factors
      
      RooAbsReal* rar_ttwj_3b2b_ratio[nBinsMET][nBinsHT];
      RooAbsReal* rar_ttwjSlSig_3b2b_ratio[nBinsMET][nBinsHT];
      RooAbsReal* rar_ttwj_4b2b_ratio[nBinsMET][nBinsHT];
      RooAbsReal* rar_ttwjSlSig_4b2b_ratio[nBinsMET][nBinsHT];


      if ( constrainBjetShape && nBinsBtag > 2 ) {

	for (int i = 0 ; i < nBinsMET ; i++) {
	  for (int j = 0 ; j < nBinsHT ; j++) {
	    if ( ignoreBin[i][j] ) continue ;

	    //--- ttwj 0lep

	    sprintf( NP_name, "ttwj_3b2b_ratio_M%d_H%d", i+1, j+1 ) ;
	    if ( ! (blindStudy && blind0lepBin[i][j]) ) {
	      if ( useLognormal ) {
                   rar_ttwj_3b2b_ratio[i][j] = makeLognormalConstraint( NP_name, ttwj_3b2b_ratio[i][j], ttwj_3b2b_ratio_err[i][j] ) ;
	      } else {
		rar_ttwj_3b2b_ratio[i][j] = makeGaussianConstraint( NP_name, ttwj_3b2b_ratio[i][j], ttwj_3b2b_ratio_err[i][j] ) ;
	      }
	    } else {
	      rar_ttwj_3b2b_ratio[i][j] = makeGaussianConstraint( NP_name, ttwj_3b2b_ratio[i][j], 0. ) ; //-- this just creates a constant.
	      workspace.import( *(rar_ttwj_3b2b_ratio[i][j]) ) ;
	    }

	    if ( nBinsBtag > 3 ) {

	      sprintf( NP_name, "ttwj_4b2b_ratio_M%d_H%d", i+1, j+1 ) ;
	      if ( ! (blindStudy && blind0lepBin[i][j]) ) {
		if ( useLognormal ) {
		  rar_ttwj_4b2b_ratio[i][j] = makeLognormalConstraint( NP_name, ttwj_4b2b_ratio[i][j], ttwj_4b2b_ratio_err[i][j] ) ;
                } else {
		  rar_ttwj_4b2b_ratio[i][j] = makeGaussianConstraint( NP_name, ttwj_4b2b_ratio[i][j], ttwj_4b2b_ratio_err[i][j] ) ;
		}
	      } else {
		rar_ttwj_4b2b_ratio[i][j] = makeGaussianConstraint( NP_name, ttwj_4b2b_ratio[i][j], 0. ) ; //-- this just creates a constant.
		workspace.import( *(rar_ttwj_4b2b_ratio[i][j]) ) ;
	      }
	    }
	    
	    //--- ttwj 1lepSig
	    
	    sprintf( NP_name, "ttwjSlSig_3b2b_ratio_M%d_H%d", i+1, j+1 ) ;
	    if ( useLognormal ) {
	      rar_ttwjSlSig_3b2b_ratio[i][j] = makeLognormalConstraint( NP_name, ttwjSlSig_3b2b_ratio[i][j], ttwjSlSig_3b2b_ratio_err[i][j] ) ;
	    } else {
	      rar_ttwjSlSig_3b2b_ratio[i][j] = makeGaussianConstraint( NP_name, ttwjSlSig_3b2b_ratio[i][j], ttwjSlSig_3b2b_ratio_err[i][j] ) ;
	    }

	    if ( nBinsBtag > 3 ) {

	      sprintf( NP_name, "ttwjSlSig_4b2b_ratio_M%d_H%d", i+1, j+1 ) ;
	      if ( useLognormal ) {
		rar_ttwjSlSig_4b2b_ratio[i][j] = makeLognormalConstraint( NP_name, ttwjSlSig_4b2b_ratio[i][j], ttwjSlSig_4b2b_ratio_err[i][j] ) ;
	      } else {
		rar_ttwjSlSig_4b2b_ratio[i][j] = makeGaussianConstraint( NP_name, ttwjSlSig_4b2b_ratio[i][j], ttwjSlSig_4b2b_ratio_err[i][j] ) ;
	      }
	    }

	    printf("++++++++++++++++\n") ;
	  } // j (ht)
	  printf("=========================\n") ;
	} // k (nbtag)

      }



      printf("\n\n QCD nB ratio factors.\n\n") ; cout << flush ;

      // QCD nB ratio 

      RooAbsReal* rar_qcd_2b1b_ratio;      
      RooAbsReal* rar_qcd_3b1b_ratio;      

      if ( qcdModelIndex == 2 ) {
	if ( useLognormal ) {
	  rar_qcd_2b1b_ratio = makeLognormalConstraint( "qcd_2b1b_ratio", qcd_2b1b_ratio, qcd_2b1b_ratio_err ) ;
	  rar_qcd_3b1b_ratio = makeLognormalConstraint( "qcd_231b_ratio", qcd_3b1b_ratio, qcd_3b1b_ratio_err ) ;
	} else {
	  rar_qcd_2b1b_ratio = makeGaussianConstraint( "qcd_2b1b_ratio", qcd_2b1b_ratio, qcd_2b1b_ratio_err ) ;
	  rar_qcd_3b1b_ratio = makeGaussianConstraint( "qcd_3b1b_ratio", qcd_3b1b_ratio, qcd_3b1b_ratio_err ) ;
	}
      }
      

      printf("\n\n Znn nB ratio factors.\n\n") ; cout << flush ;

      // Znn nB ratio 

      RooAbsReal* rar_znn_2b1b_ratio;      
      RooAbsReal* rar_znn_3b1b_ratio;      

      if ( znnModelIndex == 2 ) {
	if ( useLognormal ) {
	  rar_znn_2b1b_ratio = makeLognormalConstraint( "znn_2b1b_ratio", Znn_2b1b_ratio, Znn_2b1b_ratio_err ) ;
	  rar_znn_3b1b_ratio = makeLognormalConstraint( "znn_231b_ratio", Znn_3b1b_ratio, Znn_3b1b_ratio_err ) ;
	} else {
	  rar_znn_2b1b_ratio = makeGaussianConstraint( "znn_2b1b_ratio", Znn_2b1b_ratio, Znn_2b1b_ratio_err ) ;
	  rar_znn_3b1b_ratio = makeGaussianConstraint( "znn_3b1b_ratio", Znn_3b1b_ratio, Znn_3b1b_ratio_err ) ;
	}
      }
      


      RooAbsReal* rar_eff_sf[nBinsMET][nBinsHT][nBinsBtag] ;
      RooAbsReal* rar_eff_sf_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
      RooAbsReal* rar_eff_sf_sl[nBinsMET][nBinsHT][nBinsBtag] ;

      char effbpname[1000] ;
      sprintf( effbpname, "eff_sf" ) ;

      for (int i = 0 ; i < nBinsMET ; i++) {
         for (int j = 0 ; j < nBinsHT ; j++) {
            if ( ignoreBin[i][j] ) continue ;
            for (int k = 0 ; k < nBinsBtag ; k++) {

               sprintf( NP_name, "eff_sf_M%d_H%d_%db", i+1, j+1, k+1 ) ;
              if ( !blindStudy ) {
                  if ( useLognormal ) {
                     rar_eff_sf[i][j][k] = makeLognormalConstraint( NP_name, rv_mean_eff_sf[i][j][k]->getVal(), rv_width_eff_sf[i][j][k]->getVal() ) ;
                  } else {
                     rar_eff_sf[i][j][k] = makeGaussianConstraint( NP_name, rv_mean_eff_sf[i][j][k]->getVal(), rv_width_eff_sf[i][j][k]->getVal() ) ;
                  }
               } else {
                  rar_eff_sf[i][j][k] = makeGaussianConstraint( NP_name, rv_mean_eff_sf[i][j][k]->getVal(), 0.0 ) ; // this just makes a constant.
               }

               sprintf( NP_name, "eff_sf_slSig_M%d_H%d_%db", i+1, j+1, k+1 ) ;
              if ( !blindStudy ) {
                  if ( useLognormal ) {
                     rar_eff_sf_slSig[i][j][k] = makeLognormalConstraint( NP_name, rv_mean_eff_sf_slSig[i][j][k]->getVal(), rv_width_eff_sf_slSig[i][j][k]->getVal() ) ;
                  } else {
                     rar_eff_sf_slSig[i][j][k] = makeGaussianConstraint( NP_name, rv_mean_eff_sf_slSig[i][j][k]->getVal(), rv_width_eff_sf_slSig[i][j][k]->getVal() ) ;
                  }
               } else {
                  rar_eff_sf_slSig[i][j][k] = makeGaussianConstraint( NP_name, rv_mean_eff_sf_slSig[i][j][k]->getVal(), 0.0 ) ; // this just makes a constant.
               }


               sprintf( NP_name, "eff_sf_sl_M%d_H%d_%db", i+1, j+1, k+1 ) ;
               if ( !blindStudy ) {
                  if ( useLognormal ) {
                     rar_eff_sf_sl[i][j][k] = makeLognormalConstraint( NP_name, rv_mean_eff_sf_sl[i][j][k]->getVal(), rv_width_eff_sf_sl[i][j][k]->getVal() ) ;
                  } else {
                     rar_eff_sf_sl[i][j][k] = makeGaussianConstraint( NP_name, rv_mean_eff_sf_sl[i][j][k]->getVal(), rv_width_eff_sf_sl[i][j][k]->getVal() ) ;
                  }
               } else {
                  rar_eff_sf_sl[i][j][k] = makeGaussianConstraint( NP_name, rv_mean_eff_sf_sl[i][j][k]->getVal(), 0.0 ) ; // this just makes a constant.
               }

            } // k
         } // j
      } // i




      int nShapeSystematics(0) ;
      char shapeSystName[20][100] ;

      bool sss_return_status ;

      sprintf( shapeSystName[nShapeSystematics], "btageff_sf" ) ;
      if ( useLognormal ) {
         sss_return_status = setupShapeSyst( inputSusy_deff_dbtageff_file, "btageff_sf", 2, m0, m12, workspace ) ; // 2 = log-normal
      } else {
         sss_return_status = setupShapeSyst( inputSusy_deff_dbtageff_file, "btageff_sf", 1, m0, m12, workspace ) ; // 1 = Gaussian
      }
      if ( !sss_return_status ) { return false ; }
      nShapeSystematics++ ;

      sprintf( shapeSystName[nShapeSystematics], "btageff_lf_sf" ) ;
      if ( useLognormal ) {
         sss_return_status = setupShapeSyst( inputSusy_deff_dbtageff_lightflavor_file, "btageff_lf_sf", 2, m0, m12, workspace ) ; // 2 = log-normal
      } else {
         sss_return_status = setupShapeSyst( inputSusy_deff_dbtageff_lightflavor_file, "btageff_lf_sf", 1, m0, m12, workspace ) ; // 1 = Gaussian
      }
      if ( !sss_return_status ) { return false ; }
      nShapeSystematics++ ;

      sprintf( shapeSystName[nShapeSystematics], "JES_sf" ) ;
      if ( useLognormal ) {
         sss_return_status = setupShapeSyst( systFile1                   , "JES_sf"    , 2, m0, m12, workspace ) ; // 2 = log-normal
      } else {
         sss_return_status = setupShapeSyst( systFile1                   , "JES_sf"    , 1, m0, m12, workspace ) ; // 1 = Gaussian
      }
      if ( !sss_return_status ) { return false ; }
      nShapeSystematics++ ;

      sprintf( shapeSystName[nShapeSystematics], "pdfsyst_sf" ) ;
      if ( useLognormal ) {
         sss_return_status = setupShapeSyst( pdf_syst_file               , "pdfsyst_sf"    , 2, m0, m12, workspace ) ; // 2 = log-normal
      } else {
         sss_return_status = setupShapeSyst( pdf_syst_file               , "pdfsyst_sf"    , 1, m0, m12, workspace ) ; // 1 = Gaussian
      }
      if ( !sss_return_status ) { return false ; }
      nShapeSystematics++ ;

      sprintf( shapeSystName[nShapeSystematics], "isrsyst_sf" ) ;
      if ( useLognormal ) {
         sss_return_status = setupShapeSyst( isr_syst_file               , "isrsyst_sf"    , 2, m0, m12, workspace ) ; // 2 = log-normal
      } else {
         sss_return_status = setupShapeSyst( isr_syst_file               , "isrsyst_sf"    , 1, m0, m12, workspace ) ; // 1 = Gaussian
      }
      if ( !sss_return_status ) { return false ; }
      nShapeSystematics++ ;


      //--- Jan28, 2013: New W+jets and single top Xsec shape systematics (for ttwj closure, not SUSY signal efficiency). ------
      if ( useLognormal ) {
         sss_return_status = setupShapeSyst( wjets_xsec_shapesyst_file, "wjets_xsec", 2, -1, -1, workspace ) ; // 2 = log-normal
      } else {
         sss_return_status = setupShapeSyst( wjets_xsec_shapesyst_file, "wjets_xsec", 1, -1, -1, workspace ) ; // 1 = Gaussian
      }
      if ( !sss_return_status ) { return false ; }

      if ( useLognormal ) {
         sss_return_status = setupShapeSyst( singletop_xsec_shapesyst_file, "singletop_xsec", 2, -1, -1, workspace ) ; // 2 = log-normal
      } else {
         sss_return_status = setupShapeSyst( singletop_xsec_shapesyst_file, "singletop_xsec", 1, -1, -1, workspace ) ; // 1 = Gaussian
      }
      if ( !sss_return_status ) { return false ; }
      //-----------------------------------------------------------------------------------------------------------




      // make beta constraints for trigger efficiency corrections
      printf("\n\n trigger efficiencies.\n\n") ; cout << flush ;

      RooAbsReal* rar_trigeff[nBinsMET][nBinsHT];
      RooAbsReal* rar_trigeff_sl[nBinsMET][nBinsHT];   // this applies to both sl and slSig

      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {

          if ( ignoreBin[i][j] ) continue ;

          sprintf( NP_name, "trigeff_M%d_H%d", i+1, j+1 ) ;
          if ( useBeta ) {
	    rar_trigeff[i][j] = makeBetaConstraint( NP_name, trigeff_0L[i][j], trigefferr_0L[i][j], workspace ) ;
          } else {
             rar_trigeff[i][j] = makeGaussianConstraint( NP_name, trigeff_0L[i][j], trigefferr_0L[i][j] ) ;
          }

          sprintf( NP_name, "trigeff_sl_M%d_H%d", i+1, j+1 ) ;
          if ( useBeta ) {
	    rar_trigeff_sl[i][j] = makeBetaConstraint( NP_name, trigeff_1L[i][j], trigefferr_1L[i][j], workspace ) ;
          } else {
             rar_trigeff_sl[i][j] = makeGaussianConstraint( NP_name, trigeff_1L[i][j], trigefferr_1L[i][j] ) ;
          }
        }
      }










      //--- Oct 31, 2012: added global uncertainties.

      printf("\n\n === Global uncertainties:\n") ;
      double all_gu2(0.) ;
      for ( int gui=0; gui < n_global_uncertainties; gui++ ) {
         printf("  %20s : %5.3f\n", global_uncertainty_name[gui], global_uncertainty_val[gui] ) ;
         all_gu2 += pow( global_uncertainty_val[gui], 2 ) ;
      } // gui.
      double all_gu = sqrt( all_gu2 ) ;
      printf("  %20s : %5.3f\n\n\n", "Total", all_gu ) ;


      RooAbsReal* rar_all_gu(0x0) ;

      if ( useLognormal ) {
         rar_all_gu = makeLognormalConstraint( "all_gu", 1.0, all_gu ) ;
      } else {
         rar_all_gu = makeGaussianConstraint( "all_gu", 1.0, all_gu ) ;
      }



      //+++++++++++++++++ Relationships between parameters ++++++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining relationships between parameters.\n" ) ;
      cout << flush ;

      RooArgSet pdflist ;


      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {

          if ( ignoreBin[i][j] ) continue ;

	  for (int k = 0 ; k < nBinsBtag ; k++) {     

	    printf(" met,ht,nb : %d, %d, %d : ", i+1, j+1, k+1 ) ; cout << flush ;
	    
	    printf(" ttwj,") ; cout << flush ;
	    
	    //---- TTWJ


	    // 1lep
	    
	    TString ttwjSlString   = "mu_ttwj_sl";
	    ttwjSlString   += sMbins[i]+sHbins[j]+sBbins[k] ;
	    
	    if ( k < 2 || !constrainBjetShape ) {
	      TString ttwjSlrfvString = "@0" ;
	      rfv_mu_ttwj_sl[i][j][k] = new RooFormulaVar( ttwjSlString, ttwjSlrfvString, RooArgSet( *rv_mu_ttwj_sl[i][j][k] ) ) ;
	    }
	    else if ( k == 2 && constrainBjetShape ) {
	      TString ttwjSlrfvString = " @0 * @1" ;
	      rfv_mu_ttwj_sl[i][j][k] = new RooFormulaVar( ttwjSlString, ttwjSlrfvString,
							   RooArgSet( *rv_mu_ttwj_sl[i][j][1], *rar_ttwj_3b2b_ratio[i][j] ) ) ;
	    }
	    else if ( k == 3 && constrainBjetShape ) {
	      TString ttwjSlrfvString = " @0 * @1" ;
	      rfv_mu_ttwj_sl[i][j][k] = new RooFormulaVar( ttwjSlString, ttwjSlrfvString,
							   RooArgSet( *rv_mu_ttwj_sl[i][j][1], *rar_ttwj_4b2b_ratio[i][j] ) ) ;
	    }
	    else {
	      cout << "\nconstrainBjetShape with nB != 3,4 ??? I have no clue what I'm doing... exiting ... \n" << endl ; return false ;
	    }

	    rv_mu_ttwj_sl[i][j][k] = rfv_mu_ttwj_sl[i][j][k] ;


	    // 1lepSig
	    
	    TString ttwjSlSigString   = "mu_ttwj_slSig";
	    ttwjSlSigString   += sMbins[i]+sHbins[j]+sBbins[k] ;

	    if ( floatSLSigRatios ) { 
	      TString ttwjSlSigrfvString = " @0 * @1 * @2 * @3" ;
	      rfv_mu_ttwj_slSig[i][j][k] = new RooFormulaVar( ttwjSlSigString, ttwjSlSigrfvString,
							      RooArgSet( *rv_mu_ttwj_sl[i][j][k], *rv_ttwj_slSigsl_ratio, 
									 *rar_sf_ttwj_slSig[i][j][2], *rv_ttwj_slSigDD_ratio[i][j] ) ) ;
	    }
	    else {

	      if ( k < 2 || !constrainBjetShape ) {
		TString ttwjSlSigrfvString = " @0 * @1 * @2" ;
		rfv_mu_ttwj_slSig[i][j][k] = new RooFormulaVar( ttwjSlSigString, ttwjSlSigrfvString,
								RooArgSet( *rv_mu_ttwj_sl[i][j][k], *rv_ttwj_slSigsl_ratio, *rar_sf_ttwj_slSig[i][j][k] ) ) ;
	      }
	      else if ( k == 2 && constrainBjetShape ) {
		TString ttwjSlSigrfvString = " @0 * @1" ;
		rfv_mu_ttwj_slSig[i][j][k] = new RooFormulaVar( ttwjSlSigString, ttwjSlSigrfvString,
								RooArgSet( *rv_mu_ttwj_slSig[i][j][1], *rar_ttwjSlSig_3b2b_ratio[i][j] ) ) ;
	      }
	      else if ( k == 3 && constrainBjetShape ) {
		TString ttwjSlSigrfvString = " @0 * @1" ;
		rfv_mu_ttwj_slSig[i][j][k] = new RooFormulaVar( ttwjSlSigString, ttwjSlSigrfvString,
								RooArgSet( *rv_mu_ttwj_slSig[i][j][1], *rar_ttwjSlSig_4b2b_ratio[i][j] ) ) ;
	      }
	      else {
		cout << "\nconstrainBjetShape with nB != 3,4 ??? I have no clue what I'm doing... exiting ... \n" << endl ; return false ;
	      }

	    }

            rv_mu_ttwj_slSig[i][j][k] = rfv_mu_ttwj_slSig[i][j][k] ;


	    // 0lep

	    TString ttwjString   = "mu_ttwj";
	    ttwjString   += sMbins[i]+sHbins[j]+sBbins[k] ;

            char xsecshapesystname[1000] ;
	    
            sprintf( xsecshapesystname, "wjets_xsec_M%d_H%d_%db", i+1, j+1, k+1 ) ;
            RooAbsReal* rar_wjets_xsec_shape_syst = (RooAbsReal*) workspace.obj( xsecshapesystname ) ;
            if ( rar_wjets_xsec_shape_syst == 0x0 ) { printf("\n\n *** Missing NP : %s\n\n", xsecshapesystname ) ; return false ; }
	    
            sprintf( xsecshapesystname, "singletop_xsec_M%d_H%d_%db", i+1, j+1, k+1 ) ;
            RooAbsReal* rar_singletop_xsec_shape_syst = (RooAbsReal*) workspace.obj( xsecshapesystname ) ;
            if ( rar_singletop_xsec_shape_syst == 0x0 ) { printf("\n\n *** Missing NP : %s\n\n", xsecshapesystname ) ; return false ; }

            TString ttwjrfvString = " @0 * @1 * @2 * @3 * @4" ;
            rfv_mu_ttwj[i][j][k] = new RooFormulaVar( ttwjString, ttwjrfvString,
						      RooArgSet( *rar_wjets_xsec_shape_syst, *rar_singletop_xsec_shape_syst,
								 *rar_sf_ttwj[i][j][k], *rv_mu_ttwj_sl[i][j][k], *rv_ttwj_0lep1lep_ratio ) ) ;

            rv_mu_ttwj[i][j][k] = rfv_mu_ttwj[i][j][k] ;


	  
	    printf(" QCD,") ; cout << flush ;

	    //---- QCD

            TString qcdString    = "mu_qcd" ;
            qcdString    += sMbins[i]+sHbins[j]+sBbins[k] ;

	    if ( qcdModelIndex == 1 ) {
	      // do nothing, keep the QCD shape from MC
	    }

	    if ( qcdModelIndex == 2 ) {

	      // take the shape from 1B bins, then float 2B and 3B using loose constraints
	      
	      if ( k == 0 ) {

		TString rfvQcdString = "@0 * @1" ;

		rfv_mu_qcd[i][j][k] = new RooFormulaVar( qcdString, rfvQcdString, 
							 RooArgSet( *rv_mu_qcd[i][j][k], *rv_qcd_shape_ratio[i][j] ) ) ;
		rv_mu_qcd[i][j][k] = rfv_mu_qcd[i][j][k] ;

	      }
	      else if ( k == 1 ) {

               TString rfvQcdString = "@0 * @1" ;

               rfv_mu_qcd[i][j][k] = new RooFormulaVar( qcdString, rfvQcdString, 
                                                        RooArgSet( *rv_mu_qcd[i][j][0], *rar_qcd_2b1b_ratio ) ) ;
	       rv_mu_qcd[i][j][k] = rfv_mu_qcd[i][j][k] ;

	      }
	      else if ( k == 2 ) {

               TString rfvQcdString = "@0 * @1" ;

               rfv_mu_qcd[i][j][k] = new RooFormulaVar( qcdString, rfvQcdString, 
                                                        RooArgSet( *rv_mu_qcd[i][j][0], *rar_qcd_3b1b_ratio ) ) ;
	       rv_mu_qcd[i][j][k] = rfv_mu_qcd[i][j][k] ;

	      }

	    }



	  
	    printf(" ZNN,") ; cout << flush ;

	    //---- ZNN

            TString znnString    = "mu_znn" ;
            znnString    += sMbins[i]+sHbins[j]+sBbins[k] ;

	    if ( znnModelIndex == 1 ) {
	      // do nothing, keep the ZNN shape from MC
	    }

	    if ( znnModelIndex == 2 ) {

	      // take the shape from 1B bins, then float 2B and 3B using loose constraints
	      
	      if ( k == 0 ) {
		// nothing
	      }
	      else if ( k == 1 ) {

               TString rfvZnnString = "@0 * @1" ;

               rfv_mu_znn[i][j][k] = new RooFormulaVar( znnString, rfvZnnString, 
                                                        RooArgSet( *rv_mu_znn[i][j][0], *rar_znn_2b1b_ratio ) ) ;
	       rv_mu_znn[i][j][k] = rfv_mu_znn[i][j][k] ;

	      }
	      else if ( k == 2 ) {

               TString rfvZnnString = "@0 * @1" ;

               rfv_mu_znn[i][j][k] = new RooFormulaVar( znnString, rfvZnnString, 
                                                        RooArgSet( *rv_mu_znn[i][j][0], *rar_znn_3b1b_ratio ) ) ;
	       rv_mu_znn[i][j][k] = rfv_mu_znn[i][j][k] ;

	      }

	    }


	    printf(" susy,") ; cout << flush ;

	    //---- SUSY

	    TString susyString      = "mu_susy" ;
	    TString susySlSigString = "mu_susy_slSig" ;
	    TString susySlString    = "mu_susy_sl" ;
	    
	    susyString      += sMbins[i]+sHbins[j]+sBbins[k] ;
	    susySlSigString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    susySlString    += sMbins[i]+sHbins[j]+sBbins[k] ;



            TString rfvSusyString = "@0 * ( @1 / @2)" ;

            rv_mu_susy[i][j][k] = new RooFormulaVar( susyString, rfvSusyString,
                                                       RooArgSet( *rv_mu_susymc[i][j][k], *rv_mu_susy_all0lep, *rv_mu_susymc_all0lep ) ) ;

            TString rfvSusySlSigString  = "@0 * ( @1 / @2)" ;

            rv_mu_susy_slSig[i][j][k] = new RooFormulaVar( susySlSigString, rfvSusySlSigString,
							   RooArgSet( *rv_mu_susymc_slSig[i][j][k], *rv_mu_susy_all0lep, *rv_mu_susymc_all0lep ) ) ;

            TString rfvSusySlString  = "@0 * ( @1 / @2)" ;

            rv_mu_susy_sl[i][j][k] = new RooFormulaVar( susySlString, rfvSusySlString,
                                                        RooArgSet( *rv_mu_susymc_sl[i][j][k], *rv_mu_susy_all0lep, *rv_mu_susymc_all0lep ) ) ;


	    

            //+++++++++++++ Expected counts for observables in terms of parameters ++++++++++++++++++

	    TString nString      = "n" ;
	    TString nSlSigString = "n_slSig" ;
	    TString nSlString    = "n_sl" ;

	    nString      += sMbins[i]+sHbins[j]+sBbins[k] ;
	    nSlSigString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    nSlString    += sMbins[i]+sHbins[j]+sBbins[k] ;


            char systparname[100] ;
            char shapesystprodname[1000] ;


            char systprodeqn[1000] ;
            sprintf( systprodeqn, "@0" ) ;
            for ( int si=1; si<nShapeSystematics; si++ ) {
               char tmpstr[1000] ;
               sprintf( tmpstr, "%s * @%d", systprodeqn, si ) ;
               sprintf( systprodeqn, "%s", tmpstr ) ;
            } // si




          //--- Zero lepton : n

            RooArgSet shapeSystProdSet_zl ;
            for ( int si=0; si<nShapeSystematics; si++ ) {
	      sprintf( systparname, "%s_M%d_H%d_%db", shapeSystName[si], i+1, j+1, k+1 ) ;
	      RooAbsReal* rar_sf = (RooAbsReal*) workspace.obj( systparname ) ;
	      if ( rar_sf == 0x0 ) { printf("\n\n *** initialize: missing nuisance parameter: %s\n\n", systparname ) ; return false ; }
	      shapeSystProdSet_zl.add( *rar_sf ) ;
            } // si
            sprintf( shapesystprodname, "shapesyst_prod_M%d_H%d_%db", i+1, j+1, k+1 ) ;
            RooFormulaVar* rfv_shapeSystProd_zl = new RooFormulaVar( shapesystprodname, systprodeqn, shapeSystProdSet_zl ) ;

	    
            TString rfvNString =  "@0 * @1 + (@2 + @3 + (@4 * @5 * @6 * @7)) * @8" ;

            rv_n[i][j][k] = new RooFormulaVar( nString, rfvNString,
                                               RooArgSet( *rv_mu_qcd[i][j][k], *rar_trigeff[i][j],
                                                          *rv_mu_ttwj[i][j][k], *rv_mu_znn[i][j][k],
                                                          *rar_all_gu, *rfv_shapeSystProd_zl, *rar_eff_sf[i][j][k], *rv_mu_susy[i][j][k],
                                                          *rar_trigeff_sl[i][j] ) ) ;
							  


          //--- Sig Single lepton : n_slSig

            RooArgSet shapeSystProdSet_slSig ;
            for ( int si=0; si<nShapeSystematics; si++ ) {
	      sprintf( systparname, "%s_slSig_M%d_H%d_%db", shapeSystName[si], i+1, j+1, k+1 ) ;
	      RooAbsReal* rar_sf = (RooAbsReal*) workspace.obj( systparname ) ;
	      if ( rar_sf == 0x0 ) { printf("\n\n *** initialize: missing nuisance parameter: %s\n\n", systparname ) ; return false ; }
	      shapeSystProdSet_slSig.add( *rar_sf ) ;
            } // si
            sprintf( shapesystprodname, "shapesyst_prod_slSig_M%d_H%d_%db", i+1, j+1, k+1 ) ;
            RooFormulaVar* rfv_shapeSystProd_slSig = new RooFormulaVar( shapesystprodname, systprodeqn, shapeSystProdSet_slSig ) ;
	    
            TString rfvNSlSigString = "(@0 + (@1 * @2 * @3 * @4)) * @5" ;
	    
            rv_n_slSig[i][j][k] = new RooFormulaVar( nSlSigString, rfvNSlSigString,
						     RooArgSet( *rv_mu_ttwj_slSig[i][j][k],
								*rar_all_gu, *rfv_shapeSystProd_slSig, *rar_eff_sf_slSig[i][j][k], *rv_mu_susy_slSig[i][j][k],
								*rar_trigeff_sl[i][j] ) ) ;
	    

          //--- Single lepton : n_sl

            RooArgSet shapeSystProdSet_sl ;
            for ( int si=0; si<nShapeSystematics; si++ ) {
               sprintf( systparname, "%s_sl_M%d_H%d_%db", shapeSystName[si], i+1, j+1, k+1 ) ;
               RooAbsReal* rar_sf = (RooAbsReal*) workspace.obj( systparname ) ;
               if ( rar_sf == 0x0 ) { printf("\n\n *** initialize: missing nuisance parameter: %s\n\n", systparname ) ; return false ; }
               shapeSystProdSet_sl.add( *rar_sf ) ;
            } // si
            sprintf( shapesystprodname, "shapesyst_prod_sl_M%d_H%d_%db", i+1, j+1, k+1 ) ;
            RooFormulaVar* rfv_shapeSystProd_sl = new RooFormulaVar( shapesystprodname, systprodeqn, shapeSystProdSet_sl ) ;

            TString rfvNSlString = "(@0 + (@1 * @2 * @3 * @4)) * @5" ;

            rv_n_sl[i][j][k] = new RooFormulaVar( nSlString, rfvNSlString,
                                                  RooArgSet( *rv_mu_ttwj_sl[i][j][k], 
                                                             *rar_all_gu, *rfv_shapeSystProd_sl, *rar_eff_sf_sl[i][j][k], *rv_mu_susy_sl[i][j][k],
                                                             *rar_trigeff_sl[i][j] ) ) ;




	    //////--------------------------------


	    // pdf's

	    TString pdfN0lepString    = "pdf_N_0lep";
	    TString pdfN1lepSigString = "pdf_N_1lepSig";
	    TString pdfN1lepString    = "pdf_N_1lep";

	    pdfN0lepString    += sMbins[i]+sHbins[j]+sBbins[k] ;
	    pdfN1lepSigString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    pdfN1lepString    += sMbins[i]+sHbins[j]+sBbins[k] ;

            if ( blind0lepBin[i][j][k] ) {
               printf(" * Not including 0lep met=%d, ht=%d, nb=%d PDF in the likelihood.\n", i+1, j+1, k+1 ) ;
               pdf_N_0lep[i][j][k] = new RooConstVar( pdfN0lepString, pdfN0lepString, 1.0 ) ;
               workspace.import( *rv_n[i][j][k] ) ; //--- save the equation for the number of events.
            } else {
	       pdf_N_0lep[i][j][k] = new RooPoisson( pdfN0lepString, pdfN0lepString, *rv_0lep[i][j][k], *rv_n[i][j][k] ) ;
	       pdflist.add( *pdf_N_0lep[i][j][k] ) ;
            }

	    pdf_N_1lepSig[i][j][k] = new RooPoisson( pdfN1lepSigString, pdfN1lepSigString, *rv_1lepSig[i][j][k], *rv_n_slSig[i][j][k] ) ;
	    pdf_N_1lep[i][j][k]    = new RooPoisson( pdfN1lepString,    pdfN1lepString,    *rv_1lep[i][j][k],    *rv_n_sl[i][j][k] ) ;

	    //pdflist.add( *pdf_N_1lepSig[i][j][k] ) ;
	    pdflist.add( *pdf_N_1lep[i][j][k] ) ;

	  }

	}
      } 
     

      printf(" --- Constructing likelihood.\n" ) ; cout << flush ;

      pdflist.add( *allNuisancePdfs );
	    
      pdflist.Print() ;

      likelihood = new RooProdPdfLogSum( "likelihood", "ra2b likelihood", pdflist ) ;


      /// likelihood->printMultiline( cout, 1, kTRUE, "" ) ;


      //--- Do a simple pre-fit to tune initial values of key parameters.
      printf(" --- Performing simple pre-fit of 0-lep ttwj normalization.\n") ; cout << flush ;


      { // begin scoping bracket.

         ((RooRealVar*)rv_mu_susy_all0lep) -> setVal(0.) ;

           {
            double logL(0.) ;
            for (int i = 0 ; i < nBinsMET ; i++) {
               for (int j = 0 ; j < nBinsHT ; j++) {
                  if ( ignoreBin[i][j] ) continue ;
                  for (int k = 0 ; k < nBinsBtag ; k++) {     
                     double pdfVal = pdf_N_0lep[i][j][k] -> getVal() ;
                     if ( pdfVal > 0. ) { logL += log( pdfVal ) ; } else { printf(" *** PDF %s evaluates to %g\n", pdf_N_0lep[i][j][k] -> GetName(), pdfVal ) ; }
                     printf("       PDF %20s, val=%8.6f,  N=%5.0f,  n=%6.1f\n", pdf_N_0lep[i][j][k] -> GetName(), pdfVal, rv_0lep[i][j][k]->getVal(), rv_n[i][j][k]->getVal() ) ;
                  }
               }
            }
            printf( "\n val = %6.3f, ln L = %g\n\n\n", rv_ttwj_0lep1lep_ratio->getVal(), logL ) ;
           }

         printf("\n\n Initial guess for ttwj 0lep/1lep ratio: %6.3f\n\n", rv_ttwj_0lep1lep_ratio->getVal() ) ;

         double scanLow  = 0.7 ;
         double scanHigh = 1.5 ;
         double bestVal = rv_ttwj_0lep1lep_ratio->getVal() ;
         double bestlnL( -1.e99 ) ;
         int nScanPoints(50) ;
         if ( scanLow < 0. ) { scanLow = 1. ; }
         for ( int spi=0 ; spi<nScanPoints; spi++ ) {
            double logL(0.) ;
            double scanVal = scanLow + (scanHigh-scanLow)/(nScanPoints-1.)*spi ;
            rv_ttwj_0lep1lep_ratio -> setVal( scanVal ) ;
            for (int i = 0 ; i < nBinsMET ; i++) {
               for (int j = 0 ; j < nBinsHT ; j++) {
                  if ( ignoreBin[i][j] ) continue ;
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
            printf( " val = %6.3f, ln L = %g\n", scanVal, logL ) ;
         } // spi.

         printf("\n\n  Best val = %6.3f, ln L = %g\n\n\n", bestVal, bestlnL ) ;

         rv_ttwj_0lep1lep_ratio -> setVal(bestVal) ;

      } // end scoping bracket.



      // fix global observables to their default value

      TIterator * goIter = globalObservables->createIterator(); 

      while ( RooRealVar * gObs_rrv = (RooRealVar*)goIter->Next() ) {
	gObs_rrv->setConstant();
      }


      // parameters of interest
      RooArgSet poi(*rv_mu_susy_all0lep, "poi");
      // flat prior for POI
      RooUniform signal_prior ("signal_prior","signal_prior",*rv_mu_susy_all0lep);

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

      // set all but obs, poi and nuisance to const
      SetConstants(workspace, sbModel);
      
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
      printf("\n\n Creating output root file: %s\n\n\n", wsrootfilename) ;
      workspace.writeToFile(wsrootfilename);


      return true;

    }  // end of initialize


   //=====================================================================================================================

    bool ra2bRoostatsClass3D_3b::setSusyScanPoint( const char* inputScanFile,
						   double m0, double m12
						   ) {


      printf("\n\n Opening SUSY scan input file : %s\n", inputScanFile ) ;


      ifstream infp ;
      infp.open(inputScanFile) ;
      if ( !infp.good() ) {
	printf("\n\n *** Problem opening input file: %s.\n\n", inputScanFile ) ;
	return false ;
      } 

     //--- check that file has expected format.

      int ArraySize = 3 + 2*4*(nBinsMET*nBinsHT*nBinsBtag) ;  // 2: signal counts + associated error
                                                              // 4: number of samples (0lep, 1lepSig, 1lep, ldp)
      char command[10000] ;
      sprintf(command, "head -1 %s | awk '{print NF}' | grep -q %d", inputScanFile, ArraySize ) ;
      int returnStat = gSystem->Exec(command ) ;
      if ( returnStat !=0 ) {
         printf("\n\n\n *** setSusyScanPoint : expecting %d fields per line in input file %s.  Found ", ArraySize, inputScanFile ) ; cout << flush ;
         sprintf( command, "head -1 %s | awk '{print NF}'", inputScanFile ) ;
         gSystem->Exec(command ) ; cout << flush ;
         printf("\n\n") ;
         return false ;
      }

      
      double deltaM0(0.) ;
      double deltaM12(0.) ;
      
      deltaM0 = 25 ;
      deltaM12 = 25 ;
      
      bool found(false) ;

     
      TString sMbins[nBinsMET];
      TString sHbins[nBinsHT];
      TString sBbins[4] = {"_1b","_2b","_3b","_4b"};
      
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


      float pointM0 ;
      float pointM12 ;
      
      float n_0l_raw[nBinsMET][nBinsHT][nBinsBtag] ;
      float n_1lSig_raw[nBinsMET][nBinsHT][nBinsBtag] ;
      float n_1l_raw[nBinsMET][nBinsHT][nBinsBtag] ;
      
      float n_0l_correction[nBinsMET][nBinsHT][nBinsBtag] ;
      float n_1lSig_correction[nBinsMET][nBinsHT][nBinsBtag] ;
      float n_1l_correction[nBinsMET][nBinsHT][nBinsBtag] ;
      
      float n_0l_error[nBinsMET][nBinsHT][nBinsBtag] ;
      float n_1lSig_error[nBinsMET][nBinsHT][nBinsBtag] ;
      float n_1l_error[nBinsMET][nBinsHT][nBinsBtag] ;
      
      int nGen ;


      //--- Loop over the scan points.
      while ( infp.good() ) {

	if ( found ) break ;

	double ArrayContent[ArraySize] ;

	for (int i = 0; infp && i < ArraySize; ++ i) {
	  infp >> ArrayContent[i];
	}

	pointM0  = ArrayContent[0] ;
	pointM12 = ArrayContent[1] ;

        printf("  mgl=%g, mlsp=%g.  looking for %g, %g.\n", pointM0, pointM12, m0, m12 ) ;

	nGen = (int)ArrayContent[2] ;

	if (    fabs( pointM0 - m0 ) <= deltaM0/2 && fabs( pointM12 - m12 ) <= deltaM12/2 ) {

	  found = true ;

	  for (int i = 0 ; i < nBinsMET ; i++) {
	    for (int j = 0 ; j < nBinsHT ; j++) {
	      for (int k = 0 ; k < nBinsBtag ; k++) {     

		n_0l_raw[i][j][k]    = ArrayContent[3 + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;
		n_1lSig_raw[i][j][k] = ArrayContent[3 + (nBinsMET*nBinsHT*nBinsBtag) + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;
		n_1l_raw[i][j][k]    = ArrayContent[3 + 2*(nBinsMET*nBinsHT*nBinsBtag) + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;
		
		n_0l_correction[i][j][k]     = 1. ;
		n_1lSig_correction[i][j][k]  = 1. ;
		n_1l_correction[i][j][k]     = 1. ;
		
		n_0l_error[i][j][k]     = ArrayContent[3 + 4*(nBinsMET*nBinsHT*nBinsBtag) + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;
		n_1lSig_error[i][j][k]  = ArrayContent[3 + 5*(nBinsMET*nBinsHT*nBinsBtag) + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;
		n_1l_error[i][j][k]     = ArrayContent[3 + 6*(nBinsMET*nBinsHT*nBinsBtag) + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;

	    
		if ( n_0l_raw[i][j][k] > 0.00001 ) {
		  n_0l_error[i][j][k] = 100 * n_0l_error[i][j][k] / n_0l_raw[i][j][k] ;
		}
		else n_0l_error[i][j][k] = 0.1 ;
	    	    
		if ( n_1lSig_raw[i][j][k] > 0.00001 ) {
		  n_1lSig_error[i][j][k] = 100 * n_1lSig_error[i][j][k] / n_1lSig_raw[i][j][k] ;
		}
		else n_1lSig_error[i][j][k] = 0.1 ;
	    	    
		if ( n_1l_raw[i][j][k] > 0.00001 ) {
		  n_1l_error[i][j][k] = 100 * n_1l_error[i][j][k] / n_1l_raw[i][j][k] ;
		}
		else n_1l_error[i][j][k] = 0.1 ;

	      }
	    }
	  }

	  printf("\n\n Found event counts for point m0 = %4.0f,  m1/2 = %4.0f.\n\n", pointM0, pointM12 ) ;

	} // end check on mass point

      }	

      infp.close() ;

      if ( !found ) {
	printf("\n\n *** Point not found in scan.  Check this file %s\n  to see if this point is there: mgl=%4.0f, mlsp=%4.0f\n\n", inputScanFile, m0, m12 ) ;
	return false ;
      }
 
      double setVal_n_0l[nBinsMET][nBinsHT][nBinsBtag] ;
      double setVal_n_1lSig[nBinsMET][nBinsHT][nBinsBtag] ;
      double setVal_n_1l[nBinsMET][nBinsHT][nBinsBtag] ;
      
      double all0lep(0.) ;
      
      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  if ( ignoreBin[i][j] ) continue ;
	  for (int k = 0 ; k < nBinsBtag ; k++) {     
	    
	    setVal_n_0l[i][j][k]    = n_0l_raw[i][j][k]    * n_0l_correction[i][j][k] ;
	    setVal_n_1lSig[i][j][k] = n_1lSig_raw[i][j][k] * n_1lSig_correction[i][j][k] ;
	    setVal_n_1l[i][j][k]    = n_1l_raw[i][j][k]    * n_1l_correction[i][j][k] ;
	    
	    all0lep += setVal_n_0l[i][j][k] ;

	    rv_mu_susymc[i][j][k]      -> setVal( setVal_n_0l[i][j][k] ) ;
	    rv_mu_susymc_slSig[i][j][k]-> setVal( setVal_n_1lSig[i][j][k] ) ;
	    rv_mu_susymc_sl[i][j][k]   -> setVal( setVal_n_1l[i][j][k] ) ;

	    rv_width_eff_sf[i][j][k]      -> setVal( n_0l_error[i][j][k] / 100. ) ;
	    rv_width_eff_sf_slSig[i][j][k]-> setVal( n_1lSig_error[i][j][k] / 100. ) ;
	    rv_width_eff_sf_sl[i][j][k]   -> setVal( n_1l_error[i][j][k] / 100. ) ;
	    
	  }
	}
      }


      // print out signal values
      
      cout << "\n\nSetting susy signal: \n" << endl ;
      
      double total0lep(0.) ;
      
      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  if ( ignoreBin[i][j] ) continue ;
	  for (int k = 0 ; k < nBinsBtag ; k++) {     
	    
	    TString binString = "";
	    binString += sMbins[i]+sHbins[j]+sBbins[k] ;	  
	    
	    cout << binString + " - 0 lep    - setting susy signal to " << setVal_n_0l[i][j][k]    << " +/- " << n_0l_error[i][j][k] << " %" << endl ;
	    cout << binString + " - 1 lepSig - setting susy signal to " << setVal_n_1lSig[i][j][k] << " +/- " << n_1lSig_error[i][j][k] << " %" << endl ;
	    cout << binString + " - 1 lep    - setting susy signal to " << setVal_n_1l[i][j][k]    << " +/- " << n_1l_error[i][j][k] << " %" << endl ;
	    
	    total0lep += setVal_n_0l[i][j][k] ;
	    
	  }
	}
      }


      printf("\n\n SUSY 0lep total: %7.2f\n\n", total0lep ) ;
 


      return true ;


    }  // end of setSusyScanPoint

   //==============================================================================================================

    void ra2bRoostatsClass3D_3b::mismatchErr( char* label, TString inPar ) {

      cout << "Mismatch in input file:" << endl;
      cout << "Expecting: " << inPar << endl;
      cout << "Reading:   " << label << endl;
      cout << "\nCheck binning!\n" << endl;

      return;

    }

   //==============================================================================================================


    RooAbsReal* ra2bRoostatsClass3D_3b::makeBetaPrimeConstraint( const char* NP_name, double NP_val, double NP_err ) {
      
       if ( NP_err <= 0. ) {
          printf("  Uncertainty is zero.  Will return constant scale factor of %g.  Input val = %g, err = %g.\n", NP_val, NP_val, NP_err ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }


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
       ///// rrv_failObs = new RooRealVar( varname, varname, beta -1., 1e-5, 1e5 ) ;
       rrv_failObs = new RooRealVar( varname, varname, beta +1., 1e-5, 1e5 ) ;
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

       ///// parVal = beta-1. ;
       parVal = beta+1. ;
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



   RooAbsReal* ra2bRoostatsClass3D_3b::makeBetaConstraint( const char* NP_name, double NP_val, double NP_err, RooWorkspace& workspace ) {

       if ( NP_err <= 0. ) {
          printf("  Uncertainty is zero.  Will return constant scale factor of %g.  Input val = %g, err = %g.\n", NP_val, NP_val, NP_err ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }

       if ( NP_val >=1 ) {
          printf("\n\n *** ra2bRoostatsClass3D_3b::makeBetaConstraint:  warning.  Input efficiency value is %g.  Resetting to 0.999.\n\n", NP_val) ;
          NP_val = 0.999 ;
       }

       double alpha, beta ;
       char varname[1000] ;


       betaModeTransform( NP_val, NP_err, alpha, beta ) ;

       sprintf( varname, "alpha_%s", NP_name ) ;
       RooAbsReal* rar_alpha = new RooRealVar( varname, varname, alpha, 0., 1e6 ) ;
       ((RooRealVar*)rar_alpha) -> setConstant(kTRUE) ;

       sprintf( varname, "beta_%s", NP_name ) ;
       RooAbsReal* rar_beta = new RooRealVar( varname, varname, beta, 0., 1e6 ) ;
       ((RooRealVar*)rar_beta) -> setConstant(kTRUE) ;

       RooAbsReal* rar_np = new RooRealVar( NP_name, NP_name, NP_val, 0., 1. ) ;

       sprintf( varname, "betapdf_%s", NP_name ) ;
       RooAbsReal* rar_pdf = new RooBetaPdf( varname, varname, *rar_np, *rar_alpha, *rar_beta ) ;

       allNuisances -> add( *rar_np ) ;
       allNuisancePdfs -> add( *rar_pdf ) ;
       globalObservables -> add( *rar_alpha ) ;
       globalObservables -> add( *rar_beta ) ;

       //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

       char vname[1000] ;
       sprintf( vname, "mean_%s", NP_name ) ;
       RooConstVar* g_mean  = new RooConstVar( vname, vname, NP_val ) ;
       sprintf( vname, "sigma_%s", NP_name ) ;
       RooConstVar* g_sigma = new RooConstVar( vname, vname, NP_err ) ;
       workspace.import( *g_mean  ) ;
       workspace.import( *g_sigma ) ;


       return rar_np ;



    } // makeBetaConstraint.


   //==============================================================================================================




    RooAbsReal* ra2bRoostatsClass3D_3b::makeCorrelatedBetaPrimeConstraint( 
            const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign ) {



       if ( NP_err <= 0. ) {
          printf("  Uncertainty is zero.  Will return constant scale factor of %g.  Input val = %g, err = %g.\n", NP_val, NP_val, NP_err ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }


       double alpha, beta ;
       char varname[1000] ;

       betaPrimeModeTransform( NP_val, NP_err, alpha, beta ) ;




       sprintf( varname, "failPar_%s", NP_base_name ) ;
       RooAbsReal* baseFailPar = (RooAbsReal*) allNuisances -> find( varname ) ;
       if ( baseFailPar != 0x0 ) {
          printf("  base fail parameter : %s = %g\n", baseFailPar->GetName(), baseFailPar->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_3b::makeCorrelatedBetaPrimeConstraint: can't find base fail parameter.\n\n") ;
          return 0x0 ;
       }

       sprintf( varname, "passPar_%s", NP_base_name ) ;
       RooAbsReal* basePassPar = (RooAbsReal*) allNuisances -> find( varname ) ;
       if ( basePassPar != 0x0 ) {
          printf("  base pass parameter : %s = %g\n", basePassPar->GetName(), basePassPar->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_3b::makeCorrelatedBetaPrimeConstraint: can't find base pass parameter.\n\n") ;
          return 0x0 ;
       }





       sprintf( varname, "failObs_%s", NP_base_name ) ;
       RooAbsReal* baseFailObs = (RooAbsReal*) globalObservables -> find( varname ) ;
       if ( baseFailObs != 0x0 ) {
          printf("  base fail obs : %s = %g\n", baseFailObs->GetName(), baseFailObs->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_3b::makeCorrelatedBetaPrimeConstraint: can't find base fail obs.\n\n") ;
          return 0x0 ;
       }

       sprintf( varname, "passObs_%s", NP_base_name ) ;
       RooAbsReal* basePassObs = (RooAbsReal*) globalObservables -> find( varname ) ;
       if ( basePassObs != 0x0 ) {
          printf("  base pass obs : %s = %g\n", basePassObs->GetName(), basePassObs->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_3b::makeCorrelatedBetaPrimeConstraint: can't find base pass obs.\n\n") ;
          return 0x0 ;
       }







       if ( basePassObs->getVal() < 0. ) {
          printf("\n\n\n *** makeCorrelatedBetaPrimeConstraint : illegal base pass observable value : %g\n\n",
             basePassObs->getVal() ) ; cout << flush ;
          return 0x0 ;
       }

       if ( baseFailObs->getVal() < 0. ) {
          printf("\n\n\n *** makeCorrelatedBetaPrimeConstraint : illegal base fail observable value : %g\n\n",
             baseFailObs->getVal() ) ; cout << flush ;
          return 0x0 ;
       }



       char formula[10000] ;
       RooAbsReal *passPar(0x0), *failPar(0x0) ;

       if ( !changeSign ) {

          sprintf( formula, "%g+(@0-%g)*(%g)", (alpha-1), basePassObs->getVal(), sqrt((alpha-1)/( basePassObs->getVal() )) ) ;
          printf(" creating transformed pass parameter with formula : %s\n", formula ) ;
          sprintf( varname, "passPar_%s", NP_name ) ;
          passPar = new RooFormulaVar( varname, formula, RooArgSet( *basePassPar ) ) ;

          ////// sprintf( formula, "%g+(@0-%g)*(%g)", (beta-1), baseFailObs->getVal(), sqrt((beta-1)/( baseFailObs->getVal() )) ) ;
          sprintf( formula, "%g+(@0-%g)*(%g)", (beta-1), baseFailObs->getVal(), sqrt((beta+1)/( baseFailObs->getVal() )) ) ;
          printf(" creating transformed fail parameter with formula : %s\n", formula ) ;
          sprintf( varname, "failPar_%s", NP_name ) ;
          failPar = new RooFormulaVar( varname, formula, RooArgSet( *baseFailPar ) ) ;

       } else {

          sprintf( formula, "%g-(@0-%g)*(%g)", (alpha-1), basePassObs->getVal(), sqrt((alpha-1)/( basePassObs->getVal() )) ) ;
          printf(" creating transformed pass parameter with formula : %s\n", formula ) ;
          sprintf( varname, "passPar_%s", NP_name ) ;
          passPar = new RooFormulaVar( varname, formula, RooArgSet( *basePassPar ) ) ;

          ////// sprintf( formula, "%g-(@0-%g)*(%g)", (beta-1), baseFailObs->getVal(), sqrt((beta-1)/( baseFailObs->getVal() )) ) ;
          sprintf( formula, "%g-(@0-%g)*(%g)", (beta-1), baseFailObs->getVal(), sqrt((beta+1)/( baseFailObs->getVal() )) ) ;
          printf(" creating transformed fail parameter with formula : %s\n", formula ) ;
          sprintf( varname, "failPar_%s", NP_name ) ;
          failPar = new RooFormulaVar( varname, formula, RooArgSet( *baseFailPar ) ) ;

       }



       sprintf( varname, "%s", NP_name ) ;
       RooAbsReal* rar = new RooRatio( varname, varname, *passPar, *failPar ) ;
       printf("  Created nuisance parameter %s : val = %g\n", varname, rar -> getVal() ) ;

       return rar ;

    } // makeCorrelatedBetaPrimeConstraint.


   //==============================================================================================================


    RooAbsReal* ra2bRoostatsClass3D_3b::makeCorrelatedBetaConstraint(
        const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign ) {


       if ( NP_err <= 0. ) {
          printf("  Uncertainty is zero.  Will return constant scale factor of %g.  Input val = %g, err = %g.\n", NP_val, NP_val, NP_err ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }



       double alpha, beta ;
       char varname[1000] ;

       betaModeTransform( NP_val, NP_err, alpha, beta ) ;




       sprintf( varname, "failPar_%s", NP_base_name ) ;
       RooAbsReal* baseFailPar = (RooAbsReal*) allNuisances -> find( varname ) ;
       if ( baseFailPar != 0x0 ) {
          printf("  base fail parameter : %s = %g\n", baseFailPar->GetName(), baseFailPar->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_3b::makeCorrelatedBetaPrimeConstraint: can't find base fail parameter.\n\n") ;
          return 0x0 ;
       }

       sprintf( varname, "passPar_%s", NP_base_name ) ;
       RooAbsReal* basePassPar = (RooAbsReal*) allNuisances -> find( varname ) ;
       if ( basePassPar != 0x0 ) {
          printf("  base pass parameter : %s = %g\n", basePassPar->GetName(), basePassPar->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_3b::makeCorrelatedBetaPrimeConstraint: can't find base pass parameter.\n\n") ;
          return 0x0 ;
       }





       sprintf( varname, "failObs_%s", NP_base_name ) ;
       RooAbsReal* baseFailObs = (RooAbsReal*) globalObservables -> find( varname ) ;
       if ( baseFailObs != 0x0 ) {
          printf("  base fail obs : %s = %g\n", baseFailObs->GetName(), baseFailObs->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_3b::makeCorrelatedBetaPrimeConstraint: can't find base fail obs.\n\n") ;
          return 0x0 ;
       }

       sprintf( varname, "passObs_%s", NP_base_name ) ;
       RooAbsReal* basePassObs = (RooAbsReal*) globalObservables -> find( varname ) ;
       if ( basePassObs != 0x0 ) {
          printf("  base pass obs : %s = %g\n", basePassObs->GetName(), basePassObs->getVal() ) ;
       } else {
          printf("\n\n *** ra2bRoostatsClass3D_3b::makeCorrelatedBetaPrimeConstraint: can't find base pass obs.\n\n") ;
          return 0x0 ;
       }






       if ( basePassObs->getVal() < 0. ) {
          printf("\n\n\n *** makeCorrelatedBetaPrimeConstraint : illegal base pass observable value : %g\n\n",
             basePassObs->getVal() ) ; cout << flush ;
          return 0x0 ;
       }

       if ( baseFailObs->getVal() < 0. ) {
          printf("\n\n\n *** makeCorrelatedBetaPrimeConstraint : illegal base fail observable value : %g\n\n",
             baseFailObs->getVal() ) ; cout << flush ;
          return 0x0 ;
       }



       char formula[10000] ;
       RooAbsReal *passPar(0x0), *failPar(0x0) ;

       if ( !changeSign ) {

          sprintf( formula, "%g+(@0-%g)*(%g)", (alpha-1), basePassObs->getVal(), sqrt((alpha-1)/( basePassObs->getVal() )) ) ;
          printf(" creating transformed pass parameter with formula : %s\n", formula ) ;
          sprintf( varname, "passPar_%s", NP_name ) ;
          passPar = new RooFormulaVar( varname, formula, RooArgSet( *basePassPar ) ) ;

          sprintf( formula, "%g+(@0-%g)*(%g)", (beta-1), baseFailObs->getVal(), sqrt((beta-1)/( baseFailObs->getVal() )) ) ;
          printf(" creating transformed fail parameter with formula : %s\n", formula ) ;
          sprintf( varname, "failPar_%s", NP_name ) ;
          failPar = new RooFormulaVar( varname, formula, RooArgSet( *baseFailPar ) ) ;

       } else {

          sprintf( formula, "%g-(@0-%g)*(%g)", (alpha-1), basePassObs->getVal(), sqrt((alpha-1)/( basePassObs->getVal() )) ) ;
          printf(" creating transformed pass parameter with formula : %s\n", formula ) ;
          sprintf( varname, "passPar_%s", NP_name ) ;
          passPar = new RooFormulaVar( varname, formula, RooArgSet( *basePassPar ) ) ;

          sprintf( formula, "%g-(@0-%g)*(%g)", (beta-1), baseFailObs->getVal(), sqrt((beta-1)/( baseFailObs->getVal() )) ) ;
          printf(" creating transformed fail parameter with formula : %s\n", formula ) ;
          sprintf( varname, "failPar_%s", NP_name ) ;
          failPar = new RooFormulaVar( varname, formula, RooArgSet( *baseFailPar ) ) ;

       }







       sprintf( varname, "%s_passPlusFail", NP_name ) ;
       RooAddition* passPlusFail = new RooAddition( varname, varname, RooArgSet(*passPar, *failPar) ) ;


       sprintf( varname, "%s", NP_name ) ;
       RooAbsReal* rar = new RooRatio( varname, varname, *passPar, *passPlusFail ) ;
       printf("  Created nuisance parameter %s : val = %g\n", varname, rar -> getVal() ) ;

       return rar ;

    } // makeCorrelatedBetaPrimeConstraint.


   //==============================================================================================================


    RooAbsReal* ra2bRoostatsClass3D_3b::makeGaussianConstraint( const char* NP_name, double NP_val, double NP_err, bool allowNegative ) {

       if ( NP_err <= 0. ) {
          printf(" makeGaussianConstraint:  Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }


       double max = NP_val + 6.*NP_err ;
       double min = NP_val - 6.*NP_err ;

       if ( min < 0. && !allowNegative ) { min = 1e-5 ; }

       RooRealVar* np_rrv = new RooRealVar( NP_name, NP_name, min, max ) ;
       np_rrv -> setVal( NP_val ) ;
       np_rrv -> setConstant( kFALSE ) ;

       //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

       char vname[1000] ;
       sprintf( vname, "mean_%s", NP_name ) ;
       RooRealVar* g_mean = new RooRealVar( vname, vname, NP_val, -1000., 1000. ) ;
       g_mean->setConstant(kTRUE);
       sprintf( vname, "sigma_%s", NP_name ) ;
       RooConstVar* g_sigma = new RooConstVar( vname, vname, NP_err ) ;

       char pdfname[1000] ;
       sprintf( pdfname, "pdf_%s", NP_name ) ;
       RooGaussian* np_pdf = new RooGaussian( pdfname, pdfname, *np_rrv, *g_mean, *g_sigma ) ;

       allNuisances -> add( *np_rrv ) ;
       allNuisancePdfs -> add( *np_pdf ) ;
       globalObservables -> add( *g_mean ) ;

       printf("  makeGaussianConstraint : created nuisance parameter %s : val = %g\n", NP_name, np_rrv -> getVal() ) ;

       return np_rrv ;


    } // makeGaussianConstraint.


   //==============================================================================================================


    RooAbsReal* ra2bRoostatsClass3D_3b::makeLognormalConstraint( const char* NP_name, double NP_val, double NP_err ) {

       if ( NP_err <= 0. ) {
          printf(" makeLognormalConstraint:  Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }


       char pname[1000] ;
       sprintf( pname, "prim_%s", NP_name ) ;

       printf(" makeLognormalConstraint : creating primary log-normal variable %s\n", pname ) ;
       RooRealVar* np_prim_rrv = new RooRealVar( pname, pname, 0., -6., 6. ) ;
       np_prim_rrv -> setVal( 0. ) ;
       np_prim_rrv -> setConstant( kFALSE ) ;

       sprintf( pname, "prim_mean_%s", NP_name ) ;
       RooRealVar* np_prim_mean = new RooRealVar( pname, pname, 0., -6., 6. ) ;
       np_prim_mean->setConstant(kTRUE) ;

       sprintf( pname, "prim_sigma_%s", NP_name ) ;
       RooConstVar* np_prim_sigma = new RooConstVar( pname, pname, 1. ) ;


       char pdfname[1000] ;
       sprintf( pdfname, "pdf_prim_%s", NP_name ) ;
       RooGaussian* np_prim_pdf = new RooGaussian( pdfname, pdfname, *np_prim_rrv, *np_prim_mean, *np_prim_sigma ) ;

       allNuisances -> add( *np_prim_rrv ) ;
       allNuisancePdfs -> add( *np_prim_pdf ) ;
       globalObservables -> add( *np_prim_mean ) ;


       //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

       char vname[1000] ;
       sprintf( vname, "mean_%s", NP_name ) ;
       RooConstVar* g_mean  = new RooConstVar( vname, vname, NP_val ) ;
       sprintf( vname, "sigma_%s", NP_name ) ;
       RooConstVar* g_sigma = new RooConstVar( vname, vname, NP_err ) ;

       //-- compute the log-normal-distributed parameter from the primary parameter.

       //--- This is the old (Fedor) way. ---------------------------------------------------------
       ////// RooFormulaVar* np_rfv = new RooFormulaVar( NP_name, "@0 * pow( exp( @1/@0 ), @2)",
       //////           RooArgSet( *g_mean, *g_sigma, *np_prim_rrv ) ) ;
       //------------------------------------------------------------------------------------------

       //--- This is the new way.  RMS of lognormal is much closer to sigma when sigma is
       //    large, doing it this way.  When sigma/mean is small, they are about the same.
       //    That is, exp(sigma/mean) is close to (sigma/mean + 1).  This one is better when
       //    sigma/mean is not small.  The high-side tail is not as strong.
       //
        RooFormulaVar* np_rfv = new RooFormulaVar( NP_name, "@0 * pow( ( @1/@0 + 1. ), @2)",
                  RooArgSet( *g_mean, *g_sigma, *np_prim_rrv ) ) ;
       //------------------------------------------------------------------------------------------



       printf("  makeLognormalConstraint : created log-normal nuisance parameter %s : val = %g\n", NP_name, np_rfv -> getVal() ) ;

       return np_rfv ;


    } // makeLognormalConstraint.


   //==============================================================================================================



    RooAbsReal* ra2bRoostatsClass3D_3b::makeCorrelatedGaussianConstraint(
            const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign, bool allowNegative ) {

       if ( NP_err <= 0. ) {
          printf("  makeCorrelatedGaussianConstraint: Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }

       RooRealVar* rrv_np_base_par = (RooRealVar*) allNuisances -> find( NP_base_name ) ;

       if ( rrv_np_base_par == 0x0 ) {

          printf("\n\n makeCorrelatedGaussianConstraint : creating base nuisance parameter - %s\n\n", NP_base_name ) ;
          rrv_np_base_par = new RooRealVar( NP_base_name, NP_base_name, -6.0, 6.0 ) ;
          rrv_np_base_par -> setVal( 0. ) ;
          rrv_np_base_par -> setConstant( kFALSE ) ;
          allNuisances -> add( *rrv_np_base_par ) ;

          char vname[1000] ;
          sprintf( vname, "mean_%s", NP_base_name ) ;
          RooRealVar* g_mean = new RooRealVar( vname, vname, 0.0,-1000.,1000. ) ;
	  g_mean->setConstant(kTRUE);
          sprintf( vname, "sigma_%s", NP_base_name ) ;
          RooConstVar* g_sigma = new RooConstVar( vname, vname, 1.0 ) ;

          char pdfname[100] ;
          sprintf( pdfname, "pdf_%s", NP_base_name ) ;
          printf("\n\n makeCorrelatedGaussianConstraint : creating base nuisance parameter pdf - %s\n\n", pdfname ) ;
          ///// RooGaussian* base_np_pdf = new RooGaussian( pdfname, pdfname, *rrv_np_base_par, RooConst(0.), RooConst(1.) ) ;
          RooGaussian* base_np_pdf = new RooGaussian( pdfname, pdfname, *rrv_np_base_par, *g_mean, *g_sigma ) ;
          allNuisancePdfs -> add( *base_np_pdf ) ;
          globalObservables -> add( *g_mean ) ;

       }

       //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

       char vname[1000] ;
       sprintf( vname, "mean_%s", NP_name ) ;
       RooConstVar* g_mean = new RooConstVar( vname, vname, NP_val ) ;
       sprintf( vname, "sigma_%s", NP_name ) ;
       RooConstVar* g_sigma = new RooConstVar( vname, vname, NP_err ) ;

       RooAbsReal* rar(0x0) ;

       if ( allowNegative ) {

          char formula[1000] ;

          if ( !changeSign ) {
             sprintf( formula, "@0+@1*@2" ) ;
          } else {
             sprintf( formula, "@0-@1*@2" ) ;
          }

          rar = new RooFormulaVar( NP_name, formula, RooArgSet( *g_mean, *g_sigma, *rrv_np_base_par ) ) ;

          printf(" makeCorrelatedGaussianConstraint : creating correlated gaussian NP with formula : %s,  %s, val = %g\n", formula, NP_name, rar->getVal() ) ;

       } else {

          rar = new RooPosDefCorrGauss( NP_name, NP_name, *g_mean, *g_sigma, *rrv_np_base_par, changeSign ) ;

          printf(" makeCorrelatedGaussianConstraint : creating pos-def correlated gaussian NP  :  %s, val = %g, err = %g\n", NP_name, rar->getVal(), NP_err ) ;

       }




       return rar ;

    } // makeCorrelatedGaussianConstraint.


   //==============================================================================================================



    RooAbsReal* ra2bRoostatsClass3D_3b::makeCorrelatedLognormalConstraint(
            const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign ) {

       if ( NP_err <= 0. ) {
          printf("  makeCorrelatedLognormalConstraint: Uncertainty is zero.  Will return constant scale factor of %g for %s.  Input val = %g, err = %g.\n", NP_val, NP_name, NP_val, NP_err ) ;
          return new RooConstVar( NP_name, NP_name, NP_val ) ;
       }

       char prim_name[1000] ;
       sprintf( prim_name, "prim_%s", NP_base_name ) ;
       RooRealVar* rrv_np_base_par = (RooRealVar*) allNuisances -> find( prim_name ) ;
       
       if ( rrv_np_base_par == 0x0 ) {
	 
	 printf("\n\n makeCorrelatedLognormalConstraint : creating base nuisance parameter - %s\n\n", prim_name ) ;
	 rrv_np_base_par = new RooRealVar( prim_name, prim_name, -6.0, 6.0 ) ;
	 rrv_np_base_par -> setVal( 0. ) ;
	 rrv_np_base_par -> setConstant( kFALSE ) ;
	 allNuisances -> add( *rrv_np_base_par ) ;
	 
	 char vname[1000] ;
	 // sprintf( vname, "mean_%s", NP_base_name ) ;
	 sprintf( vname, "prim_mean_%s", NP_base_name ) ;
	 RooRealVar* g_mean = new RooRealVar( vname, vname, 0.0,-10.,10. ) ;
	 g_mean->setConstant(kTRUE);
	 // sprintf( vname, "sigma_%s", NP_base_name ) ;
	 sprintf( vname, "prim_sigma_%s", NP_base_name ) ;
	 RooConstVar* g_sigma = new RooConstVar( vname, vname, 1.0 ) ;
	 
	 char pdfname[100] ;
	 // sprintf( pdfname, "pdf_%s", NP_base_name ) ;
	 sprintf( pdfname, "pdf_prim_%s", NP_base_name ) ;
	 printf("\n\n makeCorrelatedLognormalConstraint : creating base nuisance parameter pdf - %s\n\n", pdfname ) ;
	 RooGaussian* base_np_pdf = new RooGaussian( pdfname, pdfname, *rrv_np_base_par, *g_mean, *g_sigma ) ;
	 
	 allNuisancePdfs -> add( *base_np_pdf ) ;
	 globalObservables -> add( *g_mean ) ;
	 
       }


       //-- create const variables for mean and sigma so that they can be saved and accessed from workspace later.

       char vname[1000] ;
       sprintf( vname, "mean_%s", NP_name ) ;
       RooConstVar* ln_mean  = new RooConstVar( vname, vname, NP_val ) ;
       
       sprintf( vname, "sigma_%s", NP_name ) ;
       RooConstVar* ln_sigma = new RooConstVar( vname, vname, NP_err ) ;
       
       
       RooAbsReal* rar(0x0) ;
       
       
       char formula[1000] ;
       
       if ( !changeSign ) {
	 ////// sprintf( formula, "@0 * pow( exp( @1/@0 ), @2 )" ) ;
	 sprintf( formula, "@0 * pow( ( @1/@0 + 1.), @2 )" ) ;
       } else {
	 ////// sprintf( formula, "@0 * pow( exp( @1/@0 ), -1.0 * @2 )" ) ;
	 sprintf( formula, "@0 * pow( ( @1/@0 + 1.), -1.0 * @2 )" ) ;
       }
       
       rar = new RooFormulaVar( NP_name, formula, RooArgSet( *ln_mean, *ln_sigma, *rrv_np_base_par ) ) ;
       
       printf(" makeCorrelatedLognormalConstraint : creating correlated log-normal NP with formula : %s,  %s, val = %g, mean=%g, sigma=%g\n", formula, NP_name, rar->getVal(), NP_val, NP_err ) ;
       

       return rar ;

    } // makeCorrelatedLognormalConstraint.


   //==============================================================================================================

   // copy-pasting here the helper functions from Will Reece

  void ra2bRoostatsClass3D_3b::SetConstants(RooWorkspace pWs, RooStats::ModelConfig pMc){
    //
    // Fix all variables in the PDF except observables, POI and
    // nuisance parameters. Note that global observables are fixed.
    // If you need global observables floated, you have to set them
    // to float separately.
    //

    pMc.SetWorkspace(pWs);

    RooAbsPdf * pPdf = pMc.GetPdf(); // we do not own this

    RooArgSet * pVars = pPdf->getVariables(); // we do own this

    RooArgSet * pFloated = new RooArgSet(*pMc.GetObservables());
    pFloated->add(*pMc.GetParametersOfInterest());
    pFloated->add(*pMc.GetNuisanceParameters());

    TIterator * pIter = pVars->createIterator(); // we do own this

    for(TObject * pObj = pIter->Next(); pObj; pObj = pIter->Next() ){
      std::string _name = pObj->GetName();
      RooRealVar * pFloatedObj = (RooRealVar *)pFloated->find(_name.c_str());
      if (pFloatedObj){
        ((RooRealVar *)pObj)->setConstant(kFALSE);
      }
      else{
        ((RooRealVar *)pObj)->setConstant(kTRUE);
      }
    }

    delete pIter;
    delete pVars;
    delete pFloated;

    return;

}

   //==============================================================================================================


  void ra2bRoostatsClass3D_3b::SetConstant(const RooArgSet * vars, Bool_t value ){
    //
    // Set the constant attribute for all vars in the set
    //

    TIterator * pIter = vars->createIterator(); // we do own this

    for(TObject * pObj = pIter->Next(); pObj; pObj = pIter->Next() ){
      ((RooRealVar *)pObj)->setConstant(value);
    }

    delete pIter;

    return;
  }

 //===============================================================================================================


  bool ra2bRoostatsClass3D_3b::setupShapeSyst( const char* infile,
                                               const char* systname,
                                               int constraintType,
                                               double target_mgl, double target_mlsp,
                                               RooWorkspace& workspace
                                               ) {

      printf("\n\n\n setupShapeSyst :  setting up %s systematic.  Input file %s\n\n", systname, infile ) ;

      if ( constraintType == 1 ) {
         printf(" setupShapeSyst : Constraint type for %s : Gaussian (1).\n\n", systname ) ;
      } else if ( constraintType == 2 ) {
         printf(" setupShapeSyst : Constraint type for %s : log-normal (2).\n\n", systname ) ;
      } else {
         printf("  *** setupShapeSyst : Constraint type %d not implemented.\n\n", constraintType ) ;
         return false ;
      }

      char command[1000] ;

      sprintf( command, "head -1 %s | awk '{print NF}'", infile ) ;
      const char* nfields_str = gSystem->GetFromPipe( command ) ;
      int nfields ;
      sscanf( nfields_str, "%d", &nfields ) ;
      printf(" setupShapeSyst: Nfields in %s is %d\n", infile, nfields ) ;


      bool hasSL(false) ;
      bool hasSLSig(false) ;

      int ArraySize ;
      if ( nfields == 2+2*(nBinsMET*nBinsHT*nBinsBtag) ) {
         hasSL = false ;
         printf("\n\n setupShapeSyst : %s : Format is consistent with no SL observables.\n", systname ) ;
         ArraySize = nfields ;
      } else if ( nfields == 2+3*(nBinsMET*nBinsHT*nBinsBtag) ) {
         hasSL = true ;
         printf("\n\n setupShapeSyst : %s : Format is consistent with including SL observables.\n", systname ) ;
         ArraySize = nfields ;
      } else if ( nfields == 2+4*(nBinsMET*nBinsHT*nBinsBtag) ) {
         hasSLSig = true ;
         printf("\n\n setupShapeSyst : %s : Format is consistent with including SL and SLSig observables.\n", systname ) ;
         ArraySize = nfields ;
      } else {
         printf("\n\n setupShapeSyst : %s : I don't know what to do with nfields = %d\n\n", systname, nfields ) ;
         return false ;
      }


      ifstream infq ;
      infq.open(infile) ;
      if ( !infq.good() ) {
         printf("\n\n *** setupShapeSyst: Problem opening input file: %s.\n\n", infile ) ;
         return false ;
      }


      double syst_zl[10][10][10] ;
      double syst_slSig[10][10][10] ;
      double syst_sl[10][10][10] ;
      double syst_ldp[10][10][10] ;

      double minSyst(0.) ;
      double maxSyst(0.) ;

      double matchArrayContent[ArraySize] ;
      double nearbyMatchArrayContent[ArraySize] ;

      bool foundMatch = false ;
      bool foundNearbyMatch = false ;

      int nBins = nBinsMET*nBinsHT*nBinsBtag ;

      while ( infq.good() ) {


         double readArrayContent[ArraySize] ;
         for ( int i=0; infq && i<ArraySize; ++ i) {
            infq >> readArrayContent[i] ;
         }

         double mgl  = readArrayContent[0] ;
         double mlsp = readArrayContent[1] ;

         if ( fabs( mgl-target_mgl ) < 1. && fabs( mlsp - target_mlsp ) < 1. ) {

            for ( int i=0; i<ArraySize; i++ ) { matchArrayContent[i] = readArrayContent[i] ; }
            foundMatch = true ;
            printf("\n\n setupShapeSyst : %s :Found mgl=%.0f, mlsp=%.0f\n\n", systname, mgl, mlsp ) ;
            break ;

         }

         if ( fabs( mgl-target_mgl ) < 26. && fabs( mlsp - target_mlsp ) < 26. && !foundNearbyMatch ) {

            for ( int i=0; i<ArraySize; i++ ) { nearbyMatchArrayContent[i] = readArrayContent[i] ; }
            foundNearbyMatch = true ;
            printf("\n\n setupShapeSyst : %s : Found nearby match mgl=%.0f, mlsp=%.0f  (requested mgl=%.0f, mlsp=%.0f)\n\n", systname, mgl, mlsp, target_mgl, target_mlsp ) ;

         }

      } // reading file?


      double ArrayContent[ArraySize] ;

      if ( foundMatch ) {
         for ( int i=0; i<ArraySize; i++ ) { ArrayContent[i] = matchArrayContent[i] ; }
      } else if ( foundNearbyMatch ) {
         for ( int i=0; i<ArraySize; i++ ) { ArrayContent[i] = nearbyMatchArrayContent[i] ; }
      } else {
         printf("\n\n *** setupShapeSyst : %s : Did not find match or nearby match for mgl=%.0f, mlsp=%.0f\n\n", systname, target_mgl, target_mlsp ) ;
         return false ;
      }

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {

               if ( hasSLSig ) {
                  syst_zl   [mbi][hbi][bbi]  = ArrayContent[2 + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
                  syst_slSig[mbi][hbi][bbi]  = ArrayContent[2 + nBins + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
                  syst_sl   [mbi][hbi][bbi]  = ArrayContent[2 + 2*nBins + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
                  syst_ldp  [mbi][hbi][bbi]  = ArrayContent[2 + 3*nBins + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
               } else if ( hasSL ) {
                  syst_zl   [mbi][hbi][bbi]  = ArrayContent[2 + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
		  syst_slSig[mbi][hbi][bbi]  = 0. ;
                  syst_sl   [mbi][hbi][bbi]  = ArrayContent[2 + nBins + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
                  syst_ldp  [mbi][hbi][bbi]  = ArrayContent[2 + 2*nBins + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
               } else {
                  syst_zl   [mbi][hbi][bbi]  = ArrayContent[2 + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
                  syst_slSig[mbi][hbi][bbi]  = 0. ;
                  syst_sl   [mbi][hbi][bbi]  = 0. ;
                  syst_ldp  [mbi][hbi][bbi]  = ArrayContent[2 + nBins + mbi*(nBinsHT*nBinsBtag) + hbi*(nBinsBtag) + bbi] ;
               }

               if ( syst_zl   [mbi][hbi][bbi] > maxSyst ) maxSyst = syst_zl   [mbi][hbi][bbi] ;
               if ( syst_slSig[mbi][hbi][bbi] > maxSyst ) maxSyst = syst_slSig[mbi][hbi][bbi] ;
               if ( syst_sl   [mbi][hbi][bbi] > maxSyst ) maxSyst = syst_sl   [mbi][hbi][bbi] ;
               if ( syst_ldp  [mbi][hbi][bbi] > maxSyst ) maxSyst = syst_ldp  [mbi][hbi][bbi] ;

               if ( syst_zl   [mbi][hbi][bbi] < minSyst ) minSyst = syst_zl   [mbi][hbi][bbi] ;
               if ( syst_slSig[mbi][hbi][bbi] < minSyst ) minSyst = syst_slSig[mbi][hbi][bbi] ;
               if ( syst_sl   [mbi][hbi][bbi] < minSyst ) minSyst = syst_sl   [mbi][hbi][bbi] ;
               if ( syst_ldp  [mbi][hbi][bbi] < minSyst ) minSyst = syst_ldp  [mbi][hbi][bbi] ;

            } // bbi.
         } // hbi.
      } // mbi.


      printf("\n\n setupShapeSyst: %s : Min syst = %6.2f, Max syst = %6.2f\n\n", systname, minSyst, maxSyst ) ;

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {

               char pname[100] ;
               bool changeSign ;
               RooAbsReal* rar_par ;

               sprintf( pname, "%s_M%d_H%d_%db", systname, mbi+1, hbi+1, bbi+1 ) ;
               if ( syst_zl[mbi][hbi][bbi] < 0 ) { changeSign = true ; } else { changeSign = false ; }
               if ( constraintType == 1 ) {
                  rar_par = makeCorrelatedGaussianConstraint(  pname, 1.0, fabs(syst_zl[mbi][hbi][bbi]), systname, changeSign ) ;
               } else if ( constraintType == 2 ) {
                  rar_par = makeCorrelatedLognormalConstraint( pname, 1.0, fabs(syst_zl[mbi][hbi][bbi]), systname, changeSign ) ;
               }
               cout << flush ;
               workspace.import( *rar_par ) ;

               sprintf( pname, "%s_slSig_M%d_H%d_%db", systname, mbi+1, hbi+1, bbi+1 ) ;
               if ( syst_slSig[mbi][hbi][bbi] < 0 ) { changeSign = true ; } else { changeSign = false ; }
               if ( constraintType == 1 ) {
                  rar_par = makeCorrelatedGaussianConstraint(  pname, 1.0, fabs(syst_slSig[mbi][hbi][bbi]), systname, changeSign ) ;
               } else if ( constraintType == 2 ) {
                  rar_par = makeCorrelatedLognormalConstraint( pname, 1.0, fabs(syst_slSig[mbi][hbi][bbi]), systname, changeSign ) ;
               }
               cout << flush ;
               workspace.import( *rar_par ) ;

               sprintf( pname, "%s_sl_M%d_H%d_%db", systname, mbi+1, hbi+1, bbi+1 ) ;
               if ( syst_sl[mbi][hbi][bbi] < 0 ) { changeSign = true ; } else { changeSign = false ; }
               if ( constraintType == 1 ) {
                  rar_par = makeCorrelatedGaussianConstraint(  pname, 1.0, fabs(syst_sl[mbi][hbi][bbi]), systname, changeSign ) ;
               } else if ( constraintType == 2 ) {
                  rar_par = makeCorrelatedLognormalConstraint( pname, 1.0, fabs(syst_sl[mbi][hbi][bbi]), systname, changeSign ) ;
               }
               cout << flush ;
               workspace.import( *rar_par ) ;

               sprintf( pname, "%s_ldp_M%d_H%d_%db", systname, mbi+1, hbi+1, bbi+1 ) ;
               if ( syst_ldp[mbi][hbi][bbi] < 0 ) { changeSign = true ; } else { changeSign = false ; }
               if ( constraintType == 1 ) {
                  rar_par = makeCorrelatedGaussianConstraint(  pname, 1.0, fabs(syst_ldp[mbi][hbi][bbi]), systname, changeSign ) ;
               } else if ( constraintType == 2 ) {
                  rar_par = makeCorrelatedLognormalConstraint( pname, 1.0, fabs(syst_ldp[mbi][hbi][bbi]), systname, changeSign ) ;
               }
               cout << flush ;
               workspace.import( *rar_par ) ;

            } // bbi.
         } // hbi.
      } // mbi.


      return true ;

  } // setupShapeSyst


 //===============================================================================================================















