#include "ra2bRoostatsClass3D_3b.h"

//
//  Owen: The difference between ra2bRoostatsClass3D_1 and ra2bRoostatsClass3D_3b
//        is in the implementation of the nuisance parameters.
//
//          ra2bRoostatsClass3D_1 : uses log normals
//
//          ra2bRoostatsClass3D_2b : Uses beta and beta prime distributions.
//                                   Also has total 0lep susy yield as POI
//                                   instead of 1,1,1 susy yield.
//
//          ra2bRoostatsClass3D_3b : Floating parameter for ttwj is ZL/SL ratio.
//                                  LDP content for ttwj and znn comes from
//                                  MC ZL/LDP ratios, not absolute MC prediction.
//                                  QCD model is one floating ZL/LDP ratio with
//                                  MC corrections applied.
//
//                                  NP PDFs are truncated Gaussians.
//
//                                  Added QCD Model 4.  Float 0lep/LDP ratio in
//                                  each HT bin.  Apply MET and nbjet corrections.
//
//
//            QCD Model descriptions
//
//               3 : Float one global 0lep/LDP ratio
//
//               2 : Float independent 0lep/LDP ratios for each HT bin.
//
//               4 : Float independent 0lep/LDP ratios for each HT bin.
//                   Float independent MET corrections for each MET bin > 1.
//                   Float independent nbjet corrections for each MET bin > 1.
//
//               5 : Float independent 0lep/LDP ratios for each HT bin.
//                   Float one MET correction that's common for bins with index > 1.
//                   Float one nbjet correction that's common for bins with index > 1.
//
//
//


#include <iostream>
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
					     int   qcdModelIndex,
					     const char* wsrootfilename,
					     const char* blindBinsList,
					     const char* systFile1
					     ) {








      //-- Hardwire in to ignore the highest MET bin in the lowest HT bin.
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            ignoreBin[mbi][hbi] = false ;
         } // hbi
      } // mbi
      ignoreBin[nBinsMET-1][0] = true ;

      printf("\n\n *** Ignoring these bins in the analysis.\n\n") ;
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            if ( ignoreBin[mbi][hbi] ) {
               printf("  MET %d, HT %d\n", mbi+1, hbi+1 ) ;
            }
         } // hbi
      } // mbi
      printf("\n\n\n") ;









      if ( qcdModelIndex<2 || qcdModelIndex>5 ) {
         printf("\n\n *** Unsupported QCD model index: %d.  Try again with 2,3,4,5.\n\n", qcdModelIndex) ;
         return false ;
      }




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












      double dummy ;
      dummy = t1bbbbXsec ;

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
      float NVVmc_0lep[nBinsMET][nBinsHT][nBinsBtag] ;
      float NVVmc_1lep[nBinsMET][nBinsHT][nBinsBtag] ;
      float NVVmc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

      
      // other parameters:

      float R_lsb[nBinsHT][nBinsBtag] ;
      float R_lsb_err[nBinsHT][nBinsBtag] ;

      float acc_Zee[nBinsMET] ;
      float acc_Zee_err[nBinsMET] ;
      float acc_Zmm[nBinsMET] ;
      float acc_Zmm_err[nBinsMET] ;

      float eff_Zee[nBinsHT] ;
      float eff_Zee_err[nBinsHT] ;
      float eff_Zmm[nBinsHT] ;
      float eff_Zmm_err[nBinsHT] ;
      
      float knn_1b[nBinsMET] ;
      float knn_1b_err[nBinsMET] ;
      float knn_2b ;
      float knn_2b_err ;
      float knn_3b ;
      float knn_3b_err ;

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
      float sf_ee ;
      float sf_ee_err ;       
      float sf_mm ;
      float sf_mm_err ;       
      
      float btageff_err ;

      float trigeff_0L[nBinsMET][nBinsHT];
      float trigefferr_0L[nBinsMET][nBinsHT];
      float trigeff_1L[nBinsMET][nBinsHT];
      float trigefferr_1L[nBinsMET][nBinsHT];


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

            if ( ignoreBin[i][j] ) {
               N_0lep[i][j][k] = 0. ;
               continue ;
            }

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

            if ( ignoreBin[i][j] ) {
               N_1lep[i][j][k] = 0. ;
               continue ;
            }

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

            if ( ignoreBin[i][j] ) {
               N_ldp[i][j][k] = 0. ;
               continue ;
            }

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

          if ( ignoreBin[i][j] ) {
             N_Zee[i][j] = 0. ;
             continue ;
          }

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

          if ( ignoreBin[i][j] ) {
             N_Zmm[i][j] = 0. ;
             continue ;
          }

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
      
      for (int j = 0 ; j < nBinsHT ; j++) {
	    
	TString inPar = "Z_ee_eff"+sHbins[j] ;
	fscanf( infp, "%s %g", label, &eff_Zee[j] ) ;
	
	if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	cout << inPar << " = " << eff_Zee[j] << endl ;	  
	
	inPar += "_err" ;
	fscanf( infp, "%s %g", label, &eff_Zee_err[j] ) ;
	cout << inPar << " = " << eff_Zee_err[j] << endl ;	  

      }
      
      for (int j = 0 ; j < nBinsHT ; j++) {
	    
	TString inPar = "Z_mm_eff"+sHbins[j] ;
	fscanf( infp, "%s %g", label, &eff_Zmm[j] ) ;
	
	if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	cout << inPar << " = " << eff_Zmm[j] << endl ;	  
	
	inPar += "_err" ;
	fscanf( infp, "%s %g", label, &eff_Zmm_err[j] ) ;
	cout << inPar << " = " << eff_Zmm_err[j] << endl ;	  

      }


      // VL -> nominal scale factors
      
      // 1b scale factors have a MET dependence

      for (int i = 0 ; i < nBinsMET ; i++) {
	    
	TString inPar = "knn_1b"+sMbins[i] ;
	fscanf( infp, "%s %g", label, &knn_1b[i] ) ;
	
	if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	cout << inPar << " = " << knn_1b[i] << endl ;	  
	
	inPar += "_err" ;
	fscanf( infp, "%s %g", label, &knn_1b_err[i] ) ;
	cout << inPar << " = " << knn_1b_err[i] << endl ;	  

      }
      

      // 2b and 3b scale factors:

      fscanf( infp, "%s %g", label, &knn_2b ) ; cout << "knn_2b" << " = " << knn_2b << endl ;
      fscanf( infp, "%s %g", label, &knn_2b_err ) ; cout << "knn_2b_err" << " = " << knn_2b_err << endl ;
      fscanf( infp, "%s %g", label, &knn_3b ) ; cout << "knn_3b" << " = " << knn_3b << endl ;
      fscanf( infp, "%s %g", label, &knn_3b_err ) ; cout << "knn_3b_err" << " = " << knn_3b_err << endl ;


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

      fscanf( infp, "%s %g", label, &sf_ee ) ; cout << "sf_ee" << " = " << sf_ee << endl ;
      fscanf( infp, "%s %g", label, &sf_ee_err ) ; cout << "sf_ee_err" << " = " << sf_ee_err << endl ;
      fscanf( infp, "%s %g", label, &sf_mm ) ; cout << "sf_mm" << " = " << sf_mm << endl ;
      fscanf( infp, "%s %g", label, &sf_mm_err ) ; cout << "sf_mm_err" << " = " << sf_mm_err << endl ;

      
      // btag efficiency error

      fscanf( infp, "%s %g", label, &btageff_err ) ; cout << "btageff_err" << " = " << btageff_err << endl ;
      if ( strcmp( label, "btageff_err") != 0 ) {  mismatchErr(label,"btageff_err") ; return false ; }



      // New for 3: ttwj and znn LDP/ZL MC ratios.

      float ttwj_ldp0lep_ratio[nBinsMET][nBinsHT][nBinsBtag] ;
      float  znn_ldp0lep_ratio[nBinsMET][nBinsHT][nBinsBtag] ;

      float ttwj_ldp0lep_ratio_err[nBinsMET][nBinsHT][nBinsBtag] ;
      float  znn_ldp0lep_ratio_err[nBinsMET][nBinsHT][nBinsBtag] ;

      printf("\n\n") ;
      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBtag ; k++) {
             char expectedlabel[1000] ;
             sprintf( expectedlabel, "ttwj_mc_ldpover0lep_ratio_M%d_H%d_%db", i+1, j+1, k+1 ) ;
             char vallabel[1000] ;
             char errlabel[1000] ;
             fscanf( infp, "%s %g", vallabel, &ttwj_ldp0lep_ratio[i][j][k] ) ;
             if ( strcmp( vallabel, expectedlabel ) != 0 ) { mismatchErr( label, expectedlabel ) ; return false ; }
             fscanf( infp, "%s %g", errlabel, &ttwj_ldp0lep_ratio_err[i][j][k] ) ;
             printf(" ttwj ldp/0lep ratio : %s  %6.3f +/- %5.3f\n", vallabel, ttwj_ldp0lep_ratio[i][j][k], ttwj_ldp0lep_ratio_err[i][j][k] ) ;
          } // k
        } // j
      } // i
      printf("\n\n") ;

      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          for (int k = 0 ; k < nBinsBtag ; k++) {
             char expectedlabel[1000] ;
             sprintf( expectedlabel, "znn_mc_ldpover0lep_ratio_M%d_H%d_%db", i+1, j+1, k+1 ) ;
             char vallabel[1000] ;
             char errlabel[1000] ;
             fscanf( infp, "%s %g", vallabel, &znn_ldp0lep_ratio[i][j][k] ) ;
             if ( strcmp( vallabel, expectedlabel ) != 0 ) { mismatchErr( label, expectedlabel ) ; return false ; }
             fscanf( infp, "%s %g", errlabel, &znn_ldp0lep_ratio_err[i][j][k] ) ;
             printf("  znn ldp/0lep ratio : %s  %6.3f +/- %5.3f\n", vallabel, znn_ldp0lep_ratio[i][j][k], znn_ldp0lep_ratio_err[i][j][k] ) ;
          } // k
        } // j
      } // i
      printf("\n\n") ;

      // NVVmc_0lep

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {
	    
	    TString inPar = "N_VVmc_0lep"+sMbins[i]+sHbins[j]+sBbins[k] ;
	    fscanf( infp, "%s %g", label, &NVVmc_0lep[i][j][k] ) ;
	    
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << NVVmc_0lep[i][j][k] << endl ;	  

	  }
	}
      }
      
      // NVVmc_1lep

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {
	    
	    TString inPar = "N_VVmc_1lep"+sMbins[i]+sHbins[j]+sBbins[k] ;
	    fscanf( infp, "%s %g", label, &NVVmc_1lep[i][j][k] ) ;
	    
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << NVVmc_1lep[i][j][k] << endl ;	  

	  }
	}
      }
      
      // NVVmc_ldp

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  for (int k = 0 ; k < nBinsBtag ; k++) {
	    
	    TString inPar = "N_VVmc_ldp"+sMbins[i]+sHbins[j]+sBbins[k] ;
	    fscanf( infp, "%s %g", label, &NVVmc_ldp[i][j][k] ) ;
	    
	    if ( label != inPar ) { mismatchErr(label,inPar) ; return false ; }
	    cout << inPar << " = " << NVVmc_ldp[i][j][k] << endl ;	  

	  }
	}
      }


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


      printf("\n Done reading in %s\n\n", infile ) ;
      fclose( infp ) ;
      
      
      //---- calculations for determining initial values for floating parameters.
      
      double initialval_znn_ee[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_znn_mm[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_znn[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_znn_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

      double initialval_qcd_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_qcd[nBinsMET][nBinsHT][nBinsBtag] ;

      double initialval_ttwj_sl[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_ttwj[nBinsMET][nBinsHT][nBinsBtag] ;
      double initialval_ttwj_ldp[nBinsMET][nBinsHT][nBinsBtag] ;


      //// double initialguess_ttwj_0lep1lep_ratio = 1.59 ;
      double initialguess_ttwj_0lep1lep_ratio = 1.28 ;

      printf("\n\n initial guess for ttwj 0lep/1lep ratio : %5.2f\n\n", initialguess_ttwj_0lep1lep_ratio ) ;

      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {

          if ( ignoreBin[i][j] ) continue ;

          for (int k = 0 ; k < nBinsBtag ; k++) {     

            // TTWJ stuff.

            initialval_ttwj_sl[i][j][k] = N_1lep[i][j][k] / trigeff_1L[i][j] ;
            initialval_ttwj[i][j][k] = sf_ttwj[i][j][k] * initialguess_ttwj_0lep1lep_ratio* initialval_ttwj_sl[i][j][k] ;
            initialval_ttwj_ldp[i][j][k] = ttwj_ldp0lep_ratio[i][j][k] * initialval_ttwj[i][j][k] ;

            // Z -> invis stuff :
	    
	    if ( k == 0 ) {
	      initialval_znn_ee[i][j][k] = (N_Zee[i][j]) * ( 5.95 * pur_Zee * knn_1b[i] ) / ( acc_Zee[i] * eff_Zee[j] ) ;
	      initialval_znn_mm[i][j][k] = (N_Zmm[i][j]) * ( 5.95 * pur_Zmm * knn_1b[i] ) / ( acc_Zmm[i] * eff_Zmm[j] ) ;
	    }
	    else if ( k == 1 ) {
	      initialval_znn_ee[i][j][k] = (N_Zee[i][j]) * ( 5.95 * pur_Zee * knn_2b ) / ( acc_Zee[i] * eff_Zee[j] ) ;
	      initialval_znn_mm[i][j][k] = (N_Zmm[i][j]) * ( 5.95 * pur_Zmm * knn_2b ) / ( acc_Zmm[i] * eff_Zmm[j] ) ;
	    }
	    else if ( k == 2 ) { 
	      initialval_znn_ee[i][j][k] = (N_Zee[i][j]) * ( 5.95 * pur_Zee * knn_3b ) / ( acc_Zee[i] * eff_Zee[j] ) ;
	      initialval_znn_mm[i][j][k] = (N_Zmm[i][j]) * ( 5.95 * pur_Zmm * knn_3b ) / ( acc_Zmm[i] * eff_Zmm[j] ) ;
	    }
	    else {
	      cout << "\n\n It looks like you are using more than 3 b-jet bins, you need to adjust the Z -> inv stuff! \n" << endl ;
	    }


            // simple average
            initialval_znn[i][j][k] = 0.5 * ( initialval_znn_ee[i][j][k] + initialval_znn_mm[i][j][k] ) ;
            initialval_znn_ldp[i][j][k] = znn_ldp0lep_ratio[i][j][k] * initialval_znn[i][j][k] ;


            // QCD stuff:

            //// initialval_qcd_ldp[i][j][k] = N_ldp[i][j][k] / trigeff_0L[i][j] - ( initialval_ttwj_ldp[i][j][k] + initialval_znn_ldp[i][j][k] ) ;

            initialval_qcd_ldp[i][j][k] = (1./trigeff_0L[i][j]) * ( N_ldp[i][j][k] - trigeff_1L[i][j] * ( initialval_ttwj_ldp[i][j][k] + initialval_znn_ldp[i][j][k] ) ) ;

          }
        }
      }


      printf("\n\n") ;

      double initialguess_model3_qcd_0lepLDP_ratio(0.) ;
      double initialguess_model24_qcd_0lepLDP_ratio[nBinsHT] ;
      double initialguess_model4_SFqcd_met[nBinsMET] ;
      double initialguess_model4_SFqcd_nb[nBinsBtag] ;

      //--- QCD is a bit more complicated for model 4

      if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {

         double tmp_qcd ;

         //-- Estimate the HT ratios from the first MET and nbjet bins
         for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
            ///// tmp_qcd = N_0lep[0][hbi][0] / trigeff_0L[0][hbi] - initialval_ttwj[0][hbi][0] - initialval_znn[0][hbi][0] ;
            tmp_qcd = (1./trigeff_0L[0][hbi])*( N_0lep[0][hbi][0] - trigeff_1L[0][hbi] * ( initialval_ttwj[0][hbi][0] + initialval_znn[0][hbi][0] )) ;
            if ( initialval_qcd_ldp[0][hbi][0] > 0. ) {
               initialguess_model24_qcd_0lepLDP_ratio[hbi] = tmp_qcd / ( sf_qcd[0][hbi][0] * initialval_qcd_ldp[0][hbi][0] ) ;
            } else {
               initialguess_model24_qcd_0lepLDP_ratio[hbi] = 0.0 ;
            }
            printf( "HT  bin %d : QCD 0lep/LDP ratio = (%7.1f / %7.1f) / %5.3f = %6.3f\n", hbi+1, tmp_qcd, initialval_qcd_ldp[0][hbi][0] , sf_qcd[0][hbi][0],
                 initialguess_model24_qcd_0lepLDP_ratio[hbi] ) ;
         } // hbi.

         //-- Get the nbjet scale factors.  Use lowest bins of MET and HT and compute for exact agreement.
         for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {
            if ( bbi == 0 ) {
               initialguess_model4_SFqcd_nb[0] = 1.0 ; //-- first one is 1 by definition.
            } else {
               ///// tmp_qcd = N_0lep[0][0][bbi] / trigeff_0L[0][0] - initialval_ttwj[0][0][bbi] - initialval_znn[0][0][bbi] ;
               tmp_qcd = (1./trigeff_0L[0][0])*( N_0lep[0][0][bbi] - trigeff_1L[0][0] * ( initialval_ttwj[0][0][bbi] + initialval_znn[0][0][bbi] )) ;
               double tmp_denom = sf_qcd[0][0][bbi] * initialguess_model24_qcd_0lepLDP_ratio[0] * initialval_qcd_ldp[0][0][bbi] ;
               if ( tmp_denom > 0. ) {
                  initialguess_model4_SFqcd_nb[bbi] = tmp_qcd / tmp_denom ;
               } else {
                  initialguess_model4_SFqcd_nb[bbi] = 0.0 ;
               }
               if ( bbi>1 && qcdModelIndex == 5 ) {
                  initialguess_model4_SFqcd_nb[bbi] = initialguess_model4_SFqcd_nb[1] ;
                  printf( "nb  bin %d : QCD SFnb = %6.3f (value from bin 2)\n", bbi+1, initialguess_model4_SFqcd_nb[bbi] ) ;
               } else {
                  printf( "nb  bin %d : QCD SFnb = (%7.1f / %7.1f) = %6.3f\n", bbi+1, tmp_qcd, tmp_denom,
                     initialguess_model4_SFqcd_nb[bbi] ) ;
               }
            }
         } // bbi.


         //-- Get the MET scale factors.  Use the lowest bins of HT and nbjets for exact agreement.
         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
               initialguess_model4_SFqcd_met[mbi] = 1.0 ;
            if ( mbi == 0 ) {
               initialguess_model4_SFqcd_met[0] = 1.0 ; //-- first one is 1 by definition.
            } else if ( mbi == 1 ) {
               ///// tmp_qcd = N_0lep[mbi][0][0] / trigeff_0L[mbi][0] - initialval_ttwj[mbi][0][0] - initialval_znn[mbi][0][0] ;
               tmp_qcd = (1./trigeff_0L[mbi][0])*( N_0lep[mbi][0][0] - trigeff_1L[mbi][0] * ( initialval_ttwj[mbi][0][0] + initialval_znn[mbi][0][0] )) ;
               double tmp_denom = sf_qcd[mbi][0][0] * initialguess_model24_qcd_0lepLDP_ratio[0] * initialval_qcd_ldp[mbi][0][0] ;
               if ( tmp_denom > 0. ) {
                  initialguess_model4_SFqcd_met[mbi] = tmp_qcd / tmp_denom ;
               } else {
                  initialguess_model4_SFqcd_met[mbi] = 0.0 ;
               }
               printf( "MET bin %d : QCD SFmet = (%7.1f / %7.1f) = %6.3f\n", mbi+1, tmp_qcd, tmp_denom,
                    initialguess_model4_SFqcd_met[mbi] ) ;
            } else if ( mbi == 2 ) {
               if ( qcdModelIndex == 4 ) {
                  initialguess_model4_SFqcd_met[mbi] = 1.41 ;
                  printf( "MET bin %d : QCD SFmet = %6.3f (hardwired guess)\n", mbi+1, initialguess_model4_SFqcd_met[mbi] ) ;
               } else if ( qcdModelIndex == 5 ) {
                  initialguess_model4_SFqcd_met[mbi] = initialguess_model4_SFqcd_met[1] ;
                  printf( "MET bin %d : QCD SFmet = %6.3f (value from bin 2)\n", mbi+1, initialguess_model4_SFqcd_met[mbi] ) ;
               }
            } else if ( mbi == 3 ) {
               if ( qcdModelIndex == 4 ) {
                  initialguess_model4_SFqcd_met[mbi] = 1.98 ;
                  printf( "MET bin %d : QCD SFmet = %6.3f (hardwired guess)\n", mbi+1, initialguess_model4_SFqcd_met[mbi] ) ;
               } else if ( qcdModelIndex == 5 ) {
                  initialguess_model4_SFqcd_met[mbi] = initialguess_model4_SFqcd_met[1] ;
                  printf( "MET bin %d : QCD SFmet = %6.3f (value from bin 2)\n", mbi+1, initialguess_model4_SFqcd_met[mbi] ) ;
               }
            }
         } // mbi.


         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
               for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {
                  initialval_qcd[mbi][hbi][bbi] = sf_qcd[mbi][hbi][bbi] * initialguess_model24_qcd_0lepLDP_ratio[hbi] * initialguess_model4_SFqcd_met[mbi] * initialguess_model4_SFqcd_nb[bbi] * initialval_qcd_ldp[mbi][hbi][bbi] ;
               } // bbi.
            } // hbi.
         } // mbi.


      } else if ( qcdModelIndex == 2 ) {

         double tmp_qcd ;

         //-- Estimate the HT ratios from the first MET and nbjet bins
         for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
            tmp_qcd = N_0lep[0][hbi][0] / trigeff_0L[0][hbi] - initialval_ttwj[0][hbi][0] - initialval_znn[0][hbi][0] ;
            if ( initialval_qcd_ldp[0][hbi][0] > 0. ) {
               initialguess_model24_qcd_0lepLDP_ratio[hbi] = tmp_qcd / ( sf_qcd[0][hbi][0] * initialval_qcd_ldp[0][hbi][0] ) ;
            } else {
               initialguess_model24_qcd_0lepLDP_ratio[hbi] = 0.0 ;
            }
            printf( "HT  bin %d : QCD 0lep/LDP ratio = (%7.1f / %7.1f) / %5.3f = %6.3f\n", hbi+1, tmp_qcd, initialval_qcd_ldp[0][hbi][0], sf_qcd[0][hbi][0],
                 initialguess_model24_qcd_0lepLDP_ratio[hbi] ) ;
         } // hbi.

         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
               for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {
                  initialval_qcd[mbi][hbi][bbi] = sf_qcd[mbi][hbi][bbi] * initialguess_model24_qcd_0lepLDP_ratio[hbi] * initialval_qcd_ldp[mbi][hbi][bbi] ;
               } // bbi.
            } // hbi.
         } // mbi.


      } else if ( qcdModelIndex == 3 ) {

         //-- Estimate the HT ratios from the first MET and nbjet bins
         double sum0lep(0.) ;
         double sumLDP(0.) ;
         for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
            double tmp_qcd = N_0lep[0][hbi][0] / trigeff_0L[0][hbi] - initialval_ttwj[0][hbi][0] - initialval_znn[0][hbi][0] ;
            sum0lep += tmp_qcd ;
            sumLDP += initialval_qcd_ldp[0][hbi][0] ;
         } // hbi.

         if ( sumLDP > 0. ) { initialguess_model3_qcd_0lepLDP_ratio = sum0lep / sumLDP ; }
         printf( "QCD 0lep/LDP ratio = %7.1f / %7.1f = %6.3f\n", sum0lep, sumLDP, initialguess_model3_qcd_0lepLDP_ratio ) ;

         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            for (int hbi = 0 ; hbi < nBinsHT ; hbi++) {
               for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {
                  initialval_qcd[mbi][hbi][bbi] = sf_qcd[mbi][hbi][bbi] * initialguess_model3_qcd_0lepLDP_ratio * initialval_qcd_ldp[mbi][hbi][bbi] ;
               } // bbi.
            } // hbi.
         } // mbi.

      }

      printf("\n\n") ;


      // print out initial values:


      printf("\n\n\n --------- Observables and floating parameter initial values. ------------\n\n") ;

      printf("           |   Nobs   |  Model  |    PDF        ||   ttwj   |   QCD   |   Znn   |\n") ;
      printf("-----------+----------+---------+---------------++----------+---------+---------+\n") ;

      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {

          if ( ignoreBin[i][j] ) continue ;

	  for (int k = 0 ; k < nBinsBtag ; k++) {     

            double model_0lep = trigeff_0L[i][j] * initialval_qcd[i][j][k] + trigeff_1L[i][j] * ( initialval_ttwj[i][j][k] + initialval_znn[i][j][k] ) ;
            double model_ldp  = trigeff_0L[i][j] * initialval_qcd_ldp[i][j][k] + trigeff_1L[i][j] * ( initialval_ttwj_ldp[i][j][k] + initialval_znn_ldp[i][j][k] ) ;
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
	    printf(" 0-lep     | %6d   | %7.1f | %8.6f %4s ||  %7.1f | %7.1f | %7.1f |\n", N_0lep[i][j][k], model_0lep, pdf_0lep, warning0lep,
                            trigeff_1L[i][j] * initialval_ttwj[i][j][k], trigeff_0L[i][j] * initialval_qcd[i][j][k], trigeff_1L[i][j] * initialval_znn[i][j][k] ) ;
	    printf(" ldp       | %6d   | %7.1f | %8.6f %4s ||  %7.1f | %7.1f | %7.1f |\n", N_ldp[i][j][k], model_ldp, pdf_ldp, warningldp,
                            trigeff_1L[i][j] * initialval_ttwj_ldp[i][j][k], trigeff_0L[i][j] * initialval_qcd_ldp[i][j][k], trigeff_1L[i][j] * initialval_znn_ldp[i][j][k] ) ;
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


      rv_mean_vv_sf = new RooRealVar( "mean_vv_sf", "mean_vv_sf", 0., 10. );			  
      rv_mean_vv_sf->setVal( 1.0 ) ;								  
      rv_mean_vv_sf->setConstant( kTRUE ) ;							  
      rv_width_vv_sf = new RooRealVar( "width_vv_sf", "width_vv_sf", 0., 10. );  		  
      rv_width_vv_sf->setVal( 1.00 ) ; 	   // arbitrarily set the uncertainty to 100%	  
      rv_width_vv_sf->setConstant( kTRUE ) ;							  

      RooRealVar* rv_ttwj_ldp0lep_ratio[nBinsMET][nBinsHT][nBinsBtag] ;
      RooRealVar*  rv_znn_ldp0lep_ratio[nBinsMET][nBinsHT][nBinsBtag] ;

            
      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {

          if ( ignoreBin[i][j] ) continue ;

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

            TString vv0lepString = "mu_vvmc";
            TString vv1lepString = "mu_vvmc_sl";
            TString vvldpString = "mu_vvmc_ldp";

	    TString MEffSfString     = "mean_eff_sf";
	    TString WEffSfString     = "width_eff_sf";
	    TString dEffdBtString    = "deff_dbtageff";

	    TString MEffSfSlString   = "mean_eff_sf_sl";
	    TString WEffSfSlString   = "width_eff_sf_sl";
	    TString dEffdBtSlString  = "deff_dbtageff_sl";

	    TString MEffSfLdpString  = "mean_eff_sf_ldp";
	    TString WEffSfLdpString  = "width_eff_sf_ldp";
	    TString dEffdBtLdpString = "deff_dbtageff_ldp";

            TString ttwjldp0lepString = "ttwj_ldp0lep_ratio" ;
            TString  znnldp0lepString =  "znn_ldp0lep_ratio" ;


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

            vv0lepString  += sMbins[i]+sHbins[j]+sBbins[k] ; 
            vv1lepString  += sMbins[i]+sHbins[j]+sBbins[k] ; 
            vvldpString   += sMbins[i]+sHbins[j]+sBbins[k] ; 

	    MEffSfString   += sMbins[i]+sHbins[j]+sBbins[k] ;
	    WEffSfString   += sMbins[i]+sHbins[j]+sBbins[k] ;
	    dEffdBtString  += sMbins[i]+sHbins[j]+sBbins[k] ;

            ttwjldp0lepString += sMbins[i]+sHbins[j]+sBbins[k] ;
            znnldp0lepString  += sMbins[i]+sHbins[j]+sBbins[k] ;


	    rv_mu_ttwj_sl[i][j][k] = new RooRealVar( muTtSlString, muTtSlString, 0., 100000. ) ;
	    rv_mu_ttwj_sl[i][j][k]->setVal( initialval_ttwj_sl[i][j][k] ) ;    // this is a starting value only

            //--- try allowing negative values.
	    ////// rrv_mu_qcd[i][j][k] = new RooRealVar( muQcdString, muQcdString, 0., 100000. ) ;
	    rrv_mu_qcd[i][j][k] = new RooRealVar( muQcdString, muQcdString, -1000., 100000. ) ;
	    rv_mu_qcd[i][j][k] = rrv_mu_qcd[i][j][k];
	    rrv_mu_qcd[i][j][k]->setVal( initialval_qcd[i][j][k] ) ;           // this is a starting value only

	    rrv_mu_qcd_ldp[i][j][k] = new RooRealVar( muQcdLdpString, muQcdLdpString, 0., 100000. ) ;
	    rv_mu_qcd_ldp[i][j][k] = rrv_mu_qcd_ldp[i][j][k];
	    rrv_mu_qcd_ldp[i][j][k]->setVal( initialval_qcd_ldp[i][j][k] ) ;   // this is a starting value only

            //-- owen: sept 23, 2012 : do not create these for >1 btag.
            if ( k==0 ) {
	       rrv_mu_znn[i][j][k] = new RooRealVar( muZnnString, muZnnString, 0., 100000. ) ;
	       rv_mu_znn[i][j][k] = rrv_mu_znn[i][j][k];
	       rrv_mu_znn[i][j][k]->setVal( initialval_znn[i][j][k] ) ;           // this is a starting value only
            }


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

	 // rv_mu_ttbarsingletopzjetsmc_ldp[i][j][k] = new RooRealVar( muTtLdpString, muTtLdpString, 0., 100000. ) ;
	 // rv_mu_ttbarsingletopzjetsmc_ldp[i][j][k]->setVal( Nttbarsingletopzjetsmc_ldp[i][j][k] ) ;
	 // rv_mu_ttbarsingletopzjetsmc_ldp[i][j][k]->setConstant( kTRUE ) ;

	 // rv_mu_WJmc_ldp[i][j][k] = new RooRealVar( muWjLdpString, muWjLdpString, 0., 100000. ) ;
	 // rv_mu_WJmc_ldp[i][j][k]->setVal( NWJmc_ldp[i][j][k] ) ;
	 // rv_mu_WJmc_ldp[i][j][k]->setConstant( kTRUE ) ;

	 // rv_mu_Znnmc_ldp[i][j][k] = new RooRealVar( muZnnLdpString, muZnnLdpString, 0., 100000. ) ;
	 // rv_mu_Znnmc_ldp[i][j][k]->setVal( NZnnmc_ldp[i][j][k] ) ;
	 // rv_mu_Znnmc_ldp[i][j][k]->setConstant( kTRUE ) ;

            rv_ttwj_ldp0lep_ratio[i][j][k] = new RooRealVar( ttwjldp0lepString, ttwjldp0lepString, 0., 100. ) ;
            rv_ttwj_ldp0lep_ratio[i][j][k] -> setVal( ttwj_ldp0lep_ratio[i][j][k] ) ;
            rv_ttwj_ldp0lep_ratio[i][j][k] -> setConstant( kTRUE ) ;

            rv_znn_ldp0lep_ratio[i][j][k] = new RooRealVar( znnldp0lepString, znnldp0lepString, 0., 100. ) ;
            rv_znn_ldp0lep_ratio[i][j][k] -> setVal( znn_ldp0lep_ratio[i][j][k] ) ;
            rv_znn_ldp0lep_ratio[i][j][k] -> setConstant( kTRUE ) ;


	    rv_mu_vvmc[i][j][k] = new RooRealVar( vv0lepString, vv0lepString, 0., 100. ) ;
	    rv_mu_vvmc[i][j][k]->setVal( NVVmc_0lep[i][j][k] ) ;
	    rv_mu_vvmc[i][j][k]->setConstant( kTRUE ) ;

	    rv_mu_vvmc_sl[i][j][k] = new RooRealVar( vv1lepString, vv1lepString, 0., 100. ) ;
	    rv_mu_vvmc_sl[i][j][k]->setVal( NVVmc_1lep[i][j][k] ) ;
	    rv_mu_vvmc_sl[i][j][k]->setConstant( kTRUE ) ;  

	    rv_mu_vvmc_ldp[i][j][k] = new RooRealVar( vvldpString, vvldpString, 0., 100. ) ;
	    rv_mu_vvmc_ldp[i][j][k]->setVal( NVVmc_ldp[i][j][k] ) ; 
	    rv_mu_vvmc_ldp[i][j][k]->setConstant( kTRUE ) ; 
	    
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


      //--- New in 3 : floating ttwj_0lep1lep_ratio

      RooRealVar* rv_ttwj_0lep1lep_ratio = new RooRealVar( "ttwj_0lep1lep_ratio", "ttwj_0lep1lep_ratio", 0., 5. ) ;
      rv_ttwj_0lep1lep_ratio -> setVal( initialguess_ttwj_0lep1lep_ratio ) ; // initial guess
      rv_ttwj_0lep1lep_ratio -> setConstant( kFALSE ) ;





      RooDataSet* dsObserved ;
      dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
				  *observedParametersList ) ;
      dsObserved->add( *observedParametersList ) ;

      printf(" --- Importing dataset.\n" ) ; cout << flush ;

      workspace.import(*dsObserved);



      char NP_name[1000] ;




      RooRealVar* rv_btageff_err = new RooRealVar( "btageff_err", "btageff_err", 0., 10. ) ;
      rv_btageff_err->setVal( btageff_err ) ;
      rv_btageff_err->setConstant( kTRUE ) ;


      //++++++++++

      bool ssspOk = setSusyScanPoint( inputScanFile,  m0,  m12,  isT1bbbb,  t1bbbbXsec, inputSusy_deff_dbtageff_file, systFile1 ) ;

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
      // THIS WILL NEED TO BE CAREFULLY REVISED !!!!
      // STILL USING THE LOG-NORMALS TO GET STARTED

      char formula[1024];



      //--- QCD model 1:  Use LSB pass/fail ratio for QCD 0lep/LDP ratio.
      RooRealVar* Rlsb_passfail_prim[nBinsHT][nBinsBtag];
      RooRealVar* Rlsb_passfail_nom[nBinsHT][nBinsBtag];
      RooGaussian* pdf_Rlsb_passfail[nBinsHT][nBinsBtag];
      RooFormulaVar* fv_Rlsb_passfail[nBinsHT][nBinsBtag];

      //--- QCD model 2,4:  Fit for QCD 0lep/LDP ratio in each HT bin.
      RooRealVar* rv_qcd_0lepLDP_ratio[20] ;

      //--- QCD model 4: like model 2 but includes MET and nbjet corrections.
      RooRealVar* rv_SFqcd_met[20] ;
      RooRealVar* rv_SFqcd_nb[20] ;

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
            rv_qcd_0lepLDP_ratio[htbi] = new RooRealVar( vname, vname, initialguess_model24_qcd_0lepLDP_ratio[htbi], 0., 10. ) ;
         }

      } else if ( qcdModelIndex == 3 ) {

         printf("\n\n  QCD Model 3 : One floating 0lep/LDP ratio.\n\n") ;

         char vname[1000] ;
         sprintf( vname, "qcd_0lepLDP_ratio" ) ;
         rv_qcd_0lepLDP_ratio[0] = new RooRealVar( vname, vname, initialguess_model3_qcd_0lepLDP_ratio, 0., 10. ) ;

      } else if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {

         printf("\n\n  QCD Model %d : One floating 0lep/LDP ratio for each HT bin. Includes MET and nbjet corrections.\n\n", qcdModelIndex) ;

         for ( int htbi=0; htbi < nBinsHT ; htbi++ ) {
            char vname[1000] ;
            sprintf( vname, "qcd_0lepLDP_ratio_H%d", htbi+1 ) ;
            printf("  HT bin %d : %s\n", htbi+1, vname ) ; cout << flush ;
            rv_qcd_0lepLDP_ratio[htbi] = new RooRealVar( vname, vname, initialguess_model24_qcd_0lepLDP_ratio[htbi], 0., 10. ) ;
         }

         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            if ( qcdModelIndex == 5 && mbi > 1 ) {
               rv_SFqcd_met[mbi] = rv_SFqcd_met[1] ;
            } else {
               char vname[1000] ;
               sprintf( vname, "SFqcd_met%d", mbi+1 ) ;
               printf("  MET bin %d : %s\n", mbi+1, vname ) ; cout << flush ;
               rv_SFqcd_met[mbi] = new RooRealVar( vname, vname, initialguess_model4_SFqcd_met[mbi], 0., 3. ) ;
               if ( mbi==0 ) {
                  rv_SFqcd_met[mbi] -> setConstant(kTRUE) ;
               } else {
                  rv_SFqcd_met[mbi] -> setConstant(kFALSE) ;
               }
            }
         }

         for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {
            if ( qcdModelIndex == 5 && bbi > 1 ) {
               rv_SFqcd_nb[bbi] = rv_SFqcd_nb[1] ;
            } else {
               char vname[1000] ;
               sprintf( vname, "SFqcd_nb%d", bbi+1 ) ;
               printf("  nbjet bin %d : %s\n", bbi+1, vname ) ; cout << flush ;
               rv_SFqcd_nb[bbi] = new RooRealVar( vname, vname, initialguess_model4_SFqcd_nb[bbi], 0., 3. ) ;
               if ( bbi==0 ) {
                  rv_SFqcd_nb[bbi] -> setConstant(kTRUE) ;
               } else {
                  rv_SFqcd_nb[bbi] -> setConstant(kFALSE) ;
               }
            }
         }

      } else {

         printf("\n\n *** Unknown QCD model index : %d\n\n", qcdModelIndex ) ;
         return false ;

      }





      printf("\n\n Z -> ll acceptance.\n\n") ; cout << flush ;

      // Z -> ll acceptance

      RooAbsReal* rar_acc_Zee[nBinsMET] ;
      RooAbsReal* rar_acc_Zmm[nBinsMET] ;

      for (int i = 0 ; i < nBinsMET ; i++) {

         sprintf( NP_name, "acc_Zee_M%d", i+1 ) ;
         rar_acc_Zee[i] = makeBetaConstraint( NP_name, acc_Zee[i], acc_Zee_err[i] ) ;

         sprintf( NP_name, "acc_Zmm_M%d", i+1 ) ;
         rar_acc_Zmm[i] = makeBetaConstraint( NP_name, acc_Zmm[i], acc_Zmm_err[i] ) ;

      }











      printf("\n\n Z -> ll efficiencies.\n\n") ; cout << flush ;

      // Z -> ll efficiencies

      RooAbsReal* rar_eff_Zee[nBinsHT] ;
      RooAbsReal* rar_eff_Zmm[nBinsHT] ;

      for (int j = 0 ; j < nBinsHT ; j++) {

         sprintf( NP_name, "eff_Zee_H%d", j+1 ) ;
         rar_eff_Zee[j] = makeBetaConstraint( NP_name, eff_Zee[j], eff_Zee_err[j] ) ;

         sprintf( NP_name, "eff_Zmm_H%d", j+1 ) ;
         rar_eff_Zmm[j] = makeBetaConstraint( NP_name, eff_Zmm[j], eff_Zmm_err[j] ) ;

      }




      printf("\n\n Z -> ll purities.\n\n") ; cout << flush ;

      // Z -> ll purities

       sprintf( NP_name, "pur_Zee" ) ;
       RooAbsReal* rar_pur_Zee = makeBetaConstraint( NP_name, pur_Zee, pur_Zee_err ) ;

       sprintf( NP_name, "pur_Zmm" ) ;
       RooAbsReal* rar_pur_Zmm = makeBetaConstraint( NP_name, pur_Zmm, pur_Zmm_err ) ;









      printf("\n\n Z -> ll btag scale factors.\n\n") ; cout << flush ;

      // VL -> nominal scale factors

      RooAbsReal* rar_knn_1b[nBinsMET] ;
      RooAbsReal* rar_knn_2b ;
      RooAbsReal* rar_knn_3b ;

      for (int i = 0 ; i < nBinsMET ; i++) {

         sprintf( NP_name, "knn_1b_M%d", i+1 ) ;
         rar_knn_1b[i] = makeGaussianConstraint( NP_name, knn_1b[i], knn_1b_err[i] ) ;

      }

      sprintf( NP_name, "knn_2b" ) ;
      rar_knn_2b = makeGaussianConstraint( NP_name, knn_2b, knn_2b_err ) ;

      sprintf( NP_name, "knn_3b" ) ;
      rar_knn_3b = makeGaussianConstraint( NP_name, knn_3b, knn_3b_err ) ;




      printf("\n\n Z -> ll closure systematics.\n\n") ; cout << flush ;

      // sf_ee and sf_mm derived from a common underlying parameter

      RooAbsReal* rar_sf_ee ;
      RooAbsReal* rar_sf_mm ;

      char sfllbpname[1000] ;
      sprintf( sfllbpname, "sf_ll" ) ;

      sprintf( NP_name, "sf_ee" ) ;
      rar_sf_ee = makeCorrelatedGaussianConstraint( NP_name, sf_ee, sf_ee_err, sfllbpname ) ;
      
      sprintf( NP_name, "sf_mm" ) ;
      rar_sf_mm = makeCorrelatedGaussianConstraint( NP_name, sf_mm, sf_mm_err, sfllbpname ) ;









      printf("\n\n MC scale factor.\n\n") ; cout << flush ;

      // MC scale factor

      sprintf( NP_name, "sf_mc" ) ;
      //// RooAbsReal* rar_sf_mc = makeBetaPrimeConstraint( NP_name, sf_mc, sf_mc_err ) ;
      RooAbsReal* rar_sf_mc = makeGaussianConstraint( NP_name, sf_mc, sf_mc_err ) ;







      printf("\n\n QCD and ttwj scale factors.\n\n") ; cout << flush ;

      // QCD and TTWJ scale factors

      RooAbsReal* rar_sf_qcd [nBinsMET][nBinsHT][nBinsBtag];
      RooAbsReal* rar_sf_ttwj[nBinsMET][nBinsHT][nBinsBtag];

      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {
          if ( ignoreBin[i][j] ) continue ;
          for (int k = 0 ; k < nBinsBtag ; k++) {


             //--- QCD

             sprintf( NP_name, "sf_qcd_M%d_H%d_%db", i+1, j+1, k+1 ) ;
             ///// rar_sf_qcd[i][j][k] = makeBetaPrimeConstraint( NP_name, sf_qcd[i][j][k], sf_qcd_err[i][j][k] ) ;
             rar_sf_qcd[i][j][k] = makeGaussianConstraint( NP_name, sf_qcd[i][j][k], sf_qcd_err[i][j][k] ) ;




             //--- ttwj


             sprintf( NP_name, "sf_ttwj_M%d_H%d_%db", i+1, j+1, k+1 ) ;
             ///// rar_sf_ttwj[i][j][k] = makeBetaPrimeConstraint( NP_name, sf_ttwj[i][j][k], sf_ttwj_err[i][j][k] ) ;
             rar_sf_ttwj[i][j][k] = makeGaussianConstraint( NP_name, sf_ttwj[i][j][k], sf_ttwj_err[i][j][k] ) ;


             printf("--------\n") ;

          } // i (met)
             printf("++++++++++++++++\n") ;
        } // j (ht)
             printf("=========================\n") ;
      } // k (nbtag)









      //-- If using QCD model 4, give a few floating parameters some help with a constraint PDF.

      printf("\n\n") ;
      if ( qcdModelIndex == 4 ) {

        //-- constrain MET scale factors for >= bin 3 (counting from 1).
         if ( nBinsMET==4 ) {

            char pdfname[1000] ;
            char meanname[1000] ;
            char signame[1000] ;

            sprintf( pdfname, "pdf_%s", rv_SFqcd_met[2]->GetName() ) ;
            sprintf( meanname, "pdf_mean_%s", rv_SFqcd_met[2]->GetName() ) ;
            sprintf( signame, "pdf_sigma_%s", rv_SFqcd_met[2]->GetName() ) ;
            RooConstVar* mean_mb3 = new RooConstVar( meanname, meanname, 1.34 ) ; //-- value hardwired from chi2 fit of MC ratios.
            RooConstVar* sigma_mb3 = new RooConstVar( signame, signame, 0.14 ) ; //-- value is diff between SFmet3 and SFmet2 divided by 2.
            //RooConstVar* sigma_mb3 = new RooConstVar( signame, signame, 0.34 ) ; //-- value is diff between 1-SF
	    printf(" QCD model 4 : adding constraint PDF for %s with mean %5.3f and sigma %5.3f\n", rv_SFqcd_met[2]->GetName(), mean_mb3->getVal(), sigma_mb3->getVal() ) ;
            RooGaussian* rg_mb3 = new RooGaussian( pdfname, pdfname, *rv_SFqcd_met[2], *mean_mb3, *sigma_mb3) ;
            allNuisances -> add( *rv_SFqcd_met[2] ) ;
            allNuisancePdfs -> add( *rg_mb3 ) ;

            sprintf( pdfname, "pdf_%s", rv_SFqcd_met[3]->GetName() ) ;
            sprintf( meanname, "pdf_mean_%s", rv_SFqcd_met[3]->GetName() ) ;
            sprintf( signame, "pdf_sigma_%s", rv_SFqcd_met[3]->GetName() ) ;
            RooConstVar* mean_mb4 = new RooConstVar( meanname, meanname, 1.90 ) ; //-- value hardwired from chi2 fit of MC ratios.
            RooConstVar* sigma_mb4 = new RooConstVar( signame, signame, 0.43 ) ; //-- value is diff between SFmet4 and SFmet2 divided by 2.
            //RooConstVar* sigma_mb4 = new RooConstVar( signame, signame, 0.90 ) ; //-- value is 1-SF
	    printf(" QCD model 4 : adding constraint PDF for %s with mean %5.3f and sigma %5.3f\n", rv_SFqcd_met[3]->GetName(), mean_mb4->getVal(), sigma_mb4->getVal() ) ;
            RooGaussian* rg_mb4 = new RooGaussian( pdfname, pdfname, *rv_SFqcd_met[3], *mean_mb4, *sigma_mb4) ;
            allNuisances -> add( *rv_SFqcd_met[3] ) ;
            allNuisancePdfs -> add( *rg_mb4 ) ;

            sprintf( pdfname, "pdf_%s", rv_SFqcd_nb[2]->GetName() ) ;
            sprintf( meanname, "pdf_mean_%s", rv_SFqcd_nb[2]->GetName() ) ;
            sprintf( signame, "pdf_sigma_%s", rv_SFqcd_nb[2]->GetName() ) ;
            RooConstVar* mean_nb3 = new RooConstVar( meanname, meanname, 0.69 ) ; //-- value hardwired from chi2 fit of MC ratios.
            RooConstVar* sigma_nb3 = new RooConstVar( signame, signame, 0.16 ) ; //-- guess: 1-SF / 2.  chi2 fit error is 0.07 so this is reasonable in magnitude.
            //RooConstVar* sigma_nb3 = new RooConstVar( signame, signame, 0.31 ) ; //-- guess: 1-SF 
	    printf(" QCD model 4 : adding constraint PDF for %s with mean %5.3f and sigma %5.3f\n", rv_SFqcd_nb[2]->GetName(), mean_nb3->getVal(), sigma_nb3->getVal() ) ;
            RooGaussian* rg_nb3 = new RooGaussian( pdfname, pdfname, *rv_SFqcd_nb[2], *mean_nb3, *sigma_nb3) ;
            allNuisances -> add( *rv_SFqcd_nb[2] ) ;
            allNuisancePdfs -> add( *rg_nb3 ) ;

         } else {
            printf("\n\n *** QCD model 4 : NOT adding constraint PDFs for MET SFs, since I don't know what to do for this binning.\n") ;
         }
      }
      printf("\n\n") ;










      //-- Using Gaussian here since this is a scale factor, not the actual
      //   efficiency, and the efficiency is far from 1 (much closer to zero).
      //

      RooAbsReal* rar_eff_sf[nBinsMET][nBinsHT][nBinsBtag] ;
      RooAbsReal* rar_eff_sf_sl[nBinsMET][nBinsHT][nBinsBtag] ;
      RooAbsReal* rar_eff_sf_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

      char effbpname[1000] ;
      sprintf( effbpname, "eff_sf" ) ;

      for (int i = 0 ; i < nBinsMET ; i++) {
         for (int j = 0 ; j < nBinsHT ; j++) {
            if ( ignoreBin[i][j] ) continue ;
            for (int k = 0 ; k < nBinsBtag ; k++) {

               sprintf( NP_name, "eff_sf_M%d_H%d_%db", i+1, j+1, k+1 ) ;
               if ( !blindStudy ) {
                  rar_eff_sf[i][j][k] = makeGaussianConstraint( NP_name, rv_mean_eff_sf[i][j][k]->getVal(), rv_width_eff_sf[i][j][k]->getVal() ) ;
               } else {
                  rar_eff_sf[i][j][k] = makeGaussianConstraint( NP_name, rv_mean_eff_sf[i][j][k]->getVal(), 0.0 ) ;
               }

               sprintf( NP_name, "eff_sf_sl_M%d_H%d_%db", i+1, j+1, k+1 ) ;
               if ( !blindStudy ) {
                  rar_eff_sf_sl[i][j][k] = makeGaussianConstraint( NP_name, rv_mean_eff_sf_sl[i][j][k]->getVal(), rv_width_eff_sf_sl[i][j][k]->getVal() ) ;
               } else {
                  rar_eff_sf_sl[i][j][k] = makeGaussianConstraint( NP_name, rv_mean_eff_sf_sl[i][j][k]->getVal(), 0.0 ) ;
               }

               sprintf( NP_name, "eff_sf_ldp_M%d_H%d_%db", i+1, j+1, k+1 ) ;
               if ( !blindStudy ) {
                  rar_eff_sf_ldp[i][j][k] = makeGaussianConstraint( NP_name, rv_mean_eff_sf_ldp[i][j][k]->getVal(), rv_width_eff_sf_ldp[i][j][k]->getVal() ) ;
               } else {
                  rar_eff_sf_ldp[i][j][k] = makeGaussianConstraint( NP_name, rv_mean_eff_sf_ldp[i][j][k]->getVal(), 0.0 ) ;
               }

            } // k
         } // j
      } // i







      //-- Using Gaussian here since this is a scale factor, not the actual
      //   efficiency, and the efficiency is far from 1 (much closer to zero).
      //

      RooAbsReal* rar_btageff_sf[nBinsMET][nBinsHT][nBinsBtag] ;
      RooAbsReal* rar_btageff_sf_sl[nBinsMET][nBinsHT][nBinsBtag] ;
      RooAbsReal* rar_btageff_sf_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

      char btageffbpname[1000] ;
      sprintf( btageffbpname, "btageff_sf" ) ;

      for (int i = 0 ; i < nBinsMET ; i++) {
         for (int j = 0 ; j < nBinsHT ; j++) {
            if ( ignoreBin[i][j] ) continue ;
            for (int k = 0 ; k < nBinsBtag ; k++) {

               bool changeSign ;

               sprintf( NP_name, "btageff_sf_M%d_H%d_%db", i+1, j+1, k+1 ) ;
               changeSign = false ;
               if ( rv_deff_dbtageff[i][j][k]->getVal() < 0. ) { changeSign = true ; }
               if ( !blindStudy ) {
                  rar_btageff_sf[i][j][k] = makeCorrelatedGaussianConstraint( NP_name, 1.0, fabs(rv_deff_dbtageff[i][j][k]->getVal()), btageffbpname, changeSign ) ;
               } else {
                  rar_btageff_sf[i][j][k] = makeCorrelatedGaussianConstraint( NP_name, 1.0, 0.0, btageffbpname, changeSign ) ;
               }

               sprintf( NP_name, "btageff_sf_sl_M%d_H%d_%db", i+1, j+1, k+1 ) ;
               changeSign = false ;
               if ( rv_deff_dbtageff_sl[i][j][k]->getVal() < 0. ) { changeSign = true ; }
               if ( !blindStudy ) {
                  rar_btageff_sf_sl[i][j][k] = makeCorrelatedGaussianConstraint( NP_name, 1.0, fabs(rv_deff_dbtageff_sl[i][j][k]->getVal()), btageffbpname, changeSign ) ;
               } else {
                  rar_btageff_sf_sl[i][j][k] = makeCorrelatedGaussianConstraint( NP_name, 1.0, 0.0, btageffbpname, changeSign ) ;
               }

               sprintf( NP_name, "btageff_sf_ldp_M%d_H%d_%db", i+1, j+1, k+1 ) ;
               changeSign = false ;
               if ( rv_deff_dbtageff_ldp[i][j][k]->getVal() < 0. ) { changeSign = true ; }
               if ( !blindStudy ) {
                  rar_btageff_sf_ldp[i][j][k] = makeCorrelatedGaussianConstraint( NP_name, 1.0, fabs(rv_deff_dbtageff_ldp[i][j][k]->getVal()), btageffbpname, changeSign ) ;
               } else {
                  rar_btageff_sf_ldp[i][j][k] = makeCorrelatedGaussianConstraint( NP_name, 1.0, 0.0, btageffbpname, changeSign ) ;
               }

            } // k
         } // j
      } // i

      // set the Gaussian constraint on the diboson background scale factor
      RooAbsReal* rar_vv_sf = makeGaussianConstraint( "rar_vv_sf", rv_mean_vv_sf->getVal(), rv_width_vv_sf->getVal() ) ;

      // make beta constraints for trigger efficiency corrections
      printf("\n\n trigger efficiencies.\n\n") ; cout << flush ;

      RooAbsReal* rar_trigeff[nBinsMET][nBinsHT];
      RooAbsReal* rar_trigeff_sl[nBinsMET][nBinsHT];

      for (int i = 0 ; i < nBinsMET ; i++) {
        for (int j = 0 ; j < nBinsHT ; j++) {

          if ( ignoreBin[i][j] ) continue ;

          sprintf( NP_name, "trigeff_M%d_H%d", i+1, j+1 ) ;
          rar_trigeff[i][j] = makeBetaConstraint( NP_name, trigeff_0L[i][j], trigefferr_0L[i][j] ) ;

          sprintf( NP_name, "trigeff_sl_M%d_H%d", i+1, j+1 ) ;
          rar_trigeff_sl[i][j] = makeBetaConstraint( NP_name, trigeff_1L[i][j], trigefferr_1L[i][j] ) ;
        }
      }



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

          if ( ignoreBin[i][j] ) continue ;

	  for (int k = 0 ; k < nBinsBtag ; k++) {     





             printf(" met,ht,nb : %d, %d, %d : ", i+1, j+1, k+1 ) ; cout << flush ;

             printf(" ttwj,") ; cout << flush ;

	    //---- TTWJ

	    TString ttwjString   = "mu_ttwj";
	    ttwjString   += sMbins[i]+sHbins[j]+sBbins[k] ;


            TString ttwjrfvString = " @0 * @1 * @2" ;
            rfv_mu_ttwj[i][j][k] = new RooFormulaVar( ttwjString, ttwjrfvString,
                            RooArgSet( *rar_sf_ttwj[i][j][k], *rv_mu_ttwj_sl[i][j][k], *rv_ttwj_0lep1lep_ratio ) ) ;

            rv_mu_ttwj[i][j][k] = rfv_mu_ttwj[i][j][k] ;






            TString ttwjLdpString  = "mu_ttwj_ldp" ;
            ttwjLdpString  += sMbins[i]+sHbins[j]+sBbins[k] ;

            TString rfvString = "@0 * @1" ;

            rv_mu_ttwj_ldp[i][j][k] = new RooFormulaVar( ttwjLdpString, rfvString,
                                                         RooArgSet( *rv_mu_ttwj[i][j][k], *rv_ttwj_ldp0lep_ratio[i][j][k] )) ;











             printf(" qcd,") ; cout << flush ;

           //---- QCD

            TString qcdString    = "mu_qcd" ;
            qcdString    += sMbins[i]+sHbins[j]+sBbins[k] ;


            if ( qcdModelIndex == 1 ) {

               TString rfvQcdString = "@0 * @1 * @2" ;

               rfv_mu_qcd[i][j][k] = new RooFormulaVar( qcdString, rfvQcdString, 
                                                        RooArgSet( *rv_mu_qcd_ldp[i][j][k], *rar_sf_qcd[i][j][k], *fv_Rlsb_passfail[j][k] ) ) ;


            } else if ( qcdModelIndex == 2 ) {

               TString rfvQcdString = "@0 * @1 * @2" ;

               rfv_mu_qcd[i][j][k] = new RooFormulaVar( qcdString, rfvQcdString, 
                                                        RooArgSet( *rv_mu_qcd_ldp[i][j][k], *rar_sf_qcd[i][j][k], *rv_qcd_0lepLDP_ratio[j] ) ) ;

            } else if ( qcdModelIndex == 3 ) {

               TString rfvQcdString = "@0 * @1 * @2" ;

               rfv_mu_qcd[i][j][k] = new RooFormulaVar( qcdString, rfvQcdString, 
                                                        RooArgSet( *rv_mu_qcd_ldp[i][j][k], *rar_sf_qcd[i][j][k], *rv_qcd_0lepLDP_ratio[0] ) ) ;

            } else if ( qcdModelIndex == 4 || qcdModelIndex == 5 ) {

               TString rfvQcdString = "@0 * @1 * @2 * @3 * @4" ;

               rfv_mu_qcd[i][j][k] = new RooFormulaVar( qcdString, rfvQcdString, 
                                                        RooArgSet( *rv_mu_qcd_ldp[i][j][k], *rar_sf_qcd[i][j][k], *rv_qcd_0lepLDP_ratio[j], *rv_SFqcd_met[i], *rv_SFqcd_nb[k] ) ) ;

            }

            rv_mu_qcd[i][j][k] = rfv_mu_qcd[i][j][k] ;


             printf(" susy,") ; cout << flush ;

	    //---- SUSY

	    TString susyString    = "mu_susy" ;
	    TString susySlString  = "mu_susy_sl" ;
	    TString susyLdpString = "mu_susy_ldp" ;
	    
	    susyString    += sMbins[i]+sHbins[j]+sBbins[k] ;
	    susySlString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    susyLdpString += sMbins[i]+sHbins[j]+sBbins[k] ;








            TString rfvSusyString = "@0 * ( @1 / @2)" ;

            rv_mu_susy[i][j][k] = new RooFormulaVar( susyString, rfvSusyString,
                                                       RooArgSet( *rv_mu_susymc[i][j][k], *rv_mu_susy_all0lep, *rv_mu_susymc_all0lep ) ) ;

            TString rfvSusySlString  = "@0 * ( @1 / @2)" ;

            rv_mu_susy_sl[i][j][k] = new RooFormulaVar( susySlString, rfvSusySlString,
                                                        RooArgSet( *rv_mu_susymc_sl[i][j][k], *rv_mu_susy_all0lep, *rv_mu_susymc_all0lep ) ) ;

            TString rfvSusyLdpString = "@0 * ( @1 / @2 )" ;

            rv_mu_susy_ldp[i][j][k] = new RooFormulaVar( susyLdpString, rfvSusyLdpString,
                                                         RooArgSet( *rv_mu_susymc_ldp[i][j][k], *rv_mu_susy_all0lep, *rv_mu_susymc_all0lep ) ) ;





             printf(" znn\n") ; cout << flush ;


            //---- Z -> nunu

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
						   RooArgSet( *rv_mu_znn[i][j][k], *rar_knn_1b[i], *rar_sf_ee, *rar_acc_Zee[i], 
							      *rar_eff_Zee[j], *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;


	      TString rfvZmmString = "( @0 / @1 ) * @2 * ( ( @3 * @4 ) / ( @5 * @6 ) )" ;

	      rv_mu_zmm[i][j] = new RooFormulaVar( muZmmString, rfvZmmString,
						   RooArgSet( *rv_mu_znn[i][j][k], *rar_knn_1b[i], *rar_sf_mm, *rar_acc_Zmm[i], 
							      *rar_eff_Zmm[j], *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

	    }
	    else if ( k == 1 ) {
	      
	      TString rfvZnnString = "@0 * ( @1 / @2 )" ;

	      rv_mu_znn[i][j][k] = new RooFormulaVar( muZnnString, rfvZnnString,
						      RooArgSet( *rv_mu_znn[i][j][0], *rar_knn_2b, *rar_knn_1b[i] ) ) ;
	      
	    }
	    else if ( k == 2 ) {
	      
	      TString rfvZnnString = "@0 * ( @1 / @2 )" ;

	      rv_mu_znn[i][j][k] = new RooFormulaVar( muZnnString, rfvZnnString,
						      RooArgSet( *rv_mu_znn[i][j][0], *rar_knn_3b, *rar_knn_1b[i] ) ) ;
	      
	    }
	    else {
	      cout << "\n\n It looks like you are using more than 3 b-jet bins, you need to adjust the Z -> inv stuff! \n" << endl ;
	    }
      





            TString znnLdpString   = "mu_znn_ldp" ;
            znnLdpString   += sMbins[i]+sHbins[j]+sBbins[k] ;

            TString znnMcLdpString = "@0 * @1" ;

            rv_mu_znn_ldp[i][j][k] = new RooFormulaVar( znnLdpString, znnMcLdpString,
                                                         RooArgSet( *rv_mu_znn[i][j][k], *rv_znn_ldp0lep_ratio[i][j][k] )) ;

	    TString rfvVVString = "@0 * @1" ;

            TString vvString0lep = "mu_vv";
            vvString0lep += sMbins[i]+sHbins[j]+sBbins[k] ;
	    rv_mu_vv[i][j][k] = new RooFormulaVar( vvString0lep, rfvVVString, RooArgSet( *rar_vv_sf, *rv_mu_vvmc[i][j][k] ));

            TString vvString1lep = "mu_vv";
            vvString1lep += sMbins[i]+sHbins[j]+sBbins[k] ;
	    rv_mu_vv_sl[i][j][k] = new RooFormulaVar( vvString1lep, rfvVVString, RooArgSet( *rar_vv_sf, *rv_mu_vvmc_sl[i][j][k] ));

            TString vvStringldp = "mu_vv";
            vvStringldp += sMbins[i]+sHbins[j]+sBbins[k] ;
	    rv_mu_vv_ldp[i][j][k] = new RooFormulaVar( vvStringldp, rfvVVString, RooArgSet( *rar_vv_sf, *rv_mu_vvmc_ldp[i][j][k] ));




            //+++++++++++++ Expected counts for observables in terms of parameters ++++++++++++++++++

	    TString nString    = "n" ;
	    TString nSlString  = "n_sl" ;
	    TString nLdpString = "n_ldp" ;

	    nString    += sMbins[i]+sHbins[j]+sBbins[k] ;
	    nSlString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	    nLdpString += sMbins[i]+sHbins[j]+sBbins[k] ;

            //// TString rfvNString =  "(@0 + @1 + @2 + @3 + (@4 * @5 * @6)) * @7" ;
            TString rfvNString =  "@0 * @1 + (@2 + @3 + @4 + (@5 * @6 * @7)) * @8" ;

            rv_n[i][j][k] = new RooFormulaVar( nString, rfvNString,
                                               RooArgSet( *rv_mu_qcd[i][j][k], *rar_trigeff[i][j],
                                                          *rv_mu_ttwj[i][j][k], *rv_mu_znn[i][j][k],  *rv_mu_vv[i][j][k],
                                                          *rar_btageff_sf[i][j][k], *rar_eff_sf[i][j][k], *rv_mu_susy[i][j][k], *rar_trigeff_sl[i][j] ) ) ;



            TString rfvNSlString = "(@0 + @1 + (@2 * @3 * @4)) * @5" ;

            rv_n_sl[i][j][k] = new RooFormulaVar( nSlString, rfvNSlString,
                                                  RooArgSet( *rv_mu_ttwj_sl[i][j][k], *rv_mu_vv_sl[i][j][k], *rar_btageff_sf_sl[i][j][k], 
                                                             *rar_eff_sf_sl[i][j][k], *rv_mu_susy_sl[i][j][k], *rar_trigeff_sl[i][j] ) ) ;

      //////--------------------------------
      //
      ////  TString rfvNLdpString = "@0 * @1 + @2 * ( @3 + @4 + @5 + (@6 * @7 * @8 ) ) * @9" ;
      //
      //
      ////  rv_n_ldp[i][j][k] = new RooFormulaVar( nLdpString, rfvNLdpString,
      ////                                         RooArgSet( *rar_trigeff[i][j], *rv_mu_qcd_ldp[i][j][k],
      ////                                                    *rar_sf_mc, *rv_mu_vv_ldp[i][j][k], *rv_mu_ttwj_ldp[i][j][k], *rv_mu_znn_ldp[i][j][k],
      ////                                                    *rar_btageff_sf_ldp[i][j][k], *rar_eff_sf_ldp[i][j][k], *rv_mu_susy_ldp[i][j][k],
      ////                                                    *rar_trigeff_sl[i][j] ) ) ;
      //
      ////  ///---  Above won't work because damn RooArgSet constructor only takes a maximum of 9 terms.  This stupid equation needs 10...
      //
      //////--------------------------------

            char mu_nonqcdsm_name[1000] ;
            sprintf( mu_nonqcdsm_name, "mu_nonqcdsm_ldp_M%d_H%d_%db", i+1, j+1, k+1 ) ;
            RooFormulaVar* mu_nonqcdsm_ldp = new RooFormulaVar( mu_nonqcdsm_name, "@0 + @1 + @2",
                                                                RooArgSet( *rv_mu_vv_ldp[i][j][k], *rv_mu_ttwj_ldp[i][j][k], *rv_mu_znn_ldp[i][j][k] ) ) ;

            TString rfvNLdpString = "@0 * @1 + @2 * ( @3 + (@4 * @5 * @6 ) ) * @7" ;

            rv_n_ldp[i][j][k] = new RooFormulaVar( nLdpString, rfvNLdpString,
                                                   RooArgSet( *rar_trigeff[i][j], *rv_mu_qcd_ldp[i][j][k],
                                                              *rar_sf_mc, *mu_nonqcdsm_ldp,
                                                              *rar_btageff_sf_ldp[i][j][k], *rar_eff_sf_ldp[i][j][k], *rv_mu_susy_ldp[i][j][k],
                                                              *rar_trigeff_sl[i][j] ) ) ;

      //////--------------------------------


	    // pdf's

	    TString pdfN0lepString = "pdf_N_0lep";
	    TString pdfN1lepString = "pdf_N_1lep";
	    TString pdfNldpString  = "pdf_N_ldp";

	    pdfN0lepString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    pdfN1lepString += sMbins[i]+sHbins[j]+sBbins[k] ;
	    pdfNldpString  += sMbins[i]+sHbins[j]+sBbins[k] ;

            if ( blind0lepBin[i][j][k] ) {
               printf(" * Not including 0lep met=%d, ht=%d, nb=%d PDF in the likelihood.\n", i+1, j+1, k+1 ) ;
               pdf_N_0lep[i][j][k] = new RooConstVar( pdfN0lepString, pdfN0lepString, 1.0 ) ;
            } else {
	       pdf_N_0lep[i][j][k] = new RooPoisson( pdfN0lepString, pdfN0lepString, *rv_0lep[i][j][k], *rv_n[i][j][k] ) ;
	       pdflist.add( *pdf_N_0lep[i][j][k] ) ;
            }

	    pdf_N_1lep[i][j][k] = new RooPoisson( pdfN1lepString, pdfN1lepString, *rv_1lep[i][j][k], *rv_n_sl[i][j][k] ) ;
	    pdf_N_ldp[i][j][k]  = new RooPoisson( pdfNldpString,  pdfNldpString,  *rv_ldp[i][j][k],  *rv_n_ldp[i][j][k] ) ;

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

         double scanLow  = 1.0 ;
         double scanHigh = 2.0 ;
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
						   double m0, double m12, bool isT1bbbb, double t1bbbbXsec,
						   const char* inputSusy_deff_dbtageff_file,
						   const char* systFile1
						   ) {

      double dummy ;
      dummy = t1bbbbXsec ;

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


      //--- Loop over the scan points.
      while ( infp.good() ) {

	if ( found ) break ;

	int ArraySize = 3 + 4*(nBinsMET*nBinsHT*nBinsBtag) ;
	double ArrayContent[ArraySize] ;

	for (int i = 0; infp && i < ArraySize; ++ i) {
	  infp >> ArrayContent[i];
	}

	pointM0  = ArrayContent[0] ;
	pointM12 = ArrayContent[1] ;

	nGen = (int)ArrayContent[2] ;

	if (    fabs( pointM0 - m0 ) <= deltaM0/2 && fabs( pointM12 - m12 ) <= deltaM12/2 ) {

	  found = true ;

	  for (int i = 0 ; i < nBinsMET ; i++) {
	    for (int j = 0 ; j < nBinsHT ; j++) {
	      for (int k = 0 ; k < nBinsBtag ; k++) {     
		
		n_0l_raw[i][j][k]  = ArrayContent[3 + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;
		n_1l_raw[i][j][k]  = 0. ;   // for now we are using only T1bbbb and we assume no contamination in the SL sample
		n_ldp_raw[i][j][k] = ArrayContent[3 + (nBinsMET*nBinsHT*nBinsBtag) + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;
		
		n_0l_correction[i][j][k]  = 1. ;
		n_1l_correction[i][j][k]  = 1. ;
		n_ldp_correction[i][j][k] = 1. ;

		// relative errors (in %)
		
		n_0l_error[i][j][k]  = ArrayContent[3 + 2*(nBinsMET*nBinsHT*nBinsBtag) + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;
		n_1l_error[i][j][k]  = 0.1 ;  // dummy tiny value
		n_ldp_error[i][j][k] = ArrayContent[3 + 3*(nBinsMET*nBinsHT*nBinsBtag) + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ;
	    
		if ( n_0l_raw[i][j][k] > 0.00001 ) {
		  n_0l_error[i][j][k] = 100 * n_0l_error[i][j][k] / n_0l_raw[i][j][k] ;
		}
		else n_0l_error[i][j][k] = 0.1 ;
	    
		if ( n_ldp_raw[i][j][k] > 0.00001 ) {
		  n_ldp_error[i][j][k] = 100 * n_ldp_error[i][j][k] / n_ldp_raw[i][j][k] ;
		}
		else n_ldp_error[i][j][k] = 0.1 ;

	      }
	    }
	  }

	  printf("\n\n Found event counts for point m0 = %4.0f,  m1/2 = %4.0f.", pointM0, pointM12 ) ;

	} // end check on mass point

      }	

      infp.close() ;

      if ( !found ) {
	printf("\n\n *** Point not found in scan.  Check this file %s\n  to see if this point is there: mgl=%4.0f, mlsp=%4.0f\n\n", inputScanFile, m0, m12 ) ;
	return false ;
      }

      printf("\n\n Skipping processing %s for now.\n\n", systFile1 ) ;

      // defer treatment of syst uncertainties for later
      /*
      // now open the file with the systematic uncertainties

      printf("\n\n Opening signal systematics input file : %s\n", systFile1 ) ;

      ifstream infq ;
      infq.open(systFile1) ;
      if ( !infq.good() ) {
	printf("\n\n *** Problem opening input file: %s.\n\n", systFile1 ) ;
	return false ;
      } 

      found = false ;
	  
      cout << "Systematic uncertainties (in %):" << endl;


      while ( infq.good() ) {

	if ( found ) break ;

	int ArraySize = 2 + 2*(nBinsMET*nBinsHT*nBinsBtag) ;
	double ArrayContent[ArraySize] ;

	for (int i = 0; infp && i < ArraySize; ++ i) {
	  infq >> ArrayContent[i];
	}

	pointM0  = ArrayContent[0] ;
	pointM12 = ArrayContent[1] ;
	
	//nGen = (int)ArrayContent[2] ;

	if (    fabs( pointM0 - m0 ) <= deltaM0/2 && fabs( pointM12 - m12 ) <= deltaM12/2 ) {

	  found = true ;

	  for (int i = 0 ; i < nBinsMET ; i++) {
	    for (int j = 0 ; j < nBinsHT ; j++) {
	      for (int k = 0 ; k < nBinsBtag ; k++) {     
		
		n_0l_error[i][j][k]  = fabs( ArrayContent[2 + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ) ;
		n_1l_error[i][j][k]  = 0. ;   // for now we are using only T1bbbb and we assume no contamination in the SL sample
		n_ldp_error[i][j][k] = fabs( ArrayContent[2 + (nBinsMET*nBinsHT*nBinsBtag) + i*(nBinsHT*nBinsBtag) + j*(nBinsBtag) + k] ) ;
		
	      }
	    }
	  }

	  printf("\n\n Found systematics for point m0 = %4.0f,  m1/2 = %4.0f.\n\n", pointM0, pointM12 ) ;

	}  // end check on mass point

      }

      infq.close() ;


      if ( !found ) {
	printf("\n\n *** Point not found in scan.  Check this file %s\n  to see if this point is there: mgl=%4.0f, mlsp=%4.0f\n\n", systFile1, m0, m12 ) ;
	return false ;
      }




      // print out input values

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
      */     
      

      double setVal_n_0l[nBinsMET][nBinsHT][nBinsBtag] ;
      double setVal_n_1l[nBinsMET][nBinsHT][nBinsBtag] ;
      double setVal_n_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
      
      double all0lep(0.) ;
      
      for (int i = 0 ; i < nBinsMET ; i++) {
	for (int j = 0 ; j < nBinsHT ; j++) {
	  if ( ignoreBin[i][j] ) continue ;
	  for (int k = 0 ; k < nBinsBtag ; k++) {     
		
	    setVal_n_0l[i][j][k]  = n_0l_raw[i][j][k]  * n_0l_correction[i][j][k] ;
	    setVal_n_1l[i][j][k]  = n_1l_raw[i][j][k]  * n_1l_correction[i][j][k] ;
	    setVal_n_ldp[i][j][k] = n_ldp_raw[i][j][k] * n_ldp_correction[i][j][k] ;
	    
	    all0lep += setVal_n_0l[i][j][k] ;
	    
	    rv_mu_susymc[i][j][k]     -> setVal( setVal_n_0l[i][j][k] ) ;
	    rv_mu_susymc_sl[i][j][k]  -> setVal( setVal_n_1l[i][j][k] ) ;
	    rv_mu_susymc_ldp[i][j][k] -> setVal( setVal_n_ldp[i][j][k] ) ;
	    
	    rv_width_eff_sf[i][j][k]     -> setVal( n_0l_error[i][j][k] / 100. ) ;
	    rv_width_eff_sf_sl[i][j][k]  -> setVal( n_1l_error[i][j][k] / 100. ) ;
	    rv_width_eff_sf_ldp[i][j][k] -> setVal( n_ldp_error[i][j][k] / 100. ) ;
	    
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
	    
	    cout << binString + " - 0 lep - setting susy signal to " << setVal_n_0l[i][j][k] << " +/- " << n_0l_error[i][j][k] << " %" << endl ;
	    cout << binString + " - 1 lep - setting susy signal to " << setVal_n_1l[i][j][k] << " +/- " << n_1l_error[i][j][k] << " %" << endl ;
	    cout << binString + " - ldp   - setting susy signal to " << setVal_n_ldp[i][j][k] << " +/- " << n_ldp_error[i][j][k] << " %" << endl ; 
	    
	    total0lep += setVal_n_0l[i][j][k] ;
	    
	  }
	}
      }

      printf("\n\n SUSY 0lep total: %7.2f\n\n", total0lep ) ;


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

	float pointMgl ;
	float pointMlsp ;

	float deff_dbtageff_0l[nBinsMET][nBinsHT][nBinsBtag] ;
	float deff_dbtageff_1l[nBinsMET][nBinsHT][nBinsBtag] ;
	float deff_dbtageff_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
	
	int ArraySize = 2 + 3*(nBinsMET*nBinsHT*nBinsBtag) ;
	double ArrayContent[ArraySize] ;

	for (int i = 0; infb && i < ArraySize; ++ i) {
	  infb >> ArrayContent[i];
	}

	pointMgl  = ArrayContent[0] ;
	pointMlsp = ArrayContent[1] ;

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

	//// printf("  pointMgl = %g , pointMlsp = %g\n", pointMgl, pointMlsp ) ;

	if (    fabs( pointMgl - m0 ) <= deltaM0/2. && fabs( pointMlsp - m12 ) <= deltaM12/2. ) {
	 
	  cout << "Setting susy deff_dbtag derivatives: " << endl ;

	  for (int i = 0 ; i < nBinsMET ; i++) {
	    for (int j = 0 ; j < nBinsHT ; j++) {
              if ( ignoreBin[i][j] ) continue ;
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



    RooAbsReal* ra2bRoostatsClass3D_3b::makeBetaConstraint( const char* NP_name, double NP_val, double NP_err ) {

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
       RooAbsReal* rar_alpha = new RooConstVar( varname, varname, alpha ) ;

       sprintf( varname, "beta_%s", NP_name ) ;
       RooAbsReal* rar_beta = new RooConstVar( varname, varname, beta ) ;

       RooAbsReal* rar_np = new RooRealVar( NP_name, NP_name, NP_val, 0., 1. ) ;

       sprintf( varname, "betapdf_%s", NP_name ) ;
       RooAbsReal* rar_pdf = new RooBetaPdf( varname, varname, *rar_np, *rar_alpha, *rar_beta ) ;

       allNuisances -> add( *rar_np ) ;
       allNuisancePdfs -> add( *rar_pdf ) ;

       return rar_np ;


 //  //--- Old way with two pseudo observations below here.

 //    RooRealVar *rrv_passObs, *rrv_failObs, *rrv_passPar, *rrv_failPar ;
 //    RooPoisson *passConstraint, *failConstraint ;

 //    double parVal, parErr, upperLimit, lowerLimit ;

 //    sprintf( varname, "passObs_%s", NP_name ) ;
 //    rrv_passObs = new RooRealVar( varname, varname, alpha-1., 1e-5, 1e5 ) ;
 //    rrv_passObs -> setConstant( kTRUE ) ;

 //    sprintf( varname, "failObs_%s", NP_name ) ;
 //    rrv_failObs = new RooRealVar( varname, varname, beta -1., 1e-5, 1e5 ) ;
 //    rrv_failObs -> setConstant( kTRUE ) ;

 //    globalObservables -> add( *rrv_passObs ) ;
 //    globalObservables -> add( *rrv_failObs ) ;


 //    parVal = alpha-1. ;
 //    lowerLimit = parVal - 6*sqrt(parVal) ;
 //    if ( lowerLimit <= 0. ) { lowerLimit = 1e-5 ; }
 //    upperLimit = parVal + 6*sqrt(parVal) ;
 //    if ( parVal > 0. ) { parErr = sqrt( parVal ) ; } else { parErr = 0. ; }
 //    sprintf( varname, "passPar_%s", NP_name ) ;
 //    rrv_passPar = new RooRealVar( varname, varname, parVal, lowerLimit, upperLimit ) ;
 //    rrv_passPar -> setError( parErr ) ;
 //    rrv_passPar->setConstant( kFALSE ) ;
 //    printf(" floating nuisance parameter: %s = %g +/- %g, [%g, %g]\n",
 //        varname, parVal, parErr, lowerLimit, upperLimit ) ;

 //    parVal = beta-1. ;
 //    lowerLimit = parVal - 6*sqrt(parVal) ;
 //    if ( lowerLimit <= 0. ) { lowerLimit = 1e-5 ; }
 //    upperLimit = parVal + 6*sqrt(parVal) ;
 //    if ( parVal > 0. ) { parErr = sqrt( parVal ) ; } else { parErr = 0. ; }
 //    sprintf( varname, "failPar_%s", NP_name ) ;
 //    rrv_failPar = new RooRealVar( varname, varname, parVal, lowerLimit, upperLimit ) ;
 //    rrv_failPar -> setError( parErr ) ;
 //    rrv_failPar->setConstant( kFALSE ) ;
 //    printf(" floating nuisance parameter: %s = %g +/- %g, [%g, %g]\n",
 //        varname, parVal, parErr, lowerLimit, upperLimit ) ;

 //    allNuisances -> add( *rrv_passPar ) ;
 //    allNuisances -> add( *rrv_failPar ) ;


 //    sprintf( varname, "passConstraint_%s", NP_name ) ;
 //    passConstraint = new RooPoisson( varname, varname, *rrv_passObs, *rrv_passPar ) ;
 //    printf("  Created constraint PDF : %s, val = %g, logval = %g\n",
 //      varname, passConstraint->getVal(), passConstraint->getLogVal() ) ;

 //    sprintf( varname, "failConstraint_%s", NP_name ) ;
 //    failConstraint = new RooPoisson( varname, varname, *rrv_failObs, *rrv_failPar ) ;
 //    printf("  Created constraint PDF : %s, val = %g, logval = %g\n",
 //      varname, failConstraint->getVal(), failConstraint->getLogVal() ) ;

 //    allNuisancePdfs -> add( *passConstraint ) ;
 //    allNuisancePdfs -> add( *failConstraint ) ;

 //    sprintf( varname, "%s_passPlusFail", NP_name ) ;
 //    RooAddition* passPlusFail = new RooAddition( varname, varname, RooArgSet(*rrv_passPar, *rrv_failPar) ) ;


 //    sprintf( varname, "%s", NP_name ) ;
 //    RooAbsReal* rar = new RooRatio( varname, varname, *rrv_passPar, *passPlusFail ) ;
 //    printf("  Created nuisance parameter %s : val = %g\n", varname, rar -> getVal() ) ;

 //    return rar ;

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
       RooConstVar* g_mean = new RooConstVar( vname, vname, NP_val ) ;
       sprintf( vname, "sigma_%s", NP_name ) ;
       RooConstVar* g_sigma = new RooConstVar( vname, vname, NP_err ) ;

       char pdfname[1000] ;
       sprintf( pdfname, "pdf_%s", NP_name ) ;
       RooGaussian* np_pdf = new RooGaussian( pdfname, pdfname, *np_rrv, *g_mean, *g_sigma ) ;

       allNuisances -> add( *np_rrv ) ;
       allNuisancePdfs -> add( *np_pdf ) ;

       printf("  makeGaussianConstraint : created nuisance parameter %s : val = %g\n", NP_name, np_rrv -> getVal() ) ;

       return np_rrv ;


    } // makeGaussianConstraint.


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
          RooConstVar* g_mean = new RooConstVar( vname, vname, 0.0 ) ;
          sprintf( vname, "sigma_%s", NP_base_name ) ;
          RooConstVar* g_sigma = new RooConstVar( vname, vname, 1.0 ) ;

          char pdfname[100] ;
          sprintf( pdfname, "pdf_%s", NP_base_name ) ;
          printf("\n\n makeCorrelatedGaussianConstraint : creating base nuisance parameter pdf - %s\n\n", pdfname ) ;
          ///// RooGaussian* base_np_pdf = new RooGaussian( pdfname, pdfname, *rrv_np_base_par, RooConst(0.), RooConst(1.) ) ;
          RooGaussian* base_np_pdf = new RooGaussian( pdfname, pdfname, *rrv_np_base_par, *g_mean, *g_sigma ) ;
          allNuisancePdfs -> add( *base_np_pdf ) ;

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

          printf(" makeCorrelatedGaussianConstraint : creating pos-def correlated gaussian NP  :  %s, val = %g\n", NP_name, rar->getVal() ) ;

       }




       return rar ;

    } // makeCorrelatedGaussianConstraint.


   //==============================================================================================================


