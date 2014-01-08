
#include "TMath.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TAxis.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TText.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooConstVar.h"
#include "RooStats/ModelConfig.h"
#include "TRegexp.h"

#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>

  using namespace RooFit;
  using namespace RooStats;


  void saveHist(const char* filename, const char* pat) ;

  bool getBetaPrimeModeRMS( const char* parName, RooWorkspace* ws, double &mode, double &rms, double &alpha, double &beta ) ;
  bool getBetaModeRMS( const char* parName, RooWorkspace* ws, double &mode, double &rms, double &alpha, double &beta ) ;

  double addChi2FromPullHist( TH1F* hp, double x=0.15, double y=0.85, double size=0.055 ) ;
  double addChi2FromObs( TH1F* obshist, TH1F* modelhist, double x=0.60, double y=0.85, double size=0.055 ) ;

  //------
  //
  // Note: If the 2nd argument (mu_susy_sig_val) is negative, mu_susy_sig will be floated
  //        in the fit.  If it's zero or positive, mu_susy_sig will be fixed to that value
  //        in the fit.
  //
  //------

   void ws_fitqual_plots3D_3b( const char* wsfile = "ws2.root",
			       double mu_susy_sig_val = 0.,
			       bool doNorm = false  ) {
     
     double globalChi2(0.0) ;
     double obsChi2(0.0) ;
     double npChi2(0.0) ;
     double chi2(0.0) ;
     

     // get the binning from "Binning.txt"

     int version, nBinsVar1, nBinsVar2, nBinsBtag ;
     TString label ;
     TString sVar1, sVar2, sVar3 ;

     ifstream inBinning ;
     inBinning.open("Binning.txt") ;
     
     inBinning >> label >> sVar1 ;
     inBinning >> label >> sVar2 ;
     inBinning >> label >> sVar3 ;

     inBinning >> label >> nBinsVar1 ;
     inBinning >> label >> nBinsVar2 ;
     inBinning >> label >> nBinsBtag ;

     float dummy ;
     
     for ( int i = 0 ; i < nBinsVar1 ; i++ ) {
       inBinning >> label >> dummy ;
     }
     
     for ( int i = 0 ; i < nBinsVar2 ; i++ ) {
       inBinning >> label >> dummy ;
     }
  
     inBinning >> label >> version ;
     
     inBinning.close() ;


     double hmax = 1.25 ;

     TFile* wstf = new TFile( wsfile ) ;

     RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
     ws->Print() ;


     //-- Hardwire in to ignore the highest Var1 bin in the lowest Var2 bin.
     bool ignoreBin[10][10] ;
     for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
       for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
	 ignoreBin[mbi][hbi] = false ;
       } // hbi
     } // mbi
     //ignoreBin[nBinsVar1-1][0] = true ;

     printf("\n\n *** Ignoring these bins in the analysis.\n\n") ;
     for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
       for ( int hbi=0; hbi<nBinsVar2; hbi++ ) {
	 if ( ignoreBin[mbi][hbi] ) {
	   printf("  Var1 %d, Var2 %d\n", mbi+1, hbi+1 ) ;
	 }
       } // hbi
     } // mbi
     printf("\n\n\n") ;
     
     printf("\n\n Binning: nBinsVar1 = %d, nBinsVar2 = %d\n\n", nBinsVar1, nBinsVar2 ) ;
     
     TString sMbins[nBinsVar1];
     TString sHbins[nBinsVar2];

     TString s1bins[nBinsVar1];
     TString s2bins[nBinsVar2];

     TString sBbins[4] = {"_1b","_2b","_3b","_4b"};
     
     for (int i = 0 ; i < nBinsVar1 ; i++) {
       TString base = "_M";
       TString base_1 = "_#alpha";
       stringstream sbin;
       sbin << i+1;
       base += sbin.str();
       base_1 += sbin.str();
       sMbins[i] = base;
       s1bins[i] = base_1;
     }
      
     for (int j = 0 ; j < nBinsVar2 ; j++) {
       TString base = "_H";
       TString base_2 = "_#beta";
       stringstream sbin;
       sbin << j+1;
       base += sbin.str();
       base_2 += sbin.str();
       sHbins[j] = base;
       s2bins[j] = base_2;
     }



     ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;

     printf("\n\n\n  ===== SbModel ====================\n\n") ;
     modelConfig->Print() ;



     RooDataSet* rds = (RooDataSet*) ws->obj( "ra2b_observed_rds" ) ;
     printf("\n\n\n  ===== RooDataSet ====================\n\n") ;

     rds->Print() ;
     rds->printMultiline(cout, 1, kTRUE, "") ;


     RooAbsPdf* likelihood = modelConfig->GetPdf() ;


     /////// RooRealVar* rrv_mu_susy_sig = ws->var("mu_susy_M1_H1_1b") ;
     RooRealVar* rrv_mu_susy_sig = ws->var("mu_susy_all0lep") ;
     if ( rrv_mu_susy_sig == 0x0 ) {
       printf("\n\n\n *** can't find mu_susy_all0lep in workspace.  Quitting.\n\n\n") ;
       return ;
     } else {
       printf(" current value is : %8.3f\n", rrv_mu_susy_sig->getVal() ) ; cout << flush ;
       rrv_mu_susy_sig->setConstant(kFALSE) ;
     }


     printf("\n\n\n  ===== Doing a fit with SUSY component floating ====================\n\n") ;

     RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(0) ) ;
     double logLikelihoodSusyFloat = fitResult->minNll() ;
     
     double logLikelihoodSusyFixed(0.) ;
     double testStatVal(-1.) ;
     if ( mu_susy_sig_val >= 0. ) {
       printf("\n\n\n  ===== Doing a fit with SUSY fixed ====================\n\n") ;
       printf(" fixing mu_susy_sig to %8.2f.\n", mu_susy_sig_val ) ;
       rrv_mu_susy_sig->setVal( mu_susy_sig_val ) ;
       rrv_mu_susy_sig->setConstant(kTRUE) ;
       
       fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(0) ) ;
       logLikelihoodSusyFixed = fitResult->minNll() ;
       testStatVal = 2.*(logLikelihoodSusyFixed - logLikelihoodSusyFloat) ;
       printf("\n\n\n ======= test statistic : -2 * ln (L_fixed / ln L_max) = %8.3f\n\n\n", testStatVal ) ;
     }


     printf("\n ==== Final floating parameter values\n\n") ;
     const RooArgList fitFloatVals = fitResult->floatParsFinal() ;
     {
       TIterator* parIter = fitFloatVals.createIterator() ;
       while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {
	 printf(" %20s : %8.2f\n", par->GetName(), par->getVal() ) ;
       }
     }


     printf("\n ==== Constant parameter values\n\n") ;
     const RooArgList fitConstVals = fitResult->constPars() ;
     {
       TIterator* parIter = fitConstVals.createIterator() ;
       while ( RooRealVar* par = (RooRealVar*) parIter->Next() ) {
	 printf(" %20s : %8.2f\n", par->GetName(), par->getVal() ) ;
       }
     }
     
     printf("\n ==== Function values\n\n") ;
     RooArgSet funcs = ws->allFunctions() ;
     TIterator* funcIter = funcs.createIterator() ;
     while ( RooFormulaVar* func = (RooFormulaVar*) funcIter->Next() ) {
       printf(" %20s : %8.2f\n", func->GetName(), func->getVal() ) ;
     }

     printf("\n\n") ; cout << flush ;


     double AttwjVal[nBinsVar1][nBinsVar2][nBinsBtag] ;
     double AqzoVal[nBinsVar1][nBinsVar2][nBinsBtag] ;

     for ( int i = 0 ; i < nBinsVar1 ; i++ ) {
       for ( int j = 0 ; j < nBinsVar2 ; j++ ) {

         if ( ignoreBin[i][j] ) continue ;

         char vname[1000] ;

         double trigeff(1.) ;
         sprintf( vname, "trigeff_M%d_H%d", i+1, j+1 ) ;
         RooAbsReal* rar = (RooAbsReal*) ws->obj(vname) ;
         if ( rar != 0x0 ) {
            trigeff = rar->getVal() ;
         } else {
            printf("\n\n *** %s missing\n", vname ) ;
         }

         double trigeff_sl(1.) ;
         sprintf( vname, "trigeff_sl_M%d_H%d", i+1, j+1 ) ;
         RooAbsReal* rarsl = (RooAbsReal*) ws->obj(vname) ;
         if ( rarsl != 0x0 ) {
            trigeff_sl = rarsl->getVal() ;
         } else {
            printf("\n\n *** %s missing\n", vname ) ;
         }

         for ( int k = 0 ; k < nBinsBtag ; k++ ) {

           TString ttString  = "mu_ttwj" ;
           TString qzoString = "mu_qzo" ;

           ttString  += sMbins[i]+sHbins[j]+sBbins[k] ;
           qzoString += sMbins[i]+sHbins[j]+sBbins[k] ;

           RooAbsReal* ttwj_obj = (RooAbsReal*) ws->obj(ttString) ;
           RooAbsReal* qzo_obj  = (RooAbsReal*) ws->obj(qzoString) ;

           if ( ttwj_obj == 0x0 ) {
              printf(" * %s missing\n", ttString.Data() ) ;
              AttwjVal[i][j][k] = 0. ;
           } else {
              AttwjVal[i][j][k] = trigeff_sl * ( ttwj_obj->getVal() ) ;
           }

           if ( qzo_obj == 0x0 ) {
              printf(" * %s missing\n", qzoString.Data() ) ;
              AqzoVal[i][j][k] = 0. ;
           } else {
              AqzoVal[i][j][k]  = trigeff * ( qzo_obj->getVal() ) ;
           }

         } // k
       } // j
     } // i


     //--- unpack observables.

     int dataN_0lep[nBinsVar1][nBinsVar2][nBinsBtag] ;
     int dataN_1lepSig[nBinsVar1][nBinsVar2][nBinsBtag] ;
     int dataN_1lep[nBinsVar1][nBinsVar2][nBinsBtag] ;


     for ( int i = 0 ; i < nBinsVar1 ; i++ ) {
       for ( int j = 0 ; j < nBinsVar2 ; j++ ) {
         for ( int k = 0 ; k < nBinsBtag ; k++ ) {
            dataN_0lep[i][j][k] = 0. ;
            dataN_1lepSig[i][j][k] = 0. ;
            dataN_1lep[i][j][k] = 0. ;
         }
       }
     }
     printf("\n\n Unpacking observables.\n\n") ; cout << flush ;

     const RooArgSet* dsras = rds->get() ;
     TIterator* obsIter = dsras->createIterator() ;

     while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {

       for ( int i = 0 ; i < nBinsVar1 ; i++ ) {
	 for ( int j = 0 ; j < nBinsVar2 ; j++ ) {
           if ( ignoreBin[i][j] ) continue ;
	   for ( int k = 0 ; k < nBinsBtag ; k++ ) {

	     TString String_0lep    = "N_0lep" ;
	     TString String_1lepSig = "N_1lepSig" ;
	     TString String_1lep    = "N_1lep" ;

	     String_0lep    += sMbins[i]+sHbins[j]+sBbins[k] ;
	     String_1lepSig += sMbins[i]+sHbins[j]+sBbins[k] ;
	     String_1lep    += sMbins[i]+sHbins[j]+sBbins[k] ;
	     
	     if ( strcmp( obs->GetName(), String_0lep    ) == 0 ) { dataN_0lep[i][j][k]    = obs->getVal() ; }
	     if ( strcmp( obs->GetName(), String_1lepSig ) == 0 ) { dataN_1lepSig[i][j][k] = obs->getVal() ; }
	     if ( strcmp( obs->GetName(), String_1lep    ) == 0 ) { dataN_1lep[i][j][k]    = obs->getVal() ; }

	   }
	 }
       }

     }


     gStyle->SetPadTopMargin(0.03) ;
     gStyle->SetPadBottomMargin(0.30) ;
     gStyle->SetPadRightMargin(0.10) ;



     int nbins = nBinsVar1*(nBinsVar2+1) + 1 ;

     TH1F *hfitqual_data_0lep_1b  , *hfitqual_data_0lep_2b  , *hfitqual_data_0lep_3b  , *hfitqual_data_0lep_4b ;
     TH1F *hfitqual_susy_0lep_1b  , *hfitqual_susy_0lep_2b  , *hfitqual_susy_0lep_3b  , *hfitqual_susy_0lep_4b ;
     TH1F *hfitqual_ttwj_0lep_1b  , *hfitqual_ttwj_0lep_2b  , *hfitqual_ttwj_0lep_3b  , *hfitqual_ttwj_0lep_4b ;
     TH1F *hfitqual_qzo_0lep_1b   , *hfitqual_qzo_0lep_2b   , *hfitqual_qzo_0lep_3b   , *hfitqual_qzo_0lep_4b ;
     TH1F *hfitqual_model_0lep_1b , *hfitqual_model_0lep_2b , *hfitqual_model_0lep_3b , *hfitqual_model_0lep_4b ;


     hfitqual_data_0lep_1b = new TH1F("hfitqual_data_0lep_1b", "0 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_0lep_1b = new TH1F("hfitqual_susy_0lep_1b", "0 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_0lep_1b = new TH1F("hfitqual_ttwj_0lep_1b", "0 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qzo_0lep_1b  = new TH1F("hfitqual_qzo_0lep_1b" , "0 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_0lep_1b  = new TH1F("hfitqual_model_0lep_1b" , "0 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;

     hfitqual_data_0lep_2b = new TH1F("hfitqual_data_0lep_2b", "0 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_0lep_2b = new TH1F("hfitqual_susy_0lep_2b", "0 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_0lep_2b = new TH1F("hfitqual_ttwj_0lep_2b", "0 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qzo_0lep_2b  = new TH1F("hfitqual_qzo_0lep_2b" , "0 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_0lep_2b  = new TH1F("hfitqual_model_0lep_2b" , "0 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_data_0lep_3b = new TH1F("hfitqual_data_0lep_3b", "0 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_susy_0lep_3b = new TH1F("hfitqual_susy_0lep_3b", "0 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_ttwj_0lep_3b = new TH1F("hfitqual_ttwj_0lep_3b", "0 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_qzo_0lep_3b  = new TH1F("hfitqual_qzo_0lep_3b" , "0 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       hfitqual_model_0lep_3b  = new TH1F("hfitqual_model_0lep_3b" , "0 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;

       if ( nBinsBtag > 3 ) {

	 hfitqual_data_0lep_4b = new TH1F("hfitqual_data_0lep_4b", "0 Lep, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_susy_0lep_4b = new TH1F("hfitqual_susy_0lep_4b", "0 Lep, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_ttwj_0lep_4b = new TH1F("hfitqual_ttwj_0lep_4b", "0 Lep, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_qzo_0lep_4b  = new TH1F("hfitqual_qzo_0lep_4b" , "0 Lep, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_model_0lep_4b  = new TH1F("hfitqual_model_0lep_4b" , "0 Lep, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 
       }
     }


     TH1F *hfitqual_data_1lepSig_1b  , *hfitqual_data_1lepSig_2b  , *hfitqual_data_1lepSig_3b  , *hfitqual_data_1lepSig_4b ;
     TH1F *hfitqual_susy_1lepSig_1b  , *hfitqual_susy_1lepSig_2b  , *hfitqual_susy_1lepSig_3b  , *hfitqual_susy_1lepSig_4b ;
     TH1F *hfitqual_ttwj_1lepSig_1b  , *hfitqual_ttwj_1lepSig_2b  , *hfitqual_ttwj_1lepSig_3b  , *hfitqual_ttwj_1lepSig_4b ;
     TH1F *hfitqual_qzo_1lepSig_1b   , *hfitqual_qzo_1lepSig_2b   , *hfitqual_qzo_1lepSig_3b   , *hfitqual_qzo_1lepSig_4b ;
     TH1F *hfitqual_model_1lepSig_1b , *hfitqual_model_1lepSig_2b , *hfitqual_model_1lepSig_3b , *hfitqual_model_1lepSig_4b ;


     hfitqual_data_1lepSig_1b = new TH1F("hfitqual_data_1lepSig_1b", "1 Lep Sig, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_1lepSig_1b = new TH1F("hfitqual_susy_1lepSig_1b", "1 Lep Sig, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_1lepSig_1b = new TH1F("hfitqual_ttwj_1lepSig_1b", "1 Lep Sig, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qzo_1lepSig_1b  = new TH1F("hfitqual_qzo_1lepSig_1b" , "1 Lep Sig, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_1lepSig_1b  = new TH1F("hfitqual_model_1lepSig_1b" , "1 Lep Sig, 1 btag" , nbins, 0.5, nbins+0.5 ) ;

     hfitqual_data_1lepSig_2b = new TH1F("hfitqual_data_1lepSig_2b", "1 Lep Sig, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_1lepSig_2b = new TH1F("hfitqual_susy_1lepSig_2b", "1 Lep Sig, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_1lepSig_2b = new TH1F("hfitqual_ttwj_1lepSig_2b", "1 Lep Sig, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qzo_1lepSig_2b  = new TH1F("hfitqual_qzo_1lepSig_2b" , "1 Lep Sig, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_1lepSig_2b  = new TH1F("hfitqual_model_1lepSig_2b" , "1 Lep Sig, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_data_1lepSig_3b = new TH1F("hfitqual_data_1lepSig_3b", "1 Lep Sig, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_susy_1lepSig_3b = new TH1F("hfitqual_susy_1lepSig_3b", "1 Lep Sig, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_ttwj_1lepSig_3b = new TH1F("hfitqual_ttwj_1lepSig_3b", "1 Lep Sig, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_qzo_1lepSig_3b  = new TH1F("hfitqual_qzo_1lepSig_3b" , "1 Lep Sig, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       hfitqual_model_1lepSig_3b  = new TH1F("hfitqual_model_1lepSig_3b" , "1 Lep Sig, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       
       if ( nBinsBtag > 3 ) {
       
	 hfitqual_data_1lepSig_4b = new TH1F("hfitqual_data_1lepSig_4b", "1 Lep Sig, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_susy_1lepSig_4b = new TH1F("hfitqual_susy_1lepSig_4b", "1 Lep Sig, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_ttwj_1lepSig_4b = new TH1F("hfitqual_ttwj_1lepSig_4b", "1 Lep Sig, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_qzo_1lepSig_4b  = new TH1F("hfitqual_qzo_1lepSig_4b" , "1 Lep Sig, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_model_1lepSig_4b  = new TH1F("hfitqual_model_1lepSig_4b" , "1 Lep Sig, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 
       }
     }


     TH1F *hfitqual_data_1lep_1b  , *hfitqual_data_1lep_2b  , *hfitqual_data_1lep_3b  , *hfitqual_data_1lep_4b ;
     TH1F *hfitqual_susy_1lep_1b  , *hfitqual_susy_1lep_2b  , *hfitqual_susy_1lep_3b  , *hfitqual_susy_1lep_4b ;
     TH1F *hfitqual_ttwj_1lep_1b  , *hfitqual_ttwj_1lep_2b  , *hfitqual_ttwj_1lep_3b  , *hfitqual_ttwj_1lep_4b ;
     TH1F *hfitqual_qzo_1lep_1b   , *hfitqual_qzo_1lep_2b   , *hfitqual_qzo_1lep_3b   , *hfitqual_qzo_1lep_4b ;
     TH1F *hfitqual_model_1lep_1b , *hfitqual_model_1lep_2b , *hfitqual_model_1lep_3b , *hfitqual_model_1lep_4b ;


     hfitqual_data_1lep_1b = new TH1F("hfitqual_data_1lep_1b", "1 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_1lep_1b = new TH1F("hfitqual_susy_1lep_1b", "1 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_1lep_1b = new TH1F("hfitqual_ttwj_1lep_1b", "1 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qzo_1lep_1b  = new TH1F("hfitqual_qzo_1lep_1b" , "1 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_1lep_1b  = new TH1F("hfitqual_model_1lep_1b" , "1 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;

     hfitqual_data_1lep_2b = new TH1F("hfitqual_data_1lep_2b", "1 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_1lep_2b = new TH1F("hfitqual_susy_1lep_2b", "1 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_1lep_2b = new TH1F("hfitqual_ttwj_1lep_2b", "1 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qzo_1lep_2b  = new TH1F("hfitqual_qzo_1lep_2b" , "1 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_1lep_2b  = new TH1F("hfitqual_model_1lep_2b" , "1 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_data_1lep_3b = new TH1F("hfitqual_data_1lep_3b", "1 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_susy_1lep_3b = new TH1F("hfitqual_susy_1lep_3b", "1 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_ttwj_1lep_3b = new TH1F("hfitqual_ttwj_1lep_3b", "1 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_qzo_1lep_3b  = new TH1F("hfitqual_qzo_1lep_3b" , "1 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       hfitqual_model_1lep_3b  = new TH1F("hfitqual_model_1lep_3b" , "1 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       
       if ( nBinsBtag > 3 ) {

	 hfitqual_data_1lep_4b = new TH1F("hfitqual_data_1lep_4b", "1 Lep, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_susy_1lep_4b = new TH1F("hfitqual_susy_1lep_4b", "1 Lep, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_ttwj_1lep_4b = new TH1F("hfitqual_ttwj_1lep_4b", "1 Lep, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_qzo_1lep_4b  = new TH1F("hfitqual_qzo_1lep_4b" , "1 Lep, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_model_1lep_4b  = new TH1F("hfitqual_model_1lep_4b" , "1 Lep, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 
       }
     }


     hfitqual_ttwj_0lep_1b  -> SetFillColor(kBlue-9) ;
     hfitqual_qzo_0lep_1b   -> SetFillColor(2) ;
     hfitqual_susy_0lep_1b  -> SetFillColor(6) ;
     hfitqual_data_0lep_1b->SetMarkerStyle(20) ;
     hfitqual_data_0lep_1b->SetLineWidth(2) ;
     
     hfitqual_ttwj_0lep_2b  -> SetFillColor(kBlue-9) ;
     hfitqual_qzo_0lep_2b   -> SetFillColor(2) ;
     hfitqual_susy_0lep_2b  -> SetFillColor(6) ;
     hfitqual_data_0lep_2b->SetMarkerStyle(20) ;
     hfitqual_data_0lep_2b->SetLineWidth(2) ;
     
     if ( nBinsBtag > 2 ) {
       
       hfitqual_ttwj_0lep_3b  -> SetFillColor(kBlue-9) ;
       hfitqual_qzo_0lep_3b   -> SetFillColor(2) ;
       hfitqual_susy_0lep_3b  -> SetFillColor(6) ;
       hfitqual_data_0lep_3b->SetMarkerStyle(20) ;
       hfitqual_data_0lep_3b->SetLineWidth(2) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_ttwj_0lep_4b  -> SetFillColor(kBlue-9) ;
	 hfitqual_qzo_0lep_4b   -> SetFillColor(2) ;
	 hfitqual_susy_0lep_4b  -> SetFillColor(6) ;
	 hfitqual_data_0lep_4b->SetMarkerStyle(20) ;
	 hfitqual_data_0lep_4b->SetLineWidth(2) ;

       }
     }
     
     
     hfitqual_ttwj_1lepSig_1b  -> SetFillColor(kBlue-9) ;
     hfitqual_qzo_1lepSig_1b   -> SetFillColor(2) ;
     hfitqual_susy_1lepSig_1b  -> SetFillColor(6) ;
     hfitqual_data_1lepSig_1b->SetMarkerStyle(20) ;
     hfitqual_data_1lepSig_1b->SetLineWidth(2) ;

     hfitqual_ttwj_1lepSig_2b  -> SetFillColor(kBlue-9) ;
     hfitqual_qzo_1lepSig_2b   -> SetFillColor(2) ;
     hfitqual_susy_1lepSig_2b  -> SetFillColor(6) ;
     hfitqual_data_1lepSig_2b->SetMarkerStyle(20) ;
     hfitqual_data_1lepSig_2b->SetLineWidth(2) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_ttwj_1lepSig_3b  -> SetFillColor(kBlue-9) ;
       hfitqual_qzo_1lepSig_3b   -> SetFillColor(2) ;
       hfitqual_susy_1lepSig_3b  -> SetFillColor(6) ;
       hfitqual_data_1lepSig_3b->SetMarkerStyle(20) ;
       hfitqual_data_1lepSig_3b->SetLineWidth(2) ;

       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_ttwj_1lepSig_4b  -> SetFillColor(kBlue-9) ;
	 hfitqual_qzo_1lepSig_4b   -> SetFillColor(2) ;
	 hfitqual_susy_1lepSig_4b  -> SetFillColor(6) ;
	 hfitqual_data_1lepSig_4b->SetMarkerStyle(20) ;
	 hfitqual_data_1lepSig_4b->SetLineWidth(2) ;
	 
       }
     }


     hfitqual_ttwj_1lep_1b  -> SetFillColor(kBlue-9) ;
     hfitqual_qzo_1lep_1b   -> SetFillColor(2) ;
     hfitqual_susy_1lep_1b  -> SetFillColor(6) ;
     hfitqual_data_1lep_1b->SetMarkerStyle(20) ;
     hfitqual_data_1lep_1b->SetLineWidth(2) ;

     hfitqual_ttwj_1lep_2b  -> SetFillColor(kBlue-9) ;
     hfitqual_qzo_1lep_2b   -> SetFillColor(2) ;
     hfitqual_susy_1lep_2b  -> SetFillColor(6) ;
     hfitqual_data_1lep_2b->SetMarkerStyle(20) ;
     hfitqual_data_1lep_2b->SetLineWidth(2) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_ttwj_1lep_3b  -> SetFillColor(kBlue-9) ;
       hfitqual_qzo_1lep_3b   -> SetFillColor(2) ;
       hfitqual_susy_1lep_3b  -> SetFillColor(6) ;
       hfitqual_data_1lep_3b->SetMarkerStyle(20) ;
       hfitqual_data_1lep_3b->SetLineWidth(2) ;
       
       if ( nBinsBtag > 3 ) {

	 hfitqual_ttwj_1lep_4b  -> SetFillColor(kBlue-9) ;
	 hfitqual_qzo_1lep_4b   -> SetFillColor(2) ;
	 hfitqual_susy_1lep_4b  -> SetFillColor(6) ;
	 hfitqual_data_1lep_4b->SetMarkerStyle(20) ;
	 hfitqual_data_1lep_4b->SetLineWidth(2) ;
	 
       }
     }
     

     THStack *hfitqual_fit_0lep_1b,  *hfitqual_fit_0lep_2b,  *hfitqual_fit_0lep_3b,  *hfitqual_fit_0lep_4b ;
     THStack *hfitqual_fit_1lepSig_1b,  *hfitqual_fit_1lepSig_2b,  *hfitqual_fit_1lepSig_3b,  *hfitqual_fit_1lepSig_4b ;
     THStack *hfitqual_fit_1lep_1b,  *hfitqual_fit_1lep_2b,  *hfitqual_fit_1lep_3b,  *hfitqual_fit_1lep_4b ;


     hfitqual_fit_0lep_1b = new THStack( "hfitqual_fit_0lep_1b", "RA2b likelihood fit results, fit" ) ;
     hfitqual_fit_1lepSig_1b = new THStack( "hfitqual_fit_1lepSig_1b", "RA2b likelihood fit results, fit" ) ;
     hfitqual_fit_1lep_1b = new THStack( "hfitqual_fit_1lep_1b", "RA2b likelihood fit results, fit" ) ;

     hfitqual_fit_0lep_2b = new THStack( "hfitqual_fit_0lep_2b", "RA2b likelihood fit results, fit" ) ;
     hfitqual_fit_1lepSig_2b = new THStack( "hfitqual_fit_1lepSig_2b", "RA2b likelihood fit results, fit" ) ;
     hfitqual_fit_1lep_2b = new THStack( "hfitqual_fit_1lep_2b", "RA2b likelihood fit results, fit" ) ;


     if ( nBinsBtag > 2 ) {

       hfitqual_fit_0lep_3b = new THStack( "hfitqual_fit_0lep_3b", "RA2b likelihood fit results, fit" ) ;
       hfitqual_fit_1lepSig_3b = new THStack( "hfitqual_fit_1lepSig_3b", "RA2b likelihood fit results, fit" ) ;
       hfitqual_fit_1lep_3b = new THStack( "hfitqual_fit_1lep_3b", "RA2b likelihood fit results, fit" ) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_fit_0lep_4b = new THStack( "hfitqual_fit_0lep_4b", "RA2b likelihood fit results, fit" ) ;
	 hfitqual_fit_1lepSig_4b = new THStack( "hfitqual_fit_1lepSig_4b", "RA2b likelihood fit results, fit" ) ;
	 hfitqual_fit_1lep_4b = new THStack( "hfitqual_fit_1lep_4b", "RA2b likelihood fit results, fit" ) ;
	 
       }
     }


     TAxis *xaxis_0lep_1b,    *xaxis_0lep_2b,    *xaxis_0lep_3b,    *xaxis_0lep_4b ;
     TAxis *xaxis_1lepSig_1b, *xaxis_1lepSig_2b, *xaxis_1lepSig_3b, *xaxis_1lepSig_4b ;
     TAxis *xaxis_1lep_1b,    *xaxis_1lep_2b,    *xaxis_1lep_3b,    *xaxis_1lep_4b ;


     xaxis_0lep_1b = hfitqual_data_0lep_1b->GetXaxis() ;
     xaxis_1lepSig_1b = hfitqual_data_1lepSig_1b->GetXaxis() ;
     xaxis_1lep_1b = hfitqual_data_1lep_1b->GetXaxis() ;

     xaxis_0lep_2b = hfitqual_data_0lep_2b->GetXaxis() ;
     xaxis_1lepSig_2b = hfitqual_data_1lepSig_2b->GetXaxis() ;
     xaxis_1lep_2b = hfitqual_data_1lep_2b->GetXaxis() ;

     if ( nBinsBtag > 2 ) {

       xaxis_0lep_3b = hfitqual_data_0lep_3b->GetXaxis() ;
       xaxis_1lepSig_3b = hfitqual_data_1lepSig_3b->GetXaxis() ;
       xaxis_1lep_3b = hfitqual_data_1lep_3b->GetXaxis() ;
       
       if ( nBinsBtag > 3 ) {
	 
	 xaxis_0lep_4b = hfitqual_data_0lep_4b->GetXaxis() ;
	 xaxis_1lepSig_4b = hfitqual_data_1lepSig_4b->GetXaxis() ;
	 xaxis_1lep_4b = hfitqual_data_1lep_4b->GetXaxis() ;
	 
       }
     }


     TString binLabel ;
     int  binIndex ;

     char teffvar[1000] ;

     double dataVal(0.) ;
     double dataErr(0.) ;
     double ttwjVal(0.) ;
     double qzoVal(0.) ;
     double susyVal(0.) ;
     double lhtotalVal(0.) ;
     double trigeff(1.) ;
     double trigeff_sl(1.) ;

     double eff_sf(0.) ;
     double eff_sf_slSig(0.) ;
     double eff_sf_sl(0.) ;

     double sf_mc(0.) ;


     printf("\n\n Looping through analysis bins.\n\n") ; cout << flush ;

     // loop through all the bins of the analysis

     for ( int i = 0 ; i < nBinsVar1 ; i++ ) {
       for ( int j = 0 ; j < nBinsVar2 ; j++ ) {
         if ( ignoreBin[i][j] ) continue ;
         for ( int k = 0 ; k < nBinsBtag ; k++ ) {

           printf(" here t1\n") ; cout << flush ;

           binIndex = 1 + (nBinsVar2+1)*i + j + 1 ;


           //+++++ 0 lep histograms:

           binLabel = "0 lep" ;
           binLabel += s1bins[i]+s2bins[j]+sBbins[k] ;

           TString EffSfString  = "eff_sf" ;
           TString MuSusyString = "mu_susy" ;

           EffSfString  += sMbins[i]+sHbins[j]+sBbins[k] ;
           MuSusyString += sMbins[i]+sHbins[j]+sBbins[k] ;

           printf(" here t2\n") ; cout << flush ;

           dataVal = dataN_0lep[i][j][k];
           eff_sf = 1.0 ;
           susyVal = 0. ;
           trigeff = 1. ;
           trigeff_sl = 1. ;
           sprintf( teffvar, "trigeff_M%d_H%d", i+1, j+1 ) ;
           if ( ws->obj( teffvar) != 0x0 ) {
              trigeff = ((RooAbsReal*) ws->obj(teffvar)) -> getVal() ;
           } else {
              printf(" * %s missing.\n", teffvar ) ;
              trigeff = 1. ;
           }
           sprintf( teffvar, "trigeff_sl_M%d_H%d", i+1, j+1 ) ;
           if ( ws->obj( teffvar) != 0x0 ) {
              trigeff_sl = ((RooAbsReal*) ws->obj(teffvar)) -> getVal() ;
           } else {
              printf(" * %s missing.\n", teffvar ) ;
              trigeff_sl = 1. ;
           }
           if ( ws->obj(EffSfString) != 0 ) {
              eff_sf = ((RooFormulaVar*) ws->obj(EffSfString)) -> getVal() ;
           } else {
              printf(" * %s missing.\n", EffSfString.Data() ) ;
              eff_sf = 1. ;
           }
           if ( ws->obj(MuSusyString) != 0 ) {
              susyVal = trigeff_sl * eff_sf * ( ((RooRealVar*) ws->obj(MuSusyString)) -> getVal() ) ;
           } else {
              printf(" * %s missing.\n", MuSusyString.Data() ) ;
              susyVal = 0. ;
           }

           ttwjVal = AttwjVal[i][j][k] ;
           char pname[1000] ;
           sprintf( pname, "mu_qzo_M%d_H%d_%db", i+1, j+1, k+1 ) ;
           if ( ws->obj(pname) != 0 ) {
              qzoVal = trigeff * ( ((RooAbsReal*) ws->obj(pname)) -> getVal() ) ;
           } else {
              printf(" * %s missing. \n", pname ) ;
           }

           lhtotalVal = ttwjVal + qzoVal + susyVal ;

           dataErr = sqrt(dataVal) ;


           cout << "\n" << endl ;
           cout << "binIndex = " << binIndex << endl ;
           cout << binLabel << "       susy : " << susyVal << endl ;
           cout << binLabel << "       ttwj : " << ttwjVal << endl ;
           cout << binLabel << "        qzo : " << qzoVal << endl ;
           cout << binLabel << "   LH total : " << lhtotalVal << endl ;
           cout << binLabel << "       data : " << dataVal << endl ;
           cout << flush ;

           if ( doNorm && dataVal > 0. ) {
             dataErr = dataErr / dataVal ;
             susyVal = susyVal / dataVal ;
             ttwjVal = ttwjVal / dataVal ;
             qzoVal  = qzoVal / dataVal ;
             dataVal = 1. ;
           }

           if ( k == 0 ) {
              xaxis_0lep_1b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_0lep_1b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_0lep_1b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_0lep_1b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_0lep_1b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qzo_0lep_1b  -> SetBinContent( binIndex, qzoVal ) ;
              hfitqual_model_0lep_1b  -> SetBinContent( binIndex, susyVal+ttwjVal+qzoVal ) ;
           } else if ( k == 1 ) {
              xaxis_0lep_2b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_0lep_2b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_0lep_2b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_0lep_2b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_0lep_2b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qzo_0lep_2b  -> SetBinContent( binIndex, qzoVal ) ;
              hfitqual_model_0lep_2b  -> SetBinContent( binIndex, susyVal+ttwjVal+qzoVal ) ;
           } else if ( k == 2 ) {
              xaxis_0lep_3b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_0lep_3b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_0lep_3b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_0lep_3b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_0lep_3b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qzo_0lep_3b  -> SetBinContent( binIndex, qzoVal ) ;
              hfitqual_model_0lep_3b  -> SetBinContent( binIndex, susyVal+ttwjVal+qzoVal ) ;
           } else if ( k == 3 ) {
              xaxis_0lep_4b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_0lep_4b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_0lep_4b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_0lep_4b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_0lep_4b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qzo_0lep_4b  -> SetBinContent( binIndex, qzoVal ) ;
              hfitqual_model_0lep_4b  -> SetBinContent( binIndex, susyVal+ttwjVal+qzoVal ) ;
           }




           //+++++ 1 lepSig histograms:

           binLabel = "1 lep sig" ;
           binLabel += s1bins[i]+s2bins[j]+sBbins[k] ;

           TString EffSfSlSigString  = "eff_sf_slSig" ;
           TString MuSusySlSigString = "mu_susy_slSig" ;
           TString MuTtwjSlSigString = "mu_ttwj_slSig" ;

           EffSfSlSigString  += sMbins[i]+sHbins[j]+sBbins[k] ;
           MuSusySlSigString += sMbins[i]+sHbins[j]+sBbins[k] ;
           MuTtwjSlSigString += sMbins[i]+sHbins[j]+sBbins[k] ;

           dataVal = dataN_1lepSig[i][j][k];
           eff_sf_slSig = 1. ;
           susyVal = 0. ;
           ttwjVal = 0. ;
           trigeff = 1. ;
           trigeff_sl = 1. ;

           sprintf( teffvar, "trigeff_sl_M%d_H%d", i+1, j+1 ) ;
           if ( ws->obj( teffvar) != 0x0 ) {
              trigeff_sl = ((RooAbsReal*) ws->obj(teffvar)) -> getVal() ;
           } else {
              printf(" * %s missing.\n", teffvar ) ;
              trigeff_sl = 1. ;
           }
           if ( ws->obj(EffSfSlSigString) != 0x0 ) {
              eff_sf_sl = ((RooFormulaVar*) ws->obj(EffSfSlSigString)) -> getVal() ;
           } else {
              printf(" * %s missing.\n", EffSfSlSigString.Data() ) ;
              eff_sf_sl = 1. ;
           }
           if ( ws->obj(MuSusySlSigString) != 0x0 ) {
              susyVal = trigeff_sl * eff_sf_slSig * ( ((RooRealVar*) ws->obj(MuSusySlSigString)) -> getVal() ) ;
           } else {
              printf(" * %s missing.\n", MuSusySlSigString.Data() ) ;
              susyVal = 0. ;
           }
           if ( ws->obj(MuTtwjSlSigString) != 0x0 ) {
              ttwjVal = trigeff_sl * ( ((RooRealVar*) ws->obj(MuTtwjSlSigString)) -> getVal() ) ;
           } else {
              printf(" * %s missing.\n", MuTtwjSlSigString.Data() ) ;
              ttwjVal = 0. ;
           }
           qzoVal  = 0. ;
           lhtotalVal = ttwjVal + qzoVal + susyVal ;

           dataErr = sqrt(dataVal) ;


           cout << "\n" << endl ;
           cout << binLabel << "       susy : " << susyVal << endl ;
           cout << binLabel << "       ttwj : " << ttwjVal << endl ;
           cout << binLabel << "        qzo : " << qzoVal << endl ;
           cout << binLabel << "   LH total : " << lhtotalVal << endl ;
           cout << binLabel << "       data : " << dataVal << endl ;

           if ( doNorm && dataVal > 0. ) {
             dataErr = dataErr / dataVal ;
             susyVal = susyVal / dataVal ;
             ttwjVal = ttwjVal / dataVal ;
             qzoVal  = qzoVal / dataVal ;
             dataVal = 1. ;
           }

           if ( k == 0 ) {
              xaxis_1lepSig_1b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lepSig_1b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lepSig_1b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lepSig_1b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lepSig_1b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qzo_1lepSig_1b  -> SetBinContent( binIndex, qzoVal ) ;
              hfitqual_model_1lepSig_1b  -> SetBinContent( binIndex, susyVal+ttwjVal+qzoVal ) ;
           } else if ( k == 1 ) {
              xaxis_1lepSig_2b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lepSig_2b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lepSig_2b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lepSig_2b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lepSig_2b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qzo_1lepSig_2b  -> SetBinContent( binIndex, qzoVal ) ;
              hfitqual_model_1lepSig_2b  -> SetBinContent( binIndex, susyVal+ttwjVal+qzoVal ) ;
           } else if ( k == 2 ) {
              xaxis_1lepSig_3b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lepSig_3b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lepSig_3b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lepSig_3b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lepSig_3b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qzo_1lepSig_3b  -> SetBinContent( binIndex, qzoVal ) ;
              hfitqual_model_1lepSig_3b  -> SetBinContent( binIndex, susyVal+ttwjVal+qzoVal ) ;
           } else if ( k == 3 ) {
              xaxis_1lepSig_4b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lepSig_4b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lepSig_4b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lepSig_4b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lepSig_4b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qzo_1lepSig_4b  -> SetBinContent( binIndex, qzoVal ) ;
              hfitqual_model_1lepSig_4b  -> SetBinContent( binIndex, susyVal+ttwjVal+qzoVal ) ;
           }



           //+++++ 1 lep histograms:

           binLabel = "1 lep" ;
           binLabel += s1bins[i]+s2bins[j]+sBbins[k] ;

           TString EffSfSlString  = "eff_sf_sl" ;
           TString MuSusySlString = "mu_susy_sl" ;
           TString MuTtwjSlString = "mu_ttwj_sl" ;

           EffSfSlString  += sMbins[i]+sHbins[j]+sBbins[k] ;
           MuSusySlString += sMbins[i]+sHbins[j]+sBbins[k] ;
           MuTtwjSlString += sMbins[i]+sHbins[j]+sBbins[k] ;

           dataVal = dataN_1lep[i][j][k];
           eff_sf_sl = 1. ;
           susyVal = 0. ;
           ttwjVal = 0. ;
           trigeff = 1. ;
           trigeff_sl = 1. ;

           sprintf( teffvar, "trigeff_sl_M%d_H%d", i+1, j+1 ) ;
           if ( ws->obj( teffvar) != 0x0 ) {
              trigeff_sl = ((RooAbsReal*) ws->obj(teffvar)) -> getVal() ;
           } else {
              printf(" * %s missing.\n", teffvar ) ;
              trigeff_sl = 1. ;
           }
           if ( ws->obj(EffSfSlString) != 0x0 ) {
              eff_sf_sl = ((RooFormulaVar*) ws->obj(EffSfSlString)) -> getVal() ;
           } else {
              printf(" * %s missing.\n", EffSfSlString.Data() ) ;
              eff_sf_sl = 1. ;
           }
           if ( ws->obj(MuSusySlString) != 0x0 ) {
              susyVal = trigeff_sl * eff_sf_sl * ( ((RooRealVar*) ws->obj(MuSusySlString)) -> getVal() ) ;
           } else {
              printf(" * %s missing.\n", MuSusySlString.Data() ) ;
              susyVal = 0. ;
           }
           if ( ws->obj(MuTtwjSlString) != 0x0 ) {
              ttwjVal = trigeff_sl * ( ((RooRealVar*) ws->obj(MuTtwjSlString)) -> getVal() ) ;
           } else {
              printf(" * %s missing.\n", MuTtwjSlString.Data() ) ;
              ttwjVal = 0. ;
           }
           qzoVal  = 0. ;
           lhtotalVal = ttwjVal + qzoVal + susyVal ;

           dataErr = sqrt(dataVal) ;


           cout << "\n" << endl ;
           cout << binLabel << "       susy : " << susyVal << endl ;
           cout << binLabel << "       ttwj : " << ttwjVal << endl ;
           cout << binLabel << "        qzo : " << qzoVal << endl ;
           cout << binLabel << "   LH total : " << lhtotalVal << endl ;
           cout << binLabel << "       data : " << dataVal << endl ;

           if ( doNorm && dataVal > 0. ) {
             dataErr = dataErr / dataVal ;
             susyVal = susyVal / dataVal ;
             ttwjVal = ttwjVal / dataVal ;
             qzoVal  = qzoVal / dataVal ;
             dataVal = 1. ;
           }

           if ( k == 0 ) {
              xaxis_1lep_1b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lep_1b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lep_1b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lep_1b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lep_1b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qzo_1lep_1b  -> SetBinContent( binIndex, qzoVal ) ;
              hfitqual_model_1lep_1b  -> SetBinContent( binIndex, susyVal+ttwjVal+qzoVal ) ;
           } else if ( k == 1 ) {
              xaxis_1lep_2b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lep_2b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lep_2b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lep_2b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lep_2b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qzo_1lep_2b  -> SetBinContent( binIndex, qzoVal ) ;
              hfitqual_model_1lep_2b  -> SetBinContent( binIndex, susyVal+ttwjVal+qzoVal ) ;
           } else if ( k == 2 ) {
              xaxis_1lep_3b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lep_3b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lep_3b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lep_3b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lep_3b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qzo_1lep_3b  -> SetBinContent( binIndex, qzoVal ) ;
              hfitqual_model_1lep_3b  -> SetBinContent( binIndex, susyVal+ttwjVal+qzoVal ) ;
           } else if ( k == 3 ) {
              xaxis_1lep_4b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lep_4b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lep_4b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lep_4b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lep_4b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qzo_1lep_4b  -> SetBinContent( binIndex, qzoVal ) ;
              hfitqual_model_1lep_4b  -> SetBinContent( binIndex, susyVal+ttwjVal+qzoVal ) ;
           }

	   // ????


	 }
       }
     }


     // ???


     printf(" here 6\n") ; cout << flush ;
     
     //-- Eff sf ---------
     
     /// float eff_sf_prim = ((RooRealVar*) ws->obj("eff_sf_prim")) -> getVal() ;
     
     /// hfitqual_np -> SetBinContent( 1, eff_sf_prim ) ;


     printf("\n\n\n") ;


     //--- final formatting and drawing.
     
     
     gStyle->SetOptStat(0) ;
     gStyle->SetTitleX(0.95) ;
     gStyle->SetTitleAlign(33) ;
     
     
     
     hfitqual_fit_0lep_1b->Add( hfitqual_qzo_0lep_1b ) ;
     hfitqual_fit_0lep_1b->Add( hfitqual_ttwj_0lep_1b ) ;
     hfitqual_fit_0lep_1b->Add( hfitqual_susy_0lep_1b ) ;
     
     hfitqual_fit_0lep_2b->Add( hfitqual_qzo_0lep_2b ) ;
     hfitqual_fit_0lep_2b->Add( hfitqual_ttwj_0lep_2b ) ;
     hfitqual_fit_0lep_2b->Add( hfitqual_susy_0lep_2b ) ;
     
     if ( nBinsBtag > 2 ) {
       
       hfitqual_fit_0lep_3b->Add( hfitqual_qzo_0lep_3b ) ;
       hfitqual_fit_0lep_3b->Add( hfitqual_ttwj_0lep_3b ) ;
       hfitqual_fit_0lep_3b->Add( hfitqual_susy_0lep_3b ) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_fit_0lep_4b->Add( hfitqual_qzo_0lep_4b ) ;
	 hfitqual_fit_0lep_4b->Add( hfitqual_ttwj_0lep_4b ) ;
	 hfitqual_fit_0lep_4b->Add( hfitqual_susy_0lep_4b ) ;
	 
       }
     }
     

     hfitqual_fit_1lepSig_1b->Add( hfitqual_qzo_1lepSig_1b ) ;
     hfitqual_fit_1lepSig_1b->Add( hfitqual_ttwj_1lepSig_1b ) ;
     hfitqual_fit_1lepSig_1b->Add( hfitqual_susy_1lepSig_1b ) ;
	   
     hfitqual_fit_1lepSig_2b->Add( hfitqual_qzo_1lepSig_2b ) ;
     hfitqual_fit_1lepSig_2b->Add( hfitqual_ttwj_1lepSig_2b ) ;
     hfitqual_fit_1lepSig_2b->Add( hfitqual_susy_1lepSig_2b ) ;
     
     if ( nBinsBtag > 2 ) {
       
       hfitqual_fit_1lepSig_3b->Add( hfitqual_qzo_1lepSig_3b ) ;
       hfitqual_fit_1lepSig_3b->Add( hfitqual_ttwj_1lepSig_3b ) ;
       hfitqual_fit_1lepSig_3b->Add( hfitqual_susy_1lepSig_3b ) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_fit_1lepSig_4b->Add( hfitqual_qzo_1lepSig_4b ) ;
	 hfitqual_fit_1lepSig_4b->Add( hfitqual_ttwj_1lepSig_4b ) ;
	 hfitqual_fit_1lepSig_4b->Add( hfitqual_susy_1lepSig_4b ) ;
	 
       }
     }
     

     hfitqual_fit_1lep_1b->Add( hfitqual_qzo_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_ttwj_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_susy_1lep_1b ) ;
     
     hfitqual_fit_1lep_2b->Add( hfitqual_qzo_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_ttwj_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_susy_1lep_2b ) ;
     
     if ( nBinsBtag > 2 ) {
       
       hfitqual_fit_1lep_3b->Add( hfitqual_qzo_1lep_3b ) ;
       hfitqual_fit_1lep_3b->Add( hfitqual_ttwj_1lep_3b ) ;
       hfitqual_fit_1lep_3b->Add( hfitqual_susy_1lep_3b ) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_fit_1lep_4b->Add( hfitqual_qzo_1lep_4b ) ;
	 hfitqual_fit_1lep_4b->Add( hfitqual_ttwj_1lep_4b ) ;
	 hfitqual_fit_1lep_4b->Add( hfitqual_susy_1lep_4b ) ;
	 
       }
     }
     

     TLegend* legend = new TLegend(0.6,0.45,0.9,0.95) ;
     
     legend->AddEntry( hfitqual_data_0lep_1b, "data" ) ;
     legend->AddEntry( hfitqual_susy_0lep_1b, "SUSY" ) ;
     legend->AddEntry( hfitqual_ttwj_0lep_1b, "ttwj" ) ;
     legend->AddEntry( hfitqual_qzo_0lep_1b,  "QZO" ) ;
     

     if ( doNorm ) {
       hfitqual_data_0lep_1b->SetMaximum( hmax ) ;
       hfitqual_data_1lepSig_1b->SetMaximum( hmax ) ;
       hfitqual_data_1lep_1b->SetMaximum( hmax ) ;
       
       hfitqual_data_0lep_2b->SetMaximum( hmax ) ;
       hfitqual_data_1lepSig_2b->SetMaximum( hmax ) ;
       hfitqual_data_1lep_2b->SetMaximum( hmax ) ;
       
       if ( nBinsBtag > 2 ) {
	 
	 hfitqual_data_0lep_3b->SetMaximum( hmax ) ;
	 hfitqual_data_1lepSig_3b->SetMaximum( hmax ) ;
	 hfitqual_data_1lep_3b->SetMaximum( hmax ) ;
	 
	 if ( nBinsBtag > 3 ) {
	   
	   hfitqual_data_0lep_4b->SetMaximum( hmax ) ;
	   hfitqual_data_1lepSig_4b->SetMaximum( hmax ) ;
	   hfitqual_data_1lep_4b->SetMaximum( hmax ) ;
	   
	 }
       }
       
     } else {
       
       hfitqual_data_0lep_1b->SetMaximum( hmax*(hfitqual_data_0lep_1b->GetMaximum()) ) ;
       hfitqual_data_1lepSig_1b->SetMaximum( hmax*(hfitqual_data_1lepSig_1b->GetMaximum()) ) ;
       hfitqual_data_1lep_1b->SetMaximum( hmax*(hfitqual_data_1lep_1b->GetMaximum()) ) ;
       
       hfitqual_data_0lep_2b->SetMaximum( hmax*(hfitqual_data_0lep_2b->GetMaximum()) ) ;
       hfitqual_data_1lepSig_2b->SetMaximum( hmax*(hfitqual_data_1lepSig_2b->GetMaximum()) ) ;
       hfitqual_data_1lep_2b->SetMaximum( hmax*(hfitqual_data_1lep_2b->GetMaximum()) ) ;
       
       if ( nBinsBtag > 2 ) {
	 
	 hfitqual_data_0lep_3b->SetMaximum( hmax*(hfitqual_data_0lep_3b->GetMaximum()) ) ;
	 hfitqual_data_1lepSig_3b->SetMaximum( hmax*(hfitqual_data_1lepSig_3b->GetMaximum()) ) ;
	 hfitqual_data_1lep_3b->SetMaximum( hmax*(hfitqual_data_1lep_3b->GetMaximum()) ) ;
	 
	 if ( nBinsBtag > 3 ) {
	   
	   hfitqual_data_0lep_4b->SetMaximum( hmax*(hfitqual_data_0lep_4b->GetMaximum()) ) ;
	   hfitqual_data_1lepSig_4b->SetMaximum( hmax*(hfitqual_data_1lepSig_4b->GetMaximum()) ) ;
	   hfitqual_data_1lep_4b->SetMaximum( hmax*(hfitqual_data_1lep_4b->GetMaximum()) ) ;
	   
	 }
       }
       
     }




     //====== Observables ==================================
     

     TCanvas* cfitqual = (TCanvas*) gDirectory->FindObject("cfitqual") ;
     if ( cfitqual == 0x0 ) {
       cfitqual = new TCanvas("cfitqual","RA2b fit quality", 850, 1000 ) ;
     }
     
     Int_t nColumns = max(3,nBinsBtag);
     
     cfitqual->Divide(nColumns,4);
     
     gPad->SetTicks(1,0) ;
     
     hfitqual_data_0lep_1b->SetLabelSize(0.055,"x") ;
     hfitqual_data_0lep_1b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_1lepSig_1b->SetLabelSize(0.055,"x") ;
     hfitqual_data_1lepSig_1b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_1lep_1b->SetLabelSize(0.055,"x") ;
     hfitqual_data_1lep_1b->GetXaxis()->LabelsOption("v") ;
     
     hfitqual_data_0lep_2b->SetLabelSize(0.055,"x") ;
     hfitqual_data_0lep_2b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_1lepSig_2b->SetLabelSize(0.055,"x") ;
     hfitqual_data_1lepSig_2b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_1lep_2b->SetLabelSize(0.055,"x") ;
     hfitqual_data_1lep_2b->GetXaxis()->LabelsOption("v") ;
     
     if ( nBinsBtag > 2 ) {
       
       hfitqual_data_0lep_3b->SetLabelSize(0.055,"x") ;
       hfitqual_data_0lep_3b->GetXaxis()->LabelsOption("v") ;
       hfitqual_data_1lepSig_3b->SetLabelSize(0.055,"x") ;
       hfitqual_data_1lepSig_3b->GetXaxis()->LabelsOption("v") ;
       hfitqual_data_1lep_3b->SetLabelSize(0.055,"x") ;
       hfitqual_data_1lep_3b->GetXaxis()->LabelsOption("v") ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_data_0lep_4b->SetLabelSize(0.055,"x") ;
	 hfitqual_data_0lep_4b->GetXaxis()->LabelsOption("v") ;
	 hfitqual_data_1lepSig_4b->SetLabelSize(0.055,"x") ;
	 hfitqual_data_1lepSig_4b->GetXaxis()->LabelsOption("v") ;
	 hfitqual_data_1lep_4b->SetLabelSize(0.055,"x") ;
	 hfitqual_data_1lep_4b->GetXaxis()->LabelsOption("v") ;
	 
       }
     }


     cfitqual->cd(1);
     hfitqual_data_0lep_1b->Draw("histpe") ;
     hfitqual_fit_0lep_1b->Draw("same") ;
     hfitqual_data_0lep_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_0lep_1b, hfitqual_model_0lep_1b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;
     
     cfitqual->cd(2);
     hfitqual_data_0lep_2b->Draw("histpe") ;
     hfitqual_fit_0lep_2b->Draw("same") ;
     hfitqual_data_0lep_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_0lep_2b, hfitqual_model_0lep_2b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;
     
     if ( nBinsBtag > 2 ) {
       
       cfitqual->cd(3);
       hfitqual_data_0lep_3b->Draw("histpe") ;
       hfitqual_fit_0lep_3b->Draw("same") ;
       hfitqual_data_0lep_3b->Draw("same") ;
       gPad->SetGridy(1) ;
       chi2 = addChi2FromObs( hfitqual_data_0lep_3b, hfitqual_model_0lep_3b ) ;
       globalChi2 += chi2 ;
       obsChi2 += chi2 ;
       
       if ( nBinsBtag > 3 ) {
	 
	 cfitqual->cd(4);
	 hfitqual_data_0lep_4b->Draw("histpe") ;
	 hfitqual_fit_0lep_4b->Draw("same") ;
	 hfitqual_data_0lep_4b->Draw("same") ;
	 gPad->SetGridy(1) ;
	 chi2 = addChi2FromObs( hfitqual_data_0lep_4b, hfitqual_model_0lep_4b ) ;
	 globalChi2 += chi2 ;
	 obsChi2 += chi2 ;
	 
       }
     }
     
     
     
     cfitqual->cd(nColumns+1);
     hfitqual_data_1lepSig_1b->Draw("histpe") ;
     hfitqual_fit_1lepSig_1b->Draw("same") ;
     hfitqual_data_1lepSig_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_1lepSig_1b, hfitqual_model_1lepSig_1b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;
     
     cfitqual->cd(nColumns+2);
     hfitqual_data_1lepSig_2b->Draw("histpe") ;
     hfitqual_fit_1lepSig_2b->Draw("same") ;
     hfitqual_data_1lepSig_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_1lepSig_2b, hfitqual_model_1lepSig_2b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;
     
     if ( nBinsBtag > 2 ) {
       
       cfitqual->cd(nColumns+3);
       hfitqual_data_1lepSig_3b->Draw("histpe") ;
       hfitqual_fit_1lepSig_3b->Draw("same") ;
       hfitqual_data_1lepSig_3b->Draw("same") ;
       gPad->SetGridy(1) ;
       chi2 = addChi2FromObs( hfitqual_data_1lepSig_3b, hfitqual_model_1lepSig_3b ) ;
       globalChi2 += chi2 ;
       obsChi2 += chi2 ;
	     
       if ( nBinsBtag > 3 ) {
	 
	 cfitqual->cd(nColumns+4);
	 hfitqual_data_1lepSig_4b->Draw("histpe") ;
	 hfitqual_fit_1lepSig_4b->Draw("same") ;
	 hfitqual_data_1lepSig_4b->Draw("same") ;
	 gPad->SetGridy(1) ;
	 chi2 = addChi2FromObs( hfitqual_data_1lepSig_4b, hfitqual_model_1lepSig_4b ) ;
	 globalChi2 += chi2 ;
	 obsChi2 += chi2 ;
	       
       }
     }
     


     cfitqual->cd(2*nColumns+1);
     hfitqual_data_1lep_1b->Draw("histpe") ;
     hfitqual_fit_1lep_1b->Draw("same") ;
     hfitqual_data_1lep_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_1lep_1b, hfitqual_model_1lep_1b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;
     
     cfitqual->cd(2*nColumns+2);
     hfitqual_data_1lep_2b->Draw("histpe") ;
     hfitqual_fit_1lep_2b->Draw("same") ;
     hfitqual_data_1lep_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_1lep_2b, hfitqual_model_1lep_2b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;
	   
     if ( nBinsBtag > 2 ) {
	     
       cfitqual->cd(2*nColumns+3);
       hfitqual_data_1lep_3b->Draw("histpe") ;
       hfitqual_fit_1lep_3b->Draw("same") ;
       hfitqual_data_1lep_3b->Draw("same") ;
       gPad->SetGridy(1) ;
       chi2 = addChi2FromObs( hfitqual_data_1lep_3b, hfitqual_model_1lep_3b ) ;
       globalChi2 += chi2 ;
       obsChi2 += chi2 ;
       
       if ( nBinsBtag > 3 ) {
	 
	 cfitqual->cd(2*nColumns+4);
	 hfitqual_data_1lep_4b->Draw("histpe") ;
	 hfitqual_fit_1lep_4b->Draw("same") ;
	 hfitqual_data_1lep_4b->Draw("same") ;
	 gPad->SetGridy(1) ;
	 chi2 = addChi2FromObs( hfitqual_data_1lep_4b, hfitqual_model_1lep_4b ) ;
	 globalChi2 += chi2 ;
	 obsChi2 += chi2 ;
	 
       }
     }
     


     
     cfitqual->cd(nColumns*4);
     legend->Draw() ;
     
     cfitqual->Update() ;
     


     cfitqual->cd(nColumns*4) ;
     TText* chi2text = new TText() ;
     chi2text->SetTextSize(0.055) ;
     chi2text->SetTextAlign(32) ;
     char chi2string[1000] ;
     
     sprintf( chi2string, "Overall Chi2 = %6.2f", globalChi2 ) ;
     chi2text->DrawTextNDC(0.55, 0.90, chi2string ) ;
     sprintf( chi2string, "obs Chi2 = %6.2f", obsChi2 ) ;
     chi2text->DrawTextNDC(0.55, 0.85, chi2string ) ;
     //sprintf( chi2string, "NP Chi2 = %6.2f", npChi2 ) ;
     //chi2text->DrawTextNDC(0.55, 0.80, chi2string ) ;
     
     printf("\n\n Overall Chi2 = %6.2f\n\n", globalChi2 ) ;
     printf("     Obs Chi2 = %6.2f\n", obsChi2 ) ;
     //printf("      NP Chi2 = %6.2f\n", npChi2 ) ;
     
     
     TText* frtext = new TText() ;
     frtext->SetTextSize(0.055) ;
     frtext->SetTextAlign(32) ;
     char frstring[1000] ;
     if ( mu_susy_sig_val < 0. ) {
       sprintf( frstring, "NSUSY 0lep : %6.1f +/- %5.1f", rrv_mu_susy_sig->getVal(), rrv_mu_susy_sig->getError() ) ;
       printf( "\n\n NSUSY 0lep : %6.1f +/- %5.1f\n\n", rrv_mu_susy_sig->getVal(), rrv_mu_susy_sig->getError() ) ;
     } else {
       sprintf( frstring, "NSUSY 0lep : %6.1f (fixed)", rrv_mu_susy_sig->getVal() ) ;
       printf( "\n\n NSUSY 0lep : %6.1f (fixed)\n\n", rrv_mu_susy_sig->getVal() ) ;
     }
     frtext->DrawTextNDC(0.55, 0.70, frstring ) ;
     
     

     
     TH1F* hfit_results = new TH1F("hfit_results", "Fit results", 50, 0.5, 50.5 ) ;
     
     hfit_results -> SetBinContent( 1, rrv_mu_susy_sig->getVal() ) ;
     hfit_results -> GetXaxis() -> SetBinLabel( 1, "Fit total 0lep SUSY events") ;
     
     if ( mu_susy_sig_val < 0. ) {
       hfit_results -> SetBinContent( 2, rrv_mu_susy_sig->getError() ) ;
     } else {
       hfit_results -> SetBinContent( 2, 0. ) ;
     }
     hfit_results -> GetXaxis() -> SetBinLabel( 2, "Fit total 0lep SUSY events error") ;
     
     






     TString histfile(wsfile) ;
     char newfileend[1000] ;
     if ( mu_susy_sig_val < 0 ) {
       sprintf( newfileend, "-susyFloat-fitqual.root" ) ;
     } else {
       sprintf( newfileend, "-Nsusy%.0f-fitqual.root", mu_susy_sig_val ) ;
     }
     histfile.ReplaceAll(".root", newfileend ) ;
     saveHist( histfile,"h*") ;
     
     

   }

    


//==========================================================================================

void saveHist(const char* filename, const char* pat)
{

  cout << "\n\n Saving histograms matching " << pat << " in file " << filename << "\n\n" << flush ;

  TList* list = gDirectory->GetList() ;
  TIterator* iter = list->MakeIterator();

  TRegexp re(pat,kTRUE) ;

  TFile outf(filename,"RECREATE") ;
  TObject* obj ;
  while((obj=iter->Next())) {
    if (TString(obj->GetName()).Index(re)>=0) {
      obj->Write() ;
      std::cout << "." ;
    }
  }
  std::cout << std::endl ;
  outf.Close() ;

  delete iter ;
}

//==========================================================================================


  bool getBetaPrimeModeRMS( const char* parName, RooWorkspace* ws, double &mode, double &rms, double &alpha, double &beta ) {

     mode = 1.0 ;
     rms = 0.0 ;

     char varname[1000] ;

     sprintf( varname, "passObs_%s", parName ) ;
     RooAbsReal* passObs = (RooAbsReal*) ws->obj( varname ) ;
     if ( passObs == 0x0 ) {
        printf("\n\n *** getNPModeRMS : can't find pass obs for %s\n\n", parName ) ;
        return false ;
     }
     alpha = passObs->getVal() + 1. ;

     sprintf( varname, "failObs_%s", parName ) ;
     RooAbsReal* failObs = (RooAbsReal*) ws->obj( varname ) ;
     if ( failObs == 0x0 ) {
        printf("\n\n *** getNPModeRMS : can't find fail obs for %s\n\n", parName ) ;
        return false ;
     }
     beta = failObs->getVal() + 1. ;

     mode = (alpha - 1.)/(beta + 1.) ;

     rms = sqrt( alpha * (alpha + beta - 1 ) / ( pow(beta - 1 , 2 ) * (beta - 2) ) ) ;

     return true ;

  } // getNPModeRMS.


//==========================================================================================


  bool getBetaModeRMS( const char* parName, RooWorkspace* ws, double &mode, double &rms, double &alpha, double &beta ) {

     mode = 1.0 ;
     rms = 0.0 ;

     char varname[1000] ;

     sprintf( varname, "alpha_%s", parName ) ;
     RooAbsReal* rv_alpha = (RooAbsReal*) ws->obj( varname ) ;
     if ( rv_alpha == 0x0 ) {
        printf("\n\n *** getNPModeRMS : can't find alpha for %s\n\n", parName ) ;
        return false ;
     }
     alpha = rv_alpha->getVal()  ;

     sprintf( varname, "beta_%s", parName ) ;
     RooAbsReal* rv_beta = (RooAbsReal*) ws->obj( varname ) ;
     if ( rv_beta == 0x0 ) {
        printf("\n\n *** getNPModeRMS : can't find beta for %s\n\n", parName ) ;
        return false ;
     }
     beta = rv_beta->getVal()  ;

     mode = (alpha - 1.)/(alpha+beta -2.) ;

     rms = sqrt( alpha * beta / ( pow(alpha + beta,2) * (alpha + beta + 1) ) ) ;

     return true ;

  } // getNPModeRMS.


//==========================================================================================


  double addChi2FromPullHist( TH1F* hp, double x, double y, double size ) {

     char   chi2string[100] ;
     double chi2val(0.) ;

     TText* chi2text = new TText() ;
     chi2text -> SetTextSize( size ) ;

     for ( int bi=1; bi<=hp->GetNbinsX(); bi++ ) {
        double chi = hp -> GetBinContent( bi ) ;
        chi2val += chi*chi ;
     } // bi.

     sprintf( chi2string, "Chi2 = %6.2f", chi2val ) ;
     chi2text->DrawTextNDC( x, y, chi2string ) ;

     return chi2val ;

  } // addChi2FromPullHist


//==========================================================================================


  double addChi2FromObs( TH1F* obshist, TH1F* modelhist, double x, double y, double size ) {

     char   chi2string[100] ;
     double chi2val(0.) ;

     TText* chi2text = new TText() ;
     chi2text -> SetTextSize( size ) ;

     for ( int bi=1; bi<=obshist->GetNbinsX(); bi++ ) {
        double obs = obshist -> GetBinContent( bi ) ;
        double model = modelhist -> GetBinContent( bi ) ;
        if ( obs > 0 ) {
           double chi = (obs-model)/sqrt(obs) ;
           chi2val += chi*chi ;
           printf(" %s : obs=%9.1f, model=%9.1f, chi2=%6.2f\n", obshist->GetXaxis()->GetBinLabel( bi ), obs, model, chi*chi ) ;
        }
     } // bi.

     sprintf( chi2string, "Chi2 = %6.2f", chi2val ) ;
     chi2text->DrawTextNDC( x, y, chi2string ) ;

     return chi2val ;

  } // addChi2FromPullHist


//==========================================================================================









