
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
     ignoreBin[nBinsVar1-1][0] = true ;

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
     double AqcdVal[nBinsVar1][nBinsVar2][nBinsBtag] ;

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
           TString qcdString = "mu_qcd" ;

           ttString  += sMbins[i]+sHbins[j]+sBbins[k] ;
           qcdString += sMbins[i]+sHbins[j]+sBbins[k] ;

           RooAbsReal* ttwj_obj = (RooAbsReal*) ws->obj(ttString) ;
           RooAbsReal* qcd_obj  = (RooAbsReal*) ws->obj(qcdString) ;

           if ( ttwj_obj == 0x0 ) {
              printf(" * %s missing\n", ttString.Data() ) ;
              AttwjVal[i][j][k] = 0. ;
           } else {
              AttwjVal[i][j][k] = trigeff_sl * ( ttwj_obj->getVal() ) ;
           }

           if ( qcd_obj == 0x0 ) {
              printf(" * %s missing\n", qcdString.Data() ) ;
              AqcdVal[i][j][k] = 0. ;
           } else {
              AqcdVal[i][j][k]  = trigeff * ( qcd_obj->getVal() ) ;
           }

         } // k
       } // j
     } // i


     //--- unpack observables.

     int dataN_0lep[nBinsVar1][nBinsVar2][nBinsBtag] ;
     int dataN_1lepSig[nBinsVar1][nBinsVar2][nBinsBtag] ;
     int dataN_1lep[nBinsVar1][nBinsVar2][nBinsBtag] ;
     int dataN_ldp[nBinsVar1][nBinsVar2][nBinsBtag] ;
     int dataN_Zee[nBinsVar1][nBinsVar2] ;
     int dataN_Zmm[nBinsVar1][nBinsVar2] ;

     for ( int i = 0 ; i < nBinsVar1 ; i++ ) {
       for ( int j = 0 ; j < nBinsVar2 ; j++ ) {
         for ( int k = 0 ; k < nBinsBtag ; k++ ) {
            dataN_0lep[i][j][k] = 0. ;
            dataN_1lepSig[i][j][k] = 0. ;
            dataN_1lep[i][j][k] = 0. ;
            dataN_ldp[i][j][k] = 0. ;
            dataN_Zee[i][j] = 0. ;
            dataN_Zmm[i][j] = 0. ;
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
	     TString String_ldp     = "N_ldp" ;

	     String_0lep    += sMbins[i]+sHbins[j]+sBbins[k] ;
	     String_1lepSig += sMbins[i]+sHbins[j]+sBbins[k] ;
	     String_1lep    += sMbins[i]+sHbins[j]+sBbins[k] ;
	     String_ldp     += sMbins[i]+sHbins[j]+sBbins[k] ;
	     
	     if ( strcmp( obs->GetName(), String_0lep    ) == 0 ) { dataN_0lep[i][j][k]    = obs->getVal() ; }
	     if ( strcmp( obs->GetName(), String_1lepSig ) == 0 ) { dataN_1lepSig[i][j][k] = obs->getVal() ; }
	     if ( strcmp( obs->GetName(), String_1lep    ) == 0 ) { dataN_1lep[i][j][k]    = obs->getVal() ; }
	     if ( strcmp( obs->GetName(), String_ldp     ) == 0 ) { dataN_ldp[i][j][k]     = obs->getVal() ; }

	   }
	   TString String_Zee   = "N_Zee" ;
	   TString String_Zmm   = "N_Zmm" ;
	   String_Zee   += sMbins[i]+sHbins[j] ;
	   String_Zmm   += sMbins[i]+sHbins[j] ;
	   if ( strcmp( obs->GetName(), String_Zee  ) == 0 ) { dataN_Zee[i][j]  = obs->getVal() ; }
	   if ( strcmp( obs->GetName(), String_Zmm  ) == 0 ) { dataN_Zmm[i][j]  = obs->getVal() ; }
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
     TH1F *hfitqual_qcd_0lep_1b   , *hfitqual_qcd_0lep_2b   , *hfitqual_qcd_0lep_3b   , *hfitqual_qcd_0lep_4b ;
     TH1F *hfitqual_znn_0lep_1b   , *hfitqual_znn_0lep_2b   , *hfitqual_znn_0lep_3b   , *hfitqual_znn_0lep_4b ;
     TH1F *hfitqual_model_0lep_1b , *hfitqual_model_0lep_2b , *hfitqual_model_0lep_3b , *hfitqual_model_0lep_4b ;


     hfitqual_data_0lep_1b = new TH1F("hfitqual_data_0lep_1b", "0 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_0lep_1b = new TH1F("hfitqual_susy_0lep_1b", "0 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_0lep_1b = new TH1F("hfitqual_ttwj_0lep_1b", "0 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qcd_0lep_1b  = new TH1F("hfitqual_qcd_0lep_1b" , "0 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_znn_0lep_1b  = new TH1F("hfitqual_znn_0lep_1b" , "0 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_0lep_1b  = new TH1F("hfitqual_model_0lep_1b" , "0 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;

     hfitqual_data_0lep_2b = new TH1F("hfitqual_data_0lep_2b", "0 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_0lep_2b = new TH1F("hfitqual_susy_0lep_2b", "0 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_0lep_2b = new TH1F("hfitqual_ttwj_0lep_2b", "0 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qcd_0lep_2b  = new TH1F("hfitqual_qcd_0lep_2b" , "0 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_znn_0lep_2b  = new TH1F("hfitqual_znn_0lep_2b" , "0 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_0lep_2b  = new TH1F("hfitqual_model_0lep_2b" , "0 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_data_0lep_3b = new TH1F("hfitqual_data_0lep_3b", "0 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_susy_0lep_3b = new TH1F("hfitqual_susy_0lep_3b", "0 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_ttwj_0lep_3b = new TH1F("hfitqual_ttwj_0lep_3b", "0 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_qcd_0lep_3b  = new TH1F("hfitqual_qcd_0lep_3b" , "0 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       hfitqual_znn_0lep_3b  = new TH1F("hfitqual_znn_0lep_3b" , "0 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       hfitqual_model_0lep_3b  = new TH1F("hfitqual_model_0lep_3b" , "0 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;

       if ( nBinsBtag > 3 ) {

	 hfitqual_data_0lep_4b = new TH1F("hfitqual_data_0lep_4b", "0 Lep, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_susy_0lep_4b = new TH1F("hfitqual_susy_0lep_4b", "0 Lep, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_ttwj_0lep_4b = new TH1F("hfitqual_ttwj_0lep_4b", "0 Lep, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_qcd_0lep_4b  = new TH1F("hfitqual_qcd_0lep_4b" , "0 Lep, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_znn_0lep_4b  = new TH1F("hfitqual_znn_0lep_4b" , "0 Lep, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_model_0lep_4b  = new TH1F("hfitqual_model_0lep_4b" , "0 Lep, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 
       }
     }


     TH1F *hfitqual_data_1lepSig_1b  , *hfitqual_data_1lepSig_2b  , *hfitqual_data_1lepSig_3b  , *hfitqual_data_1lepSig_4b ;
     TH1F *hfitqual_susy_1lepSig_1b  , *hfitqual_susy_1lepSig_2b  , *hfitqual_susy_1lepSig_3b  , *hfitqual_susy_1lepSig_4b ;
     TH1F *hfitqual_ttwj_1lepSig_1b  , *hfitqual_ttwj_1lepSig_2b  , *hfitqual_ttwj_1lepSig_3b  , *hfitqual_ttwj_1lepSig_4b ;
     TH1F *hfitqual_qcd_1lepSig_1b   , *hfitqual_qcd_1lepSig_2b   , *hfitqual_qcd_1lepSig_3b   , *hfitqual_qcd_1lepSig_4b ;
     TH1F *hfitqual_znn_1lepSig_1b   , *hfitqual_znn_1lepSig_2b   , *hfitqual_znn_1lepSig_3b   , *hfitqual_znn_1lepSig_4b ;
     TH1F *hfitqual_model_1lepSig_1b , *hfitqual_model_1lepSig_2b , *hfitqual_model_1lepSig_3b , *hfitqual_model_1lepSig_4b ;


     hfitqual_data_1lepSig_1b = new TH1F("hfitqual_data_1lepSig_1b", "1 Lep Sig, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_1lepSig_1b = new TH1F("hfitqual_susy_1lepSig_1b", "1 Lep Sig, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_1lepSig_1b = new TH1F("hfitqual_ttwj_1lepSig_1b", "1 Lep Sig, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qcd_1lepSig_1b  = new TH1F("hfitqual_qcd_1lepSig_1b" , "1 Lep Sig, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_znn_1lepSig_1b  = new TH1F("hfitqual_znn_1lepSig_1b" , "1 Lep Sig, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_1lepSig_1b  = new TH1F("hfitqual_model_1lepSig_1b" , "1 Lep Sig, 1 btag" , nbins, 0.5, nbins+0.5 ) ;

     hfitqual_data_1lepSig_2b = new TH1F("hfitqual_data_1lepSig_2b", "1 Lep Sig, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_1lepSig_2b = new TH1F("hfitqual_susy_1lepSig_2b", "1 Lep Sig, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_1lepSig_2b = new TH1F("hfitqual_ttwj_1lepSig_2b", "1 Lep Sig, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qcd_1lepSig_2b  = new TH1F("hfitqual_qcd_1lepSig_2b" , "1 Lep Sig, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_znn_1lepSig_2b  = new TH1F("hfitqual_znn_1lepSig_2b" , "1 Lep Sig, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_1lepSig_2b  = new TH1F("hfitqual_model_1lepSig_2b" , "1 Lep Sig, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_data_1lepSig_3b = new TH1F("hfitqual_data_1lepSig_3b", "1 Lep Sig, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_susy_1lepSig_3b = new TH1F("hfitqual_susy_1lepSig_3b", "1 Lep Sig, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_ttwj_1lepSig_3b = new TH1F("hfitqual_ttwj_1lepSig_3b", "1 Lep Sig, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_qcd_1lepSig_3b  = new TH1F("hfitqual_qcd_1lepSig_3b" , "1 Lep Sig, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       hfitqual_znn_1lepSig_3b  = new TH1F("hfitqual_znn_1lepSig_3b" , "1 Lep Sig, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       hfitqual_model_1lepSig_3b  = new TH1F("hfitqual_model_1lepSig_3b" , "1 Lep Sig, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       
       if ( nBinsBtag > 3 ) {
       
	 hfitqual_data_1lepSig_4b = new TH1F("hfitqual_data_1lepSig_4b", "1 Lep Sig, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_susy_1lepSig_4b = new TH1F("hfitqual_susy_1lepSig_4b", "1 Lep Sig, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_ttwj_1lepSig_4b = new TH1F("hfitqual_ttwj_1lepSig_4b", "1 Lep Sig, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_qcd_1lepSig_4b  = new TH1F("hfitqual_qcd_1lepSig_4b" , "1 Lep Sig, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_znn_1lepSig_4b  = new TH1F("hfitqual_znn_1lepSig_4b" , "1 Lep Sig, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_model_1lepSig_4b  = new TH1F("hfitqual_model_1lepSig_4b" , "1 Lep Sig, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 
       }
     }


     TH1F *hfitqual_data_1lep_1b  , *hfitqual_data_1lep_2b  , *hfitqual_data_1lep_3b  , *hfitqual_data_1lep_4b ;
     TH1F *hfitqual_susy_1lep_1b  , *hfitqual_susy_1lep_2b  , *hfitqual_susy_1lep_3b  , *hfitqual_susy_1lep_4b ;
     TH1F *hfitqual_ttwj_1lep_1b  , *hfitqual_ttwj_1lep_2b  , *hfitqual_ttwj_1lep_3b  , *hfitqual_ttwj_1lep_4b ;
     TH1F *hfitqual_qcd_1lep_1b   , *hfitqual_qcd_1lep_2b   , *hfitqual_qcd_1lep_3b   , *hfitqual_qcd_1lep_4b ;
     TH1F *hfitqual_znn_1lep_1b   , *hfitqual_znn_1lep_2b   , *hfitqual_znn_1lep_3b   , *hfitqual_znn_1lep_4b ;
     TH1F *hfitqual_model_1lep_1b , *hfitqual_model_1lep_2b , *hfitqual_model_1lep_3b , *hfitqual_model_1lep_4b ;


     hfitqual_data_1lep_1b = new TH1F("hfitqual_data_1lep_1b", "1 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_1lep_1b = new TH1F("hfitqual_susy_1lep_1b", "1 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_1lep_1b = new TH1F("hfitqual_ttwj_1lep_1b", "1 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qcd_1lep_1b  = new TH1F("hfitqual_qcd_1lep_1b" , "1 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_znn_1lep_1b  = new TH1F("hfitqual_znn_1lep_1b" , "1 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_1lep_1b  = new TH1F("hfitqual_model_1lep_1b" , "1 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;

     hfitqual_data_1lep_2b = new TH1F("hfitqual_data_1lep_2b", "1 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_1lep_2b = new TH1F("hfitqual_susy_1lep_2b", "1 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_1lep_2b = new TH1F("hfitqual_ttwj_1lep_2b", "1 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qcd_1lep_2b  = new TH1F("hfitqual_qcd_1lep_2b" , "1 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_znn_1lep_2b  = new TH1F("hfitqual_znn_1lep_2b" , "1 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_1lep_2b  = new TH1F("hfitqual_model_1lep_2b" , "1 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_data_1lep_3b = new TH1F("hfitqual_data_1lep_3b", "1 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_susy_1lep_3b = new TH1F("hfitqual_susy_1lep_3b", "1 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_ttwj_1lep_3b = new TH1F("hfitqual_ttwj_1lep_3b", "1 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_qcd_1lep_3b  = new TH1F("hfitqual_qcd_1lep_3b" , "1 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       hfitqual_znn_1lep_3b  = new TH1F("hfitqual_znn_1lep_3b" , "1 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       hfitqual_model_1lep_3b  = new TH1F("hfitqual_model_1lep_3b" , "1 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       
       if ( nBinsBtag > 3 ) {

	 hfitqual_data_1lep_4b = new TH1F("hfitqual_data_1lep_4b", "1 Lep, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_susy_1lep_4b = new TH1F("hfitqual_susy_1lep_4b", "1 Lep, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_ttwj_1lep_4b = new TH1F("hfitqual_ttwj_1lep_4b", "1 Lep, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_qcd_1lep_4b  = new TH1F("hfitqual_qcd_1lep_4b" , "1 Lep, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_znn_1lep_4b  = new TH1F("hfitqual_znn_1lep_4b" , "1 Lep, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_model_1lep_4b  = new TH1F("hfitqual_model_1lep_4b" , "1 Lep, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 
       }
     }


     TH1F *hfitqual_data_ldp_1b  , *hfitqual_data_ldp_2b  , *hfitqual_data_ldp_3b  , *hfitqual_data_ldp_4b ;
     TH1F *hfitqual_susy_ldp_1b  , *hfitqual_susy_ldp_2b  , *hfitqual_susy_ldp_3b  , *hfitqual_susy_ldp_4b ;
     TH1F *hfitqual_ttwj_ldp_1b  , *hfitqual_ttwj_ldp_2b  , *hfitqual_ttwj_ldp_3b  , *hfitqual_ttwj_ldp_4b ;
     TH1F *hfitqual_qcd_ldp_1b   , *hfitqual_qcd_ldp_2b   , *hfitqual_qcd_ldp_3b   , *hfitqual_qcd_ldp_4b ;
     TH1F *hfitqual_znn_ldp_1b   , *hfitqual_znn_ldp_2b   , *hfitqual_znn_ldp_3b   , *hfitqual_znn_ldp_4b ;
     TH1F *hfitqual_model_ldp_1b , *hfitqual_model_ldp_2b , *hfitqual_model_ldp_3b , *hfitqual_model_ldp_4b ;


     hfitqual_data_ldp_1b = new TH1F("hfitqual_data_ldp_1b", "LDP, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_ldp_1b = new TH1F("hfitqual_susy_ldp_1b", "LDP, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_ldp_1b = new TH1F("hfitqual_ttwj_ldp_1b", "LDP, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qcd_ldp_1b  = new TH1F("hfitqual_qcd_ldp_1b" , "LDP, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_znn_ldp_1b  = new TH1F("hfitqual_znn_ldp_1b" , "LDP, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_ldp_1b  = new TH1F("hfitqual_model_ldp_1b" , "LDP, 1 btag" , nbins, 0.5, nbins+0.5 ) ;

     hfitqual_data_ldp_2b = new TH1F("hfitqual_data_ldp_2b", "LDP, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_susy_ldp_2b = new TH1F("hfitqual_susy_ldp_2b", "LDP, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_ttwj_ldp_2b = new TH1F("hfitqual_ttwj_ldp_2b", "LDP, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     hfitqual_qcd_ldp_2b  = new TH1F("hfitqual_qcd_ldp_2b" , "LDP, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_znn_ldp_2b  = new TH1F("hfitqual_znn_ldp_2b" , "LDP, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     hfitqual_model_ldp_2b  = new TH1F("hfitqual_model_ldp_2b" , "LDP, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_data_ldp_3b = new TH1F("hfitqual_data_ldp_3b", "LDP, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_susy_ldp_3b = new TH1F("hfitqual_susy_ldp_3b", "LDP, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_ttwj_ldp_3b = new TH1F("hfitqual_ttwj_ldp_3b", "LDP, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
       hfitqual_qcd_ldp_3b  = new TH1F("hfitqual_qcd_ldp_3b" , "LDP, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       hfitqual_znn_ldp_3b  = new TH1F("hfitqual_znn_ldp_3b" , "LDP, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       hfitqual_model_ldp_3b  = new TH1F("hfitqual_model_ldp_3b" , "LDP, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_data_ldp_4b = new TH1F("hfitqual_data_ldp_4b", "LDP, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_susy_ldp_4b = new TH1F("hfitqual_susy_ldp_4b", "LDP, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_ttwj_ldp_4b = new TH1F("hfitqual_ttwj_ldp_4b", "LDP, >=4 btag", nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_qcd_ldp_4b  = new TH1F("hfitqual_qcd_ldp_4b" , "LDP, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_znn_ldp_4b  = new TH1F("hfitqual_znn_ldp_4b" , "LDP, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 hfitqual_model_ldp_4b  = new TH1F("hfitqual_model_ldp_4b" , "LDP, >=4 btag" , nbins, 0.5, nbins+0.5 ) ;
	 
       }
     }


     TH1F* hfitqual_data_zee_1b  = new TH1F("hfitqual_data_zee_1b" , "Zee" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_data_zmm_1b  = new TH1F("hfitqual_data_zmm_1b" , "Zmm" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_fit_zee_1b  = new TH1F("hfitqual_fit_zee_1b" , "Zee" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_fit_zmm_1b  = new TH1F("hfitqual_fit_zmm_1b" , "Zmm" , nbins, 0.5, nbins+0.5 ) ;

     TH1F* hfitqual_np   = new TH1F("hfitqual_np"  , "Nuisance par"  , 1, 0., 1. ) ;


     hfitqual_data_zee_1b->SetMarkerStyle(20) ;
     hfitqual_data_zee_1b->SetLineWidth(2) ;
     hfitqual_data_zmm_1b->SetMarkerStyle(20) ;
     hfitqual_data_zmm_1b->SetLineWidth(2) ;
     hfitqual_fit_zee_1b->SetFillColor(kGreen-3) ;
     hfitqual_fit_zmm_1b->SetFillColor(kGreen-3) ;



     hfitqual_ttwj_0lep_1b  -> SetFillColor(kBlue-9) ;
     hfitqual_qcd_0lep_1b   -> SetFillColor(2) ;
     hfitqual_znn_0lep_1b   -> SetFillColor(kGreen-3) ;
     hfitqual_susy_0lep_1b  -> SetFillColor(6) ;
     hfitqual_data_0lep_1b->SetMarkerStyle(20) ;
     hfitqual_data_0lep_1b->SetLineWidth(2) ;
     
     hfitqual_ttwj_0lep_2b  -> SetFillColor(kBlue-9) ;
     hfitqual_qcd_0lep_2b   -> SetFillColor(2) ;
     hfitqual_znn_0lep_2b   -> SetFillColor(kGreen-3) ;
     hfitqual_susy_0lep_2b  -> SetFillColor(6) ;
     hfitqual_data_0lep_2b->SetMarkerStyle(20) ;
     hfitqual_data_0lep_2b->SetLineWidth(2) ;
     
     if ( nBinsBtag > 2 ) {
       
       hfitqual_ttwj_0lep_3b  -> SetFillColor(kBlue-9) ;
       hfitqual_qcd_0lep_3b   -> SetFillColor(2) ;
       hfitqual_znn_0lep_3b   -> SetFillColor(kGreen-3) ;
       hfitqual_susy_0lep_3b  -> SetFillColor(6) ;
       hfitqual_data_0lep_3b->SetMarkerStyle(20) ;
       hfitqual_data_0lep_3b->SetLineWidth(2) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_ttwj_0lep_4b  -> SetFillColor(kBlue-9) ;
	 hfitqual_qcd_0lep_4b   -> SetFillColor(2) ;
	 hfitqual_znn_0lep_4b   -> SetFillColor(kGreen-3) ;
	 hfitqual_susy_0lep_4b  -> SetFillColor(6) ;
	 hfitqual_data_0lep_4b->SetMarkerStyle(20) ;
	 hfitqual_data_0lep_4b->SetLineWidth(2) ;

       }
     }
     
     
     hfitqual_ttwj_1lepSig_1b  -> SetFillColor(kBlue-9) ;
     hfitqual_qcd_1lepSig_1b   -> SetFillColor(2) ;
     hfitqual_znn_1lepSig_1b   -> SetFillColor(kGreen-3) ;
     hfitqual_susy_1lepSig_1b  -> SetFillColor(6) ;
     hfitqual_data_1lepSig_1b->SetMarkerStyle(20) ;
     hfitqual_data_1lepSig_1b->SetLineWidth(2) ;

     hfitqual_ttwj_1lepSig_2b  -> SetFillColor(kBlue-9) ;
     hfitqual_qcd_1lepSig_2b   -> SetFillColor(2) ;
     hfitqual_znn_1lepSig_2b   -> SetFillColor(kGreen-3) ;
     hfitqual_susy_1lepSig_2b  -> SetFillColor(6) ;
     hfitqual_data_1lepSig_2b->SetMarkerStyle(20) ;
     hfitqual_data_1lepSig_2b->SetLineWidth(2) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_ttwj_1lepSig_3b  -> SetFillColor(kBlue-9) ;
       hfitqual_qcd_1lepSig_3b   -> SetFillColor(2) ;
       hfitqual_znn_1lepSig_3b   -> SetFillColor(kGreen-3) ;
       hfitqual_susy_1lepSig_3b  -> SetFillColor(6) ;
       hfitqual_data_1lepSig_3b->SetMarkerStyle(20) ;
       hfitqual_data_1lepSig_3b->SetLineWidth(2) ;

       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_ttwj_1lepSig_4b  -> SetFillColor(kBlue-9) ;
	 hfitqual_qcd_1lepSig_4b   -> SetFillColor(2) ;
	 hfitqual_znn_1lepSig_4b   -> SetFillColor(kGreen-3) ;
	 hfitqual_susy_1lepSig_4b  -> SetFillColor(6) ;
	 hfitqual_data_1lepSig_4b->SetMarkerStyle(20) ;
	 hfitqual_data_1lepSig_4b->SetLineWidth(2) ;
	 
       }
     }


     hfitqual_ttwj_1lep_1b  -> SetFillColor(kBlue-9) ;
     hfitqual_qcd_1lep_1b   -> SetFillColor(2) ;
     hfitqual_znn_1lep_1b   -> SetFillColor(kGreen-3) ;
     hfitqual_susy_1lep_1b  -> SetFillColor(6) ;
     hfitqual_data_1lep_1b->SetMarkerStyle(20) ;
     hfitqual_data_1lep_1b->SetLineWidth(2) ;

     hfitqual_ttwj_1lep_2b  -> SetFillColor(kBlue-9) ;
     hfitqual_qcd_1lep_2b   -> SetFillColor(2) ;
     hfitqual_znn_1lep_2b   -> SetFillColor(kGreen-3) ;
     hfitqual_susy_1lep_2b  -> SetFillColor(6) ;
     hfitqual_data_1lep_2b->SetMarkerStyle(20) ;
     hfitqual_data_1lep_2b->SetLineWidth(2) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_ttwj_1lep_3b  -> SetFillColor(kBlue-9) ;
       hfitqual_qcd_1lep_3b   -> SetFillColor(2) ;
       hfitqual_znn_1lep_3b   -> SetFillColor(kGreen-3) ;
       hfitqual_susy_1lep_3b  -> SetFillColor(6) ;
       hfitqual_data_1lep_3b->SetMarkerStyle(20) ;
       hfitqual_data_1lep_3b->SetLineWidth(2) ;
       
       if ( nBinsBtag > 3 ) {

	 hfitqual_ttwj_1lep_4b  -> SetFillColor(kBlue-9) ;
	 hfitqual_qcd_1lep_4b   -> SetFillColor(2) ;
	 hfitqual_znn_1lep_4b   -> SetFillColor(kGreen-3) ;
	 hfitqual_susy_1lep_4b  -> SetFillColor(6) ;
	 hfitqual_data_1lep_4b->SetMarkerStyle(20) ;
	 hfitqual_data_1lep_4b->SetLineWidth(2) ;
	 
       }
     }
     

     hfitqual_ttwj_ldp_1b  -> SetFillColor(kBlue-9) ;
     hfitqual_qcd_ldp_1b   -> SetFillColor(2) ;
     hfitqual_znn_ldp_1b   -> SetFillColor(kGreen-3) ;
     hfitqual_susy_ldp_1b  -> SetFillColor(6) ;
     hfitqual_data_ldp_1b->SetMarkerStyle(20) ;
     hfitqual_data_ldp_1b->SetLineWidth(2) ;

     hfitqual_ttwj_ldp_2b  -> SetFillColor(kBlue-9) ;
     hfitqual_qcd_ldp_2b   -> SetFillColor(2) ;
     hfitqual_znn_ldp_2b   -> SetFillColor(kGreen-3) ;
     hfitqual_susy_ldp_2b  -> SetFillColor(6) ;
     hfitqual_data_ldp_2b->SetMarkerStyle(20) ;
     hfitqual_data_ldp_2b->SetLineWidth(2) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_ttwj_ldp_3b  -> SetFillColor(kBlue-9) ;
       hfitqual_qcd_ldp_3b   -> SetFillColor(2) ;
       hfitqual_znn_ldp_3b   -> SetFillColor(kGreen-3) ;
       hfitqual_susy_ldp_3b  -> SetFillColor(6) ;
       hfitqual_data_ldp_3b->SetMarkerStyle(20) ;
       hfitqual_data_ldp_3b->SetLineWidth(2) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_ttwj_ldp_4b  -> SetFillColor(kBlue-9) ;
	 hfitqual_qcd_ldp_4b   -> SetFillColor(2) ;
	 hfitqual_znn_ldp_4b   -> SetFillColor(kGreen-3) ;
	 hfitqual_susy_ldp_4b  -> SetFillColor(6) ;
	 hfitqual_data_ldp_4b->SetMarkerStyle(20) ;
	 hfitqual_data_ldp_4b->SetLineWidth(2) ;
	 
       }
     }


     hfitqual_np    -> SetFillColor(kOrange+1) ;

     THStack *hfitqual_fit_0lep_1b,  *hfitqual_fit_0lep_2b,  *hfitqual_fit_0lep_3b,  *hfitqual_fit_0lep_4b ;
     THStack *hfitqual_fit_1lepSig_1b,  *hfitqual_fit_1lepSig_2b,  *hfitqual_fit_1lepSig_3b,  *hfitqual_fit_1lepSig_4b ;
     THStack *hfitqual_fit_1lep_1b,  *hfitqual_fit_1lep_2b,  *hfitqual_fit_1lep_3b,  *hfitqual_fit_1lep_4b ;
     THStack *hfitqual_fit_ldp_1b,  *hfitqual_fit_ldp_2b,  *hfitqual_fit_ldp_3b,  *hfitqual_fit_ldp_4b ;


     hfitqual_fit_0lep_1b = new THStack( "hfitqual_fit_0lep_1b", "RA2b likelihood fit results, fit" ) ;
     hfitqual_fit_1lepSig_1b = new THStack( "hfitqual_fit_1lepSig_1b", "RA2b likelihood fit results, fit" ) ;
     hfitqual_fit_1lep_1b = new THStack( "hfitqual_fit_1lep_1b", "RA2b likelihood fit results, fit" ) ;
     hfitqual_fit_ldp_1b  = new THStack( "hfitqual_fit_ldp_1b",  "RA2b likelihood fit results, fit" ) ;

     hfitqual_fit_0lep_2b = new THStack( "hfitqual_fit_0lep_2b", "RA2b likelihood fit results, fit" ) ;
     hfitqual_fit_1lepSig_2b = new THStack( "hfitqual_fit_1lepSig_2b", "RA2b likelihood fit results, fit" ) ;
     hfitqual_fit_1lep_2b = new THStack( "hfitqual_fit_1lep_2b", "RA2b likelihood fit results, fit" ) ;
     hfitqual_fit_ldp_2b  = new THStack( "hfitqual_fit_ldp_2b",  "RA2b likelihood fit results, fit" ) ;


     if ( nBinsBtag > 2 ) {

       hfitqual_fit_0lep_3b = new THStack( "hfitqual_fit_0lep_3b", "RA2b likelihood fit results, fit" ) ;
       hfitqual_fit_1lepSig_3b = new THStack( "hfitqual_fit_1lepSig_3b", "RA2b likelihood fit results, fit" ) ;
       hfitqual_fit_1lep_3b = new THStack( "hfitqual_fit_1lep_3b", "RA2b likelihood fit results, fit" ) ;
       hfitqual_fit_ldp_3b  = new THStack( "hfitqual_fit_ldp_3b",  "RA2b likelihood fit results, fit" ) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_fit_0lep_4b = new THStack( "hfitqual_fit_0lep_4b", "RA2b likelihood fit results, fit" ) ;
	 hfitqual_fit_1lepSig_4b = new THStack( "hfitqual_fit_1lepSig_4b", "RA2b likelihood fit results, fit" ) ;
	 hfitqual_fit_1lep_4b = new THStack( "hfitqual_fit_1lep_4b", "RA2b likelihood fit results, fit" ) ;
	 hfitqual_fit_ldp_4b  = new THStack( "hfitqual_fit_ldp_4b",  "RA2b likelihood fit results, fit" ) ;
	 
       }
     }


     TAxis *xaxis_0lep_1b,    *xaxis_0lep_2b,    *xaxis_0lep_3b,    *xaxis_0lep_4b ;
     TAxis *xaxis_1lepSig_1b, *xaxis_1lepSig_2b, *xaxis_1lepSig_3b, *xaxis_1lepSig_4b ;
     TAxis *xaxis_1lep_1b,    *xaxis_1lep_2b,    *xaxis_1lep_3b,    *xaxis_1lep_4b ;
     TAxis *xaxis_ldp_1b,     *xaxis_ldp_2b,     *xaxis_ldp_3b,     *xaxis_ldp_4b ;


     xaxis_0lep_1b = hfitqual_data_0lep_1b->GetXaxis() ;
     xaxis_1lepSig_1b = hfitqual_data_1lepSig_1b->GetXaxis() ;
     xaxis_1lep_1b = hfitqual_data_1lep_1b->GetXaxis() ;
     xaxis_ldp_1b  = hfitqual_data_ldp_1b->GetXaxis() ;

     xaxis_0lep_2b = hfitqual_data_0lep_2b->GetXaxis() ;
     xaxis_1lepSig_2b = hfitqual_data_1lepSig_2b->GetXaxis() ;
     xaxis_1lep_2b = hfitqual_data_1lep_2b->GetXaxis() ;
     xaxis_ldp_2b  = hfitqual_data_ldp_2b->GetXaxis() ;

     if ( nBinsBtag > 2 ) {

       xaxis_0lep_3b = hfitqual_data_0lep_3b->GetXaxis() ;
       xaxis_1lepSig_3b = hfitqual_data_1lepSig_3b->GetXaxis() ;
       xaxis_1lep_3b = hfitqual_data_1lep_3b->GetXaxis() ;
       xaxis_ldp_3b  = hfitqual_data_ldp_3b->GetXaxis() ;
       
       if ( nBinsBtag > 3 ) {
	 
	 xaxis_0lep_4b = hfitqual_data_0lep_4b->GetXaxis() ;
	 xaxis_1lepSig_4b = hfitqual_data_1lepSig_4b->GetXaxis() ;
	 xaxis_1lep_4b = hfitqual_data_1lep_4b->GetXaxis() ;
	 xaxis_ldp_4b  = hfitqual_data_ldp_4b->GetXaxis() ;
	 
       }
     }

     TAxis* xaxis_zee_1b  = hfitqual_data_zee_1b->GetXaxis() ;
     TAxis* xaxis_zmm_1b  = hfitqual_data_zmm_1b->GetXaxis() ;

     TString binLabel ;
     int  binIndex ;

     char teffvar[1000] ;

     double dataVal(0.) ;
     double dataErr(0.) ;
     double ttwjVal(0.) ;
     double qcdVal(0.) ;
     double znnVal(0.) ;
     double susyVal(0.) ;
     double lhtotalVal(0.) ;
     double trigeff(1.) ;
     double trigeff_sl(1.) ;

     double eff_sf(0.) ;
     double eff_sf_slSig(0.) ;
     double eff_sf_sl(0.) ;
     double eff_sf_ldp(0.) ;

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
           TString ZnnString    = "mu_znn" ;

           EffSfString  += sMbins[i]+sHbins[j]+sBbins[k] ;
           MuSusyString += sMbins[i]+sHbins[j]+sBbins[k] ;
           ZnnString    += sMbins[i]+sHbins[j]+sBbins[k] ;

           printf(" here t2\n") ; cout << flush ;

           dataVal = dataN_0lep[i][j][k];
           eff_sf = 1.0 ;
           susyVal = 0. ;
           znnVal = 0. ;
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
           if ( ws->obj(ZnnString) != 0 ) {
              znnVal = trigeff_sl * ( ((RooRealVar*) ws->obj(ZnnString))  -> getVal() ) ;
           } else {
              printf(" * %s missing.\n", ZnnString.Data() ) ;
              znnVal = 0. ;
           }
           ttwjVal = AttwjVal[i][j][k] ;
           ///qcdVal  = AqcdVal[i][j][k] ;
           char pname[1000] ;
           sprintf( pname, "mu_qcd_M%d_H%d_%db", i+1, j+1, k+1 ) ;
           if ( ws->obj(pname) != 0 ) {
              qcdVal = trigeff * ( ((RooAbsReal*) ws->obj(pname)) -> getVal() ) ;
           } else {
              printf(" * %s missing.  Will try to compute it...\n", pname ) ;

              sprintf( pname, "mu_qcd_ldp_M%d_H%d_%db", i+1, j+1, k+1 ) ;
              RooAbsReal* mu_qcd_ldp = (RooAbsReal*) ws->obj(pname) ;

              sprintf( pname, "sf_qcd_M%d_H%d_%db", i+1, j+1, k+1 ) ;
              RooAbsReal* sf_qcd = (RooAbsReal*) ws->obj(pname) ;

              sprintf( pname, "qcd_0lepLDP_ratio_H%d", j+1 ) ;
              RooAbsReal* r_qcd_ht = (RooAbsReal*) ws->obj(pname) ;

              sprintf( pname, "SFqcd_1stVar%d", i+1 ) ;
              RooAbsReal* sf_qcd_1stVar = (RooAbsReal*) ws->obj(pname) ;

              sprintf( pname, "SFqcd_nb%d", k+1 ) ;
              RooAbsReal* sf_qcd_nb = (RooAbsReal*) ws->obj(pname) ;

              if ( mu_qcd_ldp != 0 && sf_qcd != 0 && r_qcd_ht != 0 && sf_qcd_1stVar != 0 && sf_qcd_nb != 0 ) {
                 qcdVal = trigeff * (mu_qcd_ldp -> getVal())
                        * (sf_qcd -> getVal())
                        * (r_qcd_ht -> getVal())
                        * (sf_qcd_1stVar -> getVal())
                        * (sf_qcd_nb -> getVal()) ;
              } else {
                 printf("\n\n *** missing one of the inputs. mu_qcd_ldp %p : sf_qcd %p : r_qcd_ht %p : sf_qcd_1stVar %p : sf_qcd_nb %p\n",
                     mu_qcd_ldp, sf_qcd, r_qcd_ht, sf_qcd_1stVar, sf_qcd_nb ) ;
              }

           }

           lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

           dataErr = sqrt(dataVal) ;


           cout << "\n" << endl ;
           cout << "binIndex = " << binIndex << endl ;
           cout << binLabel << "       susy : " << susyVal << endl ;
           cout << binLabel << "       ttwj : " << ttwjVal << endl ;
           cout << binLabel << "        qcd : " << qcdVal << endl ;
           cout << binLabel << "        znn : " << znnVal << endl ;
           cout << binLabel << "   LH total : " << lhtotalVal << endl ;
           cout << binLabel << "       data : " << dataVal << endl ;
           cout << flush ;

           if ( doNorm && dataVal > 0. ) {
             dataErr = dataErr / dataVal ;
             susyVal = susyVal / dataVal ;
             ttwjVal = ttwjVal / dataVal ;
             qcdVal  = qcdVal / dataVal ;
             znnVal  = znnVal / dataVal ;
             dataVal = 1. ;
           }

           if ( k == 0 ) {
              xaxis_0lep_1b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_0lep_1b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_0lep_1b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_0lep_1b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_0lep_1b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_0lep_1b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_0lep_1b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_0lep_1b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           } else if ( k == 1 ) {
              xaxis_0lep_2b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_0lep_2b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_0lep_2b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_0lep_2b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_0lep_2b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_0lep_2b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_0lep_2b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_0lep_2b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           } else if ( k == 2 ) {
              xaxis_0lep_3b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_0lep_3b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_0lep_3b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_0lep_3b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_0lep_3b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_0lep_3b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_0lep_3b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_0lep_3b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           } else if ( k == 3 ) {
              xaxis_0lep_4b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_0lep_4b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_0lep_4b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_0lep_4b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_0lep_4b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_0lep_4b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_0lep_4b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_0lep_4b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
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
           qcdVal  = 0. ;
           znnVal  = 0. ;
           lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

           dataErr = sqrt(dataVal) ;


           cout << "\n" << endl ;
           cout << binLabel << "       susy : " << susyVal << endl ;
           cout << binLabel << "       ttwj : " << ttwjVal << endl ;
           cout << binLabel << "        qcd : " << qcdVal << endl ;
           cout << binLabel << "        znn : " << znnVal << endl ;
           cout << binLabel << "   LH total : " << lhtotalVal << endl ;
           cout << binLabel << "       data : " << dataVal << endl ;

           if ( doNorm && dataVal > 0. ) {
             dataErr = dataErr / dataVal ;
             susyVal = susyVal / dataVal ;
             ttwjVal = ttwjVal / dataVal ;
             qcdVal  = qcdVal / dataVal ;
             znnVal  = znnVal / dataVal ;
             dataVal = 1. ;
           }

           if ( k == 0 ) {
              xaxis_1lepSig_1b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lepSig_1b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lepSig_1b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lepSig_1b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lepSig_1b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_1lepSig_1b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_1lepSig_1b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_1lepSig_1b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           } else if ( k == 1 ) {
              xaxis_1lepSig_2b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lepSig_2b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lepSig_2b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lepSig_2b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lepSig_2b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_1lepSig_2b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_1lepSig_2b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_1lepSig_2b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           } else if ( k == 2 ) {
              xaxis_1lepSig_3b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lepSig_3b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lepSig_3b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lepSig_3b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lepSig_3b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_1lepSig_3b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_1lepSig_3b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_1lepSig_3b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           } else if ( k == 3 ) {
              xaxis_1lepSig_4b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lepSig_4b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lepSig_4b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lepSig_4b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lepSig_4b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_1lepSig_4b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_1lepSig_4b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_1lepSig_4b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
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
           qcdVal  = 0. ;
           znnVal  = 0. ;
           lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

           dataErr = sqrt(dataVal) ;


           cout << "\n" << endl ;
           cout << binLabel << "       susy : " << susyVal << endl ;
           cout << binLabel << "       ttwj : " << ttwjVal << endl ;
           cout << binLabel << "        qcd : " << qcdVal << endl ;
           cout << binLabel << "        znn : " << znnVal << endl ;
           cout << binLabel << "   LH total : " << lhtotalVal << endl ;
           cout << binLabel << "       data : " << dataVal << endl ;

           if ( doNorm && dataVal > 0. ) {
             dataErr = dataErr / dataVal ;
             susyVal = susyVal / dataVal ;
             ttwjVal = ttwjVal / dataVal ;
             qcdVal  = qcdVal / dataVal ;
             znnVal  = znnVal / dataVal ;
             dataVal = 1. ;
           }

           if ( k == 0 ) {
              xaxis_1lep_1b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lep_1b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lep_1b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lep_1b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lep_1b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_1lep_1b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_1lep_1b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_1lep_1b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           } else if ( k == 1 ) {
              xaxis_1lep_2b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lep_2b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lep_2b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lep_2b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lep_2b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_1lep_2b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_1lep_2b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_1lep_2b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           } else if ( k == 2 ) {
              xaxis_1lep_3b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lep_3b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lep_3b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lep_3b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lep_3b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_1lep_3b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_1lep_3b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_1lep_3b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           } else if ( k == 3 ) {
              xaxis_1lep_4b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_1lep_4b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_1lep_4b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_1lep_4b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_1lep_4b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_1lep_4b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_1lep_4b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_1lep_4b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           }



           //+++++ ldp histograms:

           binLabel = "ldp" ;
           binLabel += s1bins[i]+s2bins[j]+sBbins[k] ;

           TString EffSfLdpString  = "eff_sf_ldp" ;
           TString MuSusyLdpString = "mu_susy_ldp" ;
           TString MuTtwjLdpString = "mu_ttwj_ldp" ;
           TString MuZnnLdpString = "mu_znn_ldp" ;
           TString MuQcdLdpString  = "mu_qcd_ldp" ;

           EffSfLdpString  += sMbins[i]+sHbins[j]+sBbins[k] ;
           MuSusyLdpString += sMbins[i]+sHbins[j]+sBbins[k] ;
           MuTtwjLdpString += sMbins[i]+sHbins[j]+sBbins[k] ;
           MuZnnLdpString  += sMbins[i]+sHbins[j]+sBbins[k] ;
           MuQcdLdpString  += sMbins[i]+sHbins[j]+sBbins[k] ;

           dataVal = dataN_ldp[i][j][k];
           eff_sf_ldp = 1. ;
           susyVal = 0. ;
           sf_mc = 1.0 ;
           ttwjVal = 0. ;
           znnVal  = 0. ;
           qcdVal  = 0. ;
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
           if ( ws->obj(EffSfLdpString) != 0 ) {
              eff_sf_ldp = ((RooFormulaVar*) ws->obj(EffSfLdpString)) -> getVal() ;
           } else {
              printf(" * %s missing\n", EffSfLdpString.Data() ) ;
              eff_sf_ldp = 1. ;
           }
           if ( ws->obj(MuSusyLdpString) != 0 ) {
              susyVal = trigeff_sl * eff_sf_ldp * ( ((RooRealVar*) ws->obj(MuSusyLdpString)) -> getVal() ) ;
           } else {
              printf(" * %s missing\n", MuSusyLdpString.Data() ) ;
              susyVal = 0. ;
           }
           if ( ws->obj("sf_mc") != 0 ) {
              sf_mc = ((RooFormulaVar*) ws->obj("sf_mc")) -> getVal() ;
           } else {
              printf(" * sf_mc missing\n" ) ;
              sf_mc = 1. ;
           }
           if ( ws->obj(MuTtwjLdpString) != 0 ) {
              ttwjVal = trigeff_sl * eff_sf_ldp * sf_mc * ((((RooRealVar*) ws->obj(MuTtwjLdpString))-> getVal() )   ) ;
           } else {
              printf(" * %s missing\n", MuTtwjLdpString.Data() ) ;
              ttwjVal = 0. ;
           }
           if ( ws->obj(MuZnnLdpString) != 0) {
              znnVal  = trigeff_sl * eff_sf_ldp * sf_mc * ((((RooRealVar*) ws->obj(MuZnnLdpString))-> getVal() )   ) ;
           } else {
              printf(" * %s missing\n", MuZnnLdpString.Data() ) ;
              znnVal = 0. ;
           }
           if ( ws->obj(MuQcdLdpString) != 0 ) {
              qcdVal  = trigeff * ( ((RooRealVar*) ws->obj(MuQcdLdpString))-> getVal() ) ;
           } else {
              printf(" * %s missing\n", MuQcdLdpString.Data() ) ;
              qcdVal = 0. ;
           }

           lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

           dataErr = sqrt(dataVal) ;

           cout << "\n" << endl ;
           cout << binLabel << "       susy : " << susyVal << endl ;
           cout << binLabel << "       ttwj : " << ttwjVal << endl ;
           cout << binLabel << "        qcd : " << qcdVal << endl ;
           cout << binLabel << "        znn : " << znnVal << endl ;
           cout << binLabel << "   LH total : " << lhtotalVal << endl ;
           cout << binLabel << "       data : " << dataVal << endl ;

           printf(" here 1\n") ; cout << flush ;

           if ( doNorm && dataVal > 0. ) {
             dataErr = dataErr / dataVal ;
             susyVal = susyVal / dataVal ;
             ttwjVal = ttwjVal / dataVal ;
             qcdVal  = qcdVal / dataVal ;
             znnVal  = znnVal / dataVal ;
             dataVal = 1. ;
           }

           printf(" here 2\n") ; cout << flush ;
           if ( k == 0 ) {
              xaxis_ldp_1b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_ldp_1b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_ldp_1b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_ldp_1b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_ldp_1b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_ldp_1b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_ldp_1b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_ldp_1b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           } else if ( k == 1 ) {
              xaxis_ldp_2b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_ldp_2b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_ldp_2b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_ldp_2b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_ldp_2b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_ldp_2b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_ldp_2b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_ldp_2b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           } else if ( k == 2 ) {
              xaxis_ldp_3b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_ldp_3b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_ldp_3b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_ldp_3b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_ldp_3b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_ldp_3b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_ldp_3b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_ldp_3b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           } else if ( k == 3 ) {
              xaxis_ldp_4b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_ldp_4b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_ldp_4b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_susy_ldp_4b -> SetBinContent( binIndex, susyVal ) ;
              hfitqual_ttwj_ldp_4b -> SetBinContent( binIndex, ttwjVal ) ;
              hfitqual_qcd_ldp_4b  -> SetBinContent( binIndex, qcdVal ) ;
              hfitqual_znn_ldp_4b  -> SetBinContent( binIndex, znnVal ) ;
              hfitqual_model_ldp_4b  -> SetBinContent( binIndex, susyVal+ttwjVal+qcdVal+znnVal ) ;
           }


           printf(" here 3\n") ; cout << flush ;


           if ( k == 0 ) {

           printf(" here 4\n") ; cout << flush ;

           //+++++ Zee histogram:

              binLabel = "Zee" ;
              binLabel += s1bins[i]+s2bins[j]+sBbins[k] ;
             

              TString nZeeModelString = "n_ee" ;
             
              nZeeModelString += sMbins[i]+sHbins[j] ;

              printf("\n Grabbing value of this from ws: %s\n", nZeeModelString.Data() ) ; cout << flush ;
              lhtotalVal = 0. ;
              if ( ws->obj(nZeeModelString) != 0 ) {
                 lhtotalVal = ((RooRealVar*) ws->obj(nZeeModelString)) -> getVal() ;
              } else {
                 printf(" * %s missing.\n", nZeeModelString.Data() ) ;
                 lhtotalVal = 0. ;
              }
            //lhtotalVal = ((RooRealVar*) ws->obj(nZeeModelString)) -> getVal() ;
              printf("\n Model value for %s is %7.1f\n", nZeeModelString.Data(), lhtotalVal ) ; cout << flush ;
             
              dataVal = dataN_Zee[i][j];
             
              dataErr = sqrt(dataVal) ;
             
              cout << "\n" << endl ;
              cout << binLabel << "      model : " << lhtotalVal << endl ;
              cout << binLabel << "       data : " << dataVal << endl ;
             
              if ( doNorm && dataVal > 0. ) {
                dataErr = dataErr / dataVal ;
        	lhtotalVal = lhtotalVal / dataVal ;
                dataVal = 1. ;
              }
             
              xaxis_zee_1b->SetBinLabel(binIndex, binLabel ) ;
              hfitqual_data_zee_1b -> SetBinContent( binIndex, dataVal ) ;
              hfitqual_data_zee_1b -> SetBinError( binIndex, dataErr ) ;
              hfitqual_fit_zee_1b  -> SetBinContent( binIndex, lhtotalVal ) ;


           printf(" here 5\n") ; cout << flush ;

           //+++++ Zmm histogram:

              binLabel = "Zmm" ;
              binLabel += s1bins[i]+s2bins[j]+sBbins[k] ;
             

              TString nZmmModelString = "n_mm" ;
             
              nZmmModelString += sMbins[i]+sHbins[j] ;

              printf("\n Grabbing value of this from ws: %s\n", nZmmModelString.Data() ) ; cout << flush ;
              lhtotalVal = 0. ;
              if ( ws->obj(nZmmModelString) != 0 ) {
                 lhtotalVal = ((RooRealVar*) ws->obj(nZmmModelString)) -> getVal() ;
              } else {
                 printf(" * %s missing.\n", nZmmModelString.Data() ) ;
                 lhtotalVal = 0. ;
              }
            //lhtotalVal = ((RooRealVar*) ws->obj(nZmmModelString)) -> getVal() ;
              printf("\n Model value for %s is %7.1f\n", nZmmModelString.Data(), lhtotalVal ) ; cout << flush ;
             
	      dataVal = dataN_Zmm[i][j];
             
	      dataErr = sqrt(dataVal) ;
             
	      cout << "\n" << endl ;
	      cout << binLabel << "      model : " << lhtotalVal << endl ;
	      cout << binLabel << "       data : " << dataVal << endl ;
             
	      if ( doNorm && dataVal > 0. ) {
	        dataErr = dataErr / dataVal ;
		lhtotalVal = lhtotalVal / dataVal ;
	        dataVal = 1. ;
	      }
             
	      xaxis_zmm_1b->SetBinLabel(binIndex, binLabel ) ;
	      hfitqual_data_zmm_1b -> SetBinContent( binIndex, dataVal ) ;
	      hfitqual_data_zmm_1b -> SetBinError( binIndex, dataErr ) ;
	      hfitqual_fit_zmm_1b  -> SetBinContent( binIndex, lhtotalVal ) ;

           } // nbjet=1?

	 } // k: btag
       } // j: Var2
     } // i: Var1


           printf(" here 6\n") ; cout << flush ;

     //-- Eff sf ---------

 /// float eff_sf_prim = ((RooRealVar*) ws->obj("eff_sf_prim")) -> getVal() ;

 /// hfitqual_np -> SetBinContent( 1, eff_sf_prim ) ;


     printf("\n\n\n") ;


     //--- final formatting and drawing.


     gStyle->SetOptStat(0) ;
     gStyle->SetTitleX(0.95) ;
     gStyle->SetTitleAlign(33) ;



     hfitqual_fit_0lep_1b->Add( hfitqual_znn_0lep_1b ) ;
     hfitqual_fit_0lep_1b->Add( hfitqual_qcd_0lep_1b ) ;
     hfitqual_fit_0lep_1b->Add( hfitqual_ttwj_0lep_1b ) ;
     hfitqual_fit_0lep_1b->Add( hfitqual_susy_0lep_1b ) ;

     hfitqual_fit_0lep_2b->Add( hfitqual_znn_0lep_2b ) ;
     hfitqual_fit_0lep_2b->Add( hfitqual_qcd_0lep_2b ) ;
     hfitqual_fit_0lep_2b->Add( hfitqual_ttwj_0lep_2b ) ;
     hfitqual_fit_0lep_2b->Add( hfitqual_susy_0lep_2b ) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_fit_0lep_3b->Add( hfitqual_znn_0lep_3b ) ;
       hfitqual_fit_0lep_3b->Add( hfitqual_qcd_0lep_3b ) ;
       hfitqual_fit_0lep_3b->Add( hfitqual_ttwj_0lep_3b ) ;
       hfitqual_fit_0lep_3b->Add( hfitqual_susy_0lep_3b ) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_fit_0lep_4b->Add( hfitqual_znn_0lep_4b ) ;
	 hfitqual_fit_0lep_4b->Add( hfitqual_qcd_0lep_4b ) ;
	 hfitqual_fit_0lep_4b->Add( hfitqual_ttwj_0lep_4b ) ;
	 hfitqual_fit_0lep_4b->Add( hfitqual_susy_0lep_4b ) ;
	 
       }
     }


     hfitqual_fit_1lepSig_1b->Add( hfitqual_znn_1lepSig_1b ) ;
     hfitqual_fit_1lepSig_1b->Add( hfitqual_qcd_1lepSig_1b ) ;
     hfitqual_fit_1lepSig_1b->Add( hfitqual_ttwj_1lepSig_1b ) ;
     hfitqual_fit_1lepSig_1b->Add( hfitqual_susy_1lepSig_1b ) ;

     hfitqual_fit_1lepSig_2b->Add( hfitqual_znn_1lepSig_2b ) ;
     hfitqual_fit_1lepSig_2b->Add( hfitqual_qcd_1lepSig_2b ) ;
     hfitqual_fit_1lepSig_2b->Add( hfitqual_ttwj_1lepSig_2b ) ;
     hfitqual_fit_1lepSig_2b->Add( hfitqual_susy_1lepSig_2b ) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_fit_1lepSig_3b->Add( hfitqual_znn_1lepSig_3b ) ;
       hfitqual_fit_1lepSig_3b->Add( hfitqual_qcd_1lepSig_3b ) ;
       hfitqual_fit_1lepSig_3b->Add( hfitqual_ttwj_1lepSig_3b ) ;
       hfitqual_fit_1lepSig_3b->Add( hfitqual_susy_1lepSig_3b ) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_fit_1lepSig_4b->Add( hfitqual_znn_1lepSig_4b ) ;
	 hfitqual_fit_1lepSig_4b->Add( hfitqual_qcd_1lepSig_4b ) ;
	 hfitqual_fit_1lepSig_4b->Add( hfitqual_ttwj_1lepSig_4b ) ;
	 hfitqual_fit_1lepSig_4b->Add( hfitqual_susy_1lepSig_4b ) ;
	 
       }
     }


     hfitqual_fit_1lep_1b->Add( hfitqual_znn_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_qcd_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_ttwj_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_susy_1lep_1b ) ;

     hfitqual_fit_1lep_2b->Add( hfitqual_znn_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_qcd_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_ttwj_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_susy_1lep_2b ) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_fit_1lep_3b->Add( hfitqual_znn_1lep_3b ) ;
       hfitqual_fit_1lep_3b->Add( hfitqual_qcd_1lep_3b ) ;
       hfitqual_fit_1lep_3b->Add( hfitqual_ttwj_1lep_3b ) ;
       hfitqual_fit_1lep_3b->Add( hfitqual_susy_1lep_3b ) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_fit_1lep_4b->Add( hfitqual_znn_1lep_4b ) ;
	 hfitqual_fit_1lep_4b->Add( hfitqual_qcd_1lep_4b ) ;
	 hfitqual_fit_1lep_4b->Add( hfitqual_ttwj_1lep_4b ) ;
	 hfitqual_fit_1lep_4b->Add( hfitqual_susy_1lep_4b ) ;
	 
       }
     }


     hfitqual_fit_ldp_1b->Add( hfitqual_znn_ldp_1b ) ;
     hfitqual_fit_ldp_1b->Add( hfitqual_qcd_ldp_1b ) ;
     hfitqual_fit_ldp_1b->Add( hfitqual_ttwj_ldp_1b ) ;
     hfitqual_fit_ldp_1b->Add( hfitqual_susy_ldp_1b ) ;

     hfitqual_fit_ldp_2b->Add( hfitqual_znn_ldp_2b ) ;
     hfitqual_fit_ldp_2b->Add( hfitqual_qcd_ldp_2b ) ;
     hfitqual_fit_ldp_2b->Add( hfitqual_ttwj_ldp_2b ) ;
     hfitqual_fit_ldp_2b->Add( hfitqual_susy_ldp_2b ) ;

     if ( nBinsBtag > 2 ) {

       hfitqual_fit_ldp_3b->Add( hfitqual_znn_ldp_3b ) ;
       hfitqual_fit_ldp_3b->Add( hfitqual_qcd_ldp_3b ) ;
       hfitqual_fit_ldp_3b->Add( hfitqual_ttwj_ldp_3b ) ;
       hfitqual_fit_ldp_3b->Add( hfitqual_susy_ldp_3b ) ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_fit_ldp_4b->Add( hfitqual_znn_ldp_4b ) ;
	 hfitqual_fit_ldp_4b->Add( hfitqual_qcd_ldp_4b ) ;
	 hfitqual_fit_ldp_4b->Add( hfitqual_ttwj_ldp_4b ) ;
	 hfitqual_fit_ldp_4b->Add( hfitqual_susy_ldp_4b ) ;
	 
       }
     }


     TLegend* legend = new TLegend(0.6,0.45,0.9,0.95) ;

     legend->AddEntry( hfitqual_data_0lep_1b, "data" ) ;
     legend->AddEntry( hfitqual_susy_0lep_1b, "SUSY" ) ;
     legend->AddEntry( hfitqual_ttwj_0lep_1b, "ttwj" ) ;
     legend->AddEntry( hfitqual_qcd_0lep_1b,  "QCD" ) ;
     legend->AddEntry( hfitqual_znn_0lep_1b,  "Znunu" ) ;
     //legend->AddEntry( hfitqual_np,           "Eff PG" ) ;


     if ( doNorm ) {
       hfitqual_data_0lep_1b->SetMaximum( hmax ) ;
       hfitqual_data_1lepSig_1b->SetMaximum( hmax ) ;
       hfitqual_data_1lep_1b->SetMaximum( hmax ) ;
       hfitqual_data_ldp_1b->SetMaximum( hmax ) ;

       hfitqual_data_0lep_2b->SetMaximum( hmax ) ;
       hfitqual_data_1lepSig_2b->SetMaximum( hmax ) ;
       hfitqual_data_1lep_2b->SetMaximum( hmax ) ;
       hfitqual_data_ldp_2b->SetMaximum( hmax ) ;

       if ( nBinsBtag > 2 ) {
	 
	 hfitqual_data_0lep_3b->SetMaximum( hmax ) ;
	 hfitqual_data_1lepSig_3b->SetMaximum( hmax ) ;
	 hfitqual_data_1lep_3b->SetMaximum( hmax ) ;
	 hfitqual_data_ldp_3b->SetMaximum( hmax ) ;
	 
	 if ( nBinsBtag > 3 ) {
	   
	   hfitqual_data_0lep_4b->SetMaximum( hmax ) ;
	   hfitqual_data_1lepSig_4b->SetMaximum( hmax ) ;
	   hfitqual_data_1lep_4b->SetMaximum( hmax ) ;
	   hfitqual_data_ldp_4b->SetMaximum( hmax ) ;
	   
	 }
       }

     } else {
       hfitqual_data_0lep_1b->SetMaximum( hmax*(hfitqual_data_0lep_1b->GetMaximum()) ) ;
       hfitqual_data_1lepSig_1b->SetMaximum( hmax*(hfitqual_data_1lepSig_1b->GetMaximum()) ) ;
       hfitqual_data_1lep_1b->SetMaximum( hmax*(hfitqual_data_1lep_1b->GetMaximum()) ) ;
       hfitqual_data_ldp_1b->SetMaximum( hmax*(hfitqual_data_ldp_1b->GetMaximum()) ) ;

       hfitqual_data_0lep_2b->SetMaximum( hmax*(hfitqual_data_0lep_2b->GetMaximum()) ) ;
       hfitqual_data_1lepSig_2b->SetMaximum( hmax*(hfitqual_data_1lepSig_2b->GetMaximum()) ) ;
       hfitqual_data_1lep_2b->SetMaximum( hmax*(hfitqual_data_1lep_2b->GetMaximum()) ) ;
       hfitqual_data_ldp_2b->SetMaximum( hmax*(hfitqual_data_ldp_2b->GetMaximum()) ) ;

       if ( nBinsBtag > 2 ) {
	 
	 hfitqual_data_0lep_3b->SetMaximum( hmax*(hfitqual_data_0lep_3b->GetMaximum()) ) ;
	 hfitqual_data_1lepSig_3b->SetMaximum( hmax*(hfitqual_data_1lepSig_3b->GetMaximum()) ) ;
	 hfitqual_data_1lep_3b->SetMaximum( hmax*(hfitqual_data_1lep_3b->GetMaximum()) ) ;
	 hfitqual_data_ldp_3b->SetMaximum( hmax*(hfitqual_data_ldp_3b->GetMaximum()) ) ;
	 
	 if ( nBinsBtag > 3 ) {
	   
	   hfitqual_data_0lep_4b->SetMaximum( hmax*(hfitqual_data_0lep_4b->GetMaximum()) ) ;
	   hfitqual_data_1lepSig_4b->SetMaximum( hmax*(hfitqual_data_1lepSig_4b->GetMaximum()) ) ;
	   hfitqual_data_1lep_4b->SetMaximum( hmax*(hfitqual_data_1lep_4b->GetMaximum()) ) ;
	   hfitqual_data_ldp_4b->SetMaximum( hmax*(hfitqual_data_ldp_4b->GetMaximum()) ) ;
	   
	 }
       }

     }
















   //=======================================================================================================


     printf("\n\n Now doing nuisance parameters.\n\n") ; cout << flush ;

    //----------  Now, do nuisance parameters.

     RooConstVar* npcheck = (RooConstVar*) ws->obj( "mean_sf_qcd_M1_H1_1b" ) ;
     if ( npcheck != 0 ) {

       TH1F *hnp_qcd_1b_val, *hnp_qcd_2b_val, *hnp_qcd_3b_val, *hnp_qcd_4b_val ;
       TH1F *hnp_qcd_1b_pull, *hnp_qcd_2b_pull, *hnp_qcd_3b_pull, *hnp_qcd_4b_pull ;
       TH1F *hnp_qcd_1b_nom, *hnp_qcd_2b_nom, *hnp_qcd_3b_nom, *hnp_qcd_4b_nom ;

       TH1F *hnp_ttwj_1b_val, *hnp_ttwj_2b_val, *hnp_ttwj_3b_val, *hnp_ttwj_4b_val ;
       TH1F *hnp_ttwj_1b_pull, *hnp_ttwj_2b_pull, *hnp_ttwj_3b_pull, *hnp_ttwj_4b_pull ;
       TH1F *hnp_ttwj_1b_nom, *hnp_ttwj_2b_nom, *hnp_ttwj_3b_nom, *hnp_ttwj_4b_nom ;


       hnp_qcd_1b_val  = new TH1F("hnp_qcd_1b_val" , "Nuisance parameters, qcd, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_qcd_1b_pull = new TH1F("hnp_qcd_1b_pull", "Nuisance parameters, qcd, 1b, pull"  , nbins, 0.5, nbins+0.5 ) ;
       hnp_qcd_1b_nom  = new TH1F("hnp_qcd_1b_nom" , "Nuisance parameters, qcd, 1b, nominal", nbins, 0.5, nbins+0.5 ) ;

       hnp_ttwj_1b_val  = new TH1F("hnp_ttwj_1b_val" , "Nuisance parameters, ttwj, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_ttwj_1b_pull = new TH1F("hnp_ttwj_1b_pull", "Nuisance parameters, ttwj, 1b, pull"  , nbins, 0.5, nbins+0.5 ) ;
       hnp_ttwj_1b_nom  = new TH1F("hnp_ttwj_1b_nom" , "Nuisance parameters, ttwj, 1b, nominal", nbins, 0.5, nbins+0.5 ) ;

       hnp_qcd_2b_val  = new TH1F("hnp_qcd_2b_val" , "Nuisance parameters, qcd, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_qcd_2b_pull = new TH1F("hnp_qcd_2b_pull", "Nuisance parameters, qcd, 2b, pull"  , nbins, 0.5, nbins+0.5 ) ;
       hnp_qcd_2b_nom  = new TH1F("hnp_qcd_2b_nom" , "Nuisance parameters, qcd, 2b, nominal", nbins, 0.5, nbins+0.5 ) ;

       hnp_ttwj_2b_val  = new TH1F("hnp_ttwj_2b_val" , "Nuisance parameters, ttwj, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_ttwj_2b_pull = new TH1F("hnp_ttwj_2b_pull", "Nuisance parameters, ttwj, 2b, pull"  , nbins, 0.5, nbins+0.5 ) ;
       hnp_ttwj_2b_nom  = new TH1F("hnp_ttwj_2b_nom" , "Nuisance parameters, ttwj, 2b, nominal", nbins, 0.5, nbins+0.5 ) ;

       if ( nBinsBtag > 2 ) {

	 hnp_qcd_3b_val  = new TH1F("hnp_qcd_3b_val" , "Nuisance parameters, qcd, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_qcd_3b_pull = new TH1F("hnp_qcd_3b_pull", "Nuisance parameters, qcd, 3b, pull"  , nbins, 0.5, nbins+0.5 ) ;
	 hnp_qcd_3b_nom  = new TH1F("hnp_qcd_3b_nom" , "Nuisance parameters, qcd, 3b, nominal", nbins, 0.5, nbins+0.5 ) ;
	 
	 hnp_ttwj_3b_val  = new TH1F("hnp_ttwj_3b_val" , "Nuisance parameters, ttwj, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_ttwj_3b_pull = new TH1F("hnp_ttwj_3b_pull", "Nuisance parameters, ttwj, 3b, pull"  , nbins, 0.5, nbins+0.5 ) ;
	 hnp_ttwj_3b_nom  = new TH1F("hnp_ttwj_3b_nom" , "Nuisance parameters, ttwj, 3b, nominal", nbins, 0.5, nbins+0.5 ) ;
	 
	 if ( nBinsBtag > 3 ) {

	   hnp_qcd_4b_val  = new TH1F("hnp_qcd_4b_val" , "Nuisance parameters, qcd, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_qcd_4b_pull = new TH1F("hnp_qcd_4b_pull", "Nuisance parameters, qcd, 4b, pull"  , nbins, 0.5, nbins+0.5 ) ;
	   hnp_qcd_4b_nom  = new TH1F("hnp_qcd_4b_nom" , "Nuisance parameters, qcd, 4b, nominal", nbins, 0.5, nbins+0.5 ) ;
	   
	   hnp_ttwj_4b_val  = new TH1F("hnp_ttwj_4b_val" , "Nuisance parameters, ttwj, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_ttwj_4b_pull = new TH1F("hnp_ttwj_4b_pull", "Nuisance parameters, ttwj, 4b, pull"  , nbins, 0.5, nbins+0.5 ) ;
	   hnp_ttwj_4b_nom  = new TH1F("hnp_ttwj_4b_nom" , "Nuisance parameters, ttwj, 4b, nominal", nbins, 0.5, nbins+0.5 ) ;
	   
	 }
       }


       TH1F *hnp_eff_sf_1b_val, *hnp_eff_sf_2b_val, *hnp_eff_sf_3b_val, *hnp_eff_sf_4b_val ;
       TH1F *hnp_eff_sf_slSig_1b_val, *hnp_eff_sf_slSig_2b_val, *hnp_eff_sf_slSig_3b_val, *hnp_eff_sf_slSig_4b_val ;
       TH1F *hnp_eff_sf_sl_1b_val, *hnp_eff_sf_sl_2b_val, *hnp_eff_sf_sl_3b_val, *hnp_eff_sf_sl_4b_val ;
       TH1F *hnp_eff_sf_ldp_1b_val, *hnp_eff_sf_ldp_2b_val, *hnp_eff_sf_ldp_3b_val, *hnp_eff_sf_ldp_4b_val ;

       TH1F *hnp_eff_sf_1b_pull, *hnp_eff_sf_2b_pull, *hnp_eff_sf_3b_pull, *hnp_eff_sf_4b_pull ;
       TH1F *hnp_eff_sf_slSig_1b_pull, *hnp_eff_sf_slSig_2b_pull, *hnp_eff_sf_slSig_3b_pull, *hnp_eff_sf_slSig_4b_pull ;
       TH1F *hnp_eff_sf_sl_1b_pull, *hnp_eff_sf_sl_2b_pull, *hnp_eff_sf_sl_3b_pull, *hnp_eff_sf_sl_4b_pull ;
       TH1F *hnp_eff_sf_ldp_1b_pull, *hnp_eff_sf_ldp_2b_pull, *hnp_eff_sf_ldp_3b_pull, *hnp_eff_sf_ldp_4b_pull ;


       hnp_eff_sf_1b_val  = new TH1F("hnp_eff_sf_1b_val" , "Nuisance parameters, eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_1b_pull  = new TH1F("hnp_eff_sf_1b_pull" , "Nuisance parameters, eff SF, 1b, pull", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_slSig_1b_val  = new TH1F("hnp_eff_sf_slSig_1b_val" , "Nuisance parameters, eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_slSig_1b_pull  = new TH1F("hnp_eff_sf_slSig_1b_pull" , "Nuisance parameters, eff SF, 1b, pull", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_sl_1b_val  = new TH1F("hnp_eff_sf_sl_1b_val" , "Nuisance parameters, eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_sl_1b_pull  = new TH1F("hnp_eff_sf_sl_1b_pull" , "Nuisance parameters, eff SF, 1b, pull", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_ldp_1b_val  = new TH1F("hnp_eff_sf_ldp_1b_val" , "Nuisance parameters, eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_ldp_1b_pull  = new TH1F("hnp_eff_sf_ldp_1b_pull" , "Nuisance parameters, eff SF, 1b, pull", nbins, 0.5, nbins+0.5 ) ;

       hnp_eff_sf_2b_val  = new TH1F("hnp_eff_sf_2b_val" , "Nuisance parameters, eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_2b_pull  = new TH1F("hnp_eff_sf_2b_pull" , "Nuisance parameters, eff SF, 2b, pull", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_slSig_2b_val  = new TH1F("hnp_eff_sf_slSig_2b_val" , "Nuisance parameters, eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_slSig_2b_pull  = new TH1F("hnp_eff_sf_slSig_2b_pull" , "Nuisance parameters, eff SF, 2b, pull", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_sl_2b_val  = new TH1F("hnp_eff_sf_sl_2b_val" , "Nuisance parameters, eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_sl_2b_pull  = new TH1F("hnp_eff_sf_sl_2b_pull" , "Nuisance parameters, eff SF, 2b, pull", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_ldp_2b_val  = new TH1F("hnp_eff_sf_ldp_2b_val" , "Nuisance parameters, eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_sf_ldp_2b_pull  = new TH1F("hnp_eff_sf_ldp_2b_pull" , "Nuisance parameters, eff SF, 2b, pull", nbins, 0.5, nbins+0.5 ) ;

       if ( nBinsBtag > 2 ) {

	 hnp_eff_sf_3b_val  = new TH1F("hnp_eff_sf_3b_val" , "Nuisance parameters, eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_eff_sf_3b_pull  = new TH1F("hnp_eff_sf_3b_pull" , "Nuisance parameters, eff SF, 3b, pull", nbins, 0.5, nbins+0.5 ) ;
	 hnp_eff_sf_slSig_3b_val  = new TH1F("hnp_eff_sf_slSig_3b_val" , "Nuisance parameters, eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_eff_sf_slSig_3b_pull  = new TH1F("hnp_eff_sf_slSig_3b_pull" , "Nuisance parameters, eff SF, 3b, pull", nbins, 0.5, nbins+0.5 ) ;
	 hnp_eff_sf_sl_3b_val  = new TH1F("hnp_eff_sf_sl_3b_val" , "Nuisance parameters, eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_eff_sf_sl_3b_pull  = new TH1F("hnp_eff_sf_sl_3b_pull" , "Nuisance parameters, eff SF, 3b, pull", nbins, 0.5, nbins+0.5 ) ;
	 hnp_eff_sf_ldp_3b_val  = new TH1F("hnp_eff_sf_ldp_3b_val" , "Nuisance parameters, eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_eff_sf_ldp_3b_pull  = new TH1F("hnp_eff_sf_ldp_3b_pull" , "Nuisance parameters, eff SF, 3b, pull", nbins, 0.5, nbins+0.5 ) ;
	 
	 if ( nBinsBtag > 3 ) {
	   
	   hnp_eff_sf_4b_val  = new TH1F("hnp_eff_sf_4b_val" , "Nuisance parameters, eff SF, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_eff_sf_4b_pull  = new TH1F("hnp_eff_sf_4b_pull" , "Nuisance parameters, eff SF, 4b, pull", nbins, 0.5, nbins+0.5 ) ;
	   hnp_eff_sf_slSig_4b_val  = new TH1F("hnp_eff_sf_slSig_4b_val" , "Nuisance parameters, eff SF, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_eff_sf_slSig_4b_pull  = new TH1F("hnp_eff_sf_slSig_4b_pull" , "Nuisance parameters, eff SF, 4b, pull", nbins, 0.5, nbins+0.5 ) ;
	   hnp_eff_sf_sl_4b_val  = new TH1F("hnp_eff_sf_sl_4b_val" , "Nuisance parameters, eff SF, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_eff_sf_sl_4b_pull  = new TH1F("hnp_eff_sf_sl_4b_pull" , "Nuisance parameters, eff SF, 4b, pull", nbins, 0.5, nbins+0.5 ) ;
	   hnp_eff_sf_ldp_4b_val  = new TH1F("hnp_eff_sf_ldp_4b_val" , "Nuisance parameters, eff SF, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_eff_sf_ldp_4b_pull  = new TH1F("hnp_eff_sf_ldp_4b_pull" , "Nuisance parameters, eff SF, 4b, pull", nbins, 0.5, nbins+0.5 ) ;
	   
	 }
       }


       TH1F *hnp_btageff_sf_1b_val, *hnp_btageff_sf_2b_val, *hnp_btageff_sf_3b_val, *hnp_btageff_sf_4b_val ;
       TH1F *hnp_btageff_sf_prod_1b_val, *hnp_btageff_sf_prod_2b_val, *hnp_btageff_sf_prod_3b_val, *hnp_btageff_sf_prod_4b_val ;
       TH1F *hnp_btageff_sf_slSig_1b_val, *hnp_btageff_sf_slSig_2b_val, *hnp_btageff_sf_slSig_3b_val, *hnp_btageff_sf_slSig_4b_val ;
       TH1F *hnp_btageff_sf_slSig_prod_1b_val, *hnp_btageff_sf_slSig_prod_2b_val, *hnp_btageff_sf_slSig_prod_3b_val, *hnp_btageff_sf_slSig_prod_4b_val ;
       TH1F *hnp_btageff_sf_sl_1b_val, *hnp_btageff_sf_sl_2b_val, *hnp_btageff_sf_sl_3b_val, *hnp_btageff_sf_sl_4b_val ;
       TH1F *hnp_btageff_sf_sl_prod_1b_val, *hnp_btageff_sf_sl_prod_2b_val, *hnp_btageff_sf_sl_prod_3b_val, *hnp_btageff_sf_sl_prod_4b_val ;
       TH1F *hnp_btageff_sf_ldp_1b_val, *hnp_btageff_sf_ldp_2b_val, *hnp_btageff_sf_ldp_3b_val, *hnp_btageff_sf_ldp_4b_val ;
       TH1F *hnp_btageff_sf_ldp_prod_1b_val, *hnp_btageff_sf_ldp_prod_2b_val, *hnp_btageff_sf_ldp_prod_3b_val, *hnp_btageff_sf_ldp_prod_4b_val ;

       TH1F *hnp_eff_btageff_sf_1b_val, *hnp_eff_btageff_sf_2b_val, *hnp_eff_btageff_sf_3b_val, *hnp_eff_btageff_sf_4b_val ;
       TH1F *hnp_eff_btageff_sf_prod_1b_val, *hnp_eff_btageff_sf_prod_2b_val, *hnp_eff_btageff_sf_prod_3b_val, *hnp_eff_btageff_sf_prod_4b_val ;
       TH1F *hnp_eff_btageff_sf_slSig_1b_val, *hnp_eff_btageff_sf_slSig_2b_val, *hnp_eff_btageff_sf_slSig_3b_val, *hnp_eff_btageff_sf_slSig_4b_val ;
       TH1F *hnp_eff_btageff_sf_slSig_prod_1b_val, *hnp_eff_btageff_sf_slSig_prod_2b_val, *hnp_eff_btageff_sf_slSig_prod_3b_val, *hnp_eff_btageff_sf_slSig_prod_4b_val ;
       TH1F *hnp_eff_btageff_sf_sl_1b_val, *hnp_eff_btageff_sf_sl_2b_val, *hnp_eff_btageff_sf_sl_3b_val, *hnp_eff_btageff_sf_sl_4b_val ;
       TH1F *hnp_eff_btageff_sf_sl_prod_1b_val, *hnp_eff_btageff_sf_sl_prod_2b_val, *hnp_eff_btageff_sf_sl_prod_3b_val, *hnp_eff_btageff_sf_sl_prod_4b_val ;
       TH1F *hnp_eff_btageff_sf_ldp_1b_val, *hnp_eff_btageff_sf_ldp_2b_val, *hnp_eff_btageff_sf_ldp_3b_val, *hnp_eff_btageff_sf_ldp_4b_val ;
       TH1F *hnp_eff_btageff_sf_ldp_prod_1b_val, *hnp_eff_btageff_sf_ldp_prod_2b_val, *hnp_eff_btageff_sf_ldp_prod_3b_val, *hnp_eff_btageff_sf_ldp_prod_4b_val ;


       hnp_btageff_sf_1b_val  = new TH1F("hnp_btageff_sf_1b_val" , "Nuisance parameters, btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_btageff_sf_prod_1b_val  = new TH1F("hnp_eff_btageff_sf_prod_1b_val" , "Nuisance parameters, eff SF * btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_btageff_sf_slSig_1b_val  = new TH1F("hnp_btageff_sf_slSig_1b_val" , "Nuisance parameters, btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_btageff_sf_slSig_prod_1b_val  = new TH1F("hnp_eff_btageff_sf_slSig_prod_1b_val" , "Nuisance parameters, eff SF * btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_btageff_sf_sl_1b_val  = new TH1F("hnp_btageff_sf_sl_1b_val" , "Nuisance parameters, btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_btageff_sf_sl_prod_1b_val  = new TH1F("hnp_eff_btageff_sf_sl_prod_1b_val" , "Nuisance parameters, eff SF * btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_btageff_sf_ldp_1b_val  = new TH1F("hnp_btageff_sf_ldp_1b_val" , "Nuisance parameters, btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_btageff_sf_ldp_prod_1b_val  = new TH1F("hnp_eff_btageff_sf_ldp_prod_1b_val" , "Nuisance parameters, eff SF * btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;

       hnp_btageff_sf_2b_val  = new TH1F("hnp_btageff_sf_2b_val" , "Nuisance parameters, btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_btageff_sf_prod_2b_val  = new TH1F("hnp_eff_btageff_sf_prod_2b_val" , "Nuisance parameters, eff SF * btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_btageff_sf_slSig_2b_val  = new TH1F("hnp_btageff_sf_slSig_2b_val" , "Nuisance parameters, btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_btageff_sf_slSig_prod_2b_val  = new TH1F("hnp_eff_btageff_sf_slSig_prod_2b_val" , "Nuisance parameters, eff SF * btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_btageff_sf_sl_2b_val  = new TH1F("hnp_btageff_sf_sl_2b_val" , "Nuisance parameters, btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_btageff_sf_sl_prod_2b_val  = new TH1F("hnp_eff_btageff_sf_sl_prod_2b_val" , "Nuisance parameters, eff SF * btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_btageff_sf_ldp_2b_val  = new TH1F("hnp_btageff_sf_ldp_2b_val" , "Nuisance parameters, btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
       hnp_eff_btageff_sf_ldp_prod_2b_val  = new TH1F("hnp_eff_btageff_sf_ldp_prod_2b_val" , "Nuisance parameters, eff SF * btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;

       if ( nBinsBtag > 2 ) {

	 hnp_btageff_sf_3b_val  = new TH1F("hnp_btageff_sf_3b_val" , "Nuisance parameters, btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_eff_btageff_sf_prod_3b_val  = new TH1F("hnp_eff_btageff_sf_prod_3b_val" , "Nuisance parameters, eff SF * btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_btageff_sf_slSig_3b_val  = new TH1F("hnp_btageff_sf_slSig_3b_val" , "Nuisance parameters, btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_eff_btageff_sf_slSig_prod_3b_val  = new TH1F("hnp_eff_btageff_sf_slSig_prod_3b_val" , "Nuisance parameters, eff SF * btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_btageff_sf_sl_3b_val  = new TH1F("hnp_btageff_sf_sl_3b_val" , "Nuisance parameters, btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_eff_btageff_sf_sl_prod_3b_val  = new TH1F("hnp_eff_btageff_sf_sl_prod_3b_val" , "Nuisance parameters, eff SF * btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_btageff_sf_ldp_3b_val  = new TH1F("hnp_btageff_sf_ldp_3b_val" , "Nuisance parameters, btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 hnp_eff_btageff_sf_ldp_prod_3b_val  = new TH1F("hnp_eff_btageff_sf_ldp_prod_3b_val" , "Nuisance parameters, eff SF * btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;
	 
	 if ( nBinsBtag > 3 ) {
	   
	   hnp_btageff_sf_4b_val  = new TH1F("hnp_btageff_sf_4b_val" , "Nuisance parameters, btag eff SF, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_eff_btageff_sf_prod_4b_val  = new TH1F("hnp_eff_btageff_sf_prod_4b_val" , "Nuisance parameters, eff SF * btag eff SF, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_btageff_sf_slSig_4b_val  = new TH1F("hnp_btageff_sf_slSig_4b_val" , "Nuisance parameters, btag eff SF, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_eff_btageff_sf_slSig_prod_4b_val  = new TH1F("hnp_eff_btageff_sf_slSig_prod_4b_val" , "Nuisance parameters, eff SF * btag eff SF, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_btageff_sf_sl_4b_val  = new TH1F("hnp_btageff_sf_sl_4b_val" , "Nuisance parameters, btag eff SF, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_eff_btageff_sf_sl_prod_4b_val  = new TH1F("hnp_eff_btageff_sf_sl_prod_4b_val" , "Nuisance parameters, eff SF * btag eff SF, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_btageff_sf_ldp_4b_val  = new TH1F("hnp_btageff_sf_ldp_4b_val" , "Nuisance parameters, btag eff SF, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   hnp_eff_btageff_sf_ldp_prod_4b_val  = new TH1F("hnp_eff_btageff_sf_ldp_prod_4b_val" , "Nuisance parameters, eff SF * btag eff SF, 4b, values", nbins, 0.5, nbins+0.5 ) ;
	   
	 }
       }


        TH1F* hnp_prim_eff = new TH1F("hnp_prim_eff", "Nuisance parameters, Efficiency primary Gaussians, pull", 8, 0.5, 8.5 ) ;

        TH1F* hnp_znn = new TH1F("hnp_znn", "Nuisance parameters, Znn, pull", 3*nBinsVar1+2*nBinsVar2+7, 0.5, 3*nBinsVar1+2*nBinsVar2+7+0.5 ) ;

        TH1F* hnp_trig_0lep = new TH1F("hnp_trig_0lep", "Nuisance parameters, trigger (0lep), pull", 2+nBinsVar1*nBinsVar2-1, 0.5, 2+nBinsVar1*nBinsVar2-1+0.5 ) ;
        TH1F* hnp_trig_1lep = new TH1F("hnp_trig_1lep", "Nuisance parameters, trigger (1lep), pull", 2+nBinsVar1*nBinsVar2-1, 0.5, 2+nBinsVar1*nBinsVar2-1+0.5 ) ;




        for ( int mbi=0 ; mbi < nBinsVar1; mbi ++ ) {
           for ( int hbi=0 ; hbi < nBinsVar2; hbi ++ ) {
              if ( ignoreBin[mbi][hbi] ) continue ;
              for ( int bbi=0 ; bbi < nBinsBtag; bbi ++ ) {

                 if ( dataN_0lep[mbi][hbi][bbi] == 0 ) continue ; //-- to do something sensible for blind fit.

                 char parName[1000] ;
                 sprintf( parName, "sf_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 RooAbsReal* np = (RooAbsReal*) ws->obj( parName ) ;
                 double npVal = 0.0 ;
                 if ( np == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
                 } else {
                    npVal = np->getVal() ;
                 }
                 char vname[1000] ;
                 sprintf( vname, "mean_%s", parName ) ;
                 RooConstVar* rcv_mean = (RooConstVar*) ws->obj( vname ) ;
                 double mean = 0. ;
                 if ( rcv_mean == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter mean? %s\n\n", vname ) ;
                 } else {
                    mean = rcv_mean->getVal() ;
                 }
                 sprintf( vname, "sigma_%s", parName ) ;
                 RooConstVar* rcv_sigma = (RooConstVar*) ws->obj( vname ) ;
                 double sigma = 1. ;
                 if ( rcv_sigma == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter sigma? %s\n\n", vname ) ;
                 } else {
                    sigma = rcv_sigma->getVal() ;
                 }

                 double pull = (npVal-mean)/sigma ;
                 printf("  %s : mean=%7.3f, sigma=%7.3f, val=%7.3f, pull=%7.3f\n", parName, mean, sigma, npVal, pull ) ;
                 binIndex = 1 + (nBinsVar2+1)*mbi + hbi + 1 ;
                 if ( bbi==0 ) hnp_qcd_1b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==1 ) hnp_qcd_2b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==2 ) hnp_qcd_3b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==3 ) hnp_qcd_4b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==0 ) hnp_qcd_1b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==1 ) hnp_qcd_2b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==2 ) hnp_qcd_3b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==3 ) hnp_qcd_4b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==0 ) hnp_qcd_1b_nom->SetBinContent( binIndex, mean ) ;
                 if ( bbi==1 ) hnp_qcd_2b_nom->SetBinContent( binIndex, mean ) ;
                 if ( bbi==2 ) hnp_qcd_3b_nom->SetBinContent( binIndex, mean ) ;
                 if ( bbi==3 ) hnp_qcd_4b_nom->SetBinContent( binIndex, mean ) ;

                 char label[1000] ;
                 sprintf( label, "#alpha%d_#beta%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 if ( bbi==0 ) hnp_qcd_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_qcd_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_qcd_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_qcd_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_qcd_1b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_qcd_2b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_qcd_3b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_qcd_4b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_qcd_1b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_qcd_2b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_qcd_3b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_qcd_4b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;


              }
           } // hbi
        } // mbi.
        printf("\n\n") ;


        for ( int mbi=0 ; mbi < nBinsVar1; mbi ++ ) {
           for ( int hbi=0 ; hbi < nBinsVar2; hbi ++ ) {
              if ( ignoreBin[mbi][hbi] ) continue ;
              for ( int bbi=0 ; bbi < nBinsBtag; bbi ++ ) {

                 if ( dataN_0lep[mbi][hbi][bbi] == 0 ) continue ; //-- to do something sensible for blind fit.

                 char parName[1000] ;
                 sprintf( parName, "sf_ttwj_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 RooAbsReal* np = (RooAbsReal*) ws->obj( parName ) ;
                 double npVal = 0.0 ;
                 if ( np == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
                 } else {
                    npVal = np->getVal() ;
                 }
                 char vname[1000] ;
                 sprintf( vname, "mean_%s", parName ) ;
                 RooConstVar* rcv_mean = (RooConstVar*) ws->obj( vname ) ;
                 double mean = 0. ;
                 if ( rcv_mean == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter mean? %s\n\n", vname ) ;
                 } else {
                    mean = rcv_mean->getVal() ;
                 }
                 sprintf( vname, "sigma_%s", parName ) ;
                 RooConstVar* rcv_sigma = (RooConstVar*) ws->obj( vname ) ;
                 double sigma = 1. ;
                 if ( rcv_sigma == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter sigma? %s\n\n", vname ) ;
                 } else {
                    sigma = rcv_sigma->getVal() ;
                 }

                 double pull = (npVal-mean)/sigma ;
                 printf("  %s : mean=%7.3f, sigma=%7.3f, val=%7.3f, pull=%7.3f\n", parName, mean, sigma, npVal, pull ) ;
                 binIndex = 1 + (nBinsVar2+1)*mbi + hbi + 1 ;
                 if ( bbi==0 ) hnp_ttwj_1b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==1 ) hnp_ttwj_2b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==2 ) hnp_ttwj_3b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==3 ) hnp_ttwj_4b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==0 ) hnp_ttwj_1b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==1 ) hnp_ttwj_2b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==2 ) hnp_ttwj_3b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==3 ) hnp_ttwj_4b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==0 ) hnp_ttwj_1b_nom->SetBinContent( binIndex, mean ) ;
                 if ( bbi==1 ) hnp_ttwj_2b_nom->SetBinContent( binIndex, mean ) ;
                 if ( bbi==2 ) hnp_ttwj_3b_nom->SetBinContent( binIndex, mean ) ;
                 if ( bbi==3 ) hnp_ttwj_4b_nom->SetBinContent( binIndex, mean ) ;

                 char label[1000] ;
                 sprintf( label, "#alpha%d_#beta%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 if ( bbi==0 ) hnp_ttwj_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_ttwj_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_ttwj_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_ttwj_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_ttwj_1b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_ttwj_2b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_ttwj_3b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_ttwj_4b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_ttwj_1b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_ttwj_2b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_ttwj_3b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_ttwj_4b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;

              }
           } // hbi
        } // mbi.
        printf("\n\n") ;


        //--- General parameters
        {
            char parName[1000] ;
            RooAbsReal* np(0x0) ;

        /////////// obsolete //////////////////////////
        /// sprintf( parName, "eff_sf" ) ;
        /// np = (RooAbsReal*) ws->obj( parName ) ;
        /// double global_eff_sf(0.) ;
        /// if ( np == 0x0 ) {
        ///    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
        ///    /// return ;
        /// } else {
        ///    global_eff_sf = np -> getVal() ;
        /// }
        /////////// obsolete //////////////////////////






	    //--- QCD SF_1stVar3

            sprintf( parName, "SFqcd_1stVar3" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double SFqcd_1stVar3_val = 0. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               SFqcd_1stVar3_val =  np -> getVal() ;
            }

            sprintf( parName, "pdf_mean_SFqcd_1stVar3" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double SFqcd_1stVar3_mean = 0. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               SFqcd_1stVar3_mean =  np -> getVal() ;
            }

            sprintf( parName, "pdf_sigma_SFqcd_1stVar3" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double SFqcd_1stVar3_sigma = 1. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               SFqcd_1stVar3_sigma =  np -> getVal() ;
            }
            double SFqcd_1stVar3_pull = (SFqcd_1stVar3_val-SFqcd_1stVar3_mean) / SFqcd_1stVar3_sigma ;

            printf("\n SFqcd_1stVar3 : val = %6.3f, mean = %6.3f, sigma = %6.3f, pull = %6.3f\n",
                 SFqcd_1stVar3_val, SFqcd_1stVar3_mean, SFqcd_1stVar3_sigma, SFqcd_1stVar3_pull ) ;




           //--- QCD SF_1stVar4

            sprintf( parName, "SFqcd_1stVar4" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double SFqcd_1stVar4_val = 0. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               SFqcd_1stVar4_val =  np -> getVal() ;
            }

            sprintf( parName, "pdf_mean_SFqcd_1stVar4" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double SFqcd_1stVar4_mean = 0. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               SFqcd_1stVar4_mean =  np -> getVal() ;
            }

            sprintf( parName, "pdf_sigma_SFqcd_1stVar4" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double SFqcd_1stVar4_sigma = 1. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               SFqcd_1stVar4_sigma =  np -> getVal() ;
            }
            double SFqcd_1stVar4_pull = (SFqcd_1stVar4_val-SFqcd_1stVar4_mean) / SFqcd_1stVar4_sigma ;

            printf("\n SFqcd_1stVar4 : val = %6.3f, mean = %6.3f, sigma = %6.3f, pull = %6.3f\n",
                 SFqcd_1stVar4_val, SFqcd_1stVar4_mean, SFqcd_1stVar4_sigma, SFqcd_1stVar4_pull ) ;






           //--- QCD SF_nb3

            sprintf( parName, "SFqcd_nb3" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double SFqcd_nb3_val = 0. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               SFqcd_nb3_val =  np -> getVal() ;
            }

            sprintf( parName, "pdf_mean_SFqcd_nb3" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double SFqcd_nb3_mean = 0. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               SFqcd_nb3_mean =  np -> getVal() ;
            }

            sprintf( parName, "pdf_sigma_SFqcd_nb3" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double SFqcd_nb3_sigma = 1. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               SFqcd_nb3_sigma =  np -> getVal() ;
            }
            double SFqcd_nb3_pull = (SFqcd_nb3_val-SFqcd_nb3_mean) / SFqcd_nb3_sigma ;

            printf("\n SFqcd_nb3 : val = %6.3f, mean = %6.3f, sigma = %6.3f, pull = %6.3f\n",
                 SFqcd_nb3_val, SFqcd_nb3_mean, SFqcd_nb3_sigma, SFqcd_nb3_pull ) ;





           //--- QCD SF_nb4

            sprintf( parName, "SFqcd_nb4" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double SFqcd_nb4_val = 0. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               SFqcd_nb4_val =  np -> getVal() ;
            }

            sprintf( parName, "pdf_mean_SFqcd_nb4" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double SFqcd_nb4_mean = 0. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               SFqcd_nb4_mean =  np -> getVal() ;
            }

            sprintf( parName, "pdf_sigma_SFqcd_nb4" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double SFqcd_nb4_sigma = 1. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               SFqcd_nb4_sigma =  np -> getVal() ;
            }
            double SFqcd_nb4_pull = (SFqcd_nb4_val-SFqcd_nb4_mean) / SFqcd_nb4_sigma ;

            printf("\n SFqcd_nb4 : val = %6.3f, mean = %6.3f, sigma = %6.3f, pull = %6.3f\n",
                 SFqcd_nb4_val, SFqcd_nb4_mean, SFqcd_nb4_sigma, SFqcd_nb4_pull ) ;







           //--- VV SF

            sprintf( parName, "rar_vv_sf" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double vv_sf_val = 0. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               vv_sf_val =  np -> getVal() ;
            }

            sprintf( parName, "mean_rar_vv_sf" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double vv_sf_mean = 0. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               vv_sf_mean =  np -> getVal() ;
            }

            sprintf( parName, "sigma_rar_vv_sf" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double vv_sf_sigma = 1. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               vv_sf_sigma =  np -> getVal() ;
            }
            double vv_sf_pull = (vv_sf_val-vv_sf_mean) / vv_sf_sigma ;

            printf("\n vv_sf : val = %6.3f, mean = %6.3f, sigma = %6.3f, pull = %6.3f\n",
                 vv_sf_val, vv_sf_mean, vv_sf_sigma, vv_sf_pull ) ;







          //-- SF MC

            sprintf( parName, "sf_mc" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double global_sf_mc = 1. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               global_sf_mc = np -> getVal() ;
            }
            sprintf( parName, "sigma_sf_mc" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double sf_mc_sigma = 1. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               sf_mc_sigma = np -> getVal() ;
            }
            double sf_mc_pull = (global_sf_mc-1.) / sf_mc_sigma ;
            printf(" Overall MC systematic val = %5.2f, sigma = %5.2f, pull = %6.3f\n", global_sf_mc, sf_mc_sigma, sf_mc_pull ) ;





            sprintf( parName, "btageff_sf" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double global_btageff_sf = 0. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               global_btageff_sf = np -> getVal() ;
            }
            printf(" Global efficiency scale factors (Gaussian with mean 0, sigma 1): btageff_sf = %6.3f\n", global_btageff_sf ) ;



            TAxis* xaxis = hnp_prim_eff -> GetXaxis() ;

            int hbin = 2 ;


            hnp_prim_eff -> SetBinContent( hbin, SFqcd_1stVar3_pull ) ;
            xaxis->SetBinLabel(hbin,"SF QCD, 1stVar3") ;
            hbin++ ;

            hnp_prim_eff -> SetBinContent( hbin, SFqcd_1stVar4_pull ) ;
            xaxis->SetBinLabel(hbin,"SF QCD, 1stVar4") ;
            hbin++ ;

            hnp_prim_eff -> SetBinContent( hbin, SFqcd_nb3_pull ) ;
            xaxis->SetBinLabel(hbin,"SF QCD, nB>=3") ;
            hbin++ ;

            hnp_prim_eff -> SetBinContent( hbin, SFqcd_nb4_pull ) ;
            xaxis->SetBinLabel(hbin,"SF QCD, nB>=4") ;
            hbin++ ;

            hnp_prim_eff -> SetBinContent( hbin, vv_sf_pull ) ;
            xaxis->SetBinLabel(hbin,"SF VV") ;
            hbin++ ;

            hnp_prim_eff -> SetBinContent(hbin, global_btageff_sf ) ;
            xaxis->SetBinLabel(hbin,"Btag Eff SF") ;
            hbin++ ;

            hnp_prim_eff -> SetBinContent(hbin, sf_mc_pull ) ;
            xaxis->SetBinLabel(hbin,"MC SF") ;
            hbin++ ;

            hnp_prim_eff -> SetFillColor(kOrange+1 ) ;

            hnp_prim_eff -> SetLabelSize(0.055,"x") ;
            hnp_prim_eff -> GetXaxis() ->LabelsOption("v") ;

            hnp_prim_eff->SetMinimum(-2.0) ;
            hnp_prim_eff->SetMaximum(2.0) ;

        }
        printf("\n\n") ;


        //--- Znn nuisance parameters.
        {

            TH1F* hp = hnp_znn ;

            hp -> SetFillColor( kGreen-3 ) ;

            TAxis* xaxis = hp -> GetXaxis() ;

            int  n_znnNP_gauss_pars(7) ;
            char znnNP_gauss_par[7][100] = { "knn_1b_M1", "knn_1b_M2", "knn_1b_M3", "knn_1b_M4",
                                             "knn_2b", "knn_3b",
                                             "sf_ll"
                                             } ;


            for ( int npi=0; npi< n_znnNP_gauss_pars; npi++ ) {

               char parName[1000] ;
               char vname[1000] ;
               RooAbsReal* np(0x0) ;
               RooAbsReal* np_mean(0x0) ;
               RooAbsReal* np_sigma(0x0) ;

               int hbin = npi+2 ;

               sprintf( parName, "%s", znnNP_gauss_par[npi] ) ;
               np = (RooAbsReal*) ws->obj( parName ) ;
               double val = 1. ;
               if ( np == 0x0 ) {
                  printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
               } else {
                 val = np->getVal() ;
               }

               sprintf( vname, "mean_%s", parName ) ;
               np_mean = (RooAbsReal*) ws->obj( vname ) ;
               double mean = 1. ;
               if ( np_mean == 0x0 ) {
                  printf("\n\n *** missing nuisance parameter mean? %s\n\n", vname ) ;
               } else {
                  mean = np_mean->getVal() ;
               }

               sprintf( vname, "sigma_%s", parName ) ;
               np_sigma = (RooAbsReal*) ws->obj( vname ) ;
               double sigma = 1. ;
               if ( np_sigma == 0x0 ) {
                  printf("\n\n *** missing nuisance parameter sigma? %s\n\n", vname ) ;
               } else {
                  sigma = np_sigma->getVal() ;
               }

               double pull = (val-mean)/sigma ;

               hp -> SetBinContent( hbin, pull ) ;
               xaxis -> SetBinLabel( hbin, parName ) ;
               printf(" Znn NP : %s = %g, mean = %g, sigma = %g, pull = %5.2f\n", parName, val, mean, sigma, pull ) ;

            } // npi

            char znnNP_beta_par[20][100] ;
            int zbnpi(0) ;

	    sprintf( znnNP_beta_par[zbnpi++], "eff_Zee" ) ;
	    sprintf( znnNP_beta_par[zbnpi++], "eff_Zmm" ) ;

            for ( int mbi=0; mbi<nBinsVar1; mbi++ ) {
               sprintf( znnNP_beta_par[zbnpi++], "acc_Zee_M%d", mbi+1 ) ;
               sprintf( znnNP_beta_par[zbnpi++], "acc_Zmm_M%d", mbi+1 ) ;
            } // i
            sprintf( znnNP_beta_par[zbnpi++], "pur_Zee" ) ;
            sprintf( znnNP_beta_par[zbnpi++], "pur_Zmm" ) ;
            int  n_znnNP_beta_pars = zbnpi ;

            for ( int npi=0; npi< n_znnNP_beta_pars; npi++ ) {

               int hbin = 2 + n_znnNP_gauss_pars + npi ;

               char parName[1000] ;
               sprintf( parName, "%s", znnNP_beta_par[npi] ) ;

               double mode, rms, alpha, beta ;
               getBetaModeRMS( parName, ws, mode, rms, alpha, beta ) ;
               RooAbsReal* np = (RooAbsReal*) ws->obj( parName ) ;
               if ( np == 0x0 ) {
                  printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
               }
               double npVal = 1. ;
               double pull = 0. ;
               if ( np != 0x0 ) {
                  npVal = np->getVal() ;
                  pull = (npVal-mode)/rms ;
               }

               hp -> SetBinContent( hbin, pull ) ;
               xaxis -> SetBinLabel( hbin, parName ) ;

               printf("  %s : mode=%7.3f, rms=%7.3f, val=%7.3f, pull=%7.3f\n", parName, mode, rms, npVal, pull ) ;

            } // npi

            hp -> SetLabelSize(0.055,"x") ;
            hp -> GetXaxis() ->LabelsOption("v") ;

            hp->SetMinimum(-2.0) ;
            hp->SetMaximum(2.0) ;

        }
        printf("\n\n") ;



       //--- Trigger, 0lep

        {

           TH1F* hp = hnp_trig_0lep ;

           hp -> SetFillColor( kOrange-3 ) ;

           TAxis* xaxis = hp -> GetXaxis() ;

           int hbin = 2 ;

           for ( int mbi=0 ; mbi < nBinsVar1; mbi ++ ) {
              for ( int hbi=0 ; hbi < nBinsVar2; hbi ++ ) {

                 if ( ignoreBin[mbi][hbi] ) continue ;

                 char parName[1000] ;

                 sprintf( parName, "trigeff_M%d_H%d", mbi+1, hbi+1 ) ;
                 double mode, rms, alpha, beta ;
                 getBetaModeRMS( parName, ws, mode, rms, alpha, beta ) ;
                 RooAbsReal* np = (RooAbsReal*) ws->obj( parName ) ;
                 if ( np == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
                 }
                 double npVal = 1. ;
                 double pull = 0. ;
                 if ( np != 0x0 ) {
                    npVal = np->getVal() ;
                    pull = (npVal-mode)/rms ;
                 }

                 hp -> SetBinContent( hbin, pull ) ;
                 xaxis -> SetBinLabel( hbin, parName ) ;

                 printf("  %s : mode=%7.3f, rms=%7.3f, val=%7.3f, pull=%7.3f\n", parName, mode, rms, npVal, pull ) ;

                 hbin++ ;

              } // hbi.
           } // mbi.

           hp -> SetLabelSize(0.055,"x") ;
           hp -> GetXaxis() ->LabelsOption("v") ;

           hp->SetMinimum(-2.0) ;
           hp->SetMaximum(2.0) ;

        }


       //--- Trigger, 1lep

        {

           TH1F* hp = hnp_trig_1lep ;

           hp -> SetFillColor( kOrange-3 ) ;

           TAxis* xaxis = hp -> GetXaxis() ;

           int hbin = 2 ;

           for ( int mbi=0 ; mbi < nBinsVar1; mbi ++ ) {
              for ( int hbi=0 ; hbi < nBinsVar2; hbi ++ ) {

                 if ( ignoreBin[mbi][hbi] ) continue ;

                 char parName[1000] ;

                 sprintf( parName, "trigeff_sl_M%d_H%d", mbi+1, hbi+1 ) ;
                 double mode, rms, alpha, beta ;
                 getBetaModeRMS( parName, ws, mode, rms, alpha, beta ) ;
                 RooAbsReal* np = (RooAbsReal*) ws->obj( parName ) ;
                 if ( np == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
                 }
                 double npVal = 1. ;
                 double pull = 0. ;
                 if ( np != 0x0 ) {
                    npVal = np->getVal() ;
                    pull = (npVal-mode)/rms ;
                 }

                 hp -> SetBinContent( hbin, pull ) ;
                 xaxis -> SetBinLabel( hbin, parName ) ;

                 printf("  %s : mode=%7.3f, rms=%7.3f, val=%7.3f, pull=%7.3f\n", parName, mode, rms, npVal, pull ) ;

                 hbin++ ;

              } // hbi.
           } // mbi.

           hp -> SetLabelSize(0.055,"x") ;
           hp -> GetXaxis() ->LabelsOption("v") ;

           hp->SetMinimum(-2.0) ;
           hp->SetMaximum(2.0) ;

        }


       //--- bin-by-bin efficiencies


        for ( int mbi=0 ; mbi < nBinsVar1; mbi ++ ) {
           for ( int hbi=0 ; hbi < nBinsVar2; hbi ++ ) {

              if ( ignoreBin[mbi][hbi] ) continue ;

              for ( int bbi=0 ; bbi < nBinsBtag; bbi ++ ) {

                 char vname[1000] ;
                 char parName[1000] ;
                 RooAbsReal* np(0x0) ;
                 RooAbsReal* np_mean(0x0) ;
                 RooAbsReal* np_sigma(0x0) ;

                 double val, mean, sigma ;


                //-- eff_sf

                 sprintf( parName, "eff_sf_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 np = (RooAbsReal*) ws->obj( parName ) ;
                 double bin_eff_sf = 1. ;
                 if ( np == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
                 } else {
                    bin_eff_sf = np -> getVal() ;
                    val = np->getVal() ;
                 }
                 double bin_eff_sf_pull = 0. ;

                 sprintf( vname, "mean_%s", parName ) ;
                 np_mean = (RooAbsReal*) ws->obj( vname ) ;
                 if ( np_mean == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter mean? %s\n\n", vname ) ;
                 } else {
                    mean = np_mean->getVal() ;

                    sprintf( vname, "sigma_%s", parName ) ;
                    np_sigma = (RooAbsReal*) ws->obj( vname ) ;
                    sigma = 1. ;
                    if ( np_sigma == 0x0 ) {
                       printf("\n\n *** missing nuisance parameter sigma? %s\n\n", vname ) ;
                    } else {
                       sigma = np_sigma->getVal() ;
                    }

                    bin_eff_sf_pull = (val-mean)/sigma ;
                 }




                //-- eff_sf_slSig

                 sprintf( parName, "eff_sf_slSig_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 np = (RooAbsReal*) ws->obj( parName ) ;
                 double bin_eff_sf_slSig = 1. ;
                 if ( np == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
                 } else {
                    bin_eff_sf_slSig = np -> getVal() ;
                    val = np->getVal() ;
                 }
                 double bin_eff_sf_slSig_pull = 0. ;

                 sprintf( vname, "mean_%s", parName ) ;
                 np_mean = (RooAbsReal*) ws->obj( vname ) ;
                 if ( np_mean == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter mean? %s\n\n", vname ) ;
                 } else {
                    mean = np_mean->getVal() ;

                    sprintf( vname, "sigma_%s", parName ) ;
                    np_sigma = (RooAbsReal*) ws->obj( vname ) ;
                    sigma = 1. ;
                    if ( np_sigma == 0x0 ) {
                       printf("\n\n *** missing nuisance parameter sigma? %s\n\n", vname ) ;
                    } else {
                       sigma = np_sigma->getVal() ;
                    }

                    bin_eff_sf_slSig_pull = (val-mean)/sigma ;
                 }




                //-- eff_sf_sl

                 sprintf( parName, "eff_sf_sl_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 np = (RooAbsReal*) ws->obj( parName ) ;
                 double bin_eff_sf_sl = 1. ;
                 if ( np == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
                 } else {
                    bin_eff_sf_sl = np -> getVal() ;
                    val = np->getVal() ;
                 }
                 double bin_eff_sf_sl_pull = 0. ;

                 sprintf( vname, "mean_%s", parName ) ;
                 np_mean = (RooAbsReal*) ws->obj( vname ) ;
                 if ( np_mean == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter mean? %s\n\n", vname ) ;
                 } else {
                    mean = np_mean->getVal() ;

                    sprintf( vname, "sigma_%s", parName ) ;
                    np_sigma = (RooAbsReal*) ws->obj( vname ) ;
                    sigma = 1. ;
                    if ( np_sigma == 0x0 ) {
                       printf("\n\n *** missing nuisance parameter sigma? %s\n\n", vname ) ;
                    } else {
                       sigma = np_sigma->getVal() ;
                    }

                    bin_eff_sf_sl_pull = (val-mean)/sigma ;
                 }




                //-- eff_sf_ldp

                 sprintf( parName, "eff_sf_ldp_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 np = (RooAbsReal*) ws->obj( parName ) ;
                 double bin_eff_sf_ldp = 1. ;
                 if ( np == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
                 } else {
                    bin_eff_sf_ldp = np -> getVal() ;
                    val = np->getVal() ;
                 }
                 double bin_eff_sf_ldp_pull = 0. ;

                 sprintf( vname, "mean_%s", parName ) ;
                 np_mean = (RooAbsReal*) ws->obj( vname ) ;
                 if ( np_mean == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter mean? %s\n\n", vname ) ;
                 } else {
                    mean = np_mean->getVal() ;

                    sprintf( vname, "sigma_%s", parName ) ;
                    np_sigma = (RooAbsReal*) ws->obj( vname ) ;
                    sigma = 1. ;
                    if ( np_sigma == 0x0 ) {
                       printf("\n\n *** missing nuisance parameter sigma? %s\n\n", vname ) ;
                    } else {
                       sigma = np_sigma->getVal() ;
                    }

                    bin_eff_sf_ldp_pull = (val-mean)/sigma ;
                 }



                 sprintf( parName, "btageff_sf_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 np = (RooAbsReal*) ws->obj( parName ) ;
                 double btageff_sf = 1. ;
                 if ( np == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
                 } else {
                    btageff_sf = np -> getVal() ;
                 }

                 double prod = bin_eff_sf * btageff_sf ;



                 sprintf( parName, "btageff_sf_slSig_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 np = (RooAbsReal*) ws->obj( parName ) ;
                 double btageff_sf_slSig = 1. ;
                 if ( np == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
                 } else {
                    btageff_sf_slSig = np -> getVal() ;
                 }

                 double prod_slSig = bin_eff_sf_slSig * btageff_sf_slSig ;



                 sprintf( parName, "btageff_sf_sl_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 np = (RooAbsReal*) ws->obj( parName ) ;
                 double btageff_sf_sl = 1. ;
                 if ( np == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
                 } else {
                    btageff_sf_sl = np -> getVal() ;
                 }

                 double prod_sl = bin_eff_sf_sl * btageff_sf_sl ;



                 sprintf( parName, "btageff_sf_ldp_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 np = (RooAbsReal*) ws->obj( parName ) ;
                 double btageff_sf_ldp = 1. ;
                 if ( np == 0x0 ) {
                    printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
                 } else {
                    btageff_sf_ldp = np -> getVal() ;
                 }

                 double prod_ldp = bin_eff_sf_ldp * btageff_sf_ldp ;




                 binIndex = 1 + (nBinsVar2+1)*mbi + hbi + 1 ;

                 if ( bbi == 0 ) { hnp_eff_sf_1b_val -> SetBinContent( binIndex, bin_eff_sf ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_2b_val -> SetBinContent( binIndex, bin_eff_sf ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_3b_val -> SetBinContent( binIndex, bin_eff_sf ) ; }
                 if ( bbi == 3 ) { hnp_eff_sf_4b_val -> SetBinContent( binIndex, bin_eff_sf ) ; }
                 if ( bbi == 0 ) { hnp_eff_sf_1b_pull -> SetBinContent( binIndex, bin_eff_sf_pull ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_2b_pull -> SetBinContent( binIndex, bin_eff_sf_pull ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_3b_pull -> SetBinContent( binIndex, bin_eff_sf_pull ) ; }
                 if ( bbi == 3 ) { hnp_eff_sf_4b_pull -> SetBinContent( binIndex, bin_eff_sf_pull ) ; }
                 if ( bbi == 0 ) { hnp_btageff_sf_1b_val -> SetBinContent( binIndex, btageff_sf ) ; }
                 if ( bbi == 1 ) { hnp_btageff_sf_2b_val -> SetBinContent( binIndex, btageff_sf ) ; }
                 if ( bbi == 2 ) { hnp_btageff_sf_3b_val -> SetBinContent( binIndex, btageff_sf ) ; }
                 if ( bbi == 3 ) { hnp_btageff_sf_4b_val -> SetBinContent( binIndex, btageff_sf ) ; }
                 if ( bbi == 0 ) { hnp_eff_btageff_sf_prod_1b_val -> SetBinContent( binIndex, prod ) ; }
                 if ( bbi == 1 ) { hnp_eff_btageff_sf_prod_2b_val -> SetBinContent( binIndex, prod ) ; }
                 if ( bbi == 2 ) { hnp_eff_btageff_sf_prod_3b_val -> SetBinContent( binIndex, prod ) ; }
                 if ( bbi == 3 ) { hnp_eff_btageff_sf_prod_4b_val -> SetBinContent( binIndex, prod ) ; }

                 if ( bbi == 0 ) { hnp_eff_sf_slSig_1b_val -> SetBinContent( binIndex, bin_eff_sf_slSig ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_slSig_2b_val -> SetBinContent( binIndex, bin_eff_sf_slSig ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_slSig_3b_val -> SetBinContent( binIndex, bin_eff_sf_slSig ) ; }
                 if ( bbi == 3 ) { hnp_eff_sf_slSig_4b_val -> SetBinContent( binIndex, bin_eff_sf_slSig ) ; }
                 if ( bbi == 0 ) { hnp_eff_sf_slSig_1b_pull -> SetBinContent( binIndex, bin_eff_sf_slSig_pull ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_slSig_2b_pull -> SetBinContent( binIndex, bin_eff_sf_slSig_pull ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_slSig_3b_pull -> SetBinContent( binIndex, bin_eff_sf_slSig_pull ) ; }
                 if ( bbi == 3 ) { hnp_eff_sf_slSig_4b_pull -> SetBinContent( binIndex, bin_eff_sf_slSig_pull ) ; }
                 if ( bbi == 0 ) { hnp_btageff_sf_slSig_1b_val -> SetBinContent( binIndex, btageff_sf_slSig ) ; }
                 if ( bbi == 1 ) { hnp_btageff_sf_slSig_2b_val -> SetBinContent( binIndex, btageff_sf_slSig ) ; }
                 if ( bbi == 2 ) { hnp_btageff_sf_slSig_3b_val -> SetBinContent( binIndex, btageff_sf_slSig ) ; }
                 if ( bbi == 3 ) { hnp_btageff_sf_slSig_4b_val -> SetBinContent( binIndex, btageff_sf_slSig ) ; }
                 if ( bbi == 0 ) { hnp_eff_btageff_sf_slSig_prod_1b_val -> SetBinContent( binIndex, prod_slSig ) ; }
                 if ( bbi == 1 ) { hnp_eff_btageff_sf_slSig_prod_2b_val -> SetBinContent( binIndex, prod_slSig ) ; }
                 if ( bbi == 2 ) { hnp_eff_btageff_sf_slSig_prod_3b_val -> SetBinContent( binIndex, prod_slSig ) ; }
                 if ( bbi == 3 ) { hnp_eff_btageff_sf_slSig_prod_4b_val -> SetBinContent( binIndex, prod_slSig ) ; }

                 if ( bbi == 0 ) { hnp_eff_sf_sl_1b_val -> SetBinContent( binIndex, bin_eff_sf_sl ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_sl_2b_val -> SetBinContent( binIndex, bin_eff_sf_sl ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_sl_3b_val -> SetBinContent( binIndex, bin_eff_sf_sl ) ; }
                 if ( bbi == 3 ) { hnp_eff_sf_sl_4b_val -> SetBinContent( binIndex, bin_eff_sf_sl ) ; }
                 if ( bbi == 0 ) { hnp_eff_sf_sl_1b_pull -> SetBinContent( binIndex, bin_eff_sf_sl_pull ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_sl_2b_pull -> SetBinContent( binIndex, bin_eff_sf_sl_pull ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_sl_3b_pull -> SetBinContent( binIndex, bin_eff_sf_sl_pull ) ; }
                 if ( bbi == 3 ) { hnp_eff_sf_sl_4b_pull -> SetBinContent( binIndex, bin_eff_sf_sl_pull ) ; }
                 if ( bbi == 0 ) { hnp_btageff_sf_sl_1b_val -> SetBinContent( binIndex, btageff_sf_sl ) ; }
                 if ( bbi == 1 ) { hnp_btageff_sf_sl_2b_val -> SetBinContent( binIndex, btageff_sf_sl ) ; }
                 if ( bbi == 2 ) { hnp_btageff_sf_sl_3b_val -> SetBinContent( binIndex, btageff_sf_sl ) ; }
                 if ( bbi == 3 ) { hnp_btageff_sf_sl_4b_val -> SetBinContent( binIndex, btageff_sf_sl ) ; }
                 if ( bbi == 0 ) { hnp_eff_btageff_sf_sl_prod_1b_val -> SetBinContent( binIndex, prod_sl ) ; }
                 if ( bbi == 1 ) { hnp_eff_btageff_sf_sl_prod_2b_val -> SetBinContent( binIndex, prod_sl ) ; }
                 if ( bbi == 2 ) { hnp_eff_btageff_sf_sl_prod_3b_val -> SetBinContent( binIndex, prod_sl ) ; }
                 if ( bbi == 3 ) { hnp_eff_btageff_sf_sl_prod_4b_val -> SetBinContent( binIndex, prod_sl ) ; }

                 if ( bbi == 0 ) { hnp_eff_sf_ldp_1b_val -> SetBinContent( binIndex, bin_eff_sf_ldp ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_ldp_2b_val -> SetBinContent( binIndex, bin_eff_sf_ldp ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_ldp_3b_val -> SetBinContent( binIndex, bin_eff_sf_ldp ) ; }
                 if ( bbi == 3 ) { hnp_eff_sf_ldp_4b_val -> SetBinContent( binIndex, bin_eff_sf_ldp ) ; }
                 if ( bbi == 0 ) { hnp_eff_sf_ldp_1b_pull -> SetBinContent( binIndex, bin_eff_sf_ldp_pull ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_ldp_2b_pull -> SetBinContent( binIndex, bin_eff_sf_ldp_pull ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_ldp_3b_pull -> SetBinContent( binIndex, bin_eff_sf_ldp_pull ) ; }
                 if ( bbi == 3 ) { hnp_eff_sf_ldp_4b_pull -> SetBinContent( binIndex, bin_eff_sf_ldp_pull ) ; }
                 if ( bbi == 0 ) { hnp_btageff_sf_ldp_1b_val -> SetBinContent( binIndex, btageff_sf_ldp ) ; }
                 if ( bbi == 1 ) { hnp_btageff_sf_ldp_2b_val -> SetBinContent( binIndex, btageff_sf_ldp ) ; }
                 if ( bbi == 2 ) { hnp_btageff_sf_ldp_3b_val -> SetBinContent( binIndex, btageff_sf_ldp ) ; }
                 if ( bbi == 3 ) { hnp_btageff_sf_ldp_4b_val -> SetBinContent( binIndex, btageff_sf_ldp ) ; }
                 if ( bbi == 0 ) { hnp_eff_btageff_sf_ldp_prod_1b_val -> SetBinContent( binIndex, prod_ldp ) ; }
                 if ( bbi == 1 ) { hnp_eff_btageff_sf_ldp_prod_2b_val -> SetBinContent( binIndex, prod_ldp ) ; }
                 if ( bbi == 2 ) { hnp_eff_btageff_sf_ldp_prod_3b_val -> SetBinContent( binIndex, prod_ldp ) ; }
                 if ( bbi == 3 ) { hnp_eff_btageff_sf_ldp_prod_4b_val -> SetBinContent( binIndex, prod_ldp ) ; }

                 printf(" m,h,b %d,%d,%d : eff_sf       = %6.3f,  btageff_sf       = %6.3f,  prod = %6.3f\n", mbi, hbi, bbi, bin_eff_sf, btageff_sf, prod ) ;
                 printf(" m,h,b %d,%d,%d : eff_sf_slSig = %6.3f,  btageff_sf_slSig = %6.3f,  prod = %6.3f\n", mbi, hbi, bbi, bin_eff_sf_slSig, btageff_sf_slSig, prod_slSig ) ;
                 printf(" m,h,b %d,%d,%d : eff_sf_sl    = %6.3f,  btageff_sf_sl    = %6.3f,  prod = %6.3f\n", mbi, hbi, bbi, bin_eff_sf_sl, btageff_sf_sl, prod_sl ) ;
                 printf(" m,h,b %d,%d,%d : eff_sf_ldp   = %6.3f,  btageff_sf_ldp   = %6.3f,  prod = %6.3f\n", mbi, hbi, bbi, bin_eff_sf_ldp, btageff_sf_ldp, prod_ldp ) ;

                 char label[1000] ;
                 sprintf( label, "#alpha%d_#beta%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 if ( bbi==0 ) hnp_eff_sf_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_eff_sf_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_sf_1b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_2b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_3b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_eff_sf_4b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_btageff_sf_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_btageff_sf_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_btageff_sf_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_btageff_sf_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_btageff_sf_prod_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_btageff_sf_prod_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_btageff_sf_prod_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_eff_btageff_sf_prod_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;

                 if ( bbi==0 ) hnp_eff_sf_slSig_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_slSig_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_slSig_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_eff_sf_slSig_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_sf_slSig_1b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_slSig_2b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_slSig_3b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_eff_sf_slSig_4b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_btageff_sf_slSig_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_btageff_sf_slSig_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_btageff_sf_slSig_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_btageff_sf_slSig_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_btageff_sf_slSig_prod_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_btageff_sf_slSig_prod_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_btageff_sf_slSig_prod_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_eff_btageff_sf_slSig_prod_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;

                 if ( bbi==0 ) hnp_eff_sf_sl_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_sl_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_sl_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_eff_sf_sl_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_sf_sl_1b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_sl_2b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_sl_3b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_eff_sf_sl_4b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_btageff_sf_sl_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_btageff_sf_sl_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_btageff_sf_sl_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_btageff_sf_sl_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_btageff_sf_sl_prod_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_btageff_sf_sl_prod_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_btageff_sf_sl_prod_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_eff_btageff_sf_sl_prod_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;

                 if ( bbi==0 ) hnp_eff_sf_ldp_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_ldp_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_ldp_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_eff_sf_ldp_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_sf_ldp_1b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_ldp_2b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_ldp_3b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_eff_sf_ldp_4b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_btageff_sf_ldp_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_btageff_sf_ldp_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_btageff_sf_ldp_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_btageff_sf_ldp_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_btageff_sf_ldp_prod_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_btageff_sf_ldp_prod_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_btageff_sf_ldp_prod_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==3 ) hnp_eff_btageff_sf_ldp_prod_4b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;

              }
           } // hbi
        } // mbi.

        hnp_eff_sf_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_1b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_1b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_1b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_prod_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_prod_1b_val->GetXaxis()->LabelsOption("v") ;

        hnp_eff_sf_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_2b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_2b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_2b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_prod_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_prod_2b_val->GetXaxis()->LabelsOption("v") ;

	if ( nBinsBtag > 2 ) {

	  hnp_eff_sf_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_eff_sf_3b_val->GetXaxis()->LabelsOption("v") ;
	  hnp_eff_sf_3b_pull->SetLabelSize(0.055,"x") ;
	  hnp_eff_sf_3b_pull->GetXaxis()->LabelsOption("v") ;
	  hnp_btageff_sf_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_btageff_sf_3b_val->GetXaxis()->LabelsOption("v") ;
	  hnp_eff_btageff_sf_prod_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_eff_btageff_sf_prod_3b_val->GetXaxis()->LabelsOption("v") ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    hnp_eff_sf_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_eff_sf_4b_val->GetXaxis()->LabelsOption("v") ;
	    hnp_eff_sf_4b_pull->SetLabelSize(0.055,"x") ;
	    hnp_eff_sf_4b_pull->GetXaxis()->LabelsOption("v") ;
	    hnp_btageff_sf_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_btageff_sf_4b_val->GetXaxis()->LabelsOption("v") ;
	    hnp_eff_btageff_sf_prod_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_eff_btageff_sf_prod_4b_val->GetXaxis()->LabelsOption("v") ;

	  }
	}


        hnp_eff_sf_slSig_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_slSig_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_slSig_1b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_slSig_1b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_slSig_1b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_slSig_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_slSig_prod_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_slSig_prod_1b_val->GetXaxis()->LabelsOption("v") ;

        hnp_eff_sf_slSig_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_slSig_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_slSig_2b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_slSig_2b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_slSig_2b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_slSig_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_slSig_prod_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_slSig_prod_2b_val->GetXaxis()->LabelsOption("v") ;

	if ( nBinsBtag > 2 ) {

	  hnp_eff_sf_slSig_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_eff_sf_slSig_3b_val->GetXaxis()->LabelsOption("v") ;
	  hnp_eff_sf_slSig_3b_pull->SetLabelSize(0.055,"x") ;
	  hnp_eff_sf_slSig_3b_pull->GetXaxis()->LabelsOption("v") ;
	  hnp_btageff_sf_slSig_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_btageff_sf_slSig_3b_val->GetXaxis()->LabelsOption("v") ;
	  hnp_eff_btageff_sf_slSig_prod_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_eff_btageff_sf_slSig_prod_3b_val->GetXaxis()->LabelsOption("v") ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    hnp_eff_sf_slSig_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_eff_sf_slSig_4b_val->GetXaxis()->LabelsOption("v") ;
	    hnp_eff_sf_slSig_4b_pull->SetLabelSize(0.055,"x") ;
	    hnp_eff_sf_slSig_4b_pull->GetXaxis()->LabelsOption("v") ;
	    hnp_btageff_sf_slSig_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_btageff_sf_slSig_4b_val->GetXaxis()->LabelsOption("v") ;
	    hnp_eff_btageff_sf_slSig_prod_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_eff_btageff_sf_slSig_prod_4b_val->GetXaxis()->LabelsOption("v") ;
	    
	  }
	}


        hnp_eff_sf_sl_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_sl_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_sl_1b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_sl_1b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_sl_1b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_sl_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_sl_prod_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_sl_prod_1b_val->GetXaxis()->LabelsOption("v") ;

        hnp_eff_sf_sl_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_sl_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_sl_2b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_sl_2b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_sl_2b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_sl_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_sl_prod_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_sl_prod_2b_val->GetXaxis()->LabelsOption("v") ;

	if ( nBinsBtag > 2 ) {

	  hnp_eff_sf_sl_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_eff_sf_sl_3b_val->GetXaxis()->LabelsOption("v") ;
	  hnp_eff_sf_sl_3b_pull->SetLabelSize(0.055,"x") ;
	  hnp_eff_sf_sl_3b_pull->GetXaxis()->LabelsOption("v") ;
	  hnp_btageff_sf_sl_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_btageff_sf_sl_3b_val->GetXaxis()->LabelsOption("v") ;
	  hnp_eff_btageff_sf_sl_prod_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_eff_btageff_sf_sl_prod_3b_val->GetXaxis()->LabelsOption("v") ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    hnp_eff_sf_sl_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_eff_sf_sl_4b_val->GetXaxis()->LabelsOption("v") ;
	    hnp_eff_sf_sl_4b_pull->SetLabelSize(0.055,"x") ;
	    hnp_eff_sf_sl_4b_pull->GetXaxis()->LabelsOption("v") ;
	    hnp_btageff_sf_sl_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_btageff_sf_sl_4b_val->GetXaxis()->LabelsOption("v") ;
	    hnp_eff_btageff_sf_sl_prod_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_eff_btageff_sf_sl_prod_4b_val->GetXaxis()->LabelsOption("v") ;
	    
	  }
	}


        hnp_eff_sf_ldp_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_ldp_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_ldp_1b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_ldp_1b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_ldp_1b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_ldp_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_ldp_prod_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_ldp_prod_1b_val->GetXaxis()->LabelsOption("v") ;

        hnp_eff_sf_ldp_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_ldp_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_ldp_2b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_ldp_2b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_ldp_2b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_ldp_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_ldp_prod_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_ldp_prod_2b_val->GetXaxis()->LabelsOption("v") ;

	if ( nBinsBtag > 2 ) {

	  hnp_eff_sf_ldp_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_eff_sf_ldp_3b_val->GetXaxis()->LabelsOption("v") ;
	  hnp_eff_sf_ldp_3b_pull->SetLabelSize(0.055,"x") ;
	  hnp_eff_sf_ldp_3b_pull->GetXaxis()->LabelsOption("v") ;
	  hnp_btageff_sf_ldp_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_btageff_sf_ldp_3b_val->GetXaxis()->LabelsOption("v") ;
	  hnp_eff_btageff_sf_ldp_prod_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_eff_btageff_sf_ldp_prod_3b_val->GetXaxis()->LabelsOption("v") ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    hnp_eff_sf_ldp_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_eff_sf_ldp_4b_val->GetXaxis()->LabelsOption("v") ;
	    hnp_eff_sf_ldp_4b_pull->SetLabelSize(0.055,"x") ;
	    hnp_eff_sf_ldp_4b_pull->GetXaxis()->LabelsOption("v") ;
	    hnp_btageff_sf_ldp_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_btageff_sf_ldp_4b_val->GetXaxis()->LabelsOption("v") ;
	    hnp_eff_btageff_sf_ldp_prod_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_eff_btageff_sf_ldp_prod_4b_val->GetXaxis()->LabelsOption("v") ;
	    
	  }
	}


        hnp_eff_sf_1b_val->SetFillColor(kOrange+1) ;        
	hnp_eff_sf_1b_pull->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_prod_1b_val->SetFillColor(kOrange+1) ;

        hnp_eff_sf_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_2b_pull->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_prod_2b_val->SetFillColor(kOrange+1) ;

	if ( nBinsBtag > 2 ) {

	  hnp_eff_sf_3b_val->SetFillColor(kOrange+1) ;
	  hnp_eff_sf_3b_pull->SetFillColor(kOrange+1) ;
	  hnp_btageff_sf_3b_val->SetFillColor(kOrange+1) ;
	  hnp_eff_btageff_sf_prod_3b_val->SetFillColor(kOrange+1) ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    hnp_eff_sf_4b_val->SetFillColor(kOrange+1) ;
	    hnp_eff_sf_4b_pull->SetFillColor(kOrange+1) ;
	    hnp_btageff_sf_4b_val->SetFillColor(kOrange+1) ;
	    hnp_eff_btageff_sf_prod_4b_val->SetFillColor(kOrange+1) ;

	  }
	}


        hnp_eff_sf_slSig_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_slSig_1b_pull->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_slSig_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_slSig_prod_1b_val->SetFillColor(kOrange+1) ;

        hnp_eff_sf_slSig_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_slSig_2b_pull->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_slSig_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_slSig_prod_2b_val->SetFillColor(kOrange+1) ;

	if ( nBinsBtag > 2 ) {

	  hnp_eff_sf_slSig_3b_val->SetFillColor(kOrange+1) ;
	  hnp_eff_sf_slSig_3b_pull->SetFillColor(kOrange+1) ;
	  hnp_btageff_sf_slSig_3b_val->SetFillColor(kOrange+1) ;
	  hnp_eff_btageff_sf_slSig_prod_3b_val->SetFillColor(kOrange+1) ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    hnp_eff_sf_slSig_4b_val->SetFillColor(kOrange+1) ;
	    hnp_eff_sf_slSig_4b_pull->SetFillColor(kOrange+1) ;
	    hnp_btageff_sf_slSig_4b_val->SetFillColor(kOrange+1) ;
	    hnp_eff_btageff_sf_slSig_prod_4b_val->SetFillColor(kOrange+1) ;
	    
	  }
	}


        hnp_eff_sf_sl_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_sl_1b_pull->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_sl_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_sl_prod_1b_val->SetFillColor(kOrange+1) ;

        hnp_eff_sf_sl_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_sl_2b_pull->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_sl_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_sl_prod_2b_val->SetFillColor(kOrange+1) ;

	if ( nBinsBtag > 2 ) {

	  hnp_eff_sf_sl_3b_val->SetFillColor(kOrange+1) ;
	  hnp_eff_sf_sl_3b_pull->SetFillColor(kOrange+1) ;
	  hnp_btageff_sf_sl_3b_val->SetFillColor(kOrange+1) ;
	  hnp_eff_btageff_sf_sl_prod_3b_val->SetFillColor(kOrange+1) ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    hnp_eff_sf_sl_4b_val->SetFillColor(kOrange+1) ;
	    hnp_eff_sf_sl_4b_pull->SetFillColor(kOrange+1) ;
	    hnp_btageff_sf_sl_4b_val->SetFillColor(kOrange+1) ;
	    hnp_eff_btageff_sf_sl_prod_4b_val->SetFillColor(kOrange+1) ;
	    
	  }
	}


        hnp_eff_sf_ldp_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_ldp_1b_pull->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_ldp_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_ldp_prod_1b_val->SetFillColor(kOrange+1) ;

        hnp_eff_sf_ldp_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_ldp_2b_pull->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_ldp_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_ldp_prod_2b_val->SetFillColor(kOrange+1) ;

	if ( nBinsBtag > 2 ) {

	  hnp_eff_sf_ldp_3b_val->SetFillColor(kOrange+1) ;
	  hnp_eff_sf_ldp_3b_pull->SetFillColor(kOrange+1) ;
	  hnp_btageff_sf_ldp_3b_val->SetFillColor(kOrange+1) ;
	  hnp_eff_btageff_sf_ldp_prod_3b_val->SetFillColor(kOrange+1) ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    hnp_eff_sf_ldp_4b_val->SetFillColor(kOrange+1) ;
	    hnp_eff_sf_ldp_4b_pull->SetFillColor(kOrange+1) ;
	    hnp_btageff_sf_ldp_4b_val->SetFillColor(kOrange+1) ;
	    hnp_eff_btageff_sf_ldp_prod_4b_val->SetFillColor(kOrange+1) ;
	    
	  }
	}


        hnp_ttwj_1b_val->SetFillColor(kBlue-9) ;
        hnp_ttwj_1b_pull->SetFillColor(kBlue-9) ;
        hnp_qcd_1b_val->SetFillColor(2) ;
        hnp_qcd_1b_pull->SetFillColor(2) ;
        hnp_ttwj_1b_val->SetMaximum(5.) ;
        hnp_ttwj_1b_pull->SetMinimum(-2.) ;
        hnp_qcd_1b_val->SetMaximum(5.) ;
        hnp_qcd_1b_pull->SetMinimum(-2.) ;
        hnp_ttwj_1b_pull->SetMaximum(2.) ;
        hnp_qcd_1b_pull->SetMaximum(2.) ;

        hnp_ttwj_2b_val->SetFillColor(kBlue-9) ;
        hnp_ttwj_2b_pull->SetFillColor(kBlue-9) ;
        hnp_qcd_2b_val->SetFillColor(2) ;
        hnp_qcd_2b_pull->SetFillColor(2) ;
        hnp_ttwj_2b_val->SetMaximum(5.) ;
        hnp_qcd_2b_val->SetMaximum(5.) ;
        hnp_ttwj_2b_pull->SetMinimum(-2.) ;
        hnp_qcd_2b_pull->SetMinimum(-2.) ;
        hnp_ttwj_2b_pull->SetMaximum(2.) ;
        hnp_qcd_2b_pull->SetMaximum(2.) ;

	if ( nBinsBtag > 2 ) {

	  hnp_ttwj_3b_val->SetFillColor(kBlue-9) ;
	  hnp_ttwj_3b_pull->SetFillColor(kBlue-9) ;
	  hnp_qcd_3b_val->SetFillColor(2) ;
	  hnp_qcd_3b_pull->SetFillColor(2) ;
	  hnp_ttwj_3b_val->SetMaximum(5.) ;
	  hnp_qcd_3b_val->SetMaximum(5.) ;
	  hnp_ttwj_3b_pull->SetMinimum(-2.) ;
	  hnp_qcd_3b_pull->SetMinimum(-2.) ;
	  hnp_ttwj_3b_pull->SetMaximum(2.) ;
	  hnp_qcd_3b_pull->SetMaximum(2.) ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    hnp_ttwj_4b_val->SetFillColor(kBlue-9) ;
	    hnp_ttwj_4b_pull->SetFillColor(kBlue-9) ;
	    hnp_qcd_4b_val->SetFillColor(2) ;
	    hnp_qcd_4b_pull->SetFillColor(2) ;
	    hnp_ttwj_4b_val->SetMaximum(5.) ;
	    hnp_qcd_4b_val->SetMaximum(5.) ;
	    hnp_ttwj_4b_pull->SetMinimum(-2.) ;
	    hnp_qcd_4b_pull->SetMinimum(-2.) ;
	    hnp_ttwj_4b_pull->SetMaximum(2.) ;
	    hnp_qcd_4b_pull->SetMaximum(2.) ;
	    
	  }
	}


        hnp_eff_sf_1b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_prod_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_slSig_1b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_slSig_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_slSig_prod_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_sl_1b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_sl_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_sl_prod_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_ldp_1b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_ldp_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_ldp_prod_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_1b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_1b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_prod_1b_val -> SetMaximum(2.0) ;
        hnp_eff_sf_slSig_1b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_slSig_1b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_slSig_prod_1b_val -> SetMaximum(2.0) ;


        hnp_eff_sf_2b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_prod_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_slSig_2b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_slSig_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_slSig_prod_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_sl_2b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_sl_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_sl_prod_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_ldp_2b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_ldp_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_ldp_prod_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_2b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_2b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_prod_2b_val -> SetMaximum(2.0) ;
        hnp_eff_sf_slSig_2b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_slSig_2b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_slSig_prod_2b_val -> SetMaximum(2.0) ;

	if ( nBinsBtag > 2 ) {

	  hnp_eff_sf_3b_val -> SetMinimum(-0.2) ;
	  hnp_btageff_sf_3b_val -> SetMinimum(-0.2) ;
	  hnp_eff_btageff_sf_prod_3b_val -> SetMinimum(-0.2) ;
	  hnp_eff_sf_slSig_3b_val -> SetMinimum(-0.2) ;
	  hnp_btageff_sf_slSig_3b_val -> SetMinimum(-0.2) ;
	  hnp_eff_btageff_sf_slSig_prod_3b_val -> SetMinimum(-0.2) ;
	  hnp_eff_sf_sl_3b_val -> SetMinimum(-0.2) ;
	  hnp_btageff_sf_sl_3b_val -> SetMinimum(-0.2) ;
	  hnp_eff_btageff_sf_sl_prod_3b_val -> SetMinimum(-0.2) ;
	  hnp_eff_sf_ldp_3b_val -> SetMinimum(-0.2) ;
	  hnp_btageff_sf_ldp_3b_val -> SetMinimum(-0.2) ;
	  hnp_eff_btageff_sf_ldp_prod_3b_val -> SetMinimum(-0.2) ;
	  hnp_eff_sf_3b_val -> SetMaximum(2.0) ;
	  hnp_btageff_sf_3b_val -> SetMaximum(2.0) ;
	  hnp_eff_btageff_sf_prod_3b_val -> SetMaximum(2.0) ;
	  hnp_eff_sf_slSig_3b_val -> SetMaximum(2.0) ;
	  hnp_btageff_sf_slSig_3b_val -> SetMaximum(2.0) ;
	  hnp_eff_btageff_sf_slSig_prod_3b_val -> SetMaximum(2.0) ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    hnp_eff_sf_4b_val -> SetMinimum(-0.2) ;
	    hnp_btageff_sf_4b_val -> SetMinimum(-0.2) ;
	    hnp_eff_btageff_sf_prod_4b_val -> SetMinimum(-0.2) ;
	    hnp_eff_sf_slSig_4b_val -> SetMinimum(-0.2) ;
	    hnp_btageff_sf_slSig_4b_val -> SetMinimum(-0.2) ;
	    hnp_eff_btageff_sf_slSig_prod_4b_val -> SetMinimum(-0.2) ;
	    hnp_eff_sf_sl_4b_val -> SetMinimum(-0.2) ;
	    hnp_btageff_sf_sl_4b_val -> SetMinimum(-0.2) ;
	    hnp_eff_btageff_sf_sl_prod_4b_val -> SetMinimum(-0.2) ;
	    hnp_eff_sf_ldp_4b_val -> SetMinimum(-0.2) ;
	    hnp_btageff_sf_ldp_4b_val -> SetMinimum(-0.2) ;
	    hnp_eff_btageff_sf_ldp_prod_4b_val -> SetMinimum(-0.2) ;
	    hnp_eff_sf_4b_val -> SetMaximum(2.0) ;
	    hnp_btageff_sf_4b_val -> SetMaximum(2.0) ;
	    hnp_eff_btageff_sf_prod_4b_val -> SetMaximum(2.0) ;
	    hnp_eff_sf_slSig_4b_val -> SetMaximum(2.0) ;
	    hnp_btageff_sf_slSig_4b_val -> SetMaximum(2.0) ;
	    hnp_eff_btageff_sf_slSig_prod_4b_val -> SetMaximum(2.0) ;
	    
	  }
	}


        hnp_eff_sf_sl_1b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_sl_1b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_sl_prod_1b_val -> SetMaximum(2.0) ;
        hnp_eff_sf_ldp_1b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_ldp_1b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_ldp_prod_1b_val -> SetMaximum(2.0) ;
        hnp_eff_sf_1b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_1b_pull -> SetMaximum( 2.0) ;
        hnp_eff_sf_slSig_1b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_slSig_1b_pull -> SetMaximum( 2.0) ;
        hnp_eff_sf_sl_1b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_sl_1b_pull -> SetMaximum( 2.0) ;
        hnp_eff_sf_ldp_1b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_ldp_1b_pull -> SetMaximum( 2.0) ;


        hnp_eff_sf_sl_2b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_sl_2b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_sl_prod_2b_val -> SetMaximum(2.0) ;
        hnp_eff_sf_ldp_2b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_ldp_2b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_ldp_prod_2b_val -> SetMaximum(2.0) ;
        hnp_eff_sf_2b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_2b_pull -> SetMaximum( 2.0) ;
        hnp_eff_sf_slSig_2b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_slSig_2b_pull -> SetMaximum( 2.0) ;
        hnp_eff_sf_sl_2b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_sl_2b_pull -> SetMaximum( 2.0) ;
        hnp_eff_sf_ldp_2b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_ldp_2b_pull -> SetMaximum( 2.0) ;

	if ( nBinsBtag > 2 ) {

	  hnp_eff_sf_sl_3b_val -> SetMaximum(2.0) ;
	  hnp_btageff_sf_sl_3b_val -> SetMaximum(2.0) ;
	  hnp_eff_btageff_sf_sl_prod_3b_val -> SetMaximum(2.0) ;
	  hnp_eff_sf_ldp_3b_val -> SetMaximum(2.0) ;
	  hnp_btageff_sf_ldp_3b_val -> SetMaximum(2.0) ;
	  hnp_eff_btageff_sf_ldp_prod_3b_val -> SetMaximum(2.0) ;
	  hnp_eff_sf_3b_pull -> SetMinimum(-2.0) ;
	  hnp_eff_sf_3b_pull -> SetMaximum( 2.0) ;
	  hnp_eff_sf_slSig_3b_pull -> SetMinimum(-2.0) ;
	  hnp_eff_sf_slSig_3b_pull -> SetMaximum( 2.0) ;
	  hnp_eff_sf_sl_3b_pull -> SetMinimum(-2.0) ;
	  hnp_eff_sf_sl_3b_pull -> SetMaximum( 2.0) ;
	  hnp_eff_sf_ldp_3b_pull -> SetMinimum(-2.0) ;
	  hnp_eff_sf_ldp_3b_pull -> SetMaximum( 2.0) ;
	  
	  if ( nBinsBtag > 3 ) {

	    hnp_eff_sf_sl_4b_val -> SetMaximum(2.0) ;
	    hnp_btageff_sf_sl_4b_val -> SetMaximum(2.0) ;
	    hnp_eff_btageff_sf_sl_prod_4b_val -> SetMaximum(2.0) ;
	    hnp_eff_sf_ldp_4b_val -> SetMaximum(2.0) ;
	    hnp_btageff_sf_ldp_4b_val -> SetMaximum(2.0) ;
	    hnp_eff_btageff_sf_ldp_prod_4b_val -> SetMaximum(2.0) ;
	    hnp_eff_sf_4b_pull -> SetMinimum(-2.0) ;
	    hnp_eff_sf_4b_pull -> SetMaximum( 2.0) ;
	    hnp_eff_sf_slSig_4b_pull -> SetMinimum(-2.0) ;
	    hnp_eff_sf_slSig_4b_pull -> SetMaximum( 2.0) ;
	    hnp_eff_sf_sl_4b_pull -> SetMinimum(-2.0) ;
	    hnp_eff_sf_sl_4b_pull -> SetMaximum( 2.0) ;
	    hnp_eff_sf_ldp_4b_pull -> SetMinimum(-2.0) ;
	    hnp_eff_sf_ldp_4b_pull -> SetMaximum( 2.0) ;

	  }
	}


        hnp_ttwj_1b_val->SetLabelSize(0.055,"x") ;
        hnp_ttwj_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_qcd_1b_val->SetLabelSize(0.055,"x") ;
        hnp_qcd_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_ttwj_1b_pull->SetLabelSize(0.055,"x") ;
        hnp_ttwj_1b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_qcd_1b_pull->SetLabelSize(0.055,"x") ;
        hnp_qcd_1b_pull->GetXaxis()->LabelsOption("v") ;

        hnp_ttwj_2b_val->SetLabelSize(0.055,"x") ;
        hnp_ttwj_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_qcd_2b_val->SetLabelSize(0.055,"x") ;
        hnp_qcd_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_ttwj_2b_pull->SetLabelSize(0.055,"x") ;
        hnp_ttwj_2b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_qcd_2b_pull->SetLabelSize(0.055,"x") ;
        hnp_qcd_2b_pull->GetXaxis()->LabelsOption("v") ;

	if ( nBinsBtag > 2 ) {

	  hnp_ttwj_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_ttwj_3b_val->GetXaxis()->LabelsOption("v") ;
	  hnp_qcd_3b_val->SetLabelSize(0.055,"x") ;
	  hnp_qcd_3b_val->GetXaxis()->LabelsOption("v") ;
	  hnp_ttwj_3b_pull->SetLabelSize(0.055,"x") ;
	  hnp_ttwj_3b_pull->GetXaxis()->LabelsOption("v") ;
	  hnp_qcd_3b_pull->SetLabelSize(0.055,"x") ;
	  hnp_qcd_3b_pull->GetXaxis()->LabelsOption("v") ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    hnp_ttwj_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_ttwj_4b_val->GetXaxis()->LabelsOption("v") ;
	    hnp_qcd_4b_val->SetLabelSize(0.055,"x") ;
	    hnp_qcd_4b_val->GetXaxis()->LabelsOption("v") ;
	    hnp_ttwj_4b_pull->SetLabelSize(0.055,"x") ;
	    hnp_ttwj_4b_pull->GetXaxis()->LabelsOption("v") ;
	    hnp_qcd_4b_pull->SetLabelSize(0.055,"x") ;
	    hnp_qcd_4b_pull->GetXaxis()->LabelsOption("v") ;

	  }
	}


        //=== Efficiency par pulls ============


        TCanvas* cnp2 = (TCanvas*) gDirectory->FindObject("cnp2") ;
        if ( cnp2 == 0x0 ) {
           cnp2 = new TCanvas("cnp2","RA2b efficiency par pull", 850, 1000 ) ;
        }
        cnp2->Divide(nBinsBtag,4) ;

        cnp2->cd(1) ;
        hnp_eff_sf_1b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_1b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp2->cd(2) ;
        hnp_eff_sf_2b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_2b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

	if ( nBinsBtag > 2 ) {
	  
	  cnp2->cd(3) ;
	  hnp_eff_sf_3b_pull->Draw() ;
	  chi2 = addChi2FromPullHist( hnp_eff_sf_3b_pull ) ;
	  globalChi2 += chi2 ;
	  npChi2 += chi2 ;

	  if ( nBinsBtag > 3 ) {

	    cnp2->cd(4) ;
	    hnp_eff_sf_4b_pull->Draw() ;
	    chi2 = addChi2FromPullHist( hnp_eff_sf_4b_pull ) ;
	    globalChi2 += chi2 ;
	    npChi2 += chi2 ;

	  }
	}


        cnp2->cd(nBinsBtag+1) ;
        hnp_eff_sf_slSig_1b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_slSig_1b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp2->cd(nBinsBtag+2) ;
        hnp_eff_sf_slSig_2b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_slSig_2b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

	if ( nBinsBtag > 2 ) {

	  cnp2->cd(nBinsBtag+3) ;
	  hnp_eff_sf_slSig_3b_pull->Draw() ;
	  chi2 = addChi2FromPullHist( hnp_eff_sf_slSig_3b_pull ) ;
	  globalChi2 += chi2 ;
	  npChi2 += chi2 ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    cnp2->cd(nBinsBtag+4) ;
	    hnp_eff_sf_slSig_4b_pull->Draw() ;
	    chi2 = addChi2FromPullHist( hnp_eff_sf_slSig_4b_pull ) ;
	    globalChi2 += chi2 ;
	    npChi2 += chi2 ;
	    
	  }
	}



        cnp2->cd(2*nBinsBtag+1) ;
        hnp_eff_sf_sl_1b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_sl_1b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp2->cd(2*nBinsBtag+2) ;
        hnp_eff_sf_sl_2b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_sl_2b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

	if ( nBinsBtag > 2 ) {

	  cnp2->cd(2*nBinsBtag+3) ;
	  hnp_eff_sf_sl_3b_pull->Draw() ;
	  chi2 = addChi2FromPullHist( hnp_eff_sf_sl_3b_pull ) ;
	  globalChi2 += chi2 ;
	  npChi2 += chi2 ;

	  if ( nBinsBtag > 3 ) {
	    
	    cnp2->cd(2*nBinsBtag+4) ;
	    hnp_eff_sf_sl_4b_pull->Draw() ;
	    chi2 = addChi2FromPullHist( hnp_eff_sf_sl_4b_pull ) ;
	    globalChi2 += chi2 ;
	    npChi2 += chi2 ;
	    
	  }
	}


        cnp2->cd(3*nBinsBtag+1) ;
        hnp_eff_sf_ldp_1b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_ldp_1b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp2->cd(3*nBinsBtag+2) ;
        hnp_eff_sf_ldp_2b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_ldp_2b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

	if ( nBinsBtag > 2 ) {

	  cnp2->cd(3*nBinsBtag+3) ;
	  hnp_eff_sf_ldp_3b_pull->Draw() ;
	  chi2 = addChi2FromPullHist( hnp_eff_sf_ldp_3b_pull ) ;
	  globalChi2 += chi2 ;
	  npChi2 += chi2 ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    cnp2->cd(3*nBinsBtag+4) ;
	    hnp_eff_sf_ldp_4b_pull->Draw() ;
	    chi2 = addChi2FromPullHist( hnp_eff_sf_ldp_4b_pull ) ;
	    globalChi2 += chi2 ;
	    npChi2 += chi2 ;
	    
	  }
	}

	cnp2->SaveAs("outputfiles/cnp2.pdf");




       //==== TTwj and QCD scale factors ========================

        gStyle->SetPadGridY(1) ;
        TCanvas* cnp = (TCanvas*) gDirectory->FindObject("cnp") ;
        if ( cnp == 0x0 ) {
           cnp = new TCanvas("cnp","RA2b nuisance pars", 850, 1000 ) ;
        }
        cnp->Divide(nBinsBtag,4) ;

      //---
        cnp->cd(1) ;
        hnp_ttwj_1b_val->Draw() ;
        hnp_ttwj_1b_nom->Draw("same") ;

        cnp->cd(2) ;
        hnp_ttwj_2b_val->Draw() ;
        hnp_ttwj_2b_nom->Draw("same") ;

	if ( nBinsBtag > 2 ) {

	  cnp->cd(3) ;
	  hnp_ttwj_3b_val->Draw() ;
	  hnp_ttwj_3b_nom->Draw("same") ;
	  
	  if ( nBinsBtag > 3 ) {

	    cnp->cd(4) ;
	    hnp_ttwj_4b_val->Draw() ;
	    hnp_ttwj_4b_nom->Draw("same") ;
	    
	  }
	}
	

      //---

        cnp->cd(nBinsBtag+1) ;
        hnp_ttwj_1b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_ttwj_1b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;


        cnp->cd(nBinsBtag+2) ;
        hnp_ttwj_2b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_ttwj_2b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

	if ( nBinsBtag > 2 ) {

	  cnp->cd(nBinsBtag+3) ;
	  hnp_ttwj_3b_pull->Draw() ;
	  chi2 = addChi2FromPullHist( hnp_ttwj_3b_pull ) ;
	  globalChi2 += chi2 ;
	  npChi2 += chi2 ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    cnp->cd(nBinsBtag+4) ;
	    hnp_ttwj_4b_pull->Draw() ;
	    chi2 = addChi2FromPullHist( hnp_ttwj_4b_pull ) ;
	    globalChi2 += chi2 ;
	    npChi2 += chi2 ;

	  }
	}


      //---

        cnp->cd(2*nBinsBtag+1) ;
        hnp_qcd_1b_val->Draw() ;
        hnp_qcd_1b_nom->Draw("same") ;

        cnp->cd(2*nBinsBtag+2) ;
        hnp_qcd_2b_val->Draw() ;
        hnp_qcd_2b_nom->Draw("same") ;

	if ( nBinsBtag > 2 ) {
	
	  cnp->cd(2*nBinsBtag+3) ;
	  hnp_qcd_3b_val->Draw() ;
	  hnp_qcd_3b_nom->Draw("same") ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    cnp->cd(2*nBinsBtag+4) ;
	    hnp_qcd_4b_val->Draw() ;
	    hnp_qcd_4b_nom->Draw("same") ;
	    
	  }
	}


      //---

        cnp->cd(3*nBinsBtag+1) ;
        hnp_qcd_1b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_qcd_1b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp->cd(3*nBinsBtag+2) ;
        hnp_qcd_2b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_qcd_2b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

	if ( nBinsBtag > 2 ) {

	  cnp->cd(3*nBinsBtag+3) ;
	  hnp_qcd_3b_pull->Draw() ;
	  chi2 = addChi2FromPullHist( hnp_qcd_3b_pull ) ;
	  globalChi2 += chi2 ;
	  npChi2 += chi2 ;
	  
	  if ( nBinsBtag > 3 ) {
	    
	    cnp->cd(3*nBinsBtag+4) ;
	    hnp_qcd_4b_pull->Draw() ;
	    chi2 = addChi2FromPullHist( hnp_qcd_4b_pull ) ;
	    globalChi2 += chi2 ;
	    npChi2 += chi2 ;
	    
	  }
	}

	cnp->SaveAs("outputfiles/cnp.pdf");



      //--- general, Znn, and trigger

        TCanvas* cnp3 = (TCanvas*) gDirectory->FindObject("cnp3") ;
        if ( cnp3 == 0x0 ) {
           cnp3 = new TCanvas("cnp3","RA2b efficiency par pull", 850, 600 ) ;
        }
        cnp3->Divide(2,2) ;


        cnp3->cd(1) ;
        hnp_prim_eff->Draw() ;
        chi2 = addChi2FromPullHist( hnp_prim_eff ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;
	
	
        cnp3->cd(2) ;
        hnp_znn->Draw() ;
        chi2 = addChi2FromPullHist( hnp_znn ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;
	

        cnp3->cd(3) ;
        hnp_trig_0lep->Draw() ;
        chi2 = addChi2FromPullHist( hnp_trig_0lep ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;
	

        cnp3->cd(4) ;
        hnp_trig_1lep->Draw() ;
        chi2 = addChi2FromPullHist( hnp_trig_1lep ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;


	cnp3->SaveAs("outputfiles/cnp3.pdf");


     } else {
        printf("\n\n *** Skipping nuisance parameters.\n\n") ;
     }


     //====== Observables ==================================


     TCanvas* cfitqual = (TCanvas*) gDirectory->FindObject("cfitqual") ;
     if ( cfitqual == 0x0 ) {
        cfitqual = new TCanvas("cfitqual","RA2b fit quality", 850, 1000 ) ;
     }

     Int_t nColumns = max(3,nBinsBtag);

     cfitqual->Divide(nColumns,5);

     gPad->SetTicks(1,0) ;

     hfitqual_data_0lep_1b->SetLabelSize(0.055,"x") ;
     hfitqual_data_0lep_1b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_1lepSig_1b->SetLabelSize(0.055,"x") ;
     hfitqual_data_1lepSig_1b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_1lep_1b->SetLabelSize(0.055,"x") ;
     hfitqual_data_1lep_1b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_ldp_1b->SetLabelSize(0.055,"x") ;
     hfitqual_data_ldp_1b->GetXaxis()->LabelsOption("v") ;
     
     hfitqual_data_0lep_2b->SetLabelSize(0.055,"x") ;
     hfitqual_data_0lep_2b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_1lepSig_2b->SetLabelSize(0.055,"x") ;
     hfitqual_data_1lepSig_2b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_1lep_2b->SetLabelSize(0.055,"x") ;
     hfitqual_data_1lep_2b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_ldp_2b->SetLabelSize(0.055,"x") ;
     hfitqual_data_ldp_2b->GetXaxis()->LabelsOption("v") ;

     if ( nBinsBtag > 2 ) {
     
       hfitqual_data_0lep_3b->SetLabelSize(0.055,"x") ;
       hfitqual_data_0lep_3b->GetXaxis()->LabelsOption("v") ;
       hfitqual_data_1lepSig_3b->SetLabelSize(0.055,"x") ;
       hfitqual_data_1lepSig_3b->GetXaxis()->LabelsOption("v") ;
       hfitqual_data_1lep_3b->SetLabelSize(0.055,"x") ;
       hfitqual_data_1lep_3b->GetXaxis()->LabelsOption("v") ;
       hfitqual_data_ldp_3b->SetLabelSize(0.055,"x") ;
       hfitqual_data_ldp_3b->GetXaxis()->LabelsOption("v") ;
       
       if ( nBinsBtag > 3 ) {
	 
	 hfitqual_data_0lep_4b->SetLabelSize(0.055,"x") ;
	 hfitqual_data_0lep_4b->GetXaxis()->LabelsOption("v") ;
	 hfitqual_data_1lepSig_4b->SetLabelSize(0.055,"x") ;
	 hfitqual_data_1lepSig_4b->GetXaxis()->LabelsOption("v") ;
	 hfitqual_data_1lep_4b->SetLabelSize(0.055,"x") ;
	 hfitqual_data_1lep_4b->GetXaxis()->LabelsOption("v") ;
	 hfitqual_data_ldp_4b->SetLabelSize(0.055,"x") ;
	 hfitqual_data_ldp_4b->GetXaxis()->LabelsOption("v") ;
	 
       }
     }

     hfitqual_data_zee_1b->SetLabelSize(0.055,"x") ;
     hfitqual_data_zee_1b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_zmm_1b->SetLabelSize(0.055,"x") ;
     hfitqual_data_zmm_1b->GetXaxis()->LabelsOption("v") ;






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



     
     cfitqual->cd(3*nColumns+1);
     hfitqual_data_ldp_1b->Draw("histpe") ;
     hfitqual_fit_ldp_1b->Draw("same") ;
     hfitqual_data_ldp_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_ldp_1b, hfitqual_model_ldp_1b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;

     cfitqual->cd(3*nColumns+2);
     hfitqual_data_ldp_2b->Draw("histpe") ;
     hfitqual_fit_ldp_2b->Draw("same") ;
     hfitqual_data_ldp_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_ldp_2b, hfitqual_model_ldp_2b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;

     if ( nBinsBtag > 2 ) {
	 
	 cfitqual->cd(3*nColumns+3);
	 hfitqual_data_ldp_3b->Draw("histpe") ;
	 hfitqual_fit_ldp_3b->Draw("same") ;
	 hfitqual_data_ldp_3b->Draw("same") ;
	 gPad->SetGridy(1) ;
	 chi2 = addChi2FromObs( hfitqual_data_ldp_3b, hfitqual_model_ldp_3b ) ;
	 globalChi2 += chi2 ;
	 obsChi2 += chi2 ;
       
       if ( nBinsBtag > 3 ) {
	 
	 cfitqual->cd(3*nColumns+4);
	 hfitqual_data_ldp_4b->Draw("histpe") ;
	 hfitqual_fit_ldp_4b->Draw("same") ;
	 hfitqual_data_ldp_4b->Draw("same") ;
	 gPad->SetGridy(1) ;
	 chi2 = addChi2FromObs( hfitqual_data_ldp_4b, hfitqual_model_ldp_4b ) ;
	 globalChi2 += chi2 ;
	 obsChi2 += chi2 ;
	 
       }
     }


     cfitqual->cd(4*nColumns+1) ;
     hfitqual_data_zee_1b->Draw("histpe") ;
     hfitqual_fit_zee_1b->Draw("same") ;
     hfitqual_data_zee_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_zee_1b, hfitqual_fit_zee_1b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;

     cfitqual->cd(4*nColumns+2) ;
     hfitqual_data_zmm_1b->Draw("histpe") ;
     hfitqual_fit_zmm_1b->Draw("same") ;
     hfitqual_data_zmm_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_zmm_1b, hfitqual_fit_zmm_1b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;


     cfitqual->cd(nColumns*5);
     legend->Draw() ;

     cfitqual->Update() ;



     cfitqual->cd(nColumns*5) ;
     TText* chi2text = new TText() ;
     chi2text->SetTextSize(0.055) ;
     chi2text->SetTextAlign(32) ;
     char chi2string[1000] ;

     sprintf( chi2string, "Overall Chi2 = %6.2f", globalChi2 ) ;
     chi2text->DrawTextNDC(0.55, 0.90, chi2string ) ;
     sprintf( chi2string, "obs Chi2 = %6.2f", obsChi2 ) ;
     chi2text->DrawTextNDC(0.55, 0.85, chi2string ) ;
     sprintf( chi2string, "NP Chi2 = %6.2f", npChi2 ) ;
     chi2text->DrawTextNDC(0.55, 0.80, chi2string ) ;

     printf("\n\n Overall Chi2 = %6.2f\n\n", globalChi2 ) ;
     printf("     Obs Chi2 = %6.2f\n", obsChi2 ) ;
     printf("      NP Chi2 = %6.2f\n", npChi2 ) ;


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

     cfitqual->SaveAs("outputfiles/cfitqual.pdf");






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









