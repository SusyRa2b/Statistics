
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


     // hardcode here the number of bins of the analysis
     int nBinsBtag = 3 ;



     double hmax = 1.25 ;

     TFile* wstf = new TFile( wsfile ) ;

     RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
     ws->Print() ;


     // figure out the binning from the variables in the likelihood
     int nBinsMET(1) ;
     for ( int bi=2; bi<20; bi++ ) {
        char vname[1000] ;
        sprintf( vname, "mu_ttwj_sl_M%d_H2_1b", bi ) ;
        RooRealVar* rrv = ws->var( vname ) ;
        if ( rrv == 0x0 ) break ;
        nBinsMET = bi ;
     }

     int nBinsHT(1) ;
     for ( int bi=2; bi<20; bi++ ) {
        char vname[1000] ;
        sprintf( vname, "mu_ttwj_sl_M1_H%d_1b", bi ) ;
        RooRealVar* rrv = ws->var( vname ) ;
        if ( rrv == 0x0 ) break ;
        nBinsHT = bi ;
     }

      //-- Hardwire in to ignore the highest MET bin in the lowest HT bin.
      bool ignoreBin[4][4] ;
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

     printf("\n\n Binning: nBinsMET = %d, nBinsHT = %d\n\n", nBinsMET, nBinsHT ) ;

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


     double AttwjVal[nBinsMET][nBinsHT][nBinsBtag] ;
     double AqcdVal[nBinsMET][nBinsHT][nBinsBtag] ;

     for ( int i = 0 ; i < nBinsMET ; i++ ) {
       for ( int j = 0 ; j < nBinsHT ; j++ ) {

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

     int dataN_0lep[nBinsMET][nBinsHT][nBinsBtag] ;
     int dataN_1lep[nBinsMET][nBinsHT][nBinsBtag] ;
     int dataN_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
     int dataN_Zee[nBinsMET][nBinsHT] ;
     int dataN_Zmm[nBinsMET][nBinsHT] ;

     for ( int i = 0 ; i < nBinsMET ; i++ ) {
       for ( int j = 0 ; j < nBinsHT ; j++ ) {
         for ( int k = 0 ; k < nBinsBtag ; k++ ) {
            dataN_0lep[i][j][k] = 0. ;
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

       for ( int i = 0 ; i < nBinsMET ; i++ ) {
	 for ( int j = 0 ; j < nBinsHT ; j++ ) {
           if ( ignoreBin[i][j] ) continue ;
	   for ( int k = 0 ; k < nBinsBtag ; k++ ) {

	     TString String_0lep  = "N_0lep" ;
	     TString String_1lep  = "N_1lep" ;
	     TString String_ldp   = "N_ldp" ;

	     String_0lep  += sMbins[i]+sHbins[j]+sBbins[k] ;
	     String_1lep  += sMbins[i]+sHbins[j]+sBbins[k] ;
	     String_ldp   += sMbins[i]+sHbins[j]+sBbins[k] ;
	     
	     if ( strcmp( obs->GetName(), String_0lep ) == 0 ) { dataN_0lep[i][j][k] = obs->getVal() ; }
	     if ( strcmp( obs->GetName(), String_1lep ) == 0 ) { dataN_1lep[i][j][k] = obs->getVal() ; }
	     if ( strcmp( obs->GetName(), String_ldp  ) == 0 ) { dataN_ldp[i][j][k]  = obs->getVal() ; }

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



     int nbins = nBinsMET*(nBinsHT+1) + 1 ;

     TH1F* hfitqual_data_0lep_1b = new TH1F("hfitqual_data_0lep_1b", "0 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_0lep_1b = new TH1F("hfitqual_susy_0lep_1b", "0 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_0lep_1b = new TH1F("hfitqual_ttwj_0lep_1b", "0 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_0lep_1b  = new TH1F("hfitqual_qcd_0lep_1b" , "0 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_0lep_1b  = new TH1F("hfitqual_znn_0lep_1b" , "0 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_model_0lep_1b  = new TH1F("hfitqual_model_0lep_1b" , "0 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;

     TH1F* hfitqual_data_0lep_2b = new TH1F("hfitqual_data_0lep_2b", "0 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_0lep_2b = new TH1F("hfitqual_susy_0lep_2b", "0 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_0lep_2b = new TH1F("hfitqual_ttwj_0lep_2b", "0 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_0lep_2b  = new TH1F("hfitqual_qcd_0lep_2b" , "0 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_0lep_2b  = new TH1F("hfitqual_znn_0lep_2b" , "0 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_model_0lep_2b  = new TH1F("hfitqual_model_0lep_2b" , "0 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     TH1F* hfitqual_data_0lep_3b = new TH1F("hfitqual_data_0lep_3b", "0 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_0lep_3b = new TH1F("hfitqual_susy_0lep_3b", "0 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_0lep_3b = new TH1F("hfitqual_ttwj_0lep_3b", "0 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_0lep_3b  = new TH1F("hfitqual_qcd_0lep_3b" , "0 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_0lep_3b  = new TH1F("hfitqual_znn_0lep_3b" , "0 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_model_0lep_3b  = new TH1F("hfitqual_model_0lep_3b" , "0 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;



     TH1F* hfitqual_data_1lep_1b = new TH1F("hfitqual_data_1lep_1b", "1 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_1lep_1b = new TH1F("hfitqual_susy_1lep_1b", "1 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_1lep_1b = new TH1F("hfitqual_ttwj_1lep_1b", "1 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_1lep_1b  = new TH1F("hfitqual_qcd_1lep_1b" , "1 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_1lep_1b  = new TH1F("hfitqual_znn_1lep_1b" , "1 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_model_1lep_1b  = new TH1F("hfitqual_model_1lep_1b" , "1 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;

     TH1F* hfitqual_data_1lep_2b = new TH1F("hfitqual_data_1lep_2b", "1 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_1lep_2b = new TH1F("hfitqual_susy_1lep_2b", "1 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_1lep_2b = new TH1F("hfitqual_ttwj_1lep_2b", "1 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_1lep_2b  = new TH1F("hfitqual_qcd_1lep_2b" , "1 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_1lep_2b  = new TH1F("hfitqual_znn_1lep_2b" , "1 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_model_1lep_2b  = new TH1F("hfitqual_model_1lep_2b" , "1 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     TH1F* hfitqual_data_1lep_3b = new TH1F("hfitqual_data_1lep_3b", "1 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_1lep_3b = new TH1F("hfitqual_susy_1lep_3b", "1 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_1lep_3b = new TH1F("hfitqual_ttwj_1lep_3b", "1 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_1lep_3b  = new TH1F("hfitqual_qcd_1lep_3b" , "1 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_1lep_3b  = new TH1F("hfitqual_znn_1lep_3b" , "1 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_model_1lep_3b  = new TH1F("hfitqual_model_1lep_3b" , "1 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;



     TH1F* hfitqual_data_ldp_1b = new TH1F("hfitqual_data_ldp_1b", "LDP, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_ldp_1b = new TH1F("hfitqual_susy_ldp_1b", "LDP, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_ldp_1b = new TH1F("hfitqual_ttwj_ldp_1b", "LDP, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_ldp_1b  = new TH1F("hfitqual_qcd_ldp_1b" , "LDP, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_ldp_1b  = new TH1F("hfitqual_znn_ldp_1b" , "LDP, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_model_ldp_1b  = new TH1F("hfitqual_model_ldp_1b" , "LDP, 1 btag" , nbins, 0.5, nbins+0.5 ) ;

     TH1F* hfitqual_data_ldp_2b = new TH1F("hfitqual_data_ldp_2b", "LDP, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_ldp_2b = new TH1F("hfitqual_susy_ldp_2b", "LDP, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_ldp_2b = new TH1F("hfitqual_ttwj_ldp_2b", "LDP, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_ldp_2b  = new TH1F("hfitqual_qcd_ldp_2b" , "LDP, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_ldp_2b  = new TH1F("hfitqual_znn_ldp_2b" , "LDP, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_model_ldp_2b  = new TH1F("hfitqual_model_ldp_2b" , "LDP, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     TH1F* hfitqual_data_ldp_3b = new TH1F("hfitqual_data_ldp_3b", "LDP, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_ldp_3b = new TH1F("hfitqual_susy_ldp_3b", "LDP, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_ldp_3b = new TH1F("hfitqual_ttwj_ldp_3b", "LDP, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_ldp_3b  = new TH1F("hfitqual_qcd_ldp_3b" , "LDP, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_ldp_3b  = new TH1F("hfitqual_znn_ldp_3b" , "LDP, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_model_ldp_3b  = new TH1F("hfitqual_model_ldp_3b" , "LDP, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;


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
     
     hfitqual_ttwj_0lep_3b  -> SetFillColor(kBlue-9) ;
     hfitqual_qcd_0lep_3b   -> SetFillColor(2) ;
     hfitqual_znn_0lep_3b   -> SetFillColor(kGreen-3) ;
     hfitqual_susy_0lep_3b  -> SetFillColor(6) ;
     hfitqual_data_0lep_3b->SetMarkerStyle(20) ;
     hfitqual_data_0lep_3b->SetLineWidth(2) ;



     
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

     hfitqual_ttwj_1lep_3b  -> SetFillColor(kBlue-9) ;
     hfitqual_qcd_1lep_3b   -> SetFillColor(2) ;
     hfitqual_znn_1lep_3b   -> SetFillColor(kGreen-3) ;
     hfitqual_susy_1lep_3b  -> SetFillColor(6) ;
     hfitqual_data_1lep_3b->SetMarkerStyle(20) ;
     hfitqual_data_1lep_3b->SetLineWidth(2) ;




     
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

     hfitqual_ttwj_ldp_3b  -> SetFillColor(kBlue-9) ;
     hfitqual_qcd_ldp_3b   -> SetFillColor(2) ;
     hfitqual_znn_ldp_3b   -> SetFillColor(kGreen-3) ;
     hfitqual_susy_ldp_3b  -> SetFillColor(6) ;
     hfitqual_data_ldp_3b->SetMarkerStyle(20) ;
     hfitqual_data_ldp_3b->SetLineWidth(2) ;




     hfitqual_np    -> SetFillColor(kOrange+1) ;

     THStack* hfitqual_fit_0lep_1b = new THStack( "hfitqual_fit_0lep_1b", "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_0lep_2b = new THStack( "hfitqual_fit_0lep_2b", "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_0lep_3b = new THStack( "hfitqual_fit_0lep_3b", "RA2b likelihood fit results, fit" ) ;

     THStack* hfitqual_fit_1lep_1b = new THStack( "hfitqual_fit_1lep_1b", "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_1lep_2b = new THStack( "hfitqual_fit_1lep_2b", "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_1lep_3b = new THStack( "hfitqual_fit_1lep_3b", "RA2b likelihood fit results, fit" ) ;

     THStack* hfitqual_fit_ldp_1b  = new THStack( "hfitqual_fit_ldp_1b",  "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_ldp_2b  = new THStack( "hfitqual_fit_ldp_2b",  "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_ldp_3b  = new THStack( "hfitqual_fit_ldp_3b",  "RA2b likelihood fit results, fit" ) ;



     TAxis* xaxis_0lep_1b = hfitqual_data_0lep_1b->GetXaxis() ;
     TAxis* xaxis_0lep_2b = hfitqual_data_0lep_2b->GetXaxis() ;
     TAxis* xaxis_0lep_3b = hfitqual_data_0lep_3b->GetXaxis() ;

     TAxis* xaxis_1lep_1b = hfitqual_data_1lep_1b->GetXaxis() ;
     TAxis* xaxis_1lep_2b = hfitqual_data_1lep_2b->GetXaxis() ;
     TAxis* xaxis_1lep_3b = hfitqual_data_1lep_3b->GetXaxis() ;

     TAxis* xaxis_ldp_1b  = hfitqual_data_ldp_1b->GetXaxis() ;
     TAxis* xaxis_ldp_2b  = hfitqual_data_ldp_2b->GetXaxis() ;
     TAxis* xaxis_ldp_3b  = hfitqual_data_ldp_3b->GetXaxis() ;

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
     double eff_sf_sl(0.) ;
     double eff_sf_ldp(0.) ;

     double sf_mc(0.) ;


     printf("\n\n Looping through analysis bins.\n\n") ; cout << flush ;

     // loop through all the bins of the analysis

     for ( int i = 0 ; i < nBinsMET ; i++ ) {
       for ( int j = 0 ; j < nBinsHT ; j++ ) {
         if ( ignoreBin[i][j] ) continue ;
         for ( int k = 0 ; k < nBinsBtag ; k++ ) {

           printf(" here t1\n") ; cout << flush ;

           binIndex = 1 + (nBinsHT+1)*i + j + 1 ;


           //+++++ 0 lep histograms:

           binLabel = "0 lep" ;
           binLabel += sMbins[i]+sHbins[j]+sBbins[k] ;

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

              sprintf( pname, "SFqcd_met%d", i+1 ) ;
              RooAbsReal* sf_qcd_met = (RooAbsReal*) ws->obj(pname) ;

              sprintf( pname, "SFqcd_nb%d", k+1 ) ;
              RooAbsReal* sf_qcd_nb = (RooAbsReal*) ws->obj(pname) ;

              if ( mu_qcd_ldp != 0 && sf_qcd != 0 && r_qcd_ht != 0 && sf_qcd_met != 0 && sf_qcd_nb != 0 ) {
                 qcdVal = trigeff * (mu_qcd_ldp -> getVal())
                        * (sf_qcd -> getVal())
                        * (r_qcd_ht -> getVal())
                        * (sf_qcd_met -> getVal())
                        * (sf_qcd_nb -> getVal()) ;
              } else {
                 printf("\n\n *** missing one of the inputs. mu_qcd_ldp %p : sf_qcd %p : r_qcd_ht %p : sf_qcd_met %p : sf_qcd_nb %p\n",
                     mu_qcd_ldp, sf_qcd, r_qcd_ht, sf_qcd_met, sf_qcd_nb ) ;
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
           }




           //+++++ 1 lep histograms:

           binLabel = "1 lep" ;
           binLabel += sMbins[i]+sHbins[j]+sBbins[k] ;

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
           }





           //+++++ ldp histograms:

           binLabel = "ldp" ;
           binLabel += sMbins[i]+sHbins[j]+sBbins[k] ;

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
           }


           printf(" here 3\n") ; cout << flush ;


           if ( k == 0 ) {

           printf(" here 4\n") ; cout << flush ;

           //+++++ Zee histogram:

              binLabel = "Zee" ;
              binLabel += sMbins[i]+sHbins[j]+sBbins[k] ;
             

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
              binLabel += sMbins[i]+sHbins[j]+sBbins[k] ;
             

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
       } // j: HT
     } // i: MET


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

     hfitqual_fit_0lep_3b->Add( hfitqual_znn_0lep_3b ) ;
     hfitqual_fit_0lep_3b->Add( hfitqual_qcd_0lep_3b ) ;
     hfitqual_fit_0lep_3b->Add( hfitqual_ttwj_0lep_3b ) ;
     hfitqual_fit_0lep_3b->Add( hfitqual_susy_0lep_3b ) ;



     hfitqual_fit_1lep_1b->Add( hfitqual_znn_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_qcd_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_ttwj_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_susy_1lep_1b ) ;

     hfitqual_fit_1lep_2b->Add( hfitqual_znn_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_qcd_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_ttwj_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_susy_1lep_2b ) ;

     hfitqual_fit_1lep_3b->Add( hfitqual_znn_1lep_3b ) ;
     hfitqual_fit_1lep_3b->Add( hfitqual_qcd_1lep_3b ) ;
     hfitqual_fit_1lep_3b->Add( hfitqual_ttwj_1lep_3b ) ;
     hfitqual_fit_1lep_3b->Add( hfitqual_susy_1lep_3b ) ;




     hfitqual_fit_ldp_1b->Add( hfitqual_znn_ldp_1b ) ;
     hfitqual_fit_ldp_1b->Add( hfitqual_qcd_ldp_1b ) ;
     hfitqual_fit_ldp_1b->Add( hfitqual_ttwj_ldp_1b ) ;
     hfitqual_fit_ldp_1b->Add( hfitqual_susy_ldp_1b ) ;

     hfitqual_fit_ldp_2b->Add( hfitqual_znn_ldp_2b ) ;
     hfitqual_fit_ldp_2b->Add( hfitqual_qcd_ldp_2b ) ;
     hfitqual_fit_ldp_2b->Add( hfitqual_ttwj_ldp_2b ) ;
     hfitqual_fit_ldp_2b->Add( hfitqual_susy_ldp_2b ) ;

     hfitqual_fit_ldp_3b->Add( hfitqual_znn_ldp_3b ) ;
     hfitqual_fit_ldp_3b->Add( hfitqual_qcd_ldp_3b ) ;
     hfitqual_fit_ldp_3b->Add( hfitqual_ttwj_ldp_3b ) ;
     hfitqual_fit_ldp_3b->Add( hfitqual_susy_ldp_3b ) ;




     TLegend* legend = new TLegend(0.6,0.45,0.9,0.95) ;

     legend->AddEntry( hfitqual_data_0lep_1b, "data" ) ;
     legend->AddEntry( hfitqual_susy_0lep_1b, "SUSY" ) ;
     legend->AddEntry( hfitqual_ttwj_0lep_1b, "ttwj" ) ;
     legend->AddEntry( hfitqual_qcd_0lep_1b,  "QCD" ) ;
     legend->AddEntry( hfitqual_znn_0lep_1b,  "Znunu" ) ;
     //legend->AddEntry( hfitqual_np,           "Eff PG" ) ;


     if ( doNorm ) {
       hfitqual_data_0lep_1b->SetMaximum( hmax ) ;
       hfitqual_data_1lep_1b->SetMaximum( hmax ) ;
       hfitqual_data_ldp_1b->SetMaximum( hmax ) ;
       hfitqual_data_0lep_2b->SetMaximum( hmax ) ;
       hfitqual_data_1lep_2b->SetMaximum( hmax ) ;
       hfitqual_data_ldp_2b->SetMaximum( hmax ) ;
       hfitqual_data_0lep_3b->SetMaximum( hmax ) ;
       hfitqual_data_1lep_3b->SetMaximum( hmax ) ;
       hfitqual_data_ldp_3b->SetMaximum( hmax ) ;
     } else {
       hfitqual_data_0lep_1b->SetMaximum( hmax*(hfitqual_data_0lep_1b->GetMaximum()) ) ;
       hfitqual_data_1lep_1b->SetMaximum( hmax*(hfitqual_data_1lep_1b->GetMaximum()) ) ;
       hfitqual_data_ldp_1b->SetMaximum( hmax*(hfitqual_data_ldp_1b->GetMaximum()) ) ;
       hfitqual_data_0lep_2b->SetMaximum( hmax*(hfitqual_data_0lep_2b->GetMaximum()) ) ;
       hfitqual_data_1lep_2b->SetMaximum( hmax*(hfitqual_data_1lep_2b->GetMaximum()) ) ;
       hfitqual_data_ldp_2b->SetMaximum( hmax*(hfitqual_data_ldp_2b->GetMaximum()) ) ;
       hfitqual_data_0lep_3b->SetMaximum( hmax*(hfitqual_data_0lep_3b->GetMaximum()) ) ;
       hfitqual_data_1lep_3b->SetMaximum( hmax*(hfitqual_data_1lep_3b->GetMaximum()) ) ;
       hfitqual_data_ldp_3b->SetMaximum( hmax*(hfitqual_data_ldp_3b->GetMaximum()) ) ;
     }



















     printf("\n\n Now doing nuisance parameters.\n\n") ; cout << flush ;

    //----------  Now, do nuisance parameters.

     RooConstVar* npcheck = (RooConstVar*) ws->obj( "mean_sf_qcd_M1_H1_1b" ) ;
     if ( npcheck != 0 ) {

        TH1F* hnp_qcd_1b_val  = new TH1F("hnp_qcd_1b_val" , "Nuisance parameters, qcd, 1b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_qcd_2b_val  = new TH1F("hnp_qcd_2b_val" , "Nuisance parameters, qcd, 2b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_qcd_3b_val  = new TH1F("hnp_qcd_3b_val" , "Nuisance parameters, qcd, 3b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_qcd_1b_pull = new TH1F("hnp_qcd_1b_pull", "Nuisance parameters, qcd, 1b, pull"  , nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_qcd_2b_pull = new TH1F("hnp_qcd_2b_pull", "Nuisance parameters, qcd, 2b, pull"  , nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_qcd_3b_pull = new TH1F("hnp_qcd_3b_pull", "Nuisance parameters, qcd, 3b, pull"  , nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_qcd_1b_nom  = new TH1F("hnp_qcd_1b_nom" , "Nuisance parameters, qcd, 1b, nominal", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_qcd_2b_nom  = new TH1F("hnp_qcd_2b_nom" , "Nuisance parameters, qcd, 2b, nominal", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_qcd_3b_nom  = new TH1F("hnp_qcd_3b_nom" , "Nuisance parameters, qcd, 3b, nominal", nbins, 0.5, nbins+0.5 ) ;

        TH1F* hnp_ttwj_1b_val  = new TH1F("hnp_ttwj_1b_val" , "Nuisance parameters, ttwj, 1b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_ttwj_2b_val  = new TH1F("hnp_ttwj_2b_val" , "Nuisance parameters, ttwj, 2b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_ttwj_3b_val  = new TH1F("hnp_ttwj_3b_val" , "Nuisance parameters, ttwj, 3b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_ttwj_1b_pull = new TH1F("hnp_ttwj_1b_pull", "Nuisance parameters, ttwj, 1b, pull"  , nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_ttwj_2b_pull = new TH1F("hnp_ttwj_2b_pull", "Nuisance parameters, ttwj, 2b, pull"  , nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_ttwj_3b_pull = new TH1F("hnp_ttwj_3b_pull", "Nuisance parameters, ttwj, 3b, pull"  , nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_ttwj_1b_nom  = new TH1F("hnp_ttwj_1b_nom" , "Nuisance parameters, ttwj, 1b, nominal", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_ttwj_2b_nom  = new TH1F("hnp_ttwj_2b_nom" , "Nuisance parameters, ttwj, 2b, nominal", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_ttwj_3b_nom  = new TH1F("hnp_ttwj_3b_nom" , "Nuisance parameters, ttwj, 3b, nominal", nbins, 0.5, nbins+0.5 ) ;



        TH1F* hnp_eff_sf_1b_val  = new TH1F("hnp_eff_sf_1b_val" , "Nuisance parameters, eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_sf_2b_val  = new TH1F("hnp_eff_sf_2b_val" , "Nuisance parameters, eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_sf_3b_val  = new TH1F("hnp_eff_sf_3b_val" , "Nuisance parameters, eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;

        TH1F* hnp_eff_sf_1b_pull  = new TH1F("hnp_eff_sf_1b_pull" , "Nuisance parameters, eff SF, 1b, pull", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_sf_2b_pull  = new TH1F("hnp_eff_sf_2b_pull" , "Nuisance parameters, eff SF, 2b, pull", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_sf_3b_pull  = new TH1F("hnp_eff_sf_3b_pull" , "Nuisance parameters, eff SF, 3b, pull", nbins, 0.5, nbins+0.5 ) ;

        TH1F* hnp_eff_sf_sl_1b_val  = new TH1F("hnp_eff_sf_sl_1b_val" , "Nuisance parameters, eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_sf_sl_2b_val  = new TH1F("hnp_eff_sf_sl_2b_val" , "Nuisance parameters, eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_sf_sl_3b_val  = new TH1F("hnp_eff_sf_sl_3b_val" , "Nuisance parameters, eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;

        TH1F* hnp_eff_sf_sl_1b_pull  = new TH1F("hnp_eff_sf_sl_1b_pull" , "Nuisance parameters, eff SF, 1b, pull", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_sf_sl_2b_pull  = new TH1F("hnp_eff_sf_sl_2b_pull" , "Nuisance parameters, eff SF, 2b, pull", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_sf_sl_3b_pull  = new TH1F("hnp_eff_sf_sl_3b_pull" , "Nuisance parameters, eff SF, 3b, pull", nbins, 0.5, nbins+0.5 ) ;

        TH1F* hnp_eff_sf_ldp_1b_val  = new TH1F("hnp_eff_sf_ldp_1b_val" , "Nuisance parameters, eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_sf_ldp_2b_val  = new TH1F("hnp_eff_sf_ldp_2b_val" , "Nuisance parameters, eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_sf_ldp_3b_val  = new TH1F("hnp_eff_sf_ldp_3b_val" , "Nuisance parameters, eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;

        TH1F* hnp_eff_sf_ldp_1b_pull  = new TH1F("hnp_eff_sf_ldp_1b_pull" , "Nuisance parameters, eff SF, 1b, pull", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_sf_ldp_2b_pull  = new TH1F("hnp_eff_sf_ldp_2b_pull" , "Nuisance parameters, eff SF, 2b, pull", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_sf_ldp_3b_pull  = new TH1F("hnp_eff_sf_ldp_3b_pull" , "Nuisance parameters, eff SF, 3b, pull", nbins, 0.5, nbins+0.5 ) ;



        TH1F* hnp_btageff_sf_1b_val  = new TH1F("hnp_btageff_sf_1b_val" , "Nuisance parameters, btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_btageff_sf_2b_val  = new TH1F("hnp_btageff_sf_2b_val" , "Nuisance parameters, btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_btageff_sf_3b_val  = new TH1F("hnp_btageff_sf_3b_val" , "Nuisance parameters, btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;

        TH1F* hnp_eff_btageff_sf_prod_1b_val  = new TH1F("hnp_eff_btageff_sf_prod_1b_val" , "Nuisance parameters, eff SF * btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_btageff_sf_prod_2b_val  = new TH1F("hnp_eff_btageff_sf_prod_2b_val" , "Nuisance parameters, eff SF * btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_btageff_sf_prod_3b_val  = new TH1F("hnp_eff_btageff_sf_prod_3b_val" , "Nuisance parameters, eff SF * btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;

        TH1F* hnp_btageff_sf_sl_1b_val  = new TH1F("hnp_btageff_sf_sl_1b_val" , "Nuisance parameters, btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_btageff_sf_sl_2b_val  = new TH1F("hnp_btageff_sf_sl_2b_val" , "Nuisance parameters, btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_btageff_sf_sl_3b_val  = new TH1F("hnp_btageff_sf_sl_3b_val" , "Nuisance parameters, btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;

        TH1F* hnp_eff_btageff_sf_sl_prod_1b_val  = new TH1F("hnp_eff_btageff_sf_sl_prod_1b_val" , "Nuisance parameters, eff SF * btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_btageff_sf_sl_prod_2b_val  = new TH1F("hnp_eff_btageff_sf_sl_prod_2b_val" , "Nuisance parameters, eff SF * btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_btageff_sf_sl_prod_3b_val  = new TH1F("hnp_eff_btageff_sf_sl_prod_3b_val" , "Nuisance parameters, eff SF * btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;

        TH1F* hnp_btageff_sf_ldp_1b_val  = new TH1F("hnp_btageff_sf_ldp_1b_val" , "Nuisance parameters, btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_btageff_sf_ldp_2b_val  = new TH1F("hnp_btageff_sf_ldp_2b_val" , "Nuisance parameters, btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_btageff_sf_ldp_3b_val  = new TH1F("hnp_btageff_sf_ldp_3b_val" , "Nuisance parameters, btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;

        TH1F* hnp_eff_btageff_sf_ldp_prod_1b_val  = new TH1F("hnp_eff_btageff_sf_ldp_prod_1b_val" , "Nuisance parameters, eff SF * btag eff SF, 1b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_btageff_sf_ldp_prod_2b_val  = new TH1F("hnp_eff_btageff_sf_ldp_prod_2b_val" , "Nuisance parameters, eff SF * btag eff SF, 2b, values", nbins, 0.5, nbins+0.5 ) ;
        TH1F* hnp_eff_btageff_sf_ldp_prod_3b_val  = new TH1F("hnp_eff_btageff_sf_ldp_prod_3b_val" , "Nuisance parameters, eff SF * btag eff SF, 3b, values", nbins, 0.5, nbins+0.5 ) ;

        TH1F* hnp_prim_eff = new TH1F("hnp_prim_eff", "Nuisance parameters, Efficiency primary Gaussians, pull", 7, 0.5, 7.5 ) ;

        TH1F* hnp_znn = new TH1F("hnp_znn", "Nuisance parameters, Znn, pull", 10+2*nBinsMET, 0.5, 10+2*nBinsMET+0.5 ) ;




        for ( int mbi=0 ; mbi < nBinsMET; mbi ++ ) {
           for ( int hbi=0 ; hbi < nBinsHT; hbi ++ ) {
              if ( ignoreBin[mbi][hbi] ) continue ;
              for ( int bbi=0 ; bbi < nBinsBtag; bbi ++ ) {
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
                 binIndex = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
                 if ( bbi==0 ) hnp_qcd_1b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==1 ) hnp_qcd_2b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==2 ) hnp_qcd_3b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==0 ) hnp_qcd_1b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==1 ) hnp_qcd_2b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==2 ) hnp_qcd_3b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==0 ) hnp_qcd_1b_nom->SetBinContent( binIndex, mean ) ;
                 if ( bbi==1 ) hnp_qcd_2b_nom->SetBinContent( binIndex, mean ) ;
                 if ( bbi==2 ) hnp_qcd_3b_nom->SetBinContent( binIndex, mean ) ;

                 char label[1000] ;
                 sprintf( label, "M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 if ( bbi==0 ) hnp_qcd_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_qcd_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_qcd_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_qcd_1b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_qcd_2b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_qcd_3b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_qcd_1b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_qcd_2b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_qcd_3b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;


              }
           } // hbi
        } // mbi.
        printf("\n\n") ;


        for ( int mbi=0 ; mbi < nBinsMET; mbi ++ ) {
           for ( int hbi=0 ; hbi < nBinsHT; hbi ++ ) {
              if ( ignoreBin[mbi][hbi] ) continue ;
              for ( int bbi=0 ; bbi < nBinsBtag; bbi ++ ) {
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
                 binIndex = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
                 if ( bbi==0 ) hnp_ttwj_1b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==1 ) hnp_ttwj_2b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==2 ) hnp_ttwj_3b_val->SetBinContent( binIndex, npVal ) ;
                 if ( bbi==0 ) hnp_ttwj_1b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==1 ) hnp_ttwj_2b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==2 ) hnp_ttwj_3b_pull->SetBinContent( binIndex, pull ) ;
                 if ( bbi==0 ) hnp_ttwj_1b_nom->SetBinContent( binIndex, mean ) ;
                 if ( bbi==1 ) hnp_ttwj_2b_nom->SetBinContent( binIndex, mean ) ;
                 if ( bbi==2 ) hnp_ttwj_3b_nom->SetBinContent( binIndex, mean ) ;

                 char label[1000] ;
                 sprintf( label, "M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 if ( bbi==0 ) hnp_ttwj_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_ttwj_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_ttwj_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_ttwj_1b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_ttwj_2b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_ttwj_3b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_ttwj_1b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_ttwj_2b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_ttwj_3b_nom-> GetXaxis() -> SetBinLabel( binIndex, label ) ;

              }
           } // hbi
        } // mbi.
        printf("\n\n") ;


        //--- Efficiency primary Gaussians.
        {
            char parName[1000] ;
            RooAbsReal* np(0x0) ;

            sprintf( parName, "eff_sf" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double global_eff_sf(0.) ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
               /// return ;
            } else {
               global_eff_sf = np -> getVal() ;
            }

            sprintf( parName, "btageff_sf" ) ;
            np = (RooAbsReal*) ws->obj( parName ) ;
            double global_btageff_sf = 1. ;
            if ( np == 0x0 ) {
               printf("\n\n *** missing nuisance parameter? %s\n\n", parName ) ;
            } else {
               global_btageff_sf = np -> getVal() ;
            }

            //--- don't know where else to put this...
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

            printf(" Global efficiency scale factors (both Gaussian with mean 0, sigma 1): eff_sf = %6.3f,  btageff_sf = %6.3f\n", global_eff_sf, global_btageff_sf ) ;

            hnp_prim_eff -> SetBinContent(2, global_eff_sf ) ;
            hnp_prim_eff -> SetBinContent(4, global_btageff_sf ) ;
            hnp_prim_eff -> SetBinContent(6, sf_mc_pull ) ;
            TAxis* xaxis = hnp_prim_eff -> GetXaxis() ;
            xaxis->SetBinLabel(2,"Eff SF") ;
            xaxis->SetBinLabel(4,"Btag Eff SF") ;
            xaxis->SetBinLabel(6,"MC SF") ;
            hnp_prim_eff -> SetFillColor(kOrange+1 ) ;

            hnp_prim_eff -> SetLabelSize(0.075,"x") ;
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

            int  n_znnNP_gauss_pars(4) ;
            char znnNP_gauss_par[4][100] = { "knn_1b", "knn_2b", "knn_3b", "sf_ll" } ;


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

            int  n_znnNP_beta_pars = 4 + 2 * nBinsMET ;
            char znnNP_beta_par[20][100] ;
            int zbnpi(0) ;
            sprintf( znnNP_beta_par[zbnpi++], "eff_Zee" ) ;
            sprintf( znnNP_beta_par[zbnpi++], "eff_Zmm" ) ;
            sprintf( znnNP_beta_par[zbnpi++], "pur_Zee" ) ;
            sprintf( znnNP_beta_par[zbnpi++], "pur_Zmm" ) ;
            for ( int i=0; i<nBinsMET; i++ ) {
               sprintf( znnNP_beta_par[zbnpi++], "acc_Zee_M%d", i+1 ) ;
               sprintf( znnNP_beta_par[zbnpi++], "acc_Zmm_M%d", i+1 ) ;
            } // i

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



       //--- bin-by-bin efficiencies

        for ( int mbi=0 ; mbi < nBinsMET; mbi ++ ) {
           for ( int hbi=0 ; hbi < nBinsHT; hbi ++ ) {
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




                 binIndex = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

                 if ( bbi == 0 ) { hnp_eff_sf_1b_val -> SetBinContent( binIndex, bin_eff_sf ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_2b_val -> SetBinContent( binIndex, bin_eff_sf ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_3b_val -> SetBinContent( binIndex, bin_eff_sf ) ; }
                 if ( bbi == 0 ) { hnp_eff_sf_1b_pull -> SetBinContent( binIndex, bin_eff_sf_pull ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_2b_pull -> SetBinContent( binIndex, bin_eff_sf_pull ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_3b_pull -> SetBinContent( binIndex, bin_eff_sf_pull ) ; }
                 if ( bbi == 0 ) { hnp_btageff_sf_1b_val -> SetBinContent( binIndex, btageff_sf ) ; }
                 if ( bbi == 1 ) { hnp_btageff_sf_2b_val -> SetBinContent( binIndex, btageff_sf ) ; }
                 if ( bbi == 2 ) { hnp_btageff_sf_3b_val -> SetBinContent( binIndex, btageff_sf ) ; }
                 if ( bbi == 0 ) { hnp_eff_btageff_sf_prod_1b_val -> SetBinContent( binIndex, prod ) ; }
                 if ( bbi == 1 ) { hnp_eff_btageff_sf_prod_2b_val -> SetBinContent( binIndex, prod ) ; }
                 if ( bbi == 2 ) { hnp_eff_btageff_sf_prod_3b_val -> SetBinContent( binIndex, prod ) ; }

                 if ( bbi == 0 ) { hnp_eff_sf_sl_1b_val -> SetBinContent( binIndex, bin_eff_sf_sl ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_sl_2b_val -> SetBinContent( binIndex, bin_eff_sf_sl ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_sl_3b_val -> SetBinContent( binIndex, bin_eff_sf_sl ) ; }
                 if ( bbi == 0 ) { hnp_eff_sf_sl_1b_pull -> SetBinContent( binIndex, bin_eff_sf_sl_pull ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_sl_2b_pull -> SetBinContent( binIndex, bin_eff_sf_sl_pull ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_sl_3b_pull -> SetBinContent( binIndex, bin_eff_sf_sl_pull ) ; }
                 if ( bbi == 0 ) { hnp_btageff_sf_sl_1b_val -> SetBinContent( binIndex, btageff_sf_sl ) ; }
                 if ( bbi == 1 ) { hnp_btageff_sf_sl_2b_val -> SetBinContent( binIndex, btageff_sf_sl ) ; }
                 if ( bbi == 2 ) { hnp_btageff_sf_sl_3b_val -> SetBinContent( binIndex, btageff_sf_sl ) ; }
                 if ( bbi == 0 ) { hnp_eff_btageff_sf_sl_prod_1b_val -> SetBinContent( binIndex, prod_sl ) ; }
                 if ( bbi == 1 ) { hnp_eff_btageff_sf_sl_prod_2b_val -> SetBinContent( binIndex, prod_sl ) ; }
                 if ( bbi == 2 ) { hnp_eff_btageff_sf_sl_prod_3b_val -> SetBinContent( binIndex, prod_sl ) ; }

                 if ( bbi == 0 ) { hnp_eff_sf_ldp_1b_val -> SetBinContent( binIndex, bin_eff_sf_ldp ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_ldp_2b_val -> SetBinContent( binIndex, bin_eff_sf_ldp ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_ldp_3b_val -> SetBinContent( binIndex, bin_eff_sf_ldp ) ; }
                 if ( bbi == 0 ) { hnp_eff_sf_ldp_1b_pull -> SetBinContent( binIndex, bin_eff_sf_ldp_pull ) ; }
                 if ( bbi == 1 ) { hnp_eff_sf_ldp_2b_pull -> SetBinContent( binIndex, bin_eff_sf_ldp_pull ) ; }
                 if ( bbi == 2 ) { hnp_eff_sf_ldp_3b_pull -> SetBinContent( binIndex, bin_eff_sf_ldp_pull ) ; }
                 if ( bbi == 0 ) { hnp_btageff_sf_ldp_1b_val -> SetBinContent( binIndex, btageff_sf_ldp ) ; }
                 if ( bbi == 1 ) { hnp_btageff_sf_ldp_2b_val -> SetBinContent( binIndex, btageff_sf_ldp ) ; }
                 if ( bbi == 2 ) { hnp_btageff_sf_ldp_3b_val -> SetBinContent( binIndex, btageff_sf_ldp ) ; }
                 if ( bbi == 0 ) { hnp_eff_btageff_sf_ldp_prod_1b_val -> SetBinContent( binIndex, prod_ldp ) ; }
                 if ( bbi == 1 ) { hnp_eff_btageff_sf_ldp_prod_2b_val -> SetBinContent( binIndex, prod_ldp ) ; }
                 if ( bbi == 2 ) { hnp_eff_btageff_sf_ldp_prod_3b_val -> SetBinContent( binIndex, prod_ldp ) ; }

                 printf(" m,h,b %d,%d,%d : eff_sf     = %6.3f,  btageff_sf     = %6.3f,  prod = %6.3f\n", mbi, hbi, bbi, bin_eff_sf, btageff_sf, prod ) ;
                 printf(" m,h,b %d,%d,%d : eff_sf_sl  = %6.3f,  btageff_sf_sl  = %6.3f,  prod = %6.3f\n", mbi, hbi, bbi, bin_eff_sf_sl, btageff_sf_sl, prod_sl ) ;
                 printf(" m,h,b %d,%d,%d : eff_sf_ldp = %6.3f,  btageff_sf_ldp = %6.3f,  prod = %6.3f\n", mbi, hbi, bbi, bin_eff_sf_ldp, btageff_sf_ldp, prod_ldp ) ;

                 char label[1000] ;
                 sprintf( label, "M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                 if ( bbi==0 ) hnp_eff_sf_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_sf_1b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_2b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_3b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_btageff_sf_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_btageff_sf_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_btageff_sf_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_btageff_sf_prod_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_btageff_sf_prod_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_btageff_sf_prod_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;

                 if ( bbi==0 ) hnp_eff_sf_sl_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_sl_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_sl_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_sf_sl_1b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_sl_2b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_sl_3b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_btageff_sf_sl_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_btageff_sf_sl_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_btageff_sf_sl_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_btageff_sf_sl_prod_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_btageff_sf_sl_prod_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_btageff_sf_sl_prod_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;

                 if ( bbi==0 ) hnp_eff_sf_ldp_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_ldp_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_ldp_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_sf_ldp_1b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_sf_ldp_2b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_sf_ldp_3b_pull-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_btageff_sf_ldp_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_btageff_sf_ldp_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_btageff_sf_ldp_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==0 ) hnp_eff_btageff_sf_ldp_prod_1b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==1 ) hnp_eff_btageff_sf_ldp_prod_2b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;
                 if ( bbi==2 ) hnp_eff_btageff_sf_ldp_prod_3b_val-> GetXaxis() -> SetBinLabel( binIndex, label ) ;

              }
           } // hbi
        } // mbi.

        hnp_eff_sf_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_3b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_3b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_1b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_1b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_2b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_2b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_3b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_3b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_1b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_2b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_3b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_3b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_prod_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_prod_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_prod_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_prod_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_prod_3b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_prod_3b_val->GetXaxis()->LabelsOption("v") ;

        hnp_eff_sf_sl_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_sl_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_sl_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_sl_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_sl_3b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_sl_3b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_sl_1b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_sl_1b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_sl_2b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_sl_2b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_sl_3b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_sl_3b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_sl_1b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_sl_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_sl_2b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_sl_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_sl_3b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_sl_3b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_sl_prod_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_sl_prod_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_sl_prod_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_sl_prod_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_sl_prod_3b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_sl_prod_3b_val->GetXaxis()->LabelsOption("v") ;

        hnp_eff_sf_ldp_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_ldp_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_ldp_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_ldp_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_ldp_3b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_ldp_3b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_ldp_1b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_ldp_1b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_ldp_2b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_ldp_2b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_eff_sf_ldp_3b_pull->SetLabelSize(0.055,"x") ;
        hnp_eff_sf_ldp_3b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_ldp_1b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_ldp_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_ldp_2b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_ldp_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_btageff_sf_ldp_3b_val->SetLabelSize(0.055,"x") ;
        hnp_btageff_sf_ldp_3b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_ldp_prod_1b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_ldp_prod_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_ldp_prod_2b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_ldp_prod_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_eff_btageff_sf_ldp_prod_3b_val->SetLabelSize(0.055,"x") ;
        hnp_eff_btageff_sf_ldp_prod_3b_val->GetXaxis()->LabelsOption("v") ;

        hnp_eff_sf_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_3b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_1b_pull->SetFillColor(kOrange+1) ;
        hnp_eff_sf_2b_pull->SetFillColor(kOrange+1) ;
        hnp_eff_sf_3b_pull->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_1b_val->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_2b_val->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_3b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_prod_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_prod_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_prod_3b_val->SetFillColor(kOrange+1) ;

        hnp_eff_sf_sl_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_sl_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_sl_3b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_sl_1b_pull->SetFillColor(kOrange+1) ;
        hnp_eff_sf_sl_2b_pull->SetFillColor(kOrange+1) ;
        hnp_eff_sf_sl_3b_pull->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_sl_1b_val->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_sl_2b_val->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_sl_3b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_sl_prod_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_sl_prod_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_sl_prod_3b_val->SetFillColor(kOrange+1) ;

        hnp_eff_sf_ldp_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_ldp_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_ldp_3b_val->SetFillColor(kOrange+1) ;
        hnp_eff_sf_ldp_1b_pull->SetFillColor(kOrange+1) ;
        hnp_eff_sf_ldp_2b_pull->SetFillColor(kOrange+1) ;
        hnp_eff_sf_ldp_3b_pull->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_ldp_1b_val->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_ldp_2b_val->SetFillColor(kOrange+1) ;
        hnp_btageff_sf_ldp_3b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_ldp_prod_1b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_ldp_prod_2b_val->SetFillColor(kOrange+1) ;
        hnp_eff_btageff_sf_ldp_prod_3b_val->SetFillColor(kOrange+1) ;


        hnp_ttwj_1b_val->SetFillColor(kBlue-9) ;
        hnp_ttwj_2b_val->SetFillColor(kBlue-9) ;
        hnp_ttwj_3b_val->SetFillColor(kBlue-9) ;

        hnp_ttwj_1b_pull->SetFillColor(kBlue-9) ;
        hnp_ttwj_2b_pull->SetFillColor(kBlue-9) ;
        hnp_ttwj_3b_pull->SetFillColor(kBlue-9) ;

        hnp_qcd_1b_val->SetFillColor(2) ;
        hnp_qcd_2b_val->SetFillColor(2) ;
        hnp_qcd_3b_val->SetFillColor(2) ;

        hnp_qcd_1b_pull->SetFillColor(2) ;
        hnp_qcd_2b_pull->SetFillColor(2) ;
        hnp_qcd_3b_pull->SetFillColor(2) ;


        hnp_ttwj_1b_val->SetMaximum(5.) ;
        hnp_ttwj_2b_val->SetMaximum(5.) ;
        hnp_ttwj_3b_val->SetMaximum(5.) ;
        hnp_qcd_1b_val->SetMaximum(5.) ;
        hnp_qcd_2b_val->SetMaximum(5.) ;
        hnp_qcd_3b_val->SetMaximum(5.) ;

        hnp_ttwj_1b_pull->SetMinimum(-2.) ;
        hnp_ttwj_2b_pull->SetMinimum(-2.) ;
        hnp_ttwj_3b_pull->SetMinimum(-2.) ;
        hnp_qcd_1b_pull->SetMinimum(-2.) ;
        hnp_qcd_2b_pull->SetMinimum(-2.) ;
        hnp_qcd_3b_pull->SetMinimum(-2.) ;

        hnp_ttwj_1b_pull->SetMaximum(2.) ;
        hnp_ttwj_2b_pull->SetMaximum(2.) ;
        hnp_ttwj_3b_pull->SetMaximum(2.) ;
        hnp_qcd_1b_pull->SetMaximum(2.) ;
        hnp_qcd_2b_pull->SetMaximum(2.) ;
        hnp_qcd_3b_pull->SetMaximum(2.) ;

        hnp_eff_sf_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_3b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_1b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_2b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_3b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_prod_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_prod_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_prod_3b_val -> SetMinimum(-0.2) ;

        hnp_eff_sf_sl_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_sl_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_sl_3b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_sl_1b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_sl_2b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_sl_3b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_sl_prod_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_sl_prod_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_sl_prod_3b_val -> SetMinimum(-0.2) ;

        hnp_eff_sf_ldp_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_ldp_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_sf_ldp_3b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_ldp_1b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_ldp_2b_val -> SetMinimum(-0.2) ;
        hnp_btageff_sf_ldp_3b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_ldp_prod_1b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_ldp_prod_2b_val -> SetMinimum(-0.2) ;
        hnp_eff_btageff_sf_ldp_prod_3b_val -> SetMinimum(-0.2) ;

        hnp_eff_sf_1b_val -> SetMaximum(2.0) ;
        hnp_eff_sf_2b_val -> SetMaximum(2.0) ;
        hnp_eff_sf_3b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_1b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_2b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_3b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_prod_1b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_prod_2b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_prod_3b_val -> SetMaximum(2.0) ;

        hnp_eff_sf_sl_1b_val -> SetMaximum(2.0) ;
        hnp_eff_sf_sl_2b_val -> SetMaximum(2.0) ;
        hnp_eff_sf_sl_3b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_sl_1b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_sl_2b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_sl_3b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_sl_prod_1b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_sl_prod_2b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_sl_prod_3b_val -> SetMaximum(2.0) ;

        hnp_eff_sf_ldp_1b_val -> SetMaximum(2.0) ;
        hnp_eff_sf_ldp_2b_val -> SetMaximum(2.0) ;
        hnp_eff_sf_ldp_3b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_ldp_1b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_ldp_2b_val -> SetMaximum(2.0) ;
        hnp_btageff_sf_ldp_3b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_ldp_prod_1b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_ldp_prod_2b_val -> SetMaximum(2.0) ;
        hnp_eff_btageff_sf_ldp_prod_3b_val -> SetMaximum(2.0) ;

        hnp_eff_sf_1b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_2b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_3b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_1b_pull -> SetMaximum( 2.0) ;
        hnp_eff_sf_2b_pull -> SetMaximum( 2.0) ;
        hnp_eff_sf_3b_pull -> SetMaximum( 2.0) ;

        hnp_eff_sf_sl_1b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_sl_2b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_sl_3b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_sl_1b_pull -> SetMaximum( 2.0) ;
        hnp_eff_sf_sl_2b_pull -> SetMaximum( 2.0) ;
        hnp_eff_sf_sl_3b_pull -> SetMaximum( 2.0) ;

        hnp_eff_sf_ldp_1b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_ldp_2b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_ldp_3b_pull -> SetMinimum(-2.0) ;
        hnp_eff_sf_ldp_1b_pull -> SetMaximum( 2.0) ;
        hnp_eff_sf_ldp_2b_pull -> SetMaximum( 2.0) ;
        hnp_eff_sf_ldp_3b_pull -> SetMaximum( 2.0) ;

        hnp_ttwj_1b_val->SetLabelSize(0.055,"x") ;
        hnp_ttwj_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_ttwj_2b_val->SetLabelSize(0.055,"x") ;
        hnp_ttwj_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_ttwj_3b_val->SetLabelSize(0.055,"x") ;
        hnp_ttwj_3b_val->GetXaxis()->LabelsOption("v") ;

        hnp_qcd_1b_val->SetLabelSize(0.055,"x") ;
        hnp_qcd_1b_val->GetXaxis()->LabelsOption("v") ;
        hnp_qcd_2b_val->SetLabelSize(0.055,"x") ;
        hnp_qcd_2b_val->GetXaxis()->LabelsOption("v") ;
        hnp_qcd_3b_val->SetLabelSize(0.055,"x") ;
        hnp_qcd_3b_val->GetXaxis()->LabelsOption("v") ;

        hnp_ttwj_1b_pull->SetLabelSize(0.055,"x") ;
        hnp_ttwj_1b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_ttwj_2b_pull->SetLabelSize(0.055,"x") ;
        hnp_ttwj_2b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_ttwj_3b_pull->SetLabelSize(0.055,"x") ;
        hnp_ttwj_3b_pull->GetXaxis()->LabelsOption("v") ;

        hnp_qcd_1b_pull->SetLabelSize(0.055,"x") ;
        hnp_qcd_1b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_qcd_2b_pull->SetLabelSize(0.055,"x") ;
        hnp_qcd_2b_pull->GetXaxis()->LabelsOption("v") ;
        hnp_qcd_3b_pull->SetLabelSize(0.055,"x") ;
        hnp_qcd_3b_pull->GetXaxis()->LabelsOption("v") ;




        //=== Efficiency and Znn nuisance par pulls ============


        TCanvas* cnp2 = (TCanvas*) gDirectory->FindObject("cnp2") ;
        if ( cnp2 == 0x0 ) {
           cnp2 = new TCanvas("cnp2","RA2b efficiency and Znn nuisance par pull", 850, 1000 ) ;
        }
        cnp2->Divide(3,4) ;

        cnp2->cd(1) ;
        hnp_prim_eff->Draw() ;
        chi2 = addChi2FromPullHist( hnp_prim_eff ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;


        cnp2->cd(2) ;
        hnp_znn->Draw() ;
        chi2 = addChi2FromPullHist( hnp_znn ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;



        cnp2->cd(4) ;
        hnp_eff_sf_1b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_1b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp2->cd(5) ;
        hnp_eff_sf_2b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_2b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp2->cd(6) ;
        hnp_eff_sf_3b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_3b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;



        cnp2->cd(7) ;
        hnp_eff_sf_sl_1b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_sl_1b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp2->cd(8) ;
        hnp_eff_sf_sl_2b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_sl_2b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp2->cd(9) ;
        hnp_eff_sf_sl_3b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_sl_3b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;



        cnp2->cd(10) ;
        hnp_eff_sf_ldp_1b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_ldp_1b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp2->cd(11) ;
        hnp_eff_sf_ldp_2b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_ldp_2b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp2->cd(12) ;
        hnp_eff_sf_ldp_3b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_eff_sf_ldp_3b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;






       //==== TTwj and QCD scale factors ========================

        gStyle->SetPadGridY(1) ;
        TCanvas* cnp = (TCanvas*) gDirectory->FindObject("cnp") ;
        if ( cnp == 0x0 ) {
           cnp = new TCanvas("cnp","RA2b nuisance pars", 850, 1000 ) ;
        }
        cnp->Divide(3,4) ;

      //---
        cnp->cd(1) ;
        hnp_ttwj_1b_val->Draw() ;
        hnp_ttwj_1b_nom->Draw("same") ;

        cnp->cd(2) ;
        hnp_ttwj_2b_val->Draw() ;
        hnp_ttwj_2b_nom->Draw("same") ;

        cnp->cd(3) ;
        hnp_ttwj_3b_val->Draw() ;
        hnp_ttwj_3b_nom->Draw("same") ;


      //---
        cnp->cd(4) ;
        hnp_ttwj_1b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_ttwj_1b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;


        cnp->cd(5) ;
        hnp_ttwj_2b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_ttwj_2b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp->cd(6) ;
        hnp_ttwj_3b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_ttwj_3b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;


      //---
        cnp->cd(7) ;
        hnp_qcd_1b_val->Draw() ;
        hnp_qcd_1b_nom->Draw("same") ;

        cnp->cd(8) ;
        hnp_qcd_2b_val->Draw() ;
        hnp_qcd_2b_nom->Draw("same") ;

        cnp->cd(9) ;
        hnp_qcd_3b_val->Draw() ;
        hnp_qcd_3b_nom->Draw("same") ;


      //---
        cnp->cd(10) ;
        hnp_qcd_1b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_qcd_1b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp->cd(11) ;
        hnp_qcd_2b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_qcd_2b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;

        cnp->cd(12) ;
        hnp_qcd_3b_pull->Draw() ;
        chi2 = addChi2FromPullHist( hnp_qcd_3b_pull ) ;
        globalChi2 += chi2 ;
        npChi2 += chi2 ;







     } else {
        printf("\n\n *** Skipping nuisance parameters.\n\n") ;
     }

   //====== Observables ==================================


     TCanvas* cfitqual = (TCanvas*) gDirectory->FindObject("cfitqual") ;
     if ( cfitqual == 0x0 ) {
        cfitqual = new TCanvas("cfitqual","RA2b fit quality", 850, 1000 ) ;
     }

     cfitqual->Divide(3,4);

     gPad->SetTicks(1,0) ;

     hfitqual_data_0lep_1b->SetLabelSize(0.055,"x") ;
     hfitqual_data_0lep_1b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_1lep_1b->SetLabelSize(0.055,"x") ;
     hfitqual_data_1lep_1b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_ldp_1b->SetLabelSize(0.055,"x") ;
     hfitqual_data_ldp_1b->GetXaxis()->LabelsOption("v") ;
     
     hfitqual_data_0lep_2b->SetLabelSize(0.055,"x") ;
     hfitqual_data_0lep_2b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_1lep_2b->SetLabelSize(0.055,"x") ;
     hfitqual_data_1lep_2b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_ldp_2b->SetLabelSize(0.055,"x") ;
     hfitqual_data_ldp_2b->GetXaxis()->LabelsOption("v") ;
     
     hfitqual_data_0lep_3b->SetLabelSize(0.055,"x") ;
     hfitqual_data_0lep_3b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_1lep_3b->SetLabelSize(0.055,"x") ;
     hfitqual_data_1lep_3b->GetXaxis()->LabelsOption("v") ;
     hfitqual_data_ldp_3b->SetLabelSize(0.055,"x") ;
     hfitqual_data_ldp_3b->GetXaxis()->LabelsOption("v") ;
     
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
     
     cfitqual->cd(3);
     hfitqual_data_0lep_3b->Draw("histpe") ;
     hfitqual_fit_0lep_3b->Draw("same") ;
     hfitqual_data_0lep_3b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_0lep_3b, hfitqual_model_0lep_3b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;
     


     cfitqual->cd(4);
     hfitqual_data_1lep_1b->Draw("histpe") ;
     hfitqual_fit_1lep_1b->Draw("same") ;
     hfitqual_data_1lep_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_1lep_1b, hfitqual_model_1lep_1b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;
     
     cfitqual->cd(5);
     hfitqual_data_1lep_2b->Draw("histpe") ;
     hfitqual_fit_1lep_2b->Draw("same") ;
     hfitqual_data_1lep_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_1lep_2b, hfitqual_model_1lep_2b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;
     
     cfitqual->cd(6);
     hfitqual_data_1lep_3b->Draw("histpe") ;
     hfitqual_fit_1lep_3b->Draw("same") ;
     hfitqual_data_1lep_3b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_1lep_3b, hfitqual_model_1lep_3b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;




     
     cfitqual->cd(7);
     hfitqual_data_ldp_1b->Draw("histpe") ;
     hfitqual_fit_ldp_1b->Draw("same") ;
     hfitqual_data_ldp_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_ldp_1b, hfitqual_model_ldp_1b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;

     cfitqual->cd(8);
     hfitqual_data_ldp_2b->Draw("histpe") ;
     hfitqual_fit_ldp_2b->Draw("same") ;
     hfitqual_data_ldp_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_ldp_2b, hfitqual_model_ldp_2b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;

     cfitqual->cd(9);
     hfitqual_data_ldp_3b->Draw("histpe") ;
     hfitqual_fit_ldp_3b->Draw("same") ;
     hfitqual_data_ldp_3b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_ldp_3b, hfitqual_model_ldp_3b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;




     cfitqual->cd(10) ;
     hfitqual_data_zee_1b->Draw("histpe") ;
     hfitqual_fit_zee_1b->Draw("same") ;
     hfitqual_data_zee_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_zee_1b, hfitqual_fit_zee_1b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;

     cfitqual->cd(11) ;
     hfitqual_data_zmm_1b->Draw("histpe") ;
     hfitqual_fit_zmm_1b->Draw("same") ;
     hfitqual_data_zmm_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     chi2 = addChi2FromObs( hfitqual_data_zmm_1b, hfitqual_fit_zmm_1b ) ;
     globalChi2 += chi2 ;
     obsChi2 += chi2 ;


     cfitqual->cd(12);
     legend->Draw() ;

     cfitqual->Update() ;



     cfitqual->cd(12) ;
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








     TString histfile(wsfile) ;
     histfile.ReplaceAll(".root","-fitqual.root") ;
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









