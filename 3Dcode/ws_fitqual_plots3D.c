
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
#include "RooStats/ModelConfig.h"
#include "TRegexp.h"

#include <string.h>
#include <sstream>

  using namespace RooFit;
  using namespace RooStats;

  void saveHist(const char* filename, const char* pat) ;

  //------
  //
  // Note: If the 2nd argument (mu_susy_sig_val) is negative, mu_susy_sig will be floated
  //        in the fit.  If it's zero or positive, mu_susy_sig will be fixed to that value
  //        in the fit.
  //
  //------

   void ws_fitqual_plots3D( const char* wsfile = "ws-met5-ht5-v2b.root", double mu_susy_sig_val = 0., bool doNorm = false, int nBinsMET=5, int nBinsHT=5 ) {


     // hardcode here the number of bins of the analysis
     int nBinsBtag = 3 ;


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


     double hmax = 1.25 ;

     TFile* wstf = new TFile( wsfile ) ;

     RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
     ws->Print() ;

     ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;

     printf("\n\n\n  ===== SbModel ====================\n\n") ;
     modelConfig->Print() ;



     RooDataSet* rds = (RooDataSet*) ws->obj( "ra2b_observed_rds" ) ;
     printf("\n\n\n  ===== RooDataSet ====================\n\n") ;

     rds->Print() ;
     rds->printMultiline(cout, 1, kTRUE, "") ;


     RooAbsPdf* likelihood = modelConfig->GetPdf() ;


     // the poi is mu_susy_M1_H1_1b, maybe later we'll want to translate it somehow to the total signal yield

     RooRealVar* rrv_mu_susy_sig = ws->var("mu_susy_M1_H1_1b") ;
     if ( rrv_mu_susy_sig == 0x0 ) {
       printf("\n\n\n *** can't find mu_susy_M1_H1_1b in workspace.  Quitting.\n\n\n") ;
       return ;
     } else {
       printf(" current value is : %8.3f\n", rrv_mu_susy_sig->getVal() ) ; cout << flush ;
       rrv_mu_susy_sig->setConstant(kFALSE) ;
     }


     printf("\n\n\n  ===== Doing a fit with SUSY component floating ====================\n\n") ;

     RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true) ) ;
     double logLikelihoodSusyFloat = fitResult->minNll() ;
     
     double logLikelihoodSusyFixed(0.) ;
     double testStatVal(-1.) ;
     if ( mu_susy_sig_val >= 0. ) {
       printf("\n\n\n  ===== Doing a fit with SUSY fixed ====================\n\n") ;
       printf(" fixing mu_susy_sig to %8.2f.\n", mu_susy_sig_val ) ;
       rrv_mu_susy_sig->setVal( mu_susy_sig_val ) ;
       rrv_mu_susy_sig->setConstant(kTRUE) ;
       
       fitResult = likelihood->fitTo( *rds, Save(true) ) ;
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

     printf("\n\n") ;


     double AttwjVal[nBinsMET][nBinsHT][nBinsBtag] ;
     double AqcdVal[nBinsMET][nBinsHT][nBinsBtag] ;

     for ( int i = 0 ; i < nBinsMET ; i++ ) {
       for ( int j = 0 ; j < nBinsHT ; j++ ) {
	 for ( int k = 0 ; k < nBinsBtag ; k++ ) {

	   TString ttString  = "mu_ttwj" ;
	   TString qcdString = "mu_qcd" ;

	   ttString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	   qcdString += sMbins[i]+sHbins[j]+sBbins[k] ;

	   TObject* ttwj_obj = ws->obj(ttString) ;
	   TObject* qcd_obj  = ws->obj(qcdString) ;
	   
	   AttwjVal[i][j][k] = ((RooRealVar*) ttwj_obj)->getVal() ;
	   AqcdVal[i][j][k]  = ((RooRealVar*) qcd_obj)->getVal() ;

	 }
       }
     }


     //--- unpack observables.

     int dataN_0lep[nBinsMET][nBinsHT][nBinsBtag] ;
     int dataN_1lep[nBinsMET][nBinsHT][nBinsBtag] ;
     int dataN_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
     int dataN_Zee[nBinsMET][nBinsHT] ;
     int dataN_Zmm[nBinsMET][nBinsHT] ;


     const RooArgSet* dsras = rds->get() ;
     TIterator* obsIter = dsras->createIterator() ;

     while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {

       for ( int i = 0 ; i < nBinsMET ; i++ ) {
	 for ( int j = 0 ; j < nBinsHT ; j++ ) {
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

     TH1F* hfitqual_data_0lep_2b = new TH1F("hfitqual_data_0lep_2b", "0 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_0lep_2b = new TH1F("hfitqual_susy_0lep_2b", "0 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_0lep_2b = new TH1F("hfitqual_ttwj_0lep_2b", "0 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_0lep_2b  = new TH1F("hfitqual_qcd_0lep_2b" , "0 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_0lep_2b  = new TH1F("hfitqual_znn_0lep_2b" , "0 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     TH1F* hfitqual_data_0lep_3b = new TH1F("hfitqual_data_0lep_3b", "0 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_0lep_3b = new TH1F("hfitqual_susy_0lep_3b", "0 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_0lep_3b = new TH1F("hfitqual_ttwj_0lep_3b", "0 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_0lep_3b  = new TH1F("hfitqual_qcd_0lep_3b" , "0 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_0lep_3b  = new TH1F("hfitqual_znn_0lep_3b" , "0 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;



     TH1F* hfitqual_data_1lep_1b = new TH1F("hfitqual_data_1lep_1b", "1 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_1lep_1b = new TH1F("hfitqual_susy_1lep_1b", "1 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_1lep_1b = new TH1F("hfitqual_ttwj_1lep_1b", "1 Lep, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_1lep_1b  = new TH1F("hfitqual_qcd_1lep_1b" , "1 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_1lep_1b  = new TH1F("hfitqual_znn_1lep_1b" , "1 Lep, 1 btag" , nbins, 0.5, nbins+0.5 ) ;

     TH1F* hfitqual_data_1lep_2b = new TH1F("hfitqual_data_1lep_2b", "1 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_1lep_2b = new TH1F("hfitqual_susy_1lep_2b", "1 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_1lep_2b = new TH1F("hfitqual_ttwj_1lep_2b", "1 Lep, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_1lep_2b  = new TH1F("hfitqual_qcd_1lep_2b" , "1 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_1lep_2b  = new TH1F("hfitqual_znn_1lep_2b" , "1 Lep, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     TH1F* hfitqual_data_1lep_3b = new TH1F("hfitqual_data_1lep_3b", "1 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_1lep_3b = new TH1F("hfitqual_susy_1lep_3b", "1 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_1lep_3b = new TH1F("hfitqual_ttwj_1lep_3b", "1 Lep, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_1lep_3b  = new TH1F("hfitqual_qcd_1lep_3b" , "1 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_1lep_3b  = new TH1F("hfitqual_znn_1lep_3b" , "1 Lep, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;



     TH1F* hfitqual_data_ldp_1b = new TH1F("hfitqual_data_ldp_1b", "LDP, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_ldp_1b = new TH1F("hfitqual_susy_ldp_1b", "LDP, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_ldp_1b = new TH1F("hfitqual_ttwj_ldp_1b", "LDP, 1 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_ldp_1b  = new TH1F("hfitqual_qcd_ldp_1b" , "LDP, 1 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_ldp_1b  = new TH1F("hfitqual_znn_ldp_1b" , "LDP, 1 btag" , nbins, 0.5, nbins+0.5 ) ;

     TH1F* hfitqual_data_ldp_2b = new TH1F("hfitqual_data_ldp_2b", "LDP, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_ldp_2b = new TH1F("hfitqual_susy_ldp_2b", "LDP, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_ldp_2b = new TH1F("hfitqual_ttwj_ldp_2b", "LDP, 2 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_ldp_2b  = new TH1F("hfitqual_qcd_ldp_2b" , "LDP, 2 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_ldp_2b  = new TH1F("hfitqual_znn_ldp_2b" , "LDP, 2 btag" , nbins, 0.5, nbins+0.5 ) ;

     TH1F* hfitqual_data_ldp_3b = new TH1F("hfitqual_data_ldp_3b", "LDP, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_susy_ldp_3b = new TH1F("hfitqual_susy_ldp_3b", "LDP, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_ttwj_ldp_3b = new TH1F("hfitqual_ttwj_ldp_3b", "LDP, >=3 btag", nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_qcd_ldp_3b  = new TH1F("hfitqual_qcd_ldp_3b" , "LDP, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;
     TH1F* hfitqual_znn_ldp_3b  = new TH1F("hfitqual_znn_ldp_3b" , "LDP, >=3 btag" , nbins, 0.5, nbins+0.5 ) ;


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

     double dataVal(0.) ;
     double dataErr(0.) ;
     double ttwjVal(0.) ;
     double qcdVal(0.) ;
     double znnVal(0.) ;
     double susyVal(0.) ;
     double lhtotalVal(0.) ;
     
     double eff_sf(0.) ;
     double eff_sf_sl(0.) ;
     double eff_sf_ldp(0.) ;

     double sf_mc(0.) ;


     // loop through all the bins of the analysis

     for ( int i = 0 ; i < nBinsMET ; i++ ) {
       for ( int j = 0 ; j < nBinsHT ; j++ ) {
	 for ( int k = 0 ; k < nBinsBtag ; k++ ) {
	   
	   /// binIndex = 2 + ( nBinsHT*i + j + k ) + k*(nBinsMET*nBinsHT) ;

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

	   dataVal = dataN_0lep[i][j][k];
	   eff_sf = ((RooFormulaVar*) ws->obj(EffSfString)) -> getVal() ;
	   susyVal = eff_sf * ( ((RooRealVar*) ws->obj(MuSusyString)) -> getVal() ) ;
	   znnVal = ((RooRealVar*) ws->obj(ZnnString))  -> getVal() ;
	   ttwjVal = AttwjVal[i][j][k] ;
	   qcdVal  = AqcdVal[i][j][k] ;
	   lhtotalVal = ttwjVal + qcdVal + znnVal + susyVal ;

	   dataErr = sqrt(dataVal) ;

	   // skip for now..
	   //-- get values for later 
	   //susyVal = ((RooRealVar*) ws->obj(MuSusyString)) -> getVal() ;
	   //susyErr = ((RooRealVar*) ws->obj(MuSusyString)) -> getError() ;
	   //znnVal  = ((RooRealVar*) ws->obj(ZnnString))  -> getVal() ;
	   //znnErr  = ((RooRealVar*) ws->obj(ZnnString))  -> getError() ;

	   cout << "\n" << endl ;
	   cout << "binIndex = " << binIndex << endl ;
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
	      xaxis_0lep_1b->SetBinLabel(binIndex, binLabel ) ;
	      hfitqual_data_0lep_1b -> SetBinContent( binIndex, dataVal ) ;
	      hfitqual_data_0lep_1b -> SetBinError( binIndex, dataErr ) ;
	      hfitqual_susy_0lep_1b -> SetBinContent( binIndex, susyVal ) ;
	      hfitqual_ttwj_0lep_1b -> SetBinContent( binIndex, ttwjVal ) ;
	      hfitqual_qcd_0lep_1b  -> SetBinContent( binIndex, qcdVal ) ;
	      hfitqual_znn_0lep_1b  -> SetBinContent( binIndex, znnVal ) ;
           } else if ( k == 1 ) {
	      xaxis_0lep_2b->SetBinLabel(binIndex, binLabel ) ;
	      hfitqual_data_0lep_2b -> SetBinContent( binIndex, dataVal ) ;
	      hfitqual_data_0lep_2b -> SetBinError( binIndex, dataErr ) ;
	      hfitqual_susy_0lep_2b -> SetBinContent( binIndex, susyVal ) ;
	      hfitqual_ttwj_0lep_2b -> SetBinContent( binIndex, ttwjVal ) ;
	      hfitqual_qcd_0lep_2b  -> SetBinContent( binIndex, qcdVal ) ;
	      hfitqual_znn_0lep_2b  -> SetBinContent( binIndex, znnVal ) ;
           } else if ( k == 2 ) {
	      xaxis_0lep_3b->SetBinLabel(binIndex, binLabel ) ;
	      hfitqual_data_0lep_3b -> SetBinContent( binIndex, dataVal ) ;
	      hfitqual_data_0lep_3b -> SetBinError( binIndex, dataErr ) ;
	      hfitqual_susy_0lep_3b -> SetBinContent( binIndex, susyVal ) ;
	      hfitqual_ttwj_0lep_3b -> SetBinContent( binIndex, ttwjVal ) ;
	      hfitqual_qcd_0lep_3b  -> SetBinContent( binIndex, qcdVal ) ;
	      hfitqual_znn_0lep_3b  -> SetBinContent( binIndex, znnVal ) ;
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
	   eff_sf_sl = ((RooFormulaVar*) ws->obj(EffSfSlString)) -> getVal() ;
	   susyVal = eff_sf_sl * ( ((RooRealVar*) ws->obj(MuSusySlString)) -> getVal() ) ;
	   ttwjVal = ((RooRealVar*) ws->obj(MuTtwjSlString)) -> getVal() ; ;
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
           } else if ( k == 1 ) {
	      xaxis_1lep_2b->SetBinLabel(binIndex, binLabel ) ;
	      hfitqual_data_1lep_2b -> SetBinContent( binIndex, dataVal ) ;
	      hfitqual_data_1lep_2b -> SetBinError( binIndex, dataErr ) ;
	      hfitqual_susy_1lep_2b -> SetBinContent( binIndex, susyVal ) ;
	      hfitqual_ttwj_1lep_2b -> SetBinContent( binIndex, ttwjVal ) ;
	      hfitqual_qcd_1lep_2b  -> SetBinContent( binIndex, qcdVal ) ;
	      hfitqual_znn_1lep_2b  -> SetBinContent( binIndex, znnVal ) ;
           } else if ( k == 2 ) {
	      xaxis_1lep_3b->SetBinLabel(binIndex, binLabel ) ;
	      hfitqual_data_1lep_3b -> SetBinContent( binIndex, dataVal ) ;
	      hfitqual_data_1lep_3b -> SetBinError( binIndex, dataErr ) ;
	      hfitqual_susy_1lep_3b -> SetBinContent( binIndex, susyVal ) ;
	      hfitqual_ttwj_1lep_3b -> SetBinContent( binIndex, ttwjVal ) ;
	      hfitqual_qcd_1lep_3b  -> SetBinContent( binIndex, qcdVal ) ;
	      hfitqual_znn_1lep_3b  -> SetBinContent( binIndex, znnVal ) ;
           }
	   




	   //+++++ ldp histograms:

	   binLabel = "ldp" ;
	   binLabel += sMbins[i]+sHbins[j]+sBbins[k] ;

	   TString EffSfLdpString  = "eff_sf_ldp" ;
	   TString MuSusyLdpString = "mu_susy_ldp" ;
	   TString MuTtwjLdpString = "mu_ttwj_ldp" ;
	   TString MuWJmcLdpString = "mu_WJmc" ;
	   TString MuQcdLdpString  = "mu_qcd_ldp" ;

	   EffSfLdpString  += sMbins[i]+sHbins[j]+sBbins[k] ;
	   MuSusyLdpString += sMbins[i]+sHbins[j]+sBbins[k] ;
	   MuTtwjLdpString += sMbins[i]+sHbins[j]+sBbins[k] ;
	   MuWJmcLdpString += sMbins[i]+sHbins[j]+sBbins[k] ;
	   MuQcdLdpString  += sMbins[i]+sHbins[j]+sBbins[k] ;

	   dataVal = dataN_ldp[i][j][k];
	   eff_sf_ldp = ((RooFormulaVar*) ws->obj(EffSfLdpString)) -> getVal() ;
	   susyVal = eff_sf_ldp * ( ((RooRealVar*) ws->obj(MuSusyLdpString)) -> getVal() ) ;
	   sf_mc = ((RooFormulaVar*) ws->obj("sf_mc")) -> getVal() ;
	   znnVal = 0. ;
	   ttwjVal = eff_sf_ldp * sf_mc * ((((RooRealVar*) ws->obj(MuTtwjLdpString))-> getVal() ) + (((RooRealVar*) ws->obj(MuWJmcLdpString))-> getVal() )  ) ; ;
	   qcdVal  = ((RooRealVar*) ws->obj(MuQcdLdpString))-> getVal() ;
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
	      xaxis_ldp_1b->SetBinLabel(binIndex, binLabel ) ;
	      hfitqual_data_ldp_1b -> SetBinContent( binIndex, dataVal ) ;
	      hfitqual_data_ldp_1b -> SetBinError( binIndex, dataErr ) ;
	      hfitqual_susy_ldp_1b -> SetBinContent( binIndex, susyVal ) ;
	      hfitqual_ttwj_ldp_1b -> SetBinContent( binIndex, ttwjVal ) ;
	      hfitqual_qcd_ldp_1b  -> SetBinContent( binIndex, qcdVal ) ;
	      hfitqual_znn_ldp_1b  -> SetBinContent( binIndex, znnVal ) ;
           } else if ( k == 1 ) {
	      xaxis_ldp_2b->SetBinLabel(binIndex, binLabel ) ;
	      hfitqual_data_ldp_2b -> SetBinContent( binIndex, dataVal ) ;
	      hfitqual_data_ldp_2b -> SetBinError( binIndex, dataErr ) ;
	      hfitqual_susy_ldp_2b -> SetBinContent( binIndex, susyVal ) ;
	      hfitqual_ttwj_ldp_2b -> SetBinContent( binIndex, ttwjVal ) ;
	      hfitqual_qcd_ldp_2b  -> SetBinContent( binIndex, qcdVal ) ;
	      hfitqual_znn_ldp_2b  -> SetBinContent( binIndex, znnVal ) ;
           } else if ( k == 2 ) {
	      xaxis_ldp_3b->SetBinLabel(binIndex, binLabel ) ;
	      hfitqual_data_ldp_3b -> SetBinContent( binIndex, dataVal ) ;
	      hfitqual_data_ldp_3b -> SetBinError( binIndex, dataErr ) ;
	      hfitqual_susy_ldp_3b -> SetBinContent( binIndex, susyVal ) ;
	      hfitqual_ttwj_ldp_3b -> SetBinContent( binIndex, ttwjVal ) ;
	      hfitqual_qcd_ldp_3b  -> SetBinContent( binIndex, qcdVal ) ;
	      hfitqual_znn_ldp_3b  -> SetBinContent( binIndex, znnVal ) ;
           }




           if ( k == 0 ) {

	   //+++++ Zee histogram:

	      binLabel = "Zee" ;
	      binLabel += sMbins[i]+sHbins[j]+sBbins[k] ;
             

              TString nZeeModelString = "n_ee" ;
             
              nZeeModelString += sMbins[i]+sHbins[j] ;

              printf("\n Grabbing value of this from ws: %s\n", nZeeModelString.Data() ) ; cout << flush ;
              lhtotalVal = ((RooRealVar*) ws->obj(nZeeModelString)) -> getVal() ;
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



	   //+++++ Zmm histogram:

	      binLabel = "Zmm" ;
	      binLabel += sMbins[i]+sHbins[j]+sBbins[k] ;
             

              TString nZmmModelString = "n_mm" ;
             
              nZmmModelString += sMbins[i]+sHbins[j] ;

              printf("\n Grabbing value of this from ws: %s\n", nZmmModelString.Data() ) ; cout << flush ;
              lhtotalVal = ((RooRealVar*) ws->obj(nZmmModelString)) -> getVal() ;
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



     //-- Eff sf ---------

     float eff_sf_prim = ((RooRealVar*) ws->obj("eff_sf_prim")) -> getVal() ;

     hfitqual_np -> SetBinContent( 1, eff_sf_prim ) ;


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




     TLegend* legend = new TLegend(0.4,0.35,0.7,0.85) ;

     legend->AddEntry( hfitqual_data_0lep_1b, "data" ) ;
     legend->AddEntry( hfitqual_susy_0lep_1b, "SUSY" ) ;
     legend->AddEntry( hfitqual_ttwj_0lep_1b, "ttwj" ) ;
     legend->AddEntry( hfitqual_qcd_0lep_1b,  "QCD" ) ;
     legend->AddEntry( hfitqual_znn_0lep_1b,  "Znunu" ) ;
     legend->AddEntry( hfitqual_np,           "Eff PG" ) ;


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


     TCanvas* cfitqual = new TCanvas("cfitqual","RA2b fit quality", 850, 1000 ) ;

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

     char histfile[10000] ;
     sprintf( histfile, "fitqual-hists-%s", wsfile ) ;
     saveHist( histfile,"h*") ;

     cfitqual->cd(1);
     hfitqual_data_0lep_1b->Draw("histpe") ;
     hfitqual_fit_0lep_1b->Draw("same") ;
     hfitqual_data_0lep_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     cfitqual->cd(2);
     hfitqual_data_0lep_2b->Draw("histpe") ;
     hfitqual_fit_0lep_2b->Draw("same") ;
     hfitqual_data_0lep_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     cfitqual->cd(3);
     hfitqual_data_0lep_3b->Draw("histpe") ;
     hfitqual_fit_0lep_3b->Draw("same") ;
     hfitqual_data_0lep_3b->Draw("same") ;
     gPad->SetGridy(1) ;
     


     cfitqual->cd(4);
     hfitqual_data_1lep_1b->Draw("histpe") ;
     hfitqual_fit_1lep_1b->Draw("same") ;
     hfitqual_data_1lep_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     cfitqual->cd(5);
     hfitqual_data_1lep_2b->Draw("histpe") ;
     hfitqual_fit_1lep_2b->Draw("same") ;
     hfitqual_data_1lep_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     cfitqual->cd(6);
     hfitqual_data_1lep_3b->Draw("histpe") ;
     hfitqual_fit_1lep_3b->Draw("same") ;
     hfitqual_data_1lep_3b->Draw("same") ;
     gPad->SetGridy(1) ;




     
     cfitqual->cd(7);
     hfitqual_data_ldp_1b->Draw("histpe") ;
     hfitqual_fit_ldp_1b->Draw("same") ;
     hfitqual_data_ldp_1b->Draw("same") ;
     gPad->SetGridy(1) ;

     cfitqual->cd(8);
     hfitqual_data_ldp_2b->Draw("histpe") ;
     hfitqual_fit_ldp_2b->Draw("same") ;
     hfitqual_data_ldp_2b->Draw("same") ;
     gPad->SetGridy(1) ;

     cfitqual->cd(9);
     hfitqual_data_ldp_3b->Draw("histpe") ;
     hfitqual_fit_ldp_3b->Draw("same") ;
     hfitqual_data_ldp_3b->Draw("same") ;
     gPad->SetGridy(1) ;




     cfitqual->cd(10) ;
     hfitqual_data_zee_1b->Draw("histpe") ;
     hfitqual_fit_zee_1b->Draw("same") ;
     hfitqual_data_zee_1b->Draw("same") ;
     gPad->SetGridy(1) ;

     cfitqual->cd(11) ;
     hfitqual_data_zmm_1b->Draw("histpe") ;
     hfitqual_fit_zmm_1b->Draw("same") ;
     hfitqual_data_zmm_1b->Draw("same") ;
     gPad->SetGridy(1) ;


     cfitqual->cd(12);
     legend->Draw() ;

     cfitqual->Update() ;



     //--- Efficiency scale factor primary Gaussian value.

     hfitqual_np->SetMinimum(-5.) ;
     hfitqual_np->SetMaximum( 5.) ;
     
     hfitqual_np->SetNdivisions(101,"x") ;
     hfitqual_np->SetNdivisions(101,"y") ;
     hfitqual_np->SetLabelOffset(99,"y") ;
     

     TPad* tp = new TPad("tp","tp",0.09,0.,0.18,1.0) ;

     tp->SetRightMargin(0.4) ;

     tp->Draw() ;
     tp->cd() ;
     hfitqual_np->SetLabelSize(0.5,"x") ;
     TAxis *xaxis ;
     xaxis = hfitqual_np->GetXaxis() ;
     xaxis->SetBinLabel(1,"Eff PG") ;
     hfitqual_np->GetXaxis()->LabelsOption("v") ;
     hfitqual_np->Draw() ;

     cfitqual->Update() ;


     TGaxis* axis = new TGaxis() ;
     axis->SetLabelOffset(0.1) ;
     axis->SetLabelSize(0.30) ;
     axis->SetTickSize(0.2) ;
     axis->DrawAxis( 1.0, -5., 1.0, 5., -5., 5., 510, "+LS") ;
     
     cfitqual->Update() ;

     cfitqual->SaveAs("fitqual.gif") ;


     // print out results

//   cout << "\n\n Results: \n" << endl ;

//   for ( int i = 0 ; i < nBinsMET ; i++ ) {
//     for ( int j = 0 ; j < nBinsHT ; j++ ) {
//       for ( int k = 0 ; k < nBinsBtag ; k++ ) {     

//         binIndex = 2 + ( nBinsHT*i + j + k ) + k*(nBinsMET*nBinsHT) ;

//         cout << "--------------------------------------------------------------------------" << endl ;
//         cout << "\t\t\t ttbar  \t QCD \t\t Znn " << endl ;

//         TString String0lep = "0-lep";
//         TString String1lep = "1-lep";
//         TString Stringldp  = "  ldp";

//         String0lep += sMbins[i]+sHbins[j]+sBbins[k] ;
//         String1lep += sMbins[i]+sHbins[j]+sBbins[k] ;
//         Stringldp  += sMbins[i]+sHbins[j]+sBbins[k] ;

//         cout << String0lep << "\t\t" << hfitqual_ttwj_0lep->GetBinContent(binIndex) 
//      	<< "\t\t" << hfitqual_qcd_0lep->GetBinContent(binIndex) 
//      	<< "\t\t" << hfitqual_znn_0lep->GetBinContent(binIndex) << endl ;

//         cout << String1lep << "\t\t" << hfitqual_ttwj_1lep->GetBinContent(binIndex) 
//      	<< "\t\t" << hfitqual_qcd_1lep->GetBinContent(binIndex) 
//      	<< "\t\t" << hfitqual_znn_1lep->GetBinContent(binIndex) << endl ;

//         cout << Stringldp << "\t\t" << hfitqual_ttwj_ldp->GetBinContent(binIndex) 
//      	<< "\t\t" << hfitqual_qcd_ldp->GetBinContent(binIndex) 
//      	<< "\t\t" << hfitqual_znn_ldp->GetBinContent(binIndex) << endl ;



//       }
//     }
//   }



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




