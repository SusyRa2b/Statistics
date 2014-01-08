#include "RooStats/ModelConfig.h"

using namespace RooFit;
using namespace RooStats;

float ComputeTestStat(TString wsfile, double mu_susy_sig_val) {

  gROOT->Reset();

  TFile* wstf = new TFile( wsfile ) ;

  RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
  ws->Print() ;
  
  ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;
  
  modelConfig->Print() ;

  RooDataSet* rds = (RooDataSet*) ws->obj( "ra2b_observed_rds" ) ;
  
  rds->Print() ;
  rds->printMultiline(cout, 1, kTRUE, "") ;
  
  RooAbsPdf* likelihood = modelConfig->GetPdf() ;
  
  RooRealVar* rrv_mu_susy_sig = ws->var("mu_susy_all0lep") ;
  if ( rrv_mu_susy_sig == 0x0 ) {
    printf("\n\n\n *** can't find mu_susy_all0lep in workspace.  Quitting.\n\n\n") ;
    return ;
  } else {
    printf(" current value is : %8.3f\n", rrv_mu_susy_sig->getVal() ) ; cout << flush ;
    rrv_mu_susy_sig->setConstant(kFALSE) ;
  }

  /*
  // check the impact of varying the qcd normalization:

  RooRealVar *rrv_qcd_0lepLDP_ratioH1 = ws->var("qcd_0lepLDP_ratio_H1");
  RooRealVar *rrv_qcd_0lepLDP_ratioH2 = ws->var("qcd_0lepLDP_ratio_H2");
  RooRealVar *rrv_qcd_0lepLDP_ratioH3 = ws->var("qcd_0lepLDP_ratio_H3");
  
  rrv_qcd_0lepLDP_ratioH1->setVal(0.3);
  rrv_qcd_0lepLDP_ratioH2->setVal(0.3);
  rrv_qcd_0lepLDP_ratioH3->setVal(0.3);
  
  rrv_qcd_0lepLDP_ratioH1->setConstant(kTRUE);
  rrv_qcd_0lepLDP_ratioH2->setConstant(kTRUE);
  rrv_qcd_0lepLDP_ratioH3->setConstant(kTRUE);
  */
  
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


  return testStatVal ;

}
