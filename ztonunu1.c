
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooArgSet.h"
#include "RooConstVar.h"

#include "RooWorkspace.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

using namespace RooFit ;
using namespace RooStats ;

    void ztonunu1() {

    // RooMsgService::instance().addStream(DEBUG,Topic(Tracing),ClassName("RooPoisson"),OutputFile("debug.log")) ;

       float acc_mm_mean( 0.98 ) ; float acc_mm_err(  0.02 ) ;
       float acc_ee_mean( 0.98 ) ; float acc_ee_err(  0.02 ) ;
       float eff_mm_mean( 0.77 ) ; float eff_mm_err(  0.08 ) ;
       float eff_ee_mean( 0.76 ) ; float eff_ee_err(  0.08 ) ;


       RooRealVar* rv_Nsbee  = new RooRealVar( "Nsbee"  ,"Nsbee"  , 0., 100. ) ;
       RooRealVar* rv_Nsbmm  = new RooRealVar( "Nsbmm"  ,"Nsbmm"  , 0., 100. ) ;
       RooRealVar* rv_Nsigee = new RooRealVar( "Nsigee" ,"Nsigee" , 0., 100. ) ;
       RooRealVar* rv_Nsigmm = new RooRealVar( "Nsigmm" ,"Nsigmm" , 0., 100. ) ;

       rv_Nsbee->setVal( 6 ) ;
       rv_Nsbmm->setVal( 0 ) ;

       rv_Nsigee->setVal( 3 ) ;
       rv_Nsigmm->setVal( 0 ) ;


       RooRealVar* rv_mu_Znnsb  = new RooRealVar( "mu_Znnsb"  , "mu_Znnsb"  , 0., 100. ) ;
       RooRealVar* rv_mu_Znnsig = new RooRealVar( "mu_Znnsig" , "mu_Znnsig" , 0., 100. ) ;


       rv_mu_Znnsb->setVal( 37 ) ; // starting value
       rv_mu_Znnsig->setVal( 17. ) ; // starting value

       RooRealVar* rv_bfRatio = new RooRealVar( "bfRatio", "bfRatio", 0., 10. ) ;

       rv_bfRatio->setVal( 5.95 ) ;
       rv_bfRatio->setConstant( kTRUE ) ;

       RooRealVar* rv_acc_mm = new RooRealVar( "acc_mm", "acc_mm", 0.001, 1.000 ) ;
       RooRealVar* rv_acc_ee = new RooRealVar( "acc_ee", "acc_ee", 0.001, 1.000 ) ;

       RooRealVar* rv_eff_mm = new RooRealVar( "eff_mm", "eff_mm", 0.001, 1.000 ) ;
       RooRealVar* rv_eff_ee = new RooRealVar( "eff_ee", "eff_ee", 0.001, 1.000 ) ;

       rv_acc_mm->setVal( acc_mm_mean ) ;
       rv_acc_ee->setVal( acc_ee_mean ) ;
       rv_eff_mm->setVal( eff_mm_mean ) ;
       rv_eff_ee->setVal( eff_ee_mean ) ;



       RooFormulaVar* rv_mu_Zeesb = new RooFormulaVar( "mu_Zeesb",
                       "mu_Znnsb * ( acc_ee * eff_ee / bfRatio )",
                       RooArgSet( *rv_mu_Znnsb, *rv_acc_ee, *rv_eff_ee, *rv_bfRatio ) ) ;

       RooFormulaVar* rv_mu_Zmmsb = new RooFormulaVar( "mu_Zmmsb",
                       "mu_Znnsb * ( acc_mm * eff_mm / bfRatio )",
                       RooArgSet( *rv_mu_Znnsb, *rv_acc_mm, *rv_eff_mm, *rv_bfRatio ) ) ;

       RooFormulaVar* rv_mu_Zeesig = new RooFormulaVar( "mu_Zeesig",
                       "mu_Znnsig * ( acc_ee * eff_ee / bfRatio )",
                       RooArgSet( *rv_mu_Znnsig, *rv_acc_ee, *rv_eff_ee, *rv_bfRatio ) ) ;

       RooFormulaVar* rv_mu_Zmmsig = new RooFormulaVar( "mu_Zmmsig",
                       "mu_Znnsig * ( acc_mm * eff_mm / bfRatio )",
                       RooArgSet( *rv_mu_Znnsig, *rv_acc_mm, *rv_eff_mm, *rv_bfRatio ) ) ;



       RooFormulaVar* rv_n_sbee  = new RooFormulaVar( "n_sbee"  , "mu_Zeesb"  , RooArgSet( *rv_mu_Zeesb  ) ) ;
       RooFormulaVar* rv_n_sbmm  = new RooFormulaVar( "n_sbmm"  , "mu_Zmmsb"  , RooArgSet( *rv_mu_Zmmsb  ) ) ;
       RooFormulaVar* rv_n_sigee = new RooFormulaVar( "n_sigee" , "mu_Zeesig" , RooArgSet( *rv_mu_Zeesig ) ) ;
       RooFormulaVar* rv_n_sigmm = new RooFormulaVar( "n_sigmm" , "mu_Zmmsig" , RooArgSet( *rv_mu_Zmmsig ) ) ;



       RooGaussian* pdf_acc_mm = new RooGaussian( "pdf_acc_mm", "Gaussian pdf for Z to mumu acceptance",
                       *rv_acc_mm, RooConst( acc_mm_mean ), RooConst( acc_mm_err ) ) ;


       RooGaussian* pdf_acc_ee = new RooGaussian( "pdf_acc_ee", "Gaussian pdf for Z to ee acceptance",
                       *rv_acc_ee, RooConst( acc_ee_mean ), RooConst( acc_ee_err ) ) ;




       RooGaussian* pdf_eff_mm = new RooGaussian( "pdf_eff_mm", "Gaussian pdf for Z to mumu efficiency",
                       *rv_eff_mm, RooConst( eff_mm_mean ), RooConst( eff_mm_err ) ) ;



       RooGaussian* pdf_eff_ee = new RooGaussian( "pdf_eff_ee", "Gaussian pdf for Z to ee efficiency",
                       *rv_eff_ee, RooConst( eff_ee_mean ), RooConst( eff_ee_err ) ) ;



       RooPoisson* pdf_Nsbee   = new RooPoisson( "pdf_Nsbee"  , "Nsb , Z to ee Poisson PDF", *rv_Nsbee  , *rv_n_sbee  ) ;
       RooPoisson* pdf_Nsbmm   = new RooPoisson( "pdf_Nsbmm"  , "Nsb , Z to mm Poisson PDF", *rv_Nsbmm  , *rv_n_sbmm  ) ;
       RooPoisson* pdf_Nsigee  = new RooPoisson( "pdf_Nsigee" , "Nsig, Z to ee Poisson PDF", *rv_Nsigee , *rv_n_sigee ) ;
       RooPoisson* pdf_Nsigmm  = new RooPoisson( "pdf_Nsigmm" , "Nsig, Z to mm Poisson PDF", *rv_Nsigmm , *rv_n_sigmm ) ;

       RooArgSet pdflist ;
       pdflist.add( *pdf_acc_mm ) ;
       pdflist.add( *pdf_acc_ee ) ;
       pdflist.add( *pdf_eff_mm ) ;
       pdflist.add( *pdf_eff_ee ) ;
       pdflist.add( *pdf_Nsbee  ) ;
       pdflist.add( *pdf_Nsbmm  ) ;
       pdflist.add( *pdf_Nsigee  ) ;
       pdflist.add( *pdf_Nsigmm  ) ;

       RooProdPdf* znnLikelihood = new RooProdPdf( "znnLikelihood", "Z to nunu likelihood", pdflist ) ;


       RooArgSet observedParametersList ;
       observedParametersList.add( *rv_Nsbee ) ;
       observedParametersList.add( *rv_Nsbmm ) ;
       observedParametersList.add( *rv_Nsigee ) ;
       observedParametersList.add( *rv_Nsigmm ) ;

       RooDataSet* dsObserved = new RooDataSet( "ztonn_rds", "Z to nunu dataset", observedParametersList ) ;
       dsObserved->add( observedParametersList ) ;
       //// RooDataSet* dsObserved = znnLikelihood->generate( observedParametersList, 1) ;

       RooWorkspace* znnWorkspace = new RooWorkspace("ztonn_ws") ;
       znnWorkspace->import( *znnLikelihood ) ;
       znnWorkspace->import( *dsObserved ) ;

       znnWorkspace->Print() ;

       dsObserved->printMultiline(cout, 1, kTRUE, "") ;
       RooFitResult* fitResult = znnLikelihood->fitTo( *dsObserved, Verbose(true), Save(true) ) ;

       printf("\n\n----- Constant parameters:\n") ;
       RooArgList constPars = fitResult->constPars() ;
       for ( int pi=0; pi<constPars.getSize(); pi++ ) {
          constPars[pi].Print() ;
       } // pi.

       printf("\n\n----- Floating parameters:\n") ;
       RooArgList floatPars = fitResult->floatParsFinal() ;
       for ( int pi=0; pi<floatPars.getSize(); pi++ ) {
          floatPars[pi].Print() ;
       } // pi.
       printf("\n\n") ;



       ProfileLikelihoodCalculator plc_znn_sb(  *dsObserved, *znnLikelihood, RooArgSet( *rv_mu_Znnsb ) ) ;
       ProfileLikelihoodCalculator plc_znn_sig( *dsObserved, *znnLikelihood, RooArgSet( *rv_mu_Znnsig ) ) ;

       plc_znn_sb.SetTestSize(0.32) ;
       plc_znn_sig.SetTestSize(0.32) ;

       ConfInterval* sb_znn_interval = plc_znn_sb.GetInterval() ;
       float sbZnnLow  = ((LikelihoodInterval*) sb_znn_interval)->LowerLimit(*rv_mu_Znnsb) ;
       float sbZnnHigh = ((LikelihoodInterval*) sb_znn_interval)->UpperLimit(*rv_mu_Znnsb) ;
       printf("\n\n znn SB interval %6.1f to %6.1f\n", sbZnnLow, sbZnnHigh ) ;

       ConfInterval* sig_znn_interval = plc_znn_sig.GetInterval() ;
       float sigZnnLow  = ((LikelihoodInterval*) sig_znn_interval)->LowerLimit(*rv_mu_Znnsig) ;
       float sigZnnHigh = ((LikelihoodInterval*) sig_znn_interval)->UpperLimit(*rv_mu_Znnsig) ;
       printf("\n\n znn SIG interval %6.1f to %6.1f\n", sigZnnLow, sigZnnHigh ) ;



       TCanvas* c_prof_sb = new TCanvas("c_prof_sb","SB Z to nunu profile") ;
       LikelihoodIntervalPlot plot_znn_sb((LikelihoodInterval*)sb_znn_interval) ;
       plot_znn_sb.Draw() ;
       c_prof_sb->SaveAs("znn_sb_profile.png") ;

       TCanvas* c_prof_sig = new TCanvas("c_prof_sig","sig Z to nunu profile") ;
       LikelihoodIntervalPlot plot_znn_sig((LikelihoodInterval*)sig_znn_interval) ;
       plot_znn_sig.Draw() ;
       c_prof_sig->SaveAs("znn_sig_profile.png") ;


    }




