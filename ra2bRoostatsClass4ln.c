
//
//   Owen Long, UCR
//   Harrison Prosper, FSU
//   Sezen Sekmen, FSU
//
//


#include "ra2bRoostatsClass4ln.h"

#include <iostream>
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

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"

#include "LikelihoodIntervalPlot.cxx"

  using namespace RooFit ;
  using namespace RooStats ;


  //=====================================================================================================


   ra2bRoostatsClass4ln::ra2bRoostatsClass4ln( bool ArgUseSigTtwjVar, bool ArgUseLdpVars, int ArgznnModel ) {

      gStyle->SetOptStat(0) ;

      useSigTtwjVar = ArgUseSigTtwjVar ;
      useLdpVars = ArgUseLdpVars ;
      znnModel = ArgznnModel ;

      if ( znnModel > 2 ) {
         printf("\n\n\n *** Unrecognized znn model number %d.  Should be 1 or 2.\n\n\n", znnModel ) ;
      }

     //--- Tell RooFit to shut up about anything less important than an ERROR.
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
  /// RooMsgService::instance().addStream(DEBUG,Topic(Tracing),OutputFile("debug.log")) ;
      printf("\n\n ==== RooFit output configuration ===============\n") ;
      RooMsgService::instance().Print("v") ;
      printf("\n\n ================================================\n") ;

      varsAtFitVals = false ;
      initialized = false ;




      //--- CLs stuff below here
      tt_bgonly = 0x0 ;
      tt_splusb = 0x0 ;

      toymean_n_sig         = 0. ;
      toymean_n_sb          = 0. ;
      toymean_n_sig_ldp     = 0. ;
      toymean_n_sb_ldp      = 0. ;
      toymean_n_sig_sl      = 0. ;
      toymean_n_sb_sl       = 0. ;
      toymean_n_lsb_0b      = 0. ;
      toymean_n_lsb_0b_ldp  = 0. ;

       //-- Znn model 1
      toymean_n_sig_ee      = 0. ;
      toymean_n_sb_ee       = 0. ;
      toymean_n_sig_mm      = 0. ;
      toymean_n_sb_mm       = 0. ;

       //-- Znn model 2
      toymean_n_sigsb_ee    = 0. ;
      toymean_n_sigsb_mm    = 0. ;





      trandom_cls = new TRandom2(12345) ;

   }




////==================================================================================================

//  bool ra2bRoostatsClass4ln::generateToyDatasetsFromLikelihood( const char* outputRootFile, int nToys ) {

//     if ( ! initialized ) {
//        printf("\n\n *** Call initialize first.\n\n") ;
//        return false ;
//     }

//     TTree* toyOutputTrueTree = new TTree("toytruetree", "ra2b toy mean true value tree" ) ;

//     toyOutputTrueTree->Branch("mu0_ttbar_sig"      , &toy_mu0_ttbar_sig, "mu0_ttbar_sig/D" ) ;
//     toyOutputTrueTree->Branch("mu0_qcd_sig"        , &toy_mu0_qcd_sig, "mu0_qcd_sig/D" ) ;
//     toyOutputTrueTree->Branch("mu0_ttbar_sb"       , &toy_mu0_ttbar_sb, "mu0_ttbar_sb/D" ) ;
//     toyOutputTrueTree->Branch("mu0_qcd_sb"         , &toy_mu0_qcd_sb, "mu0_qcd_sb/D" ) ;
//     toyOutputTrueTree->Branch("mu0_susy_sig"       , &toy_mu0_susy_sig       , "mu0_susy_sig/D" ) ;
//     toyOutputTrueTree->Branch("mu0_allbg_sig"      , &toy_mu0_allbg_sig      , "mu0_allbg_sig/D" ) ;


//    //--- Set the true values.
//     if ( useSigTtwjVar ) {
//        toy_mu0_ttbar_sig = rrv_mu_ttbar_sig->getVal() ;
//        toy_mu0_qcd_sig   = rrv_mu_qcd_sig->getVal() ;
//        toy_mu0_ttbar_sb  = rfv_mu_ttbar_sb->getVal() ;
//        toy_mu0_qcd_sb    = rfv_mu_qcd_sb->getVal() ;
//     } else {
//        toy_mu0_ttbar_sig = rfv_mu_ttbar_sig->getVal() ;
//        toy_mu0_qcd_sig   = rfv_mu_qcd_sig->getVal() ;
//        toy_mu0_ttbar_sb  = rrv_mu_ttbar_sb->getVal() ;
//        toy_mu0_qcd_sb    = rrv_mu_qcd_sb->getVal() ;
//     }
//     toy_mu0_susy_sig = rv_mu_susy_sig->getVal() ;
//     toy_mu0_allbg_sig = rv_mu_ew_sig->getVal() + toy_mu0_ttbar_sig + toy_mu0_qcd_sig ;
//     printf("\n\n") ;
//     printf("  Mean true values used in toy generation:\n") ;
//     printf("  mu0_ttbar_sig : %6.1f\n", toy_mu0_ttbar_sig ) ;
//     printf("  mu0_qcd_sig   : %6.1f\n", toy_mu0_qcd_sig ) ;
//     printf("  mu0_ttbar_sb  : %6.1f\n", toy_mu0_ttbar_sb ) ;
//     printf("  mu0_qcd_sb    : %6.1f\n", toy_mu0_qcd_sb ) ;
//     printf("  mu0_susy_sig  : %6.1f\n", toy_mu0_susy_sig ) ;
//     printf("  mu0_allbg_sig : %6.1f\n", toy_mu0_allbg_sig ) ;

//   //--- generate some toy datasets.
//     printf("\n\n Opening output toy datasets root file : %s\n", outputRootFile ) ;
//     TFile toyOutputFile( outputRootFile, "recreate" ) ;
//     RooDataSet* toyDatasets = likelihood->generate( observedParametersList, nToys ) ;
//     const TTree* toyTree = toyDatasets->tree() ;
//     toyTree->Write() ;
//     toyOutputTrueTree->Fill() ;
//     toyOutputTrueTree->Write() ;
//     printf("\n\n Closing output toy datasets root file : %s\n", outputRootFile ) ;
//     toyOutputFile.Close() ;

//     return true ;

//  } // generateToyDatasetsFromLikelihood

////==================================================================================================

//  bool ra2bRoostatsClass4ln::generateToyDatasetsFromNVals( const char* outputRootFile, int nToys ) {

//     if ( ! initialized ) {
//        printf("\n\n *** Call initialize first.\n\n") ;
//        return false ;
//     }

//     TTree* toyOutputTrueTree = new TTree("toytruetree", "ra2b toy mean true value tree" ) ;

//     toyOutputTrueTree->Branch("mu0_ttbar_sig"      , &toy_mu0_ttbar_sig, "mu0_ttbar_sig/D" ) ;
//     toyOutputTrueTree->Branch("mu0_qcd_sig"        , &toy_mu0_qcd_sig, "mu0_qcd_sig/D" ) ;
//     toyOutputTrueTree->Branch("mu0_ttbar_sb"       , &toy_mu0_ttbar_sb, "mu0_ttbar_sb/D" ) ;
//     toyOutputTrueTree->Branch("mu0_qcd_sb"         , &toy_mu0_qcd_sb, "mu0_qcd_sb/D" ) ;
//     toyOutputTrueTree->Branch("mu0_susy_sig"       , &toy_mu0_susy_sig       , "mu0_susy_sig/D" ) ;
//     toyOutputTrueTree->Branch("mu0_allbg_sig"      , &toy_mu0_allbg_sig      , "mu0_allbg_sig/D" ) ;


//    //--- Set the true values.
//     if ( useSigTtwjVar ) {
//        toy_mu0_ttbar_sig = rrv_mu_ttbar_sig->getVal() ;
//        toy_mu0_qcd_sig   = rrv_mu_qcd_sig->getVal() ;
//        toy_mu0_ttbar_sb  = rfv_mu_ttbar_sb->getVal() ;
//        toy_mu0_qcd_sb    = rfv_mu_qcd_sb->getVal() ;
//     } else {
//        toy_mu0_ttbar_sig = rfv_mu_ttbar_sig->getVal() ;
//        toy_mu0_qcd_sig   = rfv_mu_qcd_sig->getVal() ;
//        toy_mu0_ttbar_sb  = rrv_mu_ttbar_sb->getVal() ;
//        toy_mu0_qcd_sb    = rrv_mu_qcd_sb->getVal() ;
//     }
//     toy_mu0_susy_sig = rv_mu_susy_sig->getVal() ;
//     toy_mu0_allbg_sig = rv_mu_ew_sig->getVal() + toy_mu0_ttbar_sig + toy_mu0_qcd_sig ;
//     printf("\n\n") ;
//     printf("  Mean true values used in toy generation:\n") ;
//     printf("  mu0_ttbar_sig : %6.1f\n", toy_mu0_ttbar_sig ) ;
//     printf("  mu0_qcd_sig   : %6.1f\n", toy_mu0_qcd_sig ) ;
//     printf("  mu0_ttbar_sb  : %6.1f\n", toy_mu0_ttbar_sb ) ;
//     printf("  mu0_qcd_sb    : %6.1f\n", toy_mu0_qcd_sb ) ;
//     printf("  mu0_susy_sig  : %6.1f\n", toy_mu0_susy_sig ) ;
//     printf("  mu0_allbg_sig : %6.1f\n", toy_mu0_allbg_sig ) ;

//   //--- generate some toy datasets.
//     printf("\n\n Opening output toy datasets root file : %s\n", outputRootFile ) ;
//     TTree* toyTree = new TTree("likelihoodData", "ra2b observed values generated from Nvals") ;

//     double Nsig, Na, Nd,
//            Nsb1, Nsb2, Nsb3, Nsb4, Nsb5,
//            Nlsb1, Nlsb2, Nlsb3, Nlsb4, Nlsb5,
//            Nslsig1, Nslsig2, Nslsig3, Nslsig4, Nslsig5,
//            Nslsb1, Nslsb2, Nslsb3, Nslsb4, Nslsb5,
//            Nslmsb1, Nslmsb2, Nslmsb3, Nslmsb4, Nslmsb5,
//            Nqcdmca, Nqcdmcd, Nqcdmcsig, Nqcdmcsb ;

//     toyTree->Branch("Nsig"       , &Nsig      , "Nsig/D" ) ;
//     toyTree->Branch("Na"         , &Na        , "Na/D" ) ;
//     toyTree->Branch("Nd"         , &Nd        , "Nd/D" ) ;
//     toyTree->Branch("Nsb1"       , &Nsb1      , "Nsb1/D" ) ;
//     toyTree->Branch("Nsb2"       , &Nsb2      , "Nsb2/D" ) ;
//     toyTree->Branch("Nsb3"       , &Nsb3      , "Nsb3/D" ) ;
//     toyTree->Branch("Nsb4"       , &Nsb4      , "Nsb4/D" ) ;
//     toyTree->Branch("Nsb5"       , &Nsb5      , "Nsb5/D" ) ;
//     toyTree->Branch("Nlsb1"      , &Nlsb1     , "Nlsb1/D" ) ;
//     toyTree->Branch("Nlsb2"      , &Nlsb2     , "Nlsb2/D" ) ;
//     toyTree->Branch("Nlsb3"      , &Nlsb3     , "Nlsb3/D" ) ;
//     toyTree->Branch("Nlsb4"      , &Nlsb4     , "Nlsb4/D" ) ;
//     toyTree->Branch("Nlsb5"      , &Nlsb5     , "Nlsb5/D" ) ;
//     toyTree->Branch("Nslsig1"    , &Nslsig1   , "Nslsig1/D" ) ;
//     toyTree->Branch("Nslsig2"    , &Nslsig2   , "Nslsig2/D" ) ;
//     toyTree->Branch("Nslsig3"    , &Nslsig3   , "Nslsig3/D" ) ;
//     toyTree->Branch("Nslsig4"    , &Nslsig4   , "Nslsig4/D" ) ;
//     toyTree->Branch("Nslsig5"    , &Nslsig5   , "Nslsig5/D" ) ;
//     toyTree->Branch("Nslsb1"     , &Nslsb1    , "Nslsb1/D" ) ;
//     toyTree->Branch("Nslsb2"     , &Nslsb2    , "Nslsb2/D" ) ;
//     toyTree->Branch("Nslsb3"     , &Nslsb3    , "Nslsb3/D" ) ;
//     toyTree->Branch("Nslsb4"     , &Nslsb4    , "Nslsb4/D" ) ;
//     toyTree->Branch("Nslsb5"     , &Nslsb5    , "Nslsb5/D" ) ;
//     toyTree->Branch("Nslmsb1"    , &Nslmsb1   , "Nslmsb1/D" ) ;
//     toyTree->Branch("Nslmsb2"    , &Nslmsb2   , "Nslmsb2/D" ) ;
//     toyTree->Branch("Nslmsb3"    , &Nslmsb3   , "Nslmsb3/D" ) ;
//     toyTree->Branch("Nslmsb4"    , &Nslmsb4   , "Nslmsb4/D" ) ;
//     toyTree->Branch("Nslmsb5"    , &Nslmsb5   , "Nslmsb5/D" ) ;
//     toyTree->Branch("Nqcdmca"    , &Nqcdmca   , "Nqcdmca/D" ) ;
//     toyTree->Branch("Nqcdmcd"    , &Nqcdmcd   , "Nqcdmcd/D" ) ;
//     toyTree->Branch("Nqcdmcsig"  , &Nqcdmcsig , "Nqcdmcsig/D" ) ;
//     toyTree->Branch("Nqcdmcsb"   , &Nqcdmcsb  , "Nqcdmcsb/D" ) ;

//     TRandom1 rg(12345) ;

//     for ( int dsi=0; dsi<nToys; dsi++ ) {

//        Nsig      = rg.Poisson( rv_Nsig      ->getVal() ) ;
//        Na        = rg.Poisson( rv_Na        ->getVal() ) ;
//        Nd        = rg.Poisson( rv_Nd        ->getVal() ) ;
//        Nsb1      = rg.Poisson( rv_Nsb1      ->getVal() ) ;
//        Nsb2      = rg.Poisson( rv_Nsb2      ->getVal() ) ;
//        Nsb3      = rg.Poisson( rv_Nsb3      ->getVal() ) ;
//        Nsb4      = rg.Poisson( rv_Nsb4      ->getVal() ) ;
//        Nsb5      = rg.Poisson( rv_Nsb5      ->getVal() ) ;
//        Nlsb1     = rg.Poisson( rv_Nlsb1     ->getVal() ) ;
//        Nlsb2     = rg.Poisson( rv_Nlsb2     ->getVal() ) ;
//        Nlsb3     = rg.Poisson( rv_Nlsb3     ->getVal() ) ;
//        Nlsb4     = rg.Poisson( rv_Nlsb4     ->getVal() ) ;
//        Nlsb5     = rg.Poisson( rv_Nlsb5     ->getVal() ) ;
//        Nslsig1   = rg.Poisson( rv_Nslsig1   ->getVal() ) ;
//        Nslsig2   = rg.Poisson( rv_Nslsig2   ->getVal() ) ;
//        Nslsig3   = rg.Poisson( rv_Nslsig3   ->getVal() ) ;
//        Nslsig4   = rg.Poisson( rv_Nslsig4   ->getVal() ) ;
//        Nslsig5   = rg.Poisson( rv_Nslsig5   ->getVal() ) ;
//        Nslsb1    = rg.Poisson( rv_Nslsb1    ->getVal() ) ;
//        Nslsb2    = rg.Poisson( rv_Nslsb2    ->getVal() ) ;
//        Nslsb3    = rg.Poisson( rv_Nslsb3    ->getVal() ) ;
//        Nslsb4    = rg.Poisson( rv_Nslsb4    ->getVal() ) ;
//        Nslsb5    = rg.Poisson( rv_Nslsb5    ->getVal() ) ;
//        Nslmsb1   = rg.Poisson( rv_Nslmsb1   ->getVal() ) ;
//        Nslmsb2   = rg.Poisson( rv_Nslmsb2   ->getVal() ) ;
//        Nslmsb3   = rg.Poisson( rv_Nslmsb3   ->getVal() ) ;
//        Nslmsb4   = rg.Poisson( rv_Nslmsb4   ->getVal() ) ;
//        Nslmsb5   = rg.Poisson( rv_Nslmsb5   ->getVal() ) ;
//        Nqcdmca   = rg.Poisson( rv_Nqcdmca   ->getVal() ) ;
//        Nqcdmcd   = rg.Poisson( rv_Nqcdmcd   ->getVal() ) ;
//        Nqcdmcsig = rg.Poisson( rv_Nqcdmcsig ->getVal() ) ;
//        Nqcdmcsb  = rg.Poisson( rv_Nqcdmcsb  ->getVal() ) ;

//        toyTree->Fill() ;

//     } // dsi.


//     TFile toyOutputFile( outputRootFile, "recreate" ) ;
//     toyTree->Write() ;
//     toyOutputTrueTree->Fill() ;
//     toyOutputTrueTree->Write() ;
//     printf("\n\n Closing output toy datasets root file : %s\n", outputRootFile ) ;
//     toyOutputFile.Close() ;

//     return true ;

//  } // generateToyDatasetsFromNVals

  //==================================================================================================

    bool ra2bRoostatsClass4ln::doFit( ) {

       if ( ! initialized ) {
          printf("\n\n *** Call initialize first.\n\n") ;
          return false ;
       }

       printf("\n\n") ;
       printf("  Fitting with these values for the observables.\n") ;
       dsObserved->printMultiline(cout, 1, kTRUE, "") ;
       printf("\n\n") ;

    // fitResult = likelihood->fitTo(*dsObserved, Save(true), Verbose(true) );
       fitResult = likelihood->fitTo(*dsObserved, Save(true) );

       RooArgList constPars = fitResult->constPars() ;
       printf("\n\n----- Constant parameters: %d\n", constPars.getSize() ) ;
       for ( int pi=0; pi<constPars.getSize(); pi++ ) {
          printf(" par %2d : ", pi ) ;
          constPars[pi].Print() ;
       } // pi.

       RooArgList floatPars = fitResult->floatParsFinal() ;
       printf("\n\n----- Floating parameters: %d\n", floatPars.getSize() ) ;
       for ( int pi=0; pi<floatPars.getSize(); pi++ ) {
          printf(" par %2d : ", pi ) ;
          floatPars[pi].Print() ;
       } // pi.
       printf("\n\n") ;
       printf(" here 1\n") ;
       if ( znnModel == 2 ) {
          printf("\n Znn SB value : %6.2f\n\n", rfv_mu_znn_sb->getVal() ) ;
       }
       printf(" here 2\n") ;


       varsAtFitVals = true ;

       return true ;

     } // doFit .


  //==================================================================================================

     bool ra2bRoostatsClass4ln::profileSusySig( float& susySigLow, float& susySigHigh, bool makePlot, const char* plotname, double scanMax ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }

      //--- Profile likelihood for signal susy yield.

         if ( scanMax > 0 ) { rv_mu_susy_sig->setMax(scanMax) ; }

         printf("\n\n Creating ProfileLikelihoodCalculator for susy sig.\n\n") ;
      //---------------
         ProfileLikelihoodCalculator plc_susy_sig( *dsObserved, *likelihood, RooArgSet( *rv_mu_susy_sig ) ) ;
///// //---------------
/////    RooArgSet* rasNull = new RooArgSet() ;
/////    rasNull->add( *rv_mu_susy_sig ) ;
/////    ProfileLikelihoodCalculator plc_susy_sig( *dsObserved, *likelihood, RooArgSet( *rv_mu_susy_sig ),
/////         0.05, rasNull  ) ;
///// //---------------

/////    HypoTestResult* htr_susy_sig =  plc_susy_sig.GetHypoTest() ;

/////    if ( htr_susy_sig != 0x0 ) {

/////       htr_susy_sig->Print() ;

/////    } else {

/////       printf("\n\n *** GetHypoTest returned null pointer.\n\n") ;

/////    }


         plc_susy_sig.SetTestSize(0.05) ;
         ConfInterval* plinterval_susy_sig = plc_susy_sig.GetInterval() ;
         susySigLow  = ((LikelihoodInterval*) plinterval_susy_sig)->LowerLimit(*rv_mu_susy_sig) ;
         susySigHigh = ((LikelihoodInterval*) plinterval_susy_sig)->UpperLimit(*rv_mu_susy_sig) ;
         printf("\n\n") ;
         printf("    susy, SIG 95%% CL interval  [%5.1f, %5.1f]\n\n", susySigLow, susySigHigh ) ;


         if ( makePlot ) {
            TCanvas* plcplot_susy_sig = new TCanvas("plcplot_susy_sig", "susy sig, Profile likelihood", 500, 400 ) ;
            LikelihoodIntervalPlot plotInt_susy_sig((LikelihoodInterval*)plinterval_susy_sig);
            plotInt_susy_sig.Draw() ;
            gPad->SetGridy(1) ;
            plcplot_susy_sig->SaveAs( plotname ) ;
         }

         varsAtFitVals = false ;

         delete plinterval_susy_sig ; // can I safely do this???

         return true ;

     }

  //==================================================================================================

     bool ra2bRoostatsClass4ln::profileZnnSig( float& znnSigLow, float& znnSigHigh, bool makePlot, const char* plotname, double scanMax ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }

      //--- Profile likelihood for signal Z to nunu yield.

         double bestVal = rv_mu_znn_sig->getVal() ;

         if ( scanMax > 0 ) { rv_mu_znn_sig->setMax(scanMax) ; }

         ProfileLikelihoodCalculator plc_znn_sig( *dsObserved, *likelihood, RooArgSet( *rv_mu_znn_sig ) ) ;
         plc_znn_sig.SetTestSize(0.3173) ;
         ConfInterval* plinterval_znn_sig = plc_znn_sig.GetInterval() ;
         znnSigLow  = ((LikelihoodInterval*) plinterval_znn_sig)->LowerLimit(*rv_mu_znn_sig) ;
         znnSigHigh = ((LikelihoodInterval*) plinterval_znn_sig)->UpperLimit(*rv_mu_znn_sig) ;
         printf("\n\n") ;
         printf("    znn, SIG 68.3%% CL interval  [%5.1f, %5.1f]\n\n", znnSigLow, znnSigHigh ) ;
         printf("    znn SIG = %5.1f  +%5.1f  -%5.1f\n\n", bestVal, znnSigHigh-bestVal, bestVal-znnSigLow ) ;

         if ( makePlot ) {
            TCanvas* plcplot_znn_sig = new TCanvas("plcplot_znn_sig", "znn sig, Profile likelihood", 500, 400 ) ;
            LikelihoodIntervalPlot plotInt_znn_sig((LikelihoodInterval*)plinterval_znn_sig);
            plotInt_znn_sig.Draw() ;
            gPad->SetGridy(1) ;
            plcplot_znn_sig->SaveAs( plotname ) ;
         }

         varsAtFitVals = false ;

         delete plinterval_znn_sig ; // can I safely do this???

         return true ;

     }

  //==================================================================================================

     bool ra2bRoostatsClass4ln::profileqcdSig( float& qcdSigLow, float& qcdSigHigh, bool makePlot, const char* plotname, double scanMax ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }

         if ( useLdpVars ) {
            printf("\n\n *** Try again with useLdpVars set to false in the constructor.\n\n") ;
            return false ;
         }

         double bestVal = rrv_mu_qcd_sig->getVal() ;

      //--- Profile likelihood for signal qcd yield.

         if ( scanMax > 0 ) { rrv_mu_qcd_sig->setMax(scanMax) ; }

         ProfileLikelihoodCalculator plc_qcd_sig( *dsObserved, *likelihood, RooArgSet( *rrv_mu_qcd_sig ) ) ;
         plc_qcd_sig.SetTestSize(0.3173) ;
         ConfInterval* plinterval_qcd_sig = plc_qcd_sig.GetInterval() ;
         qcdSigLow  = ((LikelihoodInterval*) plinterval_qcd_sig)->LowerLimit(*rrv_mu_qcd_sig) ;
         qcdSigHigh = ((LikelihoodInterval*) plinterval_qcd_sig)->UpperLimit(*rrv_mu_qcd_sig) ;
         printf("\n\n") ;
         printf("    qcd, SIG 68.3%% CL interval  [%5.1f, %5.1f]\n\n", qcdSigLow, qcdSigHigh ) ;
         printf("    qcd SIG = %5.1f  +%5.1f  -%5.1f\n\n", bestVal, qcdSigHigh-bestVal, bestVal-qcdSigLow ) ;

         if ( makePlot ) {
            TCanvas* plcplot_qcd_sig = new TCanvas("plcplot_qcd_sig", "qcd sig, Profile likelihood", 500, 400 ) ;
            LikelihoodIntervalPlot plotInt_qcd_sig((LikelihoodInterval*)plinterval_qcd_sig);
            plotInt_qcd_sig.Draw() ;
            gPad->SetGridy(1) ;
            plcplot_qcd_sig->SaveAs( plotname ) ;
         }

         varsAtFitVals = false ;

         delete plinterval_qcd_sig ; // can I safely do this???

         return true ;

     }

  //==================================================================================================
     bool ra2bRoostatsClass4ln::profileqcdSb( float& qcdSbLow, float& qcdSbHigh, bool makePlot, const char* plotname, double scanMax ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }

         if ( useLdpVars ) {
            printf("\n\n *** Try again with useLdpVars set to false in the constructor.\n\n") ;
            return false ;
         }

         double bestVal = rrv_mu_qcd_sb->getVal() ;

      //--- Profile likelihood for SB QCD yield.

         if ( scanMax > 0 ) { rrv_mu_qcd_sb->setMax(scanMax) ; }

         ProfileLikelihoodCalculator plc_qcd_sb( *dsObserved, *likelihood, RooArgSet( *rrv_mu_qcd_sb ) ) ;
         plc_qcd_sb.SetTestSize(0.3173) ;
         ConfInterval* plinterval_qcd_sb = plc_qcd_sb.GetInterval() ;
         qcdSbLow  = ((LikelihoodInterval*) plinterval_qcd_sb)->LowerLimit(*rrv_mu_qcd_sb) ;
         qcdSbHigh = ((LikelihoodInterval*) plinterval_qcd_sb)->UpperLimit(*rrv_mu_qcd_sb) ;
         printf("\n\n") ;
         printf("    qcd, sb 68.3%% CL interval  [%5.1f, %5.1f]\n\n", qcdSbLow, qcdSbHigh ) ;
         printf("    qcd SB = %5.1f  +%5.1f  -%5.1f\n\n", bestVal, qcdSbHigh-bestVal, bestVal-qcdSbLow ) ;

         if ( makePlot ) {
            TCanvas* plcplot_qcd_sb = new TCanvas("plcplot_qcd_sb", "qcd sb, Profile likelihood", 500, 400 ) ;
            LikelihoodIntervalPlot plotInt_qcd_sb((LikelihoodInterval*)plinterval_qcd_sb);
            plotInt_qcd_sb.Draw() ;
            gPad->SetGridy(1) ;
            plcplot_qcd_sb->SaveAs( plotname ) ;
         }

         varsAtFitVals = false ;

         delete plinterval_qcd_sb ; // can I safely do this???

         return true ;

     }

  //==================================================================================================

     bool ra2bRoostatsClass4ln::profilettwjSig( float& ttwjSigLow, float& ttwjSigHigh, bool makePlot, const char* plotname, double scanMax ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }

         if ( !useSigTtwjVar ) {
            printf("\n\n *** Try again with useSigTtwjVar set to true in the constructor.\n\n") ;
            return false ;
         }

         double bestVal = rrv_mu_ttwj_sig->getVal() ;

      //--- Profile likelihood for signal ttwj yield.

         if ( scanMax > 0 ) { rrv_mu_ttwj_sig->setMax(scanMax) ; }

         ProfileLikelihoodCalculator plc_ttwj_sig( *dsObserved, *likelihood, RooArgSet( *rrv_mu_ttwj_sig ) ) ;
         plc_ttwj_sig.SetTestSize(0.3173) ;
         ConfInterval* plinterval_ttwj_sig = plc_ttwj_sig.GetInterval() ;
         ttwjSigLow  = ((LikelihoodInterval*) plinterval_ttwj_sig)->LowerLimit(*rrv_mu_ttwj_sig) ;
         ttwjSigHigh = ((LikelihoodInterval*) plinterval_ttwj_sig)->UpperLimit(*rrv_mu_ttwj_sig) ;
         printf("\n\n") ;
         printf("    ttwj, SIG 68.3%% CL interval  [%5.1f, %5.1f]\n\n", ttwjSigLow, ttwjSigHigh ) ;
         printf("    ttwj SIG = %5.1f  +%5.1f  -%5.1f\n\n", bestVal, ttwjSigHigh-bestVal, bestVal-ttwjSigLow ) ;

         if ( makePlot ) {
            TCanvas* plcplot_ttwj_sig = new TCanvas("plcplot_ttwj_sig", "ttwj sig, Profile likelihood", 500, 400 ) ;
            LikelihoodIntervalPlot plotInt_ttwj_sig((LikelihoodInterval*)plinterval_ttwj_sig);
            plotInt_ttwj_sig.Draw() ;
            gPad->SetGridy(1) ;
            plcplot_ttwj_sig->SaveAs( plotname) ;
         }

         varsAtFitVals = false ;

         delete plinterval_ttwj_sig ; // can I safely do this???

         return true ;

     }

  //==================================================================================================
     bool ra2bRoostatsClass4ln::profilettwjSb( float& ttwjSbLow, float& ttwjSbHigh, bool makePlot, const char* plotname, double scanMax ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }

         if ( useSigTtwjVar ) {
            printf("\n\n *** Try again with useSigTtwjVar set to false in the constructor.\n\n") ;
            return false ;
         }

         double bestVal = rrv_mu_ttwj_sb->getVal() ;

      //--- Profile likelihood for ttwj SB yield.

         if ( scanMax > 0 ) { rrv_mu_ttwj_sb->setMax(scanMax) ; }

         ProfileLikelihoodCalculator plc_ttwj_sb( *dsObserved, *likelihood, RooArgSet( *rrv_mu_ttwj_sb ) ) ;
         plc_ttwj_sb.SetTestSize(0.3173) ;
         ConfInterval* plinterval_ttwj_sb = plc_ttwj_sb.GetInterval() ;
         ttwjSbLow  = ((LikelihoodInterval*) plinterval_ttwj_sb)->LowerLimit(*rrv_mu_ttwj_sb) ;
         ttwjSbHigh = ((LikelihoodInterval*) plinterval_ttwj_sb)->UpperLimit(*rrv_mu_ttwj_sb) ;
         printf("\n\n") ;
         printf("    ttwj, sb 68.3%% CL interval  [%5.1f, %5.1f]\n\n", ttwjSbLow, ttwjSbHigh ) ;
         printf("    ttwj SB = %5.1f  +%5.1f  -%5.1f\n\n", bestVal, ttwjSbHigh-bestVal, bestVal-ttwjSbLow ) ;

         if ( makePlot ) {
            TCanvas* plcplot_ttwj_sb = new TCanvas("plcplot_ttwj_sb", "ttwj sb, Profile likelihood", 500, 400 ) ;
            LikelihoodIntervalPlot plotInt_ttwj_sb((LikelihoodInterval*)plinterval_ttwj_sb);
            plotInt_ttwj_sb.Draw() ;
            gPad->SetGridy(1) ;
            plcplot_ttwj_sb->SaveAs( plotname ) ;
         }

         varsAtFitVals = false ;

         delete plinterval_ttwj_sb ; // can I safely do this???

         return true ;

     }

  //==================================================================================================

     bool ra2bRoostatsClass4ln::profileZnnSb( float& znnSbLow, float& znnSbHigh, bool makePlot, const char* plotname, double scanMax ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }
         if ( znnModel != 1 ) {
            printf("\n\n *** can't do this in Znn model 2.\n\n") ;
            return false ;
         }

         double bestVal = rrv_mu_znn_sb->getVal() ;

      //--- Profile likelihood for Z to nunu SB yield.

         if ( scanMax > 0 ) { rrv_mu_znn_sb->setMax(scanMax) ; }

         ProfileLikelihoodCalculator plc_znn_sb( *dsObserved, *likelihood, RooArgSet( *rrv_mu_znn_sb ) ) ;
         plc_znn_sb.SetTestSize(0.3173) ;
         ConfInterval* plinterval_znn_sb = plc_znn_sb.GetInterval() ;
         znnSbLow  = ((LikelihoodInterval*) plinterval_znn_sb)->LowerLimit(*rrv_mu_znn_sb) ;
         znnSbHigh = ((LikelihoodInterval*) plinterval_znn_sb)->UpperLimit(*rrv_mu_znn_sb) ;
         printf("\n\n") ;
         printf("    znn, SB 68.3%% CL interval  [%5.1f, %5.1f]\n\n", znnSbLow, znnSbHigh ) ;
         printf("    znn SB = %5.1f  +%5.1f  -%5.1f\n\n", bestVal, znnSbHigh-bestVal, bestVal-znnSbLow ) ;

         if ( makePlot ) {
            TCanvas* plcplot_znn_sb = new TCanvas("plcplot_znn_sb", "znn sb, Profile likelihood", 500, 400 ) ;
            LikelihoodIntervalPlot plotInt_znn_sb((LikelihoodInterval*)plinterval_znn_sb);
            plotInt_znn_sb.Draw() ;
            gPad->SetGridy(1) ;
            plcplot_znn_sb->SaveAs( plotname ) ;
         }

         varsAtFitVals = false ;

         delete plinterval_znn_sb ; // can I safely do this???

         return true ;

     }


  //===================================================================

   ra2bRoostatsClass4ln::~ra2bRoostatsClass4ln() {

      if ( !initialized ) return ;

//    delete dsObserved ;
//    //delete plinterval_susy_sig ;

//    delete rv_Nsig ;
//    delete rv_Na ;
//    delete rv_Nd ;
//    delete rv_Nsb1 ;
//    delete rv_Nsb2 ;
//    delete rv_Nsb3 ;
//    delete rv_Nsb4 ;
//    delete rv_Nsb5 ;
//    delete rv_Nlsb1 ;
//    delete rv_Nlsb2 ;
//    delete rv_Nlsb3 ;
//    delete rv_Nlsb4 ;
//    delete rv_Nlsb5 ;
//    delete rv_Nslsig1 ;
//    delete rv_Nslsig2 ;
//    delete rv_Nslsig3 ;
//    delete rv_Nslsig4 ;
//    delete rv_Nslsig5 ;
//    delete rv_Nslsb1 ;
//    delete rv_Nslsb2 ;
//    delete rv_Nslsb3 ;
//    delete rv_Nslsb4 ;
//    delete rv_Nslsb5 ;
//    delete rv_Nslmsb1 ;
//    delete rv_Nslmsb2 ;
//    delete rv_Nslmsb3 ;
//    delete rv_Nslmsb4 ;
//    delete rv_Nslmsb5 ;
//    delete rv_mu_ttbar_sig ;
//    delete rv_mu_qcd_sig ;
//    delete rv_mu_ttbar_sb ;
//    delete rv_mu_qcd_sb ;
//    delete rv_mu_ew_sig    ;
//    delete rv_mu_susy_sig  ;
//    delete rv_mu_ew_sb1 ;
//    delete rv_mu_ew_sb2 ;
//    delete rv_mu_ew_sb3 ;
//    delete rv_mu_ew_sb4 ;
//    delete rv_mu_ew_sb5 ;

//    delete rv_mu_susy_sb1 ;
//    delete rv_mu_susy_sb2 ;
//    delete rv_mu_susy_sb3 ;
//    delete rv_mu_susy_sb4 ;
//    delete rv_mu_susy_sb5 ;
//    delete rv_mu_ttbar_a ;
//    delete rv_mu_qcd_a   ;
//    delete rv_mu_ew_a    ;
//    delete rv_mu_susy_a  ;
//    delete rv_mu_ttbar_d ;
//    delete rv_mu_qcd_d   ;
//    delete rv_mu_ew_d    ;
//    delete rv_mu_susy_d  ;
//    delete rv_mu_qcd_lsb1 ;
//    delete rv_mu_qcd_lsb2 ;
//    delete rv_mu_qcd_lsb3 ;
//    delete rv_mu_qcd_lsb4 ;
//    delete rv_mu_qcd_lsb5 ;
//    delete rv_mu_qcdmc_sig ;
//    delete rv_mu_qcdmc_sb  ;
//    delete rv_mu_qcdmc_a   ;
//    delete rv_mu_qcdmc_d   ;

//    delete rv_mu_sl_ttbar_sig1 ;
//    delete rv_mu_sl_ttbar_sig2 ;
//    delete rv_mu_sl_ttbar_sig3 ;
//    delete rv_mu_sl_ttbar_sig4 ;
//    delete rv_mu_sl_ttbar_sig5 ;
//    delete rv_mu_sl_ttbar_sb1 ;
//    delete rv_mu_sl_ttbar_sb2 ;
//    delete rv_mu_sl_ttbar_sb3 ;
//    delete rv_mu_sl_ttbar_sb4 ;
//    delete rv_mu_sl_ttbar_sb5 ;
//    delete rv_mu_sl_ttbar_msb1 ;
//    delete rv_mu_sl_ttbar_msb2 ;
//    delete rv_mu_sl_ttbar_msb3 ;
//    delete rv_mu_sl_ttbar_msb4 ;
//    delete rv_mu_sl_ttbar_msb5 ;
//    delete rv_mu_sl_ew_sig1 ;
//    delete rv_mu_sl_ew_sig2 ;
//    delete rv_mu_sl_ew_sig3 ;
//    delete rv_mu_sl_ew_sig4 ;
//    delete rv_mu_sl_ew_sig5 ;
//    delete rv_mu_sl_ew_sb1 ;
//    delete rv_mu_sl_ew_sb2 ;
//    delete rv_mu_sl_ew_sb3 ;
//    delete rv_mu_sl_ew_sb4 ;
//    delete rv_mu_sl_ew_sb5 ;
//    delete rv_mu_sl_ew_msb1 ;
//    delete rv_mu_sl_ew_msb2 ;
//    delete rv_mu_sl_ew_msb3 ;
//    delete rv_mu_sl_ew_msb4 ;
//    delete rv_mu_sl_ew_msb5 ;

//    delete rv_mu_sl_susy_sig1 ;
//    delete rv_mu_sl_susy_sig2 ;
//    delete rv_mu_sl_susy_sig3;
//    delete rv_mu_sl_susy_sig4 ;
//    delete rv_mu_sl_susy_sig5 ;
//    delete rv_mu_sl_susy_sb1 ;
//    delete rv_mu_sl_susy_sb2 ;
//    delete rv_mu_sl_susy_sb3 ;
//    delete rv_mu_sl_susy_sb4 ;
//    delete rv_mu_sl_susy_sb5 ;
//    delete rv_mu_sl_susy_msb1 ;
//    delete rv_mu_sl_susy_msb2 ;
//    delete rv_mu_sl_susy_msb3 ;
//    delete rv_mu_sl_susy_msb4 ;
//    delete rv_mu_sl_susy_msb5 ;

//    delete rv_mu_sl_ttbar_sb ;
//    delete rv_mu_sl_ttbar_sig ;
//    delete rv_mu_qcd_lsb ;

//    delete rv_f_qcd_lsb1 ;
//    delete rv_f_qcd_lsb2 ;
//    delete rv_f_qcd_lsb3 ;
//    delete rv_f_qcd_lsb4 ;
//    delete rv_f_qcd_lsb5 ;
//    delete rv_mu_qcd_sb1 ;
//    delete rv_mu_qcd_sb2 ;
//    delete rv_mu_qcd_sb3 ;
//    delete rv_mu_qcd_sb4 ;
//    delete rv_mu_qcd_sb5 ;
//    delete rv_mu_sl_ttbar1 ;
//    delete rv_mu_sl_ttbar2 ;
//    delete rv_mu_sl_ttbar3 ;
//    delete rv_mu_sl_ttbar4 ;
//    delete rv_mu_sl_ttbar5 ;
//    delete rv_mu_sl_ttbar ;
//    delete rv_f_sl_ttbar1 ;
//    delete rv_f_sl_ttbar2 ;
//    delete rv_f_sl_ttbar3 ;
//    delete rv_f_sl_ttbar4 ;
//    delete rv_f_sl_ttbar5 ;
//    delete rv_mu_ttbar_sb1 ;
//    delete rv_mu_ttbar_sb2 ;
//    delete rv_mu_ttbar_sb3 ;
//    delete rv_mu_ttbar_sb4 ;
//    delete rv_mu_ttbar_sb5 ;
//    delete rv_n_sig ;
//    delete rv_n_sb1 ;
//    delete rv_n_sb2 ;
//    delete rv_n_sb3 ;
//    delete rv_n_sb4 ;
//    delete rv_n_sb5 ;
//    delete rv_n_a ;
//    delete rv_n_d ;
//    delete rv_n_lsb1 ;
//    delete rv_n_lsb2 ;
//    delete rv_n_lsb3 ;
//    delete rv_n_lsb4 ;
//    delete rv_n_lsb5 ;
//    delete rv_n_sl_sig1 ;
//    delete rv_n_sl_sig2 ;
//    delete rv_n_sl_sig3 ;
//    delete rv_n_sl_sig4 ;
//    delete rv_n_sl_sig5 ;
//    delete rv_n_sl_sb1 ;
//    delete rv_n_sl_sb2 ;
//    delete rv_n_sl_sb3 ;
//    delete rv_n_sl_sb4 ;
//    delete rv_n_sl_sb5 ;
//    delete rv_n_sl_msb1 ;
//    delete rv_n_sl_msb2 ;
//    delete rv_n_sl_msb3 ;
//    delete rv_n_sl_msb4 ;
//    delete rv_n_sl_msb5 ;

//    delete pdf_Nsig ;
//    delete pdf_Na ;
//    delete pdf_Nd ;
//    delete pdf_Nsb1 ;
//    delete pdf_Nsb2 ;
//    delete pdf_Nsb3 ;
//    delete pdf_Nsb4 ;
//    delete pdf_Nsb5 ;
//    delete pdf_Nlsb1 ;
//    delete pdf_Nlsb2 ;
//    delete pdf_Nlsb3 ;
//    delete pdf_Nlsb4 ;
//    delete pdf_Nlsb5 ;
//    delete pdf_Nsl_sig1 ;
//    delete pdf_Nsl_sig2 ;
//    delete pdf_Nsl_sig3 ;
//    delete pdf_Nsl_sig4 ;
//    delete pdf_Nsl_sig5 ;
//    delete pdf_Nsl_sb1 ;
//    delete pdf_Nsl_sb2 ;
//    delete pdf_Nsl_sb3 ;
//    delete pdf_Nsl_sb4 ;
//    delete pdf_Nsl_sb5 ;
//    delete pdf_Nsl_msb1 ;
//    delete pdf_Nsl_msb2 ;
//    delete pdf_Nsl_msb3 ;
//    delete pdf_Nsl_msb4 ;
//    delete pdf_Nsl_msb5 ;
//    delete pdf_Nqcdmc_sig ;
//    delete pdf_Nqcdmc_sb  ;
//    delete pdf_Nqcdmc_a   ;
//    delete pdf_Nqcdmc_d   ;

//    delete likelihood ;

//    delete workspace ;



   }



  //===================================================================

    bool ra2bRoostatsClass4ln::initialize( const char* infile ) {


       printf( "\n\n Opening input file : %s\n\n", infile ) ;

       FILE* infp ;
       if ( (infp=fopen( infile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
          return false ;
       }

       sprintf( initializeFile, "%s", infile ) ;


       int    Nsig                  ;
       int    Nsb                   ;
       int    Nsig_sl               ;
       int    Nsb_sl                ;
       int    Nsig_ldp              ;
       int    Nsb_ldp               ;
       int    Nlsb                  ;
       int    Nlsb_ldp              ;
       int    Nlsb_0b               ;
       int    Nlsb_0b_ldp           ;
       float  Nqcdmc_sig            ;
       float  Nqcdmc_sig_err        ;
       float  Nqcdmc_sb             ;
       float  Nqcdmc_sb_err         ;
       float  Nqcdmc_sig_sl         ;
       float  Nqcdmc_sig_sl_err     ;
       float  Nqcdmc_sb_sl          ;
       float  Nqcdmc_sb_sl_err      ;
       float  Nqcdmc_sig_ldp        ;
       float  Nqcdmc_sig_ldp_err    ;
       float  Nqcdmc_sb_ldp         ;
       float  Nqcdmc_sb_ldp_err     ;
       float  Nqcdmc_lsb            ;
       float  Nqcdmc_lsb_err        ;
       float  Nqcdmc_lsb_ldp        ;
       float  Nqcdmc_lsb_ldp_err    ;
       float  Nqcdmc_lsb_0b         ;
       float  Nqcdmc_lsb_0b_err     ;
       float  Nqcdmc_lsb_0b_ldp     ;
       float  Nqcdmc_lsb_0b_ldp_err ;
       float  Nttbarmc_sig          ;
       float  Nttbarmc_sb           ;
       float  Nttbarmc_sig_sl       ;
       float  Nttbarmc_sb_sl        ;
       float  Nttbarmc_sig_ldp      ;
       float  Nttbarmc_sb_ldp       ;
       float  Nttbarmc_lsb          ;
       float  Nttbarmc_lsb_ldp      ;
       float  NWJmc_sig             ;
       float  NWJmc_sb              ;
       float  NWJmc_sig_sl          ;
       float  NWJmc_sb_sl           ;
       float  NWJmc_sig_ldp         ;
       float  NWJmc_sb_ldp          ;
       float  NWJmc_lsb             ;
       float  NWJmc_lsb_ldp         ;
       float  NZnnmc_sig            ;
       float  NZnnmc_sb             ;
       float  NZnnmc_sig_sl         ;
       float  NZnnmc_sb_sl          ;
       float  NZnnmc_sig_ldp        ;
       float  NZnnmc_sb_ldp         ;
       float  NZnnmc_lsb            ;
       float  NZnnmc_lsb_ldp        ;
       float  NEwomc_sig            ;
       float  NEwomc_sb             ;
       float  NEwomc_sig_sl         ;
       float  NEwomc_sb_sl          ;
       float  NEwomc_sig_ldp        ;
       float  NEwomc_sb_ldp         ;
       float  NEwomc_lsb            ;
       float  NEwomc_lsb_ldp        ;
       float  Nsusymc_sig           ;
       float  Nsusymc_sb            ;
       float  Nsusymc_sig_sl        ;
       float  Nsusymc_sb_sl         ;
       float  Nsusymc_sig_ldp       ;
       float  Nsusymc_sb_ldp        ;
       float  Nsusymc_lsb           ;
       float  Nsusymc_lsb_ldp       ;
       float  Nsusymc_lsb_0b        ;
       float  Nsusymc_lsb_0b_ldp    ;
       int    Nhtonlytrig_lsb_0b      ;
       int    Nhtonlytrig_lsb_0b_ldp  ;
       int    Nsb_ee                  ;
       int    Nsig_ee                 ;
       int    Nsb_mm                  ;
       int    Nsig_mm                 ;

       float  sf_mc            ;
       float  sf_mc_err        ;
       float  sf_qcd_sb        ;
       float  sf_qcd_sb_err    ;
       float  sf_qcd_sig       ;
       float  sf_qcd_sig_err   ;
       float  sf_ttwj_sig      ;
       float  sf_ttwj_sig_err  ;
       float  sf_ee            ;
       float  sf_ee_err        ;
       float  sf_mm            ;
       float  sf_mm_err        ;











       //--- read in description line.
       printf("\n\n") ;
       char c(0) ;
       while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
       printf("\n\n") ;


       char label[1000] ;

      //--- Inputs generated with gen_roostats_input4.c
      //    The order here must be consistent with the order there!

       fscanf( infp, "%s %g", label, &EffScaleFactor        ) ;   printf( "%s %g\n", label, EffScaleFactor        ) ;
       fscanf( infp, "%s %g", label, &EffScaleFactorErr     ) ;   printf( "%s %g\n", label, EffScaleFactorErr     ) ;
       fscanf( infp, "%s %d", label, &Nsig                  ) ;   printf( "%s %d\n", label, Nsig                  ) ;
       fscanf( infp, "%s %d", label, &Nsb                   ) ;   printf( "%s %d\n", label, Nsb                   ) ;
       fscanf( infp, "%s %d", label, &Nsig_sl               ) ;   printf( "%s %d\n", label, Nsig_sl               ) ;
       fscanf( infp, "%s %d", label, &Nsb_sl                ) ;   printf( "%s %d\n", label, Nsb_sl                ) ;
       fscanf( infp, "%s %d", label, &Nsig_ldp              ) ;   printf( "%s %d\n", label, Nsig_ldp              ) ;
       fscanf( infp, "%s %d", label, &Nsb_ldp               ) ;   printf( "%s %d\n", label, Nsb_ldp               ) ;
       fscanf( infp, "%s %d", label, &Nlsb                  ) ;   printf( "%s %d\n", label, Nlsb                  ) ;
       fscanf( infp, "%s %d", label, &Nlsb_ldp              ) ;   printf( "%s %d\n", label, Nlsb_ldp              ) ;
       fscanf( infp, "%s %d", label, &Nlsb_0b               ) ;   printf( "%s %d\n", label, Nlsb_0b               ) ;
       fscanf( infp, "%s %d", label, &Nlsb_0b_ldp           ) ;   printf( "%s %d\n", label, Nlsb_0b_ldp           ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig            ) ;   printf( "%s %g\n", label, Nqcdmc_sig            ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_err        ) ;   printf( "%s %g\n", label, Nqcdmc_sig_err        ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb             ) ;   printf( "%s %g\n", label, Nqcdmc_sb             ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_err         ) ;   printf( "%s %g\n", label, Nqcdmc_sb_err         ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_sl         ) ;   printf( "%s %g\n", label, Nqcdmc_sig_sl         ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_sl_err     ) ;   printf( "%s %g\n", label, Nqcdmc_sig_sl_err     ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_sl          ) ;   printf( "%s %g\n", label, Nqcdmc_sb_sl          ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_sl_err      ) ;   printf( "%s %g\n", label, Nqcdmc_sb_sl_err      ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_ldp        ) ;   printf( "%s %g\n", label, Nqcdmc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_ldp_err    ) ;   printf( "%s %g\n", label, Nqcdmc_sig_ldp_err    ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_ldp         ) ;   printf( "%s %g\n", label, Nqcdmc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_ldp_err     ) ;   printf( "%s %g\n", label, Nqcdmc_sb_ldp_err     ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb            ) ;   printf( "%s %g\n", label, Nqcdmc_lsb            ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_err        ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_err        ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_ldp        ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_ldp        ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_ldp_err    ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_ldp_err    ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_0b         ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_0b         ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_0b_err     ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_0b_err     ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_0b_ldp     ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_0b_ldp     ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_0b_ldp_err ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_0b_ldp_err ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sig          ) ;   printf( "%s %g\n", label, Nttbarmc_sig          ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sb           ) ;   printf( "%s %g\n", label, Nttbarmc_sb           ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sig_sl       ) ;   printf( "%s %g\n", label, Nttbarmc_sig_sl       ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sb_sl        ) ;   printf( "%s %g\n", label, Nttbarmc_sb_sl        ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sig_ldp      ) ;   printf( "%s %g\n", label, Nttbarmc_sig_ldp      ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sb_ldp       ) ;   printf( "%s %g\n", label, Nttbarmc_sb_ldp       ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_lsb          ) ;   printf( "%s %g\n", label, Nttbarmc_lsb          ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_lsb_ldp      ) ;   printf( "%s %g\n", label, Nttbarmc_lsb_ldp      ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig             ) ;   printf( "%s %g\n", label, NWJmc_sig             ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb              ) ;   printf( "%s %g\n", label, NWJmc_sb              ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig_sl          ) ;   printf( "%s %g\n", label, NWJmc_sig_sl          ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb_sl           ) ;   printf( "%s %g\n", label, NWJmc_sb_sl           ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig_ldp         ) ;   printf( "%s %g\n", label, NWJmc_sig_ldp         ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb_ldp          ) ;   printf( "%s %g\n", label, NWJmc_sb_ldp          ) ;
       fscanf( infp, "%s %g", label, &NWJmc_lsb             ) ;   printf( "%s %g\n", label, NWJmc_lsb             ) ;
       fscanf( infp, "%s %g", label, &NWJmc_lsb_ldp         ) ;   printf( "%s %g\n", label, NWJmc_lsb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig            ) ;   printf( "%s %g\n", label, NZnnmc_sig            ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb             ) ;   printf( "%s %g\n", label, NZnnmc_sb             ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig_sl         ) ;   printf( "%s %g\n", label, NZnnmc_sig_sl         ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb_sl          ) ;   printf( "%s %g\n", label, NZnnmc_sb_sl          ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig_ldp        ) ;   printf( "%s %g\n", label, NZnnmc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb_ldp         ) ;   printf( "%s %g\n", label, NZnnmc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_lsb            ) ;   printf( "%s %g\n", label, NZnnmc_lsb            ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_lsb_ldp        ) ;   printf( "%s %g\n", label, NZnnmc_lsb_ldp        ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig            ) ;   printf( "%s %g\n", label, NEwomc_sig            ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb             ) ;   printf( "%s %g\n", label, NEwomc_sb             ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig_sl         ) ;   printf( "%s %g\n", label, NEwomc_sig_sl         ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb_sl          ) ;   printf( "%s %g\n", label, NEwomc_sb_sl          ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig_ldp        ) ;   printf( "%s %g\n", label, NEwomc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb_ldp         ) ;   printf( "%s %g\n", label, NEwomc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NEwomc_lsb            ) ;   printf( "%s %g\n", label, NEwomc_lsb            ) ;
       fscanf( infp, "%s %g", label, &NEwomc_lsb_ldp        ) ;   printf( "%s %g\n", label, NEwomc_lsb_ldp        ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sig           ) ;   printf( "%s %g\n", label, Nsusymc_sig           ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sb            ) ;   printf( "%s %g\n", label, Nsusymc_sb            ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sig_sl        ) ;   printf( "%s %g\n", label, Nsusymc_sig_sl        ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sb_sl         ) ;   printf( "%s %g\n", label, Nsusymc_sb_sl         ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sig_ldp       ) ;   printf( "%s %g\n", label, Nsusymc_sig_ldp       ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sb_ldp        ) ;   printf( "%s %g\n", label, Nsusymc_sb_ldp        ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_lsb           ) ;   printf( "%s %g\n", label, Nsusymc_lsb           ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_lsb_ldp       ) ;   printf( "%s %g\n", label, Nsusymc_lsb_ldp       ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_lsb_0b        ) ;   printf( "%s %g\n", label, Nsusymc_lsb_0b        ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_lsb_0b_ldp    ) ;   printf( "%s %g\n", label, Nsusymc_lsb_0b_ldp    ) ;
       fscanf( infp, "%s %d", label, &Nhtonlytrig_lsb_0b        ) ;   printf( "%s %d\n", label, Nhtonlytrig_lsb_0b        ) ;
       fscanf( infp, "%s %d", label, &Nhtonlytrig_lsb_0b_ldp    ) ;   printf( "%s %d\n", label, Nhtonlytrig_lsb_0b_ldp    ) ;
       fscanf( infp, "%s %g", label, &DataLumi                  ) ;   printf( "%s %g\n", label, DataLumi                  ) ;
       fscanf( infp, "%s %d", label, &Nsb_ee                    ) ;   printf( "%s %d\n", label, Nsb_ee                    ) ;
       fscanf( infp, "%s %d", label, &Nsig_ee                   ) ;   printf( "%s %d\n", label, Nsig_ee                   ) ;
       fscanf( infp, "%s %d", label, &Nsb_mm                    ) ;   printf( "%s %d\n", label, Nsb_mm                    ) ;
       fscanf( infp, "%s %d", label, &Nsig_mm                   ) ;   printf( "%s %d\n", label, Nsig_mm                   ) ;
       fscanf( infp, "%s %g", label, &acc_ee_mean               ) ;   printf( "%s %g\n", label, acc_ee_mean               ) ;
       fscanf( infp, "%s %g", label, &acc_ee_err                ) ;   printf( "%s %g\n", label, acc_ee_err                ) ;
       fscanf( infp, "%s %g", label, &acc_mm_mean               ) ;   printf( "%s %g\n", label, acc_mm_mean               ) ;
       fscanf( infp, "%s %g", label, &acc_mm_err                ) ;   printf( "%s %g\n", label, acc_mm_err                ) ;
       fscanf( infp, "%s %g", label, &eff_ee_mean               ) ;   printf( "%s %g\n", label, eff_ee_mean               ) ;
       fscanf( infp, "%s %g", label, &eff_ee_err                ) ;   printf( "%s %g\n", label, eff_ee_err                ) ;
       fscanf( infp, "%s %g", label, &eff_mm_mean               ) ;   printf( "%s %g\n", label, eff_mm_mean               ) ;
       fscanf( infp, "%s %g", label, &eff_mm_err                ) ;   printf( "%s %g\n", label, eff_mm_err                ) ;
       fscanf( infp, "%s %g", label, &Ztoll_lumi                ) ;   printf( "%s %g\n", label, Ztoll_lumi                ) ;
       fscanf( infp, "%s %g", label, &knn_sig_mean              ) ;   printf( "%s %g\n", label, knn_sig_mean              ) ;
       fscanf( infp, "%s %g", label, &knn_sig_err               ) ;   printf( "%s %g\n", label, knn_sig_err               ) ;
       fscanf( infp, "%s %g", label, &knn_sb_mean               ) ;   printf( "%s %g\n", label, knn_sb_mean               ) ;
       fscanf( infp, "%s %g", label, &knn_sb_err                ) ;   printf( "%s %g\n", label, knn_sb_err                ) ;
       fscanf( infp, "%s %g", label, &fsig_ee_mean              ) ;   printf( "%s %g\n", label, fsig_ee_mean              ) ;
       fscanf( infp, "%s %g", label, &fsig_ee_err               ) ;   printf( "%s %g\n", label, fsig_ee_err               ) ;
       fscanf( infp, "%s %g", label, &fsig_mm_mean              ) ;   printf( "%s %g\n", label, fsig_mm_mean              ) ;
       fscanf( infp, "%s %g", label, &fsig_mm_err               ) ;   printf( "%s %g\n", label, fsig_mm_err               ) ;

       fscanf( infp, "%s %g", label, &sf_mc                     ) ;   printf( "%s %g\n", label, sf_mc                     ) ;
       fscanf( infp, "%s %g", label, &sf_mc_err                 ) ;   printf( "%s %g\n", label, sf_mc_err                 ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb                 ) ;   printf( "%s %g\n", label, sf_qcd_sb                 ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb_err             ) ;   printf( "%s %g\n", label, sf_qcd_sb_err             ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig                ) ;   printf( "%s %g\n", label, sf_qcd_sig                ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig_err            ) ;   printf( "%s %g\n", label, sf_qcd_sig_err            ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig               ) ;   printf( "%s %g\n", label, sf_ttwj_sig               ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig_err           ) ;   printf( "%s %g\n", label, sf_ttwj_sig_err           ) ;
       fscanf( infp, "%s %g", label, &sf_ee                     ) ;   printf( "%s %g\n", label, sf_ee                     ) ;
       fscanf( infp, "%s %g", label, &sf_ee_err                 ) ;   printf( "%s %g\n", label, sf_ee_err                 ) ;
       fscanf( infp, "%s %g", label, &sf_mm                     ) ;   printf( "%s %g\n", label, sf_mm                     ) ;
       fscanf( infp, "%s %g", label, &sf_mm_err                 ) ;   printf( "%s %g\n", label, sf_mm_err                 ) ;



       printf("\n Done reading in %s\n\n", infile ) ;
       fclose( infp ) ;


       //+++++ Owen, Aug 6.  set all lsf_ factors to 1.  Input MC counts are now all weighted.
       //                    Calculations of stat errors below are now wrong.  Fix later.
       lsf_WJmc = 1.0 ;
       lsf_Znnmc = 1.0 ;
       lsf_Ewomc = 0.0 ; // don't use ewo.
       sf_ttbarmc = 1.0 ;

       //--- Print out a nice summary of the inputs.

       float Nsm_sig         = Nttbarmc_sig          +  lsf_WJmc*NWJmc_sig          +  Nqcdmc_sig          +  lsf_Znnmc*NZnnmc_sig          +  lsf_Ewomc*NEwomc_sig         ;
       float Nsm_sb          = Nttbarmc_sb           +  lsf_WJmc*NWJmc_sb           +  Nqcdmc_sb           +  lsf_Znnmc*NZnnmc_sb           +  lsf_Ewomc*NEwomc_sb          ;
       float Nsm_sig_sl      = Nttbarmc_sig_sl       +  lsf_WJmc*NWJmc_sig_sl       +  Nqcdmc_sig_sl       +  lsf_Znnmc*NZnnmc_sig_sl       +  lsf_Ewomc*NEwomc_sig_sl      ;
       float Nsm_sb_sl       = Nttbarmc_sb_sl        +  lsf_WJmc*NWJmc_sb_sl        +  Nqcdmc_sb_sl        +  lsf_Znnmc*NZnnmc_sb_sl        +  lsf_Ewomc*NEwomc_sb_sl       ;
       float Nsm_sig_ldp     = Nttbarmc_sig_ldp      +  lsf_WJmc*NWJmc_sig_ldp      +  Nqcdmc_sig_ldp      +  lsf_Znnmc*NZnnmc_sig_ldp      +  lsf_Ewomc*NEwomc_sig_ldp     ;
       float Nsm_sb_ldp      = Nttbarmc_sb_ldp       +  lsf_WJmc*NWJmc_sb_ldp       +  Nqcdmc_sb_ldp       +  lsf_Znnmc*NZnnmc_sb_ldp       +  lsf_Ewomc*NEwomc_sb_ldp      ;
       float Nsm_lsb         = Nttbarmc_lsb          +  lsf_WJmc*NWJmc_lsb          +  Nqcdmc_lsb          +  lsf_Znnmc*NZnnmc_lsb          +  lsf_Ewomc*NEwomc_lsb         ;
       float Nsm_lsb_ldp     = Nttbarmc_lsb_ldp      +  lsf_WJmc*NWJmc_lsb_ldp      +  Nqcdmc_lsb_ldp      +  lsf_Znnmc*NZnnmc_lsb_ldp      +  lsf_Ewomc*NEwomc_lsb_ldp     ;
       float Nsm_lsb_0b      =                                                         Nqcdmc_lsb_0b                                                                        ;
       float Nsm_lsb_0b_ldp  =                                                         Nqcdmc_lsb_0b_ldp                                                                    ;



       //--- QCD min Delta phi N  pass / fail ratios.

       float Rqcd_sig        =  0. ;
       if ( Nqcdmc_sig_ldp > 0. ) { Rqcd_sig        =  Nqcdmc_sig        / Nqcdmc_sig_ldp ;  }
       float Rqcd_sig_err    =  0. ;
       if ( Nqcdmc_sig > 0. && Nqcdmc_sig_ldp > 0. ) { Rqcd_sig_err  =  Rqcd_sig * sqrt( pow(Nqcdmc_sig_err/Nqcdmc_sig,2) + pow(Nqcdmc_sig_ldp_err/Nqcdmc_sig_ldp,2) ) ; }


       float Rqcd_sb        =  0. ;
       if ( Nqcdmc_sb_ldp > 0. ) { Rqcd_sb        =  Nqcdmc_sb        / Nqcdmc_sb_ldp ;  }
       float Rqcd_sb_err    =  0. ;
       if ( Nqcdmc_sb > 0. && Nqcdmc_sb_ldp > 0. ) { Rqcd_sb_err  =  Rqcd_sb * sqrt( pow(Nqcdmc_sb_err/Nqcdmc_sb,2) + pow(Nqcdmc_sb_ldp_err/Nqcdmc_sb_ldp,2) ) ; }


       float Rqcd_lsb        =  0. ;
       if ( Nqcdmc_lsb_ldp > 0. ) { Rqcd_lsb        =  Nqcdmc_lsb        / Nqcdmc_lsb_ldp ;  }
       float Rqcd_lsb_err    =  0. ;
       if ( Nqcdmc_lsb > 0. && Nqcdmc_lsb_ldp > 0. ) { Rqcd_lsb_err  =  Rqcd_lsb * sqrt( pow(Nqcdmc_lsb_err/Nqcdmc_lsb,2) + pow(Nqcdmc_lsb_ldp_err/Nqcdmc_lsb_ldp,2) ) ; }


       float Nlsb_corrected     = Nlsb     - Nttbarmc_lsb     - lsf_WJmc*NWJmc_lsb     - lsf_Znnmc*NZnnmc_lsb     - lsf_Ewomc*NEwomc_lsb ;
       float Nlsb_ldp_corrected = Nlsb_ldp - Nttbarmc_lsb_ldp - lsf_WJmc*NWJmc_lsb_ldp - lsf_Znnmc*NZnnmc_lsb_ldp - lsf_Ewomc*NEwomc_lsb_ldp ;
       float Rdata_lsb = 0. ;
       if ( Nlsb_ldp_corrected > 0 ) { Rdata_lsb = (Nlsb_corrected) / (Nlsb_ldp_corrected) ; }
       float Rdata_lsb_err = 0. ;
       if ( Nlsb_corrected > 0 && Nlsb_ldp_corrected > 0 ) { Rdata_lsb_err = Rdata_lsb * sqrt( 1.0/(Nlsb_corrected) + 1.0/(Nlsb_ldp_corrected) ) ; }

       float RdataHTOT_lsb_0b = 0. ;
       if ( Nhtonlytrig_lsb_0b_ldp > 0 ) { RdataHTOT_lsb_0b = (1.0*Nhtonlytrig_lsb_0b) / (1.0*Nhtonlytrig_lsb_0b_ldp) ; }
       float RdataHTOT_lsb_0b_err = 0. ;
       if ( Nhtonlytrig_lsb_0b_ldp > 0 && Nhtonlytrig_lsb_0b > 0 ) { RdataHTOT_lsb_0b_err = RdataHTOT_lsb_0b * sqrt( 1.0/(1.0*Nhtonlytrig_lsb_0b) + 1.0/(1.0*Nhtonlytrig_lsb_0b_ldp) ) ; }


       float Rqcd_lsb_0b        =  0. ;
       if ( Nqcdmc_lsb_0b_ldp > 0. ) { Rqcd_lsb_0b        =  Nqcdmc_lsb_0b        / Nqcdmc_lsb_0b_ldp ;  }
       float Rqcd_lsb_0b_err    =  0. ;
       if ( Nqcdmc_lsb_0b > 0. && Nqcdmc_lsb_0b_ldp > 0. ) { Rqcd_lsb_0b_err  =  Rqcd_lsb_0b * sqrt( pow(Nqcdmc_lsb_0b_err/Nqcdmc_lsb_0b,2) + pow(Nqcdmc_lsb_0b_ldp_err/Nqcdmc_lsb_0b_ldp,2) ) ; }


       float Rdata_lsb_0b = 0. ;
       if ( Nlsb_0b_ldp > 0 ) { Rdata_lsb_0b = (1.0*Nlsb_0b) / (1.0*Nlsb_0b_ldp) ; }
       float Rdata_lsb_0b_err = 0. ;
       if ( Nlsb_0b > 0 && Nlsb_0b_ldp > 0 ) { Rdata_lsb_0b_err = Rdata_lsb_0b * sqrt( 1.0/(1.0*Nlsb_0b) + 1.0/(1.0*Nlsb_0b_ldp) ) ; }


       //--- ttbar + Wjets  MET  SIG / SB ratios

       float Rttwj = 0. ;
       if ( (Nttbarmc_sb+lsf_WJmc*NWJmc_sb) > 0. ) { Rttwj    =   (Nttbarmc_sig+lsf_WJmc*NWJmc_sig) / (Nttbarmc_sb+lsf_WJmc*NWJmc_sb) ; }
       float Rttwj_err2 = pow( lsf_WJmc/ (Nttbarmc_sb+lsf_WJmc*NWJmc_sb),2)*NWJmc_sig
                        + pow( (Nttbarmc_sig+lsf_WJmc*NWJmc_sig) * lsf_WJmc / pow((Nttbarmc_sb+lsf_WJmc*NWJmc_sb),2), 2)*NWJmc_sb ;
       float Rttwj_err = sqrt(Rttwj_err2) ;

       float Rttwj_sl = 0. ;
       if ( (Nttbarmc_sb_sl+lsf_WJmc*NWJmc_sb_sl) > 0. ) { Rttwj_sl    =   (Nttbarmc_sig_sl+lsf_WJmc*NWJmc_sig_sl) / (Nttbarmc_sb_sl+lsf_WJmc*NWJmc_sb_sl) ; }
       float Rttwj_sl_err2 = pow( lsf_WJmc / (Nttbarmc_sb_sl+lsf_WJmc*NWJmc_sb_sl), 2) * NWJmc_sig_sl
                           + pow( lsf_WJmc * (Nttbarmc_sig_sl+lsf_WJmc*NWJmc_sig_sl) / pow((Nttbarmc_sb_sl+lsf_WJmc*NWJmc_sb_sl),2), 2)*NWJmc_sb_sl ;
       float Rttwj_sl_err = sqrt(Rttwj_sl_err2) ;

       float Rdata_sl = 0. ;
       if ( Nsb_sl > 0 ) { Rdata_sl = (1.0*Nsig_sl)/(1.0*Nsb_sl) ; }
       float Rdata_sl_err = 0. ;
       if ( Nsig_sl>0 && Nsb_sl>0 ) { Rdata_sl_err = Rdata_sl * sqrt( 1.0/(1.0*Nsig_sl) + 1.0/(1.0*Nsb_sl) ) ; }


       //--- Simple MC closure tests.

       float comp_mc_qcd_sb = Nqcdmc_sb_ldp * ( Nqcdmc_lsb_0b / Nqcdmc_lsb_0b_ldp ) ;
       float comp_mc_qcd_sb_err = Nqcdmc_sb_ldp_err * ( Nqcdmc_lsb_0b / Nqcdmc_lsb_0b_ldp ) ;

       float comp_mc_qcd_sig = Nqcdmc_sig_ldp * ( Nqcdmc_lsb_0b / Nqcdmc_lsb_0b_ldp ) ;
       float comp_mc_qcd_sig_err = Nqcdmc_sig_ldp_err * ( Nqcdmc_lsb_0b / Nqcdmc_lsb_0b_ldp ) ;

       float comp_mc_ttwj_sig = (Nttbarmc_sb + lsf_WJmc*NWJmc_sb) * ( (Nttbarmc_sig_sl + lsf_WJmc*NWJmc_sig_sl) / (Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl) ) ;

       //--- below ignores ttbar errors, uses sqrt(N) on raw WJ counts.
       float comp_mc_ttwj_sig_err2 = pow( lsf_WJmc*( (Nttbarmc_sig_sl + lsf_WJmc*NWJmc_sig_sl) / (Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl) ),2)*NWJmc_sb
                                   + pow( (Nttbarmc_sb + lsf_WJmc*NWJmc_sb) *(lsf_WJmc/ (Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl) ), 2)*NWJmc_sig_sl
                                   + pow( (Nttbarmc_sb + lsf_WJmc*NWJmc_sb) *  (Nttbarmc_sig_sl + lsf_WJmc*NWJmc_sig_sl) / pow((Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl) ,2)*lsf_WJmc, 2 )*NWJmc_sb_sl ;
       float comp_mc_ttwj_sig_err = sqrt(comp_mc_ttwj_sig_err2) ;

       //--- Simple data calculations.


       float Nsb_ldp_corrected  = Nsb_ldp  - (Nttbarmc_sb_ldp  + lsf_WJmc*NWJmc_sb_ldp  + lsf_Znnmc*NZnnmc_sb_ldp  + lsf_Ewomc*NEwomc_sb_ldp  ) ;
       float comp_data_qcd_sb = Nsb_ldp_corrected * ( (1.0*Nlsb_0b)/(1.0*Nlsb_0b_ldp) ) ;
       float comp_data_qcd_sb_err = sqrt( Nsb_ldp + Nttbarmc_sb_ldp  + lsf_WJmc*NWJmc_sb_ldp  + lsf_Znnmc*NZnnmc_sb_ldp  + lsf_Ewomc*NEwomc_sb_ldp ) * ( (1.0*Nlsb_0b)/(1.0*Nlsb_0b_ldp) ) ;

       float Nsig_ldp_corrected = Nsig_ldp - (Nttbarmc_sig_ldp + lsf_WJmc*NWJmc_sig_ldp + lsf_Znnmc*NZnnmc_sig_ldp + lsf_Ewomc*NEwomc_sig_ldp ) ;
       float comp_data_qcd_sig = Nsig_ldp_corrected * ( (1.0*Nlsb_0b)/(1.0*Nlsb_0b_ldp) ) ;
       float comp_data_qcd_sig_err = sqrt( Nsig_ldp + Nttbarmc_sig_ldp  + lsf_WJmc*NWJmc_sig_ldp  + lsf_Znnmc*NZnnmc_sig_ldp  + lsf_Ewomc*NEwomc_sig_ldp ) * ( (1.0*Nlsb_0b)/(1.0*Nlsb_0b_ldp) ) ;

       float comp_data_ttwj_sig = (Nsb - (comp_data_qcd_sb + lsf_Znnmc*NZnnmc_sb + lsf_Ewomc*NEwomc_sb)) * ( (1.0*Nsig_sl) / (1.0*Nsb_sl) ) ;
       float comp_data_ttwj_sig_err2 = pow( ( (1.0*Nsig_sl) / (1.0*Nsb_sl) ), 2)*Nsb
                                     + pow( ( (1.0*Nsig_sl) / (1.0*Nsb_sl) ), 2)*pow(comp_data_qcd_sb_err,2)
                                     + pow( (Nsb - (comp_data_qcd_sb + lsf_Znnmc*NZnnmc_sb + lsf_Ewomc*NEwomc_sb)) / (1.0*Nsb_sl), 2 )*Nsig_sl
                                     + pow( (Nsb - (comp_data_qcd_sb + lsf_Znnmc*NZnnmc_sb + lsf_Ewomc*NEwomc_sb)) * (1.0*Nsig_sl) / pow(1.0*Nsb_sl,2), 2 ) * Nsb_sl ;
       float comp_data_ttwj_sig_err = sqrt(comp_data_ttwj_sig_err2) ;


       float comp_znn_sig_ee(2.) ;
       float comp_znn_sig_mm(2.) ;
       float comp_znn_sig(2.) ;
       float comp_znn_sb_ee(2.) ;
       float comp_znn_sb_mm(2.) ;
       float comp_znn_sb(2.) ;

       float comp_znn_sig_ee_err(0.) ;
       float comp_znn_sig_mm_err(0.) ;
       float comp_znn_sig_err(0.) ;
       float comp_znn_sb_ee_err(0.) ;
       float comp_znn_sb_mm_err(0.) ;
       float comp_znn_sb_err(0.) ;


       if ( znnModel == 1 ) {

         comp_znn_sig_ee = Nsig_ee * ( 5.95 * DataLumi * fsig_ee_mean ) / ( acc_ee_mean * eff_ee_mean * Ztoll_lumi ) ;
         comp_znn_sb_ee  = Nsb_ee  * ( 5.95 * DataLumi * fsig_ee_mean ) / ( acc_ee_mean * eff_ee_mean * Ztoll_lumi ) ;

         comp_znn_sig_mm = Nsig_mm * ( 5.95 * DataLumi * fsig_mm_mean ) / ( acc_mm_mean * eff_mm_mean * Ztoll_lumi ) ;
         comp_znn_sb_mm  = Nsb_mm  * ( 5.95 * DataLumi * fsig_mm_mean ) / ( acc_mm_mean * eff_mm_mean * Ztoll_lumi ) ;

         if ( Nsig_ee > 0 ) { comp_znn_sig_ee_err = comp_znn_sig_ee * sqrt( 1.0/(1.0*Nsig_ee) + pow(fsig_ee_err/fsig_ee_mean,2) + pow(acc_ee_err/acc_ee_mean,2) + pow(eff_ee_err/eff_ee_mean,2) ) ; }
         if ( Nsb_ee  > 0 ) { comp_znn_sb_ee_err  = comp_znn_sb_ee  * sqrt( 1.0/(1.0*Nsb_ee)  + pow(fsig_ee_err/fsig_ee_mean,2) + pow(acc_ee_err/acc_ee_mean,2) + pow(eff_ee_err/eff_ee_mean,2) ) ; }

         if ( Nsig_mm > 0 ) { comp_znn_sig_mm_err = comp_znn_sig_mm * sqrt( 1.0/(1.0*Nsig_mm) + pow(fsig_mm_err/fsig_mm_mean,2) + pow(acc_mm_err/acc_mm_mean,2) + pow(eff_mm_err/eff_mm_mean,2) ) ; }
         if ( Nsb_mm  > 0 ) { comp_znn_sb_mm_err  = comp_znn_sb_mm  * sqrt( 1.0/(1.0*Nsb_mm)  + pow(fsig_mm_err/fsig_mm_mean,2) + pow(acc_mm_err/acc_mm_mean,2) + pow(eff_mm_err/eff_mm_mean,2) ) ; }



       } else if ( znnModel == 2 ) {

         comp_znn_sig_ee = (Nsig_ee + Nsb_ee) * ( 5.95 * DataLumi * fsig_ee_mean * knn_sig_mean ) / ( acc_ee_mean * eff_ee_mean * Ztoll_lumi ) ;
         comp_znn_sb_ee  = (Nsig_ee + Nsb_ee) * ( 5.95 * DataLumi * fsig_ee_mean * knn_sb_mean  ) / ( acc_ee_mean * eff_ee_mean * Ztoll_lumi ) ;

         comp_znn_sig_mm = (Nsig_mm + Nsb_mm) * ( 5.95 * DataLumi * fsig_mm_mean * knn_sig_mean ) / ( acc_mm_mean * eff_mm_mean * Ztoll_lumi ) ;
         comp_znn_sb_mm  = (Nsig_mm + Nsb_mm) * ( 5.95 * DataLumi * fsig_mm_mean * knn_sb_mean  ) / ( acc_mm_mean * eff_mm_mean * Ztoll_lumi ) ;

         if ( (Nsig_ee+Nsb_ee) > 0 ) { comp_znn_sig_ee_err = comp_znn_sig_ee * sqrt( 1.0/(1.0*(Nsig_ee+Nsb_ee)) + pow(fsig_ee_err/fsig_ee_mean,2) + pow(acc_ee_err/acc_ee_mean,2) + pow(eff_ee_err/eff_ee_mean,2) + pow(knn_sig_err/knn_sig_mean,2) ) ; }
         if ( (Nsig_ee+Nsb_ee) > 0 ) { comp_znn_sb_ee_err  = comp_znn_sb_ee  * sqrt( 1.0/(1.0*(Nsig_ee+Nsb_ee)) + pow(fsig_ee_err/fsig_ee_mean,2) + pow(acc_ee_err/acc_ee_mean,2) + pow(eff_ee_err/eff_ee_mean,2) + pow(knn_sig_err/knn_sig_mean,2) ) ; }

         if ( (Nsig_mm+Nsb_mm) > 0 ) { comp_znn_sig_mm_err = comp_znn_sig_mm * sqrt( 1.0/(1.0*(Nsig_mm+Nsb_mm)) + pow(fsig_mm_err/fsig_mm_mean,2) + pow(acc_mm_err/acc_mm_mean,2) + pow(eff_mm_err/eff_mm_mean,2) + pow(knn_sig_err/knn_sig_mean,2) ) ; }
         if ( (Nsig_mm+Nsb_mm) > 0 ) { comp_znn_sb_mm_err  = comp_znn_sb_mm  * sqrt( 1.0/(1.0*(Nsig_mm+Nsb_mm)) + pow(fsig_mm_err/fsig_mm_mean,2) + pow(acc_mm_err/acc_mm_mean,2) + pow(eff_mm_err/eff_mm_mean,2) + pow(knn_sig_err/knn_sig_mean,2) ) ; }


       }

       //-- really dumb ave.
       comp_znn_sig = 0.5 * ( comp_znn_sig_ee + comp_znn_sig_mm ) ;
       comp_znn_sig_err = comp_znn_sig_ee_err ;
       if ( comp_znn_sig_mm_err > comp_znn_sig_ee_err ) comp_znn_sig_err = comp_znn_sig_mm_err ;

       comp_znn_sb = 0.5 * ( comp_znn_sb_ee + comp_znn_sb_mm ) ;
       comp_znn_sb_err = comp_znn_sb_ee_err ;
       if ( comp_znn_sb_mm_err > comp_znn_sb_ee_err ) comp_znn_sb_err = comp_znn_sb_mm_err ;


       printf("\n\n\n") ;

       printf("------------+----------+----------+------------------------+----------+----------+------------+------------+---------------\n") ;
       printf("   Sample   |  ttbar   |  W+jets  |          QCD           |  Z to nn |   other  |   All SM   |    Data    |     SUSY      \n") ;
       printf("------------+----------+----------+------------------------+----------+----------+------------+------------+---------------\n") ;
       printf("            |          |          |                        |          |          |            |            |               \n") ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "sig", Nttbarmc_sig, (lsf_WJmc*NWJmc_sig), Nqcdmc_sig, Nqcdmc_sig_err, lsf_Znnmc*NZnnmc_sig, lsf_Ewomc*NEwomc_sig, Nsm_sig, Nsig, Nsusymc_sig ) ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "sb", Nttbarmc_sb, (lsf_WJmc*NWJmc_sb), Nqcdmc_sb, Nqcdmc_sb_err, lsf_Znnmc*NZnnmc_sb, lsf_Ewomc*NEwomc_sb, Nsm_sb, Nsb, Nsusymc_sb ) ;
       printf("            |          |          |                        |          |          |            |            |               \n") ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "sig_sl", Nttbarmc_sig_sl, (lsf_WJmc*NWJmc_sig_sl), Nqcdmc_sig_sl, Nqcdmc_sig_sl_err, lsf_Znnmc*NZnnmc_sig_sl, lsf_Ewomc*NEwomc_sig_sl, Nsm_sig_sl, Nsig_sl, Nsusymc_sig_sl ) ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "sb_sl", Nttbarmc_sb_sl, (lsf_WJmc*NWJmc_sb_sl), Nqcdmc_sb_sl, Nqcdmc_sb_sl_err, lsf_Znnmc*NZnnmc_sb_sl, lsf_Ewomc*NEwomc_sb_sl, Nsm_sb_sl, Nsb_sl, Nsusymc_sb_sl ) ;
       printf("            |          |          |                        |          |          |            |            |               \n") ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "sig_ldp", Nttbarmc_sig_ldp, (lsf_WJmc*NWJmc_sig_ldp), Nqcdmc_sig_ldp, Nqcdmc_sig_ldp_err, lsf_Znnmc*NZnnmc_sig_ldp, lsf_Ewomc*NEwomc_sig_ldp, Nsm_sig_ldp, Nsig_ldp, Nsusymc_sig_ldp ) ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "sb_ldp", Nttbarmc_sb_ldp, (lsf_WJmc*NWJmc_sb_ldp), Nqcdmc_sb_ldp, Nqcdmc_sb_ldp_err, lsf_Znnmc*NZnnmc_sb_ldp, lsf_Ewomc*NEwomc_sb_ldp, Nsm_sb_ldp, Nsb_ldp, Nsusymc_sb_ldp ) ;
       printf("            |          |          |                        |          |          |            |            |               \n") ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "lsb", Nttbarmc_lsb, (lsf_WJmc*NWJmc_lsb), Nqcdmc_lsb, Nqcdmc_lsb_err, lsf_Znnmc*NZnnmc_lsb, lsf_Ewomc*NEwomc_lsb, Nsm_lsb, Nlsb, Nsusymc_lsb ) ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "lsb_ldp", Nttbarmc_lsb_ldp, (lsf_WJmc*NWJmc_lsb_ldp), Nqcdmc_lsb_ldp, Nqcdmc_lsb_ldp_err, lsf_Znnmc*NZnnmc_lsb_ldp, lsf_Ewomc*NEwomc_lsb_ldp, Nsm_lsb_ldp, Nlsb_ldp, Nsusymc_lsb_ldp ) ;
       printf("            |          |          |                        |          |          |            |            |               \n") ;
       printf(" %10s | -------- | -------- | %10.1f +/- %7.1f | -------- | -------- | %10.1f | %10d | %8.1f\n",
            "lsb_0b",  Nqcdmc_lsb_0b, Nqcdmc_lsb_0b_err, Nsm_lsb_0b, Nlsb_0b, Nsusymc_lsb_0b ) ;
       printf(" %10s | -------- | -------- | %10.1f +/- %7.1f | -------- | -------- | %10.1f | %10d | %8.1f\n",
            "lsb_0b_ldp",  Nqcdmc_lsb_0b_ldp, Nqcdmc_lsb_0b_ldp_err, Nsm_lsb_0b_ldp, Nlsb_0b_ldp, Nsusymc_lsb_0b_ldp ) ;
       printf("            |          |          |                        |          |          |            |            |               \n") ;
       printf("------------+----------+----------+------------------------+----------+----------+------------+------------+---------------\n") ;


       printf("\n\n\n") ;

       printf(" R QCD  :  sig    /    sig_ldp  : ( %10.1f +/- %6.1f ) / ( %10.1f +/- %7.1f )    =    %5.3f +/- %5.3f\n",
            Nqcdmc_sig, Nqcdmc_sig_err, Nqcdmc_sig_ldp, Nqcdmc_sig_ldp_err, Rqcd_sig, Rqcd_sig_err ) ;

       printf(" R QCD  :  sb     /     sb_ldp  : ( %10.1f +/- %6.1f ) / ( %10.1f +/- %7.1f )    =    %5.3f +/- %5.3f\n",
            Nqcdmc_sb, Nqcdmc_sb_err, Nqcdmc_sb_ldp, Nqcdmc_sb_ldp_err, Rqcd_sb, Rqcd_sb_err ) ;

       printf(" R QCD  :  lsb    /    lsb_ldp  : ( %10.1f +/- %6.1f ) / ( %10.1f +/- %7.1f )    =    %5.3f +/- %5.3f\n",
            Nqcdmc_lsb, Nqcdmc_lsb_err, Nqcdmc_lsb_ldp, Nqcdmc_lsb_ldp_err, Rqcd_lsb, Rqcd_lsb_err ) ;

       printf(" R data :  lsb    /    lsb_ldp  : ( %10.1f            ) / ( %10.1f             )    =    %5.3f +/- %5.3f\n",
            Nlsb_corrected, Nlsb_ldp_corrected, Rdata_lsb, Rdata_lsb_err ) ;

       printf(" R QCD  :  lsb_0b / lsb_0b_ldp  : ( %10.1f +/- %6.1f ) / ( %10.1f +/- %7.1f )    =    %5.3f +/- %5.3f\n",
            Nqcdmc_lsb_0b, Nqcdmc_lsb_0b_err, Nqcdmc_lsb_0b_ldp, Nqcdmc_lsb_0b_ldp_err, Rqcd_lsb_0b, Rqcd_lsb_0b_err ) ;

       printf(" R data :  lsb_0b / lsb_0b_ldp  : ( %8d              ) / ( %8d               )    =    %5.3f +/- %5.3f\n",
            Nlsb_0b, Nlsb_0b_ldp, Rdata_lsb_0b, Rdata_lsb_0b_err ) ;

       printf(" R data :  lsb_0b / lsb_0b_ldp  : ( %8d              ) / ( %8d               )    =    %5.3f +/- %5.3f (HT-only trigger)\n",
            Nhtonlytrig_lsb_0b, Nhtonlytrig_lsb_0b_ldp, RdataHTOT_lsb_0b, RdataHTOT_lsb_0b_err ) ;



       printf("\n\n\n") ;


       printf(" R ttwj :  sig    /    sb    :  ( %5.1f +/- %4.1f ) / ( %5.1f +/- %4.1f )   =   %5.3f +/- %5.3f\n",
                 (Nttbarmc_sig+lsf_WJmc*NWJmc_sig), lsf_WJmc*sqrt(NWJmc_sig),
                 (Nttbarmc_sb+lsf_WJmc*NWJmc_sb), lsf_WJmc*sqrt(NWJmc_sb),
                 Rttwj, Rttwj_err ) ;

       printf(" R ttwj :  sig_sl /    sb_sl :  ( %5.1f +/- %4.1f ) / ( %5.1f +/- %4.1f )   =   %5.3f +/- %5.3f\n",
                 (Nttbarmc_sig_sl+lsf_WJmc*NWJmc_sig_sl), lsf_WJmc*sqrt(NWJmc_sig_sl),
                 (Nttbarmc_sb_sl+lsf_WJmc*NWJmc_sb_sl), lsf_WJmc*sqrt(NWJmc_sb_sl),
                 Rttwj_sl, Rttwj_sl_err ) ;

       printf(" R data :  sig_sl /    sb_sl :    %3d              /   %3d                =   %5.3f +/- %5.3f\n",
                 Nsig_sl, Nsb_sl, Rdata_sl, Rdata_sl_err ) ;



       printf("\n\n\n") ;


       printf(" ----- Simple MC closure tests\n\n") ;

       printf("  QCD :  sb =  sb_ldp * (lsb_0b/lsb_0b_ldp) : ( %6.1f +/- %4.1f ) * ( %8.1f / %10.1f ) = %8.1f +/- %6.1f\n",
             Nqcdmc_sb_ldp, Nqcdmc_sb_ldp_err, Nqcdmc_lsb_0b, Nqcdmc_lsb_0b_ldp, comp_mc_qcd_sb, comp_mc_qcd_sb_err ) ;
       printf("                                                                               MC truth value : %8.1f +/- %6.1f\n",
                   Nqcdmc_sb, Nqcdmc_sb_err ) ;

       printf("\n") ;
       printf("  QCD : sig = sig_ldp * (lsb_0b/lsb_0b_ldp) : ( %6.1f +/- %4.1f ) * ( %8.1f / %10.1f ) = %8.1f +/- %6.1f\n",
             Nqcdmc_sig_ldp, Nqcdmc_sig_ldp_err, Nqcdmc_lsb_0b, Nqcdmc_lsb_0b_ldp, comp_mc_qcd_sig, comp_mc_qcd_sig_err ) ;
       printf("                                                                               MC truth value : %8.1f +/- %6.1f\n",
                   Nqcdmc_sig, Nqcdmc_sig_err ) ;


       printf("\n") ;
       printf("  ttwj : sig = ttwj_sb * ( ttwj_sig_sl / ttwj_sb_sl ) :  %6.1f * ( %6.1f / %6.1f ) = %6.1f +/- %4.1f\n",
              (Nttbarmc_sb + lsf_WJmc*NWJmc_sb), (Nttbarmc_sig_sl + lsf_WJmc*NWJmc_sig_sl), (Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl),
              comp_mc_ttwj_sig, comp_mc_ttwj_sig_err ) ;
       printf("                                                                       MC truth value : %6.1f +/- %4.1f\n",
                (Nttbarmc_sig + lsf_WJmc*NWJmc_sig), lsf_WJmc*sqrt(NWJmc_sig) ) ;


       printf("\n\n\n") ;

       printf(" ----- Simple Data calculations\n\n") ;


       printf(" QCD :  sb  = (  sb_ldp - (ttwj+znn+other)_ldp_mc ) * (lsb_0b/lsb_0b_ldp) : ( %4d - %5.1f ) * ( %5d / %7d ) = %6.1f +/- %4.1f\n",
           Nsb_ldp, (Nttbarmc_sb_ldp  + lsf_WJmc*NWJmc_sb_ldp  + lsf_Znnmc*NZnnmc_sb_ldp  + lsf_Ewomc*NEwomc_sb_ldp  ),
           Nlsb_0b, Nlsb_0b_ldp, comp_data_qcd_sb, comp_data_qcd_sb_err ) ;
       printf("                                                                                                  MC truth value : %8.1f +/- %4.1f\n",
                   Nqcdmc_sb, Nqcdmc_sb_err ) ;

       printf("\n") ;
       printf(" QCD :  sig = ( sig_ldp - (ttwj+znn+other)_ldp_mc ) * (lsb_0b/lsb_0b_ldp) : ( %4d - %5.1f ) * ( %5d / %7d ) = %6.1f +/- %4.1f\n",
           Nsig_ldp, (Nttbarmc_sig_ldp  + lsf_WJmc*NWJmc_sig_ldp  + lsf_Znnmc*NZnnmc_sig_ldp  + lsf_Ewomc*NEwomc_sig_ldp  ),
           Nlsb_0b, Nlsb_0b_ldp, comp_data_qcd_sig, comp_data_qcd_sig_err ) ;
       printf("                                                                                                  MC truth value : %8.1f +/- %4.1f\n",
                   Nqcdmc_sig, Nqcdmc_sig_err ) ;


       printf("\n") ;
       printf(" ttwj : sig = ( sb - (qcd_sb + (znn+other)_sb_mc) ) * (sig_sl/sb_sl) : ( %5d - ( %5.1f + %5.1f))*(%3d/%3d) = %5.1f +/- %4.1f\n",
              Nsb, comp_data_qcd_sb, (lsf_Znnmc*NZnnmc_sb + lsf_Ewomc*NEwomc_sb), Nsig_sl, Nsb_sl, comp_data_ttwj_sig, comp_data_ttwj_sig_err ) ;
       printf("                                                                                              MC truth value : %5.1f +/- %4.1f\n",
                 (Nttbarmc_sig + lsf_WJmc*NWJmc_sig), lsf_WJmc*sqrt(NWJmc_sig) ) ;


       printf("\n\n\n") ;

       printf(" ---- Z to nunu calculations\n\n" ) ;

       printf("  Znn SIG : %6.1f +/- %5.1f\n", comp_znn_sig, comp_znn_sig_err ) ;
       printf("  Znn SB  : %6.1f +/- %5.1f\n", comp_znn_sb , comp_znn_sb_err ) ;


       printf("\n\n\n") ;



     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


       printf(" --- Defining observables.\n" ) ;


      rv_Nsig        = new RooRealVar( "Nsig"        , "Nsig"        , 0.0, 1000000. ) ;
      rv_Nsb         = new RooRealVar( "Nsb"         , "Nsb"         , 0.0, 1000000. ) ;

      rv_Nsig_sl     = new RooRealVar( "Nsig_sl"     , "Nsig_sl"     , 0.0, 1000000. ) ;
      rv_Nsb_sl      = new RooRealVar( "Nsb_sl"      , "Nsb_sl"      , 0.0, 1000000. ) ;

      rv_Nsig_ldp    = new RooRealVar( "Nsig_ldp"    , "Nsig_ldp"    , 0.0, 1000000. ) ;
      rv_Nsb_ldp     = new RooRealVar( "Nsb_ldp"     , "Nsb_ldp"     , 0.0, 1000000. ) ;

      rv_Nlsb_0b     = new RooRealVar( "Nlsb_0b"     , "Nlsb_0b"     , 0.0, 1000000. ) ;
      rv_Nlsb_0b_ldp = new RooRealVar( "Nlsb_0b_ldp" , "Nlsb_0b_ldp" , 0.0, 1000000. ) ;

      if ( znnModel == 1 ) {

         rv_Nsb_ee      = new RooRealVar( "Nsb_ee"      ,"Nsb_ee"       , 0., 100. ) ;
         rv_Nsig_ee     = new RooRealVar( "Nsig_ee"     ,"Nsig_ee"      , 0., 100. ) ;

         rv_Nsb_mm      = new RooRealVar( "Nsb_mm"      ,"Nsb_mm"       , 0., 100. ) ;
         rv_Nsig_mm     = new RooRealVar( "Nsig_mm"     ,"Nsig_mm"      , 0., 100. ) ;

      } else if ( znnModel == 2 ) {

         rv_Nsigsb_ee     = new RooRealVar( "Nsigsb_ee"     ,"Nsigsb_ee"      , 0., 100. ) ;
         rv_Nsigsb_mm     = new RooRealVar( "Nsigsb_mm"     ,"Nsigsb_mm"      , 0., 100. ) ;

      }





      if ( Nsig < 0 ) {
         printf("\n\n *** Negative value for Nsig in input file.  Will set Nsig to MC expectation, which is %d.\n\n",
             TMath::Nint( Nsm_sig ) ) ;
         Nsig = TMath::Nint( Nsm_sig ) ;
      }

      rv_Nsig        -> setVal( Nsig ) ;
      rv_Nsb         -> setVal( Nsb ) ;

      rv_Nsig_sl     -> setVal( Nsig_sl ) ;
      rv_Nsb_sl      -> setVal( Nsb_sl ) ;

      rv_Nsig_ldp    -> setVal( Nsig_ldp ) ;
      rv_Nsb_ldp     -> setVal( Nsb_ldp ) ;

      rv_Nlsb_0b     -> setVal( Nhtonlytrig_lsb_0b ) ;
      rv_Nlsb_0b_ldp -> setVal( Nhtonlytrig_lsb_0b_ldp ) ;

      if ( znnModel == 1 ) {

         rv_Nsb_ee      -> setVal( Nsb_ee ) ;
         rv_Nsig_ee     -> setVal( Nsig_ee ) ;

         rv_Nsb_mm      -> setVal( Nsb_mm ) ;
         rv_Nsig_mm     -> setVal( Nsig_mm ) ;

      } else if ( znnModel == 2 ) {

         rv_Nsigsb_ee     -> setVal( Nsig_ee + Nsb_ee ) ;
         rv_Nsigsb_mm     -> setVal( Nsig_mm + Nsb_mm ) ;

      }







    //++++++++ Parameters of the likelihood +++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining parameters.\n" ) ;



    //____ Counts in SIG ______________________

      if ( useSigTtwjVar ) {
         rrv_mu_ttwj_sig = new RooRealVar( "mu_ttwj_sig"   , "mu_ttwj_sig"   , 0.0, 10000. ) ;
         rv_mu_ttwj_sig = rrv_mu_ttwj_sig ;
         rrv_mu_ttwj_sig   -> setVal( Nttbarmc_sig + lsf_WJmc*NWJmc_sig ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sig  = new RooRealVar( "mu_qcd_sig"    , "mu_qcd_sig"    , 0.0, 50. ) ;
         rv_mu_qcd_sig = rrv_mu_qcd_sig ;
         rrv_mu_qcd_sig  -> setVal( Nqcdmc_sig ) ; //-- this is a starting value only.
      }


      //-- Note: Ewo is rfv

      rv_mu_znn_sig      = new RooRealVar( "mu_znn_sig"    , "mu_znn_sig"    , 0.0, 80. ) ;

      float maxSusySig = 4.0*Nsig ;
      rv_mu_susy_sig     = new RooRealVar( "mu_susy_sig"   , "mu_susy_sig"   , 0.0, maxSusySig ) ;


      rv_mu_znn_sig   -> setVal( comp_znn_sig ) ;  //-- this is a starting value only.
      rv_mu_susy_sig    -> setVal( 0. ) ;  //-- this is a starting value only.





    //____ Counts in SB  ______________________


      if ( !useSigTtwjVar ) {
         rrv_mu_ttwj_sb  = new RooRealVar( "mu_ttwj_sb"    , "mu_ttwj_sb"    , 0.0, 10000. ) ;
         rv_mu_ttwj_sb = rrv_mu_ttwj_sb ;
         rrv_mu_ttwj_sb   -> setVal( Nttbarmc_sb + lsf_WJmc*NWJmc_sb ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sb  = new RooRealVar( "mu_qcd_sb"    , "mu_qcd_sb"    , 0.0, 100. ) ;
         rv_mu_qcd_sb = rrv_mu_qcd_sb ;
         rrv_mu_qcd_sb  -> setVal( Nqcdmc_sb ) ; //-- this is a starting value only.
      }

      //-- Note: QCD is rfv
      //-- Note: Ewo is rfv
      //-- Note: SUSY is rfv

      if ( znnModel == 1 ) {

         rrv_mu_znn_sb       = new RooRealVar( "mu_znn_sb"     , "mu_znn_sb"     , 0.0, 150. ) ;

         rrv_mu_znn_sb   -> setVal( comp_znn_sb ) ;  //-- this is a starting value only.

         rv_mu_znn_sb = rrv_mu_znn_sb ;
      }
      //-- Note: Znn is rfv in Znn model 2.






    //____ Counts in SIG, SL  ______________________

      rv_mu_ttwj_sig_sl  = new RooRealVar( "mu_ttwj_sig_sl"    , "mu_ttwj_sig_sl"    , 0.0, 500. ) ;

      //-- Note: QCD, Ewo, and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      rv_mu_ttwj_sig_sl  -> setVal( Nttbarmc_sig_sl + lsf_WJmc*NWJmc_sig_sl ) ;  //-- this is a starting value only.







    //____ Counts in SB, SL  ______________________

      rv_mu_ttwj_sb_sl  = new RooRealVar( "mu_ttwj_sb_sl"    , "mu_ttwj_sb_sl"    , 0.0, 3000. ) ;

      //-- Note: QCD, Ewo, and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      rv_mu_ttwj_sb_sl  -> setVal( Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl ) ;  //-- this is a starting value only.







    //____ Counts in SIG, LDP  ______________________

      if ( useLdpVars ) {
         rrv_mu_qcd_sig_ldp  = new RooRealVar( "mu_qcd_sig_ldp"    , "mu_qcd_sig_ldp"    , 0.0, 500. ) ;
         rv_mu_qcd_sig_ldp = rrv_mu_qcd_sig_ldp ;
         rrv_mu_qcd_sig_ldp  -> setVal( Nqcdmc_sig_ldp ) ; //-- this is a starting value only.
      }

      //-- Note: Ewo is assumed to be negligible and is not explicitly included.
      //-- Note: Znn is rfv (MC)
      //-- Note: ttwj is rfv (MC)
      //-- Note: SUSY is rfv








    //____ Counts in SB, LDP  ______________________

      if ( useLdpVars ) {
         rrv_mu_qcd_sb_ldp  = new RooRealVar( "mu_qcd_sb_ldp"    , "mu_qcd_sb_ldp"    , 0.0, 1000. ) ;
         rv_mu_qcd_sb_ldp = rrv_mu_qcd_sb_ldp ;
         rrv_mu_qcd_sb_ldp  -> setVal( Nqcdmc_sb_ldp ) ; //-- this is a starting value only.
      }


      //-- Note: Ewo is assumed to be negligible and is not explicitly included.
      //-- Note: Znn is rfv (MC)
      //-- Note: ttwj is rfv (MC)
      //-- Note: SUSY is rfv








    //____ Counts in LSB, 0b  ______________________

      rv_mu_qcd_lsb_0b  = new RooRealVar( "mu_qcd_lsb_0b"    ,  "mu_qcd_lsb_0b" ,  0.    ,  10000. ) ;

      //-- Note: The 0btag LSB is assumed to be 100% QCD.

      rv_mu_qcd_lsb_0b  -> setVal( Nhtonlytrig_lsb_0b ) ;  //-- this is a starting value only.






    //____ Counts in LSB, 0b, LDP  ______________________

      rv_mu_qcd_lsb_0b_ldp  = new RooRealVar( "mu_qcd_lsb_0b_ldp"    ,  "mu_qcd_lsb_0b_ldp",   0. ,  10000. ) ;

      //-- Note: The 0btag LSB, LDP is assumed to be 100% QCD.

      rv_mu_qcd_lsb_0b_ldp  -> setVal( Nhtonlytrig_lsb_0b_ldp ) ;  //-- this is a starting value only.






    //____ Counts in SB, ee  ______________________


      //-- Note: The Z to ee sample is assumed to be 100% Z to ee.
      //-- Note: zee is rfv






    //____ Counts in SIG, ee  ______________________


      //-- Note: The Z to ee sample is assumed to be 100% Z to ee.
      //-- Note: zee is rfv








    //____ Counts in SB, mm  ______________________


      //-- Note: The Z to mm sample is assumed to be 100% Z to mm.
      //-- Note: zmm is rfv






    //____ Counts in SIG, mm  ______________________


      //-- Note: The Z to mm sample is assumed to be 100% Z to mm.
      //-- Note: zmm is rfv






    //____ MC inputs _______________________________

      printf(" --- Defining MC parameters.\n" ) ;

     //-- SUSY

      rv_mu_susymc_sig      = new RooRealVar( "mu_susymc_sig"     , "mu_susymc_sig"     , 0.0, 100000. ) ;
      rv_mu_susymc_sb       = new RooRealVar( "mu_susymc_sb"      , "mu_susymc_sb"      , 0.0, 100000. ) ;
      rv_mu_susymc_sig_sl   = new RooRealVar( "mu_susymc_sig_sl"  , "mu_susymc_sig_sl"  , 0.0, 100000. ) ;
      rv_mu_susymc_sb_sl    = new RooRealVar( "mu_susymc_sb_sl"   , "mu_susymc_sb_sl"   , 0.0, 100000. ) ;
      rv_mu_susymc_sig_ldp  = new RooRealVar( "mu_susymc_sig_ldp" , "mu_susymc_sig_ldp" , 0.0, 100000. ) ;
      rv_mu_susymc_sb_ldp   = new RooRealVar( "mu_susymc_sb_ldp"  , "mu_susymc_sb_ldp"  , 0.0, 100000. ) ;

      rv_mu_susymc_sig     -> setVal( 0.1 ) ;
      rv_mu_susymc_sb      -> setVal( 0. ) ;
      rv_mu_susymc_sig_sl  -> setVal( 0. ) ;
      rv_mu_susymc_sb_sl   -> setVal( 0. ) ;
      rv_mu_susymc_sig_ldp -> setVal( 0. ) ;
      rv_mu_susymc_sb_ldp  -> setVal( 0. ) ;

      rv_mu_susymc_sig     -> setConstant(kTRUE) ;
      rv_mu_susymc_sb      -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_sl  -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_sl   -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_ldp -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_ldp  -> setConstant(kTRUE) ;



     //-- SIG, LDP

      rv_mu_ttbarmc_sig_ldp   = new RooRealVar( "mu_ttbarmc_sig_ldp" ,"mu_ttbarmc_sig_ldp" , 0., 1000. ) ;
      rv_mu_WJmc_sig_ldp      = new RooRealVar( "mu_WJmc_sig_ldp"    ,"mu_WJmc_sig_ldp"    , 0., 1000. ) ;
      rv_mu_Znnmc_sig_ldp     = new RooRealVar( "mu_Znnmc_sig_ldp"   ,"mu_Znnmc_sig_ldp"   , 0., 1000. ) ;
      rv_mu_Ewomc_sig_ldp     = new RooRealVar( "mu_Ewomc_sig_ldp"   ,"mu_Ewomc_sig_ldp"   , 0., 1000. ) ;

      rv_mu_ttbarmc_sig_ldp  -> setVal( Nttbarmc_sig_ldp ) ;
      rv_mu_WJmc_sig_ldp     -> setVal( NWJmc_sig_ldp ) ;
      rv_mu_Znnmc_sig_ldp    -> setVal( NZnnmc_sig_ldp ) ;
      rv_mu_Ewomc_sig_ldp    -> setVal( NEwomc_sig_ldp ) ;

      rv_mu_ttbarmc_sig_ldp  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sig_ldp     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sig_ldp    -> setConstant( kTRUE ) ;
      rv_mu_Ewomc_sig_ldp    -> setConstant( kTRUE ) ;


     //-- SB, LDP

      rv_mu_ttbarmc_sb_ldp   = new RooRealVar( "mu_ttbarmc_sb_ldp" ,"mu_ttbarmc_sb_ldp" , 0., 1000. ) ;
      rv_mu_WJmc_sb_ldp      = new RooRealVar( "mu_WJmc_sb_ldp"    ,"mu_WJmc_sb_ldp"    , 0., 1000. ) ;
      rv_mu_Znnmc_sb_ldp     = new RooRealVar( "mu_Znnmc_sb_ldp"   ,"mu_Znnmc_sb_ldp"   , 0., 1000. ) ;
      rv_mu_Ewomc_sb_ldp     = new RooRealVar( "mu_Ewomc_sb_ldp"   ,"mu_Ewomc_sb_ldp"   , 0., 1000. ) ;

      rv_mu_ttbarmc_sb_ldp  -> setVal( Nttbarmc_sb_ldp ) ;
      rv_mu_WJmc_sb_ldp     -> setVal( NWJmc_sb_ldp ) ;
      rv_mu_Znnmc_sb_ldp    -> setVal( NZnnmc_sb_ldp ) ;
      rv_mu_Ewomc_sb_ldp    -> setVal( NEwomc_sb_ldp ) ;

      rv_mu_ttbarmc_sb_ldp  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sb_ldp     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sb_ldp    -> setConstant( kTRUE ) ;
      rv_mu_Ewomc_sb_ldp    -> setConstant( kTRUE ) ;


     //-- SIG

      rv_mu_Ewomc_sig     = new RooRealVar( "mu_Ewomc_sig"   ,"mu_Ewomc_sig"   , 0., 1000. ) ;
      rv_mu_Ewomc_sig    -> setVal( NEwomc_sig ) ;
      rv_mu_Ewomc_sig    -> setConstant( kTRUE ) ;

     //-- SB

      rv_mu_Ewomc_sb     = new RooRealVar( "mu_Ewomc_sb"   ,"mu_Ewomc_sb"   , 0., 1000. ) ;
      rv_mu_Ewomc_sb    -> setVal( NEwomc_sb ) ;
      rv_mu_Ewomc_sb    -> setConstant( kTRUE ) ;





    //+++++++ Gaussian constraints ++++++++++++++++++++++++++++++++

      printf(" --- Defining Gaussian constraint and constant parameters.\n" ) ;

    //_______ Efficiency scale factor.  Applied to SUSY and all MC inputs _______________

    //   August 23, 2011:  switching to correlated log-normal PDFs
    //
    //
      double pmin, pmax ;


    //--- mean parameters.
      rv_mean_eff_sf_sig     = new RooRealVar( "mean_eff_sf_sig"    , "mean_eff_sf_sig", 0., 10. ) ;
      rv_mean_eff_sf_sb      = new RooRealVar( "mean_eff_sf_sb"     , "mean_eff_sf_sb", 0., 10. ) ;
      rv_mean_eff_sf_sig_sl  = new RooRealVar( "mean_eff_sf_sig_sl" , "mean_eff_sf_sig_sl", 0., 10. ) ;
      rv_mean_eff_sf_sb_sl   = new RooRealVar( "mean_eff_sf_sb_sl"  , "mean_eff_sf_sb_sl", 0., 10. ) ;
      rv_mean_eff_sf_sig_ldp = new RooRealVar( "mean_eff_sf_sig_ldp", "mean_eff_sf_sig_ldp", 0., 10. ) ;
      rv_mean_eff_sf_sb_ldp  = new RooRealVar( "mean_eff_sf_sb_ldp" , "mean_eff_sf_sb_ldp", 0., 10. ) ;

    //--- width parameters.
      rv_width_eff_sf_sig     = new RooRealVar( "width_eff_sf_sig"    , "width_eff_sf_sig", 0., 10. ) ;
      rv_width_eff_sf_sb      = new RooRealVar( "width_eff_sf_sb"     , "width_eff_sf_sb", 0., 10. ) ;
      rv_width_eff_sf_sig_sl  = new RooRealVar( "width_eff_sf_sig_sl" , "width_eff_sf_sig_sl", 0., 10. ) ;
      rv_width_eff_sf_sb_sl   = new RooRealVar( "width_eff_sf_sb_sl"  , "width_eff_sf_sb_sl", 0., 10. ) ;
      rv_width_eff_sf_sig_ldp = new RooRealVar( "width_eff_sf_sig_ldp", "width_eff_sf_sig_ldp", 0., 10. ) ;
      rv_width_eff_sf_sb_ldp  = new RooRealVar( "width_eff_sf_sb_ldp" , "width_eff_sf_sb_ldp", 0., 10. ) ;


      rv_mean_eff_sf_sig     -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb      -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_sl  -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_sl   -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_ldp -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_ldp  -> setVal( 1.00 ) ;

      rv_mean_eff_sf_sig     -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb      -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_sl  -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_sl   -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_ldp -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_ldp  -> setConstant( kTRUE ) ;


      //--- Initialize all width parameters to 15%.  They will be reset in the susy scan.
      rv_width_eff_sf_sig     -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb      -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_sl  -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_sl   -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_ldp -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_ldp  -> setVal( 0.15 ) ;

      rv_width_eff_sf_sig     -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb      -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_sl  -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_sl   -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_ldp -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_ldp  -> setConstant( kTRUE ) ;


    //--- Underlying Gaussian variable for log-normal.

      rv_eff_sf_prim = new RooRealVar("eff_sf_prim", "eff_sf_prim", 0., -5., 5. ) ;
      rv_eff_sf_nom  = new RooRealVar("eff_sf_nom" , "eff_sf_nom" , 0., -5., 5. ) ;
      rv_eff_sf_nom->setConstant() ;


    //--- Systematics

      pmin = (sf_mc-4*sf_mc_err) ;
      pmax = (sf_mc+4*sf_mc_err) ;
      if ( pmin < 0 ) pmin = 0.0 ;
      rv_sf_mc = new RooRealVar( "sf_mc", "sf_mc", pmin, pmax ) ;
      rv_sf_mc -> setVal( sf_mc ) ;


      pmin = (sf_qcd_sb-4*sf_qcd_sb_err) ;
      pmax = (sf_qcd_sb+4*sf_qcd_sb_err) ;
      if ( pmin < 0 ) pmin = 0.0 ;
      rv_sf_qcd_sb = new RooRealVar( "sf_qcd_sb", "sf_qcd_sb", pmin, pmax ) ;
      rv_sf_qcd_sb -> setVal( sf_qcd_sb ) ;


      pmin = (sf_qcd_sig-4*sf_qcd_sig_err) ;
      pmax = (sf_qcd_sig+4*sf_qcd_sig_err) ;
      if ( pmin < 0 ) pmin = 0.0 ;
      rv_sf_qcd_sig = new RooRealVar( "sf_qcd_sig", "sf_qcd_sig", pmin, pmax ) ;
      rv_sf_qcd_sig -> setVal( sf_qcd_sig ) ;


      pmin = (sf_ttwj_sig-4*sf_ttwj_sig_err) ;
      pmax = (sf_ttwj_sig+4*sf_ttwj_sig_err) ;
      if ( pmin < 0 ) pmin = 0.0 ;
      rv_sf_ttwj_sig = new RooRealVar( "sf_ttwj_sig", "sf_ttwj_sig", pmin, pmax ) ;
      rv_sf_ttwj_sig -> setVal( sf_ttwj_sig ) ;



    //-- Z to nunu stuff

      pmin = (acc_ee_mean-4.*acc_ee_err) ;
      pmax = (acc_ee_mean+4.*acc_ee_err) ;
      if ( pmin < 0 ) pmin = 0.1 ;
      rv_acc_ee  = new RooRealVar( "acc_ee", "acc_ee", pmin, pmax ) ;
      rv_acc_ee  -> setVal( acc_ee_mean ) ;

      pmin = (acc_mm_mean-4.*acc_mm_err) ;
      pmax = (acc_mm_mean+4.*acc_mm_err) ;
      if ( pmin < 0 ) pmin = 0.1 ;
      rv_acc_mm  = new RooRealVar( "acc_mm", "acc_mm", pmin, pmax ) ;
      rv_acc_mm  -> setVal( acc_mm_mean ) ;


      pmin = (eff_ee_mean-4.*eff_ee_err) ;
      pmax = (eff_ee_mean+4.*eff_ee_err) ;
      if ( pmin < 0 ) pmin = 0.1 ;
      rv_eff_ee  = new RooRealVar( "eff_ee", "eff_ee", pmin, pmax ) ;
      rv_eff_ee  -> setVal( eff_ee_mean ) ;

      pmin = (eff_mm_mean-4.*eff_mm_err) ;
      pmax = (eff_mm_mean+4.*eff_mm_err) ;
      if ( pmin < 0 ) pmin = 0.1 ;
      rv_eff_mm  = new RooRealVar( "eff_mm", "eff_mm", pmin, pmax ) ;
      rv_eff_mm  -> setVal( eff_mm_mean ) ;


      rv_znnoverll_bfratio = new RooRealVar( "znnoverll_bfratio", "znnoverll_bfratio", 0.1, 10. ) ;
      rv_znnoverll_bfratio -> setVal( 5.95 ) ;
      rv_znnoverll_bfratio -> setConstant( kTRUE ) ;

      rv_dataoverll_lumiratio = new RooRealVar( "dataoverll_lumiratio", "dataoverll_lumiratio", 0.1, 10.0 ) ;
      rv_dataoverll_lumiratio  -> setVal( DataLumi / Ztoll_lumi ) ;
      rv_dataoverll_lumiratio  -> setConstant( kTRUE ) ;

      if ( znnModel == 2 ) {
         rv_knn_sig = new RooRealVar( "knn_sig" , "knn_sig" , 0.01, 1.0 ) ;
         rv_knn_sb  = new RooRealVar( "knn_sb"  , "knn_sb"  , 0.01, 1.0 ) ;
         rv_knn_sig  -> setVal( knn_sig_mean ) ;
         rv_knn_sb   -> setVal( knn_sb_mean  ) ;
      }


      pmin = (fsig_ee_mean-4.*fsig_ee_err) ;
      pmax = (fsig_ee_mean+4.*fsig_ee_err) ;
      if ( pmin < 0 ) pmin = 0.1 ;
      rv_fsig_ee  = new RooRealVar( "fsig_ee", "fsig_ee", pmin, 1.0 ) ;
      rv_fsig_ee  -> setVal( fsig_ee_mean ) ;

      pmin = (fsig_mm_mean-4.*fsig_mm_err) ;
      pmax = (fsig_mm_mean+4.*fsig_mm_err) ;
      if ( pmin < 0 ) pmin = 0.1 ;
      rv_fsig_mm  = new RooRealVar( "fsig_mm", "fsig_mm", pmin, 1.0 ) ;
      rv_fsig_mm  -> setVal( fsig_mm_mean ) ;




      pmin = (sf_ee-4*sf_ee_err) ;
      pmax = (sf_ee+4*sf_ee_err) ;
      if ( pmin < 0 ) pmin = 0.0 ;
      rv_sf_ee = new RooRealVar( "sf_ee", "sf_ee", pmin, pmax ) ;
      rv_sf_ee -> setVal( sf_ee ) ;


      pmin = (sf_mm-4*sf_mm_err) ;
      pmax = (sf_mm+4*sf_mm_err) ;
      if ( pmin < 0 ) pmin = 0.0 ;
      rv_sf_mm = new RooRealVar( "sf_mm", "sf_mm", pmin, pmax ) ;
      rv_sf_mm -> setVal( sf_mm ) ;




     //+++++++++++++++++ Relationships between parameters ++++++++++++++++++++++++++++++++++++++++++++

       printf(" --- Defining relationships between parameters.\n" ) ;




    //-- ttwj

      if ( useSigTtwjVar ) {
         rfv_mu_ttwj_sb = new RooFormulaVar( "mu_ttwj_sb",
                                                       "mu_ttwj_sig * (1.0/sf_ttwj_sig) * (mu_ttwj_sb_sl/mu_ttwj_sig_sl)",
                                                       RooArgSet( *rv_mu_ttwj_sig, *rv_sf_ttwj_sig, *rv_mu_ttwj_sb_sl, *rv_mu_ttwj_sig_sl ) ) ;
         rv_mu_ttwj_sb = rfv_mu_ttwj_sb ;
      } else {
         rfv_mu_ttwj_sig = new RooFormulaVar( "mu_ttwj_sig",
                                                       "mu_ttwj_sb * sf_ttwj_sig * (mu_ttwj_sig_sl/mu_ttwj_sb_sl)",
                                                       RooArgSet( *rv_mu_ttwj_sb, *rv_sf_ttwj_sig, *rv_mu_ttwj_sig_sl, *rv_mu_ttwj_sb_sl ) ) ;
         rv_mu_ttwj_sig = rfv_mu_ttwj_sig ;
      }




      rv_mu_ttwj_sig_ldp = new RooFormulaVar( "mu_ttwj_sig_ldp",
                              "mu_ttbarmc_sig_ldp + mu_WJmc_sig_ldp",
                              RooArgSet( *rv_mu_ttbarmc_sig_ldp,
                                         *rv_mu_WJmc_sig_ldp ) ) ;


      rv_mu_ttwj_sb_ldp = new RooFormulaVar( "mu_ttwj_sb_ldp",
                              "mu_ttbarmc_sb_ldp + mu_WJmc_sb_ldp",
                              RooArgSet( *rv_mu_ttbarmc_sb_ldp,
                                         *rv_mu_WJmc_sb_ldp ) ) ;






    //-- QCD

      if ( useLdpVars ) {

         rfv_mu_qcd_sig = new RooFormulaVar( "mu_qcd_sig",
                                     "mu_qcd_sig_ldp * sf_qcd_sig * ( mu_qcd_lsb_0b / mu_qcd_lsb_0b_ldp )",
                                     RooArgSet( *rv_mu_qcd_sig_ldp, *rv_sf_qcd_sig, *rv_mu_qcd_lsb_0b, *rv_mu_qcd_lsb_0b_ldp ) ) ;
         rv_mu_qcd_sig = rfv_mu_qcd_sig ;

         rfv_mu_qcd_sb = new RooFormulaVar( "mu_qcd_sb",
                                     "mu_qcd_sb_ldp * sf_qcd_sb * ( mu_qcd_lsb_0b / mu_qcd_lsb_0b_ldp )",
                                     RooArgSet( *rv_mu_qcd_sb_ldp, *rv_sf_qcd_sb, *rv_mu_qcd_lsb_0b, *rv_mu_qcd_lsb_0b_ldp ) ) ;
         rv_mu_qcd_sb = rfv_mu_qcd_sb ;

      } else {

         rfv_mu_qcd_sig_ldp = new RooFormulaVar( "mu_qcd_sig_ldp",
                                     "mu_qcd_sig * (1.0/sf_qcd_sig) * ( mu_qcd_lsb_0b_ldp / mu_qcd_lsb_0b )",
                                     RooArgSet( *rv_mu_qcd_sig, *rv_sf_qcd_sig, *rv_mu_qcd_lsb_0b_ldp, *rv_mu_qcd_lsb_0b ) ) ;
         rv_mu_qcd_sig_ldp = rfv_mu_qcd_sig_ldp ;

         rfv_mu_qcd_sb_ldp = new RooFormulaVar( "mu_qcd_sb_ldp",
                                     "mu_qcd_sb * (1.0/sf_qcd_sb) * ( mu_qcd_lsb_0b_ldp / mu_qcd_lsb_0b )",
                                     RooArgSet( *rv_mu_qcd_sb, *rv_sf_qcd_sb, *rv_mu_qcd_lsb_0b_ldp, *rv_mu_qcd_lsb_0b ) ) ;
         rv_mu_qcd_sb_ldp = rfv_mu_qcd_sb_ldp ;

      }




    //-- SUSY

      rv_mu_susy_sb = new RooFormulaVar( "mu_susy_sb",
                                        "mu_susymc_sb * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sb, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sig_sl = new RooFormulaVar( "mu_susy_sig_sl",
                                        "mu_susymc_sig_sl * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sig_sl, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sb_sl = new RooFormulaVar( "mu_susy_sb_sl",
                                        "mu_susymc_sb_sl * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sb_sl, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sig_ldp = new RooFormulaVar( "mu_susy_sig_ldp",
                                        "mu_susymc_sig_ldp * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sig_ldp, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sb_ldp = new RooFormulaVar( "mu_susy_sb_ldp",
                                        "mu_susymc_sb_ldp * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sb_ldp, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;





    //-- Z to nu nu

      rv_mu_znn_sig_ldp = new RooFormulaVar( "mu_znn_sig_ldp",
                              "mu_Znnmc_sig_ldp",
                              RooArgSet( *rv_mu_Znnmc_sig_ldp ) ) ;

      rv_mu_znn_sb_ldp = new RooFormulaVar( "mu_znn_sb_ldp",
                              "mu_Znnmc_sb_ldp",
                              RooArgSet( *rv_mu_Znnmc_sb_ldp ) ) ;


      if ( znnModel == 1 ) {

         rv_mu_zee_sb_ee = new RooFormulaVar( "mu_zee_sb_ee",
                                         "mu_znn_sb * sf_ee * ( acc_ee * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
                                         RooArgSet( *rv_mu_znn_sb, *rv_sf_ee, *rv_acc_ee, *rv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

         rv_mu_zee_sig_ee = new RooFormulaVar( "mu_zee_sig_ee",
                                         "mu_znn_sig * sf_ee * ( acc_ee * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
                                         RooArgSet( *rv_mu_znn_sig, *rv_sf_ee, *rv_acc_ee, *rv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

         rv_mu_zmm_sb_mm = new RooFormulaVar( "mu_zmm_sb_mm",
                                         "mu_znn_sb * sf_mm * ( acc_mm * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
                                         RooArgSet( *rv_mu_znn_sb, *rv_sf_mm, *rv_acc_mm, *rv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

         rv_mu_zmm_sig_mm = new RooFormulaVar( "mu_zmm_sig_mm",
                                         "mu_znn_sig * sf_mm * ( acc_mm * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
                                         RooArgSet( *rv_mu_znn_sig, *rv_sf_mm, *rv_acc_mm, *rv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

      } else if ( znnModel == 2 ) {

         rfv_mu_znn_sb = new RooFormulaVar( "mu_znn_sb",
                                          "mu_znn_sig * ( knn_sb / knn_sig )",
                                          RooArgSet( *rv_mu_znn_sig, *rv_knn_sb, *rv_knn_sig ) ) ;

         rv_mu_znn_sb = rfv_mu_znn_sb ;

         rv_mu_zee_sigsb_ee = new RooFormulaVar( "mu_zee_sigsb_ee",
                                       "( mu_znn_sig / knn_sig ) * sf_ee * ( (acc_ee * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                         RooArgSet( *rv_mu_znn_sig, *rv_knn_sig, *rv_sf_ee, *rv_acc_ee, *rv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

         rv_mu_zmm_sigsb_mm = new RooFormulaVar( "mu_zmm_sigsb_mm",
                                       "( mu_znn_sig / knn_sig ) * sf_mm * ( (acc_mm * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                         RooArgSet( *rv_mu_znn_sig, *rv_knn_sig, *rv_sf_mm, *rv_acc_mm, *rv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

      }

    //-- EWO

      rv_mu_ewo_sig     = new RooFormulaVar( "mu_ewo_sig"     , "mu_Ewomc_sig"     , RooArgSet( *rv_mu_Ewomc_sig     ) ) ;
      rv_mu_ewo_sb      = new RooFormulaVar( "mu_ewo_sb"      , "mu_Ewomc_sb"      , RooArgSet( *rv_mu_Ewomc_sb      ) ) ;
      rv_mu_ewo_sig_ldp = new RooFormulaVar( "mu_ewo_sig_ldp" , "mu_Ewomc_sig_ldp" , RooArgSet( *rv_mu_Ewomc_sig_ldp ) ) ;
      rv_mu_ewo_sb_ldp  = new RooFormulaVar( "mu_ewo_sb_ldp"  , "mu_Ewomc_sb_ldp"  , RooArgSet( *rv_mu_Ewomc_sb_ldp  ) ) ;




    //-- Parametric relations between correlated log-normal efficiency scale factors.


      rv_eff_sf_sig = new RooFormulaVar( "eff_sf_sig",
                                         "mean_eff_sf_sig * pow( exp( width_eff_sf_sig/mean_eff_sf_sig ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig, *rv_width_eff_sf_sig, *rv_mean_eff_sf_sig, *rv_eff_sf_prim ) ) ;

      rv_eff_sf_sb = new RooFormulaVar( "eff_sf_sb",
                                         "mean_eff_sf_sb * pow( exp( width_eff_sf_sb/mean_eff_sf_sb ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb, *rv_width_eff_sf_sb, *rv_mean_eff_sf_sb, *rv_eff_sf_prim ) ) ;

      rv_eff_sf_sig_sl = new RooFormulaVar( "eff_sf_sig_sl",
                                         "mean_eff_sf_sig_sl * pow( exp( width_eff_sf_sig_sl/mean_eff_sf_sig_sl ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_sl, *rv_width_eff_sf_sig_sl, *rv_mean_eff_sf_sig_sl, *rv_eff_sf_prim ) ) ;

      rv_eff_sf_sb_sl = new RooFormulaVar( "eff_sf_sb_sl",
                                         "mean_eff_sf_sb_sl * pow( exp( width_eff_sf_sb_sl/mean_eff_sf_sb_sl ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_sl, *rv_width_eff_sf_sb_sl, *rv_mean_eff_sf_sb_sl, *rv_eff_sf_prim ) ) ;

      rv_eff_sf_sig_ldp = new RooFormulaVar( "eff_sf_sig_ldp",
                                         "mean_eff_sf_sig_ldp * pow( exp( width_eff_sf_sig_ldp/mean_eff_sf_sig_ldp ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_ldp, *rv_width_eff_sf_sig_ldp, *rv_mean_eff_sf_sig_ldp, *rv_eff_sf_prim ) ) ;

      rv_eff_sf_sb_ldp = new RooFormulaVar( "eff_sf_sb_ldp",
                                         "mean_eff_sf_sb_ldp * pow( exp( width_eff_sf_sb_ldp/mean_eff_sf_sb_ldp ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_ldp, *rv_width_eff_sf_sb_ldp, *rv_mean_eff_sf_sb_ldp, *rv_eff_sf_prim ) ) ;



    //+++++++++++++ Expected counts for observables in terms of parameters ++++++++++++++++++

       printf(" --- Defining expected counts in terms of parameters.\n" ) ;


      rv_n_sig         = new RooFormulaVar( "n_sig",
                                     "mu_ttwj_sig + mu_qcd_sig + mu_znn_sig + eff_sf_sig*( mu_ewo_sig + mu_susy_sig)",
                                     RooArgSet( *rv_mu_ttwj_sig, *rv_mu_qcd_sig, *rv_mu_znn_sig, *rv_eff_sf_sig, *rv_mu_ewo_sig, *rv_mu_susy_sig ) ) ;

      rv_n_sb          = new RooFormulaVar( "n_sb",
                                     "mu_ttwj_sb  + mu_qcd_sb  + mu_znn_sb  + eff_sf_sb*( mu_ewo_sb  + mu_susy_sb )",
                                     RooArgSet( *rv_mu_ttwj_sb , *rv_mu_qcd_sb , *rv_mu_znn_sb , *rv_eff_sf_sb, *rv_mu_ewo_sb , *rv_mu_susy_sb  ) ) ;

      rv_n_sig_ldp     = new RooFormulaVar( "n_sig_ldp",
                                     "mu_qcd_sig_ldp + eff_sf_sig_ldp*( sf_mc * (mu_ttwj_sig_ldp + mu_znn_sig_ldp + mu_ewo_sig_ldp) + mu_susy_sig_ldp)",
                                     RooArgSet( *rv_mu_qcd_sig_ldp, *rv_eff_sf_sig_ldp, *rv_sf_mc, *rv_mu_ttwj_sig_ldp, *rv_mu_znn_sig_ldp, *rv_mu_ewo_sig_ldp, *rv_mu_susy_sig_ldp ) ) ;

      rv_n_sb_ldp      = new RooFormulaVar( "n_sb_ldp",
                                     "mu_qcd_sb_ldp + eff_sf_sb_ldp*( sf_mc * (mu_ttwj_sb_ldp + mu_znn_sb_ldp + mu_ewo_sb_ldp) + mu_susy_sb_ldp)",
                                     RooArgSet( *rv_mu_qcd_sb_ldp, *rv_eff_sf_sb_ldp, *rv_sf_mc, *rv_mu_ttwj_sb_ldp, *rv_mu_znn_sb_ldp, *rv_mu_ewo_sb_ldp, *rv_mu_susy_sb_ldp ) ) ;

      rv_n_sig_sl      = new RooFormulaVar( "n_sig_sl",
                                     "mu_ttwj_sig_sl + eff_sf_sig_sl*mu_susy_sig_sl",
                                     RooArgSet( *rv_mu_ttwj_sig_sl, *rv_eff_sf_sig_sl, *rv_mu_susy_sig_sl ) ) ;

      rv_n_sb_sl       = new RooFormulaVar( "n_sb_sl",
                                     "mu_ttwj_sb_sl + eff_sf_sb_sl*mu_susy_sb_sl",
                                     RooArgSet( *rv_mu_ttwj_sb_sl, *rv_eff_sf_sb_sl, *rv_mu_susy_sb_sl ) ) ;

      rv_n_lsb_0b      = new RooFormulaVar( "n_lsb_0b",
                                     "mu_qcd_lsb_0b",
                                     RooArgSet( *rv_mu_qcd_lsb_0b ) ) ;

      rv_n_lsb_0b_ldp  = new RooFormulaVar( "n_lsb_0b_ldp",
                                     "mu_qcd_lsb_0b_ldp",
                                     RooArgSet( *rv_mu_qcd_lsb_0b_ldp ) ) ;

      if ( znnModel == 1 ) {

         rv_n_sig_ee      = new RooFormulaVar( "n_sig_ee",
                                        "mu_zee_sig_ee / fsig_ee",
                                        RooArgSet( *rv_mu_zee_sig_ee, *rv_fsig_ee ) ) ;

         rv_n_sb_ee       = new RooFormulaVar( "n_sb_ee",
                                        "mu_zee_sb_ee / fsig_ee",
                                        RooArgSet( *rv_mu_zee_sb_ee, *rv_fsig_ee ) ) ;

         rv_n_sig_mm      = new RooFormulaVar( "n_sig_mm",
                                        "mu_zmm_sig_mm / fsig_mm",
                                        RooArgSet( *rv_mu_zmm_sig_mm, *rv_fsig_mm ) ) ;

         rv_n_sb_mm       = new RooFormulaVar( "n_sb_mm",
                                        "mu_zmm_sb_mm / fsig_mm",
                                        RooArgSet( *rv_mu_zmm_sb_mm, *rv_fsig_mm ) ) ;

      } else if ( znnModel == 2 ) {

         rv_n_sigsb_ee      = new RooFormulaVar( "n_sigsb_ee",
                                        "mu_zee_sigsb_ee / fsig_ee",
                                        RooArgSet( *rv_mu_zee_sigsb_ee, *rv_fsig_ee ) ) ;

         rv_n_sigsb_mm      = new RooFormulaVar( "n_sigsb_mm",
                                        "mu_zmm_sigsb_mm / fsig_mm",
                                        RooArgSet( *rv_mu_zmm_sigsb_mm, *rv_fsig_mm ) ) ;

      }

   //++++++++++++ PDFs for the likelihood +++++++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining PDFs of the likelihood.\n" ) ;

      pdf_Nsig        = new RooPoisson( "pdf_Nsig"        , "Nsig Poisson PDF"        , *rv_Nsig        , *rv_n_sig ) ;
      pdf_Nsb         = new RooPoisson( "pdf_Nsb"         , "Nsb Poisson PDF"         , *rv_Nsb         , *rv_n_sb ) ;
      pdf_Nsig_ldp    = new RooPoisson( "pdf_Nsig_ldp"    , "Nsig_ldp Poisson PDF"    , *rv_Nsig_ldp    , *rv_n_sig_ldp ) ;
      pdf_Nsb_ldp     = new RooPoisson( "pdf_Nsb_ldp"     , "Nsb_ldp Poisson PDF"     , *rv_Nsb_ldp     , *rv_n_sb_ldp ) ;
      pdf_Nsig_sl     = new RooPoisson( "pdf_Nsig_sl"     , "Nsig_sl Poisson PDF"     , *rv_Nsig_sl     , *rv_n_sig_sl ) ;
      pdf_Nsb_sl      = new RooPoisson( "pdf_Nsb_sl"      , "Nsb_sl Poisson PDF"      , *rv_Nsb_sl      , *rv_n_sb_sl ) ;
      pdf_Nlsb_0b     = new RooPoisson( "pdf_Nlsb_0b"     , "Nlsb_0b Poisson PDF"     , *rv_Nlsb_0b     , *rv_n_lsb_0b ) ;
      pdf_Nlsb_0b_ldp = new RooPoisson( "pdf_Nlsb_0b_ldp" , "Nlsb_0b_ldp Poisson PDF" , *rv_Nlsb_0b_ldp , *rv_n_lsb_0b_ldp ) ;

      if ( znnModel == 1 ) {
         pdf_Nsig_ee     = new RooPoisson( "pdf_Nsig_ee"     , "Nsig_ee Poisson PDF"     , *rv_Nsig_ee     , *rv_n_sig_ee ) ;
         pdf_Nsb_ee      = new RooPoisson( "pdf_Nsb_ee"      , "Nsb_ee Poisson PDF"      , *rv_Nsb_ee      , *rv_n_sb_ee ) ;
         pdf_Nsig_mm     = new RooPoisson( "pdf_Nsig_mm"     , "Nsig_mm Poisson PDF"     , *rv_Nsig_mm     , *rv_n_sig_mm ) ;
         pdf_Nsb_mm      = new RooPoisson( "pdf_Nsb_mm"      , "Nsb_mm Poisson PDF"      , *rv_Nsb_mm      , *rv_n_sb_mm ) ;
      } else if ( znnModel == 2 ) {
         pdf_Nsigsb_ee     = new RooPoisson( "pdf_Nsigsb_ee"     , "Nsigsb_ee Poisson PDF"     , *rv_Nsigsb_ee     , *rv_n_sigsb_ee ) ;
         pdf_Nsigsb_mm     = new RooPoisson( "pdf_Nsigsb_mm"     , "Nsigsb_mm Poisson PDF"     , *rv_Nsigsb_mm     , *rv_n_sigsb_mm ) ;
      }


      np_m_sf_mc = new RooRealVar("np_m_sf_mc","Mean sf_mc", sf_mc, 0.001, 10000. ) ;
      np_m_sf_mc -> setConstant( kTRUE ) ;
      np_s_sf_mc = new RooRealVar("np_s_sf_mc","Sigma sf_mc", sf_mc_err, 0.001, 10000. ) ;
      np_s_sf_mc -> setConstant( kTRUE ) ;
      pdf_sf_mc    = new RooGaussian( "pdf_sf_mc" , "Gaussian pdf for MC scale factor",
                                       *rv_sf_mc , *np_m_sf_mc , *np_s_sf_mc ) ;

      np_m_sf_qcd_sb = new RooRealVar("np_m_sf_qcd_sb", "Mean sf_qcd_sb", sf_qcd_sb, 0.001, 10000. ) ;
      np_m_sf_qcd_sb -> setConstant( kTRUE ) ;
      np_s_sf_qcd_sb = new RooRealVar("np_s_sf_qcd_sb", "Sigma sf_qcd_sb", sf_qcd_sb_err, 0.001, 10000. ) ;
      np_s_sf_qcd_sb -> setConstant( kTRUE ) ;
      pdf_sf_qcd_sb    = new RooGaussian( "pdf_sf_qcd_sb" , "Gaussian pdf for QCD pass/fail ratio SB scale factor",
                                       *rv_sf_qcd_sb , *np_m_sf_qcd_sb , *np_s_sf_qcd_sb  ) ;

      np_m_sf_qcd_sig = new RooRealVar("np_m_sf_qcd_sig", "Mean sf_qcd_sig", sf_qcd_sig, 0.001, 10000. ) ;
      np_m_sf_qcd_sig -> setConstant( kTRUE ) ;
      np_s_sf_qcd_sig = new RooRealVar("np_s_sf_qcd_sig", "Sigma sf_qcd_sig", sf_qcd_sig_err, 0.001, 10000. ) ;
      np_s_sf_qcd_sig -> setConstant( kTRUE ) ;
      pdf_sf_qcd_sig    = new RooGaussian( "pdf_sf_qcd_sig" , "Gaussian pdf for QCD pass/fail ratio SIG scale factor",
                                       *rv_sf_qcd_sig , *np_m_sf_qcd_sig , *np_s_sf_qcd_sig ) ;

      np_m_sf_ttwj_sig = new RooRealVar("np_m_sf_ttwj_sig", "Mean sf_ttwj_sig", sf_ttwj_sig, 0.001, 10000. ) ;
      np_m_sf_ttwj_sig -> setConstant( kTRUE ) ;
      np_s_sf_ttwj_sig = new RooRealVar("np_s_sf_ttwj_sig", "Sigma sf_ttwj_sig", sf_ttwj_sig_err, 0.001, 10000. ) ;
      np_s_sf_ttwj_sig -> setConstant( kTRUE ) ;
      pdf_sf_ttwj_sig    = new RooGaussian( "pdf_sf_ttwj_sig" , "Gaussian pdf for ttwj SIG/SB scale factor",
                                       *rv_sf_ttwj_sig , *np_m_sf_ttwj_sig , *np_s_sf_ttwj_sig  ) ;




      np_m_acc_ee = new RooRealVar("np_m_acc_ee","Mean acc_ee", acc_ee_mean, 0.001, 10000. ) ;
      np_m_acc_ee -> setConstant( kTRUE ) ;
      np_s_acc_ee = new RooRealVar("np_s_acc_ee","Sigma acc_ee", acc_ee_err, 0.001, 10000. ) ;
      np_s_acc_ee -> setConstant( kTRUE ) ;
      pdf_acc_ee   = new RooGaussian( "pdf_acc_ee" , "Gaussian pdf for Z to ee acceptance",
                                          *rv_acc_ee , *np_m_acc_ee  , *np_s_acc_ee  ) ;

      np_m_acc_mm = new RooRealVar("np_m_acc_mm","Mean acc_mm", acc_mm_mean, 0.001, 10000. ) ;
      np_m_acc_mm -> setConstant( kTRUE ) ;
      np_s_acc_mm = new RooRealVar("np_s_acc_mm","Sigma acc_mm", acc_mm_err, 0.001, 10000. ) ;
      np_s_acc_mm -> setConstant( kTRUE ) ;
      pdf_acc_mm   = new RooGaussian( "pdf_acc_mm" , "Gaussian pdf for Z to mm acceptance",
                                          *rv_acc_mm , *np_m_acc_mm  , *np_s_acc_mm  ) ;

      np_m_eff_ee = new RooRealVar("np_m_eff_ee","Mean eff_ee", eff_ee_mean, 0.001, 10000. ) ;
      np_m_eff_ee -> setConstant( kTRUE ) ;
      np_s_eff_ee = new RooRealVar("np_s_eff_ee","Sigma eff_ee", eff_ee_err, 0.001, 10000. ) ;
      np_s_eff_ee -> setConstant( kTRUE ) ;
      pdf_eff_ee   = new RooGaussian( "pdf_eff_ee" , "Gaussian pdf for Z to ee efficiency",
                                          *rv_eff_ee , *np_m_eff_ee  , *np_s_eff_ee  ) ;

      np_m_eff_mm = new RooRealVar("np_m_eff_mm","Mean eff_mm", eff_mm_mean, 0.001, 10000. ) ;
      np_m_eff_mm -> setConstant( kTRUE ) ;
      np_s_eff_mm = new RooRealVar("np_s_eff_mm","Sigma eff_mm", eff_mm_err, 0.001, 10000. ) ;
      np_s_eff_mm -> setConstant( kTRUE ) ;
      pdf_eff_mm   = new RooGaussian( "pdf_eff_mm" , "Gaussian pdf for Z to mm efficiency",
                                          *rv_eff_mm , *np_m_eff_mm  , *np_s_eff_mm  ) ;

      np_m_fsig_ee = new RooRealVar("np_m_fsig_ee","Mean fsig_ee", fsig_ee_mean, 0.001, 10000. ) ;
      np_m_fsig_ee -> setConstant( kTRUE ) ;
      np_s_fsig_ee = new RooRealVar("np_s_fsig_ee","Sigma fsig_ee", fsig_ee_err, 0.001, 10000. ) ;
      np_s_fsig_ee -> setConstant( kTRUE ) ;
      pdf_fsig_ee   = new RooGaussian( "pdf_fsig_ee" , "Gaussian pdf for Z to ee purity",
                                          *rv_fsig_ee , *np_m_fsig_ee  , *np_s_fsig_ee  ) ;

      np_m_fsig_mm = new RooRealVar("np_m_fsig_mm","Mean fsig_mm", fsig_mm_mean, 0.001, 10000. ) ;
      np_m_fsig_mm -> setConstant( kTRUE ) ;
      np_s_fsig_mm = new RooRealVar("np_s_fsig_mm","Sigma fsig_mm", fsig_mm_err, 0.001, 10000. ) ;
      np_s_fsig_mm -> setConstant( kTRUE ) ;
      pdf_fsig_mm   = new RooGaussian( "pdf_fsig_mm" , "Gaussian pdf for Z to mm purity",
                                          *rv_fsig_mm , *np_m_fsig_mm  , *np_s_fsig_mm  ) ;

      np_m_sf_ee = new RooRealVar("np_m_sf_ee","Mean sf_ee", sf_ee, 0.001, 10000. ) ;
      np_m_sf_ee -> setConstant( kTRUE ) ;
      np_s_sf_ee = new RooRealVar("np_s_sf_ee","Sigma sf_ee", sf_ee_err, 0.001, 10000. ) ;
      np_s_sf_ee -> setConstant( kTRUE ) ;
      pdf_sf_ee    = new RooGaussian( "pdf_sf_ee" , "Gaussian pdf for Z to ee scale factor",
                                       *rv_sf_ee , *np_m_sf_ee , *np_s_sf_ee  ) ;

      np_m_sf_mm = new RooRealVar("np_m_sf_mm","Mean sf_mm", sf_mm, 0.001, 10000. ) ;
      np_m_sf_mm -> setConstant( kTRUE ) ;
      np_s_sf_mm = new RooRealVar("np_s_sf_mm","Sigma sf_mm", sf_mm_err, 0.001, 10000. ) ;
      np_s_sf_mm -> setConstant( kTRUE ) ;
      pdf_sf_mm    = new RooGaussian( "pdf_sf_mm" , "Gaussian pdf for Z to mm scale factor",
                                       *rv_sf_mm , *np_m_sf_mm , *np_s_sf_mm  ) ;



      if ( znnModel == 2 ) {

         np_m_knn_sig = new RooRealVar("np_m_knn_sig","Mean knn_sig", knn_sig_mean, 0.001, 10000. ) ;
         np_m_knn_sig -> setConstant( kTRUE ) ;
         np_s_knn_sig = new RooRealVar("np_s_knn_sig","Sigma knn_sig", knn_sig_err, 0.001, 10000. ) ;
         np_s_knn_sig -> setConstant( kTRUE ) ;
         pdf_knn_sig = new RooGaussian( "pdf_knn_sig", "Gaussian pdf for Z to nunu loose to tight signal scale factor",
                                          *rv_knn_sig , *np_m_knn_sig  , *np_s_knn_sig  ) ;

         np_m_knn_sb = new RooRealVar("np_m_knn_sb","Mean knn_sb", knn_sb_mean, 0.001, 10000. ) ;
         np_m_knn_sb -> setConstant( kTRUE ) ;
         np_s_knn_sb = new RooRealVar("np_s_knn_sb","sbma knn_sb", knn_sb_err, 0.001, 10000. ) ;
         np_s_knn_sb -> setConstant( kTRUE ) ;
         pdf_knn_sb = new RooGaussian( "pdf_knn_sb", "Gaussian pdf for Z to nunu loose to tight SB scale factor",
                                          *rv_knn_sb , *np_m_knn_sb  , *np_s_knn_sb  ) ;

      }



      pdf_eff_sf_prim     = new RooGaussian( "pdf_eff_sf_prim", "Master Gaussian pdf for log-normal Efficiency scale factors",
                                          *rv_eff_sf_prim, *rv_eff_sf_nom , RooConst( 1.0 ) ) ;





      {
         RooArgSet pdflist ;
         pdflist.add( *pdf_Nsig        ) ;
         pdflist.add( *pdf_Nsb         ) ;
         pdflist.add( *pdf_Nsig_ldp    ) ;
         pdflist.add( *pdf_Nsb_ldp     ) ;
         pdflist.add( *pdf_Nsig_sl     ) ;
         pdflist.add( *pdf_Nsb_sl      ) ;
         pdflist.add( *pdf_Nlsb_0b     ) ;
         pdflist.add( *pdf_Nlsb_0b_ldp ) ;
         if ( znnModel == 1 ) {
            pdflist.add( *pdf_Nsig_ee     ) ;
            pdflist.add( *pdf_Nsb_ee      ) ;
            pdflist.add( *pdf_Nsig_mm     ) ;
            pdflist.add( *pdf_Nsb_mm      ) ;
         } else if ( znnModel == 2 ) {
            pdflist.add( *pdf_Nsigsb_ee   ) ;
            pdflist.add( *pdf_Nsigsb_mm   ) ;
            pdflist.add( *pdf_knn_sig     ) ;
            pdflist.add( *pdf_knn_sb      ) ;
         }
         pdflist.add( *pdf_sf_mc      ) ;
         pdflist.add( *pdf_sf_qcd_sb      ) ;
         pdflist.add( *pdf_sf_qcd_sig      ) ;
         pdflist.add( *pdf_sf_ttwj_sig      ) ;
         pdflist.add( *pdf_sf_ee      ) ;
         pdflist.add( *pdf_sf_mm      ) ;

         pdflist.add( *pdf_acc_ee      ) ;
         pdflist.add( *pdf_acc_mm      ) ;
         pdflist.add( *pdf_eff_ee      ) ;
         pdflist.add( *pdf_eff_mm      ) ;
         pdflist.add( *pdf_fsig_ee      ) ;
         pdflist.add( *pdf_fsig_mm      ) ;

         pdflist.add( *pdf_eff_sf_prim      ) ;

         likelihood = new RooProdPdf( "likelihood", "ra2b likelihood", pdflist ) ;
      }


     //---- Define the list of observables.

       observedParametersList.add( *rv_Nsig        ) ;
       observedParametersList.add( *rv_Nsb         ) ;
       observedParametersList.add( *rv_Nsig_sl     ) ;
       observedParametersList.add( *rv_Nsb_sl      ) ;
       observedParametersList.add( *rv_Nsig_ldp    ) ;
       observedParametersList.add( *rv_Nsb_ldp     ) ;
       observedParametersList.add( *rv_Nlsb_0b     ) ;
       observedParametersList.add( *rv_Nlsb_0b_ldp ) ;
       if ( znnModel == 1 ) {
          observedParametersList.add( *rv_Nsb_ee      ) ;
          observedParametersList.add( *rv_Nsig_ee     ) ;
          observedParametersList.add( *rv_Nsb_mm      ) ;
          observedParametersList.add( *rv_Nsig_mm     ) ;
       } else if ( znnModel == 2 ) {
          observedParametersList.add( *rv_Nsigsb_ee      ) ;
          observedParametersList.add( *rv_Nsigsb_mm      ) ;
       }



       dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
                                      observedParametersList ) ;
       dsObserved->add( observedParametersList ) ;


       workspace = new RooWorkspace("ra2bv4ws") ;
       workspace->import( *likelihood ) ;
       workspace->import( *dsObserved ) ;
       printf("\n\n ========== Likelihood configuration:\n\n") ;
       workspace->Print() ;

       initialized = true ;

       return true ;


    } // initialize.

  //===================================================================

    bool ra2bRoostatsClass4ln::reinitialize( ) {


       printf( "\n\n Opening input file : %s\n\n", initializeFile ) ;

       FILE* infp ;
       if ( (infp=fopen( initializeFile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", initializeFile ) ;
          return false ;
       }

       int    Nsig                  ;
       int    Nsb                   ;
       int    Nsig_sl               ;
       int    Nsb_sl                ;
       int    Nsig_ldp              ;
       int    Nsb_ldp               ;
       int    Nlsb                  ;
       int    Nlsb_ldp              ;
       int    Nlsb_0b               ;
       int    Nlsb_0b_ldp           ;
       float  Nqcdmc_sig            ;
       float  Nqcdmc_sig_err        ;
       float  Nqcdmc_sb             ;
       float  Nqcdmc_sb_err         ;
       float  Nqcdmc_sig_sl         ;
       float  Nqcdmc_sig_sl_err     ;
       float  Nqcdmc_sb_sl          ;
       float  Nqcdmc_sb_sl_err      ;
       float  Nqcdmc_sig_ldp        ;
       float  Nqcdmc_sig_ldp_err    ;
       float  Nqcdmc_sb_ldp         ;
       float  Nqcdmc_sb_ldp_err     ;
       float  Nqcdmc_lsb            ;
       float  Nqcdmc_lsb_err        ;
       float  Nqcdmc_lsb_ldp        ;
       float  Nqcdmc_lsb_ldp_err    ;
       float  Nqcdmc_lsb_0b         ;
       float  Nqcdmc_lsb_0b_err     ;
       float  Nqcdmc_lsb_0b_ldp     ;
       float  Nqcdmc_lsb_0b_ldp_err ;
       float  Nttbarmc_sig          ;
       float  Nttbarmc_sb           ;
       float  Nttbarmc_sig_sl       ;
       float  Nttbarmc_sb_sl        ;
       float  Nttbarmc_sig_ldp      ;
       float  Nttbarmc_sb_ldp       ;
       float  Nttbarmc_lsb          ;
       float  Nttbarmc_lsb_ldp      ;
       float  NWJmc_sig             ;
       float  NWJmc_sb              ;
       float  NWJmc_sig_sl          ;
       float  NWJmc_sb_sl           ;
       float  NWJmc_sig_ldp         ;
       float  NWJmc_sb_ldp          ;
       float  NWJmc_lsb             ;
       float  NWJmc_lsb_ldp         ;
       float  NZnnmc_sig            ;
       float  NZnnmc_sb             ;
       float  NZnnmc_sig_sl         ;
       float  NZnnmc_sb_sl          ;
       float  NZnnmc_sig_ldp        ;
       float  NZnnmc_sb_ldp         ;
       float  NZnnmc_lsb            ;
       float  NZnnmc_lsb_ldp        ;
       float  NEwomc_sig            ;
       float  NEwomc_sb             ;
       float  NEwomc_sig_sl         ;
       float  NEwomc_sb_sl          ;
       float  NEwomc_sig_ldp        ;
       float  NEwomc_sb_ldp         ;
       float  NEwomc_lsb            ;
       float  NEwomc_lsb_ldp        ;
       float  Nsusymc_sig           ;
       float  Nsusymc_sb            ;
       float  Nsusymc_sig_sl        ;
       float  Nsusymc_sb_sl         ;
       float  Nsusymc_sig_ldp       ;
       float  Nsusymc_sb_ldp        ;
       float  Nsusymc_lsb           ;
       float  Nsusymc_lsb_ldp       ;
       float  Nsusymc_lsb_0b        ;
       float  Nsusymc_lsb_0b_ldp    ;
       int    Nhtonlytrig_lsb_0b      ;
       int    Nhtonlytrig_lsb_0b_ldp  ;
       int    Nsb_ee                  ;
       int    Nsig_ee                 ;
       int    Nsb_mm                  ;
       int    Nsig_mm                 ;

       float  sf_mc            ;
       float  sf_mc_err        ;
       float  sf_qcd_sb        ;
       float  sf_qcd_sb_err    ;
       float  sf_qcd_sig       ;
       float  sf_qcd_sig_err   ;
       float  sf_ttwj_sig      ;
       float  sf_ttwj_sig_err  ;
       float  sf_ee            ;
       float  sf_ee_err        ;
       float  sf_mm            ;
       float  sf_mm_err        ;



       //--- read in description line.
       printf("\n\n") ;
       char c(0) ;
       while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
       printf("\n\n") ;


       char label[1000] ;

      //--- Inputs generated with gen_roostats_input4.c
      //    The order here must be consistent with the order there!

       fscanf( infp, "%s %g", label, &EffScaleFactor        ) ;// printf( "%s %g\n", label, EffScaleFactor        ) ;
       fscanf( infp, "%s %g", label, &EffScaleFactorErr     ) ;// printf( "%s %g\n", label, EffScaleFactorErr     ) ;
       fscanf( infp, "%s %d", label, &Nsig                  ) ;// printf( "%s %d\n", label, Nsig                  ) ;
       fscanf( infp, "%s %d", label, &Nsb                   ) ;// printf( "%s %d\n", label, Nsb                   ) ;
       fscanf( infp, "%s %d", label, &Nsig_sl               ) ;// printf( "%s %d\n", label, Nsig_sl               ) ;
       fscanf( infp, "%s %d", label, &Nsb_sl                ) ;// printf( "%s %d\n", label, Nsb_sl                ) ;
       fscanf( infp, "%s %d", label, &Nsig_ldp              ) ;// printf( "%s %d\n", label, Nsig_ldp              ) ;
       fscanf( infp, "%s %d", label, &Nsb_ldp               ) ;// printf( "%s %d\n", label, Nsb_ldp               ) ;
       fscanf( infp, "%s %d", label, &Nlsb                  ) ;// printf( "%s %d\n", label, Nlsb                  ) ;
       fscanf( infp, "%s %d", label, &Nlsb_ldp              ) ;// printf( "%s %d\n", label, Nlsb_ldp              ) ;
       fscanf( infp, "%s %d", label, &Nlsb_0b               ) ;// printf( "%s %d\n", label, Nlsb_0b               ) ;
       fscanf( infp, "%s %d", label, &Nlsb_0b_ldp           ) ;// printf( "%s %d\n", label, Nlsb_0b_ldp           ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig            ) ;// printf( "%s %g\n", label, Nqcdmc_sig            ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_err        ) ;// printf( "%s %g\n", label, Nqcdmc_sig_err        ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb             ) ;// printf( "%s %g\n", label, Nqcdmc_sb             ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_err         ) ;// printf( "%s %g\n", label, Nqcdmc_sb_err         ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_sl         ) ;// printf( "%s %g\n", label, Nqcdmc_sig_sl         ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_sl_err     ) ;// printf( "%s %g\n", label, Nqcdmc_sig_sl_err     ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_sl          ) ;// printf( "%s %g\n", label, Nqcdmc_sb_sl          ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_sl_err      ) ;// printf( "%s %g\n", label, Nqcdmc_sb_sl_err      ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_ldp        ) ;// printf( "%s %g\n", label, Nqcdmc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_ldp_err    ) ;// printf( "%s %g\n", label, Nqcdmc_sig_ldp_err    ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_ldp         ) ;// printf( "%s %g\n", label, Nqcdmc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_ldp_err     ) ;// printf( "%s %g\n", label, Nqcdmc_sb_ldp_err     ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb            ) ;// printf( "%s %g\n", label, Nqcdmc_lsb            ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_err        ) ;// printf( "%s %g\n", label, Nqcdmc_lsb_err        ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_ldp        ) ;// printf( "%s %g\n", label, Nqcdmc_lsb_ldp        ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_ldp_err    ) ;// printf( "%s %g\n", label, Nqcdmc_lsb_ldp_err    ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_0b         ) ;// printf( "%s %g\n", label, Nqcdmc_lsb_0b         ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_0b_err     ) ;// printf( "%s %g\n", label, Nqcdmc_lsb_0b_err     ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_0b_ldp     ) ;// printf( "%s %g\n", label, Nqcdmc_lsb_0b_ldp     ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_0b_ldp_err ) ;// printf( "%s %g\n", label, Nqcdmc_lsb_0b_ldp_err ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sig          ) ;// printf( "%s %g\n", label, Nttbarmc_sig          ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sb           ) ;// printf( "%s %g\n", label, Nttbarmc_sb           ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sig_sl       ) ;// printf( "%s %g\n", label, Nttbarmc_sig_sl       ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sb_sl        ) ;// printf( "%s %g\n", label, Nttbarmc_sb_sl        ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sig_ldp      ) ;// printf( "%s %g\n", label, Nttbarmc_sig_ldp      ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sb_ldp       ) ;// printf( "%s %g\n", label, Nttbarmc_sb_ldp       ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_lsb          ) ;// printf( "%s %g\n", label, Nttbarmc_lsb          ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_lsb_ldp      ) ;// printf( "%s %g\n", label, Nttbarmc_lsb_ldp      ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig             ) ;// printf( "%s %g\n", label, NWJmc_sig             ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb              ) ;// printf( "%s %g\n", label, NWJmc_sb              ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig_sl          ) ;// printf( "%s %g\n", label, NWJmc_sig_sl          ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb_sl           ) ;// printf( "%s %g\n", label, NWJmc_sb_sl           ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig_ldp         ) ;// printf( "%s %g\n", label, NWJmc_sig_ldp         ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb_ldp          ) ;// printf( "%s %g\n", label, NWJmc_sb_ldp          ) ;
       fscanf( infp, "%s %g", label, &NWJmc_lsb             ) ;// printf( "%s %g\n", label, NWJmc_lsb             ) ;
       fscanf( infp, "%s %g", label, &NWJmc_lsb_ldp         ) ;// printf( "%s %g\n", label, NWJmc_lsb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig            ) ;// printf( "%s %g\n", label, NZnnmc_sig            ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb             ) ;// printf( "%s %g\n", label, NZnnmc_sb             ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig_sl         ) ;// printf( "%s %g\n", label, NZnnmc_sig_sl         ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb_sl          ) ;// printf( "%s %g\n", label, NZnnmc_sb_sl          ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig_ldp        ) ;// printf( "%s %g\n", label, NZnnmc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb_ldp         ) ;// printf( "%s %g\n", label, NZnnmc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_lsb            ) ;// printf( "%s %g\n", label, NZnnmc_lsb            ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_lsb_ldp        ) ;// printf( "%s %g\n", label, NZnnmc_lsb_ldp        ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig            ) ;// printf( "%s %g\n", label, NEwomc_sig            ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb             ) ;// printf( "%s %g\n", label, NEwomc_sb             ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig_sl         ) ;// printf( "%s %g\n", label, NEwomc_sig_sl         ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb_sl          ) ;// printf( "%s %g\n", label, NEwomc_sb_sl          ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig_ldp        ) ;// printf( "%s %g\n", label, NEwomc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb_ldp         ) ;// printf( "%s %g\n", label, NEwomc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NEwomc_lsb            ) ;// printf( "%s %g\n", label, NEwomc_lsb            ) ;
       fscanf( infp, "%s %g", label, &NEwomc_lsb_ldp        ) ;// printf( "%s %g\n", label, NEwomc_lsb_ldp        ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sig           ) ;// printf( "%s %g\n", label, Nsusymc_sig           ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sb            ) ;// printf( "%s %g\n", label, Nsusymc_sb            ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sig_sl        ) ;// printf( "%s %g\n", label, Nsusymc_sig_sl        ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sb_sl         ) ;// printf( "%s %g\n", label, Nsusymc_sb_sl         ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sig_ldp       ) ;// printf( "%s %g\n", label, Nsusymc_sig_ldp       ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sb_ldp        ) ;// printf( "%s %g\n", label, Nsusymc_sb_ldp        ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_lsb           ) ;// printf( "%s %g\n", label, Nsusymc_lsb           ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_lsb_ldp       ) ;// printf( "%s %g\n", label, Nsusymc_lsb_ldp       ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_lsb_0b        ) ;// printf( "%s %g\n", label, Nsusymc_lsb_0b        ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_lsb_0b_ldp    ) ;// printf( "%s %g\n", label, Nsusymc_lsb_0b_ldp    ) ;
       fscanf( infp, "%s %d", label, &Nhtonlytrig_lsb_0b        ) ;// printf( "%s %d\n", label, Nhtonlytrig_lsb_0b        ) ;
       fscanf( infp, "%s %d", label, &Nhtonlytrig_lsb_0b_ldp    ) ;// printf( "%s %d\n", label, Nhtonlytrig_lsb_0b_ldp    ) ;
       fscanf( infp, "%s %g", label, &DataLumi                  ) ;// printf( "%s %g\n", label, DataLumi                  ) ;
       fscanf( infp, "%s %d", label, &Nsb_ee                    ) ;// printf( "%s %d\n", label, Nsb_ee                    ) ;
       fscanf( infp, "%s %d", label, &Nsig_ee                   ) ;// printf( "%s %d\n", label, Nsig_ee                   ) ;
       fscanf( infp, "%s %d", label, &Nsb_mm                    ) ;// printf( "%s %d\n", label, Nsb_mm                    ) ;
       fscanf( infp, "%s %d", label, &Nsig_mm                   ) ;// printf( "%s %d\n", label, Nsig_mm                   ) ;
       fscanf( infp, "%s %g", label, &acc_ee_mean               ) ;// printf( "%s %g\n", label, acc_ee_mean               ) ;
       fscanf( infp, "%s %g", label, &acc_ee_err                ) ;// printf( "%s %g\n", label, acc_ee_err                ) ;
       fscanf( infp, "%s %g", label, &acc_mm_mean               ) ;// printf( "%s %g\n", label, acc_mm_mean               ) ;
       fscanf( infp, "%s %g", label, &acc_mm_err                ) ;// printf( "%s %g\n", label, acc_mm_err                ) ;
       fscanf( infp, "%s %g", label, &eff_ee_mean               ) ;// printf( "%s %g\n", label, eff_ee_mean               ) ;
       fscanf( infp, "%s %g", label, &eff_ee_err                ) ;// printf( "%s %g\n", label, eff_ee_err                ) ;
       fscanf( infp, "%s %g", label, &eff_mm_mean               ) ;// printf( "%s %g\n", label, eff_mm_mean               ) ;
       fscanf( infp, "%s %g", label, &eff_mm_err                ) ;// printf( "%s %g\n", label, eff_mm_err                ) ;
       fscanf( infp, "%s %g", label, &Ztoll_lumi                ) ;// printf( "%s %g\n", label, Ztoll_lumi                ) ;
       fscanf( infp, "%s %g", label, &knn_sig_mean              ) ;// printf( "%s %g\n", label, knn_sig_mean              ) ;
       fscanf( infp, "%s %g", label, &knn_sig_err               ) ;// printf( "%s %g\n", label, knn_sig_err               ) ;
       fscanf( infp, "%s %g", label, &knn_sb_mean               ) ;// printf( "%s %g\n", label, knn_sb_mean               ) ;
       fscanf( infp, "%s %g", label, &knn_sb_err                ) ;// printf( "%s %g\n", label, knn_sb_err                ) ;
       fscanf( infp, "%s %g", label, &fsig_ee_mean              ) ;// printf( "%s %g\n", label, fsig_ee_mean              ) ;
       fscanf( infp, "%s %g", label, &fsig_ee_err               ) ;// printf( "%s %g\n", label, fsig_ee_err               ) ;
       fscanf( infp, "%s %g", label, &fsig_mm_mean              ) ;// printf( "%s %g\n", label, fsig_mm_mean              ) ;
       fscanf( infp, "%s %g", label, &fsig_mm_err               ) ;// printf( "%s %g\n", label, fsig_mm_err               ) ;

       fscanf( infp, "%s %g", label, &sf_mc                     ) ;// printf( "%s %g\n", label, sf_mc                     ) ;
       fscanf( infp, "%s %g", label, &sf_mc_err                 ) ;// printf( "%s %g\n", label, sf_mc_err                 ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb                 ) ;// printf( "%s %g\n", label, sf_qcd_sb                 ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb_err             ) ;// printf( "%s %g\n", label, sf_qcd_sb_err             ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig                ) ;// printf( "%s %g\n", label, sf_qcd_sig                ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig_err            ) ;// printf( "%s %g\n", label, sf_qcd_sig_err            ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig               ) ;// printf( "%s %g\n", label, sf_ttwj_sig               ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig_err           ) ;// printf( "%s %g\n", label, sf_ttwj_sig_err           ) ;
       fscanf( infp, "%s %g", label, &sf_ee                     ) ;// printf( "%s %g\n", label, sf_ee                     ) ;
       fscanf( infp, "%s %g", label, &sf_ee_err                 ) ;// printf( "%s %g\n", label, sf_ee_err                 ) ;
       fscanf( infp, "%s %g", label, &sf_mm                     ) ;// printf( "%s %g\n", label, sf_mm                     ) ;
       fscanf( infp, "%s %g", label, &sf_mm_err                 ) ;// printf( "%s %g\n", label, sf_mm_err                 ) ;

  //   printf("\n Done reading in %s\n\n", initializeFile ) ;
       fclose( infp ) ;

       //+++++ Owen, Aug 6.  set all lsf_ factors to 1.  Input MC counts are now all weighted.
       //                    Calculations of stat errors below are now wrong.  Fix later.
       lsf_WJmc = 1.0 ;
       lsf_Znnmc = 1.0 ;
       lsf_Ewomc = 0.0 ; // don't use ewo.
       sf_ttbarmc = 1.0 ;


       //--- Print out a nice summary of the inputs.

       float Nsm_sig         = Nttbarmc_sig          +  lsf_WJmc*NWJmc_sig          +  Nqcdmc_sig          +  lsf_Znnmc*NZnnmc_sig          +  lsf_Ewomc*NEwomc_sig         ;


 //    printf("\n\n\n") ;


       float comp_znn_sig_ee(2.) ;
       float comp_znn_sig_mm(2.) ;
       float comp_znn_sig(2.) ;
       float comp_znn_sb_ee(2.) ;
       float comp_znn_sb_mm(2.) ;
       float comp_znn_sb(2.) ;

       float comp_znn_sig_ee_err(0.) ;
       float comp_znn_sig_mm_err(0.) ;
       float comp_znn_sig_err(0.) ;
       float comp_znn_sb_ee_err(0.) ;
       float comp_znn_sb_mm_err(0.) ;
       float comp_znn_sb_err(0.) ;


       if ( znnModel == 1 ) {

         comp_znn_sig_ee = Nsig_ee * ( 5.95 * DataLumi * fsig_ee_mean ) / ( acc_ee_mean * eff_ee_mean * Ztoll_lumi ) ;
         comp_znn_sb_ee  = Nsb_ee  * ( 5.95 * DataLumi * fsig_ee_mean ) / ( acc_ee_mean * eff_ee_mean * Ztoll_lumi ) ;

         comp_znn_sig_mm = Nsig_mm * ( 5.95 * DataLumi * fsig_mm_mean ) / ( acc_mm_mean * eff_mm_mean * Ztoll_lumi ) ;
         comp_znn_sb_mm  = Nsb_mm  * ( 5.95 * DataLumi * fsig_mm_mean ) / ( acc_mm_mean * eff_mm_mean * Ztoll_lumi ) ;

         if ( Nsig_ee > 0 ) { comp_znn_sig_ee_err = comp_znn_sig_ee * sqrt( 1.0/(1.0*Nsig_ee) + pow(fsig_ee_err/fsig_ee_mean,2) + pow(acc_ee_err/acc_ee_mean,2) + pow(eff_ee_err/eff_ee_mean,2) ) ; }
         if ( Nsb_ee  > 0 ) { comp_znn_sb_ee_err  = comp_znn_sb_ee  * sqrt( 1.0/(1.0*Nsb_ee)  + pow(fsig_ee_err/fsig_ee_mean,2) + pow(acc_ee_err/acc_ee_mean,2) + pow(eff_ee_err/eff_ee_mean,2) ) ; }

         if ( Nsig_mm > 0 ) { comp_znn_sig_mm_err = comp_znn_sig_mm * sqrt( 1.0/(1.0*Nsig_mm) + pow(fsig_mm_err/fsig_mm_mean,2) + pow(acc_mm_err/acc_mm_mean,2) + pow(eff_mm_err/eff_mm_mean,2) ) ; }
         if ( Nsb_mm  > 0 ) { comp_znn_sb_mm_err  = comp_znn_sb_mm  * sqrt( 1.0/(1.0*Nsb_mm)  + pow(fsig_mm_err/fsig_mm_mean,2) + pow(acc_mm_err/acc_mm_mean,2) + pow(eff_mm_err/eff_mm_mean,2) ) ; }



       } else if ( znnModel == 2 ) {

         comp_znn_sig_ee = (Nsig_ee + Nsb_ee) * ( 5.95 * DataLumi * fsig_ee_mean * knn_sig_mean ) / ( acc_ee_mean * eff_ee_mean * Ztoll_lumi ) ;
         comp_znn_sb_ee  = (Nsig_ee + Nsb_ee) * ( 5.95 * DataLumi * fsig_ee_mean * knn_sb_mean  ) / ( acc_ee_mean * eff_ee_mean * Ztoll_lumi ) ;

         comp_znn_sig_mm = (Nsig_mm + Nsb_mm) * ( 5.95 * DataLumi * fsig_mm_mean * knn_sig_mean ) / ( acc_mm_mean * eff_mm_mean * Ztoll_lumi ) ;
         comp_znn_sb_mm  = (Nsig_mm + Nsb_mm) * ( 5.95 * DataLumi * fsig_mm_mean * knn_sb_mean  ) / ( acc_mm_mean * eff_mm_mean * Ztoll_lumi ) ;

         if ( (Nsig_ee+Nsb_ee) > 0 ) { comp_znn_sig_ee_err = comp_znn_sig_ee * sqrt( 1.0/(1.0*(Nsig_ee+Nsb_ee)) + pow(fsig_ee_err/fsig_ee_mean,2) + pow(acc_ee_err/acc_ee_mean,2) + pow(eff_ee_err/eff_ee_mean,2) + pow(knn_sig_err/knn_sig_mean,2) ) ; }
         if ( (Nsig_ee+Nsb_ee) > 0 ) { comp_znn_sb_ee_err  = comp_znn_sb_ee  * sqrt( 1.0/(1.0*(Nsig_ee+Nsb_ee)) + pow(fsig_ee_err/fsig_ee_mean,2) + pow(acc_ee_err/acc_ee_mean,2) + pow(eff_ee_err/eff_ee_mean,2) + pow(knn_sig_err/knn_sig_mean,2) ) ; }

         if ( (Nsig_mm+Nsb_mm) > 0 ) { comp_znn_sig_mm_err = comp_znn_sig_mm * sqrt( 1.0/(1.0*(Nsig_mm+Nsb_mm)) + pow(fsig_mm_err/fsig_mm_mean,2) + pow(acc_mm_err/acc_mm_mean,2) + pow(eff_mm_err/eff_mm_mean,2) + pow(knn_sig_err/knn_sig_mean,2) ) ; }
         if ( (Nsig_mm+Nsb_mm) > 0 ) { comp_znn_sb_mm_err  = comp_znn_sb_mm  * sqrt( 1.0/(1.0*(Nsig_mm+Nsb_mm)) + pow(fsig_mm_err/fsig_mm_mean,2) + pow(acc_mm_err/acc_mm_mean,2) + pow(eff_mm_err/eff_mm_mean,2) + pow(knn_sig_err/knn_sig_mean,2) ) ; }


       }

       //-- really dumb ave.
       comp_znn_sig = 0.5 * ( comp_znn_sig_ee + comp_znn_sig_mm ) ;
       comp_znn_sig_err = comp_znn_sig_ee_err ;
       if ( comp_znn_sig_mm_err > comp_znn_sig_ee_err ) comp_znn_sig_err = comp_znn_sig_mm_err ;

       comp_znn_sb = 0.5 * ( comp_znn_sb_ee + comp_znn_sb_mm ) ;
       comp_znn_sb_err = comp_znn_sb_ee_err ;
       if ( comp_znn_sb_mm_err > comp_znn_sb_ee_err ) comp_znn_sb_err = comp_znn_sb_mm_err ;





     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




      if ( Nsig < 0 ) {
      // printf("\n\n *** Negative value for Nsig in input file.  Will set Nsig to MC expectation, which is %d.\n\n", TMath::Nint( Nsm_sig ) ) ;
         Nsig = TMath::Nint( Nsm_sig ) ;
      }

      rv_Nsig        -> setVal( Nsig ) ;
      rv_Nsb         -> setVal( Nsb ) ;

      rv_Nsig_sl     -> setVal( Nsig_sl ) ;
      rv_Nsb_sl      -> setVal( Nsb_sl ) ;

      rv_Nsig_ldp    -> setVal( Nsig_ldp ) ;
      rv_Nsb_ldp     -> setVal( Nsb_ldp ) ;

      rv_Nlsb_0b     -> setVal( Nhtonlytrig_lsb_0b ) ;
      rv_Nlsb_0b_ldp -> setVal( Nhtonlytrig_lsb_0b_ldp ) ;

      if ( znnModel == 1 ) {

         rv_Nsb_ee      -> setVal( Nsb_ee ) ;
         rv_Nsig_ee     -> setVal( Nsig_ee ) ;

         rv_Nsb_mm      -> setVal( Nsb_mm ) ;
         rv_Nsig_mm     -> setVal( Nsig_mm ) ;

      } else if ( znnModel == 2 ) {

         rv_Nsigsb_ee     -> setVal( Nsig_ee + Nsb_ee ) ;
         rv_Nsigsb_mm     -> setVal( Nsig_mm + Nsb_mm ) ;

      }







      if ( useSigTtwjVar ) {
         rrv_mu_ttwj_sig   -> setVal( Nttbarmc_sig + lsf_WJmc*NWJmc_sig ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sig  -> setVal( Nqcdmc_sig ) ; //-- this is a starting value only.
      }



      rv_mu_znn_sig   -> setVal( comp_znn_sig ) ;  //-- this is a starting value only.
      rv_mu_susy_sig    -> setVal( 0. ) ;  //-- this is a starting value only.


      if ( !useSigTtwjVar ) {
         rrv_mu_ttwj_sb   -> setVal( Nttbarmc_sb + lsf_WJmc*NWJmc_sb ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sb  -> setVal( Nqcdmc_sb ) ; //-- this is a starting value only.
      }

      if ( znnModel == 1 ) {
         rrv_mu_znn_sb   -> setVal( comp_znn_sb ) ;  //-- this is a starting value only.
      }
      rv_mu_ttwj_sig_sl  -> setVal( Nttbarmc_sig_sl + lsf_WJmc*NWJmc_sig_sl ) ;  //-- this is a starting value only.
      rv_mu_ttwj_sb_sl  -> setVal( Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl ) ;  //-- this is a starting value only.


      if ( useLdpVars ) {
         rrv_mu_qcd_sig_ldp  -> setVal( Nqcdmc_sig_ldp ) ; //-- this is a starting value only.
      }

      if ( useLdpVars ) {
         rrv_mu_qcd_sb_ldp  -> setVal( Nqcdmc_sb_ldp ) ; //-- this is a starting value only.
      }

      rv_mu_qcd_lsb_0b  -> setVal( Nhtonlytrig_lsb_0b ) ;  //-- this is a starting value only.

      rv_mu_qcd_lsb_0b_ldp  -> setVal( Nhtonlytrig_lsb_0b_ldp ) ;  //-- this is a starting value only.

      rv_mu_ttbarmc_sig_ldp  -> setVal( Nttbarmc_sig_ldp ) ;
      rv_mu_WJmc_sig_ldp     -> setVal( NWJmc_sig_ldp ) ;
      rv_mu_Znnmc_sig_ldp    -> setVal( NZnnmc_sig_ldp ) ;
      rv_mu_Ewomc_sig_ldp    -> setVal( NEwomc_sig_ldp ) ;

      rv_mu_ttbarmc_sig_ldp  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sig_ldp     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sig_ldp    -> setConstant( kTRUE ) ;
      rv_mu_Ewomc_sig_ldp    -> setConstant( kTRUE ) ;


      rv_mu_ttbarmc_sb_ldp  -> setVal( Nttbarmc_sb_ldp ) ;
      rv_mu_WJmc_sb_ldp     -> setVal( NWJmc_sb_ldp ) ;
      rv_mu_Znnmc_sb_ldp    -> setVal( NZnnmc_sb_ldp ) ;
      rv_mu_Ewomc_sb_ldp    -> setVal( NEwomc_sb_ldp ) ;

      rv_mu_ttbarmc_sb_ldp  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sb_ldp     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sb_ldp    -> setConstant( kTRUE ) ;
      rv_mu_Ewomc_sb_ldp    -> setConstant( kTRUE ) ;


      rv_mu_Ewomc_sig    -> setVal( NEwomc_sig ) ;
      rv_mu_Ewomc_sig    -> setConstant( kTRUE ) ;

      rv_mu_Ewomc_sb    -> setVal( NEwomc_sb ) ;
      rv_mu_Ewomc_sb    -> setConstant( kTRUE ) ;



      rv_eff_sf_prim  -> setVal( 0.0 ) ;

      rv_sf_mc -> setVal( sf_mc ) ;
      rv_sf_qcd_sb -> setVal( sf_qcd_sb ) ;
      rv_sf_qcd_sig -> setVal( sf_qcd_sig ) ;
      rv_sf_ttwj_sig -> setVal( sf_ttwj_sig ) ;
      rv_sf_ee -> setVal( sf_ee ) ;
      rv_sf_mm -> setVal( sf_mm ) ;

      rv_acc_ee  -> setVal( acc_ee_mean ) ;
      rv_acc_mm  -> setVal( acc_mm_mean ) ;
      rv_eff_ee  -> setVal( eff_ee_mean ) ;
      rv_eff_mm  -> setVal( eff_mm_mean ) ;
      rv_fsig_ee  -> setVal( fsig_ee_mean ) ;
      rv_fsig_mm  -> setVal( fsig_mm_mean ) ;
      if ( znnModel == 2 ) {
         rv_knn_sig -> setVal( knn_sig_mean ) ;
         rv_knn_sb  -> setVal( knn_sb_mean ) ;
      }
      rv_znnoverll_bfratio -> setVal( 5.95 ) ;
      rv_znnoverll_bfratio -> setConstant( kTRUE ) ;
      rv_dataoverll_lumiratio  -> setVal( DataLumi / Ztoll_lumi ) ;
      rv_dataoverll_lumiratio  -> setConstant( kTRUE ) ;



       initialized = true ;

       return true ;


    } // reinitialize.


  //===================================================================

    bool ra2bRoostatsClass4ln::susyScanNoContam( const char* inputScanFile, const char* outputFilebase ) {


       //-- First, (re)do the fit and susy signal profile scan.

       reinitialize() ;
       doFit() ;
       float susySigLow, susySigHigh ;
       profileSusySig( susySigLow, susySigHigh, false ) ;

       printf("  Upper limit on SUSY SIG yield : %6.1f\n\n", susySigHigh ) ;




       printf("\n\n Opening SUSY scan input file : %s\n", inputScanFile ) ;

       FILE* infp ;
       if ( (infp=fopen( inputScanFile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", inputScanFile ) ;
          return false ;
       }

       int nM0bins ;
       float minM0 ;
       float maxM0 ;
       float deltaM0 ;
       int nM12bins ;
       float minM12 ;
       float maxM12 ;
       float deltaM12 ;
       int nScanPoints ;

       char label[1000] ;

       fscanf( infp, "%s %d %f %f %f", label, &nM0bins, &minM0, &maxM0, &deltaM0 ) ;
       fscanf( infp, "%s %d %f %f %f", label, &nM12bins, &minM12, &maxM12, &deltaM12 ) ;
       fscanf( infp, "%s %d", label, &nScanPoints ) ;

       printf( "\n\n" ) ;
       printf( "  M0   :  Npoints = %4d,  min=%4.0f, max=%4.0f\n", nM0bins, minM0, maxM0 ) ;
       printf( "  M1/2 :  Npoints = %4d,  min=%4.0f, max=%4.0f\n", nM12bins, minM12, maxM12 ) ;
       printf( "\n\n" ) ;

       TH2F* hsusyscanExcluded = new TH2F("hsusyscanExcluded", "SUSY m1/2 vs m0 parameter scan",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;


       //--- read in the column headers line.
       char c(0) ;
       c = fgetc( infp ) ;
       c = 0 ;
       while ( c!=10  ) { c = fgetc( infp ) ; }

       //--- Loop over the scan points.
       for ( int pi = 0 ; pi < nScanPoints ; pi++ ) {


          float pointM0 ;
          float pointM12 ;
          float pointXsec ;
          int    n_sig ;
          int    n_sb ;
          int    n_sig_sl ;
          int    n_sb_sl ;
          int    n_sig_ldp ;
          int    n_sb_ldp ;

          fscanf( infp, "%f %f %f   %d %d %d   %d  %d  %d",
            &pointM0, &pointM12, &pointXsec,
            &n_sig, &n_sb, &n_sig_sl,
            &n_sb_sl, &n_sig_ldp, &n_sb_ldp ) ;

          int nGenPerPoint(10000) ;
          float nselWeighted = ( n_sig * pointXsec * DataLumi ) / ( 1.0*nGenPerPoint ) ;


          int m0bin  = hsusyscanExcluded->GetXaxis()->FindBin( pointM0 ) ;
          int m12bin = hsusyscanExcluded->GetYaxis()->FindBin( pointM12 ) ;

          printf(" m0 = %4.0f (%3d),  m1/2 = %4.0f (%3d),  Npred = %7.1f", pointM0, m0bin, pointM12, m12bin, nselWeighted ) ;

          //--- Owen : Do a sanity check.
          //           If Nobs is 2 or less, give the signal UL assuming zero background.
          //           From PDG
          //           Nobs  UL(95%)
          //           0     3.00
          //           1     4.74
          //           2     6.30

          if ( nselWeighted > susySigHigh ) {
             printf(" Excluded\n") ;
             hsusyscanExcluded->SetBinContent( m0bin, m12bin, 1. ) ;
          } else {
             printf("\n") ;
          }


       } // pi .

       fclose( infp ) ;


 //---------------------------------
 //    TStringLong infilestr( inputScanFile ) ;
 //    TStringLong pngoutputfilestr = infilestr ;
 //    pngoutputfilestr.ReplaceAll("input","output") ;
 //    pngoutputfilestr.ReplaceAll(".txt", outputEndname ) ;
 //    printf("\n\n png output file : %s\n\n", pngoutputfilestr.Data() ) ;

 //    TStringLong rootoutputfilestr = pngoutputfilestr ;
 //    rootoutputfilestr.ReplaceAll("png","root") ;
 //    printf("\n\n root output file : %s\n\n", rootoutputfilestr.Data() ) ;
 //---------------------------------


       char pngoutputfile[10000] ;
       sprintf( pngoutputfile, "%s.png", outputFilebase ) ;
       printf("\n\n png output file : %s\n\n", pngoutputfile ) ;

       char rootoutputfile[10000] ;
       sprintf( rootoutputfile, "%s.root", outputFilebase ) ;
       printf("\n\n root output file : %s\n\n", rootoutputfile ) ;


       gStyle->SetPadGridX(1) ;
       gStyle->SetPadGridY(1) ;
       TCanvas* csusy = new TCanvas("csusy","SUSY m1/2 vs m0 scan") ;
       hsusyscanExcluded->Draw("col") ;
       csusy->SaveAs( pngoutputfile ) ;
       TFile* f = new TFile( rootoutputfile ,"recreate") ;
       hsusyscanExcluded->Write() ;
       f->Write() ;
       f->Close() ;

       return true ;

    } // susyScanNoContam.


  //===================================================================

    bool ra2bRoostatsClass4ln::susyScanWithContam( const char* inputScanFile, const char* outputFilebase, bool isT1bbbb ) {


       if ( isT1bbbb ) {
          printf("\n\n +++++ Input is t1bbbb MC.\n\n") ;
       } else {
          printf("\n\n +++++ Input is tanb40 MC.\n\n") ;
       }



    //=== Owen, Aug 14, 2011: Adapting to Josh's input format.

       printf("\n\n Opening SUSY scan input file : %s\n", inputScanFile ) ;

       FILE* infp ;
       if ( (infp=fopen( inputScanFile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", inputScanFile ) ;
          return false ;
       }

       int nM0bins ;
       float minM0 ;
       float maxM0 ;
       float deltaM0 ;
       int nM12bins ;
       float minM12 ;
       float maxM12 ;
       float deltaM12 ;



       if ( !isT1bbbb ) {
          nM0bins = 93 ;
          minM0 = 160 ;
          maxM0 = 2000 ;
          deltaM0 = 20 ;

          nM12bins = 38 ;
          minM12 = 20 ;
          maxM12 = 760 ;
          deltaM12 = 20 ;
       } else {
          nM0bins = 60 ;
          minM0 = 0 ;
          maxM0 = 1500 ;
          deltaM0 = 25 ;

          nM12bins = 60 ;
          minM12 = 0 ;
          maxM12 = 1500 ;
          deltaM12 = 25 ;
       }


       printf( "\n\n" ) ;
       printf( "  M0   :  Npoints = %4d,  min=%4.0f, max=%4.0f\n", nM0bins, minM0, maxM0 ) ;
       printf( "  M1/2 :  Npoints = %4d,  min=%4.0f, max=%4.0f\n", nM12bins, minM12, maxM12 ) ;
       printf( "\n\n" ) ;

       TH2F* hsusyscanNgen = new TH2F("hsusyscanNgen", "SUSY m1/2 vs m0 parameter scan, Ngen",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;

       TH2F* hsusyscanEffCorr = new TH2F("hsusyscanEffCorr", "SUSY m1/2 vs m0 parameter scan, EffCorr",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;

       TH2F* hsusyscanExcluded = new TH2F("hsusyscanExcluded", "SUSY m1/2 vs m0 parameter scan",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;

       TH2F* hsusyscanNsigul = new TH2F("hsusyscanNsigul", "SUSY m1/2 vs m0 parameter scan, upper limit on Nsig",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;

       TH2F* hsusyscanNsigpred = new TH2F("hsusyscanNsigpred", "SUSY m1/2 vs m0 parameter scan, predicted Nsig",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;

       TH2F* hsusyscanEffError_sig = new TH2F("hsusyscanEffError_sig", "SUSY m1/2 vs m0 parameter scan, EffError, sig",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;

       TH2F* hsusyscanEffError_sb = new TH2F("hsusyscanEffError_sb", "SUSY m1/2 vs m0 parameter scan, EffError, sb",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;

       TH2F* hsusyscanEffError_sig_sl = new TH2F("hsusyscanEffError_sig_sl", "SUSY m1/2 vs m0 parameter scan, EffError, sig_sl",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;

       TH2F* hsusyscanEffError_sb_sl = new TH2F("hsusyscanEffError_sb_sl", "SUSY m1/2 vs m0 parameter scan, EffError, sb_sl",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;

       TH2F* hsusyscanEffError_sig_ldp = new TH2F("hsusyscanEffError_sig_ldp", "SUSY m1/2 vs m0 parameter scan, EffError, sig_ldp",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;

       TH2F* hsusyscanEffError_sb_ldp = new TH2F("hsusyscanEffError_sb_ldp", "SUSY m1/2 vs m0 parameter scan, EffError, sb_ldp",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;

       TH2F* hsusyscanXsecul(0x0) ;
       TH2F* hsusyscanEfficiency(0x0) ;

       if ( isT1bbbb ) {
          hsusyscanXsecul = new TH2F( "hsusyscanXsecul","T1bbbb Xsec upper limit",
               nM0bins, minM0, maxM0,
               nM12bins, minM12, maxM12 ) ;

          hsusyscanEfficiency = new TH2F( "hsusyscanEfficiency", "T1bbbb efficiency",
               nM0bins, minM0, maxM0,
               nM12bins, minM12, maxM12 ) ;
       }

       printf(" t1bbbb scan point:  mGluino   mLSP     Eff     NsigUL    XsecUL\n" ) ;

       //--- Loop over the scan points.
       while ( !feof( infp ) ) {

          float pointM0 ;
          float pointM12 ;



   //+++ ORL: Aug 14, 2011 ++++++++++++++++++++++++++
          float n_sig_raw ;
          float n_sb_raw ;
          float n_sig_sl_raw ;
          float n_sb_sl_raw ;
          float n_sig_ldp_raw ;
          float n_sb_ldp_raw ;

          float n_sig_correction ;
          float n_sb_correction ;
          float n_sig_sl_correction ;
          float n_sb_sl_correction ;
          float n_sig_ldp_correction ;
          float n_sb_ldp_correction ;

          float n_sig_error ;
          float n_sb_error ;
          float n_sig_sl_error ;
          float n_sb_sl_error ;
          float n_sig_ldp_error ;
          float n_sb_ldp_error ;

          int nGen ;

          fscanf( infp, "%f %f %d  %f %f %f %f %f %f  %f %f %f %f %f %f %f %f %f %f %f %f",
            &pointM0, &pointM12, &nGen,
            &n_sig_raw, &n_sb_raw, &n_sig_sl_raw, &n_sb_sl_raw, &n_sig_ldp_raw, &n_sb_ldp_raw,
            &n_sig_correction, &n_sb_correction, &n_sig_sl_correction, &n_sb_sl_correction, &n_sig_ldp_correction, &n_sb_ldp_correction,
            &n_sig_error, &n_sb_error, &n_sig_sl_error, &n_sb_sl_error, &n_sig_ldp_error, &n_sb_ldp_error ) ;

          if ( feof(infp) ) break ;
          if ( n_sig_raw < 0.00001 ) continue ;
          if ( nGen != 10000 ) continue ; // get rid of bad scan points.

           double nGenPerPoint = 10000 ; // for t1bbbb

          printf("\n\n") ;
          printf("  SUSY efficiency  systematic uncertainty,   n_sig_error     = %6.1f %%\n", n_sig_error     ) ;
          printf("  SUSY efficiency  systematic uncertainty,   n_sb_error      = %6.1f %%\n", n_sb_error      ) ;
          printf("  SUSY efficiency  systematic uncertainty,   n_sig_sl_error  = %6.1f %%\n", n_sig_sl_error  ) ;
          printf("  SUSY efficiency  systematic uncertainty,   n_sb_sl_error   = %6.1f %%\n", n_sb_sl_error   ) ;
          printf("  SUSY efficiency  systematic uncertainty,   n_sig_ldp_error = %6.1f %%\n", n_sig_ldp_error ) ;
          printf("  SUSY efficiency  systematic uncertainty,   n_sb_ldp_error  = %6.1f %%\n", n_sb_ldp_error  ) ;
          printf("\n\n") ;

   //--Include the stat error on the efficiency for t1bbbb.
          if ( isT1bbbb ) {

           //-- absolute raw eff
              float n_sig_raw_eff     = n_sig_raw     / nGenPerPoint ;
              float n_sb_raw_eff      = n_sb_raw      / nGenPerPoint ;
              float n_sig_sl_raw_eff  = n_sig_sl_raw  / nGenPerPoint ;
              float n_sb_sl_raw_eff   = n_sb_sl_raw   / nGenPerPoint ;
              float n_sig_ldp_raw_eff = n_sig_ldp_raw / nGenPerPoint ;
              float n_sb_ldp_raw_eff  = n_sb_ldp_raw  / nGenPerPoint ;


            //-- absolute stat err.
              float n_sig_stat_error     =  sqrt(  n_sig_raw_eff    * ( 1.0 - n_sig_raw_eff     ) / nGenPerPoint ) ;
              float n_sb_stat_error      =  sqrt(  n_sb_raw_eff     * ( 1.0 - n_sb_raw_eff      ) / nGenPerPoint ) ;
              float n_sig_sl_stat_error  =  sqrt(  n_sig_sl_raw_eff * ( 1.0 - n_sig_sl_raw_eff  ) / nGenPerPoint ) ;
              float n_sb_sl_stat_error   =  sqrt(  n_sb_sl_raw_eff  * ( 1.0 - n_sb_sl_raw_eff   ) / nGenPerPoint ) ;
              float n_sig_ldp_stat_error =  sqrt(  n_sig_ldp_raw_eff* ( 1.0 - n_sig_ldp_raw_eff ) / nGenPerPoint ) ;
              float n_sb_ldp_stat_error  =  sqrt(  n_sb_ldp_raw_eff * ( 1.0 - n_sb_ldp_raw_eff  ) / nGenPerPoint ) ;

            //-- relative stat err in percent.
              if ( n_sig_raw_eff     > 0 ) { n_sig_stat_error     = 100.* n_sig_stat_error     / n_sig_raw_eff     ; } else { n_sig_stat_error     = 0. ; }
              if ( n_sb_raw_eff      > 0 ) { n_sb_stat_error      = 100.* n_sb_stat_error      / n_sb_raw_eff      ; } else { n_sb_stat_error      = 0. ; }
              if ( n_sig_sl_raw_eff  > 0 ) { n_sig_sl_stat_error  = 100.* n_sig_sl_stat_error  / n_sig_sl_raw_eff  ; } else { n_sig_sl_stat_error  = 0. ; }
              if ( n_sb_sl_raw_eff   > 0 ) { n_sb_sl_stat_error   = 100.* n_sb_sl_stat_error   / n_sb_sl_raw_eff   ; } else { n_sb_sl_stat_error   = 0. ; }
              if ( n_sig_ldp_raw_eff > 0 ) { n_sig_ldp_stat_error = 100.* n_sig_ldp_stat_error / n_sig_ldp_raw_eff ; } else { n_sig_ldp_stat_error = 0. ; }
              if ( n_sb_ldp_raw_eff  > 0 ) { n_sb_ldp_stat_error  = 100.* n_sb_ldp_stat_error  / n_sb_ldp_raw_eff  ; } else { n_sb_ldp_stat_error  = 0. ; }

              printf("  SUSY efficiency  statistical uncertainty,   n_sig_stat_error     = %6.1f %%\n", n_sig_stat_error     ) ;
              printf("  SUSY efficiency  statistical uncertainty,   n_sb_stat_error      = %6.1f %%\n", n_sb_stat_error      ) ;
              printf("  SUSY efficiency  statistical uncertainty,   n_sig_sl_stat_error  = %6.1f %%\n", n_sig_sl_stat_error  ) ;
              printf("  SUSY efficiency  statistical uncertainty,   n_sb_sl_stat_error   = %6.1f %%\n", n_sb_sl_stat_error   ) ;
              printf("  SUSY efficiency  statistical uncertainty,   n_sig_ldp_stat_error = %6.1f %%\n", n_sig_ldp_stat_error ) ;
              printf("  SUSY efficiency  statistical uncertainty,   n_sb_ldp_stat_error  = %6.1f %%\n", n_sb_ldp_stat_error  ) ;

            //-- total err in percent.
              n_sig_error     = sqrt( pow( n_sig_error    , 2) + pow( n_sig_stat_error    , 2) ) ;
              n_sb_error      = sqrt( pow( n_sb_error     , 2) + pow( n_sb_stat_error     , 2) ) ;
              n_sig_sl_error  = sqrt( pow( n_sig_sl_error , 2) + pow( n_sig_sl_stat_error , 2) ) ;
              n_sb_sl_error   = sqrt( pow( n_sb_sl_error  , 2) + pow( n_sb_sl_stat_error  , 2) ) ;
              n_sig_ldp_error = sqrt( pow( n_sig_ldp_error, 2) + pow( n_sig_ldp_stat_error, 2) ) ;
              n_sb_ldp_error  = sqrt( pow( n_sb_ldp_error , 2) + pow( n_sb_ldp_stat_error , 2) ) ;

              printf("\n\n") ;
              printf("  SUSY efficiency  total uncertainty,   n_sig_error     = %6.1f %%\n", n_sig_error     ) ;
              printf("  SUSY efficiency  total uncertainty,   n_sb_error      = %6.1f %%\n", n_sb_error      ) ;
              printf("  SUSY efficiency  total uncertainty,   n_sig_sl_error  = %6.1f %%\n", n_sig_sl_error  ) ;
              printf("  SUSY efficiency  total uncertainty,   n_sb_sl_error   = %6.1f %%\n", n_sb_sl_error   ) ;
              printf("  SUSY efficiency  total uncertainty,   n_sig_ldp_error = %6.1f %%\n", n_sig_ldp_error ) ;
              printf("  SUSY efficiency  total uncertainty,   n_sb_ldp_error  = %6.1f %%\n", n_sb_ldp_error  ) ;
              printf("\n\n") ;

          }

    //-- no longer needed with log-normal efficiency.
    //
          //-- enforce a maximum efficiency uncertainty (to avoid negative scale factors).
     ///  if ( n_sig_error     > 35. ) { n_sig_error     = 35. ; }
     ///  if ( n_sb_error      > 35. ) { n_sb_error      = 35. ; }
     ///  if ( n_sig_sl_error  > 35. ) { n_sig_sl_error  = 35. ; }
     ///  if ( n_sb_sl_error   > 35. ) { n_sb_sl_error   = 35. ; }
     ///  if ( n_sig_ldp_error > 35. ) { n_sig_ldp_error = 35. ; }
     ///  if ( n_sb_ldp_error  > 35. ) { n_sb_ldp_error  = 35. ; }

   //++++++++++++++++++++++++++++++++++++++++++++++++




          int m0bin  = hsusyscanExcluded->GetXaxis()->FindBin( pointM0 ) ;
          int m12bin = hsusyscanExcluded->GetYaxis()->FindBin( pointM12 ) ;


       //--- Set up the likelihood to include the SUSY contributions to the non-SIG regions.


          reinitialize() ;


          rv_mu_susymc_sig       -> setVal( n_sig_raw     * n_sig_correction     ) ;
          rv_mu_susymc_sb        -> setVal( n_sb_raw      * n_sb_correction      ) ;
          rv_mu_susymc_sig_sl    -> setVal( n_sig_sl_raw  * n_sig_sl_correction  ) ;
          rv_mu_susymc_sb_sl     -> setVal( n_sb_sl_raw   * n_sb_sl_correction   ) ;
          rv_mu_susymc_sig_ldp   -> setVal( n_sig_ldp_raw * n_sig_ldp_correction ) ;
          rv_mu_susymc_sb_ldp    -> setVal( n_sb_ldp_raw  * n_sb_ldp_correction  ) ;

          rv_width_eff_sf_sig     -> setVal( n_sig_error     / 100. ) ;
          rv_width_eff_sf_sb      -> setVal( n_sb_error      / 100. ) ;
          rv_width_eff_sf_sig_sl  -> setVal( n_sig_sl_error  / 100. ) ;
          rv_width_eff_sf_sb_sl   -> setVal( n_sb_sl_error   / 100. ) ;
          rv_width_eff_sf_sig_ldp -> setVal( n_sig_ldp_error / 100. ) ;
          rv_width_eff_sf_sb_ldp  -> setVal( n_sb_ldp_error  / 100. ) ;

          printf("\n\n") ;
          printf(" Setting susy N_sig     to  %7.1f\n", n_sig_raw     * n_sig_correction      ) ;
          printf(" Setting susy N_sb      to  %7.1f\n", n_sb_raw      * n_sb_correction       ) ;
          printf(" Setting susy N_sig_sl  to  %7.1f\n", n_sig_sl_raw  * n_sig_sl_correction   ) ;
          printf(" Setting susy N_sb_sl   to  %7.1f\n", n_sb_sl_raw   * n_sb_sl_correction    ) ;
          printf(" Setting susy N_sig_ldp to  %7.1f\n", n_sig_ldp_raw * n_sig_ldp_correction  ) ;
          printf(" Setting susy N_sb_ldp  to  %7.1f\n", n_sb_ldp_raw  * n_sb_ldp_correction   ) ;
          printf("\n\n") ;

          parameterSnapshot() ;

          doFit() ;

          parameterSnapshot() ;

          float susySigLow, susySigHigh ;
          susySigLow = 0. ;
          susySigHigh = 9999. ;
          profileSusySig( susySigLow, susySigHigh, false ) ;



          float nselWeighted =  n_sig_raw     * n_sig_correction ;


          printf("  Upper limit on SUSY SIG yield : %6.1f\n\n", susySigHigh ) ;
          printf(" m0 = %4.0f (%d),  m1/2 = %4.0f (%d),  Npred = %7.1f", pointM0, m0bin, pointM12, m12bin, nselWeighted ) ;

          hsusyscanNsigul->SetBinContent( m0bin, m12bin, susySigHigh ) ;
          hsusyscanNsigpred->SetBinContent( m0bin, m12bin, nselWeighted ) ;
          hsusyscanNgen->SetBinContent( m0bin, m12bin, 1.0*nGen ) ;
          hsusyscanEffCorr->SetBinContent( m0bin, m12bin, n_sig_correction ) ;

          hsusyscanEffError_sig->SetBinContent( m0bin, m12bin, n_sig_error ) ;
          hsusyscanEffError_sb->SetBinContent( m0bin, m12bin, n_sb_error ) ;
          hsusyscanEffError_sig_sl->SetBinContent( m0bin, m12bin, n_sig_sl_error ) ;
          hsusyscanEffError_sb_sl->SetBinContent( m0bin, m12bin, n_sb_sl_error ) ;
          hsusyscanEffError_sig_ldp->SetBinContent( m0bin, m12bin, n_sig_ldp_error ) ;
          hsusyscanEffError_sb_ldp->SetBinContent( m0bin, m12bin, n_sb_ldp_error ) ;

          float t1bbbbEff = 0. ;
          float t1bbbbXsecUL = 999. ;

          if ( isT1bbbb ) {
             if(n_sig_raw > 0) {
                t1bbbbEff = (n_sig_raw*n_sig_correction)/10000. ;
                t1bbbbXsecUL = susySigHigh/(DataLumi*t1bbbbEff) ;
            //---------
            // hsusyscanXsecul->SetBinContent( pointM0/25.+1., pointM12/25.+1., susySigHigh*10000./n_sig_raw/1143.);
            //---------
            // hsusyscanXsecul->SetBinContent( pointM0/25.+1., pointM12/25.+1., susySigHigh/(1143.*(n_sig_raw/10000.)) );
            // hsusyscanEfficiency->SetBinContent( pointM0/25.+1., pointM12/25.+1., n_sig_raw/10000.);
            //---------
            // hsusyscanXsecul->SetBinContent( pointM0/25.+1., pointM12/25.+1., susySigHigh/(DataLumi*((n_sig_raw*n_sig_correction)/10000.)) );
            // hsusyscanEfficiency->SetBinContent( pointM0/25.+1., pointM12/25.+1., (n_sig_raw*n_sig_correction)/10000.);
            //---------
               hsusyscanXsecul->SetBinContent( pointM0/25.+1., pointM12/25.+1., t1bbbbXsecUL );
               hsusyscanEfficiency->SetBinContent( pointM0/25.+1., pointM12/25.+1., t1bbbbEff );
            //---------
             }
          }


          if ( nselWeighted > susySigHigh ) {
             if ( ! isT1bbbb ) {
                printf(" Excluded\n") ;
                hsusyscanExcluded->SetBinContent( m0bin, m12bin, 1. ) ;
             } else {
                printf("\n") ;
             }
          } else {
             printf("\n") ;
          }

          if ( isT1bbbb ) {
             //--- note that pointM0 is really mGluino and pointM12 is really mLSP for t1bbbb.
             printf(" t1bbbb scan point: %6.0f  %6.0f  %8.4f  %8.1f  %8.2f\n",
                    pointM0, pointM12, t1bbbbEff, susySigHigh, t1bbbbXsecUL ) ;
          }


       } // while not end of input file.

       fclose( infp ) ;








       char pngoutputfile[10000] ;
       sprintf( pngoutputfile, "%s.png", outputFilebase ) ;
       printf("\n\n png output file : %s\n\n", pngoutputfile ) ;

       char rootoutputfile[10000] ;
       sprintf( rootoutputfile, "%s.root", outputFilebase ) ;
       printf("\n\n root output file : %s\n\n", rootoutputfile ) ;



       gStyle->SetPadGridX(1) ;
       gStyle->SetPadGridY(1) ;
       if ( !isT1bbbb ) {
           TCanvas* csusy = new TCanvas("csusy","SUSY m1/2 vs m0 scan") ;
           hsusyscanExcluded->Draw("col") ;
           gPad->SetGridx(1) ;
           gPad->SetGridy(1) ;
           csusy->SaveAs( pngoutputfile ) ;
       } else {
           TCanvas* csusy = new TCanvas("csusy","SUSY mLSP vs mGluino scan, Cross section upper limit") ;
           hsusyscanXsecul->Draw("colz") ;
           gPad->SetGridx(1) ;
           gPad->SetGridy(1) ;
           csusy->SaveAs( pngoutputfile ) ;
       }


       TFile* f = new TFile( rootoutputfile,"recreate") ;

       hsusyscanNsigul->Write() ;
       hsusyscanNsigpred->Write() ;
       hsusyscanNgen->Write() ;
       hsusyscanEffCorr->Write() ;
       hsusyscanEffError_sig->Write() ;
       hsusyscanEffError_sb->Write() ;
       hsusyscanEffError_sig_sl->Write() ;
       hsusyscanEffError_sb_sl->Write() ;
       hsusyscanEffError_sig_ldp->Write() ;
       hsusyscanEffError_sb_ldp->Write() ;
       if ( isT1bbbb ) {
          hsusyscanXsecul->Write() ;
          hsusyscanEfficiency->Write() ;
       } else {
          hsusyscanExcluded->Write() ;
       }

       f->Write() ;
       f->Close() ;

       return true ;

    } // susyScanWithContam.

  //===================================================================

    bool ra2bRoostatsClass4ln::discoveryScanWithContam( const char* inputScanFile ) {


       printf("\n\n Opening SUSY scan input file : %s\n", inputScanFile ) ;

       FILE* infp ;
       if ( (infp=fopen( inputScanFile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", inputScanFile ) ;
          return false ;
       }

       int nM0bins ;
       float minM0 ;
       float maxM0 ;
       float deltaM0 ;
       int nM12bins ;
       float minM12 ;
       float maxM12 ;
       float deltaM12 ;
       int nScanPoints ;

       char label[1000] ;

       fscanf( infp, "%s %d %f %f %f", label, &nM0bins, &minM0, &maxM0, &deltaM0 ) ;
       fscanf( infp, "%s %d %f %f %f", label, &nM12bins, &minM12, &maxM12, &deltaM12 ) ;
       fscanf( infp, "%s %d", label, &nScanPoints ) ;

       printf( "\n\n" ) ;
       printf( "  M0   :  Npoints = %4d,  min=%4.0f, max=%4.0f\n", nM0bins, minM0, maxM0 ) ;
       printf( "  M1/2 :  Npoints = %4d,  min=%4.0f, max=%4.0f\n", nM12bins, minM12, maxM12 ) ;
       printf( "\n\n" ) ;

       TH2F* hsusyscanDeltalogl = new TH2F("hsusyscanDeltalogl", "SUSY m1/2 vs m0 parameter scan : Delta log likelihood",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;


       TH2F* hsusyscanNsigma = new TH2F("hsusyscanNsigma", "SUSY m1/2 vs m0 parameter scan : sqrt(2Delta log likelihood)",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;


       //--- read in the column headers line.
       char c(0) ;
       c = fgetc( infp ) ;
       c = 0 ;
       while ( c!=10  ) { c = fgetc( infp ) ; }

       //--- Loop over the scan points.
       for ( int pi = 0 ; pi < nScanPoints ; pi++ ) {

          float pointM0 ;
          float pointM12 ;
          float pointXsec ;
          int    n_sig ;
          int    n_sb ;
          int    n_sig_sl ;
          int    n_sb_sl ;
          int    n_sig_ldp ;
          int    n_sb_ldp ;

          fscanf( infp, "%f %f %f   %d %d %d   %d  %d  %d",
            &pointM0, &pointM12, &pointXsec,
            &n_sig, &n_sb, &n_sig_sl,
            &n_sb_sl, &n_sig_ldp, &n_sb_ldp ) ;



          int nGenPerPoint(10000) ;


          int m0bin  = hsusyscanDeltalogl->GetXaxis()->FindBin( pointM0 ) ;
          int m12bin = hsusyscanDeltalogl->GetYaxis()->FindBin( pointM12 ) ;


       //--- Set up the likelihood to include the SUSY contributions to the non-SIG regions.

          float weight = ( pointXsec * DataLumi ) / ( 1.0*nGenPerPoint ) ;

          reinitialize() ;

          rv_mu_susymc_sig     -> setVal( n_sig * weight ) ;
          rv_mu_susymc_sb      -> setVal( n_sb * weight ) ;
          rv_mu_susymc_sig_sl  -> setVal( n_sig_sl * weight ) ;
          rv_mu_susymc_sb_sl   -> setVal( n_sb_sl * weight ) ;
          rv_mu_susymc_sig_ldp -> setVal( n_sig_ldp * weight ) ;
          rv_mu_susymc_sb_ldp  -> setVal( n_sb_ldp * weight ) ;

       //--- Add susy contributions to the observables.

          rv_Nsig     -> setVal(  rv_Nsig     ->getVal() + n_sig     * weight ) ;
          rv_Nsb      -> setVal(  rv_Nsb      ->getVal() + n_sb      * weight ) ;
          rv_Nsig_sl  -> setVal(  rv_Nsig_sl  ->getVal() + n_sig_sl  * weight ) ;
          rv_Nsb_sl   -> setVal(  rv_Nsb_sl   ->getVal() + n_sb_sl   * weight ) ;
          rv_Nsig_ldp -> setVal(  rv_Nsig_ldp ->getVal() + n_sig_ldp * weight ) ;
          rv_Nsb_ldp  -> setVal(  rv_Nsb_ldp  ->getVal() + n_sb_ldp  * weight ) ;

       //--- start the susy parameter close to the true value.

          rv_mu_susy_sig -> setVal( n_sig     * weight ) ;



          parameterSnapshot() ;

          RooArgSet discFitobservedParametersList ;
          discFitobservedParametersList.add( *rv_Nsig        ) ;
          discFitobservedParametersList.add( *rv_Nsb         ) ;
          discFitobservedParametersList.add( *rv_Nsig_sl     ) ;
          discFitobservedParametersList.add( *rv_Nsb_sl      ) ;
          discFitobservedParametersList.add( *rv_Nsig_ldp    ) ;
          discFitobservedParametersList.add( *rv_Nsb_ldp     ) ;
          discFitobservedParametersList.add( *rv_Nlsb_0b     ) ;
          discFitobservedParametersList.add( *rv_Nlsb_0b_ldp ) ;
          if ( znnModel == 1 ) {
             discFitobservedParametersList.add( *rv_Nsb_ee      ) ;
             discFitobservedParametersList.add( *rv_Nsig_ee     ) ;
             discFitobservedParametersList.add( *rv_Nsb_mm      ) ;
             discFitobservedParametersList.add( *rv_Nsig_mm     ) ;
          } else if ( znnModel == 2 ) {
             discFitobservedParametersList.add( *rv_Nsigsb_ee      ) ;
             discFitobservedParametersList.add( *rv_Nsigsb_mm      ) ;
          }


          RooDataSet* discFitdsObserved = new RooDataSet("discfit_ra2b_observed_rds", "RA2b observed data values",
                                         discFitobservedParametersList ) ;
          discFitdsObserved->add( discFitobservedParametersList ) ;

          rv_mu_susy_sig->setConstant( kFALSE ) ;

          printf("\n\n") ;
          printf("  Fitting with these values for the observables.\n") ;
          discFitdsObserved->printMultiline(cout, 1, kTRUE, "") ;
          printf("\n\n") ;
          fitResult = likelihood->fitTo(*discFitdsObserved, Save(true));

          double bestLogLikelihood = fitResult->minNll() ;
          int bestFitStatus = fitResult->status() ;
          int bestFitCovQual = fitResult->covQual() ;

          printf("\n\n ++++++++ Minimum log likelihood, susy signal floating: %g, status %d, covQual %d.\n\n", bestLogLikelihood, bestFitStatus, bestFitCovQual ) ;

          parameterSnapshot() ;


          rv_mu_susy_sig->setVal(0.) ;
          rv_mu_susy_sig->setConstant( kTRUE ) ;

          fitResult = likelihood->fitTo(*discFitdsObserved, Save(true));


          double zerosusyLogLikelihood = fitResult->minNll() ;
          int nosusyFitStatus = fitResult->status() ;
          int nosusyCovQual = fitResult->covQual() ;

          printf("\n\n ++++++++ No SUSY log likelihood. : %g, fit status %d, covQual %d\n\n", zerosusyLogLikelihood, nosusyFitStatus, nosusyCovQual ) ;
          printf("  +++++++ Delta log likelihood : %g\n\n", zerosusyLogLikelihood - bestLogLikelihood ) ;

          parameterSnapshot() ;

          double deltaLogLikelihood = zerosusyLogLikelihood - bestLogLikelihood ;

          if ( TMath::IsNaN( deltaLogLikelihood ) == 0 && bestFitCovQual==3 && nosusyCovQual==3 ) {
             hsusyscanDeltalogl -> SetBinContent( m0bin, m12bin, deltaLogLikelihood ) ;
             if ( deltaLogLikelihood>0 ) { hsusyscanNsigma -> SetBinContent( m0bin, m12bin, sqrt(2*deltaLogLikelihood) ) ; }
          }

          delete discFitdsObserved ;

       } // pi .

       fclose( infp ) ;




       TStringLong infilestr( inputScanFile ) ;
       TStringLong pngoutputfilestr = infilestr ;
       pngoutputfilestr.ReplaceAll("input","output") ;
       pngoutputfilestr.ReplaceAll(".txt","-discovery-deltaLogL-withcontam.png") ;
       printf("\n\n png output file : %s\n\n", pngoutputfilestr.Data() ) ;

       TCanvas* cdsusy = new TCanvas("cdsusy","SUSY m1/2 vs m0 scan") ;
       hsusyscanDeltalogl->SetMinimum(0.) ;
       hsusyscanDeltalogl->SetMaximum(8.) ;
       hsusyscanDeltalogl->Draw("colz") ;
       cdsusy->SaveAs( pngoutputfilestr.Data() ) ;



       TStringLong pngoutputfilestr2 = infilestr ;
       pngoutputfilestr2.ReplaceAll("input","output") ;
       pngoutputfilestr2.ReplaceAll(".txt","-discovery-Nsigma-withcontam.png") ;
       printf("\n\n png output file : %s\n\n", pngoutputfilestr2.Data() ) ;

       hsusyscanNsigma->SetMinimum(0.) ;
       hsusyscanNsigma->SetMaximum(4.) ;
       hsusyscanNsigma->Draw("colz") ;
       cdsusy->SaveAs( pngoutputfilestr2.Data() ) ;



       TStringLong rootoutputfilestr = pngoutputfilestr ;
       rootoutputfilestr.ReplaceAll("png","root") ;
       printf("\n\n root output file : %s\n\n", rootoutputfilestr.Data() ) ;
       TFile* f = new TFile( rootoutputfilestr.Data(),"recreate") ;
       hsusyscanDeltalogl->Write() ;
       hsusyscanNsigma->Write() ;
       f->Write() ;
       f->Close() ;

       return true ;

    } // discoveryScanWithContam.
  //===================================================================

    bool ra2bRoostatsClass4ln::nosusysignifScanWithContam( const char* inputScanFile ) {


       printf("\n\n Opening SUSY scan input file : %s\n", inputScanFile ) ;

       FILE* infp ;
       if ( (infp=fopen( inputScanFile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", inputScanFile ) ;
          return false ;
       }

       int nM0bins ;
       float minM0 ;
       float maxM0 ;
       float deltaM0 ;
       int nM12bins ;
       float minM12 ;
       float maxM12 ;
       float deltaM12 ;
       int nScanPoints ;

       char label[1000] ;

       fscanf( infp, "%s %d %f %f %f", label, &nM0bins, &minM0, &maxM0, &deltaM0 ) ;
       fscanf( infp, "%s %d %f %f %f", label, &nM12bins, &minM12, &maxM12, &deltaM12 ) ;
       fscanf( infp, "%s %d", label, &nScanPoints ) ;

       printf( "\n\n" ) ;
       printf( "  M0   :  Npoints = %4d,  min=%4.0f, max=%4.0f\n", nM0bins, minM0, maxM0 ) ;
       printf( "  M1/2 :  Npoints = %4d,  min=%4.0f, max=%4.0f\n", nM12bins, minM12, maxM12 ) ;
       printf( "\n\n" ) ;

       TH2F* hsusyscanNosusysignif = new TH2F("hsusyscanNosusysignif", "SUSY m1/2 vs m0 parameter scan : Delta log likelihood",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;


       TH2F* hsusyscanNosusyNsigma = new TH2F("hsusyscanNosusyNsigma", "SUSY m1/2 vs m0 parameter scan : sqrt(2Delta log likelihood)",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;


       //--- read in the column headers line.
       char c(0) ;
       c = fgetc( infp ) ;
       c = 0 ;
       while ( c!=10  ) { c = fgetc( infp ) ; }

       //--- Loop over the scan points.
       for ( int pi = 0 ; pi < nScanPoints ; pi++ ) {

          float pointM0 ;
          float pointM12 ;
          float pointXsec ;
          int    n_sig ;
          int    n_sb ;
          int    n_sig_sl ;
          int    n_sb_sl ;
          int    n_sig_ldp ;
          int    n_sb_ldp ;

          fscanf( infp, "%f %f %f   %d %d %d   %d  %d  %d",
            &pointM0, &pointM12, &pointXsec,
            &n_sig, &n_sb, &n_sig_sl,
            &n_sb_sl, &n_sig_ldp, &n_sb_ldp ) ;



          int nGenPerPoint(10000) ;


          int m0bin  = hsusyscanNosusysignif->GetXaxis()->FindBin( pointM0 ) ;
          int m12bin = hsusyscanNosusysignif->GetYaxis()->FindBin( pointM12 ) ;


       //--- Set up the likelihood to include the SUSY contributions to the non-SIG regions.

          float weight = ( pointXsec * DataLumi ) / ( 1.0*nGenPerPoint ) ;

          reinitialize() ;

          rv_mu_susymc_sig     -> setVal( n_sig * weight ) ;
          rv_mu_susymc_sb      -> setVal( n_sb * weight ) ;
          rv_mu_susymc_sig_sl  -> setVal( n_sig_sl * weight ) ;
          rv_mu_susymc_sb_sl   -> setVal( n_sb_sl * weight ) ;
          rv_mu_susymc_sig_ldp -> setVal( n_sig_ldp * weight ) ;
          rv_mu_susymc_sb_ldp  -> setVal( n_sb_ldp * weight ) ;

          rv_mu_susy_sig -> setVal( 0.0 ) ;



          parameterSnapshot() ;

          RooArgSet discFitobservedParametersList ;
          discFitobservedParametersList.add( *rv_Nsig        ) ;
          discFitobservedParametersList.add( *rv_Nsb         ) ;
          discFitobservedParametersList.add( *rv_Nsig_sl     ) ;
          discFitobservedParametersList.add( *rv_Nsb_sl      ) ;
          discFitobservedParametersList.add( *rv_Nsig_ldp    ) ;
          discFitobservedParametersList.add( *rv_Nsb_ldp     ) ;
          discFitobservedParametersList.add( *rv_Nlsb_0b     ) ;
          discFitobservedParametersList.add( *rv_Nlsb_0b_ldp ) ;
          if ( znnModel == 1 ) {
             discFitobservedParametersList.add( *rv_Nsb_ee      ) ;
             discFitobservedParametersList.add( *rv_Nsig_ee     ) ;
             discFitobservedParametersList.add( *rv_Nsb_mm      ) ;
             discFitobservedParametersList.add( *rv_Nsig_mm     ) ;
          } else if ( znnModel == 2 ) {
             discFitobservedParametersList.add( *rv_Nsigsb_ee      ) ;
             discFitobservedParametersList.add( *rv_Nsigsb_mm      ) ;
          }


          RooDataSet* discFitdsObserved = new RooDataSet("discfit_ra2b_observed_rds", "RA2b observed data values",
                                         discFitobservedParametersList ) ;
          discFitdsObserved->add( discFitobservedParametersList ) ;

          rv_mu_susy_sig->setConstant( kFALSE ) ;

          printf("\n\n") ;
          printf("  Fitting with these values for the observables.\n") ;
          discFitdsObserved->printMultiline(cout, 1, kTRUE, "") ;
          printf("\n\n") ;
          fitResult = likelihood->fitTo(*discFitdsObserved, Save(true));

          double bestLogLikelihood = fitResult->minNll() ;
          int bestFitStatus = fitResult->status() ;
          int bestFitCovQual = fitResult->covQual() ;

          printf("\n\n ++++++++ Minimum log likelihood, susy signal floating: %g, status %d, covQual %d.\n\n", bestLogLikelihood, bestFitStatus, bestFitCovQual ) ;

          parameterSnapshot() ;


          rv_mu_susy_sig->setVal(0.) ;
          rv_mu_susy_sig->setConstant( kTRUE ) ;

          fitResult = likelihood->fitTo(*discFitdsObserved, Save(true));


          double zerosusyLogLikelihood = fitResult->minNll() ;
          int nosusyFitStatus = fitResult->status() ;
          int nosusyCovQual = fitResult->covQual() ;

          printf("\n\n ++++++++ No SUSY log likelihood. : %g, fit status %d, covQual %d\n\n", zerosusyLogLikelihood, nosusyFitStatus, nosusyCovQual ) ;
          printf("  +++++++ Delta log likelihood : %g\n\n", zerosusyLogLikelihood - bestLogLikelihood ) ;

          parameterSnapshot() ;

          double deltaLogLikelihood = zerosusyLogLikelihood - bestLogLikelihood ;

          if ( TMath::IsNaN( deltaLogLikelihood ) == 0 && bestFitCovQual==3 && nosusyCovQual==3 ) {
             hsusyscanNosusysignif -> SetBinContent( m0bin, m12bin, deltaLogLikelihood ) ;
             if ( deltaLogLikelihood>0 ) { hsusyscanNosusyNsigma -> SetBinContent( m0bin, m12bin, sqrt(2*deltaLogLikelihood) ) ; }
          }

          delete discFitdsObserved ;

       } // pi .

       fclose( infp ) ;




       TStringLong infilestr( inputScanFile ) ;
       TStringLong pngoutputfilestr = infilestr ;
       pngoutputfilestr.ReplaceAll("input","output") ;
       pngoutputfilestr.ReplaceAll(".txt","-nosusy-deltaLogL-withcontam.png") ;
       printf("\n\n png output file : %s\n\n", pngoutputfilestr.Data() ) ;

       TCanvas* cdsusy = new TCanvas("cdsusy","SUSY m1/2 vs m0 scan") ;
       hsusyscanNosusysignif->SetMinimum(0.) ;
       hsusyscanNosusysignif->SetMaximum(8.) ;
       hsusyscanNosusysignif->Draw("colz") ;
       cdsusy->SaveAs( pngoutputfilestr.Data() ) ;



       TStringLong pngoutputfilestr2 = infilestr ;
       pngoutputfilestr2.ReplaceAll("input","output") ;
       pngoutputfilestr2.ReplaceAll(".txt","-nosusy-Nsigma-withcontam.png") ;
       printf("\n\n png output file : %s\n\n", pngoutputfilestr2.Data() ) ;

       hsusyscanNosusyNsigma->SetMinimum(0.) ;
       hsusyscanNosusyNsigma->SetMaximum(4.) ;
       hsusyscanNosusyNsigma->Draw("colz") ;
       cdsusy->SaveAs( pngoutputfilestr2.Data() ) ;



       TStringLong rootoutputfilestr = pngoutputfilestr ;
       rootoutputfilestr.ReplaceAll("png","root") ;
       printf("\n\n root output file : %s\n\n", rootoutputfilestr.Data() ) ;
       TFile* f = new TFile( rootoutputfilestr.Data(),"recreate") ;
       hsusyscanNosusysignif->Write() ;
       hsusyscanNosusyNsigma->Write() ;
       f->Write() ;
       f->Close() ;

       return true ;

    } // discoveryScanWithContam.

  //===================================================================

    bool ra2bRoostatsClass4ln::setSusyScanPoint( const char* inputScanFile, double m0, double m12, bool isT1bbbb, double t1bbbbXsec ) {


       //--- Aug 15, 2011: updated to new format for AN, v3.


       printf("\n\n Opening SUSY scan input file : %s\n", inputScanFile ) ;

       FILE* infp ;
       if ( (infp=fopen( inputScanFile,"r"))==NULL ) {
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

       //--- Loop over the scan points.
       while ( !feof( infp ) ) {

          float pointM0 ;
          float pointM12 ;



   //+++ ORL: Aug 14, 2011 ++++++++++++++++++++++++++
          float n_sig_raw ;
          float n_sb_raw ;
          float n_sig_sl_raw ;
          float n_sb_sl_raw ;
          float n_sig_ldp_raw ;
          float n_sb_ldp_raw ;

          float n_sig_correction ;
          float n_sb_correction ;
          float n_sig_sl_correction ;
          float n_sb_sl_correction ;
          float n_sig_ldp_correction ;
          float n_sb_ldp_correction ;

          float n_sig_error ;
          float n_sb_error ;
          float n_sig_sl_error ;
          float n_sb_sl_error ;
          float n_sig_ldp_error ;
          float n_sb_ldp_error ;

          int nGen ;

          fscanf( infp, "%f %f %d  %f %f %f %f %f %f  %f %f %f %f %f %f %f %f %f %f %f %f",
            &pointM0, &pointM12, &nGen,
            &n_sig_raw, &n_sb_raw, &n_sig_sl_raw, &n_sb_sl_raw, &n_sig_ldp_raw, &n_sb_ldp_raw,
            &n_sig_correction, &n_sb_correction, &n_sig_sl_correction, &n_sb_sl_correction, &n_sig_ldp_correction, &n_sb_ldp_correction,
            &n_sig_error, &n_sb_error, &n_sig_sl_error, &n_sb_sl_error, &n_sig_ldp_error, &n_sb_ldp_error ) ;

          if ( feof(infp) ) break ;
          if ( n_sig_raw < 0.00001 ) continue ;
   //--If you are asking for it, I'll assume it's good.  Josh is using 0 for ngen dummy in LM9.
   /////  if ( nGen != 10000 ) continue ; // get rid of bad scan points.

          if (    fabs( pointM0 - m0 ) <= deltaM0/2.
               && fabs( pointM12 - m12 ) <= deltaM12/2. ) {

              double nGenPerPoint = 10000 ; // for t1bbbb

             printf("\n\n") ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_error     = %6.1f %%\n", n_sig_error     ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_error      = %6.1f %%\n", n_sb_error      ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_sl_error  = %6.1f %%\n", n_sig_sl_error  ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_sl_error   = %6.1f %%\n", n_sb_sl_error   ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_ldp_error = %6.1f %%\n", n_sig_ldp_error ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_ldp_error  = %6.1f %%\n", n_sb_ldp_error  ) ;
             printf("\n\n") ;

      //--Include the stat error on the efficiency for t1bbbb.
             if ( isT1bbbb ) {

              //-- absolute raw eff
                 float n_sig_raw_eff     = n_sig_raw     / nGenPerPoint ;
                 float n_sb_raw_eff      = n_sb_raw      / nGenPerPoint ;
                 float n_sig_sl_raw_eff  = n_sig_sl_raw  / nGenPerPoint ;
                 float n_sb_sl_raw_eff   = n_sb_sl_raw   / nGenPerPoint ;
                 float n_sig_ldp_raw_eff = n_sig_ldp_raw / nGenPerPoint ;
                 float n_sb_ldp_raw_eff  = n_sb_ldp_raw  / nGenPerPoint ;


               //-- absolute stat err.
                 float n_sig_stat_error     =  sqrt(  n_sig_raw_eff    * ( 1.0 - n_sig_raw_eff     ) / nGenPerPoint ) ;
                 float n_sb_stat_error      =  sqrt(  n_sb_raw_eff     * ( 1.0 - n_sb_raw_eff      ) / nGenPerPoint ) ;
                 float n_sig_sl_stat_error  =  sqrt(  n_sig_sl_raw_eff * ( 1.0 - n_sig_sl_raw_eff  ) / nGenPerPoint ) ;
                 float n_sb_sl_stat_error   =  sqrt(  n_sb_sl_raw_eff  * ( 1.0 - n_sb_sl_raw_eff   ) / nGenPerPoint ) ;
                 float n_sig_ldp_stat_error =  sqrt(  n_sig_ldp_raw_eff* ( 1.0 - n_sig_ldp_raw_eff ) / nGenPerPoint ) ;
                 float n_sb_ldp_stat_error  =  sqrt(  n_sb_ldp_raw_eff * ( 1.0 - n_sb_ldp_raw_eff  ) / nGenPerPoint ) ;

               //-- relative stat err in percent.
                 if ( n_sig_raw_eff     > 0 ) { n_sig_stat_error     = 100.* n_sig_stat_error     / n_sig_raw_eff     ; } else { n_sig_stat_error     = 0. ; }
                 if ( n_sb_raw_eff      > 0 ) { n_sb_stat_error      = 100.* n_sb_stat_error      / n_sb_raw_eff      ; } else { n_sb_stat_error      = 0. ; }
                 if ( n_sig_sl_raw_eff  > 0 ) { n_sig_sl_stat_error  = 100.* n_sig_sl_stat_error  / n_sig_sl_raw_eff  ; } else { n_sig_sl_stat_error  = 0. ; }
                 if ( n_sb_sl_raw_eff   > 0 ) { n_sb_sl_stat_error   = 100.* n_sb_sl_stat_error   / n_sb_sl_raw_eff   ; } else { n_sb_sl_stat_error   = 0. ; }
                 if ( n_sig_ldp_raw_eff > 0 ) { n_sig_ldp_stat_error = 100.* n_sig_ldp_stat_error / n_sig_ldp_raw_eff ; } else { n_sig_ldp_stat_error = 0. ; }
                 if ( n_sb_ldp_raw_eff  > 0 ) { n_sb_ldp_stat_error  = 100.* n_sb_ldp_stat_error  / n_sb_ldp_raw_eff  ; } else { n_sb_ldp_stat_error  = 0. ; }

                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_stat_error     = %6.1f %%\n", n_sig_stat_error     ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_stat_error      = %6.1f %%\n", n_sb_stat_error      ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_sl_stat_error  = %6.1f %%\n", n_sig_sl_stat_error  ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_sl_stat_error   = %6.1f %%\n", n_sb_sl_stat_error   ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_ldp_stat_error = %6.1f %%\n", n_sig_ldp_stat_error ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_ldp_stat_error  = %6.1f %%\n", n_sb_ldp_stat_error  ) ;

               //-- total err in percent.
                 n_sig_error     = sqrt( pow( n_sig_error    , 2) + pow( n_sig_stat_error    , 2) ) ;
                 n_sb_error      = sqrt( pow( n_sb_error     , 2) + pow( n_sb_stat_error     , 2) ) ;
                 n_sig_sl_error  = sqrt( pow( n_sig_sl_error , 2) + pow( n_sig_sl_stat_error , 2) ) ;
                 n_sb_sl_error   = sqrt( pow( n_sb_sl_error  , 2) + pow( n_sb_sl_stat_error  , 2) ) ;
                 n_sig_ldp_error = sqrt( pow( n_sig_ldp_error, 2) + pow( n_sig_ldp_stat_error, 2) ) ;
                 n_sb_ldp_error  = sqrt( pow( n_sb_ldp_error , 2) + pow( n_sb_ldp_stat_error , 2) ) ;

                 printf("\n\n") ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_error     = %6.1f %%\n", n_sig_error     ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_error      = %6.1f %%\n", n_sb_error      ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_sl_error  = %6.1f %%\n", n_sig_sl_error  ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_sl_error   = %6.1f %%\n", n_sb_sl_error   ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_ldp_error = %6.1f %%\n", n_sig_ldp_error ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_ldp_error  = %6.1f %%\n", n_sb_ldp_error  ) ;
                 printf("\n\n") ;

             }

       //--- Not needed with log-normal
        ///  //-- enforce a maximum efficiency uncertainty (to avoid negative scale factors).
        ///  if ( n_sig_error     > 35. ) { n_sig_error     = 35. ; }
        ///  if ( n_sb_error      > 35. ) { n_sb_error      = 35. ; }
        ///  if ( n_sig_sl_error  > 35. ) { n_sig_sl_error  = 35. ; }
        ///  if ( n_sb_sl_error   > 35. ) { n_sb_sl_error   = 35. ; }
        ///  if ( n_sig_ldp_error > 35. ) { n_sig_ldp_error = 35. ; }
        ///  if ( n_sb_ldp_error  > 35. ) { n_sb_ldp_error  = 35. ; }

      //++++++++++++++++++++++++++++++++++++++++++++++++



             double setVal_n_sig(0.) ;
             double setVal_n_sb(0.) ;
             double setVal_n_sig_sl(0.) ;
             double setVal_n_sb_sl(0.) ;
             double setVal_n_sig_ldp(0.) ;
             double setVal_n_sb_ldp(0.) ;


             if ( !isT1bbbb ) {
                //-- tanb40
                setVal_n_sig     = n_sig_raw     * n_sig_correction     ;
                setVal_n_sb      = n_sb_raw      * n_sb_correction      ;
                setVal_n_sig_sl  = n_sig_sl_raw  * n_sig_sl_correction  ;
                setVal_n_sb_sl   = n_sb_sl_raw   * n_sb_sl_correction   ;
                setVal_n_sig_ldp = n_sig_ldp_raw * n_sig_ldp_correction ;
                setVal_n_sb_ldp  = n_sb_ldp_raw  * n_sb_ldp_correction  ;
             } else {
                //-- t1bbbb
                setVal_n_sig     = DataLumi * t1bbbbXsec * (( n_sig_raw     * n_sig_correction     )/ nGenPerPoint ) ;
                setVal_n_sb      = DataLumi * t1bbbbXsec * (( n_sb_raw      * n_sb_correction      )/ nGenPerPoint ) ;
                setVal_n_sig_sl  = DataLumi * t1bbbbXsec * (( n_sig_sl_raw  * n_sig_sl_correction  )/ nGenPerPoint ) ;
                setVal_n_sb_sl   = DataLumi * t1bbbbXsec * (( n_sb_sl_raw   * n_sb_sl_correction   )/ nGenPerPoint ) ;
                setVal_n_sig_ldp = DataLumi * t1bbbbXsec * (( n_sig_ldp_raw * n_sig_ldp_correction )/ nGenPerPoint ) ;
                setVal_n_sb_ldp  = DataLumi * t1bbbbXsec * (( n_sb_ldp_raw  * n_sb_ldp_correction  )/ nGenPerPoint ) ;
             }

             rv_mu_susymc_sig       -> setVal( setVal_n_sig      ) ;
             rv_mu_susymc_sb        -> setVal( setVal_n_sb       ) ;
             rv_mu_susymc_sig_sl    -> setVal( setVal_n_sig_sl   ) ;
             rv_mu_susymc_sb_sl     -> setVal( setVal_n_sb_sl    ) ;
             rv_mu_susymc_sig_ldp   -> setVal( setVal_n_sig_ldp  ) ;
             rv_mu_susymc_sb_ldp    -> setVal( setVal_n_sb_ldp   ) ;

             rv_width_eff_sf_sig     -> setVal( n_sig_error     / 100. ) ;
             rv_width_eff_sf_sb      -> setVal( n_sb_error      / 100. ) ;
             rv_width_eff_sf_sig_sl  -> setVal( n_sig_sl_error  / 100. ) ;
             rv_width_eff_sf_sb_sl   -> setVal( n_sb_sl_error   / 100. ) ;
             rv_width_eff_sf_sig_ldp -> setVal( n_sig_ldp_error / 100. ) ;
             rv_width_eff_sf_sb_ldp  -> setVal( n_sb_ldp_error  / 100. ) ;



             if ( !isT1bbbb ) {
                printf("\n\n Found point m0 = %4.0f,  m1/2 = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, setVal_n_sig ) ;
             } else {
                printf("\n\n Found point mGluino = %4.0f,  mLSP = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, setVal_n_sig ) ;
             }


             printf("\n\n") ;
             printf(" Setting susy N_sig     to  %7.1f\n", setVal_n_sig       ) ;
             printf(" Setting susy N_sb      to  %7.1f\n", setVal_n_sb        ) ;
             printf(" Setting susy N_sig_sl  to  %7.1f\n", setVal_n_sig_sl    ) ;
             printf(" Setting susy N_sb_sl   to  %7.1f\n", setVal_n_sb_sl     ) ;
             printf(" Setting susy N_sig_ldp to  %7.1f\n", setVal_n_sig_ldp   ) ;
             printf(" Setting susy N_sb_ldp  to  %7.1f\n", setVal_n_sb_ldp    ) ;
             printf("\n\n") ;

             found = true ;

             break ;

          } // point match?

       } // not eof ?

       fclose( infp ) ;

       if ( found ) {
          return true ;
       } else {
          printf("\n\n *** Point not found in scan.\n\n" ) ;
          return false ;
       }

    } // setSusyScanPoint.


  //===================================================================

  //===================================================================

    bool ra2bRoostatsClass4ln::fitQualityPlot( bool doNorm, const char* plotname, double hmax ) {

     //  Shows all of the Likelihood inputs, normalized to the data values,
     //  compared to the fit result.

     //
     //
     //            S
     //            L S
     //      S     S L                   E
     //      I S   I S                   S
     //   | |G B| |G B| |D A| | QCDMC | |F| |
     //    x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x
     //    1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
     //                      1                   2                   3
     //

      int nbins(10) ;

      TH1F* hfitqual_data = new TH1F("hfitqual_data","RA2b likelihood fit results, data",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_ttwj = new TH1F("hfitqual_ttwj","RA2b likelihood fit results, ttbar",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_qcd = new TH1F("hfitqual_qcd","RA2b likelihood fit results, QCD",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_znn = new TH1F("hfitqual_znn","RA2b likelihood fit results, Ztonunu",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_ewo = new TH1F("hfitqual_ewo","RA2b likelihood fit results, EW-other",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_susy = new TH1F("hfitqual_susy","RA2b likelihood fit results, SUSY",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_gaus = new TH1F("hfitqual_gaus","RA2b likelihood fit results, Gaussian constraints",
                            nbins, 0.5, nbins+0.5 ) ;



      hfitqual_ttwj  -> SetFillColor(kBlue-9) ;
      hfitqual_qcd   -> SetFillColor(2) ;
      hfitqual_znn   -> SetFillColor(kGreen-3) ;
      hfitqual_ewo   -> SetFillColor(kGreen-5) ;
      hfitqual_susy  -> SetFillColor(6) ;
      hfitqual_gaus  -> SetFillColor(kOrange+1) ;

      THStack* hfitqual_fit = new THStack( "hfitqual_fit", "RA2b likelihood fit results, fit" ) ;

      hfitqual_data->SetMarkerStyle(20) ;
      hfitqual_data->SetLineWidth(2) ;

      TAxis* xaxis = hfitqual_data->GetXaxis() ;

      int binIndex ;

      binIndex = 2 ;

      double dataVal ;
      double dataErr ;
      double dataErrNorm ;

      double ttwjVal ;
      double qcdVal ;
      double znnVal ;
      double ewoVal ;
      double susyVal ;

      double ttwjValNorm ;
      double qcdValNorm ;
      double znnValNorm ;
      double ewoValNorm ;
      double susyValNorm ;

      double eff_sf_sig ;
      double eff_sf_sb ;
      double eff_sf_sig_sl ;
      double eff_sf_sb_sl ;
      double eff_sf_sig_ldp ;
      double eff_sf_sb_ldp ;

      char   binLabel[1000] ;



      eff_sf_sig     = rv_eff_sf_sig->getVal() ;
      eff_sf_sb      = rv_eff_sf_sb->getVal() ;
      eff_sf_sig_sl  = rv_eff_sf_sig_sl->getVal() ;
      eff_sf_sb_sl   = rv_eff_sf_sb_sl->getVal() ;
      eff_sf_sig_ldp = rv_eff_sf_sig_ldp->getVal() ;
      eff_sf_sb_ldp  = rv_eff_sf_sb_ldp->getVal() ;


     //-- SIG -------------------------------------------------------

      sprintf( binLabel, "SIG" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_Nsig->getVal() ;
      if ( useSigTtwjVar ) {
         ttwjVal = rrv_mu_ttwj_sig->getVal() ;
      } else {
         ttwjVal = rfv_mu_ttwj_sig->getVal() ;
      }
      if ( useLdpVars ) {
         qcdVal   = rfv_mu_qcd_sig->getVal() ;
      } else {
         qcdVal   = rrv_mu_qcd_sig->getVal() ;
      }
      znnVal   = (rv_mu_znn_sig->getVal()) ;
      ewoVal   = (rv_mu_ewo_sig->getVal()) ;
      susyVal  = rv_mu_susy_sig->getVal() ;
      ewoVal   = eff_sf_sig * ewoVal ;
      susyVal = eff_sf_sig * susyVal ;


      dataErr = sqrt(dataVal) ;

      dataErrNorm = dataErr ;
      ttwjValNorm = ttwjVal ;
      qcdValNorm = qcdVal ;
      znnValNorm = znnVal ;
      ewoValNorm = ewoVal ;
      susyValNorm = susyVal ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         ttwjValNorm = ttwjVal / dataVal ;
         qcdValNorm   = qcdVal / dataVal ;
         znnValNorm    = znnVal / dataVal ;
         ewoValNorm    = ewoVal / dataVal ;
         susyValNorm  = susyVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_ttwj -> SetBinContent( binIndex, ttwjValNorm ) ;
      hfitqual_qcd   -> SetBinContent( binIndex, qcdValNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, znnValNorm ) ;
      hfitqual_ewo    -> SetBinContent( binIndex, ewoValNorm ) ;
      hfitqual_susy  -> SetBinContent( binIndex, susyValNorm ) ;


      binIndex++ ;





     //-- SB -------------------------------------------------------

      sprintf( binLabel, "SB" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_Nsb->getVal() ;
      if ( useSigTtwjVar ) {
         ttwjVal = rfv_mu_ttwj_sb->getVal() ;
      } else {
         ttwjVal = rrv_mu_ttwj_sb->getVal() ;
      }
      if ( useLdpVars ) {
         qcdVal   = rfv_mu_qcd_sb->getVal() ;
      } else {
         qcdVal   = rrv_mu_qcd_sb->getVal() ;
      }
      if ( znnModel == 1 ) {
         znnVal   = (rrv_mu_znn_sb->getVal()) ;
      } else {
         znnVal   = (rfv_mu_znn_sb->getVal()) ;
      }
      ewoVal   = (rv_mu_ewo_sb->getVal()) ;
      susyVal  = rv_mu_susy_sb->getVal() ;
      ewoVal   = eff_sf_sb * ewoVal ;
      susyVal = eff_sf_sb * susyVal ;


      dataErr = sqrt(dataVal) ;

      dataErrNorm = dataErr ;
      ttwjValNorm = ttwjVal ;
      qcdValNorm = qcdVal ;
      znnValNorm = znnVal ;
      ewoValNorm = ewoVal ;
      susyValNorm = susyVal ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         ttwjValNorm = ttwjVal / dataVal ;
         qcdValNorm   = qcdVal / dataVal ;
         znnValNorm    = znnVal / dataVal ;
         ewoValNorm    = ewoVal / dataVal ;
         susyValNorm  = susyVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_ttwj -> SetBinContent( binIndex, ttwjValNorm ) ;
      hfitqual_qcd   -> SetBinContent( binIndex, qcdValNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, znnValNorm ) ;
      hfitqual_ewo    -> SetBinContent( binIndex, ewoValNorm ) ;
      hfitqual_susy  -> SetBinContent( binIndex, susyValNorm ) ;


      binIndex++ ;





     //-- SIG,SL -------------------------------------------------------

      sprintf( binLabel, "SIG,SL" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_Nsig_sl->getVal() ;
      ttwjVal = rv_mu_ttwj_sig_sl->getVal() ;
      qcdVal   = 0. ;
      znnVal   = 0. ;
      ewoVal   = 0. ;
      susyVal  = rv_mu_susy_sig_sl->getVal() ;
      susyVal = eff_sf_sig_sl * susyVal ;


      dataErr = sqrt(dataVal) ;

      dataErrNorm = dataErr ;
      ttwjValNorm = ttwjVal ;
      qcdValNorm = qcdVal ;
      znnValNorm = znnVal ;
      ewoValNorm = ewoVal ;
      susyValNorm = susyVal ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         ttwjValNorm = ttwjVal / dataVal ;
         qcdValNorm   = qcdVal / dataVal ;
         znnValNorm    = znnVal / dataVal ;
         ewoValNorm    = ewoVal / dataVal ;
         susyValNorm  = susyVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_ttwj -> SetBinContent( binIndex, ttwjValNorm ) ;
      hfitqual_qcd   -> SetBinContent( binIndex, qcdValNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, znnValNorm ) ;
      hfitqual_ewo    -> SetBinContent( binIndex, ewoValNorm ) ;
      hfitqual_susy  -> SetBinContent( binIndex, susyValNorm ) ;


      binIndex++ ;





     //-- SB,SL -------------------------------------------------------

      sprintf( binLabel, "SB,SL" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_Nsb_sl->getVal() ;
      ttwjVal = rv_mu_ttwj_sb_sl->getVal() ;
      qcdVal   = 0. ;
      znnVal   = 0. ;
      ewoVal   = 0. ;
      susyVal  = rv_mu_susy_sb_sl->getVal() ;
      susyVal = eff_sf_sb_sl * susyVal ;


      dataErr = sqrt(dataVal) ;

      dataErrNorm = dataErr ;
      ttwjValNorm = ttwjVal ;
      qcdValNorm = qcdVal ;
      znnValNorm = znnVal ;
      ewoValNorm = ewoVal ;
      susyValNorm = susyVal ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         ttwjValNorm = ttwjVal / dataVal ;
         qcdValNorm   = qcdVal / dataVal ;
         znnValNorm    = znnVal / dataVal ;
         ewoValNorm    = ewoVal / dataVal ;
         susyValNorm  = susyVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_ttwj -> SetBinContent( binIndex, ttwjValNorm ) ;
      hfitqual_qcd   -> SetBinContent( binIndex, qcdValNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, znnValNorm ) ;
      hfitqual_ewo    -> SetBinContent( binIndex, ewoValNorm ) ;
      hfitqual_susy  -> SetBinContent( binIndex, susyValNorm ) ;


      binIndex++ ;





     //-- SIG,LDP -------------------------------------------------------

      sprintf( binLabel, "SIG,LDP" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_Nsig_ldp->getVal() ;
      ttwjVal = rv_mu_ttwj_sig_ldp->getVal() ;
      if ( useLdpVars ) {
         qcdVal   = rrv_mu_qcd_sig_ldp->getVal() ;
      } else {
         qcdVal   = rfv_mu_qcd_sig_ldp->getVal() ;
      }
      znnVal   = rv_mu_znn_sig_ldp->getVal() ;
      ewoVal   = (rv_mu_ewo_sig_ldp->getVal())  ;
      susyVal  = rv_mu_susy_sig_ldp->getVal() ;
      ttwjVal  = eff_sf_sig_ldp * ttwjVal ;
      znnVal   = eff_sf_sig_ldp * znnVal ;
      ewoVal   = eff_sf_sig_ldp * ewoVal ;
      susyVal = eff_sf_sig_ldp * susyVal ;


      dataErr = sqrt(dataVal) ;

      dataErrNorm = dataErr ;
      ttwjValNorm = ttwjVal ;
      qcdValNorm = qcdVal ;
      znnValNorm = znnVal ;
      ewoValNorm = ewoVal ;
      susyValNorm = susyVal ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         ttwjValNorm = ttwjVal / dataVal ;
         qcdValNorm   = qcdVal / dataVal ;
         znnValNorm    = znnVal / dataVal ;
         ewoValNorm    = ewoVal / dataVal ;
         susyValNorm  = susyVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_ttwj -> SetBinContent( binIndex, ttwjValNorm ) ;
      hfitqual_qcd   -> SetBinContent( binIndex, qcdValNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, znnValNorm ) ;
      hfitqual_ewo    -> SetBinContent( binIndex, ewoValNorm ) ;
      hfitqual_susy  -> SetBinContent( binIndex, susyValNorm ) ;


      binIndex++ ;





     //-- SB,LDP -------------------------------------------------------

      sprintf( binLabel, "SB,LDP" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_Nsb_ldp->getVal() ;
      ttwjVal = rv_mu_ttwj_sb_ldp->getVal() ;
      if ( useLdpVars ) {
         qcdVal   = rrv_mu_qcd_sb_ldp->getVal() ;
      } else {
         qcdVal   = rfv_mu_qcd_sb_ldp->getVal() ;
      }
      znnVal   = rv_mu_znn_sb_ldp->getVal() ;
      ewoVal   = (rv_mu_ewo_sb_ldp->getVal())  ;
      susyVal  = rv_mu_susy_sb_ldp->getVal() ;
      ttwjVal  = eff_sf_sb_ldp * ttwjVal ;
      znnVal   = eff_sf_sb_ldp * znnVal ;
      ewoVal   = eff_sf_sb_ldp * ewoVal ;
      susyVal = eff_sf_sb_ldp * susyVal ;


      dataErr = sqrt(dataVal) ;

      dataErrNorm = dataErr ;
      ttwjValNorm = ttwjVal ;
      qcdValNorm = qcdVal ;
      znnValNorm = znnVal ;
      ewoValNorm = ewoVal ;
      susyValNorm = susyVal ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         ttwjValNorm = ttwjVal / dataVal ;
         qcdValNorm   = qcdVal / dataVal ;
         znnValNorm    = znnVal / dataVal ;
         ewoValNorm    = ewoVal / dataVal ;
         susyValNorm  = susyVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_ttwj -> SetBinContent( binIndex, ttwjValNorm ) ;
      hfitqual_qcd   -> SetBinContent( binIndex, qcdValNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, znnValNorm ) ;
      hfitqual_ewo    -> SetBinContent( binIndex, ewoValNorm ) ;
      hfitqual_susy  -> SetBinContent( binIndex, susyValNorm ) ;


      binIndex++ ;





      binIndex++ ;

    //-- Eff SF

      double effsfScale(1.5) ;

      printf(" hist max: %g\n", hfitqual_data->GetMaximum()) ;
      double effsfsf = hmax*(hfitqual_data->GetMaximum())/(effsfScale)  ;

      sprintf( binLabel, "Eff SF" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal = 1.0 ;
      float mcVal   = rv_eff_sf_sig->getVal() ;

      dataErr = rv_width_eff_sf_sig->getVal() ;

      dataErrNorm = dataErr ;
      float mcValNorm = mcVal ;


      hfitqual_data->SetBinContent( binIndex, effsfsf*dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         mcValNorm = mcVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, effsfsf*dataErrNorm ) ;

      hfitqual_gaus -> SetBinContent( binIndex, effsfsf*mcValNorm ) ;

      printf(" %10s : err=%4.2f ;   MC = %4.2f\n", binLabel, dataErrNorm, mcValNorm ) ;

      binIndex++ ;




     //-- final formatting

      hfitqual_fit->Add( hfitqual_znn ) ;
      hfitqual_fit->Add( hfitqual_qcd ) ;
      hfitqual_fit->Add( hfitqual_ttwj ) ;
      hfitqual_fit->Add( hfitqual_susy ) ;
      hfitqual_fit->Add( hfitqual_gaus ) ;

      TLegend* legend = new TLegend(0.65,0.2,0.77,0.9) ;

      legend->AddEntry( hfitqual_data,  "data" ) ;
      legend->AddEntry( hfitqual_susy,  "SUSY" ) ;
      legend->AddEntry( hfitqual_ttwj, "ttwj" ) ;
      legend->AddEntry( hfitqual_qcd,   "QCD" ) ;
      legend->AddEntry( hfitqual_znn,   "Znunu" ) ;
      legend->AddEntry( hfitqual_gaus,  "Eff SF" ) ;

      if ( doNorm ) {
         hfitqual_data->SetMaximum( hmax ) ;
      } else {
         hfitqual_data->SetMaximum( hmax*(hfitqual_data->GetMaximum()) ) ;
      }


      double oldBM = gStyle->GetPadBottomMargin() ;
      double oldRM = gStyle->GetPadRightMargin() ;


      gStyle->SetPadBottomMargin(0.20) ;
      gStyle->SetPadRightMargin(0.40) ;


      TCanvas* cfitqual = new TCanvas("cfitqual","RA2b fit quality", 780, 550 ) ;

      gPad->SetTicks(1,0) ;


      hfitqual_data->SetLabelSize(0.055,"x") ;
      hfitqual_data->GetXaxis()->LabelsOption("v") ;

      hfitqual_data->Draw() ;
      hfitqual_fit->Draw("same") ;
      hfitqual_data->Draw("same") ;
      gPad->SetGridy(1) ;
      legend->Draw() ;



      TGaxis* axis = new TGaxis() ;
      axis->DrawAxis( nbins+0.5, 0, nbins+0.5, hfitqual_data->GetMaximum(), 0., effsfScale, 510, "+L") ;

      TLine* line = new TLine() ;
      line->SetLineWidth(2) ;

      line->DrawLine( nbins-2., 0., nbins-2., hfitqual_data->GetMaximum() ) ;



      //-- add fit results.

      TText* fittext = new TText() ;
      fittext->SetTextSize(0.045) ;
      double fitval ;
      double fiterr ;
      double tx = 0.80 ;
      double ty = 0.71 ;
      double dy = 0.115 ;
      char fitvalchars[1000] ;



      fitval = rv_mu_susy_sig->getVal() ;
      fiterr = rv_mu_susy_sig->getError() ;
      if ( fitval<10.) {
         sprintf( fitvalchars, "%3.1f +/- %3.1f", fitval, fiterr ) ;
      } else {
         sprintf( fitvalchars, "%3.0f +/- %3.0f", fitval, fiterr ) ;
      }
      fittext->DrawTextNDC( tx, ty, fitvalchars ) ;



      if ( useSigTtwjVar ) {
         fitval = rrv_mu_ttwj_sig->getVal() ;
         fiterr = rrv_mu_ttwj_sig->getError() ;
         if ( fitval<10.) {
            sprintf( fitvalchars, "%3.1f +/- %3.1f", fitval, fiterr ) ;
         } else {
            sprintf( fitvalchars, "%3.0f +/- %3.0f", fitval, fiterr ) ;
         }
      } else {
         fitval = rfv_mu_ttwj_sig->getVal() ;
         sprintf( fitvalchars, "%3.0f", fitval ) ;
         if ( fitval<10.) {
            sprintf( fitvalchars, "%3.1f", fitval ) ;
         } else {
            sprintf( fitvalchars, "%3.0f", fitval ) ;
         }
      }
      ty = ty - dy ;
      fittext->DrawTextNDC( tx, ty, fitvalchars ) ;




      if ( !useLdpVars ) {
         fitval = rrv_mu_qcd_sig->getVal() ;
         fiterr = rrv_mu_qcd_sig->getError() ;
         if ( fitval<10.) {
            sprintf( fitvalchars, "%3.1f +/- %3.1f", fitval, fiterr ) ;
         } else {
            sprintf( fitvalchars, "%3.0f +/- %3.0f", fitval, fiterr ) ;
         }
      } else {
         fitval = rfv_mu_qcd_sig->getVal() ;
         if ( fitval<10.) {
            sprintf( fitvalchars, "%3.1f", fitval ) ;
         } else {
            sprintf( fitvalchars, "%3.0f", fitval ) ;
         }
      }
      ty = ty - dy ;
      fittext->DrawTextNDC( tx, ty, fitvalchars ) ;




      fitval = rv_mu_znn_sig->getVal() ;
      fiterr = rv_mu_znn_sig->getError() ;
      if ( fitval<10.) {
         sprintf( fitvalchars, "%3.1f +/- %3.1f", fitval, fiterr ) ;
      } else {
         sprintf( fitvalchars, "%3.0f +/- %3.0f", fitval, fiterr ) ;
      }
      ty = ty - dy ;
      fittext->DrawTextNDC( tx, ty, fitvalchars ) ;




      fitval = rv_eff_sf_sig->getVal() ;
      sprintf( fitvalchars, "%4.2f", fitval ) ;
      ty = ty - dy ;
      fittext->DrawTextNDC( tx, ty, fitvalchars ) ;







      cfitqual->SaveAs( plotname ) ;


      gStyle->SetPadBottomMargin(oldBM) ;
      gStyle->SetPadRightMargin(oldRM) ;

      return true ;

    } // fitQualityPlot

  //===================================================================


    void ra2bRoostatsClass4ln::setAndFixSusySig( double setVal ) {

       rv_mu_susy_sig->setVal( setVal ) ;
       rv_mu_susy_sig->setConstant( kTRUE ) ;

       printf("\n\n Set SUSY SIG yield to %.1f and fixed it.\n\n\n", setVal ) ;

    }

  //===================================================================


    void ra2bRoostatsClass4ln::setAndFixSusySigToPredictedValue( ) {

       rv_mu_susy_sig->setVal( rv_mu_susymc_sig->getVal() ) ;
       rv_mu_susy_sig->setConstant( kTRUE ) ;

       printf("\n\n Set SUSY SIG yield to %.1f and fixed it.\n\n\n", rv_mu_susymc_sig->getVal() ) ;

    }

  //===================================================================


    void ra2bRoostatsClass4ln::freeSusySig( ) {

       rv_mu_susy_sig->setConstant( kFALSE ) ;

       printf("\n\n Freed SUSY SIG yield fit parameter.\n\n\n" ) ;

    }

  //===================================================================




    void ra2bRoostatsClass4ln::parameterSnapshot() {

       printf("\n\n====================================================================\n\n") ;

       printf("       Snapshot of observables and likelihood parameters.\n\n") ;


       printf("\n\n") ;


       printf("                     Nobs      all  |     ttwj      qcd     Znn       Ewo     SUSY \n") ;
       printf("------------------------------------+--------------------------------------------------------\n") ;

       float ttwjSig ;
       float ttwjSb ;

       if ( useSigTtwjVar ) {
          ttwjSig = rrv_mu_ttwj_sig->getVal() ;
          ttwjSb  = rfv_mu_ttwj_sb->getVal() ;
       } else {
          ttwjSig = rfv_mu_ttwj_sig->getVal() ;
          ttwjSb  = rrv_mu_ttwj_sb->getVal() ;
       }

       float qcdSig, qcdSb, qcdSigLdp, qcdSbLdp ;
       if ( useLdpVars ) {
          qcdSig    = rfv_mu_qcd_sig->getVal() ;
          qcdSb     = rfv_mu_qcd_sb->getVal() ;
          qcdSigLdp = rrv_mu_qcd_sig_ldp->getVal() ;
          qcdSbLdp  = rrv_mu_qcd_sb_ldp->getVal() ;
       } else {
          qcdSig    = rrv_mu_qcd_sig->getVal() ;
          qcdSb     = rrv_mu_qcd_sb->getVal() ;
          qcdSigLdp = rfv_mu_qcd_sig_ldp->getVal() ;
          qcdSbLdp  = rfv_mu_qcd_sb_ldp->getVal() ;
       }

       printf("                                    |\n") ;
       printf("     Nsig        :  %5.0f   %7.1f |  %7.1f  %7.1f  %7.1f  %7.1f  %7.1f (%7.1f)  \n",
              rv_Nsig->getVal(),
              rv_n_sig->getVal(),
              ttwjSig,
              qcdSig,
              rv_mu_znn_sig->getVal(),
              (rv_mu_ewo_sig->getVal())*(rv_eff_sf_sig->getVal()),
              (rv_mu_susy_sig->getVal())*(rv_eff_sf_sig->getVal()),
              rv_mu_susy_sig->getVal()
              ) ;


       double znnsbval(0.) ;
       if ( znnModel == 1 ) {
          znnsbval = rrv_mu_znn_sb->getVal() ;
       } else if (znnModel == 2 ) {
          znnsbval = rfv_mu_znn_sb->getVal() ;
       }
       printf("     Nsb         :  %5.0f   %7.1f |  %7.1f  %7.1f  %7.1f  %7.1f  %7.1f (%7.1f)  \n",
              rv_Nsb->getVal(),
              rv_n_sb->getVal(),
              ttwjSb,
              qcdSb,
              znnsbval,
              (rv_mu_ewo_sb->getVal())*(rv_eff_sf_sb->getVal()),
              (rv_mu_susy_sb->getVal())*(rv_eff_sf_sb->getVal()),
              rv_mu_susy_sb->getVal()
              ) ;


       printf("                                    |\n") ;
       printf("     Nsig,SL     :  %5.0f   %7.1f |  %7.1f  -------  -------  -------  %7.1f (%7.1f)  \n",
              rv_Nsig_sl->getVal(),
              rv_n_sig_sl->getVal(),
              rv_mu_ttwj_sig_sl->getVal(),
              (rv_mu_susy_sig_sl->getVal())*(rv_eff_sf_sig_sl->getVal()),
              rv_mu_susy_sig_sl->getVal()
              ) ;


       printf("     Nsb,SL      :  %5.0f   %7.1f |  %7.1f  -------  -------  -------  %7.1f (%7.1f)  \n",
              rv_Nsb_sl->getVal(),
              rv_n_sb_sl->getVal(),
              rv_mu_ttwj_sb_sl->getVal(),
              (rv_mu_susy_sb_sl->getVal())*(rv_eff_sf_sb_sl->getVal()),
              rv_mu_susy_sb_sl->getVal()
              ) ;


       printf("                                    |\n") ;
       printf("     Nsig,ldp    :  %5.0f   %7.1f |  %7.1f  %7.1f  %7.1f  %7.1f  %7.1f (%7.1f)  \n",
              rv_Nsig_ldp->getVal(),
              rv_n_sig_ldp->getVal(),
              rv_mu_ttwj_sig_ldp->getVal(),
              qcdSigLdp,
              rv_mu_znn_sig_ldp->getVal(),
              (rv_mu_ewo_sig_ldp->getVal())*(rv_eff_sf_sig_ldp->getVal()),
              (rv_mu_susy_sig_ldp->getVal())*(rv_eff_sf_sig_ldp->getVal()),
              rv_mu_susy_sig_ldp->getVal()
              ) ;

       printf("     Nsb,ldp     :  %5.0f   %7.1f |  %7.1f  %7.1f  %7.1f  %7.1f  %7.1f (%7.1f)  \n",
              rv_Nsb_ldp->getVal(),
              rv_n_sb_ldp->getVal(),
              rv_mu_ttwj_sb_ldp->getVal(),
              qcdSbLdp,
              rv_mu_znn_sb_ldp->getVal(),
              (rv_mu_ewo_sb_ldp->getVal())*(rv_eff_sf_sb_ldp->getVal()),
              (rv_mu_susy_sb_ldp->getVal())*(rv_eff_sf_sb_ldp->getVal()),
              rv_mu_susy_sb_ldp->getVal()
              ) ;

       printf("                                    |\n") ;
       printf("------------------------------------+--------------------------------------------------------\n") ;

       printf("\n") ;

       float delta ;

       printf("\n") ;
       delta = 0. ;
       if ( EffScaleFactorErr>0 ) { delta = (rv_eff_sf_sig->getVal()-1.0)/(rv_width_eff_sf_sig->getVal()) ; }
       printf("  Eff scale fac. :  input = %4.2f +/- %4.2f ,   fit = %4.2f, delta = %4.2f sigma.  \n",
              EffScaleFactor,
              EffScaleFactorErr,
              rv_eff_sf_sig->getVal(),
              delta ) ;

    // printf("\n") ;
    // delta = 0. ;
    // if ( lsf_WJmc_err>0 ) { delta = (rv_lsf_WJmc->getVal()-lsf_WJmc)/lsf_WJmc_err ; }
    // printf("  LSF W+jets     :  input = %4.2f +/- %4.2f ,   fit = %4.2f, delta = %4.2f sigma.  \n",
    //        lsf_WJmc,
    //        lsf_WJmc_err,
    //        rv_lsf_WJmc->getVal(),
    //        delta ) ;

    // delta = 0. ;
    // if ( lsf_Znnmc_err>0 ) { delta = (rv_lsf_Znnmc->getVal()-lsf_Znnmc)/lsf_Znnmc_err ; }
    // printf("  LSF Z->invis   :  input = %4.2f +/- %4.2f ,   fit = %4.2f, delta = %4.2f sigma.  \n",
    //        lsf_Znnmc,
    //        lsf_Znnmc_err,
    //        rv_lsf_Znnmc->getVal(),
    //        delta ) ;


       printf("\n\n====================================================================\n\n") ;

    } // parameterSnapshot


  //===================================================================


    void ra2bRoostatsClass4ln::saveToymeanSnapshot( ) {

       printf("\n\n Saving model predictions for observables.\n\n\n" ) ;

       toymean_n_sig        = rv_n_sig        ->getVal() ;
       toymean_n_sb         = rv_n_sb         ->getVal() ;
       toymean_n_sig_sl     = rv_n_sig_sl     ->getVal() ;
       toymean_n_sb_sl      = rv_n_sb_sl      ->getVal() ;
       toymean_n_sig_ldp    = rv_n_sig_ldp    ->getVal() ;
       toymean_n_sb_ldp     = rv_n_sb_ldp     ->getVal() ;
       toymean_n_lsb_0b     = rv_n_lsb_0b     ->getVal() ;
       toymean_n_lsb_0b_ldp = rv_n_lsb_0b_ldp ->getVal() ;

       if ( znnModel == 1 ) {
          toymean_n_sig_ee = rv_n_sig_ee ->getVal() ;
          toymean_n_sb_ee  = rv_n_sb_ee  ->getVal() ;
          toymean_n_sig_mm = rv_n_sig_mm ->getVal() ;
          toymean_n_sb_mm  = rv_n_sb_mm  ->getVal() ;
       } else if ( znnModel == 2 ) {
          toymean_n_sigsb_ee = rv_n_sigsb_ee ->getVal() ;
          toymean_n_sigsb_mm = rv_n_sigsb_mm ->getVal() ;
       }

       printf( " toymean_n_sig        %8.2f\n", toymean_n_sig        ) ;
       printf( " toymean_n_sb         %8.2f\n", toymean_n_sb         ) ;
       printf( " toymean_n_sig_sl     %8.2f\n", toymean_n_sig_sl     ) ;
       printf( " toymean_n_sb_sl      %8.2f\n", toymean_n_sb_sl      ) ;
       printf( " toymean_n_sig_ldp    %8.2f\n", toymean_n_sig_ldp    ) ;
       printf( " toymean_n_sb_ldp     %8.2f\n", toymean_n_sb_ldp     ) ;
       printf( " toymean_n_lsb_0b     %8.2f\n", toymean_n_lsb_0b     ) ;
       printf( " toymean_n_lsb_0b_ldp %8.2f\n", toymean_n_lsb_0b_ldp ) ;
       if ( znnModel == 1 ) {
          printf( " toymean_n_sig_ee        %8.2f\n", toymean_n_sig_ee        ) ;
          printf( " toymean_n_sb_ee         %8.2f\n", toymean_n_sb_ee         ) ;
          printf( " toymean_n_sig_mm        %8.2f\n", toymean_n_sig_mm        ) ;
          printf( " toymean_n_sb_mm         %8.2f\n", toymean_n_sb_mm         ) ;
       } else if ( znnModel == 2 ) {
          printf( " toymean_n_sigsb_ee        %8.2f\n", toymean_n_sigsb_ee        ) ;
          printf( " toymean_n_sigsb_mm        %8.2f\n", toymean_n_sigsb_mm        ) ;
       }
       printf( "\n\n\n") ;

    } // saveToymeanSnapshot

  //===================================================================


    void ra2bRoostatsClass4ln::saveDataObservableVals( ) {

       printf("\n\n Saving data values for observables.\n\n\n" ) ;

       dataval_Nsig        = rv_Nsig        ->getVal() ;
       dataval_Nsb         = rv_Nsb         ->getVal() ;
       dataval_Nsig_sl     = rv_Nsig_sl     ->getVal() ;
       dataval_Nsb_sl      = rv_Nsb_sl      ->getVal() ;
       dataval_Nsig_ldp    = rv_Nsig_ldp    ->getVal() ;
       dataval_Nsb_ldp     = rv_Nsb_ldp     ->getVal() ;
       dataval_Nlsb_0b     = rv_Nlsb_0b     ->getVal() ;
       dataval_Nlsb_0b_ldp = rv_Nlsb_0b_ldp ->getVal() ;

       if ( znnModel == 1 ) {
          dataval_Nsig_ee = rv_Nsig_ee ->getVal() ;
          dataval_Nsb_ee  = rv_Nsb_ee  ->getVal() ;
          dataval_Nsig_mm = rv_Nsig_mm ->getVal() ;
          dataval_Nsb_mm  = rv_Nsb_mm  ->getVal() ;
       } else if ( znnModel == 2 ) {
          dataval_Nsigsb_ee = rv_Nsigsb_ee ->getVal() ;
          dataval_Nsigsb_mm = rv_Nsigsb_mm ->getVal() ;
       }

       printf( " dataval_Nsig        %8.2f\n", dataval_Nsig        ) ;
       printf( " dataval_Nsb         %8.2f\n", dataval_Nsb         ) ;
       printf( " dataval_Nsig_sl     %8.2f\n", dataval_Nsig_sl     ) ;
       printf( " dataval_Nsb_sl      %8.2f\n", dataval_Nsb_sl      ) ;
       printf( " dataval_Nsig_ldp    %8.2f\n", dataval_Nsig_ldp    ) ;
       printf( " dataval_Nsb_ldp     %8.2f\n", dataval_Nsb_ldp     ) ;
       printf( " dataval_Nlsb_0b     %8.2f\n", dataval_Nlsb_0b     ) ;
       printf( " dataval_Nlsb_0b_ldp %8.2f\n", dataval_Nlsb_0b_ldp ) ;
       if ( znnModel == 1 ) {
          printf( " dataval_Nsig_ee        %8.2f\n", dataval_Nsig_ee        ) ;
          printf( " dataval_Nsb_ee         %8.2f\n", dataval_Nsb_ee         ) ;
          printf( " dataval_Nsig_mm        %8.2f\n", dataval_Nsig_mm        ) ;
          printf( " dataval_Nsb_mm         %8.2f\n", dataval_Nsb_mm         ) ;
       } else if ( znnModel == 2 ) {
          printf( " dataval_Nsigsb_ee        %8.2f\n", dataval_Nsigsb_ee        ) ;
          printf( " dataval_Nsigsb_mm        %8.2f\n", dataval_Nsigsb_mm        ) ;
       }
       printf( "\n\n\n") ;

    } // saveDataObservableVals

  //===================================================================


    void ra2bRoostatsClass4ln::genToyExperiment() {

       rv_Nsig        -> setVal( trandom_cls->Poisson( toymean_n_sig        ) ) ;
       rv_Nsb         -> setVal( trandom_cls->Poisson( toymean_n_sb         ) ) ;
       rv_Nsig_ldp    -> setVal( trandom_cls->Poisson( toymean_n_sig_ldp    ) ) ;
       rv_Nsb_ldp     -> setVal( trandom_cls->Poisson( toymean_n_sb_ldp     ) ) ;
       rv_Nsig_sl     -> setVal( trandom_cls->Poisson( toymean_n_sig_sl     ) ) ;
       rv_Nsb_sl      -> setVal( trandom_cls->Poisson( toymean_n_sb_sl      ) ) ;
       rv_Nlsb_0b     -> setVal( trandom_cls->Poisson( toymean_n_lsb_0b     ) ) ;
       rv_Nlsb_0b_ldp -> setVal( trandom_cls->Poisson( toymean_n_lsb_0b_ldp ) ) ;

       if ( znnModel == 1 ) {
          rv_Nsig_ee     -> setVal( trandom_cls->Poisson( toymean_n_sig_ee ) ) ;
          rv_Nsb_ee      -> setVal( trandom_cls->Poisson( toymean_n_sb_ee  ) ) ;
          rv_Nsig_mm     -> setVal( trandom_cls->Poisson( toymean_n_sig_mm ) ) ;
          rv_Nsb_mm      -> setVal( trandom_cls->Poisson( toymean_n_sb_mm  ) ) ;
       } else if ( znnModel == 2 ) {
          rv_Nsigsb_ee   -> setVal( trandom_cls->Poisson( toymean_n_sigsb_ee  ) ) ;
          rv_Nsigsb_mm   -> setVal( trandom_cls->Poisson( toymean_n_sigsb_mm  ) ) ;
       }


    } // genToyExperiment

  //===================================================================

    void ra2bRoostatsClass4ln::setObservablesToDataVals() {

       rv_Nsig        -> setVal(  dataval_Nsig         ) ;
       rv_Nsb         -> setVal(  dataval_Nsb          ) ;
       rv_Nsig_ldp    -> setVal(  dataval_Nsig_ldp     ) ;
       rv_Nsb_ldp     -> setVal(  dataval_Nsb_ldp      ) ;
       rv_Nsig_sl     -> setVal(  dataval_Nsig_sl      ) ;
       rv_Nsb_sl      -> setVal(  dataval_Nsb_sl       ) ;
       rv_Nlsb_0b     -> setVal(  dataval_Nlsb_0b      ) ;
       rv_Nlsb_0b_ldp -> setVal(  dataval_Nlsb_0b_ldp  ) ;

       if ( znnModel == 1 ) {
          rv_Nsig_ee     -> setVal(  dataval_Nsig_ee  ) ;
          rv_Nsb_ee      -> setVal(  dataval_Nsb_ee   ) ;
          rv_Nsig_mm     -> setVal(  dataval_Nsig_mm  ) ;
          rv_Nsb_mm      -> setVal(  dataval_Nsb_mm   ) ;
       } else if ( znnModel == 2 ) {
          rv_Nsigsb_ee   -> setVal(  dataval_Nsigsb_ee   ) ;
          rv_Nsigsb_mm   -> setVal(  dataval_Nsigsb_mm   ) ;
       }


    } // setObservablesToDataVals

  //===================================================================


    //-- returns p-value of data for this hypothesis, given the data-fit value of the test statistic q.

    double ra2bRoostatsClass4ln::doToyStudy( int nToys, bool isBgonlyStudy, double data_q ) {

       double retVal = 1.0 ;

       saveDataObservableVals() ;


       int n_sig         ;
       int n_sb          ;
       int n_sig_ldp     ;
       int n_sb_ldp      ;
       int n_sig_sl      ;
       int n_sb_sl       ;
       int n_lsb_0b      ;
       int n_lsb_0b_ldp  ;

       int n_sig_ee      ;
       int n_sb_ee       ;
       int n_sig_mm      ;
       int n_sb_mm       ;
       int n_sigsb_ee      ;
       int n_sigsb_mm      ;

       double tm_sig         ;
       double tm_sb          ;
       double tm_sig_ldp     ;
       double tm_sb_ldp      ;
       double tm_sig_sl      ;
       double tm_sb_sl       ;
       double tm_lsb_0b      ;
       double tm_lsb_0b_ldp  ;

       double tm_sig_ee      ;
       double tm_sb_ee       ;
       double tm_sig_mm      ;
       double tm_sb_mm       ;
       double tm_sigsb_ee      ;
       double tm_sigsb_mm      ;

       double maxLogL ;
       double sfixedLogL ;
       double testStat ;

       int maxCovQual ;
       int sfixedCovQual ;


       TTree* tt(0x0) ;

       if ( isBgonlyStudy ) {
          printf(" Creating bgonly ttree.\n") ;
          tt = new TTree("tt_cls_bgonly", "CLs BG-only toy study") ;
       } else {
          printf(" Creating splusb ttree.\n") ;
          tt = new TTree("tt_cls_splusb", "CLs SIG+BG toy study") ;
       }

   //-- To save disk quota space, need to strip the ntuple WAY down.

       tt->Branch( "maxLogL"        , &maxLogL           , "maxLogL/D"            ) ;
       tt->Branch( "sfixedLogL"        , &sfixedLogL           , "sfixedLogL/D"            ) ;
       tt->Branch( "testStat"        , &testStat           , "testStat/D"            ) ;

       tt->Branch( "maxCovQual"        , &maxCovQual           , "maxCovQual/I"            ) ;
       tt->Branch( "sfixedCovQual"        , &sfixedCovQual           , "sfixedCovQual/I"            ) ;

       tt->Branch( "n_sig"        , &n_sig           , "n_sig/I"            ) ;
       tt->Branch( "n_sb"         , &n_sb            , "n_sb/I"             ) ;
       tt->Branch( "n_sig_ldp"    , &n_sig_ldp       , "n_sig_ldp/I"        ) ;
       tt->Branch( "n_sb_ldp"     , &n_sb_ldp        , "n_sb_ldp/I"         ) ;
       tt->Branch( "n_sig_sl"     , &n_sig_sl        , "n_sig_sl/I"         ) ;
       tt->Branch( "n_sb_sl"      , &n_sb_sl         , "n_sb_sl/I"          ) ;
       tt->Branch( "n_lsb_0b"     , &n_lsb_0b        , "n_lsb_0b/I"         ) ;
       tt->Branch( "n_lsb_0b_ldp" , &n_lsb_0b_ldp    , "n_lsb_0b_ldp/I"     ) ;
       if ( znnModel==1 ) {
          tt->Branch( "n_sig_ee"        , &n_sig_ee           , "n_sig_ee/I"            ) ;
          tt->Branch( "n_sb_ee"         , &n_sb_ee            , "n_sb_ee/I"             ) ;
          tt->Branch( "n_sig_mm"        , &n_sig_mm           , "n_sig_mm/I"            ) ;
          tt->Branch( "n_sb_mm"         , &n_sb_mm            , "n_sb_mm/I"             ) ;
       } else if (znnModel==2) {
          tt->Branch( "n_sigsb_ee"        , &n_sigsb_ee           , "n_sigsb_ee/I"            ) ;
          tt->Branch( "n_sigsb_mm"        , &n_sigsb_mm           , "n_sigsb_mm/I"            ) ;
       }

       tt->Branch( "tm_sig"        , &tm_sig           , "tm_sig/D"            ) ;
       tt->Branch( "tm_sb"         , &tm_sb            , "tm_sb/D"             ) ;
       tt->Branch( "tm_sig_ldp"    , &tm_sig_ldp       , "tm_sig_ldp/D"        ) ;
       tt->Branch( "tm_sb_ldp"     , &tm_sb_ldp        , "tm_sb_ldp/D"         ) ;
       tt->Branch( "tm_sig_sl"     , &tm_sig_sl        , "tm_sig_sl/D"         ) ;
       tt->Branch( "tm_sb_sl"      , &tm_sb_sl         , "tm_sb_sl/D"          ) ;
       tt->Branch( "tm_lsb_0b"     , &tm_lsb_0b        , "tm_lsb_0b/D"         ) ;
       tt->Branch( "tm_lsb_0b_ldp" , &tm_lsb_0b_ldp    , "tm_lsb_0b_ldp/D"     ) ;
       if ( znnModel==1 ) {
          tt->Branch( "tm_sig_ee"        , &tm_sig_ee           , "tm_sig_ee/D"            ) ;
          tt->Branch( "tm_sb_ee"         , &tm_sb_ee            , "tm_sb_ee/D"             ) ;
          tt->Branch( "tm_sig_mm"        , &tm_sig_mm           , "tm_sig_mm/D"            ) ;
          tt->Branch( "tm_sb_mm"         , &tm_sb_mm            , "tm_sb_mm/D"             ) ;
       } else if (znnModel==2) {
          tt->Branch( "tm_sigsb_ee"        , &tm_sigsb_ee           , "tm_sigsb_ee/D"            ) ;
          tt->Branch( "tm_sigsb_mm"        , &tm_sigsb_mm           , "tm_sigsb_mm/D"            ) ;
       }


       double original_mean_eff_sf_sig     = rv_mean_eff_sf_sig -> getVal() ;
       double original_mean_eff_sf_sb      = rv_mean_eff_sf_sb -> getVal() ;
       double original_mean_eff_sf_sig_sl  = rv_mean_eff_sf_sig_sl -> getVal() ;
       double original_mean_eff_sf_sb_sl   = rv_mean_eff_sf_sb_sl -> getVal() ;
       double original_mean_eff_sf_sig_ldp = rv_mean_eff_sf_sig_ldp -> getVal() ;
       double original_mean_eff_sf_sb_ldp  = rv_mean_eff_sf_sb_ldp -> getVal() ;

       double original_np_m_sf_mc          = np_m_sf_mc  -> getVal() ;
       double original_np_m_sf_qcd_sb      = np_m_sf_qcd_sb  -> getVal() ;
       double original_np_m_sf_qcd_sig     = np_m_sf_qcd_sig  -> getVal() ;
       double original_np_m_sf_ttwj_sig    = np_m_sf_ttwj_sig  -> getVal() ;
       double original_np_m_sf_ee          = np_m_sf_ee  -> getVal() ;
       double original_np_m_sf_mm          = np_m_sf_mm  -> getVal() ;
       double original_np_m_acc_ee         = np_m_acc_ee  -> getVal() ;
       double original_np_m_acc_mm         = np_m_acc_mm  -> getVal() ;
       double original_np_m_eff_ee         = np_m_eff_ee  -> getVal() ;
       double original_np_m_eff_mm         = np_m_eff_mm  -> getVal() ;
       double original_np_m_fsig_ee        = np_m_fsig_ee  -> getVal() ;
       double original_np_m_fsig_mm        = np_m_fsig_mm  -> getVal() ;
       double original_np_m_knn_sig(0.) ;
       double original_np_m_knn_sb(0.) ;
       if ( znnModel==2 ) {
          original_np_m_knn_sig        = np_m_knn_sig  -> getVal() ;
          original_np_m_knn_sb         = np_m_knn_sb  -> getVal() ;
       }


       //-- memory management debugging
       ////// RooTrace::active(kTRUE) ;

       int nWorse(0) ;

       for ( int ti=0; ti<nToys; ti++ ) {

          RooFitResult* toyFitResult(0x0) ;

          printf("\n\n\n\n\n +++++++++++++++++ Beginning of toy experiment %6d +++++++++++++++++\n\n\n", ti ) ;

        //-- memory management debugging
          ////// RooTrace::mark() ;

         //-- Steps 1,2) Generate value of mean for each Gaussian PDF (nuisance parameter "measurements").
         //              Build a new likelihood with these Gaussian PDFs.




          double grn ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_mean_eff_sf_sig, rv_width_eff_sf_sig->getVal() ) ; }
          rv_mean_eff_sf_sig -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_mean_eff_sf_sb, rv_width_eff_sf_sb->getVal() ) ; }
          rv_mean_eff_sf_sb -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_mean_eff_sf_sig_sl, rv_width_eff_sf_sig_sl->getVal() ) ; }
          rv_mean_eff_sf_sig_sl -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_mean_eff_sf_sb_sl, rv_width_eff_sf_sb_sl->getVal() ) ; }
          rv_mean_eff_sf_sb_sl -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_mean_eff_sf_sig_ldp, rv_width_eff_sf_sig_ldp->getVal() ) ; }
          rv_mean_eff_sf_sig_ldp -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_mean_eff_sf_sb_ldp, rv_width_eff_sf_sb_ldp->getVal() ) ; }
          rv_mean_eff_sf_sb_ldp -> setVal( grn ) ;




          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_sf_mc, np_s_sf_mc->getVal() ) ; }
          np_m_sf_mc -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_sf_qcd_sb, np_s_sf_qcd_sb->getVal() ) ; }
          np_m_sf_qcd_sb -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_sf_qcd_sig, np_s_sf_qcd_sig->getVal() ) ; }
          np_m_sf_qcd_sig -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_sf_ttwj_sig, np_s_sf_ttwj_sig->getVal() ) ; }
          np_m_sf_ttwj_sig -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_sf_ee, np_s_sf_ee->getVal() ) ; }
          np_m_sf_ee -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_sf_mm, np_s_sf_mm->getVal() ) ; }
          np_m_sf_mm -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_acc_ee, np_s_acc_ee->getVal() ) ; }
          np_m_acc_ee -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_acc_mm, np_s_acc_mm->getVal() ) ; }
          np_m_acc_mm -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_eff_ee, np_s_eff_ee->getVal() ) ; }
          np_m_eff_ee -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_eff_mm, np_s_eff_mm->getVal() ) ; }
          np_m_eff_mm -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_fsig_ee, np_s_fsig_ee->getVal() ) ; }
          np_m_fsig_ee -> setVal( grn ) ;

          grn = -1. ;
          while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_fsig_mm, np_s_fsig_mm->getVal() ) ; }
          np_m_fsig_mm -> setVal( grn ) ;

          if ( znnModel == 2 ) {

             grn = -1. ;
             while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_knn_sig, np_s_knn_sig->getVal() ) ; }
             np_m_knn_sig -> setVal( grn ) ;

             grn = -1. ;
             while ( grn < 0. ) { grn = trandom_cls->Gaus( original_np_m_knn_sb, np_s_knn_sb->getVal() ) ; }
             np_m_knn_sb -> setVal( grn ) ;

          }









       //--- Step 3) Fit the actual data with this likelihood function.

          reinitialize() ;

          setObservablesToDataVals() ;

          RooArgSet toyGenFitobservedParametersList ;
          toyGenFitobservedParametersList.add( *rv_Nsig        ) ;
          toyGenFitobservedParametersList.add( *rv_Nsb         ) ;
          toyGenFitobservedParametersList.add( *rv_Nsig_sl     ) ;
          toyGenFitobservedParametersList.add( *rv_Nsb_sl      ) ;
          toyGenFitobservedParametersList.add( *rv_Nsig_ldp    ) ;
          toyGenFitobservedParametersList.add( *rv_Nsb_ldp     ) ;
          toyGenFitobservedParametersList.add( *rv_Nlsb_0b     ) ;
          toyGenFitobservedParametersList.add( *rv_Nlsb_0b_ldp ) ;
          if ( znnModel == 1 ) {
             toyGenFitobservedParametersList.add( *rv_Nsb_ee      ) ;
             toyGenFitobservedParametersList.add( *rv_Nsig_ee     ) ;
             toyGenFitobservedParametersList.add( *rv_Nsb_mm      ) ;
             toyGenFitobservedParametersList.add( *rv_Nsig_mm     ) ;
          } else if ( znnModel == 2 ) {
             toyGenFitobservedParametersList.add( *rv_Nsigsb_ee      ) ;
             toyGenFitobservedParametersList.add( *rv_Nsigsb_mm      ) ;
          }



          RooDataSet* toyGenFitdsObserved = new RooDataSet("toyGenFit_ra2b_observed_rds", "ToyGenL RA2b observed data values",
                                         toyGenFitobservedParametersList ) ;
          toyGenFitdsObserved->add( toyGenFitobservedParametersList ) ;


          if ( isBgonlyStudy ) {
            //-- fit with susy yield fixed to zero.
             rv_mu_susy_sig -> setVal( 0. ) ;
          } else {
            //-- fit with susy yield fixed to predicted value.
             rv_mu_susy_sig -> setVal( rv_mu_susymc_sig->getVal() ) ;
          }
          rv_mu_susy_sig->setConstant( kTRUE ) ;

          printf("\n\n") ;
          printf("  Fitting with these values for the observables.\n") ;
          toyGenFitdsObserved->printMultiline(cout, 1, kTRUE, "") ;
          printf("\n\n") ;
       // toyFitResult = tgl_likelihood->fitTo(*toyGenFitdsObserved, Save(true));
  /////// toyFitResult = tgl_likelihood->fitTo(*toyGenFitdsObserved, Save(true), PrintLevel(-1));
          toyFitResult = likelihood->fitTo(*toyGenFitdsObserved, Save(true), PrintLevel(-1));
          maxLogL = toyFitResult->minNll() ;
          maxCovQual = toyFitResult->covQual() ;
          if ( isBgonlyStudy ) {
             printf("\n  Fit result with susy fixed to zero for toy %d : %8.3f \n\n", ti, maxLogL ) ;
          } else {
             printf("\n  Fit result with susy fixed to predicted value for toy %d : %8.3f \n\n", ti, maxLogL ) ;
          }
          parameterSnapshot() ;
          delete toyFitResult ;



       //--- Step 4) Use the results of this fit to generate a toy dataset.

          saveToymeanSnapshot() ;
          genToyExperiment() ;


       //--- Step 5) Do the analysis (i.e. evaluate the test statistic) for this toy dataset using
       //            the nominal likelihood (NOT the one for the fit above).

          rv_mean_eff_sf_sig     -> setVal(  original_mean_eff_sf_sig     ) ;
          rv_mean_eff_sf_sb      -> setVal(  original_mean_eff_sf_sb      ) ;
          rv_mean_eff_sf_sig_sl  -> setVal(  original_mean_eff_sf_sig_sl  ) ;
          rv_mean_eff_sf_sb_sl   -> setVal(  original_mean_eff_sf_sb_sl   ) ;
          rv_mean_eff_sf_sig_ldp -> setVal(  original_mean_eff_sf_sig_ldp ) ;
          rv_mean_eff_sf_sb_ldp  -> setVal(  original_mean_eff_sf_sb_ldp  ) ;


          np_m_sf_mc        -> setVal( original_np_m_sf_mc       ) ;
          np_m_sf_qcd_sb    -> setVal( original_np_m_sf_qcd_sb   ) ;
          np_m_sf_qcd_sig   -> setVal( original_np_m_sf_qcd_sig  ) ;
          np_m_sf_ttwj_sig  -> setVal( original_np_m_sf_ttwj_sig ) ;
          np_m_sf_ee        -> setVal( original_np_m_sf_ee       ) ;
          np_m_sf_mm        -> setVal( original_np_m_sf_mm       ) ;
          np_m_acc_ee       -> setVal( original_np_m_acc_ee      ) ;
          np_m_acc_mm       -> setVal( original_np_m_acc_mm      ) ;
          np_m_eff_ee       -> setVal( original_np_m_eff_ee      ) ;
          np_m_eff_mm       -> setVal( original_np_m_eff_mm      ) ;
          np_m_fsig_ee      -> setVal( original_np_m_fsig_ee     ) ;
          np_m_fsig_mm      -> setVal( original_np_m_fsig_mm     ) ;
          if ( znnModel == 2 ) {
             np_m_knn_sig      -> setVal( original_np_m_knn_sig     ) ;
             np_m_knn_sb       -> setVal( original_np_m_knn_sb      ) ;
          }



          RooArgSet toyFitobservedParametersList ;
          toyFitobservedParametersList.add( *rv_Nsig        ) ;
          toyFitobservedParametersList.add( *rv_Nsb         ) ;
          toyFitobservedParametersList.add( *rv_Nsig_sl     ) ;
          toyFitobservedParametersList.add( *rv_Nsb_sl      ) ;
          toyFitobservedParametersList.add( *rv_Nsig_ldp    ) ;
          toyFitobservedParametersList.add( *rv_Nsb_ldp     ) ;
          toyFitobservedParametersList.add( *rv_Nlsb_0b     ) ;
          toyFitobservedParametersList.add( *rv_Nlsb_0b_ldp ) ;
          if ( znnModel == 1 ) {
             toyFitobservedParametersList.add( *rv_Nsb_ee      ) ;
             toyFitobservedParametersList.add( *rv_Nsig_ee     ) ;
             toyFitobservedParametersList.add( *rv_Nsb_mm      ) ;
             toyFitobservedParametersList.add( *rv_Nsig_mm     ) ;
          } else if ( znnModel == 2 ) {
             toyFitobservedParametersList.add( *rv_Nsigsb_ee      ) ;
             toyFitobservedParametersList.add( *rv_Nsigsb_mm      ) ;
          }


          RooDataSet* toyFitdsObserved = new RooDataSet("toyfit_ra2b_observed_rds", "RA2b toy observed data values",
                                         toyFitobservedParametersList ) ;
          toyFitdsObserved->add( toyFitobservedParametersList ) ;


         //-- fit with susy yield floating to get the absolute maximum log likelihood.

          rv_mu_susy_sig->setConstant( kFALSE ) ;

          printf("\n\n") ;
          printf("  Fitting with these values for the observables.\n") ;
          toyFitdsObserved->printMultiline(cout, 1, kTRUE, "") ;
          printf("\n\n") ;
       // toyFitResult = likelihood->fitTo(*toyFitdsObserved, Save(true));
          toyFitResult = likelihood->fitTo(*toyFitdsObserved, Save(true), PrintLevel(-1));
          double maxL_mu_susy_sig = rv_mu_susy_sig->getVal() ;
          maxLogL = toyFitResult->minNll() ;
          maxCovQual = toyFitResult->covQual() ;
          printf("\n  Fit result with susy floating for toy %d : %8.3f \n\n", ti, maxLogL ) ;
          parameterSnapshot() ;
          delete toyFitResult ;


        //--- Only bother doing the next step if the maxL value of the susy yield is
        //    less than the predicted value.
        //    If it the maxL value is greater than the predicted value, set the test statistic to zero.
        //    This is to get a one-sided limit.

          if ( maxL_mu_susy_sig < (rv_mu_susymc_sig->getVal()) ) {
            //-- fit with susy yield fixed to susy model prediction.
             rv_mu_susy_sig -> setVal( rv_mu_susymc_sig->getVal() ) ;
             rv_mu_susy_sig->setConstant( kTRUE ) ;
             printf("\n  Fitting with susy SIG yield fixed to %8.2f\n\n", rv_mu_susy_sig->getVal() ) ;
          // toyFitResult = likelihood->fitTo(*toyFitdsObserved, Save(true));
             toyFitResult = likelihood->fitTo(*toyFitdsObserved, Save(true), PrintLevel(-1));
             sfixedLogL = toyFitResult->minNll() ;
             sfixedCovQual = toyFitResult->covQual() ;
             testStat = 2.*(sfixedLogL-maxLogL) ;
             printf("\n  Fit result with susy fixed to prediction for toy %d : %8.3f, %8.3f \n\n", ti, sfixedLogL, testStat ) ;
             parameterSnapshot() ;
             delete toyFitResult ;
          } else {
             sfixedLogL = 0.0 ;
             sfixedCovQual = 0 ;
             testStat = 0.0 ;
             printf("\n  SUSY yield from susy floating is greater than predicted value: %9.2f > %9.2f\n",
                          maxL_mu_susy_sig, (rv_mu_susymc_sig->getVal()) ) ;
             printf("  +++++ Setting test statistic to zero.\n\n") ;
          }








         //--- save this experiment in the TTree.

          n_sig        = rv_Nsig        ->getVal() ;
          n_sb         = rv_Nsb         ->getVal() ;
          n_sig_sl     = rv_Nsig_sl     ->getVal() ;
          n_sb_sl      = rv_Nsb_sl      ->getVal() ;
          n_sig_ldp    = rv_Nsig_ldp    ->getVal() ;
          n_sb_ldp     = rv_Nsb_ldp     ->getVal() ;
          n_lsb_0b     = rv_Nlsb_0b     ->getVal() ;
          n_lsb_0b_ldp = rv_Nlsb_0b_ldp ->getVal() ;

          if ( znnModel == 1 ) {
             n_sig_ee = rv_Nsig_ee ->getVal() ;
             n_sb_ee  = rv_Nsb_ee  ->getVal() ;
             n_sig_mm = rv_Nsig_mm ->getVal() ;
             n_sb_mm  = rv_Nsb_mm  ->getVal() ;
          } else if ( znnModel == 2 ) {
             n_sigsb_ee = rv_Nsigsb_ee ->getVal() ;
             n_sigsb_mm = rv_Nsigsb_mm ->getVal() ;
          }

          tm_sig        = toymean_n_sig ;
          tm_sb         = toymean_n_sb ;
          tm_sig_sl     = toymean_n_sig_sl ;
          tm_sb_sl      = toymean_n_sb_sl ;
          tm_sig_ldp    = toymean_n_sig_ldp ;
          tm_sb_ldp     = toymean_n_sb_ldp ;
          tm_lsb_0b     = toymean_n_lsb_0b ;
          tm_lsb_0b_ldp = toymean_n_lsb_0b_ldp ;
          if ( znnModel == 1 ) {
             tm_sig_ee   = toymean_n_sig_ee ;
             tm_sb_ee    = toymean_n_sb_ee ;
             tm_sig_mm   = toymean_n_sig_mm ;
             tm_sb_mm    = toymean_n_sb_mm ;
          } else if ( znnModel == 2 ) {
             tm_sigsb_ee = toymean_n_sigsb_ee ;
             tm_sigsb_mm = toymean_n_sigsb_mm ;
          }

          //-- check for Nans
          if ( TMath::IsNaN( testStat ) != 0 ) { testStat = -1. ; }
          if ( TMath::IsNaN( maxLogL ) != 0 ) { maxLogL = -1. ; }
          if ( TMath::IsNaN( sfixedLogL ) != 0 ) { sfixedLogL = -1. ; }

          tt->Fill() ;

          if ( testStat >= data_q ) { nWorse++ ; }

          if ( isBgonlyStudy ) {
             printf("\n\n ===== Value of test statistic for BG-only  toy %5d is %9.3f\n\n", ti, testStat) ;
          } else {
             printf("\n\n ===== Value of test statistic for S-plus-B toy %5d is %9.3f\n\n", ti, testStat) ;
          }


          //-- clean up. (This doesn't seem to help with the memory leak...)

          printf("\n\n  Before deletes.\n\n") ;
          ////// RooTrace::dump(cout,kTRUE) ;


          printf(" Cleaning up...\n") ;
          toyGenFitobservedParametersList.removeAll() ;
          toyFitobservedParametersList.removeAll() ;
          delete toyGenFitdsObserved ;
          delete toyFitdsObserved ;
          printf(" Done cleaning up.\n") ;

        //--- memory management debugging
          printf("\n\n  After deletes.\n\n") ;
          ////// RooTrace::dump(cout,kTRUE) ;

       } // ti.


  //-- turn off root output to save disk quota.
       tt->Write() ;

       retVal = (1.0*nWorse)/(1.0*nToys) ;

       printf("\n\n\n ra2bRoostatsClass4ln :  p-value for data_q = %8.3f is  %d / %d = %8.3f\n\n\n",
            data_q, nWorse, nToys, retVal ) ;

       return retVal ;


    } // doToyStudy


  //===================================================================


    double ra2bRoostatsClass4ln::doToyStudyNoSusyInFit( int nToys, const char* trueValsInputFile ) {

       float true_ttwj_sig(0.) ;
       float true_qcd_sig(0.) ;
       float true_znn_sig(0.) ;

       if ( strlen(trueValsInputFile) > 0 ) {
          FILE* infp ;
          if ( (infp=fopen( trueValsInputFile,"r"))==NULL ) {
             printf("\n\n *** Problem opening true values input file: %s.\n\n", trueValsInputFile ) ;
             return 0.0 ;
          }
          char label[1000] ;
          //--- order matters!
          fscanf( infp, "%s %g", label, &true_ttwj_sig ) ;
          fscanf( infp, "%s %g", label, &true_qcd_sig ) ;
          fscanf( infp, "%s %g", label, &true_znn_sig ) ;
          printf("\n\n True values of SIG yields: ttwj =%5.1f, QCD =%5.1f, Znn =%5.1f\n\n", true_ttwj_sig, true_qcd_sig, true_znn_sig ) ;
          fclose( infp ) ;
       }

       double retVal = 1.0 ;

       int n_sig         ;
       int n_sb          ;
       int n_sig_ldp     ;
       int n_sb_ldp      ;
       int n_sig_sl      ;
       int n_sb_sl       ;
       int n_lsb_0b      ;
       int n_lsb_0b_ldp  ;

       int n_sig_ee      ;
       int n_sb_ee       ;
       int n_sig_mm      ;
       int n_sb_mm       ;
       int n_sigsb_ee      ;
       int n_sigsb_mm      ;

       double ttwj_sig_true ;
       double ttwj_sig_fit ;
       double ttwj_sig_err ;
       double ttwj_sig_pull ;

       double qcd_sig_true ;
       double qcd_sig_fit ;
       double qcd_sig_err ;
       double qcd_sig_pull ;

       double znn_sig_true ;
       double znn_sig_fit ;
       double znn_sig_err ;
       double znn_sig_pull ;

       int maxCovQual ;


       TTree* tt(0x0) ;

       tt = new TTree("tt_toy_nosusyfit", "toy study for fit with SUSY fixed to zero") ;

       tt->Branch( "maxCovQual"        , &maxCovQual           , "maxCovQual/I"            ) ;

       tt->Branch( "n_sig"        , &n_sig           , "n_sig/I"            ) ;
       tt->Branch( "n_sb"         , &n_sb            , "n_sb/I"             ) ;
       tt->Branch( "n_sig_ldp"    , &n_sig_ldp       , "n_sig_ldp/I"        ) ;
       tt->Branch( "n_sb_ldp"     , &n_sb_ldp        , "n_sb_ldp/I"         ) ;
       tt->Branch( "n_sig_sl"     , &n_sig_sl        , "n_sig_sl/I"         ) ;
       tt->Branch( "n_sb_sl"      , &n_sb_sl         , "n_sb_sl/I"          ) ;
       tt->Branch( "n_lsb_0b"     , &n_lsb_0b        , "n_lsb_0b/I"         ) ;
       tt->Branch( "n_lsb_0b_ldp" , &n_lsb_0b_ldp    , "n_lsb_0b_ldp/I"     ) ;
       if ( znnModel==1 ) {
          tt->Branch( "n_sig_ee"        , &n_sig_ee           , "n_sig_ee/I"            ) ;
          tt->Branch( "n_sb_ee"         , &n_sb_ee            , "n_sb_ee/I"             ) ;
          tt->Branch( "n_sig_mm"        , &n_sig_mm           , "n_sig_mm/I"            ) ;
          tt->Branch( "n_sb_mm"         , &n_sb_mm            , "n_sb_mm/I"             ) ;
       } else if (znnModel==2) {
          tt->Branch( "n_sigsb_ee"        , &n_sigsb_ee           , "n_sigsb_ee/I"            ) ;
          tt->Branch( "n_sigsb_mm"        , &n_sigsb_mm           , "n_sigsb_mm/I"            ) ;
       }

       tt->Branch( "ttwj_sig_true"        , &ttwj_sig_true           , "ttwj_sig_true/D"            ) ;
       tt->Branch( "ttwj_sig_fit"         , &ttwj_sig_fit            , "ttwj_sig_fit/D"             ) ;
       tt->Branch( "ttwj_sig_err"         , &ttwj_sig_err            , "ttwj_sig_err/D"             ) ;
       tt->Branch( "ttwj_sig_pull"        , &ttwj_sig_pull           , "ttwj_sig_pull/D"            ) ;

       tt->Branch( "qcd_sig_true"        , &qcd_sig_true           , "qcd_sig_true/D"            ) ;
       tt->Branch( "qcd_sig_fit"         , &qcd_sig_fit            , "qcd_sig_fit/D"             ) ;
       tt->Branch( "qcd_sig_err"         , &qcd_sig_err            , "qcd_sig_err/D"             ) ;
       tt->Branch( "qcd_sig_pull"        , &qcd_sig_pull           , "qcd_sig_pull/D"            ) ;

       tt->Branch( "znn_sig_true"        , &znn_sig_true           , "znn_sig_true/D"            ) ;
       tt->Branch( "znn_sig_fit"         , &znn_sig_fit            , "znn_sig_fit/D"             ) ;
       tt->Branch( "znn_sig_err"         , &znn_sig_err            , "znn_sig_err/D"             ) ;
       tt->Branch( "znn_sig_pull"        , &znn_sig_pull           , "znn_sig_pull/D"            ) ;


       ttwj_sig_true = true_ttwj_sig ;
       qcd_sig_true  = true_qcd_sig ;
       znn_sig_true  = true_znn_sig ;


       for ( int ti=0; ti<nToys; ti++ ) {

          reinitialize() ;

          genToyExperiment() ;

          RooArgSet toyFitobservedParametersList ;
          toyFitobservedParametersList.add( *rv_Nsig        ) ;
          toyFitobservedParametersList.add( *rv_Nsb         ) ;
          toyFitobservedParametersList.add( *rv_Nsig_sl     ) ;
          toyFitobservedParametersList.add( *rv_Nsb_sl      ) ;
          toyFitobservedParametersList.add( *rv_Nsig_ldp    ) ;
          toyFitobservedParametersList.add( *rv_Nsb_ldp     ) ;
          toyFitobservedParametersList.add( *rv_Nlsb_0b     ) ;
          toyFitobservedParametersList.add( *rv_Nlsb_0b_ldp ) ;
          if ( znnModel == 1 ) {
             toyFitobservedParametersList.add( *rv_Nsb_ee      ) ;
             toyFitobservedParametersList.add( *rv_Nsig_ee     ) ;
             toyFitobservedParametersList.add( *rv_Nsb_mm      ) ;
             toyFitobservedParametersList.add( *rv_Nsig_mm     ) ;
          } else if ( znnModel == 2 ) {
             toyFitobservedParametersList.add( *rv_Nsigsb_ee      ) ;
             toyFitobservedParametersList.add( *rv_Nsigsb_mm      ) ;
          }


          RooDataSet* toyFitdsObserved = new RooDataSet("toyfit_ra2b_observed_rds", "RA2b toy observed data values",
                                         toyFitobservedParametersList ) ;
          toyFitdsObserved->add( toyFitobservedParametersList ) ;


         //-- fit with susy yield fixed to zero.

          rv_mu_susy_sig -> setVal( 0.0 ) ;
          rv_mu_susy_sig->setConstant( kTRUE ) ;

          printf("\n\n") ;
          printf("  Fitting with these values for the observables.\n") ;
          toyFitdsObserved->printMultiline(cout, 1, kTRUE, "") ;
          printf("\n\n") ;
       // fitResult = likelihood->fitTo(*toyFitdsObserved, Save(true));
          fitResult = likelihood->fitTo(*toyFitdsObserved, Save(true), PrintLevel(-1));
          maxCovQual = fitResult->covQual() ;
          parameterSnapshot() ;


         //--- save this experiment in the TTree.

          n_sig        = rv_Nsig        ->getVal() ;
          n_sb         = rv_Nsb         ->getVal() ;
          n_sig_sl     = rv_Nsig_sl     ->getVal() ;
          n_sb_sl      = rv_Nsb_sl      ->getVal() ;
          n_sig_ldp    = rv_Nsig_ldp    ->getVal() ;
          n_sb_ldp     = rv_Nsb_ldp     ->getVal() ;
          n_lsb_0b     = rv_Nlsb_0b     ->getVal() ;
          n_lsb_0b_ldp = rv_Nlsb_0b_ldp ->getVal() ;

          if ( znnModel == 1 ) {
             n_sig_ee = rv_Nsig_ee ->getVal() ;
             n_sb_ee  = rv_Nsb_ee  ->getVal() ;
             n_sig_mm = rv_Nsig_mm ->getVal() ;
             n_sb_mm  = rv_Nsb_mm  ->getVal() ;
          } else if ( znnModel == 2 ) {
             n_sig_ee = rv_Nsigsb_ee ->getVal() ;
             n_sig_ee = rv_Nsigsb_mm ->getVal() ;
          }


          if ( useSigTtwjVar ) {
             ttwj_sig_fit = rrv_mu_ttwj_sig->getVal() ;
             ttwj_sig_err = rrv_mu_ttwj_sig->getError() ;
             ttwj_sig_pull = -12. ;
             if ( ttwj_sig_err > 0. ) { ttwj_sig_pull = (ttwj_sig_fit - ttwj_sig_true)/ttwj_sig_err ; }
          } else {
             ttwj_sig_fit = rfv_mu_ttwj_sig->getVal() ;
             ttwj_sig_err = -1. ;
             ttwj_sig_pull = -12. ;
          }

          if ( !useLdpVars ) {
             qcd_sig_fit = rrv_mu_qcd_sig->getVal() ;
             qcd_sig_err = rrv_mu_qcd_sig->getError() ;
             qcd_sig_pull = -12. ;
             if ( qcd_sig_err > 0. ) { qcd_sig_pull = (qcd_sig_fit - qcd_sig_true)/qcd_sig_err ; }
          } else {
             qcd_sig_fit = rfv_mu_qcd_sig->getVal() ;
             qcd_sig_err = -1. ;
             qcd_sig_pull = -12. ;
          }

          znn_sig_fit = rv_mu_znn_sig->getVal() ;
          znn_sig_err = rv_mu_znn_sig->getError() ;
          znn_sig_pull = -12. ;
          if ( znn_sig_err > 0. ) { znn_sig_pull = (znn_sig_fit - znn_sig_true)/znn_sig_err ; }


          tt->Fill() ;

       } // ti.


       tt->Write() ;

       return retVal ;


    } // doToyStudyNoSusyInFit


  //===================================================================


    double ra2bRoostatsClass4ln::getLogLikelihoodValue( ) {

        return fitResult->minNll() ;

    }

  //===================================================================














