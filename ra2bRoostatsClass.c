
//
//   Owen Long, UCR
//   Harrison Prosper, FSU
//   Sezen Sekmen, FSU
//
//


#include "ra2bRoostatsClass.h"

#include <iostream>


#include "TCanvas.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TText.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"

#include "RooArgSet.h"
#include "RooConstVar.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

  using namespace RooFit ;
  using namespace RooStats ;


  //=====================================================================================================


   ra2bRoostatsClass::ra2bRoostatsClass( bool ArgUseSigBgVars ) {

      gStyle->SetOptStat(0) ;

      useSigBgVars = ArgUseSigBgVars ;

     //--- Tell RooFit to shut up about anything less important than an ERROR.
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
      printf("\n\n ==== RooFit output configuration ===============\n") ;
      RooMsgService::instance().Print() ;
      printf("\n\n ================================================\n") ;

      varsAtFitVals = false ;
      initialized = false ;

   }




  //==================================================================================================

    bool ra2bRoostatsClass::generateToyDatasetsFromLikelihood( const char* outputRootFile, int nToys ) {

       if ( ! initialized ) {
          printf("\n\n *** Call initialize first.\n\n") ;
          return false ;
       }

       TTree* toyOutputTrueTree = new TTree("toytruetree", "ra2b toy mean true value tree" ) ;

       toyOutputTrueTree->Branch("mu0_ttbar_sig"      , &toy_mu0_ttbar_sig, "mu0_ttbar_sig/D" ) ;
       toyOutputTrueTree->Branch("mu0_qcd_sig"        , &toy_mu0_qcd_sig, "mu0_qcd_sig/D" ) ;
       toyOutputTrueTree->Branch("mu0_ttbar_sb"       , &toy_mu0_ttbar_sb, "mu0_ttbar_sb/D" ) ;
       toyOutputTrueTree->Branch("mu0_qcd_sb"         , &toy_mu0_qcd_sb, "mu0_qcd_sb/D" ) ;
       toyOutputTrueTree->Branch("mu0_susy_sig"       , &toy_mu0_susy_sig       , "mu0_susy_sig/D" ) ;
       toyOutputTrueTree->Branch("mu0_allbg_sig"      , &toy_mu0_allbg_sig      , "mu0_allbg_sig/D" ) ;


      //--- Set the true values.
       if ( useSigBgVars ) {
          toy_mu0_ttbar_sig = rrv_mu_ttbar_sig->getVal() ;
          toy_mu0_qcd_sig   = rrv_mu_qcd_sig->getVal() ;
          toy_mu0_ttbar_sb  = rfv_mu_ttbar_sb->getVal() ;
          toy_mu0_qcd_sb    = rfv_mu_qcd_sb->getVal() ;
       } else {
          toy_mu0_ttbar_sig = rfv_mu_ttbar_sig->getVal() ;
          toy_mu0_qcd_sig   = rfv_mu_qcd_sig->getVal() ;
          toy_mu0_ttbar_sb  = rrv_mu_ttbar_sb->getVal() ;
          toy_mu0_qcd_sb    = rrv_mu_qcd_sb->getVal() ;
       }
       toy_mu0_susy_sig = rv_mu_susy_sig->getVal() ;
       toy_mu0_allbg_sig = rv_mu_ew_sig->getVal() + toy_mu0_ttbar_sig + toy_mu0_qcd_sig ;
       printf("\n\n") ;
       printf("  Mean true values used in toy generation:\n") ;
       printf("  mu0_ttbar_sig : %6.1f\n", toy_mu0_ttbar_sig ) ;
       printf("  mu0_qcd_sig   : %6.1f\n", toy_mu0_qcd_sig ) ;
       printf("  mu0_ttbar_sb  : %6.1f\n", toy_mu0_ttbar_sb ) ;
       printf("  mu0_qcd_sb    : %6.1f\n", toy_mu0_qcd_sb ) ;
       printf("  mu0_susy_sig  : %6.1f\n", toy_mu0_susy_sig ) ;
       printf("  mu0_allbg_sig : %6.1f\n", toy_mu0_allbg_sig ) ;

     //--- generate some toy datasets.
       printf("\n\n Opening output toy datasets root file : %s\n", outputRootFile ) ;
       TFile toyOutputFile( outputRootFile, "recreate" ) ;
       RooDataSet* toyDatasets = likelihood->generate( observedParametersList, nToys ) ;
       const TTree* toyTree = toyDatasets->tree() ;
       toyTree->Write() ;
       toyOutputTrueTree->Fill() ;
       toyOutputTrueTree->Write() ;
       printf("\n\n Closing output toy datasets root file : %s\n", outputRootFile ) ;
       toyOutputFile.Close() ;

       return true ;

    } // generateToyDatasetsFromLikelihood

  //==================================================================================================

    bool ra2bRoostatsClass::generateToyDatasetsFromNVals( const char* outputRootFile, int nToys ) {

       if ( ! initialized ) {
          printf("\n\n *** Call initialize first.\n\n") ;
          return false ;
       }

       TTree* toyOutputTrueTree = new TTree("toytruetree", "ra2b toy mean true value tree" ) ;

       toyOutputTrueTree->Branch("mu0_ttbar_sig"      , &toy_mu0_ttbar_sig, "mu0_ttbar_sig/D" ) ;
       toyOutputTrueTree->Branch("mu0_qcd_sig"        , &toy_mu0_qcd_sig, "mu0_qcd_sig/D" ) ;
       toyOutputTrueTree->Branch("mu0_ttbar_sb"       , &toy_mu0_ttbar_sb, "mu0_ttbar_sb/D" ) ;
       toyOutputTrueTree->Branch("mu0_qcd_sb"         , &toy_mu0_qcd_sb, "mu0_qcd_sb/D" ) ;
       toyOutputTrueTree->Branch("mu0_susy_sig"       , &toy_mu0_susy_sig       , "mu0_susy_sig/D" ) ;
       toyOutputTrueTree->Branch("mu0_allbg_sig"      , &toy_mu0_allbg_sig      , "mu0_allbg_sig/D" ) ;


      //--- Set the true values.
       if ( useSigBgVars ) {
          toy_mu0_ttbar_sig = rrv_mu_ttbar_sig->getVal() ;
          toy_mu0_qcd_sig   = rrv_mu_qcd_sig->getVal() ;
          toy_mu0_ttbar_sb  = rfv_mu_ttbar_sb->getVal() ;
          toy_mu0_qcd_sb    = rfv_mu_qcd_sb->getVal() ;
       } else {
          toy_mu0_ttbar_sig = rfv_mu_ttbar_sig->getVal() ;
          toy_mu0_qcd_sig   = rfv_mu_qcd_sig->getVal() ;
          toy_mu0_ttbar_sb  = rrv_mu_ttbar_sb->getVal() ;
          toy_mu0_qcd_sb    = rrv_mu_qcd_sb->getVal() ;
       }
       toy_mu0_susy_sig = rv_mu_susy_sig->getVal() ;
       toy_mu0_allbg_sig = rv_mu_ew_sig->getVal() + toy_mu0_ttbar_sig + toy_mu0_qcd_sig ;
       printf("\n\n") ;
       printf("  Mean true values used in toy generation:\n") ;
       printf("  mu0_ttbar_sig : %6.1f\n", toy_mu0_ttbar_sig ) ;
       printf("  mu0_qcd_sig   : %6.1f\n", toy_mu0_qcd_sig ) ;
       printf("  mu0_ttbar_sb  : %6.1f\n", toy_mu0_ttbar_sb ) ;
       printf("  mu0_qcd_sb    : %6.1f\n", toy_mu0_qcd_sb ) ;
       printf("  mu0_susy_sig  : %6.1f\n", toy_mu0_susy_sig ) ;
       printf("  mu0_allbg_sig : %6.1f\n", toy_mu0_allbg_sig ) ;

     //--- generate some toy datasets.
       printf("\n\n Opening output toy datasets root file : %s\n", outputRootFile ) ;
       TTree* toyTree = new TTree("likelihoodData", "ra2b observed values generated from Nvals") ;

       double Nsig, Na, Nd,
              Nsb1, Nsb2, Nsb3, Nsb4, Nsb5,
              Nlsb1, Nlsb2, Nlsb3, Nlsb4, Nlsb5,
              Nslsig1, Nslsig2, Nslsig3, Nslsig4, Nslsig5,
              Nslsb1, Nslsb2, Nslsb3, Nslsb4, Nslsb5,
              Nslmsb1, Nslmsb2, Nslmsb3, Nslmsb4, Nslmsb5,
              Nqcdmca, Nqcdmcd, Nqcdmcsig, Nqcdmcsb ;

       toyTree->Branch("Nsig"       , &Nsig      , "Nsig/D" ) ;
       toyTree->Branch("Na"         , &Na        , "Na/D" ) ;
       toyTree->Branch("Nd"         , &Nd        , "Nd/D" ) ;
       toyTree->Branch("Nsb1"       , &Nsb1      , "Nsb1/D" ) ;
       toyTree->Branch("Nsb2"       , &Nsb2      , "Nsb2/D" ) ;
       toyTree->Branch("Nsb3"       , &Nsb3      , "Nsb3/D" ) ;
       toyTree->Branch("Nsb4"       , &Nsb4      , "Nsb4/D" ) ;
       toyTree->Branch("Nsb5"       , &Nsb5      , "Nsb5/D" ) ;
       toyTree->Branch("Nlsb1"      , &Nlsb1     , "Nlsb1/D" ) ;
       toyTree->Branch("Nlsb2"      , &Nlsb2     , "Nlsb2/D" ) ;
       toyTree->Branch("Nlsb3"      , &Nlsb3     , "Nlsb3/D" ) ;
       toyTree->Branch("Nlsb4"      , &Nlsb4     , "Nlsb4/D" ) ;
       toyTree->Branch("Nlsb5"      , &Nlsb5     , "Nlsb5/D" ) ;
       toyTree->Branch("Nslsig1"    , &Nslsig1   , "Nslsig1/D" ) ;
       toyTree->Branch("Nslsig2"    , &Nslsig2   , "Nslsig2/D" ) ;
       toyTree->Branch("Nslsig3"    , &Nslsig3   , "Nslsig3/D" ) ;
       toyTree->Branch("Nslsig4"    , &Nslsig4   , "Nslsig4/D" ) ;
       toyTree->Branch("Nslsig5"    , &Nslsig5   , "Nslsig5/D" ) ;
       toyTree->Branch("Nslsb1"     , &Nslsb1    , "Nslsb1/D" ) ;
       toyTree->Branch("Nslsb2"     , &Nslsb2    , "Nslsb2/D" ) ;
       toyTree->Branch("Nslsb3"     , &Nslsb3    , "Nslsb3/D" ) ;
       toyTree->Branch("Nslsb4"     , &Nslsb4    , "Nslsb4/D" ) ;
       toyTree->Branch("Nslsb5"     , &Nslsb5    , "Nslsb5/D" ) ;
       toyTree->Branch("Nslmsb1"    , &Nslmsb1   , "Nslmsb1/D" ) ;
       toyTree->Branch("Nslmsb2"    , &Nslmsb2   , "Nslmsb2/D" ) ;
       toyTree->Branch("Nslmsb3"    , &Nslmsb3   , "Nslmsb3/D" ) ;
       toyTree->Branch("Nslmsb4"    , &Nslmsb4   , "Nslmsb4/D" ) ;
       toyTree->Branch("Nslmsb5"    , &Nslmsb5   , "Nslmsb5/D" ) ;
       toyTree->Branch("Nqcdmca"    , &Nqcdmca   , "Nqcdmca/D" ) ;
       toyTree->Branch("Nqcdmcd"    , &Nqcdmcd   , "Nqcdmcd/D" ) ;
       toyTree->Branch("Nqcdmcsig"  , &Nqcdmcsig , "Nqcdmcsig/D" ) ;
       toyTree->Branch("Nqcdmcsb"   , &Nqcdmcsb  , "Nqcdmcsb/D" ) ;

       TRandom1 rg(12345) ;

       for ( int dsi=0; dsi<nToys; dsi++ ) {

          Nsig      = rg.Poisson( rv_Nsig      ->getVal() ) ;
          Na        = rg.Poisson( rv_Na        ->getVal() ) ;
          Nd        = rg.Poisson( rv_Nd        ->getVal() ) ;
          Nsb1      = rg.Poisson( rv_Nsb1      ->getVal() ) ;
          Nsb2      = rg.Poisson( rv_Nsb2      ->getVal() ) ;
          Nsb3      = rg.Poisson( rv_Nsb3      ->getVal() ) ;
          Nsb4      = rg.Poisson( rv_Nsb4      ->getVal() ) ;
          Nsb5      = rg.Poisson( rv_Nsb5      ->getVal() ) ;
          Nlsb1     = rg.Poisson( rv_Nlsb1     ->getVal() ) ;
          Nlsb2     = rg.Poisson( rv_Nlsb2     ->getVal() ) ;
          Nlsb3     = rg.Poisson( rv_Nlsb3     ->getVal() ) ;
          Nlsb4     = rg.Poisson( rv_Nlsb4     ->getVal() ) ;
          Nlsb5     = rg.Poisson( rv_Nlsb5     ->getVal() ) ;
          Nslsig1   = rg.Poisson( rv_Nslsig1   ->getVal() ) ;
          Nslsig2   = rg.Poisson( rv_Nslsig2   ->getVal() ) ;
          Nslsig3   = rg.Poisson( rv_Nslsig3   ->getVal() ) ;
          Nslsig4   = rg.Poisson( rv_Nslsig4   ->getVal() ) ;
          Nslsig5   = rg.Poisson( rv_Nslsig5   ->getVal() ) ;
          Nslsb1    = rg.Poisson( rv_Nslsb1    ->getVal() ) ;
          Nslsb2    = rg.Poisson( rv_Nslsb2    ->getVal() ) ;
          Nslsb3    = rg.Poisson( rv_Nslsb3    ->getVal() ) ;
          Nslsb4    = rg.Poisson( rv_Nslsb4    ->getVal() ) ;
          Nslsb5    = rg.Poisson( rv_Nslsb5    ->getVal() ) ;
          Nslmsb1   = rg.Poisson( rv_Nslmsb1   ->getVal() ) ;
          Nslmsb2   = rg.Poisson( rv_Nslmsb2   ->getVal() ) ;
          Nslmsb3   = rg.Poisson( rv_Nslmsb3   ->getVal() ) ;
          Nslmsb4   = rg.Poisson( rv_Nslmsb4   ->getVal() ) ;
          Nslmsb5   = rg.Poisson( rv_Nslmsb5   ->getVal() ) ;
          Nqcdmca   = rg.Poisson( rv_Nqcdmca   ->getVal() ) ;
          Nqcdmcd   = rg.Poisson( rv_Nqcdmcd   ->getVal() ) ;
          Nqcdmcsig = rg.Poisson( rv_Nqcdmcsig ->getVal() ) ;
          Nqcdmcsb  = rg.Poisson( rv_Nqcdmcsb  ->getVal() ) ;

          toyTree->Fill() ;

       } // dsi.


       TFile toyOutputFile( outputRootFile, "recreate" ) ;
       toyTree->Write() ;
       toyOutputTrueTree->Fill() ;
       toyOutputTrueTree->Write() ;
       printf("\n\n Closing output toy datasets root file : %s\n", outputRootFile ) ;
       toyOutputFile.Close() ;

       return true ;

    } // generateToyDatasetsFromNVals

  //==================================================================================================

    bool ra2bRoostatsClass::doFit( ) {

       if ( ! initialized ) {
          printf("\n\n *** Call initialize first.\n\n") ;
          return false ;
       }

       printf("\n\n") ;
       printf("  Fitting with these values for the observables.\n") ;
       dsObserved->printMultiline(cout, 1, kTRUE, "") ;
       printf("\n\n") ;

       fitResult = likelihood->fitTo(*dsObserved, Save(true));

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



       float Amc = rv_mu_qcdmc_a->getVal() ;
       float Dmc = rv_mu_qcdmc_d->getVal() ;
       float Bmc = rv_mu_qcdmc_sb->getVal() ;
       float Cmc = rv_mu_qcdmc_sig->getVal() ;

       float K = Amc * Cmc / ( Bmc * Dmc ) ;

       qcdCorrection = K ;

       float covAA = pow( rv_mu_qcdmc_a->getError(), 2 ) ;
       float covDD = pow( rv_mu_qcdmc_d->getError(), 2 ) ;
       float covBB = pow( rv_mu_qcdmc_sb->getError(), 2 ) ;
       float covCC = pow( rv_mu_qcdmc_sig->getError(), 2 ) ;

       float rhoAB = fitResult->correlation( "mu_qcdmc_a", "mu_qcdmc_sb" ) ;
       float rhoAC = fitResult->correlation( "mu_qcdmc_a", "mu_qcdmc_sig" ) ;
       float rhoAD = fitResult->correlation( "mu_qcdmc_a", "mu_qcdmc_d" ) ;

       float rhoBC = fitResult->correlation( "mu_qcdmc_sb", "mu_qcdmc_sig" ) ;
       float rhoBD = fitResult->correlation( "mu_qcdmc_sb", "mu_qcdmc_d" ) ;

       float rhoCD = fitResult->correlation( "mu_qcdmc_sig", "mu_qcdmc_d" ) ;

       float covAB = rhoAB * sqrt( covAA * covBB ) ;
       float covAC = rhoAC * sqrt( covAA * covCC ) ;
       float covAD = rhoAD * sqrt( covAA * covDD ) ;

       float covBC = rhoBC * sqrt( covBB * covCC ) ;
       float covBD = rhoBD * sqrt( covBB * covDD ) ;

       float covCD = rhoCD * sqrt( covCC * covDD ) ;

       qcdCorrectionErr = K * sqrt( 
              covAA/(Amc*Amc) + covBB/(Bmc*Bmc) + covCC/(Cmc*Cmc) + covDD/(Dmc*Dmc)
              -2 * covAB / (Amc*Bmc) +2 * covAC / (Amc*Cmc) -2 * covAD / (Amc*Dmc)
              +2 * covBC / (Bmc*Cmc) -2 * covBD / (Bmc*Dmc) 
              -2 * covCD / (Cmc*Dmc)
       ) ;


       printf("  QCD bias correction: %4.2f +/- %4.2f\n", qcdCorrection, qcdCorrectionErr ) ;

       varsAtFitVals = true ;

       return true ;

     } // doFit .


  //==================================================================================================

     bool ra2bRoostatsClass::profileSusySig( float& susySigLow, float& susySigHigh ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }

      //--- Profile likelihood for signal susy yield.

         ProfileLikelihoodCalculator plc_susy_sig( *dsObserved, *likelihood, RooArgSet( *rv_mu_susy_sig ) ) ;
         plc_susy_sig.SetTestSize(0.05) ;
         ConfInterval* plinterval_susy_sig = plc_susy_sig.GetInterval() ;
         susySigLow  = ((LikelihoodInterval*) plinterval_susy_sig)->LowerLimit(*rv_mu_susy_sig) ;
         susySigHigh = ((LikelihoodInterval*) plinterval_susy_sig)->UpperLimit(*rv_mu_susy_sig) ;
         printf("\n\n") ;
         printf("    susy, SIG 95%% CL interval  [%5.1f, %5.1f]\n\n", susySigLow, susySigHigh ) ;
         TCanvas* plcplot_susy_sig = new TCanvas("plcplot_susy_sig", "susy sig, Profile likelihood", 500, 400 ) ;
         LikelihoodIntervalPlot plotInt_susy_sig((LikelihoodInterval*)plinterval_susy_sig);
         plotInt_susy_sig.Draw() ;
         plcplot_susy_sig->SaveAs("plscan_susy_sig.pdf") ;
         plcplot_susy_sig->SaveAs("plscan_susy_sig.png") ;

         varsAtFitVals = false ;

         return true ;

     }

  //==================================================================================================

     bool ra2bRoostatsClass::profileTtbarSig( float& ttbarSigLow, float& ttbarSigHigh ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }

         if ( ! useSigBgVars ) {
            printf("\n\n\n *** Can't do it.  Need to use ra2bRoostatsClass with constructor arg useSigBgVars set to true.\n\n") ;
            return false ;
         }

      //--- Profile likelihood for SIG ttbar yield.

         ProfileLikelihoodCalculator plc_ttbar_sig( *dsObserved, *likelihood, RooArgSet( *rv_mu_ttbar_sig ) ) ;
         plc_ttbar_sig.SetTestSize(0.32) ;
         ConfInterval* plinterval_ttbar_sig = plc_ttbar_sig.GetInterval() ;
         ttbarSigLow  = ((LikelihoodInterval*) plinterval_ttbar_sig)->LowerLimit(*rrv_mu_ttbar_sig) ;
         ttbarSigHigh = ((LikelihoodInterval*) plinterval_ttbar_sig)->UpperLimit(*rrv_mu_ttbar_sig) ;
         printf("\n\n") ;
         printf("    ttbar, SIG 68%% CL interval  [%5.1f, %5.1f]\n\n", ttbarSigLow, ttbarSigHigh ) ;
         TCanvas* plcplot_ttbar_sig = new TCanvas("plcplot_ttbar_sig", "ttbar sig, Profile likelihood", 500, 400 ) ;
         LikelihoodIntervalPlot plotInt_ttbar_sig((LikelihoodInterval*)plinterval_ttbar_sig);
         plotInt_ttbar_sig.Draw() ;
         plcplot_ttbar_sig->SaveAs("plscan_ttbar_sig.pdf") ;
         plcplot_ttbar_sig->SaveAs("plscan_ttbar_sig.png") ;

         varsAtFitVals = false ;

         return true ;

     }

  //==================================================================================================

     bool ra2bRoostatsClass::profileQcdSig( float& qcdSigLow, float& qcdSigHigh ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }

         if ( ! useSigBgVars ) {
            printf("\n\n\n *** Can't do it.  Need to use ra2bRoostatsClass with constructor arg useSigBgVars set to true.\n\n") ;
            return false ;
         }

      //--- Profile likelihood for SIG qcd yield.

         ProfileLikelihoodCalculator plc_qcd_sig( *dsObserved, *likelihood, RooArgSet( *rv_mu_qcd_sig ) ) ;
         plc_qcd_sig.SetTestSize(0.32) ;
         ConfInterval* plinterval_qcd_sig = plc_qcd_sig.GetInterval() ;
         qcdSigLow  = ((LikelihoodInterval*) plinterval_qcd_sig)->LowerLimit(*rrv_mu_qcd_sig) ;
         qcdSigHigh = ((LikelihoodInterval*) plinterval_qcd_sig)->UpperLimit(*rrv_mu_qcd_sig) ;
         printf("\n\n") ;
         printf("    qcd, SIG 68%% CL interval  [%5.1f, %5.1f]\n\n", qcdSigLow, qcdSigHigh ) ;
         TCanvas* plcplot_qcd_sig = new TCanvas("plcplot_qcd_sig", "qcd sig, Profile likelihood", 500, 400 ) ;
         LikelihoodIntervalPlot plotInt_qcd_sig((LikelihoodInterval*)plinterval_qcd_sig);
         plotInt_qcd_sig.Draw() ;
         plcplot_qcd_sig->SaveAs("plscan_qcd_sig.pdf") ;
         plcplot_qcd_sig->SaveAs("plscan_qcd_sig.png") ;

         varsAtFitVals = false ;

         return true ;

     }

  //==================================================================================================


     bool ra2bRoostatsClass::profileTtbarSb( float& ttbarSbLow, float& ttbarSbHigh ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }

         if ( useSigBgVars ) {
            printf("\n\n\n *** Can't do it.  Need to use ra2bRoostatsClass with constructor arg useSigBgVars set to false.\n\n") ;
            return false ;
         }

      //--- Profile likelihood for SB ttbar yield.

         ProfileLikelihoodCalculator plc_ttbar_sb( *dsObserved, *likelihood, RooArgSet( *rv_mu_ttbar_sb ) ) ;
         plc_ttbar_sb.SetTestSize(0.32) ;
         ConfInterval* plinterval_ttbar_sb = plc_ttbar_sb.GetInterval() ;
         ttbarSbLow  = ((LikelihoodInterval*) plinterval_ttbar_sb)->LowerLimit(*rrv_mu_ttbar_sb) ;
         ttbarSbHigh = ((LikelihoodInterval*) plinterval_ttbar_sb)->UpperLimit(*rrv_mu_ttbar_sb) ;
         printf("\n\n") ;
         printf("    ttbar, SB 68%% CL interval  [%5.1f, %5.1f]\n\n", ttbarSbLow, ttbarSbHigh ) ;
         TCanvas* plcplot_ttbar_sb = new TCanvas("plcplot_ttbar_sb", "ttbar sb, Profile likelihood", 500, 400 ) ;
         LikelihoodIntervalPlot plotInt_ttbar_sb((LikelihoodInterval*)plinterval_ttbar_sb);
         plotInt_ttbar_sb.Draw() ;
         plcplot_ttbar_sb->SaveAs("plscan_ttbar_sb.pdf") ;
         plcplot_ttbar_sb->SaveAs("plscan_ttbar_sb.png") ;

         varsAtFitVals = false ;

         return true ;

     }

  //==================================================================================================


     bool ra2bRoostatsClass::profileQcdSb( float& qcdSbLow, float& qcdSbHigh ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }

         if ( useSigBgVars ) {
            printf("\n\n\n *** Can't do it.  Need to use ra2bRoostatsClass with constructor arg useSigBgVars set to false.\n\n") ;
            return false ;
         }

      //--- Profile likelihood for SB qcd yield.

         ProfileLikelihoodCalculator plc_qcd_sb( *dsObserved, *likelihood, RooArgSet( *rv_mu_qcd_sb ) ) ;
         plc_qcd_sb.SetTestSize(0.32) ;
         ConfInterval* plinterval_qcd_sb = plc_qcd_sb.GetInterval() ;
         qcdSbLow  = ((LikelihoodInterval*) plinterval_qcd_sb)->LowerLimit(*rrv_mu_qcd_sb) ;
         qcdSbHigh = ((LikelihoodInterval*) plinterval_qcd_sb)->UpperLimit(*rrv_mu_qcd_sb) ;
         printf("\n\n") ;
         printf("    qcd, SB 68%% CL interval  [%5.1f, %5.1f]\n\n", qcdSbLow, qcdSbHigh ) ;
         TCanvas* plcplot_qcd_sb = new TCanvas("plcplot_qcd_sb", "qcd sb, Profile likelihood", 500, 400 ) ;
         LikelihoodIntervalPlot plotInt_qcd_sb((LikelihoodInterval*)plinterval_qcd_sb);
         plotInt_qcd_sb.Draw() ;
         plcplot_qcd_sb->SaveAs("plscan_qcd_sb.pdf") ;
         plcplot_qcd_sb->SaveAs("plscan_qcd_sb.png") ;

         varsAtFitVals = false ;

         return true ;

     }

  //==================================================================================================



///    float ttbar_sb_val(0.) ;
///    float qcd_sb_val(0.) ;
///    float ttbar_sb_err(0.) ;
///    float qcd_sb_err(0.) ;

///    float ttbar_sig_val(0.) ;
///    float qcd_sig_val(0.) ;
///    float ttbar_sig_err(0.) ;
///    float qcd_sig_err(0.) ;

///    if ( useSigBgVars ) {
///       ttbar_sb_val = ((RooFormulaVar*)rv_mu_ttbar_sb) ->getVal() ;
///       qcd_sb_val = ((RooFormulaVar*)rv_mu_qcd_sb) ->getVal() ;
///       ttbar_sig_val = ((RooRealVar*)rv_mu_ttbar_sig) ->getVal() ;
///       qcd_sig_val = ((RooRealVar*)rv_mu_qcd_sig) ->getVal() ;
///       ttbar_sig_err = ((RooRealVar*)rv_mu_ttbar_sig) ->getError() ;
///       qcd_sig_err = ((RooRealVar*)rv_mu_qcd_sig) ->getError() ;
///    } else {
///       ttbar_sig_val = ((RooFormulaVar*)rv_mu_ttbar_sig) ->getVal() ;
///       qcd_sig_val = ((RooFormulaVar*)rv_mu_qcd_sig) ->getVal() ;
///       ttbar_sb_val = ((RooRealVar*)rv_mu_ttbar_sb) ->getVal() ;
///       qcd_sb_val = ((RooRealVar*)rv_mu_qcd_sb) ->getVal() ;
///       ttbar_sb_err = ((RooRealVar*)rv_mu_ttbar_sb) ->getError() ;
///       qcd_sb_err = ((RooRealVar*)rv_mu_qcd_sb) ->getError() ;
///    }


///    if ( useSigBgVars ) {

///    //--- Profile likelihood for signal ttbar yield.

///        ProfileLikelihoodCalculator plc_ttbar_sig( *dsObserved, *likelihood, RooArgSet( *rv_mu_ttbar_sig ) ) ;
///        plc_ttbar_sig.SetTestSize(0.32) ;
///        ConfInterval* plinterval_ttbar_sig = plc_ttbar_sig.GetInterval() ;
///        float ttbar_sig_p1sig = ((LikelihoodInterval*) plinterval_ttbar_sig)->UpperLimit(*((RooRealVar*)rv_mu_ttbar_sig)) ;
///        float ttbar_sig_m1sig = ((LikelihoodInterval*) plinterval_ttbar_sig)->LowerLimit(*((RooRealVar*)rv_mu_ttbar_sig)) ;
///        printf("\n\n") ;
///        printf("    ttbar, SIG 68%% CL interval  [%5.1f, %5.1f]\n\n", ttbar_sig_m1sig, ttbar_sig_p1sig) ;
///        TCanvas* plcplot_ttbar_sig = new TCanvas("plcplot_ttbar_sig", "ttbar sig, Profile likelihood", 500, 400 ) ;
///        LikelihoodIntervalPlot plotInt_ttbar_sig((LikelihoodInterval*)plinterval_ttbar_sig);
///        plotInt_ttbar_sig.Draw() ;
///        plcplot_ttbar_sig->SaveAs("plscan_ttbar_sig.pdf") ;
///        plcplot_ttbar_sig->SaveAs("plscan_ttbar_sig.png") ;

///    //--- Profile likelihood for signal qcd yield.

///        ProfileLikelihoodCalculator plc_qcd_sig( *dsObserved, *likelihood, RooArgSet( *rv_mu_qcd_sig ) ) ;
///        plc_qcd_sig.SetTestSize(0.32) ;
///        ConfInterval* plinterval_qcd_sig = plc_qcd_sig.GetInterval() ;
///        float qcd_sig_p1sig = ((LikelihoodInterval*) plinterval_qcd_sig)->UpperLimit(*((RooRealVar*)rv_mu_qcd_sig)) ;
///        float qcd_sig_m1sig = ((LikelihoodInterval*) plinterval_qcd_sig)->LowerLimit(*((RooRealVar*)rv_mu_qcd_sig)) ;
///        printf("\n\n") ;
///        printf("    qcd, SIG 68%% CL interval  [%5.1f, %5.1f]\n\n", qcd_sig_m1sig, qcd_sig_p1sig) ;
///        TCanvas* plcplot_qcd_sig = new TCanvas("plcplot_qcd_sig", "qcd sig, Profile likelihood", 500, 400 ) ;
///        LikelihoodIntervalPlot plotInt_qcd_sig((LikelihoodInterval*)plinterval_qcd_sig);
///        plotInt_qcd_sig.Draw() ;
///        plcplot_qcd_sig->SaveAs("plscan_qcd_sig.pdf") ;
///        plcplot_qcd_sig->SaveAs("plscan_qcd_sig.png") ;

///    } else {

///    //--- Profile likelihood for sideband ttbar yield.

///        ProfileLikelihoodCalculator plc_ttbar_sb( *dsObserved, *likelihood, RooArgSet( *rv_mu_ttbar_sb ) ) ;
///        plc_ttbar_sb.SetTestSize(0.32) ;
///        ConfInterval* plinterval_ttbar_sb = plc_ttbar_sb.GetInterval() ;
///        float ttbar_sb_p1sb = ((LikelihoodInterval*) plinterval_ttbar_sb)->UpperLimit(*((RooRealVar*)rv_mu_ttbar_sb)) ;
///        float ttbar_sb_m1sb = ((LikelihoodInterval*) plinterval_ttbar_sb)->LowerLimit(*((RooRealVar*)rv_mu_ttbar_sb)) ;
///        printf("\n\n") ;
///        printf("    ttbar, SB 68%% CL interval  [%5.1f, %5.1f]\n\n", ttbar_sb_m1sb, ttbar_sb_p1sb) ;
///        TCanvas* plcplot_ttbar_sb = new TCanvas("plcplot_ttbar_sb", "ttbar sb, Profile likelihood", 500, 400 ) ;
///        LikelihoodIntervalPlot plotInt_ttbar_sb((LikelihoodInterval*)plinterval_ttbar_sb);
///        plotInt_ttbar_sb.Draw() ;
///        plcplot_ttbar_sb->SaveAs("plscan_ttbar_sb.pdf") ;
///        plcplot_ttbar_sb->SaveAs("plscan_ttbar_sb.png") ;

///    //--- Profile likelihood for sideband qcd yield.

///        ProfileLikelihoodCalculator plc_qcd_sb( *dsObserved, *likelihood, RooArgSet( *rv_mu_qcd_sb ) ) ;
///        plc_qcd_sb.SetTestSize(0.32) ;
///        ConfInterval* plinterval_qcd_sb = plc_qcd_sb.GetInterval() ;
///        float qcd_sb_p1sb = ((LikelihoodInterval*) plinterval_qcd_sb)->UpperLimit(*((RooRealVar*)rv_mu_qcd_sb)) ;
///        float qcd_sb_m1sb = ((LikelihoodInterval*) plinterval_qcd_sb)->LowerLimit(*((RooRealVar*)rv_mu_qcd_sb)) ;
///        printf("\n\n") ;
///        printf("    sb, SB 68%% CL interval  [%5.1f, %5.1f]\n\n", qcd_sb_m1sb, qcd_sb_p1sb) ;
///        TCanvas* plcplot_qcd_sb = new TCanvas("plcplot_qcd_sb", "qcd sb, Profile likelihood", 500, 400 ) ;
///        LikelihoodIntervalPlot plotInt_qcd_sb((LikelihoodInterval*)plinterval_qcd_sb);
///        plotInt_qcd_sb.Draw() ;
///        plcplot_qcd_sb->SaveAs("plscan_qcd_sb.pdf") ;
///        plcplot_qcd_sb->SaveAs("plscan_qcd_sb.png") ;

///    }



    //====================================================================================================================

       bool ra2bRoostatsClass::sbPlotsUniformBins( const char* plotBaseName ) {

          if ( ! initialized ) {
             printf("\n\n *** Call initialize first.\n\n") ;
             return false ;
          }


          if ( ! varsAtFitVals ) {
             printf("\n\n *** Try this right after calling doFit.\n\n") ;
             return false ;
          }

       //--  Drawn with same-width bins.

          TH1F* hm3j_sb_data  = new TH1F("hm3j_sb_data" ,"3-jet mass, SB data" , 5, 0.5, 5.5 ) ;
          TH1F* hm3j_sb_ttbar = new TH1F("hm3j_sb_ttbar","3-jet mass, SB ttbar", 5, 0.5, 5.5 ) ;
          TH1F* hm3j_sb_qcd   = new TH1F("hm3j_sb_qcd"  ,"3-jet mass, SB qcd"  , 5, 0.5, 5.5 ) ;
          TH1F* hm3j_sb_ew    = new TH1F("hm3j_sb_ew"   ,"3-jet mass, SB ew"   , 5, 0.5, 5.5 ) ;

          hm3j_sb_data->SetLineWidth(2) ;
          hm3j_sb_data->SetMarkerStyle(20) ;
          hm3j_sb_ew->SetFillColor(49) ;
          hm3j_sb_qcd->SetFillColor(46) ;
          hm3j_sb_ttbar->SetFillColor(42) ;

          hm3j_sb_data->SetBinContent( 1, rv_Nsb1->getVal() ) ;
          hm3j_sb_data->SetBinContent( 2, rv_Nsb2->getVal() ) ;
          hm3j_sb_data->SetBinContent( 3, rv_Nsb3->getVal() ) ;
          hm3j_sb_data->SetBinContent( 4, rv_Nsb4->getVal() ) ;
          hm3j_sb_data->SetBinContent( 5, rv_Nsb5->getVal() ) ;

          hm3j_sb_ttbar->SetBinContent( 1, rv_mu_ttbar_sb1->getVal() ) ;
          hm3j_sb_ttbar->SetBinContent( 2, rv_mu_ttbar_sb2->getVal() ) ;
          hm3j_sb_ttbar->SetBinContent( 3, rv_mu_ttbar_sb3->getVal() ) ;
          hm3j_sb_ttbar->SetBinContent( 4, rv_mu_ttbar_sb4->getVal() ) ;
          hm3j_sb_ttbar->SetBinContent( 5, rv_mu_ttbar_sb5->getVal() ) ;

          hm3j_sb_qcd->SetBinContent( 1, rv_mu_qcd_sb1->getVal() ) ;
          hm3j_sb_qcd->SetBinContent( 2, rv_mu_qcd_sb2->getVal() ) ;
          hm3j_sb_qcd->SetBinContent( 3, rv_mu_qcd_sb3->getVal() ) ;
          hm3j_sb_qcd->SetBinContent( 4, rv_mu_qcd_sb4->getVal() ) ;
          hm3j_sb_qcd->SetBinContent( 5, rv_mu_qcd_sb5->getVal() ) ;

          hm3j_sb_ew->SetBinContent( 1, rv_mu_ew_sb1->getVal() ) ;
          hm3j_sb_ew->SetBinContent( 2, rv_mu_ew_sb2->getVal() ) ;
          hm3j_sb_ew->SetBinContent( 3, rv_mu_ew_sb3->getVal() ) ;
          hm3j_sb_ew->SetBinContent( 4, rv_mu_ew_sb4->getVal() ) ;
          hm3j_sb_ew->SetBinContent( 5, rv_mu_ew_sb5->getVal() ) ;

          THStack* hstack_m3j_sb_fit = new THStack( "hstack_m3j_sb_fit", "SB fit, 3-jet mass" ) ;
          hstack_m3j_sb_fit->Add( hm3j_sb_ew ) ;
          hstack_m3j_sb_fit->Add( hm3j_sb_qcd ) ;
          hstack_m3j_sb_fit->Add( hm3j_sb_ttbar ) ;


          TCanvas* can_sbfit = new TCanvas("can_sbfit", "SB 3-jet mass fit", 700, 500 ) ;

          hm3j_sb_data->SetMaximum( 1.4*(hm3j_sb_data->GetMaximum()) ) ;
          hm3j_sb_data->SetLabelSize(0.06,"x") ;
          hm3j_sb_data->GetXaxis()->SetBinLabel(1, "Bin 1") ;
          hm3j_sb_data->GetXaxis()->SetBinLabel(2, "Bin 2") ;
          hm3j_sb_data->GetXaxis()->SetBinLabel(3, "Bin 3") ;
          hm3j_sb_data->GetXaxis()->SetBinLabel(4, "Bin 4") ;
          hm3j_sb_data->GetXaxis()->SetBinLabel(5, "Bin 5") ;


          hm3j_sb_data->Draw("histpe") ;
          hstack_m3j_sb_fit->Draw( "same" ) ;
          hm3j_sb_data->Draw("samehistpe") ;

          TLegend* m3j_legend = new TLegend(0.62,0.7,0.97,0.95) ;
          m3j_legend->SetFillColor( kWhite ) ;
          m3j_legend->AddEntry( hm3j_sb_data, "Data" ) ;
          m3j_legend->AddEntry( hm3j_sb_ttbar, "ttbar" ) ;
          m3j_legend->AddEntry( hm3j_sb_qcd, "QCD" ) ;
          m3j_legend->AddEntry( hm3j_sb_ew, "EW" ) ;
          m3j_legend->Draw() ;

          TText* fittext = new TText() ;
          fittext->SetTextSize(0.04) ;
          char fitlabel[1000] ;

          float ew_sb(0.) ;
          ew_sb += rv_mu_ew_sb1->getVal() ;
          ew_sb += rv_mu_ew_sb2->getVal() ;
          ew_sb += rv_mu_ew_sb3->getVal() ;
          ew_sb += rv_mu_ew_sb4->getVal() ;
          ew_sb += rv_mu_ew_sb5->getVal() ;
          int Nsb = rv_Nsb1->getVal() +rv_Nsb2->getVal() +rv_Nsb3->getVal() +rv_Nsb4->getVal() +rv_Nsb5->getVal()  ;

          float ttbar_sb_val(0.) ;
          float qcd_sb_val(0.) ;
          float ttbar_sb_err(0.) ;
          float qcd_sb_err(0.) ;

          float ttbar_sig_val(0.) ;
          float qcd_sig_val(0.) ;
          float ttbar_sig_err(0.) ;
          float qcd_sig_err(0.) ;

          if ( useSigBgVars ) {
             ttbar_sb_val = ((RooFormulaVar*)rv_mu_ttbar_sb) ->getVal() ;
             qcd_sb_val = ((RooFormulaVar*)rv_mu_qcd_sb) ->getVal() ;
             ttbar_sig_val = ((RooRealVar*)rv_mu_ttbar_sig) ->getVal() ;
             qcd_sig_val = ((RooRealVar*)rv_mu_qcd_sig) ->getVal() ;
             ttbar_sig_err = ((RooRealVar*)rv_mu_ttbar_sig) ->getError() ;
             qcd_sig_err = ((RooRealVar*)rv_mu_qcd_sig) ->getError() ;
          } else {
             ttbar_sig_val = ((RooFormulaVar*)rv_mu_ttbar_sig) ->getVal() ;
             qcd_sig_val = ((RooFormulaVar*)rv_mu_qcd_sig) ->getVal() ;
             ttbar_sb_val = ((RooRealVar*)rv_mu_ttbar_sb) ->getVal() ;
             qcd_sb_val = ((RooRealVar*)rv_mu_qcd_sb) ->getVal() ;
             ttbar_sb_err = ((RooRealVar*)rv_mu_ttbar_sb) ->getError() ;
             qcd_sb_err = ((RooRealVar*)rv_mu_qcd_sb) ->getError() ;
          }

          float ltop = 0.90 ;
          float lx = 0.78 ;
          float dy = 0.06 ;
          if ( useSigBgVars ) {
             sprintf( fitlabel, "%5d", Nsb ) ;
             fittext->DrawTextNDC( lx, ltop, fitlabel ) ;
             sprintf( fitlabel, "%4.0f", ttbar_sb_val ) ;
             fittext->DrawTextNDC( lx, ltop-dy, fitlabel ) ;
             sprintf( fitlabel, "%4.0f", qcd_sb_val ) ;
             fittext->DrawTextNDC( lx, ltop-2*dy, fitlabel ) ;
             sprintf( fitlabel, "%4.0f", ew_sb ) ;
             fittext->DrawTextNDC( lx, ltop-3*dy, fitlabel ) ;
          } else {
             sprintf( fitlabel, "%5d", Nsb ) ;
             fittext->DrawTextNDC( lx, ltop, fitlabel ) ;
             sprintf( fitlabel, "%4.0f +/- %4.0f", ttbar_sb_val, ttbar_sb_err ) ;
             fittext->DrawTextNDC( lx, ltop-dy, fitlabel ) ;
             sprintf( fitlabel, "%4.0f +/- %4.0f", qcd_sb_val, qcd_sb_err ) ;
             fittext->DrawTextNDC( lx, ltop-2*dy, fitlabel ) ;
             sprintf( fitlabel, "%4.0f", ew_sb ) ;
             fittext->DrawTextNDC( lx, ltop-3*dy, fitlabel ) ;
          }

          char outfile[1000] ;
          sprintf( outfile, "%s-sb.png", plotBaseName ) ;
          can_sbfit->SaveAs( outfile ) ;

          return true ;

       } // sbPlotsUniformBins

  //==================================================================================================================

       bool ra2bRoostatsClass::sbPlotsVariableBins( const char* plotBaseName ) {

          if ( ! initialized ) {
             printf("\n\n *** Call initialize first.\n\n") ;
             return false ;
          }

          if ( ! varsAtFitVals ) {
             printf("\n\n *** Try this right after calling doFit.\n\n") ;
             return false ;
          }

       //--  Drawn with variable-width bins.

          int nxbins = 5 ;
          float xbinedges[6] = { 0., 160., 180., 260., 400., 800. } ;
          float xbinwid[5] ;
          for ( int bi=0; bi<nxbins; bi++ ) { xbinwid[bi] = xbinedges[bi+1] - xbinedges[bi] ; }

          TH1F* hvbm3j_sb_data  = new TH1F("hvbm3j_sb_data" ,"3-jet mass, SB data" , nxbins, xbinedges ) ;
          TH1F* hvbm3j_sb_ttbar = new TH1F("hvbm3j_sb_ttbar","3-jet mass, SB ttbar", nxbins, xbinedges ) ;
          TH1F* hvbm3j_sb_qcd   = new TH1F("hvbm3j_sb_qcd"  ,"3-jet mass, SB qcd"  , nxbins, xbinedges ) ;
          TH1F* hvbm3j_sb_ew    = new TH1F("hvbm3j_sb_ew"   ,"3-jet mass, SB ew"   , nxbins, xbinedges ) ;

          hvbm3j_sb_data->SetLineWidth(2) ;
          hvbm3j_sb_data->SetMarkerStyle(20) ;
          hvbm3j_sb_ew->SetFillColor(49) ;
          hvbm3j_sb_qcd->SetFillColor(46) ;
          hvbm3j_sb_ttbar->SetFillColor(42) ;

          hvbm3j_sb_data->SetBinContent( 1, (rv_Nsb1->getVal()) * (xbinwid[1]/xbinwid[0]) ) ;
          hvbm3j_sb_data->SetBinContent( 2, (rv_Nsb2->getVal()) * (xbinwid[1]/xbinwid[1]) ) ;
          hvbm3j_sb_data->SetBinContent( 3, (rv_Nsb3->getVal()) * (xbinwid[1]/xbinwid[2]) ) ;
          hvbm3j_sb_data->SetBinContent( 4, (rv_Nsb4->getVal()) * (xbinwid[1]/xbinwid[3]) ) ;
          hvbm3j_sb_data->SetBinContent( 5, (rv_Nsb5->getVal()) * (xbinwid[1]/xbinwid[4]) ) ;

          hvbm3j_sb_data->SetBinError( 1, sqrt(rv_Nsb1->getVal()) * (xbinwid[1]/xbinwid[0]) ) ;
          hvbm3j_sb_data->SetBinError( 2, sqrt(rv_Nsb2->getVal()) * (xbinwid[1]/xbinwid[1]) ) ;
          hvbm3j_sb_data->SetBinError( 3, sqrt(rv_Nsb3->getVal()) * (xbinwid[1]/xbinwid[2]) ) ;
          hvbm3j_sb_data->SetBinError( 4, sqrt(rv_Nsb4->getVal()) * (xbinwid[1]/xbinwid[3]) ) ;
          hvbm3j_sb_data->SetBinError( 5, sqrt(rv_Nsb5->getVal()) * (xbinwid[1]/xbinwid[4]) ) ;

          hvbm3j_sb_ttbar->SetBinContent( 1, (rv_mu_ttbar_sb1->getVal()) * (xbinwid[1]/xbinwid[0]) ) ;
          hvbm3j_sb_ttbar->SetBinContent( 2, (rv_mu_ttbar_sb2->getVal()) * (xbinwid[1]/xbinwid[1]) ) ;
          hvbm3j_sb_ttbar->SetBinContent( 3, (rv_mu_ttbar_sb3->getVal()) * (xbinwid[1]/xbinwid[2]) ) ;
          hvbm3j_sb_ttbar->SetBinContent( 4, (rv_mu_ttbar_sb4->getVal()) * (xbinwid[1]/xbinwid[3]) ) ;
          hvbm3j_sb_ttbar->SetBinContent( 5, (rv_mu_ttbar_sb5->getVal()) * (xbinwid[1]/xbinwid[4]) ) ;

          hvbm3j_sb_qcd->SetBinContent( 1, (rv_mu_qcd_sb1->getVal()) * (xbinwid[1]/xbinwid[0]) ) ;
          hvbm3j_sb_qcd->SetBinContent( 2, (rv_mu_qcd_sb2->getVal()) * (xbinwid[1]/xbinwid[1]) ) ;
          hvbm3j_sb_qcd->SetBinContent( 3, (rv_mu_qcd_sb3->getVal()) * (xbinwid[1]/xbinwid[2]) ) ;
          hvbm3j_sb_qcd->SetBinContent( 4, (rv_mu_qcd_sb4->getVal()) * (xbinwid[1]/xbinwid[3]) ) ;
          hvbm3j_sb_qcd->SetBinContent( 5, (rv_mu_qcd_sb5->getVal()) * (xbinwid[1]/xbinwid[4]) ) ;

          hvbm3j_sb_ew->SetBinContent( 1, (rv_mu_ew_sb1->getVal()) * (xbinwid[1]/xbinwid[0]) ) ;
          hvbm3j_sb_ew->SetBinContent( 2, (rv_mu_ew_sb2->getVal()) * (xbinwid[1]/xbinwid[1]) ) ;
          hvbm3j_sb_ew->SetBinContent( 3, (rv_mu_ew_sb3->getVal()) * (xbinwid[1]/xbinwid[2]) ) ;
          hvbm3j_sb_ew->SetBinContent( 4, (rv_mu_ew_sb4->getVal()) * (xbinwid[1]/xbinwid[3]) ) ;
          hvbm3j_sb_ew->SetBinContent( 5, (rv_mu_ew_sb5->getVal()) * (xbinwid[1]/xbinwid[4]) ) ;

          THStack* hvbstack_m3j_sb_fit = new THStack( "hvbstack_m3j_sb_fit", "SB fit, 3-jet mass" ) ;
          hvbstack_m3j_sb_fit->Add( hvbm3j_sb_ew ) ;
          hvbstack_m3j_sb_fit->Add( hvbm3j_sb_qcd ) ;
          hvbstack_m3j_sb_fit->Add( hvbm3j_sb_ttbar ) ;

          TCanvas* can_vbsbfit = new TCanvas("can_vbsbfit", "SB 3-jet mass fit", 700, 500 ) ;

          //hvbm3j_sb_data->SetMaximum( 1.4*(hvbm3j_sb_data->GetMaximum()) ) ;


          hvbm3j_sb_data->Draw("histpe") ;
          hvbstack_m3j_sb_fit->Draw( "same" ) ;
          hvbm3j_sb_data->Draw("samehistpe") ;

          TLegend* vbm3j_legend = new TLegend(0.62,0.7,0.97,0.95) ;
          vbm3j_legend->SetFillColor( kWhite ) ;
          vbm3j_legend->AddEntry( hvbm3j_sb_data, "Data" ) ;
          vbm3j_legend->AddEntry( hvbm3j_sb_ttbar, "ttbar" ) ;
          vbm3j_legend->AddEntry( hvbm3j_sb_qcd, "QCD" ) ;
          vbm3j_legend->AddEntry( hvbm3j_sb_ew, "EW" ) ;
          vbm3j_legend->Draw() ;

          TText* fittext = new TText() ;
          fittext->SetTextSize(0.04) ;
          char fitlabel[1000] ;

          float ew_sb(0.) ;
          ew_sb += rv_mu_ew_sb1->getVal() ;
          ew_sb += rv_mu_ew_sb2->getVal() ;
          ew_sb += rv_mu_ew_sb3->getVal() ;
          ew_sb += rv_mu_ew_sb4->getVal() ;
          ew_sb += rv_mu_ew_sb5->getVal() ;
          int Nsb = rv_Nsb1->getVal() +rv_Nsb2->getVal() +rv_Nsb3->getVal() +rv_Nsb4->getVal() +rv_Nsb5->getVal()  ;
          float ltop = 0.90 ;
          float lx = 0.78 ;
          float dy = 0.06 ;

          float ttbar_sb_val(0.) ;
          float qcd_sb_val(0.) ;
          float ttbar_sb_err(0.) ;
          float qcd_sb_err(0.) ;

          float ttbar_sig_val(0.) ;
          float qcd_sig_val(0.) ;
          float ttbar_sig_err(0.) ;
          float qcd_sig_err(0.) ;

          if ( useSigBgVars ) {
             ttbar_sb_val = ((RooFormulaVar*)rv_mu_ttbar_sb) ->getVal() ;
             qcd_sb_val = ((RooFormulaVar*)rv_mu_qcd_sb) ->getVal() ;
             ttbar_sig_val = ((RooRealVar*)rv_mu_ttbar_sig) ->getVal() ;
             qcd_sig_val = ((RooRealVar*)rv_mu_qcd_sig) ->getVal() ;
             ttbar_sig_err = ((RooRealVar*)rv_mu_ttbar_sig) ->getError() ;
             qcd_sig_err = ((RooRealVar*)rv_mu_qcd_sig) ->getError() ;
          } else {
             ttbar_sig_val = ((RooFormulaVar*)rv_mu_ttbar_sig) ->getVal() ;
             qcd_sig_val = ((RooFormulaVar*)rv_mu_qcd_sig) ->getVal() ;
             ttbar_sb_val = ((RooRealVar*)rv_mu_ttbar_sb) ->getVal() ;
             qcd_sb_val = ((RooRealVar*)rv_mu_qcd_sb) ->getVal() ;
             ttbar_sb_err = ((RooRealVar*)rv_mu_ttbar_sb) ->getError() ;
             qcd_sb_err = ((RooRealVar*)rv_mu_qcd_sb) ->getError() ;
          }

          if ( useSigBgVars ) {
             sprintf( fitlabel, "%5d", Nsb ) ;
             fittext->DrawTextNDC( lx, ltop, fitlabel ) ;
             sprintf( fitlabel, "%4.0f", ttbar_sb_val ) ;
             fittext->DrawTextNDC( lx, ltop-dy, fitlabel ) ;
             sprintf( fitlabel, "%4.0f", qcd_sb_val ) ;
             fittext->DrawTextNDC( lx, ltop-2*dy, fitlabel ) ;
             sprintf( fitlabel, "%4.0f", ew_sb ) ;
             fittext->DrawTextNDC( lx, ltop-3*dy, fitlabel ) ;
          } else {
             sprintf( fitlabel, "%5d", Nsb ) ;
             fittext->DrawTextNDC( lx, ltop, fitlabel ) ;
             sprintf( fitlabel, "%4.0f +/- %4.0f", ttbar_sb_val, ttbar_sb_err ) ;
             fittext->DrawTextNDC( lx, ltop-dy, fitlabel ) ;
             sprintf( fitlabel, "%4.0f +/- %4.0f", qcd_sb_val, qcd_sb_err ) ;
             fittext->DrawTextNDC( lx, ltop-2*dy, fitlabel ) ;
             sprintf( fitlabel, "%4.0f", ew_sb ) ;
             fittext->DrawTextNDC( lx, ltop-3*dy, fitlabel ) ;
          }

          char outfile[1000] ;
          sprintf( outfile, "%s-sb-vb.png", plotBaseName ) ;
          can_vbsbfit->SaveAs( outfile ) ;

          return true ;

       } // sbPlotsVariableBins

  //===================================================================

   ra2bRoostatsClass::~ra2bRoostatsClass() {

      if ( !initialized ) return ;

      delete dsObserved ;
      //delete plinterval_susy_sig ;

      delete rv_Nsig ;
      delete rv_Na ;
      delete rv_Nd ;
      delete rv_Nsb1 ;
      delete rv_Nsb2 ;
      delete rv_Nsb3 ;
      delete rv_Nsb4 ;
      delete rv_Nsb5 ;
      delete rv_Nlsb1 ;
      delete rv_Nlsb2 ;
      delete rv_Nlsb3 ;
      delete rv_Nlsb4 ;
      delete rv_Nlsb5 ;
      delete rv_Nslsig1 ;
      delete rv_Nslsig2 ;
      delete rv_Nslsig3 ;
      delete rv_Nslsig4 ;
      delete rv_Nslsig5 ;
      delete rv_Nslsb1 ;
      delete rv_Nslsb2 ;
      delete rv_Nslsb3 ;
      delete rv_Nslsb4 ;
      delete rv_Nslsb5 ;
      delete rv_Nslmsb1 ;
      delete rv_Nslmsb2 ;
      delete rv_Nslmsb3 ;
      delete rv_Nslmsb4 ;
      delete rv_Nslmsb5 ;
      delete rv_mu_ttbar_sig ;
      delete rv_mu_qcd_sig ;
      delete rv_mu_ttbar_sb ;
      delete rv_mu_qcd_sb ;
      delete rv_mu_ew_sig    ;
      delete rv_mu_susy_sig  ;
      delete rv_mu_ew_sb1 ;
      delete rv_mu_ew_sb2 ;
      delete rv_mu_ew_sb3 ;
      delete rv_mu_ew_sb4 ;
      delete rv_mu_ew_sb5 ;

      delete rv_mu_susy_sb1 ;
      delete rv_mu_susy_sb2 ;
      delete rv_mu_susy_sb3 ;
      delete rv_mu_susy_sb4 ;
      delete rv_mu_susy_sb5 ;
      delete rv_mu_ttbar_a ;
      delete rv_mu_qcd_a   ;
      delete rv_mu_ew_a    ;
      delete rv_mu_susy_a  ;
      delete rv_mu_ttbar_d ;
      delete rv_mu_qcd_d   ;
      delete rv_mu_ew_d    ;
      delete rv_mu_susy_d  ;
      delete rv_mu_qcd_lsb1 ;
      delete rv_mu_qcd_lsb2 ;
      delete rv_mu_qcd_lsb3 ;
      delete rv_mu_qcd_lsb4 ;
      delete rv_mu_qcd_lsb5 ;
      delete rv_mu_qcdmc_sig ;
      delete rv_mu_qcdmc_sb  ;
      delete rv_mu_qcdmc_a   ;
      delete rv_mu_qcdmc_d   ;

      delete rv_mu_sl_ttbar_sig1 ;
      delete rv_mu_sl_ttbar_sig2 ;
      delete rv_mu_sl_ttbar_sig3 ;
      delete rv_mu_sl_ttbar_sig4 ;
      delete rv_mu_sl_ttbar_sig5 ;
      delete rv_mu_sl_ttbar_sb1 ;
      delete rv_mu_sl_ttbar_sb2 ;
      delete rv_mu_sl_ttbar_sb3 ;
      delete rv_mu_sl_ttbar_sb4 ;
      delete rv_mu_sl_ttbar_sb5 ;
      delete rv_mu_sl_ttbar_msb1 ;
      delete rv_mu_sl_ttbar_msb2 ;
      delete rv_mu_sl_ttbar_msb3 ;
      delete rv_mu_sl_ttbar_msb4 ;
      delete rv_mu_sl_ttbar_msb5 ;
      delete rv_mu_sl_ew_sig1 ;
      delete rv_mu_sl_ew_sig2 ;
      delete rv_mu_sl_ew_sig3 ;
      delete rv_mu_sl_ew_sig4 ;
      delete rv_mu_sl_ew_sig5 ;
      delete rv_mu_sl_ew_sb1 ;
      delete rv_mu_sl_ew_sb2 ;
      delete rv_mu_sl_ew_sb3 ;
      delete rv_mu_sl_ew_sb4 ;
      delete rv_mu_sl_ew_sb5 ;
      delete rv_mu_sl_ew_msb1 ;
      delete rv_mu_sl_ew_msb2 ;
      delete rv_mu_sl_ew_msb3 ;
      delete rv_mu_sl_ew_msb4 ;
      delete rv_mu_sl_ew_msb5 ;

      delete rv_mu_sl_susy_sig1 ;
      delete rv_mu_sl_susy_sig2 ;
      delete rv_mu_sl_susy_sig3;
      delete rv_mu_sl_susy_sig4 ;
      delete rv_mu_sl_susy_sig5 ;
      delete rv_mu_sl_susy_sb1 ;
      delete rv_mu_sl_susy_sb2 ;
      delete rv_mu_sl_susy_sb3 ;
      delete rv_mu_sl_susy_sb4 ;
      delete rv_mu_sl_susy_sb5 ;
      delete rv_mu_sl_susy_msb1 ;
      delete rv_mu_sl_susy_msb2 ;
      delete rv_mu_sl_susy_msb3 ;
      delete rv_mu_sl_susy_msb4 ;
      delete rv_mu_sl_susy_msb5 ;

      delete rv_mu_sl_ttbar_sb ;
      delete rv_mu_sl_ttbar_sig ;
      delete rv_mu_qcd_lsb ;

      delete rv_f_qcd_lsb1 ;
      delete rv_f_qcd_lsb2 ;
      delete rv_f_qcd_lsb3 ;
      delete rv_f_qcd_lsb4 ;
      delete rv_f_qcd_lsb5 ;
      delete rv_mu_qcd_sb1 ;
      delete rv_mu_qcd_sb2 ;
      delete rv_mu_qcd_sb3 ;
      delete rv_mu_qcd_sb4 ;
      delete rv_mu_qcd_sb5 ;
      delete rv_mu_sl_ttbar1 ;
      delete rv_mu_sl_ttbar2 ;
      delete rv_mu_sl_ttbar3 ;
      delete rv_mu_sl_ttbar4 ;
      delete rv_mu_sl_ttbar5 ;
      delete rv_mu_sl_ttbar ;
      delete rv_f_sl_ttbar1 ;
      delete rv_f_sl_ttbar2 ;
      delete rv_f_sl_ttbar3 ;
      delete rv_f_sl_ttbar4 ;
      delete rv_f_sl_ttbar5 ;
      delete rv_mu_ttbar_sb1 ;
      delete rv_mu_ttbar_sb2 ;
      delete rv_mu_ttbar_sb3 ;
      delete rv_mu_ttbar_sb4 ;
      delete rv_mu_ttbar_sb5 ;
      delete rv_n_sig ;
      delete rv_n_sb1 ;
      delete rv_n_sb2 ;
      delete rv_n_sb3 ;
      delete rv_n_sb4 ;
      delete rv_n_sb5 ;
      delete rv_n_a ;
      delete rv_n_d ;
      delete rv_n_lsb1 ;
      delete rv_n_lsb2 ;
      delete rv_n_lsb3 ;
      delete rv_n_lsb4 ;
      delete rv_n_lsb5 ;
      delete rv_n_sl_sig1 ;
      delete rv_n_sl_sig2 ;
      delete rv_n_sl_sig3 ;
      delete rv_n_sl_sig4 ;
      delete rv_n_sl_sig5 ;
      delete rv_n_sl_sb1 ;
      delete rv_n_sl_sb2 ;
      delete rv_n_sl_sb3 ;
      delete rv_n_sl_sb4 ;
      delete rv_n_sl_sb5 ;
      delete rv_n_sl_msb1 ;
      delete rv_n_sl_msb2 ;
      delete rv_n_sl_msb3 ;
      delete rv_n_sl_msb4 ;
      delete rv_n_sl_msb5 ;

      delete pdf_Nsig ;
      delete pdf_Na ;
      delete pdf_Nd ;
      delete pdf_Nsb1 ;
      delete pdf_Nsb2 ;
      delete pdf_Nsb3 ;
      delete pdf_Nsb4 ;
      delete pdf_Nsb5 ;
      delete pdf_Nlsb1 ;
      delete pdf_Nlsb2 ;
      delete pdf_Nlsb3 ;
      delete pdf_Nlsb4 ;
      delete pdf_Nlsb5 ;
      delete pdf_Nsl_sig1 ;
      delete pdf_Nsl_sig2 ;
      delete pdf_Nsl_sig3 ;
      delete pdf_Nsl_sig4 ;
      delete pdf_Nsl_sig5 ;
      delete pdf_Nsl_sb1 ;
      delete pdf_Nsl_sb2 ;
      delete pdf_Nsl_sb3 ;
      delete pdf_Nsl_sb4 ;
      delete pdf_Nsl_sb5 ;
      delete pdf_Nsl_msb1 ;
      delete pdf_Nsl_msb2 ;
      delete pdf_Nsl_msb3 ;
      delete pdf_Nsl_msb4 ;
      delete pdf_Nsl_msb5 ;
      delete pdf_Nqcdmc_sig ;
      delete pdf_Nqcdmc_sb  ;
      delete pdf_Nqcdmc_a   ;
      delete pdf_Nqcdmc_d   ;

      delete likelihood ;

      delete workspace ;



   }



  //===================================================================

    bool ra2bRoostatsClass::initialize( const char* infile ) {


       printf( "\n\n Opening input file : %s\n\n", infile ) ;

       FILE* infp ;
       if ( (infp=fopen( infile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
          return false ;
       }

       sprintf( initializeFile, "%s", infile ) ;

       int N3jmBins(0) ; //-- number of 3-jet mass bins.

       int Nsig(0), Na(0), Nd(0) ; //-- data counts in signal region, A, and D.
       int Nsb1(0), Nsb2(0), Nsb3(0), Nsb4(0), Nsb5(0) ; //-- data counts in 5 3-jet mass bins of SB.
       int Nlsb1(0), Nlsb2(0), Nlsb3(0), Nlsb4(0), Nlsb5(0) ; //-- data counts in 5 3-jet mass bins of LSB.
       int Nslsig1(0), Nslsig2(0), Nslsig3(0), Nslsig4(0), Nslsig5(0) ; //-- data counts in 5 3-jet mass bins of SL, SIG.
       int Nslsb1(0), Nslsb2(0), Nslsb3(0), Nslsb4(0), Nslsb5(0) ; ; //-- data counts in 5 3-jet mass bins of SL, SB.
       int Nslmsb1(0), Nslmsb2(0), Nslmsb3(0), Nslmsb4(0), Nslmsb5(0) ; //-- data counts in 5 3-jet mass bins of SL, MSB.

       float Nqcdmcsig(0.), Nqcdmcsb(0.), Nqcdmca(0.), Nqcdmcd(0.) ; //-- QCD MC counts in SIG, SB, A, and D.
       float Nqcdmcsigerr(0.), Nqcdmcsberr(0.), Nqcdmcaerr(0.), Nqcdmcderr(0.) ; //-- QCD MC uncertainties in SIG, SB, A, and D.
       float Nqcdmcslsig(0.), Nqcdmcslsb(0.), Nqcdmcslmsb(0.) ; //-- QCD MC counts in SL; SIG, LSB, and MSB.
       float Nttbarmcsig(0.), Nttbarmcsb(0.), Nttbarmca(0.), Nttbarmcd(0.) ; //-- ttbar MC counts in SIG, SB, A, and D.
       float Nttbarmcslsig(0.), Nttbarmcslsb(0.), Nttbarmcslmsb(0.) ; //-- ttbar MC counts in SL; SIG, SB, and MSB.
       float Newmcsig(0.), Newmca(0.), Newmcd(0.) ; //-- EW MC counts in SIG, A, and D.
       float Newmcsb1(0.), Newmcsb2(0.), Newmcsb3(0.), Newmcsb4(0.), Newmcsb5(0.) ; //-- EW MC counts in 5 3-jet mass bins of SB.
       float Newmcslsig1(0.), Newmcslsig2(0.), Newmcslsig3(0.), Newmcslsig4(0.), Newmcslsig5(0.) ; //-- EW MC, SL counts in SIG, bins of 3-jet mass
       float Newmcslsb1(0.), Newmcslsb2(0.), Newmcslsb3(0.), Newmcslsb4(0.), Newmcslsb5(0.) ; //-- EW MC, SL counts in SB, bins of 3-jet mass
       float Newmcslmsb1(0.), Newmcslmsb2(0.), Newmcslmsb3(0.), Newmcslmsb4(0.), Newmcslmsb5(0.) ; //-- EW MC, SL counts in MSB, bins of 3-jet mass

       float Nsusymcsig(0.), Nsusymca(0.), Nsusymcd(0.) ; //-- SUSY MC counts in SIG, A, and D.
       float Nsusymcsb1(0.), Nsusymcsb2(0.), Nsusymcsb3(0.), Nsusymcsb4(0.), Nsusymcsb5(0.) ; //-- SUSY MC counts in 5 3-jet mass bins of SB.
       float Nsusymcslsig1(0.), Nsusymcslsig2(0.), Nsusymcslsig3(0.), Nsusymcslsig4(0.), Nsusymcslsig5(0.) ; //-- SUSY MC counts in 5 3-jet mass bins of SL, SIG.
       float Nsusymcslsb1(0.), Nsusymcslsb2(0.), Nsusymcslsb3(0.), Nsusymcslsb4(0.), Nsusymcslsb5(0.) ; //-- SUSY MC counts in 5 3-jet mass bins of SL, SB.
       float Nsusymcslmsb1(0.), Nsusymcslmsb2(0.), Nsusymcslmsb3(0.), Nsusymcslmsb4(0.), Nsusymcslmsb5(0.) ; //-- SUSY MC counts in 5 3-jet mass bins of SL, MSB.

       //--- read in description line.
       printf("\n\n") ;
       char c(0) ;
       while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
       printf("\n\n") ;


       char label[1000] ;

      //--- Inputs generated with gen_roostats_input.c
      //    The order here must be consistent with the order there!

       fscanf( infp, "%s %d", label, &N3jmBins ) ;             printf( "%s %d\n", label, N3jmBins ) ;                
       fscanf( infp, "%s %d", label, &Nsig ) ;                 printf( "%s %d\n", label, Nsig ) ;         
       fscanf( infp, "%s %d", label, &Na ) ;                   printf( "%s %d\n", label, Na ) ;           
       fscanf( infp, "%s %d", label, &Nd ) ;                   printf( "%s %d\n", label, Nd ) ;           
       fscanf( infp, "%s %d", label, &Nsb1 ) ;                 printf( "%s %d\n", label, Nsb1 ) ;         
       fscanf( infp, "%s %d", label, &Nsb2 ) ;                 printf( "%s %d\n", label, Nsb2 ) ;         
       fscanf( infp, "%s %d", label, &Nsb3 ) ;                 printf( "%s %d\n", label, Nsb3 ) ;         
       fscanf( infp, "%s %d", label, &Nsb4 ) ;                 printf( "%s %d\n", label, Nsb4 ) ;         
       fscanf( infp, "%s %d", label, &Nsb5 ) ;                 printf( "%s %d\n", label, Nsb5 ) ;         
       fscanf( infp, "%s %d", label, &Nlsb1 ) ;                printf( "%s %d\n", label, Nlsb1 ) ;        
       fscanf( infp, "%s %d", label, &Nlsb2 ) ;                printf( "%s %d\n", label, Nlsb2 ) ;        
       fscanf( infp, "%s %d", label, &Nlsb3 ) ;                printf( "%s %d\n", label, Nlsb3 ) ;        
       fscanf( infp, "%s %d", label, &Nlsb4 ) ;                printf( "%s %d\n", label, Nlsb4 ) ;        
       fscanf( infp, "%s %d", label, &Nlsb5 ) ;                printf( "%s %d\n", label, Nlsb5 ) ;        
       fscanf( infp, "%s %d", label, &Nslsig1 ) ;              printf( "%s %d\n", label, Nslsig1 ) ;      
       fscanf( infp, "%s %d", label, &Nslsig2 ) ;              printf( "%s %d\n", label, Nslsig2 ) ;      
       fscanf( infp, "%s %d", label, &Nslsig3 ) ;              printf( "%s %d\n", label, Nslsig3 ) ;      
       fscanf( infp, "%s %d", label, &Nslsig4 ) ;              printf( "%s %d\n", label, Nslsig4 ) ;      
       fscanf( infp, "%s %d", label, &Nslsig5 ) ;              printf( "%s %d\n", label, Nslsig5 ) ;      
       fscanf( infp, "%s %d", label, &Nslsb1 ) ;               printf( "%s %d\n", label, Nslsb1 ) ;       
       fscanf( infp, "%s %d", label, &Nslsb2 ) ;               printf( "%s %d\n", label, Nslsb2 ) ;       
       fscanf( infp, "%s %d", label, &Nslsb3 ) ;               printf( "%s %d\n", label, Nslsb3 ) ;       
       fscanf( infp, "%s %d", label, &Nslsb4 ) ;               printf( "%s %d\n", label, Nslsb4 ) ;       
       fscanf( infp, "%s %d", label, &Nslsb5 ) ;               printf( "%s %d\n", label, Nslsb5 ) ;       
       fscanf( infp, "%s %d", label, &Nslmsb1 ) ;              printf( "%s %d\n", label, Nslmsb1 ) ;      
       fscanf( infp, "%s %d", label, &Nslmsb2 ) ;              printf( "%s %d\n", label, Nslmsb2 ) ;      
       fscanf( infp, "%s %d", label, &Nslmsb3 ) ;              printf( "%s %d\n", label, Nslmsb3 ) ;      
       fscanf( infp, "%s %d", label, &Nslmsb4 ) ;              printf( "%s %d\n", label, Nslmsb4 ) ;      
       fscanf( infp, "%s %d", label, &Nslmsb5 ) ;              printf( "%s %d\n", label, Nslmsb5 ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcsig ) ;            printf( "%s %g\n", label, Nqcdmcsig ) ;    
       fscanf( infp, "%s %g", label, &Nqcdmcsigerr ) ;         printf( "%s %g\n", label, Nqcdmcsigerr ) ;    
       fscanf( infp, "%s %g", label, &Nqcdmcsb ) ;             printf( "%s %g\n", label, Nqcdmcsb ) ;     
       fscanf( infp, "%s %g", label, &Nqcdmcsberr ) ;          printf( "%s %g\n", label, Nqcdmcsberr ) ;     
       fscanf( infp, "%s %g", label, &Nqcdmca ) ;              printf( "%s %g\n", label, Nqcdmca ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcaerr ) ;           printf( "%s %g\n", label, Nqcdmcaerr ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcd ) ;              printf( "%s %g\n", label, Nqcdmcd ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcderr ) ;           printf( "%s %g\n", label, Nqcdmcderr ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcslsig ) ;          printf( "%s %g\n", label, Nqcdmcslsig ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcslsb ) ;           printf( "%s %g\n", label, Nqcdmcslsb ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcslmsb ) ;          printf( "%s %g\n", label, Nqcdmcslmsb ) ;      
       fscanf( infp, "%s %g", label, &Nttbarmcsig ) ;          printf( "%s %g\n", label, Nttbarmcsig ) ;  
       fscanf( infp, "%s %g", label, &Nttbarmcsb ) ;           printf( "%s %g\n", label, Nttbarmcsb ) ;   
       fscanf( infp, "%s %g", label, &Nttbarmca ) ;            printf( "%s %g\n", label, Nttbarmca ) ;    
       fscanf( infp, "%s %g", label, &Nttbarmcd ) ;            printf( "%s %g\n", label, Nttbarmcd ) ;    
       fscanf( infp, "%s %g", label, &Nttbarmcslsig ) ;        printf( "%s %g\n", label, Nttbarmcslsig ) ;      
       fscanf( infp, "%s %g", label, &Nttbarmcslsb ) ;         printf( "%s %g\n", label, Nttbarmcslsb ) ;      
       fscanf( infp, "%s %g", label, &Nttbarmcslmsb ) ;        printf( "%s %g\n", label, Nttbarmcslmsb ) ;      
       fscanf( infp, "%s %g", label, &Newmcsig ) ;             printf( "%s %g\n", label, Newmcsig ) ;     
       fscanf( infp, "%s %g", label, &Newmca ) ;               printf( "%s %g\n", label, Newmca ) ;       
       fscanf( infp, "%s %g", label, &Newmcd ) ;               printf( "%s %g\n", label, Newmcd ) ;       
       fscanf( infp, "%s %g", label, &Newmcsb1 ) ;             printf( "%s %g\n", label, Newmcsb1 ) ;     
       fscanf( infp, "%s %g", label, &Newmcsb2 ) ;             printf( "%s %g\n", label, Newmcsb2 ) ;     
       fscanf( infp, "%s %g", label, &Newmcsb3 ) ;             printf( "%s %g\n", label, Newmcsb3 ) ;     
       fscanf( infp, "%s %g", label, &Newmcsb4 ) ;             printf( "%s %g\n", label, Newmcsb4 ) ;     
       fscanf( infp, "%s %g", label, &Newmcsb5 ) ;             printf( "%s %g\n", label, Newmcsb5 ) ;     
       fscanf( infp, "%s %g", label, &Newmcslsig1 ) ;          printf( "%s %g\n", label, Newmcslsig1 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsig2 ) ;          printf( "%s %g\n", label, Newmcslsig2 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsig3 ) ;          printf( "%s %g\n", label, Newmcslsig3 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsig4 ) ;          printf( "%s %g\n", label, Newmcslsig4 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsig5 ) ;          printf( "%s %g\n", label, Newmcslsig5 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsb1 ) ;          printf( "%s %g\n", label, Newmcslsb1 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsb2 ) ;          printf( "%s %g\n", label, Newmcslsb2 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsb3 ) ;          printf( "%s %g\n", label, Newmcslsb3 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsb4 ) ;          printf( "%s %g\n", label, Newmcslsb4 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsb5 ) ;          printf( "%s %g\n", label, Newmcslsb5 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslmsb1 ) ;          printf( "%s %g\n", label, Newmcslmsb1 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslmsb2 ) ;          printf( "%s %g\n", label, Newmcslmsb2 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslmsb3 ) ;          printf( "%s %g\n", label, Newmcslmsb3 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslmsb4 ) ;          printf( "%s %g\n", label, Newmcslmsb4 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslmsb5 ) ;          printf( "%s %g\n", label, Newmcslmsb5 ) ;      
       fscanf( infp, "%s %g", label, &Nsusymcsig ) ;           printf( "%s %g\n", label, Nsusymcsig ) ;   
       fscanf( infp, "%s %g", label, &Nsusymca ) ;             printf( "%s %g\n", label, Nsusymca ) ;     
       fscanf( infp, "%s %g", label, &Nsusymcd ) ;             printf( "%s %g\n", label, Nsusymcd ) ;     
       fscanf( infp, "%s %g", label, &Nsusymcsb1 ) ;           printf( "%s %g\n", label, Nsusymcsb1 ) ;   
       fscanf( infp, "%s %g", label, &Nsusymcsb2 ) ;           printf( "%s %g\n", label, Nsusymcsb2 ) ;   
       fscanf( infp, "%s %g", label, &Nsusymcsb3 ) ;           printf( "%s %g\n", label, Nsusymcsb3 ) ;   
       fscanf( infp, "%s %g", label, &Nsusymcsb4 ) ;           printf( "%s %g\n", label, Nsusymcsb4 ) ;   
       fscanf( infp, "%s %g", label, &Nsusymcsb5 ) ;           printf( "%s %g\n", label, Nsusymcsb5 ) ;   
       fscanf( infp, "%s %g", label, &Nsusymcslsig1 ) ;        printf( "%s %g\n", label, Nsusymcslsig1 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslsig2 ) ;        printf( "%s %g\n", label, Nsusymcslsig2 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslsig3 ) ;        printf( "%s %g\n", label, Nsusymcslsig3 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslsig4 ) ;        printf( "%s %g\n", label, Nsusymcslsig4 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslsig5 ) ;        printf( "%s %g\n", label, Nsusymcslsig5 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslsb1 ) ;         printf( "%s %g\n", label, Nsusymcslsb1 ) ; 
       fscanf( infp, "%s %g", label, &Nsusymcslsb2 ) ;         printf( "%s %g\n", label, Nsusymcslsb2 ) ; 
       fscanf( infp, "%s %g", label, &Nsusymcslsb3 ) ;         printf( "%s %g\n", label, Nsusymcslsb3 ) ; 
       fscanf( infp, "%s %g", label, &Nsusymcslsb4 ) ;         printf( "%s %g\n", label, Nsusymcslsb4 ) ; 
       fscanf( infp, "%s %g", label, &Nsusymcslsb5 ) ;         printf( "%s %g\n", label, Nsusymcslsb5 ) ; 
       fscanf( infp, "%s %g", label, &Nsusymcslmsb1 ) ;        printf( "%s %g\n", label, Nsusymcslmsb1 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslmsb2 ) ;        printf( "%s %g\n", label, Nsusymcslmsb2 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslmsb3 ) ;        printf( "%s %g\n", label, Nsusymcslmsb3 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslmsb4 ) ;        printf( "%s %g\n", label, Nsusymcslmsb4 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslmsb5 ) ;        printf( "%s %g\n", label, Nsusymcslmsb5 ) ;


       printf("\n Done reading in %s\n\n", infile ) ;
       fclose( infp ) ;






       //--- Print out a nice summary of the inputs.


       float Nsmsig = Nttbarmcsig + Nqcdmcsig + Newmcsig ;

       int   Nsb = Nsb1+Nsb2+Nsb3+Nsb4+Nsb5 ;
       float Newmcsb =Newmcsb1+Newmcsb2+Newmcsb3+Newmcsb4+Newmcsb5;
       float Nsusymcsb =Nsusymcsb1+Nsusymcsb2+Nsusymcsb3+Nsusymcsb4+Nsusymcsb5;
       float Nsmsb = Nttbarmcsb + Nqcdmcsb + Newmcsb ;

       float Nsma = Nttbarmca + Nqcdmca + Newmca ;
       float Nsmd = Nttbarmcd + Nqcdmcd + Newmcd ;

       float Newmcslsig = Newmcslsig1 + Newmcslsig2 + Newmcslsig3 + Newmcslsig4 + Newmcslsig5 ;
       float Nsmslsig = Nttbarmcslsig + Nqcdmcslsig + Newmcslsig ;
       int   Nslsig = Nslsig1+Nslsig2+Nslsig3+Nslsig4+Nslsig5 ;
       float Nsusymcslsig =Nsusymcslsig1+Nsusymcslsig2+Nsusymcslsig3+Nsusymcslsig4+Nsusymcslsig5;

       float Newmcslsb = Newmcslsb1 + Newmcslsb2 + Newmcslsb3 + Newmcslsb4 + Newmcslsb5 ;
       float Nsmslsb = Nttbarmcslsb + Nqcdmcslsb + Newmcslsb ;
       int   Nslsb = Nslsb1+Nslsb2+Nslsb3+Nslsb4+Nslsb5 ;
       float Nsusymcslsb =Nsusymcslsb1+Nsusymcslsb2+Nsusymcslsb3+Nsusymcslsb4+Nsusymcslsb5;

       int   Nslmsb = Nslmsb1+Nslmsb2+Nslmsb3+Nslmsb4+Nslmsb5 ;
       float Nsusymcslmsb =Nsusymcslmsb1+Nsusymcslmsb2+Nsusymcslmsb3+Nsusymcslmsb4+Nsusymcslmsb5;

       float Nttbarmcslmsbsbsig = Nttbarmcslmsb+Nttbarmcslsb+Nttbarmcslsig ;
       float Nqcdmcslmsbsbsig = Nqcdmcslmsb+Nqcdmcslsb+Nqcdmcslsig ;
       float Newmcslmsb = Newmcslmsb1 + Newmcslmsb2 + Newmcslmsb3 + Newmcslmsb4 + Newmcslmsb5 ;
       float Newmcslmsbsbsig = Newmcslmsb+Newmcslsb+Newmcslsig ;
       float Nsmmcslmsbsbsig = Nttbarmcslmsbsbsig + Nqcdmcslmsbsbsig + Newmcslmsbsbsig ;
       int   Nslmsbsbsig = Nslmsb+Nslsb+Nslsig ;
       float Nsusymcslmsbsbsig = Nsusymcslmsb+Nsusymcslsb+Nsusymcslsig ;

       float Rttmc = Nttbarmcsig / Nttbarmcsb ;
       float Rslmc = Nsmslsig / Nsmslsb ;
       float Rsldata = (1.0*Nslsig) / (1.0*Nslsb) ;
       float Rsldataerr = Rsldata*sqrt(1.0/Nslsig + 1.0/Nslsb) ;

       float Rqcdmctruth = Nqcdmcd / Nqcdmca ;
       float Rqcddata    = (Nd-Nttbarmcd-Newmcd)/(Na-Nttbarmca-Newmca) ;
       float Rqcddataerr = sqrt( pow(1./(Na-Nttbarmca-Newmca),2)*Nd + pow( ((Nd-Nttbarmcd-Newmcd)/pow((Na-Nttbarmca-Newmca),2)),2)*Na ) ;

       float frsigttbar = Nttbarmcsig / Nsmsig ;
       float frsigqcd   = Nqcdmcsig / Nsmsig ;
       float frsigew    = Newmcsig / Nsmsig ;

       float frsbttbar = Nttbarmcsb / Nsmsb ;
       float frsbqcd   = Nqcdmcsb / Nsmsb ;
       float frsbew    = Newmcsb / Nsmsb ;

       float frslttbar = Nttbarmcslmsbsbsig / Nsmmcslmsbsbsig ;
       float frslqcd   = Nqcdmcslmsbsbsig / Nsmmcslmsbsbsig ;
       float frslew    = Newmcslmsbsbsig / Nsmmcslmsbsbsig ;



       printf("\n\n\n") ;

       printf("                        ttbar    qcd       EW    all SM      data   SUSY\n") ;
       printf("----------------------------------------------------------------------------\n") ;
       printf("  Signal box     :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmcsig, Nqcdmcsig, Newmcsig, Nsmsig, Nsig, Nsusymcsig ) ;
       printf("  Sideband       :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmcsb, Nqcdmcsb, Newmcsb, Nsmsb, Nsb, Nsusymcsb ) ;
       printf("\n") ;
       printf("  A              :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmca, Nqcdmca, Newmca, Nsma, Na, Nsusymca ) ;
       printf("  D              :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmcd, Nqcdmcd, Newmcd, Nsmd, Nd, Nsusymcd ) ;
       printf("\n") ;
       printf("  SL, SIG        :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmcslsig, Nqcdmcslsig, Newmcslsig, Nsmslsig, Nslsig, Nsusymcslsig ) ;
       printf("  SL, SB         :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmcslsb, Nqcdmcslsb, Newmcslsb, Nsmslsb, Nslsb, Nsusymcslsb ) ;
       printf("\n") ;
       printf("  SL, MSB+SB+SIG :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmcslmsbsbsig, Nqcdmcslmsbsbsig, Newmcslmsbsbsig, Nsmmcslmsbsbsig, Nslmsbsbsig, Nsusymcslmsbsbsig ) ;
       printf("----------------------------------------------------------------------------\n") ;
       printf("\n") ;
       printf("    Rttmc     = %5.1f / %5.1f = %5.3f\n", Nttbarmcsig, Nttbarmcsb, Rttmc ) ;
       printf("    Rsl,mc    = %5.1f / %5.1f = %5.3f\n", Nsmslsig, Nsmslsb, Rslmc ) ;
       printf("    Rsl,data  = %5d / %5d = %5.3f +/- %5.3f\n", Nslsig, Nslsb, Rsldata, Rsldataerr ) ;
       printf("\n") ;
       printf("    Rqcd,mc   = %5.1f / %6.1f = %5.3f\n", Nqcdmcd, Nqcdmca, Rqcdmctruth ) ;
       printf("    Rqcd,data = %5.1f / %6.1f = %5.3f +/- %5.3f\n", (Nd-Nttbarmcd-Newmcd), (Na-Nttbarmca-Newmca), Rqcddata, Rqcddataerr ) ;
       printf("\n") ;
       printf("                              ttbar    qcd     EW\n") ;
       printf("-----------------------------------------------------\n") ;
       printf("   SIG region  fractions :   %5.2f   %5.2f   %5.2f\n", frsigttbar, frsigqcd, frsigew ) ;
       printf("   SB  region  fractions :   %5.2f   %5.2f   %5.2f\n", frsbttbar, frsbqcd, frsbew ) ;
       printf("   SL (MET>50) fractions :   %5.2f   %5.2f   %5.2f\n", frslttbar, frslqcd, frslew ) ;
       printf("\n\n\n") ;



       printf(" --- Defining observables.\n" ) ;

      //-- data counts in signal region, A, and D.

      rv_Nsig = new RooRealVar( "Nsig", "Nsig", 0.5, 400. ) ;
      rv_Na = new RooRealVar( "Na", "Na", 0.5, 1000000. ) ;
      rv_Nd = new RooRealVar( "Nd", "Nd", 0.5, 1000000. ) ;

      rv_Nsig -> setVal( Nsig ) ;
      rv_Na -> setVal( Na ) ;
      rv_Nd -> setVal( Nd ) ;

      //-- data counts in 5 3-jet mass bins of SB.

      rv_Nsb1 = new RooRealVar( "Nsb1", "Nsb1", 0.5, 1000000. ) ;
      rv_Nsb2 = new RooRealVar( "Nsb2", "Nsb2", 0.5, 1000000. ) ;
      rv_Nsb3 = new RooRealVar( "Nsb3", "Nsb3", 0.5, 1000000. ) ;
      rv_Nsb4 = new RooRealVar( "Nsb4", "Nsb4", 0.5, 1000000. ) ;
      rv_Nsb5 = new RooRealVar( "Nsb5", "Nsb5", 0.5, 1000000. ) ;

      rv_Nsb1 -> setVal( Nsb1 ) ;
      rv_Nsb2 -> setVal( Nsb2 ) ;
      rv_Nsb3 -> setVal( Nsb3 ) ;
      rv_Nsb4 -> setVal( Nsb4 ) ;
      rv_Nsb5 -> setVal( Nsb5 ) ;

      //-- data counts in 5 3-jet mass bins of LSB.

      rv_Nlsb1 = new RooRealVar( "Nlsb1", "Nlsb1", 0.5, 10000000. ) ;
      rv_Nlsb2 = new RooRealVar( "Nlsb2", "Nlsb2", 0.5, 10000000. ) ;
      rv_Nlsb3 = new RooRealVar( "Nlsb3", "Nlsb3", 0.5, 10000000. ) ;
      rv_Nlsb4 = new RooRealVar( "Nlsb4", "Nlsb4", 0.5, 10000000. ) ;
      rv_Nlsb5 = new RooRealVar( "Nlsb5", "Nlsb5", 0.5, 10000000. ) ;

      rv_Nlsb1 -> setVal( Nlsb1 ) ;
      rv_Nlsb2 -> setVal( Nlsb2 ) ;
      rv_Nlsb3 -> setVal( Nlsb3 ) ;
      rv_Nlsb4 -> setVal( Nlsb4 ) ;
      rv_Nlsb5 -> setVal( Nlsb5 ) ;


      //-- data counts in 5 3-jet mass bins of SL, SIG.

      rv_Nslsig1 = new RooRealVar( "Nslsig1", "Nslsig1", 0.5, 1000000. ) ;
      rv_Nslsig2 = new RooRealVar( "Nslsig2", "Nslsig2", 0.5, 1000000. ) ;
      rv_Nslsig3 = new RooRealVar( "Nslsig3", "Nslsig3", 0.5, 1000000. ) ;
      rv_Nslsig4 = new RooRealVar( "Nslsig4", "Nslsig4", 0.5, 1000000. ) ;
      rv_Nslsig5 = new RooRealVar( "Nslsig5", "Nslsig5", 0.5, 1000000. ) ;

      rv_Nslsig1 -> setVal( Nslsig1 ) ;
      rv_Nslsig2 -> setVal( Nslsig2 ) ;
      rv_Nslsig3 -> setVal( Nslsig3 ) ;
      rv_Nslsig4 -> setVal( Nslsig4 ) ;
      rv_Nslsig5 -> setVal( Nslsig5 ) ;


      //-- data counts in 5 3-jet mass bins of SL, SB.

      rv_Nslsb1 = new RooRealVar( "Nslsb1", "Nslsb1", 0.5, 1000000. ) ;
      rv_Nslsb2 = new RooRealVar( "Nslsb2", "Nslsb2", 0.5, 1000000. ) ;
      rv_Nslsb3 = new RooRealVar( "Nslsb3", "Nslsb3", 0.5, 1000000. ) ;
      rv_Nslsb4 = new RooRealVar( "Nslsb4", "Nslsb4", 0.5, 1000000. ) ;
      rv_Nslsb5 = new RooRealVar( "Nslsb5", "Nslsb5", 0.5, 1000000. ) ;

      rv_Nslsb1 -> setVal( Nslsb1 ) ;
      rv_Nslsb2 -> setVal( Nslsb2 ) ;
      rv_Nslsb3 -> setVal( Nslsb3 ) ;
      rv_Nslsb4 -> setVal( Nslsb4 ) ;
      rv_Nslsb5 -> setVal( Nslsb5 ) ;


      //-- data counts in 5 3-jet mass bins of SL, MSB.

      rv_Nslmsb1 = new RooRealVar( "Nslmsb1", "Nslmsb1", 0.5, 1000000. ) ;
      rv_Nslmsb2 = new RooRealVar( "Nslmsb2", "Nslmsb2", 0.5, 1000000. ) ;
      rv_Nslmsb3 = new RooRealVar( "Nslmsb3", "Nslmsb3", 0.5, 1000000. ) ;
      rv_Nslmsb4 = new RooRealVar( "Nslmsb4", "Nslmsb4", 0.5, 1000000. ) ;
      rv_Nslmsb5 = new RooRealVar( "Nslmsb5", "Nslmsb5", 0.5, 1000000. ) ;

      rv_Nslmsb1 -> setVal( Nslmsb1 ) ;
      rv_Nslmsb2 -> setVal( Nslmsb2 ) ;
      rv_Nslmsb3 -> setVal( Nslmsb3 ) ;
      rv_Nslmsb4 -> setVal( Nslmsb4 ) ;
      rv_Nslmsb5 -> setVal( Nslmsb5 ) ;





      //-- QCD MC counts

      rv_Nqcdmca   = new RooRealVar( "Nqcdmca", "Nqcdmca", 0.5, 1000000. ) ;
      rv_Nqcdmcd   = new RooRealVar( "Nqcdmcd", "Nqcdmcd", 0.5, 1000000. ) ;
      rv_Nqcdmcsig = new RooRealVar( "Nqcdmcsig", "Nqcdmcsig", 0.5, 1000000. ) ;
      rv_Nqcdmcsb  = new RooRealVar( "Nqcdmcsb", "Nqcdmcsb", 0.5, 1000000. ) ;

      rv_Nqcdmca   -> setVal( Nqcdmca ) ;
      rv_Nqcdmcd   -> setVal( Nqcdmcd ) ;
      rv_Nqcdmcsig -> setVal( Nqcdmcsig ) ;
      rv_Nqcdmcsb  -> setVal( Nqcdmcsb ) ;







      //++++++++ Parameters of the likelihood

       printf(" --- Defining parameters.\n" ) ;

    //===================================================================================================
     //-- Use these lines to define things in terms of the SIG parameters.
      if ( useSigBgVars ) {

         rv_mu_ttbar_sig = new RooRealVar( "mu_ttbar_sig", "mu_ttbar_sig", 0.0, 100. ) ;
         rv_mu_qcd_sig   = new RooRealVar( "mu_qcd_sig"  , "mu_qcd_sig"  , 0.0, 100. ) ;

         rrv_mu_ttbar_sig = ((RooRealVar*) rv_mu_ttbar_sig) ;
         rrv_mu_qcd_sig = ((RooRealVar*) rv_mu_qcd_sig) ;

         rrv_mu_ttbar_sig -> setVal( Nttbarmcsig ) ;  //-- this is a starting value only.
         rrv_mu_qcd_sig   -> setVal( Nqcdmcsig ) ;  //-- this is a starting value only.

      } else {
    //===================================================================================================
     //-- Use these lines to define things in terms of the SB parameters.

         rv_mu_ttbar_sb = new RooRealVar( "mu_ttbar_sb", "mu_ttbar_sb", 0.0, 300. ) ;
         rv_mu_qcd_sb   = new RooRealVar( "mu_qcd_sb"  , "mu_qcd_sb"  , 0.0, 300. ) ;

         rrv_mu_ttbar_sb = ((RooRealVar*) rv_mu_ttbar_sb) ;
         rrv_mu_qcd_sb = ((RooRealVar*) rv_mu_qcd_sb) ;

         rrv_mu_ttbar_sb -> setVal( Nttbarmcsb ) ;  //-- this is a starting value only.
         rrv_mu_qcd_sb   -> setVal( Nqcdmcsb ) ;  //-- this is a starting value only.

      }
    //===================================================================================================


     //-- Counts in SIG, signal selection

      rv_mu_ew_sig    = new RooRealVar( "mu_ew_sig"   , "mu_ew_sig"   , 0.0, 10000. ) ;
      rv_mu_susy_sig  = new RooRealVar( "mu_susy_sig" , "mu_susy_sig" , 0.0, 100. ) ;

      rv_mu_ew_sig    -> setVal( Newmcsig ) ;
      rv_mu_susy_sig  -> setVal( 0. ) ;  //-- this is a starting value only.

      rv_mu_ew_sig -> setConstant(kTRUE) ;

     //-- Counts in SB, bins of 3-jet mass, EW, signal selection

      rv_mu_ew_sb1 = new RooRealVar( "mu_ew_sb1", "mu_ew_sb1", 0.0, 10000. ) ;
      rv_mu_ew_sb2 = new RooRealVar( "mu_ew_sb2", "mu_ew_sb2", 0.0, 10000. ) ;
      rv_mu_ew_sb3 = new RooRealVar( "mu_ew_sb3", "mu_ew_sb3", 0.0, 10000. ) ;
      rv_mu_ew_sb4 = new RooRealVar( "mu_ew_sb4", "mu_ew_sb4", 0.0, 10000. ) ;
      rv_mu_ew_sb5 = new RooRealVar( "mu_ew_sb5", "mu_ew_sb5", 0.0, 10000. ) ;

      rv_mu_ew_sb1 -> setVal( Newmcsb1 ) ;
      rv_mu_ew_sb2 -> setVal( Newmcsb2 ) ;
      rv_mu_ew_sb3 -> setVal( Newmcsb3 ) ;
      rv_mu_ew_sb4 -> setVal( Newmcsb4 ) ;
      rv_mu_ew_sb5 -> setVal( Newmcsb5 ) ;

      rv_mu_ew_sb1 -> setConstant(kTRUE) ;
      rv_mu_ew_sb2 -> setConstant(kTRUE) ;
      rv_mu_ew_sb3 -> setConstant(kTRUE) ;
      rv_mu_ew_sb4 -> setConstant(kTRUE) ;
      rv_mu_ew_sb5 -> setConstant(kTRUE) ;

     //-- Counts in SB, bins of 3-jet mass, SUSY, signal selection

      rv_mu_susy_sb1 = new RooRealVar( "mu_susy_sb1", "mu_susy_sb1", 0.0, 10000. ) ;
      rv_mu_susy_sb2 = new RooRealVar( "mu_susy_sb2", "mu_susy_sb2", 0.0, 10000. ) ;
      rv_mu_susy_sb3 = new RooRealVar( "mu_susy_sb3", "mu_susy_sb3", 0.0, 10000. ) ;
      rv_mu_susy_sb4 = new RooRealVar( "mu_susy_sb4", "mu_susy_sb4", 0.0, 10000. ) ;
      rv_mu_susy_sb5 = new RooRealVar( "mu_susy_sb5", "mu_susy_sb5", 0.0, 10000. ) ;

      rv_mu_susy_sb1 -> setVal( 0. ) ;
      rv_mu_susy_sb2 -> setVal( 0. ) ;
      rv_mu_susy_sb3 -> setVal( 0. ) ;
      rv_mu_susy_sb4 -> setVal( 0. ) ;
      rv_mu_susy_sb5 -> setVal( 0. ) ;

      rv_mu_susy_sb1 -> setConstant(kTRUE) ;
      rv_mu_susy_sb2 -> setConstant(kTRUE) ;
      rv_mu_susy_sb3 -> setConstant(kTRUE) ;
      rv_mu_susy_sb4 -> setConstant(kTRUE) ;
      rv_mu_susy_sb5 -> setConstant(kTRUE) ;

     //-- Counts in A, signal selection

      rv_mu_ttbar_a = new RooRealVar( "mu_ttbar_a", "mu_ttbar_a", 0.0, 10000. ) ;
      rv_mu_qcd_a   = new RooRealVar( "mu_qcd_a"  , "mu_qcd_a"  , 0.0, 10000. ) ;
      rv_mu_ew_a    = new RooRealVar( "mu_ew_a"   , "mu_ew_a"   , 0.0, 10000. ) ;
      rv_mu_susy_a  = new RooRealVar( "mu_susy_a" , "mu_susy_a" , 0.0, 10000. ) ;

      rv_mu_ttbar_a -> setVal( Nttbarmca ) ;
      rv_mu_qcd_a   -> setVal( Na - Nttbarmca - Newmca ) ; //-- this is a starting value only.
      rv_mu_ew_a    -> setVal( Newmca ) ;
      rv_mu_susy_a  -> setVal( 0. ) ;

      rv_mu_susy_a  ->setConstant(kTRUE) ;
      rv_mu_ew_a    ->setConstant(kTRUE) ;
      rv_mu_ttbar_a ->setConstant(kTRUE) ;


     //-- Counts in D, signal selection

      rv_mu_ttbar_d = new RooRealVar( "mu_ttbar_d", "mu_ttbar_d", 0.0, 10000. ) ;
      rv_mu_qcd_d   = new RooRealVar( "mu_qcd_d"  , "mu_qcd_d"  , 0.0, 10000. ) ;
      rv_mu_ew_d    = new RooRealVar( "mu_ew_d"   , "mu_ew_d"   , 0.0, 10000. ) ;
      rv_mu_susy_d  = new RooRealVar( "mu_susy_d" , "mu_susy_d" , 0.0, 10000. ) ;

      rv_mu_ttbar_d -> setVal( Nttbarmcd ) ;
      rv_mu_qcd_d   -> setVal( Nd - Nttbarmcd - Newmcd ) ; //-- this is a starting value only.
      rv_mu_ew_d    -> setVal( Newmcd ) ;
      rv_mu_susy_d  -> setVal( 0. ) ;

      rv_mu_ew_d    -> setConstant(kTRUE) ;
      rv_mu_ttbar_d -> setConstant(kTRUE) ;
      rv_mu_susy_d  -> setConstant(kTRUE) ;

     //-- Counts in LSB, bins of 3-jet mass, qcd, signal selection

      rv_mu_qcd_lsb1 = new RooRealVar( "mu_qcd_lsb1", "mu_qcd_lsb1", 0.0, 1000000. ) ;
      rv_mu_qcd_lsb2 = new RooRealVar( "mu_qcd_lsb2", "mu_qcd_lsb2", 0.0, 1000000. ) ;
      rv_mu_qcd_lsb3 = new RooRealVar( "mu_qcd_lsb3", "mu_qcd_lsb3", 0.0, 1000000. ) ;
      rv_mu_qcd_lsb4 = new RooRealVar( "mu_qcd_lsb4", "mu_qcd_lsb4", 0.0, 1000000. ) ;
      rv_mu_qcd_lsb5 = new RooRealVar( "mu_qcd_lsb5", "mu_qcd_lsb5", 0.0, 1000000. ) ;

      rv_mu_qcd_lsb1 -> setVal( Nlsb1 ) ;
      rv_mu_qcd_lsb2 -> setVal( Nlsb2 ) ;
      rv_mu_qcd_lsb3 -> setVal( Nlsb3 ) ;
      rv_mu_qcd_lsb4 -> setVal( Nlsb4 ) ;
      rv_mu_qcd_lsb5 -> setVal( Nlsb5 ) ;


     //-- QCD MC, SIG,SB,A,D counts, signal selection

      rv_mu_qcdmc_sig = new RooRealVar( "mu_qcdmc_sig", "mu_qcdmc_sig", 0.0, 10000. ) ;
      rv_mu_qcdmc_sb  = new RooRealVar( "mu_qcdmc_sb" , "mu_qcdmc_sb" , 0.0, 10000. ) ;
      rv_mu_qcdmc_a   = new RooRealVar( "mu_qcdmc_a"  , "mu_qcdmc_a"  , 0.0, 10000. ) ;
      rv_mu_qcdmc_d   = new RooRealVar( "mu_qcdmc_d"  , "mu_qcdmc_d"  , 0.0, 10000. ) ;

      rv_mu_qcdmc_sig -> setVal( Nqcdmcsig ) ;
      rv_mu_qcdmc_sb  -> setVal( Nqcdmcsb ) ;
      rv_mu_qcdmc_a   -> setVal( Nqcdmca ) ;
      rv_mu_qcdmc_d   -> setVal( Nqcdmcd ) ;

     //-- Single Lepton (SL) counts, ttbar, SIG region, bins of 3-jet mass

      rv_mu_sl_ttbar_sig1 = new RooRealVar( "mu_sl_ttbar_sig1", "mu_sl_ttbar_sig1", 0.0, 10000. ) ;
      rv_mu_sl_ttbar_sig2 = new RooRealVar( "mu_sl_ttbar_sig2", "mu_sl_ttbar_sig2", 0.0, 10000. ) ;
      rv_mu_sl_ttbar_sig3 = new RooRealVar( "mu_sl_ttbar_sig3", "mu_sl_ttbar_sig3", 0.0, 10000. ) ;
      rv_mu_sl_ttbar_sig4 = new RooRealVar( "mu_sl_ttbar_sig4", "mu_sl_ttbar_sig4", 0.0, 10000. ) ;
      rv_mu_sl_ttbar_sig5 = new RooRealVar( "mu_sl_ttbar_sig5", "mu_sl_ttbar_sig5", 0.0, 10000. ) ;

      rv_mu_sl_ttbar_sig1 -> setVal( Nslsig1 - Newmcslsig1 ) ;
      rv_mu_sl_ttbar_sig2 -> setVal( Nslsig2 - Newmcslsig2 ) ;
      rv_mu_sl_ttbar_sig3 -> setVal( Nslsig3 - Newmcslsig3 ) ;
      rv_mu_sl_ttbar_sig4 -> setVal( Nslsig4 - Newmcslsig4 ) ;
      rv_mu_sl_ttbar_sig5 -> setVal( Nslsig5 - Newmcslsig5 ) ;


     //-- Single Lepton (SL) counts, ttbar, SB region, bins of 3-jet mass

      rv_mu_sl_ttbar_sb1 = new RooRealVar( "mu_sl_ttbar_sb1", "mu_sl_ttbar_sb1", 0.0, 10000. ) ;
      rv_mu_sl_ttbar_sb2 = new RooRealVar( "mu_sl_ttbar_sb2", "mu_sl_ttbar_sb2", 0.0, 10000. ) ;
      rv_mu_sl_ttbar_sb3 = new RooRealVar( "mu_sl_ttbar_sb3", "mu_sl_ttbar_sb3", 0.0, 10000. ) ;
      rv_mu_sl_ttbar_sb4 = new RooRealVar( "mu_sl_ttbar_sb4", "mu_sl_ttbar_sb4", 0.0, 10000. ) ;
      rv_mu_sl_ttbar_sb5 = new RooRealVar( "mu_sl_ttbar_sb5", "mu_sl_ttbar_sb5", 0.0, 10000. ) ;

      rv_mu_sl_ttbar_sb1 -> setVal( Nslsb1 - Newmcslsb1 ) ;
      rv_mu_sl_ttbar_sb2 -> setVal( Nslsb2 - Newmcslsb2 ) ;
      rv_mu_sl_ttbar_sb3 -> setVal( Nslsb3 - Newmcslsb3 ) ;
      rv_mu_sl_ttbar_sb4 -> setVal( Nslsb4 - Newmcslsb4 ) ;
      rv_mu_sl_ttbar_sb5 -> setVal( Nslsb5 - Newmcslsb5 ) ;


     //-- Single Lepton (SL) counts, ttbar, MSB region, bins of 3-jet mass

      rv_mu_sl_ttbar_msb1 = new RooRealVar( "mu_sl_ttbar_msb1", "mu_sl_ttbar_msb1", 0.0, 10000. ) ;
      rv_mu_sl_ttbar_msb2 = new RooRealVar( "mu_sl_ttbar_msb2", "mu_sl_ttbar_msb2", 0.0, 10000. ) ;
      rv_mu_sl_ttbar_msb3 = new RooRealVar( "mu_sl_ttbar_msb3", "mu_sl_ttbar_msb3", 0.0, 10000. ) ;
      rv_mu_sl_ttbar_msb4 = new RooRealVar( "mu_sl_ttbar_msb4", "mu_sl_ttbar_msb4", 0.0, 10000. ) ;
      rv_mu_sl_ttbar_msb5 = new RooRealVar( "mu_sl_ttbar_msb5", "mu_sl_ttbar_msb5", 0.0, 10000. ) ;

      rv_mu_sl_ttbar_msb1 -> setVal( Nslmsb1 - Newmcslmsb1 ) ;
      rv_mu_sl_ttbar_msb2 -> setVal( Nslmsb2 - Newmcslmsb2 ) ;
      rv_mu_sl_ttbar_msb3 -> setVal( Nslmsb3 - Newmcslmsb3 ) ;
      rv_mu_sl_ttbar_msb4 -> setVal( Nslmsb4 - Newmcslmsb4 ) ;
      rv_mu_sl_ttbar_msb5 -> setVal( Nslmsb5 - Newmcslmsb5 ) ;







     //-- Single Lepton (SL) counts, EW, SIG region, bins of 3-jet mass

      rv_mu_sl_ew_sig1 = new RooRealVar( "mu_sl_ew_sig1", "mu_sl_ew_sig1", 0.0, 10000. ) ;
      rv_mu_sl_ew_sig2 = new RooRealVar( "mu_sl_ew_sig2", "mu_sl_ew_sig2", 0.0, 10000. ) ;
      rv_mu_sl_ew_sig3 = new RooRealVar( "mu_sl_ew_sig3", "mu_sl_ew_sig3", 0.0, 10000. ) ;
      rv_mu_sl_ew_sig4 = new RooRealVar( "mu_sl_ew_sig4", "mu_sl_ew_sig4", 0.0, 10000. ) ;
      rv_mu_sl_ew_sig5 = new RooRealVar( "mu_sl_ew_sig5", "mu_sl_ew_sig5", 0.0, 10000. ) ;

      rv_mu_sl_ew_sig1 -> setVal( Newmcslsig1 ) ;
      rv_mu_sl_ew_sig2 -> setVal( Newmcslsig2 ) ;
      rv_mu_sl_ew_sig3 -> setVal( Newmcslsig3 ) ;
      rv_mu_sl_ew_sig4 -> setVal( Newmcslsig4 ) ;
      rv_mu_sl_ew_sig5 -> setVal( Newmcslsig5 ) ;

      rv_mu_sl_ew_sig1 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sig2 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sig3 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sig4 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sig5 -> setConstant(kTRUE) ;



     //-- Single Lepton (SL) counts, EW, SB region, bins of 3-jet mass

      rv_mu_sl_ew_sb1 = new RooRealVar( "mu_sl_ew_sb1", "mu_sl_ew_sb1", 0.0, 10000. ) ;
      rv_mu_sl_ew_sb2 = new RooRealVar( "mu_sl_ew_sb2", "mu_sl_ew_sb2", 0.0, 10000. ) ;
      rv_mu_sl_ew_sb3 = new RooRealVar( "mu_sl_ew_sb3", "mu_sl_ew_sb3", 0.0, 10000. ) ;
      rv_mu_sl_ew_sb4 = new RooRealVar( "mu_sl_ew_sb4", "mu_sl_ew_sb4", 0.0, 10000. ) ;
      rv_mu_sl_ew_sb5 = new RooRealVar( "mu_sl_ew_sb5", "mu_sl_ew_sb5", 0.0, 10000. ) ;

      rv_mu_sl_ew_sb1 -> setVal( Newmcslsb1 ) ;
      rv_mu_sl_ew_sb2 -> setVal( Newmcslsb2 ) ;
      rv_mu_sl_ew_sb3 -> setVal( Newmcslsb3 ) ;
      rv_mu_sl_ew_sb4 -> setVal( Newmcslsb4 ) ;
      rv_mu_sl_ew_sb5 -> setVal( Newmcslsb5 ) ;

      rv_mu_sl_ew_sb1 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sb2 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sb3 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sb4 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sb5 -> setConstant(kTRUE) ;



     //-- Single Lepton (SL) counts, EW, MSB region, bins of 3-jet mass

      rv_mu_sl_ew_msb1 = new RooRealVar( "mu_sl_ew_msb1", "mu_sl_ew_msb1", 0.0, 10000. ) ;
      rv_mu_sl_ew_msb2 = new RooRealVar( "mu_sl_ew_msb2", "mu_sl_ew_msb2", 0.0, 10000. ) ;
      rv_mu_sl_ew_msb3 = new RooRealVar( "mu_sl_ew_msb3", "mu_sl_ew_msb3", 0.0, 10000. ) ;
      rv_mu_sl_ew_msb4 = new RooRealVar( "mu_sl_ew_msb4", "mu_sl_ew_msb4", 0.0, 10000. ) ;
      rv_mu_sl_ew_msb5 = new RooRealVar( "mu_sl_ew_msb5", "mu_sl_ew_msb5", 0.0, 10000. ) ;

      rv_mu_sl_ew_msb1 -> setVal( Newmcslmsb1 ) ;
      rv_mu_sl_ew_msb2 -> setVal( Newmcslmsb2 ) ;
      rv_mu_sl_ew_msb3 -> setVal( Newmcslmsb3 ) ;
      rv_mu_sl_ew_msb4 -> setVal( Newmcslmsb4 ) ;
      rv_mu_sl_ew_msb5 -> setVal( Newmcslmsb5 ) ;

      rv_mu_sl_ew_msb1 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_msb2 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_msb3 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_msb4 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_msb5 -> setConstant(kTRUE) ;








     //-- Single Lepton (SL) counts, susy, SIG region, bins of 3-jet mass

      rv_mu_sl_susy_sig1 = new RooRealVar( "mu_sl_susy_sig1", "mu_sl_susy_sig1", 0.0, 10000. ) ;
      rv_mu_sl_susy_sig2 = new RooRealVar( "mu_sl_susy_sig2", "mu_sl_susy_sig2", 0.0, 10000. ) ;
      rv_mu_sl_susy_sig3 = new RooRealVar( "mu_sl_susy_sig3", "mu_sl_susy_sig3", 0.0, 10000. ) ;
      rv_mu_sl_susy_sig4 = new RooRealVar( "mu_sl_susy_sig4", "mu_sl_susy_sig4", 0.0, 10000. ) ;
      rv_mu_sl_susy_sig5 = new RooRealVar( "mu_sl_susy_sig5", "mu_sl_susy_sig5", 0.0, 10000. ) ;

      rv_mu_sl_susy_sig1 -> setVal( 0. ) ;
      rv_mu_sl_susy_sig2 -> setVal( 0. ) ;
      rv_mu_sl_susy_sig3 -> setVal( 0. ) ;
      rv_mu_sl_susy_sig4 -> setVal( 0. ) ;
      rv_mu_sl_susy_sig5 -> setVal( 0. ) ;

      rv_mu_sl_susy_sig1 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sig2 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sig3 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sig4 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sig5 -> setConstant( kTRUE ) ;


     //-- Single Lepton (SL) counts, susy, SB region, bins of 3-jet mass

      rv_mu_sl_susy_sb1 = new RooRealVar( "mu_sl_susy_sb1", "mu_sl_susy_sb1", 0.0, 10000. ) ;
      rv_mu_sl_susy_sb2 = new RooRealVar( "mu_sl_susy_sb2", "mu_sl_susy_sb2", 0.0, 10000. ) ;
      rv_mu_sl_susy_sb3 = new RooRealVar( "mu_sl_susy_sb3", "mu_sl_susy_sb3", 0.0, 10000. ) ;
      rv_mu_sl_susy_sb4 = new RooRealVar( "mu_sl_susy_sb4", "mu_sl_susy_sb4", 0.0, 10000. ) ;
      rv_mu_sl_susy_sb5 = new RooRealVar( "mu_sl_susy_sb5", "mu_sl_susy_sb5", 0.0, 10000. ) ;

      rv_mu_sl_susy_sb1 -> setVal( 0. ) ;
      rv_mu_sl_susy_sb2 -> setVal( 0. ) ;
      rv_mu_sl_susy_sb3 -> setVal( 0. ) ;
      rv_mu_sl_susy_sb4 -> setVal( 0. ) ;
      rv_mu_sl_susy_sb5 -> setVal( 0. ) ;

      rv_mu_sl_susy_sb1 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sb2 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sb3 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sb4 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sb5 -> setConstant( kTRUE ) ;


     //-- Single Lepton (SL) counts, susy, MSB region, bins of 3-jet mass

      rv_mu_sl_susy_msb1 = new RooRealVar( "mu_sl_susy_msb1", "mu_sl_susy_msb1", 0.0, 10000. ) ;
      rv_mu_sl_susy_msb2 = new RooRealVar( "mu_sl_susy_msb2", "mu_sl_susy_msb2", 0.0, 10000. ) ;
      rv_mu_sl_susy_msb3 = new RooRealVar( "mu_sl_susy_msb3", "mu_sl_susy_msb3", 0.0, 10000. ) ;
      rv_mu_sl_susy_msb4 = new RooRealVar( "mu_sl_susy_msb4", "mu_sl_susy_msb4", 0.0, 10000. ) ;
      rv_mu_sl_susy_msb5 = new RooRealVar( "mu_sl_susy_msb5", "mu_sl_susy_msb5", 0.0, 10000. ) ;

      rv_mu_sl_susy_msb1 -> setVal( 0. ) ;
      rv_mu_sl_susy_msb2 -> setVal( 0. ) ;
      rv_mu_sl_susy_msb3 -> setVal( 0. ) ;
      rv_mu_sl_susy_msb4 -> setVal( 0. ) ;
      rv_mu_sl_susy_msb5 -> setVal( 0. ) ;

      rv_mu_sl_susy_msb1 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_msb2 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_msb3 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_msb4 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_msb5 -> setConstant( kTRUE ) ;




     //+++++++++++++++++ Relationships between parameters

       printf(" --- Defining relationships between parameters.\n" ) ;

      rv_mu_sl_ttbar_sb = new RooFormulaVar( "mu_sl_ttbar_sb",
                                                             "mu_sl_ttbar_sb1 + mu_sl_ttbar_sb2 + mu_sl_ttbar_sb3 + mu_sl_ttbar_sb4 + mu_sl_ttbar_sb5",
                                                             RooArgSet( *rv_mu_sl_ttbar_sb1, *rv_mu_sl_ttbar_sb2, *rv_mu_sl_ttbar_sb3, *rv_mu_sl_ttbar_sb4, *rv_mu_sl_ttbar_sb5 ) ) ;

      rv_mu_sl_ttbar_sig = new RooFormulaVar( "mu_sl_ttbar_sig",
                                                             "mu_sl_ttbar_sig1 + mu_sl_ttbar_sig2 + mu_sl_ttbar_sig3 + mu_sl_ttbar_sig4 + mu_sl_ttbar_sig5",
                                                             RooArgSet( *rv_mu_sl_ttbar_sig1, *rv_mu_sl_ttbar_sig2, *rv_mu_sl_ttbar_sig3, *rv_mu_sl_ttbar_sig4, *rv_mu_sl_ttbar_sig5 ) ) ;




    //===================================================================================================
     //-- Use these lines to define things in terms of the SIG parameters.
      if ( useSigBgVars ) {

         rv_mu_ttbar_sb = new RooFormulaVar( "mu_ttbar_sb",
                                                          "mu_ttbar_sig*(mu_sl_ttbar_sb/mu_sl_ttbar_sig)",
                                                          RooArgSet( *rv_mu_ttbar_sig, *rv_mu_sl_ttbar_sb, *rv_mu_sl_ttbar_sig ) ) ;


         rv_mu_qcd_sb = new RooFormulaVar( "mu_qcd_sb",
                                                        "mu_qcdmc_sb * (mu_qcd_sig/mu_qcdmc_sig) * (mu_qcd_a/mu_qcdmc_a) * (mu_qcdmc_d/mu_qcd_d)",
                                                        RooArgSet( *rv_mu_qcdmc_sb,
                                                                   *rv_mu_qcd_sig, *rv_mu_qcdmc_sig,
                                                                   *rv_mu_qcd_a,  *rv_mu_qcdmc_a,
                                                                   *rv_mu_qcdmc_d, *rv_mu_qcd_d ) ) ;

         rfv_mu_ttbar_sb = ((RooFormulaVar*) rv_mu_ttbar_sb ) ;
         rfv_mu_qcd_sb   = ((RooFormulaVar*) rv_mu_qcd_sb ) ;

      } else {
    //===================================================================================================
     //-- Use these lines to define things in terms of the SB parameters.

         rv_mu_ttbar_sig = new RooFormulaVar( "mu_ttbar_sig",
                                                          "mu_ttbar_sb*(mu_sl_ttbar_sig/mu_sl_ttbar_sb)",
                                                          RooArgSet( *rv_mu_ttbar_sb, *rv_mu_sl_ttbar_sig, *rv_mu_sl_ttbar_sb ) ) ;


         rv_mu_qcd_sig = new RooFormulaVar( "mu_qcd_sig",
                                                        "mu_qcdmc_sig * (mu_qcd_sb/mu_qcdmc_sb) * (mu_qcdmc_a/mu_qcd_a) * (mu_qcd_d/mu_qcdmc_d)",
                                                        RooArgSet( *rv_mu_qcdmc_sig,
                                                                   *rv_mu_qcd_sb, *rv_mu_qcdmc_sb,
                                                                   *rv_mu_qcdmc_a,  *rv_mu_qcd_a,
                                                                   *rv_mu_qcd_d, *rv_mu_qcdmc_d ) ) ;

         rfv_mu_ttbar_sig = ((RooFormulaVar*) rv_mu_ttbar_sig ) ;
         rfv_mu_qcd_sig   = ((RooFormulaVar*) rv_mu_qcd_sig ) ;

      }
    //===================================================================================================






      rv_mu_qcd_lsb = new RooFormulaVar( "mu_qcd_lsb",
                                          "mu_qcd_lsb1 + mu_qcd_lsb2 + mu_qcd_lsb3 + mu_qcd_lsb4 + mu_qcd_lsb5",
                                          RooArgSet( *rv_mu_qcd_lsb1, *rv_mu_qcd_lsb2, *rv_mu_qcd_lsb3, *rv_mu_qcd_lsb4, *rv_mu_qcd_lsb5 ) ) ;

      rv_f_qcd_lsb1 = new RooFormulaVar( "f_qcd_lsb1",
                                          "mu_qcd_lsb1 / mu_qcd_lsb",
                                          RooArgSet( *rv_mu_qcd_lsb1, *rv_mu_qcd_lsb ) ) ;

      rv_f_qcd_lsb2 = new RooFormulaVar( "f_qcd_lsb2",
                                          "mu_qcd_lsb2 / mu_qcd_lsb",
                                          RooArgSet( *rv_mu_qcd_lsb2, *rv_mu_qcd_lsb ) ) ;

      rv_f_qcd_lsb3 = new RooFormulaVar( "f_qcd_lsb3",
                                          "mu_qcd_lsb3 / mu_qcd_lsb",
                                          RooArgSet( *rv_mu_qcd_lsb3, *rv_mu_qcd_lsb ) ) ;

      rv_f_qcd_lsb4 = new RooFormulaVar( "f_qcd_lsb4",
                                          "mu_qcd_lsb4 / mu_qcd_lsb",
                                          RooArgSet( *rv_mu_qcd_lsb4, *rv_mu_qcd_lsb ) ) ;

      rv_f_qcd_lsb5 = new RooFormulaVar( "f_qcd_lsb5",
                                          "mu_qcd_lsb5 / mu_qcd_lsb",
                                          RooArgSet( *rv_mu_qcd_lsb5, *rv_mu_qcd_lsb ) ) ;



      rv_mu_qcd_sb1 = new RooFormulaVar( "mu_qcd_sb1",
                                          "f_qcd_lsb1 * mu_qcd_sb",
                                          RooArgSet( *rv_f_qcd_lsb1, *rv_mu_qcd_sb ) ) ;

      rv_mu_qcd_sb2 = new RooFormulaVar( "mu_qcd_sb2",
                                          "f_qcd_lsb2 * mu_qcd_sb",
                                          RooArgSet( *rv_f_qcd_lsb2, *rv_mu_qcd_sb ) ) ;

      rv_mu_qcd_sb3 = new RooFormulaVar( "mu_qcd_sb3",
                                          "f_qcd_lsb3 * mu_qcd_sb",
                                          RooArgSet( *rv_f_qcd_lsb3, *rv_mu_qcd_sb ) ) ;

      rv_mu_qcd_sb4 = new RooFormulaVar( "mu_qcd_sb4",
                                          "f_qcd_lsb4 * mu_qcd_sb",
                                          RooArgSet( *rv_f_qcd_lsb4, *rv_mu_qcd_sb ) ) ;

      rv_mu_qcd_sb5 = new RooFormulaVar( "mu_qcd_sb5",
                                          "f_qcd_lsb5 * mu_qcd_sb",
                                          RooArgSet( *rv_f_qcd_lsb5, *rv_mu_qcd_sb ) ) ;





      rv_mu_sl_ttbar1 = new RooFormulaVar( "mu_sl_ttbar1",
                                            "mu_sl_ttbar_msb1 + mu_sl_ttbar_sb1 + mu_sl_ttbar_sig1",
                                            RooArgSet( *rv_mu_sl_ttbar_msb1, *rv_mu_sl_ttbar_sb1, *rv_mu_sl_ttbar_sig1 ) ) ;

      rv_mu_sl_ttbar2 = new RooFormulaVar( "mu_sl_ttbar2",
                                            "mu_sl_ttbar_msb2 + mu_sl_ttbar_sb2 + mu_sl_ttbar_sig2",
                                            RooArgSet( *rv_mu_sl_ttbar_msb2, *rv_mu_sl_ttbar_sb2, *rv_mu_sl_ttbar_sig2 ) ) ;

      rv_mu_sl_ttbar3 = new RooFormulaVar( "mu_sl_ttbar3",
                                            "mu_sl_ttbar_msb3 + mu_sl_ttbar_sb3 + mu_sl_ttbar_sig3",
                                            RooArgSet( *rv_mu_sl_ttbar_msb3, *rv_mu_sl_ttbar_sb3, *rv_mu_sl_ttbar_sig3 ) ) ;

      rv_mu_sl_ttbar4 = new RooFormulaVar( "mu_sl_ttbar4",
                                            "mu_sl_ttbar_msb4 + mu_sl_ttbar_sb4 + mu_sl_ttbar_sig4",
                                            RooArgSet( *rv_mu_sl_ttbar_msb4, *rv_mu_sl_ttbar_sb4, *rv_mu_sl_ttbar_sig4 ) ) ;

      rv_mu_sl_ttbar5 = new RooFormulaVar( "mu_sl_ttbar5",
                                            "mu_sl_ttbar_msb5 + mu_sl_ttbar_sb5 + mu_sl_ttbar_sig5",
                                            RooArgSet( *rv_mu_sl_ttbar_msb5, *rv_mu_sl_ttbar_sb5, *rv_mu_sl_ttbar_sig5 ) ) ;




      rv_mu_sl_ttbar = new RooFormulaVar( "mu_sl_ttbar",
                                           "mu_sl_ttbar1 + mu_sl_ttbar2 + mu_sl_ttbar3 + mu_sl_ttbar4 + mu_sl_ttbar5" ,
                                           RooArgSet( *rv_mu_sl_ttbar1, *rv_mu_sl_ttbar2, *rv_mu_sl_ttbar3, *rv_mu_sl_ttbar4, *rv_mu_sl_ttbar5  )  ) ;





      rv_f_sl_ttbar1 = new RooFormulaVar( "f_sl_ttbar1",
                                           "mu_sl_ttbar1 / mu_sl_ttbar",
                                           RooArgSet( *rv_mu_sl_ttbar1, *rv_mu_sl_ttbar ) ) ;

      rv_f_sl_ttbar2 = new RooFormulaVar( "f_sl_ttbar2",
                                           "mu_sl_ttbar2 / mu_sl_ttbar",
                                           RooArgSet( *rv_mu_sl_ttbar2, *rv_mu_sl_ttbar ) ) ;

      rv_f_sl_ttbar3 = new RooFormulaVar( "f_sl_ttbar3",
                                           "mu_sl_ttbar3 / mu_sl_ttbar",
                                           RooArgSet( *rv_mu_sl_ttbar3, *rv_mu_sl_ttbar ) ) ;

      rv_f_sl_ttbar4 = new RooFormulaVar( "f_sl_ttbar4",
                                           "mu_sl_ttbar4 / mu_sl_ttbar",
                                           RooArgSet( *rv_mu_sl_ttbar4, *rv_mu_sl_ttbar ) ) ;

      rv_f_sl_ttbar5 = new RooFormulaVar( "f_sl_ttbar5",
                                           "mu_sl_ttbar5 / mu_sl_ttbar",
                                           RooArgSet( *rv_mu_sl_ttbar5, *rv_mu_sl_ttbar ) ) ;




      rv_mu_ttbar_sb1 = new RooFormulaVar( "mu_ttbar_sb1",
                                            "f_sl_ttbar1 * mu_ttbar_sb",
                                            RooArgSet( *rv_f_sl_ttbar1, *rv_mu_ttbar_sb ) ) ;

      rv_mu_ttbar_sb2 = new RooFormulaVar( "mu_ttbar_sb2",
                                            "f_sl_ttbar2 * mu_ttbar_sb",
                                            RooArgSet( *rv_f_sl_ttbar2, *rv_mu_ttbar_sb ) ) ;

      rv_mu_ttbar_sb3 = new RooFormulaVar( "mu_ttbar_sb3",
                                            "f_sl_ttbar3 * mu_ttbar_sb",
                                            RooArgSet( *rv_f_sl_ttbar3, *rv_mu_ttbar_sb ) ) ;

      rv_mu_ttbar_sb4 = new RooFormulaVar( "mu_ttbar_sb4",
                                            "f_sl_ttbar4 * mu_ttbar_sb",
                                            RooArgSet( *rv_f_sl_ttbar4, *rv_mu_ttbar_sb ) ) ;

      rv_mu_ttbar_sb5 = new RooFormulaVar( "mu_ttbar_sb5",
                                            "f_sl_ttbar5 * mu_ttbar_sb",
                                            RooArgSet( *rv_f_sl_ttbar5, *rv_mu_ttbar_sb ) ) ;





    //+++++++++++++ Expected counts for observables in terms of parameters ++++++++++++++++++

       printf(" --- Defining expected counts in terms of parameters.\n" ) ;


      rv_n_sig = new RooFormulaVar( "n_sig",
                                     "mu_ttbar_sig + mu_qcd_sig + mu_ew_sig + mu_susy_sig",
                                     RooArgSet( *rv_mu_ttbar_sig, *rv_mu_qcd_sig,
                                     *rv_mu_ew_sig, *rv_mu_susy_sig ) ) ;





      rv_n_sb1 = new RooFormulaVar( "n_sb1",
                                     "mu_ttbar_sb1 + mu_qcd_sb1 + mu_ew_sb1 + mu_susy_sb1",
                                     RooArgSet( *rv_mu_ttbar_sb1, *rv_mu_qcd_sb1,
                                     *rv_mu_ew_sb1, *rv_mu_susy_sb1 ) ) ;

      rv_n_sb2 = new RooFormulaVar( "n_sb2",
                                     "mu_ttbar_sb2 + mu_qcd_sb2 + mu_ew_sb2 + mu_susy_sb2",
                                     RooArgSet( *rv_mu_ttbar_sb2, *rv_mu_qcd_sb2,
                                     *rv_mu_ew_sb2, *rv_mu_susy_sb2 ) ) ;

      rv_n_sb3 = new RooFormulaVar( "n_sb3",
                                     "mu_ttbar_sb3 + mu_qcd_sb3 + mu_ew_sb3 + mu_susy_sb3",
                                     RooArgSet( *rv_mu_ttbar_sb3, *rv_mu_qcd_sb3,
                                     *rv_mu_ew_sb3, *rv_mu_susy_sb3 ) ) ;

      rv_n_sb4 = new RooFormulaVar( "n_sb4",
                                     "mu_ttbar_sb4 + mu_qcd_sb4 + mu_ew_sb4 + mu_susy_sb4",
                                     RooArgSet( *rv_mu_ttbar_sb4, *rv_mu_qcd_sb4,
                                     *rv_mu_ew_sb4, *rv_mu_susy_sb4 ) ) ;

      rv_n_sb5 = new RooFormulaVar( "n_sb5",
                                     "mu_ttbar_sb5 + mu_qcd_sb5 + mu_ew_sb5 + mu_susy_sb5",
                                     RooArgSet( *rv_mu_ttbar_sb5, *rv_mu_qcd_sb5,
                                     *rv_mu_ew_sb5, *rv_mu_susy_sb5 ) ) ;





      rv_n_a = new RooFormulaVar( "n_a",
                                   "mu_qcd_a + mu_ttbar_a + mu_ew_a + mu_susy_a",
                                   RooArgSet( *rv_mu_qcd_a, *rv_mu_ttbar_a,
                                   *rv_mu_ew_a, *rv_mu_susy_a ) ) ;

      rv_n_d = new RooFormulaVar( "n_d",
                                   "mu_qcd_d + mu_ttbar_d + mu_ew_d + mu_susy_d",
                                   RooArgSet( *rv_mu_qcd_d, *rv_mu_ttbar_d,
                                   *rv_mu_ew_d, *rv_mu_susy_d ) ) ;




      rv_n_lsb1 = new RooFormulaVar( "n_lsb1",
                                      "mu_qcd_lsb1",
                                      RooArgSet( *rv_mu_qcd_lsb1 ) ) ;

      rv_n_lsb2 = new RooFormulaVar( "n_lsb2",
                                      "mu_qcd_lsb2",
                                      RooArgSet( *rv_mu_qcd_lsb2 ) ) ;

      rv_n_lsb3 = new RooFormulaVar( "n_lsb3",
                                      "mu_qcd_lsb3",
                                      RooArgSet( *rv_mu_qcd_lsb3 ) ) ;

      rv_n_lsb4 = new RooFormulaVar( "n_lsb4",
                                      "mu_qcd_lsb4",
                                      RooArgSet( *rv_mu_qcd_lsb4 ) ) ;

      rv_n_lsb5 = new RooFormulaVar( "n_lsb5",
                                      "mu_qcd_lsb5",
                                      RooArgSet( *rv_mu_qcd_lsb5 ) ) ;



      rv_n_sl_sig1 = new RooFormulaVar( "n_sl_sig1",
                                         "mu_sl_ttbar_sig1 + mu_sl_ew_sig1 + mu_sl_susy_sig1",
                                         RooArgSet( *rv_mu_sl_ttbar_sig1, *rv_mu_sl_ew_sig1, *rv_mu_sl_susy_sig1 ) ) ;

      rv_n_sl_sig2 = new RooFormulaVar( "n_sl_sig2",
                                         "mu_sl_ttbar_sig2 + mu_sl_ew_sig2 + mu_sl_susy_sig2",
                                         RooArgSet( *rv_mu_sl_ttbar_sig2, *rv_mu_sl_ew_sig2, *rv_mu_sl_susy_sig2 ) ) ;

      rv_n_sl_sig3 = new RooFormulaVar( "n_sl_sig3",
                                         "mu_sl_ttbar_sig3 + mu_sl_ew_sig3 + mu_sl_susy_sig3",
                                         RooArgSet( *rv_mu_sl_ttbar_sig3, *rv_mu_sl_ew_sig3, *rv_mu_sl_susy_sig3 ) ) ;

      rv_n_sl_sig4 = new RooFormulaVar( "n_sl_sig4",
                                         "mu_sl_ttbar_sig4 + mu_sl_ew_sig4 + mu_sl_susy_sig4",
                                         RooArgSet( *rv_mu_sl_ttbar_sig4, *rv_mu_sl_ew_sig4, *rv_mu_sl_susy_sig4 ) ) ;

      rv_n_sl_sig5 = new RooFormulaVar( "n_sl_sig5",
                                         "mu_sl_ttbar_sig5 + mu_sl_ew_sig5 + mu_sl_susy_sig5",
                                         RooArgSet( *rv_mu_sl_ttbar_sig5, *rv_mu_sl_ew_sig5, *rv_mu_sl_susy_sig5 ) ) ;




      rv_n_sl_sb1 = new RooFormulaVar( "n_sl_sb1",
                                         "mu_sl_ttbar_sb1 + mu_sl_ew_sb1 + mu_sl_susy_sb1",
                                         RooArgSet( *rv_mu_sl_ttbar_sb1, *rv_mu_sl_ew_sb1, *rv_mu_sl_susy_sb1 ) ) ;

      rv_n_sl_sb2 = new RooFormulaVar( "n_sl_sb2",
                                         "mu_sl_ttbar_sb2 + mu_sl_ew_sb2 + mu_sl_susy_sb2",
                                         RooArgSet( *rv_mu_sl_ttbar_sb2, *rv_mu_sl_ew_sb2, *rv_mu_sl_susy_sb2 ) ) ;

      rv_n_sl_sb3 = new RooFormulaVar( "n_sl_sb3",
                                         "mu_sl_ttbar_sb3 + mu_sl_ew_sb3 + mu_sl_susy_sb3",
                                         RooArgSet( *rv_mu_sl_ttbar_sb3, *rv_mu_sl_ew_sb3, *rv_mu_sl_susy_sb3 ) ) ;

      rv_n_sl_sb4 = new RooFormulaVar( "n_sl_sb4",
                                         "mu_sl_ttbar_sb4 + mu_sl_ew_sb4 + mu_sl_susy_sb4",
                                         RooArgSet( *rv_mu_sl_ttbar_sb4, *rv_mu_sl_ew_sb4, *rv_mu_sl_susy_sb4 ) ) ;

      rv_n_sl_sb5 = new RooFormulaVar( "n_sl_sb5",
                                         "mu_sl_ttbar_sb5 + mu_sl_ew_sb5 + mu_sl_susy_sb5",
                                         RooArgSet( *rv_mu_sl_ttbar_sb5, *rv_mu_sl_ew_sb5, *rv_mu_sl_susy_sb5 ) ) ;




      rv_n_sl_msb1 = new RooFormulaVar( "n_sl_msb1",
                                         "mu_sl_ttbar_msb1 + mu_sl_ew_msb1 + mu_sl_susy_msb1",
                                         RooArgSet( *rv_mu_sl_ttbar_msb1, *rv_mu_sl_ew_msb1, *rv_mu_sl_susy_msb1 ) ) ;

      rv_n_sl_msb2 = new RooFormulaVar( "n_sl_msb2",
                                         "mu_sl_ttbar_msb2 + mu_sl_ew_msb2 + mu_sl_susy_msb2",
                                         RooArgSet( *rv_mu_sl_ttbar_msb2, *rv_mu_sl_ew_msb2, *rv_mu_sl_susy_msb2 ) ) ;

      rv_n_sl_msb3 = new RooFormulaVar( "n_sl_msb3",
                                         "mu_sl_ttbar_msb3 + mu_sl_ew_msb3 + mu_sl_susy_msb3",
                                         RooArgSet( *rv_mu_sl_ttbar_msb3, *rv_mu_sl_ew_msb3, *rv_mu_sl_susy_msb3 ) ) ;

      rv_n_sl_msb4 = new RooFormulaVar( "n_sl_msb4",
                                         "mu_sl_ttbar_msb4 + mu_sl_ew_msb4 + mu_sl_susy_msb4",
                                         RooArgSet( *rv_mu_sl_ttbar_msb4, *rv_mu_sl_ew_msb4, *rv_mu_sl_susy_msb4 ) ) ;

      rv_n_sl_msb5 = new RooFormulaVar( "n_sl_msb5",
                                         "mu_sl_ttbar_msb5 + mu_sl_ew_msb5 + mu_sl_susy_msb5",
                                         RooArgSet( *rv_mu_sl_ttbar_msb5, *rv_mu_sl_ew_msb5, *rv_mu_sl_susy_msb5 ) ) ;













   //++++++++++++ PDFs for the likelihood +++++++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining PDFs of the likelihood.\n" ) ;

      pdf_Nsig = new RooPoisson( "pdf_Nsig", "Nsig Poisson PDF", *rv_Nsig, *rv_n_sig ) ;

      pdf_Na = new RooPoisson( "pdf_Na", "Na Poisson PDF", *rv_Na, *rv_n_a ) ;
      pdf_Nd = new RooPoisson( "pdf_Nd", "Nd Poisson PDF", *rv_Nd, *rv_n_d ) ;

      pdf_Nsb1 = new RooPoisson( "pdf_Nsb1", "Nsb1 Poisson PDF", *rv_Nsb1, *rv_n_sb1 ) ;
      pdf_Nsb2 = new RooPoisson( "pdf_Nsb2", "Nsb2 Poisson PDF", *rv_Nsb2, *rv_n_sb2 ) ;
      pdf_Nsb3 = new RooPoisson( "pdf_Nsb3", "Nsb3 Poisson PDF", *rv_Nsb3, *rv_n_sb3 ) ;
      pdf_Nsb4 = new RooPoisson( "pdf_Nsb4", "Nsb4 Poisson PDF", *rv_Nsb4, *rv_n_sb4 ) ;
      pdf_Nsb5 = new RooPoisson( "pdf_Nsb5", "Nsb5 Poisson PDF", *rv_Nsb5, *rv_n_sb5 ) ;

      pdf_Nlsb1 = new RooPoisson( "pdf_Nlsb1", "Nlsb1 Poisson PDF", *rv_Nlsb1, *rv_n_lsb1 ) ;
      pdf_Nlsb2 = new RooPoisson( "pdf_Nlsb2", "Nlsb2 Poisson PDF", *rv_Nlsb2, *rv_n_lsb2 ) ;
      pdf_Nlsb3 = new RooPoisson( "pdf_Nlsb3", "Nlsb3 Poisson PDF", *rv_Nlsb3, *rv_n_lsb3 ) ;
      pdf_Nlsb4 = new RooPoisson( "pdf_Nlsb4", "Nlsb4 Poisson PDF", *rv_Nlsb4, *rv_n_lsb4 ) ;
      pdf_Nlsb5 = new RooPoisson( "pdf_Nlsb5", "Nlsb5 Poisson PDF", *rv_Nlsb5, *rv_n_lsb5 ) ;

      pdf_Nsl_sig1 = new RooPoisson( "pdf_Nsl_sig1", "Nsl,sig1 Poisson PDF", *rv_Nslsig1, *rv_n_sl_sig1 ) ;
      pdf_Nsl_sig2 = new RooPoisson( "pdf_Nsl_sig2", "Nsl,sig2 Poisson PDF", *rv_Nslsig2, *rv_n_sl_sig2 ) ;
      pdf_Nsl_sig3 = new RooPoisson( "pdf_Nsl_sig3", "Nsl,sig3 Poisson PDF", *rv_Nslsig3, *rv_n_sl_sig3 ) ;
      pdf_Nsl_sig4 = new RooPoisson( "pdf_Nsl_sig4", "Nsl,sig4 Poisson PDF", *rv_Nslsig4, *rv_n_sl_sig4 ) ;
      pdf_Nsl_sig5 = new RooPoisson( "pdf_Nsl_sig5", "Nsl,sig5 Poisson PDF", *rv_Nslsig5, *rv_n_sl_sig5 ) ;

      pdf_Nsl_sb1 = new RooPoisson( "pdf_Nsl_sb1", "Nsl,sb1 Poisson PDF", *rv_Nslsb1, *rv_n_sl_sb1 ) ;
      pdf_Nsl_sb2 = new RooPoisson( "pdf_Nsl_sb2", "Nsl,sb2 Poisson PDF", *rv_Nslsb2, *rv_n_sl_sb2 ) ;
      pdf_Nsl_sb3 = new RooPoisson( "pdf_Nsl_sb3", "Nsl,sb3 Poisson PDF", *rv_Nslsb3, *rv_n_sl_sb3 ) ;
      pdf_Nsl_sb4 = new RooPoisson( "pdf_Nsl_sb4", "Nsl,sb4 Poisson PDF", *rv_Nslsb4, *rv_n_sl_sb4 ) ;
      pdf_Nsl_sb5 = new RooPoisson( "pdf_Nsl_sb5", "Nsl,sb5 Poisson PDF", *rv_Nslsb5, *rv_n_sl_sb5 ) ;

      pdf_Nsl_msb1 = new RooPoisson( "pdf_Nsl_msb1", "Nsl,msb1 Poisson PDF", *rv_Nslmsb1, *rv_n_sl_msb1 ) ;
      pdf_Nsl_msb2 = new RooPoisson( "pdf_Nsl_msb2", "Nsl,msb2 Poisson PDF", *rv_Nslmsb2, *rv_n_sl_msb2 ) ;
      pdf_Nsl_msb3 = new RooPoisson( "pdf_Nsl_msb3", "Nsl,msb3 Poisson PDF", *rv_Nslmsb3, *rv_n_sl_msb3 ) ;
      pdf_Nsl_msb4 = new RooPoisson( "pdf_Nsl_msb4", "Nsl,msb4 Poisson PDF", *rv_Nslmsb4, *rv_n_sl_msb4 ) ;
      pdf_Nsl_msb5 = new RooPoisson( "pdf_Nsl_msb5", "Nsl,msb5 Poisson PDF", *rv_Nslmsb5, *rv_n_sl_msb5 ) ;

////  pdf_Nqcdmc_sig  = new RooGaussian( "pdf_Nqcdmc_sig", "Gaussian pdf for Nqcdmc,sig",
////                                      *rv_mu_qcdmc_sig, RooConst( Nqcdmcsig ), RooConst( Nqcdmcsigerr ) ) ;
////  pdf_Nqcdmc_sb  = new RooGaussian( "pdf_Nqcdmc_sb", "Gaussian pdf for Nqcdmc,sb",
////                                      *rv_mu_qcdmc_sb, RooConst( Nqcdmcsb ), RooConst( Nqcdmcsberr ) ) ;
////  pdf_Nqcdmc_a   = new RooGaussian( "pdf_Nqcdmc_a", "Gaussian pdf for Nqcdmc,a",
////                                      *rv_mu_qcdmc_a, RooConst( Nqcdmca ), RooConst( Nqcdmcaerr ) ) ;
////  pdf_Nqcdmc_d   = new RooGaussian( "pdf_Nqcdmc_d", "Gaussian pdf for Nqcdmc,d",
////                                      *rv_mu_qcdmc_d, RooConst( Nqcdmcd ), RooConst( Nqcdmcderr ) ) ;

      pdf_Nqcdmc_sig  = new RooGaussian( "pdf_Nqcdmc_sig", "Gaussian pdf for Nqcdmc,sig",
                                          *rv_mu_qcdmc_sig, *rv_Nqcdmcsig, RooConst( Nqcdmcsigerr ) ) ;
      pdf_Nqcdmc_sb  = new RooGaussian( "pdf_Nqcdmc_sb", "Gaussian pdf for Nqcdmc,sb",
                                          *rv_mu_qcdmc_sb, *rv_Nqcdmcsb, RooConst( Nqcdmcsberr ) ) ;
      pdf_Nqcdmc_a   = new RooGaussian( "pdf_Nqcdmc_a", "Gaussian pdf for Nqcdmc,a",
                                          *rv_mu_qcdmc_a, *rv_Nqcdmca, RooConst( Nqcdmcaerr ) ) ;
      pdf_Nqcdmc_d   = new RooGaussian( "pdf_Nqcdmc_d", "Gaussian pdf for Nqcdmc,d",
                                          *rv_mu_qcdmc_d, *rv_Nqcdmcd, RooConst( Nqcdmcderr ) ) ;

      {
         RooArgSet pdflist ;
         pdflist.add( *pdf_Nsig ) ;
         pdflist.add( *pdf_Na ) ;
         pdflist.add( *pdf_Nd ) ;
         pdflist.add( *pdf_Nsb1 ) ;
         pdflist.add( *pdf_Nsb2 ) ;
         pdflist.add( *pdf_Nsb3 ) ;
         pdflist.add( *pdf_Nsb4 ) ;
         pdflist.add( *pdf_Nsb5 ) ;
         pdflist.add( *pdf_Nlsb1 ) ;
         pdflist.add( *pdf_Nlsb2 ) ;
         pdflist.add( *pdf_Nlsb3 ) ;
         pdflist.add( *pdf_Nlsb4 ) ;
         pdflist.add( *pdf_Nlsb5 ) ;
         pdflist.add( *pdf_Nsl_sig1 ) ;
         pdflist.add( *pdf_Nsl_sig2 ) ;
         pdflist.add( *pdf_Nsl_sig3 ) ;
         pdflist.add( *pdf_Nsl_sig4 ) ;
         pdflist.add( *pdf_Nsl_sig5 ) ;
         pdflist.add( *pdf_Nsl_sb1 ) ;
         pdflist.add( *pdf_Nsl_sb2 ) ;
         pdflist.add( *pdf_Nsl_sb3 ) ;
         pdflist.add( *pdf_Nsl_sb4 ) ;
         pdflist.add( *pdf_Nsl_sb5 ) ;
         pdflist.add( *pdf_Nsl_msb1 ) ;
         pdflist.add( *pdf_Nsl_msb2 ) ;
         pdflist.add( *pdf_Nsl_msb3 ) ;
         pdflist.add( *pdf_Nsl_msb4 ) ;
         pdflist.add( *pdf_Nsl_msb5 ) ;
         pdflist.add( *pdf_Nqcdmc_sig ) ;
         pdflist.add( *pdf_Nqcdmc_sb ) ;
         pdflist.add( *pdf_Nqcdmc_a ) ;
         pdflist.add( *pdf_Nqcdmc_d ) ;
         likelihood = new RooProdPdf( "likelihood", "ra2b likelihood", pdflist ) ;
      }


     //---- Define the list of observables.

       observedParametersList.add( *rv_Nsig ) ;
       observedParametersList.add( *rv_Na ) ;
       observedParametersList.add( *rv_Nd ) ;
       observedParametersList.add( *rv_Nsb1 ) ;
       observedParametersList.add( *rv_Nsb2 ) ;
       observedParametersList.add( *rv_Nsb3 ) ;
       observedParametersList.add( *rv_Nsb4 ) ;
       observedParametersList.add( *rv_Nsb5 ) ;
       observedParametersList.add( *rv_Nlsb1 ) ;
       observedParametersList.add( *rv_Nlsb2 ) ;
       observedParametersList.add( *rv_Nlsb3 ) ;
       observedParametersList.add( *rv_Nlsb4 ) ;
       observedParametersList.add( *rv_Nlsb5 ) ;
       observedParametersList.add( *rv_Nslsig1 ) ;
       observedParametersList.add( *rv_Nslsig2 ) ;
       observedParametersList.add( *rv_Nslsig3 ) ;
       observedParametersList.add( *rv_Nslsig4 ) ;
       observedParametersList.add( *rv_Nslsig5 ) ;
       observedParametersList.add( *rv_Nslsb1 ) ;
       observedParametersList.add( *rv_Nslsb2 ) ;
       observedParametersList.add( *rv_Nslsb3 ) ;
       observedParametersList.add( *rv_Nslsb4 ) ;
       observedParametersList.add( *rv_Nslsb5 ) ;
       observedParametersList.add( *rv_Nslmsb1 ) ;
       observedParametersList.add( *rv_Nslmsb2 ) ;
       observedParametersList.add( *rv_Nslmsb3 ) ;
       observedParametersList.add( *rv_Nslmsb4 ) ;
       observedParametersList.add( *rv_Nslmsb5 ) ;

       observedParametersList.add( *rv_Nqcdmca ) ;
       observedParametersList.add( *rv_Nqcdmcd ) ;
       observedParametersList.add( *rv_Nqcdmcsig ) ;
       observedParametersList.add( *rv_Nqcdmcsb ) ;


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

    bool ra2bRoostatsClass::reinitialize( ) {

       printf( "\n\n Opening input file : %s\n\n", initializeFile ) ;

       FILE* infp ;
       if ( (infp=fopen( initializeFile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", initializeFile ) ;
          return false ;
       }

       int N3jmBins(0) ; //-- number of 3-jet mass bins.

       int Nsig(0), Na(0), Nd(0) ; //-- data counts in signal region, A, and D.
       int Nsb1(0), Nsb2(0), Nsb3(0), Nsb4(0), Nsb5(0) ; //-- data counts in 5 3-jet mass bins of SB.
       int Nlsb1(0), Nlsb2(0), Nlsb3(0), Nlsb4(0), Nlsb5(0) ; //-- data counts in 5 3-jet mass bins of LSB.
       int Nslsig1(0), Nslsig2(0), Nslsig3(0), Nslsig4(0), Nslsig5(0) ; //-- data counts in 5 3-jet mass bins of SL, SIG.
       int Nslsb1(0), Nslsb2(0), Nslsb3(0), Nslsb4(0), Nslsb5(0) ; ; //-- data counts in 5 3-jet mass bins of SL, SB.
       int Nslmsb1(0), Nslmsb2(0), Nslmsb3(0), Nslmsb4(0), Nslmsb5(0) ; //-- data counts in 5 3-jet mass bins of SL, MSB.

       float Nqcdmcsig(0.), Nqcdmcsb(0.), Nqcdmca(0.), Nqcdmcd(0.) ; //-- QCD MC counts in SIG, SB, A, and D.
       float Nqcdmcsigerr(0.), Nqcdmcsberr(0.), Nqcdmcaerr(0.), Nqcdmcderr(0.) ; //-- QCD MC uncertainties in SIG, SB, A, and D.
       float Nqcdmcslsig(0.), Nqcdmcslsb(0.), Nqcdmcslmsb(0.) ; //-- QCD MC counts in SL; SIG, LSB, and MSB.
       float Nttbarmcsig(0.), Nttbarmcsb(0.), Nttbarmca(0.), Nttbarmcd(0.) ; //-- ttbar MC counts in SIG, SB, A, and D.
       float Nttbarmcslsig(0.), Nttbarmcslsb(0.), Nttbarmcslmsb(0.) ; //-- ttbar MC counts in SL; SIG, SB, and MSB.
       float Newmcsig(0.), Newmca(0.), Newmcd(0.) ; //-- EW MC counts in SIG, A, and D.
       float Newmcsb1(0.), Newmcsb2(0.), Newmcsb3(0.), Newmcsb4(0.), Newmcsb5(0.) ; //-- EW MC counts in 5 3-jet mass bins of SB.
       float Newmcslsig1(0.), Newmcslsig2(0.), Newmcslsig3(0.), Newmcslsig4(0.), Newmcslsig5(0.) ; //-- EW MC, SL counts in SIG, bins of 3-jet mass
       float Newmcslsb1(0.), Newmcslsb2(0.), Newmcslsb3(0.), Newmcslsb4(0.), Newmcslsb5(0.) ; //-- EW MC, SL counts in SB, bins of 3-jet mass
       float Newmcslmsb1(0.), Newmcslmsb2(0.), Newmcslmsb3(0.), Newmcslmsb4(0.), Newmcslmsb5(0.) ; //-- EW MC, SL counts in MSB, bins of 3-jet mass

       float Nsusymcsig(0.), Nsusymca(0.), Nsusymcd(0.) ; //-- SUSY MC counts in SIG, A, and D.
       float Nsusymcsb1(0.), Nsusymcsb2(0.), Nsusymcsb3(0.), Nsusymcsb4(0.), Nsusymcsb5(0.) ; //-- SUSY MC counts in 5 3-jet mass bins of SB.
       float Nsusymcslsig1(0.), Nsusymcslsig2(0.), Nsusymcslsig3(0.), Nsusymcslsig4(0.), Nsusymcslsig5(0.) ; //-- SUSY MC counts in 5 3-jet mass bins of SL, SIG.
       float Nsusymcslsb1(0.), Nsusymcslsb2(0.), Nsusymcslsb3(0.), Nsusymcslsb4(0.), Nsusymcslsb5(0.) ; //-- SUSY MC counts in 5 3-jet mass bins of SL, SB.
       float Nsusymcslmsb1(0.), Nsusymcslmsb2(0.), Nsusymcslmsb3(0.), Nsusymcslmsb4(0.), Nsusymcslmsb5(0.) ; //-- SUSY MC counts in 5 3-jet mass bins of SL, MSB.

       //--- read in description line.
       printf("\n\n") ;
       char c(0) ;
       while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
       printf("\n\n") ;


       char label[1000] ;

      //--- Inputs generated with gen_roostats_input.c
      //    The order here must be consistent with the order there!

       fscanf( infp, "%s %d", label, &N3jmBins ) ;             printf( "%s %d\n", label, N3jmBins ) ;                
       fscanf( infp, "%s %d", label, &Nsig ) ;                 printf( "%s %d\n", label, Nsig ) ;         
       fscanf( infp, "%s %d", label, &Na ) ;                   printf( "%s %d\n", label, Na ) ;           
       fscanf( infp, "%s %d", label, &Nd ) ;                   printf( "%s %d\n", label, Nd ) ;           
       fscanf( infp, "%s %d", label, &Nsb1 ) ;                 printf( "%s %d\n", label, Nsb1 ) ;         
       fscanf( infp, "%s %d", label, &Nsb2 ) ;                 printf( "%s %d\n", label, Nsb2 ) ;         
       fscanf( infp, "%s %d", label, &Nsb3 ) ;                 printf( "%s %d\n", label, Nsb3 ) ;         
       fscanf( infp, "%s %d", label, &Nsb4 ) ;                 printf( "%s %d\n", label, Nsb4 ) ;         
       fscanf( infp, "%s %d", label, &Nsb5 ) ;                 printf( "%s %d\n", label, Nsb5 ) ;         
       fscanf( infp, "%s %d", label, &Nlsb1 ) ;                printf( "%s %d\n", label, Nlsb1 ) ;        
       fscanf( infp, "%s %d", label, &Nlsb2 ) ;                printf( "%s %d\n", label, Nlsb2 ) ;        
       fscanf( infp, "%s %d", label, &Nlsb3 ) ;                printf( "%s %d\n", label, Nlsb3 ) ;        
       fscanf( infp, "%s %d", label, &Nlsb4 ) ;                printf( "%s %d\n", label, Nlsb4 ) ;        
       fscanf( infp, "%s %d", label, &Nlsb5 ) ;                printf( "%s %d\n", label, Nlsb5 ) ;        
       fscanf( infp, "%s %d", label, &Nslsig1 ) ;              printf( "%s %d\n", label, Nslsig1 ) ;      
       fscanf( infp, "%s %d", label, &Nslsig2 ) ;              printf( "%s %d\n", label, Nslsig2 ) ;      
       fscanf( infp, "%s %d", label, &Nslsig3 ) ;              printf( "%s %d\n", label, Nslsig3 ) ;      
       fscanf( infp, "%s %d", label, &Nslsig4 ) ;              printf( "%s %d\n", label, Nslsig4 ) ;      
       fscanf( infp, "%s %d", label, &Nslsig5 ) ;              printf( "%s %d\n", label, Nslsig5 ) ;      
       fscanf( infp, "%s %d", label, &Nslsb1 ) ;               printf( "%s %d\n", label, Nslsb1 ) ;       
       fscanf( infp, "%s %d", label, &Nslsb2 ) ;               printf( "%s %d\n", label, Nslsb2 ) ;       
       fscanf( infp, "%s %d", label, &Nslsb3 ) ;               printf( "%s %d\n", label, Nslsb3 ) ;       
       fscanf( infp, "%s %d", label, &Nslsb4 ) ;               printf( "%s %d\n", label, Nslsb4 ) ;       
       fscanf( infp, "%s %d", label, &Nslsb5 ) ;               printf( "%s %d\n", label, Nslsb5 ) ;       
       fscanf( infp, "%s %d", label, &Nslmsb1 ) ;              printf( "%s %d\n", label, Nslmsb1 ) ;      
       fscanf( infp, "%s %d", label, &Nslmsb2 ) ;              printf( "%s %d\n", label, Nslmsb2 ) ;      
       fscanf( infp, "%s %d", label, &Nslmsb3 ) ;              printf( "%s %d\n", label, Nslmsb3 ) ;      
       fscanf( infp, "%s %d", label, &Nslmsb4 ) ;              printf( "%s %d\n", label, Nslmsb4 ) ;      
       fscanf( infp, "%s %d", label, &Nslmsb5 ) ;              printf( "%s %d\n", label, Nslmsb5 ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcsig ) ;            printf( "%s %g\n", label, Nqcdmcsig ) ;    
       fscanf( infp, "%s %g", label, &Nqcdmcsigerr ) ;         printf( "%s %g\n", label, Nqcdmcsigerr ) ;    
       fscanf( infp, "%s %g", label, &Nqcdmcsb ) ;             printf( "%s %g\n", label, Nqcdmcsb ) ;     
       fscanf( infp, "%s %g", label, &Nqcdmcsberr ) ;          printf( "%s %g\n", label, Nqcdmcsberr ) ;     
       fscanf( infp, "%s %g", label, &Nqcdmca ) ;              printf( "%s %g\n", label, Nqcdmca ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcaerr ) ;           printf( "%s %g\n", label, Nqcdmcaerr ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcd ) ;              printf( "%s %g\n", label, Nqcdmcd ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcderr ) ;           printf( "%s %g\n", label, Nqcdmcderr ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcslsig ) ;          printf( "%s %g\n", label, Nqcdmcslsig ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcslsb ) ;           printf( "%s %g\n", label, Nqcdmcslsb ) ;      
       fscanf( infp, "%s %g", label, &Nqcdmcslmsb ) ;          printf( "%s %g\n", label, Nqcdmcslmsb ) ;      
       fscanf( infp, "%s %g", label, &Nttbarmcsig ) ;          printf( "%s %g\n", label, Nttbarmcsig ) ;  
       fscanf( infp, "%s %g", label, &Nttbarmcsb ) ;           printf( "%s %g\n", label, Nttbarmcsb ) ;   
       fscanf( infp, "%s %g", label, &Nttbarmca ) ;            printf( "%s %g\n", label, Nttbarmca ) ;    
       fscanf( infp, "%s %g", label, &Nttbarmcd ) ;            printf( "%s %g\n", label, Nttbarmcd ) ;    
       fscanf( infp, "%s %g", label, &Nttbarmcslsig ) ;        printf( "%s %g\n", label, Nttbarmcslsig ) ;      
       fscanf( infp, "%s %g", label, &Nttbarmcslsb ) ;         printf( "%s %g\n", label, Nttbarmcslsb ) ;      
       fscanf( infp, "%s %g", label, &Nttbarmcslmsb ) ;        printf( "%s %g\n", label, Nttbarmcslmsb ) ;      
       fscanf( infp, "%s %g", label, &Newmcsig ) ;             printf( "%s %g\n", label, Newmcsig ) ;     
       fscanf( infp, "%s %g", label, &Newmca ) ;               printf( "%s %g\n", label, Newmca ) ;       
       fscanf( infp, "%s %g", label, &Newmcd ) ;               printf( "%s %g\n", label, Newmcd ) ;       
       fscanf( infp, "%s %g", label, &Newmcsb1 ) ;             printf( "%s %g\n", label, Newmcsb1 ) ;     
       fscanf( infp, "%s %g", label, &Newmcsb2 ) ;             printf( "%s %g\n", label, Newmcsb2 ) ;     
       fscanf( infp, "%s %g", label, &Newmcsb3 ) ;             printf( "%s %g\n", label, Newmcsb3 ) ;     
       fscanf( infp, "%s %g", label, &Newmcsb4 ) ;             printf( "%s %g\n", label, Newmcsb4 ) ;     
       fscanf( infp, "%s %g", label, &Newmcsb5 ) ;             printf( "%s %g\n", label, Newmcsb5 ) ;     
       fscanf( infp, "%s %g", label, &Newmcslsig1 ) ;          printf( "%s %g\n", label, Newmcslsig1 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsig2 ) ;          printf( "%s %g\n", label, Newmcslsig2 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsig3 ) ;          printf( "%s %g\n", label, Newmcslsig3 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsig4 ) ;          printf( "%s %g\n", label, Newmcslsig4 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsig5 ) ;          printf( "%s %g\n", label, Newmcslsig5 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsb1 ) ;          printf( "%s %g\n", label, Newmcslsb1 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsb2 ) ;          printf( "%s %g\n", label, Newmcslsb2 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsb3 ) ;          printf( "%s %g\n", label, Newmcslsb3 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsb4 ) ;          printf( "%s %g\n", label, Newmcslsb4 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslsb5 ) ;          printf( "%s %g\n", label, Newmcslsb5 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslmsb1 ) ;          printf( "%s %g\n", label, Newmcslmsb1 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslmsb2 ) ;          printf( "%s %g\n", label, Newmcslmsb2 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslmsb3 ) ;          printf( "%s %g\n", label, Newmcslmsb3 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslmsb4 ) ;          printf( "%s %g\n", label, Newmcslmsb4 ) ;      
       fscanf( infp, "%s %g", label, &Newmcslmsb5 ) ;          printf( "%s %g\n", label, Newmcslmsb5 ) ;      
       fscanf( infp, "%s %g", label, &Nsusymcsig ) ;           printf( "%s %g\n", label, Nsusymcsig ) ;   
       fscanf( infp, "%s %g", label, &Nsusymca ) ;             printf( "%s %g\n", label, Nsusymca ) ;     
       fscanf( infp, "%s %g", label, &Nsusymcd ) ;             printf( "%s %g\n", label, Nsusymcd ) ;     
       fscanf( infp, "%s %g", label, &Nsusymcsb1 ) ;           printf( "%s %g\n", label, Nsusymcsb1 ) ;   
       fscanf( infp, "%s %g", label, &Nsusymcsb2 ) ;           printf( "%s %g\n", label, Nsusymcsb2 ) ;   
       fscanf( infp, "%s %g", label, &Nsusymcsb3 ) ;           printf( "%s %g\n", label, Nsusymcsb3 ) ;   
       fscanf( infp, "%s %g", label, &Nsusymcsb4 ) ;           printf( "%s %g\n", label, Nsusymcsb4 ) ;   
       fscanf( infp, "%s %g", label, &Nsusymcsb5 ) ;           printf( "%s %g\n", label, Nsusymcsb5 ) ;   
       fscanf( infp, "%s %g", label, &Nsusymcslsig1 ) ;        printf( "%s %g\n", label, Nsusymcslsig1 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslsig2 ) ;        printf( "%s %g\n", label, Nsusymcslsig2 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslsig3 ) ;        printf( "%s %g\n", label, Nsusymcslsig3 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslsig4 ) ;        printf( "%s %g\n", label, Nsusymcslsig4 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslsig5 ) ;        printf( "%s %g\n", label, Nsusymcslsig5 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslsb1 ) ;         printf( "%s %g\n", label, Nsusymcslsb1 ) ; 
       fscanf( infp, "%s %g", label, &Nsusymcslsb2 ) ;         printf( "%s %g\n", label, Nsusymcslsb2 ) ; 
       fscanf( infp, "%s %g", label, &Nsusymcslsb3 ) ;         printf( "%s %g\n", label, Nsusymcslsb3 ) ; 
       fscanf( infp, "%s %g", label, &Nsusymcslsb4 ) ;         printf( "%s %g\n", label, Nsusymcslsb4 ) ; 
       fscanf( infp, "%s %g", label, &Nsusymcslsb5 ) ;         printf( "%s %g\n", label, Nsusymcslsb5 ) ; 
       fscanf( infp, "%s %g", label, &Nsusymcslmsb1 ) ;        printf( "%s %g\n", label, Nsusymcslmsb1 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslmsb2 ) ;        printf( "%s %g\n", label, Nsusymcslmsb2 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslmsb3 ) ;        printf( "%s %g\n", label, Nsusymcslmsb3 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslmsb4 ) ;        printf( "%s %g\n", label, Nsusymcslmsb4 ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslmsb5 ) ;        printf( "%s %g\n", label, Nsusymcslmsb5 ) ;


       printf("\n Done reading in %s\n\n", initializeFile ) ;
       fclose( infp ) ;






       //--- Print out a nice summary of the inputs.


       float Nsmsig = Nttbarmcsig + Nqcdmcsig + Newmcsig ;

       int   Nsb = Nsb1+Nsb2+Nsb3+Nsb4+Nsb5 ;
       float Newmcsb =Newmcsb1+Newmcsb2+Newmcsb3+Newmcsb4+Newmcsb5;
       float Nsusymcsb =Nsusymcsb1+Nsusymcsb2+Nsusymcsb3+Nsusymcsb4+Nsusymcsb5;
       float Nsmsb = Nttbarmcsb + Nqcdmcsb + Newmcsb ;

       float Nsma = Nttbarmca + Nqcdmca + Newmca ;
       float Nsmd = Nttbarmcd + Nqcdmcd + Newmcd ;

       float Newmcslsig = Newmcslsig1 + Newmcslsig2 + Newmcslsig3 + Newmcslsig4 + Newmcslsig5 ;
       float Nsmslsig = Nttbarmcslsig + Nqcdmcslsig + Newmcslsig ;
       int   Nslsig = Nslsig1+Nslsig2+Nslsig3+Nslsig4+Nslsig5 ;
       float Nsusymcslsig =Nsusymcslsig1+Nsusymcslsig2+Nsusymcslsig3+Nsusymcslsig4+Nsusymcslsig5;

       float Newmcslsb = Newmcslsb1 + Newmcslsb2 + Newmcslsb3 + Newmcslsb4 + Newmcslsb5 ;
       float Nsmslsb = Nttbarmcslsb + Nqcdmcslsb + Newmcslsb ;
       int   Nslsb = Nslsb1+Nslsb2+Nslsb3+Nslsb4+Nslsb5 ;
       float Nsusymcslsb =Nsusymcslsb1+Nsusymcslsb2+Nsusymcslsb3+Nsusymcslsb4+Nsusymcslsb5;

       int   Nslmsb = Nslmsb1+Nslmsb2+Nslmsb3+Nslmsb4+Nslmsb5 ;
       float Nsusymcslmsb =Nsusymcslmsb1+Nsusymcslmsb2+Nsusymcslmsb3+Nsusymcslmsb4+Nsusymcslmsb5;

       float Nttbarmcslmsbsbsig = Nttbarmcslmsb+Nttbarmcslsb+Nttbarmcslsig ;
       float Nqcdmcslmsbsbsig = Nqcdmcslmsb+Nqcdmcslsb+Nqcdmcslsig ;
       float Newmcslmsb = Newmcslmsb1 + Newmcslmsb2 + Newmcslmsb3 + Newmcslmsb4 + Newmcslmsb5 ;
       float Newmcslmsbsbsig = Newmcslmsb+Newmcslsb+Newmcslsig ;
       float Nsmmcslmsbsbsig = Nttbarmcslmsbsbsig + Nqcdmcslmsbsbsig + Newmcslmsbsbsig ;
       int   Nslmsbsbsig = Nslmsb+Nslsb+Nslsig ;
       float Nsusymcslmsbsbsig = Nsusymcslmsb+Nsusymcslsb+Nsusymcslsig ;

       float Rttmc = Nttbarmcsig / Nttbarmcsb ;
       float Rslmc = Nsmslsig / Nsmslsb ;
       float Rsldata = (1.0*Nslsig) / (1.0*Nslsb) ;
       float Rsldataerr = Rsldata*sqrt(1.0/Nslsig + 1.0/Nslsb) ;

       float Rqcdmctruth = Nqcdmcd / Nqcdmca ;
       float Rqcddata    = (Nd-Nttbarmcd-Newmcd)/(Na-Nttbarmca-Newmca) ;
       float Rqcddataerr = sqrt( pow(1./(Na-Nttbarmca-Newmca),2)*Nd + pow( ((Nd-Nttbarmcd-Newmcd)/pow((Na-Nttbarmca-Newmca),2)),2)*Na ) ;

       float frsigttbar = Nttbarmcsig / Nsmsig ;
       float frsigqcd   = Nqcdmcsig / Nsmsig ;
       float frsigew    = Newmcsig / Nsmsig ;

       float frsbttbar = Nttbarmcsb / Nsmsb ;
       float frsbqcd   = Nqcdmcsb / Nsmsb ;
       float frsbew    = Newmcsb / Nsmsb ;

       float frslttbar = Nttbarmcslmsbsbsig / Nsmmcslmsbsbsig ;
       float frslqcd   = Nqcdmcslmsbsbsig / Nsmmcslmsbsbsig ;
       float frslew    = Newmcslmsbsbsig / Nsmmcslmsbsbsig ;



       printf("\n\n\n") ;

       printf("                        ttbar    qcd       EW    all SM      data   SUSY\n") ;
       printf("----------------------------------------------------------------------------\n") ;
       printf("  Signal box     :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmcsig, Nqcdmcsig, Newmcsig, Nsmsig, Nsig, Nsusymcsig ) ;
       printf("  Sideband       :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmcsb, Nqcdmcsb, Newmcsb, Nsmsb, Nsb, Nsusymcsb ) ;
       printf("\n") ;
       printf("  A              :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmca, Nqcdmca, Newmca, Nsma, Na, Nsusymca ) ;
       printf("  D              :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmcd, Nqcdmcd, Newmcd, Nsmd, Nd, Nsusymcd ) ;
       printf("\n") ;
       printf("  SL, SIG        :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmcslsig, Nqcdmcslsig, Newmcslsig, Nsmslsig, Nslsig, Nsusymcslsig ) ;
       printf("  SL, SB         :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmcslsb, Nqcdmcslsb, Newmcslsb, Nsmslsb, Nslsb, Nsusymcslsb ) ;
       printf("\n") ;
       printf("  SL, MSB+SB+SIG :  %8.1f %8.1f %8.1f %8.1f %8d %8.1f\n",
             Nttbarmcslmsbsbsig, Nqcdmcslmsbsbsig, Newmcslmsbsbsig, Nsmmcslmsbsbsig, Nslmsbsbsig, Nsusymcslmsbsbsig ) ;
       printf("----------------------------------------------------------------------------\n") ;
       printf("\n") ;
       printf("    Rttmc     = %5.1f / %5.1f = %5.3f\n", Nttbarmcsig, Nttbarmcsb, Rttmc ) ;
       printf("    Rsl,mc    = %5.1f / %5.1f = %5.3f\n", Nsmslsig, Nsmslsb, Rslmc ) ;
       printf("    Rsl,data  = %5d / %5d = %5.3f +/- %5.3f\n", Nslsig, Nslsb, Rsldata, Rsldataerr ) ;
       printf("\n") ;
       printf("    Rqcd,mc   = %5.1f / %6.1f = %5.3f\n", Nqcdmcd, Nqcdmca, Rqcdmctruth ) ;
       printf("    Rqcd,data = %5.1f / %6.1f = %5.3f +/- %5.3f\n", (Nd-Nttbarmcd-Newmcd), (Na-Nttbarmca-Newmca), Rqcddata, Rqcddataerr ) ;
       printf("\n") ;
       printf("                              ttbar    qcd     EW\n") ;
       printf("-----------------------------------------------------\n") ;
       printf("   SIG region  fractions :   %5.2f   %5.2f   %5.2f\n", frsigttbar, frsigqcd, frsigew ) ;
       printf("   SB  region  fractions :   %5.2f   %5.2f   %5.2f\n", frsbttbar, frsbqcd, frsbew ) ;
       printf("   SL (MET>50) fractions :   %5.2f   %5.2f   %5.2f\n", frslttbar, frslqcd, frslew ) ;
       printf("\n\n\n") ;



       printf(" --- Defining observables.\n" ) ;

      //-- data counts in signal region, A, and D.

      rv_Nsig -> setVal( Nsig ) ;
      rv_Na -> setVal( Na ) ;
      rv_Nd -> setVal( Nd ) ;

      //-- data counts in 5 3-jet mass bins of SB.

      rv_Nsb1 -> setVal( Nsb1 ) ;
      rv_Nsb2 -> setVal( Nsb2 ) ;
      rv_Nsb3 -> setVal( Nsb3 ) ;
      rv_Nsb4 -> setVal( Nsb4 ) ;
      rv_Nsb5 -> setVal( Nsb5 ) ;

      //-- data counts in 5 3-jet mass bins of LSB.

      rv_Nlsb1 -> setVal( Nlsb1 ) ;
      rv_Nlsb2 -> setVal( Nlsb2 ) ;
      rv_Nlsb3 -> setVal( Nlsb3 ) ;
      rv_Nlsb4 -> setVal( Nlsb4 ) ;
      rv_Nlsb5 -> setVal( Nlsb5 ) ;


      //-- data counts in 5 3-jet mass bins of SL, SIG.

      rv_Nslsig1 -> setVal( Nslsig1 ) ;
      rv_Nslsig2 -> setVal( Nslsig2 ) ;
      rv_Nslsig3 -> setVal( Nslsig3 ) ;
      rv_Nslsig4 -> setVal( Nslsig4 ) ;
      rv_Nslsig5 -> setVal( Nslsig5 ) ;


      //-- data counts in 5 3-jet mass bins of SL, SB.

      rv_Nslsb1 -> setVal( Nslsb1 ) ;
      rv_Nslsb2 -> setVal( Nslsb2 ) ;
      rv_Nslsb3 -> setVal( Nslsb3 ) ;
      rv_Nslsb4 -> setVal( Nslsb4 ) ;
      rv_Nslsb5 -> setVal( Nslsb5 ) ;


      //-- data counts in 5 3-jet mass bins of SL, MSB.

      rv_Nslmsb1 -> setVal( Nslmsb1 ) ;
      rv_Nslmsb2 -> setVal( Nslmsb2 ) ;
      rv_Nslmsb3 -> setVal( Nslmsb3 ) ;
      rv_Nslmsb4 -> setVal( Nslmsb4 ) ;
      rv_Nslmsb5 -> setVal( Nslmsb5 ) ;





      //-- QCD MC counts

      rv_Nqcdmca   -> setVal( Nqcdmca ) ;
      rv_Nqcdmcd   -> setVal( Nqcdmcd ) ;
      rv_Nqcdmcsig -> setVal( Nqcdmcsig ) ;
      rv_Nqcdmcsb  -> setVal( Nqcdmcsb ) ;







      //++++++++ Parameters of the likelihood

       printf(" --- Resetting parameters.\n" ) ;

    //===================================================================================================
     //-- Use these lines to define things in terms of the SIG parameters.
      if ( useSigBgVars ) {

         rrv_mu_ttbar_sig -> setVal( Nttbarmcsig ) ;  //-- this is a starting value only.
         rrv_mu_qcd_sig   -> setVal( Nqcdmcsig ) ;  //-- this is a starting value only.

      } else {
    //===================================================================================================
     //-- Use these lines to define things in terms of the SB parameters.

         rrv_mu_ttbar_sb -> setVal( Nttbarmcsb ) ;  //-- this is a starting value only.
         rrv_mu_qcd_sb   -> setVal( Nqcdmcsb ) ;  //-- this is a starting value only.

      }
    //===================================================================================================


     //-- Counts in SIG, signal selection

      rv_mu_ew_sig    -> setVal( Newmcsig ) ;
      rv_mu_susy_sig  -> setVal( 0. ) ;  //-- this is a starting value only.

      rv_mu_ew_sig -> setConstant(kTRUE) ;

     //-- Counts in SB, bins of 3-jet mass, EW, signal selection

      rv_mu_ew_sb1 -> setVal( Newmcsb1 ) ;
      rv_mu_ew_sb2 -> setVal( Newmcsb2 ) ;
      rv_mu_ew_sb3 -> setVal( Newmcsb3 ) ;
      rv_mu_ew_sb4 -> setVal( Newmcsb4 ) ;
      rv_mu_ew_sb5 -> setVal( Newmcsb5 ) ;

      rv_mu_ew_sb1 -> setConstant(kTRUE) ;
      rv_mu_ew_sb2 -> setConstant(kTRUE) ;
      rv_mu_ew_sb3 -> setConstant(kTRUE) ;
      rv_mu_ew_sb4 -> setConstant(kTRUE) ;
      rv_mu_ew_sb5 -> setConstant(kTRUE) ;

     //-- Counts in SB, bins of 3-jet mass, SUSY, signal selection

      rv_mu_susy_sb1 -> setVal( 0. ) ;
      rv_mu_susy_sb2 -> setVal( 0. ) ;
      rv_mu_susy_sb3 -> setVal( 0. ) ;
      rv_mu_susy_sb4 -> setVal( 0. ) ;
      rv_mu_susy_sb5 -> setVal( 0. ) ;

      rv_mu_susy_sb1 -> setConstant(kTRUE) ;
      rv_mu_susy_sb2 -> setConstant(kTRUE) ;
      rv_mu_susy_sb3 -> setConstant(kTRUE) ;
      rv_mu_susy_sb4 -> setConstant(kTRUE) ;
      rv_mu_susy_sb5 -> setConstant(kTRUE) ;

     //-- Counts in A, signal selection

      rv_mu_ttbar_a -> setVal( Nttbarmca ) ;
      rv_mu_qcd_a   -> setVal( Na - Nttbarmca - Newmca ) ; //-- this is a starting value only.
      rv_mu_ew_a    -> setVal( Newmca ) ;
      rv_mu_susy_a  -> setVal( 0. ) ;

      rv_mu_susy_a  ->setConstant(kTRUE) ;
      rv_mu_ew_a    ->setConstant(kTRUE) ;
      rv_mu_ttbar_a ->setConstant(kTRUE) ;


     //-- Counts in D, signal selection

      rv_mu_ttbar_d -> setVal( Nttbarmcd ) ;
      rv_mu_qcd_d   -> setVal( Nd - Nttbarmcd - Newmcd ) ; //-- this is a starting value only.
      rv_mu_ew_d    -> setVal( Newmcd ) ;
      rv_mu_susy_d  -> setVal( 0. ) ;

      rv_mu_ew_d    -> setConstant(kTRUE) ;
      rv_mu_ttbar_d -> setConstant(kTRUE) ;
      rv_mu_susy_d  -> setConstant(kTRUE) ;

     //-- Counts in LSB, bins of 3-jet mass, qcd, signal selection

      rv_mu_qcd_lsb1 -> setVal( Nlsb1 ) ;
      rv_mu_qcd_lsb2 -> setVal( Nlsb2 ) ;
      rv_mu_qcd_lsb3 -> setVal( Nlsb3 ) ;
      rv_mu_qcd_lsb4 -> setVal( Nlsb4 ) ;
      rv_mu_qcd_lsb5 -> setVal( Nlsb5 ) ;


     //-- QCD MC, SIG,SB,A,D counts, signal selection

      rv_mu_qcdmc_sig -> setVal( Nqcdmcsig ) ;
      rv_mu_qcdmc_sb  -> setVal( Nqcdmcsb ) ;
      rv_mu_qcdmc_a   -> setVal( Nqcdmca ) ;
      rv_mu_qcdmc_d   -> setVal( Nqcdmcd ) ;

     //-- Single Lepton (SL) counts, ttbar, SIG region, bins of 3-jet mass

      rv_mu_sl_ttbar_sig1 -> setVal( Nslsig1 - Newmcslsig1 ) ;
      rv_mu_sl_ttbar_sig2 -> setVal( Nslsig2 - Newmcslsig2 ) ;
      rv_mu_sl_ttbar_sig3 -> setVal( Nslsig3 - Newmcslsig3 ) ;
      rv_mu_sl_ttbar_sig4 -> setVal( Nslsig4 - Newmcslsig4 ) ;
      rv_mu_sl_ttbar_sig5 -> setVal( Nslsig5 - Newmcslsig5 ) ;


     //-- Single Lepton (SL) counts, ttbar, SB region, bins of 3-jet mass

      rv_mu_sl_ttbar_sb1 -> setVal( Nslsb1 - Newmcslsb1 ) ;
      rv_mu_sl_ttbar_sb2 -> setVal( Nslsb2 - Newmcslsb2 ) ;
      rv_mu_sl_ttbar_sb3 -> setVal( Nslsb3 - Newmcslsb3 ) ;
      rv_mu_sl_ttbar_sb4 -> setVal( Nslsb4 - Newmcslsb4 ) ;
      rv_mu_sl_ttbar_sb5 -> setVal( Nslsb5 - Newmcslsb5 ) ;


     //-- Single Lepton (SL) counts, ttbar, MSB region, bins of 3-jet mass

      rv_mu_sl_ttbar_msb1 -> setVal( Nslmsb1 - Newmcslmsb1 ) ;
      rv_mu_sl_ttbar_msb2 -> setVal( Nslmsb2 - Newmcslmsb2 ) ;
      rv_mu_sl_ttbar_msb3 -> setVal( Nslmsb3 - Newmcslmsb3 ) ;
      rv_mu_sl_ttbar_msb4 -> setVal( Nslmsb4 - Newmcslmsb4 ) ;
      rv_mu_sl_ttbar_msb5 -> setVal( Nslmsb5 - Newmcslmsb5 ) ;







     //-- Single Lepton (SL) counts, EW, SIG region, bins of 3-jet mass

      rv_mu_sl_ew_sig1 -> setVal( Newmcslsig1 ) ;
      rv_mu_sl_ew_sig2 -> setVal( Newmcslsig2 ) ;
      rv_mu_sl_ew_sig3 -> setVal( Newmcslsig3 ) ;
      rv_mu_sl_ew_sig4 -> setVal( Newmcslsig4 ) ;
      rv_mu_sl_ew_sig5 -> setVal( Newmcslsig5 ) ;

      rv_mu_sl_ew_sig1 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sig2 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sig3 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sig4 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sig5 -> setConstant(kTRUE) ;



     //-- Single Lepton (SL) counts, EW, SB region, bins of 3-jet mass

      rv_mu_sl_ew_sb1 -> setVal( Newmcslsb1 ) ;
      rv_mu_sl_ew_sb2 -> setVal( Newmcslsb2 ) ;
      rv_mu_sl_ew_sb3 -> setVal( Newmcslsb3 ) ;
      rv_mu_sl_ew_sb4 -> setVal( Newmcslsb4 ) ;
      rv_mu_sl_ew_sb5 -> setVal( Newmcslsb5 ) ;

      rv_mu_sl_ew_sb1 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sb2 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sb3 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sb4 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_sb5 -> setConstant(kTRUE) ;



     //-- Single Lepton (SL) counts, EW, MSB region, bins of 3-jet mass

      rv_mu_sl_ew_msb1 -> setVal( Newmcslmsb1 ) ;
      rv_mu_sl_ew_msb2 -> setVal( Newmcslmsb2 ) ;
      rv_mu_sl_ew_msb3 -> setVal( Newmcslmsb3 ) ;
      rv_mu_sl_ew_msb4 -> setVal( Newmcslmsb4 ) ;
      rv_mu_sl_ew_msb5 -> setVal( Newmcslmsb5 ) ;

      rv_mu_sl_ew_msb1 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_msb2 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_msb3 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_msb4 -> setConstant(kTRUE) ;
      rv_mu_sl_ew_msb5 -> setConstant(kTRUE) ;








     //-- Single Lepton (SL) counts, susy, SIG region, bins of 3-jet mass

      rv_mu_sl_susy_sig1 -> setVal( 0. ) ;
      rv_mu_sl_susy_sig2 -> setVal( 0. ) ;
      rv_mu_sl_susy_sig3 -> setVal( 0. ) ;
      rv_mu_sl_susy_sig4 -> setVal( 0. ) ;
      rv_mu_sl_susy_sig5 -> setVal( 0. ) ;

      rv_mu_sl_susy_sig1 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sig2 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sig3 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sig4 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sig5 -> setConstant( kTRUE ) ;


     //-- Single Lepton (SL) counts, susy, SB region, bins of 3-jet mass

      rv_mu_sl_susy_sb1 -> setVal( 0. ) ;
      rv_mu_sl_susy_sb2 -> setVal( 0. ) ;
      rv_mu_sl_susy_sb3 -> setVal( 0. ) ;
      rv_mu_sl_susy_sb4 -> setVal( 0. ) ;
      rv_mu_sl_susy_sb5 -> setVal( 0. ) ;

      rv_mu_sl_susy_sb1 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sb2 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sb3 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sb4 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_sb5 -> setConstant( kTRUE ) ;


     //-- Single Lepton (SL) counts, susy, MSB region, bins of 3-jet mass

      rv_mu_sl_susy_msb1 -> setVal( 0. ) ;
      rv_mu_sl_susy_msb2 -> setVal( 0. ) ;
      rv_mu_sl_susy_msb3 -> setVal( 0. ) ;
      rv_mu_sl_susy_msb4 -> setVal( 0. ) ;
      rv_mu_sl_susy_msb5 -> setVal( 0. ) ;

      rv_mu_sl_susy_msb1 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_msb2 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_msb3 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_msb4 -> setConstant( kTRUE ) ;
      rv_mu_sl_susy_msb5 -> setConstant( kTRUE ) ;


       initialized = true ;

       return true ;


    } // reinitialize.

  //===================================================================

    bool ra2bRoostatsClass::readToyDataset( const char* inputRootFile, int dsIndex ) {

          if ( ! initialized ) {
             printf("\n\n *** Call initialize first.\n\n") ;
             return false ;
          }

         TFile intoyfile( inputRootFile, "READ" ) ;
         gDirectory->ls() ;
         TTree* toytree = (TTree*) gDirectory->FindObjectAny("likelihoodData") ;
         TTree* toytruetree = (TTree*) gDirectory->FindObjectAny("toytruetree") ;

         if ( toytree == 0 ) {
            printf(" \n\n *** Can't find TTree likelihoodData in file %s\n", inputRootFile ) ;
         } else {
            toytree->Print("toponly") ;
         }

         if ( toytruetree == 0 ) {
            printf(" \n\n *** Can't find TTree toytruetree in file %s\n", inputRootFile ) ;
         } else {
            toytruetree->Print("toponly") ;
         }


         Double_t toyNsig;
         Double_t toyNa;
         Double_t toyNd;
         Double_t toyNsb1;
         Double_t toyNsb2;
         Double_t toyNsb3;
         Double_t toyNsb4;
         Double_t toyNsb5;
         Double_t toyNlsb1;
         Double_t toyNlsb2;
         Double_t toyNlsb3;
         Double_t toyNlsb4;
         Double_t toyNlsb5;
         Double_t toyNslsig1;
         Double_t toyNslsig2;
         Double_t toyNslsig3;
         Double_t toyNslsig4;
         Double_t toyNslsig5;
         Double_t toyNslsb1;
         Double_t toyNslsb2;
         Double_t toyNslsb3;
         Double_t toyNslsb4;
         Double_t toyNslsb5;
         Double_t toyNslmsb1;
         Double_t toyNslmsb2;
         Double_t toyNslmsb3;
         Double_t toyNslmsb4;
         Double_t toyNslmsb5;

         Double_t toyNqcdmca ;
         Double_t toyNqcdmcd ;
         Double_t toyNqcdmcsig ;
         Double_t toyNqcdmcsb ;

         // List of branches
         TBranch *b_toyNsig;   //!
         TBranch *b_toyNa;   //!
         TBranch *b_toyNd;   //!
         TBranch *b_toyNsb1;   //!
         TBranch *b_toyNsb2;   //!
         TBranch *b_toyNsb3;   //!
         TBranch *b_toyNsb4;   //!
         TBranch *b_toyNsb5;   //!
         TBranch *b_toyNlsb1;   //!
         TBranch *b_toyNlsb2;   //!
         TBranch *b_toyNlsb3;   //!
         TBranch *b_toyNlsb4;   //!
         TBranch *b_toyNlsb5;   //!
         TBranch *b_toyNslsig1;   //!
         TBranch *b_toyNslsig2;   //!
         TBranch *b_toyNslsig3;   //!
         TBranch *b_toyNslsig4;   //!
         TBranch *b_toyNslsig5;   //!
         TBranch *b_toyNslsb1;   //!
         TBranch *b_toyNslsb2;   //!
         TBranch *b_toyNslsb3;   //!
         TBranch *b_toyNslsb4;   //!
         TBranch *b_toyNslsb5;   //!
         TBranch *b_toyNslmsb1;   //!
         TBranch *b_toyNslmsb2;   //!
         TBranch *b_toyNslmsb3;   //!
         TBranch *b_toyNslmsb4;   //!
         TBranch *b_toyNslmsb5;   //!

         TBranch *b_toyNqcdmca ;
         TBranch *b_toyNqcdmcd ;
         TBranch *b_toyNqcdmcsig ;
         TBranch *b_toyNqcdmcsb ;

         TBranch *b_toy_mu0_ttbar_sig ;
         TBranch *b_toy_mu0_qcd_sig ;
         TBranch *b_toy_mu0_ttbar_sb ;
         TBranch *b_toy_mu0_qcd_sb ;
         TBranch *b_toy_mu0_susy_sig ;
         TBranch *b_toy_mu0_allbg_sig ;

         printf("\n\n Setting branch addresses...\n") ;
         cout << " Toy tree pointer: " << toytree << endl ;
         toytree->SetBranchAddress("Nsig", &toyNsig, &b_toyNsig ) ;
         toytree->SetBranchAddress("Na", &toyNa, &b_toyNa);
         toytree->SetBranchAddress("Nd", &toyNd, &b_toyNd);
         toytree->SetBranchAddress("Nsb1", &toyNsb1, &b_toyNsb1);
         toytree->SetBranchAddress("Nsb2", &toyNsb2, &b_toyNsb2);
         toytree->SetBranchAddress("Nsb3", &toyNsb3, &b_toyNsb3);
         toytree->SetBranchAddress("Nsb4", &toyNsb4, &b_toyNsb4);
         toytree->SetBranchAddress("Nsb5", &toyNsb5, &b_toyNsb5);
         toytree->SetBranchAddress("Nlsb1", &toyNlsb1, &b_toyNlsb1);
         toytree->SetBranchAddress("Nlsb2", &toyNlsb2, &b_toyNlsb2);
         toytree->SetBranchAddress("Nlsb3", &toyNlsb3, &b_toyNlsb3);
         toytree->SetBranchAddress("Nlsb4", &toyNlsb4, &b_toyNlsb4);
         toytree->SetBranchAddress("Nlsb5", &toyNlsb5, &b_toyNlsb5);
         toytree->SetBranchAddress("Nslsig1", &toyNslsig1, &b_toyNslsig1);
         toytree->SetBranchAddress("Nslsig2", &toyNslsig2, &b_toyNslsig2);
         toytree->SetBranchAddress("Nslsig3", &toyNslsig3, &b_toyNslsig3);
         toytree->SetBranchAddress("Nslsig4", &toyNslsig4, &b_toyNslsig4);
         toytree->SetBranchAddress("Nslsig5", &toyNslsig5, &b_toyNslsig5);
         toytree->SetBranchAddress("Nslsb1", &toyNslsb1, &b_toyNslsb1);
         toytree->SetBranchAddress("Nslsb2", &toyNslsb2, &b_toyNslsb2);
         toytree->SetBranchAddress("Nslsb3", &toyNslsb3, &b_toyNslsb3);
         toytree->SetBranchAddress("Nslsb4", &toyNslsb4, &b_toyNslsb4);
         toytree->SetBranchAddress("Nslsb5", &toyNslsb5, &b_toyNslsb5);
         toytree->SetBranchAddress("Nslmsb1", &toyNslmsb1, &b_toyNslmsb1);
         toytree->SetBranchAddress("Nslmsb2", &toyNslmsb2, &b_toyNslmsb2);
         toytree->SetBranchAddress("Nslmsb3", &toyNslmsb3, &b_toyNslmsb3);
         toytree->SetBranchAddress("Nslmsb4", &toyNslmsb4, &b_toyNslmsb4);
         toytree->SetBranchAddress("Nslmsb5", &toyNslmsb5, &b_toyNslmsb5);

         toytree->SetBranchAddress("Nqcdmca"  , &toyNqcdmca  , &b_toyNqcdmca);
         toytree->SetBranchAddress("Nqcdmcd"  , &toyNqcdmcd  , &b_toyNqcdmcd);
         toytree->SetBranchAddress("Nqcdmcsig", &toyNqcdmcsig, &b_toyNqcdmcsig);
         toytree->SetBranchAddress("Nqcdmcsb" , &toyNqcdmcsb , &b_toyNqcdmcsb);

         toytruetree->SetBranchAddress("mu0_ttbar_sig" , &toy_mu0_ttbar_sig , &b_toy_mu0_ttbar_sig);
         toytruetree->SetBranchAddress("mu0_qcd_sig"   , &toy_mu0_qcd_sig   , &b_toy_mu0_qcd_sig);
         toytruetree->SetBranchAddress("mu0_ttbar_sb"  , &toy_mu0_ttbar_sb  , &b_toy_mu0_ttbar_sb);
         toytruetree->SetBranchAddress("mu0_qcd_sb"    , &toy_mu0_qcd_sb    , &b_toy_mu0_qcd_sb);
         toytruetree->SetBranchAddress("mu0_susy_sig"  , &toy_mu0_susy_sig  , &b_toy_mu0_susy_sig);
         toytruetree->SetBranchAddress("mu0_allbg_sig" , &toy_mu0_allbg_sig , &b_toy_mu0_allbg_sig);

         printf("\n Done.\n\n") ;


         toytruetree->GetEntry(0) ;

         toytree->GetEntry( dsIndex ) ;

         rv_Nsig -> setVal( toyNsig ) ;
         rv_Na -> setVal( toyNa ) ;
         rv_Nd -> setVal( toyNd ) ;
         rv_Nsb1 -> setVal( toyNsb1 ) ;
         rv_Nsb2 -> setVal( toyNsb2 ) ;
         rv_Nsb3 -> setVal( toyNsb3 ) ;
         rv_Nsb4 -> setVal( toyNsb4 ) ;
         rv_Nsb5 -> setVal( toyNsb5 ) ;
         rv_Nlsb1 -> setVal( toyNlsb1 ) ;
         rv_Nlsb2 -> setVal( toyNlsb2 ) ;
         rv_Nlsb3 -> setVal( toyNlsb3 ) ;
         rv_Nlsb4 -> setVal( toyNlsb4 ) ;
         rv_Nlsb5 -> setVal( toyNlsb5 ) ;
         rv_Nslsig1 -> setVal( toyNslsig1 ) ;
         rv_Nslsig2 -> setVal( toyNslsig2 ) ;
         rv_Nslsig3 -> setVal( toyNslsig3 ) ;
         rv_Nslsig4 -> setVal( toyNslsig4 ) ;
         rv_Nslsig5 -> setVal( toyNslsig5 ) ;
         rv_Nslsb1 -> setVal( toyNslsb1 ) ;
         rv_Nslsb2 -> setVal( toyNslsb2 ) ;
         rv_Nslsb3 -> setVal( toyNslsb3 ) ;
         rv_Nslsb4 -> setVal( toyNslsb4 ) ;
         rv_Nslsb5 -> setVal( toyNslsb5 ) ;
         rv_Nslmsb1 -> setVal( toyNslmsb1 ) ;
         rv_Nslmsb2 -> setVal( toyNslmsb2 ) ;
         rv_Nslmsb3 -> setVal( toyNslmsb3 ) ;
         rv_Nslmsb4 -> setVal( toyNslmsb4 ) ;
         rv_Nslmsb5 -> setVal( toyNslmsb5 ) ;

         rv_Nqcdmca   -> setVal( toyNqcdmca ) ;
         rv_Nqcdmcd   -> setVal( toyNqcdmcd ) ;
         rv_Nqcdmcsig -> setVal( toyNqcdmcsig ) ;
         rv_Nqcdmcsb  -> setVal( toyNqcdmcsb ) ;


         if ( dsObserved != 0x0 ) delete dsObserved ;

         dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
                                      observedParametersList ) ;
         dsObserved->add( observedParametersList ) ;
         printf("\n\n") ;
         dsObserved->printMultiline(cout, 1, kTRUE, "") ;
         printf("\n\n") ;


         return true ;

    } // readToyDataset

  //===================================================================



    bool ra2bRoostatsClass::readTextDataset( const char* infile ) {

          if ( ! initialized ) {
             printf("\n\n *** Call initialize first.\n\n") ;
             return false ;
          }


       printf( "\n\n Opening input file : %s\n\n", infile ) ;

       FILE* infp ;
       if ( (infp=fopen( infile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
          return false ;
       }

       int N3jmBins(0) ; //-- number of 3-jet mass bins.

       int Nsig(0), Na(0), Nd(0) ; //-- data counts in signal region, A, and D.
       int Nsb1(0), Nsb2(0), Nsb3(0), Nsb4(0), Nsb5(0) ; //-- data counts in 5 3-jet mass bins of SB.
       int Nlsb1(0), Nlsb2(0), Nlsb3(0), Nlsb4(0), Nlsb5(0) ; //-- data counts in 5 3-jet mass bins of LSB.
       int Nslsig1(0), Nslsig2(0), Nslsig3(0), Nslsig4(0), Nslsig5(0) ; //-- data counts in 5 3-jet mass bins of SL, SIG.
       int Nslsb1(0), Nslsb2(0), Nslsb3(0), Nslsb4(0), Nslsb5(0) ; ; //-- data counts in 5 3-jet mass bins of SL, SB.
       int Nslmsb1(0), Nslmsb2(0), Nslmsb3(0), Nslmsb4(0), Nslmsb5(0) ; //-- data counts in 5 3-jet mass bins of SL, MSB.


       //--- read in description line.
       printf("\n\n") ;
       char c(0) ;
       while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
       printf("\n\n") ;


       char label[1000] ;

      //--- Inputs generated with gen_roostats_input.c
      //    The order here must be consistent with the order there!

       fscanf( infp, "%s %d", label, &N3jmBins ) ;             printf( "%s %d\n", label, N3jmBins ) ;                
       fscanf( infp, "%s %d", label, &Nsig ) ;                 printf( "%s %d\n", label, Nsig ) ;         
       fscanf( infp, "%s %d", label, &Na ) ;                   printf( "%s %d\n", label, Na ) ;           
       fscanf( infp, "%s %d", label, &Nd ) ;                   printf( "%s %d\n", label, Nd ) ;           
       fscanf( infp, "%s %d", label, &Nsb1 ) ;                 printf( "%s %d\n", label, Nsb1 ) ;         
       fscanf( infp, "%s %d", label, &Nsb2 ) ;                 printf( "%s %d\n", label, Nsb2 ) ;         
       fscanf( infp, "%s %d", label, &Nsb3 ) ;                 printf( "%s %d\n", label, Nsb3 ) ;         
       fscanf( infp, "%s %d", label, &Nsb4 ) ;                 printf( "%s %d\n", label, Nsb4 ) ;         
       fscanf( infp, "%s %d", label, &Nsb5 ) ;                 printf( "%s %d\n", label, Nsb5 ) ;         
       fscanf( infp, "%s %d", label, &Nlsb1 ) ;                printf( "%s %d\n", label, Nlsb1 ) ;        
       fscanf( infp, "%s %d", label, &Nlsb2 ) ;                printf( "%s %d\n", label, Nlsb2 ) ;        
       fscanf( infp, "%s %d", label, &Nlsb3 ) ;                printf( "%s %d\n", label, Nlsb3 ) ;        
       fscanf( infp, "%s %d", label, &Nlsb4 ) ;                printf( "%s %d\n", label, Nlsb4 ) ;        
       fscanf( infp, "%s %d", label, &Nlsb5 ) ;                printf( "%s %d\n", label, Nlsb5 ) ;        
       fscanf( infp, "%s %d", label, &Nslsig1 ) ;              printf( "%s %d\n", label, Nslsig1 ) ;      
       fscanf( infp, "%s %d", label, &Nslsig2 ) ;              printf( "%s %d\n", label, Nslsig2 ) ;      
       fscanf( infp, "%s %d", label, &Nslsig3 ) ;              printf( "%s %d\n", label, Nslsig3 ) ;      
       fscanf( infp, "%s %d", label, &Nslsig4 ) ;              printf( "%s %d\n", label, Nslsig4 ) ;      
       fscanf( infp, "%s %d", label, &Nslsig5 ) ;              printf( "%s %d\n", label, Nslsig5 ) ;      
       fscanf( infp, "%s %d", label, &Nslsb1 ) ;               printf( "%s %d\n", label, Nslsb1 ) ;       
       fscanf( infp, "%s %d", label, &Nslsb2 ) ;               printf( "%s %d\n", label, Nslsb2 ) ;       
       fscanf( infp, "%s %d", label, &Nslsb3 ) ;               printf( "%s %d\n", label, Nslsb3 ) ;       
       fscanf( infp, "%s %d", label, &Nslsb4 ) ;               printf( "%s %d\n", label, Nslsb4 ) ;       
       fscanf( infp, "%s %d", label, &Nslsb5 ) ;               printf( "%s %d\n", label, Nslsb5 ) ;       
       fscanf( infp, "%s %d", label, &Nslmsb1 ) ;              printf( "%s %d\n", label, Nslmsb1 ) ;      
       fscanf( infp, "%s %d", label, &Nslmsb2 ) ;              printf( "%s %d\n", label, Nslmsb2 ) ;      
       fscanf( infp, "%s %d", label, &Nslmsb3 ) ;              printf( "%s %d\n", label, Nslmsb3 ) ;      
       fscanf( infp, "%s %d", label, &Nslmsb4 ) ;              printf( "%s %d\n", label, Nslmsb4 ) ;      
       fscanf( infp, "%s %d", label, &Nslmsb5 ) ;              printf( "%s %d\n", label, Nslmsb5 ) ;      

       rv_Nsig -> setVal( Nsig ) ;
       rv_Na -> setVal( Na ) ;
       rv_Nd -> setVal( Nd ) ;
       rv_Nsb1 -> setVal( Nsb1 ) ;
       rv_Nsb2 -> setVal( Nsb2 ) ;
       rv_Nsb3 -> setVal( Nsb3 ) ;
       rv_Nsb4 -> setVal( Nsb4 ) ;
       rv_Nsb5 -> setVal( Nsb5 ) ;
       rv_Nlsb1 -> setVal( Nlsb1 ) ;
       rv_Nlsb2 -> setVal( Nlsb2 ) ;
       rv_Nlsb3 -> setVal( Nlsb3 ) ;
       rv_Nlsb4 -> setVal( Nlsb4 ) ;
       rv_Nlsb5 -> setVal( Nlsb5 ) ;
       rv_Nslsig1 -> setVal( Nslsig1 ) ;
       rv_Nslsig2 -> setVal( Nslsig2 ) ;
       rv_Nslsig3 -> setVal( Nslsig3 ) ;
       rv_Nslsig4 -> setVal( Nslsig4 ) ;
       rv_Nslsig5 -> setVal( Nslsig5 ) ;
       rv_Nslsb1 -> setVal( Nslsb1 ) ;
       rv_Nslsb2 -> setVal( Nslsb2 ) ;
       rv_Nslsb3 -> setVal( Nslsb3 ) ;
       rv_Nslsb4 -> setVal( Nslsb4 ) ;
       rv_Nslsb5 -> setVal( Nslsb5 ) ;
       rv_Nslmsb1 -> setVal( Nslmsb1 ) ;
       rv_Nslmsb2 -> setVal( Nslmsb2 ) ;
       rv_Nslmsb3 -> setVal( Nslmsb3 ) ;
       rv_Nslmsb4 -> setVal( Nslmsb4 ) ;
       rv_Nslmsb5 -> setVal( Nslmsb5 ) ;


       if ( dsObserved != 0x0 ) delete dsObserved ;

       dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
                                    observedParametersList ) ;
       dsObserved->add( observedParametersList ) ;
       printf("\n\n") ;
       dsObserved->printMultiline(cout, 1, kTRUE, "") ;
       printf("\n\n") ;


       return true ;



    } // readTextDataset


  //===================================================================


    bool ra2bRoostatsClass::doToyStudy( const char* inputRootFile, const char* outputRootFile, int dsFirst, int nToys ) {


         double toyNsig ;
         double toy_mu_ttbar_sig ;
         double toy_mu_qcd_sig ;
         double toy_mu_ttbar_sb ;
         double toy_mu_qcd_sb ;

         double toy_mu_ttbar_sig_err ;
         double toy_mu_qcd_sig_err ;
         double toy_mu_ttbar_sb_err ;
         double toy_mu_qcd_sb_err ;

         double toy_rho_ttbar_qcd_sig ;
         double toy_rho_ttbar_qcd_sb ;

         double toy_mu_susy_sig ;
         double toy_mu_susy_sig_err ;
         double toy_mu_susy_sig_ul ;

         double toy_mu_allbg_sig ;
         double toy_mu_allbg_sig_err ;


         TTree* toyfittree = new TTree("toyfittree","ra2b toy fit tree") ;

         toyfittree->Branch("Nsig"              , &toyNsig, "toyNsig/D" ) ;
         toyfittree->Branch("mu_ttbar_sig"      , &toy_mu_ttbar_sig, "mu_ttbar_sig/D" ) ;
         toyfittree->Branch("mu_qcd_sig"        , &toy_mu_qcd_sig, "mu_qcd_sig/D" ) ;
         toyfittree->Branch("mu_ttbar_sb"       , &toy_mu_ttbar_sb, "mu_ttbar_sb/D" ) ;
         toyfittree->Branch("mu_qcd_sb"         , &toy_mu_qcd_sb, "mu_qcd_sb/D" ) ;
         toyfittree->Branch("mu_ttbar_sig_err"  , &toy_mu_ttbar_sig_err, "mu_ttbar_sig_err/D" ) ;
         toyfittree->Branch("mu_qcd_sig_err"    , &toy_mu_qcd_sig_err, "mu_qcd_sig_err/D" ) ;
         toyfittree->Branch("mu_ttbar_sb_err"   , &toy_mu_ttbar_sb_err, "mu_ttbar_sb_err/D" ) ;
         toyfittree->Branch("mu_qcd_sb_err"     , &toy_mu_qcd_sb_err, "mu_qcd_sb_err/D" ) ;

         toyfittree->Branch("rho_ttbar_qcd_sig" , &toy_rho_ttbar_qcd_sig     , "rho_ttbar_qcd_sig/D" ) ;
         toyfittree->Branch("rho_ttbar_qcd_sb"  , &toy_rho_ttbar_qcd_sb      , "rho_ttbar_qcd_sb/D" ) ;

         toyfittree->Branch("mu_susy_sig"       , &toy_mu_susy_sig       , "mu_susy_sig/D" ) ;
         toyfittree->Branch("mu_susy_sig_err"   , &toy_mu_susy_sig_err   , "mu_susy_sig_err/D" ) ;
         toyfittree->Branch("mu_susy_sig_ul"    , &toy_mu_susy_sig_ul    , "mu_susy_sig_ul/D" ) ;
         toyfittree->Branch("mu_allbg_sig"      , &toy_mu_allbg_sig      , "mu_allbg_sig/D" ) ;
         toyfittree->Branch("mu_allbg_sig_err"  , &toy_mu_allbg_sig_err  , "mu_allbg_sig_err/D" ) ;

         toyfittree->Branch("mu0_ttbar_sig"      , &toy_mu0_ttbar_sig, "mu0_ttbar_sig/D" ) ;
         toyfittree->Branch("mu0_qcd_sig"        , &toy_mu0_qcd_sig, "mu0_qcd_sig/D" ) ;
         toyfittree->Branch("mu0_ttbar_sb"       , &toy_mu0_ttbar_sb, "mu0_ttbar_sb/D" ) ;
         toyfittree->Branch("mu0_qcd_sb"         , &toy_mu0_qcd_sb, "mu0_qcd_sb/D" ) ;
         toyfittree->Branch("mu0_susy_sig"       , &toy_mu0_susy_sig       , "mu0_susy_sig/D" ) ;
         toyfittree->Branch("mu0_allbg_sig"      , &toy_mu0_allbg_sig      , "mu0_allbg_sig/D" ) ;


         for ( int dsi=dsFirst; dsi<(dsFirst+nToys); dsi++ ) {

            printf("\n\n\n ============ Begin fit of dataset %d ====================== \n\n", dsi ) ;

            reinitialize() ;

            readToyDataset( inputRootFile, dsi ) ;

            doFit() ;

            toyNsig = rv_Nsig->getVal() ;

            if ( useSigBgVars ) {
               toy_mu_ttbar_sig = ((RooRealVar*)rv_mu_ttbar_sig)->getVal() ;
               toy_mu_qcd_sig = ((RooRealVar*)rv_mu_qcd_sig)->getVal() ;
               toy_mu_ttbar_sb = ((RooFormulaVar*)rv_mu_ttbar_sb)->getVal() ;
               toy_mu_qcd_sb = ((RooFormulaVar*)rv_mu_qcd_sb)->getVal() ;
               toy_mu_ttbar_sig_err = ((RooRealVar*)rv_mu_ttbar_sig)->getError() ;
               toy_mu_qcd_sig_err = ((RooRealVar*)rv_mu_qcd_sig)->getError() ;
               toy_mu_ttbar_sb_err = 0. ;
               toy_mu_qcd_sb_err = 0. ;
            } else {
               toy_mu_ttbar_sig = ((RooFormulaVar*)rv_mu_ttbar_sig)->getVal() ;
               toy_mu_qcd_sig = ((RooFormulaVar*)rv_mu_qcd_sig)->getVal() ;
               toy_mu_ttbar_sb = ((RooRealVar*)rv_mu_ttbar_sb)->getVal() ;
               toy_mu_qcd_sb = ((RooRealVar*)rv_mu_qcd_sb)->getVal() ;
               toy_mu_ttbar_sb_err = ((RooRealVar*)rv_mu_ttbar_sb)->getError() ;
               toy_mu_qcd_sb_err = ((RooRealVar*)rv_mu_qcd_sb)->getError() ;
               toy_mu_ttbar_sig_err = 0. ;
               toy_mu_qcd_sig_err = 0. ;
            }

            toy_mu_susy_sig = rv_mu_susy_sig->getVal() ;
            toy_mu_susy_sig_err = rv_mu_susy_sig->getError() ;

            toy_rho_ttbar_qcd_sig = fitResult->correlation( "mu_ttbar_sig", "mu_qcd_sig" ) ;
            toy_rho_ttbar_qcd_sb  = fitResult->correlation( "mu_ttbar_sb" , "mu_qcd_sb"  ) ;

       //--- Profile likelihood for signal susy yield.

            ProfileLikelihoodCalculator plc_susy_sig( *dsObserved, *likelihood, RooArgSet( *rv_mu_susy_sig ) ) ;
            plc_susy_sig.SetTestSize(0.05) ;
            ConfInterval* plinterval_susy_sig = plc_susy_sig.GetInterval() ;
            float susy_sig_ul = ((LikelihoodInterval*) plinterval_susy_sig)->UpperLimit(*((RooRealVar*)rv_mu_susy_sig)) ;
            float susy_sig_ll = ((LikelihoodInterval*) plinterval_susy_sig)->LowerLimit(*((RooRealVar*)rv_mu_susy_sig)) ;
            printf("\n\n") ;
            printf("    susy, SIG 95%% CL interval  [%5.1f, %5.1f]\n\n", susy_sig_ll, susy_sig_ul) ;
            printf("\n\n") ;
            toy_mu_susy_sig_ul = susy_sig_ul ;
            delete plinterval_susy_sig ;


            toy_mu_allbg_sig = toy_mu_ttbar_sig + toy_mu_qcd_sig + rv_mu_ew_sig->getVal() ;
            { double err2 = pow(toy_mu_ttbar_sig_err,2) + pow(toy_mu_qcd_sig_err,2)
                         + 2.*toy_rho_ttbar_qcd_sig*toy_mu_ttbar_sig_err*toy_mu_qcd_sig_err ;
               toy_mu_allbg_sig_err =  0. ;
               if ( err2 > 0. ) toy_mu_allbg_sig_err = sqrt( err2 ) ;
            }

            toyfittree->Fill() ;

            printf("\n\n\n ============ End fit of dataset %d ====================== \n\n", dsi ) ;

         } // dsi.

         TFile outputfile(outputRootFile,"recreate") ;
         toyfittree->Write() ;
         outputfile.Close() ;

         return true ;


    } // doToyStudy

  //===================================================================





















