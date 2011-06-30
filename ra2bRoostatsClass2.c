
//
//   Owen Long, UCR
//   Harrison Prosper, FSU
//   Sezen Sekmen, FSU
//
//


#include "ra2bRoostatsClass2.h"

#include <iostream>


#include "TCanvas.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TText.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"
#include "TH2F.h"
#include "TGaxis.h"
#include "TLine.h"

#include "RooArgSet.h"
#include "RooConstVar.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "LikelihoodIntervalPlot.cxx"

  using namespace RooFit ;
  using namespace RooStats ;


  //=====================================================================================================


   ra2bRoostatsClass2::ra2bRoostatsClass2( bool ArgUseSigBgVars ) {

      gStyle->SetOptStat(0) ;

      useSigBgVars = ArgUseSigBgVars ;

     //--- Tell RooFit to shut up about anything less important than an ERROR.
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
  /// RooMsgService::instance().addStream(DEBUG,Topic(Tracing),OutputFile("debug.log")) ;
      printf("\n\n ==== RooFit output configuration ===============\n") ;
      RooMsgService::instance().Print("v") ;
      printf("\n\n ================================================\n") ;

      varsAtFitVals = false ;
      initialized = false ;

   }




////==================================================================================================

//  bool ra2bRoostatsClass2::generateToyDatasetsFromLikelihood( const char* outputRootFile, int nToys ) {

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
//     if ( useSigBgVars ) {
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

//  bool ra2bRoostatsClass2::generateToyDatasetsFromNVals( const char* outputRootFile, int nToys ) {

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
//     if ( useSigBgVars ) {
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

    bool ra2bRoostatsClass2::doFit( ) {

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



//     float Amc = rv_mu_qcdmc_a->getVal() ;
//     float Dmc = rv_mu_qcdmc_d->getVal() ;
//     float Bmc = rv_mu_qcdmc_sb->getVal() ;
//     float Cmc = rv_mu_qcdmc_sig->getVal() ;

//     float K = Amc * Cmc / ( Bmc * Dmc ) ;

//     qcdCorrection = K ;

//     float covAA = pow( rv_mu_qcdmc_a->getError(), 2 ) ;
//     float covDD = pow( rv_mu_qcdmc_d->getError(), 2 ) ;
//     float covBB = pow( rv_mu_qcdmc_sb->getError(), 2 ) ;
//     float covCC = pow( rv_mu_qcdmc_sig->getError(), 2 ) ;

//     float rhoAB = fitResult->correlation( "mu_qcdmc_a", "mu_qcdmc_sb" ) ;
//     float rhoAC = fitResult->correlation( "mu_qcdmc_a", "mu_qcdmc_sig" ) ;
//     float rhoAD = fitResult->correlation( "mu_qcdmc_a", "mu_qcdmc_d" ) ;

//     float rhoBC = fitResult->correlation( "mu_qcdmc_sb", "mu_qcdmc_sig" ) ;
//     float rhoBD = fitResult->correlation( "mu_qcdmc_sb", "mu_qcdmc_d" ) ;

//     float rhoCD = fitResult->correlation( "mu_qcdmc_sig", "mu_qcdmc_d" ) ;

//     float covAB = rhoAB * sqrt( covAA * covBB ) ;
//     float covAC = rhoAC * sqrt( covAA * covCC ) ;
//     float covAD = rhoAD * sqrt( covAA * covDD ) ;

//     float covBC = rhoBC * sqrt( covBB * covCC ) ;
//     float covBD = rhoBD * sqrt( covBB * covDD ) ;

//     float covCD = rhoCD * sqrt( covCC * covDD ) ;

//     qcdCorrectionErr = K * sqrt( 
//            covAA/(Amc*Amc) + covBB/(Bmc*Bmc) + covCC/(Cmc*Cmc) + covDD/(Dmc*Dmc)
//            -2 * covAB / (Amc*Bmc) +2 * covAC / (Amc*Cmc) -2 * covAD / (Amc*Dmc)
//            +2 * covBC / (Bmc*Cmc) -2 * covBD / (Bmc*Dmc) 
//            -2 * covCD / (Cmc*Dmc)
//     ) ;


//     printf("  QCD bias correction: %4.2f +/- %4.2f\n", qcdCorrection, qcdCorrectionErr ) ;

       varsAtFitVals = true ;

       return true ;

     } // doFit .


  //==================================================================================================

     bool ra2bRoostatsClass2::profileSusySig( float& susySigLow, float& susySigHigh, bool makePlot ) {

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

         if ( makePlot ) {
            TCanvas* plcplot_susy_sig = new TCanvas("plcplot_susy_sig", "susy sig, Profile likelihood", 500, 400 ) ;
            LikelihoodIntervalPlot plotInt_susy_sig((LikelihoodInterval*)plinterval_susy_sig);
            plotInt_susy_sig.Draw() ;
            plcplot_susy_sig->SaveAs("plscan_susy_sig.pdf") ;
            plcplot_susy_sig->SaveAs("plscan_susy_sig.png") ;
         }

         varsAtFitVals = false ;

         delete plinterval_susy_sig ; // can I safely do this???

         return true ;

     }

  //==================================================================================================

     bool ra2bRoostatsClass2::profileTtbarSig( float& ttbarSigLow, float& ttbarSigHigh ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }

         if ( ! useSigBgVars ) {
            printf("\n\n\n *** Can't do it.  Need to use ra2bRoostatsClass2 with constructor arg useSigBgVars set to true.\n\n") ;
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

     bool ra2bRoostatsClass2::profileQcdSig( float& qcdSigLow, float& qcdSigHigh ) {

         if ( ! initialized ) {
            printf("\n\n *** Call initialize first.\n\n") ;
            return false ;
         }


      //--- Profile likelihood for SIG qcd yield.

         ProfileLikelihoodCalculator plc_qcd_sig( *dsObserved, *likelihood, RooArgSet( *rv_mu_qcd_sig ) ) ;
         plc_qcd_sig.SetTestSize(0.32) ;
         ConfInterval* plinterval_qcd_sig = plc_qcd_sig.GetInterval() ;
         qcdSigLow  = ((LikelihoodInterval*) plinterval_qcd_sig)->LowerLimit(*rv_mu_qcd_sig) ;
         qcdSigHigh = ((LikelihoodInterval*) plinterval_qcd_sig)->UpperLimit(*rv_mu_qcd_sig) ;
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


//   bool ra2bRoostatsClass2::profileTtbarSb( float& ttbarSbLow, float& ttbarSbHigh ) {

//       if ( ! initialized ) {
//          printf("\n\n *** Call initialize first.\n\n") ;
//          return false ;
//       }

//       if ( useSigBgVars ) {
//          printf("\n\n\n *** Can't do it.  Need to use ra2bRoostatsClass2 with constructor arg useSigBgVars set to false.\n\n") ;
//          return false ;
//       }

//    //--- Profile likelihood for SB ttbar yield.

//       ProfileLikelihoodCalculator plc_ttbar_sb( *dsObserved, *likelihood, RooArgSet( *rv_mu_ttbar_sb ) ) ;
//       plc_ttbar_sb.SetTestSize(0.32) ;
//       ConfInterval* plinterval_ttbar_sb = plc_ttbar_sb.GetInterval() ;
//       ttbarSbLow  = ((LikelihoodInterval*) plinterval_ttbar_sb)->LowerLimit(*rrv_mu_ttbar_sb) ;
//       ttbarSbHigh = ((LikelihoodInterval*) plinterval_ttbar_sb)->UpperLimit(*rrv_mu_ttbar_sb) ;
//       printf("\n\n") ;
//       printf("    ttbar, SB 68%% CL interval  [%5.1f, %5.1f]\n\n", ttbarSbLow, ttbarSbHigh ) ;
//       TCanvas* plcplot_ttbar_sb = new TCanvas("plcplot_ttbar_sb", "ttbar sb, Profile likelihood", 500, 400 ) ;
//       LikelihoodIntervalPlot plotInt_ttbar_sb((LikelihoodInterval*)plinterval_ttbar_sb);
//       plotInt_ttbar_sb.Draw() ;
//       plcplot_ttbar_sb->SaveAs("plscan_ttbar_sb.pdf") ;
//       plcplot_ttbar_sb->SaveAs("plscan_ttbar_sb.png") ;

//       varsAtFitVals = false ;

//       return true ;

//   }

////==================================================================================================


//   bool ra2bRoostatsClass2::profileQcdSb( float& qcdSbLow, float& qcdSbHigh ) {

//       if ( ! initialized ) {
//          printf("\n\n *** Call initialize first.\n\n") ;
//          return false ;
//       }

//       if ( useSigBgVars ) {
//          printf("\n\n\n *** Can't do it.  Need to use ra2bRoostatsClass2 with constructor arg useSigBgVars set to false.\n\n") ;
//          return false ;
//       }

//    //--- Profile likelihood for SB qcd yield.

//       ProfileLikelihoodCalculator plc_qcd_sb( *dsObserved, *likelihood, RooArgSet( *rv_mu_qcd_sb ) ) ;
//       plc_qcd_sb.SetTestSize(0.32) ;
//       ConfInterval* plinterval_qcd_sb = plc_qcd_sb.GetInterval() ;
//       qcdSbLow  = ((LikelihoodInterval*) plinterval_qcd_sb)->LowerLimit(*rrv_mu_qcd_sb) ;
//       qcdSbHigh = ((LikelihoodInterval*) plinterval_qcd_sb)->UpperLimit(*rrv_mu_qcd_sb) ;
//       printf("\n\n") ;
//       printf("    qcd, SB 68%% CL interval  [%5.1f, %5.1f]\n\n", qcdSbLow, qcdSbHigh ) ;
//       TCanvas* plcplot_qcd_sb = new TCanvas("plcplot_qcd_sb", "qcd sb, Profile likelihood", 500, 400 ) ;
//       LikelihoodIntervalPlot plotInt_qcd_sb((LikelihoodInterval*)plinterval_qcd_sb);
//       plotInt_qcd_sb.Draw() ;
//       plcplot_qcd_sb->SaveAs("plscan_qcd_sb.pdf") ;
//       plcplot_qcd_sb->SaveAs("plscan_qcd_sb.png") ;

//       varsAtFitVals = false ;

//       return true ;

//   }

////==================================================================================================



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



//  //====================================================================================================================

//     bool ra2bRoostatsClass2::sbPlotsUniformBins( const char* plotBaseName ) {

//        if ( ! initialized ) {
//           printf("\n\n *** Call initialize first.\n\n") ;
//           return false ;
//        }


//        if ( ! varsAtFitVals ) {
//           printf("\n\n *** Try this right after calling doFit.\n\n") ;
//           return false ;
//        }

//     //--  Drawn with same-width bins.

//        TH1F* hm3j_sb_data  = new TH1F("hm3j_sb_data" ,"3-jet mass, SB data" , 5, 0.5, 5.5 ) ;
//        TH1F* hm3j_sb_ttbar = new TH1F("hm3j_sb_ttbar","3-jet mass, SB ttbar", 5, 0.5, 5.5 ) ;
//        TH1F* hm3j_sb_qcd   = new TH1F("hm3j_sb_qcd"  ,"3-jet mass, SB qcd"  , 5, 0.5, 5.5 ) ;
//        TH1F* hm3j_sb_ew    = new TH1F("hm3j_sb_ew"   ,"3-jet mass, SB ew"   , 5, 0.5, 5.5 ) ;

//        hm3j_sb_data->SetLineWidth(2) ;
//        hm3j_sb_data->SetMarkerStyle(20) ;
//        hm3j_sb_ew->SetFillColor(49) ;
//        hm3j_sb_qcd->SetFillColor(46) ;
//        hm3j_sb_ttbar->SetFillColor(42) ;

//        hm3j_sb_data->SetBinContent( 1, rv_Nsb1->getVal() ) ;
//        hm3j_sb_data->SetBinContent( 2, rv_Nsb2->getVal() ) ;
//        hm3j_sb_data->SetBinContent( 3, rv_Nsb3->getVal() ) ;
//        hm3j_sb_data->SetBinContent( 4, rv_Nsb4->getVal() ) ;
//        hm3j_sb_data->SetBinContent( 5, rv_Nsb5->getVal() ) ;

//        hm3j_sb_ttbar->SetBinContent( 1, rv_mu_ttbar_sb1->getVal() ) ;
//        hm3j_sb_ttbar->SetBinContent( 2, rv_mu_ttbar_sb2->getVal() ) ;
//        hm3j_sb_ttbar->SetBinContent( 3, rv_mu_ttbar_sb3->getVal() ) ;
//        hm3j_sb_ttbar->SetBinContent( 4, rv_mu_ttbar_sb4->getVal() ) ;
//        hm3j_sb_ttbar->SetBinContent( 5, rv_mu_ttbar_sb5->getVal() ) ;

//        hm3j_sb_qcd->SetBinContent( 1, rv_mu_qcd_sb1->getVal() ) ;
//        hm3j_sb_qcd->SetBinContent( 2, rv_mu_qcd_sb2->getVal() ) ;
//        hm3j_sb_qcd->SetBinContent( 3, rv_mu_qcd_sb3->getVal() ) ;
//        hm3j_sb_qcd->SetBinContent( 4, rv_mu_qcd_sb4->getVal() ) ;
//        hm3j_sb_qcd->SetBinContent( 5, rv_mu_qcd_sb5->getVal() ) ;

//        hm3j_sb_ew->SetBinContent( 1, rv_mu_ew_sb1->getVal() ) ;
//        hm3j_sb_ew->SetBinContent( 2, rv_mu_ew_sb2->getVal() ) ;
//        hm3j_sb_ew->SetBinContent( 3, rv_mu_ew_sb3->getVal() ) ;
//        hm3j_sb_ew->SetBinContent( 4, rv_mu_ew_sb4->getVal() ) ;
//        hm3j_sb_ew->SetBinContent( 5, rv_mu_ew_sb5->getVal() ) ;

//        THStack* hstack_m3j_sb_fit = new THStack( "hstack_m3j_sb_fit", "SB fit, 3-jet mass" ) ;
//        hstack_m3j_sb_fit->Add( hm3j_sb_ew ) ;
//        hstack_m3j_sb_fit->Add( hm3j_sb_qcd ) ;
//        hstack_m3j_sb_fit->Add( hm3j_sb_ttbar ) ;


//        TCanvas* can_sbfit = new TCanvas("can_sbfit", "SB 3-jet mass fit", 700, 500 ) ;

//        hm3j_sb_data->SetMaximum( 1.4*(hm3j_sb_data->GetMaximum()) ) ;
//        hm3j_sb_data->SetLabelSize(0.06,"x") ;
//        hm3j_sb_data->GetXaxis()->SetBinLabel(1, "Bin 1") ;
//        hm3j_sb_data->GetXaxis()->SetBinLabel(2, "Bin 2") ;
//        hm3j_sb_data->GetXaxis()->SetBinLabel(3, "Bin 3") ;
//        hm3j_sb_data->GetXaxis()->SetBinLabel(4, "Bin 4") ;
//        hm3j_sb_data->GetXaxis()->SetBinLabel(5, "Bin 5") ;


//        hm3j_sb_data->Draw("histpe") ;
//        hstack_m3j_sb_fit->Draw( "same" ) ;
//        hm3j_sb_data->Draw("samehistpe") ;

//        TLegend* m3j_legend = new TLegend(0.62,0.7,0.97,0.95) ;
//        m3j_legend->SetFillColor( kWhite ) ;
//        m3j_legend->AddEntry( hm3j_sb_data, "Data" ) ;
//        m3j_legend->AddEntry( hm3j_sb_ttbar, "ttbar" ) ;
//        m3j_legend->AddEntry( hm3j_sb_qcd, "QCD" ) ;
//        m3j_legend->AddEntry( hm3j_sb_ew, "EW" ) ;
//        m3j_legend->Draw() ;

//        TText* fittext = new TText() ;
//        fittext->SetTextSize(0.04) ;
//        char fitlabel[1000] ;

//        float ew_sb(0.) ;
//        ew_sb += rv_mu_ew_sb1->getVal() ;
//        ew_sb += rv_mu_ew_sb2->getVal() ;
//        ew_sb += rv_mu_ew_sb3->getVal() ;
//        ew_sb += rv_mu_ew_sb4->getVal() ;
//        ew_sb += rv_mu_ew_sb5->getVal() ;
//        int Nsb = rv_Nsb1->getVal() +rv_Nsb2->getVal() +rv_Nsb3->getVal() +rv_Nsb4->getVal() +rv_Nsb5->getVal()  ;

//        float ttbar_sb_val(0.) ;
//        float qcd_sb_val(0.) ;
//        float ttbar_sb_err(0.) ;
//        float qcd_sb_err(0.) ;

//        float ttbar_sig_val(0.) ;
//        float qcd_sig_val(0.) ;
//        float ttbar_sig_err(0.) ;
//        float qcd_sig_err(0.) ;

//        if ( useSigBgVars ) {
//           ttbar_sb_val = ((RooFormulaVar*)rv_mu_ttbar_sb) ->getVal() ;
//           qcd_sb_val = ((RooFormulaVar*)rv_mu_qcd_sb) ->getVal() ;
//           ttbar_sig_val = ((RooRealVar*)rv_mu_ttbar_sig) ->getVal() ;
//           qcd_sig_val = ((RooRealVar*)rv_mu_qcd_sig) ->getVal() ;
//           ttbar_sig_err = ((RooRealVar*)rv_mu_ttbar_sig) ->getError() ;
//           qcd_sig_err = ((RooRealVar*)rv_mu_qcd_sig) ->getError() ;
//        } else {
//           ttbar_sig_val = ((RooFormulaVar*)rv_mu_ttbar_sig) ->getVal() ;
//           qcd_sig_val = ((RooFormulaVar*)rv_mu_qcd_sig) ->getVal() ;
//           ttbar_sb_val = ((RooRealVar*)rv_mu_ttbar_sb) ->getVal() ;
//           qcd_sb_val = ((RooRealVar*)rv_mu_qcd_sb) ->getVal() ;
//           ttbar_sb_err = ((RooRealVar*)rv_mu_ttbar_sb) ->getError() ;
//           qcd_sb_err = ((RooRealVar*)rv_mu_qcd_sb) ->getError() ;
//        }

//        float ltop = 0.90 ;
//        float lx = 0.78 ;
//        float dy = 0.06 ;
//        if ( useSigBgVars ) {
//           sprintf( fitlabel, "%5d", Nsb ) ;
//           fittext->DrawTextNDC( lx, ltop, fitlabel ) ;
//           sprintf( fitlabel, "%4.0f", ttbar_sb_val ) ;
//           fittext->DrawTextNDC( lx, ltop-dy, fitlabel ) ;
//           sprintf( fitlabel, "%4.0f", qcd_sb_val ) ;
//           fittext->DrawTextNDC( lx, ltop-2*dy, fitlabel ) ;
//           sprintf( fitlabel, "%4.0f", ew_sb ) ;
//           fittext->DrawTextNDC( lx, ltop-3*dy, fitlabel ) ;
//        } else {
//           sprintf( fitlabel, "%5d", Nsb ) ;
//           fittext->DrawTextNDC( lx, ltop, fitlabel ) ;
//           sprintf( fitlabel, "%4.0f +/- %4.0f", ttbar_sb_val, ttbar_sb_err ) ;
//           fittext->DrawTextNDC( lx, ltop-dy, fitlabel ) ;
//           sprintf( fitlabel, "%4.0f +/- %4.0f", qcd_sb_val, qcd_sb_err ) ;
//           fittext->DrawTextNDC( lx, ltop-2*dy, fitlabel ) ;
//           sprintf( fitlabel, "%4.0f", ew_sb ) ;
//           fittext->DrawTextNDC( lx, ltop-3*dy, fitlabel ) ;
//        }

//        char outfile[1000] ;
//        sprintf( outfile, "%s-sb.png", plotBaseName ) ;
//        can_sbfit->SaveAs( outfile ) ;

//        return true ;

//     } // sbPlotsUniformBins

////==================================================================================================================

//     bool ra2bRoostatsClass2::sbPlotsVariableBins( const char* plotBaseName ) {

//        if ( ! initialized ) {
//           printf("\n\n *** Call initialize first.\n\n") ;
//           return false ;
//        }

//        if ( ! varsAtFitVals ) {
//           printf("\n\n *** Try this right after calling doFit.\n\n") ;
//           return false ;
//        }

//     //--  Drawn with variable-width bins.

//        int nxbins = 5 ;
//        float xbinedges[6] = { 0., 160., 180., 260., 400., 800. } ;
//        float xbinwid[5] ;
//        for ( int bi=0; bi<nxbins; bi++ ) { xbinwid[bi] = xbinedges[bi+1] - xbinedges[bi] ; }

//        TH1F* hvbm3j_sb_data  = new TH1F("hvbm3j_sb_data" ,"3-jet mass, SB data" , nxbins, xbinedges ) ;
//        TH1F* hvbm3j_sb_ttbar = new TH1F("hvbm3j_sb_ttbar","3-jet mass, SB ttbar", nxbins, xbinedges ) ;
//        TH1F* hvbm3j_sb_qcd   = new TH1F("hvbm3j_sb_qcd"  ,"3-jet mass, SB qcd"  , nxbins, xbinedges ) ;
//        TH1F* hvbm3j_sb_ew    = new TH1F("hvbm3j_sb_ew"   ,"3-jet mass, SB ew"   , nxbins, xbinedges ) ;

//        hvbm3j_sb_data->SetLineWidth(2) ;
//        hvbm3j_sb_data->SetMarkerStyle(20) ;
//        hvbm3j_sb_ew->SetFillColor(49) ;
//        hvbm3j_sb_qcd->SetFillColor(46) ;
//        hvbm3j_sb_ttbar->SetFillColor(42) ;

//        hvbm3j_sb_data->SetBinContent( 1, (rv_Nsb1->getVal()) * (xbinwid[1]/xbinwid[0]) ) ;
//        hvbm3j_sb_data->SetBinContent( 2, (rv_Nsb2->getVal()) * (xbinwid[1]/xbinwid[1]) ) ;
//        hvbm3j_sb_data->SetBinContent( 3, (rv_Nsb3->getVal()) * (xbinwid[1]/xbinwid[2]) ) ;
//        hvbm3j_sb_data->SetBinContent( 4, (rv_Nsb4->getVal()) * (xbinwid[1]/xbinwid[3]) ) ;
//        hvbm3j_sb_data->SetBinContent( 5, (rv_Nsb5->getVal()) * (xbinwid[1]/xbinwid[4]) ) ;

//        hvbm3j_sb_data->SetBinError( 1, sqrt(rv_Nsb1->getVal()) * (xbinwid[1]/xbinwid[0]) ) ;
//        hvbm3j_sb_data->SetBinError( 2, sqrt(rv_Nsb2->getVal()) * (xbinwid[1]/xbinwid[1]) ) ;
//        hvbm3j_sb_data->SetBinError( 3, sqrt(rv_Nsb3->getVal()) * (xbinwid[1]/xbinwid[2]) ) ;
//        hvbm3j_sb_data->SetBinError( 4, sqrt(rv_Nsb4->getVal()) * (xbinwid[1]/xbinwid[3]) ) ;
//        hvbm3j_sb_data->SetBinError( 5, sqrt(rv_Nsb5->getVal()) * (xbinwid[1]/xbinwid[4]) ) ;

//        hvbm3j_sb_ttbar->SetBinContent( 1, (rv_mu_ttbar_sb1->getVal()) * (xbinwid[1]/xbinwid[0]) ) ;
//        hvbm3j_sb_ttbar->SetBinContent( 2, (rv_mu_ttbar_sb2->getVal()) * (xbinwid[1]/xbinwid[1]) ) ;
//        hvbm3j_sb_ttbar->SetBinContent( 3, (rv_mu_ttbar_sb3->getVal()) * (xbinwid[1]/xbinwid[2]) ) ;
//        hvbm3j_sb_ttbar->SetBinContent( 4, (rv_mu_ttbar_sb4->getVal()) * (xbinwid[1]/xbinwid[3]) ) ;
//        hvbm3j_sb_ttbar->SetBinContent( 5, (rv_mu_ttbar_sb5->getVal()) * (xbinwid[1]/xbinwid[4]) ) ;

//        hvbm3j_sb_qcd->SetBinContent( 1, (rv_mu_qcd_sb1->getVal()) * (xbinwid[1]/xbinwid[0]) ) ;
//        hvbm3j_sb_qcd->SetBinContent( 2, (rv_mu_qcd_sb2->getVal()) * (xbinwid[1]/xbinwid[1]) ) ;
//        hvbm3j_sb_qcd->SetBinContent( 3, (rv_mu_qcd_sb3->getVal()) * (xbinwid[1]/xbinwid[2]) ) ;
//        hvbm3j_sb_qcd->SetBinContent( 4, (rv_mu_qcd_sb4->getVal()) * (xbinwid[1]/xbinwid[3]) ) ;
//        hvbm3j_sb_qcd->SetBinContent( 5, (rv_mu_qcd_sb5->getVal()) * (xbinwid[1]/xbinwid[4]) ) ;

//        hvbm3j_sb_ew->SetBinContent( 1, (rv_mu_ew_sb1->getVal()) * (xbinwid[1]/xbinwid[0]) ) ;
//        hvbm3j_sb_ew->SetBinContent( 2, (rv_mu_ew_sb2->getVal()) * (xbinwid[1]/xbinwid[1]) ) ;
//        hvbm3j_sb_ew->SetBinContent( 3, (rv_mu_ew_sb3->getVal()) * (xbinwid[1]/xbinwid[2]) ) ;
//        hvbm3j_sb_ew->SetBinContent( 4, (rv_mu_ew_sb4->getVal()) * (xbinwid[1]/xbinwid[3]) ) ;
//        hvbm3j_sb_ew->SetBinContent( 5, (rv_mu_ew_sb5->getVal()) * (xbinwid[1]/xbinwid[4]) ) ;

//        THStack* hvbstack_m3j_sb_fit = new THStack( "hvbstack_m3j_sb_fit", "SB fit, 3-jet mass" ) ;
//        hvbstack_m3j_sb_fit->Add( hvbm3j_sb_ew ) ;
//        hvbstack_m3j_sb_fit->Add( hvbm3j_sb_qcd ) ;
//        hvbstack_m3j_sb_fit->Add( hvbm3j_sb_ttbar ) ;

//        TCanvas* can_vbsbfit = new TCanvas("can_vbsbfit", "SB 3-jet mass fit", 700, 500 ) ;

//        //hvbm3j_sb_data->SetMaximum( 1.4*(hvbm3j_sb_data->GetMaximum()) ) ;


//        hvbm3j_sb_data->Draw("histpe") ;
//        hvbstack_m3j_sb_fit->Draw( "same" ) ;
//        hvbm3j_sb_data->Draw("samehistpe") ;

//        TLegend* vbm3j_legend = new TLegend(0.62,0.7,0.97,0.95) ;
//        vbm3j_legend->SetFillColor( kWhite ) ;
//        vbm3j_legend->AddEntry( hvbm3j_sb_data, "Data" ) ;
//        vbm3j_legend->AddEntry( hvbm3j_sb_ttbar, "ttbar" ) ;
//        vbm3j_legend->AddEntry( hvbm3j_sb_qcd, "QCD" ) ;
//        vbm3j_legend->AddEntry( hvbm3j_sb_ew, "EW" ) ;
//        vbm3j_legend->Draw() ;

//        TText* fittext = new TText() ;
//        fittext->SetTextSize(0.04) ;
//        char fitlabel[1000] ;

//        float ew_sb(0.) ;
//        ew_sb += rv_mu_ew_sb1->getVal() ;
//        ew_sb += rv_mu_ew_sb2->getVal() ;
//        ew_sb += rv_mu_ew_sb3->getVal() ;
//        ew_sb += rv_mu_ew_sb4->getVal() ;
//        ew_sb += rv_mu_ew_sb5->getVal() ;
//        int Nsb = rv_Nsb1->getVal() +rv_Nsb2->getVal() +rv_Nsb3->getVal() +rv_Nsb4->getVal() +rv_Nsb5->getVal()  ;
//        float ltop = 0.90 ;
//        float lx = 0.78 ;
//        float dy = 0.06 ;

//        float ttbar_sb_val(0.) ;
//        float qcd_sb_val(0.) ;
//        float ttbar_sb_err(0.) ;
//        float qcd_sb_err(0.) ;

//        float ttbar_sig_val(0.) ;
//        float qcd_sig_val(0.) ;
//        float ttbar_sig_err(0.) ;
//        float qcd_sig_err(0.) ;

//        if ( useSigBgVars ) {
//           ttbar_sb_val = ((RooFormulaVar*)rv_mu_ttbar_sb) ->getVal() ;
//           qcd_sb_val = ((RooFormulaVar*)rv_mu_qcd_sb) ->getVal() ;
//           ttbar_sig_val = ((RooRealVar*)rv_mu_ttbar_sig) ->getVal() ;
//           qcd_sig_val = ((RooRealVar*)rv_mu_qcd_sig) ->getVal() ;
//           ttbar_sig_err = ((RooRealVar*)rv_mu_ttbar_sig) ->getError() ;
//           qcd_sig_err = ((RooRealVar*)rv_mu_qcd_sig) ->getError() ;
//        } else {
//           ttbar_sig_val = ((RooFormulaVar*)rv_mu_ttbar_sig) ->getVal() ;
//           qcd_sig_val = ((RooFormulaVar*)rv_mu_qcd_sig) ->getVal() ;
//           ttbar_sb_val = ((RooRealVar*)rv_mu_ttbar_sb) ->getVal() ;
//           qcd_sb_val = ((RooRealVar*)rv_mu_qcd_sb) ->getVal() ;
//           ttbar_sb_err = ((RooRealVar*)rv_mu_ttbar_sb) ->getError() ;
//           qcd_sb_err = ((RooRealVar*)rv_mu_qcd_sb) ->getError() ;
//        }

//        if ( useSigBgVars ) {
//           sprintf( fitlabel, "%5d", Nsb ) ;
//           fittext->DrawTextNDC( lx, ltop, fitlabel ) ;
//           sprintf( fitlabel, "%4.0f", ttbar_sb_val ) ;
//           fittext->DrawTextNDC( lx, ltop-dy, fitlabel ) ;
//           sprintf( fitlabel, "%4.0f", qcd_sb_val ) ;
//           fittext->DrawTextNDC( lx, ltop-2*dy, fitlabel ) ;
//           sprintf( fitlabel, "%4.0f", ew_sb ) ;
//           fittext->DrawTextNDC( lx, ltop-3*dy, fitlabel ) ;
//        } else {
//           sprintf( fitlabel, "%5d", Nsb ) ;
//           fittext->DrawTextNDC( lx, ltop, fitlabel ) ;
//           sprintf( fitlabel, "%4.0f +/- %4.0f", ttbar_sb_val, ttbar_sb_err ) ;
//           fittext->DrawTextNDC( lx, ltop-dy, fitlabel ) ;
//           sprintf( fitlabel, "%4.0f +/- %4.0f", qcd_sb_val, qcd_sb_err ) ;
//           fittext->DrawTextNDC( lx, ltop-2*dy, fitlabel ) ;
//           sprintf( fitlabel, "%4.0f", ew_sb ) ;
//           fittext->DrawTextNDC( lx, ltop-3*dy, fitlabel ) ;
//        }

//        char outfile[1000] ;
//        sprintf( outfile, "%s-sb-vb.png", plotBaseName ) ;
//        can_vbsbfit->SaveAs( outfile ) ;

//        return true ;

//     } // sbPlotsVariableBins

  //===================================================================

   ra2bRoostatsClass2::~ra2bRoostatsClass2() {

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

    bool ra2bRoostatsClass2::initialize( const char* infile ) {


       printf( "\n\n Opening input file : %s\n\n", infile ) ;

       FILE* infp ;
       if ( (infp=fopen( infile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
          return false ;
       }

       sprintf( initializeFile, "%s", infile ) ;

       int Nsig(0), Na(0), Nd(0) ; //-- data counts in signal region, A, and D.
       int Nsb(0) ; //-- data counts in SB.
       int Nslsig(0) ; //-- data counts in SL, SIG.
       int Nslsb(0) ; //-- data counts in SL, SB.


       float Nqcdmcsig(0.), Nqcdmcsb(0.), Nqcdmca(0.), Nqcdmcd(0.) ; //-- QCD MC counts in SIG, SB, A, and D.
       float Nqcdmcslsig(0.), Nqcdmcslsb(0.) ; //-- QCD MC counts in SL; SIG and SB.
       float Nttbarmcsig(0.), Nttbarmcsb(0.), Nttbarmca(0.), Nttbarmcd(0.) ; //-- ttbar MC counts in SIG, SB, A, and D.
       float Nttbarmcslsig(0.), Nttbarmcslsb(0.) ; //-- ttbar MC counts in SL; SIG and SB.

       int NWJmcsig, NWJmca, NWJmcd, NWJmcsb, NWJmcslsig, NWJmcslsb ;
       int NZnnmcsig, NZnnmca, NZnnmcd, NZnnmcsb, NZnnmcslsig, NZnnmcslsb ;
       int NEwomcsig, NEwomca, NEwomcd, NEwomcsb, NEwomcslsig, NEwomcslsb ;

       float Nsusymcsig(0.), Nsusymca(0.), Nsusymcd(0.) ; //-- SUSY MC counts in SIG, A, and D.
       float Nsusymcsb(0.) ; //-- SUSY MC counts in SB.
       float Nsusymcslsig(0.) ; //-- SUSY MC counts in SL, SIG.
       float Nsusymcslsb(0.) ; //-- SUSY MC counts in SL, SB.

       //--- read in description line.
       printf("\n\n") ;
       char c(0) ;
       while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
       printf("\n\n") ;


       char label[1000] ;

      //--- Inputs generated with gen_roostats_input.c
      //    The order here must be consistent with the order there!

       fscanf( infp, "%s %g", label, &EffScaleFactor ) ;       printf( "%s %g\n", label, EffScaleFactor ) ;         
       fscanf( infp, "%s %g", label, &EffScaleFactorErr ) ;    printf( "%s %g\n", label, EffScaleFactorErr ) ;         
       fscanf( infp, "%s %d", label, &Nsig ) ;                 printf( "%s %d\n", label, Nsig ) ;         
       fscanf( infp, "%s %d", label, &Na ) ;                   printf( "%s %d\n", label, Na ) ;           
       fscanf( infp, "%s %d", label, &Nd ) ;                   printf( "%s %d\n", label, Nd ) ;           
       fscanf( infp, "%s %d", label, &Nsb ) ;                 printf( "%s %d\n", label, Nsb ) ;         
       fscanf( infp, "%s %d", label, &Nslsig ) ;              printf( "%s %d\n", label, Nslsig ) ;      
       fscanf( infp, "%s %d", label, &Nslsb ) ;               printf( "%s %d\n", label, Nslsb ) ;       
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
       fscanf( infp, "%s %g", label, &Nttbarmcsig ) ;          printf( "%s %g\n", label, Nttbarmcsig ) ;  
       fscanf( infp, "%s %g", label, &Nttbarmcsb ) ;           printf( "%s %g\n", label, Nttbarmcsb ) ;   
       fscanf( infp, "%s %g", label, &Nttbarmca ) ;            printf( "%s %g\n", label, Nttbarmca ) ;    
       fscanf( infp, "%s %g", label, &Nttbarmcd ) ;            printf( "%s %g\n", label, Nttbarmcd ) ;    
       fscanf( infp, "%s %g", label, &Nttbarmcslsig ) ;        printf( "%s %g\n", label, Nttbarmcslsig ) ;      
       fscanf( infp, "%s %g", label, &Nttbarmcslsb ) ;        printf( "%s %g\n", label, Nttbarmcslsb ) ;      

       fscanf( infp, "%s %g", label, &lsf_WJmc ) ;        printf( "%s %g\n", label, lsf_WJmc ) ;      
       fscanf( infp, "%s %g", label, &lsf_WJmc_err ) ;        printf( "%s %g\n", label, lsf_WJmc_err ) ;      
       fscanf( infp, "%s %d", label, &NWJmcsig ) ;        printf( "%s %d\n", label, NWJmcsig ) ;      
       fscanf( infp, "%s %d", label, &NWJmca ) ;        printf( "%s %d\n", label, NWJmca ) ;      
       fscanf( infp, "%s %d", label, &NWJmcd ) ;        printf( "%s %d\n", label, NWJmcd ) ;      
       fscanf( infp, "%s %d", label, &NWJmcsb ) ;        printf( "%s %d\n", label, NWJmcsb ) ;      
       fscanf( infp, "%s %d", label, &NWJmcslsig ) ;        printf( "%s %d\n", label, NWJmcslsig ) ;      
       fscanf( infp, "%s %d", label, &NWJmcslsb ) ;        printf( "%s %d\n", label, NWJmcslsb ) ;      

       fscanf( infp, "%s %g", label, &lsf_Znnmc ) ;        printf( "%s %g\n", label, lsf_Znnmc ) ;      
       fscanf( infp, "%s %g", label, &lsf_Znnmc_err ) ;        printf( "%s %g\n", label, lsf_Znnmc_err ) ;      
       fscanf( infp, "%s %d", label, &NZnnmcsig ) ;        printf( "%s %d\n", label, NZnnmcsig ) ;      
       fscanf( infp, "%s %d", label, &NZnnmca ) ;        printf( "%s %d\n", label, NZnnmca ) ;      
       fscanf( infp, "%s %d", label, &NZnnmcd ) ;        printf( "%s %d\n", label, NZnnmcd ) ;      
       fscanf( infp, "%s %d", label, &NZnnmcsb ) ;        printf( "%s %d\n", label, NZnnmcsb ) ;      
       fscanf( infp, "%s %d", label, &NZnnmcslsig ) ;        printf( "%s %d\n", label, NZnnmcslsig ) ;      
       fscanf( infp, "%s %d", label, &NZnnmcslsb ) ;        printf( "%s %d\n", label, NZnnmcslsb ) ;      

       fscanf( infp, "%s %g", label, &lsf_Ewomc ) ;        printf( "%s %g\n", label, lsf_Ewomc ) ;      
       fscanf( infp, "%s %g", label, &lsf_Ewomc_err ) ;        printf( "%s %g\n", label, lsf_Ewomc_err ) ;      
       fscanf( infp, "%s %d", label, &NEwomcsig ) ;        printf( "%s %d\n", label, NEwomcsig ) ;      
       fscanf( infp, "%s %d", label, &NEwomca ) ;        printf( "%s %d\n", label, NEwomca ) ;      
       fscanf( infp, "%s %d", label, &NEwomcd ) ;        printf( "%s %d\n", label, NEwomcd ) ;      
       fscanf( infp, "%s %d", label, &NEwomcsb ) ;        printf( "%s %d\n", label, NEwomcsb ) ;      
       fscanf( infp, "%s %d", label, &NEwomcslsig ) ;        printf( "%s %d\n", label, NEwomcslsig ) ;      
       fscanf( infp, "%s %d", label, &NEwomcslsb ) ;        printf( "%s %d\n", label, NEwomcslsb ) ;      

       fscanf( infp, "%s %g", label, &Nsusymcsig ) ;           printf( "%s %g\n", label, Nsusymcsig ) ;   
       fscanf( infp, "%s %g", label, &Nsusymca ) ;             printf( "%s %g\n", label, Nsusymca ) ;     
       fscanf( infp, "%s %g", label, &Nsusymcd ) ;             printf( "%s %g\n", label, Nsusymcd ) ;     
       fscanf( infp, "%s %g", label, &Nsusymcsb ) ;           printf( "%s %g\n", label, Nsusymcsb ) ;   
       fscanf( infp, "%s %g", label, &Nsusymcslsig ) ;        printf( "%s %g\n", label, Nsusymcslsig ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslsb ) ;        printf( "%s %g\n", label, Nsusymcslsb ) ;


       printf("\n Done reading in %s\n\n", infile ) ;
       fclose( infp ) ;






       //--- Print out a nice summary of the inputs.


       float Newmcsig = lsf_WJmc * NWJmcsig
                      + lsf_Znnmc * NZnnmcsig
                      + lsf_Ewomc * NEwomcsig ;

       float Newmcsb = lsf_WJmc * NWJmcsb
                      + lsf_Znnmc * NZnnmcsb
                      + lsf_Ewomc * NEwomcsb ;

       float Newmca = lsf_WJmc * NWJmca
                      + lsf_Znnmc * NZnnmca
                      + lsf_Ewomc * NEwomca ;

       float Newmcd = lsf_WJmc * NWJmcd
                      + lsf_Znnmc * NZnnmcd
                      + lsf_Ewomc * NEwomcd ;

       float Newmcslsig = lsf_WJmc * NWJmcslsig
                      + lsf_Znnmc * NZnnmcslsig
                      + lsf_Ewomc * NEwomcslsig ;

       float Newmcslsb = lsf_WJmc * NWJmcslsb
                      + lsf_Znnmc * NZnnmcslsb
                      + lsf_Ewomc * NEwomcslsb ;



       float Nsmsig = Nttbarmcsig + Nqcdmcsig + Newmcsig ;

       float Nsmsb = Nttbarmcsb + Nqcdmcsb + Newmcsb ;

       float Nsma = Nttbarmca + Nqcdmca + Newmca ;
       float Nsmd = Nttbarmcd + Nqcdmcd + Newmcd ;

       float Nsmslsig = Nttbarmcslsig + Nqcdmcslsig + Newmcslsig ;

       float Nsmslsb = Nttbarmcslsb + Nqcdmcslsb + Newmcslsb ;



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
       printf("\n\n\n") ;



       printf(" --- Defining observables.\n" ) ;

      //-- data counts in signal region, A, and D, SB, SL SB, SL SIG.

      rv_Nsig = new RooRealVar( "Nsig", "Nsig", 0.5, 400. ) ;
      rv_Na = new RooRealVar( "Na", "Na", 0.5, 1000000. ) ;
      rv_Nd = new RooRealVar( "Nd", "Nd", 0.5, 1000000. ) ;
      rv_Nsb = new RooRealVar( "Nsb", "Nsb", 0.5, 1000000. ) ;
      rv_Nslsig = new RooRealVar( "Nslsig", "Nslsig", 0.5, 1000000. ) ;
      rv_Nslsb = new RooRealVar( "Nslsb", "Nslsb", 0.5, 1000000. ) ;

      rv_Nsig -> setVal( Nsig ) ;
      rv_Na -> setVal( Na ) ;
      rv_Nd -> setVal( Nd ) ;
      rv_Nsb -> setVal( Nsb ) ;
      rv_Nslsig -> setVal( Nslsig ) ;
      rv_Nslsb -> setVal( Nslsb ) ;





      //-- QCD MC counts

      rv_Nqcdmca   = new RooRealVar( "Nqcdmca", "Nqcdmca", 0.5, 1000000. ) ;
      rv_Nqcdmcd   = new RooRealVar( "Nqcdmcd", "Nqcdmcd", 0.5, 1000000. ) ;
      rv_Nqcdmcsig = new RooRealVar( "Nqcdmcsig", "Nqcdmcsig", 0.5, 1000000. ) ;
      rv_Nqcdmcsb  = new RooRealVar( "Nqcdmcsb", "Nqcdmcsb", 0.5, 1000000. ) ;

      rv_Nqcdmca   -> setVal( Nqcdmca ) ;
      rv_Nqcdmcd   -> setVal( Nqcdmcd ) ;
      rv_Nqcdmcsig -> setVal( Nqcdmcsig ) ;
      rv_Nqcdmcsb  -> setVal( Nqcdmcsb ) ;






      //-- EW MC counts

      rv_NWJmcsig   = new RooRealVar( "NWJmcsig"  , "NWJmcsig"  , 0.0, 100. ) ;
      rv_NWJmca     = new RooRealVar( "NWJmca"    , "NWJmca"    , 0.0, 100. ) ;
      rv_NWJmcd     = new RooRealVar( "NWJmcd"    , "NWJmcd"    , 0.0, 100. ) ;
      rv_NWJmcsb    = new RooRealVar( "NWJmcsb"   , "NWJmcsb"   , 0.0, 100. ) ;
      rv_NWJmcslsig = new RooRealVar( "NWJmcslsig", "NWJmcslsig", 0.0, 100. ) ;
      rv_NWJmcslsb  = new RooRealVar( "NWJmcslsb" , "NWJmcslsb" , 0.0, 100. ) ;

      rv_NZnnmcsig   = new RooRealVar( "NZnnmcsig"  , "NZnnmcsig"  , 0.0, 100. ) ;
      rv_NZnnmca     = new RooRealVar( "NZnnmca"    , "NZnnmca"    , 0.0, 100. ) ;
      rv_NZnnmcd     = new RooRealVar( "NZnnmcd"    , "NZnnmcd"    , 0.0, 100. ) ;
      rv_NZnnmcsb    = new RooRealVar( "NZnnmcsb"   , "NZnnmcsb"   , 0.0, 100. ) ;
      rv_NZnnmcslsig = new RooRealVar( "NZnnmcslsig", "NZnnmcslsig", 0.0, 100. ) ;
      rv_NZnnmcslsb  = new RooRealVar( "NZnnmcslsb" , "NZnnmcslsb" , 0.0, 100. ) ;

      rv_NEwomcsig   = new RooRealVar( "NEwomcsig"  , "NEwomcsig"  , 0.0, 100. ) ;
      rv_NEwomca     = new RooRealVar( "NEwomca"    , "NEwomca"    , 0.0, 100. ) ;
      rv_NEwomcd     = new RooRealVar( "NEwomcd"    , "NEwomcd"    , 0.0, 100. ) ;
      rv_NEwomcsb    = new RooRealVar( "NEwomcsb"   , "NEwomcsb"   , 0.0, 100. ) ;
      rv_NEwomcslsig = new RooRealVar( "NEwomcslsig", "NEwomcslsig", 0.0, 100. ) ;
      rv_NEwomcslsb  = new RooRealVar( "NEwomcslsb" , "NEwomcslsb" , 0.0, 100. ) ;


      rv_NWJmcsig   -> setVal( NWJmcsig ) ;
      rv_NWJmca     -> setVal( NWJmca ) ;
      rv_NWJmcd     -> setVal( NWJmcd ) ;
      rv_NWJmcsb    -> setVal( NWJmcsb ) ;
      rv_NWJmcslsig -> setVal( NWJmcslsig ) ;
      rv_NWJmcslsb  -> setVal( NWJmcslsb ) ;

      rv_NZnnmcsig   -> setVal( NZnnmcsig ) ;
      rv_NZnnmca     -> setVal( NZnnmca ) ;
      rv_NZnnmcd     -> setVal( NZnnmcd ) ;
      rv_NZnnmcsb    -> setVal( NZnnmcsb ) ;
      rv_NZnnmcslsig -> setVal( NZnnmcslsig ) ;
      rv_NZnnmcslsb  -> setVal( NZnnmcslsb ) ;

      rv_NEwomcsig   -> setVal( NEwomcsig ) ;
      rv_NEwomca     -> setVal( NEwomca ) ;
      rv_NEwomcd     -> setVal( NEwomcd ) ;
      rv_NEwomcsb    -> setVal( NEwomcsb ) ;
      rv_NEwomcslsig -> setVal( NEwomcslsig ) ;
      rv_NEwomcslsb  -> setVal( NEwomcslsb ) ;











      //++++++++ Parameters of the likelihood

      printf(" --- Defining parameters.\n" ) ;




     //-- Counts in SIG

      if ( useSigBgVars ) {
         rrv_mu_ttbar_sig   = new RooRealVar( "mu_ttbar_sig", "mu_ttbar_sig", 0.0, 100. ) ;
         rv_mu_ttbar_sig = rrv_mu_ttbar_sig ;
      }
      rv_mu_qcd_sig     = new RooRealVar( "mu_qcd_sig"  , "mu_qcd_sig"  , 0.0,  50. ) ;
      //-- ew is rfv
      rv_mu_susy_sig    = new RooRealVar( "mu_susy_sig" , "mu_susy_sig" , 0.0, 150. ) ;
      rv_mu_susymc_sig  = new RooRealVar( "mu_susymc_sig" , "mu_susymc_sig" , 0.0, 100000. ) ;


      if ( useSigBgVars ) {
         rrv_mu_ttbar_sig   -> setVal( Nttbarmcsig ) ;  //-- this is a starting value only.
      }
      rv_mu_qcd_sig     -> setVal( Nqcdmcsig ) ;  //-- this is a starting value only.
      rv_mu_susy_sig    -> setVal( 0. ) ;  //-- this is a starting value only.
      rv_mu_susymc_sig  -> setVal( 0.1 ) ;

      rv_mu_susymc_sig -> setConstant(kTRUE) ;





     //-- Counts in SB

      if ( !useSigBgVars ) {
         rrv_mu_ttbar_sb   = new RooRealVar( "mu_ttbar_sb", "mu_ttbar_sb", 0.0, 500. ) ;
         rv_mu_ttbar_sb = rrv_mu_ttbar_sb ;
      }
      rv_mu_qcd_sb    = new RooRealVar( "mu_qcd_sb"  , "mu_qcd_sb"  , 0.0, 200. ) ;
      //-- ew is rfv
      rv_mu_susymc_sb = new RooRealVar( "mu_susymc_sb", "mu_susymc_sb", 0.0, 10000. ) ;
      //-- susy is rfv


      if ( !useSigBgVars ) {
         rrv_mu_ttbar_sb   -> setVal( Nttbarmcsb ) ;  //-- this is a starting value only.
      }
      rv_mu_qcd_sb    -> setVal( Nqcdmcsb ) ;  //-- this is a starting value only.
      rv_mu_susymc_sb -> setVal( 0. ) ;

      rv_mu_susymc_sb -> setConstant(kTRUE) ;






     //-- Counts in A, signal selection

      rv_mu_ttbar_a   = new RooRealVar( "mu_ttbar_a", "mu_ttbar_a", 0.0, 10000. ) ;
      //-- qcd is rfv
      //-- ew is rfv
      rv_mu_susymc_a  = new RooRealVar( "mu_susymc_a" , "mu_susymc_a" , 0.0, 10000. ) ;
      //-- susy is rfv

      rv_mu_ttbar_a   -> setVal( Nttbarmca ) ;
      rv_mu_susymc_a  -> setVal( 0. ) ;

      rv_mu_susymc_a  ->setConstant(kTRUE) ;
      rv_mu_ttbar_a   ->setConstant(kTRUE) ;




     //-- Counts in D, signal selection

      rv_mu_ttbar_d   = new RooRealVar( "mu_ttbar_d", "mu_ttbar_d", 0.0, 10000. ) ;
      //-- qcd is rfv
      //-- ew is rfv
      rv_mu_susymc_d  = new RooRealVar( "mu_susymc_d" , "mu_susymc_d" , 0.0, 10000. ) ;
      //-- susy is rfv

      rv_mu_ttbar_d   -> setVal( Nttbarmcd ) ;
      rv_mu_susymc_d  -> setVal( 0. ) ;

      rv_mu_ttbar_d   -> setConstant(kTRUE) ;
      rv_mu_susymc_d  -> setConstant(kTRUE) ;




     //-- Single Lepton (SL) counts, SIG region

      rv_mu_sl_ttbar_sig  = new RooRealVar( "mu_sl_ttbar_sig", "mu_sl_ttbar_sig", 1.0, 250. ) ;
      //-- ignore QCD in sl
      //-- ew is rfv
      rv_mu_sl_susymc_sig = new RooRealVar( "mu_sl_susymc_sig", "mu_sl_susymc_sig", 0.0, 10000. ) ;
      //-- susy is rfv

      rv_mu_sl_ttbar_sig  -> setVal( Nslsig - Newmcslsig ) ;
      rv_mu_sl_susymc_sig -> setVal( 0. ) ;

      rv_mu_sl_susymc_sig -> setConstant( kTRUE ) ;





     //-- Single Lepton (SL) counts, SB region

      rv_mu_sl_ttbar_sb  = new RooRealVar( "mu_sl_ttbar_sb", "mu_sl_ttbar_sb", 0.0, 250. ) ;
      //-- ignore QCD in sl
      //-- ew is rfv
      rv_mu_sl_susymc_sb = new RooRealVar( "mu_sl_susymc_sb", "mu_sl_susymc_sb", 0.0, 10000. ) ;
      //-- susy is rfv

      rv_mu_sl_ttbar_sb  -> setVal( Nslsb - Newmcslsb ) ;
      rv_mu_sl_susymc_sb -> setVal( 0. ) ;

      rv_mu_sl_susymc_sb -> setConstant( kTRUE ) ;






     //-- QCD MC, SIG,SB,A,D counts, signal selection

      rv_mu_qcdmc_sig = new RooRealVar( "mu_qcdmc_sig", "mu_qcdmc_sig", 0.1, 10000. ) ;
      rv_mu_qcdmc_sb  = new RooRealVar( "mu_qcdmc_sb" , "mu_qcdmc_sb" , 0.1, 10000. ) ;
      rv_mu_qcdmc_a   = new RooRealVar( "mu_qcdmc_a"  , "mu_qcdmc_a"  , 0.0, 10000. ) ;
      rv_mu_qcdmc_d   = new RooRealVar( "mu_qcdmc_d"  , "mu_qcdmc_d"  , 0.0, 10000. ) ;

      rv_mu_qcdmc_sig -> setVal( Nqcdmcsig ) ;
      rv_mu_qcdmc_sb  -> setVal( Nqcdmcsb ) ;
      rv_mu_qcdmc_a   -> setVal( Nqcdmca ) ;
      rv_mu_qcdmc_d   -> setVal( Nqcdmcd ) ;

 //++++++++ try fixing (!)

///// rv_mu_qcdmc_sig -> setConstant( kTRUE ) ;
///// rv_mu_qcdmc_sb  -> setConstant( kTRUE ) ;
///// rv_mu_qcdmc_a   -> setConstant( kTRUE ) ;
///// rv_mu_qcdmc_d   -> setConstant( kTRUE ) ;





     //-- EW MC counts

      rv_mu_wjmc_sig   = new RooRealVar( "mu_wjmc_sig"  , "mu_wjmc_sig"  , 0.00, (3+NWJmcsig+5*sqrt(1.*NWJmcsig)) ) ;
      rv_mu_wjmc_a     = new RooRealVar( "mu_wjmc_a"    , "mu_wjmc_a"    , 0.00, 50. ) ;
      rv_mu_wjmc_d     = new RooRealVar( "mu_wjmc_d"    , "mu_wjmc_d"    , 0.00, 50. ) ;
      rv_mu_wjmc_sb    = new RooRealVar( "mu_wjmc_sb"   , "mu_wjmc_sb"   , 0.00, 50. ) ;
      rv_mu_wjmc_slsig = new RooRealVar( "mu_wjmc_slsig", "mu_wjmc_slsig", 0.00, 50. ) ;
      rv_mu_wjmc_slsb  = new RooRealVar( "mu_wjmc_slsb" , "mu_wjmc_slsb" , 0.00, 50. ) ;

      rv_mu_znnmc_sig   = new RooRealVar( "mu_znnmc_sig"  , "mu_znnmc_sig"  , 0.00, (3+NZnnmcsig+5*sqrt(1.*NZnnmcsig)) ) ;
      rv_mu_znnmc_a     = new RooRealVar( "mu_znnmc_a"    , "mu_znnmc_a"    , 0.00, 50. ) ;
      rv_mu_znnmc_d     = new RooRealVar( "mu_znnmc_d"    , "mu_znnmc_d"    , 0.00, 50. ) ;
      rv_mu_znnmc_sb    = new RooRealVar( "mu_znnmc_sb"   , "mu_znnmc_sb"   , 0.00, 50. ) ;
      rv_mu_znnmc_slsig = new RooRealVar( "mu_znnmc_slsig", "mu_znnmc_slsig", 0.00, 50. ) ;
      rv_mu_znnmc_slsb  = new RooRealVar( "mu_znnmc_slsb" , "mu_znnmc_slsb" , 0.00, 50. ) ;

      rv_mu_ewomc_sig   = new RooRealVar( "mu_ewomc_sig"  , "mu_ewomc_sig"  , 0.00, 50. ) ;
      rv_mu_ewomc_a     = new RooRealVar( "mu_ewomc_a"    , "mu_ewomc_a"    , 0.00, 50. ) ;
      rv_mu_ewomc_d     = new RooRealVar( "mu_ewomc_d"    , "mu_ewomc_d"    , 0.00, 50. ) ;
      rv_mu_ewomc_sb    = new RooRealVar( "mu_ewomc_sb"   , "mu_ewomc_sb"   , 0.00, 50. ) ;
      rv_mu_ewomc_slsig = new RooRealVar( "mu_ewomc_slsig", "mu_ewomc_slsig", 0.00, 50. ) ;
      rv_mu_ewomc_slsb  = new RooRealVar( "mu_ewomc_slsb" , "mu_ewomc_slsb" , 0.00, 50. ) ;


      rv_mu_wjmc_sig   -> setVal( NWJmcsig ) ;
      rv_mu_wjmc_a     -> setVal( NWJmca ) ;
      rv_mu_wjmc_d     -> setVal( NWJmcd ) ;
      rv_mu_wjmc_sb    -> setVal( NWJmcsb ) ;
      rv_mu_wjmc_slsig -> setVal( NWJmcslsig ) ;
      rv_mu_wjmc_slsb  -> setVal( NWJmcslsb ) ;

      rv_mu_znnmc_sig   -> setVal( NZnnmcsig ) ;
      rv_mu_znnmc_a     -> setVal( NZnnmca ) ;
      rv_mu_znnmc_d     -> setVal( NZnnmcd ) ;
      rv_mu_znnmc_sb    -> setVal( NZnnmcsb ) ;
      rv_mu_znnmc_slsig -> setVal( NZnnmcslsig ) ;
      rv_mu_znnmc_slsb  -> setVal( NZnnmcslsb ) ;

      rv_mu_ewomc_sig   -> setVal( NEwomcsig ) ;
      rv_mu_ewomc_a     -> setVal( NEwomca ) ;
      rv_mu_ewomc_d     -> setVal( NEwomcd ) ;
      rv_mu_ewomc_sb    -> setVal( NEwomcsb ) ;
      rv_mu_ewomc_slsig -> setVal( NEwomcslsig ) ;
      rv_mu_ewomc_slsb  -> setVal( NEwomcslsb ) ;


   //+++ Dont use EW other for now.
      rv_mu_ewomc_sig   -> setConstant( kTRUE ) ;
      rv_mu_ewomc_a     -> setConstant( kTRUE ) ;
      rv_mu_ewomc_d     -> setConstant( kTRUE ) ;
      rv_mu_ewomc_sb    -> setConstant( kTRUE ) ;
      rv_mu_ewomc_slsig -> setConstant( kTRUE ) ;
      rv_mu_ewomc_slsb  -> setConstant( kTRUE ) ;

      rv_mu_znnmc_slsig -> setConstant( kTRUE ) ;
      rv_mu_znnmc_slsb  -> setConstant( kTRUE ) ;

   //++++++ fix all EW for now.
////  rv_mu_wjmc_sig   -> setConstant( kTRUE ) ;
////  rv_mu_wjmc_a     -> setConstant( kTRUE ) ;
////  rv_mu_wjmc_d     -> setConstant( kTRUE ) ;
////  rv_mu_wjmc_sb    -> setConstant( kTRUE ) ;
////  rv_mu_wjmc_slsig -> setConstant( kTRUE ) ;
////  rv_mu_wjmc_slsb  -> setConstant( kTRUE ) ;

////  rv_mu_znnmc_sig   -> setConstant( kTRUE ) ;
////  rv_mu_znnmc_a     -> setConstant( kTRUE ) ;
////  rv_mu_znnmc_d     -> setConstant( kTRUE ) ;
////  rv_mu_znnmc_sb    -> setConstant( kTRUE ) ;











    //-- Efficiency scale factor.

      double pmin, pmax ;

      pmin = (EffScaleFactor-4.*EffScaleFactorErr) ;
      pmax = (EffScaleFactor+4.*EffScaleFactorErr) ;
      if ( pmin < 0 ) pmin = 0.1 ;
      rv_eff_sf  = new RooRealVar( "eff_sf"      , "eff_sf"      , pmin, pmax ) ;
      rv_eff_sf  -> setVal( EffScaleFactor ) ;


      pmin = (lsf_WJmc-4.*lsf_WJmc_err) ;
      pmax = (lsf_WJmc+4.*lsf_WJmc_err) ;
      if ( pmin < 0 ) pmin = 0.1 ;
      rv_lsf_wjmc  = new RooRealVar( "lsf_wjmc", "lsf_wjmc", pmin, pmax ) ;
      rv_lsf_wjmc  -> setVal( lsf_WJmc ) ;

      pmin = (lsf_Znnmc-4.*lsf_Znnmc_err) ;
      pmax = (lsf_Znnmc+4.*lsf_Znnmc_err) ;
      if ( pmin < 0 ) pmin = 0.1 ;
      rv_lsf_znnmc  = new RooRealVar( "lsf_znnmc", "lsf_znnmc", pmin, pmax ) ;
      rv_lsf_znnmc  -> setVal( lsf_Znnmc ) ;

      pmin = (lsf_Ewomc-4.*lsf_Ewomc_err) ;
      pmax = (lsf_Ewomc+4.*lsf_Ewomc_err) ;
      if ( pmin < 0 ) pmin = 0.1 ;
      rv_lsf_ewomc  = new RooRealVar( "lsf_ewomc", "lsf_ewomc", pmin, pmax ) ;
      rv_lsf_ewomc  -> setVal( lsf_Ewomc ) ;

   //+++ Dont use EW other for now.
      rv_lsf_ewomc -> setConstant( kTRUE ) ;

/////++++++ try fixed lsf.
///   rv_lsf_wjmc  -> setConstant( kTRUE ) ;
///   rv_lsf_znnmc -> setConstant( kTRUE ) ;

/////+++++++++ try fixed eff_sf
///   rv_eff_sf    -> setConstant( kTRUE ) ;





     //+++++++++++++++++ Relationships between parameters

       printf(" --- Defining relationships between parameters.\n" ) ;


      if ( useSigBgVars ) {
         rfv_mu_ttbar_sb = new RooFormulaVar( "mu_ttbar_sb",
                                                       "mu_ttbar_sig*(mu_sl_ttbar_sb/mu_sl_ttbar_sig)",
                                                       RooArgSet( *rv_mu_ttbar_sig, *rv_mu_sl_ttbar_sb, *rv_mu_sl_ttbar_sig ) ) ;
         rv_mu_ttbar_sb = rfv_mu_ttbar_sb ;
      } else {
         rfv_mu_ttbar_sig = new RooFormulaVar( "mu_ttbar_sig",
                                                       "mu_ttbar_sb*(mu_sl_ttbar_sig/mu_sl_ttbar_sb)",
                                                       RooArgSet( *rv_mu_ttbar_sb, *rv_mu_sl_ttbar_sig, *rv_mu_sl_ttbar_sb ) ) ;
         rv_mu_ttbar_sig = rfv_mu_ttbar_sig ;
      }




      rv_mu_qcd_a = new RooFormulaVar( "mu_qcd_a",
                                         "mu_qcd_sb*(mu_qcdmc_a/mu_qcdmc_sb)",
                                         RooArgSet( *rv_mu_qcd_sb, *rv_mu_qcdmc_a, *rv_mu_qcdmc_sb ) ) ;


      rv_mu_qcd_d = new RooFormulaVar( "mu_qcd_d",
                                         "mu_qcd_sig*(mu_qcdmc_d/mu_qcdmc_sig)",
                                         RooArgSet( *rv_mu_qcd_sig, *rv_mu_qcdmc_d, *rv_mu_qcdmc_sig ) ) ;







      rv_mu_susy_a = new RooFormulaVar( "mu_susy_a",
                                        "mu_susymc_a * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_a, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_d = new RooFormulaVar( "mu_susy_d",
                                        "mu_susymc_d * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_d, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;


      rv_mu_susy_sb = new RooFormulaVar( "mu_susy_sb",
                                        "mu_susymc_sb * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sb, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;



      rv_mu_sl_susy_sb = new RooFormulaVar( "mu_sl_susy_sb",
                                        "mu_sl_susymc_sb * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_sl_susymc_sb, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;


      rv_mu_sl_susy_sig = new RooFormulaVar( "mu_sl_susy_sig",
                                        "mu_sl_susymc_sig * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_sl_susymc_sig, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;






      rv_mu_ew_sig = new RooFormulaVar( "mu_ew_sig",
                                        "lsf_wjmc*mu_wjmc_sig  + lsf_znnmc*mu_znnmc_sig + lsf_ewomc*mu_ewomc_sig",
                                        RooArgSet( *rv_lsf_wjmc, *rv_mu_wjmc_sig,
                                                   *rv_lsf_znnmc, *rv_mu_znnmc_sig,
                                                   *rv_lsf_ewomc, *rv_mu_ewomc_sig ) ) ;

      rv_mu_ew_a = new RooFormulaVar( "mu_ew_a",
                                        "lsf_wjmc*mu_wjmc_a  + lsf_znnmc*mu_znnmc_a + lsf_ewomc*mu_ewomc_a",
                                        RooArgSet( *rv_lsf_wjmc, *rv_mu_wjmc_a,
                                                   *rv_lsf_znnmc, *rv_mu_znnmc_a,
                                                   *rv_lsf_ewomc, *rv_mu_ewomc_a ) ) ;

      rv_mu_ew_d = new RooFormulaVar( "mu_ew_d",
                                        "lsf_wjmc*mu_wjmc_d  + lsf_znnmc*mu_znnmc_d + lsf_ewomc*mu_ewomc_d",
                                        RooArgSet( *rv_lsf_wjmc, *rv_mu_wjmc_d,
                                                   *rv_lsf_znnmc, *rv_mu_znnmc_d,
                                                   *rv_lsf_ewomc, *rv_mu_ewomc_d ) ) ;

      rv_mu_ew_sb = new RooFormulaVar( "mu_ew_sb",
                                        "lsf_wjmc*mu_wjmc_sb  + lsf_znnmc*mu_znnmc_sb + lsf_ewomc*mu_ewomc_sb",
                                        RooArgSet( *rv_lsf_wjmc, *rv_mu_wjmc_sb,
                                                   *rv_lsf_znnmc, *rv_mu_znnmc_sb,
                                                   *rv_lsf_ewomc, *rv_mu_ewomc_sb ) ) ;

      rv_mu_sl_ew_sig = new RooFormulaVar( "mu_sl_ew_sig",
                                        "lsf_wjmc*mu_wjmc_slsig  + lsf_znnmc*mu_znnmc_slsig + lsf_ewomc*mu_ewomc_slsig",
                                        RooArgSet( *rv_lsf_wjmc, *rv_mu_wjmc_slsig,
                                                   *rv_lsf_znnmc, *rv_mu_znnmc_slsig,
                                                   *rv_lsf_ewomc, *rv_mu_ewomc_slsig ) ) ;

      rv_mu_sl_ew_sb = new RooFormulaVar( "mu_sl_ew_sb",
                                        "lsf_wjmc*mu_wjmc_slsb  + lsf_znnmc*mu_znnmc_slsb + lsf_ewomc*mu_ewomc_slsb",
                                        RooArgSet( *rv_lsf_wjmc, *rv_mu_wjmc_slsb,
                                                   *rv_lsf_znnmc, *rv_mu_znnmc_slsb,
                                                   *rv_lsf_ewomc, *rv_mu_ewomc_slsb ) ) ;



    //+++++++++++++ Expected counts for observables in terms of parameters ++++++++++++++++++

       printf(" --- Defining expected counts in terms of parameters.\n" ) ;



      rv_n_sig = new RooFormulaVar( "n_sig",
                                     "mu_ttbar_sig + mu_qcd_sig + mu_ew_sig + eff_sf*mu_susy_sig",
                                     RooArgSet( *rv_mu_ttbar_sig, *rv_mu_qcd_sig,
                                     *rv_mu_ew_sig, *rv_eff_sf, *rv_mu_susy_sig ) ) ;

      rv_n_sb = new RooFormulaVar( "n_sb",
                                     "mu_ttbar_sb + mu_qcd_sb + mu_ew_sb + eff_sf*mu_susy_sb",
                                     RooArgSet( *rv_mu_ttbar_sb, *rv_mu_qcd_sb,
                                     *rv_mu_ew_sb, *rv_eff_sf, *rv_mu_susy_sb ) ) ;

      rv_n_a = new RooFormulaVar( "n_a",
                                   "mu_qcd_a + mu_ttbar_a + mu_ew_a + eff_sf*mu_susy_a",
                                   RooArgSet( *rv_mu_qcd_a, *rv_mu_ttbar_a,
                                   *rv_mu_ew_a, *rv_eff_sf, *rv_mu_susy_a ) ) ;

      rv_n_d = new RooFormulaVar( "n_d",
                                   "mu_qcd_d + mu_ttbar_d + mu_ew_d + eff_sf*mu_susy_d",
                                   RooArgSet( *rv_mu_qcd_d, *rv_mu_ttbar_d,
                                   *rv_mu_ew_d, *rv_eff_sf, *rv_mu_susy_d ) ) ;

      rv_n_sl_sig = new RooFormulaVar( "n_sl_sig",
                                         "mu_sl_ttbar_sig + mu_sl_ew_sig + eff_sf*mu_sl_susy_sig",
                                         RooArgSet( *rv_mu_sl_ttbar_sig, *rv_mu_sl_ew_sig, *rv_eff_sf, *rv_mu_sl_susy_sig ) ) ;

      rv_n_sl_sb = new RooFormulaVar( "n_sl_sb",
                                         "mu_sl_ttbar_sb + mu_sl_ew_sb + eff_sf*mu_sl_susy_sb",
                                         RooArgSet( *rv_mu_sl_ttbar_sb, *rv_mu_sl_ew_sb, *rv_eff_sf, *rv_mu_sl_susy_sb ) ) ;



      rv_n_wjmc_sig    = new RooFormulaVar( "n_wjmc_sig"   , "mu_wjmc_sig"  , RooArgSet( *rv_mu_wjmc_sig ) ) ;
      rv_n_wjmc_a      = new RooFormulaVar( "n_wjmc_a"     , "mu_wjmc_a"    , RooArgSet( *rv_mu_wjmc_a ) ) ;
      rv_n_wjmc_d      = new RooFormulaVar( "n_wjmc_d"     , "mu_wjmc_d"    , RooArgSet( *rv_mu_wjmc_d ) ) ;
      rv_n_wjmc_sb     = new RooFormulaVar( "n_wjmc_sb"    , "mu_wjmc_sb"   , RooArgSet( *rv_mu_wjmc_sb ) ) ;
      rv_n_wjmc_sl_sig = new RooFormulaVar( "n_wjmc_sl_sig", "mu_wjmc_slsig", RooArgSet( *rv_mu_wjmc_slsig ) ) ;
      rv_n_wjmc_sl_sb  = new RooFormulaVar( "n_wjmc_sl_sb" , "mu_wjmc_slsb" , RooArgSet( *rv_mu_wjmc_slsb ) ) ;

      rv_n_znnmc_sig    = new RooFormulaVar( "n_znnmc_sig"   , "mu_znnmc_sig"  , RooArgSet( *rv_mu_znnmc_sig ) ) ;
      rv_n_znnmc_a      = new RooFormulaVar( "n_znnmc_a"     , "mu_znnmc_a"    , RooArgSet( *rv_mu_znnmc_a ) ) ;
      rv_n_znnmc_d      = new RooFormulaVar( "n_znnmc_d"     , "mu_znnmc_d"    , RooArgSet( *rv_mu_znnmc_d ) ) ;
      rv_n_znnmc_sb     = new RooFormulaVar( "n_znnmc_sb"    , "mu_znnmc_sb"   , RooArgSet( *rv_mu_znnmc_sb ) ) ;
      rv_n_znnmc_sl_sig = new RooFormulaVar( "n_znnmc_sl_sig", "mu_znnmc_slsig", RooArgSet( *rv_mu_znnmc_slsig ) ) ;
      rv_n_znnmc_sl_sb  = new RooFormulaVar( "n_znnmc_sl_sb" , "mu_znnmc_slsb" , RooArgSet( *rv_mu_znnmc_slsb ) ) ;

      rv_n_ewomc_sig    = new RooFormulaVar( "n_ewomc_sig"   , "mu_ewomc_sig"  , RooArgSet( *rv_mu_ewomc_sig ) ) ;
      rv_n_ewomc_a      = new RooFormulaVar( "n_ewomc_a"     , "mu_ewomc_a"    , RooArgSet( *rv_mu_ewomc_a ) ) ;
      rv_n_ewomc_d      = new RooFormulaVar( "n_ewomc_d"     , "mu_ewomc_d"    , RooArgSet( *rv_mu_ewomc_d ) ) ;
      rv_n_ewomc_sb     = new RooFormulaVar( "n_ewomc_sb"    , "mu_ewomc_sb"   , RooArgSet( *rv_mu_ewomc_sb ) ) ;
      rv_n_ewomc_sl_sig = new RooFormulaVar( "n_ewomc_sl_sig", "mu_ewomc_slsig", RooArgSet( *rv_mu_ewomc_slsig ) ) ;
      rv_n_ewomc_sl_sb  = new RooFormulaVar( "n_ewomc_sl_sb" , "mu_ewomc_slsb" , RooArgSet( *rv_mu_ewomc_slsb ) ) ;













   //++++++++++++ PDFs for the likelihood +++++++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining PDFs of the likelihood.\n" ) ;

      pdf_Nsig = new RooPoisson( "pdf_Nsig", "Nsig Poisson PDF", *rv_Nsig, *rv_n_sig ) ;
      pdf_Na = new RooPoisson( "pdf_Na", "Na Poisson PDF", *rv_Na, *rv_n_a ) ;
      pdf_Nd = new RooPoisson( "pdf_Nd", "Nd Poisson PDF", *rv_Nd, *rv_n_d ) ;
      pdf_Nsb = new RooPoisson( "pdf_Nsb", "Nsb Poisson PDF", *rv_Nsb, *rv_n_sb ) ;
      pdf_Nsl_sig = new RooPoisson( "pdf_Nsl_sig", "Nsl,sig Poisson PDF", *rv_Nslsig, *rv_n_sl_sig ) ;
      pdf_Nsl_sb = new RooPoisson( "pdf_Nsl_sb", "Nsl,sb Poisson PDF", *rv_Nslsb, *rv_n_sl_sb ) ;

      pdf_NWJmcsig = new RooPoisson( "pdf_NWJmcsig", "Nsig Poisson PDF", *rv_NWJmcsig, *rv_n_wjmc_sig ) ;
      pdf_NWJmca = new RooPoisson( "pdf_NWJmca", "Na Poisson PDF", *rv_NWJmca, *rv_n_wjmc_a ) ;
      pdf_NWJmcd = new RooPoisson( "pdf_NWJmcd", "Nd Poisson PDF", *rv_NWJmcd, *rv_n_wjmc_d ) ;
      pdf_NWJmcsb = new RooPoisson( "pdf_NWJmcsb", "Nsb Poisson PDF", *rv_NWJmcsb, *rv_n_wjmc_sb ) ;
      pdf_NWJmcsl_sig = new RooPoisson( "pdf_NWJmcsl_sig", "Nsl,sig Poisson PDF", *rv_NWJmcslsig, *rv_n_wjmc_sl_sig ) ;
      pdf_NWJmcsl_sb = new RooPoisson( "pdf_NWJmcsl_sb", "Nsl,sb Poisson PDF", *rv_NWJmcslsb, *rv_n_wjmc_sl_sb ) ;

      pdf_NZnnmcsig = new RooPoisson( "pdf_NZnnmcsig", "Nsig Poisson PDF", *rv_NZnnmcsig, *rv_n_znnmc_sig ) ;
      pdf_NZnnmca = new RooPoisson( "pdf_NZnnmca", "Na Poisson PDF", *rv_NZnnmca, *rv_n_znnmc_a ) ;
      pdf_NZnnmcd = new RooPoisson( "pdf_NZnnmcd", "Nd Poisson PDF", *rv_NZnnmcd, *rv_n_znnmc_d ) ;
      pdf_NZnnmcsb = new RooPoisson( "pdf_NZnnmcsb", "Nsb Poisson PDF", *rv_NZnnmcsb, *rv_n_znnmc_sb ) ;
 ///  pdf_NZnnmcsl_sig = new RooPoisson( "pdf_NZnnmcsl_sig", "Nsl,sig Poisson PDF", *rv_NZnnmcslsig, *rv_n_znnmc_sl_sig ) ;
 ///  pdf_NZnnmcsl_sb = new RooPoisson( "pdf_NZnnmcsl_sb", "Nsl,sb Poisson PDF", *rv_NZnnmcslsb, *rv_n_znnmc_sl_sb ) ;

 ///  pdf_NEwomcsig = new RooPoisson( "pdf_NEwomcsig", "Nsig Poisson PDF", *rv_NEwomcsig, *rv_n_ewomc_sig ) ;
 ///  pdf_NEwomca = new RooPoisson( "pdf_NEwomca", "Na Poisson PDF", *rv_NEwomca, *rv_n_ewomc_a ) ;
 ///  pdf_NEwomcd = new RooPoisson( "pdf_NEwomcd", "Nd Poisson PDF", *rv_NEwomcd, *rv_n_ewomc_d ) ;
 ///  pdf_NEwomcsb = new RooPoisson( "pdf_NEwomcsb", "Nsb Poisson PDF", *rv_NEwomcsb, *rv_n_ewomc_sb ) ;
 ///  pdf_NEwomcsl_sig = new RooPoisson( "pdf_NEwomcsl_sig", "Nsl,sig Poisson PDF", *rv_NEwomcslsig, *rv_n_ewomc_sl_sig ) ;
 ///  pdf_NEwomcsl_sb = new RooPoisson( "pdf_NEwomcsl_sb", "Nsl,sb Poisson PDF", *rv_NEwomcslsb, *rv_n_ewomc_sl_sb ) ;



      pdf_Nqcdmc_sig  = new RooGaussian( "pdf_Nqcdmc_sig", "Gaussian pdf for Nqcdmc,sig",
                                          *rv_mu_qcdmc_sig, *rv_Nqcdmcsig, RooConst( Nqcdmcsigerr ) ) ;
      pdf_Nqcdmc_sb  = new RooGaussian( "pdf_Nqcdmc_sb", "Gaussian pdf for Nqcdmc,sb",
                                          *rv_mu_qcdmc_sb, *rv_Nqcdmcsb, RooConst( Nqcdmcsberr ) ) ;
      pdf_Nqcdmc_a   = new RooGaussian( "pdf_Nqcdmc_a", "Gaussian pdf for Nqcdmc,a",
                                          *rv_mu_qcdmc_a, *rv_Nqcdmca, RooConst( Nqcdmcaerr ) ) ;
      pdf_Nqcdmc_d   = new RooGaussian( "pdf_Nqcdmc_d", "Gaussian pdf for Nqcdmc,d",
                                          *rv_mu_qcdmc_d, *rv_Nqcdmcd, RooConst( Nqcdmcderr ) ) ;



      pdf_Eff_sf     = new RooGaussian( "pdf_Eff_sf", "Gaussian pdf for Efficiency scale factor",
                                          *rv_eff_sf, RooConst( EffScaleFactor ) , RooConst( EffScaleFactorErr ) ) ;


      pdf_lsf_WJmc  = new RooGaussian( "pdf_lsf_Wjmc", "Gaussian pdf for lsf, WJmc",
                                          *rv_lsf_wjmc, RooConst( lsf_WJmc ), RooConst( lsf_WJmc_err ) ) ;
      pdf_lsf_Znnmc  = new RooGaussian( "pdf_lsf_Znnmc", "Gaussian pdf for lsf, Znnmc",
                                          *rv_lsf_znnmc, RooConst( lsf_Znnmc ), RooConst( lsf_Znnmc_err ) ) ;
 ///  pdf_lsf_Ewomc  = new RooGaussian( "pdf_lsf_Ewomc", "Gaussian pdf for lsf, Ewomc",
 ///                                      *rv_lsf_ewomc, RooConst( lsf_Ewomc ), RooConst( lsf_Ewomc_err ) ) ;

      {
         RooArgSet pdflist ;
         pdflist.add( *pdf_Nsig ) ;
         pdflist.add( *pdf_Na ) ;
         pdflist.add( *pdf_Nd ) ;
         pdflist.add( *pdf_Nsb ) ;
         pdflist.add( *pdf_Nsl_sig ) ;
         pdflist.add( *pdf_Nsl_sb ) ;
         pdflist.add( *pdf_NWJmcsig ) ;
         pdflist.add( *pdf_NWJmca ) ;
         pdflist.add( *pdf_NWJmcd ) ;
         pdflist.add( *pdf_NWJmcsb ) ;
         pdflist.add( *pdf_NWJmcsl_sig ) ;
         pdflist.add( *pdf_NWJmcsl_sb ) ;
         pdflist.add( *pdf_NZnnmcsig ) ;
         pdflist.add( *pdf_NZnnmca ) ;
         pdflist.add( *pdf_NZnnmcd ) ;
         pdflist.add( *pdf_NZnnmcsb ) ;
 ///     pdflist.add( *pdf_NZnnmcsl_sig ) ;
 ///     pdflist.add( *pdf_NZnnmcsl_sb ) ;
 ///     pdflist.add( *pdf_NEwomcsig ) ;
 ///     pdflist.add( *pdf_NEwomca ) ;
 ///     pdflist.add( *pdf_NEwomcd ) ;
 ///     pdflist.add( *pdf_NEwomcsb ) ;
 ///     pdflist.add( *pdf_NEwomcsl_sig ) ;
 ///     pdflist.add( *pdf_NEwomcsl_sb ) ;
         pdflist.add( *pdf_Nqcdmc_sig ) ;
         pdflist.add( *pdf_Nqcdmc_sb ) ;
         pdflist.add( *pdf_Nqcdmc_a ) ;
         pdflist.add( *pdf_Nqcdmc_d ) ;
         pdflist.add( *pdf_Eff_sf ) ;
         pdflist.add( *pdf_lsf_WJmc ) ;
         pdflist.add( *pdf_lsf_Znnmc ) ;
 ////    pdflist.add( *pdf_lsf_Ewomc ) ;
         likelihood = new RooProdPdf( "likelihood", "ra2b likelihood", pdflist ) ;
      }


     //---- Define the list of observables.

       observedParametersList.add( *rv_Nsig ) ;
       observedParametersList.add( *rv_Na ) ;
       observedParametersList.add( *rv_Nd ) ;
       observedParametersList.add( *rv_Nsb ) ;
       observedParametersList.add( *rv_Nslsig ) ;
       observedParametersList.add( *rv_Nslsb ) ;

       observedParametersList.add( *rv_NWJmcsig ) ;
       observedParametersList.add( *rv_NWJmca ) ;
       observedParametersList.add( *rv_NWJmcd ) ;
       observedParametersList.add( *rv_NWJmcsb ) ;
       observedParametersList.add( *rv_NWJmcslsig ) ;
       observedParametersList.add( *rv_NWJmcslsb ) ;

       observedParametersList.add( *rv_NZnnmcsig ) ;
       observedParametersList.add( *rv_NZnnmca ) ;
       observedParametersList.add( *rv_NZnnmcd ) ;
       observedParametersList.add( *rv_NZnnmcsb ) ;
 //    observedParametersList.add( *rv_NZnnmcslsig ) ;
 //    observedParametersList.add( *rv_NZnnmcslsb ) ;

 //    observedParametersList.add( *rv_NEwomcsig ) ;
 //    observedParametersList.add( *rv_NEwomca ) ;
 //    observedParametersList.add( *rv_NEwomcd ) ;
 //    observedParametersList.add( *rv_NEwomcsb ) ;
 //    observedParametersList.add( *rv_NEwomcslsig ) ;
 //    observedParametersList.add( *rv_NEwomcslsb ) ;

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

    bool ra2bRoostatsClass2::reinitialize( ) {


       printf( "\n\n Opening input file : %s\n\n", initializeFile ) ;

       FILE* infp ;
       if ( (infp=fopen( initializeFile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", initializeFile ) ;
          return false ;
       }

       int Nsig(0), Na(0), Nd(0) ; //-- data counts in signal region, A, and D.
       int Nsb(0) ; //-- data counts in SB.
       int Nslsig(0) ; //-- data counts in SL, SIG.
       int Nslsb(0) ; //-- data counts in SL, SB.


       float Nqcdmcsig(0.), Nqcdmcsb(0.), Nqcdmca(0.), Nqcdmcd(0.) ; //-- QCD MC counts in SIG, SB, A, and D.
       float Nqcdmcslsig(0.), Nqcdmcslsb(0.) ; //-- QCD MC counts in SL; SIG and SB.
       float Nttbarmcsig(0.), Nttbarmcsb(0.), Nttbarmca(0.), Nttbarmcd(0.) ; //-- ttbar MC counts in SIG, SB, A, and D.
       float Nttbarmcslsig(0.), Nttbarmcslsb(0.) ; //-- ttbar MC counts in SL; SIG and SB.

       int NWJmcsig, NWJmca, NWJmcd, NWJmcsb, NWJmcslsig, NWJmcslsb ;
       int NZnnmcsig, NZnnmca, NZnnmcd, NZnnmcsb, NZnnmcslsig, NZnnmcslsb ;
       int NEwomcsig, NEwomca, NEwomcd, NEwomcsb, NEwomcslsig, NEwomcslsb ;

       float Nsusymcsig(0.), Nsusymca(0.), Nsusymcd(0.) ; //-- SUSY MC counts in SIG, A, and D.
       float Nsusymcsb(0.) ; //-- SUSY MC counts in SB.
       float Nsusymcslsig(0.) ; //-- SUSY MC counts in SL, SIG.
       float Nsusymcslsb(0.) ; //-- SUSY MC counts in SL, SB.

       //--- read in description line.
       printf("\n\n") ;
       char c(0) ;
       while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
       printf("\n\n") ;


       char label[1000] ;

      //--- Inputs generated with gen_roostats_input.c
      //    The order here must be consistent with the order there!

       fscanf( infp, "%s %g", label, &EffScaleFactor ) ;       printf( "%s %g\n", label, EffScaleFactor ) ;         
       fscanf( infp, "%s %g", label, &EffScaleFactorErr ) ;    printf( "%s %g\n", label, EffScaleFactorErr ) ;         
       fscanf( infp, "%s %d", label, &Nsig ) ;                 printf( "%s %d\n", label, Nsig ) ;         
       fscanf( infp, "%s %d", label, &Na ) ;                   printf( "%s %d\n", label, Na ) ;           
       fscanf( infp, "%s %d", label, &Nd ) ;                   printf( "%s %d\n", label, Nd ) ;           
       fscanf( infp, "%s %d", label, &Nsb ) ;                 printf( "%s %d\n", label, Nsb ) ;         
       fscanf( infp, "%s %d", label, &Nslsig ) ;              printf( "%s %d\n", label, Nslsig ) ;      
       fscanf( infp, "%s %d", label, &Nslsb ) ;               printf( "%s %d\n", label, Nslsb ) ;       
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
       fscanf( infp, "%s %g", label, &Nttbarmcsig ) ;          printf( "%s %g\n", label, Nttbarmcsig ) ;  
       fscanf( infp, "%s %g", label, &Nttbarmcsb ) ;           printf( "%s %g\n", label, Nttbarmcsb ) ;   
       fscanf( infp, "%s %g", label, &Nttbarmca ) ;            printf( "%s %g\n", label, Nttbarmca ) ;    
       fscanf( infp, "%s %g", label, &Nttbarmcd ) ;            printf( "%s %g\n", label, Nttbarmcd ) ;    
       fscanf( infp, "%s %g", label, &Nttbarmcslsig ) ;        printf( "%s %g\n", label, Nttbarmcslsig ) ;      
       fscanf( infp, "%s %g", label, &Nttbarmcslsb ) ;        printf( "%s %g\n", label, Nttbarmcslsb ) ;      

       fscanf( infp, "%s %g", label, &lsf_WJmc ) ;        printf( "%s %g\n", label, lsf_WJmc ) ;      
       fscanf( infp, "%s %g", label, &lsf_WJmc_err ) ;        printf( "%s %g\n", label, lsf_WJmc_err ) ;      
       fscanf( infp, "%s %d", label, &NWJmcsig ) ;        printf( "%s %d\n", label, NWJmcsig ) ;      
       fscanf( infp, "%s %d", label, &NWJmca ) ;        printf( "%s %d\n", label, NWJmca ) ;      
       fscanf( infp, "%s %d", label, &NWJmcd ) ;        printf( "%s %d\n", label, NWJmcd ) ;      
       fscanf( infp, "%s %d", label, &NWJmcsb ) ;        printf( "%s %d\n", label, NWJmcsb ) ;      
       fscanf( infp, "%s %d", label, &NWJmcslsig ) ;        printf( "%s %d\n", label, NWJmcslsig ) ;      
       fscanf( infp, "%s %d", label, &NWJmcslsb ) ;        printf( "%s %d\n", label, NWJmcslsb ) ;      

       fscanf( infp, "%s %g", label, &lsf_Znnmc ) ;        printf( "%s %g\n", label, lsf_Znnmc ) ;      
       fscanf( infp, "%s %g", label, &lsf_Znnmc_err ) ;        printf( "%s %g\n", label, lsf_Znnmc_err ) ;      
       fscanf( infp, "%s %d", label, &NZnnmcsig ) ;        printf( "%s %d\n", label, NZnnmcsig ) ;      
       fscanf( infp, "%s %d", label, &NZnnmca ) ;        printf( "%s %d\n", label, NZnnmca ) ;      
       fscanf( infp, "%s %d", label, &NZnnmcd ) ;        printf( "%s %d\n", label, NZnnmcd ) ;      
       fscanf( infp, "%s %d", label, &NZnnmcsb ) ;        printf( "%s %d\n", label, NZnnmcsb ) ;      
       fscanf( infp, "%s %d", label, &NZnnmcslsig ) ;        printf( "%s %d\n", label, NZnnmcslsig ) ;      
       fscanf( infp, "%s %d", label, &NZnnmcslsb ) ;        printf( "%s %d\n", label, NZnnmcslsb ) ;      

       fscanf( infp, "%s %g", label, &lsf_Ewomc ) ;        printf( "%s %g\n", label, lsf_Ewomc ) ;      
       fscanf( infp, "%s %g", label, &lsf_Ewomc_err ) ;        printf( "%s %g\n", label, lsf_Ewomc_err ) ;      
       fscanf( infp, "%s %d", label, &NEwomcsig ) ;        printf( "%s %d\n", label, NEwomcsig ) ;      
       fscanf( infp, "%s %d", label, &NEwomca ) ;        printf( "%s %d\n", label, NEwomca ) ;      
       fscanf( infp, "%s %d", label, &NEwomcd ) ;        printf( "%s %d\n", label, NEwomcd ) ;      
       fscanf( infp, "%s %d", label, &NEwomcsb ) ;        printf( "%s %d\n", label, NEwomcsb ) ;      
       fscanf( infp, "%s %d", label, &NEwomcslsig ) ;        printf( "%s %d\n", label, NEwomcslsig ) ;      
       fscanf( infp, "%s %d", label, &NEwomcslsb ) ;        printf( "%s %d\n", label, NEwomcslsb ) ;      

       fscanf( infp, "%s %g", label, &Nsusymcsig ) ;           printf( "%s %g\n", label, Nsusymcsig ) ;   
       fscanf( infp, "%s %g", label, &Nsusymca ) ;             printf( "%s %g\n", label, Nsusymca ) ;     
       fscanf( infp, "%s %g", label, &Nsusymcd ) ;             printf( "%s %g\n", label, Nsusymcd ) ;     
       fscanf( infp, "%s %g", label, &Nsusymcsb ) ;           printf( "%s %g\n", label, Nsusymcsb ) ;   
       fscanf( infp, "%s %g", label, &Nsusymcslsig ) ;        printf( "%s %g\n", label, Nsusymcslsig ) ;
       fscanf( infp, "%s %g", label, &Nsusymcslsb ) ;        printf( "%s %g\n", label, Nsusymcslsb ) ;


       printf("\n Done reading in %s\n\n", initializeFile ) ;
       fclose( infp ) ;






       //--- Print out a nice summary of the inputs.


       float Newmcsig = lsf_WJmc * NWJmcsig
                      + lsf_Znnmc * NZnnmcsig
                      + lsf_Ewomc * NEwomcsig ;

       float Newmcsb = lsf_WJmc * NWJmcsb
                      + lsf_Znnmc * NZnnmcsb
                      + lsf_Ewomc * NEwomcsb ;

       float Newmca = lsf_WJmc * NWJmca
                      + lsf_Znnmc * NZnnmca
                      + lsf_Ewomc * NEwomca ;

       float Newmcd = lsf_WJmc * NWJmcd
                      + lsf_Znnmc * NZnnmcd
                      + lsf_Ewomc * NEwomcd ;

       float Newmcslsig = lsf_WJmc * NWJmcslsig
                      + lsf_Znnmc * NZnnmcslsig
                      + lsf_Ewomc * NEwomcslsig ;

       float Newmcslsb = lsf_WJmc * NWJmcslsb
                      + lsf_Znnmc * NZnnmcslsb
                      + lsf_Ewomc * NEwomcslsb ;



       float Nsmsig = Nttbarmcsig + Nqcdmcsig + Newmcsig ;

       float Nsmsb = Nttbarmcsb + Nqcdmcsb + Newmcsb ;

       float Nsma = Nttbarmca + Nqcdmca + Newmca ;
       float Nsmd = Nttbarmcd + Nqcdmcd + Newmcd ;

       float Nsmslsig = Nttbarmcslsig + Nqcdmcslsig + Newmcslsig ;

       float Nsmslsb = Nttbarmcslsb + Nqcdmcslsb + Newmcslsb ;



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
       printf("\n\n\n") ;



       printf(" --- Defining observables.\n" ) ;

      //-- data counts in signal region, A, and D, SB, SL SB, SL SIG.

      rv_Nsig -> setVal( Nsig ) ;
      rv_Na -> setVal( Na ) ;
      rv_Nd -> setVal( Nd ) ;
      rv_Nsb -> setVal( Nsb ) ;
      rv_Nslsig -> setVal( Nslsig ) ;
      rv_Nslsb -> setVal( Nslsb ) ;





      //-- QCD MC counts

      rv_Nqcdmca   -> setVal( Nqcdmca ) ;
      rv_Nqcdmcd   -> setVal( Nqcdmcd ) ;
      rv_Nqcdmcsig -> setVal( Nqcdmcsig ) ;
      rv_Nqcdmcsb  -> setVal( Nqcdmcsb ) ;






      //-- EW MC counts

      rv_NWJmcsig   -> setVal( NWJmcsig ) ;
      rv_NWJmca     -> setVal( NWJmca ) ;
      rv_NWJmcd     -> setVal( NWJmcd ) ;
      rv_NWJmcsb    -> setVal( NWJmcsb ) ;
      rv_NWJmcslsig -> setVal( NWJmcslsig ) ;
      rv_NWJmcslsb  -> setVal( NWJmcslsb ) ;

      rv_NZnnmcsig   -> setVal( NZnnmcsig ) ;
      rv_NZnnmca     -> setVal( NZnnmca ) ;
      rv_NZnnmcd     -> setVal( NZnnmcd ) ;
      rv_NZnnmcsb    -> setVal( NZnnmcsb ) ;
      rv_NZnnmcslsig -> setVal( NZnnmcslsig ) ;
      rv_NZnnmcslsb  -> setVal( NZnnmcslsb ) ;

      rv_NEwomcsig   -> setVal( NEwomcsig ) ;
      rv_NEwomca     -> setVal( NEwomca ) ;
      rv_NEwomcd     -> setVal( NEwomcd ) ;
      rv_NEwomcsb    -> setVal( NEwomcsb ) ;
      rv_NEwomcslsig -> setVal( NEwomcslsig ) ;
      rv_NEwomcslsb  -> setVal( NEwomcslsb ) ;











      //++++++++ Parameters of the likelihood

      printf(" --- Defining parameters.\n" ) ;




     //-- Counts in SIG

      if ( useSigBgVars ) {
         rrv_mu_ttbar_sig   -> setVal( Nttbarmcsig ) ;  //-- this is a starting value only.
      }
      rv_mu_qcd_sig     -> setVal( Nqcdmcsig ) ;  //-- this is a starting value only.
      rv_mu_susy_sig    -> setVal( 0. ) ;  //-- this is a starting value only.
      rv_mu_susymc_sig  -> setVal( 0. ) ;

      rv_mu_susymc_sig -> setConstant(kTRUE) ;





     //-- Counts in SB

      if ( !useSigBgVars ) {
         rrv_mu_ttbar_sb   -> setVal( Nttbarmcsb ) ;  //-- this is a starting value only.
      }
      rv_mu_qcd_sb    -> setVal( Nqcdmcsb ) ;  //-- this is a starting value only.
      rv_mu_susymc_sb -> setVal( 0. ) ;

      rv_mu_susymc_sb -> setConstant(kTRUE) ;






     //-- Counts in A, signal selection

      rv_mu_ttbar_a   -> setVal( Nttbarmca ) ;
      rv_mu_susymc_a  -> setVal( 0. ) ;

      rv_mu_susymc_a  ->setConstant(kTRUE) ;
      rv_mu_ttbar_a   ->setConstant(kTRUE) ;




     //-- Counts in D, signal selection

      rv_mu_ttbar_d   -> setVal( Nttbarmcd ) ;
      rv_mu_susymc_d  -> setVal( 0. ) ;

      rv_mu_ttbar_d   -> setConstant(kTRUE) ;
      rv_mu_susymc_d  -> setConstant(kTRUE) ;




     //-- Single Lepton (SL) counts, SIG region

      rv_mu_sl_ttbar_sig  -> setVal( Nslsig - Newmcslsig ) ;
      rv_mu_sl_susymc_sig -> setVal( 0. ) ;

      rv_mu_sl_susymc_sig -> setConstant( kTRUE ) ;





     //-- Single Lepton (SL) counts, SB region

      rv_mu_sl_ttbar_sb  -> setVal( Nslsb - Newmcslsb ) ;
      rv_mu_sl_susymc_sb -> setVal( 0. ) ;

      rv_mu_sl_susymc_sb -> setConstant( kTRUE ) ;






     //-- QCD MC, SIG,SB,A,D counts, signal selection

      rv_mu_qcdmc_sig -> setVal( Nqcdmcsig ) ;
      rv_mu_qcdmc_sb  -> setVal( Nqcdmcsb ) ;
      rv_mu_qcdmc_a   -> setVal( Nqcdmca ) ;
      rv_mu_qcdmc_d   -> setVal( Nqcdmcd ) ;



     //-- EW MC counts

      rv_mu_wjmc_sig   -> setVal( NWJmcsig ) ;
      rv_mu_wjmc_a     -> setVal( NWJmca ) ;
      rv_mu_wjmc_d     -> setVal( NWJmcd ) ;
      rv_mu_wjmc_sb    -> setVal( NWJmcsb ) ;
      rv_mu_wjmc_slsig -> setVal( NWJmcslsig ) ;
      rv_mu_wjmc_slsb  -> setVal( NWJmcslsb ) ;

      rv_mu_znnmc_sig   -> setVal( NZnnmcsig ) ;
      rv_mu_znnmc_a     -> setVal( NZnnmca ) ;
      rv_mu_znnmc_d     -> setVal( NZnnmcd ) ;
      rv_mu_znnmc_sb    -> setVal( NZnnmcsb ) ;
      rv_mu_znnmc_slsig -> setVal( NZnnmcslsig ) ;
      rv_mu_znnmc_slsb  -> setVal( NZnnmcslsb ) ;

      rv_mu_ewomc_sig   -> setVal( NEwomcsig ) ;
      rv_mu_ewomc_a     -> setVal( NEwomca ) ;
      rv_mu_ewomc_d     -> setVal( NEwomcd ) ;
      rv_mu_ewomc_sb    -> setVal( NEwomcsb ) ;
      rv_mu_ewomc_slsig -> setVal( NEwomcslsig ) ;
      rv_mu_ewomc_slsb  -> setVal( NEwomcslsb ) ;








    //-- Efficiency scale factor.

      rv_eff_sf  -> setVal( EffScaleFactor ) ;

      rv_lsf_wjmc  -> setVal( lsf_WJmc ) ;
      rv_lsf_znnmc  -> setVal( lsf_Znnmc ) ;
      rv_lsf_ewomc  -> setVal( lsf_Ewomc ) ;




       initialized = true ;

       return true ;


    } // reinitialize.

  //===================================================================

//  bool ra2bRoostatsClass2::readToyDataset( const char* inputRootFile, int dsIndex ) {

//        if ( ! initialized ) {
//           printf("\n\n *** Call initialize first.\n\n") ;
//           return false ;
//        }

//       TFile intoyfile( inputRootFile, "READ" ) ;
//       gDirectory->ls() ;
//       TTree* toytree = (TTree*) gDirectory->FindObjectAny("likelihoodData") ;
//       TTree* toytruetree = (TTree*) gDirectory->FindObjectAny("toytruetree") ;

//       if ( toytree == 0 ) {
//          printf(" \n\n *** Can't find TTree likelihoodData in file %s\n", inputRootFile ) ;
//       } else {
//          toytree->Print("toponly") ;
//       }

//       if ( toytruetree == 0 ) {
//          printf(" \n\n *** Can't find TTree toytruetree in file %s\n", inputRootFile ) ;
//       } else {
//          toytruetree->Print("toponly") ;
//       }


//       Double_t toyNsig;
//       Double_t toyNa;
//       Double_t toyNd;
//       Double_t toyNsb1;
//       Double_t toyNsb2;
//       Double_t toyNsb3;
//       Double_t toyNsb4;
//       Double_t toyNsb5;
//       Double_t toyNlsb1;
//       Double_t toyNlsb2;
//       Double_t toyNlsb3;
//       Double_t toyNlsb4;
//       Double_t toyNlsb5;
//       Double_t toyNslsig1;
//       Double_t toyNslsig2;
//       Double_t toyNslsig3;
//       Double_t toyNslsig4;
//       Double_t toyNslsig5;
//       Double_t toyNslsb1;
//       Double_t toyNslsb2;
//       Double_t toyNslsb3;
//       Double_t toyNslsb4;
//       Double_t toyNslsb5;
//       Double_t toyNslmsb1;
//       Double_t toyNslmsb2;
//       Double_t toyNslmsb3;
//       Double_t toyNslmsb4;
//       Double_t toyNslmsb5;

//       Double_t toyNqcdmca ;
//       Double_t toyNqcdmcd ;
//       Double_t toyNqcdmcsig ;
//       Double_t toyNqcdmcsb ;

//       // List of branches
//       TBranch *b_toyNsig;   //!
//       TBranch *b_toyNa;   //!
//       TBranch *b_toyNd;   //!
//       TBranch *b_toyNsb1;   //!
//       TBranch *b_toyNsb2;   //!
//       TBranch *b_toyNsb3;   //!
//       TBranch *b_toyNsb4;   //!
//       TBranch *b_toyNsb5;   //!
//       TBranch *b_toyNlsb1;   //!
//       TBranch *b_toyNlsb2;   //!
//       TBranch *b_toyNlsb3;   //!
//       TBranch *b_toyNlsb4;   //!
//       TBranch *b_toyNlsb5;   //!
//       TBranch *b_toyNslsig1;   //!
//       TBranch *b_toyNslsig2;   //!
//       TBranch *b_toyNslsig3;   //!
//       TBranch *b_toyNslsig4;   //!
//       TBranch *b_toyNslsig5;   //!
//       TBranch *b_toyNslsb1;   //!
//       TBranch *b_toyNslsb2;   //!
//       TBranch *b_toyNslsb3;   //!
//       TBranch *b_toyNslsb4;   //!
//       TBranch *b_toyNslsb5;   //!
//       TBranch *b_toyNslmsb1;   //!
//       TBranch *b_toyNslmsb2;   //!
//       TBranch *b_toyNslmsb3;   //!
//       TBranch *b_toyNslmsb4;   //!
//       TBranch *b_toyNslmsb5;   //!

//       TBranch *b_toyNqcdmca ;
//       TBranch *b_toyNqcdmcd ;
//       TBranch *b_toyNqcdmcsig ;
//       TBranch *b_toyNqcdmcsb ;

//       TBranch *b_toy_mu0_ttbar_sig ;
//       TBranch *b_toy_mu0_qcd_sig ;
//       TBranch *b_toy_mu0_ttbar_sb ;
//       TBranch *b_toy_mu0_qcd_sb ;
//       TBranch *b_toy_mu0_susy_sig ;
//       TBranch *b_toy_mu0_allbg_sig ;

//       printf("\n\n Setting branch addresses...\n") ;
//       cout << " Toy tree pointer: " << toytree << endl ;
//       toytree->SetBranchAddress("Nsig", &toyNsig, &b_toyNsig ) ;
//       toytree->SetBranchAddress("Na", &toyNa, &b_toyNa);
//       toytree->SetBranchAddress("Nd", &toyNd, &b_toyNd);
//       toytree->SetBranchAddress("Nsb1", &toyNsb1, &b_toyNsb1);
//       toytree->SetBranchAddress("Nsb2", &toyNsb2, &b_toyNsb2);
//       toytree->SetBranchAddress("Nsb3", &toyNsb3, &b_toyNsb3);
//       toytree->SetBranchAddress("Nsb4", &toyNsb4, &b_toyNsb4);
//       toytree->SetBranchAddress("Nsb5", &toyNsb5, &b_toyNsb5);
//       toytree->SetBranchAddress("Nlsb1", &toyNlsb1, &b_toyNlsb1);
//       toytree->SetBranchAddress("Nlsb2", &toyNlsb2, &b_toyNlsb2);
//       toytree->SetBranchAddress("Nlsb3", &toyNlsb3, &b_toyNlsb3);
//       toytree->SetBranchAddress("Nlsb4", &toyNlsb4, &b_toyNlsb4);
//       toytree->SetBranchAddress("Nlsb5", &toyNlsb5, &b_toyNlsb5);
//       toytree->SetBranchAddress("Nslsig1", &toyNslsig1, &b_toyNslsig1);
//       toytree->SetBranchAddress("Nslsig2", &toyNslsig2, &b_toyNslsig2);
//       toytree->SetBranchAddress("Nslsig3", &toyNslsig3, &b_toyNslsig3);
//       toytree->SetBranchAddress("Nslsig4", &toyNslsig4, &b_toyNslsig4);
//       toytree->SetBranchAddress("Nslsig5", &toyNslsig5, &b_toyNslsig5);
//       toytree->SetBranchAddress("Nslsb1", &toyNslsb1, &b_toyNslsb1);
//       toytree->SetBranchAddress("Nslsb2", &toyNslsb2, &b_toyNslsb2);
//       toytree->SetBranchAddress("Nslsb3", &toyNslsb3, &b_toyNslsb3);
//       toytree->SetBranchAddress("Nslsb4", &toyNslsb4, &b_toyNslsb4);
//       toytree->SetBranchAddress("Nslsb5", &toyNslsb5, &b_toyNslsb5);
//       toytree->SetBranchAddress("Nslmsb1", &toyNslmsb1, &b_toyNslmsb1);
//       toytree->SetBranchAddress("Nslmsb2", &toyNslmsb2, &b_toyNslmsb2);
//       toytree->SetBranchAddress("Nslmsb3", &toyNslmsb3, &b_toyNslmsb3);
//       toytree->SetBranchAddress("Nslmsb4", &toyNslmsb4, &b_toyNslmsb4);
//       toytree->SetBranchAddress("Nslmsb5", &toyNslmsb5, &b_toyNslmsb5);

//       toytree->SetBranchAddress("Nqcdmca"  , &toyNqcdmca  , &b_toyNqcdmca);
//       toytree->SetBranchAddress("Nqcdmcd"  , &toyNqcdmcd  , &b_toyNqcdmcd);
//       toytree->SetBranchAddress("Nqcdmcsig", &toyNqcdmcsig, &b_toyNqcdmcsig);
//       toytree->SetBranchAddress("Nqcdmcsb" , &toyNqcdmcsb , &b_toyNqcdmcsb);

//       toytruetree->SetBranchAddress("mu0_ttbar_sig" , &toy_mu0_ttbar_sig , &b_toy_mu0_ttbar_sig);
//       toytruetree->SetBranchAddress("mu0_qcd_sig"   , &toy_mu0_qcd_sig   , &b_toy_mu0_qcd_sig);
//       toytruetree->SetBranchAddress("mu0_ttbar_sb"  , &toy_mu0_ttbar_sb  , &b_toy_mu0_ttbar_sb);
//       toytruetree->SetBranchAddress("mu0_qcd_sb"    , &toy_mu0_qcd_sb    , &b_toy_mu0_qcd_sb);
//       toytruetree->SetBranchAddress("mu0_susy_sig"  , &toy_mu0_susy_sig  , &b_toy_mu0_susy_sig);
//       toytruetree->SetBranchAddress("mu0_allbg_sig" , &toy_mu0_allbg_sig , &b_toy_mu0_allbg_sig);

//       printf("\n Done.\n\n") ;


//       toytruetree->GetEntry(0) ;

//       toytree->GetEntry( dsIndex ) ;

//       rv_Nsig -> setVal( toyNsig ) ;
//       rv_Na -> setVal( toyNa ) ;
//       rv_Nd -> setVal( toyNd ) ;
//       rv_Nsb1 -> setVal( toyNsb1 ) ;
//       rv_Nsb2 -> setVal( toyNsb2 ) ;
//       rv_Nsb3 -> setVal( toyNsb3 ) ;
//       rv_Nsb4 -> setVal( toyNsb4 ) ;
//       rv_Nsb5 -> setVal( toyNsb5 ) ;
//       rv_Nlsb1 -> setVal( toyNlsb1 ) ;
//       rv_Nlsb2 -> setVal( toyNlsb2 ) ;
//       rv_Nlsb3 -> setVal( toyNlsb3 ) ;
//       rv_Nlsb4 -> setVal( toyNlsb4 ) ;
//       rv_Nlsb5 -> setVal( toyNlsb5 ) ;
//       rv_Nslsig1 -> setVal( toyNslsig1 ) ;
//       rv_Nslsig2 -> setVal( toyNslsig2 ) ;
//       rv_Nslsig3 -> setVal( toyNslsig3 ) ;
//       rv_Nslsig4 -> setVal( toyNslsig4 ) ;
//       rv_Nslsig5 -> setVal( toyNslsig5 ) ;
//       rv_Nslsb1 -> setVal( toyNslsb1 ) ;
//       rv_Nslsb2 -> setVal( toyNslsb2 ) ;
//       rv_Nslsb3 -> setVal( toyNslsb3 ) ;
//       rv_Nslsb4 -> setVal( toyNslsb4 ) ;
//       rv_Nslsb5 -> setVal( toyNslsb5 ) ;
//       rv_Nslmsb1 -> setVal( toyNslmsb1 ) ;
//       rv_Nslmsb2 -> setVal( toyNslmsb2 ) ;
//       rv_Nslmsb3 -> setVal( toyNslmsb3 ) ;
//       rv_Nslmsb4 -> setVal( toyNslmsb4 ) ;
//       rv_Nslmsb5 -> setVal( toyNslmsb5 ) ;

//       rv_Nqcdmca   -> setVal( toyNqcdmca ) ;
//       rv_Nqcdmcd   -> setVal( toyNqcdmcd ) ;
//       rv_Nqcdmcsig -> setVal( toyNqcdmcsig ) ;
//       rv_Nqcdmcsb  -> setVal( toyNqcdmcsb ) ;


//       if ( dsObserved != 0x0 ) delete dsObserved ;

//       dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
//                                    observedParametersList ) ;
//       dsObserved->add( observedParametersList ) ;
//       printf("\n\n") ;
//       dsObserved->printMultiline(cout, 1, kTRUE, "") ;
//       printf("\n\n") ;


//       return true ;

//  } // readToyDataset

  //===================================================================



//  bool ra2bRoostatsClass2::readTextDataset( const char* infile ) {

//        if ( ! initialized ) {
//           printf("\n\n *** Call initialize first.\n\n") ;
//           return false ;
//        }


//     printf( "\n\n Opening input file : %s\n\n", infile ) ;

//     FILE* infp ;
//     if ( (infp=fopen( infile,"r"))==NULL ) {
//        printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
//        return false ;
//     }

//     int N3jmBins(0) ; //-- number of 3-jet mass bins.

//     int Nsig(0), Na(0), Nd(0) ; //-- data counts in signal region, A, and D.
//     int Nsb1(0), Nsb2(0), Nsb3(0), Nsb4(0), Nsb5(0) ; //-- data counts in 5 3-jet mass bins of SB.
//     int Nlsb1(0), Nlsb2(0), Nlsb3(0), Nlsb4(0), Nlsb5(0) ; //-- data counts in 5 3-jet mass bins of LSB.
//     int Nslsig1(0), Nslsig2(0), Nslsig3(0), Nslsig4(0), Nslsig5(0) ; //-- data counts in 5 3-jet mass bins of SL, SIG.
//     int Nslsb1(0), Nslsb2(0), Nslsb3(0), Nslsb4(0), Nslsb5(0) ; ; //-- data counts in 5 3-jet mass bins of SL, SB.
//     int Nslmsb1(0), Nslmsb2(0), Nslmsb3(0), Nslmsb4(0), Nslmsb5(0) ; //-- data counts in 5 3-jet mass bins of SL, MSB.


//     //--- read in description line.
//     printf("\n\n") ;
//     char c(0) ;
//     while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
//     printf("\n\n") ;


//     char label[1000] ;

//    //--- Inputs generated with gen_roostats_input.c
//    //    The order here must be consistent with the order there!

//     fscanf( infp, "%s %d", label, &N3jmBins ) ;             printf( "%s %d\n", label, N3jmBins ) ;                
//     fscanf( infp, "%s %d", label, &Nsig ) ;                 printf( "%s %d\n", label, Nsig ) ;         
//     fscanf( infp, "%s %d", label, &Na ) ;                   printf( "%s %d\n", label, Na ) ;           
//     fscanf( infp, "%s %d", label, &Nd ) ;                   printf( "%s %d\n", label, Nd ) ;           
//     fscanf( infp, "%s %d", label, &Nsb1 ) ;                 printf( "%s %d\n", label, Nsb1 ) ;         
//     fscanf( infp, "%s %d", label, &Nsb2 ) ;                 printf( "%s %d\n", label, Nsb2 ) ;         
//     fscanf( infp, "%s %d", label, &Nsb3 ) ;                 printf( "%s %d\n", label, Nsb3 ) ;         
//     fscanf( infp, "%s %d", label, &Nsb4 ) ;                 printf( "%s %d\n", label, Nsb4 ) ;         
//     fscanf( infp, "%s %d", label, &Nsb5 ) ;                 printf( "%s %d\n", label, Nsb5 ) ;         
//     fscanf( infp, "%s %d", label, &Nlsb1 ) ;                printf( "%s %d\n", label, Nlsb1 ) ;        
//     fscanf( infp, "%s %d", label, &Nlsb2 ) ;                printf( "%s %d\n", label, Nlsb2 ) ;        
//     fscanf( infp, "%s %d", label, &Nlsb3 ) ;                printf( "%s %d\n", label, Nlsb3 ) ;        
//     fscanf( infp, "%s %d", label, &Nlsb4 ) ;                printf( "%s %d\n", label, Nlsb4 ) ;        
//     fscanf( infp, "%s %d", label, &Nlsb5 ) ;                printf( "%s %d\n", label, Nlsb5 ) ;        
//     fscanf( infp, "%s %d", label, &Nslsig1 ) ;              printf( "%s %d\n", label, Nslsig1 ) ;      
//     fscanf( infp, "%s %d", label, &Nslsig2 ) ;              printf( "%s %d\n", label, Nslsig2 ) ;      
//     fscanf( infp, "%s %d", label, &Nslsig3 ) ;              printf( "%s %d\n", label, Nslsig3 ) ;      
//     fscanf( infp, "%s %d", label, &Nslsig4 ) ;              printf( "%s %d\n", label, Nslsig4 ) ;      
//     fscanf( infp, "%s %d", label, &Nslsig5 ) ;              printf( "%s %d\n", label, Nslsig5 ) ;      
//     fscanf( infp, "%s %d", label, &Nslsb1 ) ;               printf( "%s %d\n", label, Nslsb1 ) ;       
//     fscanf( infp, "%s %d", label, &Nslsb2 ) ;               printf( "%s %d\n", label, Nslsb2 ) ;       
//     fscanf( infp, "%s %d", label, &Nslsb3 ) ;               printf( "%s %d\n", label, Nslsb3 ) ;       
//     fscanf( infp, "%s %d", label, &Nslsb4 ) ;               printf( "%s %d\n", label, Nslsb4 ) ;       
//     fscanf( infp, "%s %d", label, &Nslsb5 ) ;               printf( "%s %d\n", label, Nslsb5 ) ;       
//     fscanf( infp, "%s %d", label, &Nslmsb1 ) ;              printf( "%s %d\n", label, Nslmsb1 ) ;      
//     fscanf( infp, "%s %d", label, &Nslmsb2 ) ;              printf( "%s %d\n", label, Nslmsb2 ) ;      
//     fscanf( infp, "%s %d", label, &Nslmsb3 ) ;              printf( "%s %d\n", label, Nslmsb3 ) ;      
//     fscanf( infp, "%s %d", label, &Nslmsb4 ) ;              printf( "%s %d\n", label, Nslmsb4 ) ;      
//     fscanf( infp, "%s %d", label, &Nslmsb5 ) ;              printf( "%s %d\n", label, Nslmsb5 ) ;      

//     rv_Nsig -> setVal( Nsig ) ;
//     rv_Na -> setVal( Na ) ;
//     rv_Nd -> setVal( Nd ) ;
//     rv_Nsb1 -> setVal( Nsb1 ) ;
//     rv_Nsb2 -> setVal( Nsb2 ) ;
//     rv_Nsb3 -> setVal( Nsb3 ) ;
//     rv_Nsb4 -> setVal( Nsb4 ) ;
//     rv_Nsb5 -> setVal( Nsb5 ) ;
//     rv_Nlsb1 -> setVal( Nlsb1 ) ;
//     rv_Nlsb2 -> setVal( Nlsb2 ) ;
//     rv_Nlsb3 -> setVal( Nlsb3 ) ;
//     rv_Nlsb4 -> setVal( Nlsb4 ) ;
//     rv_Nlsb5 -> setVal( Nlsb5 ) ;
//     rv_Nslsig1 -> setVal( Nslsig1 ) ;
//     rv_Nslsig2 -> setVal( Nslsig2 ) ;
//     rv_Nslsig3 -> setVal( Nslsig3 ) ;
//     rv_Nslsig4 -> setVal( Nslsig4 ) ;
//     rv_Nslsig5 -> setVal( Nslsig5 ) ;
//     rv_Nslsb1 -> setVal( Nslsb1 ) ;
//     rv_Nslsb2 -> setVal( Nslsb2 ) ;
//     rv_Nslsb3 -> setVal( Nslsb3 ) ;
//     rv_Nslsb4 -> setVal( Nslsb4 ) ;
//     rv_Nslsb5 -> setVal( Nslsb5 ) ;
//     rv_Nslmsb1 -> setVal( Nslmsb1 ) ;
//     rv_Nslmsb2 -> setVal( Nslmsb2 ) ;
//     rv_Nslmsb3 -> setVal( Nslmsb3 ) ;
//     rv_Nslmsb4 -> setVal( Nslmsb4 ) ;
//     rv_Nslmsb5 -> setVal( Nslmsb5 ) ;


//     if ( dsObserved != 0x0 ) delete dsObserved ;

//     dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
//                                  observedParametersList ) ;
//     dsObserved->add( observedParametersList ) ;
//     printf("\n\n") ;
//     dsObserved->printMultiline(cout, 1, kTRUE, "") ;
//     printf("\n\n") ;


//     return true ;



//  } // readTextDataset


  //===================================================================


//  bool ra2bRoostatsClass2::doToyStudy( const char* inputRootFile, const char* outputRootFile, int dsFirst, int nToys ) {


//       double toyNsig ;
//       double toy_mu_ttbar_sig ;
//       double toy_mu_qcd_sig ;
//       double toy_mu_ttbar_sb ;
//       double toy_mu_qcd_sb ;

//       double toy_mu_ttbar_sig_err ;
//       double toy_mu_qcd_sig_err ;
//       double toy_mu_ttbar_sb_err ;
//       double toy_mu_qcd_sb_err ;

//       double toy_rho_ttbar_qcd_sig ;
//       double toy_rho_ttbar_qcd_sb ;

//       double toy_mu_susy_sig ;
//       double toy_mu_susy_sig_err ;
//       double toy_mu_susy_sig_ul ;

//       double toy_mu_allbg_sig ;
//       double toy_mu_allbg_sig_err ;


//       TTree* toyfittree = new TTree("toyfittree","ra2b toy fit tree") ;

//       toyfittree->Branch("Nsig"              , &toyNsig, "toyNsig/D" ) ;
//       toyfittree->Branch("mu_ttbar_sig"      , &toy_mu_ttbar_sig, "mu_ttbar_sig/D" ) ;
//       toyfittree->Branch("mu_qcd_sig"        , &toy_mu_qcd_sig, "mu_qcd_sig/D" ) ;
//       toyfittree->Branch("mu_ttbar_sb"       , &toy_mu_ttbar_sb, "mu_ttbar_sb/D" ) ;
//       toyfittree->Branch("mu_qcd_sb"         , &toy_mu_qcd_sb, "mu_qcd_sb/D" ) ;
//       toyfittree->Branch("mu_ttbar_sig_err"  , &toy_mu_ttbar_sig_err, "mu_ttbar_sig_err/D" ) ;
//       toyfittree->Branch("mu_qcd_sig_err"    , &toy_mu_qcd_sig_err, "mu_qcd_sig_err/D" ) ;
//       toyfittree->Branch("mu_ttbar_sb_err"   , &toy_mu_ttbar_sb_err, "mu_ttbar_sb_err/D" ) ;
//       toyfittree->Branch("mu_qcd_sb_err"     , &toy_mu_qcd_sb_err, "mu_qcd_sb_err/D" ) ;

//       toyfittree->Branch("rho_ttbar_qcd_sig" , &toy_rho_ttbar_qcd_sig     , "rho_ttbar_qcd_sig/D" ) ;
//       toyfittree->Branch("rho_ttbar_qcd_sb"  , &toy_rho_ttbar_qcd_sb      , "rho_ttbar_qcd_sb/D" ) ;

//       toyfittree->Branch("mu_susy_sig"       , &toy_mu_susy_sig       , "mu_susy_sig/D" ) ;
//       toyfittree->Branch("mu_susy_sig_err"   , &toy_mu_susy_sig_err   , "mu_susy_sig_err/D" ) ;
//       toyfittree->Branch("mu_susy_sig_ul"    , &toy_mu_susy_sig_ul    , "mu_susy_sig_ul/D" ) ;
//       toyfittree->Branch("mu_allbg_sig"      , &toy_mu_allbg_sig      , "mu_allbg_sig/D" ) ;
//       toyfittree->Branch("mu_allbg_sig_err"  , &toy_mu_allbg_sig_err  , "mu_allbg_sig_err/D" ) ;

//       toyfittree->Branch("mu0_ttbar_sig"      , &toy_mu0_ttbar_sig, "mu0_ttbar_sig/D" ) ;
//       toyfittree->Branch("mu0_qcd_sig"        , &toy_mu0_qcd_sig, "mu0_qcd_sig/D" ) ;
//       toyfittree->Branch("mu0_ttbar_sb"       , &toy_mu0_ttbar_sb, "mu0_ttbar_sb/D" ) ;
//       toyfittree->Branch("mu0_qcd_sb"         , &toy_mu0_qcd_sb, "mu0_qcd_sb/D" ) ;
//       toyfittree->Branch("mu0_susy_sig"       , &toy_mu0_susy_sig       , "mu0_susy_sig/D" ) ;
//       toyfittree->Branch("mu0_allbg_sig"      , &toy_mu0_allbg_sig      , "mu0_allbg_sig/D" ) ;


//       for ( int dsi=dsFirst; dsi<(dsFirst+nToys); dsi++ ) {

//          printf("\n\n\n ============ Begin fit of dataset %d ====================== \n\n", dsi ) ;

//          reinitialize() ;

//          readToyDataset( inputRootFile, dsi ) ;

//          doFit() ;

//          toyNsig = rv_Nsig->getVal() ;

//          if ( useSigBgVars ) {
//             toy_mu_ttbar_sig = ((RooRealVar*)rv_mu_ttbar_sig)->getVal() ;
//             toy_mu_qcd_sig = ((RooRealVar*)rv_mu_qcd_sig)->getVal() ;
//             toy_mu_ttbar_sb = ((RooFormulaVar*)rv_mu_ttbar_sb)->getVal() ;
//             toy_mu_qcd_sb = ((RooFormulaVar*)rv_mu_qcd_sb)->getVal() ;
//             toy_mu_ttbar_sig_err = ((RooRealVar*)rv_mu_ttbar_sig)->getError() ;
//             toy_mu_qcd_sig_err = ((RooRealVar*)rv_mu_qcd_sig)->getError() ;
//             toy_mu_ttbar_sb_err = 0. ;
//             toy_mu_qcd_sb_err = 0. ;
//          } else {
//             toy_mu_ttbar_sig = ((RooFormulaVar*)rv_mu_ttbar_sig)->getVal() ;
//             toy_mu_qcd_sig = ((RooFormulaVar*)rv_mu_qcd_sig)->getVal() ;
//             toy_mu_ttbar_sb = ((RooRealVar*)rv_mu_ttbar_sb)->getVal() ;
//             toy_mu_qcd_sb = ((RooRealVar*)rv_mu_qcd_sb)->getVal() ;
//             toy_mu_ttbar_sb_err = ((RooRealVar*)rv_mu_ttbar_sb)->getError() ;
//             toy_mu_qcd_sb_err = ((RooRealVar*)rv_mu_qcd_sb)->getError() ;
//             toy_mu_ttbar_sig_err = 0. ;
//             toy_mu_qcd_sig_err = 0. ;
//          }

//          toy_mu_susy_sig = rv_mu_susy_sig->getVal() ;
//          toy_mu_susy_sig_err = rv_mu_susy_sig->getError() ;

//          toy_rho_ttbar_qcd_sig = fitResult->correlation( "mu_ttbar_sig", "mu_qcd_sig" ) ;
//          toy_rho_ttbar_qcd_sb  = fitResult->correlation( "mu_ttbar_sb" , "mu_qcd_sb"  ) ;

//     //--- Profile likelihood for signal susy yield.

//          ProfileLikelihoodCalculator plc_susy_sig( *dsObserved, *likelihood, RooArgSet( *rv_mu_susy_sig ) ) ;
//          plc_susy_sig.SetTestSize(0.05) ;
//          ConfInterval* plinterval_susy_sig = plc_susy_sig.GetInterval() ;
//          float susy_sig_ul = ((LikelihoodInterval*) plinterval_susy_sig)->UpperLimit(*((RooRealVar*)rv_mu_susy_sig)) ;
//          float susy_sig_ll = ((LikelihoodInterval*) plinterval_susy_sig)->LowerLimit(*((RooRealVar*)rv_mu_susy_sig)) ;
//          printf("\n\n") ;
//          printf("    susy, SIG 95%% CL interval  [%5.1f, %5.1f]\n\n", susy_sig_ll, susy_sig_ul) ;
//          printf("\n\n") ;
//          toy_mu_susy_sig_ul = susy_sig_ul ;
//          delete plinterval_susy_sig ;


//          toy_mu_allbg_sig = toy_mu_ttbar_sig + toy_mu_qcd_sig + rv_mu_ew_sig->getVal() ;
//          { double err2 = pow(toy_mu_ttbar_sig_err,2) + pow(toy_mu_qcd_sig_err,2)
//                       + 2.*toy_rho_ttbar_qcd_sig*toy_mu_ttbar_sig_err*toy_mu_qcd_sig_err ;
//             toy_mu_allbg_sig_err =  0. ;
//             if ( err2 > 0. ) toy_mu_allbg_sig_err = sqrt( err2 ) ;
//          }

//          toyfittree->Fill() ;

//          printf("\n\n\n ============ End fit of dataset %d ====================== \n\n", dsi ) ;

//       } // dsi.

//       TFile outputfile(outputRootFile,"recreate") ;
//       toyfittree->Write() ;
//       outputfile.Close() ;

//       return true ;


//  } // doToyStudy












  //===================================================================

    bool ra2bRoostatsClass2::susyScanNoContam( const char* inputScanFile, double dataLumi ) {


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
          int    nsel ;
          int    n_a ;
          int    n_d ;
          int    n_sb ;
          int    n_sl_sb ;
          int    n_sl_sig ;

          fscanf( infp, "%f %f %f   %d %d %d   %d  %d  %d",
            &pointM0, &pointM12, &pointXsec,
            &nsel, &n_a, &n_d,
            &n_sb,
            &n_sl_sb,
            &n_sl_sig ) ;

          int nGenPerPoint(10000) ;
          float nselWeighted = ( nsel * pointXsec * dataLumi ) / ( 1.0*nGenPerPoint ) ;


          int m0bin  = hsusyscanExcluded->GetXaxis()->FindBin( pointM0 ) ;
          int m12bin = hsusyscanExcluded->GetYaxis()->FindBin( pointM12 ) ;

    //    printf(" m0 = %4.0f (%d),  m1/2 = %4.0f (%d),  Npred = %7.1f", pointM0, m0bin, pointM12, m12bin, nselWeighted ) ;

          if ( nselWeighted > susySigHigh ) {
    //       printf(" Excluded\n") ;
             hsusyscanExcluded->SetBinContent( m0bin, m12bin, 1. ) ;
          } else {
    //       printf("\n") ;
          }


       } // pi .

       fclose( infp ) ;

       gStyle->SetPadGridX(1) ;
       gStyle->SetPadGridY(1) ;
       TCanvas* csusy = new TCanvas("csusy","SUSY m1/2 vs m0 scan") ;
       hsusyscanExcluded->Draw("col") ;
       csusy->SaveAs("susyScan.png") ;
       TFile* f = new TFile("scan.root","recreate") ;
       hsusyscanExcluded->Write() ;
       f->Write() ;
       f->Close() ;

       return true ;

    } // susyScanNoContam.


  //===================================================================

    bool ra2bRoostatsClass2::susyScanWithContam( const char* inputScanFile, double dataLumi ) {


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

       TH2F* hsusyscanNsigul = new TH2F("hsusyscanNsigul", "SUSY m1/2 vs m0 parameter scan, upper limit on Nsig",
            nM0bins, minM0-deltaM0/2., maxM0+deltaM0/2.,
            nM12bins, minM12-deltaM12/2., maxM12+deltaM12/2. ) ;

       TH2F* hsusyscanNsigpred = new TH2F("hsusyscanNsigpred", "SUSY m1/2 vs m0 parameter scan, predicted Nsig",
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
          int    nsel ;
          int    n_a ;
          int    n_d ;
          int    n_sb ;
          int    n_sl_sb ;
          int    n_sl_sig ;

          fscanf( infp, "%f %f %f   %d %d %d   %d  %d  %d",
            &pointM0, &pointM12, &pointXsec,
            &nsel, &n_a, &n_d,
            &n_sb,
            &n_sl_sb,
            &n_sl_sig ) ;


          int nGenPerPoint(10000) ;


          int m0bin  = hsusyscanExcluded->GetXaxis()->FindBin( pointM0 ) ;
          int m12bin = hsusyscanExcluded->GetYaxis()->FindBin( pointM12 ) ;


       //--- Set up the likelihood to include the SUSY contributions to the non-SIG regions.

          float weight = ( pointXsec * dataLumi ) / ( 1.0*nGenPerPoint ) ;

          reinitialize() ;

          rv_mu_susymc_sig -> setVal( nsel * weight ) ;
          rv_mu_susymc_a  -> setVal( n_a * weight ) ;
          rv_mu_susymc_d  -> setVal( n_d * weight ) ;
          rv_mu_susymc_sb -> setVal( n_sb * weight ) ;
          rv_mu_sl_susymc_sb -> setVal( n_sl_sb * weight ) ;
          rv_mu_sl_susymc_sig -> setVal( n_sl_sig * weight ) ;

          parameterSnapshot() ;

          doFit() ;

          parameterSnapshot() ;

          float susySigLow, susySigHigh ;
          profileSusySig( susySigLow, susySigHigh, false ) ;


          float nselWeighted =  nsel * weight ;


          printf("  Upper limit on SUSY SIG yield : %6.1f\n\n", susySigHigh ) ;
          printf(" m0 = %4.0f (%d),  m1/2 = %4.0f (%d),  Npred = %7.1f", pointM0, m0bin, pointM12, m12bin, nselWeighted ) ;

          hsusyscanNsigul->SetBinContent( m0bin, m12bin, susySigHigh ) ;
          hsusyscanNsigpred->SetBinContent( m0bin, m12bin, nselWeighted ) ;

          if ( nselWeighted > susySigHigh ) {
             printf(" Excluded\n") ;
             hsusyscanExcluded->SetBinContent( m0bin, m12bin, 1. ) ;
          } else {
             printf("\n") ;
          }


       } // pi .

       fclose( infp ) ;

       TCanvas* csusy = new TCanvas("csusy","SUSY m1/2 vs m0 scan") ;
       gStyle->SetPadGridX(1) ;
       gStyle->SetPadGridY(1) ;
       hsusyscanExcluded->Draw("col") ;
       csusy->SaveAs("susyScan.png") ;
       TFile* f = new TFile("scan.root","recreate") ;
       hsusyscanExcluded->Write() ;
       hsusyscanNsigul->Write() ;
       hsusyscanNsigpred->Write() ;
       f->Write() ;
       f->Close() ;

       return true ;

    } // susyScanWithContam.

  //===================================================================

    bool ra2bRoostatsClass2::setSusyScanPoint( const char* inputScanFile, double dataLumi, double m0, double m12 ) {



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


       //--- read in the column headers line.
       char c(0) ;
       c = fgetc( infp ) ;
       c = 0 ;
       while ( c!=10  ) { c = fgetc( infp ) ; }

       bool found(false) ;

       //--- Loop over the scan points.
       for ( int pi = 0 ; pi < nScanPoints ; pi++ ) {

          float pointM0 ;
          float pointM12 ;
          float pointXsec ;
          int    nsel ;
          int    n_a ;
          int    n_d ;
          int    n_sb ;
          int    n_sl_sb ;
          int    n_sl_sig ;

          fscanf( infp, "%f %f %f   %d %d %d   %d  %d  %d",
            &pointM0, &pointM12, &pointXsec,
            &nsel, &n_a, &n_d,
            &n_sb,
            &n_sl_sb,
            &n_sl_sig ) ;


          if (    fabs( pointM0 - m0 ) <= deltaM0/2.
               && fabs( pointM12 - m12 ) <= deltaM12/2. ) {

             int nGenPerPoint(10000) ;
             float nselWeighted = ( nsel * pointXsec * dataLumi ) / ( 1.0*nGenPerPoint ) ;


             printf("\n\n Found point m0 = %4.0f,  m1/2 = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, nselWeighted ) ;

          //--- Set up the likelihood to include the SUSY contributions to the non-SIG regions.

             float weight = ( pointXsec * dataLumi ) / ( 1.0*nGenPerPoint ) ;

             rv_mu_susymc_sig    -> setVal( nsel * weight ) ;
             rv_mu_susymc_a      -> setVal( n_a * weight ) ;
             rv_mu_susymc_d      -> setVal( n_d * weight ) ;
             rv_mu_susymc_sb     -> setVal( n_sb * weight ) ;
             rv_mu_sl_susymc_sb  -> setVal( n_sl_sb * weight ) ;
             rv_mu_sl_susymc_sig -> setVal( n_sl_sig * weight ) ;

             found = true ;

             break ;

          } // point match?

       } // pi .

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

    bool ra2bRoostatsClass2::fitQualityPlot( bool doNorm, double hmax ) {

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

      int nbins(32) ;

      TH1F* hfitqual_data = new TH1F("hfitqual_data","RA2b likelihood fit results, data",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_ttbar = new TH1F("hfitqual_ttbar","RA2b likelihood fit results, ttbar",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_qcd = new TH1F("hfitqual_qcd","RA2b likelihood fit results, QCD",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_wj = new TH1F("hfitqual_wj","RA2b likelihood fit results, W+jets",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_znn = new TH1F("hfitqual_znn","RA2b likelihood fit results, Ztonunu",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_ewo = new TH1F("hfitqual_ewo","RA2b likelihood fit results, EW-other",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_susy = new TH1F("hfitqual_susy","RA2b likelihood fit results, SUSY",
                            nbins, 0.5, nbins+0.5 ) ;
      TH1F* hfitqual_gaus = new TH1F("hfitqual_gaus","RA2b likelihood fit results, Gaussian constraints",
                            nbins, 0.5, nbins+0.5 ) ;



      hfitqual_ttbar -> SetFillColor(kBlue-9) ;
      hfitqual_qcd   -> SetFillColor(2) ;
      hfitqual_wj    -> SetFillColor(kGreen) ;
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

      double ttbarVal ;
      double qcdVal ;
      double wjVal ;
      double znnVal ;
      double ewoVal ;
      double susyVal ;

      double ttbarValNorm ;
      double qcdValNorm ;
      double wjValNorm ;
      double znnValNorm ;
      double ewoValNorm ;
      double susyValNorm ;

      double eff_sf ;

      char   binLabel[1000] ;



      eff_sf = rv_eff_sf->getVal() ;


     //-- SIG -------------------------------------------------------

      sprintf( binLabel, "SIG" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_Nsig->getVal() ;
      if ( useSigBgVars ) {
         ttbarVal = rrv_mu_ttbar_sig->getVal() ;
      } else {
         ttbarVal = rfv_mu_ttbar_sig->getVal() ;
      }
      qcdVal   = rv_mu_qcd_sig->getVal() ;
      wjVal    = (rv_mu_wjmc_sig->getVal()) * (rv_lsf_wjmc->getVal()) ;
      znnVal   = (rv_mu_znnmc_sig->getVal()) * (rv_lsf_znnmc->getVal()) ;
      ewoVal   = (rv_mu_ewomc_sig->getVal()) * (rv_lsf_ewomc->getVal()) ;
      susyVal  = rv_mu_susy_sig->getVal() ;
      susyVal = eff_sf * susyVal ;


      dataErr = sqrt(dataVal) ;

      dataErrNorm = dataErr ;
      ttbarValNorm = ttbarVal ;
      qcdValNorm = qcdVal ;
      wjValNorm = wjVal ;
      znnValNorm = znnVal ;
      ewoValNorm = ewoVal ;
      susyValNorm = susyVal ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         ttbarValNorm = ttbarVal / dataVal ;
         qcdValNorm   = qcdVal / dataVal ;
         wjValNorm    = wjVal / dataVal ;
         znnValNorm    = znnVal / dataVal ;
         ewoValNorm    = ewoVal / dataVal ;
         susyValNorm  = susyVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_ttbar -> SetBinContent( binIndex, ttbarValNorm ) ;
      hfitqual_qcd   -> SetBinContent( binIndex, qcdValNorm ) ;
      hfitqual_wj    -> SetBinContent( binIndex, wjValNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, znnValNorm ) ;
      hfitqual_ewo    -> SetBinContent( binIndex, ewoValNorm ) ;
      hfitqual_susy  -> SetBinContent( binIndex, susyValNorm ) ;


      binIndex++ ;





     //-- SB  -------------------------------------------------------

      sprintf( binLabel, "SB" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_Nsb->getVal() ;
      if ( useSigBgVars ) {
         ttbarVal = rfv_mu_ttbar_sb->getVal() ;
      } else {
         ttbarVal = rrv_mu_ttbar_sb->getVal() ;
      }
      qcdVal   = rv_mu_qcd_sb->getVal() ;
      wjVal    = (rv_mu_wjmc_sb->getVal()) * (rv_lsf_wjmc->getVal()) ;
      znnVal   = (rv_mu_znnmc_sb->getVal()) * (rv_lsf_znnmc->getVal()) ;
      ewoVal   = (rv_mu_ewomc_sb->getVal()) * (rv_lsf_ewomc->getVal()) ;
      susyVal  = rv_mu_susy_sb->getVal() ;
      susyVal = eff_sf * susyVal ;


      dataErr = sqrt(dataVal) ;

      dataErrNorm = dataErr ;
      ttbarValNorm = ttbarVal ;
      qcdValNorm = qcdVal ;
      wjValNorm = wjVal ;
      znnValNorm = znnVal ;
      ewoValNorm = ewoVal ;
      susyValNorm = susyVal ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         ttbarValNorm = ttbarVal / dataVal ;
         qcdValNorm   = qcdVal / dataVal ;
         wjValNorm    = wjVal / dataVal ;
         znnValNorm    = znnVal / dataVal ;
         ewoValNorm    = ewoVal / dataVal ;
         susyValNorm  = susyVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_ttbar -> SetBinContent( binIndex, ttbarValNorm ) ;
      hfitqual_qcd   -> SetBinContent( binIndex, qcdValNorm ) ;
      hfitqual_wj    -> SetBinContent( binIndex, wjValNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, znnValNorm ) ;
      hfitqual_ewo    -> SetBinContent( binIndex, ewoValNorm ) ;
      hfitqual_susy  -> SetBinContent( binIndex, susyValNorm ) ;


      binIndex++ ;






     //-- SL, SIG  -------------------------------------------------------

      sprintf( binLabel, "SL SIG" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_Nslsig->getVal() ;
      ttbarVal = rv_mu_sl_ttbar_sig->getVal() ;
      wjVal    = (rv_mu_wjmc_slsig->getVal()) * (rv_lsf_wjmc->getVal()) ;
      znnVal   = (rv_mu_znnmc_slsig->getVal()) * (rv_lsf_znnmc->getVal()) ;
      ewoVal   = (rv_mu_ewomc_slsig->getVal()) * (rv_lsf_ewomc->getVal()) ;
      susyVal  = rv_mu_sl_susy_sig->getVal() ;
      susyVal = eff_sf * susyVal ;


      dataErr = sqrt(dataVal) ;

      dataErrNorm = dataErr ;
      ttbarValNorm = ttbarVal ;
      wjValNorm = wjVal ;
      znnValNorm = znnVal ;
      ewoValNorm = ewoVal ;
      susyValNorm = susyVal ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         ttbarValNorm = ttbarVal / dataVal ;
         wjValNorm    = wjVal / dataVal ;
         znnValNorm    = znnVal / dataVal ;
         ewoValNorm    = ewoVal / dataVal ;
         susyValNorm  = susyVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_ttbar -> SetBinContent( binIndex, ttbarValNorm ) ;
      hfitqual_wj    -> SetBinContent( binIndex, wjValNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, znnValNorm ) ;
      hfitqual_ewo    -> SetBinContent( binIndex, ewoValNorm ) ;
      hfitqual_susy  -> SetBinContent( binIndex, susyValNorm ) ;


      binIndex++ ;





     //-- SL, SB -------------------------------------------------------

      sprintf( binLabel, "SL SB" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_Nslsb->getVal() ;
      ttbarVal = rv_mu_sl_ttbar_sb->getVal() ;
      wjVal    = (rv_mu_wjmc_slsb->getVal()) * (rv_lsf_wjmc->getVal()) ;
      znnVal   = (rv_mu_znnmc_slsb->getVal()) * (rv_lsf_znnmc->getVal()) ;
      ewoVal   = (rv_mu_ewomc_slsb->getVal()) * (rv_lsf_ewomc->getVal()) ;
      susyVal  = rv_mu_sl_susy_sb->getVal() ;
      susyVal = eff_sf * susyVal ;


      dataErr = sqrt(dataVal) ;

      dataErrNorm = dataErr ;
      ttbarValNorm = ttbarVal ;
      wjValNorm = wjVal ;
      znnValNorm = znnVal ;
      ewoValNorm = ewoVal ;
      susyValNorm = susyVal ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         ttbarValNorm = ttbarVal / dataVal ;
         wjValNorm    = wjVal / dataVal ;
         znnValNorm    = znnVal / dataVal ;
         ewoValNorm    = ewoVal / dataVal ;
         susyValNorm  = susyVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_ttbar -> SetBinContent( binIndex, ttbarValNorm ) ;
      hfitqual_wj    -> SetBinContent( binIndex, wjValNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, znnValNorm ) ;
      hfitqual_ewo    -> SetBinContent( binIndex, ewoValNorm ) ;
      hfitqual_susy  -> SetBinContent( binIndex, susyValNorm ) ;


      binIndex++ ;




     //-- D -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "D" ) ;
      } else {
         sprintf( binLabel, "D (x0.1)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_Nd->getVal() ;
      ttbarVal = rv_mu_ttbar_d->getVal() ;
      qcdVal   = rv_mu_qcd_d->getVal() ;
      wjVal    = (rv_mu_wjmc_d->getVal()) * (rv_lsf_wjmc->getVal()) ;
      znnVal   = (rv_mu_znnmc_d->getVal()) * (rv_lsf_znnmc->getVal()) ;
      ewoVal   = (rv_mu_ewomc_d->getVal()) * (rv_lsf_ewomc->getVal()) ;
      susyVal  = rv_mu_susy_d->getVal() ;
      susyVal = eff_sf * susyVal ;


      dataErr = sqrt(dataVal) ;

      dataErrNorm = dataErr ;
      ttbarValNorm = ttbarVal ;
      qcdValNorm = qcdVal ;
      wjValNorm = wjVal ;
      znnValNorm = znnVal ;
      ewoValNorm = ewoVal ;
      susyValNorm = susyVal ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         ttbarValNorm = ttbarVal / dataVal ;
         qcdValNorm   = qcdVal / dataVal ;
         wjValNorm    = wjVal / dataVal ;
         znnValNorm    = znnVal / dataVal ;
         ewoValNorm    = ewoVal / dataVal ;
         susyValNorm  = susyVal / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*0.1 ) ;
         dataErrNorm = dataErrNorm*0.1 ;
         ttbarValNorm = ttbarValNorm*0.1 ;
         wjValNorm = wjValNorm*0.1 ;
         znnValNorm = znnValNorm*0.1 ;
         ewoValNorm = ewoValNorm*0.1 ;
         qcdValNorm = qcdValNorm*0.1 ;
         susyValNorm = susyValNorm*0.1 ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_ttbar -> SetBinContent( binIndex, ttbarValNorm ) ;
      hfitqual_qcd   -> SetBinContent( binIndex, qcdValNorm ) ;
      hfitqual_wj    -> SetBinContent( binIndex, wjValNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, znnValNorm ) ;
      hfitqual_ewo    -> SetBinContent( binIndex, ewoValNorm ) ;
      hfitqual_susy  -> SetBinContent( binIndex, susyValNorm ) ;


      binIndex++ ;



     //-- A -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "A" ) ;
      } else {
         sprintf( binLabel, "A (x0.1)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_Na->getVal() ;
      ttbarVal = rv_mu_ttbar_a->getVal() ;
      qcdVal   = rv_mu_qcd_a->getVal() ;
      wjVal    = (rv_mu_wjmc_a->getVal()) * (rv_lsf_wjmc->getVal()) ;
      znnVal   = (rv_mu_znnmc_a->getVal()) * (rv_lsf_znnmc->getVal()) ;
      ewoVal   = (rv_mu_ewomc_a->getVal()) * (rv_lsf_ewomc->getVal()) ;
      susyVal  = rv_mu_susy_a->getVal() ;
      susyVal = eff_sf * susyVal ;


      dataErr = sqrt(dataVal) ;

      dataErrNorm = dataErr ;
      ttbarValNorm = ttbarVal ;
      qcdValNorm = qcdVal ;
      wjValNorm = wjVal ;
      znnValNorm = znnVal ;
      ewoValNorm = ewoVal ;
      susyValNorm = susyVal ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         ttbarValNorm = ttbarVal / dataVal ;
         qcdValNorm   = qcdVal / dataVal ;
         wjValNorm    = wjVal / dataVal ;
         znnValNorm    = znnVal / dataVal ;
         ewoValNorm    = ewoVal / dataVal ;
         susyValNorm  = susyVal / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*0.1 ) ;
         dataErrNorm = dataErrNorm*0.1 ;
         ttbarValNorm = ttbarValNorm*0.1 ;
         wjValNorm = wjValNorm*0.1 ;
         znnValNorm = znnValNorm*0.1 ;
         ewoValNorm = ewoValNorm*0.1 ;
         qcdValNorm = qcdValNorm*0.1 ;
         susyValNorm = susyValNorm*0.1 ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_ttbar -> SetBinContent( binIndex, ttbarValNorm ) ;
      hfitqual_qcd   -> SetBinContent( binIndex, qcdValNorm ) ;
      hfitqual_wj    -> SetBinContent( binIndex, wjValNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, znnValNorm ) ;
      hfitqual_ewo    -> SetBinContent( binIndex, ewoValNorm ) ;
      hfitqual_susy  -> SetBinContent( binIndex, susyValNorm ) ;


      binIndex++ ;















      double val, valNorm ;


      binIndex++ ;

    //++++++++ W+jets MC ++++++++++++++++++++++++++++++++++++++++

     //-- W+jets, SIG -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "Wjets, SIG" ) ;
      } else {
         sprintf( binLabel, "Wjets, SIG (x10)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_NWJmcsig->getVal() ;
      val      = (rv_mu_wjmc_sig->getVal())*(rv_lsf_wjmc->getVal()) ;

      dataErr = sqrt(dataVal) ;

      dataVal = dataVal*(rv_lsf_wjmc->getVal()) ;
      dataErr = dataErr*(rv_lsf_wjmc->getVal()) ;

      dataErrNorm = dataErr ;
      valNorm = val ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         valNorm = val / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*10. ) ;
         dataErrNorm = dataErrNorm*10. ;
         valNorm = valNorm*10. ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;
      hfitqual_wj    -> SetBinContent( binIndex, valNorm ) ;

      printf(" %10s : err=%4.2f ;   val = %4.2f\n", binLabel, dataErrNorm, valNorm ) ;

      binIndex++ ;



     //-- W+jets, SB -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "Wjets, SB" ) ;
      } else {
         sprintf( binLabel, "Wjets, SB (x10)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_NWJmcsb->getVal() ;
      val      = (rv_mu_wjmc_sb->getVal())*(rv_lsf_wjmc->getVal()) ;

      dataErr = sqrt(dataVal) ;

      dataVal = dataVal*(rv_lsf_wjmc->getVal()) ;
      dataErr = dataErr*(rv_lsf_wjmc->getVal()) ;

      dataErrNorm = dataErr ;
      valNorm = val ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         valNorm = val / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*10. ) ;
         dataErrNorm = dataErrNorm*10. ;
         valNorm = valNorm*10. ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;
      hfitqual_wj    -> SetBinContent( binIndex, valNorm ) ;

      printf(" %10s : err=%4.2f ;   val = %4.2f\n", binLabel, dataErrNorm, valNorm ) ;

      binIndex++ ;



     //-- W+jets, SL, SIG -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "Wjets, SL,SIG" ) ;
      } else {
         sprintf( binLabel, "Wjets, SL,SIG (x10)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_NWJmcslsig->getVal() ;
      val      = (rv_mu_wjmc_slsig->getVal())*(rv_lsf_wjmc->getVal()) ;

      dataErr = sqrt(dataVal) ;

      dataVal = dataVal*(rv_lsf_wjmc->getVal()) ;
      dataErr = dataErr*(rv_lsf_wjmc->getVal()) ;

      dataErrNorm = dataErr ;
      valNorm = val ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         valNorm = val / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*10. ) ;
         dataErrNorm = dataErrNorm*10. ;
         valNorm = valNorm*10. ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;
      hfitqual_wj    -> SetBinContent( binIndex, valNorm ) ;

      printf(" %10s : err=%4.2f ;   val = %4.2f\n", binLabel, dataErrNorm, valNorm ) ;

      binIndex++ ;



     //-- W+jets, SL, SB -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "Wjets, SL,SB" ) ;
      } else {
         sprintf( binLabel, "Wjets, SL,SB (x10)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_NWJmcslsb->getVal() ;
      val      = (rv_mu_wjmc_slsb->getVal())*(rv_lsf_wjmc->getVal()) ;

      dataErr = sqrt(dataVal) ;

      dataVal = dataVal*(rv_lsf_wjmc->getVal()) ;
      dataErr = dataErr*(rv_lsf_wjmc->getVal()) ;

      dataErrNorm = dataErr ;
      valNorm = val ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         valNorm = val / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*10. ) ;
         dataErrNorm = dataErrNorm*10. ;
         valNorm = valNorm*10. ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;
      hfitqual_wj    -> SetBinContent( binIndex, valNorm ) ;

      printf(" %10s : err=%4.2f ;   val = %4.2f\n", binLabel, dataErrNorm, valNorm ) ;

      binIndex++ ;



     //-- W+jets, D -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "Wjets, D" ) ;
      } else {
         sprintf( binLabel, "Wjets, D (x10)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_NWJmcd->getVal() ;
      val      = (rv_mu_wjmc_d->getVal())*(rv_lsf_wjmc->getVal()) ;

      dataErr = sqrt(dataVal) ;

      dataVal = dataVal*(rv_lsf_wjmc->getVal()) ;
      dataErr = dataErr*(rv_lsf_wjmc->getVal()) ;

      dataErrNorm = dataErr ;
      valNorm = val ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         valNorm = val / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*10. ) ;
         dataErrNorm = dataErrNorm*10. ;
         valNorm = valNorm*10. ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;
      hfitqual_wj    -> SetBinContent( binIndex, valNorm ) ;

      printf(" %10s : err=%4.2f ;   val = %4.2f\n", binLabel, dataErrNorm, valNorm ) ;

      binIndex++ ;



     //-- W+jets, A -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "Wjets, A" ) ;
      } else {
         sprintf( binLabel, "Wjets, A (x10)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_NWJmca->getVal() ;
      val      = (rv_mu_wjmc_a->getVal())*(rv_lsf_wjmc->getVal()) ;

      dataErr = sqrt(dataVal) ;

      dataVal = dataVal*(rv_lsf_wjmc->getVal()) ;
      dataErr = dataErr*(rv_lsf_wjmc->getVal()) ;

      dataErrNorm = dataErr ;
      valNorm = val ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         valNorm = val / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*10. ) ;
         dataErrNorm = dataErrNorm*10. ;
         valNorm = valNorm*10. ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;
      hfitqual_wj    -> SetBinContent( binIndex, valNorm ) ;

      printf(" %10s : err=%4.2f ;   val = %4.2f\n", binLabel, dataErrNorm, valNorm ) ;

      binIndex++ ;













      binIndex++ ;

    //++++++++ Z invis MC ++++++++++++++++++++++++++++++++++++++++

     //-- Z invis, SIG -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "Zinvis, SIG" ) ;
      } else {
         sprintf( binLabel, "Zinvis, SIG (x10)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_NZnnmcsig->getVal() ;
      val      = (rv_mu_znnmc_sig->getVal())*(rv_lsf_znnmc->getVal()) ;

      dataErr = sqrt(dataVal) ;

      dataVal = dataVal*(rv_lsf_znnmc->getVal()) ;
      dataErr = dataErr*(rv_lsf_znnmc->getVal()) ;

      dataErrNorm = dataErr ;
      valNorm = val ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         valNorm = val / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*10. ) ;
         dataErrNorm = dataErrNorm*10. ;
         valNorm = valNorm*10. ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, valNorm ) ;

      printf(" %10s : err=%4.2f ;   val = %4.2f\n", binLabel, dataErrNorm, valNorm ) ;

      binIndex++ ;



     //-- Z invis, SB -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "Zinvis, SB" ) ;
      } else {
         sprintf( binLabel, "Zinvis, SB (x10)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_NZnnmcsb->getVal() ;
      val      = (rv_mu_znnmc_sb->getVal())*(rv_lsf_znnmc->getVal()) ;

      dataErr = sqrt(dataVal) ;

      dataVal = dataVal*(rv_lsf_znnmc->getVal()) ;
      dataErr = dataErr*(rv_lsf_znnmc->getVal()) ;

      dataErrNorm = dataErr ;
      valNorm = val ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         valNorm = val / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*10. ) ;
         dataErrNorm = dataErrNorm*10. ;
         valNorm = valNorm*10. ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, valNorm ) ;

      printf(" %10s : err=%4.2f ;   val = %4.2f\n", binLabel, dataErrNorm, valNorm ) ;

      binIndex++ ;



     //-- Z invis, SL, SIG -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "Zinvis, SL,SIG" ) ;
      } else {
         sprintf( binLabel, "Zinvis, SL,SIG (x10)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;


      dataVal  = rv_NZnnmcslsig->getVal() ;
      val      = (rv_mu_znnmc_slsig->getVal())*(rv_lsf_znnmc->getVal()) ;

      if ( dataVal > 0.11 ) {

         dataErr = sqrt(dataVal) ;

      dataVal = dataVal*(rv_lsf_znnmc->getVal()) ;
      dataErr = dataErr*(rv_lsf_znnmc->getVal()) ;


         dataErrNorm = dataErr ;
         valNorm = val ;
         hfitqual_data->SetBinContent( binIndex, dataVal ) ;
         if ( dataVal > 0. && doNorm ) { 
            hfitqual_data->SetBinContent( binIndex, 1. ) ;
            dataErrNorm  = dataErr / dataVal ; 
            valNorm = val / dataVal ;
         }
         if ( !doNorm ) {
            hfitqual_data->SetBinContent( binIndex, dataVal*10. ) ;
            dataErrNorm = dataErrNorm*10. ;
            valNorm = valNorm*10. ;
         }

         hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;
         hfitqual_znn    -> SetBinContent( binIndex, valNorm ) ;

      } else {

         hfitqual_data->SetBinContent( binIndex, 0. ) ;
         hfitqual_znn    -> SetBinContent( binIndex, 0. ) ;

      }

      printf(" %10s : err=%4.2f ;   val = %4.2f\n", binLabel, dataErrNorm, valNorm ) ;

      binIndex++ ;



     //-- Z invis, SL, SB -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "Zinvis, SL,SB" ) ;
      } else {
         sprintf( binLabel, "Zinvis, SL,SB (x10)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_NZnnmcslsb->getVal() ;
      val      = (rv_mu_znnmc_slsb->getVal())*(rv_lsf_znnmc->getVal()) ;

      if ( dataVal > 0.11 ) {

         dataErr = sqrt(dataVal) ;

      dataVal = dataVal*(rv_lsf_znnmc->getVal()) ;
      dataErr = dataErr*(rv_lsf_znnmc->getVal()) ;

         dataErrNorm = dataErr ;
         valNorm = val ;
         hfitqual_data->SetBinContent( binIndex, dataVal ) ;
         if ( dataVal > 0. && doNorm ) { 
            hfitqual_data->SetBinContent( binIndex, 1. ) ;
            dataErrNorm  = dataErr / dataVal ; 
            valNorm = val / dataVal ;
         }
         if ( !doNorm ) {
            hfitqual_data->SetBinContent( binIndex, dataVal*10. ) ;
            dataErrNorm = dataErrNorm*10. ;
            valNorm = valNorm*10. ;
         }

         hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;
         hfitqual_znn    -> SetBinContent( binIndex, valNorm ) ;

      } else {

         hfitqual_data->SetBinContent( binIndex, 0. ) ;
         hfitqual_znn    -> SetBinContent( binIndex, 0. ) ;

      }

      printf(" %10s : err=%4.2f ;   val = %4.2f\n", binLabel, dataErrNorm, valNorm ) ;

      binIndex++ ;



     //-- Z invis, D -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "Zinvis, D" ) ;
      } else {
         sprintf( binLabel, "Zinvis, D (x10)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_NZnnmcd->getVal() ;
      val      = (rv_mu_znnmc_d->getVal())*(rv_lsf_znnmc->getVal()) ;

      dataErr = sqrt(dataVal) ;

      dataVal = dataVal*(rv_lsf_znnmc->getVal()) ;
      dataErr = dataErr*(rv_lsf_znnmc->getVal()) ;

      dataErrNorm = dataErr ;
      valNorm = val ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         valNorm = val / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*10. ) ;
         dataErrNorm = dataErrNorm*10. ;
         valNorm = valNorm*10. ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, valNorm ) ;

      printf(" %10s : err=%4.2f ;   val = %4.2f\n", binLabel, dataErrNorm, valNorm ) ;

      binIndex++ ;



     //-- Z invis, A -------------------------------------------------------

      if ( doNorm ) {
         sprintf( binLabel, "Zinvis, A" ) ;
      } else {
         sprintf( binLabel, "Zinvis, A (x10)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal  = rv_NZnnmca->getVal() ;
      val      = (rv_mu_znnmc_a->getVal())*(rv_lsf_znnmc->getVal()) ;

      dataErr = sqrt(dataVal) ;

      dataVal = dataVal*(rv_lsf_znnmc->getVal()) ;
      dataErr = dataErr*(rv_lsf_znnmc->getVal()) ;


      dataErrNorm = dataErr ;
      valNorm = val ;
      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         valNorm = val / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*10. ) ;
         dataErrNorm = dataErrNorm*10. ;
         valNorm = valNorm*10. ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;
      hfitqual_znn    -> SetBinContent( binIndex, valNorm ) ;

      printf(" %10s : err=%4.2f ;   val = %4.2f\n", binLabel, dataErrNorm, valNorm ) ;

      binIndex++ ;


















      binIndex++ ;

      double mcVal ;
      double mcValNorm ;

    //-- Gaussian constraint parameters.


    //-- MC SIG

      sprintf( binLabel, "QCD MC SIG" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal = rv_Nqcdmcsig->getVal() ;
      mcVal   = rv_mu_qcdmc_sig->getVal() ;

      dataErr = Nqcdmcsigerr ;

      dataErrNorm = dataErr ;
      mcValNorm = mcVal ;

      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         mcValNorm = mcVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_gaus -> SetBinContent( binIndex, mcValNorm ) ;

      printf(" %10s : err=%4.2f ;   MC = %4.2f\n", binLabel, dataErrNorm, mcValNorm ) ;

      binIndex++ ;


    //-- MC SB

      sprintf( binLabel, "QCD MC SB" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal = rv_Nqcdmcsb->getVal() ;
      mcVal   = rv_mu_qcdmc_sb->getVal() ;

      dataErr = Nqcdmcsberr ;

      dataErrNorm = dataErr ;
      mcValNorm = mcVal ;

      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         mcValNorm = mcVal / dataVal ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_gaus -> SetBinContent( binIndex, mcValNorm ) ;

      printf(" %10s : err=%4.2f ;   MC = %4.2f\n", binLabel, dataErrNorm, mcValNorm ) ;

      binIndex++ ;



    //-- MC D

      if ( doNorm ) {
         sprintf( binLabel, "QCD MC D" ) ;
      } else {
         sprintf( binLabel, "QCD MC D (x0.1)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal = rv_Nqcdmcd->getVal() ;
      mcVal   = rv_mu_qcdmc_d->getVal() ;

      dataErr = Nqcdmcderr ;

      dataErrNorm = dataErr ;
      mcValNorm = mcVal ;

      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         mcValNorm = mcVal / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*0.1 ) ;
         dataErrNorm = dataErrNorm*0.1 ;
         mcValNorm = mcValNorm*0.1 ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_gaus -> SetBinContent( binIndex, mcValNorm ) ;

      printf(" %10s : err=%4.2f ;   MC = %4.2f\n", binLabel, dataErrNorm, mcValNorm ) ;

      binIndex++ ;





    //-- MC A

      if ( doNorm ) {
         sprintf( binLabel, "QCD MC A" ) ;
      } else {
         sprintf( binLabel, "QCD MC A (x0.1)" ) ;
      }
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal = rv_Nqcdmca->getVal() ;
      mcVal   = rv_mu_qcdmc_a->getVal() ;

      dataErr = Nqcdmcaerr ;

      dataErrNorm = dataErr ;
      mcValNorm = mcVal ;

      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         mcValNorm = mcVal / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*0.1 ) ;
         dataErrNorm = dataErrNorm*0.1 ;
         mcValNorm = mcValNorm*0.1 ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_gaus -> SetBinContent( binIndex, mcValNorm ) ;

      printf(" %10s : err=%4.2f ;   MC = %4.2f\n", binLabel, dataErrNorm, mcValNorm ) ;

      binIndex++ ;






      binIndex++ ;

    //-- lsf, W+jets

      sprintf( binLabel, "LSF, W+jets" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal = lsf_WJmc ;
      mcVal   = rv_lsf_wjmc->getVal() ;

      dataErr = lsf_WJmc_err ;

      dataErrNorm = dataErr ;
      mcValNorm = mcVal ;

      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         mcValNorm = mcVal / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*0.1 ) ;
         dataErrNorm = dataErrNorm*0.1 ;
         mcValNorm = mcValNorm*0.1 ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_gaus -> SetBinContent( binIndex, mcValNorm ) ;

      printf(" %10s : err=%4.2f ;   MC = %4.2f\n", binLabel, dataErrNorm, mcValNorm ) ;

      binIndex++ ;




    //-- lsf, Z to nunu

      sprintf( binLabel, "LSF, Z invis" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal = lsf_Znnmc ;
      mcVal   = rv_lsf_znnmc->getVal() ;

      dataErr = lsf_Znnmc_err ;

      dataErrNorm = dataErr ;
      mcValNorm = mcVal ;

      hfitqual_data->SetBinContent( binIndex, dataVal ) ;
      if ( dataVal > 0. && doNorm ) { 
         hfitqual_data->SetBinContent( binIndex, 1. ) ;
         dataErrNorm  = dataErr / dataVal ; 
         mcValNorm = mcVal / dataVal ;
      }
      if ( !doNorm ) {
         hfitqual_data->SetBinContent( binIndex, dataVal*0.1 ) ;
         dataErrNorm = dataErrNorm*0.1 ;
         mcValNorm = mcValNorm*0.1 ;
      }


      hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

      hfitqual_gaus -> SetBinContent( binIndex, mcValNorm ) ;

      printf(" %10s : err=%4.2f ;   MC = %4.2f\n", binLabel, dataErrNorm, mcValNorm ) ;

      binIndex++ ;




//  //-- lsf, EW other

//    sprintf( binLabel, "LSF, W+jets" ) ;
//    xaxis->SetBinLabel(binIndex, binLabel ) ;

//    dataVal = lsf_WJmc ;
//    mcVal   = rv_lsf_wjmc->getVal() ;

//    dataErr = lsf_WJmc_err ;

//    dataErrNorm = dataErr ;
//    mcValNorm = mcVal ;

//    hfitqual_data->SetBinContent( binIndex, dataVal ) ;
//    if ( dataVal > 0. && doNorm ) { 
//       hfitqual_data->SetBinContent( binIndex, 1. ) ;
//       dataErrNorm  = dataErr / dataVal ; 
//       mcValNorm = mcVal / dataVal ;
//    }
//    if ( !doNorm ) {
//       hfitqual_data->SetBinContent( binIndex, dataVal*0.1 ) ;
//       dataErrNorm = dataErrNorm*0.1 ;
//       mcValNorm = mcValNorm*0.1 ;
//    }


//    hfitqual_data->SetBinError( binIndex, dataErrNorm ) ;

//    hfitqual_gaus -> SetBinContent( binIndex, mcValNorm ) ;

//    printf(" %10s : err=%4.2f ;   MC = %4.2f\n", binLabel, dataErrNorm, mcValNorm ) ;

//    binIndex++ ;










      binIndex++ ;

    //-- Eff SF

      double effsfScale(1.5) ;

      printf(" hist max: %g\n", hfitqual_data->GetMaximum()) ;
      double effsfsf = hmax*(hfitqual_data->GetMaximum())/(effsfScale)  ;

      sprintf( binLabel, "Eff SF (right scale)" ) ;
      xaxis->SetBinLabel(binIndex, binLabel ) ;

      dataVal = EffScaleFactor ;
      mcVal   = rv_eff_sf->getVal() ;

      dataErr = EffScaleFactorErr ;

      dataErrNorm = dataErr ;
      mcValNorm = mcVal ;


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

      hfitqual_fit->Add( hfitqual_wj ) ;
      hfitqual_fit->Add( hfitqual_znn ) ;
      hfitqual_fit->Add( hfitqual_qcd ) ;
      hfitqual_fit->Add( hfitqual_ttbar ) ;
      hfitqual_fit->Add( hfitqual_susy ) ;
      hfitqual_fit->Add( hfitqual_gaus ) ;

      TLegend* legend = new TLegend(0.88,0.3,0.98,0.9) ;

      legend->AddEntry( hfitqual_data,  "data" ) ;
      legend->AddEntry( hfitqual_susy,  "SUSY" ) ;
      legend->AddEntry( hfitqual_ttbar, "ttbar" ) ;
      legend->AddEntry( hfitqual_qcd,   "QCD" ) ;
      legend->AddEntry( hfitqual_znn,   "Znunu" ) ;
      legend->AddEntry( hfitqual_wj,    "W+jets" ) ;
      legend->AddEntry( hfitqual_gaus,  "Gaus.C." ) ;

      if ( doNorm ) {
         hfitqual_data->SetMaximum( hmax ) ;
      } else {
         hfitqual_data->SetMaximum( hmax*(hfitqual_data->GetMaximum()) ) ;
      }


      double oldBM = gStyle->GetPadBottomMargin() ;
      double oldRM = gStyle->GetPadRightMargin() ;


      gStyle->SetPadBottomMargin(0.3) ;
      gStyle->SetPadRightMargin(0.17) ;


      TCanvas* cfitqual = new TCanvas("cfitqual","RA2b fit quality", 1200, 550 ) ;

      gPad->SetTicks(1,0) ;


      hfitqual_data->SetLabelSize(0.05,"x") ;
      hfitqual_data->GetXaxis()->LabelsOption("v") ;

      hfitqual_data->Draw() ;
      hfitqual_fit->Draw("same") ;
      hfitqual_data->Draw("same") ;
      legend->Draw() ;



      TGaxis* axis = new TGaxis() ;
      axis->DrawAxis( nbins+0.5, 0, nbins+0.5, hfitqual_data->GetMaximum(), 0., effsfScale, 510, "+L") ;

      TLine* line = new TLine() ;
      line->SetLineWidth(2) ;

      line->DrawLine( nbins-2., 0., nbins-2., hfitqual_data->GetMaximum() ) ;


      cfitqual->SaveAs("fitqual.png") ;


      gStyle->SetPadBottomMargin(oldBM) ;
      gStyle->SetPadRightMargin(oldRM) ;

      return true ;

    } // fitQualityPlot

  //===================================================================


    void ra2bRoostatsClass2::setAndFixSusySig( double setVal ) {

       rv_mu_susy_sig->setVal( setVal ) ;
       rv_mu_susy_sig->setConstant( kTRUE ) ;

       printf("\n\n Set SUSY SIG yield to %.1f and fixed it.\n\n\n", setVal ) ;

    }

  //===================================================================


    void ra2bRoostatsClass2::freeSusySig( ) {

       rv_mu_susy_sig->setConstant( kFALSE ) ;

       printf("\n\n Freed SUSY SIG yield fit parameter.\n\n\n" ) ;

    }

  //===================================================================



    void ra2bRoostatsClass2::parameterSnapshot() {

       printf("\n\n====================================================================\n\n") ;

       printf("       Snapshot of observables and likelihood parameters.\n\n") ;


       printf("\n\n") ;


       printf("                     Nobs      all    ttbar      qcd       EW     SUSY x Eff SF.\n") ;
       printf("---------------------------------------------------------------------------------------------\n") ;

       float ttbarSig ;
       float ttbarSb ;

       if ( useSigBgVars ) {
          ttbarSig = rrv_mu_ttbar_sig->getVal() ;
          ttbarSb  = rfv_mu_ttbar_sb->getVal() ;
       } else {
          ttbarSig = rfv_mu_ttbar_sig->getVal() ;
          ttbarSb  = rrv_mu_ttbar_sb->getVal() ;
       }

       printf("\n") ;
       printf("     Nsig        :  %5.0f   %7.1f  %7.1f  %7.1f  %7.1f  %7.1f x %4.2f = %7.1f \n",
              rv_Nsig->getVal(),
              rv_n_sig->getVal(),
              ttbarSig,
              rv_mu_qcd_sig->getVal(),
              rv_mu_ew_sig->getVal(),
              rv_mu_susy_sig->getVal(),
              rv_eff_sf->getVal(),
              (rv_mu_susy_sig->getVal())*(rv_eff_sf->getVal())
              ) ;




       printf("     Nsb         :  %5.0f   %7.1f  %7.1f  %7.1f  %7.1f  %7.1f x %4.2f = %7.1f\n",
              rv_Nsb->getVal(),
              rv_n_sb->getVal(),
              ttbarSb,
              rv_mu_qcd_sb->getVal(),
              rv_mu_ew_sb->getVal(),
              rv_mu_susy_sb->getVal(),
              rv_eff_sf->getVal(),
              (rv_eff_sf->getVal())*(rv_mu_susy_sb->getVal()) ) ;




       printf(" SL Nsb          :  %5.0f   %7.1f  %7.1f           %7.1f  %7.1f x %4.2f = %7.1f\n",
              rv_Nslsb->getVal(),
              rv_n_sl_sb->getVal(),
              rv_mu_sl_ttbar_sb->getVal(),
              rv_mu_sl_ew_sb->getVal(),
              rv_mu_sl_susy_sb->getVal(),
              rv_eff_sf->getVal(),
              (rv_eff_sf->getVal())*(rv_mu_sl_susy_sb->getVal()) ) ;



       printf(" SL Nsig         :  %5.0f   %7.1f  %7.1f           %7.1f  %7.1f x %4.2f = %7.1f\n",
              rv_Nslsig->getVal(),
              rv_n_sl_sig->getVal(),
              rv_mu_sl_ttbar_sig->getVal(),
              rv_mu_sl_ew_sig->getVal(),
              rv_mu_sl_susy_sig->getVal(),
              rv_eff_sf->getVal(),
              (rv_eff_sf->getVal())*(rv_mu_sl_susy_sig->getVal()) ) ;



       printf("     N_A         :  %5.0f   %7.1f  %7.1f  %7.1f  %7.1f  %7.1f x %4.2f = %7.1f\n",
              rv_Na->getVal(),
              rv_n_a->getVal(),
              rv_mu_ttbar_a->getVal(),
              rv_mu_qcd_a->getVal(),
              rv_mu_ew_a->getVal(),
              rv_mu_susy_a->getVal(),
              rv_eff_sf->getVal(),
              (rv_eff_sf->getVal())*(rv_mu_susy_a->getVal()) ) ;

       printf("     N_D         :  %5.0f   %7.1f  %7.1f  %7.1f  %7.1f  %7.1f x %4.2f = %7.1f\n",
              rv_Nd->getVal(),
              rv_n_d->getVal(),
              rv_mu_ttbar_d->getVal(),
              rv_mu_qcd_d->getVal(),
              rv_mu_ew_d->getVal(),
              rv_mu_susy_d->getVal(),
              rv_eff_sf->getVal(),
              (rv_eff_sf->getVal())*(rv_mu_susy_d->getVal()) ) ;



       printf("\n") ;
       printf("  QCD MC A       :  input = %7.1f +/- %5.1f ,   fit = %7.1f  \n",
              rv_Nqcdmca->getVal(),
              Nqcdmcaerr,
              rv_mu_qcdmc_a->getVal() ) ;

       printf("  QCD MC D       :  input = %7.1f +/- %5.1f ,   fit = %7.1f  \n",
              rv_Nqcdmcd->getVal(),
              Nqcdmcderr,
              rv_mu_qcdmc_d->getVal() ) ;

       printf("  QCD MC SB      :  input = %7.1f +/- %5.1f ,   fit = %7.1f  \n",
              rv_Nqcdmcsb->getVal(),
              Nqcdmcsberr,
              rv_mu_qcdmc_sb->getVal() ) ;

       printf("  QCD MC SIG     :  input = %7.1f +/- %5.1f ,   fit = %7.1f  \n",
              rv_Nqcdmcsig->getVal(),
              Nqcdmcsigerr,
              rv_mu_qcdmc_sig->getVal() ) ;


       printf("\n") ;
       printf("  Eff scale fac. :  input = %4.2f +/- %4.2f ,   fit = %4.2f  \n",
              EffScaleFactor,
              EffScaleFactorErr,
              rv_eff_sf->getVal() ) ;

       printf("\n") ;
       printf("  LSF W+jets     :  input = %4.2f +/- %4.2f ,   fit = %4.2f  \n",
              lsf_WJmc,
              lsf_WJmc_err,
              rv_lsf_wjmc->getVal() ) ;

       printf("  LSF Z->invis   :  input = %4.2f +/- %4.2f ,   fit = %4.2f  \n",
              lsf_Znnmc,
              lsf_Znnmc_err,
              rv_lsf_znnmc->getVal() ) ;


       printf("\n\n====================================================================\n\n") ;

    }

  //===================================================================
