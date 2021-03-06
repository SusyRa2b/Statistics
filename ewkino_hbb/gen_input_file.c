
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TSystem.h"
#include "histio.c"


   int  nbgcomps(2) ;
   char bgcompname[2][100] = { "tt", "znn" } ;
   TChain* bgcompchain[2] ;

   TChain* sigchain ;

   float dataIntLumiIPB(20000.) ;




   float met_bin_edges_1bins[2] = { 50., 10000. } ;

   float met_bin_edges_4bins[5] = { 30., 50., 100., 150., 10000. } ;

   float met_bin_edges[100] ;




  //================================================================================================

   void gen_input_file( const char* outfilename = "outputfiles/input-file.txt", int sigmass = 250, float sig_strength = 0., int bins_of_met = 4, float min_met = 50. ) {

      if ( bins_of_met == 1 ) {
         for ( int bi=0; bi<=bins_of_met; bi++ ) { met_bin_edges[bi] = met_bin_edges_1bins[bi] ; }
         met_bin_edges[0] = min_met ;
      } else if ( bins_of_met == 4 ) {
         for ( int bi=0; bi<=bins_of_met; bi++ ) { met_bin_edges[bi] = met_bin_edges_4bins[bi] ; }
      } else {
         printf("\n\n *** Don't know how to do %d met bins.\n\n", bins_of_met ) ;
      }

      char metvarname[100] ;
      sprintf( metvarname, "METsig" ) ;


      gDirectory -> Delete( "h*" ) ;




    //--- setup chains of reduced trees for the comps.

      printf("\n\n Setting up reduced tree chains.\n\n" ) ;

      for ( int si=0; si<nbgcomps; si++ ) { bgcompchain[si] = new TChain("reducedTree") ; }

      char rtdir[10000] = "/data/cms/hadronic-susy-bjets/hbb/reduced-trees-may23-2013" ;
      /////char rtdir[10000] = "/Users/owen/work/cms/hadronic-susy-bjets/hbb/reduced-trees-may23-2013" ;

      int compIndex(0) ;

      char pathandfile[10000] ;

     //--- ttbar, 1 and 2 lepton
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1_AODSIM_UCSB1606_v66-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2_AODSIM_UCSB1596_v66-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      compIndex++ ;

     //--- Znn
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1525_v66-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7C-v1_AODSIM_UCSB1607_v66-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1524_v66-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1594_v66-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1523_v66-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_UCSB1602_v66-slimskim.root", rtdir ) ;
      bgcompchain[compIndex] -> Add( pathandfile ) ;
      compIndex++ ;

      float signal_weight = 1. ;

      if ( sigmass == 250 ) {

     // //--- signal, 250, original MC (10k gen)
     //  sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.test_TChihh_250_v68-slimskim.root", rtdir ) ;
     //  sigchain = new TChain("reducedTree") ;
     //  sigchain -> Add( pathandfile ) ;
     //  signal_weight = 9.0e-6 ;

        //--- signal, 250, new MC (100k gen)
         sprintf( pathandfile, "/data/cms/hadronic-susy-bjets/hbb/reduced-trees-june14-2013/TChihh_250_Jim-slimskim.root" ) ;
         sigchain = new TChain("reducedTree") ;
         sigchain -> Add( pathandfile ) ;
         signal_weight = 9.0e-7 ;

      } else if ( sigmass == 400 ) {

     // //--- signal, 400, original MC (10k gen)
     //  sprintf( pathandfile, "%s/reducedTree.CSVM_PF2PATjets_JES0_JER0_PFMETTypeI_METunc0_PUunc0_BTagEff05_HLTEff0.test_TChihh_400_v68-slimskim.root", rtdir ) ;
     //  sigchain = new TChain("reducedTree") ;
     //  sigchain -> Add( pathandfile ) ;
     //  signal_weight = 1.0e-6 ;

        //--- signal, 400, original MC (100k gen)
         sprintf( pathandfile, "/data/cms/hadronic-susy-bjets/hbb/reduced-trees-june14-2013/TChihh_400_Jim-slimskim.root" ) ;
         sigchain = new TChain("reducedTree") ;
         sigchain -> Add( pathandfile ) ;
         signal_weight = 1.0e-7 ;

      }





     //--- Define cuts.

      printf("\n\n Setting up cuts.\n\n") ;

     //--- These are included in the skim definition of doSlimSkim.c.  Doesn't hurt to apply them here.

      char basiccuts[10000] ;
      sprintf( basiccuts, "cutPV==1&&passCleaning==1&&buggyEvent==0" ) ;

      char triggercuts[10000] ;
      sprintf( triggercuts, "(passMC_DiCentralPFJet30_PFMET80_BTagCSV07==1||passMC_DiPFJet80_DiPFJet30_BTagCSVd07d05==1||passMC_PFMET150==1)" ) ;

      char njetcuts[10000] ;
      sprintf( njetcuts, "njets20>=4&&njets20<=5" ) ;

      char skimcuts[10000] ;
      sprintf( skimcuts, "((%s)&&(%s)&&(%s))", basiccuts, triggercuts, njetcuts ) ;


     //--- These are beyond the skim selection.

      char masssigcuts[10000] ;
      sprintf( masssigcuts, "%s", "abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<20&&((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>100)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<140))" ) ;

      char masssbcuts[10000] ;
      sprintf( masssbcuts, "%s", "!(abs(higgsMbb1MassDiff-higgsMbb2MassDiff)<30&&((0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)>90)&&(0.5*(higgsMbb1MassDiff+higgsMbb2MassDiff)<150)))" ) ;

      char btag4cuts[10000] ;
      sprintf( btag4cuts, "%s", "CSVbest2>0.898&&CSVbest3>0.679&&CSVbest4>0.244" ) ;
      char btag3cuts[10000] ;
      sprintf( btag3cuts, "%s", "CSVbest2>0.898&&CSVbest3>0.679&&CSVbest4<0.244" ) ;
      char btag2cuts[10000] ;
      sprintf( btag2cuts, "%s", "CSVbest2>0.898&&CSVbest3<0.244" ) ;

      char leptonveto[10000] ;
      sprintf( leptonveto, "%s", "nMuons==0&&nElectrons==0&&nIsoTracks15_005_03==0&&nIsoTracks5_005_03<2&&nTausLoose==0" ) ;

      char drmaxcut[10000] ;
      sprintf( drmaxcut, "%s", "deltaRmax_hh<2.2" ) ;







     //--- prepare histograms.

      printf("\n\n Booking histograms.\n\n") ;

      TH1F* h_4b_msig_bg[10] ;
      TH1F* h_3b_msig_bg[10] ;
      TH1F* h_2b_msig_bg[10] ;

      TH1F* h_4b_msb_bg[10] ;
      TH1F* h_3b_msb_bg[10] ;
      TH1F* h_2b_msb_bg[10] ;

      char hname[1000] ;
      char htitle[1000] ;

      for ( int si=0; si<nbgcomps; si++ ) {

         sprintf( hname, "h_%s_4b_msig_bg_%s", metvarname, bgcompname[si] ) ;
         sprintf( htitle, "%s, 4b, mass signal box, background, %s", metvarname, bgcompname[si] ) ;
         h_4b_msig_bg[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_4b_msig_bg[si] -> Sumw2() ;

         sprintf( hname, "h_%s_3b_msig_bg_%s", metvarname, bgcompname[si] ) ;
         sprintf( htitle, "%s, 3b, mass signal box, background, %s", metvarname, bgcompname[si] ) ;
         h_3b_msig_bg[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_3b_msig_bg[si] -> Sumw2() ;

         sprintf( hname, "h_%s_2b_msig_bg_%s", metvarname, bgcompname[si] ) ;
         sprintf( htitle, "%s, 2b, mass signal box, background, %s", metvarname, bgcompname[si] ) ;
         h_2b_msig_bg[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_2b_msig_bg[si] -> Sumw2() ;

         sprintf( hname, "h_%s_4b_msb_bg_%s", metvarname, bgcompname[si] ) ;
         sprintf( htitle, "%s, 4b, mass sideband, background, %s", metvarname, bgcompname[si] ) ;
         h_4b_msb_bg[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_4b_msb_bg[si] -> Sumw2() ;

         sprintf( hname, "h_%s_3b_msb_bg_%s", metvarname, bgcompname[si] ) ;
         sprintf( htitle, "%s, 3b, mass sideband, background, %s", metvarname, bgcompname[si] ) ;
         h_3b_msb_bg[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_3b_msb_bg[si] -> Sumw2() ;

         sprintf( hname, "h_%s_2b_msb_bg_%s", metvarname, bgcompname[si] ) ;
         sprintf( htitle, "%s, 2b, mass sideband, background, %s", metvarname, bgcompname[si] ) ;
         h_2b_msb_bg[si] = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
         h_2b_msb_bg[si] -> Sumw2() ;

      } // si.

      sprintf( hname, "h_%s_4b_msig_smc", metvarname ) ;
      sprintf( htitle, "%s, 4b, mass signal box, signal MC", metvarname ) ;
      TH1F* h_4b_msig_smc = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
      h_4b_msig_smc -> Sumw2() ;

      sprintf( hname, "h_%s_3b_msig_smc", metvarname ) ;
      sprintf( htitle, "%s, 3b, mass signal box, signal MC", metvarname ) ;
      TH1F* h_3b_msig_smc = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
      h_3b_msig_smc -> Sumw2() ;

      sprintf( hname, "h_%s_2b_msig_smc", metvarname ) ;
      sprintf( htitle, "%s, 2b, mass signal box, signal MC", metvarname ) ;
      TH1F* h_2b_msig_smc = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
      h_2b_msig_smc -> Sumw2() ;

      sprintf( hname, "h_%s_4b_msb_smc", metvarname ) ;
      sprintf( htitle, "%s, 4b, mass sideband, signal MC", metvarname ) ;
      TH1F* h_4b_msb_smc = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
      h_4b_msb_smc -> Sumw2() ;

      sprintf( hname, "h_%s_3b_msb_smc", metvarname ) ;
      sprintf( htitle, "%s, 3b, mass sideband, signal MC", metvarname ) ;
      TH1F* h_3b_msb_smc = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
      h_3b_msb_smc -> Sumw2() ;

      sprintf( hname, "h_%s_2b_msb_smc", metvarname ) ;
      sprintf( htitle, "%s, 2b, mass sideband, signal MC", metvarname ) ;
      TH1F* h_2b_msb_smc = new TH1F( hname, htitle, bins_of_met, met_bin_edges ) ;
      h_2b_msb_smc -> Sumw2() ;







     //--- fill the histograms.

      TCanvas* can = new TCanvas("can","plots") ;

      printf("\n\n Filling the histograms.\n\n" ) ; fflush(stdout) ;

      char arg1[1000] ;
      char allcuts[10000] ;

      for ( int si=0; si<nbgcomps; si++ ) {

         printf("\n\n ++++++++ %s component\n\n", bgcompname[si] ) ;

         sprintf( arg1, "%s>>h_%s_4b_msig_bg_%s", metvarname, metvarname, bgcompname[si] ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s)&&(%s)&&(%s))*weight3*%.0f", skimcuts, leptonveto, drmaxcut, masssigcuts, btag4cuts, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         bgcompchain[si] -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_4b_msig_bg[si] -> Print("all") ;

         sprintf( arg1, "%s>>h_%s_3b_msig_bg_%s", metvarname, metvarname, bgcompname[si] ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s)&&(%s)&&(%s))*weight3*%.0f", skimcuts, leptonveto, drmaxcut, masssigcuts, btag3cuts, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         bgcompchain[si] -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_3b_msig_bg[si] -> Print("all") ;

         sprintf( arg1, "%s>>h_%s_2b_msig_bg_%s", metvarname, metvarname, bgcompname[si] ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s)&&(%s)&&(%s))*weight3*%.0f", skimcuts, leptonveto, drmaxcut, masssigcuts, btag2cuts, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         bgcompchain[si] -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_2b_msig_bg[si] -> Print("all") ;


         sprintf( arg1, "%s>>h_%s_4b_msb_bg_%s", metvarname, metvarname, bgcompname[si] ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s)&&(%s)&&(%s))*weight3*%.0f", skimcuts, leptonveto, drmaxcut, masssbcuts, btag4cuts, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         bgcompchain[si] -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_4b_msb_bg[si] -> Print("all") ;

         sprintf( arg1, "%s>>h_%s_3b_msb_bg_%s", metvarname, metvarname, bgcompname[si] ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s)&&(%s)&&(%s))*weight3*%.0f", skimcuts, leptonveto, drmaxcut, masssbcuts, btag3cuts, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         bgcompchain[si] -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_3b_msb_bg[si] -> Print("all") ;

         sprintf( arg1, "%s>>h_%s_2b_msb_bg_%s", metvarname, metvarname, bgcompname[si] ) ;
         sprintf( allcuts, "((%s)&&(%s)&&(%s)&&(%s)&&(%s))*weight3*%.0f", skimcuts, leptonveto, drmaxcut, masssbcuts, btag2cuts, dataIntLumiIPB ) ;
         printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
         bgcompchain[si] -> Draw( arg1, allcuts ) ;
         can->Update() ; can->Draw() ;
         h_2b_msb_bg[si] -> Print("all") ;


      } // si.

      printf("\n\n ++++++++ signal\n\n" ) ;

      printf("\n\n 4b, mass sig\n\n") ; fflush(stdout) ;
      sprintf( arg1, "%s>>h_%s_4b_msig_smc", metvarname, metvarname ) ;
      sprintf( allcuts, "((%s)&&(%s)&&(%s)&&(%s)&&(%s))*%g*%.0f", skimcuts, leptonveto, drmaxcut, masssigcuts, btag4cuts, signal_weight, dataIntLumiIPB ) ;
      printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
      sigchain -> Draw( arg1, allcuts ) ;
      can->Update() ; can->Draw() ;
      h_4b_msig_smc -> Print("all") ;

      printf("\n\n 3b, mass sig\n\n") ; fflush(stdout) ;
      sprintf( arg1, "%s>>h_%s_3b_msig_smc", metvarname, metvarname ) ;
      sprintf( allcuts, "((%s)&&(%s)&&(%s)&&(%s)&&(%s))*%g*%.0f", skimcuts, leptonveto, drmaxcut, masssigcuts, btag3cuts, signal_weight, dataIntLumiIPB ) ;
      printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
      sigchain -> Draw( arg1, allcuts ) ;
      can->Update() ; can->Draw() ;
      h_3b_msig_smc -> Print("all") ;

      printf("\n\n 2b, mass sig\n\n") ; fflush(stdout) ;
      sprintf( arg1, "%s>>h_%s_2b_msig_smc", metvarname, metvarname ) ;
      sprintf( allcuts, "((%s)&&(%s)&&(%s)&&(%s)&&(%s))*%g*%.0f", skimcuts, leptonveto, drmaxcut, masssigcuts, btag2cuts, signal_weight, dataIntLumiIPB ) ;
      printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
      sigchain -> Draw( arg1, allcuts ) ;
      can->Update() ; can->Draw() ;
      h_2b_msig_smc -> Print("all") ;



      printf("\n\n 4b, mass sb\n\n") ; fflush(stdout) ;
      sprintf( arg1, "%s>>h_%s_4b_msb_smc", metvarname, metvarname ) ;
      sprintf( allcuts, "((%s)&&(%s)&&(%s)&&(%s)&&(%s))*%g*%.0f", skimcuts, leptonveto, drmaxcut, masssbcuts, btag4cuts, signal_weight, dataIntLumiIPB ) ;
      printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
      sigchain -> Draw( arg1, allcuts ) ;
      can->Update() ; can->Draw() ;
      h_4b_msb_smc -> Print("all") ;

      printf("\n\n 3b, mass sb\n\n") ; fflush(stdout) ;
      sprintf( arg1, "%s>>h_%s_3b_msb_smc", metvarname, metvarname ) ;
      sprintf( allcuts, "((%s)&&(%s)&&(%s)&&(%s)&&(%s))*%g*%.0f", skimcuts, leptonveto, drmaxcut, masssbcuts, btag3cuts, signal_weight, dataIntLumiIPB ) ;
      printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
      sigchain -> Draw( arg1, allcuts ) ;
      can->Update() ; can->Draw() ;
      h_3b_msb_smc -> Print("all") ;

      printf("\n\n 2b, mass sb\n\n") ; fflush(stdout) ;
      sprintf( arg1, "%s>>h_%s_2b_msb_smc", metvarname, metvarname ) ;
      sprintf( allcuts, "((%s)&&(%s)&&(%s)&&(%s)&&(%s))*%g*%.0f", skimcuts, leptonveto, drmaxcut, masssbcuts, btag2cuts, signal_weight, dataIntLumiIPB ) ;
      printf("\n %s : %s\n", arg1, allcuts ) ; fflush(stdout) ;
      sigchain -> Draw( arg1, allcuts ) ;
      can->Update() ; can->Draw() ;
      h_2b_msb_smc -> Print("all") ;





     //--- Print a nice table.


      float n_4b_msig[50] ;
      float n_3b_msig[50] ;
      float n_2b_msig[50] ;

      float n_4b_msb[50] ;
      float n_3b_msb[50] ;
      float n_2b_msb[50] ;

      float smc_4b_msig[50] ;
      float smc_3b_msig[50] ;
      float smc_2b_msig[50] ;

      float smc_4b_msb[50] ;
      float smc_3b_msb[50] ;
      float smc_2b_msb[50] ;

      float Rsigsb_4b_val[50] ;
      float Rsigsb_3b_val[50] ;
      float Rsigsb_2b_val[50] ;

      float Rsigsb_4b_err[50] ;
      float Rsigsb_3b_err[50] ;
      float Rsigsb_2b_err[50] ;

      float Rsigsb_4b_staterr2[50] ;
      float Rsigsb_3b_staterr2[50] ;
      float Rsigsb_2b_staterr2[50] ;


      printf("\n\n\n") ;
      printf("============================================================================================================================================================================\n") ;
      printf(" METsig    comp  |   4bSB   4bSIG     4bSIG/4bSB    |   3bSB   3bSIG     3bSIG/3bSB     |   2bSB   2bSIG     2bSIG/2bSB     |    BG val.        3b pred.       2b pred.     |\n") ;
      printf("============================================================================================================================================================================\n") ;
      fflush(stdout) ;
      for ( int hbi=1; hbi<=bins_of_met; hbi++ ) {

         float metsiglow  = h_4b_msig_smc->GetXaxis()->GetBinLowEdge( hbi ) ;
         float metsighigh = h_4b_msig_smc->GetXaxis()->GetBinLowEdge( hbi+1 ) ;
         if ( hbi==bins_of_met ) metsighigh = 999. ;

         float bgsum_nmsig_4b_val(0.) ; float bgsum_nmsig_4b_err2(0.) ;
         float bgsum_nmsig_3b_val(0.) ; float bgsum_nmsig_3b_err2(0.) ;
         float bgsum_nmsig_2b_val(0.) ; float bgsum_nmsig_2b_err2(0.) ;

         float bgsum_nmsb_4b_val(0.) ; float bgsum_nmsb_4b_err2(0.) ;
         float bgsum_nmsb_3b_val(0.) ; float bgsum_nmsb_3b_err2(0.) ;
         float bgsum_nmsb_2b_val(0.) ; float bgsum_nmsb_2b_err2(0.) ;

         for ( int si=0; si<nbgcomps; si++ ) {

            float nmsig_4b_val = h_4b_msig_bg[si]->GetBinContent(hbi) ; float nmsig_4b_err = h_4b_msig_bg[si]->GetBinError(hbi) ;
            float nmsig_3b_val = h_3b_msig_bg[si]->GetBinContent(hbi) ; float nmsig_3b_err = h_3b_msig_bg[si]->GetBinError(hbi) ;
            float nmsig_2b_val = h_2b_msig_bg[si]->GetBinContent(hbi) ; float nmsig_2b_err = h_2b_msig_bg[si]->GetBinError(hbi) ;

            float nmsb_4b_val = h_4b_msb_bg[si]->GetBinContent(hbi) ; float nmsb_4b_err = h_4b_msb_bg[si]->GetBinError(hbi) ;
            float nmsb_3b_val = h_3b_msb_bg[si]->GetBinContent(hbi) ; float nmsb_3b_err = h_3b_msb_bg[si]->GetBinError(hbi) ;
            float nmsb_2b_val = h_2b_msb_bg[si]->GetBinContent(hbi) ; float nmsb_2b_err = h_2b_msb_bg[si]->GetBinError(hbi) ;


            printf( "[%3.0f,%3.0f] %6s |  %6.1f %6.1f                   |  %6.1f %6.1f                    |  %6.1f %6.1f                    |                                               |\n",
               metsiglow, metsighigh,
               bgcompname[si],
               nmsb_4b_val, nmsig_4b_val,
               nmsb_3b_val, nmsig_3b_val,
               nmsb_2b_val, nmsig_2b_val
             ) ;

            if ( si < nbgcomps ) {
               bgsum_nmsig_4b_val += nmsig_4b_val ;
               bgsum_nmsig_3b_val += nmsig_3b_val ;
               bgsum_nmsig_2b_val += nmsig_2b_val ;
               bgsum_nmsig_4b_err2 += pow( nmsig_4b_err, 2.) ;
               bgsum_nmsig_3b_err2 += pow( nmsig_3b_err, 2.) ;
               bgsum_nmsig_2b_err2 += pow( nmsig_2b_err, 2.) ;
               bgsum_nmsb_4b_val += nmsb_4b_val ;
               bgsum_nmsb_3b_val += nmsb_3b_val ;
               bgsum_nmsb_2b_val += nmsb_2b_val ;
               bgsum_nmsb_4b_err2 += pow( nmsb_4b_err, 2.) ;
               bgsum_nmsb_3b_err2 += pow( nmsb_3b_err, 2.) ;
               bgsum_nmsb_2b_err2 += pow( nmsb_2b_err, 2.) ;
            }

            if ( si == (nbgcomps-1) ) {

               Rsigsb_4b_val[hbi] =  bgsum_nmsig_4b_val / bgsum_nmsb_4b_val ;
               Rsigsb_4b_err[hbi] =  Rsigsb_4b_val[hbi] * sqrt( bgsum_nmsig_4b_err2 / pow(bgsum_nmsig_4b_val,2) + bgsum_nmsb_4b_err2 / pow(bgsum_nmsb_4b_val,2) ) ;
               Rsigsb_3b_val[hbi] =  bgsum_nmsig_3b_val / bgsum_nmsb_3b_val ;
               Rsigsb_3b_err[hbi] =  Rsigsb_3b_val[hbi] * sqrt( bgsum_nmsig_3b_err2 / pow(bgsum_nmsig_3b_val,2) + bgsum_nmsb_3b_err2 / pow(bgsum_nmsb_3b_val,2) ) ;
               Rsigsb_2b_val[hbi] =  bgsum_nmsig_2b_val / bgsum_nmsb_2b_val ;
               Rsigsb_2b_err[hbi] =  Rsigsb_2b_val[hbi] * sqrt( bgsum_nmsig_2b_err2 / pow(bgsum_nmsig_2b_val,2) + bgsum_nmsb_2b_err2 / pow(bgsum_nmsb_2b_val,2) ) ;

               Rsigsb_4b_staterr2[hbi] = 0. ;
               if ( bgsum_nmsig_4b_val > 0. ) { Rsigsb_4b_staterr2[hbi] += Rsigsb_4b_val[hbi] * Rsigsb_4b_val[hbi] / bgsum_nmsig_4b_val ; }
               if ( bgsum_nmsb_4b_val  > 0. ) { Rsigsb_4b_staterr2[hbi] += Rsigsb_4b_val[hbi] * Rsigsb_4b_val[hbi] / bgsum_nmsb_4b_val ; }
               if ( Rsigsb_4b_staterr2[hbi] < 0.0001 ) { Rsigsb_4b_staterr2[hbi] = 0.15 ; }

               Rsigsb_3b_staterr2[hbi] = 0. ;
               if ( bgsum_nmsig_3b_val > 0. ) { Rsigsb_3b_staterr2[hbi] += Rsigsb_3b_val[hbi] * Rsigsb_3b_val[hbi] / bgsum_nmsig_3b_val ; }
               if ( bgsum_nmsb_3b_val  > 0. ) { Rsigsb_3b_staterr2[hbi] += Rsigsb_3b_val[hbi] * Rsigsb_3b_val[hbi] / bgsum_nmsb_3b_val ; }
               if ( Rsigsb_3b_staterr2[hbi] < 0.0001 ) { Rsigsb_3b_staterr2[hbi] = 0.15 ; }

               Rsigsb_2b_staterr2[hbi] = 0. ;
               if ( bgsum_nmsig_2b_val > 0. ) { Rsigsb_2b_staterr2[hbi] += Rsigsb_2b_val[hbi] * Rsigsb_2b_val[hbi] / bgsum_nmsig_2b_val ; }
               if ( bgsum_nmsb_2b_val  > 0. ) { Rsigsb_2b_staterr2[hbi] += Rsigsb_2b_val[hbi] * Rsigsb_2b_val[hbi] / bgsum_nmsb_2b_val ; }
               if ( Rsigsb_2b_staterr2[hbi] < 0.0001 ) { Rsigsb_2b_staterr2[hbi] = 0.15 ; }

               float bgpred_from3b_val = bgsum_nmsb_4b_val * Rsigsb_3b_val[hbi] ;
               float bgpred_from3b_err = bgpred_from3b_val * sqrt( bgsum_nmsb_4b_err2/pow(bgsum_nmsb_4b_val,2.) + pow( Rsigsb_3b_err[hbi]/Rsigsb_3b_val[hbi],2.) ) ;

               float bgpred_from2b_val = bgsum_nmsb_4b_val * Rsigsb_2b_val[hbi] ;
               float bgpred_from2b_err = bgpred_from2b_val * sqrt( bgsum_nmsb_4b_err2/pow(bgsum_nmsb_4b_val,2.) + pow( Rsigsb_2b_err[hbi]/Rsigsb_2b_val[hbi],2.) ) ;

               printf( "[%3.0f,%3.0f] %6s |  %6.1f %6.1f   %5.3f +/- %5.3f |  %6.1f %6.1f   %5.3f +/- %5.3f  |  %6.1f %6.1f   %5.3f +/- %5.3f  | %4.1f +/- %4.1f   %4.1f +/- %4.1f  %4.1f +/- %4.1f  |\n",
                  metsiglow, metsighigh,
                  "bg sum",
                  bgsum_nmsb_4b_val, bgsum_nmsig_4b_val,   Rsigsb_4b_val[hbi], Rsigsb_4b_err[hbi],
                  bgsum_nmsb_3b_val, bgsum_nmsig_3b_val,   Rsigsb_3b_val[hbi], Rsigsb_3b_err[hbi],
                  bgsum_nmsb_2b_val, bgsum_nmsig_2b_val,   Rsigsb_2b_val[hbi], Rsigsb_2b_err[hbi],
                  bgsum_nmsig_4b_val, sqrt(bgsum_nmsig_4b_err2),
                  bgpred_from3b_val, bgpred_from3b_err,
                  bgpred_from2b_val, bgpred_from2b_err
                ) ;
                printf("----------------------------------------------------------------------------------------------------------------------------|                                               |\n") ;
            }

         } // si.

         float nmsig_4b_val_smc = h_4b_msig_smc->GetBinContent(hbi) ;
         float nmsig_3b_val_smc = h_3b_msig_smc->GetBinContent(hbi) ;
         float nmsig_2b_val_smc = h_2b_msig_smc->GetBinContent(hbi) ;

         float nmsb_4b_val_smc = h_4b_msb_smc->GetBinContent(hbi) ;
         float nmsb_3b_val_smc = h_3b_msb_smc->GetBinContent(hbi) ;
         float nmsb_2b_val_smc = h_2b_msb_smc->GetBinContent(hbi) ;


         printf( "[%3.0f,%3.0f] %6s |  %6.1f %6.1f                   |  %6.1f %6.1f                    |  %6.1f %6.1f                    |                                               |\n",
            metsiglow, metsighigh,
            "sig MC",
            nmsb_4b_val_smc, nmsig_4b_val_smc,
            nmsb_3b_val_smc, nmsig_3b_val_smc,
            nmsb_2b_val_smc, nmsig_2b_val_smc
          ) ;


         printf("============================================================================================================================================================================\n") ;

         n_4b_msig[hbi] = bgsum_nmsig_4b_val ;
         n_3b_msig[hbi] = bgsum_nmsig_3b_val ;
         n_2b_msig[hbi] = bgsum_nmsig_2b_val ;
         n_4b_msb[hbi]  = bgsum_nmsb_4b_val ;
         n_3b_msb[hbi]  = bgsum_nmsb_3b_val ;
         n_2b_msb[hbi]  = bgsum_nmsb_2b_val ;

         smc_4b_msig[hbi] =  nmsig_4b_val_smc ;
         smc_3b_msig[hbi] =  nmsig_3b_val_smc ;
         smc_2b_msig[hbi] =  nmsig_2b_val_smc ;
         smc_4b_msb[hbi]  =  nmsb_4b_val_smc ;
         smc_3b_msb[hbi]  =  nmsb_3b_val_smc ;
         smc_2b_msb[hbi]  =  nmsb_2b_val_smc ;

      } // hbi.

      printf("\n\n\n") ;



      //-- Compute SIG/SB corrections and uncertainties.
      //
      //  The fit will do a weighted average of the 2b, 3b, and 4b samples, so corrections
      //  should be made with respect to that average.  In computing the weighted average
      //  below, use sqrt(N) for uncertainty, not MC stats error, since that's what the
      //  fit will do.
      //
      //  Taking the uncertainty on the correction as the (MC stat error / value) on the difference
      //  with respect to the weighted ave combined with half of the difference with respect
      //  to the weighted ave in quadrature.
      //
      //  If bins_of_met is set to 1, only use 2b and 4b (ignore 3b).

      printf("\n\n") ;
      float wave_Rsigsb_val[50] ;
      for ( int hbi=1; hbi<=bins_of_met; hbi++ ) {
         if ( bins_of_met > 1 ) {
            wave_Rsigsb_val[hbi] =  ( Rsigsb_2b_val[hbi] / Rsigsb_2b_staterr2[hbi] + Rsigsb_3b_val[hbi] / Rsigsb_3b_staterr2[hbi] + Rsigsb_4b_val[hbi] / Rsigsb_4b_staterr2[hbi] ) /
                                    (                1.  / Rsigsb_2b_staterr2[hbi] +                1.  / Rsigsb_3b_staterr2[hbi] +                1.  / Rsigsb_4b_staterr2[hbi] ) ;
         } else {
            wave_Rsigsb_val[hbi] =  ( Rsigsb_2b_val[hbi] / Rsigsb_2b_staterr2[hbi] + Rsigsb_4b_val[hbi] / Rsigsb_4b_staterr2[hbi] ) /
                                    (                1.  / Rsigsb_2b_staterr2[hbi] +                1.  / Rsigsb_4b_staterr2[hbi] ) ;
         }
         printf(" expected sqrt(N) weighted Rsig/sb ave for METsig bin %d :  %5.3f\n", hbi, wave_Rsigsb_val[hbi] ) ;
      } // hbi

      float correction_Rsigsb_4b[50] ;
      float correction_Rsigsb_3b[50] ;
      float correction_Rsigsb_2b[50] ;

      float syst_Rsigsb_4b[50] ;
      float syst_Rsigsb_3b[50] ;
      float syst_Rsigsb_2b[50] ;

      printf("\n\n") ;

      for ( int hbi=1; hbi<=bins_of_met; hbi++ ) {

         correction_Rsigsb_4b[hbi] = Rsigsb_4b_val[hbi] / wave_Rsigsb_val[hbi] ;
         correction_Rsigsb_3b[hbi] = Rsigsb_3b_val[hbi] / wave_Rsigsb_val[hbi] ;
         correction_Rsigsb_2b[hbi] = Rsigsb_2b_val[hbi] / wave_Rsigsb_val[hbi] ;

         syst_Rsigsb_4b[hbi] = sqrt( pow( (Rsigsb_4b_err[hbi]/Rsigsb_4b_val[hbi]), 2. )  +  pow( (correction_Rsigsb_4b[hbi]-1.)/2.,2. ) ) ;
         syst_Rsigsb_3b[hbi] = sqrt( pow( (Rsigsb_3b_err[hbi]/Rsigsb_4b_val[hbi]), 2. )  +  pow( (correction_Rsigsb_3b[hbi]-1.)/2.,2. ) ) ;
         syst_Rsigsb_2b[hbi] = sqrt( pow( (Rsigsb_2b_err[hbi]/Rsigsb_4b_val[hbi]), 2. )  +  pow( (correction_Rsigsb_2b[hbi]-1.)/2.,2. ) ) ;

         if ( bins_of_met > 1 ) {
            printf("  Rsig/sb corrections for METsig bin %d :     4b = %5.3f +/- %5.3f,       3b =  %5.3f +/- %5.3f,      2b =  %5.3f +/- %5.3f\n",
                  hbi,
                  correction_Rsigsb_4b[hbi], syst_Rsigsb_4b[hbi],
                  correction_Rsigsb_3b[hbi], syst_Rsigsb_3b[hbi],
                  correction_Rsigsb_2b[hbi], syst_Rsigsb_2b[hbi] ) ;
         } else {
            printf("  Rsig/sb corrections for METsig bin %d :     4b = %5.3f +/- %5.3f,       2b =  %5.3f +/- %5.3f\n",
                  hbi,
                  correction_Rsigsb_4b[hbi], syst_Rsigsb_4b[hbi],
                  correction_Rsigsb_2b[hbi], syst_Rsigsb_2b[hbi] ) ;
         }

      } // hbi.





     //--- Generate the likelihood input file.

      gSystem -> Exec("mkdir -p outputfiles") ;

      printf("\n\n output file with path: %s\n\n", outfilename ) ; fflush(stdout) ;




      char command[10000] ;
      sprintf( command, "ls %s >& /dev/null", outfilename ) ;
      int returnstat = gSystem->Exec( command ) ;
      if ( returnstat == 0 ) {
         char mvfile[10000] ;
         sprintf( mvfile, "%s-old", outfilename ) ;
         printf("\n\n *** Output file already exists.  Moving it to %s\n\n", mvfile ) ;
         sprintf( command, "mv %s %s", outfilename, mvfile ) ;
         gSystem->Exec( command ) ;
      }
      FILE* outfile = fopen( outfilename, "w" ) ;

      fprintf( outfile, "met_variable  %s\n", metvarname ) ;
      fprintf( outfile, "bins_of_met  %d\n", bins_of_met ) ;

      //-- all observables, BG totals.  Added in signal if requested (sig_strength>0).
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "N_4b_msig_met%d   %8.2f\n", mbi, (n_4b_msig[mbi] + sig_strength * smc_4b_msig[mbi] )  ) ; }
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "N_4b_msb_met%d    %8.2f\n", mbi, (n_4b_msb[mbi]  + sig_strength * smc_4b_msb[mbi]  )  ) ; }
      if ( bins_of_met > 1 ) { for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "N_3b_msig_met%d   %8.2f\n", mbi, (n_3b_msig[mbi] + sig_strength * smc_3b_msig[mbi] )  ) ; } }
      if ( bins_of_met > 1 ) { for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "N_3b_msb_met%d    %8.2f\n", mbi, (n_3b_msb[mbi]  + sig_strength * smc_3b_msb[mbi]  )  ) ; } }
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "N_2b_msig_met%d   %8.2f\n", mbi, (n_2b_msig[mbi] + sig_strength * smc_2b_msig[mbi] )  ) ; }
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "N_2b_msb_met%d    %8.2f\n", mbi, (n_2b_msb[mbi]  + sig_strength * smc_2b_msb[mbi]  )  ) ; }


      //-- signal MC values.
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "smc_4b_msig_met%d   %8.2f\n", mbi, smc_4b_msig[mbi] ) ; }
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "smc_4b_msb_met%d    %8.2f\n", mbi, smc_4b_msb[mbi]  ) ; }
      if ( bins_of_met > 1 ) { for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "smc_3b_msig_met%d   %8.2f\n", mbi, smc_3b_msig[mbi] ) ; } }
      if ( bins_of_met > 1 ) { for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "smc_3b_msb_met%d    %8.2f\n", mbi, smc_3b_msb[mbi]  ) ; } }
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "smc_2b_msig_met%d   %8.2f\n", mbi, smc_2b_msig[mbi] ) ; }
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) { fprintf( outfile, "smc_2b_msb_met%d    %8.2f\n", mbi, smc_2b_msb[mbi]  ) ; }


      //-- Rsig/sb correction factors and systematics.
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {
         fprintf( outfile, "Rsigsb_corr_4b_met%d   %5.3f\n", mbi, correction_Rsigsb_4b[mbi] ) ;
         fprintf( outfile, "Rsigsb_syst_4b_met%d   %5.3f\n", mbi, syst_Rsigsb_4b[mbi] ) ;
      }
      if ( bins_of_met > 1 ) {
         for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {
            fprintf( outfile, "Rsigsb_corr_3b_met%d   %5.3f\n", mbi, correction_Rsigsb_3b[mbi] ) ;
            fprintf( outfile, "Rsigsb_syst_3b_met%d   %5.3f\n", mbi, syst_Rsigsb_3b[mbi] ) ;
         }
      }
      for ( int mbi=1; mbi<=bins_of_met; mbi++ ) {
         fprintf( outfile, "Rsigsb_corr_2b_met%d   %5.3f\n", mbi, correction_Rsigsb_2b[mbi] ) ;
         fprintf( outfile, "Rsigsb_syst_2b_met%d   %5.3f\n", mbi, syst_Rsigsb_2b[mbi] ) ;
      }



      fclose( outfile ) ;


      printf("\n\n Created likelihood input file: %s\n\n", outfilename ) ;




      printf("\n\n Saving histograms to outputfiles/gen_input.root\n\n") ;

      saveHist( "outputfiles/gen_input.root", "h*" ) ;




   } // gen_input_file







