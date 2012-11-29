
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>

  //-- met4-ht4-v15
      const int nBinsMET   = 4 ;
      const int nBinsHT    = 4 ;
      float Mbins[nBinsMET+1] = {125.,150.,250.,350.,99999.};
      float Hbins[nBinsHT+1] = {400.,500.,800.,1000.,99999.};

    int nComps(7) ;
    char compname[7][100] ;
    int compcolor[7] ;
    double complumi[7] ;
    TChain* compchain[7] ;
    double compscale[7] ;

    TCanvas* dcan ;

   //----------------------------
   // prototypes

    void bookSet( const char* hname_base, const char* htitle_base, int nbins, double xmin, double xmax ) ;
    void bookSetLHB( const char* hname_base, const char* htitle_base ) ;
    void fillSet( const char* hname_base, const char* varname, const char* cuts, bool doTrigCorr = true ) ;
    void fillSetLHB( const char* hname_base, const char* cuts, bool doTrigCorr = true ) ;
    void saveHist(const char* filename, const char* pat) ;

   //----------------------------

    void dmc_fillplots2( const char* dataset_string = "all" ) {

       gDirectory -> Delete("h*") ;

       gStyle -> SetOptStat(0) ;
       gStyle -> SetTitleH(0.035 ) ;

       for ( int ci=0; ci<nComps; ci++ ) { compchain[ci] = new TChain("tree") ; }

       int compind ;

    //-- Data
       compind = 0 ;
       sprintf( compname[compind], "data" ) ;
       if ( strcmp(dataset_string, "all" ) == 0 ) {
          printf("\n\n Loading all data.\n\n") ;
          compchain[compind] -> Add( "filesHCP_53_v6/MET_2012A.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/MET_2012B.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/MET_2012C_pr.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/MET_2012C_rr.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/HT_2012A.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/HTMHT_2012B.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/HTMHT_2012C_pr.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/HTMHT_2012C_rr.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/JetHT_2012B.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/JetHT_2012C_pr.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/JetHT_2012C_rr.root" ) ;
          complumi[compind] = 12.03 ;
       } else if ( strcmp(dataset_string, "RunsAB" ) == 0 ) {
          printf("\n\n Loading runs A and B.\n\n") ;
          compchain[compind] -> Add( "filesHCP_53_v6/MET_2012A.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/MET_2012B.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/HT_2012A.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/HTMHT_2012B.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/JetHT_2012B.root" ) ;
          complumi[compind] = 5.23 ;
       } else if ( strcmp(dataset_string, "RunC" ) == 0 ) {
          printf("\n\n Loading run C.\n\n") ;
          compchain[compind] -> Add( "filesHCP_53_v6/MET_2012C_pr.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/MET_2012C_rr.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/HTMHT_2012C_pr.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/HTMHT_2012C_rr.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/JetHT_2012C_pr.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v6/JetHT_2012C_rr.root" ) ;
          complumi[compind] = 6.81 ;
       } else {
          printf("\n\n *** Unknown dataset: %s\n\n", dataset_string ) ;
          return ;
       }
  //------------
       compcolor[compind] = 0 ;
       compscale[compind] = 1. ;


       double MClumi = 12.03 ;

    //-- Diboson.
       compind = 1 ;
       sprintf( compname[compind], "diboson" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/WW.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/WZ.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/ZZ.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 800+7 ; // kOrange = 800
       compscale[compind] = complumi[0] / complumi[compind] ;





    //-- Zinvisible
       compind = 2 ;
       sprintf( compname[compind], "znn" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/Zinv-100to200.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/Zinv-200to400.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/Zinv-400.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 416-3 ; // kGreen = 416
       compscale[compind] = complumi[0] / complumi[compind] ;




    //-- QCD
       compind = 3 ;
       sprintf( compname[compind], "qcd" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/QCD-120to170.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/QCD-170to300.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/QCD-300to470.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/QCD-470to600.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/QCD-600to800.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/QCD-800to1000.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/QCD-1000to1400.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/QCD-1400to1800.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/QCD-1800.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 2 ;
       compscale[compind] = complumi[0] / complumi[compind] ;




    //-- single top
       compind = 4 ;
       sprintf( compname[compind], "singlet" ) ;
       ////// compchain[compind] -> Add( "filesHCP_53_v6/T-t.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/Tbar-t.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/T-s.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/Tbar-s.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/T-tW.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/Tbar-tW.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 600-9 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;



    //-- W+jets
       compind = 5 ;
       sprintf( compname[compind], "wjets" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/WJets-250to300.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/WJets-300to400.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/WJets-400.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 600-4 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;



    //-- ttbar
       compind = 6 ;
       sprintf( compname[compind], "ttbar" ) ;
       ///////compchain[compind] -> Add( "filesHCP_53_v6/TT.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/TT_FullLept.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/TT_SemiLept.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/TT_FullHad.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 600-7 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;




    //--- cuts

       char basecuts_1lep[10000] ;
       sprintf( basecuts_1lep, "MT<100&&maxChNMultDiff<40&&pfOcaloMET<2&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nometratiocut[10000] ;
       sprintf( basecuts_1lep_nometratiocut, "MT<100&&maxChNMultDiff<40&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nonb_nometratiocut[10000] ;
       sprintf( basecuts_1lep_nonb_nometratiocut, "MT<100&&maxChNMultDiff<40&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nomindphin[10000] ;
       sprintf( basecuts_1lep_nomindphin, "MT<100&&maxChNMultDiff<40&&pfOcaloMET<2&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nonb[10000] ;
       sprintf( basecuts_1lep_nonb, "MT<100&&maxChNMultDiff<40&&pfOcaloMET<2&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nonjet[10000] ;
       sprintf( basecuts_1lep_nonjet, "MT<100&&maxChNMultDiff<40&&pfOcaloMET<2&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nB>0&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nonjet_nonb[10000] ;
       sprintf( basecuts_1lep_nonjet_nonb, "MT<100&&maxChNMultDiff<40&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nonb_noisotrk[10000] ;
       sprintf( basecuts_1lep_nonb_noisotrk, "MT<100&&maxChNMultDiff<40&&pfOcaloMET<2&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nonb_nomaxchmultdiff[10000] ;
       sprintf( basecuts_1lep_nonb_nomaxchmultdiff, "MT<100&&pfOcaloMET<2&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nomt[10000] ;
       sprintf( basecuts_1lep_nomt, "maxChNMultDiff<40&&pfOcaloMET<2&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nonb_nomt[10000] ;
       sprintf( basecuts_1lep_nonb_nomt, "maxChNMultDiff<40&&pfOcaloMET<2&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;




       char basecuts_ldp[10000] ;
       sprintf( basecuts_ldp, "maxChNMultDiff<40&&nIsoTrk==0&&pfOcaloMET<2&&minDelPhiN<=4&&(nMu==0&&nEl==0)&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_ldp_nometratiocut[10000] ;
       sprintf( basecuts_ldp_nometratiocut, "maxChNMultDiff<40&&nIsoTrk==0&&minDelPhiN<=4&&(nMu==0&&nEl==0)&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_ldp_nonb_nometratiocut[10000] ;
       sprintf( basecuts_ldp_nonb_nometratiocut, "maxChNMultDiff<40&&nIsoTrk==0&&minDelPhiN<=4&&(nMu==0&&nEl==0)&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_ldp_nonb[10000] ;
       sprintf( basecuts_ldp_nonb, "maxChNMultDiff<40&&nIsoTrk==0&&pfOcaloMET<2&&minDelPhiN<=4&&(nMu==0&&nEl==0)&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_ldp_nonjet[10000] ;
       sprintf( basecuts_ldp_nonjet, "maxChNMultDiff<40&&nIsoTrk==0&&pfOcaloMET<2&&minDelPhiN<=4&&(nMu==0&&nEl==0)&&nB>0&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_ldp_nonjet_nonb[10000] ;
       sprintf( basecuts_ldp_nonjet_nonb, "maxChNMultDiff<40&&nIsoTrk==0&&minDelPhiN<=4&&(nMu==0&&nEl==0)&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_ldp_nonb_noisotrk[10000] ;
       sprintf( basecuts_ldp_nonb_noisotrk, "maxChNMultDiff<40&&pfOcaloMET<2&&minDelPhiN<=4&&(nMu==0&&nEl==0)&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_ldp_nonb_nomaxchmultdiff[10000] ;
       sprintf( basecuts_ldp_nonb_nomaxchmultdiff, "nIsoTrk==0&&pfOcaloMET<2&&minDelPhiN<=4&&(nMu==0&&nEl==0)&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;




       char basecuts_0lep_nomindphin[10000] ;
       sprintf( basecuts_0lep_nomindphin, "maxChNMultDiff<40&&nIsoTrk==0&&pfOcaloMET<2&&(nMu==0&&nEl==0)&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_0lep_nomindphin_nonb[10000] ;
       sprintf( basecuts_0lep_nomindphin_nonb, "maxChNMultDiff<40&&nIsoTrk==0&&pfOcaloMET<2&&(nMu==0&&nEl==0)&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_0lep[10000] ;
       sprintf( basecuts_0lep, "maxChNMultDiff<40&&nIsoTrk==0&&pfOcaloMET<2&&minDelPhiN>4&&(nMu==0&&nEl==0)&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_0lep_nometratiocut[10000] ;
       sprintf( basecuts_0lep_nometratiocut, "maxChNMultDiff<40&&nIsoTrk==0&&minDelPhiN>4&&(nMu==0&&nEl==0)&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_0lep_nonb_nometratiocut[10000] ;
       sprintf( basecuts_0lep_nonb_nometratiocut, "maxChNMultDiff<40&&nIsoTrk==0&&minDelPhiN>4&&(nMu==0&&nEl==0)&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_0lep_nonb[10000] ;
       sprintf( basecuts_0lep_nonb, "maxChNMultDiff<40&&nIsoTrk==0&&pfOcaloMET<2&&minDelPhiN>4&&(nMu==0&&nEl==0)&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_0lep_nonb_noisotrk[10000] ;
       sprintf( basecuts_0lep_nonb_noisotrk, "maxChNMultDiff<40&&pfOcaloMET<2&&minDelPhiN>4&&(nMu==0&&nEl==0)&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_0lep_nonb_nomaxchmultdiff[10000] ;
       sprintf( basecuts_0lep_nonb_nomaxchmultdiff, "nIsoTrk==0&&pfOcaloMET<2&&minDelPhiN>4&&(nMu==0&&nEl==0)&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_0lep_nonjet[10000] ;
       sprintf( basecuts_0lep, "maxChNMultDiff<40&&nIsoTrk==0&&pfOcaloMET<2&&minDelPhiN>4&&(nMu==0&&nEl==0)&&nB>0&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_0lep_nonjet_nonb[10000] ;
       sprintf( basecuts_0lep, "maxChNMultDiff<40&&nIsoTrk==0&&pfOcaloMET<2&&minDelPhiN>4&&(nMu==0&&nEl==0)&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;



       dcan = (TCanvas*) gDirectory->FindObject("dcan") ;
       if ( dcan == 0x0 ) {
          dcan = new TCanvas("dcan","") ;
       }

       char htitle[1000] ;

    //-- MET

       sprintf( htitle, "MET, SL, %s", dataset_string ) ;
       bookSet( "h_met_sl_all", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_sl_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_sl_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, SL, HT<800, %s", dataset_string ) ;
       bookSet( "h_met_sl_htlt800_all", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, HT<800, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_htlt800_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, HT<800, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_sl_htlt800_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, HT<800, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_sl_htlt800_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, SL, HT>800, %s", dataset_string ) ;
       bookSet( "h_met_sl_htgt800_all", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, HT>800, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_htgt800_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, HT>800, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_sl_htgt800_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, HT>800, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_sl_htgt800_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, SL, no MC trig correction, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_all", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, SL, no MC trig correction, HT[400,500], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht1_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, HT[500,800], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht2_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, HT[800,1000], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht3_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, HT>1000, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht4_nb1", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, SL, no MC trig correction, HT[400,500], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht1_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, HT[500,800], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht2_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, HT[800,1000], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht3_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, HT>1000, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht4_nb1_wide1", htitle, 30, 125., 1500. ) ;

       sprintf( htitle, "MET, SL, HT[400,500], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht1_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, HT[500,800], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht2_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, HT[800,1000], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht3_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, HT>1000, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht4_nb1", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, SL, HT[400,500], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht1_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, SL, HT[500,800], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht2_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, SL, HT[800,1000], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht3_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, SL, HT>1000, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht4_nb1_wide1", htitle, 30, 125., 1500. ) ;

       sprintf( htitle, "MET, SL, no MC trig correction, HT[400,500], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht1_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, HT[500,800], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht2_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, HT[800,1000], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht3_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, HT>1000, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht4_nb0", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, SL, no MC trig correction, HT[400,500], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht1_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, HT[500,800], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht2_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, HT[800,1000], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht3_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, HT>1000, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_ht4_nb0_wide1", htitle, 30, 125., 1500. ) ;

       sprintf( htitle, "MET, SL, HT[400,500], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht1_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, HT[500,800], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht2_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, HT[800,1000], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht3_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, HT>1000, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht4_nb0", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, SL, HT[400,500], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht1_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, SL, HT[500,800], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht2_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, SL, HT[800,1000], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht3_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, SL, HT>1000, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_sl_ht4_nb0_wide1", htitle, 30, 125., 1500. ) ;


       sprintf( htitle, "MET, LDP, %s", dataset_string ) ;
       bookSet( "h_met_ldp_all", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, LDP, HT<800, %s", dataset_string ) ;
       bookSet( "h_met_ldp_htlt800_all", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, HT<800, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_htlt800_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, HT<800, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_ldp_htlt800_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, HT<800, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_ldp_htlt800_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, LDP, HT>800, %s", dataset_string ) ;
       bookSet( "h_met_ldp_htgt800_all", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, HT>800, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_htgt800_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, HT>800, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_ldp_htgt800_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, HT>800, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_ldp_htgt800_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, LDP, no MC trig correction, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_all", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, LDP, no MC trig correction, HT[400,500], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht1_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, HT[500,800], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht2_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, HT[800,1000], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht3_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, HT>1000, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht4_nb1", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, LDP, no MC trig correction, HT[400,500], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht1_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, HT[500,800], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht2_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, HT[800,1000], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht3_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, HT>1000, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht4_nb1_wide1", htitle, 30, 125., 1500. ) ;

       sprintf( htitle, "MET, LDP, HT[400,500], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht1_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, HT[500,800], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht2_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, HT[800,1000], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht3_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, HT>1000, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht4_nb1", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, LDP, HT[400,500], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht1_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, LDP, HT[500,800], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht2_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, LDP, HT[800,1000], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht3_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, LDP, HT>1000, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht4_nb1_wide1", htitle, 30, 125., 1500. ) ;

       sprintf( htitle, "MET, LDP, no MC trig correction, HT[400,500], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht1_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, HT[500,800], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht2_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, HT[800,1000], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht3_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, HT>1000, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht4_nb0", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, LDP, no MC trig correction, HT[400,500], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht1_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, HT[500,800], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht2_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, HT[800,1000], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht3_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, HT>1000, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_ht4_nb0_wide1", htitle, 30, 125., 1500. ) ;

       sprintf( htitle, "MET, LDP, HT[400,500], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht1_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, HT[500,800], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht2_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, HT[800,1000], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht3_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, HT>1000, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht4_nb0", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, LDP, HT[400,500], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht1_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, LDP, HT[500,800], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht2_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, LDP, HT[800,1000], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht3_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, LDP, HT>1000, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_ldp_ht4_nb0_wide1", htitle, 30, 125., 1500. ) ;




       sprintf( htitle, "MET, ZL, no MC trig correction, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, ZL, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_zl_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_zl_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, ZL, no MC trig correction, HT[400,500], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht1_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, HT[500,800], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht2_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, HT[800,1000], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht3_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, HT>1000, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht4_nb1", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, ZL, HT[400,500], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht1_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, HT[500,800], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht2_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, HT[800,1000], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht3_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, HT>1000, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht4_nb1", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, ZL, HT[400,500], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht1_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, ZL, HT[500,800], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht2_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, ZL, HT[800,1000], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht3_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, ZL, HT>1000, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht4_nb1_wide1", htitle, 30, 125., 1500. ) ;

       sprintf( htitle, "MET, ZL, no MC trig correction, HT[400,500], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht1_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, HT[500,800], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht2_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, HT[800,1000], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht3_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, HT>1000, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht4_nb0", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, ZL, HT[400,500], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht1_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, HT[500,800], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht2_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, HT[800,1000], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht3_nb0", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, ZL, HT>1000, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht4_nb0", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, ZL, HT[400,500], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht1_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, ZL, HT[500,800], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht2_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, ZL, HT[800,1000], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht3_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, ZL, HT>1000, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_ht4_nb0_wide1", htitle, 30, 125., 1500. ) ;


       sprintf( htitle, "MET, ZL, no MC trig correction, HT[400,500], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht1_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, HT[500,800], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht2_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, HT[800,1000], nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht3_nb0_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, HT>1000, nB=0, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht4_nb0_wide1", htitle, 30, 125., 1500. ) ;

       sprintf( htitle, "MET, ZL, no MC trig correction, HT[400,500], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht1_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, HT[500,800], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht2_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, HT[800,1000], nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht3_nb1_wide1", htitle, 30, 125., 1500. ) ;
       sprintf( htitle, "MET, ZL, no MC trig correction, HT>1000, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_zl_nomctrigcorr_ht4_nb1_wide1", htitle, 30, 125., 1500. ) ;



     //-- HT

       sprintf( htitle, "HT, SL, %s", dataset_string ) ;
       bookSet( "h_ht_sl_all", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB=0, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nb3", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, SL, MET<250, %s", dataset_string ) ;
       bookSet( "h_ht_sl_metlt250_all", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, MET<250, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_sl_metlt250_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, MET<250, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_sl_metlt250_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, MET<250, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_sl_metlt250_nb3", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, SL, MET>250, %s", dataset_string ) ;
       bookSet( "h_ht_sl_metgt250_all", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, MET>250, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_sl_metgt250_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, MET>250, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_sl_metgt250_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, MET>250, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_sl_metgt250_nb3", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, SL, no MC trig correction, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_all", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, no MC trig correction, nB=0, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, no MC trig correction, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_nb3", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, SL, nB=1, no MC trig correction, MET[125,150], %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_met1_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB=1, no MC trig correction, MET[150,250], %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_met2_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB=1, no MC trig correction, MET[250,350], %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_met3_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB=1, no MC trig correction, MET>350, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_met4_nb1", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, SL, nB=1, MET[125,150], %s", dataset_string ) ;
       bookSet( "h_ht_sl_met1_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB=1, MET[150,250], %s", dataset_string ) ;
       bookSet( "h_ht_sl_met2_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB=1, MET[250,350], %s", dataset_string ) ;
       bookSet( "h_ht_sl_met3_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB=1, MET>350, %s", dataset_string ) ;
       bookSet( "h_ht_sl_met4_nb1", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, SL, nb=0, no MC trig correction, MET[125,150], %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_met1_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nb=0, no MC trig correction, MET[150,250], %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_met2_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nb=0, no MC trig correction, MET[250,350], %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_met3_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nb=0, no MC trig correction, MET>350, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_met4_nb0", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, SL, nb=0, MET[125,150], %s", dataset_string ) ;
       bookSet( "h_ht_sl_met1_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nb=0, MET[150,250], %s", dataset_string ) ;
       bookSet( "h_ht_sl_met2_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nb=0, MET[250,350], %s", dataset_string ) ;
       bookSet( "h_ht_sl_met3_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nb=0, MET>350, %s", dataset_string ) ;
       bookSet( "h_ht_sl_met4_nb0", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, LDP, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_all", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB=0, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nb3", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, LDP, MET<250, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_metlt250_all", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, MET<250, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_metlt250_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, MET<250, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_metlt250_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, MET<250, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_metlt250_nb3", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, LDP, MET>250, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_metgt250_all", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, MET>250, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_metgt250_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, MET>250, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_metgt250_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, MET>250, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_metgt250_nb3", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, LDP, no MC trig correction, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_all", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, no MC trig correction, nB=0, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, no MC trig correction, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_nb3", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, LDP, nB=1, no MC trig correction, MET[125,150], %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_met1_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB=1, no MC trig correction, MET[150,250], %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_met2_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB=1, no MC trig correction, MET[250,350], %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_met3_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB=1, no MC trig correction, MET>350, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_met4_nb1", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, LDP, nB=1, MET[125,150], %s", dataset_string ) ;
       bookSet( "h_ht_ldp_met1_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB=1, MET[150,250], %s", dataset_string ) ;
       bookSet( "h_ht_ldp_met2_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB=1, MET[250,350], %s", dataset_string ) ;
       bookSet( "h_ht_ldp_met3_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB=1, MET>350, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_met4_nb1", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, LDP, nb=0, no MC trig correction, MET[125,150], %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_met1_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nb=0, no MC trig correction, MET[150,250], %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_met2_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nb=0, no MC trig correction, MET[250,350], %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_met3_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nb=0, no MC trig correction, MET>350, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_met4_nb0", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, LDP, nb=0, MET[125,150], %s", dataset_string ) ;
       bookSet( "h_ht_ldp_met1_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nb=0, MET[150,250], %s", dataset_string ) ;
       bookSet( "h_ht_ldp_met2_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nb=0, MET[250,350], %s", dataset_string ) ;
       bookSet( "h_ht_ldp_met3_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nb=0, MET>350, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_met4_nb0", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, ZL, no MC trig correction, nB=0, %s", dataset_string ) ;
       bookSet( "h_ht_zl_nomctrigcorr_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_zl_nomctrigcorr_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, no MC trig correction, MET<350, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_zl_nomctrigcorr_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, no MC trig correction, MET<350, nB=3, %s", dataset_string ) ;
       bookSet( "h_ht_zl_nomctrigcorr_nb3", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, ZL, nB=0, %s", dataset_string ) ;
       bookSet( "h_ht_zl_nb0", htitle, 30, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_zl_nb1", htitle, 30, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_zl_nb2", htitle, 30, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_zl_nb3", htitle, 30, 400., 2000. ) ;

       sprintf( htitle, "HT, ZL, nB=1, no MC trig correction, MET[125,150], %s", dataset_string ) ;
       bookSet( "h_ht_zl_nomctrigcorr_met1_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nB=1, no MC trig correction, MET[150,250], %s", dataset_string ) ;
       bookSet( "h_ht_zl_nomctrigcorr_met2_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nB=1, no MC trig correction, MET[250,350], %s", dataset_string ) ;
       bookSet( "h_ht_zl_nomctrigcorr_met3_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nB=1, no MC trig correction, MET>350, %s", dataset_string ) ;
       bookSet( "h_ht_zl_nomctrigcorr_met4_nb1", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, ZL, nB=1, MET[125,150], %s", dataset_string ) ;
       bookSet( "h_ht_zl_met1_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nB=1, MET[150,250], %s", dataset_string ) ;
       bookSet( "h_ht_zl_met2_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nB=1, MET[250,350], %s", dataset_string ) ;
       bookSet( "h_ht_zl_met3_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nB=1, MET>350, %s", dataset_string ) ;
       bookSet( "h_ht_zl_met4_nb1", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, ZL, nb=0, no MC trig correction, MET[125,150], %s", dataset_string ) ;
       bookSet( "h_ht_zl_nomctrigcorr_met1_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nb=0, no MC trig correction, MET[150,250], %s", dataset_string ) ;
       bookSet( "h_ht_zl_nomctrigcorr_met2_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nb=0, no MC trig correction, MET[250,350], %s", dataset_string ) ;
       bookSet( "h_ht_zl_nomctrigcorr_met3_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nb=0, no MC trig correction, MET>350, %s", dataset_string ) ;
       bookSet( "h_ht_zl_nomctrigcorr_met4_nb0", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "HT, ZL, nb=0, MET[125,150], %s", dataset_string ) ;
       bookSet( "h_ht_zl_met1_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nb=0, MET[150,250], %s", dataset_string ) ;
       bookSet( "h_ht_zl_met2_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nb=0, MET[250,350], %s", dataset_string ) ;
       bookSet( "h_ht_zl_met3_nb0", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, ZL, nb=0, MET>350, %s", dataset_string ) ;
       bookSet( "h_ht_zl_met4_nb0", htitle, 32, 400., 2000. ) ;








   //-- likelihood bins

       sprintf( htitle, "HT vs MET, SL, nB=0, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nb0", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, nB=1, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nb1", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, nB=2, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nb2", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, nB=3, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nb3", htitle ) ;

       sprintf( htitle, "HT vs MET, LDP, nB=0, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nb0", htitle ) ;
       sprintf( htitle, "HT vs MET, LDP, nB=1, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nb1", htitle ) ;
       sprintf( htitle, "HT vs MET, LDP, nB=2, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nb2", htitle ) ;
       sprintf( htitle, "HT vs MET, LDP, nB=3, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nb3", htitle ) ;


       sprintf( htitle, "HT vs MET, SL, no MC trig correction, nB=0, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nomctrigcorr_nb0", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nomctrigcorr_nb1", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nomctrigcorr_nb2", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, no MC trig correction, nB=3, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nomctrigcorr_nb3", htitle ) ;

       sprintf( htitle, "HT vs MET, LDP, no MC trig correction, nB=0, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nomctrigcorr_nb0", htitle ) ;
       sprintf( htitle, "HT vs MET, LDP, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nomctrigcorr_nb1", htitle ) ;
       sprintf( htitle, "HT vs MET, LDP, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nomctrigcorr_nb2", htitle ) ;
       sprintf( htitle, "HT vs MET, LDP, no MC trig correction, nB=3, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nomctrigcorr_nb3", htitle ) ;


       sprintf( htitle, "HT vs MET, ZL, nB=0, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_zl_nb0", htitle ) ;
       sprintf( htitle, "HT vs MET, ZL, nB=1, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_zl_nb1", htitle ) ;
       sprintf( htitle, "HT vs MET, ZL, nB=2, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_zl_nb2", htitle ) ;
       sprintf( htitle, "HT vs MET, ZL, nB>=3, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_zl_nb3", htitle ) ;

       sprintf( htitle, "HT vs MET, ZL, no MC trig correction, nB=0, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_zl_nomctrigcorr_nb0", htitle ) ;
       sprintf( htitle, "HT vs MET, ZL, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_zl_nomctrigcorr_nb1", htitle ) ;
       sprintf( htitle, "HT vs MET, ZL, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_zl_nomctrigcorr_nb2", htitle ) ;
       sprintf( htitle, "HT vs MET, ZL, no MC trig correction, nB>=3, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_zl_nomctrigcorr_nb3", htitle ) ;



     //--- n b jets

       sprintf( htitle, "nbjets, SL, %s", dataset_string ) ;
       bookSet( "h_nb_sl_all", htitle, 6, -0.5, 5.5 ) ;

       sprintf( htitle, "nbjets, LDP, %s", dataset_string ) ;
       bookSet( "h_nb_ldp_all", htitle, 6, -0.5, 5.5 ) ;

       sprintf( htitle, "nbjets, ZL, %s", dataset_string ) ;
       bookSet( "h_nb_zl_all", htitle, 6, -0.5, 5.5 ) ;


       sprintf( htitle, "nbjets, SL, %s", dataset_string ) ;
       bookSet( "h_nbgt0_sl_all", htitle, 6, -0.5, 5.5 ) ;

       sprintf( htitle, "nbjets, LDP, %s", dataset_string ) ;
       bookSet( "h_nbgt0_ldp_all", htitle, 6, -0.5, 5.5 ) ;

       sprintf( htitle, "nbjets, ZL, %s", dataset_string ) ;
       bookSet( "h_nbgt0_zl_all", htitle, 6, -0.5, 5.5 ) ;




     //--- N jets

       sprintf( htitle, "N jets, SL, %s", dataset_string ) ;
       bookSet( "h_njets_sl_all", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, SL, nB=0, %s", dataset_string ) ;
       bookSet( "h_njets_sl_nb0", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, SL, nB=1, %s", dataset_string ) ;
       bookSet( "h_njets_sl_nb1", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, SL, nB=2, %s", dataset_string ) ;
       bookSet( "h_njets_sl_nb2", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, SL, nB=3, %s", dataset_string ) ;
       bookSet( "h_njets_sl_nb3", htitle, 11, -0.5, 10.5 ) ;

       sprintf( htitle, "N jets, LDP, %s", dataset_string ) ;
       bookSet( "h_njets_ldp_all", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, LDP, nB=0, %s", dataset_string ) ;
       bookSet( "h_njets_ldp_nb0", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, LDP, nB=1, %s", dataset_string ) ;
       bookSet( "h_njets_ldp_nb1", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, LDP, nB=2, %s", dataset_string ) ;
       bookSet( "h_njets_ldp_nb2", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, LDP, nB=3, %s", dataset_string ) ;
       bookSet( "h_njets_ldp_nb3", htitle, 11, -0.5, 10.5 ) ;

       sprintf( htitle, "N jets, ZL, %s", dataset_string ) ;
       bookSet( "h_njets_zl_all", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, ZL, nB=0, %s", dataset_string ) ;
       bookSet( "h_njets_zl_nb0", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, ZL, nB=1, %s", dataset_string ) ;
       bookSet( "h_njets_zl_nb1", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, ZL, nB=2, %s", dataset_string ) ;
       bookSet( "h_njets_zl_nb2", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, ZL, nB=3, %s", dataset_string ) ;
       bookSet( "h_njets_zl_nb3", htitle, 11, -0.5, 10.5 ) ;


      //--- MT

       sprintf( htitle, "MT, SL, %s", dataset_string ) ;
       bookSet( "h_mt_sl_all", htitle, 25, 0., 250. ) ;
       sprintf( htitle, "MT, SL, nB=0, %s", dataset_string ) ;
       bookSet( "h_mt_sl_nb0", htitle, 25, 0., 250. ) ;
       sprintf( htitle, "MT, SL, nB=1, %s", dataset_string ) ;
       bookSet( "h_mt_sl_nb1", htitle, 25, 0., 250. ) ;
       sprintf( htitle, "MT, SL, nB=2, %s", dataset_string ) ;
       bookSet( "h_mt_sl_nb2", htitle, 25, 0., 250. ) ;
       sprintf( htitle, "MT, SL, nB=3, %s", dataset_string ) ;
       bookSet( "h_mt_sl_nb3", htitle, 25, 0., 250. ) ;


      //--- npv

       sprintf( htitle, "N PV, SL, nB=0, %s", dataset_string ) ;
       bookSet( "h_npv_sl_nb0", htitle, 36, -0.5, 35.5 ) ;
       sprintf( htitle, "N PV, LDP, nB=0, %s", dataset_string ) ;
       bookSet( "h_npv_ldp_nb0", htitle, 36, -0.5, 35.5 ) ;
       sprintf( htitle, "N PV, ZL, nB=0, %s", dataset_string ) ;
       bookSet( "h_npv_zl_nb0", htitle, 36, -0.5, 35.5 ) ;

       sprintf( htitle, "N PV, SL, nB=1, %s", dataset_string ) ;
       bookSet( "h_npv_sl_nb1", htitle, 36, -0.5, 35.5 ) ;
       sprintf( htitle, "N PV, LDP, nB=1, %s", dataset_string ) ;
       bookSet( "h_npv_ldp_nb1", htitle, 36, -0.5, 35.5 ) ;
       sprintf( htitle, "N PV, ZL, nB=1, %s", dataset_string ) ;
       bookSet( "h_npv_zl_nb1", htitle, 36, -0.5, 35.5 ) ;




      //--- min delta phi

       sprintf( htitle, "min Delta Phi, 0 leptons, nB=0 %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_nb0", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=0, MET1, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_met1_nb0", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=0, MET2, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_met2_nb0", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=0, MET3, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_met3_nb0", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=0, MET4, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_met4_nb0", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=0, HT1, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_ht1_nb0", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=0, HT2, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_ht2_nb0", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=0, HT3, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_ht3_nb0", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=0, HT4, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_ht4_nb0", htitle, 64, 0., 3.2 ) ;

       sprintf( htitle, "min Delta Phi, 0 leptons, nB=1 %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_nb1", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=1, MET1, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_met1_nb1", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=1, MET2, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_met2_nb1", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=1, MET3, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_met3_nb1", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=1, MET4, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_met4_nb1", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=1, HT1, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_ht1_nb1", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=1, HT2, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_ht2_nb1", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=1, HT3, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_ht3_nb1", htitle, 64, 0., 3.2 ) ;
       sprintf( htitle, "min Delta Phi, 0 leptons, nB=1, HT4, %s", dataset_string ) ;
       bookSet( "h_mindphi_0lep_ht4_nb1", htitle, 64, 0., 3.2 ) ;





     //--- min Delta Phi N

       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=0 %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_nb0", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=0, MET1, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_met1_nb0", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=0, MET2, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_met2_nb0", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=0, MET3, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_met3_nb0", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=0, MET4, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_met4_nb0", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=0, HT1, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_ht1_nb0", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=0, HT2, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_ht2_nb0", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=0, HT3, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_ht3_nb0", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=0, HT4, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_ht4_nb0", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=1 %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=1, MET1, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_met1_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=1, MET2, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_met2_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=1, MET3, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_met3_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=1, MET4, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_met4_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=1, HT1, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_ht1_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=1, HT2, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_ht2_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=1, HT3, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_ht3_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 0 leptons, nB=1, HT4, %s", dataset_string ) ;
       bookSet( "h_mindphin_0lep_ht4_nb1", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "min Delta Phi N, 1 lepton, nB=1 %s", dataset_string ) ;
       bookSet( "h_mindphin_1lep_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 1 lepton, nB=1, MET1, %s", dataset_string ) ;
       bookSet( "h_mindphin_1lep_met1_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 1 lepton, nB=1, MET2, %s", dataset_string ) ;
       bookSet( "h_mindphin_1lep_met2_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 1 lepton, nB=1, MET3, %s", dataset_string ) ;
       bookSet( "h_mindphin_1lep_met3_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 1 lepton, nB=1, MET4, %s", dataset_string ) ;
       bookSet( "h_mindphin_1lep_met4_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 1 lepton, nB=1, HT1, %s", dataset_string ) ;
       bookSet( "h_mindphin_1lep_ht1_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 1 lepton, nB=1, HT2, %s", dataset_string ) ;
       bookSet( "h_mindphin_1lep_ht2_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 1 lepton, nB=1, HT3, %s", dataset_string ) ;
       bookSet( "h_mindphin_1lep_ht3_nb1", htitle, 30, 0., 30. ) ;
       sprintf( htitle, "min Delta Phi N, 1 lepton, nB=1, HT4, %s", dataset_string ) ;
       bookSet( "h_mindphin_1lep_ht4_nb1", htitle, 30, 0., 30. ) ;


      //--- met ratio

       sprintf( htitle, "pfMET/caloMET, LDP, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_ldp", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_ldp_wide", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "pfMET/caloMET, LDP, nB=0, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_ldp_nb0", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_ldp_nb0_wide", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "pfMET/caloMET, LDP, nB=1, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_ldp_nb1", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_ldp_nb1_wide", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "pfMET/caloMET, LDP, nB=2, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_ldp_nb2", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_ldp_nb2_wide", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "pfMET/caloMET, LDP, nB>=3, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_ldp_nb3", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_ldp_nb3_wide", htitle, 30, 0., 30. ) ;


       sprintf( htitle, "pfMET/caloMET, SL, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_sl", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_sl_wide", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "pfMET/caloMET, SL, nB=0, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_sl_nb0", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_sl_nb0_wide", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "pfMET/caloMET, SL, nB=1, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_sl_nb1", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_sl_nb1_wide", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "pfMET/caloMET, SL, nB=2, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_sl_nb2", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_sl_nb2_wide", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "pfMET/caloMET, SL, nB>=3, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_sl_nb3", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_sl_nb3_wide", htitle, 30, 0., 30. ) ;




       sprintf( htitle, "pfMET/caloMET, ZL, nB=0, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_zl_nb0", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_zl_nb0_wide", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "pfMET/caloMET, ZL, nB=1, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_zl_nb1", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_zl_nb1_wide", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "pfMET/caloMET, ZL, nB=2, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_zl_nb2", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_zl_nb2_wide", htitle, 30, 0., 30. ) ;

       sprintf( htitle, "pfMET/caloMET, ZL, nB>=3, %s", dataset_string ) ;
       bookSet( "h_pfOcaloMET_zl_nb3", htitle, 30, 0., 3. ) ;
       bookSet( "h_pfOcaloMET_zl_nb3_wide", htitle, 30, 0., 30. ) ;





       sprintf( htitle, "N iso track, ZL, nB=0, %s", dataset_string ) ;
       bookSet( "h_nisotrk_zl_nb0", htitle, 6, -0.5, 5.5 ) ;
       sprintf( htitle, "N iso track, ZL, nB=1, %s", dataset_string ) ;
       bookSet( "h_nisotrk_zl_nb1", htitle, 6, -0.5, 5.5 ) ;
       sprintf( htitle, "N iso track, ZL, nB=2, %s", dataset_string ) ;
       bookSet( "h_nisotrk_zl_nb2", htitle, 6, -0.5, 5.5 ) ;
       sprintf( htitle, "N iso track, ZL, nB>=3, %s", dataset_string ) ;
       bookSet( "h_nisotrk_zl_nb3", htitle, 6, -0.5, 5.5 ) ;

       sprintf( htitle, "N iso track, LDP, nB=0, %s", dataset_string ) ;
       bookSet( "h_nisotrk_ldp_nb0", htitle, 6, -0.5, 5.5 ) ;
       sprintf( htitle, "N iso track, LDP, nB=1, %s", dataset_string ) ;
       bookSet( "h_nisotrk_ldp_nb1", htitle, 6, -0.5, 5.5 ) ;
       sprintf( htitle, "N iso track, LDP, nB=2, %s", dataset_string ) ;
       bookSet( "h_nisotrk_ldp_nb2", htitle, 6, -0.5, 5.5 ) ;
       sprintf( htitle, "N iso track, LDP, nB>=3, %s", dataset_string ) ;
       bookSet( "h_nisotrk_ldp_nb3", htitle, 6, -0.5, 5.5 ) ;

       sprintf( htitle, "N iso track, SL, nB=0, %s", dataset_string ) ;
       bookSet( "h_nisotrk_sl_nb0", htitle, 6, -0.5, 5.5 ) ;
       sprintf( htitle, "N iso track, SL, nB=1, %s", dataset_string ) ;
       bookSet( "h_nisotrk_sl_nb1", htitle, 6, -0.5, 5.5 ) ;
       sprintf( htitle, "N iso track, SL, nB=2, %s", dataset_string ) ;
       bookSet( "h_nisotrk_sl_nb2", htitle, 6, -0.5, 5.5 ) ;
       sprintf( htitle, "N iso track, SL, nB>=3, %s", dataset_string ) ;
       bookSet( "h_nisotrk_sl_nb3", htitle, 6, -0.5, 5.5 ) ;






       sprintf( htitle, "Max ch mult diff, ZL, nB=0, %s", dataset_string ) ;
       bookSet( "h_maxchnmultdiff_zl_nb0", htitle, 53, -106, 106 ) ;
       sprintf( htitle, "Max ch mult diff, ZL, nB=1, %s", dataset_string ) ;
       bookSet( "h_maxchnmultdiff_zl_nb1", htitle, 53, -106, 106 ) ;
       sprintf( htitle, "Max ch mult diff, ZL, nB=2, %s", dataset_string ) ;
       bookSet( "h_maxchnmultdiff_zl_nb2", htitle, 53, -106, 106 ) ;
       sprintf( htitle, "Max ch mult diff, ZL, nB>=3, %s", dataset_string ) ;
       bookSet( "h_maxchnmultdiff_zl_nb3", htitle, 53, -106, 106 ) ;

       sprintf( htitle, "Max ch mult diff, LDP, nB=0, %s", dataset_string ) ;
       bookSet( "h_maxchnmultdiff_ldp_nb0", htitle, 53, -106, 106 ) ;
       sprintf( htitle, "Max ch mult diff, LDP, nB=1, %s", dataset_string ) ;
       bookSet( "h_maxchnmultdiff_ldp_nb1", htitle, 53, -106, 106 ) ;
       sprintf( htitle, "Max ch mult diff, LDP, nB=2, %s", dataset_string ) ;
       bookSet( "h_maxchnmultdiff_ldp_nb2", htitle, 53, -106, 106 ) ;
       sprintf( htitle, "Max ch mult diff, LDP, nB>=3, %s", dataset_string ) ;
       bookSet( "h_maxchnmultdiff_ldp_nb3", htitle, 53, -106, 106 ) ;

       sprintf( htitle, "Max ch mult diff, SL, nB=0, %s", dataset_string ) ;
       bookSet( "h_maxchnmultdiff_sl_nb0", htitle, 53, -106, 106 ) ;
       sprintf( htitle, "Max ch mult diff, SL, nB=1, %s", dataset_string ) ;
       bookSet( "h_maxchnmultdiff_sl_nb1", htitle, 53, -106, 106 ) ;
       sprintf( htitle, "Max ch mult diff, SL, nB=2, %s", dataset_string ) ;
       bookSet( "h_maxchnmultdiff_sl_nb2", htitle, 53, -106, 106 ) ;
       sprintf( htitle, "Max ch mult diff, SL, nB>=3, %s", dataset_string ) ;
       bookSet( "h_maxchnmultdiff_sl_nb3", htitle, 53, -106, 106 ) ;











     //--- Fill the histograms. -----------------------------------------------

       char cuts[10000] ;

       sprintf( cuts, "%s", basecuts_1lep_nonb ) ;
       fillSet( "h_nb_sl_all", "nB", cuts ) ;

       sprintf( cuts, "%s", basecuts_ldp_nonb ) ;
       fillSet( "h_nb_ldp_all", "nB", cuts ) ;

       sprintf( cuts, "%s", basecuts_0lep_nonb ) ;
       fillSet( "h_nb_zl_all", "nB", cuts ) ;


       sprintf( cuts, "%s&&nB>0", basecuts_1lep_nonb ) ;
       fillSet( "h_nbgt0_sl_all", "nB", cuts ) ;

       sprintf( cuts, "%s&&nB>0", basecuts_ldp_nonb ) ;
       fillSet( "h_nbgt0_ldp_all", "nB", cuts ) ;

       sprintf( cuts, "%s&&nB>0", basecuts_0lep_nonb ) ;
       fillSet( "h_nbgt0_zl_all", "nB", cuts ) ;


       sprintf( cuts, "%s", basecuts_1lep ) ;
       fillSet( "h_met_sl_all", "MET", cuts ) ;
       fillSet( "h_ht_sl_all", "HT", cuts ) ;
       fillSet( "h_met_sl_nomctrigcorr_all", "MET", cuts, false ) ;
       fillSet( "h_ht_sl_nomctrigcorr_all", "HT", cuts, false ) ;

       sprintf( cuts, "%s&&nB==0", basecuts_1lep_nonb ) ;
       fillSet( "h_met_sl_nb0", "MET", cuts ) ;
       fillSet( "h_ht_sl_nb0", "HT", cuts ) ;
       fillSet( "h_met_sl_nomctrigcorr_nb0", "MET", cuts, false ) ;
       fillSet( "h_ht_sl_nomctrigcorr_nb0", "HT", cuts, false ) ;
       fillSetLHB( "h_lhb_sl_nb0", cuts ) ;
       fillSetLHB( "h_lhb_sl_nomctrigcorr_nb0", cuts, false ) ;

       sprintf( cuts, "%s&&nB==1", basecuts_1lep ) ;
       fillSet( "h_met_sl_nb1", "MET", cuts ) ;
       fillSet( "h_ht_sl_nb1", "HT", cuts ) ;
       fillSet( "h_met_sl_nomctrigcorr_nb1", "MET", cuts, false ) ;
       fillSet( "h_ht_sl_nomctrigcorr_nb1", "HT", cuts, false ) ;
       fillSetLHB( "h_lhb_sl_nb1", cuts ) ;
       fillSetLHB( "h_lhb_sl_nomctrigcorr_nb1", cuts, false ) ;

       sprintf( cuts, "%s&&nB==2", basecuts_1lep ) ;
       fillSet( "h_met_sl_nb2", "MET", cuts ) ;
       fillSet( "h_ht_sl_nb2", "HT", cuts ) ;
       fillSet( "h_met_sl_nomctrigcorr_nb2", "MET", cuts, false ) ;
       fillSet( "h_ht_sl_nomctrigcorr_nb2", "HT", cuts, false ) ;
       fillSetLHB( "h_lhb_sl_nb2", cuts ) ;
       fillSetLHB( "h_lhb_sl_nomctrigcorr_nb2", cuts, false ) ;

       sprintf( cuts, "%s&&nB>=3", basecuts_1lep ) ;
       fillSet( "h_met_sl_nb3", "MET", cuts ) ;
       fillSet( "h_ht_sl_nb3", "HT", cuts ) ;
       fillSet( "h_met_sl_nomctrigcorr_nb3", "MET", cuts, false ) ;
       fillSet( "h_ht_sl_nomctrigcorr_nb3", "HT", cuts, false ) ;
       fillSetLHB( "h_lhb_sl_nb3", cuts ) ;
       fillSetLHB( "h_lhb_sl_nomctrigcorr_nb3", cuts, false ) ;


       sprintf( cuts, "%s&&HT<800", basecuts_1lep ) ;
       fillSet( "h_met_sl_htlt800_all", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT<800&&nB==1", basecuts_1lep ) ;
       fillSet( "h_met_sl_htlt800_nb1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT<800&&nB==2", basecuts_1lep ) ;
       fillSet( "h_met_sl_htlt800_nb2", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT<800&&nB>=3", basecuts_1lep ) ;
       fillSet( "h_met_sl_htlt800_nb3", "MET", cuts ) ;

       sprintf( cuts, "%s&&HT>800", basecuts_1lep ) ;
       fillSet( "h_met_sl_htgt800_all", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>800&&nB==1", basecuts_1lep ) ;
       fillSet( "h_met_sl_htgt800_nb1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>800&&nB==2", basecuts_1lep ) ;
       fillSet( "h_met_sl_htgt800_nb2", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>800&&nB>=3", basecuts_1lep ) ;
       fillSet( "h_met_sl_htgt800_nb3", "MET", cuts ) ;


       sprintf( cuts, "%s&&MET<250", basecuts_1lep ) ;
       fillSet( "h_ht_sl_metlt250_all", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET<250&&nB==1", basecuts_1lep ) ;
       fillSet( "h_ht_sl_metlt250_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET<250&&nB==2", basecuts_1lep ) ;
       fillSet( "h_ht_sl_metlt250_nb2", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET<250&&nB>=3", basecuts_1lep ) ;
       fillSet( "h_ht_sl_metlt250_nb3", "HT", cuts ) ;

       sprintf( cuts, "%s&&MET>250", basecuts_1lep ) ;
       fillSet( "h_ht_sl_metgt250_all", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>250&&nB==1", basecuts_1lep ) ;
       fillSet( "h_ht_sl_metgt250_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>250&&nB==2", basecuts_1lep ) ;
       fillSet( "h_ht_sl_metgt250_nb2", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>250&&nB>=3", basecuts_1lep ) ;
       fillSet( "h_ht_sl_metgt250_nb3", "HT", cuts ) ;




       sprintf( cuts, "%s", basecuts_ldp ) ;
       fillSet( "h_met_ldp_all", "MET", cuts ) ;
       fillSet( "h_ht_ldp_all", "HT", cuts ) ;
       fillSet( "h_met_ldp_nomctrigcorr_all", "MET", cuts, false ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_all", "HT", cuts, false ) ;

       sprintf( cuts, "%s&&nB==0", basecuts_ldp_nonb ) ;
       fillSet( "h_met_ldp_nb0", "MET", cuts ) ;
       fillSet( "h_ht_ldp_nb0", "HT", cuts ) ;
       fillSet( "h_met_ldp_nomctrigcorr_nb0", "MET", cuts, false ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_nb0", "HT", cuts, false ) ;
       fillSetLHB( "h_lhb_ldp_nb0", cuts ) ;
       fillSetLHB( "h_lhb_ldp_nomctrigcorr_nb0", cuts, false ) ;

       sprintf( cuts, "%s&&nB==1", basecuts_ldp ) ;
       fillSet( "h_met_ldp_nb1", "MET", cuts ) ;
       fillSet( "h_ht_ldp_nb1", "HT", cuts ) ;
       fillSet( "h_met_ldp_nomctrigcorr_nb1", "MET", cuts, false ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_nb1", "HT", cuts, false ) ;
       fillSetLHB( "h_lhb_ldp_nb1", cuts ) ;
       fillSetLHB( "h_lhb_ldp_nomctrigcorr_nb1", cuts, false ) ;

       sprintf( cuts, "%s&&nB==2", basecuts_ldp ) ;
       fillSet( "h_met_ldp_nb2", "MET", cuts ) ;
       fillSet( "h_ht_ldp_nb2", "HT", cuts ) ;
       fillSet( "h_met_ldp_nomctrigcorr_nb2", "MET", cuts, false ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_nb2", "HT", cuts, false ) ;
       fillSetLHB( "h_lhb_ldp_nb2", cuts ) ;
       fillSetLHB( "h_lhb_ldp_nomctrigcorr_nb2", cuts, false ) ;

       sprintf( cuts, "%s&&nB>=3", basecuts_ldp ) ;
       fillSet( "h_met_ldp_nb3", "MET", cuts ) ;
       fillSet( "h_ht_ldp_nb3", "HT", cuts ) ;
       fillSet( "h_met_ldp_nomctrigcorr_nb3", "MET", cuts, false ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_nb3", "HT", cuts, false ) ;
       fillSetLHB( "h_lhb_ldp_nb3", cuts ) ;
       fillSetLHB( "h_lhb_ldp_nomctrigcorr_nb3", cuts, false ) ;


       sprintf( cuts, "%s&&HT<800", basecuts_ldp ) ;
       fillSet( "h_met_ldp_htlt800_all", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT<800&&nB==1", basecuts_ldp ) ;
       fillSet( "h_met_ldp_htlt800_nb1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT<800&&nB==2", basecuts_ldp ) ;
       fillSet( "h_met_ldp_htlt800_nb2", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT<800&&nB>=3", basecuts_ldp ) ;
       fillSet( "h_met_ldp_htlt800_nb3", "MET", cuts ) ;

       sprintf( cuts, "%s&&HT>800", basecuts_ldp ) ;
       fillSet( "h_met_ldp_htgt800_all", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>800&&nB==1", basecuts_ldp ) ;
       fillSet( "h_met_ldp_htgt800_nb1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>800&&nB==2", basecuts_ldp ) ;
       fillSet( "h_met_ldp_htgt800_nb2", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>800&&nB>=3", basecuts_ldp ) ;
       fillSet( "h_met_ldp_htgt800_nb3", "MET", cuts ) ;


       sprintf( cuts, "%s&&MET<250", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_metlt250_all", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET<250&&nB==1", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_metlt250_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET<250&&nB==2", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_metlt250_nb2", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET<250&&nB>=3", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_metlt250_nb3", "HT", cuts ) ;

       sprintf( cuts, "%s&&MET>250", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_metgt250_all", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>250&&nB==1", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_metgt250_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>250&&nB==2", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_metgt250_nb2", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>250&&nB>=3", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_metgt250_nb3", "HT", cuts ) ;


       sprintf( cuts, "%s&&HT>400&&HT<=500&&nB==1", basecuts_1lep ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht1_nb1", "MET", cuts, false ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht1_nb1_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_sl_ht1_nb1", "MET", cuts ) ;
       fillSet( "h_met_sl_ht1_nb1_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>500&&HT<=800&&nB==1", basecuts_1lep ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht2_nb1", "MET", cuts, false ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht2_nb1_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_sl_ht2_nb1", "MET", cuts ) ;
       fillSet( "h_met_sl_ht2_nb1_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>800&&HT<=1000&&nB==1", basecuts_1lep ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht3_nb1", "MET", cuts, false ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht3_nb1_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_sl_ht3_nb1", "MET", cuts ) ;
       fillSet( "h_met_sl_ht3_nb1_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>1000&&nB==1", basecuts_1lep ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht4_nb1", "MET", cuts, false ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht4_nb1_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_sl_ht4_nb1", "MET", cuts ) ;
       fillSet( "h_met_sl_ht4_nb1_wide1", "MET", cuts ) ;

       sprintf( cuts, "%s&&HT>400&&HT<=500&&nB==0", basecuts_1lep_nonb ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht1_nb0", "MET", cuts, false ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht1_nb0_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_sl_ht1_nb0", "MET", cuts ) ;
       fillSet( "h_met_sl_ht1_nb0_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>500&&HT<=800&&nB==0", basecuts_1lep_nonb ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht2_nb0", "MET", cuts, false ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht2_nb0_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_sl_ht2_nb0", "MET", cuts ) ;
       fillSet( "h_met_sl_ht2_nb0_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>800&&HT<=1000&&nB==0", basecuts_1lep_nonb ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht3_nb0", "MET", cuts, false ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht3_nb0_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_sl_ht3_nb0", "MET", cuts ) ;
       fillSet( "h_met_sl_ht3_nb0_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>1000&&nB==0", basecuts_1lep_nonb ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht4_nb0", "MET", cuts, false ) ;
       fillSet( "h_met_sl_nomctrigcorr_ht4_nb0_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_sl_ht4_nb0", "MET", cuts ) ;
       fillSet( "h_met_sl_ht4_nb0_wide1", "MET", cuts ) ;

       sprintf( cuts, "%s&&HT>400&&HT<=500&&nB==1", basecuts_ldp ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht1_nb1", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht1_nb1_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_ht1_nb1", "MET", cuts ) ;
       fillSet( "h_met_ldp_ht1_nb1_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>500&&HT<=800&&nB==1", basecuts_ldp ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht2_nb1", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht2_nb1_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_ht2_nb1", "MET", cuts ) ;
       fillSet( "h_met_ldp_ht2_nb1_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>800&&HT<=1000&&nB==1", basecuts_ldp ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht3_nb1", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht3_nb1_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_ht3_nb1", "MET", cuts ) ;
       fillSet( "h_met_ldp_ht3_nb1_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>1000&&nB==1", basecuts_ldp ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht4_nb1", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht4_nb1_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_ht4_nb1", "MET", cuts ) ;
       fillSet( "h_met_ldp_ht4_nb1_wide1", "MET", cuts ) ;

       sprintf( cuts, "%s&&HT>400&&HT<=500&&nB==0", basecuts_ldp_nonb ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht1_nb0", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht1_nb0_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_ht1_nb0", "MET", cuts ) ;
       fillSet( "h_met_ldp_ht1_nb0_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>500&&HT<=800&&nB==0", basecuts_ldp_nonb ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht2_nb0", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht2_nb0_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_ht2_nb0", "MET", cuts ) ;
       fillSet( "h_met_ldp_ht2_nb0_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>800&&HT<=1000&&nB==0", basecuts_ldp_nonb ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht3_nb0", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht3_nb0_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_ht3_nb0", "MET", cuts ) ;
       fillSet( "h_met_ldp_ht3_nb0_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>1000&&nB==0", basecuts_ldp_nonb ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht4_nb0", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_nomctrigcorr_ht4_nb0_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_ldp_ht4_nb0", "MET", cuts ) ;
       fillSet( "h_met_ldp_ht4_nb0_wide1", "MET", cuts ) ;

       sprintf( cuts, "%s&&HT>400&&HT<=500&&nB==1", basecuts_0lep ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht1_nb1", "MET", cuts, false ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht1_nb1_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_zl_ht1_nb1", "MET", cuts ) ;
       fillSet( "h_met_zl_ht1_nb1_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>500&&HT<=800&&nB==1", basecuts_0lep ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht2_nb1", "MET", cuts, false ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht2_nb1_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_zl_ht2_nb1", "MET", cuts ) ;
       fillSet( "h_met_zl_ht2_nb1_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>800&&HT<=1000&&nB==1", basecuts_0lep ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht3_nb1", "MET", cuts, false ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht3_nb1_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_zl_ht3_nb1", "MET", cuts ) ;
       fillSet( "h_met_zl_ht3_nb1_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>1000&&nB==1", basecuts_0lep ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht4_nb1", "MET", cuts, false ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht4_nb1_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_zl_ht4_nb1", "MET", cuts ) ;
       fillSet( "h_met_zl_ht4_nb1_wide1", "MET", cuts ) ;

       sprintf( cuts, "%s&&HT>400&&HT<=500&&nB==0", basecuts_0lep_nonb ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht1_nb0", "MET", cuts, false ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht1_nb0_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_zl_ht1_nb0", "MET", cuts ) ;
       fillSet( "h_met_zl_ht1_nb0_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>500&&HT<=800&&nB==0", basecuts_0lep_nonb ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht2_nb0", "MET", cuts, false ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht2_nb0_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_zl_ht2_nb0", "MET", cuts ) ;
       fillSet( "h_met_zl_ht2_nb0_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>800&&HT<=1000&&nB==0", basecuts_0lep_nonb ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht3_nb0", "MET", cuts, false ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht3_nb0_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_zl_ht3_nb0", "MET", cuts ) ;
       fillSet( "h_met_zl_ht3_nb0_wide1", "MET", cuts ) ;
       sprintf( cuts, "%s&&HT>1000&&nB==0", basecuts_0lep_nonb ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht4_nb0", "MET", cuts, false ) ;
       fillSet( "h_met_zl_nomctrigcorr_ht4_nb0_wide1", "MET", cuts, false ) ;
       fillSet( "h_met_zl_ht4_nb0", "MET", cuts ) ;
       fillSet( "h_met_zl_ht4_nb0_wide1", "MET", cuts ) ;





       sprintf( cuts, "%s&&MET>125&&MET<=150&&nB==1", basecuts_1lep ) ;
       fillSet( "h_ht_sl_nomctrigcorr_met1_nb1", "HT", cuts, false ) ;
       fillSet( "h_ht_sl_met1_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>150&&MET<=250&&nB==1", basecuts_1lep ) ;
       fillSet( "h_ht_sl_nomctrigcorr_met2_nb1", "HT", cuts, false ) ;
       fillSet( "h_ht_sl_met2_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>250&&MET<=350&&nB==1", basecuts_1lep ) ;
       fillSet( "h_ht_sl_nomctrigcorr_met3_nb1", "HT", cuts, false ) ;
       fillSet( "h_ht_sl_met3_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>350&&nB==1", basecuts_1lep ) ;
       fillSet( "h_ht_sl_nomctrigcorr_met4_nb1", "HT", cuts, false ) ;
       fillSet( "h_ht_sl_met4_nb1", "HT", cuts ) ;

       sprintf( cuts, "%s&&MET>125&&MET<=150&&nB==0", basecuts_1lep_nonb ) ;
       fillSet( "h_ht_sl_nomctrigcorr_met1_nb0", "HT", cuts, false ) ;
       fillSet( "h_ht_sl_met1_nb0", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>150&&MET<=250&&nB==0", basecuts_1lep_nonb ) ;
       fillSet( "h_ht_sl_nomctrigcorr_met2_nb0", "HT", cuts, false ) ;
       fillSet( "h_ht_sl_met2_nb0", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>250&&MET<=350&&nB==0", basecuts_1lep_nonb ) ;
       fillSet( "h_ht_sl_nomctrigcorr_met3_nb0", "HT", cuts, false ) ;
       fillSet( "h_ht_sl_met3_nb0", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>350&&nB==0", basecuts_1lep_nonb ) ;
       fillSet( "h_ht_sl_nomctrigcorr_met4_nb0", "HT", cuts, false ) ;
       fillSet( "h_ht_sl_met4_nb0", "HT", cuts ) ;


       sprintf( cuts, "%s&&MET>125&&MET<=150&&nB==1", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_met1_nb1", "HT", cuts, false ) ;
       fillSet( "h_ht_ldp_met1_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>150&&MET<=250&&nB==1", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_met2_nb1", "HT", cuts, false ) ;
       fillSet( "h_ht_ldp_met2_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>250&&MET<=350&&nB==1", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_met3_nb1", "HT", cuts, false ) ;
       fillSet( "h_ht_ldp_met3_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>350&&nB==1", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_met4_nb1", "HT", cuts, false ) ;
       fillSet( "h_ht_ldp_met4_nb1", "HT", cuts ) ;

       sprintf( cuts, "%s&&MET>125&&MET<=150&&nB==0", basecuts_ldp_nonb ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_met1_nb0", "HT", cuts, false ) ;
       fillSet( "h_ht_ldp_met1_nb0", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>150&&MET<=250&&nB==0", basecuts_ldp_nonb ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_met2_nb0", "HT", cuts, false ) ;
       fillSet( "h_ht_ldp_met2_nb0", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>250&&MET<=350&&nB==0", basecuts_ldp_nonb ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_met3_nb0", "HT", cuts, false ) ;
       fillSet( "h_ht_ldp_met3_nb0", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>350&&nB==0", basecuts_ldp_nonb ) ;
       fillSet( "h_ht_ldp_nomctrigcorr_met4_nb0", "HT", cuts, false ) ;
       fillSet( "h_ht_ldp_met4_nb0", "HT", cuts ) ;


       sprintf( cuts, "%s&&MET>125&&MET<=150&&nB==1", basecuts_0lep ) ;
       fillSet( "h_ht_zl_nomctrigcorr_met1_nb1", "HT", cuts, false ) ;
       fillSet( "h_ht_zl_met1_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>150&&MET<=250&&nB==1", basecuts_0lep ) ;
       fillSet( "h_ht_zl_nomctrigcorr_met2_nb1", "HT", cuts, false ) ;
       fillSet( "h_ht_zl_met2_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>250&&MET<=350&&nB==1", basecuts_0lep ) ;
       fillSet( "h_ht_zl_nomctrigcorr_met3_nb1", "HT", cuts, false ) ;
       fillSet( "h_ht_zl_met3_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>350&&nB==1", basecuts_0lep ) ;
       fillSet( "h_ht_zl_nomctrigcorr_met4_nb1", "HT", cuts, false ) ;
       fillSet( "h_ht_zl_met4_nb1", "HT", cuts ) ;

       sprintf( cuts, "%s&&MET>125&&MET<=150&&nB==0", basecuts_0lep_nonb ) ;
       fillSet( "h_ht_zl_nomctrigcorr_met1_nb0", "HT", cuts, false ) ;
       fillSet( "h_ht_zl_met1_nb0", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>150&&MET<=250&&nB==0", basecuts_0lep_nonb ) ;
       fillSet( "h_ht_zl_nomctrigcorr_met2_nb0", "HT", cuts, false ) ;
       fillSet( "h_ht_zl_met2_nb0", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>250&&MET<=350&&nB==0", basecuts_0lep_nonb ) ;
       fillSet( "h_ht_zl_nomctrigcorr_met3_nb0", "HT", cuts, false ) ;
       fillSet( "h_ht_zl_met3_nb0", "HT", cuts ) ;
       sprintf( cuts, "%s&&MET>350&&nB==0", basecuts_0lep_nonb ) ;
       fillSet( "h_ht_zl_nomctrigcorr_met4_nb0", "HT", cuts, false ) ;
       fillSet( "h_ht_zl_met4_nb0", "HT", cuts ) ;






       sprintf( cuts, "%s&&nB==0", basecuts_0lep_nonb ) ;
       fillSet( "h_met_zl_nb0", "MET", cuts ) ;
       fillSet( "h_ht_zl_nb0", "HT", cuts ) ;
       fillSet( "h_met_zl_nomctrigcorr_nb0", "MET", cuts, false ) ;
       fillSet( "h_ht_zl_nomctrigcorr_nb0", "HT", cuts, false ) ;
       fillSetLHB( "h_lhb_zl_nb0", cuts ) ;
       fillSetLHB( "h_lhb_zl_nomctrigcorr_nb0", cuts, false ) ;

       sprintf( cuts, "%s&&nB==1", basecuts_0lep ) ;
       fillSet( "h_met_zl_nb1", "MET", cuts ) ;
       fillSet( "h_ht_zl_nb1", "HT", cuts ) ;
       fillSet( "h_met_zl_nomctrigcorr_nb1", "MET", cuts, false ) ;
       fillSet( "h_ht_zl_nomctrigcorr_nb1", "HT", cuts, false ) ;
       fillSetLHB( "h_lhb_zl_nb1", cuts ) ;
       fillSetLHB( "h_lhb_zl_nomctrigcorr_nb1", cuts, false ) ;

       sprintf( cuts, "%s&&nB==2", basecuts_0lep ) ;
       fillSet( "h_met_zl_nb2", "MET", cuts ) ;
       fillSet( "h_ht_zl_nb2", "HT", cuts ) ;
       fillSet( "h_met_zl_nomctrigcorr_nb2", "MET", cuts, false ) ;
       fillSet( "h_ht_zl_nomctrigcorr_nb2", "HT", cuts, false ) ;
       fillSetLHB( "h_lhb_zl_nb2", cuts ) ;
       fillSetLHB( "h_lhb_zl_nomctrigcorr_nb2", cuts, false ) ;

       sprintf( cuts, "%s&&nB>=3", basecuts_0lep ) ;
       fillSet( "h_met_zl_nb3", "MET", cuts ) ;
       fillSet( "h_ht_zl_nb3", "HT", cuts ) ;
       fillSet( "h_met_zl_nomctrigcorr_nb3", "MET", cuts, false ) ;
       fillSet( "h_ht_zl_nomctrigcorr_nb3", "HT", cuts, false ) ;
       fillSetLHB( "h_lhb_zl_nb3", cuts ) ;
       fillSetLHB( "h_lhb_zl_nomctrigcorr_nb3", cuts, false ) ;






       sprintf( cuts, "%s", basecuts_1lep_nonjet ) ;
       fillSet( "h_njets_sl_all", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==0", basecuts_1lep_nonjet_nonb ) ;
       fillSet( "h_njets_sl_nb0", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_1lep_nonjet ) ;
       fillSet( "h_njets_sl_nb1", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_1lep_nonjet ) ;
       fillSet( "h_njets_sl_nb2", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==3", basecuts_1lep_nonjet ) ;
       fillSet( "h_njets_sl_nb3", "nJets", cuts ) ;

       sprintf( cuts, "%s", basecuts_ldp_nonjet ) ;
       fillSet( "h_njets_ldp_all", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==0", basecuts_ldp_nonjet_nonb ) ;
       fillSet( "h_njets_ldp_nb0", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_ldp_nonjet ) ;
       fillSet( "h_njets_ldp_nb1", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_ldp_nonjet ) ;
       fillSet( "h_njets_ldp_nb2", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==3", basecuts_ldp_nonjet ) ;
       fillSet( "h_njets_ldp_nb3", "nJets", cuts ) ;

       sprintf( cuts, "%s", basecuts_0lep_nonjet ) ;
       fillSet( "h_njets_zl_all", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==0", basecuts_0lep_nonjet_nonb ) ;
       fillSet( "h_njets_zl_nb0", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_0lep_nonjet ) ;
       fillSet( "h_njets_zl_nb1", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_0lep_nonjet ) ;
       fillSet( "h_njets_zl_nb2", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==3", basecuts_0lep_nonjet ) ;
       fillSet( "h_njets_zl_nb3", "nJets", cuts ) ;


     //-- MT

       sprintf( cuts, "%s", basecuts_1lep_nomt ) ;
       fillSet( "h_mt_sl_all", "MT", cuts ) ;
       sprintf( cuts, "%s&&nB==0", basecuts_1lep_nonb_nomt ) ;
       fillSet( "h_mt_sl_nb0", "MT", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_1lep_nonb_nomt ) ;
       fillSet( "h_mt_sl_nb1", "MT", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_1lep_nonb_nomt ) ;
       fillSet( "h_mt_sl_nb2", "MT", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_1lep_nonb_nomt ) ;
       fillSet( "h_mt_sl_nb3", "MT", cuts ) ;


     //-- NPV

       sprintf( cuts, "%s&&nB==0", basecuts_1lep_nonb ) ;
       fillSet( "h_npv_sl_nb0", "nPV", cuts ) ;
       sprintf( cuts, "%s&&nB==0", basecuts_ldp_nonb ) ;
       fillSet( "h_npv_ldp_nb0", "nPV", cuts ) ;
       sprintf( cuts, "%s&&nB==0", basecuts_0lep_nonb ) ;
       fillSet( "h_npv_zl_nb0", "nPV", cuts ) ;

       sprintf( cuts, "%s&&nB==1", basecuts_1lep_nonb ) ;
       fillSet( "h_npv_sl_nb1", "nPV", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_ldp_nonb ) ;
       fillSet( "h_npv_ldp_nb1", "nPV", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_0lep_nonb ) ;
       fillSet( "h_npv_zl_nb1", "nPV", cuts ) ;


      //-- min Delta phi

       sprintf( cuts, "%s&&nB==0", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphi_0lep_nb0", "minDelPhi", cuts ) ;

       sprintf( cuts, "%s&&nB==0&&MET>125&&MET<150", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphi_0lep_met1_nb0", "minDelPhi", cuts ) ;
       sprintf( cuts, "%s&&nB==0&&MET>150&&MET<250", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphi_0lep_met2_nb0", "minDelPhi", cuts ) ;
       sprintf( cuts, "%s&&nB==0&&MET>250&&MET<350", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphi_0lep_met3_nb0", "minDelPhi", cuts ) ;
       sprintf( cuts, "%s&&nB==0&&MET>350", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphi_0lep_met4_nb0", "minDelPhi", cuts ) ;

       sprintf( cuts, "%s&&nB==0&&HT>400&&HT<500", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphi_0lep_ht1_nb0", "minDelPhi", cuts ) ;
       sprintf( cuts, "%s&&nB==0&&HT>500&&HT<800", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphi_0lep_ht2_nb0", "minDelPhi", cuts ) ;
       sprintf( cuts, "%s&&nB==0&&HT>800&&HT<1000", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphi_0lep_ht3_nb0", "minDelPhi", cuts ) ;
       sprintf( cuts, "%s&&nB==0&&HT>1000", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphi_0lep_ht4_nb0", "minDelPhi", cuts ) ;


       sprintf( cuts, "%s&&nB==1", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphi_0lep_nb1", "minDelPhi", cuts ) ;

       sprintf( cuts, "%s&&nB==1&&MET>125&&MET<150", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphi_0lep_met1_nb1", "minDelPhi", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&MET>150&&MET<250", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphi_0lep_met2_nb1", "minDelPhi", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&MET>250&&MET<350", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphi_0lep_met3_nb1", "minDelPhi", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&MET>350", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphi_0lep_met4_nb1", "minDelPhi", cuts ) ;

       sprintf( cuts, "%s&&nB==1&&HT>400&&HT<500", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphi_0lep_ht1_nb1", "minDelPhi", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&HT>500&&HT<800", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphi_0lep_ht2_nb1", "minDelPhi", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&HT>800&&HT<1000", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphi_0lep_ht3_nb1", "minDelPhi", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&HT>1000", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphi_0lep_ht4_nb1", "minDelPhi", cuts ) ;





      //-- min Delta phi N

       sprintf( cuts, "%s&&nB==0", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphin_0lep_nb0", "minDelPhiN", cuts ) ;

       sprintf( cuts, "%s&&nB==0&&MET>125&&MET<150", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphin_0lep_met1_nb0", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==0&&MET>150&&MET<250", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphin_0lep_met2_nb0", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==0&&MET>250&&MET<350", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphin_0lep_met3_nb0", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==0&&MET>350", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphin_0lep_met4_nb0", "minDelPhiN", cuts ) ;

       sprintf( cuts, "%s&&nB==0&&HT>400&&HT<500", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphin_0lep_ht1_nb0", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==0&&HT>500&&HT<800", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphin_0lep_ht2_nb0", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==0&&HT>800&&HT<1000", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphin_0lep_ht3_nb0", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==0&&HT>1000", basecuts_0lep_nomindphin_nonb ) ;
       fillSet( "h_mindphin_0lep_ht4_nb0", "minDelPhiN", cuts ) ;


       sprintf( cuts, "%s&&nB==1", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_0lep_nb1", "minDelPhiN", cuts ) ;

       sprintf( cuts, "%s&&nB==1&&MET>125&&MET<150", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_0lep_met1_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&MET>150&&MET<250", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_0lep_met2_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&MET>250&&MET<350", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_0lep_met3_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&MET>350", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_0lep_met4_nb1", "minDelPhiN", cuts ) ;

       sprintf( cuts, "%s&&nB==1&&HT>400&&HT<500", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_0lep_ht1_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&HT>500&&HT<800", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_0lep_ht2_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&HT>800&&HT<1000", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_0lep_ht3_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&HT>1000", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_0lep_ht4_nb1", "minDelPhiN", cuts ) ;



       sprintf( cuts, "%s&&nB==1", basecuts_1lep_nomindphin ) ;
       fillSet( "h_mindphin_1lep_nb1", "minDelPhiN", cuts ) ;

       sprintf( cuts, "%s&&nB==1&&MET>125&&MET<150", basecuts_1lep_nomindphin ) ;
       fillSet( "h_mindphin_1lep_met1_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&MET>150&&MET<250", basecuts_1lep_nomindphin ) ;
       fillSet( "h_mindphin_1lep_met2_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&MET>250&&MET<350", basecuts_1lep_nomindphin ) ;
       fillSet( "h_mindphin_1lep_met3_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&MET>350", basecuts_1lep_nomindphin ) ;
       fillSet( "h_mindphin_1lep_met4_nb1", "minDelPhiN", cuts ) ;

       sprintf( cuts, "%s&&nB==1&&HT>400&&HT<500", basecuts_1lep_nomindphin ) ;
       fillSet( "h_mindphin_1lep_ht1_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&HT>500&&HT<800", basecuts_1lep_nomindphin ) ;
       fillSet( "h_mindphin_1lep_ht2_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&HT>800&&HT<1000", basecuts_1lep_nomindphin ) ;
       fillSet( "h_mindphin_1lep_ht3_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1&&HT>1000", basecuts_1lep_nomindphin ) ;
       fillSet( "h_mindphin_1lep_ht4_nb1", "minDelPhiN", cuts ) ;


     //--- met ratio

       sprintf( cuts, "%s", basecuts_ldp_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_ldp", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_ldp_wide", "pfOcaloMET", cuts ) ;

       sprintf( cuts, "%s&&nB==0", basecuts_ldp_nonb_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_ldp_nb0", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_ldp_nb0_wide", "pfOcaloMET", cuts ) ;

       sprintf( cuts, "%s&&nB==1", basecuts_ldp_nonb_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_ldp_nb1", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_ldp_nb1_wide", "pfOcaloMET", cuts ) ;

       sprintf( cuts, "%s&&nB==2", basecuts_ldp_nonb_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_ldp_nb2", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_ldp_nb2_wide", "pfOcaloMET", cuts ) ;

       sprintf( cuts, "%s&&nB>=3", basecuts_ldp_nonb_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_ldp_nb3", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_ldp_nb3_wide", "pfOcaloMET", cuts ) ;




       sprintf( cuts, "%s", basecuts_1lep_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_sl", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_sl_wide", "pfOcaloMET", cuts ) ;

       sprintf( cuts, "%s&&nB==0", basecuts_1lep_nonb_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_sl_nb0", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_sl_nb0_wide", "pfOcaloMET", cuts ) ;

       sprintf( cuts, "%s&&nB==1", basecuts_1lep_nonb_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_sl_nb1", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_sl_nb1_wide", "pfOcaloMET", cuts ) ;

       sprintf( cuts, "%s&&nB==2", basecuts_1lep_nonb_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_sl_nb2", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_sl_nb2_wide", "pfOcaloMET", cuts ) ;

       sprintf( cuts, "%s&&nB>=3", basecuts_1lep_nonb_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_sl_nb3", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_sl_nb3_wide", "pfOcaloMET", cuts ) ;





       sprintf( cuts, "%s&&nB==0", basecuts_0lep_nonb_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_zl_nb0", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_zl_nb0_wide", "pfOcaloMET", cuts ) ;

       sprintf( cuts, "%s&&nB==1", basecuts_0lep_nonb_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_zl_nb1", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_zl_nb1_wide", "pfOcaloMET", cuts ) ;

       sprintf( cuts, "%s&&nB==2", basecuts_0lep_nonb_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_zl_nb2", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_zl_nb2_wide", "pfOcaloMET", cuts ) ;

       sprintf( cuts, "%s&&nB>=3", basecuts_0lep_nonb_nometratiocut ) ;
       fillSet( "h_pfOcaloMET_zl_nb3", "pfOcaloMET", cuts ) ;
       fillSet( "h_pfOcaloMET_zl_nb3_wide", "pfOcaloMET", cuts ) ;






     //--- number of isolated tracks.

       sprintf( cuts, "%s&&nB==0", basecuts_0lep_nonb_noisotrk ) ;
       fillSet( "h_nisotrk_zl_nb0", "nIsoTrk", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_0lep_nonb_noisotrk ) ;
       fillSet( "h_nisotrk_zl_nb1", "nIsoTrk", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_0lep_nonb_noisotrk ) ;
       fillSet( "h_nisotrk_zl_nb2", "nIsoTrk", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_0lep_nonb_noisotrk ) ;
       fillSet( "h_nisotrk_zl_nb3", "nIsoTrk", cuts ) ;

       sprintf( cuts, "%s&&nB==0", basecuts_ldp_nonb_noisotrk ) ;
       fillSet( "h_nisotrk_ldp_nb0", "nIsoTrk", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_ldp_nonb_noisotrk ) ;
       fillSet( "h_nisotrk_ldp_nb1", "nIsoTrk", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_ldp_nonb_noisotrk ) ;
       fillSet( "h_nisotrk_ldp_nb2", "nIsoTrk", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_ldp_nonb_noisotrk ) ;
       fillSet( "h_nisotrk_ldp_nb3", "nIsoTrk", cuts ) ;

       sprintf( cuts, "%s&&nB==0", basecuts_1lep_nonb_noisotrk ) ;
       fillSet( "h_nisotrk_sl_nb0", "nIsoTrk", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_1lep_nonb_noisotrk ) ;
       fillSet( "h_nisotrk_sl_nb1", "nIsoTrk", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_1lep_nonb_noisotrk ) ;
       fillSet( "h_nisotrk_sl_nb2", "nIsoTrk", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_1lep_nonb_noisotrk ) ;
       fillSet( "h_nisotrk_sl_nb3", "nIsoTrk", cuts ) ;



     //--- max Ncharged - Nneutral for jets in |eta| from 0.9 to 1.9 with Pt>50.

       sprintf( cuts, "%s&&nB==0", basecuts_0lep_nonb_nomaxchmultdiff ) ;
       fillSet( "h_maxchnmultdiff_zl_nb0", "maxChNMultDiff", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_0lep_nonb_nomaxchmultdiff ) ;
       fillSet( "h_maxchnmultdiff_zl_nb1", "maxChNMultDiff", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_0lep_nonb_nomaxchmultdiff ) ;
       fillSet( "h_maxchnmultdiff_zl_nb2", "maxChNMultDiff", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_0lep_nonb_nomaxchmultdiff ) ;
       fillSet( "h_maxchnmultdiff_zl_nb3", "maxChNMultDiff", cuts ) ;

       sprintf( cuts, "%s&&nB==0", basecuts_ldp_nonb_nomaxchmultdiff ) ;
       fillSet( "h_maxchnmultdiff_ldp_nb0", "maxChNMultDiff", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_ldp_nonb_nomaxchmultdiff ) ;
       fillSet( "h_maxchnmultdiff_ldp_nb1", "maxChNMultDiff", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_ldp_nonb_nomaxchmultdiff ) ;
       fillSet( "h_maxchnmultdiff_ldp_nb2", "maxChNMultDiff", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_ldp_nonb_nomaxchmultdiff ) ;
       fillSet( "h_maxchnmultdiff_ldp_nb3", "maxChNMultDiff", cuts ) ;

       sprintf( cuts, "%s&&nB==0", basecuts_1lep_nonb_nomaxchmultdiff ) ;
       fillSet( "h_maxchnmultdiff_sl_nb0", "maxChNMultDiff", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_1lep_nonb_nomaxchmultdiff ) ;
       fillSet( "h_maxchnmultdiff_sl_nb1", "maxChNMultDiff", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_1lep_nonb_nomaxchmultdiff ) ;
       fillSet( "h_maxchnmultdiff_sl_nb2", "maxChNMultDiff", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_1lep_nonb_nomaxchmultdiff ) ;
       fillSet( "h_maxchnmultdiff_sl_nb3", "maxChNMultDiff", cuts ) ;





       char outfile[10000] ;
       sprintf( outfile, "rootfiles/dmc_plots2_%s.root", dataset_string ) ;
       saveHist( outfile, "h*" ) ;

    } // dmc_fillplots.

   //----------------------------


    void bookSetLHB( const char* hname_base, const char* htitle_base ) {

       printf("\n\n") ;
       for ( int ci=0; ci<nComps; ci++ ) {

          char hname[1000] ;
          char htitle[1000] ;

          sprintf( hname, "%s_%s", hname_base, compname[ci] ) ;
          sprintf( htitle, "%s, %s", htitle_base, compname[ci] ) ;

          printf(" bookSet: %s : %s\n", hname, htitle ) ;
          TH2F* hp = new TH2F( hname, htitle, nBinsMET, Mbins, nBinsHT, Hbins ) ;
          if ( compcolor[ci] > 0 ) { hp -> SetFillColor( compcolor[ci] ) ; }
          hp -> Sumw2() ;

       } // ci
       printf("\n\n") ;


    } // bookSetLHB

   //----------------------------


    void bookSet( const char* hname_base, const char* htitle_base, int nbins, double xmin, double xmax ) {

       printf("\n\n") ;
       for ( int ci=0; ci<nComps; ci++ ) {

          char hname[1000] ;
          char htitle[1000] ;

          sprintf( hname, "%s_%s", hname_base, compname[ci] ) ;
          sprintf( htitle, "%s, %s", htitle_base, compname[ci] ) ;

          printf(" bookSet: %s : %s\n", hname, htitle ) ;
          TH1F* hp = new TH1F( hname, htitle, nbins, xmin, xmax ) ;
          if ( compcolor[ci] > 0 ) { hp -> SetFillColor( compcolor[ci] ) ; }
          hp -> Sumw2() ;

       } // ci
       printf("\n\n") ;


    } // bookSet

   //----------------------------



    void fillSetLHB( const char* hname_base, const char* cuts, bool doTrigCorr ) {

       printf("\n\n") ;
       for ( int ci=0; ci<nComps; ci++ ) {

          char hname[1000] ;
          char arg1[10000] ;
          char arg2[10000] ;

          sprintf( hname, "%s_%s", hname_base, compname[ci] ) ;
          sprintf( arg1, "HT:MET>>%s", hname ) ;
          if ( doTrigCorr ) {
             sprintf( arg2, "trigWeight*(%s)", cuts ) ;
          } else {
             sprintf( arg2, "(%s)", cuts ) ;
          }
          printf(" %s : %s\n", arg1, arg2 ) ;
          compchain[ci] -> Draw( arg1, arg2 ) ;

          TH2F* hp = (TH2F*) gDirectory->FindObject( hname ) ;
          if ( hp == 0x0 ) { printf("\n *** missing hist! %s \n\n", hname ) ; return ; }

          hp -> Scale( compscale[ci] ) ;

          hp -> Draw("box") ;

          dcan->Update() ;

       } // ci
       printf("\n\n") ;


    } // fillSetLHB

   //----------------------------

    void fillSet( const char* hname_base, const char* varname, const char* cuts, bool doTrigCorr ) {

       printf("\n\n") ;
       for ( int ci=0; ci<nComps; ci++ ) {

          char hname[1000] ;
          char arg1[10000] ;
          char arg2[10000] ;

          sprintf( hname, "%s_%s", hname_base, compname[ci] ) ;
          sprintf( arg1, "%s>>%s", varname, hname ) ;
          if ( doTrigCorr ) {
             sprintf( arg2, "trigWeight*(%s)", cuts ) ;
          } else {
             sprintf( arg2, "(%s)", cuts ) ;
          }
          printf(" %s : %s\n", arg1, arg2 ) ;
          compchain[ci] -> Draw( arg1, arg2 ) ;

          TH1F* hp = (TH1F*) gDirectory->FindObject( hname ) ;
          if ( hp == 0x0 ) { printf("\n *** missing hist! %s \n\n", hname ) ; return ; }

          hp -> Scale( compscale[ci] ) ;

          int nbinsx = hp->GetNbinsX() ;
          double last_val     = hp->GetBinContent( nbinsx  ) ;
          double last_err     = hp->GetBinError( nbinsx  ) ;
          double overflow_val = hp->GetBinContent( nbinsx + 1 ) ;
          double overflow_err = hp->GetBinError( nbinsx + 1 ) ;

          hp -> SetBinContent( nbinsx, last_val + overflow_val ) ;
          hp -> SetBinError( nbinsx, sqrt( pow(last_err,2) + pow(overflow_err,2) ) ) ;


          hp -> Draw("hist") ;
          hp -> Draw("esame") ;

          dcan->Update() ;

       } // ci
       printf("\n\n") ;


    } // fillSet

   //----------------------------

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

