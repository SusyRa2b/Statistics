
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

       for ( int ci=0; ci<nComps; ci++ ) { compchain[ci] = new TChain("tree") ; }

       int compind ;

    //-- Data
       compind = 0 ;
       sprintf( compname[compind], "data" ) ;
       if ( strcmp(dataset_string, "all" ) == 0 ) {
          printf("\n\n Loading all data.\n\n") ;
          compchain[compind] -> Add( "filesHCP_53_v1/MET_2012A_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/MET_2012B_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/MET_2012C_pr_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/MET_2012C_rr_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/HT_2012A_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/HTMHT_2012B_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/HTMHT_2012C_pr_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/HTMHT_2012C_rr_BLIND.root" ) ;
          complumi[compind] = 12.03 ;
       } else if ( strcmp(dataset_string, "RunsAB" ) == 0 ) {
          printf("\n\n Loading runs A and B.\n\n") ;
          compchain[compind] -> Add( "filesHCP_53_v1/MET_2012A_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/MET_2012B_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/HT_2012A_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/HTMHT_2012B_BLIND.root" ) ;
          complumi[compind] = 5.23 ;
       } else if ( strcmp(dataset_string, "RunC" ) == 0 ) {
          printf("\n\n Loading run C.\n\n") ;
          compchain[compind] -> Add( "filesHCP_53_v1/MET_2012C_pr_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/MET_2012C_rr_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/HTMHT_2012C_pr_BLIND.root" ) ;
          compchain[compind] -> Add( "filesHCP_53_v1/HTMHT_2012C_rr_BLIND.root" ) ;
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
       compchain[compind] -> Add( "filesHCP_53_v1/WW.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/WZ.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/ZZ.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 800+7 ; // kOrange = 800
       compscale[compind] = complumi[0] / complumi[compind] ;





    //-- Zinvisible
       compind = 2 ;
       sprintf( compname[compind], "znn" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/Zinv-100to200.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/Zinv-200to400.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/Zinv-400.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 416-3 ; // kGreen = 416
       compscale[compind] = complumi[0] / complumi[compind] ;




    //-- QCD
       compind = 3 ;
       sprintf( compname[compind], "qcd" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-120to170.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-170to300.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-300to470.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-470to600.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-600to800.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-800to1000.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-1000to1400.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-1400to1800.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-1800.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 2 ;
       compscale[compind] = complumi[0] / complumi[compind] ;




    //-- single top
       compind = 4 ;
       sprintf( compname[compind], "singlet" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/T-t.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/Tbar-t.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/T-s.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/Tbar-s.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/T-tW.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/Tbar-tW.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 600-9 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;



    //-- W+jets
       compind = 5 ;
       sprintf( compname[compind], "wjets" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/WJets-250to300.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/WJets-300to400.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/WJets-400.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 600-4 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;



    //-- ttbar
       compind = 6 ;
       sprintf( compname[compind], "ttbar" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/TT.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 600-7 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;




    //--- cuts

       char basecuts_1lep[10000] ;
       sprintf( basecuts_1lep, "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nonb[10000] ;
       sprintf( basecuts_1lep_nonb, "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_1lep_nonjet[10000] ;
       sprintf( basecuts_1lep_nonjet, "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nB>0&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_ldp[10000] ;
       sprintf( basecuts_ldp, "minDelPhiN<=4&&(nMu==0&&nEl==0)&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_ldp_nonb[10000] ;
       sprintf( basecuts_ldp_nonb, "minDelPhiN<=4&&(nMu==0&&nEl==0)&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_ldp_nonjet[10000] ;
       sprintf( basecuts_ldp_nonjet, "minDelPhiN<=4&&(nMu==0&&nEl==0)&&nB>0&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;

       char basecuts_0lep_nomindphin[10000] ;
       sprintf( basecuts_0lep_nomindphin, "(nMu==0&&nEl==0)&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;



       dcan = (TCanvas*) gDirectory->FindObject("dcan") ;
       if ( dcan == 0x0 ) {
          dcan = new TCanvas("dcan","") ;
       }

       char htitle[1000] ;

    //-- MET

       sprintf( htitle, "MET, SL, %s", dataset_string ) ;
       bookSet( "h_met_sl_all", htitle, 30, 125., 500. ) ;
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
       sprintf( htitle, "MET, SL, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, no MC trig correction, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_sl_nomctrigcorr_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "MET, LDP, %s", dataset_string ) ;
       bookSet( "h_met_ldp_all", htitle, 30, 125., 500. ) ;
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
       sprintf( htitle, "MET, LDP, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, no MC trig correction, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nomctrigcorr_nb3", htitle, 30, 125., 500. ) ;



     //-- HT

       sprintf( htitle, "HT, SL, %s", dataset_string ) ;
       bookSet( "h_ht_sl_all", htitle, 32, 400., 2000. ) ;
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
       sprintf( htitle, "HT, SL, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, no MC trig correction, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nomctrigcorr_nb3", htitle, 32, 400., 2000. ) ;



       sprintf( htitle, "HT, LDP, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_all", htitle, 32, 400., 2000. ) ;
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
       sprintf( htitle, "HT, LDP, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, no MC trig correction, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nomctrigcorr_nb3", htitle, 32, 400., 2000. ) ;



       sprintf( htitle, "HT vs MET, SL, nB=1, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nb1", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, nB=2, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nb2", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, nB=3, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nb3", htitle ) ;

       sprintf( htitle, "HT vs MET, LDP, nB=1, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nb1", htitle ) ;
       sprintf( htitle, "HT vs MET, LDP, nB=2, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nb2", htitle ) ;
       sprintf( htitle, "HT vs MET, LDP, nB=3, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nb3", htitle ) ;


       sprintf( htitle, "HT vs MET, SL, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nomctrigcorr_nb1", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nomctrigcorr_nb2", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, no MC trig correction, nB=3, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_nomctrigcorr_nb3", htitle ) ;

       sprintf( htitle, "HT vs MET, LDP, no MC trig correction, nB=1, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nomctrigcorr_nb1", htitle ) ;
       sprintf( htitle, "HT vs MET, LDP, no MC trig correction, nB=2, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nomctrigcorr_nb2", htitle ) ;
       sprintf( htitle, "HT vs MET, LDP, no MC trig correction, nB=3, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_ldp_nomctrigcorr_nb3", htitle ) ;





       sprintf( htitle, "nbjets, SL, %s", dataset_string ) ;
       bookSet( "h_nb_sl_all", htitle, 6, -0.5, 5.5 ) ;

       sprintf( htitle, "nbjets, LDP, %s", dataset_string ) ;
       bookSet( "h_nb_ldp_all", htitle, 6, -0.5, 5.5 ) ;


       sprintf( htitle, "minDphiN, no leptons, %s", dataset_string ) ;
       bookSet( "h_mindphin_nolep_all", htitle, 30, 0., 15. ) ;
       sprintf( htitle, "minDphiN, no leptons, nB=1, %s", dataset_string ) ;
       bookSet( "h_mindphin_nolep_nb1", htitle, 30, 0., 15. ) ;
       sprintf( htitle, "minDphiN, no leptons, nB=2, %s", dataset_string ) ;
       bookSet( "h_mindphin_nolep_nb2", htitle, 30, 0., 15. ) ;
       sprintf( htitle, "minDphiN, no leptons, nB=3, %s", dataset_string ) ;
       bookSet( "h_mindphin_nolep_nb3", htitle, 30, 0., 15. ) ;


       sprintf( htitle, "N jets, SL, %s", dataset_string ) ;
       bookSet( "h_njets_sl_all", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, SL, nB=1, %s", dataset_string ) ;
       bookSet( "h_njets_sl_nb1", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, SL, nB=2, %s", dataset_string ) ;
       bookSet( "h_njets_sl_nb2", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, SL, nB=3, %s", dataset_string ) ;
       bookSet( "h_njets_sl_nb3", htitle, 11, -0.5, 10.5 ) ;

       sprintf( htitle, "N jets, LDP, %s", dataset_string ) ;
       bookSet( "h_njets_ldp_all", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, LDP, nB=1, %s", dataset_string ) ;
       bookSet( "h_njets_ldp_nb1", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, LDP, nB=2, %s", dataset_string ) ;
       bookSet( "h_njets_ldp_nb2", htitle, 11, -0.5, 10.5 ) ;
       sprintf( htitle, "N jets, LDP, nB=3, %s", dataset_string ) ;
       bookSet( "h_njets_ldp_nb3", htitle, 11, -0.5, 10.5 ) ;



     //--- Fill the histograms. -----------------------------------------------

       char cuts[10000] ;

       sprintf( cuts, "%s", basecuts_1lep_nonb ) ;
       fillSet( "h_nb_sl_all", "nB", cuts ) ;

       sprintf( cuts, "%s", basecuts_ldp_nonb ) ;
       fillSet( "h_nb_ldp_all", "nB", cuts ) ;


       sprintf( cuts, "%s", basecuts_1lep ) ;
       fillSet( "h_met_sl_all", "MET", cuts ) ;
       fillSet( "h_ht_sl_all", "HT", cuts ) ;
       fillSet( "h_met_sl_nomctrigcorr_all", "MET", cuts, false ) ;
       fillSet( "h_ht_sl_nomctrigcorr_all", "HT", cuts, false ) ;

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














       sprintf( cuts, "%s", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_nolep_all", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_nolep_nb1", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_nolep_nb2", "minDelPhiN", cuts ) ;
       sprintf( cuts, "%s&&nB==3", basecuts_0lep_nomindphin ) ;
       fillSet( "h_mindphin_nolep_nb3", "minDelPhiN", cuts ) ;

       sprintf( cuts, "%s", basecuts_1lep_nonjet ) ;
       fillSet( "h_njets_sl_all", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_1lep_nonjet ) ;
       fillSet( "h_njets_sl_nb1", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_1lep_nonjet ) ;
       fillSet( "h_njets_sl_nb2", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==3", basecuts_1lep_nonjet ) ;
       fillSet( "h_njets_sl_nb3", "nJets", cuts ) ;

       sprintf( cuts, "%s", basecuts_ldp_nonjet ) ;
       fillSet( "h_njets_ldp_all", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_ldp_nonjet ) ;
       fillSet( "h_njets_ldp_nb1", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_ldp_nonjet ) ;
       fillSet( "h_njets_ldp_nb2", "nJets", cuts ) ;
       sprintf( cuts, "%s&&nB==3", basecuts_ldp_nonjet ) ;
       fillSet( "h_njets_ldp_nb3", "nJets", cuts ) ;




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

