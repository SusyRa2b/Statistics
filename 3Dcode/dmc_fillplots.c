
#include "TStyle.h"
#include "TH1F.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TFile.h"
#include <iostream>


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
    void fillSet( const char* hname_base, const char* varname, const char* cuts ) ;
    void saveHist(const char* filename, const char* pat) ;

   //----------------------------

    void dmc_fillplots( const char* dataset_string = "all" ) {

       gDirectory -> Delete("h*") ;

       gStyle -> SetOptStat(0) ;

       for ( int ci=0; ci<nComps; ci++ ) { compchain[ci] = new TChain("tree") ; }

       int compind ;

    //-- Data
       compind = 0 ;
       sprintf( compname[compind], "data" ) ;
  //------------
  //   compchain[compind] -> Add( "filesHCP_53_v1/MET_2012A_BLIND.root" ) ;
  //   compchain[compind] -> Add( "filesHCP_53_v1/MET_2012B_BLIND.root" ) ;
  //   compchain[compind] -> Add( "filesHCP_53_v1/HT_2012A_BLIND.root" ) ;
  //   compchain[compind] -> Add( "filesHCP_53_v1/HTMHT_2012B_BLIND.root" ) ;
  //   complumi[compind] = 5.295 ; // AB
  //   complumi[compind] = 4.421 ; // B
  //------------
  //   compchain[compind] -> Add( "filesHCP_53_v1/MET_2012C_pr_BLIND.root" ) ;
  //   compchain[compind] -> Add( "filesHCP_53_v1/MET_2012C_rr_BLIND.root" ) ;
  //   compchain[compind] -> Add( "filesHCP_53_v1/HTMHT_2012C_pr_BLIND.root" ) ;
  //   compchain[compind] -> Add( "filesHCP_53_v1/HTMHT_2012C_rr_BLIND.root" ) ;
  //   complumi[compind] = 6.806 ;
  //------------
       compchain[compind] -> Add( "filesHCP_53_v1/MET_2012A_BLIND.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/MET_2012B_BLIND.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/MET_2012C_pr_BLIND.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/MET_2012C_rr_BLIND.root" ) ;
       complumi[compind] = 12.117 ;
  //------------
       compcolor[compind] = 0 ;
       compscale[compind] = 1. ;


    //-- Diboson.
       compind = 1 ;
       sprintf( compname[compind], "diboson" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/WW.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/WZ.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/ZZ.root" ) ;
       complumi[compind] = 11.9 ;
       compcolor[compind] = 800+7 ; // kOrange = 800
       compscale[compind] = complumi[0] / complumi[compind] ;





    //-- Zinvisible
       compind = 2 ;
       sprintf( compname[compind], "znn" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/Zinv-400.root" ) ;
       complumi[compind] = 11.9 ;
       compcolor[compind] = 416-3 ; // kGreen = 416
       compscale[compind] = complumi[0] / complumi[compind] ;




    //-- QCD
       compind = 3 ;
       sprintf( compname[compind], "qcd" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-170to300.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-300to470.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-470to600.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-600to800.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-800to1000.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-1000to1400.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-1400to1800.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/QCD-1800.root" ) ;
       complumi[compind] = 11.9 ;
       compcolor[compind] = 2 ;
       compscale[compind] = complumi[0] / complumi[compind] ;




    //-- single top
       compind = 4 ;
       sprintf( compname[compind], "singlet" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/T-s.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/Tbar-s.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/T-tW.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/Tbar-tW.root" ) ;
       complumi[compind] = 11.9 ;
       compcolor[compind] = 600-9 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;



    //-- W+jets
       compind = 5 ;
       sprintf( compname[compind], "wjets" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/WJets-250to300.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/WJets-300to400.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/WJets-400.root" ) ;
       complumi[compind] = 11.9 ;
       compcolor[compind] = 600-4 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;



    //-- ttbar
       compind = 6 ;
       sprintf( compname[compind], "ttbar" ) ;
       compchain[compind] -> Add( "filesHCP_53_v1/TT.root" ) ;
       complumi[compind] = 11.9 ;
       compcolor[compind] = 600-7 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;






       char basecuts_1lep[10000] ;
    // sprintf( basecuts_1lep, "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;
       sprintf( basecuts_1lep, "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400") ;

       char basecuts_1lep_nonb[10000] ;
    // sprintf( basecuts_1lep_nonb, "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;
       sprintf( basecuts_1lep_nonb, "minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400") ;

       char basecuts_ldp[10000] ;
    // sprintf( basecuts_ldp, "minDelPhiN<=4&&(nMu==0&&nEl==0)&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;
       sprintf( basecuts_ldp, "minDelPhiN<=4&&(nMu==0&&nEl==0)&&nB>0&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400") ;

       char basecuts_ldp_nonb[10000] ;
    // sprintf( basecuts_ldp_nonb, "minDelPhiN<=4&&(nMu==0&&nEl==0)&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;
       sprintf( basecuts_ldp_nonb, "minDelPhiN<=4&&(nMu==0&&nEl==0)&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400") ;



       dcan = (TCanvas*) gDirectory->FindObject("dcan") ;
       if ( dcan == 0x0 ) {
          dcan = new TCanvas("dcan","") ;
       }

       char htitle[1000] ;


       sprintf( htitle, "MET, SL, %s", dataset_string ) ;
       bookSet( "h_met_sl_all", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_sl_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_sl_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, SL, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_sl_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "HT, SL, %s", dataset_string ) ;
       bookSet( "h_ht_sl_all", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, SL, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_sl_nb3", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "nbjets, SL, %s", dataset_string ) ;
       bookSet( "h_nb_sl_all", htitle, 6, -0.5, 5.5 ) ;


       sprintf( htitle, "MET, LDP, %s", dataset_string ) ;
       bookSet( "h_met_ldp_all", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, nB=1, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nb1", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, nB=2, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nb2", htitle, 30, 125., 500. ) ;
       sprintf( htitle, "MET, LDP, nB>=3, %s", dataset_string ) ;
       bookSet( "h_met_ldp_nb3", htitle, 30, 125., 500. ) ;

       sprintf( htitle, "HT, LDP, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_all", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB=1, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nb1", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB=2, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nb2", htitle, 32, 400., 2000. ) ;
       sprintf( htitle, "HT, LDP, nB>=3, %s", dataset_string ) ;
       bookSet( "h_ht_ldp_nb3", htitle, 32, 400., 2000. ) ;

       sprintf( htitle, "nbjets, LDP, %s", dataset_string ) ;
       bookSet( "h_nb_ldp_all", htitle, 6, -0.5, 5.5 ) ;






       char cuts[10000] ;

       sprintf( cuts, "%s", basecuts_1lep ) ;
       fillSet( "h_met_sl_all", "MET", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_1lep ) ;
       fillSet( "h_met_sl_nb1", "MET", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_1lep ) ;
       fillSet( "h_met_sl_nb2", "MET", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_1lep ) ;
       fillSet( "h_met_sl_nb3", "MET", cuts ) ;

       sprintf( cuts, "%s", basecuts_1lep ) ;
       fillSet( "h_ht_sl_all", "HT", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_1lep ) ;
       fillSet( "h_ht_sl_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_1lep ) ;
       fillSet( "h_ht_sl_nb2", "HT", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_1lep ) ;
       fillSet( "h_ht_sl_nb3", "HT", cuts ) ;

       sprintf( cuts, "%s", basecuts_1lep_nonb ) ;
       fillSet( "h_nb_sl_all", "nB", cuts ) ;



       sprintf( cuts, "%s", basecuts_ldp ) ;
       fillSet( "h_met_ldp_all", "MET", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_ldp ) ;
       fillSet( "h_met_ldp_nb1", "MET", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_ldp ) ;
       fillSet( "h_met_ldp_nb2", "MET", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_ldp ) ;
       fillSet( "h_met_ldp_nb3", "MET", cuts ) ;

       sprintf( cuts, "%s", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_all", "HT", cuts ) ;
       sprintf( cuts, "%s&&nB==1", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_nb1", "HT", cuts ) ;
       sprintf( cuts, "%s&&nB==2", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_nb2", "HT", cuts ) ;
       sprintf( cuts, "%s&&nB>=3", basecuts_ldp ) ;
       fillSet( "h_ht_ldp_nb3", "HT", cuts ) ;

       sprintf( cuts, "%s", basecuts_ldp_nonb ) ;
       fillSet( "h_nb_ldp_all", "nB", cuts ) ;






       char outfile[10000] ;
       sprintf( outfile, "rootfiles/dmc_plots1_%s.root", dataset_string ) ;
       saveHist( outfile, "h*" ) ;

    } // dmc_fillplots.

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



    void fillSet( const char* hname_base, const char* varname, const char* cuts ) {

       printf("\n\n") ;
       for ( int ci=0; ci<nComps; ci++ ) {

          char hname[1000] ;
          char arg1[10000] ;

          sprintf( hname, "%s_%s", hname_base, compname[ci] ) ;
          sprintf( arg1, "%s>>%s", varname, hname ) ;
          printf(" %s : %s\n", arg1, cuts ) ;
          compchain[ci] -> Draw( arg1, cuts ) ;

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

