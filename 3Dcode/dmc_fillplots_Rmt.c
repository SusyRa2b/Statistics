
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

    int nComps(4) ;
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
    TH2F* evenBinLHB( TH2F* h2 ) ;

   //----------------------------

    void dmc_fillplots_Rmt( const char* dataset_string = "all" ) {

       gDirectory -> Delete("h*") ;

       gStyle -> SetOptStat(0) ;
       gStyle -> SetTitleH(0.035 ) ;
       gStyle -> SetMarkerSize(2) ;
       gStyle -> SetPaintTextFormat(".3f") ;

       for ( int ci=0; ci<nComps; ci++ ) { compchain[ci] = new TChain("tree") ; }

       int compind ;

  ////-- Data
  //   compind = 0 ;
  //   sprintf( compname[compind], "data" ) ;
  //   if ( strcmp(dataset_string, "all" ) == 0 ) {
  //      printf("\n\n Loading all data.\n\n") ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/MET_2012A.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/MET_2012B.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/MET_2012C_pr.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/MET_2012C_rr.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/HT_2012A.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/HTMHT_2012B.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/HTMHT_2012C_pr.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/HTMHT_2012C_rr.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/JetHT_2012B.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/JetHT_2012C_pr.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/JetHT_2012C_rr.root" ) ;
  //      complumi[compind] = 12.03 ;
  //   } else if ( strcmp(dataset_string, "RunsAB" ) == 0 ) {
  //      printf("\n\n Loading runs A and B.\n\n") ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/MET_2012A.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/MET_2012B.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/HT_2012A.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/HTMHT_2012B.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/JetHT_2012B.root" ) ;
  //      complumi[compind] = 5.23 ;
  //   } else if ( strcmp(dataset_string, "RunC" ) == 0 ) {
  //      printf("\n\n Loading run C.\n\n") ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/MET_2012C_pr.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/MET_2012C_rr.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/HTMHT_2012C_pr.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/HTMHT_2012C_rr.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/JetHT_2012C_pr.root" ) ;
  //      compchain[compind] -> Add( "filesHCP_53_v6/JetHT_2012C_rr.root" ) ;
  //      complumi[compind] = 6.81 ;
  //   } else {
  //      printf("\n\n *** Unknown dataset: %s\n\n", dataset_string ) ;
  //      return ;
  //   }
  // //------------
  //   compcolor[compind] = 0 ;
  //   compscale[compind] = 1. ;


       double MClumi = 12.03 ;

    //-- ttbar, single semilep
       compind = 0 ;
       sprintf( compname[compind], "ttsl" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/TT_SemiLept.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 600-7 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;


    //-- ttbar, double semilep
       compind = 1 ;
       sprintf( compname[compind], "ttdl" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/TT_FullLept.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 600-10 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;


    //-- single top
       compind = 2 ;
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
       compind = 3 ;
       sprintf( compname[compind], "wjets" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/WJets-250to300.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/WJets-300to400.root" ) ;
       compchain[compind] -> Add( "filesHCP_53_v6/WJets-400.root" ) ;
       complumi[compind] = MClumi ;
       compcolor[compind] = 600-4 ; // kBlue = 600
       compscale[compind] = complumi[0] / complumi[compind] ;






    //--- cuts

       char basecuts_1lep_nomt_nonb[10000] ;
       sprintf( basecuts_1lep_nomt_nonb, "maxChNMultDiff<40&&pfOcaloMET<2&&minDelPhiN>4&&( (nMu==1&&nEl==0) || (nMu==0&&nEl==1) )&&nJets>=3&&(pt_1st_leadJet>70&&pt_2nd_leadJet>70&&pt_3rd_leadJet>50)&&MET>125&&HT>400&&passedTrigger==1") ;




       dcan = (TCanvas*) gDirectory->FindObject("dcan") ;
       if ( dcan == 0x0 ) {
          dcan = new TCanvas("dcan","") ;
       }

       char htitle[1000] ;









   //-- likelihood bins

       sprintf( htitle, "HT vs MET, SL, MT<100, nB=0, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_mtl_nb0", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, MT<100, nB>0, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_mtl_nbgt0", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, MT<100, nB=1, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_mtl_nb1", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, MT<100, nB=2, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_mtl_nb2", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, MT<100, nB>=3, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_mtl_nb3", htitle ) ;


       sprintf( htitle, "HT vs MET, SL, MT>100, nB=0, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_mth_nb0", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, MT>100, nB>0, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_mth_nbgt0", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, MT>100, nB=1, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_mth_nb1", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, MT>100, nB=2, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_mth_nb2", htitle ) ;
       sprintf( htitle, "HT vs MET, SL, MT>100, nB>=3, %s", dataset_string ) ;
       bookSetLHB( "h_lhb_sl_mth_nb3", htitle ) ;






     //--- Fill the histograms. -----------------------------------------------

       char cuts[10000] ;


       sprintf( cuts, "%s&&MT<100&&nB==0", basecuts_1lep_nomt_nonb ) ;
       fillSetLHB( "h_lhb_sl_mtl_nb0", cuts ) ;

       sprintf( cuts, "%s&&MT<100&&nB>0", basecuts_1lep_nomt_nonb ) ;
       fillSetLHB( "h_lhb_sl_mtl_nbgt0", cuts ) ;

       sprintf( cuts, "%s&&MT<100&&nB==1", basecuts_1lep_nomt_nonb ) ;
       fillSetLHB( "h_lhb_sl_mtl_nb1", cuts ) ;

       sprintf( cuts, "%s&&MT<100&&nB==2", basecuts_1lep_nomt_nonb ) ;
       fillSetLHB( "h_lhb_sl_mtl_nb2", cuts ) ;

       sprintf( cuts, "%s&&MT<100&&nB>=3", basecuts_1lep_nomt_nonb ) ;
       fillSetLHB( "h_lhb_sl_mtl_nb3", cuts ) ;




       sprintf( cuts, "%s&&MT>100&&nB==0", basecuts_1lep_nomt_nonb ) ;
       fillSetLHB( "h_lhb_sl_mth_nb0", cuts ) ;

       sprintf( cuts, "%s&&MT>100&&nB>0", basecuts_1lep_nomt_nonb ) ;
       fillSetLHB( "h_lhb_sl_mth_nbgt0", cuts ) ;

       sprintf( cuts, "%s&&MT>100&&nB==1", basecuts_1lep_nomt_nonb ) ;
       fillSetLHB( "h_lhb_sl_mth_nb1", cuts ) ;

       sprintf( cuts, "%s&&MT>100&&nB==2", basecuts_1lep_nomt_nonb ) ;
       fillSetLHB( "h_lhb_sl_mth_nb2", cuts ) ;

       sprintf( cuts, "%s&&MT>100&&nB>=3", basecuts_1lep_nomt_nonb ) ;
       fillSetLHB( "h_lhb_sl_mth_nb3", cuts ) ;






       char outfile[10000] ;
       sprintf( outfile, "rootfiles/dmc_plots_Rmt_%s.root", dataset_string ) ;
       saveHist( outfile, "h*" ) ;

    } // dmc_fillplots_Rmt.

   //==========================================================================================================


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

   //==========================================================================================================


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

   //==========================================================================================================



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

          TH2F* heb = evenBinLHB( hp ) ;

          heb -> Draw("box") ;

          dcan->Update() ;

       } // ci
       printf("\n\n") ;


    } // fillSetLHB

   //==========================================================================================================


    TH2F* evenBinLHB( TH2F* h2 ) {

       char hname_eb[1000] ;
       sprintf( hname_eb, "%s_eb", h2->GetName() ) ;

       TH2F* heb = new TH2F( hname_eb, h2->GetTitle(), nBinsMET, 0.5, 0.5+nBinsMET,  nBinsHT, 0.5, 0.5+nBinsHT ) ;

       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          char binlabel[100] ;
          sprintf( binlabel, "MET%d", mbi+1 ) ;
          heb->GetXaxis()->SetBinLabel( mbi+1, binlabel ) ;
          for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
             sprintf( binlabel, "HT%d", hbi+1 ) ;
             heb->GetYaxis()->SetBinLabel( hbi+1, binlabel ) ;
             heb -> SetBinContent( mbi+1, hbi+1,  h2->GetBinContent( mbi+1, hbi+1 ) ) ;
             heb -> SetBinError  ( mbi+1, hbi+1,  h2->GetBinError  ( mbi+1, hbi+1 ) ) ;
          } // hbi
       } // mbi

       heb->SetLabelSize(0.04,"x") ;
       heb->SetFillColor( h2->GetFillColor() ) ;

       return heb ;

    } // evenBinLHB

   //==========================================================================================================


    TH2F* evenBinLHB( TH2F* h2 ) {

       char hname_eb[1000] ;
       sprintf( hname_eb, "%s_eb", h2->GetName() ) ;

       TH2F* heb = new TH2F( hname_eb, h2->GetTitle(), nBinsMET, 0.5, 0.5+nBinsMET,  nBinsHT, 0.5, 0.5+nBinsHT ) ;

       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          char binlabel[100] ;
          sprintf( binlabel, "MET%d", mbi+1 ) ;
          heb->GetXaxis()->SetBinLabel( mbi+1, binlabel ) ;
          for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
             sprintf( binlabel, "HT%d", hbi+1 ) ;
             heb->GetYaxis()->SetBinLabel( hbi+1, binlabel ) ;
             heb -> SetBinContent( mbi+1, hbi+1,  h2->GetBinContent( mbi+1, hbi+1 ) ) ;
             heb -> SetBinError  ( mbi+1, hbi+1,  h2->GetBinError  ( mbi+1, hbi+1 ) ) ;
          } // hbi
       } // mbi

       heb->SetLabelSize(0.04,"x") ;
       heb->SetFillColor( h2->GetFillColor() ) ;

       return heb ;

    } // evenBinLHB


   //==========================================================================================================
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

