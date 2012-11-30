#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TText.h"
#include "TString.h"
#include "TRegexp.h"

#include <iostream>


  //-- met4-ht4-v15
      const int nBinsMET   = 4 ;
      const int nBinsHT    = 4 ;
      float Mbins[nBinsMET+1] = {125.,150.,250.,350.,99999.};
      float Hbins[nBinsHT+1] = {400.,500.,800.,1000.,99999.};

    bool recycleCanvas(false) ;
    bool savePdf ;

    int nComps(4) ;
    char compname[4][100] = { "ttsl", "ttdl", "singlet", "wjets" } ;

    TH2F* combineCompsLHB( const char* hname_base ) ;
    TH2F* combineToponlyLHB( const char* hname_base ) ;

   TH2F* makeSFplotLHB( const char* selname, const char* compstring, float& ave_Rmt ) ;


   void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;


   //------------

   void dmc_drawplots_Rmc( const char* infile = "rootfiles/dmc_plots_Rmt_all.root",
                           bool arg_savePdf = false,
                           bool arg_recycleCanvas = false ) {

       loadHist( infile ) ;

       savePdf = arg_savePdf ;
       recycleCanvas = arg_recycleCanvas ;

       gStyle -> SetOptStat(0) ;
       gStyle -> SetTitleH(0.045 ) ;
       gStyle -> SetTitleW(0.85 ) ;
       gStyle -> SetMarkerSize(2.5) ;
       gStyle -> SetPadRightMargin(0.15) ;
       gStyle -> SetPaintTextFormat(".2f") ;
       gStyle -> SetLabelSize(0.06,"x") ;
       gStyle -> SetLabelSize(0.06,"y") ;


  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       float ave_Rmt(-1.) ;

       TH2F* h_sf_Rmt_nbgt0 = makeSFplotLHB( "nbgt0", "comb", ave_Rmt ) ; if ( h_sf_Rmt_nbgt0 == 0x0 ) return ;
       TH2F* h_sf_Rmt_nb0   = makeSFplotLHB( "nb0"  , "comb", ave_Rmt ) ; if ( h_sf_Rmt_nb0   == 0x0 ) return ;
       TH2F* h_sf_Rmt_nb1   = makeSFplotLHB( "nb1"  , "comb", ave_Rmt ) ; if ( h_sf_Rmt_nb1   == 0x0 ) return ;
       TH2F* h_sf_Rmt_nb2   = makeSFplotLHB( "nb2"  , "comb", ave_Rmt ) ; if ( h_sf_Rmt_nb2   == 0x0 ) return ;
       TH2F* h_sf_Rmt_nb3   = makeSFplotLHB( "nb3"  , "comb", ave_Rmt ) ; if ( h_sf_Rmt_nb3   == 0x0 ) return ;

       TH2F* h_Rmt_nbgt0 = (TH2F*) gDirectory->FindObject( "h_Rmt_nbgt0_comb"  ) ;
       TH2F* h_Rmt_nb0   = (TH2F*) gDirectory->FindObject( "h_Rmt_nb0_comb"    ) ;
       TH2F* h_Rmt_nb1   = (TH2F*) gDirectory->FindObject( "h_Rmt_nb1_comb"    ) ;
       TH2F* h_Rmt_nb2   = (TH2F*) gDirectory->FindObject( "h_Rmt_nb2_comb"    ) ;
       TH2F* h_Rmt_nb3   = (TH2F*) gDirectory->FindObject( "h_Rmt_nb3_comb"    ) ;

       h_sf_Rmt_nbgt0 -> SetMaximum(2.5) ;
       h_sf_Rmt_nb0   -> SetMaximum(2.5) ;
       h_sf_Rmt_nb1   -> SetMaximum(2.5) ;
       h_sf_Rmt_nb2   -> SetMaximum(2.5) ;
       h_sf_Rmt_nb3   -> SetMaximum(2.5) ;

       h_sf_Rmt_nbgt0 -> SetMinimum(0.3) ;
       h_sf_Rmt_nb0   -> SetMinimum(0.3) ;
       h_sf_Rmt_nb1   -> SetMinimum(0.3) ;
       h_sf_Rmt_nb2   -> SetMinimum(0.3) ;
       h_sf_Rmt_nb3   -> SetMinimum(0.3) ;

       h_Rmt_nbgt0 -> SetMaximum(0.3) ;
       h_Rmt_nb0   -> SetMaximum(0.3) ;
       h_Rmt_nb1   -> SetMaximum(0.3) ;
       h_Rmt_nb2   -> SetMaximum(0.3) ;
       h_Rmt_nb3   -> SetMaximum(0.3) ;

       h_Rmt_nbgt0 -> SetMinimum(0.00) ;
       h_Rmt_nb0   -> SetMinimum(0.00) ;
       h_Rmt_nb1   -> SetMinimum(0.00) ;
       h_Rmt_nb2   -> SetMinimum(0.00) ;
       h_Rmt_nb3   -> SetMinimum(0.00) ;

     //---------

       gStyle -> SetPaintTextFormat(".2f") ;

       TCanvas* c_sf_Rmt = new TCanvas( "c_sf_Rmt","c_sf_Rmt", 400, 1000 ) ;
       c_sf_Rmt->Clear() ;
     //c_sf_Rmt->Divide(1,4) ;
       c_sf_Rmt->Divide(1,3) ;


     //c_sf_Rmt->cd(1) ;
     //h_sf_Rmt_nb0->Draw("colz") ;
     //h_sf_Rmt_nb0->Draw("texte same") ;

       c_sf_Rmt->cd(1) ;
       h_sf_Rmt_nb1->Draw("colz") ;
       h_sf_Rmt_nb1->Draw("texte same") ;

       c_sf_Rmt->cd(2) ;
       h_sf_Rmt_nb2->Draw("colz") ;
       h_sf_Rmt_nb2->Draw("texte same") ;

       c_sf_Rmt->cd(3) ;
       h_sf_Rmt_nb3->Draw("colz") ;
       h_sf_Rmt_nb3->Draw("texte same") ;

       c_sf_Rmt->Update() ; c_sf_Rmt->Draw() ;


     //---------

       gStyle -> SetPaintTextFormat(".3f") ;

       TCanvas* c_Rmt = new TCanvas( "c_Rmt","c_Rmt", 400, 1000 ) ;
       c_Rmt->Clear() ;
     //c_Rmt->Divide(1,4) ;
       c_Rmt->Divide(1,3) ;


     //c_Rmt->cd(1) ;
     //h_Rmt_nb0->Draw("colz") ;
     //h_Rmt_nb0->Draw("texte same") ;

       c_Rmt->cd(1) ;
       h_Rmt_nb1->Draw("colz") ;
       h_Rmt_nb1->Draw("texte same") ;

       c_Rmt->cd(2) ;
       h_Rmt_nb2->Draw("colz") ;
       h_Rmt_nb2->Draw("texte same") ;

       c_Rmt->cd(3) ;
       h_Rmt_nb3->Draw("colz") ;
       h_Rmt_nb3->Draw("texte same") ;


  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       float ave_Rmt_toponly(-1.) ;

       TH2F* h_sf_Rmt_nbgt0_toponly = makeSFplotLHB( "nbgt0", "toponly", ave_Rmt_toponly ) ; if ( h_sf_Rmt_nbgt0_toponly == 0x0 ) return ;
       TH2F* h_sf_Rmt_nb0_toponly   = makeSFplotLHB( "nb0"  , "toponly", ave_Rmt_toponly ) ; if ( h_sf_Rmt_nb0_toponly   == 0x0 ) return ;
       TH2F* h_sf_Rmt_nb1_toponly   = makeSFplotLHB( "nb1"  , "toponly", ave_Rmt_toponly ) ; if ( h_sf_Rmt_nb1_toponly   == 0x0 ) return ;
       TH2F* h_sf_Rmt_nb2_toponly   = makeSFplotLHB( "nb2"  , "toponly", ave_Rmt_toponly ) ; if ( h_sf_Rmt_nb2_toponly   == 0x0 ) return ;
       TH2F* h_sf_Rmt_nb3_toponly   = makeSFplotLHB( "nb3"  , "toponly", ave_Rmt_toponly ) ; if ( h_sf_Rmt_nb3_toponly   == 0x0 ) return ;

       TH2F* h_Rmt_nbgt0_toponly = (TH2F*) gDirectory->FindObject( "h_Rmt_nbgt0_toponly"  ) ;
       TH2F* h_Rmt_nb0_toponly   = (TH2F*) gDirectory->FindObject( "h_Rmt_nb0_toponly"    ) ;
       TH2F* h_Rmt_nb1_toponly   = (TH2F*) gDirectory->FindObject( "h_Rmt_nb1_toponly"    ) ;
       TH2F* h_Rmt_nb2_toponly   = (TH2F*) gDirectory->FindObject( "h_Rmt_nb2_toponly"    ) ;
       TH2F* h_Rmt_nb3_toponly   = (TH2F*) gDirectory->FindObject( "h_Rmt_nb3_toponly"    ) ;

       h_sf_Rmt_nbgt0_toponly -> SetMaximum(2.5) ;
       h_sf_Rmt_nb0_toponly   -> SetMaximum(2.5) ;
       h_sf_Rmt_nb1_toponly   -> SetMaximum(2.5) ;
       h_sf_Rmt_nb2_toponly   -> SetMaximum(2.5) ;
       h_sf_Rmt_nb3_toponly   -> SetMaximum(2.5) ;

       h_sf_Rmt_nbgt0_toponly -> SetMinimum(0.3) ;
       h_sf_Rmt_nb0_toponly   -> SetMinimum(0.3) ;
       h_sf_Rmt_nb1_toponly   -> SetMinimum(0.3) ;
       h_sf_Rmt_nb2_toponly   -> SetMinimum(0.3) ;
       h_sf_Rmt_nb3_toponly   -> SetMinimum(0.3) ;

       h_Rmt_nbgt0_toponly -> SetMaximum(0.3) ;
       h_Rmt_nb0_toponly   -> SetMaximum(0.3) ;
       h_Rmt_nb1_toponly   -> SetMaximum(0.3) ;
       h_Rmt_nb2_toponly   -> SetMaximum(0.3) ;
       h_Rmt_nb3_toponly   -> SetMaximum(0.3) ;

       h_Rmt_nbgt0_toponly -> SetMinimum(0.00) ;
       h_Rmt_nb0_toponly   -> SetMinimum(0.00) ;
       h_Rmt_nb1_toponly   -> SetMinimum(0.00) ;
       h_Rmt_nb2_toponly   -> SetMinimum(0.00) ;
       h_Rmt_nb3_toponly   -> SetMinimum(0.00) ;

     //---------

       gStyle -> SetPaintTextFormat(".2f") ;

       TCanvas* c_sf_Rmt_toponly = new TCanvas( "c_sf_Rmt_toponly","c_sf_Rmt_toponly", 400, 1000 ) ;
       c_sf_Rmt_toponly->Clear() ;
     //c_sf_Rmt_toponly->Divide(1,4) ;
       c_sf_Rmt_toponly->Divide(1,3) ;


     //c_sf_Rmt_toponly->cd(1) ;
     //h_sf_Rmt_nb0_toponly->Draw("colz") ;
     //h_sf_Rmt_nb0_toponly->Draw("texte same") ;

       c_sf_Rmt_toponly->cd(1) ;
       h_sf_Rmt_nb1_toponly->Draw("colz") ;
       h_sf_Rmt_nb1_toponly->Draw("texte same") ;

       c_sf_Rmt_toponly->cd(2) ;
       h_sf_Rmt_nb2_toponly->Draw("colz") ;
       h_sf_Rmt_nb2_toponly->Draw("texte same") ;

       c_sf_Rmt_toponly->cd(3) ;
       h_sf_Rmt_nb3_toponly->Draw("colz") ;
       h_sf_Rmt_nb3_toponly->Draw("texte same") ;

       c_sf_Rmt_toponly->Update() ; c_sf_Rmt_toponly->Draw() ;


     //---------

       gStyle -> SetPaintTextFormat(".3f") ;

       TCanvas* c_Rmt_toponly = new TCanvas( "c_Rmt_toponly","c_Rmt_toponly", 400, 1000 ) ;
       c_Rmt_toponly->Clear() ;
     //c_Rmt_toponly->Divide(1,4) ;
       c_Rmt_toponly->Divide(1,3) ;


     //c_Rmt_toponly->cd(1) ;
     //h_Rmt_nb0_toponly->Draw("colz") ;
     //h_Rmt_nb0_toponly->Draw("texte same") ;

       c_Rmt_toponly->cd(1) ;
       h_Rmt_nb1_toponly->Draw("colz") ;
       h_Rmt_nb1_toponly->Draw("texte same") ;

       c_Rmt_toponly->cd(2) ;
       h_Rmt_nb2_toponly->Draw("colz") ;
       h_Rmt_nb2_toponly->Draw("texte same") ;

       c_Rmt_toponly->cd(3) ;
       h_Rmt_nb3_toponly->Draw("colz") ;
       h_Rmt_nb3_toponly->Draw("texte same") ;


  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       float ave_Rmt_nbgt0_ttsl(-1.) ;
       float ave_Rmt_nbgt0_ttdl(-1.) ;
       TH2F* h_sf_Rmt_nbgt0_ttsl   = makeSFplotLHB( "nbgt0"  , "ttsl", ave_Rmt_nbgt0_ttsl ) ; if ( h_sf_Rmt_nbgt0_ttsl   == 0x0 ) return ;
       TH2F* h_sf_Rmt_nbgt0_ttdl   = makeSFplotLHB( "nbgt0"  , "ttdl", ave_Rmt_nbgt0_ttdl ) ; if ( h_sf_Rmt_nbgt0_ttdl   == 0x0 ) return ;

       TH2F* h_Rmt_nbgt0_ttsl = (TH2F*) gDirectory->FindObject( "h_Rmt_nbgt0_ttsl" ) ; if ( h_Rmt_nbgt0_ttsl == 0x0 ) return ;
       TH2F* h_Rmt_nbgt0_ttdl = (TH2F*) gDirectory->FindObject( "h_Rmt_nbgt0_ttdl" ) ; if ( h_Rmt_nbgt0_ttdl == 0x0 ) return ;

       h_Rmt_nbgt0_ttsl -> SetMinimum(0.0) ;
       h_Rmt_nbgt0_ttdl -> SetMinimum(0.0) ;

       h_sf_Rmt_nbgt0_ttsl -> SetMaximum(2.5) ;
       h_sf_Rmt_nbgt0_ttdl -> SetMaximum(2.5) ;
       h_sf_Rmt_nbgt0_ttsl -> SetMinimum(0.0) ;
       h_sf_Rmt_nbgt0_ttdl -> SetMinimum(0.0) ;

       TCanvas* c_sl_dl_gt0b = new TCanvas("c_sl_dl_gt0b","c_sl_dl_gt0b", 1000, 800 ) ;

       c_sl_dl_gt0b -> Clear() ;
       c_sl_dl_gt0b -> Divide(2,2) ;

       c_sl_dl_gt0b->cd(1) ;
       h_Rmt_nbgt0_ttsl -> Draw("colz") ;
       h_Rmt_nbgt0_ttsl -> Draw("texte same") ;

       c_sl_dl_gt0b->cd(2) ;
       h_Rmt_nbgt0_ttdl -> Draw("colz") ;
       h_Rmt_nbgt0_ttdl -> Draw("texte same") ;

       c_sl_dl_gt0b->cd(3) ;
       h_sf_Rmt_nbgt0_ttsl -> Draw("colz") ;
       h_sf_Rmt_nbgt0_ttsl -> Draw("texte same") ;

       c_sl_dl_gt0b->cd(4) ;
       h_sf_Rmt_nbgt0_ttdl -> Draw("colz") ;
       h_sf_Rmt_nbgt0_ttdl -> Draw("texte same") ;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       float ave_Rmt_nb1_ttsl(-1.) ;
       float ave_Rmt_nb1_ttdl(-1.) ;
       TH2F* h_sf_Rmt_nb1_ttsl   = makeSFplotLHB( "nb1"  , "ttsl", ave_Rmt_nb1_ttsl ) ; if ( h_sf_Rmt_nb1_ttsl   == 0x0 ) return ;
       TH2F* h_sf_Rmt_nb1_ttdl   = makeSFplotLHB( "nb1"  , "ttdl", ave_Rmt_nb1_ttdl ) ; if ( h_sf_Rmt_nb1_ttdl   == 0x0 ) return ;

       TH2F* h_Rmt_nb1_ttsl = (TH2F*) gDirectory->FindObject( "h_Rmt_nb1_ttsl" ) ; if ( h_Rmt_nb1_ttsl == 0x0 ) return ;
       TH2F* h_Rmt_nb1_ttdl = (TH2F*) gDirectory->FindObject( "h_Rmt_nb1_ttdl" ) ; if ( h_Rmt_nb1_ttdl == 0x0 ) return ;

       h_Rmt_nb1_ttsl -> SetMinimum(0.0) ;
       h_Rmt_nb1_ttdl -> SetMinimum(0.0) ;

       h_sf_Rmt_nb1_ttsl -> SetMaximum(2.5) ;
       h_sf_Rmt_nb1_ttdl -> SetMaximum(2.5) ;
       h_sf_Rmt_nb1_ttsl -> SetMinimum(0.0) ;
       h_sf_Rmt_nb1_ttdl -> SetMinimum(0.0) ;

       TCanvas* c_sl_dl_1b = new TCanvas("c_sl_dl_1b","c_sl_dl_1b", 1000, 800 ) ;

       c_sl_dl_1b -> Clear() ;
       c_sl_dl_1b -> Divide(2,2) ;

       c_sl_dl_1b->cd(1) ;
       h_Rmt_nb1_ttsl -> Draw("colz") ;
       h_Rmt_nb1_ttsl -> Draw("texte same") ;

       c_sl_dl_1b->cd(2) ;
       h_Rmt_nb1_ttdl -> Draw("colz") ;
       h_Rmt_nb1_ttdl -> Draw("texte same") ;

       c_sl_dl_1b->cd(3) ;
       h_sf_Rmt_nb1_ttsl -> Draw("colz") ;
       h_sf_Rmt_nb1_ttsl -> Draw("texte same") ;

       c_sl_dl_1b->cd(4) ;
       h_sf_Rmt_nb1_ttdl -> Draw("colz") ;
       h_sf_Rmt_nb1_ttdl -> Draw("texte same") ;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       float ave_Rmt_nb2_ttsl(-1.) ;
       float ave_Rmt_nb2_ttdl(-1.) ;
       TH2F* h_sf_Rmt_nb2_ttsl   = makeSFplotLHB( "nb2"  , "ttsl", ave_Rmt_nb2_ttsl ) ; if ( h_sf_Rmt_nb2_ttsl   == 0x0 ) return ;
       TH2F* h_sf_Rmt_nb2_ttdl   = makeSFplotLHB( "nb2"  , "ttdl", ave_Rmt_nb2_ttdl ) ; if ( h_sf_Rmt_nb2_ttdl   == 0x0 ) return ;

       TH2F* h_Rmt_nb2_ttsl = (TH2F*) gDirectory->FindObject( "h_Rmt_nb2_ttsl" ) ; if ( h_Rmt_nb2_ttsl == 0x0 ) return ;
       TH2F* h_Rmt_nb2_ttdl = (TH2F*) gDirectory->FindObject( "h_Rmt_nb2_ttdl" ) ; if ( h_Rmt_nb2_ttdl == 0x0 ) return ;

       h_Rmt_nb2_ttsl -> SetMinimum(0.0) ;
       h_Rmt_nb2_ttdl -> SetMinimum(0.0) ;

       h_sf_Rmt_nb2_ttsl -> SetMaximum(2.5) ;
       h_sf_Rmt_nb2_ttdl -> SetMaximum(2.5) ;
       h_sf_Rmt_nb2_ttsl -> SetMinimum(0.0) ;
       h_sf_Rmt_nb2_ttdl -> SetMinimum(0.0) ;

       TCanvas* c_sl_dl_2b = new TCanvas("c_sl_dl_2b","c_sl_dl_2b", 1000, 800 ) ;

       c_sl_dl_2b -> Clear() ;
       c_sl_dl_2b -> Divide(2,2) ;

       c_sl_dl_2b->cd(1) ;
       h_Rmt_nb2_ttsl -> Draw("colz") ;
       h_Rmt_nb2_ttsl -> Draw("texte same") ;

       c_sl_dl_2b->cd(2) ;
       h_Rmt_nb2_ttdl -> Draw("colz") ;
       h_Rmt_nb2_ttdl -> Draw("texte same") ;

       c_sl_dl_2b->cd(3) ;
       h_sf_Rmt_nb2_ttsl -> Draw("colz") ;
       h_sf_Rmt_nb2_ttsl -> Draw("texte same") ;

       c_sl_dl_2b->cd(4) ;
       h_sf_Rmt_nb2_ttdl -> Draw("colz") ;
       h_sf_Rmt_nb2_ttdl -> Draw("texte same") ;

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       float ave_Rmt_nb3_ttsl(-1.) ;
       float ave_Rmt_nb3_ttdl(-1.) ;
       TH2F* h_sf_Rmt_nb3_ttsl   = makeSFplotLHB( "nb3"  , "ttsl", ave_Rmt_nb3_ttsl ) ; if ( h_sf_Rmt_nb3_ttsl   == 0x0 ) return ;
       TH2F* h_sf_Rmt_nb3_ttdl   = makeSFplotLHB( "nb3"  , "ttdl", ave_Rmt_nb3_ttdl ) ; if ( h_sf_Rmt_nb3_ttdl   == 0x0 ) return ;

       TH2F* h_Rmt_nb3_ttsl = (TH2F*) gDirectory->FindObject( "h_Rmt_nb3_ttsl" ) ; if ( h_Rmt_nb3_ttsl == 0x0 ) return ;
       TH2F* h_Rmt_nb3_ttdl = (TH2F*) gDirectory->FindObject( "h_Rmt_nb3_ttdl" ) ; if ( h_Rmt_nb3_ttdl == 0x0 ) return ;

       h_Rmt_nb3_ttsl -> SetMinimum(0.0) ;
       h_Rmt_nb3_ttdl -> SetMinimum(0.0) ;

       h_sf_Rmt_nb3_ttsl -> SetMaximum(2.5) ;
       h_sf_Rmt_nb3_ttdl -> SetMaximum(2.5) ;
       h_sf_Rmt_nb3_ttsl -> SetMinimum(0.0) ;
       h_sf_Rmt_nb3_ttdl -> SetMinimum(0.0) ;

       TCanvas* c_sl_dl_3b = new TCanvas("c_sl_dl_3b","c_sl_dl_3b", 1000, 800 ) ;

       c_sl_dl_3b -> Clear() ;
       c_sl_dl_3b -> Divide(2,2) ;

       c_sl_dl_3b->cd(1) ;
       h_Rmt_nb3_ttsl -> Draw("colz") ;
       h_Rmt_nb3_ttsl -> Draw("texte same") ;

       c_sl_dl_3b->cd(2) ;
       h_Rmt_nb3_ttdl -> Draw("colz") ;
       h_Rmt_nb3_ttdl -> Draw("texte same") ;

       c_sl_dl_3b->cd(3) ;
       h_sf_Rmt_nb3_ttsl -> Draw("colz") ;
       h_sf_Rmt_nb3_ttsl -> Draw("texte same") ;

       c_sl_dl_3b->cd(4) ;
       h_sf_Rmt_nb3_ttdl -> Draw("colz") ;
       h_sf_Rmt_nb3_ttdl -> Draw("texte same") ;

   } // dmc_drawplots_Rmc

   //=======================================================================================

   TH2F* makeSFplotLHB( const char* selname, const char* compstring, float& ave_Rmt ) {

       TH2F* h_mtl(0x0) ;
       TH2F* h_mth(0x0) ;

       char hname[1000] ;

       if ( strcmp( compstring, "comb" ) == 0 ) {
          sprintf( hname, "h_lhb_sl_mtl_%s", selname ) ;
          h_mtl = combineCompsLHB( hname ) ; if ( h_mtl == 0x0 ) return 0x0 ;
          sprintf( hname, "h_lhb_sl_mth_%s", selname ) ;
          h_mth = combineCompsLHB( hname ) ; if ( h_mth == 0x0 ) return 0x0 ;
       } else if ( strcmp( compstring, "toponly" ) == 0 ) {
          sprintf( hname, "h_lhb_sl_mtl_%s", selname ) ;
          h_mtl = combineToponlyLHB( hname ) ; if ( h_mtl == 0x0 ) return 0x0 ;
          sprintf( hname, "h_lhb_sl_mth_%s", selname ) ;
          h_mth = combineToponlyLHB( hname ) ; if ( h_mth == 0x0 ) return 0x0 ;
       } else {
          sprintf( hname, "h_lhb_sl_mtl_%s_%s_eb", selname, compstring ) ;
          h_mtl = (TH2F*) gDirectory->FindObject( hname ) ; if ( h_mtl == 0x0 ) return 0x0 ;
          sprintf( hname, "h_lhb_sl_mth_%s_%s_eb", selname, compstring ) ;
          h_mth = (TH2F*) gDirectory->FindObject( hname ) ; if ( h_mth == 0x0 ) return 0x0 ;
       }


       if ( ave_Rmt < 0. ) {
          double integral_mth = h_mth -> Integral() ;
          double integral_mtl = h_mtl -> Integral() ;
          ave_Rmt = integral_mth / integral_mtl ;
          printf("\n Ave R_MT for %s, %s is %3.5f\n\n", selname, compstring, ave_Rmt ) ;
       }

       sprintf( hname, "h_Rmt_%s_%s", selname, compstring ) ;
       TH2F* h_Rmt = (TH2F*) h_mth -> Clone( hname ) ;
       TString newtitle = h_Rmt->GetTitle() ;
       newtitle.ReplaceAll("MT>100","R_MT") ;
       h_Rmt->SetTitle( newtitle ) ;

       h_Rmt -> Divide( h_mtl ) ;

       sprintf( hname, "h_sf_Rmt_%s_%s", selname, compstring ) ;
       TH2F* h_sf_Rmt = (TH2F*) h_Rmt -> Clone( hname ) ;
       h_sf_Rmt -> Scale( 1./ave_Rmt ) ;
       newtitle.ReplaceAll("R_MT", "SF=R_MT/R_MT_ave") ;
       h_sf_Rmt -> SetTitle( newtitle ) ;

       h_Rmt -> UseCurrentStyle() ;
       h_sf_Rmt -> UseCurrentStyle() ;

       h_Rmt -> SetBinContent(4,1, -1.) ;
       h_sf_Rmt -> SetBinContent(4,1, -1.) ;

       return h_sf_Rmt ;

   } // makeSFplotLHB.



   //=======================================================================================

   TH2F* combineCompsLHB( const char* hname_base ) {

      TH2F* hret(0x0) ;

      for ( int ci=0; ci<nComps; ci++ ) {
         char hname[100] ;
         sprintf( hname, "%s_%s_eb", hname_base, compname[ci] ) ;
         TH2F* hcomp = (TH2F*) gDirectory->FindObject( hname ) ;
         if ( hcomp == 0x0 ) { printf("\n\n *** combineCompsLHB : can't find %s\n\n", hname ) ; cout << flush ; return 0x0 ; }
         if ( ci == 0 ) {
            sprintf( hname, "%s_comb_eb", hname_base ) ;
            hret = (TH2F*) hcomp -> Clone( hname ) ;
            TString combtitle = hret -> GetTitle() ;
            combtitle.ReplaceAll( compname[ci], "combined" ) ;
            hret->SetTitle( combtitle ) ;
         } else {
            hret -> Add( hcomp ) ;
         }
      } // ci.
      return hret ;
   }


   //=======================================================================================

   TH2F* combineToponlyLHB( const char* hname_base ) {

      TH2F* hret(0x0) ;

      for ( int ci=0; ci<(nComps-1); ci++ ) {
         char hname[100] ;
         sprintf( hname, "%s_%s_eb", hname_base, compname[ci] ) ;
         TH2F* hcomp = (TH2F*) gDirectory->FindObject( hname ) ;
         if ( hcomp == 0x0 ) { printf("\n\n *** combineToponlyLHB : can't find %s\n\n", hname ) ; cout << flush ; return 0x0 ; }
         if ( ci == 0 ) {
            sprintf( hname, "%s_toponly_eb", hname_base ) ;
            hret = (TH2F*) hcomp -> Clone( hname ) ;
            TString combtitle = hret -> GetTitle() ;
            combtitle.ReplaceAll( compname[ci], "combined" ) ;
            hret->SetTitle( combtitle ) ;
         } else {
            hret -> Add( hcomp ) ;
         }
      } // ci.
      return hret ;
   }


   //=======================================================================================
void loadHist(const char* filename, const char* pfx, const char* pat, Bool_t doAdd, Double_t scaleFactor)
{
  cout << " Reading histograms from file: " << filename << endl << flush ;
  TFile inf(filename) ;
  //inf.ReadAll() ;
  TList* list = inf.GetListOfKeys() ;
  TIterator* iter = list->MakeIterator();

  TRegexp re(pat,kTRUE) ;
  std::cout << "pat = " << pat << std::endl ;

  gDirectory->cd("Rint:") ;

  TObject* obj ;
  TKey* key ;
  std::cout << "doAdd = " << (doAdd?"T":"F") << std::endl ;
  std::cout << "loadHist: reading." ;
  while((key=(TKey*)iter->Next())) {
   
    Int_t ridx = TString(key->GetName()).Index(re) ;    
    if (ridx==-1) {
      continue ;
    }

    obj = inf.Get(key->GetName()) ;
    TObject* clone ;
    if (pfx) {

      // Find existing TH1-derived objects
      TObject* oldObj = 0 ;
      if (doAdd){
	oldObj = gDirectory->Get(Form("%s_%s",pfx,obj->GetName())) ;
	if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	  oldObj = 0 ;
	}
      }
      if (oldObj) {
	clone = oldObj ;
        if ( scaleFactor > 0 ) {
           ((TH1*)clone)->Sumw2() ;
           ((TH1*)clone)->Add((TH1*)obj, scaleFactor) ;
        } else {
           ((TH1*)clone)->Add((TH1*)obj) ;
        }
      } else {
	clone = obj->Clone(Form("%s_%s",pfx,obj->GetName())) ;
      }


    } else {

      // Find existing TH1-derived objects
      TObject* oldObj = 0 ;
      if (doAdd){
	oldObj = gDirectory->Get(key->GetName()) ;
	if (oldObj && !oldObj->IsA()->InheritsFrom(TH1::Class())) {
	  oldObj = 0 ;
	}
      }

      if (oldObj) {
	clone = oldObj ;
        if ( scaleFactor > 0 ) {
           ((TH1*)clone)->Sumw2() ;
           ((TH1*)clone)->Add((TH1*)obj, scaleFactor) ;
        } else {
           ((TH1*)clone)->Add((TH1*)obj) ;
        }
      } else {
	clone = obj->Clone() ;
      }
    }
    if ( scaleFactor > 0 && !doAdd ) {
       ((TH1*) clone)->Sumw2() ;
       ((TH1*) clone)->Scale(scaleFactor) ;
    }
    if (!gDirectory->GetList()->FindObject(clone)) {
      gDirectory->Append(clone) ;
    }
    std::cout << "." ;
    std::cout.flush() ;
  }
  std::cout << std::endl;
  inf.Close() ;
  delete iter ;
}

