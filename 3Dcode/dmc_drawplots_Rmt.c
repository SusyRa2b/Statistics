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
#include "TLatex.h"

#include <iostream>

    double kfactor[7] ;

  //-- met4-ht4-v15
      const int nBinsBtag  = 3 ;
      const int nBinsMET   = 4 ;
      const int nBinsHT    = 4 ;
      float Mbins[nBinsMET+1] = {125.,150.,250.,350.,99999.};
      float Hbins[nBinsHT+1] = {400.,500.,800.,1000.,99999.};

    bool recycleCanvas(false) ;
    bool savePdf ;

    bool kfactorAlreadyApplied(false) ;

    char inrootfile[10000] ;

    int nComps(4) ;
    char compname[4][100] ;

    TH2F* combineCompsLHB( const char* hname_base ) ;
    TH2F* combineToponlyLHB( const char* hname_base ) ;
    TH1F* flattenLHB( TH2F* h2 ) ;

   TH2F* makeSFplotLHB( const char* selname, const char* compstring, float& ave_Rmt ) ;

   void drawSet( const char* hname_base, const char* xtitle ) ;


   void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;


   //------------

   void dmc_drawplots_Rmc( const char* infile = "rootfiles/dmc_plots_Rmt_all.root",
                           bool arg_savePdf = false,
                           bool arg_recycleCanvas = false ) {

       gDirectory->Delete("h*") ;

       loadHist( infile ) ;

       sprintf( inrootfile, "%s", infile ) ;

       savePdf = arg_savePdf ;
       recycleCanvas = arg_recycleCanvas ;

       gStyle -> SetOptStat(0) ;
       gStyle -> SetTitleH(0.06 ) ;
       gStyle -> SetTitleW(0.85 ) ;
       gStyle -> SetMarkerSize(2.5) ;
       gStyle -> SetPadRightMargin(0.15) ;
       gStyle -> SetPaintTextFormat(".2f") ;
       gStyle -> SetLabelSize(0.07,"x") ;
       gStyle -> SetLabelSize(0.07,"y") ;
       gStyle -> SetPadBottomMargin(0.12) ;

       gStyle -> SetOptStat(0) ;

       sprintf( compname[0], "ttdl" ) ;        kfactor[1] = 1.00 ;
       sprintf( compname[1], "singlet" ) ;     kfactor[2] = 1.00 ;
       sprintf( compname[2], "wjets" ) ;       kfactor[3] = 1.00 ;
       sprintf( compname[3], "ttsl" ) ;        kfactor[0] = 1.00 ;

  //   for ( int bbi=0; bbi<nBinsBtag; bbi++ ) {
  //      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
  //         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
  //            char hname[1000] ;
  //            sprintf( hname, "h_mt_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
  //            drawSet( hname, "MT" ) ;
  //         } // hbi.
  //      } // mbi.
  //   } // bbi.


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

       TCanvas* c_sf_Rmt = (TCanvas*) gDirectory->FindObject("c_sf_Rmt") ;
       if ( c_sf_Rmt == 0x0 ) {
          c_sf_Rmt = new TCanvas( "c_sf_Rmt","c_sf_Rmt", 1000, 300 ) ;
       }
       c_sf_Rmt->Clear() ;
       c_sf_Rmt->Divide(3,1) ;

       h_sf_Rmt_nb1->UseCurrentStyle() ;
       h_sf_Rmt_nb2->UseCurrentStyle() ;
       h_sf_Rmt_nb3->UseCurrentStyle() ;

       h_sf_Rmt_nb1 -> SetTitle("SF = R_MT / (ave over MET,HT,nB R_MT),   nB=1") ;
       h_sf_Rmt_nb2 -> SetTitle("SF = R_MT / (ave over MET,HT,nB R_MT),   nB=2") ;
       h_sf_Rmt_nb3 -> SetTitle("SF = R_MT / (ave over MET,HT,nB R_MT),   nB=3") ;

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

       c_sf_Rmt->Update() ; c_sf_Rmt->Draw() ;
       if ( savePdf ) {
          c_sf_Rmt->SaveAs("outputfiles/rmt_sf_ht_vs_met.pdf") ;
       }


     //---------

       gStyle -> SetPaintTextFormat(".3f") ;

       TCanvas* c_Rmt = (TCanvas*) gDirectory->FindObject("c_Rmt") ;
       if ( c_Rmt == 0x0 ) {
          c_Rmt = new TCanvas( "c_Rmt","c_Rmt", 1000, 300 ) ;
       }
       c_Rmt->Clear() ;
       c_Rmt->Divide(3,1) ;

       h_Rmt_nb1->UseCurrentStyle() ;
       h_Rmt_nb2->UseCurrentStyle() ;
       h_Rmt_nb3->UseCurrentStyle() ;

       h_Rmt_nb1 -> SetTitle("R_MT = N_SL(MT>100)/N_SL(MT<100),   nB=1") ;
       h_Rmt_nb2 -> SetTitle("R_MT = N_SL(MT>100)/N_SL(MT<100),   nB=2") ;
       h_Rmt_nb3 -> SetTitle("R_MT = N_SL(MT>100)/N_SL(MT<100),   nB=3") ;



       c_Rmt->cd(1) ;
       h_Rmt_nb1->Draw("colz") ;
       h_Rmt_nb1->Draw("texte same") ;

       c_Rmt->cd(2) ;
       h_Rmt_nb2->Draw("colz") ;
       h_Rmt_nb2->Draw("texte same") ;

       c_Rmt->cd(3) ;
       h_Rmt_nb3->Draw("colz") ;
       h_Rmt_nb3->Draw("texte same") ;

       c_Rmt->Update() ; c_Rmt->Draw() ;
       if ( savePdf ) {
          c_Rmt->SaveAs("outputfiles/rmt_ht_vs_met.pdf") ;
       }

     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       TH2F* h_sf2_Rmt_nb1 = (TH2F*) h_Rmt_nb1 -> Clone( "h_sf2_Rmt_nb1" ) ;
       TH2F* h_sf2_Rmt_nb2 = (TH2F*) h_Rmt_nb2 -> Clone( "h_sf2_Rmt_nb2" ) ;
       TH2F* h_sf2_Rmt_nb3 = (TH2F*) h_Rmt_nb3 -> Clone( "h_sf2_Rmt_nb3" ) ;

       h_sf2_Rmt_nb1 -> Divide( h_Rmt_nbgt0 ) ;
       h_sf2_Rmt_nb2 -> Divide( h_Rmt_nbgt0 ) ;
       h_sf2_Rmt_nb3 -> Divide( h_Rmt_nbgt0 ) ;

       gStyle -> SetPaintTextFormat(".2f") ;

       TCanvas* c_sf2_Rmt = (TCanvas*) gDirectory->FindObject("c_sf2_Rmt") ;
       if ( c_sf2_Rmt == 0x0 ) {
          c_sf2_Rmt = new TCanvas( "c_sf2_Rmt","c_sf2_Rmt", 1000, 300 ) ;
       }
       c_sf2_Rmt->Clear() ;
       c_sf2_Rmt->Divide(3,1) ;

       h_sf2_Rmt_nb1->UseCurrentStyle() ;
       h_sf2_Rmt_nb2->UseCurrentStyle() ;
       h_sf2_Rmt_nb3->UseCurrentStyle() ;

       h_sf2_Rmt_nb1 -> SetTitle("SF' = R_MT / (ave over nB R_MT),   nB=1") ;
       h_sf2_Rmt_nb2 -> SetTitle("SF' = R_MT / (ave over nB R_MT),   nB=2") ;
       h_sf2_Rmt_nb3 -> SetTitle("SF' = R_MT / (ave over nB R_MT),   nB=3") ;

       h_sf2_Rmt_nb1 -> SetMinimum( 0.30 ) ;
       h_sf2_Rmt_nb2 -> SetMinimum( 0.30 ) ;
       h_sf2_Rmt_nb3 -> SetMinimum( 0.30 ) ;

       h_sf2_Rmt_nb1 -> SetMaximum( 2.5 ) ;
       h_sf2_Rmt_nb2 -> SetMaximum( 2.5 ) ;
       h_sf2_Rmt_nb3 -> SetMaximum( 2.5 ) ;

       c_sf2_Rmt->cd(1) ;
       h_sf2_Rmt_nb1->Draw("colz") ;
       h_sf2_Rmt_nb1->Draw("texte same") ;

       c_sf2_Rmt->cd(2) ;
       h_sf2_Rmt_nb2->Draw("colz") ;
       h_sf2_Rmt_nb2->Draw("texte same") ;

       c_sf2_Rmt->cd(3) ;
       h_sf2_Rmt_nb3->Draw("colz") ;
       h_sf2_Rmt_nb3->Draw("texte same") ;

       c_sf2_Rmt->Update() ; c_sf2_Rmt->Draw() ;

       c_sf2_Rmt->Update() ; c_sf2_Rmt->Draw() ;
       if ( savePdf ) {
          c_sf2_Rmt->SaveAs("outputfiles/rmt_sf2_ht_vs_met.pdf") ;
       }



     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       gStyle->SetMarkerSize(1) ;

       TH1F* h_sf_Rmt_nb1_flat = flattenLHB( h_sf_Rmt_nb1 ) ;
       TH1F* h_sf_Rmt_nb2_flat = flattenLHB( h_sf_Rmt_nb2 ) ;
       TH1F* h_sf_Rmt_nb3_flat = flattenLHB( h_sf_Rmt_nb3 ) ;

       h_sf_Rmt_nb1_flat->SetMinimum(0.) ;
       h_sf_Rmt_nb2_flat->SetMinimum(0.) ;
       h_sf_Rmt_nb3_flat->SetMinimum(0.) ;

       h_sf_Rmt_nb1_flat->SetMaximum(3.0) ;
       h_sf_Rmt_nb2_flat->SetMaximum(3.0) ;
       h_sf_Rmt_nb3_flat->SetMaximum(3.0) ;

       h_sf_Rmt_nb1_flat->SetMarkerStyle(20) ;
       h_sf_Rmt_nb2_flat->SetMarkerStyle(20) ;
       h_sf_Rmt_nb3_flat->SetMarkerStyle(20) ;



       TCanvas* c_sf_Rmt_flat = (TCanvas*) gDirectory->FindObject("c_sf_Rmt_flat") ;
       if ( c_sf_Rmt_flat == 0x0 ) {
          c_sf_Rmt_flat = new TCanvas( "c_sf_Rmt_flat","c_sf_Rmt_flat", 1000, 300 ) ;
       }
       c_sf_Rmt_flat->Clear() ;
       c_sf_Rmt_flat->Divide(3,1) ;

       c_sf_Rmt_flat->cd(1) ;
       h_sf_Rmt_nb1_flat->Draw() ;
       gPad->SetGridy(1) ;

       c_sf_Rmt_flat->cd(2) ;
       h_sf_Rmt_nb2_flat->Draw() ;
       gPad->SetGridy(1) ;

       c_sf_Rmt_flat->cd(3) ;
       h_sf_Rmt_nb3_flat->Draw() ;
       gPad->SetGridy(1) ;

       c_sf_Rmt_flat->Update() ; c_sf_Rmt_flat->Draw() ;
       if ( savePdf ) {
          c_sf_Rmt_flat->SaveAs("outputfiles/rmt_sf_flat.pdf") ;
       }


     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       gStyle->SetMarkerSize(1) ;

       TH1F* h_sf2_Rmt_nb1_flat = flattenLHB( h_sf2_Rmt_nb1 ) ;
       TH1F* h_sf2_Rmt_nb2_flat = flattenLHB( h_sf2_Rmt_nb2 ) ;
       TH1F* h_sf2_Rmt_nb3_flat = flattenLHB( h_sf2_Rmt_nb3 ) ;

       h_sf2_Rmt_nb1_flat->SetMinimum(0.) ;
       h_sf2_Rmt_nb2_flat->SetMinimum(0.) ;
       h_sf2_Rmt_nb3_flat->SetMinimum(0.) ;

       h_sf2_Rmt_nb1_flat->SetMaximum(3.0) ;
       h_sf2_Rmt_nb2_flat->SetMaximum(3.0) ;
       h_sf2_Rmt_nb3_flat->SetMaximum(3.0) ;

       h_sf2_Rmt_nb1_flat->SetMarkerStyle(20) ;
       h_sf2_Rmt_nb2_flat->SetMarkerStyle(20) ;
       h_sf2_Rmt_nb3_flat->SetMarkerStyle(20) ;



       TCanvas* c_sf2_Rmt_flat = (TCanvas*) gDirectory->FindObject("c_sf2_Rmt_flat") ;
       if ( c_sf2_Rmt_flat == 0x0 ) {
          c_sf2_Rmt_flat = new TCanvas( "c_sf2_Rmt_flat","c_sf2_Rmt_flat", 1000, 300 ) ;
       }
       c_sf2_Rmt_flat->Clear() ;
       c_sf2_Rmt_flat->Divide(3,1) ;

       c_sf2_Rmt_flat->cd(1) ;
       h_sf2_Rmt_nb1_flat->Draw() ;
       gPad->SetGridy(1) ;

       c_sf2_Rmt_flat->cd(2) ;
       h_sf2_Rmt_nb2_flat->Draw() ;
       gPad->SetGridy(1) ;

       c_sf2_Rmt_flat->cd(3) ;
       h_sf2_Rmt_nb3_flat->Draw() ;
       gPad->SetGridy(1) ;

       c_sf2_Rmt_flat->Update() ; c_sf2_Rmt_flat->Draw() ;
       if ( savePdf ) {
          c_sf2_Rmt_flat->SaveAs("outputfiles/rmt_sf2_flat.pdf") ;
       }





     //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //  float ave_Rmt_toponly(-1.) ;

   //  TH2F* h_sf_Rmt_nbgt0_toponly = makeSFplotLHB( "nbgt0", "toponly", ave_Rmt_toponly ) ; if ( h_sf_Rmt_nbgt0_toponly == 0x0 ) return ;
   //  TH2F* h_sf_Rmt_nb0_toponly   = makeSFplotLHB( "nb0"  , "toponly", ave_Rmt_toponly ) ; if ( h_sf_Rmt_nb0_toponly   == 0x0 ) return ;
   //  TH2F* h_sf_Rmt_nb1_toponly   = makeSFplotLHB( "nb1"  , "toponly", ave_Rmt_toponly ) ; if ( h_sf_Rmt_nb1_toponly   == 0x0 ) return ;
   //  TH2F* h_sf_Rmt_nb2_toponly   = makeSFplotLHB( "nb2"  , "toponly", ave_Rmt_toponly ) ; if ( h_sf_Rmt_nb2_toponly   == 0x0 ) return ;
   //  TH2F* h_sf_Rmt_nb3_toponly   = makeSFplotLHB( "nb3"  , "toponly", ave_Rmt_toponly ) ; if ( h_sf_Rmt_nb3_toponly   == 0x0 ) return ;

   //  TH2F* h_Rmt_nbgt0_toponly = (TH2F*) gDirectory->FindObject( "h_Rmt_nbgt0_toponly"  ) ;
   //  TH2F* h_Rmt_nb0_toponly   = (TH2F*) gDirectory->FindObject( "h_Rmt_nb0_toponly"    ) ;
   //  TH2F* h_Rmt_nb1_toponly   = (TH2F*) gDirectory->FindObject( "h_Rmt_nb1_toponly"    ) ;
   //  TH2F* h_Rmt_nb2_toponly   = (TH2F*) gDirectory->FindObject( "h_Rmt_nb2_toponly"    ) ;
   //  TH2F* h_Rmt_nb3_toponly   = (TH2F*) gDirectory->FindObject( "h_Rmt_nb3_toponly"    ) ;

   //  h_sf_Rmt_nbgt0_toponly -> SetMaximum(2.5) ;
   //  h_sf_Rmt_nb0_toponly   -> SetMaximum(2.5) ;
   //  h_sf_Rmt_nb1_toponly   -> SetMaximum(2.5) ;
   //  h_sf_Rmt_nb2_toponly   -> SetMaximum(2.5) ;
   //  h_sf_Rmt_nb3_toponly   -> SetMaximum(2.5) ;

   //  h_sf_Rmt_nbgt0_toponly -> SetMinimum(0.3) ;
   //  h_sf_Rmt_nb0_toponly   -> SetMinimum(0.3) ;
   //  h_sf_Rmt_nb1_toponly   -> SetMinimum(0.3) ;
   //  h_sf_Rmt_nb2_toponly   -> SetMinimum(0.3) ;
   //  h_sf_Rmt_nb3_toponly   -> SetMinimum(0.3) ;

   //  h_Rmt_nbgt0_toponly -> SetMaximum(0.3) ;
   //  h_Rmt_nb0_toponly   -> SetMaximum(0.3) ;
   //  h_Rmt_nb1_toponly   -> SetMaximum(0.3) ;
   //  h_Rmt_nb2_toponly   -> SetMaximum(0.3) ;
   //  h_Rmt_nb3_toponly   -> SetMaximum(0.3) ;

   //  h_Rmt_nbgt0_toponly -> SetMinimum(0.00) ;
   //  h_Rmt_nb0_toponly   -> SetMinimum(0.00) ;
   //  h_Rmt_nb1_toponly   -> SetMinimum(0.00) ;
   //  h_Rmt_nb2_toponly   -> SetMinimum(0.00) ;
   //  h_Rmt_nb3_toponly   -> SetMinimum(0.00) ;

   ////---------

   //  gStyle -> SetPaintTextFormat(".2f") ;

   //  TCanvas* c_sf_Rmt_toponly = new TCanvas( "c_sf_Rmt_toponly","c_sf_Rmt_toponly", 400, 1000 ) ;
   //  c_sf_Rmt_toponly->Clear() ;
   //  c_sf_Rmt_toponly->Divide(1,3) ;


   //  c_sf_Rmt_toponly->cd(1) ;
   //  h_sf_Rmt_nb1_toponly->Draw("colz") ;
   //  h_sf_Rmt_nb1_toponly->Draw("texte same") ;

   //  c_sf_Rmt_toponly->cd(2) ;
   //  h_sf_Rmt_nb2_toponly->Draw("colz") ;
   //  h_sf_Rmt_nb2_toponly->Draw("texte same") ;

   //  c_sf_Rmt_toponly->cd(3) ;
   //  h_sf_Rmt_nb3_toponly->Draw("colz") ;
   //  h_sf_Rmt_nb3_toponly->Draw("texte same") ;

   //  c_sf_Rmt_toponly->Update() ; c_sf_Rmt_toponly->Draw() ;


   ////---------

   //  gStyle -> SetPaintTextFormat(".3f") ;

   //  TCanvas* c_Rmt_toponly = new TCanvas( "c_Rmt_toponly","c_Rmt_toponly", 400, 1000 ) ;
   //  c_Rmt_toponly->Clear() ;
   //  c_Rmt_toponly->Divide(1,3) ;


   //  c_Rmt_toponly->cd(1) ;
   //  h_Rmt_nb1_toponly->Draw("colz") ;
   //  h_Rmt_nb1_toponly->Draw("texte same") ;

   //  c_Rmt_toponly->cd(2) ;
   //  h_Rmt_nb2_toponly->Draw("colz") ;
   //  h_Rmt_nb2_toponly->Draw("texte same") ;

   //  c_Rmt_toponly->cd(3) ;
   //  h_Rmt_nb3_toponly->Draw("colz") ;
   //  h_Rmt_nb3_toponly->Draw("texte same") ;


   ////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   //  float ave_Rmt_nbgt0_ttsl(-1.) ;
   //  float ave_Rmt_nbgt0_ttdl(-1.) ;
   //  TH2F* h_sf_Rmt_nbgt0_ttsl   = makeSFplotLHB( "nbgt0"  , "ttsl", ave_Rmt_nbgt0_ttsl ) ; if ( h_sf_Rmt_nbgt0_ttsl   == 0x0 ) return ;
   //  TH2F* h_sf_Rmt_nbgt0_ttdl   = makeSFplotLHB( "nbgt0"  , "ttdl", ave_Rmt_nbgt0_ttdl ) ; if ( h_sf_Rmt_nbgt0_ttdl   == 0x0 ) return ;

   //  TH2F* h_Rmt_nbgt0_ttsl = (TH2F*) gDirectory->FindObject( "h_Rmt_nbgt0_ttsl" ) ; if ( h_Rmt_nbgt0_ttsl == 0x0 ) return ;
   //  TH2F* h_Rmt_nbgt0_ttdl = (TH2F*) gDirectory->FindObject( "h_Rmt_nbgt0_ttdl" ) ; if ( h_Rmt_nbgt0_ttdl == 0x0 ) return ;

   //  h_Rmt_nbgt0_ttsl -> SetMinimum(0.0) ;
   //  h_Rmt_nbgt0_ttdl -> SetMinimum(0.0) ;

   //  h_sf_Rmt_nbgt0_ttsl -> SetMaximum(2.5) ;
   //  h_sf_Rmt_nbgt0_ttdl -> SetMaximum(2.5) ;
   //  h_sf_Rmt_nbgt0_ttsl -> SetMinimum(0.0) ;
   //  h_sf_Rmt_nbgt0_ttdl -> SetMinimum(0.0) ;

   //  TCanvas* c_sl_dl_gt0b = new TCanvas("c_sl_dl_gt0b","c_sl_dl_gt0b", 1000, 800 ) ;

   //  c_sl_dl_gt0b -> Clear() ;
   //  c_sl_dl_gt0b -> Divide(2,2) ;

   //  c_sl_dl_gt0b->cd(1) ;
   //  h_Rmt_nbgt0_ttsl -> Draw("colz") ;
   //  h_Rmt_nbgt0_ttsl -> Draw("texte same") ;

   //  c_sl_dl_gt0b->cd(2) ;
   //  h_Rmt_nbgt0_ttdl -> Draw("colz") ;
   //  h_Rmt_nbgt0_ttdl -> Draw("texte same") ;

   //  c_sl_dl_gt0b->cd(3) ;
   //  h_sf_Rmt_nbgt0_ttsl -> Draw("colz") ;
   //  h_sf_Rmt_nbgt0_ttsl -> Draw("texte same") ;

   //  c_sl_dl_gt0b->cd(4) ;
   //  h_sf_Rmt_nbgt0_ttdl -> Draw("colz") ;
   //  h_sf_Rmt_nbgt0_ttdl -> Draw("texte same") ;

   ////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   //  float ave_Rmt_nb1_ttsl(-1.) ;
   //  float ave_Rmt_nb1_ttdl(-1.) ;
   //  TH2F* h_sf_Rmt_nb1_ttsl   = makeSFplotLHB( "nb1"  , "ttsl", ave_Rmt_nb1_ttsl ) ; if ( h_sf_Rmt_nb1_ttsl   == 0x0 ) return ;
   //  TH2F* h_sf_Rmt_nb1_ttdl   = makeSFplotLHB( "nb1"  , "ttdl", ave_Rmt_nb1_ttdl ) ; if ( h_sf_Rmt_nb1_ttdl   == 0x0 ) return ;

   //  TH2F* h_Rmt_nb1_ttsl = (TH2F*) gDirectory->FindObject( "h_Rmt_nb1_ttsl" ) ; if ( h_Rmt_nb1_ttsl == 0x0 ) return ;
   //  TH2F* h_Rmt_nb1_ttdl = (TH2F*) gDirectory->FindObject( "h_Rmt_nb1_ttdl" ) ; if ( h_Rmt_nb1_ttdl == 0x0 ) return ;

   //  h_Rmt_nb1_ttsl -> SetMinimum(0.0) ;
   //  h_Rmt_nb1_ttdl -> SetMinimum(0.0) ;

   //  h_sf_Rmt_nb1_ttsl -> SetMaximum(2.5) ;
   //  h_sf_Rmt_nb1_ttdl -> SetMaximum(2.5) ;
   //  h_sf_Rmt_nb1_ttsl -> SetMinimum(0.0) ;
   //  h_sf_Rmt_nb1_ttdl -> SetMinimum(0.0) ;

   //  TCanvas* c_sl_dl_1b = new TCanvas("c_sl_dl_1b","c_sl_dl_1b", 1000, 800 ) ;

   //  c_sl_dl_1b -> Clear() ;
   //  c_sl_dl_1b -> Divide(2,2) ;

   //  c_sl_dl_1b->cd(1) ;
   //  h_Rmt_nb1_ttsl -> Draw("colz") ;
   //  h_Rmt_nb1_ttsl -> Draw("texte same") ;

   //  c_sl_dl_1b->cd(2) ;
   //  h_Rmt_nb1_ttdl -> Draw("colz") ;
   //  h_Rmt_nb1_ttdl -> Draw("texte same") ;

   //  c_sl_dl_1b->cd(3) ;
   //  h_sf_Rmt_nb1_ttsl -> Draw("colz") ;
   //  h_sf_Rmt_nb1_ttsl -> Draw("texte same") ;

   //  c_sl_dl_1b->cd(4) ;
   //  h_sf_Rmt_nb1_ttdl -> Draw("colz") ;
   //  h_sf_Rmt_nb1_ttdl -> Draw("texte same") ;

   ////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   //  float ave_Rmt_nb2_ttsl(-1.) ;
   //  float ave_Rmt_nb2_ttdl(-1.) ;
   //  TH2F* h_sf_Rmt_nb2_ttsl   = makeSFplotLHB( "nb2"  , "ttsl", ave_Rmt_nb2_ttsl ) ; if ( h_sf_Rmt_nb2_ttsl   == 0x0 ) return ;
   //  TH2F* h_sf_Rmt_nb2_ttdl   = makeSFplotLHB( "nb2"  , "ttdl", ave_Rmt_nb2_ttdl ) ; if ( h_sf_Rmt_nb2_ttdl   == 0x0 ) return ;

   //  TH2F* h_Rmt_nb2_ttsl = (TH2F*) gDirectory->FindObject( "h_Rmt_nb2_ttsl" ) ; if ( h_Rmt_nb2_ttsl == 0x0 ) return ;
   //  TH2F* h_Rmt_nb2_ttdl = (TH2F*) gDirectory->FindObject( "h_Rmt_nb2_ttdl" ) ; if ( h_Rmt_nb2_ttdl == 0x0 ) return ;

   //  h_Rmt_nb2_ttsl -> SetMinimum(0.0) ;
   //  h_Rmt_nb2_ttdl -> SetMinimum(0.0) ;

   //  h_sf_Rmt_nb2_ttsl -> SetMaximum(2.5) ;
   //  h_sf_Rmt_nb2_ttdl -> SetMaximum(2.5) ;
   //  h_sf_Rmt_nb2_ttsl -> SetMinimum(0.0) ;
   //  h_sf_Rmt_nb2_ttdl -> SetMinimum(0.0) ;

   //  TCanvas* c_sl_dl_2b = new TCanvas("c_sl_dl_2b","c_sl_dl_2b", 1000, 800 ) ;

   //  c_sl_dl_2b -> Clear() ;
   //  c_sl_dl_2b -> Divide(2,2) ;

   //  c_sl_dl_2b->cd(1) ;
   //  h_Rmt_nb2_ttsl -> Draw("colz") ;
   //  h_Rmt_nb2_ttsl -> Draw("texte same") ;

   //  c_sl_dl_2b->cd(2) ;
   //  h_Rmt_nb2_ttdl -> Draw("colz") ;
   //  h_Rmt_nb2_ttdl -> Draw("texte same") ;

   //  c_sl_dl_2b->cd(3) ;
   //  h_sf_Rmt_nb2_ttsl -> Draw("colz") ;
   //  h_sf_Rmt_nb2_ttsl -> Draw("texte same") ;

   //  c_sl_dl_2b->cd(4) ;
   //  h_sf_Rmt_nb2_ttdl -> Draw("colz") ;
   //  h_sf_Rmt_nb2_ttdl -> Draw("texte same") ;

   ////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   //  float ave_Rmt_nb3_ttsl(-1.) ;
   //  float ave_Rmt_nb3_ttdl(-1.) ;
   //  TH2F* h_sf_Rmt_nb3_ttsl   = makeSFplotLHB( "nb3"  , "ttsl", ave_Rmt_nb3_ttsl ) ; if ( h_sf_Rmt_nb3_ttsl   == 0x0 ) return ;
   //  TH2F* h_sf_Rmt_nb3_ttdl   = makeSFplotLHB( "nb3"  , "ttdl", ave_Rmt_nb3_ttdl ) ; if ( h_sf_Rmt_nb3_ttdl   == 0x0 ) return ;

   //  TH2F* h_Rmt_nb3_ttsl = (TH2F*) gDirectory->FindObject( "h_Rmt_nb3_ttsl" ) ; if ( h_Rmt_nb3_ttsl == 0x0 ) return ;
   //  TH2F* h_Rmt_nb3_ttdl = (TH2F*) gDirectory->FindObject( "h_Rmt_nb3_ttdl" ) ; if ( h_Rmt_nb3_ttdl == 0x0 ) return ;

   //  h_Rmt_nb3_ttsl -> SetMinimum(0.0) ;
   //  h_Rmt_nb3_ttdl -> SetMinimum(0.0) ;

   //  h_sf_Rmt_nb3_ttsl -> SetMaximum(2.5) ;
   //  h_sf_Rmt_nb3_ttdl -> SetMaximum(2.5) ;
   //  h_sf_Rmt_nb3_ttsl -> SetMinimum(0.0) ;
   //  h_sf_Rmt_nb3_ttdl -> SetMinimum(0.0) ;

   //  TCanvas* c_sl_dl_3b = new TCanvas("c_sl_dl_3b","c_sl_dl_3b", 1000, 800 ) ;

   //  c_sl_dl_3b -> Clear() ;
   //  c_sl_dl_3b -> Divide(2,2) ;

   //  c_sl_dl_3b->cd(1) ;
   //  h_Rmt_nb3_ttsl -> Draw("colz") ;
   //  h_Rmt_nb3_ttsl -> Draw("texte same") ;

   //  c_sl_dl_3b->cd(2) ;
   //  h_Rmt_nb3_ttdl -> Draw("colz") ;
   //  h_Rmt_nb3_ttdl -> Draw("texte same") ;

   //  c_sl_dl_3b->cd(3) ;
   //  h_sf_Rmt_nb3_ttsl -> Draw("colz") ;
   //  h_sf_Rmt_nb3_ttsl -> Draw("texte same") ;

   //  c_sl_dl_3b->cd(4) ;
   //  h_sf_Rmt_nb3_ttdl -> Draw("colz") ;
   //  h_sf_Rmt_nb3_ttdl -> Draw("texte same") ;

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

      for ( int ci=0; ci<nComps; ci++ ) {
         if ( strcmp( compname[ci], "wjets" ) == 0 ) continue ;
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

   //=======================================================================================

   void drawSet( const char* hname_base, const char* xtitle ) {

      bool isMTPlot(false) ;
      TString hnstr( hname_base ) ;

      if ( hnstr.Contains("h_mt") ) {
         isMTPlot = true ;
      }

      printf(" drawSet : %s\n", hname_base ) ;

      bool islogy = gStyle->GetOptLogy() ;

      char cname[1000] ;
      if ( islogy ) {
         if ( recycleCanvas ) {
            sprintf( cname, "can_logy" ) ;
         } else {
            sprintf( cname, "can_logy_%s", hname_base ) ;
         }
      } else {
         if ( recycleCanvas ) {
            sprintf( cname, "can" ) ;
         } else {
            sprintf( cname, "can_%s", hname_base ) ;
         }
      }
      TCanvas* dmccan = (TCanvas*) gDirectory->FindObject( cname ) ;
      if ( dmccan == 0x0 ) {
         dmccan = new TCanvas( cname, hname_base, 600, 450 ) ;
      }
      dmccan->Clear() ;

      char hname[1000] ;

      sprintf( hname, "%s_mcstack", hname_base ) ;
      THStack* hmcstack = new THStack() ;

      sprintf( hname, "%s_mcsum", hname_base ) ;
      TH1F* hmcsum(0x0) ;

      TLegend* legend = new TLegend( 0.80, 0.67, 0.95, 0.92 ) ;
      legend->SetFillColor(kWhite) ;

      for ( int ci=0; ci<nComps; ci++ ) {


         sprintf( hname, "%s_%s", hname_base, compname[ci] ) ;
         TH1F* hmc = (TH1F*) gDirectory->FindObject( hname ) ;
         if ( ci==0 ) {
            sprintf( hname, "%s_mcsum", hname_base ) ;
            hmcsum = (TH1F*) hmc->Clone( hname ) ;
            hmcsum -> Reset() ;
         }
         if ( !kfactorAlreadyApplied ) { hmc->Scale( kfactor[ci] ) ; }
         if ( hmc == 0x0 ) { printf("\n\n *** drawSet: missing MC hist %s\n", hname ) ; return ; }
         hmcsum -> Add( hmc ) ;
         hmcstack -> Add( hmc ) ;

      }

      for ( int ci=nComps-1; ci>=0; ci-- ) {
         sprintf( hname, "%s_%s", hname_base, compname[ci] ) ;
         TH1F* hmc = (TH1F*) gDirectory->FindObject( hname ) ;
         legend -> AddEntry( hmc, compname[ci] ) ;
      }


      hmcsum->SetXTitle( xtitle ) ;

      hmcsum->SetMarkerStyle(0) ;

      gStyle->SetOptTitle(0) ;
      hmcsum->UseCurrentStyle() ;
      hmcsum->Draw("e") ;
      hmcstack->Draw("hist same") ;
      hmcsum->Draw("esame") ;
      hmcsum->Draw("axis same") ;
      legend->Draw() ;
      TText* title = new TText() ;
      title->SetTextSize(0.040) ;
      TString newtitle( hmcsum->GetTitle() ) ;
      char erasestring[100] ;
      sprintf( erasestring, ", %s", compname[0] ) ;
      newtitle.ReplaceAll( erasestring, "" ) ;
      title->DrawTextNDC( 0.05, 0.95, newtitle ) ;


      dmccan->Update() ;
      dmccan->Draw() ;



      if ( isMTPlot ) {
         TLine* line = new TLine() ;
         line->SetLineStyle(2) ;
         line->DrawLine(100., 0., 100., (gPad->GetUymax()) ) ;
         double int_mtlt100err(0.), int_mtgt100err(0.) ;
         double int_mtlt100 = hmcsum->IntegralAndError(  1, 10, int_mtlt100err) ;
         double int_mtgt100 = hmcsum->IntegralAndError( 11, 25, int_mtgt100err) ;
         double rmt(0.) ;
         if ( int_mtlt100 > 0 ) {
            rmt = int_mtgt100 / int_mtlt100 ;
         }
         double rmt_err(0.) ;
         if ( int_mtgt100>0 && int_mtlt100>0 ) {
            rmt_err = rmt * sqrt( pow(int_mtlt100err/int_mtlt100,2) + pow(int_mtgt100err/int_mtgt100,2) ) ;
         }
         char rmtstring[1000] ;
         sprintf( rmtstring, "R_{MT} = %5.3f #pm %5.3f", rmt, rmt_err ) ;
         TLatex* tlatex = new TLatex() ;
         tlatex -> DrawLatex( 115., 0.9*(gPad->GetUymax()), rmtstring ) ;
      }


      if ( savePdf ) {
         TString dataset( inrootfile ) ;
         dataset.ReplaceAll("rootfiles/dmc_plots_Rmt_","") ;
         dataset.ReplaceAll(".root","") ;
         char filename[10000] ;
         if ( islogy ) {
            sprintf( filename, "outputfiles/%s_%s_logy.pdf", hname_base, dataset.Data() ) ;
         } else {
            sprintf( filename, "outputfiles/%s_%s.pdf", hname_base, dataset.Data() ) ;
         }
         dmccan->SaveAs( filename ) ;
      }

   } // drawSet

   //=======================================================================================


   TH1F* flattenLHB( TH2F* h2 ) {

      if ( h2 == 0x0 ) { printf("\n\n *** drawSet: missing lhb hist\n" ) ; return 0x0 ; }

      int nflatbins(0) ;

      nflatbins = 1 + (nBinsHT+1)*nBinsMET ;

      char flatname[1000] ;
      sprintf( flatname, "%s_flat", h2->GetName() ) ;

      gStyle->SetLabelSize(0.05,"x") ;
      gStyle->SetLabelSize(0.05,"y") ;
      TH1F* hpf = new TH1F( flatname, h2->GetTitle(), nflatbins, 0.5, nflatbins+0.5 ) ;
      hpf->Sumw2() ;
      hpf->SetFillColor( h2->GetFillColor() ) ;

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            int fbi = 1 + mbi*(nBinsHT+1) + hbi + 1 ;
            hpf -> SetBinContent( fbi, h2->GetBinContent( mbi+1, hbi+1 ) ) ;
            hpf -> SetBinError(   fbi, h2->GetBinError  ( mbi+1, hbi+1 ) ) ;
            char binlabel[100] ;
            sprintf(binlabel, "M%d_H%d", mbi+1, hbi+1) ;
            hpf->GetXaxis()->SetBinLabel( fbi, binlabel ) ;
         } // hbi
      } // mbi

      hpf->SetLabelSize(0.04,"x") ;
      hpf->GetXaxis()->LabelsOption("v") ;

      return hpf ;


   } // flattenLHB

   //=======================================================================================
