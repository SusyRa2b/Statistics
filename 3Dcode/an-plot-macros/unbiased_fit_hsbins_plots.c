
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TDirectory.h"

#include "fill_arrays.c"

   void unbiased_fit_hsbins_plots( bool doLogy=false ) {

       fill_arrays() ;

       int fillStyle1sig = 1 ;
       int fillStyle2sig = 1 ;
       int fillColor1sig = kBlue-9 ;
       int fillColor2sig = kBlue-10 ;

       gStyle->SetOptStat(0) ;
       gStyle->SetLabelSize(0.12,"x") ;
       gStyle->SetLabelSize(0.08,"y") ;
       gStyle->SetTitleW(0.8) ;
       gStyle->SetTitleH(0.08) ;
       gStyle->SetTitleAlign(13) ;
       gStyle->SetPadLeftMargin(0.14) ;
       gStyle->SetPadRightMargin(0.02) ;
       gStyle->SetLabelOffset(0.02,"y") ;
       gStyle->SetPadGridY(1) ;

       if ( doLogy ) {
          gStyle->SetOptLogy(1) ;
       } else {
          gStyle->SetOptLogy(0) ;
       }

       gDirectory->Delete("h*") ;

       double xvals[4] = { 1.05, 2.05, 3.05, 4.05 } ;
       double xel[4] = {0.5, 0.5, 0.5, 0.5} ;
       double xeh[4] = {0.5, 0.5, 0.5, 0.5} ;


       TCanvas* cmip = (TCanvas*) gDirectory->FindObject("cmip") ;
       if ( cmip == 0x0 ) {
          cmip = new TCanvas("cmip", "MIP", 800, 600) ;
       }
       cmip->Clear() ;

       cmip->Divide(4,1) ;



    //--------

       cmip->cd(1) ;

       TH1F* h_data_nb2_m4 = new TH1F( "h_data_nb2_m4", "Data, nB=2, MET4", 4, 0.5, 4.5 ) ;

       h_data_nb2_m4 -> GetXaxis() -> SetBinLabel( 1, "HT1" ) ;
       h_data_nb2_m4 -> GetXaxis() -> SetBinLabel( 2, "HT2" ) ;
       h_data_nb2_m4 -> GetXaxis() -> SetBinLabel( 3, "HT3" ) ;
       h_data_nb2_m4 -> GetXaxis() -> SetBinLabel( 4, "HT4" ) ;

       h_data_nb2_m4 -> SetBinContent( 1, 0 ) ;
       h_data_nb2_m4 -> SetBinContent( 2, data_2b_obs[0] ) ;
       h_data_nb2_m4 -> SetBinContent( 3, data_2b_obs[1] ) ;
       h_data_nb2_m4 -> SetBinContent( 4, data_2b_obs[2] ) ;

       h_data_nb2_m4 -> SetMarkerStyle(20) ;
       h_data_nb2_m4 -> SetMarkerSize(1.2) ;
       h_data_nb2_m4 -> SetLineWidth(2) ;
       if ( doLogy) h_data_nb2_m4 -> SetMinimum( 10 ) ;
       h_data_nb2_m4 -> SetMaximum( 110 ) ;


       h_data_nb2_m4 -> Draw("e") ;


       double nb2_m4_yvals[4]    = { -99., ub_2b_val[0], ub_2b_val[1], ub_2b_val[2] } ;
       double nb2_m4_yerrh[4]    = {  0.0, ub_2b_p1s[0], ub_2b_p1s[1], ub_2b_p1s[2] } ;
       double nb2_m4_yerrl[4]    = {  0.0, ub_2b_m1s[0], ub_2b_m1s[1], ub_2b_m1s[2] } ;
       double nb2_m4_yerrh_2s[4] = {  0.0, ub_2b_p2s[0], ub_2b_p2s[1], ub_2b_p2s[2] } ;
       double nb2_m4_yerrl_2s[4] = {  0.0, ub_2b_m2s[0], ub_2b_m2s[1], ub_2b_m2s[2] } ;

       TGraphAsymmErrors* gr_nb2_m4 = new TGraphAsymmErrors( 4, xvals, nb2_m4_yvals, xel, xeh,
            nb2_m4_yerrl, nb2_m4_yerrh ) ;

       TGraphAsymmErrors* gr_nb2_m4_2s = new TGraphAsymmErrors( 4, xvals, nb2_m4_yvals, xel, xeh,
            nb2_m4_yerrl_2s, nb2_m4_yerrh_2s ) ;

       gr_nb2_m4 -> SetFillStyle(fillStyle1sig) ;
       gr_nb2_m4 -> SetFillColor(fillColor1sig) ;
       gr_nb2_m4 -> SetLineColor(4) ;
       gr_nb2_m4 -> SetLineWidth(1) ;
       gr_nb2_m4 -> SetMarkerColor(4) ;
       //gr_nb2_m4 -> SetMarkerStyle(24) ;

       gr_nb2_m4_2s -> SetFillStyle(fillStyle2sig) ;
       gr_nb2_m4_2s -> SetFillColor(fillColor2sig) ;
       gr_nb2_m4_2s -> SetLineColor(4) ;
       gr_nb2_m4_2s -> SetMarkerColor(4) ;
       //gr_nb2_m4_2s -> SetMarkerStyle(24) ;

       gr_nb2_m4_2s -> Draw("p2") ;
       gr_nb2_m4_2s -> Draw("p") ;
       gr_nb2_m4 -> Draw("p2") ;
       gr_nb2_m4 -> Draw("p") ;
       h_data_nb2_m4 -> Draw("same e") ;



    //--------

       cmip->cd(2) ;

       TH1F* h_data_nb3_m2 = new TH1F( "h_data_nb3_m2", "Data, nB=3, MET2", 4, 0.5, 4.5 ) ;

       h_data_nb3_m2 -> GetXaxis() -> SetBinLabel( 1, "HT1" ) ;
       h_data_nb3_m2 -> GetXaxis() -> SetBinLabel( 2, "HT2" ) ;
       h_data_nb3_m2 -> GetXaxis() -> SetBinLabel( 3, "HT3" ) ;
       h_data_nb3_m2 -> GetXaxis() -> SetBinLabel( 4, "HT4" ) ;

       h_data_nb3_m2 -> SetBinContent( 1, data_3b_obs[0][0] ) ;
       h_data_nb3_m2 -> SetBinContent( 2, data_3b_obs[0][1] ) ;
       h_data_nb3_m2 -> SetBinContent( 3, data_3b_obs[0][2] ) ;
       h_data_nb3_m2 -> SetBinContent( 4, data_3b_obs[0][3] ) ;

       h_data_nb3_m2 -> SetMarkerStyle(20) ;
       h_data_nb3_m2 -> SetMarkerSize(1.2) ;
       h_data_nb3_m2 -> SetLineWidth(2) ;
       if ( doLogy) h_data_nb3_m2 -> SetMinimum( 4 ) ;
       h_data_nb3_m2 -> SetMaximum( 400 ) ;


       h_data_nb3_m2 -> Draw("e") ;


       double nb3_m2_yvals[4]    = { ub_3b_val[0][0], ub_3b_val[0][1], ub_3b_val[0][2], ub_3b_val[0][3] } ;
       double nb3_m2_yerrh[4]    = { ub_3b_p1s[0][0], ub_3b_p1s[0][1], ub_3b_p1s[0][2], ub_3b_p1s[0][3] } ;
       double nb3_m2_yerrl[4]    = { ub_3b_m1s[0][0], ub_3b_m1s[0][1], ub_3b_m1s[0][2], ub_3b_m1s[0][3] } ;
       double nb3_m2_yerrh_2s[4] = { ub_3b_p2s[0][0], ub_3b_p2s[0][1], ub_3b_p2s[0][2], ub_3b_p2s[0][3] } ;
       double nb3_m2_yerrl_2s[4] = { ub_3b_m2s[0][0], ub_3b_m2s[0][1], ub_3b_m2s[0][2], ub_3b_m2s[0][3] } ;


       TGraphAsymmErrors* gr_nb3_m2 = new TGraphAsymmErrors( 4, xvals, nb3_m2_yvals, xel, xeh,
            nb3_m2_yerrl, nb3_m2_yerrh ) ;

       TGraphAsymmErrors* gr_nb3_m2_2s = new TGraphAsymmErrors( 4, xvals, nb3_m2_yvals, xel, xeh,
            nb3_m2_yerrl_2s, nb3_m2_yerrh_2s ) ;


       gr_nb3_m2 -> SetFillStyle(fillStyle1sig) ;
       gr_nb3_m2 -> SetFillColor(fillColor1sig) ;
       gr_nb3_m2 -> SetLineColor(4) ;
       gr_nb3_m2 -> SetLineWidth(1) ;
       gr_nb3_m2 -> SetMarkerColor(4) ;
       //gr_nb3_m2 -> SetMarkerStyle(24) ;

       gr_nb3_m2_2s -> SetFillStyle(fillStyle2sig) ;
       gr_nb3_m2_2s -> SetFillColor(fillColor2sig) ;
       gr_nb3_m2_2s -> SetLineColor(4) ;
       gr_nb3_m2_2s -> SetMarkerColor(4) ;
       //gr_nb3_m2_2s -> SetMarkerStyle(24) ;


       gr_nb3_m2_2s -> Draw("p2") ;
       gr_nb3_m2_2s -> Draw("p") ;
       gr_nb3_m2 -> Draw("p2") ;
       gr_nb3_m2 -> Draw("p") ;
       h_data_nb3_m2 -> Draw("same e") ;



    //--------

       cmip->cd(3) ;

       TH1F* h_data_nb3_m3 = new TH1F( "h_data_nb3_m3", "Data, nB=3, MET3", 4, 0.5, 4.5 ) ;

       h_data_nb3_m3 -> GetXaxis() -> SetBinLabel( 1, "HT1" ) ;
       h_data_nb3_m3 -> GetXaxis() -> SetBinLabel( 2, "HT2" ) ;
       h_data_nb3_m3 -> GetXaxis() -> SetBinLabel( 3, "HT3" ) ;
       h_data_nb3_m3 -> GetXaxis() -> SetBinLabel( 4, "HT4" ) ;

       h_data_nb3_m3 -> SetBinContent( 1, data_3b_obs[1][0] ) ;
       h_data_nb3_m3 -> SetBinContent( 2, data_3b_obs[1][1] ) ;
       h_data_nb3_m3 -> SetBinContent( 3, data_3b_obs[1][2] ) ;
       h_data_nb3_m3 -> SetBinContent( 4, data_3b_obs[1][3] ) ;

       h_data_nb3_m3 -> SetMarkerStyle(20) ;
       h_data_nb3_m3 -> SetMarkerSize(1.2) ;
       h_data_nb3_m3 -> SetLineWidth(2) ;
       if ( doLogy) h_data_nb3_m3 -> SetMinimum( 0.3 ) ;
       h_data_nb3_m3 -> SetMaximum( 50 ) ;


       h_data_nb3_m3 -> Draw("e") ;


       double nb3_m3_yvals[4]    = { ub_3b_val[1][0], ub_3b_val[1][1], ub_3b_val[1][2], ub_3b_val[1][3] } ;
       double nb3_m3_yerrh[4]    = { ub_3b_p1s[1][0], ub_3b_p1s[1][1], ub_3b_p1s[1][2], ub_3b_p1s[1][3] } ;
       double nb3_m3_yerrl[4]    = { ub_3b_m1s[1][0], ub_3b_m1s[1][1], ub_3b_m1s[1][2], ub_3b_m1s[1][3] } ;
       double nb3_m3_yerrh_2s[4] = { ub_3b_p2s[1][0], ub_3b_p2s[1][1], ub_3b_p2s[1][2], ub_3b_p2s[1][3] } ;
       double nb3_m3_yerrl_2s[4] = { ub_3b_m2s[1][0], ub_3b_m2s[1][1], ub_3b_m2s[1][2], ub_3b_m2s[1][3] } ;

       TGraphAsymmErrors* gr_nb3_m3 = new TGraphAsymmErrors( 4, xvals, nb3_m3_yvals, xel, xeh,
            nb3_m3_yerrl, nb3_m3_yerrh ) ;

       TGraphAsymmErrors* gr_nb3_m3_2s = new TGraphAsymmErrors( 4, xvals, nb3_m3_yvals, xel, xeh,
            nb3_m3_yerrl_2s, nb3_m3_yerrh_2s ) ;

       gr_nb3_m3 -> SetFillStyle(fillStyle1sig) ;
       gr_nb3_m3 -> SetFillColor(fillColor1sig) ;
       gr_nb3_m3 -> SetLineColor(4) ;
       gr_nb3_m3 -> SetLineWidth(1) ;
       gr_nb3_m3 -> SetMarkerColor(4) ;
       //gr_nb3_m3 -> SetMarkerStyle(24) ;

       gr_nb3_m3_2s -> SetFillStyle(fillStyle2sig) ;
       gr_nb3_m3_2s -> SetFillColor(fillColor2sig) ;
       gr_nb3_m3_2s -> SetLineColor(4) ;
       gr_nb3_m3_2s -> SetMarkerColor(4) ;
       //gr_nb3_m3_2s -> SetMarkerStyle(24) ;

       gr_nb3_m3_2s -> Draw("p2") ;
       gr_nb3_m3_2s -> Draw("p") ;
       gr_nb3_m3 -> Draw("p2") ;
       gr_nb3_m3 -> Draw("p") ;
       h_data_nb3_m3 -> Draw("same e") ;




    //--------

       cmip->cd(4) ;

       TH1F* h_data_nb3_m4 = new TH1F( "h_data_nb3_m4", "Data, nB=3, MET4", 4, 0.5, 4.5 ) ;

       h_data_nb3_m4 -> GetXaxis() -> SetBinLabel( 1, "HT1" ) ;
       h_data_nb3_m4 -> GetXaxis() -> SetBinLabel( 2, "HT2" ) ;
       h_data_nb3_m4 -> GetXaxis() -> SetBinLabel( 3, "HT3" ) ;
       h_data_nb3_m4 -> GetXaxis() -> SetBinLabel( 4, "HT4" ) ;

       h_data_nb3_m4 -> SetBinContent( 1, 0 ) ;
       h_data_nb3_m4 -> SetBinContent( 2, data_3b_obs[2][1] ) ;
       h_data_nb3_m4 -> SetBinContent( 3, data_3b_obs[2][2] ) ;
       h_data_nb3_m4 -> SetBinContent( 4, data_3b_obs[2][3] ) ;

       h_data_nb3_m4 -> SetMarkerStyle(20) ;
       h_data_nb3_m4 -> SetMarkerSize(1.2) ;
       h_data_nb3_m4 -> SetLineWidth(2) ;
       if ( doLogy) h_data_nb3_m4 -> SetMinimum( 0.1 ) ;
       h_data_nb3_m4 -> SetMaximum( 20 ) ;


       h_data_nb3_m4 -> Draw("e") ;


       double nb3_m4_yvals[4]    = { -99., ub_3b_val[2][1], ub_3b_val[2][2], ub_3b_val[2][3] } ;
       double nb3_m4_yerrh[4]    = {  0.0, ub_3b_p1s[2][1], ub_3b_p1s[2][2], ub_3b_p1s[2][3] } ;
       double nb3_m4_yerrl[4]    = {  0.0, ub_3b_m1s[2][1], ub_3b_m1s[2][2], ub_3b_m1s[2][3] } ;
       double nb3_m4_yerrh_2s[4] = {  0.0, ub_3b_p2s[2][1], ub_3b_p2s[2][2], ub_3b_p2s[2][3] } ;
       double nb3_m4_yerrl_2s[4] = {  0.0, ub_3b_m2s[2][1], ub_3b_m2s[2][2], ub_3b_m2s[2][3] } ;


       TGraphAsymmErrors* gr_nb3_m4 = new TGraphAsymmErrors( 4, xvals, nb3_m4_yvals, xel, xeh,
            nb3_m4_yerrl, nb3_m4_yerrh ) ;

       TGraphAsymmErrors* gr_nb3_m4_2s = new TGraphAsymmErrors( 4, xvals, nb3_m4_yvals, xel, xeh,
            nb3_m4_yerrl_2s, nb3_m4_yerrh_2s ) ;

       gr_nb3_m4 -> SetFillStyle(fillStyle1sig) ;
       gr_nb3_m4 -> SetFillColor(fillColor1sig) ;
       gr_nb3_m4 -> SetLineColor(4) ;
       gr_nb3_m4 -> SetLineWidth(1) ;
       gr_nb3_m4 -> SetMarkerColor(4) ;
       //gr_nb3_m4 -> SetMarkerStyle(24) ;

       gr_nb3_m4_2s -> SetFillStyle(fillStyle2sig) ;
       gr_nb3_m4_2s -> SetFillColor(fillColor2sig) ;
       gr_nb3_m4_2s -> SetLineColor(4) ;
       gr_nb3_m4_2s -> SetMarkerColor(4) ;
       //gr_nb3_m4_2s -> SetMarkerStyle(24) ;

       gr_nb3_m4_2s -> Draw("p2") ;
       gr_nb3_m4_2s -> Draw("p") ;
       gr_nb3_m4 -> Draw("p2") ;
       gr_nb3_m4 -> Draw("p") ;
       h_data_nb3_m4 -> Draw("same e") ;

       cmip->Update() ;
       cmip->Draw() ;

       if ( doLogy ) {
          cmip->SaveAs("fitresult-unbiased-hsbins-18fb-logy.pdf" ) ;
       } else {
          cmip->SaveAs("fitresult-unbiased-hsbins-18fb.pdf" ) ;
       }


   }
