
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TDirectory.h"

#include "fill_arrays.c"


   void binsum_hsbins_plots( bool doLogy=false ) {

       fill_arrays() ;

       int fillStyle1sig = 1001 ;
       int fillStyle2sig = 1001 ;
       int fillColor1sig = kBlue-9 ;
       int fillColor2sig = kBlue-10 ;
       int fillColor1sig_ff = kGreen-7 ;
       int fillColor2sig_ff = kGreen-10 ;

       gStyle->SetOptStat(0) ;
       gStyle->SetLabelSize(0.12,"x") ;
       gStyle->SetLabelSize(0.08,"y") ;
       gStyle->SetTitleW(0.8) ;
       gStyle->SetTitleH(0.08) ;
       gStyle->SetTitleAlign(13) ;
       gStyle->SetPadLeftMargin(0.14) ;
       gStyle->SetPadRightMargin(0.02) ;
       gStyle->SetPadBottomMargin(0.30) ;
       gStyle->SetLabelOffset(0.02,"y") ;
       gStyle->SetPadGridY(1) ;

       if ( doLogy ) {
          gStyle->SetOptLogy(1) ;
       } else {
          gStyle->SetOptLogy(0) ;
       }

       double xvals[4] = { 1.05, 2.05, 3.05, 4.05 } ;
       double xel[4] = {0.5, 0.5, 0.5, 0.5} ;
       double xeh[4] = {0.5, 0.5, 0.5, 0.5} ;


       TCanvas* cmipbs = (TCanvas*) gDirectory->FindObject("cmipbs") ;
       if ( cmipbs == 0x0 ) {
          cmipbs = new TCanvas("cmipbs", "MIP", 800, 600) ;
       }
       cmipbs->Clear() ;

       cmipbs->Divide(4,1) ;



    //--------

       cmipbs->cd(1) ;

       TH1F* h_data_nb3_hs = new TH1F( "h_data_nb3_hs", "Data, nB=3, HT bins", 4, 0.5, 4.5 ) ;

       h_data_nb3_hs -> GetXaxis() -> SetBinLabel( 1, "HT1, sum of MET2,3" ) ;
       h_data_nb3_hs -> GetXaxis() -> SetBinLabel( 2, "HT2, sum of MET2,3,4" ) ;
       h_data_nb3_hs -> GetXaxis() -> SetBinLabel( 3, "HT3, sum of MET2,3,4" ) ;
       h_data_nb3_hs -> GetXaxis() -> SetBinLabel( 4, "HT4, sum of MET2,3,4" ) ;

       h_data_nb3_hs -> GetXaxis() -> LabelsOption("v") ;

       h_data_nb3_hs -> SetBinContent( 1, data_3b_obs[3][0] ) ;
       h_data_nb3_hs -> SetBinContent( 2, data_3b_obs[3][1] ) ;
       h_data_nb3_hs -> SetBinContent( 3, data_3b_obs[3][2] ) ;
       h_data_nb3_hs -> SetBinContent( 4, data_3b_obs[3][3] ) ;

       h_data_nb3_hs -> SetMarkerStyle(20) ;
       h_data_nb3_hs -> SetMarkerSize(1.2) ;
       h_data_nb3_hs -> SetLineWidth(2) ;
       if ( doLogy ) h_data_nb3_hs -> SetMinimum( 4 ) ;
       h_data_nb3_hs -> SetMaximum( 440 ) ;


       h_data_nb3_hs -> Draw("e") ;


       double nb3_hs_yvals[4]    = { ub_3b_val[3][0], ub_3b_val[3][1], ub_3b_val[3][2], ub_3b_val[3][3] } ;
       double nb3_hs_yerrh[4]    = { ub_3b_p1s[3][0], ub_3b_p1s[3][1], ub_3b_p1s[3][2], ub_3b_p1s[3][3] } ;
       double nb3_hs_yerrl[4]    = { ub_3b_m1s[3][0], ub_3b_m1s[3][1], ub_3b_m1s[3][2], ub_3b_m1s[3][3] } ;
       double nb3_hs_yerrh_2s[4] = { ub_3b_p2s[3][0], ub_3b_p2s[3][1], ub_3b_p2s[3][2], ub_3b_p2s[3][3] } ;
       double nb3_hs_yerrl_2s[4] = { ub_3b_m2s[3][0], ub_3b_m2s[3][1], ub_3b_m2s[3][2], ub_3b_m2s[3][3] } ;


       TGraphAsymmErrors* gr_nb3_hs = new TGraphAsymmErrors( 4, xvals, nb3_hs_yvals, xel, xeh,
            nb3_hs_yerrl, nb3_hs_yerrh ) ;

       TGraphAsymmErrors* gr_nb3_hs_2s = new TGraphAsymmErrors( 4, xvals, nb3_hs_yvals, xel, xeh,
            nb3_hs_yerrl_2s, nb3_hs_yerrh_2s ) ;

       gr_nb3_hs -> SetFillStyle(fillStyle1sig) ;
       gr_nb3_hs -> SetFillColor(fillColor1sig) ;
       gr_nb3_hs -> SetLineColor(4) ;
       gr_nb3_hs -> SetLineWidth(1) ;
       gr_nb3_hs -> SetMarkerColor(4) ;

       gr_nb3_hs_2s -> SetFillStyle(fillStyle2sig) ;
       gr_nb3_hs_2s -> SetFillColor(fillColor2sig) ;
       gr_nb3_hs_2s -> SetLineColor(4) ;
       gr_nb3_hs_2s -> SetMarkerColor(4) ;

       gr_nb3_hs_2s -> Draw("p2") ;
       gr_nb3_hs_2s -> Draw("p") ;
       gr_nb3_hs -> Draw("p2") ;
       gr_nb3_hs -> Draw("p") ;
       h_data_nb3_hs -> Draw("same e") ;



    //--------

       cmipbs->cd(2) ;

       TH1F* h_data_nb3_ms = new TH1F( "h_data_nb3_ms", "Data, nB=3, MET bins", 3, 0.5, 3.5 ) ;

       h_data_nb3_ms -> GetXaxis() -> SetBinLabel( 1, "MET2, sum of HT1,2,3,4" ) ;
       h_data_nb3_ms -> GetXaxis() -> SetBinLabel( 2, "MET3, sum of HT1,2,3,4" ) ;
       h_data_nb3_ms -> GetXaxis() -> SetBinLabel( 3, "MET4, sum of HT2,3,4" ) ;

       h_data_nb3_ms -> GetXaxis() -> LabelsOption("v") ;

       h_data_nb3_ms -> SetBinContent( 1, data_3b_obs[0][4] ) ;
       h_data_nb3_ms -> SetBinContent( 2, data_3b_obs[1][4] ) ;
       h_data_nb3_ms -> SetBinContent( 3, data_3b_obs[2][4] ) ;

       h_data_nb3_ms -> SetMarkerStyle(20) ;
       h_data_nb3_ms -> SetMarkerSize(1.2) ;
       h_data_nb3_ms -> SetLineWidth(2) ;
       if ( doLogy ) h_data_nb3_ms -> SetMinimum( 4 ) ;
       h_data_nb3_ms -> SetMaximum( 580 ) ;


       h_data_nb3_ms -> Draw("e") ;


       double nb3_ms_yvals[3]    = { ub_3b_val[0][4], ub_3b_val[1][4], ub_3b_val[2][4] } ;
       double nb3_ms_yerrh[3]    = { ub_3b_p1s[0][4], ub_3b_p1s[1][4], ub_3b_p1s[2][4] } ;
       double nb3_ms_yerrl[3]    = { ub_3b_m1s[0][4], ub_3b_m1s[1][4], ub_3b_m1s[2][4] } ;
       double nb3_ms_yerrh_2s[3] = { ub_3b_p2s[0][4], ub_3b_p2s[1][4], ub_3b_p2s[2][4] } ;
       double nb3_ms_yerrl_2s[3] = { ub_3b_m2s[0][4], ub_3b_m2s[1][4], ub_3b_m2s[2][4] } ;

       TGraphAsymmErrors* gr_nb3_ms = new TGraphAsymmErrors( 4, xvals, nb3_ms_yvals, xel, xeh,
            nb3_ms_yerrl, nb3_ms_yerrh ) ;

       TGraphAsymmErrors* gr_nb3_ms_2s = new TGraphAsymmErrors( 4, xvals, nb3_ms_yvals, xel, xeh,
            nb3_ms_yerrl_2s, nb3_ms_yerrh_2s ) ;

       gr_nb3_ms -> SetFillStyle(fillStyle1sig) ;
       gr_nb3_ms -> SetFillColor(fillColor1sig) ;
       gr_nb3_ms -> SetLineColor(4) ;
       gr_nb3_ms -> SetLineWidth(1) ;
       gr_nb3_ms -> SetMarkerColor(4) ;

       gr_nb3_ms_2s -> SetFillStyle(fillStyle2sig) ;
       gr_nb3_ms_2s -> SetFillColor(fillColor2sig) ;
       gr_nb3_ms_2s -> SetLineColor(4) ;
       gr_nb3_ms_2s -> SetMarkerColor(4) ;

       gr_nb3_ms_2s -> Draw("p2") ;
       gr_nb3_ms_2s -> Draw("p") ;
       gr_nb3_ms -> Draw("p2") ;
       gr_nb3_ms -> Draw("p") ;
       h_data_nb3_ms -> Draw("same e") ;



    //++++++++++++++++

       cmipbs->cd(3) ;



       h_data_nb3_hs -> Draw("e") ;


       double nb3_hs_ff_yvals[4]    = { ff_3b_val[3][0], ff_3b_val[3][1], ff_3b_val[3][2], ff_3b_val[3][3] } ;
       double nb3_hs_ff_yerrh[4]    = { ff_3b_p1s[3][0], ff_3b_p1s[3][1], ff_3b_p1s[3][2], ff_3b_p1s[3][3] } ;
       double nb3_hs_ff_yerrl[4]    = { ff_3b_m1s[3][0], ff_3b_m1s[3][1], ff_3b_m1s[3][2], ff_3b_m1s[3][3] } ;
       double nb3_hs_ff_yerrh_2s[4] = { ff_3b_p2s[3][0], ff_3b_p2s[3][1], ff_3b_p2s[3][2], ff_3b_p2s[3][3] } ;
       double nb3_hs_ff_yerrl_2s[4] = { ff_3b_m2s[3][0], ff_3b_m2s[3][1], ff_3b_m2s[3][2], ff_3b_m2s[3][3] } ;



       TGraphAsymmErrors* gr_nb3_hs_ff = new TGraphAsymmErrors( 4, xvals, nb3_hs_ff_yvals, xel, xeh,
            nb3_hs_ff_yerrl, nb3_hs_ff_yerrh ) ;

       TGraphAsymmErrors* gr_nb3_hs_ff_2s = new TGraphAsymmErrors( 4, xvals, nb3_hs_ff_yvals, xel, xeh,
            nb3_hs_ff_yerrl_2s, nb3_hs_ff_yerrh_2s ) ;

       gr_nb3_hs_ff -> SetFillStyle(fillStyle1sig) ;
       gr_nb3_hs_ff -> SetFillColor(fillColor1sig_ff) ;
       gr_nb3_hs_ff -> SetLineColor(kGreen+2) ;
       gr_nb3_hs_ff -> SetLineWidth(1) ;
       gr_nb3_hs_ff -> SetMarkerColor(4) ;

       gr_nb3_hs_ff_2s -> SetFillStyle(fillStyle2sig) ;
       gr_nb3_hs_ff_2s -> SetFillColor(fillColor2sig_ff) ;
       gr_nb3_hs_ff_2s -> SetLineColor(kGreen+2) ;
       gr_nb3_hs_ff_2s -> SetMarkerColor(4) ;

       gr_nb3_hs_ff_2s -> Draw("p2") ;
       gr_nb3_hs_ff_2s -> Draw("p") ;
       gr_nb3_hs_ff -> Draw("p2") ;
       gr_nb3_hs_ff -> Draw("p") ;
       h_data_nb3_hs -> Draw("same e") ;



    //--------

       cmipbs->cd(4) ;

       h_data_nb3_ms -> Draw("e") ;


       double nb3_ms_ff_yvals[3]    = { ff_3b_val[0][4], ff_3b_val[1][4], ff_3b_val[2][4] } ;
       double nb3_ms_ff_yerrh[3]    = { ff_3b_p1s[0][4], ff_3b_p1s[1][4], ff_3b_p1s[2][4] } ;
       double nb3_ms_ff_yerrl[3]    = { ff_3b_m1s[0][4], ff_3b_m1s[1][4], ff_3b_m1s[2][4] } ;
       double nb3_ms_ff_yerrh_2s[3] = { ff_3b_p2s[0][4], ff_3b_p2s[1][4], ff_3b_p2s[2][4] } ;
       double nb3_ms_ff_yerrl_2s[3] = { ff_3b_m2s[0][4], ff_3b_m2s[1][4], ff_3b_m2s[2][4] } ;

       TGraphAsymmErrors* gr_nb3_ms_ff = new TGraphAsymmErrors( 4, xvals, nb3_ms_ff_yvals, xel, xeh,
            nb3_ms_ff_yerrl, nb3_ms_ff_yerrh ) ;

       TGraphAsymmErrors* gr_nb3_ms_ff_2s = new TGraphAsymmErrors( 4, xvals, nb3_ms_ff_yvals, xel, xeh,
            nb3_ms_ff_yerrl_2s, nb3_ms_ff_yerrh_2s ) ;

       gr_nb3_ms_ff -> SetFillStyle(fillStyle1sig) ;
       gr_nb3_ms_ff -> SetFillColor(fillColor1sig_ff) ;
       gr_nb3_ms_ff -> SetLineColor(kGreen+2) ;
       gr_nb3_ms_ff -> SetLineWidth(1) ;
       gr_nb3_ms_ff -> SetMarkerColor(4) ;

       gr_nb3_ms_ff_2s -> SetFillStyle(fillStyle2sig) ;
       gr_nb3_ms_ff_2s -> SetFillColor(fillColor2sig_ff) ;
       gr_nb3_ms_ff_2s -> SetLineColor(kGreen+2) ;
       gr_nb3_ms_ff_2s -> SetMarkerColor(4) ;

       gr_nb3_ms_ff_2s -> Draw("p2") ;
       gr_nb3_ms_ff_2s -> Draw("p") ;
       gr_nb3_ms_ff -> Draw("p2") ;
       gr_nb3_ms_ff -> Draw("p") ;
       h_data_nb3_ms -> Draw("same e") ;


       if ( doLogy ) {
          cmipbs->SaveAs("fitresult-binsums-hsbins-2sig-19fb-logy.pdf") ;
       } else {
          cmipbs->SaveAs("fitresult-binsums-hsbins-2sig-19fb.pdf") ;
       }



   }
