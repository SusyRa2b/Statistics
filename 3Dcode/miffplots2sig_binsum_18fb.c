
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TDirectory.h"



   void miffplots2sig_binsum_18fb( bool doLogy=false ) {

       int fillStyle1sig = 1 ;
       int fillStyle2sig = 1 ;
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

       h_data_nb3_hs -> SetBinContent( 1, 162 ) ;
       h_data_nb3_hs -> SetBinContent( 2, 210 ) ;
       h_data_nb3_hs -> SetBinContent( 3, 23 ) ;
       h_data_nb3_hs -> SetBinContent( 4, 22 ) ;

       h_data_nb3_hs -> SetMarkerStyle(20) ;
       h_data_nb3_hs -> SetMarkerSize(1.2) ;
       h_data_nb3_hs -> SetLineWidth(2) ;
       if ( doLogy ) h_data_nb3_hs -> SetMinimum( 4 ) ;
       h_data_nb3_hs -> SetMaximum( 400 ) ;


       h_data_nb3_hs -> Draw("e") ;


       double nb3_hs_yvals[4] = { 118.0, 178.9, 32.7, 11.9 } ;
       double nb3_hs_yerrh[4] = { 28, 59, 7.9, 3.8 } ;
       double nb3_hs_yerrl[4] = { 19, 34, 6.5, 3.0 } ;
       double nb3_hs_yerrh_2s[4] = { 82, 193, 17.8, 8.5 } ;
       double nb3_hs_yerrl_2s[4] = { 35, 61, 11.8, 5.4 } ;


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

       h_data_nb3_ms -> SetBinContent( 1, 347 ) ;
       h_data_nb3_ms -> SetBinContent( 2, 56 ) ;
       h_data_nb3_ms -> SetBinContent( 3, 14 ) ;

       h_data_nb3_ms -> SetMarkerStyle(20) ;
       h_data_nb3_ms -> SetMarkerSize(1.2) ;
       h_data_nb3_ms -> SetLineWidth(2) ;
       if ( doLogy ) h_data_nb3_ms -> SetMinimum( 4 ) ;
       h_data_nb3_ms -> SetMaximum( 530 ) ;


       h_data_nb3_ms -> Draw("e") ;


       double nb3_ms_yvals[3] = { 289, 42.2, 10.0 } ;
       double nb3_ms_yerrh[3] = { 64, 8.7, 3.1 } ;
       double nb3_ms_yerrl[3] = { 43, 7.1, 2.5 } ;
       double nb3_ms_yerrh_2s[3] = { 199, 19, 7 } ;
       double nb3_ms_yerrl_2s[3] = { 79, 13, 4.5 } ;

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


       double nb3_hs_ff_yvals[4] = { 157, 206, 27.6, 16.7 } ;
       double nb3_hs_ff_yerrh[4] = { 12.6, 13.9, 4.2, 3.1 } ;
       double nb3_hs_ff_yerrl[4] = { 11.6, 13.3, 3.7, 2.7 } ;
       double nb3_hs_ff_yerrh_2s[4] = { 25.7, 28.3, 8.9, 6.6 } ;
       double nb3_hs_ff_yerrl_2s[4] = { 22.6, 25.8, 7.1, 5.1 } ;


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


       double nb3_ms_ff_yvals[3] = { 344, 51.3, 11.8 } ;
       double nb3_ms_ff_yerrh[3] = { 18.6, 6.0, 2.5 } ;
       double nb3_ms_ff_yerrl[3] = { 17.4, 5.4, 2.1 } ;
       double nb3_ms_ff_yerrh_2s[3] = { 37, 12.3, 5.3 } ;
       double nb3_ms_ff_yerrl_2s[3] = { 34, 10.4, 3.9 } ;

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





   }
