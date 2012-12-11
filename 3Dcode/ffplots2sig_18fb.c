
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TDirectory.h"



   void ffplots2sig_18fb( bool doLogy=false ) {

       int fillStyle1sig = 1 ;
       int fillStyle2sig = 1 ;
       int fillColor1sig = kGreen-7 ;
       int fillColor2sig = kGreen-10 ;

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

       gDirectory->Delete("h*") ;

       if ( doLogy ) {
          gStyle->SetOptLogy(1) ;
       } else {
          gStyle->SetOptLogy(0) ;
       }

       double xvals[4] = { 1.05, 2.05, 3.05, 4.05 } ;
       double xel[4] = {0.5, 0.5, 0.5, 0.5} ;
       double xeh[4] = {0.5, 0.5, 0.5, 0.5} ;


       TCanvas* cffp = (TCanvas*) gDirectory->FindObject("cffp") ;
       if ( cffp == 0x0 ) {
          cffp = new TCanvas("cffp", "MIP", 800, 600) ;
       }
       cffp->Clear() ;

       cffp->Divide(4,1) ;



    //--------

       cffp->cd(1) ;

       TH1F* h_data_nb2_m4 = new TH1F( "h_data_nb2_m4", "Data, nB=2, MET4", 4, 0.5, 4.5 ) ;

       h_data_nb2_m4 -> GetXaxis() -> SetBinLabel( 1, "HT1" ) ;
       h_data_nb2_m4 -> GetXaxis() -> SetBinLabel( 2, "HT2" ) ;
       h_data_nb2_m4 -> GetXaxis() -> SetBinLabel( 3, "HT3" ) ;
       h_data_nb2_m4 -> GetXaxis() -> SetBinLabel( 4, "HT4" ) ;

       h_data_nb2_m4 -> SetBinContent( 1, 0 ) ;
       h_data_nb2_m4 -> SetBinContent( 2, 57 ) ;
       h_data_nb2_m4 -> SetBinContent( 3, 17 ) ;
       h_data_nb2_m4 -> SetBinContent( 4, 19 ) ;

       h_data_nb2_m4 -> SetMarkerStyle(20) ;
       h_data_nb2_m4 -> SetMarkerSize(1.2) ;
       h_data_nb2_m4 -> SetLineWidth(2) ;
       if ( doLogy ) h_data_nb2_m4 -> SetMinimum( 10 ) ;
       h_data_nb2_m4 -> SetMaximum( 100 ) ;


       h_data_nb2_m4 -> Draw("e") ;


       double nb2_m4_yvals[4] = { -99., 63.6, 18.4, 18.7 } ;
       double nb2_m4_yerrh[4] = { 0.0, 6.1, 3.1, 3.3 } ;
       double nb2_m4_yerrl[4] = { 0.0, 5.6, 3.1, 3.3 } ;
       double nb2_m4_yerrh_2s[4] = { 0.0, 12.5, 6.5, 7.1 } ;
       double nb2_m4_yerrl_2s[4] = { 0.0, 10.9, 5.2, 5.3 } ;

       TGraphAsymmErrors* gr_nb2_m4 = new TGraphAsymmErrors( 4, xvals, nb2_m4_yvals, xel, xeh,
            nb2_m4_yerrl, nb2_m4_yerrh ) ;

       TGraphAsymmErrors* gr_nb2_m4_2s = new TGraphAsymmErrors( 4, xvals, nb2_m4_yvals, xel, xeh,
            nb2_m4_yerrl_2s, nb2_m4_yerrh_2s ) ;

       gr_nb2_m4 -> SetFillStyle(fillStyle1sig) ;
       gr_nb2_m4 -> SetFillColor(fillColor1sig) ;
       gr_nb2_m4 -> SetLineColor(kGreen+2) ;
       gr_nb2_m4 -> SetLineWidth(1) ;
       gr_nb2_m4 -> SetMarkerColor(4) ;
       //gr_nb2_m4 -> SetMarkerStyle(24) ;

       gr_nb2_m4_2s -> SetFillStyle(fillStyle2sig) ;
       gr_nb2_m4_2s -> SetFillColor(fillColor2sig) ;
       gr_nb2_m4_2s -> SetLineColor(kGreen+2) ;
       gr_nb2_m4_2s -> SetMarkerColor(4) ;
       //gr_nb2_m4_2s -> SetMarkerStyle(24) ;

       gr_nb2_m4_2s -> Draw("p2") ;
       gr_nb2_m4_2s -> Draw("p") ;
       gr_nb2_m4 -> Draw("p2") ;
       gr_nb2_m4 -> Draw("p") ;
       h_data_nb2_m4 -> Draw("same e") ;



    //--------

       cffp->cd(2) ;

       TH1F* h_data_nb3_m2 = new TH1F( "h_data_nb3_m2", "Data, nB=3, MET2", 4, 0.5, 4.5 ) ;

       h_data_nb3_m2 -> GetXaxis() -> SetBinLabel( 1, "HT1" ) ;
       h_data_nb3_m2 -> GetXaxis() -> SetBinLabel( 2, "HT2" ) ;
       h_data_nb3_m2 -> GetXaxis() -> SetBinLabel( 3, "HT3" ) ;
       h_data_nb3_m2 -> GetXaxis() -> SetBinLabel( 4, "HT4" ) ;

       h_data_nb3_m2 -> SetBinContent( 1, 149 ) ;
       h_data_nb3_m2 -> SetBinContent( 2, 168 ) ;
       h_data_nb3_m2 -> SetBinContent( 3, 16 ) ;
       h_data_nb3_m2 -> SetBinContent( 4, 14 ) ;

       h_data_nb3_m2 -> SetMarkerStyle(20) ;
       h_data_nb3_m2 -> SetMarkerSize(1.2) ;
       h_data_nb3_m2 -> SetLineWidth(2) ;
       if ( doLogy ) h_data_nb3_m2 -> SetMinimum( 4 ) ;
       h_data_nb3_m2 -> SetMaximum( 370 ) ;


       h_data_nb3_m2 -> Draw("e") ;


       double nb3_m2_yvals[4] = { 144.8, 167.2, 20.3, 12.1 } ;
       double nb3_m2_yerrh[4] = { 12.3, 12.8, 3.7, 2.7 } ;
       double nb3_m2_yerrl[4] = { 11.3, 12.0, 3.2, 2.3} ;
       double nb3_m2_yerrh_2s[4] = { 25, 27, 7.9, 5.8 } ;
       double nb3_m2_yerrl_2s[4] = { 22, 23, 6.1, 4.3 } ;

       TGraphAsymmErrors* gr_nb3_m2 = new TGraphAsymmErrors( 4, xvals, nb3_m2_yvals, xel, xeh,
            nb3_m2_yerrl, nb3_m2_yerrh ) ;

       TGraphAsymmErrors* gr_nb3_m2_2s = new TGraphAsymmErrors( 4, xvals, nb3_m2_yvals, xel, xeh,
            nb3_m2_yerrl_2s, nb3_m2_yerrh_2s ) ;


       gr_nb3_m2 -> SetFillStyle(fillStyle1sig) ;
       gr_nb3_m2 -> SetFillColor(fillColor1sig) ;
       gr_nb3_m2 -> SetLineColor(kGreen+2) ;
       gr_nb3_m2 -> SetLineWidth(1) ;
       gr_nb3_m2 -> SetMarkerColor(4) ;
       //gr_nb3_m2 -> SetMarkerStyle(24) ;

       gr_nb3_m2_2s -> SetFillStyle(fillStyle2sig) ;
       gr_nb3_m2_2s -> SetFillColor(fillColor2sig) ;
       gr_nb3_m2_2s -> SetLineColor(kGreen+2) ;
       gr_nb3_m2_2s -> SetMarkerColor(4) ;
       //gr_nb3_m2_2s -> SetMarkerStyle(24) ;


       gr_nb3_m2_2s -> Draw("p2") ;
       gr_nb3_m2_2s -> Draw("p") ;
       gr_nb3_m2 -> Draw("p2") ;
       gr_nb3_m2 -> Draw("p") ;
       h_data_nb3_m2 -> Draw("same e") ;



    //--------

       cffp->cd(3) ;

       TH1F* h_data_nb3_m3 = new TH1F( "h_data_nb3_m3", "Data, nB=3, MET3", 4, 0.5, 4.5 ) ;

       h_data_nb3_m3 -> GetXaxis() -> SetBinLabel( 1, "HT1" ) ;
       h_data_nb3_m3 -> GetXaxis() -> SetBinLabel( 2, "HT2" ) ;
       h_data_nb3_m3 -> GetXaxis() -> SetBinLabel( 3, "HT3" ) ;
       h_data_nb3_m3 -> GetXaxis() -> SetBinLabel( 4, "HT4" ) ;

       h_data_nb3_m3 -> SetBinContent( 1, 13 ) ;
       h_data_nb3_m3 -> SetBinContent( 2, 34 ) ;
       h_data_nb3_m3 -> SetBinContent( 3, 5 ) ;
       h_data_nb3_m3 -> SetBinContent( 4, 4 ) ;

       h_data_nb3_m3 -> SetMarkerStyle(20) ;
       h_data_nb3_m3 -> SetMarkerSize(1.2) ;
       h_data_nb3_m3 -> SetLineWidth(2) ;
       if ( doLogy ) h_data_nb3_m3 -> SetMinimum( 0.3 ) ;
       h_data_nb3_m3 -> SetMaximum( 45 ) ;


       h_data_nb3_m3 -> Draw("e") ;


       double nb3_m3_yvals[4] = { 12.5, 30.4, 5.6, 2.9 } ;
       double nb3_m3_yerrh[4] = {  2.7, 4.8, 1.8, 1.4 } ;
       double nb3_m3_yerrl[4] = {  2.3, 4.3, 1.4, 1.0 } ;
       double nb3_m3_yerrh_2s[4] = { 5.8, 10.2, 4.0, 3.2 } ;
       double nb3_m3_yerrl_2s[4] = { 4.2, 8.0, 2.6, 1.7 } ;

       TGraphAsymmErrors* gr_nb3_m3 = new TGraphAsymmErrors( 4, xvals, nb3_m3_yvals, xel, xeh,
            nb3_m3_yerrl, nb3_m3_yerrh ) ;

       TGraphAsymmErrors* gr_nb3_m3_2s = new TGraphAsymmErrors( 4, xvals, nb3_m3_yvals, xel, xeh,
            nb3_m3_yerrl_2s, nb3_m3_yerrh_2s ) ;

       gr_nb3_m3 -> SetFillStyle(fillStyle1sig) ;
       gr_nb3_m3 -> SetFillColor(fillColor1sig) ;
       gr_nb3_m3 -> SetLineColor(kGreen+2) ;
       gr_nb3_m3 -> SetLineWidth(1) ;
       gr_nb3_m3 -> SetMarkerColor(4) ;
       //gr_nb3_m3 -> SetMarkerStyle(24) ;

       gr_nb3_m3_2s -> SetFillStyle(fillStyle2sig) ;
       gr_nb3_m3_2s -> SetFillColor(fillColor2sig) ;
       gr_nb3_m3_2s -> SetLineColor(kGreen+2) ;
       gr_nb3_m3_2s -> SetMarkerColor(4) ;
       //gr_nb3_m3_2s -> SetMarkerStyle(24) ;

       gr_nb3_m3_2s -> Draw("p2") ;
       gr_nb3_m3_2s -> Draw("p") ;
       gr_nb3_m3 -> Draw("p2") ;
       gr_nb3_m3 -> Draw("p") ;
       h_data_nb3_m3 -> Draw("same e") ;




    //--------

       cffp->cd(4) ;

       TH1F* h_data_nb3_m4 = new TH1F( "h_data_nb3_m4", "Data, nB=3, MET4", 4, 0.5, 4.5 ) ;

       h_data_nb3_m4 -> GetXaxis() -> SetBinLabel( 1, "HT1" ) ;
       h_data_nb3_m4 -> GetXaxis() -> SetBinLabel( 2, "HT2" ) ;
       h_data_nb3_m4 -> GetXaxis() -> SetBinLabel( 3, "HT3" ) ;
       h_data_nb3_m4 -> GetXaxis() -> SetBinLabel( 4, "HT4" ) ;

       h_data_nb3_m4 -> SetBinContent( 1, 0 ) ;
       h_data_nb3_m4 -> SetBinContent( 2, 8 ) ;
       h_data_nb3_m4 -> SetBinContent( 3, 2 ) ;
       h_data_nb3_m4 -> SetBinContent( 4, 4 ) ;

       h_data_nb3_m4 -> SetMarkerStyle(20) ;
       h_data_nb3_m4 -> SetMarkerSize(1.2) ;
       h_data_nb3_m4 -> SetLineWidth(2) ;
       if ( doLogy ) h_data_nb3_m4 -> SetMinimum( 0.1 ) ;
       h_data_nb3_m4 -> SetMaximum( 16 ) ;


       h_data_nb3_m4 -> Draw("e") ;


       double nb3_m4_yvals[4] = { -99., 8.4, 1.8, 1.7 } ;
       double nb3_m4_yerrh[4] = { 0.0, 2.1, 1.0, 1.0 } ;
       double nb3_m4_yerrl[4] = { 0.0, 1.8, 0.7, 0.7 } ;
       double nb3_m4_yerrh_2s[4] = { 0.0, 4.6, 2.3, 2.3 } ;
       double nb3_m4_yerrl_2s[4] = { 0.0, 3.2, 1.1, 1.2 } ;

       TGraphAsymmErrors* gr_nb3_m4 = new TGraphAsymmErrors( 4, xvals, nb3_m4_yvals, xel, xeh,
            nb3_m4_yerrl, nb3_m4_yerrh ) ;

       TGraphAsymmErrors* gr_nb3_m4_2s = new TGraphAsymmErrors( 4, xvals, nb3_m4_yvals, xel, xeh,
            nb3_m4_yerrl_2s, nb3_m4_yerrh_2s ) ;

       gr_nb3_m4 -> SetFillStyle(fillStyle1sig) ;
       gr_nb3_m4 -> SetFillColor(fillColor1sig) ;
       gr_nb3_m4 -> SetLineColor(kGreen+2) ;
       gr_nb3_m4 -> SetLineWidth(1) ;
       gr_nb3_m4 -> SetMarkerColor(4) ;
       //gr_nb3_m4 -> SetMarkerStyle(24) ;

       gr_nb3_m4_2s -> SetFillStyle(fillStyle2sig) ;
       gr_nb3_m4_2s -> SetFillColor(fillColor2sig) ;
       gr_nb3_m4_2s -> SetLineColor(kGreen+2) ;
       gr_nb3_m4_2s -> SetMarkerColor(4) ;
       //gr_nb3_m4_2s -> SetMarkerStyle(24) ;

       gr_nb3_m4_2s -> Draw("p2") ;
       gr_nb3_m4_2s -> Draw("p") ;
       gr_nb3_m4 -> Draw("p2") ;
       gr_nb3_m4 -> Draw("p") ;
       h_data_nb3_m4 -> Draw("same e") ;

       cffp->Update() ;
       cffp->Draw() ;

       cffp->SaveAs("outputfiles/fitresult-fullfit-hsbins-18fb.pdf" ) ;


   }
