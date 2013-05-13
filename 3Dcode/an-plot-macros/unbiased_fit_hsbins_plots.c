#include "TROOT.h"
#include "TLatex.h"

#include "TPad.h"
#include "TLegend.h"
#include "TFrame.h"

#include "TStyle.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TDirectory.h"

#include "fill_arrays.c"
bool usePublicStyle_=true;
#include "../defineCmsStyle.c"


void formatClone(TGraphAsymmErrors * gg) {
  //what i want is no y err
  for (int kk=0; kk< gg->GetN();kk++) {
    gg->SetPointEYlow(kk,0);
    gg->SetPointEYhigh(kk,0);
  }
  gg->SetMarkerSize(0);
  gg->SetLineWidth(2);
}

   void unbiased_fit_hsbins_plots( bool doLogy=false ) {

     TLatex* nblabel=new TLatex();
     float header_x=0.37; float header_y=0.96;
     if (usePublicStyle_) {
       initCmsStyle();
       gROOT->SetStyle("CMS");
       gROOT->ForceStyle();
       gStyle->SetPadTopMargin(0.06) ; //was 0.02

       nblabel->SetTextSize(0.095);
       nblabel->SetTextFont(42);
       nblabel->SetNDC();
     }

       fill_arrays() ;

       int fillStyle1sig = 1001 ;
       int fillStyle2sig = 1001 ;
       int fillColor1sig = kBlue-9 ;
       int fillColor2sig = kBlue-10 ;

       gStyle->SetOptStat(0) ;
       gStyle->SetLabelSize(usePublicStyle_? 0.17 : 0.12,"x") ;
       gStyle->SetLabelSize(usePublicStyle_? 0.12 : 0.08,"y") ;
       gStyle->SetTitleW(0.8) ;
       gStyle->SetTitleH(0.08) ;
       gStyle->SetTitleAlign(13) ;
       gStyle->SetPadLeftMargin(usePublicStyle_ ? 0.35 : 0.14) ;
       gStyle->SetPadRightMargin(0.02) ;
       gStyle->SetLabelOffset(0.02,"y") ;
       if (usePublicStyle_) gStyle->SetTitleSize(0.12,"y");
       if (usePublicStyle_) gStyle->SetTitleOffset(1.5,"y");
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

       if (usePublicStyle_) {
	 for (int k=0;k<4;k++) xvals[k]-=0.05; //cancel out owen's fancy offset
       }

       TCanvas* cmip = (TCanvas*) gDirectory->FindObject("cmip") ;
       if ( cmip == 0x0 ) {
          cmip = new TCanvas("cmip", "MIP", 800, 600) ;
       }
       cmip->Clear() ;

       //make pads by hand
	 /* original coordinates using divide command
	    pad 0.01 0.01 0.24 0.99
	    pad 0.26 0.01 0.49 0.99
	    pad 0.51 0.01 0.74 0.99
	    pad 0.76 0.01 0.99 0.99
	 */
       float ymax=usePublicStyle_ ? 0.9 : 0.99;
       TPad* 	 pad1 = new TPad("pad1","pad1",	0.01,0.01,0.24,ymax);
       TPad*	 pad2 = new TPad("pad2","pad2",	0.26,0.01,0.49,ymax);
       TPad*	 pad3 = new TPad("pad3","pad3",	0.51,0.01,0.74,ymax);
       TPad*	 pad4 = new TPad("pad4","pad4",	0.76,0.01,0.99,ymax);
       TPad* padh = usePublicStyle_ ? new TPad("padh","padh",0.01,ymax+0.02,0.99,0.99) : 0;
       pad1->Draw(); pad2->Draw(); pad3->Draw(); pad4->Draw();
       if (usePublicStyle_) padh->Draw();

       TLatex* plotheader=0;
      if (usePublicStyle_) {  //copied from drawReducedTrees.h
	padh->cd();
	TString astring;
	//astring.Form("CMS Preliminary, L_{int} = %.1f fb^{-1}, #sqrt{s} = 8 TeV",19.39);
	astring.Form("CMS, L_{int} = %.1f fb^{-1}, #sqrt{s} = 8 TeV",19.39);
	plotheader = new TLatex(3.570061,23.08044,astring);
	plotheader->SetNDC();
	//	plotheader->SetTextAlign(13);
	plotheader->SetX(0.05);
	plotheader->SetY(0.4);
	plotheader->SetTextFont(42);
	plotheader->SetTextSizePixels(24);
	plotheader->SetTextSize(0.6);
	plotheader->Draw();


      }

    //--------

      //      double xlow,ylow,xup,yup;
      //      cmip->cd(1)->GetPadPar(xlow,ylow,xup,yup) ;
      //      cout<<"pad "<<xlow<<" "<<ylow<<" "<<xup<<" "<<yup<<endl;
      pad1->cd();

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

       if (usePublicStyle_) {
	h_data_nb2_m4->SetYTitle("Events / bin");
	h_data_nb2_m4->GetXaxis()->LabelsOption("v");
       }


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

       TString drawopt = usePublicStyle_? "2" : "p2";
       gr_nb2_m4_2s -> Draw(drawopt) ;
       if (!usePublicStyle_)  gr_nb2_m4_2s -> Draw("p") ;
       gr_nb2_m4 -> Draw(drawopt) ;

       TGraphAsymmErrors* gr_nb2_m4_noyerr=(TGraphAsymmErrors*) gr_nb2_m4 ->Clone("gr_nb2_m4_noyerr");
       if (usePublicStyle_ ) {
	 formatClone(gr_nb2_m4_noyerr);
	 gr_nb2_m4_noyerr->Draw("pz");
       }
       else gr_nb2_m4 -> Draw("p") ;
       h_data_nb2_m4 -> Draw("same e") ;

    //--------
       pad2->cd();
       //      cmip->cd(2)->GetPadPar(xlow,ylow,xup,yup) ;
       //      cout<<"pad "<<xlow<<" "<<ylow<<" "<<xup<<" "<<yup<<endl;
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


       if (usePublicStyle_) {
	h_data_nb3_m2->SetYTitle("Events / bin");
	h_data_nb3_m2->GetXaxis()->LabelsOption("v");
       }

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

       gr_nb3_m2_2s -> Draw(drawopt) ;
       if (!usePublicStyle_)       gr_nb3_m2_2s -> Draw("p") ;
       gr_nb3_m2 -> Draw(drawopt) ;

       TGraphAsymmErrors* gr_nb3_m2_noyerr=(TGraphAsymmErrors*) gr_nb3_m2 ->Clone("gr_nb3_m2_noyerr");
       if (usePublicStyle_ ) {
	 formatClone(gr_nb3_m2_noyerr);
	 gr_nb3_m2_noyerr->Draw("pz");
       }
       else  gr_nb3_m2 -> Draw("p") ;
       h_data_nb3_m2 -> Draw("same e") ;



    //--------

       pad3->cd();
       //      cmip->cd(3)->GetPadPar(xlow,ylow,xup,yup) ;
       //      cout<<"pad "<<xlow<<" "<<ylow<<" "<<xup<<" "<<yup<<endl;
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


       if (usePublicStyle_) {
	h_data_nb3_m3->SetYTitle("Events / bin");
	h_data_nb3_m3->GetXaxis()->LabelsOption("v");
       }
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

       gr_nb3_m3_2s -> Draw(drawopt) ;
       if (!usePublicStyle_) gr_nb3_m3_2s -> Draw("p") ;
       gr_nb3_m3 -> Draw(drawopt) ;

       TGraphAsymmErrors* gr_nb3_m3_noyerr=(TGraphAsymmErrors*) gr_nb3_m3 ->Clone("gr_nb3_m3_noyerr");
       if (usePublicStyle_ ) {
	 formatClone(gr_nb3_m3_noyerr);
	 gr_nb3_m3_noyerr->Draw("pz");
       }
       else  gr_nb3_m3 -> Draw("p") ;
       h_data_nb3_m3 -> Draw("same e") ;




    //--------

       pad4->cd();
       //      cmip->cd(4)->GetPadPar(xlow,ylow,xup,yup) ;
       //      cout<<"pad "<<xlow<<" "<<ylow<<" "<<xup<<" "<<yup<<endl;
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


 
       if (usePublicStyle_) {
	h_data_nb3_m4->SetYTitle("Events / bin");
	h_data_nb3_m4->GetXaxis()->LabelsOption("v");
       }
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

       gr_nb3_m4_2s -> Draw(drawopt) ;
       if (!usePublicStyle_) gr_nb3_m4_2s -> Draw("p") ;
       gr_nb3_m4 -> Draw(drawopt) ;
       TGraphAsymmErrors* gr_nb3_m4_noyerr=(TGraphAsymmErrors*) gr_nb3_m4 ->Clone("gr_nb3_m4_noyerr");
       if (usePublicStyle_ ) {
	 formatClone(gr_nb3_m4_noyerr);
	 gr_nb3_m4_noyerr->Draw("pz");
       }
       else gr_nb3_m4 -> Draw("p") ;
       h_data_nb3_m4 -> Draw("same e") ;


       TLegend * leg1=0;
       if (usePublicStyle_) {
	 pad1->cd();
	 nblabel->DrawLatex(header_x,header_y,"N_{b-jet} = 2, MET4");
	 pad2->cd();
	 nblabel->DrawLatex(header_x,header_y,"N_{b-jet} #geq 3, MET2");
	 pad3->cd();
	 nblabel->DrawLatex(header_x,header_y,"N_{b-jet} #geq 3, MET3");
	 pad4->cd();
	 nblabel->DrawLatex(header_x,header_y,"N_{b-jet} #geq 3, MET4");

	 pad1->GetFrame()->Draw();
	 pad1->RedrawAxis();
	 pad2->GetFrame()->Draw();
	 pad2->RedrawAxis();
	 pad3->GetFrame()->Draw();
	 pad3->RedrawAxis();
	 pad4->GetFrame()->Draw();
	 pad4->RedrawAxis();

	 padh->cd();
	 float leg_x1=0.67,leg_y1=0.04,leg_x2=0.98,leg_y2=0.92;
	 leg1 = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2);
	 leg1->SetBorderSize(0);
	 leg1->SetLineStyle(0);
	 leg1->SetTextFont(42);
	 leg1->SetFillStyle(0);
	 leg1->SetNColumns(2);
	 leg1->AddEntry(gr_nb2_m4_noyerr,"Sideband fit","lf");
	 leg1->AddEntry(h_data_nb2_m4,"Data");
	 leg1->Draw();



       }

       cmip->Update() ;
       cmip->Draw() ;


       if ( doLogy ) {
          cmip->SaveAs("fitresult-unbiased-hsbins-2sig-19fb-logy.pdf" ) ;
       } else {
          cmip->SaveAs("fitresult-unbiased-hsbins-2sig-19fb.eps" ) ;
          cmip->SaveAs("fitresult-unbiased-hsbins-2sig-19fb.png" ) ;
          cmip->SaveAs("fitresult-unbiased-hsbins-2sig-19fb.pdf" ) ;
       }


   }
