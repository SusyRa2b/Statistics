
#include "TStyle.h"
#include "TLegend.h"
#include "TLine.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TDirectory.h"
#include "shifthist.c"

#include <iostream>

 using std::cout ;



   void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;

//----------------


    void calc_wjets_stop_shapesyst( const char* infile = "rootfiles/gi-plots-met4-ht4-v15.root" ) {

       gStyle->SetOptStat(0) ;
       gStyle->SetPadGridX(0) ;
       gStyle->SetPadGridY(1) ;
       gStyle->SetPadBottomMargin(0.18) ;


       loadHist( infile ) ;

     //--- inputs

       TH1F* hmctruth_ttwj_0lep_1b = (TH1F*) gDirectory->FindObject( "hmctruth_ttwj_0lep_1b" ) ;
       TH1F* hmctruth_ttwj_0lep_2b = (TH1F*) gDirectory->FindObject( "hmctruth_ttwj_0lep_2b" ) ;
       TH1F* hmctruth_ttwj_0lep_3b = (TH1F*) gDirectory->FindObject( "hmctruth_ttwj_0lep_3b" ) ;
       TH1F* hmctruth_ttwj_1lep_1b = (TH1F*) gDirectory->FindObject( "hmctruth_ttwj_1lep_1b" ) ;
       TH1F* hmctruth_ttwj_1lep_2b = (TH1F*) gDirectory->FindObject( "hmctruth_ttwj_1lep_2b" ) ;
       TH1F* hmctruth_ttwj_1lep_3b = (TH1F*) gDirectory->FindObject( "hmctruth_ttwj_1lep_3b" ) ;

       TH1F* hmctruth_ttbar_0lep_1b = (TH1F*) gDirectory->FindObject( "hmctruth_ttbar_0lep_1b" ) ;
       TH1F* hmctruth_ttbar_0lep_2b = (TH1F*) gDirectory->FindObject( "hmctruth_ttbar_0lep_2b" ) ;
       TH1F* hmctruth_ttbar_0lep_3b = (TH1F*) gDirectory->FindObject( "hmctruth_ttbar_0lep_3b" ) ;
       TH1F* hmctruth_ttbar_1lep_1b = (TH1F*) gDirectory->FindObject( "hmctruth_ttbar_1lep_1b" ) ;
       TH1F* hmctruth_ttbar_1lep_2b = (TH1F*) gDirectory->FindObject( "hmctruth_ttbar_1lep_2b" ) ;
       TH1F* hmctruth_ttbar_1lep_3b = (TH1F*) gDirectory->FindObject( "hmctruth_ttbar_1lep_3b" ) ;

       TH1F* hmctruth_wjetsonly_0lep_1b = (TH1F*) gDirectory->FindObject( "hmctruth_wjetsonly_0lep_1b" ) ;
       TH1F* hmctruth_wjetsonly_0lep_2b = (TH1F*) gDirectory->FindObject( "hmctruth_wjetsonly_0lep_2b" ) ;
       TH1F* hmctruth_wjetsonly_0lep_3b = (TH1F*) gDirectory->FindObject( "hmctruth_wjetsonly_0lep_3b" ) ;
       TH1F* hmctruth_wjetsonly_1lep_1b = (TH1F*) gDirectory->FindObject( "hmctruth_wjetsonly_1lep_1b" ) ;
       TH1F* hmctruth_wjetsonly_1lep_2b = (TH1F*) gDirectory->FindObject( "hmctruth_wjetsonly_1lep_2b" ) ;
       TH1F* hmctruth_wjetsonly_1lep_3b = (TH1F*) gDirectory->FindObject( "hmctruth_wjetsonly_1lep_3b" ) ;

       TH1F* hmctruth_singletop_0lep_1b = (TH1F*) gDirectory->FindObject( "hmctruth_singletop_0lep_1b" ) ;
       TH1F* hmctruth_singletop_0lep_2b = (TH1F*) gDirectory->FindObject( "hmctruth_singletop_0lep_2b" ) ;
       TH1F* hmctruth_singletop_0lep_3b = (TH1F*) gDirectory->FindObject( "hmctruth_singletop_0lep_3b" ) ;
       TH1F* hmctruth_singletop_1lep_1b = (TH1F*) gDirectory->FindObject( "hmctruth_singletop_1lep_1b" ) ;
       TH1F* hmctruth_singletop_1lep_2b = (TH1F*) gDirectory->FindObject( "hmctruth_singletop_1lep_2b" ) ;
       TH1F* hmctruth_singletop_1lep_3b = (TH1F*) gDirectory->FindObject( "hmctruth_singletop_1lep_3b" ) ;

      //--- zero out the ignored max MET, min HT bin (MET4,HT1)

       hmctruth_ttwj_0lep_1b->SetBinContent( 17, 0. ) ;
       hmctruth_ttwj_0lep_2b->SetBinContent( 17, 0. ) ;
       hmctruth_ttwj_0lep_3b->SetBinContent( 17, 0. ) ;
       hmctruth_ttwj_1lep_1b->SetBinContent( 17, 0. ) ;
       hmctruth_ttwj_1lep_2b->SetBinContent( 17, 0. ) ;
       hmctruth_ttwj_1lep_3b->SetBinContent( 17, 0. ) ;

       hmctruth_ttbar_0lep_1b->SetBinContent( 17, 0. ) ;
       hmctruth_ttbar_0lep_2b->SetBinContent( 17, 0. ) ;
       hmctruth_ttbar_0lep_3b->SetBinContent( 17, 0. ) ;
       hmctruth_ttbar_1lep_1b->SetBinContent( 17, 0. ) ;
       hmctruth_ttbar_1lep_2b->SetBinContent( 17, 0. ) ;
       hmctruth_ttbar_1lep_3b->SetBinContent( 17, 0. ) ;

       hmctruth_wjetsonly_0lep_1b->SetBinContent( 17, 0. ) ;
       hmctruth_wjetsonly_0lep_2b->SetBinContent( 17, 0. ) ;
       hmctruth_wjetsonly_0lep_3b->SetBinContent( 17, 0. ) ;
       hmctruth_wjetsonly_1lep_1b->SetBinContent( 17, 0. ) ;
       hmctruth_wjetsonly_1lep_2b->SetBinContent( 17, 0. ) ;
       hmctruth_wjetsonly_1lep_3b->SetBinContent( 17, 0. ) ;

       hmctruth_singletop_0lep_1b->SetBinContent( 17, 0. ) ;
       hmctruth_singletop_0lep_2b->SetBinContent( 17, 0. ) ;
       hmctruth_singletop_0lep_3b->SetBinContent( 17, 0. ) ;
       hmctruth_singletop_1lep_1b->SetBinContent( 17, 0. ) ;
       hmctruth_singletop_1lep_2b->SetBinContent( 17, 0. ) ;
       hmctruth_singletop_1lep_3b->SetBinContent( 17, 0. ) ;










     //--- fraction plots

       TH1F* h_frac_wjets_0lep_1b  = (TH1F*) hmctruth_wjetsonly_0lep_1b->Clone("h_frac_wjets_0lep_1b") ;
       TH1F* h_frac_wjets_0lep_2b  = (TH1F*) hmctruth_wjetsonly_0lep_2b->Clone("h_frac_wjets_0lep_2b") ;
       TH1F* h_frac_wjets_0lep_3b  = (TH1F*) hmctruth_wjetsonly_0lep_3b->Clone("h_frac_wjets_0lep_3b") ;

       h_frac_wjets_0lep_1b->Divide( hmctruth_ttwj_0lep_1b ) ;
       h_frac_wjets_0lep_2b->Divide( hmctruth_ttwj_0lep_2b ) ;
       h_frac_wjets_0lep_3b->Divide( hmctruth_ttwj_0lep_3b ) ;

       TH1F* h_frac_wjets_1lep_1b  = (TH1F*) hmctruth_wjetsonly_1lep_1b->Clone("h_frac_wjets_1lep_1b") ;
       TH1F* h_frac_wjets_1lep_2b  = (TH1F*) hmctruth_wjetsonly_1lep_2b->Clone("h_frac_wjets_1lep_2b") ;
       TH1F* h_frac_wjets_1lep_3b  = (TH1F*) hmctruth_wjetsonly_1lep_3b->Clone("h_frac_wjets_1lep_3b") ;

       h_frac_wjets_1lep_1b->Divide( hmctruth_ttwj_1lep_1b ) ;
       h_frac_wjets_1lep_2b->Divide( hmctruth_ttwj_1lep_2b ) ;
       h_frac_wjets_1lep_3b->Divide( hmctruth_ttwj_1lep_3b ) ;


     //---

       TH1F* h_frac_singletop_0lep_1b_ns  = (TH1F*) hmctruth_singletop_0lep_1b->Clone("h_frac_singletop_0lep_1b") ;
       TH1F* h_frac_singletop_0lep_2b_ns  = (TH1F*) hmctruth_singletop_0lep_2b->Clone("h_frac_singletop_0lep_2b") ;
       TH1F* h_frac_singletop_0lep_3b_ns  = (TH1F*) hmctruth_singletop_0lep_3b->Clone("h_frac_singletop_0lep_3b") ;

       h_frac_singletop_0lep_1b_ns->Divide( hmctruth_ttwj_0lep_1b ) ;
       h_frac_singletop_0lep_2b_ns->Divide( hmctruth_ttwj_0lep_2b ) ;
       h_frac_singletop_0lep_3b_ns->Divide( hmctruth_ttwj_0lep_3b ) ;

       TH1F* h_frac_singletop_1lep_1b_ns  = (TH1F*) hmctruth_singletop_1lep_1b->Clone("h_frac_singletop_1lep_1b") ;
       TH1F* h_frac_singletop_1lep_2b_ns  = (TH1F*) hmctruth_singletop_1lep_2b->Clone("h_frac_singletop_1lep_2b") ;
       TH1F* h_frac_singletop_1lep_3b_ns  = (TH1F*) hmctruth_singletop_1lep_3b->Clone("h_frac_singletop_1lep_3b") ;

       h_frac_singletop_1lep_1b_ns->Divide( hmctruth_ttwj_1lep_1b ) ;
       h_frac_singletop_1lep_2b_ns->Divide( hmctruth_ttwj_1lep_2b ) ;
       h_frac_singletop_1lep_3b_ns->Divide( hmctruth_ttwj_1lep_3b ) ;

       TH1F* h_frac_singletop_0lep_1b = shifthist( h_frac_singletop_0lep_1b_ns, 0.15 ) ;
       TH1F* h_frac_singletop_0lep_2b = shifthist( h_frac_singletop_0lep_2b_ns, 0.15 ) ;
       TH1F* h_frac_singletop_0lep_3b = shifthist( h_frac_singletop_0lep_3b_ns, 0.15 ) ;
       TH1F* h_frac_singletop_1lep_1b = shifthist( h_frac_singletop_1lep_1b_ns, 0.15 ) ;
       TH1F* h_frac_singletop_1lep_2b = shifthist( h_frac_singletop_1lep_2b_ns, 0.15 ) ;
       TH1F* h_frac_singletop_1lep_3b = shifthist( h_frac_singletop_1lep_3b_ns, 0.15 ) ;


       double fracmax(0.45) ;
       double fracmin(-0.05) ;

       h_frac_wjets_0lep_1b->SetMinimum( fracmin ) ;
       h_frac_wjets_0lep_2b->SetMinimum( fracmin ) ;
       h_frac_wjets_0lep_3b->SetMinimum( fracmin ) ;
       h_frac_wjets_1lep_1b->SetMinimum( fracmin ) ;
       h_frac_wjets_1lep_2b->SetMinimum( fracmin ) ;
       h_frac_wjets_1lep_3b->SetMinimum( fracmin ) ;

       h_frac_wjets_0lep_1b->SetMaximum( fracmax ) ;
       h_frac_wjets_0lep_2b->SetMaximum( fracmax ) ;
       h_frac_wjets_0lep_3b->SetMaximum( fracmax ) ;
       h_frac_wjets_1lep_1b->SetMaximum( fracmax ) ;
       h_frac_wjets_1lep_2b->SetMaximum( fracmax ) ;
       h_frac_wjets_1lep_3b->SetMaximum( fracmax ) ;


       h_frac_wjets_0lep_1b->SetMarkerStyle(20) ;
       h_frac_wjets_0lep_2b->SetMarkerStyle(20) ;
       h_frac_wjets_0lep_3b->SetMarkerStyle(20) ;
       h_frac_wjets_1lep_1b->SetMarkerStyle(20) ;
       h_frac_wjets_1lep_2b->SetMarkerStyle(20) ;
       h_frac_wjets_1lep_3b->SetMarkerStyle(20) ;


       h_frac_wjets_0lep_1b->SetLineColor(4) ;
       h_frac_wjets_0lep_2b->SetLineColor(4) ;
       h_frac_wjets_0lep_3b->SetLineColor(4) ;
       h_frac_wjets_1lep_1b->SetLineColor(4) ;
       h_frac_wjets_1lep_2b->SetLineColor(4) ;
       h_frac_wjets_1lep_3b->SetLineColor(4) ;

       h_frac_wjets_0lep_1b->SetMarkerColor(4) ;
       h_frac_wjets_0lep_2b->SetMarkerColor(4) ;
       h_frac_wjets_0lep_3b->SetMarkerColor(4) ;
       h_frac_wjets_1lep_1b->SetMarkerColor(4) ;
       h_frac_wjets_1lep_2b->SetMarkerColor(4) ;
       h_frac_wjets_1lep_3b->SetMarkerColor(4) ;


       h_frac_singletop_0lep_1b->SetMarkerStyle(22) ;
       h_frac_singletop_0lep_2b->SetMarkerStyle(22) ;
       h_frac_singletop_0lep_3b->SetMarkerStyle(22) ;
       h_frac_singletop_1lep_1b->SetMarkerStyle(22) ;
       h_frac_singletop_1lep_2b->SetMarkerStyle(22) ;
       h_frac_singletop_1lep_3b->SetMarkerStyle(22) ;


       h_frac_singletop_0lep_1b->SetLineColor(2) ;
       h_frac_singletop_0lep_2b->SetLineColor(2) ;
       h_frac_singletop_0lep_3b->SetLineColor(2) ;
       h_frac_singletop_1lep_1b->SetLineColor(2) ;
       h_frac_singletop_1lep_2b->SetLineColor(2) ;
       h_frac_singletop_1lep_3b->SetLineColor(2) ;

       h_frac_singletop_0lep_1b->SetMarkerColor(2) ;
       h_frac_singletop_0lep_2b->SetMarkerColor(2) ;
       h_frac_singletop_0lep_3b->SetMarkerColor(2) ;
       h_frac_singletop_1lep_1b->SetMarkerColor(2) ;
       h_frac_singletop_1lep_2b->SetMarkerColor(2) ;
       h_frac_singletop_1lep_3b->SetMarkerColor(2) ;


       TLine* line = new TLine() ;
       line->SetLineColor(1) ;


       TLegend* legend = new TLegend( 0.15, 0.74, 0.68, 0.88 ) ;
       legend->SetFillColor(kWhite) ;
       legend->AddEntry( h_frac_wjets_0lep_1b, "W+jets fraction of ttwj" ) ;
       legend->AddEntry( h_frac_singletop_0lep_1b, "single t fraction of ttwj" ) ;


       TCanvas* can1 = new TCanvas( "can1", "Wjets fraction of ttwj", 1000, 700 ) ;

       TH1F* hist1(0x0) ;
       TH1F* hist2(0x0) ;

       can1->Divide(3,2) ;



       can1->cd(1) ;
       hist1 = h_frac_wjets_0lep_1b ;
       hist2 = h_frac_singletop_0lep_1b ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       line->DrawLine( hist1 -> GetXaxis() -> GetBinLowEdge( 1 ), 0, hist1 -> GetXaxis() -> GetBinUpEdge( hist1  -> GetNbinsX() ), 0 ) ;
       legend->Draw() ;


       can1->cd(2) ;
       hist1 = h_frac_wjets_0lep_2b ;
       hist2 = h_frac_singletop_0lep_2b ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       line->DrawLine( hist1 -> GetXaxis() -> GetBinLowEdge( 1 ), 0, hist1 -> GetXaxis() -> GetBinUpEdge( hist1  -> GetNbinsX() ), 0 ) ;
       legend->Draw() ;


       can1->cd(3) ;
       hist1 = h_frac_wjets_0lep_3b ;
       hist2 = h_frac_singletop_0lep_3b ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       line->DrawLine( hist1 -> GetXaxis() -> GetBinLowEdge( 1 ), 0, hist1 -> GetXaxis() -> GetBinUpEdge( hist1  -> GetNbinsX() ), 0 ) ;
       legend->Draw() ;



       can1->cd(4) ;
       hist1 = h_frac_wjets_1lep_1b ;
       hist2 = h_frac_singletop_1lep_1b ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       line->DrawLine( hist1 -> GetXaxis() -> GetBinLowEdge( 1 ), 0, hist1 -> GetXaxis() -> GetBinUpEdge( hist1  -> GetNbinsX() ), 0 ) ;
       legend->Draw() ;


       can1->cd(5) ;
       hist1 = h_frac_wjets_1lep_2b ;
       hist2 = h_frac_singletop_1lep_2b ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       line->DrawLine( hist1 -> GetXaxis() -> GetBinLowEdge( 1 ), 0, hist1 -> GetXaxis() -> GetBinUpEdge( hist1  -> GetNbinsX() ), 0 ) ;
       legend->Draw() ;


       can1->cd(6) ;
       hist1 = h_frac_wjets_1lep_3b ;
       hist2 = h_frac_singletop_1lep_3b ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       line->DrawLine( hist1 -> GetXaxis() -> GetBinLowEdge( 1 ), 0, hist1 -> GetXaxis() -> GetBinUpEdge( hist1  -> GetNbinsX() ), 0 ) ;
       legend->Draw() ;





     //--- Event count plots.

       TH1F* hist3(0x0) ;


       legend = new TLegend(0.63, 0.70, 0.88, 0.88 ) ;
       legend->SetFillColor(kWhite) ;
       legend->AddEntry( hmctruth_ttbar_0lep_1b, "ttbar" ) ;
       legend->AddEntry( hmctruth_singletop_0lep_1b, "single t" ) ;
       legend->AddEntry( hmctruth_wjetsonly_0lep_1b, "Wjets" ) ;


       TCanvas* can2 = new TCanvas("can2","Wjets event counts", 1000, 700 ) ;

       THStack* stack1(0x0) ;
       THStack* stack2(0x0) ;

       double maxfac(1.5) ;

       can2->Divide(3,2) ;


       can2->cd(1) ;
       hist1 = hmctruth_wjetsonly_0lep_1b ;
       hist2 = hmctruth_singletop_0lep_1b ;
       hist3 = hmctruth_ttbar_0lep_1b ;
       stack1 = new THStack( "hs1_0lep_1b", "0lep, 1b" ) ;
       stack2 = new THStack( "hs2_0lep_1b", "0lep, 1b" ) ;
       hist1 -> SetFillColor( 600-4 ) ;
       hist3 -> SetFillColor( 600-9 ) ;
       hist2 -> SetFillColor( 600-7 ) ;
       stack1->Add( hist1 ) ;
       stack1->Add( hist2 ) ;
       stack2->Add( hist1 ) ;
       stack2->Add( hist2 ) ;
       stack2->Add( hist3 ) ;
       hist1->SetMaximum( maxfac*( stack1->GetMaximum()) ) ;
       hist1->Draw() ;
       stack2->Draw("hist same") ;
       stack2->Draw("same") ;
       legend->Draw() ;


       can2->cd(2) ;
       hist1 = hmctruth_wjetsonly_0lep_2b ;
       hist2 = hmctruth_singletop_0lep_2b ;
       hist3 = hmctruth_ttbar_0lep_2b ;
       stack1 = new THStack( "hs1_0lep_2b", "0lep, 2b" ) ;
       stack2 = new THStack( "hs2_0lep_2b", "0lep, 2b" ) ;
       hist1 -> SetFillColor( 600-4 ) ;
       hist3 -> SetFillColor( 600-9 ) ;
       hist2 -> SetFillColor( 600-7 ) ;
       stack1->Add( hist1 ) ;
       stack1->Add( hist2 ) ;
       stack2->Add( hist1 ) ;
       stack2->Add( hist2 ) ;
       stack2->Add( hist3 ) ;
       hist1->SetMaximum( maxfac*( stack1->GetMaximum()) ) ;
       hist1->Draw() ;
       stack2->Draw("hist same") ;
       stack2->Draw("same") ;


       can2->cd(3) ;
       hist1 = hmctruth_wjetsonly_0lep_3b ;
       hist2 = hmctruth_singletop_0lep_3b ;
       hist3 = hmctruth_ttbar_0lep_3b ;
       stack1 = new THStack( "hs1_0lep_3b", "0lep, 3b" ) ;
       stack2 = new THStack( "hs2_0lep_3b", "0lep, 3b" ) ;
       hist1 -> SetFillColor( 600-4 ) ;
       hist3 -> SetFillColor( 600-9 ) ;
       hist2 -> SetFillColor( 600-7 ) ;
       stack1->Add( hist1 ) ;
       stack1->Add( hist2 ) ;
       stack2->Add( hist1 ) ;
       stack2->Add( hist2 ) ;
       stack2->Add( hist3 ) ;
       hist1->SetMaximum( maxfac*( stack1->GetMaximum()) ) ;
       hist1->Draw() ;
       stack2->Draw("hist same") ;
       stack2->Draw("same") ;





       can2->cd(4) ;
       hist1 = hmctruth_wjetsonly_1lep_1b ;
       hist2 = hmctruth_singletop_1lep_1b ;
       hist3 = hmctruth_ttbar_1lep_1b ;
       stack1 = new THStack( "hs1_1lep_1b", "1lep, 1b" ) ;
       stack2 = new THStack( "hs2_1lep_1b", "1lep, 1b" ) ;
       hist1 -> SetFillColor( 600-4 ) ;
       hist3 -> SetFillColor( 600-9 ) ;
       hist2 -> SetFillColor( 600-7 ) ;
       stack1->Add( hist1 ) ;
       stack1->Add( hist2 ) ;
       stack2->Add( hist1 ) ;
       stack2->Add( hist2 ) ;
       stack2->Add( hist3 ) ;
       hist1->SetMaximum( maxfac*( stack1->GetMaximum()) ) ;
       hist1->Draw() ;
       stack2->Draw("hist same") ;
       stack2->Draw("same") ;
       legend->Draw() ;


       can2->cd(5) ;
       hist1 = hmctruth_wjetsonly_1lep_2b ;
       hist2 = hmctruth_singletop_1lep_2b ;
       hist3 = hmctruth_ttbar_1lep_2b ;
       stack1 = new THStack( "hs1_1lep_2b", "1lep, 2b" ) ;
       stack2 = new THStack( "hs2_1lep_2b", "1lep, 2b" ) ;
       hist1 -> SetFillColor( 600-4 ) ;
       hist3 -> SetFillColor( 600-9 ) ;
       hist2 -> SetFillColor( 600-7 ) ;
       stack1->Add( hist1 ) ;
       stack1->Add( hist2 ) ;
       stack2->Add( hist1 ) ;
       stack2->Add( hist2 ) ;
       stack2->Add( hist3 ) ;
       hist1->SetMaximum( maxfac*( stack1->GetMaximum()) ) ;
       hist1->Draw() ;
       stack2->Draw("hist same") ;
       stack2->Draw("same") ;


       can2->cd(6) ;
       hist1 = hmctruth_wjetsonly_1lep_3b ;
       hist2 = hmctruth_singletop_1lep_3b ;
       hist3 = hmctruth_ttbar_1lep_3b ;
       stack1 = new THStack( "hs1_1lep_3b", "1lep, 3b" ) ;
       stack2 = new THStack( "hs2_1lep_3b", "1lep, 3b" ) ;
       hist1 -> SetFillColor( 600-4 ) ;
       hist3 -> SetFillColor( 600-9 ) ;
       hist2 -> SetFillColor( 600-7 ) ;
       stack1->Add( hist1 ) ;
       stack1->Add( hist2 ) ;
       stack2->Add( hist1 ) ;
       stack2->Add( hist2 ) ;
       stack2->Add( hist3 ) ;
       hist1->SetMaximum( maxfac*( stack1->GetMaximum()) ) ;
       hist1->Draw() ;
       stack2->Draw("hist same") ;
       stack2->Draw("same") ;







    //--- Zl/SL Ratios.

       TH1F* h_zlsl_ratio_ttbar_1b = (TH1F*) hmctruth_ttbar_0lep_1b->Clone( "h_zlsl_ratio_ttbar_1b" ) ;
       TH1F* h_zlsl_ratio_ttbar_2b = (TH1F*) hmctruth_ttbar_0lep_2b->Clone( "h_zlsl_ratio_ttbar_2b" ) ;
       TH1F* h_zlsl_ratio_ttbar_3b = (TH1F*) hmctruth_ttbar_0lep_3b->Clone( "h_zlsl_ratio_ttbar_3b" ) ;

       h_zlsl_ratio_ttbar_1b->Divide( hmctruth_ttbar_1lep_1b ) ;
       h_zlsl_ratio_ttbar_2b->Divide( hmctruth_ttbar_1lep_2b ) ;
       h_zlsl_ratio_ttbar_3b->Divide( hmctruth_ttbar_1lep_3b ) ;


       TH1F* h_zlsl_ratio_singletop_1b = (TH1F*) hmctruth_singletop_0lep_1b->Clone( "h_zlsl_ratio_singletop_1b" ) ;
       TH1F* h_zlsl_ratio_singletop_2b = (TH1F*) hmctruth_singletop_0lep_2b->Clone( "h_zlsl_ratio_singletop_2b" ) ;
       TH1F* h_zlsl_ratio_singletop_3b = (TH1F*) hmctruth_singletop_0lep_3b->Clone( "h_zlsl_ratio_singletop_3b" ) ;

       h_zlsl_ratio_singletop_1b->Divide( hmctruth_singletop_1lep_1b ) ;
       h_zlsl_ratio_singletop_2b->Divide( hmctruth_singletop_1lep_2b ) ;
       h_zlsl_ratio_singletop_3b->Divide( hmctruth_singletop_1lep_3b ) ;

       TH1F* h_zlsl_ratio_singletop_1b_s = shifthist( h_zlsl_ratio_singletop_1b, 0.1 ) ;
       TH1F* h_zlsl_ratio_singletop_2b_s = shifthist( h_zlsl_ratio_singletop_2b, 0.1 ) ;
       TH1F* h_zlsl_ratio_singletop_3b_s = shifthist( h_zlsl_ratio_singletop_3b, 0.1 ) ;


       TH1F* h_zlsl_ratio_wjetsonly_1b = (TH1F*) hmctruth_wjetsonly_0lep_1b->Clone( "h_zlsl_ratio_wjetsonly_1b" ) ;
       TH1F* h_zlsl_ratio_wjetsonly_2b = (TH1F*) hmctruth_wjetsonly_0lep_2b->Clone( "h_zlsl_ratio_wjetsonly_2b" ) ;
       TH1F* h_zlsl_ratio_wjetsonly_3b = (TH1F*) hmctruth_wjetsonly_0lep_3b->Clone( "h_zlsl_ratio_wjetsonly_3b" ) ;

       h_zlsl_ratio_wjetsonly_1b->Divide( hmctruth_wjetsonly_1lep_1b ) ;
       h_zlsl_ratio_wjetsonly_2b->Divide( hmctruth_wjetsonly_1lep_2b ) ;
       h_zlsl_ratio_wjetsonly_3b->Divide( hmctruth_wjetsonly_1lep_3b ) ;

       TH1F* h_zlsl_ratio_wjetsonly_1b_s = shifthist( h_zlsl_ratio_wjetsonly_1b, -0.1 ) ;
       TH1F* h_zlsl_ratio_wjetsonly_2b_s = shifthist( h_zlsl_ratio_wjetsonly_2b, -0.1 ) ;
       TH1F* h_zlsl_ratio_wjetsonly_3b_s = shifthist( h_zlsl_ratio_wjetsonly_3b, -0.1 ) ;


       legend = new TLegend( 0.12, 0.71, 0.63, 0.88 ) ;
       legend->SetFillColor(kWhite) ;
       h_zlsl_ratio_ttbar_1b->SetFillColor(0) ;
       legend->AddEntry( h_zlsl_ratio_ttbar_1b, "ZL/SL ratio, ttbar" ) ;
       legend->AddEntry( h_zlsl_ratio_singletop_1b_s, "ZL/SL ratio, single t" ) ;
       legend->AddEntry( h_zlsl_ratio_wjetsonly_1b_s, "ZL/SL ratio, W+jets" ) ;



       TCanvas* can3 = new TCanvas("can3","ZL/SL ratios", 1000, 350 ) ;

       double hmax(3.5) ;
       double hmin(0.0) ;

       can3->Divide(3,1) ;

       can3->cd(1) ;
       hist1 = h_zlsl_ratio_ttbar_1b ;
       hist2 = h_zlsl_ratio_singletop_1b_s ;
       hist3 = h_zlsl_ratio_wjetsonly_1b_s ;
       hist1->SetTitle("1btag") ;
       hist1->SetLineColor(2) ;
       hist2->SetLineColor(4) ;
       hist3->SetLineColor(1) ;
       hist1->SetMarkerColor(2) ;
       hist2->SetMarkerColor(4) ;
       hist3->SetMarkerColor(1) ;
       hist1->SetMarkerStyle(20) ;
       hist2->SetMarkerStyle(22) ;
       hist3->SetMarkerStyle(32) ;
       hist1->SetMaximum(hmax) ;
       hist1->SetMinimum(hmin) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;
       legend->Draw() ;



       can3->cd(2) ;
       hist1 = h_zlsl_ratio_ttbar_2b ;
       hist2 = h_zlsl_ratio_singletop_2b_s ;
       hist3 = h_zlsl_ratio_wjetsonly_2b_s ;
       hist1->SetTitle("2btag") ;
       hist1->SetLineColor(2) ;
       hist2->SetLineColor(4) ;
       hist3->SetLineColor(1) ;
       hist1->SetMarkerColor(2) ;
       hist2->SetMarkerColor(4) ;
       hist3->SetMarkerColor(1) ;
       hist1->SetMarkerStyle(20) ;
       hist2->SetMarkerStyle(22) ;
       hist3->SetMarkerStyle(32) ;
       hist1->SetMaximum(hmax) ;
       hist1->SetMinimum(hmin) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;



       can3->cd(3) ;
       hist1 = h_zlsl_ratio_ttbar_3b ;
       hist2 = h_zlsl_ratio_singletop_3b_s ;
       hist3 = h_zlsl_ratio_wjetsonly_3b_s ;
       hist1->SetTitle("ge3btag") ;
       hist1->SetLineColor(2) ;
       hist2->SetLineColor(4) ;
       hist3->SetLineColor(1) ;
       hist1->SetMarkerColor(2) ;
       hist2->SetMarkerColor(4) ;
       hist3->SetMarkerColor(1) ;
       hist1->SetMarkerStyle(20) ;
       hist2->SetMarkerStyle(22) ;
       hist3->SetMarkerStyle(32) ;
       hist1->SetMaximum(hmax) ;
       hist1->SetMinimum(hmin) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;









   //--- ZL/SL ratios and scale factors, nominal


       TH1F* h_zl_ttwj_nom_1b = (TH1F*) hmctruth_ttbar_0lep_1b->Clone( "h_zl_ttwj_nom_1b" ) ;
       TH1F* h_zl_ttwj_nom_2b = (TH1F*) hmctruth_ttbar_0lep_2b->Clone( "h_zl_ttwj_nom_2b" ) ;
       TH1F* h_zl_ttwj_nom_3b = (TH1F*) hmctruth_ttbar_0lep_3b->Clone( "h_zl_ttwj_nom_3b" ) ;

       h_zl_ttwj_nom_1b->Add( hmctruth_wjetsonly_0lep_1b ) ;
       h_zl_ttwj_nom_2b->Add( hmctruth_wjetsonly_0lep_2b ) ;
       h_zl_ttwj_nom_3b->Add( hmctruth_wjetsonly_0lep_3b ) ;

       h_zl_ttwj_nom_1b->Add( hmctruth_singletop_0lep_1b ) ;
       h_zl_ttwj_nom_2b->Add( hmctruth_singletop_0lep_2b ) ;
       h_zl_ttwj_nom_3b->Add( hmctruth_singletop_0lep_3b ) ;


       TH1F* h_sl_ttwj_nom_1b = (TH1F*) hmctruth_ttbar_1lep_1b->Clone( "h_sl_ttwj_nom_1b" ) ;
       TH1F* h_sl_ttwj_nom_2b = (TH1F*) hmctruth_ttbar_1lep_2b->Clone( "h_sl_ttwj_nom_2b" ) ;
       TH1F* h_sl_ttwj_nom_3b = (TH1F*) hmctruth_ttbar_1lep_3b->Clone( "h_sl_ttwj_nom_3b" ) ;

       h_sl_ttwj_nom_1b->Add( hmctruth_wjetsonly_1lep_1b ) ;
       h_sl_ttwj_nom_2b->Add( hmctruth_wjetsonly_1lep_2b ) ;
       h_sl_ttwj_nom_3b->Add( hmctruth_wjetsonly_1lep_3b ) ;

       h_sl_ttwj_nom_1b->Add( hmctruth_singletop_1lep_1b ) ;
       h_sl_ttwj_nom_2b->Add( hmctruth_singletop_1lep_2b ) ;
       h_sl_ttwj_nom_3b->Add( hmctruth_singletop_1lep_3b ) ;


       double sum_zl_1b = h_zl_ttwj_nom_1b->Integral() ;
       double sum_zl_2b = h_zl_ttwj_nom_2b->Integral() ;
       double sum_zl_3b = h_zl_ttwj_nom_3b->Integral() ;

       double sum_sl_1b = h_sl_ttwj_nom_1b->Integral() ;
       double sum_sl_2b = h_sl_ttwj_nom_2b->Integral() ;
       double sum_sl_3b = h_sl_ttwj_nom_3b->Integral() ;

       double ave_zlsl_ratio_all = ( sum_zl_1b + sum_zl_2b + sum_zl_3b ) / ( sum_sl_1b + sum_sl_2b + sum_sl_3b ) ;
       double ave_zlsl_ratio_2b = sum_zl_2b / sum_sl_2b ;
       double ave_zlsl_ratio_3b = sum_zl_3b / sum_sl_3b ;

       printf("\n\n") ;
       printf(" Average ZL/SL, all = %6.3f\n", ave_zlsl_ratio_all ) ;
       printf(" Average ZL/SL, 2b  = %6.3f\n", ave_zlsl_ratio_2b ) ;
       printf(" Average ZL/SL, 3b  = %6.3f\n", ave_zlsl_ratio_3b ) ;
       printf("\n\n") ;


       TH1F* h_zlsl_ratio_ttwj_nom_1b = (TH1F*) h_zl_ttwj_nom_1b->Clone( "h_zlsl_ratio_ttwj_nom_1b" ) ;
       TH1F* h_zlsl_ratio_ttwj_nom_2b = (TH1F*) h_zl_ttwj_nom_2b->Clone( "h_zlsl_ratio_ttwj_nom_2b" ) ;
       TH1F* h_zlsl_ratio_ttwj_nom_3b = (TH1F*) h_zl_ttwj_nom_3b->Clone( "h_zlsl_ratio_ttwj_nom_3b" ) ;

       h_zlsl_ratio_ttwj_nom_1b->Divide( h_sl_ttwj_nom_1b ) ;
       h_zlsl_ratio_ttwj_nom_2b->Divide( h_sl_ttwj_nom_2b ) ;
       h_zlsl_ratio_ttwj_nom_3b->Divide( h_sl_ttwj_nom_3b ) ;


       TH1F* h_zlsl_sf_ttwj_nom_1b = (TH1F*) h_zlsl_ratio_ttwj_nom_1b->Clone( "h_zlsl_ratio_ttwj_nom_1b" ) ;
       TH1F* h_zlsl_sf_ttwj_nom_2b = (TH1F*) h_zlsl_ratio_ttwj_nom_2b->Clone( "h_zlsl_ratio_ttwj_nom_2b" ) ;
       TH1F* h_zlsl_sf_ttwj_nom_3b = (TH1F*) h_zlsl_ratio_ttwj_nom_3b->Clone( "h_zlsl_ratio_ttwj_nom_3b" ) ;

       h_zlsl_sf_ttwj_nom_1b->Scale( 1./ave_zlsl_ratio_all ) ;
       h_zlsl_sf_ttwj_nom_2b->Scale( 1./ave_zlsl_ratio_all ) ;
       h_zlsl_sf_ttwj_nom_3b->Scale( 1./ave_zlsl_ratio_all ) ;




   //--- ZL/SL ratios and scale factors, single top +30%.


       TH1F* hmctruth_singletop_0lep_1b_up = (TH1F*) hmctruth_singletop_0lep_1b->Clone( "hmctruth_singletop_0lep_1b_up" ) ;
       TH1F* hmctruth_singletop_0lep_2b_up = (TH1F*) hmctruth_singletop_0lep_2b->Clone( "hmctruth_singletop_0lep_2b_up" ) ;
       TH1F* hmctruth_singletop_0lep_3b_up = (TH1F*) hmctruth_singletop_0lep_3b->Clone( "hmctruth_singletop_0lep_3b_up" ) ;
       TH1F* hmctruth_singletop_1lep_1b_up = (TH1F*) hmctruth_singletop_1lep_1b->Clone( "hmctruth_singletop_1lep_1b_up" ) ;
       TH1F* hmctruth_singletop_1lep_2b_up = (TH1F*) hmctruth_singletop_1lep_2b->Clone( "hmctruth_singletop_1lep_2b_up" ) ;
       TH1F* hmctruth_singletop_1lep_3b_up = (TH1F*) hmctruth_singletop_1lep_3b->Clone( "hmctruth_singletop_1lep_3b_up" ) ;

       hmctruth_singletop_0lep_1b_up->Scale(1.3) ;
       hmctruth_singletop_0lep_2b_up->Scale(1.3) ;
       hmctruth_singletop_0lep_3b_up->Scale(1.3) ;
       hmctruth_singletop_1lep_1b_up->Scale(1.3) ;
       hmctruth_singletop_1lep_2b_up->Scale(1.3) ;
       hmctruth_singletop_1lep_3b_up->Scale(1.3) ;

       TH1F* h_zl_ttwj_stup_1b = (TH1F*) hmctruth_ttbar_0lep_1b->Clone( "h_zl_ttwj_stup_1b" ) ;
       TH1F* h_zl_ttwj_stup_2b = (TH1F*) hmctruth_ttbar_0lep_2b->Clone( "h_zl_ttwj_stup_2b" ) ;
       TH1F* h_zl_ttwj_stup_3b = (TH1F*) hmctruth_ttbar_0lep_3b->Clone( "h_zl_ttwj_stup_3b" ) ;

       h_zl_ttwj_stup_1b->Add( hmctruth_wjetsonly_0lep_1b ) ;
       h_zl_ttwj_stup_2b->Add( hmctruth_wjetsonly_0lep_2b ) ;
       h_zl_ttwj_stup_3b->Add( hmctruth_wjetsonly_0lep_3b ) ;

       h_zl_ttwj_stup_1b->Add( hmctruth_singletop_0lep_1b_up ) ;
       h_zl_ttwj_stup_2b->Add( hmctruth_singletop_0lep_2b_up ) ;
       h_zl_ttwj_stup_3b->Add( hmctruth_singletop_0lep_3b_up ) ;


       TH1F* h_sl_ttwj_stup_1b = (TH1F*) hmctruth_ttbar_1lep_1b->Clone( "h_sl_ttwj_stup_1b" ) ;
       TH1F* h_sl_ttwj_stup_2b = (TH1F*) hmctruth_ttbar_1lep_2b->Clone( "h_sl_ttwj_stup_2b" ) ;
       TH1F* h_sl_ttwj_stup_3b = (TH1F*) hmctruth_ttbar_1lep_3b->Clone( "h_sl_ttwj_stup_3b" ) ;

       h_sl_ttwj_stup_1b->Add( hmctruth_wjetsonly_1lep_1b ) ;
       h_sl_ttwj_stup_2b->Add( hmctruth_wjetsonly_1lep_2b ) ;
       h_sl_ttwj_stup_3b->Add( hmctruth_wjetsonly_1lep_3b ) ;

       h_sl_ttwj_stup_1b->Add( hmctruth_singletop_1lep_1b_up ) ;
       h_sl_ttwj_stup_2b->Add( hmctruth_singletop_1lep_2b_up ) ;
       h_sl_ttwj_stup_3b->Add( hmctruth_singletop_1lep_3b_up ) ;


       TH1F* h_zlsl_ratio_ttwj_stup_1b = (TH1F*) h_zl_ttwj_stup_1b->Clone( "h_zlsl_ratio_ttwj_stup_1b" ) ;
       TH1F* h_zlsl_ratio_ttwj_stup_2b = (TH1F*) h_zl_ttwj_stup_2b->Clone( "h_zlsl_ratio_ttwj_stup_2b" ) ;
       TH1F* h_zlsl_ratio_ttwj_stup_3b = (TH1F*) h_zl_ttwj_stup_3b->Clone( "h_zlsl_ratio_ttwj_stup_3b" ) ;

       h_zlsl_ratio_ttwj_stup_1b->Divide( h_sl_ttwj_stup_1b ) ;
       h_zlsl_ratio_ttwj_stup_2b->Divide( h_sl_ttwj_stup_2b ) ;
       h_zlsl_ratio_ttwj_stup_3b->Divide( h_sl_ttwj_stup_3b ) ;


       TH1F* h_zlsl_sf_ttwj_stup_1b = (TH1F*) h_zlsl_ratio_ttwj_stup_1b->Clone( "h_zlsl_ratio_ttwj_stup_1b" ) ;
       TH1F* h_zlsl_sf_ttwj_stup_2b = (TH1F*) h_zlsl_ratio_ttwj_stup_2b->Clone( "h_zlsl_ratio_ttwj_stup_2b" ) ;
       TH1F* h_zlsl_sf_ttwj_stup_3b = (TH1F*) h_zlsl_ratio_ttwj_stup_3b->Clone( "h_zlsl_ratio_ttwj_stup_3b" ) ;

       h_zlsl_sf_ttwj_stup_1b->Scale( 1./ave_zlsl_ratio_all ) ;
       h_zlsl_sf_ttwj_stup_2b->Scale( 1./ave_zlsl_ratio_all ) ;
       h_zlsl_sf_ttwj_stup_3b->Scale( 1./ave_zlsl_ratio_all ) ;



   //--- ZL/SL ratios and scale factors, single top -30%.


       TH1F* hmctruth_singletop_0lep_1b_dn = (TH1F*) hmctruth_singletop_0lep_1b->Clone( "hmctruth_singletop_0lep_1b_dn" ) ;
       TH1F* hmctruth_singletop_0lep_2b_dn = (TH1F*) hmctruth_singletop_0lep_2b->Clone( "hmctruth_singletop_0lep_2b_dn" ) ;
       TH1F* hmctruth_singletop_0lep_3b_dn = (TH1F*) hmctruth_singletop_0lep_3b->Clone( "hmctruth_singletop_0lep_3b_dn" ) ;
       TH1F* hmctruth_singletop_1lep_1b_dn = (TH1F*) hmctruth_singletop_1lep_1b->Clone( "hmctruth_singletop_1lep_1b_dn" ) ;
       TH1F* hmctruth_singletop_1lep_2b_dn = (TH1F*) hmctruth_singletop_1lep_2b->Clone( "hmctruth_singletop_1lep_2b_dn" ) ;
       TH1F* hmctruth_singletop_1lep_3b_dn = (TH1F*) hmctruth_singletop_1lep_3b->Clone( "hmctruth_singletop_1lep_3b_dn" ) ;

       hmctruth_singletop_0lep_1b_dn->Scale(0.7) ;
       hmctruth_singletop_0lep_2b_dn->Scale(0.7) ;
       hmctruth_singletop_0lep_3b_dn->Scale(0.7) ;
       hmctruth_singletop_1lep_1b_dn->Scale(0.7) ;
       hmctruth_singletop_1lep_2b_dn->Scale(0.7) ;
       hmctruth_singletop_1lep_3b_dn->Scale(0.7) ;

       TH1F* h_zl_ttwj_stdn_1b = (TH1F*) hmctruth_ttbar_0lep_1b->Clone( "h_zl_ttwj_stdn_1b" ) ;
       TH1F* h_zl_ttwj_stdn_2b = (TH1F*) hmctruth_ttbar_0lep_2b->Clone( "h_zl_ttwj_stdn_2b" ) ;
       TH1F* h_zl_ttwj_stdn_3b = (TH1F*) hmctruth_ttbar_0lep_3b->Clone( "h_zl_ttwj_stdn_3b" ) ;

       h_zl_ttwj_stdn_1b->Add( hmctruth_wjetsonly_0lep_1b ) ;
       h_zl_ttwj_stdn_2b->Add( hmctruth_wjetsonly_0lep_2b ) ;
       h_zl_ttwj_stdn_3b->Add( hmctruth_wjetsonly_0lep_3b ) ;

       h_zl_ttwj_stdn_1b->Add( hmctruth_singletop_0lep_1b_dn ) ;
       h_zl_ttwj_stdn_2b->Add( hmctruth_singletop_0lep_2b_dn ) ;
       h_zl_ttwj_stdn_3b->Add( hmctruth_singletop_0lep_3b_dn ) ;


       TH1F* h_sl_ttwj_stdn_1b = (TH1F*) hmctruth_ttbar_1lep_1b->Clone( "h_sl_ttwj_stdn_1b" ) ;
       TH1F* h_sl_ttwj_stdn_2b = (TH1F*) hmctruth_ttbar_1lep_2b->Clone( "h_sl_ttwj_stdn_2b" ) ;
       TH1F* h_sl_ttwj_stdn_3b = (TH1F*) hmctruth_ttbar_1lep_3b->Clone( "h_sl_ttwj_stdn_3b" ) ;

       h_sl_ttwj_stdn_1b->Add( hmctruth_wjetsonly_1lep_1b ) ;
       h_sl_ttwj_stdn_2b->Add( hmctruth_wjetsonly_1lep_2b ) ;
       h_sl_ttwj_stdn_3b->Add( hmctruth_wjetsonly_1lep_3b ) ;

       h_sl_ttwj_stdn_1b->Add( hmctruth_singletop_1lep_1b_dn ) ;
       h_sl_ttwj_stdn_2b->Add( hmctruth_singletop_1lep_2b_dn ) ;
       h_sl_ttwj_stdn_3b->Add( hmctruth_singletop_1lep_3b_dn ) ;


       TH1F* h_zlsl_ratio_ttwj_stdn_1b = (TH1F*) h_zl_ttwj_stdn_1b->Clone( "h_zlsl_ratio_ttwj_stdn_1b" ) ;
       TH1F* h_zlsl_ratio_ttwj_stdn_2b = (TH1F*) h_zl_ttwj_stdn_2b->Clone( "h_zlsl_ratio_ttwj_stdn_2b" ) ;
       TH1F* h_zlsl_ratio_ttwj_stdn_3b = (TH1F*) h_zl_ttwj_stdn_3b->Clone( "h_zlsl_ratio_ttwj_stdn_3b" ) ;

       h_zlsl_ratio_ttwj_stdn_1b->Divide( h_sl_ttwj_stdn_1b ) ;
       h_zlsl_ratio_ttwj_stdn_2b->Divide( h_sl_ttwj_stdn_2b ) ;
       h_zlsl_ratio_ttwj_stdn_3b->Divide( h_sl_ttwj_stdn_3b ) ;


       TH1F* h_zlsl_sf_ttwj_stdn_1b = (TH1F*) h_zlsl_ratio_ttwj_stdn_1b->Clone( "h_zlsl_ratio_ttwj_stdn_1b" ) ;
       TH1F* h_zlsl_sf_ttwj_stdn_2b = (TH1F*) h_zlsl_ratio_ttwj_stdn_2b->Clone( "h_zlsl_ratio_ttwj_stdn_2b" ) ;
       TH1F* h_zlsl_sf_ttwj_stdn_3b = (TH1F*) h_zlsl_ratio_ttwj_stdn_3b->Clone( "h_zlsl_ratio_ttwj_stdn_3b" ) ;

       h_zlsl_sf_ttwj_stdn_1b->Scale( 1./ave_zlsl_ratio_all ) ;
       h_zlsl_sf_ttwj_stdn_2b->Scale( 1./ave_zlsl_ratio_all ) ;
       h_zlsl_sf_ttwj_stdn_3b->Scale( 1./ave_zlsl_ratio_all ) ;


   //--- ZL/SL ratios and scale factors, Wjets +100%.


       TH1F* hmctruth_wjetsonly_0lep_1b_up = (TH1F*) hmctruth_wjetsonly_0lep_1b->Clone( "hmctruth_wjetsonly_0lep_1b_up" ) ;
       TH1F* hmctruth_wjetsonly_0lep_2b_up = (TH1F*) hmctruth_wjetsonly_0lep_2b->Clone( "hmctruth_wjetsonly_0lep_2b_up" ) ;
       TH1F* hmctruth_wjetsonly_0lep_3b_up = (TH1F*) hmctruth_wjetsonly_0lep_3b->Clone( "hmctruth_wjetsonly_0lep_3b_up" ) ;
       TH1F* hmctruth_wjetsonly_1lep_1b_up = (TH1F*) hmctruth_wjetsonly_1lep_1b->Clone( "hmctruth_wjetsonly_1lep_1b_up" ) ;
       TH1F* hmctruth_wjetsonly_1lep_2b_up = (TH1F*) hmctruth_wjetsonly_1lep_2b->Clone( "hmctruth_wjetsonly_1lep_2b_up" ) ;
       TH1F* hmctruth_wjetsonly_1lep_3b_up = (TH1F*) hmctruth_wjetsonly_1lep_3b->Clone( "hmctruth_wjetsonly_1lep_3b_up" ) ;

       hmctruth_wjetsonly_0lep_1b_up->Scale(2.0) ;
       hmctruth_wjetsonly_0lep_2b_up->Scale(2.0) ;
       hmctruth_wjetsonly_0lep_3b_up->Scale(2.0) ;
       hmctruth_wjetsonly_1lep_1b_up->Scale(2.0) ;
       hmctruth_wjetsonly_1lep_2b_up->Scale(2.0) ;
       hmctruth_wjetsonly_1lep_3b_up->Scale(2.0) ;

       TH1F* h_zl_ttwj_wjup_1b = (TH1F*) hmctruth_ttbar_0lep_1b->Clone( "h_zl_ttwj_wjup_1b" ) ;
       TH1F* h_zl_ttwj_wjup_2b = (TH1F*) hmctruth_ttbar_0lep_2b->Clone( "h_zl_ttwj_wjup_2b" ) ;
       TH1F* h_zl_ttwj_wjup_3b = (TH1F*) hmctruth_ttbar_0lep_3b->Clone( "h_zl_ttwj_wjup_3b" ) ;

       h_zl_ttwj_wjup_1b->Add( hmctruth_wjetsonly_0lep_1b_up ) ;
       h_zl_ttwj_wjup_2b->Add( hmctruth_wjetsonly_0lep_2b_up ) ;
       h_zl_ttwj_wjup_3b->Add( hmctruth_wjetsonly_0lep_3b_up ) ;

       h_zl_ttwj_wjup_1b->Add( hmctruth_singletop_0lep_1b ) ;
       h_zl_ttwj_wjup_2b->Add( hmctruth_singletop_0lep_2b ) ;
       h_zl_ttwj_wjup_3b->Add( hmctruth_singletop_0lep_3b ) ;


       TH1F* h_sl_ttwj_wjup_1b = (TH1F*) hmctruth_ttbar_1lep_1b->Clone( "h_sl_ttwj_wjup_1b" ) ;
       TH1F* h_sl_ttwj_wjup_2b = (TH1F*) hmctruth_ttbar_1lep_2b->Clone( "h_sl_ttwj_wjup_2b" ) ;
       TH1F* h_sl_ttwj_wjup_3b = (TH1F*) hmctruth_ttbar_1lep_3b->Clone( "h_sl_ttwj_wjup_3b" ) ;

       h_sl_ttwj_wjup_1b->Add( hmctruth_wjetsonly_1lep_1b_up ) ;
       h_sl_ttwj_wjup_2b->Add( hmctruth_wjetsonly_1lep_2b_up ) ;
       h_sl_ttwj_wjup_3b->Add( hmctruth_wjetsonly_1lep_3b_up ) ;

       h_sl_ttwj_wjup_1b->Add( hmctruth_singletop_1lep_1b ) ;
       h_sl_ttwj_wjup_2b->Add( hmctruth_singletop_1lep_2b ) ;
       h_sl_ttwj_wjup_3b->Add( hmctruth_singletop_1lep_3b ) ;


       TH1F* h_zlsl_ratio_ttwj_wjup_1b = (TH1F*) h_zl_ttwj_wjup_1b->Clone( "h_zlsl_ratio_ttwj_wjup_1b" ) ;
       TH1F* h_zlsl_ratio_ttwj_wjup_2b = (TH1F*) h_zl_ttwj_wjup_2b->Clone( "h_zlsl_ratio_ttwj_wjup_2b" ) ;
       TH1F* h_zlsl_ratio_ttwj_wjup_3b = (TH1F*) h_zl_ttwj_wjup_3b->Clone( "h_zlsl_ratio_ttwj_wjup_3b" ) ;

       h_zlsl_ratio_ttwj_wjup_1b->Divide( h_sl_ttwj_wjup_1b ) ;
       h_zlsl_ratio_ttwj_wjup_2b->Divide( h_sl_ttwj_wjup_2b ) ;
       h_zlsl_ratio_ttwj_wjup_3b->Divide( h_sl_ttwj_wjup_3b ) ;


       TH1F* h_zlsl_sf_ttwj_wjup_1b = (TH1F*) h_zlsl_ratio_ttwj_wjup_1b->Clone( "h_zlsl_ratio_ttwj_wjup_1b" ) ;
       TH1F* h_zlsl_sf_ttwj_wjup_2b = (TH1F*) h_zlsl_ratio_ttwj_wjup_2b->Clone( "h_zlsl_ratio_ttwj_wjup_2b" ) ;
       TH1F* h_zlsl_sf_ttwj_wjup_3b = (TH1F*) h_zlsl_ratio_ttwj_wjup_3b->Clone( "h_zlsl_ratio_ttwj_wjup_3b" ) ;

       h_zlsl_sf_ttwj_wjup_1b->Scale( 1./ave_zlsl_ratio_all ) ;
       h_zlsl_sf_ttwj_wjup_2b->Scale( 1./ave_zlsl_ratio_all ) ;
       h_zlsl_sf_ttwj_wjup_3b->Scale( 1./ave_zlsl_ratio_all ) ;




   //--- ZL/SL ratios and scale factors, Wjets -100%.


       TH1F* hmctruth_wjetsonly_0lep_1b_dn = (TH1F*) hmctruth_wjetsonly_0lep_1b->Clone( "hmctruth_wjetsonly_0lep_1b_dn" ) ;
       TH1F* hmctruth_wjetsonly_0lep_2b_dn = (TH1F*) hmctruth_wjetsonly_0lep_2b->Clone( "hmctruth_wjetsonly_0lep_2b_dn" ) ;
       TH1F* hmctruth_wjetsonly_0lep_3b_dn = (TH1F*) hmctruth_wjetsonly_0lep_3b->Clone( "hmctruth_wjetsonly_0lep_3b_dn" ) ;
       TH1F* hmctruth_wjetsonly_1lep_1b_dn = (TH1F*) hmctruth_wjetsonly_1lep_1b->Clone( "hmctruth_wjetsonly_1lep_1b_dn" ) ;
       TH1F* hmctruth_wjetsonly_1lep_2b_dn = (TH1F*) hmctruth_wjetsonly_1lep_2b->Clone( "hmctruth_wjetsonly_1lep_2b_dn" ) ;
       TH1F* hmctruth_wjetsonly_1lep_3b_dn = (TH1F*) hmctruth_wjetsonly_1lep_3b->Clone( "hmctruth_wjetsonly_1lep_3b_dn" ) ;

       hmctruth_wjetsonly_0lep_1b_dn->Scale(0.0) ;
       hmctruth_wjetsonly_0lep_2b_dn->Scale(0.0) ;
       hmctruth_wjetsonly_0lep_3b_dn->Scale(0.0) ;
       hmctruth_wjetsonly_1lep_1b_dn->Scale(0.0) ;
       hmctruth_wjetsonly_1lep_2b_dn->Scale(0.0) ;
       hmctruth_wjetsonly_1lep_3b_dn->Scale(0.0) ;

       TH1F* h_zl_ttwj_wjdn_1b = (TH1F*) hmctruth_ttbar_0lep_1b->Clone( "h_zl_ttwj_wjdn_1b" ) ;
       TH1F* h_zl_ttwj_wjdn_2b = (TH1F*) hmctruth_ttbar_0lep_2b->Clone( "h_zl_ttwj_wjdn_2b" ) ;
       TH1F* h_zl_ttwj_wjdn_3b = (TH1F*) hmctruth_ttbar_0lep_3b->Clone( "h_zl_ttwj_wjdn_3b" ) ;

       h_zl_ttwj_wjdn_1b->Add( hmctruth_wjetsonly_0lep_1b_dn ) ;
       h_zl_ttwj_wjdn_2b->Add( hmctruth_wjetsonly_0lep_2b_dn ) ;
       h_zl_ttwj_wjdn_3b->Add( hmctruth_wjetsonly_0lep_3b_dn ) ;

       h_zl_ttwj_wjdn_1b->Add( hmctruth_singletop_0lep_1b ) ;
       h_zl_ttwj_wjdn_2b->Add( hmctruth_singletop_0lep_2b ) ;
       h_zl_ttwj_wjdn_3b->Add( hmctruth_singletop_0lep_3b ) ;


       TH1F* h_sl_ttwj_wjdn_1b = (TH1F*) hmctruth_ttbar_1lep_1b->Clone( "h_sl_ttwj_wjdn_1b" ) ;
       TH1F* h_sl_ttwj_wjdn_2b = (TH1F*) hmctruth_ttbar_1lep_2b->Clone( "h_sl_ttwj_wjdn_2b" ) ;
       TH1F* h_sl_ttwj_wjdn_3b = (TH1F*) hmctruth_ttbar_1lep_3b->Clone( "h_sl_ttwj_wjdn_3b" ) ;

       h_sl_ttwj_wjdn_1b->Add( hmctruth_wjetsonly_1lep_1b_dn ) ;
       h_sl_ttwj_wjdn_2b->Add( hmctruth_wjetsonly_1lep_2b_dn ) ;
       h_sl_ttwj_wjdn_3b->Add( hmctruth_wjetsonly_1lep_3b_dn ) ;

       h_sl_ttwj_wjdn_1b->Add( hmctruth_singletop_1lep_1b ) ;
       h_sl_ttwj_wjdn_2b->Add( hmctruth_singletop_1lep_2b ) ;
       h_sl_ttwj_wjdn_3b->Add( hmctruth_singletop_1lep_3b ) ;


       TH1F* h_zlsl_ratio_ttwj_wjdn_1b = (TH1F*) h_zl_ttwj_wjdn_1b->Clone( "h_zlsl_ratio_ttwj_wjdn_1b" ) ;
       TH1F* h_zlsl_ratio_ttwj_wjdn_2b = (TH1F*) h_zl_ttwj_wjdn_2b->Clone( "h_zlsl_ratio_ttwj_wjdn_2b" ) ;
       TH1F* h_zlsl_ratio_ttwj_wjdn_3b = (TH1F*) h_zl_ttwj_wjdn_3b->Clone( "h_zlsl_ratio_ttwj_wjdn_3b" ) ;

       h_zlsl_ratio_ttwj_wjdn_1b->Divide( h_sl_ttwj_wjdn_1b ) ;
       h_zlsl_ratio_ttwj_wjdn_2b->Divide( h_sl_ttwj_wjdn_2b ) ;
       h_zlsl_ratio_ttwj_wjdn_3b->Divide( h_sl_ttwj_wjdn_3b ) ;


       TH1F* h_zlsl_sf_ttwj_wjdn_1b = (TH1F*) h_zlsl_ratio_ttwj_wjdn_1b->Clone( "h_zlsl_ratio_ttwj_wjdn_1b" ) ;
       TH1F* h_zlsl_sf_ttwj_wjdn_2b = (TH1F*) h_zlsl_ratio_ttwj_wjdn_2b->Clone( "h_zlsl_ratio_ttwj_wjdn_2b" ) ;
       TH1F* h_zlsl_sf_ttwj_wjdn_3b = (TH1F*) h_zlsl_ratio_ttwj_wjdn_3b->Clone( "h_zlsl_ratio_ttwj_wjdn_3b" ) ;

       h_zlsl_sf_ttwj_wjdn_1b->Scale( 1./ave_zlsl_ratio_all ) ;
       h_zlsl_sf_ttwj_wjdn_2b->Scale( 1./ave_zlsl_ratio_all ) ;
       h_zlsl_sf_ttwj_wjdn_3b->Scale( 1./ave_zlsl_ratio_all ) ;







    //--- Ratios


       TCanvas* can4 = new TCanvas("can4", "ZL/SL ratio, single top variation (top), Wjets variation (bottom)", 1000, 700 ) ;

       can4->Divide(3,2) ;

       can4->cd(1) ;
       hist1 = h_zlsl_ratio_ttwj_nom_1b ;
       hist2 = h_zlsl_ratio_ttwj_stup_1b ;
       hist3 = h_zlsl_ratio_ttwj_stdn_1b ;
       hist2->SetLineColor(2) ;
       hist3->SetLineColor(4) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;


       can4->cd(2) ;
       hist1 = h_zlsl_ratio_ttwj_nom_2b ;
       hist2 = h_zlsl_ratio_ttwj_stup_2b ;
       hist3 = h_zlsl_ratio_ttwj_stdn_2b ;
       hist2->SetLineColor(2) ;
       hist3->SetLineColor(4) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;


       can4->cd(3) ;
       hist1 = h_zlsl_ratio_ttwj_nom_3b ;
       hist2 = h_zlsl_ratio_ttwj_stup_3b ;
       hist3 = h_zlsl_ratio_ttwj_stdn_3b ;
       hist2->SetLineColor(2) ;
       hist3->SetLineColor(4) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;


       can4->cd(4) ;
       hist1 = h_zlsl_ratio_ttwj_nom_1b ;
       hist2 = h_zlsl_ratio_ttwj_wjup_1b ;
       hist3 = h_zlsl_ratio_ttwj_wjdn_1b ;
       hist2->SetLineColor(2) ;
       hist3->SetLineColor(4) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;


       can4->cd(5) ;
       hist1 = h_zlsl_ratio_ttwj_nom_2b ;
       hist2 = h_zlsl_ratio_ttwj_wjup_2b ;
       hist3 = h_zlsl_ratio_ttwj_wjdn_2b ;
       hist2->SetLineColor(2) ;
       hist3->SetLineColor(4) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;


       can4->cd(6) ;
       hist1 = h_zlsl_ratio_ttwj_nom_3b ;
       hist2 = h_zlsl_ratio_ttwj_wjup_3b ;
       hist3 = h_zlsl_ratio_ttwj_wjdn_3b ;
       hist2->SetLineColor(2) ;
       hist3->SetLineColor(4) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;



    //--- SFs



       TCanvas* can5 = new TCanvas("can5", "ZL/SL SF, single top variation (top), Wjets variation (bottom)", 1000, 700 ) ;

       can5->Divide(3,2) ;

       can5->cd(1) ;
       hist1 = h_zlsl_sf_ttwj_nom_1b ;
       hist2 = h_zlsl_sf_ttwj_stup_1b ;
       hist3 = h_zlsl_sf_ttwj_stdn_1b ;
       hist2->SetLineColor(2) ;
       hist3->SetLineColor(4) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;


       can5->cd(2) ;
       hist1 = h_zlsl_sf_ttwj_nom_2b ;
       hist2 = h_zlsl_sf_ttwj_stup_2b ;
       hist3 = h_zlsl_sf_ttwj_stdn_2b ;
       hist2->SetLineColor(2) ;
       hist3->SetLineColor(4) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;


       can5->cd(3) ;
       hist1 = h_zlsl_sf_ttwj_nom_3b ;
       hist2 = h_zlsl_sf_ttwj_stup_3b ;
       hist3 = h_zlsl_sf_ttwj_stdn_3b ;
       hist2->SetLineColor(2) ;
       hist3->SetLineColor(4) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;


       can5->cd(4) ;
       hist1 = h_zlsl_sf_ttwj_nom_1b ;
       hist2 = h_zlsl_sf_ttwj_wjup_1b ;
       hist3 = h_zlsl_sf_ttwj_wjdn_1b ;
       hist2->SetLineColor(2) ;
       hist3->SetLineColor(4) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;


       can5->cd(5) ;
       hist1 = h_zlsl_sf_ttwj_nom_2b ;
       hist2 = h_zlsl_sf_ttwj_wjup_2b ;
       hist3 = h_zlsl_sf_ttwj_wjdn_2b ;
       hist2->SetLineColor(2) ;
       hist3->SetLineColor(4) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;


       can5->cd(6) ;
       hist1 = h_zlsl_sf_ttwj_nom_3b ;
       hist2 = h_zlsl_sf_ttwj_wjup_3b ;
       hist3 = h_zlsl_sf_ttwj_wjdn_3b ;
       hist2->SetLineColor(2) ;
       hist3->SetLineColor(4) ;
       hist1->Draw() ;
       hist2->Draw("same") ;
       hist3->Draw("same") ;



   //--- Shape syst
   //
   //  delta_+ = var+ / nom - 1
   //
   //  delta_- = var- / nom - 1
   //
   //  <delta>_+ = (delta_+ - delta_-)/2 = (1/2)(var+ - var-)/nom
   //
   //

       TH1F* h_zlsl_sf_shapesyst_singletop_1b = (TH1F*) h_zlsl_sf_ttwj_stup_1b->Clone( "h_zlsl_sf_shapesyst_singletop_1b" ) ;
       h_zlsl_sf_shapesyst_singletop_1b -> Add( h_zlsl_sf_ttwj_stdn_1b, -1. ) ;
       h_zlsl_sf_shapesyst_singletop_1b -> Divide( h_zlsl_sf_ttwj_nom_1b ) ;
       h_zlsl_sf_shapesyst_singletop_1b -> Scale( 0.5 ) ;

       TH1F* h_zlsl_sf_shapesyst_singletop_2b = (TH1F*) h_zlsl_sf_ttwj_stup_2b->Clone( "h_zlsl_sf_shapesyst_singletop_2b" ) ;
       h_zlsl_sf_shapesyst_singletop_2b -> Add( h_zlsl_sf_ttwj_stdn_2b, -1. ) ;
       h_zlsl_sf_shapesyst_singletop_2b -> Divide( h_zlsl_sf_ttwj_nom_2b ) ;
       h_zlsl_sf_shapesyst_singletop_2b -> Scale( 0.5 ) ;

       TH1F* h_zlsl_sf_shapesyst_singletop_3b = (TH1F*) h_zlsl_sf_ttwj_stup_3b->Clone( "h_zlsl_sf_shapesyst_singletop_3b" ) ;
       h_zlsl_sf_shapesyst_singletop_3b -> Add( h_zlsl_sf_ttwj_stdn_3b, -1. ) ;
       h_zlsl_sf_shapesyst_singletop_3b -> Divide( h_zlsl_sf_ttwj_nom_3b ) ;
       h_zlsl_sf_shapesyst_singletop_3b -> Scale( 0.5 ) ;


       TH1F* h_zlsl_sf_shapesyst_wjets_1b = (TH1F*) h_zlsl_sf_ttwj_wjup_1b->Clone( "h_zlsl_sf_shapesyst_wjets_1b" ) ;
       h_zlsl_sf_shapesyst_wjets_1b -> Add( h_zlsl_sf_ttwj_wjdn_1b, -1. ) ;
       h_zlsl_sf_shapesyst_wjets_1b -> Divide( h_zlsl_sf_ttwj_nom_1b ) ;
       h_zlsl_sf_shapesyst_wjets_1b -> Scale( 0.5 ) ;

       TH1F* h_zlsl_sf_shapesyst_wjets_2b = (TH1F*) h_zlsl_sf_ttwj_wjup_2b->Clone( "h_zlsl_sf_shapesyst_wjets_2b" ) ;
       h_zlsl_sf_shapesyst_wjets_2b -> Add( h_zlsl_sf_ttwj_wjdn_2b, -1. ) ;
       h_zlsl_sf_shapesyst_wjets_2b -> Divide( h_zlsl_sf_ttwj_nom_2b ) ;
       h_zlsl_sf_shapesyst_wjets_2b -> Scale( 0.5 ) ;

       TH1F* h_zlsl_sf_shapesyst_wjets_3b = (TH1F*) h_zlsl_sf_ttwj_wjup_3b->Clone( "h_zlsl_sf_shapesyst_wjets_3b" ) ;
       h_zlsl_sf_shapesyst_wjets_3b -> Add( h_zlsl_sf_ttwj_wjdn_3b, -1. ) ;
       h_zlsl_sf_shapesyst_wjets_3b -> Divide( h_zlsl_sf_ttwj_nom_3b ) ;
       h_zlsl_sf_shapesyst_wjets_3b -> Scale( 0.5 ) ;


       TCanvas* can6 = new TCanvas("can6", "ZL/SL SF, single top shapesyst (top), Wjets shapesyst (bottom)", 1000, 700 ) ;

       can6->Divide(3,2) ;

       can6->cd(1) ;
       hist1 = h_zlsl_sf_shapesyst_singletop_1b ;
       hist1 -> SetFillColor(11) ;
       hist1 -> SetMinimum(-0.1) ;
       hist1 -> SetMaximum(+0.1) ;
       hist1->Draw("hist") ;

       can6->cd(2) ;
       hist1 = h_zlsl_sf_shapesyst_singletop_2b ;
       hist1 -> SetFillColor(11) ;
       hist1 -> SetMinimum(-0.1) ;
       hist1 -> SetMaximum(+0.1) ;
       hist1->Draw("hist") ;

       can6->cd(3) ;
       hist1 = h_zlsl_sf_shapesyst_singletop_3b ;
       hist1 -> SetFillColor(11) ;
       hist1 -> SetMinimum(-0.1) ;
       hist1 -> SetMaximum(+0.1) ;
       hist1->Draw("hist") ;

       can6->cd(4) ;
       hist1 = h_zlsl_sf_shapesyst_wjets_1b ;
       hist1 -> SetFillColor(11) ;
       hist1 -> SetMinimum(-0.1) ;
       hist1 -> SetMaximum(+0.1) ;
       hist1->Draw("hist") ;

       can6->cd(5) ;
       hist1 = h_zlsl_sf_shapesyst_wjets_2b ;
       hist1 -> SetFillColor(11) ;
       hist1 -> SetMinimum(-0.1) ;
       hist1 -> SetMaximum(+0.1) ;
       hist1->Draw("hist") ;

       can6->cd(6) ;
       hist1 = h_zlsl_sf_shapesyst_wjets_3b ;
       hist1 -> SetFillColor(11) ;
       hist1 -> SetMinimum(-0.1) ;
       hist1 -> SetMaximum(+0.1) ;
       hist1->Draw("hist") ;



    } // draw_ttwjfracs


//================================================================

void loadHist(const char* filename, const char* pfx, const char* pat, Bool_t doAdd, Double_t scaleFactor)
{
  cout << " Reading histograms from file: " << filename << endl ;
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
   //// ((TH1*)clone)->Add((TH1*)obj) ;
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
 /////  ((TH1*)clone)->Add((TH1*)obj) ;
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

//------------------------------------------------------------------------



