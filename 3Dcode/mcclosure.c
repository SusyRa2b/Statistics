
#include "TH1F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TKey.h"


#include <iostream>

  using std::cout ;
  using std::endl ;

   void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;

  //----------------------------------------

   void mcclosure( const char* infile = "rootfiles/gi-plots-met3-ht3.root" ) {

      TString infileStr( infile ) ;

      gSystem->Exec("mkdir -p outputfiles") ;

      gStyle->SetOptStat(0) ;
      gStyle->SetPadTopMargin(0.03) ;
      gStyle->SetPadBottomMargin(0.30) ;
      gStyle->SetPadRightMargin(0.10) ;
      gStyle->SetTitleX(0.98) ;
      gStyle->SetTitleAlign(33) ;

      gStyle->SetPadGridY(1) ;

      TLine* line = new TLine() ;
      line->SetLineStyle(1) ;
      line->SetLineWidth(2) ;
      line->SetLineColor(2) ;

      gDirectory->Delete("h*") ;

      loadHist( infile ) ;

      TH1F* hmctruth_ttwj_0lep_1b = (TH1F*) gDirectory->FindObject("hmctruth_ttwj_0lep_1b") ;
      if ( hmctruth_ttwj_0lep_1b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_ttwj_0lep_1b.\n\n") ; return ; }
      TH1F* hmctruth_ttwj_1lep_1b = (TH1F*) gDirectory->FindObject("hmctruth_ttwj_1lep_1b") ;
      if ( hmctruth_ttwj_1lep_1b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_ttwj_1lep_1b.\n\n") ; return ; }

      TH1F* hmctruth_ttwj_0over1ratio_1b = (TH1F*) hmctruth_ttwj_0lep_1b->Clone("hmctruth_ttwj_0over1ratio_1b") ;
      hmctruth_ttwj_0over1ratio_1b->Divide( hmctruth_ttwj_1lep_1b ) ;


      TH1F* hmctruth_ttwj_0lep_2b = (TH1F*) gDirectory->FindObject("hmctruth_ttwj_0lep_2b") ;
      if ( hmctruth_ttwj_0lep_2b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_ttwj_0lep_2b.\n\n") ; return ; }
      TH1F* hmctruth_ttwj_1lep_2b = (TH1F*) gDirectory->FindObject("hmctruth_ttwj_1lep_2b") ;
      if ( hmctruth_ttwj_1lep_2b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_ttwj_1lep_2b.\n\n") ; return ; }

      TH1F* hmctruth_ttwj_0over1ratio_2b = (TH1F*) hmctruth_ttwj_0lep_2b->Clone("hmctruth_ttwj_0over1ratio_2b") ;
      hmctruth_ttwj_0over1ratio_2b->Divide( hmctruth_ttwj_1lep_2b ) ;


      TH1F* hmctruth_ttwj_0lep_3b = (TH1F*) gDirectory->FindObject("hmctruth_ttwj_0lep_3b") ;
      if ( hmctruth_ttwj_0lep_3b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_ttwj_0lep_3b.\n\n") ; return ; }
      TH1F* hmctruth_ttwj_1lep_3b = (TH1F*) gDirectory->FindObject("hmctruth_ttwj_1lep_3b") ;
      if ( hmctruth_ttwj_1lep_3b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_ttwj_1lep_3b.\n\n") ; return ; }

      TH1F* hmctruth_ttwj_0over1ratio_3b = (TH1F*) hmctruth_ttwj_0lep_3b->Clone("hmctruth_ttwj_0over1ratio_3b") ;
      hmctruth_ttwj_0over1ratio_3b->Divide( hmctruth_ttwj_1lep_3b ) ;

      hmctruth_ttwj_0over1ratio_1b->SetLineColor(2) ;
      hmctruth_ttwj_0over1ratio_2b->SetLineColor(6) ;
      hmctruth_ttwj_0over1ratio_3b->SetLineColor(4) ;

      hmctruth_ttwj_0over1ratio_1b->SetMarkerStyle(20) ;
      hmctruth_ttwj_0over1ratio_2b->SetMarkerStyle(25) ;
      hmctruth_ttwj_0over1ratio_3b->SetMarkerStyle(30) ;




      //-- compute dumb ave ratio

      int nBins = hmctruth_ttwj_0over1ratio_1b->GetNbinsX() ;

      double total0lep(0.) ;
      double total1lep(0.) ;
      double sumw20lep(0.) ;
      double sumw21lep(0.) ;

      for ( int bi=1; bi<=nBins; bi++ ) {

         total0lep += hmctruth_ttwj_0lep_1b->GetBinContent( bi ) ;
         total0lep += hmctruth_ttwj_0lep_2b->GetBinContent( bi ) ;
         total0lep += hmctruth_ttwj_0lep_3b->GetBinContent( bi ) ;

         total1lep += hmctruth_ttwj_1lep_1b->GetBinContent( bi ) ;
         total1lep += hmctruth_ttwj_1lep_2b->GetBinContent( bi ) ;
         total1lep += hmctruth_ttwj_1lep_3b->GetBinContent( bi ) ;

         sumw20lep += pow( hmctruth_ttwj_0lep_1b->GetBinError( bi ), 2 ) ;
         sumw20lep += pow( hmctruth_ttwj_0lep_2b->GetBinError( bi ), 2 ) ;
         sumw20lep += pow( hmctruth_ttwj_0lep_3b->GetBinError( bi ), 2 ) ;

         sumw21lep += pow( hmctruth_ttwj_1lep_1b->GetBinError( bi ), 2 ) ;
         sumw21lep += pow( hmctruth_ttwj_1lep_2b->GetBinError( bi ), 2 ) ;
         sumw21lep += pow( hmctruth_ttwj_1lep_3b->GetBinError( bi ), 2 ) ;

      } // bi.

      double simpleAveR_0over1 = total0lep / total1lep ;
      double simpleAveR_0over1_err = simpleAveR_0over1 * sqrt( sumw20lep/(total0lep*total0lep) + sumw21lep/(total1lep*total1lep)  ) ;

      printf("\n\n Simple average 0lep/1lep = %5.3f +/- %5.3f\n\n", simpleAveR_0over1, simpleAveR_0over1_err ) ;


      TH1F* hscalefactor_ttwj_0over1ratio_1b = (TH1F*) hmctruth_ttwj_0over1ratio_1b->Clone("hscalefactor_ttwj_0lep_1b") ;
      TH1F* hscalefactor_ttwj_0over1ratio_2b = (TH1F*) hmctruth_ttwj_0over1ratio_2b->Clone("hscalefactor_ttwj_0lep_2b") ;
      TH1F* hscalefactor_ttwj_0over1ratio_3b = (TH1F*) hmctruth_ttwj_0over1ratio_3b->Clone("hscalefactor_ttwj_0lep_3b") ;

      hscalefactor_ttwj_0over1ratio_1b->Scale(1./simpleAveR_0over1) ;
      hscalefactor_ttwj_0over1ratio_2b->Scale(1./simpleAveR_0over1) ;
      hscalefactor_ttwj_0over1ratio_3b->Scale(1./simpleAveR_0over1) ;


      hmctruth_ttwj_0over1ratio_1b->SetMinimum(0.) ;
      hmctruth_ttwj_0over1ratio_2b->SetMinimum(0.) ;
      hmctruth_ttwj_0over1ratio_3b->SetMinimum(0.) ;
      hmctruth_ttwj_0over1ratio_1b->SetMaximum(5.) ;
      hmctruth_ttwj_0over1ratio_2b->SetMaximum(5.) ;
      hmctruth_ttwj_0over1ratio_3b->SetMaximum(5.) ;

      hscalefactor_ttwj_0over1ratio_1b->SetMinimum(0.) ;
      hscalefactor_ttwj_0over1ratio_2b->SetMinimum(0.) ;
      hscalefactor_ttwj_0over1ratio_3b->SetMinimum(0.) ;
      hscalefactor_ttwj_0over1ratio_1b->SetMaximum(5.) ;
      hscalefactor_ttwj_0over1ratio_2b->SetMaximum(5.) ;
      hscalefactor_ttwj_0over1ratio_3b->SetMaximum(5.) ;




      TCanvas* cttwj = (TCanvas*) gDirectory->FindObject("cttwj") ;
      if ( cttwj == 0x0 ) {
         cttwj = new TCanvas("cttwj","ttwj closure", 700, 950) ;
      }
      cttwj->Clear() ;
      cttwj->Divide(1,2) ;

      TLegend* legend_ttwj = new TLegend( 0.88, 0.77,  0.97, 0.93 ) ;
      legend_ttwj->SetFillColor(kWhite) ;
      legend_ttwj->AddEntry( hmctruth_ttwj_0over1ratio_1b, "=1b") ;
      legend_ttwj->AddEntry( hmctruth_ttwj_0over1ratio_2b, "=2b") ;
      legend_ttwj->AddEntry( hmctruth_ttwj_0over1ratio_3b, ">=3b") ;

      cttwj->cd(1) ;
      hmctruth_ttwj_0over1ratio_1b->SetTitle("ttwj: 0 Lepton / 1 Lepton, Ratio") ;
      hmctruth_ttwj_0over1ratio_1b->Draw() ;
      hmctruth_ttwj_0over1ratio_2b->Draw("same") ;
      hmctruth_ttwj_0over1ratio_3b->Draw("same") ;
      legend_ttwj->Draw() ;

      cttwj->cd(2) ;
      hscalefactor_ttwj_0over1ratio_1b->SetTitle("ttwj: 0 Lepton / 1 Lepton, Scale Factor") ;
      hscalefactor_ttwj_0over1ratio_1b->Draw() ;
      hscalefactor_ttwj_0over1ratio_2b->Draw("same") ;
      hscalefactor_ttwj_0over1ratio_3b->Draw("same") ;
      line->DrawLine(0.5,1.,hscalefactor_ttwj_0over1ratio_1b->GetNbinsX()+0.5,1.) ;
      legend_ttwj->Draw() ;

      TString outttwj( infileStr ) ;
      outttwj.ReplaceAll("rootfiles","outputfiles") ;
      outttwj.ReplaceAll(".root","-mcclosure-ttwj1.pdf") ;
      cttwj->SaveAs( outttwj ) ;






      TCanvas* cttwj3 = (TCanvas*) gDirectory->FindObject("cttwj3") ;
      if ( cttwj3 == 0x0 ) {
         cttwj3 = new TCanvas("cttwj3","ttwj closure", 1200, 950) ;
      }
      cttwj3->Clear() ;
      cttwj3->Divide(3,2) ;

      cttwj3->cd(1) ;
      hmctruth_ttwj_0over1ratio_1b->SetTitle("ttwj =1b: 0 Lepton / 1 Lepton, Ratio") ;
      hmctruth_ttwj_0over1ratio_1b->Draw() ;
      cttwj3->cd(2) ;
      hmctruth_ttwj_0over1ratio_2b->SetTitle("ttwj =2b: 0 Lepton / 1 Lepton, Ratio") ;
      hmctruth_ttwj_0over1ratio_2b->Draw() ;
      cttwj3->cd(3) ;
      hmctruth_ttwj_0over1ratio_3b->SetTitle("ttwj >=3b: 0 Lepton / 1 Lepton, Ratio") ;
      hmctruth_ttwj_0over1ratio_3b->Draw() ;

      cttwj3->cd(4) ;
      hscalefactor_ttwj_0over1ratio_1b->SetTitle("ttwj =1b: 0 Lepton / 1 Lepton, Scale Factor") ;
      hscalefactor_ttwj_0over1ratio_1b->Draw() ;
      line->DrawLine(0.5,1.,hscalefactor_ttwj_0over1ratio_1b->GetNbinsX()+0.5,1.) ;
      hscalefactor_ttwj_0over1ratio_1b->Draw("same") ;
      cttwj3->cd(5) ;
      hscalefactor_ttwj_0over1ratio_2b->SetTitle("ttwj =2b: 0 Lepton / 1 Lepton, Scale Factor") ;
      hscalefactor_ttwj_0over1ratio_2b->Draw() ;
      line->DrawLine(0.5,1.,hscalefactor_ttwj_0over1ratio_1b->GetNbinsX()+0.5,1.) ;
      hscalefactor_ttwj_0over1ratio_2b->Draw("same") ;
      cttwj3->cd(6) ;
      hscalefactor_ttwj_0over1ratio_3b->SetTitle("ttwj >=3b: 0 Lepton / 1 Lepton, Scale Factor") ;
      hscalefactor_ttwj_0over1ratio_3b->Draw() ;
      line->DrawLine(0.5,1.,hscalefactor_ttwj_0over1ratio_1b->GetNbinsX()+0.5,1.) ;
      hscalefactor_ttwj_0over1ratio_3b->Draw("same") ;


      TString outttwj3( infileStr ) ;
      outttwj3.ReplaceAll("rootfiles","outputfiles") ;
      outttwj3.ReplaceAll(".root","-mcclosure-ttwj3.pdf") ;
      cttwj3->SaveAs( outttwj3 ) ;



     //---   Q C D  Part   --------------



      //--- get LSB ratios

      TH1F* hmctruth_qcd_lsb_pass = (TH1F*) gDirectory->FindObject("hmctruth_qcd_lsb_pass") ;
      if ( hmctruth_qcd_lsb_pass == 0x0 ) { printf("\n\n *** Can't find hmctruth_qcd_lsb_pass.\n\n") ; return ; }
      TH1F* hmctruth_qcd_lsb_fail = (TH1F*) gDirectory->FindObject("hmctruth_qcd_lsb_fail") ;
      if ( hmctruth_qcd_lsb_fail == 0x0 ) { printf("\n\n *** Can't find hmctruth_qcd_lsb_fail.\n\n") ; return ; }

      int nBinsBjets(3) ;
      int nBinsHT = (hmctruth_qcd_lsb_pass->GetNbinsX() - 1) / nBinsBjets - 1 ;
      printf("\n\n Number of HT bins : %d\n\n", nBinsHT ) ;

      double R_lsb[3][20] ;
      double R_lsb_err[3][20] ;

      for ( int nbi=0; nbi<nBinsBjets; nbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            int hbin = 1 + (nBinsHT+1)*nbi + hbi + 1 ;
            double pass = hmctruth_qcd_lsb_pass->GetBinContent( hbin ) ;
            double fail = hmctruth_qcd_lsb_fail->GetBinContent( hbin ) ;
            double pass_err = hmctruth_qcd_lsb_pass->GetBinError( hbin ) ;
            double fail_err = hmctruth_qcd_lsb_fail->GetBinError( hbin ) ;
            double ratio(0.) ;
            double err(0.) ;
            if ( fail > 0. ) { ratio = pass / fail ; }
            if ( fail > 0. && fail_err > 0. && pass > 0. && pass_err > 0. ) {
               err = ratio * sqrt( pow(pass_err/pass,2) + pow(fail_err/fail,2) ) ;
            }
            R_lsb[nbi][hbi] = ratio ;
            R_lsb_err[nbi][hbi] = err ;
         } // hbi.
      } // nbi.

      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         printf("  HT bin %d: LSB pass/fail ratio R_lsb :  %6.4f +/- %6.4f\n", hbi+1, R_lsb[0][hbi], R_lsb_err[0][hbi] ) ;
      }
      printf("\n\n") ;



      TH1F* hmctruth_qcd_0lep_1b = (TH1F*) gDirectory->FindObject("hmctruth_qcd_0lep_1b") ;
      if ( hmctruth_qcd_0lep_1b == 0x0 ) { printf("\n\n *** Can't find hmctruth_qcd_0lep_1b.\n\n") ; return ; }
      TH1F* hmctruth_qcd_ldp_1b = (TH1F*) gDirectory->FindObject("hmctruth_qcd_ldp_1b") ;
      if ( hmctruth_qcd_ldp_1b == 0x0 ) { printf("\n\n *** Can't find hmctruth_qcd_ldp_1b.\n\n") ; return ; }

      TH1F* hmctruth_qcd_0lep_2b = (TH1F*) gDirectory->FindObject("hmctruth_qcd_0lep_2b") ;
      if ( hmctruth_qcd_0lep_2b == 0x0 ) { printf("\n\n *** Can't find hmctruth_qcd_0lep_2b.\n\n") ; return ; }
      TH1F* hmctruth_qcd_ldp_2b = (TH1F*) gDirectory->FindObject("hmctruth_qcd_ldp_2b") ;
      if ( hmctruth_qcd_ldp_2b == 0x0 ) { printf("\n\n *** Can't find hmctruth_qcd_ldp_2b.\n\n") ; return ; }

      TH1F* hmctruth_qcd_0lep_3b = (TH1F*) gDirectory->FindObject("hmctruth_qcd_0lep_3b") ;
      if ( hmctruth_qcd_0lep_3b == 0x0 ) { printf("\n\n *** Can't find hmctruth_qcd_0lep_3b.\n\n") ; return ; }
      TH1F* hmctruth_qcd_ldp_3b = (TH1F*) gDirectory->FindObject("hmctruth_qcd_ldp_3b") ;
      if ( hmctruth_qcd_ldp_3b == 0x0 ) { printf("\n\n *** Can't find hmctruth_qcd_ldp_3b.\n\n") ; return ; }

      TH1F* hmctruth_qcd_0lepoverldpratio_1b = (TH1F*) hmctruth_qcd_0lep_1b->Clone("hmctruth_qcd_0lepoverldpratio_1b") ;
      TH1F* hmctruth_qcd_0lepoverldpratio_2b = (TH1F*) hmctruth_qcd_0lep_2b->Clone("hmctruth_qcd_0lepoverldpratio_2b") ;
      TH1F* hmctruth_qcd_0lepoverldpratio_3b = (TH1F*) hmctruth_qcd_0lep_3b->Clone("hmctruth_qcd_0lepoverldpratio_3b") ;

      hmctruth_qcd_0lepoverldpratio_1b->Divide( hmctruth_qcd_ldp_1b ) ;
      hmctruth_qcd_0lepoverldpratio_2b->Divide( hmctruth_qcd_ldp_2b ) ;
      hmctruth_qcd_0lepoverldpratio_3b->Divide( hmctruth_qcd_ldp_3b ) ;

      hmctruth_qcd_0lepoverldpratio_1b->SetLineColor(2) ;
      hmctruth_qcd_0lepoverldpratio_2b->SetLineColor(6) ;
      hmctruth_qcd_0lepoverldpratio_3b->SetLineColor(4) ;

      hmctruth_qcd_0lepoverldpratio_1b->SetMarkerStyle(20) ;
      hmctruth_qcd_0lepoverldpratio_2b->SetMarkerStyle(25) ;
      hmctruth_qcd_0lepoverldpratio_3b->SetMarkerStyle(30) ;


      int nBinsMET = (hmctruth_qcd_0lep_1b->GetNbinsX()-1)/(nBinsHT+1) ;
      printf("\n\n Number of MET bins : %d\n\n", nBinsMET ) ;

      TH1F* hscalefactor_qcd_0lepoverldpratio_1b = (TH1F*) hmctruth_qcd_0lepoverldpratio_1b->Clone("hscalefactor_qcd_0lep_1b") ;
      TH1F* hscalefactor_qcd_0lepoverldpratio_2b = (TH1F*) hmctruth_qcd_0lepoverldpratio_2b->Clone("hscalefactor_qcd_0lep_2b") ;
      TH1F* hscalefactor_qcd_0lepoverldpratio_3b = (TH1F*) hmctruth_qcd_0lepoverldpratio_3b->Clone("hscalefactor_qcd_0lep_3b") ;

      for ( int mbi=0; mbi<nBinsMET ; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT ; hbi++ ) {
            int hbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
            if ( R_lsb[0][hbi] > 0. ) {
               printf("  MET,HT bin %d,%d  scale factor : %6.3f / %6.3f\n", mbi+1, hbi+1, hmctruth_qcd_0lepoverldpratio_1b->GetBinContent( hbin ),
                    R_lsb[0][hbi] ) ;
               hscalefactor_qcd_0lepoverldpratio_1b -> SetBinContent( hbin, hmctruth_qcd_0lepoverldpratio_1b->GetBinContent( hbin ) / R_lsb[0][hbi] ) ;
               hscalefactor_qcd_0lepoverldpratio_1b -> SetBinError(   hbin, hmctruth_qcd_0lepoverldpratio_1b->GetBinError( hbin ) / R_lsb[0][hbi] ) ;
            }
            if ( R_lsb[1][hbi] > 0. ) {
               hscalefactor_qcd_0lepoverldpratio_2b -> SetBinContent( hbin, hmctruth_qcd_0lepoverldpratio_2b->GetBinContent( hbin ) / R_lsb[1][hbi] ) ;
               hscalefactor_qcd_0lepoverldpratio_2b -> SetBinError(   hbin, hmctruth_qcd_0lepoverldpratio_2b->GetBinError( hbin ) / R_lsb[1][hbi] ) ;
            }
            if ( R_lsb[2][hbi] > 0. ) {
               hscalefactor_qcd_0lepoverldpratio_3b -> SetBinContent( hbin, hmctruth_qcd_0lepoverldpratio_3b->GetBinContent( hbin ) / R_lsb[2][hbi] ) ;
               hscalefactor_qcd_0lepoverldpratio_3b -> SetBinError(   hbin, hmctruth_qcd_0lepoverldpratio_3b->GetBinError( hbin ) / R_lsb[2][hbi] ) ;
            }
         } // hbi
      } // mbi


      hmctruth_qcd_0lepoverldpratio_1b->SetMinimum(-0.1) ;
      hmctruth_qcd_0lepoverldpratio_2b->SetMinimum(-0.1) ;
      hmctruth_qcd_0lepoverldpratio_3b->SetMinimum(-0.1) ;
      hmctruth_qcd_0lepoverldpratio_1b->SetMaximum(0.5) ;
      hmctruth_qcd_0lepoverldpratio_2b->SetMaximum(0.5) ;
      hmctruth_qcd_0lepoverldpratio_3b->SetMaximum(0.5) ;

      hscalefactor_qcd_0lepoverldpratio_1b->SetMinimum(-2.) ;
      hscalefactor_qcd_0lepoverldpratio_2b->SetMinimum(-2.) ;
      hscalefactor_qcd_0lepoverldpratio_3b->SetMinimum(-2.) ;
      hscalefactor_qcd_0lepoverldpratio_1b->SetMaximum(8.) ;
      hscalefactor_qcd_0lepoverldpratio_2b->SetMaximum(8.) ;
      hscalefactor_qcd_0lepoverldpratio_3b->SetMaximum(8.) ;




      TCanvas* cqcd = (TCanvas*) gDirectory->FindObject("cqcd") ;
      if ( cqcd == 0x0 ) {
         cqcd = new TCanvas("cqcd","qcd closure",700,950) ;
      }
      cqcd->Clear() ;
      cqcd->Divide(1,2) ;

      TLegend* legend_qcd = new TLegend( 0.88, 0.77,  0.97, 0.93 ) ;
      legend_qcd->SetFillColor(kWhite) ;
      legend_qcd->AddEntry( hmctruth_qcd_0lepoverldpratio_1b, "=1b") ;
      legend_qcd->AddEntry( hmctruth_qcd_0lepoverldpratio_2b, "=2b") ;
      legend_qcd->AddEntry( hmctruth_qcd_0lepoverldpratio_3b, ">=3b") ;

      cqcd->cd(1) ;
      hmctruth_qcd_0lepoverldpratio_1b->SetTitle("QCD: 0 Lepton / LDP, Ratio") ;
      hmctruth_qcd_0lepoverldpratio_1b->Draw() ;
      hmctruth_qcd_0lepoverldpratio_2b->Draw("same") ;
      hmctruth_qcd_0lepoverldpratio_3b->Draw("same") ;
      legend_qcd->Draw() ;

      cqcd->cd(2) ;
      hscalefactor_qcd_0lepoverldpratio_1b->SetTitle("QCD: 0 Lepton / LDP, Scale Factor") ;
      hscalefactor_qcd_0lepoverldpratio_1b->Draw() ;
      hscalefactor_qcd_0lepoverldpratio_2b->Draw("same") ;
      hscalefactor_qcd_0lepoverldpratio_3b->Draw("same") ;
      line->DrawLine(0.5,1.,hscalefactor_qcd_0lepoverldpratio_1b->GetNbinsX()+0.5,1.) ;
      legend_qcd->Draw() ;

      TString outqcd( infileStr ) ;
      outqcd.ReplaceAll("rootfiles","outputfiles") ;
      outqcd.ReplaceAll(".root","-mcclosure-qcd1.pdf") ;
      cqcd->SaveAs( outqcd ) ;






      TCanvas* cqcd3 = (TCanvas*) gDirectory->FindObject("cqcd3") ;
      if ( cqcd3 == 0x0 ) {
         cqcd3 = new TCanvas("cqcd3","qcd closure",1200,950) ;
      }
      cqcd3->Clear() ;
      cqcd3->Divide(3,2) ;

      cqcd3->cd(1) ;
      hmctruth_qcd_0lepoverldpratio_1b->SetTitle("QCD =1b: 0 Lepton / LDP, Ratio") ;
      hmctruth_qcd_0lepoverldpratio_1b->Draw() ;
      cqcd3->cd(2) ;
      hmctruth_qcd_0lepoverldpratio_2b->SetTitle("QCD =2b: 0 Lepton / LDP, Ratio") ;
      hmctruth_qcd_0lepoverldpratio_2b->Draw() ;
      cqcd3->cd(3) ;
      hmctruth_qcd_0lepoverldpratio_3b->SetTitle("QCD >=3b: 0 Lepton / LDP, Ratio") ;
      hmctruth_qcd_0lepoverldpratio_3b->Draw() ;

      cqcd3->cd(4) ;
      hscalefactor_qcd_0lepoverldpratio_1b->SetTitle("QCD =1b: 0 Lepton / LDP, Scale Factor") ;
      hscalefactor_qcd_0lepoverldpratio_1b->Draw() ;
      line->DrawLine(0.5,1.,hscalefactor_qcd_0lepoverldpratio_1b->GetNbinsX()+0.5,1.) ;
      hscalefactor_qcd_0lepoverldpratio_1b->Draw("same") ;
      cqcd3->cd(5) ;
      hscalefactor_qcd_0lepoverldpratio_2b->SetTitle("QCD =2b: 0 Lepton / LDP, Scale Factor") ;
      hscalefactor_qcd_0lepoverldpratio_2b->Draw() ;
      line->DrawLine(0.5,1.,hscalefactor_qcd_0lepoverldpratio_1b->GetNbinsX()+0.5,1.) ;
      hscalefactor_qcd_0lepoverldpratio_2b->Draw("same") ;
      cqcd3->cd(6) ;
      hscalefactor_qcd_0lepoverldpratio_3b->SetTitle("QCD >=3b: 0 Lepton / LDP, Scale Factor") ;
      hscalefactor_qcd_0lepoverldpratio_3b->Draw() ;
      line->DrawLine(0.5,1.,hscalefactor_qcd_0lepoverldpratio_1b->GetNbinsX()+0.5,1.) ;
      hscalefactor_qcd_0lepoverldpratio_3b->Draw("same") ;

      TString outqcd3( infileStr ) ;
      outqcd3.ReplaceAll("rootfiles","outputfiles") ;
      outqcd3.ReplaceAll(".root","-mcclosure-qcd3.pdf") ;
      cqcd3->SaveAs( outqcd3 ) ;


   }



//==========================================================================================


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

//==========================================================================================

