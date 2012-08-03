
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

#include "updateFileValue.c"

#include <iostream>

  using std::cout ;
  using std::endl ;

   void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;

  //----------------------------------------

   void mcclosure2( const char* infile = "rootfiles/gi-plots-met4-ht4.root", const char* datfile = "Input-met4-ht4-wsyst1.dat" ) {

      int nBinsMET(0), nBinsHT(0) ;

      TString infileStr( infile ) ;

      gSystem->Exec("mkdir -p outputfiles") ;

      gStyle->SetOptStat(0) ;
      gStyle->SetPadTopMargin(0.03) ;
      gStyle->SetPadBottomMargin(0.30) ;
      gStyle->SetPadRightMargin(0.05) ;
      gStyle->SetPadLeftMargin(0.15) ;
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


      char binlabel[1000] ;
      sprintf( binlabel, "%s", hmctruth_ttwj_0lep_1b -> GetXaxis() -> GetBinLabel( hmctruth_ttwj_0lep_1b->GetNbinsX() - 1 ) ) ;
      sscanf( binlabel, "0lep_M%d_H%d_1b", &nBinsMET, &nBinsHT ) ;
      printf("\n\n Bin label: %s,  nmet=%d, nht=%d\n\n", binlabel, nBinsMET, nBinsHT ) ;




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
      hmctruth_ttwj_0over1ratio_1b->SetMaximum(3.) ;
      hmctruth_ttwj_0over1ratio_2b->SetMaximum(3.) ;
      hmctruth_ttwj_0over1ratio_3b->SetMaximum(3.) ;

      hscalefactor_ttwj_0over1ratio_1b->SetMinimum(0.) ;
      hscalefactor_ttwj_0over1ratio_2b->SetMinimum(0.) ;
      hscalefactor_ttwj_0over1ratio_3b->SetMinimum(0.) ;
      hscalefactor_ttwj_0over1ratio_1b->SetMaximum(2.5) ;
      hscalefactor_ttwj_0over1ratio_2b->SetMaximum(2.5) ;
      hscalefactor_ttwj_0over1ratio_3b->SetMaximum(2.5) ;



      gStyle->SetPadRightMargin(0.08) ;
      TCanvas* cttwj = (TCanvas*) gDirectory->FindObject("cttwj") ;
      if ( cttwj == 0x0 ) {
         cttwj = new TCanvas("cttwj","ttwj closure", 700, 950) ;
      }
      cttwj->Clear() ;
      cttwj->Divide(1,2) ;

      TLegend* legend_ttwj = new TLegend( 0.91, 0.77,  0.99, 0.93 ) ;
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




   //----


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

   //----
   

      hmctruth_ttwj_0lep_1b->SetFillColor(kBlue-9) ;
      hmctruth_ttwj_0lep_2b->SetFillColor(kBlue-9) ;
      hmctruth_ttwj_0lep_3b->SetFillColor(kBlue-9) ;

      hmctruth_ttwj_1lep_1b->SetFillColor(kBlue-9) ;
      hmctruth_ttwj_1lep_2b->SetFillColor(kBlue-9) ;
      hmctruth_ttwj_1lep_3b->SetFillColor(kBlue-9) ;

      hmctruth_ttwj_0lep_1b->SetLineWidth(2) ;
      hmctruth_ttwj_0lep_2b->SetLineWidth(2) ;
      hmctruth_ttwj_0lep_3b->SetLineWidth(2) ;

      hmctruth_ttwj_1lep_1b->SetLineWidth(2) ;
      hmctruth_ttwj_1lep_2b->SetLineWidth(2) ;
      hmctruth_ttwj_1lep_3b->SetLineWidth(2) ;

      TCanvas* cttwj3b = (TCanvas*) gDirectory->FindObject("cttwj3b") ;
      if ( cttwj3b == 0x0 ) {
         cttwj3b = new TCanvas("cttwj3b","ttwj closure", 1200, 950) ;
      }
      cttwj3b->Clear() ;
      cttwj3b->Divide(3,2) ;

      cttwj3b->cd(1) ;
      hmctruth_ttwj_0lep_1b->SetTitle("ttwj =1b: 0 Lepton") ;
      hmctruth_ttwj_0lep_1b->Draw() ;
      hmctruth_ttwj_0lep_1b->Draw("histsame") ;
      hmctruth_ttwj_0lep_1b->Draw("esame") ;
      cttwj3b->cd(2) ;
      hmctruth_ttwj_0lep_2b->SetTitle("ttwj =2b: 0 Lepton") ;
      hmctruth_ttwj_0lep_2b->Draw() ;
      hmctruth_ttwj_0lep_2b->Draw("histsame") ;
      hmctruth_ttwj_0lep_2b->Draw("esame") ;
      cttwj3b->cd(3) ;
      hmctruth_ttwj_0lep_3b->SetTitle("ttwj >=3b: 0 Lepton") ;
      hmctruth_ttwj_0lep_3b->Draw() ;
      hmctruth_ttwj_0lep_3b->Draw("histsame") ;
      hmctruth_ttwj_0lep_3b->Draw("esame") ;

      cttwj3b->cd(4) ;
      hmctruth_ttwj_1lep_1b->SetTitle("ttwj =1b: 1 Lepton") ;
      hmctruth_ttwj_1lep_1b->Draw() ;
      hmctruth_ttwj_1lep_1b->Draw("histsame") ;
      hmctruth_ttwj_1lep_1b->Draw("esame") ;
      cttwj3b->cd(5) ;
      hmctruth_ttwj_1lep_2b->SetTitle("ttwj =2b: 1 Lepton") ;
      hmctruth_ttwj_1lep_2b->Draw() ;
      hmctruth_ttwj_1lep_2b->Draw("histsame") ;
      hmctruth_ttwj_1lep_2b->Draw("esame") ;
      cttwj3b->cd(6) ;
      hmctruth_ttwj_1lep_3b->SetTitle("ttwj >=3b: 1 Lepton") ;
      hmctruth_ttwj_1lep_3b->Draw() ;
      hmctruth_ttwj_1lep_3b->Draw("histsame") ;
      hmctruth_ttwj_1lep_3b->Draw("esame") ;



      TString outttwj3b( infileStr ) ;
      outttwj3b.ReplaceAll("rootfiles","outputfiles") ;
      outttwj3b.ReplaceAll(".root","-mcclosure-ttwj3b.pdf") ;
      cttwj3b->SaveAs( outttwj3b ) ;

   //----



    //--- insert ttwj numbers into file.

      for ( int mbi=0; mbi<nBinsMET ; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT ; hbi++ ) {

            int hbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

            char parameterName[1000] ;
            double err, diff, systValue, correction ;


            sprintf( parameterName, "sf_ttwj_M%d_H%d_1b_err", mbi+1, hbi+1 ) ;
            printf( "met=%d, ht=%d : %s %s\n", mbi+1, hbi+1,
                  parameterName,
                  hmctruth_ttwj_0lep_1b -> GetXaxis() -> GetBinLabel( hbin ) ) ;
            err  = hscalefactor_ttwj_0over1ratio_1b->GetBinError(hbin)  ;
            correction = hscalefactor_ttwj_0over1ratio_1b->GetBinContent(hbin)  ;
            diff = correction - 1. ;
            if ( err > 0. ) {
               systValue = sqrt( pow(err,2) + pow(0.5*diff,2) ) ;
            } else {
               systValue = 3.0 ;
            }
            updateFileValue( datfile, parameterName, systValue ) ;
            sprintf( parameterName, "sf_ttwj_M%d_H%d_1b", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, correction ) ;

            sprintf( parameterName, "sf_ttwj_M%d_H%d_2b_err", mbi+1, hbi+1 ) ;
            err  = hscalefactor_ttwj_0over1ratio_2b->GetBinError(hbin)  ;
            correction = hscalefactor_ttwj_0over1ratio_2b->GetBinContent(hbin) ;
            diff = correction - 1. ;
            if ( err > 0. ) {
               systValue = sqrt( pow(err,2) + pow(0.5*diff,2) ) ;
            } else {
               systValue = 3.0 ;
            }
            updateFileValue( datfile, parameterName, systValue ) ;
            sprintf( parameterName, "sf_ttwj_M%d_H%d_2b", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, correction ) ;

            sprintf( parameterName, "sf_ttwj_M%d_H%d_3b_err", mbi+1, hbi+1 ) ;
            err  = hscalefactor_ttwj_0over1ratio_3b->GetBinError(hbin)  ;
            correction = hscalefactor_ttwj_0over1ratio_3b->GetBinContent(hbin)  ;
            diff = correction - 1. ;
            if ( err > 0. ) {
               systValue = sqrt( pow(err,2) + pow(0.5*diff,2) ) ;
            } else {
               systValue = 3.0 ;
            }
            updateFileValue( datfile, parameterName, systValue ) ;
            sprintf( parameterName, "sf_ttwj_M%d_H%d_3b", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, correction ) ;

         } // hbi.
      } // mbi.



















     //---   Q C D  Part   --------------






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

      TH1F* hmctruth_qcd_0lepoverldpratio_1b_htgrouping = (TH1F*) hmctruth_qcd_0lepoverldpratio_1b->Clone("hmctruth_qcd_0lepoverldpratio_1b_htgrouping") ;
      TH1F* hmctruth_qcd_0lepoverldpratio_2b_htgrouping = (TH1F*) hmctruth_qcd_0lepoverldpratio_2b->Clone("hmctruth_qcd_0lepoverldpratio_2b_htgrouping") ;
      TH1F* hmctruth_qcd_0lepoverldpratio_3b_htgrouping = (TH1F*) hmctruth_qcd_0lepoverldpratio_3b->Clone("hmctruth_qcd_0lepoverldpratio_3b_htgrouping") ;

      TH1F* hmctruth_qcd_0lep_1b_htgrouping = (TH1F*) hmctruth_qcd_0lep_1b->Clone("hmctruth_qcd_0lep_1b_htgrouping") ;
      TH1F* hmctruth_qcd_0lep_2b_htgrouping = (TH1F*) hmctruth_qcd_0lep_2b->Clone("hmctruth_qcd_0lep_2b_htgrouping") ;
      TH1F* hmctruth_qcd_0lep_3b_htgrouping = (TH1F*) hmctruth_qcd_0lep_3b->Clone("hmctruth_qcd_0lep_3b_htgrouping") ;

      TH1F* hmctruth_qcd_ldp_1b_htgrouping = (TH1F*) hmctruth_qcd_ldp_1b->Clone("hmctruth_qcd_ldp_1b_htgrouping") ;
      TH1F* hmctruth_qcd_ldp_2b_htgrouping = (TH1F*) hmctruth_qcd_ldp_2b->Clone("hmctruth_qcd_ldp_2b_htgrouping") ;
      TH1F* hmctruth_qcd_ldp_3b_htgrouping = (TH1F*) hmctruth_qcd_ldp_3b->Clone("hmctruth_qcd_ldp_3b_htgrouping") ;









      double zlsum[100] ;
      double ldpsum[100] ;

      double zlsumw2[100] ;
      double ldpsumw2[100] ;

      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         zlsum[hbi] = 0. ;
         zlsumw2[hbi] = 0. ;
         ldpsum[hbi] = 0. ;
         ldpsumw2[hbi] = 0. ;
      }

      for ( int mbi=0; mbi<nBinsMET ; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT ; hbi++ ) {

            int hbinin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

            zlsum[hbi] += hmctruth_qcd_0lep_1b->GetBinContent( hbinin ) ;
            zlsum[hbi] += hmctruth_qcd_0lep_2b->GetBinContent( hbinin ) ;
            zlsum[hbi] += hmctruth_qcd_0lep_3b->GetBinContent( hbinin ) ;

            zlsumw2[hbi] += pow( hmctruth_qcd_0lep_1b->GetBinError( hbinin ), 2 ) ;
            zlsumw2[hbi] += pow( hmctruth_qcd_0lep_2b->GetBinError( hbinin ), 2 ) ;
            zlsumw2[hbi] += pow( hmctruth_qcd_0lep_3b->GetBinError( hbinin ), 2 ) ;


            ldpsum[hbi] += hmctruth_qcd_ldp_1b->GetBinContent( hbinin ) ;
            ldpsum[hbi] += hmctruth_qcd_ldp_2b->GetBinContent( hbinin ) ;
            ldpsum[hbi] += hmctruth_qcd_ldp_3b->GetBinContent( hbinin ) ;

            ldpsumw2[hbi] += pow( hmctruth_qcd_ldp_1b->GetBinError( hbinin ), 2 ) ;
            ldpsumw2[hbi] += pow( hmctruth_qcd_ldp_2b->GetBinError( hbinin ), 2 ) ;
            ldpsumw2[hbi] += pow( hmctruth_qcd_ldp_3b->GetBinError( hbinin ), 2 ) ;


         } // hbi
      } // mbi




      printf("\n\n") ;

      double qcd_0lepoverldpratio[100] ;
      double qcd_0lepoverldpratio_err[100] ;

      for ( int hbi=0; hbi<nBinsHT ; hbi++ ) {

          qcd_0lepoverldpratio[hbi] = 0. ;
          if ( ldpsum[hbi] > 0. ) {
             qcd_0lepoverldpratio[hbi] =  zlsum[hbi] / ldpsum[hbi] ;
          }
          qcd_0lepoverldpratio_err[hbi] = 0. ;
          if ( ldpsum[hbi] > 0. && zlsum[hbi] > 0. ) {
             qcd_0lepoverldpratio_err[hbi] = qcd_0lepoverldpratio[hbi] * sqrt(
                    zlsumw2[hbi] / pow( zlsum[hbi],2)
                  + ldpsumw2[hbi] / pow( ldpsum[hbi],2)
                    ) ;
          }

          printf(" HT bin %d : 0lep/LDP ratio = %5.3f +/- %5.3f\n", hbi+1, qcd_0lepoverldpratio[hbi], qcd_0lepoverldpratio_err[hbi] ) ;

      } // hbi.

      printf("\n\n") ;



      for ( int mbi=0; mbi<nBinsMET ; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT ; hbi++ ) {

            int hbinin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
            int hbinout = 1 + (nBinsMET+1)*hbi + mbi + 1 ;

            hmctruth_qcd_0lepoverldpratio_1b_htgrouping->SetBinContent( hbinout, hmctruth_qcd_0lepoverldpratio_1b->GetBinContent( hbinin ) ) ;
            hmctruth_qcd_0lepoverldpratio_1b_htgrouping->SetBinError(   hbinout, hmctruth_qcd_0lepoverldpratio_1b->GetBinError(   hbinin ) ) ;
            hmctruth_qcd_0lepoverldpratio_1b_htgrouping->GetXaxis()->SetBinLabel( hbinout,
                 hmctruth_qcd_0lepoverldpratio_1b->GetXaxis()->GetBinLabel( hbinin ) ) ;

            hmctruth_qcd_0lepoverldpratio_2b_htgrouping->SetBinContent( hbinout, hmctruth_qcd_0lepoverldpratio_2b->GetBinContent( hbinin ) ) ;
            hmctruth_qcd_0lepoverldpratio_2b_htgrouping->SetBinError(   hbinout, hmctruth_qcd_0lepoverldpratio_2b->GetBinError(   hbinin ) ) ;
            hmctruth_qcd_0lepoverldpratio_2b_htgrouping->GetXaxis()->SetBinLabel( hbinout,
                 hmctruth_qcd_0lepoverldpratio_2b->GetXaxis()->GetBinLabel( hbinin ) ) ;

            hmctruth_qcd_0lepoverldpratio_3b_htgrouping->SetBinContent( hbinout, hmctruth_qcd_0lepoverldpratio_3b->GetBinContent( hbinin ) ) ;
            hmctruth_qcd_0lepoverldpratio_3b_htgrouping->SetBinError(   hbinout, hmctruth_qcd_0lepoverldpratio_3b->GetBinError(   hbinin ) ) ;
            hmctruth_qcd_0lepoverldpratio_3b_htgrouping->GetXaxis()->SetBinLabel( hbinout,
                 hmctruth_qcd_0lepoverldpratio_3b->GetXaxis()->GetBinLabel( hbinin ) ) ;



            hmctruth_qcd_0lep_1b_htgrouping->SetBinContent( hbinout, hmctruth_qcd_0lep_1b->GetBinContent( hbinin ) ) ;
            hmctruth_qcd_0lep_1b_htgrouping->SetBinError(   hbinout, hmctruth_qcd_0lep_1b->GetBinError(   hbinin ) ) ;
            hmctruth_qcd_0lep_1b_htgrouping->GetXaxis()->SetBinLabel( hbinout,
                 hmctruth_qcd_0lep_1b->GetXaxis()->GetBinLabel( hbinin ) ) ;

            hmctruth_qcd_0lep_2b_htgrouping->SetBinContent( hbinout, hmctruth_qcd_0lep_2b->GetBinContent( hbinin ) ) ;
            hmctruth_qcd_0lep_2b_htgrouping->SetBinError(   hbinout, hmctruth_qcd_0lep_2b->GetBinError(   hbinin ) ) ;
            hmctruth_qcd_0lep_2b_htgrouping->GetXaxis()->SetBinLabel( hbinout,
                 hmctruth_qcd_0lep_2b->GetXaxis()->GetBinLabel( hbinin ) ) ;

            hmctruth_qcd_0lep_3b_htgrouping->SetBinContent( hbinout, hmctruth_qcd_0lep_3b->GetBinContent( hbinin ) ) ;
            hmctruth_qcd_0lep_3b_htgrouping->SetBinError(   hbinout, hmctruth_qcd_0lep_3b->GetBinError(   hbinin ) ) ;
            hmctruth_qcd_0lep_3b_htgrouping->GetXaxis()->SetBinLabel( hbinout,
                 hmctruth_qcd_0lep_3b->GetXaxis()->GetBinLabel( hbinin ) ) ;



            hmctruth_qcd_ldp_1b_htgrouping->SetBinContent( hbinout, hmctruth_qcd_ldp_1b->GetBinContent( hbinin ) ) ;
            hmctruth_qcd_ldp_1b_htgrouping->SetBinError(   hbinout, hmctruth_qcd_ldp_1b->GetBinError(   hbinin ) ) ;
            hmctruth_qcd_ldp_1b_htgrouping->GetXaxis()->SetBinLabel( hbinout,
                 hmctruth_qcd_ldp_1b->GetXaxis()->GetBinLabel( hbinin ) ) ;

            hmctruth_qcd_ldp_2b_htgrouping->SetBinContent( hbinout, hmctruth_qcd_ldp_2b->GetBinContent( hbinin ) ) ;
            hmctruth_qcd_ldp_2b_htgrouping->SetBinError(   hbinout, hmctruth_qcd_ldp_2b->GetBinError(   hbinin ) ) ;
            hmctruth_qcd_ldp_2b_htgrouping->GetXaxis()->SetBinLabel( hbinout,
                 hmctruth_qcd_ldp_2b->GetXaxis()->GetBinLabel( hbinin ) ) ;

            hmctruth_qcd_ldp_3b_htgrouping->SetBinContent( hbinout, hmctruth_qcd_ldp_3b->GetBinContent( hbinin ) ) ;
            hmctruth_qcd_ldp_3b_htgrouping->SetBinError(   hbinout, hmctruth_qcd_ldp_3b->GetBinError(   hbinin ) ) ;
            hmctruth_qcd_ldp_3b_htgrouping->GetXaxis()->SetBinLabel( hbinout,
                 hmctruth_qcd_ldp_3b->GetXaxis()->GetBinLabel( hbinin ) ) ;



         } // hbi
      } // mbi








      TH1F* hscalefactor_qcd_0lepoverldpratio_1b_htgrouping = (TH1F*) hmctruth_qcd_0lepoverldpratio_1b_htgrouping->Clone("hscalefactor_qcd_0lep_1b_htgrouping") ;
      TH1F* hscalefactor_qcd_0lepoverldpratio_2b_htgrouping = (TH1F*) hmctruth_qcd_0lepoverldpratio_2b_htgrouping->Clone("hscalefactor_qcd_0lep_2b_htgrouping") ;
      TH1F* hscalefactor_qcd_0lepoverldpratio_3b_htgrouping = (TH1F*) hmctruth_qcd_0lepoverldpratio_3b_htgrouping->Clone("hscalefactor_qcd_0lep_3b_htgrouping") ;


      for ( int mbi=0; mbi<nBinsMET ; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT ; hbi++ ) {

            int hbin = 1 + (nBinsMET+1)*hbi + mbi + 1 ;

            double sf ;
            if ( qcd_0lepoverldpratio[hbi] > 0. ) {

               sf = hmctruth_qcd_0lepoverldpratio_1b_htgrouping->GetBinContent( hbin ) / qcd_0lepoverldpratio[hbi] ;
               printf("  1b MET,HT bin %d,%d  scale factor : %6.3f / %6.3f = %6.3f\n", mbi+1, hbi+1, hmctruth_qcd_0lepoverldpratio_1b_htgrouping->GetBinContent( hbin ),
                    qcd_0lepoverldpratio[hbi], sf ) ;
               hscalefactor_qcd_0lepoverldpratio_1b_htgrouping -> SetBinContent( hbin, hmctruth_qcd_0lepoverldpratio_1b_htgrouping->GetBinContent( hbin ) /
                                                                                       qcd_0lepoverldpratio[hbi] ) ;
               hscalefactor_qcd_0lepoverldpratio_1b_htgrouping -> SetBinError(   hbin, hmctruth_qcd_0lepoverldpratio_1b_htgrouping->GetBinError( hbin ) /
                                                                                       qcd_0lepoverldpratio[hbi] ) ;

               sf = hmctruth_qcd_0lepoverldpratio_2b_htgrouping->GetBinContent( hbin ) / qcd_0lepoverldpratio[hbi] ;
               printf("  2b MET,HT bin %d,%d  scale factor : %6.3f / %6.3f = %6.3f\n", mbi+1, hbi+1, hmctruth_qcd_0lepoverldpratio_2b_htgrouping->GetBinContent( hbin ),
                    qcd_0lepoverldpratio[hbi], sf ) ;
               hscalefactor_qcd_0lepoverldpratio_2b_htgrouping -> SetBinContent( hbin, hmctruth_qcd_0lepoverldpratio_2b_htgrouping->GetBinContent( hbin ) /
                                                                                       qcd_0lepoverldpratio[hbi] ) ;
               hscalefactor_qcd_0lepoverldpratio_2b_htgrouping -> SetBinError(   hbin, hmctruth_qcd_0lepoverldpratio_2b_htgrouping->GetBinError( hbin ) /
                                                                                       qcd_0lepoverldpratio[hbi] ) ;

               sf = hmctruth_qcd_0lepoverldpratio_3b_htgrouping->GetBinContent( hbin ) / qcd_0lepoverldpratio[hbi] ;
               printf("  3b MET,HT bin %d,%d  scale factor : %6.3f / %6.3f = %6.3f\n", mbi+1, hbi+1, hmctruth_qcd_0lepoverldpratio_3b_htgrouping->GetBinContent( hbin ),
                    qcd_0lepoverldpratio[hbi], sf ) ;
               hscalefactor_qcd_0lepoverldpratio_3b_htgrouping -> SetBinContent( hbin, hmctruth_qcd_0lepoverldpratio_3b_htgrouping->GetBinContent( hbin ) /
                                                                                       qcd_0lepoverldpratio[hbi] ) ;
               hscalefactor_qcd_0lepoverldpratio_3b_htgrouping -> SetBinError(   hbin, hmctruth_qcd_0lepoverldpratio_3b_htgrouping->GetBinError( hbin ) /
                                                                                       qcd_0lepoverldpratio[hbi] ) ;

            }

         } // hbi
      } // mbi


      hmctruth_qcd_0lepoverldpratio_1b->SetMinimum(-0.1) ;
      hmctruth_qcd_0lepoverldpratio_2b->SetMinimum(-0.1) ;
      hmctruth_qcd_0lepoverldpratio_3b->SetMinimum(-0.1) ;
      hmctruth_qcd_0lepoverldpratio_1b->SetMaximum(0.5) ;
      hmctruth_qcd_0lepoverldpratio_2b->SetMaximum(0.5) ;
      hmctruth_qcd_0lepoverldpratio_3b->SetMaximum(0.5) ;

      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->SetMinimum(-0.1) ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->SetMinimum(-0.1) ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->SetMinimum(-0.1) ;
      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->SetMaximum(0.5) ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->SetMaximum(0.5) ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->SetMaximum(0.5) ;

      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->SetMinimum(-2.) ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->SetMinimum(-2.) ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->SetMinimum(-2.) ;
      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->SetMaximum(8.) ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->SetMaximum(8.) ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->SetMaximum(8.) ;

      hmctruth_qcd_0lep_1b_htgrouping->GetYaxis()->SetNdivisions(505) ;
      hmctruth_qcd_0lep_2b_htgrouping->GetYaxis()->SetNdivisions(505) ;
      hmctruth_qcd_0lep_3b_htgrouping->GetYaxis()->SetNdivisions(505) ;

      hmctruth_qcd_ldp_1b_htgrouping->GetYaxis()->SetNdivisions(505) ;
      hmctruth_qcd_ldp_2b_htgrouping->GetYaxis()->SetNdivisions(505) ;
      hmctruth_qcd_ldp_3b_htgrouping->GetYaxis()->SetNdivisions(505) ;

      hmctruth_qcd_0lep_1b_htgrouping->GetYaxis()->SetLabelSize(0.07) ;
      hmctruth_qcd_0lep_2b_htgrouping->GetYaxis()->SetLabelSize(0.07) ;
      hmctruth_qcd_0lep_3b_htgrouping->GetYaxis()->SetLabelSize(0.07) ;

      hmctruth_qcd_ldp_1b_htgrouping->GetYaxis()->SetLabelSize(0.07) ;
      hmctruth_qcd_ldp_2b_htgrouping->GetYaxis()->SetLabelSize(0.07) ;
      hmctruth_qcd_ldp_3b_htgrouping->GetYaxis()->SetLabelSize(0.07) ;

      hmctruth_qcd_0lep_1b_htgrouping->GetXaxis()->SetLabelSize(0.06) ;
      hmctruth_qcd_0lep_2b_htgrouping->GetXaxis()->SetLabelSize(0.06) ;
      hmctruth_qcd_0lep_3b_htgrouping->GetXaxis()->SetLabelSize(0.06) ;

      hmctruth_qcd_ldp_1b_htgrouping->GetXaxis()->SetLabelSize(0.06) ;
      hmctruth_qcd_ldp_2b_htgrouping->GetXaxis()->SetLabelSize(0.06) ;
      hmctruth_qcd_ldp_3b_htgrouping->GetXaxis()->SetLabelSize(0.06) ;

      hmctruth_qcd_0lep_1b_htgrouping->GetYaxis()->SetLabelOffset(0.01) ;
      hmctruth_qcd_0lep_2b_htgrouping->GetYaxis()->SetLabelOffset(0.01) ;
      hmctruth_qcd_0lep_3b_htgrouping->GetYaxis()->SetLabelOffset(0.01) ;

      hmctruth_qcd_ldp_1b_htgrouping->GetYaxis()->SetLabelOffset(0.01) ;
      hmctruth_qcd_ldp_2b_htgrouping->GetYaxis()->SetLabelOffset(0.01) ;
      hmctruth_qcd_ldp_3b_htgrouping->GetYaxis()->SetLabelOffset(0.01) ;

      hmctruth_qcd_0lep_1b_htgrouping->SetFillColor(2) ;
      hmctruth_qcd_0lep_2b_htgrouping->SetFillColor(2) ;
      hmctruth_qcd_0lep_3b_htgrouping->SetFillColor(2) ;

      hmctruth_qcd_ldp_1b_htgrouping->SetFillColor(2) ;
      hmctruth_qcd_ldp_2b_htgrouping->SetFillColor(2) ;
      hmctruth_qcd_ldp_3b_htgrouping->SetFillColor(2) ;


      hmctruth_qcd_0lep_1b_htgrouping->SetLineWidth(2) ;
      hmctruth_qcd_0lep_2b_htgrouping->SetLineWidth(2) ;
      hmctruth_qcd_0lep_3b_htgrouping->SetLineWidth(2) ;

      hmctruth_qcd_ldp_1b_htgrouping->SetLineWidth(2) ;
      hmctruth_qcd_ldp_2b_htgrouping->SetLineWidth(2) ;
      hmctruth_qcd_ldp_3b_htgrouping->SetLineWidth(2) ;


      hmctruth_qcd_ldp_1b_htgrouping->SetMaximum( 1.15*hmctruth_qcd_ldp_1b_htgrouping->GetMaximum() ) ;
      hmctruth_qcd_ldp_2b_htgrouping->SetMaximum( 1.15*hmctruth_qcd_ldp_2b_htgrouping->GetMaximum() ) ;
      hmctruth_qcd_ldp_3b_htgrouping->SetMaximum( 1.15*hmctruth_qcd_ldp_3b_htgrouping->GetMaximum() ) ;

      hmctruth_qcd_0lep_1b_htgrouping->SetMaximum( 0.10*hmctruth_qcd_ldp_1b_htgrouping->GetMaximum() ) ;
      hmctruth_qcd_0lep_2b_htgrouping->SetMaximum( 0.10*hmctruth_qcd_ldp_2b_htgrouping->GetMaximum() ) ;
      hmctruth_qcd_0lep_3b_htgrouping->SetMaximum( 0.10*hmctruth_qcd_ldp_3b_htgrouping->GetMaximum() ) ;

  /// double qcdldpmax(700) ;
  /// hmctruth_qcd_ldp_1b_htgrouping->SetMaximum(qcdldpmax) ;
  /// hmctruth_qcd_ldp_2b_htgrouping->SetMaximum(qcdldpmax) ;
  /// hmctruth_qcd_ldp_3b_htgrouping->SetMaximum(qcdldpmax) ;




    //---------

      TCanvas* cqcd3 = (TCanvas*) gDirectory->FindObject("cqcd3") ;
      if ( cqcd3 == 0x0 ) {
         cqcd3 = new TCanvas("cqcd3","qcd closure",1200,950) ;
      }
      cqcd3->Clear() ;
      cqcd3->Divide(3,2) ;

      cqcd3->cd(1) ;
      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->SetTitle("QCD =1b: 0 Lepton / LDP, Ratio") ;
      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->DrawCopy() ;
      cqcd3->cd(2) ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->SetTitle("QCD =2b: 0 Lepton / LDP, Ratio") ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->DrawCopy() ;
      cqcd3->cd(3) ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->SetTitle("QCD >=3b: 0 Lepton / LDP, Ratio") ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->DrawCopy() ;

      cqcd3->cd(4) ;
      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->SetTitle("QCD =1b: 0 Lepton / LDP, Scale Factor") ;
      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->DrawCopy() ;
      line->DrawLine(0.5,1.,hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->GetNbinsX()+0.5,1.) ;
      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->Draw("same") ;
      cqcd3->cd(5) ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->SetTitle("QCD =2b: 0 Lepton / LDP, Scale Factor") ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->DrawCopy() ;
      line->DrawLine(0.5,1.,hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->GetNbinsX()+0.5,1.) ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->Draw("same") ;
      cqcd3->cd(6) ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->SetTitle("QCD >=3b: 0 Lepton / LDP, Scale Factor") ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->DrawCopy() ;
      line->DrawLine(0.5,1.,hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->GetNbinsX()+0.5,1.) ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->Draw("same") ;

      TString outqcd3( infileStr ) ;
      outqcd3.ReplaceAll("rootfiles","outputfiles") ;
      outqcd3.ReplaceAll(".root","-mcclosure-qcd3.pdf") ;
      cqcd3->SaveAs( outqcd3 ) ;

    //---------



      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->SetLineWidth(2) ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->SetLineWidth(2) ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->SetLineWidth(2) ;

      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->GetYaxis()->SetNdivisions(505) ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->GetYaxis()->SetNdivisions(505) ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->GetYaxis()->SetNdivisions(505) ;

      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->GetYaxis()->SetLabelSize(0.07) ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->GetYaxis()->SetLabelSize(0.07) ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->GetYaxis()->SetLabelSize(0.07) ;

      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->GetYaxis()->SetLabelOffset(0.01) ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->GetYaxis()->SetLabelOffset(0.01) ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->GetYaxis()->SetLabelOffset(0.01) ;


      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->SetLineWidth(2) ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->SetLineWidth(2) ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->SetLineWidth(2) ;

      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->GetYaxis()->SetNdivisions(505) ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->GetYaxis()->SetNdivisions(505) ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->GetYaxis()->SetNdivisions(505) ;

      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->GetYaxis()->SetLabelSize(0.07) ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->GetYaxis()->SetLabelSize(0.07) ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->GetYaxis()->SetLabelSize(0.07) ;

      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->GetYaxis()->SetLabelOffset(0.01) ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->GetYaxis()->SetLabelOffset(0.01) ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->GetYaxis()->SetLabelOffset(0.01) ;

      hmctruth_qcd_0lep_1b_htgrouping->SetYTitle("Events") ;
      hmctruth_qcd_0lep_2b_htgrouping->SetYTitle("Events") ;
      hmctruth_qcd_0lep_3b_htgrouping->SetYTitle("Events") ;

      hmctruth_qcd_ldp_1b_htgrouping->SetYTitle("Events") ;
      hmctruth_qcd_ldp_2b_htgrouping->SetYTitle("Events") ;
      hmctruth_qcd_ldp_3b_htgrouping->SetYTitle("Events") ;

      hmctruth_qcd_0lep_1b_htgrouping->SetTitleOffset(1.11,"y") ;
      hmctruth_qcd_0lep_2b_htgrouping->SetTitleOffset(1.11,"y") ;
      hmctruth_qcd_0lep_3b_htgrouping->SetTitleOffset(1.11,"y") ;

      hmctruth_qcd_0lep_1b_htgrouping->SetTitleSize(0.07,"y") ;
      hmctruth_qcd_0lep_2b_htgrouping->SetTitleSize(0.07,"y") ;
      hmctruth_qcd_0lep_3b_htgrouping->SetTitleSize(0.07,"y") ;

      hmctruth_qcd_ldp_1b_htgrouping->SetTitleOffset(1.11,"y") ;
      hmctruth_qcd_ldp_2b_htgrouping->SetTitleOffset(1.11,"y") ;
      hmctruth_qcd_ldp_3b_htgrouping->SetTitleOffset(1.11,"y") ;

      hmctruth_qcd_ldp_1b_htgrouping->SetTitleSize(0.07,"y") ;
      hmctruth_qcd_ldp_2b_htgrouping->SetTitleSize(0.07,"y") ;
      hmctruth_qcd_ldp_3b_htgrouping->SetTitleSize(0.07,"y") ;

      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->SetYTitle("0Lep / LDP") ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->SetYTitle("0Lep / LDP") ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->SetYTitle("0Lep / LDP") ;

      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->SetTitleOffset(1.11,"y") ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->SetTitleOffset(1.11,"y") ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->SetTitleOffset(1.11,"y") ;

      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->SetTitleSize(0.07,"y") ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->SetTitleSize(0.07,"y") ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->SetTitleSize(0.07,"y") ;


      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->SetYTitle("Closure") ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->SetYTitle("Closure") ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->SetYTitle("Closure") ;

      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->SetTitleOffset(1.11,"y") ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->SetTitleOffset(1.11,"y") ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->SetTitleOffset(1.11,"y") ;

      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->SetTitleSize(0.07,"y") ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->SetTitleSize(0.07,"y") ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->SetTitleSize(0.07,"y") ;

      TCanvas* cqcd3b = (TCanvas*) gDirectory->FindObject("cqcd3b") ;
      if ( cqcd3b == 0x0 ) {
         cqcd3b = new TCanvas("cqcd3b","qcd closure",1000,1000) ;
      }
      cqcd3b->Clear() ;
      cqcd3b->Divide(3,4) ;

      int pad(1) ;

      cqcd3b->cd(pad++) ;
      hmctruth_qcd_0lep_1b_htgrouping->SetTitle("QCD =1b: 0 Lepton") ;
      hmctruth_qcd_0lep_1b_htgrouping->Draw() ;
      hmctruth_qcd_0lep_1b_htgrouping->Draw("histsame") ;
      hmctruth_qcd_0lep_1b_htgrouping->Draw("esame") ;
      cqcd3b->cd(pad++) ;
      hmctruth_qcd_0lep_2b_htgrouping->SetTitle("QCD =2b: 0 Lepton") ;
      hmctruth_qcd_0lep_2b_htgrouping->Draw() ;
      hmctruth_qcd_0lep_2b_htgrouping->Draw("histsame") ;
      hmctruth_qcd_0lep_2b_htgrouping->Draw("esame") ;
      cqcd3b->cd(pad++) ;
      hmctruth_qcd_0lep_3b_htgrouping->SetTitle("QCD >=3b: 0 Lepton") ;
      hmctruth_qcd_0lep_3b_htgrouping->Draw() ;
      hmctruth_qcd_0lep_3b_htgrouping->Draw("histsame") ;
      hmctruth_qcd_0lep_3b_htgrouping->Draw("esame") ;

      cqcd3b->cd(pad++) ;
      hmctruth_qcd_ldp_1b_htgrouping->SetTitle("QCD =1b: Low minDPhiN") ;
      hmctruth_qcd_ldp_1b_htgrouping->Draw() ;
      hmctruth_qcd_ldp_1b_htgrouping->Draw("histsame") ;
      hmctruth_qcd_ldp_1b_htgrouping->Draw("esame") ;
      cqcd3b->cd(pad++) ;
      hmctruth_qcd_ldp_2b_htgrouping->SetTitle("QCD =2b: Low minDPhiN") ;
      hmctruth_qcd_ldp_2b_htgrouping->Draw() ;
      hmctruth_qcd_ldp_2b_htgrouping->Draw("histsame") ;
      hmctruth_qcd_ldp_2b_htgrouping->Draw("esame") ;
      cqcd3b->cd(pad++) ;
      hmctruth_qcd_ldp_3b_htgrouping->SetTitle("QCD >=3b: Low minDPhiN") ;
      hmctruth_qcd_ldp_3b_htgrouping->Draw() ;
      hmctruth_qcd_ldp_3b_htgrouping->Draw("histsame") ;
      hmctruth_qcd_ldp_3b_htgrouping->Draw("esame") ;

      cqcd3b->cd(pad++) ;
      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->SetTitle("QCD =1b: 0 Lepton / LDP, Ratio") ;
      hmctruth_qcd_0lepoverldpratio_1b_htgrouping->Draw() ;
      cqcd3b->cd(pad++) ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->SetTitle("QCD =2b: 0 Lepton / LDP, Ratio") ;
      hmctruth_qcd_0lepoverldpratio_2b_htgrouping->Draw() ;
      cqcd3b->cd(pad++) ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->SetTitle("QCD >=3b: 0 Lepton / LDP, Ratio") ;
      hmctruth_qcd_0lepoverldpratio_3b_htgrouping->Draw() ;

      cqcd3b->cd(pad++) ;
      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->SetTitle("QCD =1b: 0 Lepton / LDP, Scale Factor") ;
      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->Draw() ;
      line->DrawLine(0.5,1.,hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->GetNbinsX()+0.5,1.) ;
      hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->Draw("same") ;
      cqcd3b->cd(pad++) ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->SetTitle("QCD =2b: 0 Lepton / LDP, Scale Factor") ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->Draw() ;
      line->DrawLine(0.5,1.,hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->GetNbinsX()+0.5,1.) ;
      hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->Draw("same") ;
      cqcd3b->cd(pad++) ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->SetTitle("QCD >=3b: 0 Lepton / LDP, Scale Factor") ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->Draw() ;
      line->DrawLine(0.5,1.,hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->GetNbinsX()+0.5,1.) ;
      hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->Draw("same") ;

      TString outqcd3b( infileStr ) ;
      outqcd3b.ReplaceAll("rootfiles","outputfiles") ;
      outqcd3b.ReplaceAll(".root","-mcclosure-qcd3b.pdf") ;
      cqcd3b->SaveAs( outqcd3b ) ;




    //--- insert numbers into file.

      for ( int mbi=0; mbi<nBinsMET ; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT ; hbi++ ) {

            int hbin = 1 + (nBinsMET+1)*hbi + mbi + 1 ;

            char parameterName[1000] ;
            double err, diff, systValue ;

            sprintf( parameterName, "sf_qcd_M%d_H%d_1b_err", mbi+1, hbi+1 ) ;
            err  = hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->GetBinError(hbin)  ;
            diff = hscalefactor_qcd_0lepoverldpratio_1b_htgrouping->GetBinContent(hbin) - 1. ;
            if ( err > 0. ) {
               systValue = sqrt( pow(err,2) + pow(diff,2) ) ;
            } else {
               systValue = 3.0 ;
            }
            updateFileValue( datfile, parameterName, systValue ) ;

            sprintf( parameterName, "sf_qcd_M%d_H%d_2b_err", mbi+1, hbi+1 ) ;
            err  = hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->GetBinError(hbin)  ;
            diff = hscalefactor_qcd_0lepoverldpratio_2b_htgrouping->GetBinContent(hbin) - 1. ;
            if ( err > 0. ) {
               systValue = sqrt( pow(err,2) + pow(diff,2) ) ;
            } else {
               systValue = 3.0 ;
            }
            updateFileValue( datfile, parameterName, systValue ) ;

            sprintf( parameterName, "sf_qcd_M%d_H%d_3b_err", mbi+1, hbi+1 ) ;
            err  = hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->GetBinError(hbin)  ;
            diff = hscalefactor_qcd_0lepoverldpratio_3b_htgrouping->GetBinContent(hbin) - 1. ;
            if ( err > 0. ) {
               systValue = sqrt( pow(err,2) + pow(diff,2) ) ;
            } else {
               systValue = 3.0 ;
            }
            updateFileValue( datfile, parameterName, systValue ) ;

         } // hbi.
      } // mbi.












     //--- ttwj MC LDP / 0lep ratios

      TH1F* hmctruth_ttwj_ldp_1b = (TH1F*) gDirectory->FindObject("hmctruth_ttwj_ldp_1b") ;
      if ( hmctruth_ttwj_ldp_1b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_ttwj_ldp_1b.\n\n") ; return ; }

      TH1F* hmctruth_ttwj_ldpover0lep_ratio_1b = (TH1F*) hmctruth_ttwj_ldp_1b->Clone("hmctruth_ttwj_ldpover0lep_ratio_1b") ;
      hmctruth_ttwj_ldpover0lep_ratio_1b->Divide( hmctruth_ttwj_0lep_1b ) ;


      TH1F* hmctruth_ttwj_ldp_2b = (TH1F*) gDirectory->FindObject("hmctruth_ttwj_ldp_2b") ;
      if ( hmctruth_ttwj_ldp_2b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_ttwj_ldp_2b.\n\n") ; return ; }

      TH1F* hmctruth_ttwj_ldpover0lep_ratio_2b = (TH1F*) hmctruth_ttwj_ldp_2b->Clone("hmctruth_ttwj_ldpover0lep_ratio_2b") ;
      hmctruth_ttwj_ldpover0lep_ratio_2b->Divide( hmctruth_ttwj_0lep_2b ) ;


      TH1F* hmctruth_ttwj_ldp_3b = (TH1F*) gDirectory->FindObject("hmctruth_ttwj_ldp_3b") ;
      if ( hmctruth_ttwj_ldp_3b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_ttwj_ldp_3b.\n\n") ; return ; }

      TH1F* hmctruth_ttwj_ldpover0lep_ratio_3b = (TH1F*) hmctruth_ttwj_ldp_3b->Clone("hmctruth_ttwj_ldpover0lep_ratio_3b") ;
      hmctruth_ttwj_ldpover0lep_ratio_3b->Divide( hmctruth_ttwj_0lep_3b ) ;

      hmctruth_ttwj_ldpover0lep_ratio_1b->SetLineColor(2) ;
      hmctruth_ttwj_ldpover0lep_ratio_2b->SetLineColor(6) ;
      hmctruth_ttwj_ldpover0lep_ratio_3b->SetLineColor(4) ;

      hmctruth_ttwj_ldpover0lep_ratio_1b->SetMarkerStyle(20) ;
      hmctruth_ttwj_ldpover0lep_ratio_2b->SetMarkerStyle(25) ;
      hmctruth_ttwj_ldpover0lep_ratio_3b->SetMarkerStyle(30) ;



      gStyle->SetPadRightMargin(0.08) ;
      TCanvas* cttwjmcr = (TCanvas*) gDirectory->FindObject("cttwjmcr") ;
      if ( cttwjmcr == 0x0 ) {
         cttwjmcr = new TCanvas("cttwjmcr","ttwjmcr closure", 700, 475) ;
      }
      cttwjmcr->Clear() ;

      hmctruth_ttwj_ldpover0lep_ratio_1b->Draw() ;
      hmctruth_ttwj_ldpover0lep_ratio_2b->Draw("same") ;
      hmctruth_ttwj_ldpover0lep_ratio_3b->Draw("same") ;


      //--- write to file

      for ( int mbi=0; mbi<nBinsMET ; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT ; hbi++ ) {

            int hbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

            char parameterName[1000] ;
            double err, value ;


            err   = hmctruth_ttwj_ldpover0lep_ratio_1b->GetBinError(hbin)  ;
            value = hmctruth_ttwj_ldpover0lep_ratio_1b->GetBinContent(hbin)  ;
            sprintf( parameterName, "ttwj_mc_ldpover1lep_ratio_M%d_H%d_1b", mbi+1, hbi+1 ) ;
            printf( "met=%d, ht=%d : %s %s  %6.3f +/- %5.3f\n", mbi+1, hbi+1,
                  parameterName,
                  hmctruth_ttwj_ldpover0lep_ratio_1b -> GetXaxis() -> GetBinLabel( hbin ), value, err ) ;
            updateFileValue( datfile, parameterName, value ) ;
            sprintf( parameterName, "ttwj_mc_ldpover1lep_ratio_M%d_H%d_1b_err", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, err ) ;

            err   = hmctruth_ttwj_ldpover0lep_ratio_2b->GetBinError(hbin)  ;
            value = hmctruth_ttwj_ldpover0lep_ratio_2b->GetBinContent(hbin)  ;
            sprintf( parameterName, "ttwj_mc_ldpover1lep_ratio_M%d_H%d_2b", mbi+1, hbi+1 ) ;
            printf( "met=%d, ht=%d : %s %s  %6.3f +/- %5.3f\n", mbi+1, hbi+1,
                  parameterName,
                  hmctruth_ttwj_ldpover0lep_ratio_2b -> GetXaxis() -> GetBinLabel( hbin ), value, err ) ;
            updateFileValue( datfile, parameterName, value ) ;
            sprintf( parameterName, "ttwj_mc_ldpover1lep_ratio_M%d_H%d_2b_err", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, err ) ;

            err   = hmctruth_ttwj_ldpover0lep_ratio_3b->GetBinError(hbin)  ;
            value = hmctruth_ttwj_ldpover0lep_ratio_3b->GetBinContent(hbin)  ;
            sprintf( parameterName, "ttwj_mc_ldpover1lep_ratio_M%d_H%d_3b", mbi+1, hbi+1 ) ;
            printf( "met=%d, ht=%d : %s %s  %6.3f +/- %5.3f\n", mbi+1, hbi+1,
                  parameterName,
                  hmctruth_ttwj_ldpover0lep_ratio_3b -> GetXaxis() -> GetBinLabel( hbin ), value, err ) ;
            updateFileValue( datfile, parameterName, value ) ;
            sprintf( parameterName, "ttwj_mc_ldpover1lep_ratio_M%d_H%d_3b_err", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, err ) ;

            printf("\n") ;


         } // hbi.
      } // mbi.










     //--- znn MC LDP / 0lep ratios

      TH1F* hmctruth_znn_0lep_1b = (TH1F*) gDirectory->FindObject("hmctruth_znn_0lep_1b") ;
      if ( hmctruth_znn_0lep_1b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_znn_0lep_1b.\n\n") ; return ; }
      TH1F* hmctruth_znn_ldp_1b = (TH1F*) gDirectory->FindObject("hmctruth_znn_ldp_1b") ;
      if ( hmctruth_znn_ldp_1b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_znn_ldp_1b.\n\n") ; return ; }

      TH1F* hmctruth_znn_ldpover0lep_ratio_1b = (TH1F*) hmctruth_znn_ldp_1b->Clone("hmctruth_znn_ldpover0lep_ratio_1b") ;
      hmctruth_znn_ldpover0lep_ratio_1b->Divide( hmctruth_znn_0lep_1b ) ;


      TH1F* hmctruth_znn_0lep_2b = (TH1F*) gDirectory->FindObject("hmctruth_znn_0lep_2b") ;
      if ( hmctruth_znn_0lep_2b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_znn_0lep_2b.\n\n") ; return ; }
      TH1F* hmctruth_znn_ldp_2b = (TH1F*) gDirectory->FindObject("hmctruth_znn_ldp_2b") ;
      if ( hmctruth_znn_ldp_2b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_znn_ldp_2b.\n\n") ; return ; }

      TH1F* hmctruth_znn_ldpover0lep_ratio_2b = (TH1F*) hmctruth_znn_ldp_2b->Clone("hmctruth_znn_ldpover0lep_ratio_2b") ;
      hmctruth_znn_ldpover0lep_ratio_2b->Divide( hmctruth_znn_0lep_2b ) ;


      TH1F* hmctruth_znn_0lep_3b = (TH1F*) gDirectory->FindObject("hmctruth_znn_0lep_3b") ;
      if ( hmctruth_znn_0lep_3b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_znn_0lep_3b.\n\n") ; return ; }
      TH1F* hmctruth_znn_ldp_3b = (TH1F*) gDirectory->FindObject("hmctruth_znn_ldp_3b") ;
      if ( hmctruth_znn_ldp_3b == 0x0 ) { printf("\n\n\n *** can't find hmctruth_znn_ldp_3b.\n\n") ; return ; }

      TH1F* hmctruth_znn_ldpover0lep_ratio_3b = (TH1F*) hmctruth_znn_ldp_3b->Clone("hmctruth_znn_ldpover0lep_ratio_3b") ;
      hmctruth_znn_ldpover0lep_ratio_3b->Divide( hmctruth_znn_0lep_3b ) ;

      hmctruth_znn_ldpover0lep_ratio_1b->SetLineColor(2) ;
      hmctruth_znn_ldpover0lep_ratio_2b->SetLineColor(6) ;
      hmctruth_znn_ldpover0lep_ratio_3b->SetLineColor(4) ;

      hmctruth_znn_ldpover0lep_ratio_1b->SetMarkerStyle(20) ;
      hmctruth_znn_ldpover0lep_ratio_2b->SetMarkerStyle(25) ;
      hmctruth_znn_ldpover0lep_ratio_3b->SetMarkerStyle(30) ;



      gStyle->SetPadRightMargin(0.08) ;
      TCanvas* cznnmcr = (TCanvas*) gDirectory->FindObject("cznnmcr") ;
      if ( cznnmcr == 0x0 ) {
         cznnmcr = new TCanvas("cznnmcr","znnmcr closure", 700, 475) ;
      }
      cznnmcr->Clear() ;

      hmctruth_znn_ldpover0lep_ratio_1b->Draw() ;
      hmctruth_znn_ldpover0lep_ratio_2b->Draw("same") ;
      hmctruth_znn_ldpover0lep_ratio_3b->Draw("same") ;


      //--- MC stats are too low for 2b and >=3b cases.  Use 1b values for all three.

      //--- write to file

      for ( int mbi=0; mbi<nBinsMET ; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT ; hbi++ ) {

            int hbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

            char parameterName[1000] ;
            double err, value ;


            err   = hmctruth_znn_ldpover0lep_ratio_1b->GetBinError(hbin)  ;
            value = hmctruth_znn_ldpover0lep_ratio_1b->GetBinContent(hbin)  ;
            sprintf( parameterName, "znn_mc_ldpover1lep_ratio_M%d_H%d_1b", mbi+1, hbi+1 ) ;
            printf( "met=%d, ht=%d : %s %s  %6.3f +/- %5.3f\n", mbi+1, hbi+1,
                  parameterName,
                  hmctruth_znn_ldpover0lep_ratio_1b -> GetXaxis() -> GetBinLabel( hbin ), value, err ) ;
            updateFileValue( datfile, parameterName, value ) ;
            sprintf( parameterName, "znn_mc_ldpover1lep_ratio_M%d_H%d_1b_err", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, err ) ;
            sprintf( parameterName, "znn_mc_ldpover1lep_ratio_M%d_H%d_2b", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, value ) ;
            sprintf( parameterName, "znn_mc_ldpover1lep_ratio_M%d_H%d_2b_err", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, err ) ;
            sprintf( parameterName, "znn_mc_ldpover1lep_ratio_M%d_H%d_3b", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, value ) ;
            sprintf( parameterName, "znn_mc_ldpover1lep_ratio_M%d_H%d_3b_err", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, err ) ;

      ///   err   = hmctruth_znn_ldpover0lep_ratio_2b->GetBinError(hbin)  ;
      ///   value = hmctruth_znn_ldpover0lep_ratio_2b->GetBinContent(hbin)  ;
      ///   sprintf( parameterName, "znn_mc_ldpover1lep_ratio_M%d_H%d_2b", mbi+1, hbi+1 ) ;
      ///   printf( "met=%d, ht=%d : %s %s  %6.3f +/- %5.3f\n", mbi+1, hbi+1,
      ///         parameterName,
      ///         hmctruth_znn_ldpover0lep_ratio_2b -> GetXaxis() -> GetBinLabel( hbin ), value, err ) ;
      ///   updateFileValue( datfile, parameterName, value ) ;
      ///   sprintf( parameterName, "znn_mc_ldpover1lep_ratio_M%d_H%d_2b_err", mbi+1, hbi+1 ) ;
      ///   updateFileValue( datfile, parameterName, err ) ;

      ///   err   = hmctruth_znn_ldpover0lep_ratio_3b->GetBinError(hbin)  ;
      ///   value = hmctruth_znn_ldpover0lep_ratio_3b->GetBinContent(hbin)  ;
      ///   sprintf( parameterName, "znn_mc_ldpover1lep_ratio_M%d_H%d_3b", mbi+1, hbi+1 ) ;
      ///   printf( "met=%d, ht=%d : %s %s  %6.3f +/- %5.3f\n", mbi+1, hbi+1,
      ///         parameterName,
      ///         hmctruth_znn_ldpover0lep_ratio_3b -> GetXaxis() -> GetBinLabel( hbin ), value, err ) ;
      ///   updateFileValue( datfile, parameterName, value ) ;
      ///   sprintf( parameterName, "znn_mc_ldpover1lep_ratio_M%d_H%d_3b_err", mbi+1, hbi+1 ) ;
      ///   updateFileValue( datfile, parameterName, err ) ;

            printf("\n") ;


         } // hbi.
      } // mbi.









   } // mcclosure2



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

