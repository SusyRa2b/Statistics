
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

   void mcclosure3( const char* infile = "rootfiles/gi-plots-met4-ht4.root",
                    const char* datfile = "Input-met4-ht4-wsyst1.dat",
                    bool doQCDCorrection = true,
                    bool doQCDSyst = true,
                    bool doTTWJCorrection = true,
                    bool doTTWJSyst = true
                    ) {

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
            if ( doTTWJSyst ) {
               updateFileValue( datfile, parameterName, systValue ) ;
            } else {
               updateFileValue( datfile, parameterName, 0.0 ) ;
            }
            sprintf( parameterName, "sf_ttwj_M%d_H%d_1b", mbi+1, hbi+1 ) ;
            if ( doTTWJCorrection ) {
               updateFileValue( datfile, parameterName, correction ) ;
            } else {
               updateFileValue( datfile, parameterName, 1.0 ) ;
            }

            sprintf( parameterName, "sf_ttwj_M%d_H%d_2b_err", mbi+1, hbi+1 ) ;
            err  = hscalefactor_ttwj_0over1ratio_2b->GetBinError(hbin)  ;
            correction = hscalefactor_ttwj_0over1ratio_2b->GetBinContent(hbin) ;
            diff = correction - 1. ;
            if ( err > 0. ) {
               systValue = sqrt( pow(err,2) + pow(0.5*diff,2) ) ;
            } else {
               systValue = 3.0 ;
            }
            if ( doTTWJSyst ) {
               updateFileValue( datfile, parameterName, systValue ) ;
            } else {
               updateFileValue( datfile, parameterName, 0.0 ) ;
            }
            sprintf( parameterName, "sf_ttwj_M%d_H%d_2b", mbi+1, hbi+1 ) ;
            if ( doTTWJCorrection ) {
               updateFileValue( datfile, parameterName, correction ) ;
            } else {
               updateFileValue( datfile, parameterName, 1.0 ) ;
            }

            sprintf( parameterName, "sf_ttwj_M%d_H%d_3b_err", mbi+1, hbi+1 ) ;
            err  = hscalefactor_ttwj_0over1ratio_3b->GetBinError(hbin)  ;
            correction = hscalefactor_ttwj_0over1ratio_3b->GetBinContent(hbin)  ;
            diff = correction - 1. ;
            if ( err > 0. ) {
               systValue = sqrt( pow(err,2) + pow(0.5*diff,2) ) ;
            } else {
               systValue = 3.0 ;
            }
            if ( doTTWJSyst ) {
               updateFileValue( datfile, parameterName, systValue ) ;
            } else {
               updateFileValue( datfile, parameterName, 0.0 ) ;
            }
            sprintf( parameterName, "sf_ttwj_M%d_H%d_3b", mbi+1, hbi+1 ) ;
            if ( doTTWJCorrection ) {
               updateFileValue( datfile, parameterName, correction ) ;
            } else {
               updateFileValue( datfile, parameterName, 1.0 ) ;
            }

         } // hbi.
      } // mbi.



















     //---   Q C D  Part   --------------

     //---  This is written to go with QCD model 3, which is a single ZL/LDP floating scale factor.




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





      float zlsum_all(0.) ;
      float ldpsum_all(0.) ;

      float zlsumw2_all(0.) ;
      float ldpsumw2_all(0.) ;

      for ( int hbi=0; hbi<nBinsHT ; hbi++ ) {

         zlsum_all += zlsum[hbi] ;
         ldpsum_all += ldpsum[hbi] ;

         zlsumw2_all += zlsumw2[hbi] ;
         ldpsumw2_all += ldpsumw2[hbi] ;


      } // hbi.

      float qcd_0lepoverldpratio_all(0.) ;
      float qcd_0lepoverldpratio_all_err(0.) ;

      if ( ldpsum_all > 0 && zlsum_all > 0 ) {
         qcd_0lepoverldpratio_all = zlsum_all / ldpsum_all ;
         qcd_0lepoverldpratio_all_err = qcd_0lepoverldpratio_all * sqrt( zlsumw2_all / pow( zlsum_all,2) + ldpsumw2_all / pow( ldpsum_all, 2) ) ;
      }
      printf(" Average QCD 0lep/LDP ratio = %5.3f +/- %5.3f\n", qcd_0lepoverldpratio_all, qcd_0lepoverldpratio_all_err ) ;

      printf("\n\n") ;



     //----------

      TH1F* hscalefactor_qcd_0lepoverldpratio_1b = (TH1F*)  hmctruth_qcd_0lepoverldpratio_1b->Clone( "hscalefactor_qcd_0lepoverldpratio_1b" ) ;
      TH1F* hscalefactor_qcd_0lepoverldpratio_2b = (TH1F*)  hmctruth_qcd_0lepoverldpratio_2b->Clone( "hscalefactor_qcd_0lepoverldpratio_2b" ) ;
      TH1F* hscalefactor_qcd_0lepoverldpratio_3b = (TH1F*)  hmctruth_qcd_0lepoverldpratio_3b->Clone( "hscalefactor_qcd_0lepoverldpratio_3b" ) ;

      hscalefactor_qcd_0lepoverldpratio_1b -> Scale( 1./ qcd_0lepoverldpratio_all ) ;
      hscalefactor_qcd_0lepoverldpratio_2b -> Scale( 1./ qcd_0lepoverldpratio_all ) ;
      hscalefactor_qcd_0lepoverldpratio_3b -> Scale( 1./ qcd_0lepoverldpratio_all ) ;





      TCanvas* cqcd = (TCanvas*) gDirectory->FindObject("cqcd") ;
      if ( cqcd == 0x0 ) {
         cqcd = new TCanvas("cqcd","qcd closure", 700, 950) ;
      }
      cqcd->Clear() ;
      cqcd->Divide(1,2) ;

      TLegend* legend_qcd = new TLegend( 0.91, 0.77,  0.99, 0.93 ) ;
      legend_qcd->SetFillColor(kWhite) ;
      legend_qcd->AddEntry( hmctruth_qcd_0lepoverldpratio_1b, "=1b") ;
      legend_qcd->AddEntry( hmctruth_qcd_0lepoverldpratio_2b, "=2b") ;
      legend_qcd->AddEntry( hmctruth_qcd_0lepoverldpratio_3b, ">=3b") ;

      cqcd->cd(1) ;
      hmctruth_qcd_0lepoverldpratio_1b->SetTitle("qcd: 0 Lepton / LDP, Ratio") ;
      hmctruth_qcd_0lepoverldpratio_1b->Draw() ;
      hmctruth_qcd_0lepoverldpratio_2b->Draw("same") ;
      hmctruth_qcd_0lepoverldpratio_3b->Draw("same") ;
      legend_qcd->Draw() ;

      cqcd->cd(2) ;
      hscalefactor_qcd_0lepoverldpratio_1b->SetTitle("qcd: 0 Lepton / LDP, Scale Factor") ;
      hscalefactor_qcd_0lepoverldpratio_1b->Draw() ;
      hscalefactor_qcd_0lepoverldpratio_2b->Draw("same") ;
      hscalefactor_qcd_0lepoverldpratio_3b->Draw("same") ;
      line->DrawLine(0.5,1.,hscalefactor_qcd_0lepoverldpratio_1b->GetNbinsX()+0.5,1.) ;
      legend_qcd->Draw() ;

      TString outqcd( infileStr ) ;
      outqcd.ReplaceAll("rootfiles","outputfiles") ;
      outqcd.ReplaceAll(".root","-mcclosure-qcd1.pdf") ;
      cqcd->SaveAs( outqcd ) ;







    //--- insert numbers into file.

      for ( int mbi=0; mbi<nBinsMET ; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT ; hbi++ ) {

            int hbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

            char parameterName[1000] ;
            double sf_val, sf_err ;

            if ( doQCDCorrection ) {

               sprintf( parameterName, "sf_qcd_M%d_H%d_1b", mbi+1, hbi+1 ) ;
               sf_val = hscalefactor_qcd_0lepoverldpratio_1b -> GetBinContent( hbin ) ;
               if ( sf_val == 0. ) { sf_val = 1.0 ; }
               updateFileValue( datfile, parameterName, sf_val ) ;

               sprintf( parameterName, "sf_qcd_M%d_H%d_2b", mbi+1, hbi+1 ) ;
               sf_val = hscalefactor_qcd_0lepoverldpratio_2b -> GetBinContent( hbin ) ;
               if ( sf_val == 0. ) { sf_val = 1.0 ; }
               updateFileValue( datfile, parameterName, sf_val ) ;

               sprintf( parameterName, "sf_qcd_M%d_H%d_3b", mbi+1, hbi+1 ) ;
               sf_val = hscalefactor_qcd_0lepoverldpratio_3b -> GetBinContent( hbin ) ;
               if ( sf_val == 0. ) { sf_val = 1.0 ; }
               updateFileValue( datfile, parameterName, sf_val ) ;

            } else {

               sprintf( parameterName, "sf_qcd_M%d_H%d_1b", mbi+1, hbi+1 ) ;
               updateFileValue( datfile, parameterName, 1.0 ) ;

               sprintf( parameterName, "sf_qcd_M%d_H%d_2b", mbi+1, hbi+1 ) ;
               updateFileValue( datfile, parameterName, 1.0 ) ;

               sprintf( parameterName, "sf_qcd_M%d_H%d_3b", mbi+1, hbi+1 ) ;
               updateFileValue( datfile, parameterName, 1.0 ) ;


            }

            if ( doQCDSyst ) {

               sprintf( parameterName, "sf_qcd_M%d_H%d_1b_err", mbi+1, hbi+1 ) ;
               sf_err = hscalefactor_qcd_0lepoverldpratio_1b -> GetBinError( hbin ) ;
               printf(" bin check : %s  ,  %s\n", parameterName, hscalefactor_qcd_0lepoverldpratio_1b -> GetXaxis() -> GetBinLabel( hbin ) ) ;
               if ( sf_err == 0. ) { sf_err = 3.0 ; }
               updateFileValue( datfile, parameterName, sf_err ) ;

               sprintf( parameterName, "sf_qcd_M%d_H%d_2b_err", mbi+1, hbi+1 ) ;
               sf_err = hscalefactor_qcd_0lepoverldpratio_2b -> GetBinError( hbin ) ;
               printf(" bin check : %s  ,  %s\n", parameterName, hscalefactor_qcd_0lepoverldpratio_2b -> GetXaxis() -> GetBinLabel( hbin ) ) ;
               if ( sf_err == 0. ) { sf_err = 3.0 ; }
               updateFileValue( datfile, parameterName, sf_err ) ;

               sprintf( parameterName, "sf_qcd_M%d_H%d_3b_err", mbi+1, hbi+1 ) ;
               sf_err = hscalefactor_qcd_0lepoverldpratio_3b -> GetBinError( hbin ) ;
               printf(" bin check : %s  ,  %s\n", parameterName, hscalefactor_qcd_0lepoverldpratio_3b -> GetXaxis() -> GetBinLabel( hbin ) ) ;
               if ( sf_err == 0. ) { sf_err = 3.0 ; }
               updateFileValue( datfile, parameterName, sf_err ) ;

            } else {

               sprintf( parameterName, "sf_qcd_M%d_H%d_1b_err", mbi+1, hbi+1 ) ;
               updateFileValue( datfile, parameterName, 0.0 ) ;

               sprintf( parameterName, "sf_qcd_M%d_H%d_2b_err", mbi+1, hbi+1 ) ;
               updateFileValue( datfile, parameterName, 0.0 ) ;

               sprintf( parameterName, "sf_qcd_M%d_H%d_3b_err", mbi+1, hbi+1 ) ;
               updateFileValue( datfile, parameterName, 0.0 ) ;

            }


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
            sprintf( parameterName, "ttwj_mc_ldpover0lep_ratio_M%d_H%d_1b", mbi+1, hbi+1 ) ;
            printf( "met=%d, ht=%d : %s %s  %6.3f +/- %5.3f\n", mbi+1, hbi+1,
                  parameterName,
                  hmctruth_ttwj_ldpover0lep_ratio_1b -> GetXaxis() -> GetBinLabel( hbin ), value, err ) ;
            updateFileValue( datfile, parameterName, value ) ;
            sprintf( parameterName, "ttwj_mc_ldpover0lep_ratio_M%d_H%d_1b_err", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, err ) ;

            err   = hmctruth_ttwj_ldpover0lep_ratio_2b->GetBinError(hbin)  ;
            value = hmctruth_ttwj_ldpover0lep_ratio_2b->GetBinContent(hbin)  ;
            sprintf( parameterName, "ttwj_mc_ldpover0lep_ratio_M%d_H%d_2b", mbi+1, hbi+1 ) ;
            printf( "met=%d, ht=%d : %s %s  %6.3f +/- %5.3f\n", mbi+1, hbi+1,
                  parameterName,
                  hmctruth_ttwj_ldpover0lep_ratio_2b -> GetXaxis() -> GetBinLabel( hbin ), value, err ) ;
            updateFileValue( datfile, parameterName, value ) ;
            sprintf( parameterName, "ttwj_mc_ldpover0lep_ratio_M%d_H%d_2b_err", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, err ) ;

            err   = hmctruth_ttwj_ldpover0lep_ratio_3b->GetBinError(hbin)  ;
            value = hmctruth_ttwj_ldpover0lep_ratio_3b->GetBinContent(hbin)  ;
            sprintf( parameterName, "ttwj_mc_ldpover0lep_ratio_M%d_H%d_3b", mbi+1, hbi+1 ) ;
            printf( "met=%d, ht=%d : %s %s  %6.3f +/- %5.3f\n", mbi+1, hbi+1,
                  parameterName,
                  hmctruth_ttwj_ldpover0lep_ratio_3b -> GetXaxis() -> GetBinLabel( hbin ), value, err ) ;
            updateFileValue( datfile, parameterName, value ) ;
            sprintf( parameterName, "ttwj_mc_ldpover0lep_ratio_M%d_H%d_3b_err", mbi+1, hbi+1 ) ;
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
            sprintf( parameterName, "znn_mc_ldpover0lep_ratio_M%d_H%d_1b", mbi+1, hbi+1 ) ;
            printf( "met=%d, ht=%d : %s %s  %6.3f +/- %5.3f\n", mbi+1, hbi+1,
                  parameterName,
                  hmctruth_znn_ldpover0lep_ratio_1b -> GetXaxis() -> GetBinLabel( hbin ), value, err ) ;
            updateFileValue( datfile, parameterName, value ) ;
            sprintf( parameterName, "znn_mc_ldpover0lep_ratio_M%d_H%d_1b_err", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, err ) ;
            sprintf( parameterName, "znn_mc_ldpover0lep_ratio_M%d_H%d_2b", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, value ) ;
            sprintf( parameterName, "znn_mc_ldpover0lep_ratio_M%d_H%d_2b_err", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, err ) ;
            sprintf( parameterName, "znn_mc_ldpover0lep_ratio_M%d_H%d_3b", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, value ) ;
            sprintf( parameterName, "znn_mc_ldpover0lep_ratio_M%d_H%d_3b_err", mbi+1, hbi+1 ) ;
            updateFileValue( datfile, parameterName, err ) ;

      ///   err   = hmctruth_znn_ldpover0lep_ratio_2b->GetBinError(hbin)  ;
      ///   value = hmctruth_znn_ldpover0lep_ratio_2b->GetBinContent(hbin)  ;
      ///   sprintf( parameterName, "znn_mc_ldpover0lep_ratio_M%d_H%d_2b", mbi+1, hbi+1 ) ;
      ///   printf( "met=%d, ht=%d : %s %s  %6.3f +/- %5.3f\n", mbi+1, hbi+1,
      ///         parameterName,
      ///         hmctruth_znn_ldpover0lep_ratio_2b -> GetXaxis() -> GetBinLabel( hbin ), value, err ) ;
      ///   updateFileValue( datfile, parameterName, value ) ;
      ///   sprintf( parameterName, "znn_mc_ldpover0lep_ratio_M%d_H%d_2b_err", mbi+1, hbi+1 ) ;
      ///   updateFileValue( datfile, parameterName, err ) ;

      ///   err   = hmctruth_znn_ldpover0lep_ratio_3b->GetBinError(hbin)  ;
      ///   value = hmctruth_znn_ldpover0lep_ratio_3b->GetBinContent(hbin)  ;
      ///   sprintf( parameterName, "znn_mc_ldpover0lep_ratio_M%d_H%d_3b", mbi+1, hbi+1 ) ;
      ///   printf( "met=%d, ht=%d : %s %s  %6.3f +/- %5.3f\n", mbi+1, hbi+1,
      ///         parameterName,
      ///         hmctruth_znn_ldpover0lep_ratio_3b -> GetXaxis() -> GetBinLabel( hbin ), value, err ) ;
      ///   updateFileValue( datfile, parameterName, value ) ;
      ///   sprintf( parameterName, "znn_mc_ldpover0lep_ratio_M%d_H%d_3b_err", mbi+1, hbi+1 ) ;
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

