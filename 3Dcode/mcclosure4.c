
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TMinuit.h"

#include "updateFileValue.c"
#include "getFileValue.c"

#include <iostream>
#include <fstream>

  using std::cout ;
  using std::endl ;

   void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;
   TH1F* bookHist(const char* hname, const char* htitle, int sampleIndex ) ;

   double data_Rqcd[20][20][10] ;
   double data_Rqcd_err[20][20][10] ;

   double fit_Rqcd_HT[20] ;
   double fit_SFqcd_MET[20] ;
   double fit_SFqcd_nb[10] ;

   int nBinsMET ;
   int nBinsHT ;
   const int nBinsBjets(3) ;
   const int nQcdSamples(9) ;

  //-----------

   void minuit_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {

      int idummy = npar ;
      double fdummy = gin[0] ;
      fdummy = 0. ;
      idummy = iflag ;

      f = 0. ;

      //--- unpack the stupid par vector.
      int parind(0) ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         fit_Rqcd_HT[hbi] = par[parind] ;
         parind ++ ;
      } // hbi.
      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         if ( mbi == 0 ) {
            fit_SFqcd_MET[mbi] = 1.0 ;
         } else {
            fit_SFqcd_MET[mbi] = par[parind] ;
            parind++ ;
         }
      } // mbi.
      for ( int bbi=0; bbi<nBinsMET; bbi++ ) {
         if ( bbi == 0 ) {
            fit_SFqcd_nb[bbi] = 1.0 ;
         } else {
            fit_SFqcd_nb[bbi] = par[parind] ;
            parind++ ;
         }
      } // bbi.

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
               if ( !( data_Rqcd[mbi][hbi][bbi] > 0. && data_Rqcd_err[mbi][hbi][bbi] > 0. ) ) { continue ; }
               double delta = data_Rqcd[mbi][hbi][bbi] - fit_Rqcd_HT[hbi] * fit_SFqcd_MET[mbi] * fit_SFqcd_nb[bbi] ;
               f += delta*delta / (data_Rqcd_err[mbi][hbi][bbi] * data_Rqcd_err[mbi][hbi][bbi] ) ;
            } // bbi.
         } // hbi.
      } // mbi.

   } // minuit_fcn.

  //-----------

  //=============================================================================================================================

   void mcclosure4( const char* infile = "rootfiles/gi-plots-met4-ht4.root",
                    const char* datfile = "Input-met4-ht4-wsyst1.dat",
                    int qcdModelIndex = 4,
                    bool setExpected0lepObs = true
                    ) {

      if ( qcdModelIndex < 2 || qcdModelIndex > 4 ) {
         printf("\n\n *** Unsupported qcdModelIndex (%d).  Try 2, 3, or 4.\n\n", qcdModelIndex ) ;
      }

      bool doTTWJSyst = true ;
      bool doTTWJCorrection = true ;
      bool savePlots = true ;

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

     double htbin_0lepldp_ratio[nBinsHT] ;
     double htbin_0lepldp_ratio_global(0.) ;

     { //-- scoping bracket for QCD part.


      char samplename[9][100] = {
        "qcd_0120_to_0170"
       ,"qcd_0170_to_0300"
       ,"qcd_0300_to_0470"
       ,"qcd_0470_to_0600"
       ,"qcd_0600_to_0800"
       ,"qcd_0800_to_1000"
       ,"qcd_1000_to_1400"
       ,"qcd_1400_to_1800"
       ,"qcd_1800_to_9999"
      } ;

      int samplecolor[9] = {
         800+8,
         632,
         616+2,
         600+1,
         860+2,
         416+2,
         900-2,
         880+2,
         800+2
      } ;


      TH2F*   h0lep[nQcdSamples][nBinsBjets] ;
      TH2F*   hldp [nQcdSamples][nBinsBjets] ;

      TCanvas* cqcd = (TCanvas*) gDirectory->FindObject("cqcd") ;
      if ( cqcd == 0x0 ) {
         cqcd = new TCanvas("cqcd", "qcd study", 700, 900 ) ;
      }

      char hname[1000] ;
      char htitle[1000] ;

      loadHist( "rootfiles/qcd-study1.root" ) ;

      for ( int si=0; si<nQcdSamples; si++ ) {
         for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

            sprintf( hname, "h_0lep_%db_%s", bbi+1, samplename[si] ) ;
            printf("  loading %s\n", hname ) ;
            h0lep[si][bbi] = (TH2F*) gDirectory->FindObject( hname ) ;
            if ( h0lep[si][bbi] == 0x0 ) { printf("\n\n *** %s missing.\n\n", hname ) ; return ; }
            sprintf( hname, "h_ldp_%db_%s", bbi+1, samplename[si] ) ;
            printf("  loading %s\n", hname ) ;
            hldp [si][bbi] = (TH2F*) gDirectory->FindObject( hname ) ;
            if ( hldp [si][bbi] == 0x0 ) { printf("\n\n *** %s missing.\n\n", hname ) ; return ; }

         } // bbi.
      } // si.






     //--- Flatten 2D histos and compute 0lep/LDP ratios for each sample.

      int nbins = nBinsMET*(nBinsHT+1) + 1 ;

      TH1F* hflat_0lep[nQcdSamples][nBinsBjets] ;
      TH1F* hflat_ldp [nQcdSamples][nBinsBjets] ;
      TH1F* hflat_0lepldp_ratio[nQcdSamples][nBinsBjets] ;

      sprintf( hname, "hflat_0lep_all" ) ;
      sprintf( htitle, "All QCD 0lep events" ) ;
      TH1F* hflat_0lep_all = bookHist( hname, htitle, 5 ) ;
      sprintf( hname, "hflat_ldp_all" ) ;
      sprintf( htitle, "All QCD ldp events" ) ;
      TH1F* hflat_ldp_all = bookHist( hname, htitle, 5 ) ;


      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
         for ( int si=0; si<nQcdSamples; si++ ) {

            sprintf( hname, "hflat_0lep_%db_%s", bbi+1, samplename[si] ) ;
            sprintf( htitle, "QCD 0lep events, nb=%d, %s", bbi+1, samplename[si] ) ;
            hflat_0lep[si][bbi] = bookHist( hname, htitle, si ) ;
            hflat_0lep[si][bbi]->SetFillColor(11) ;
            hldp[si][bbi]->Print("all") ;

            sprintf( hname, "hflat_ldp_%db_%s", bbi+1, samplename[si] ) ;
            sprintf( htitle, "QCD LDP events, nb=%d, %s", bbi+1, samplename[si] ) ;
            hflat_ldp[si][bbi] = bookHist( hname, htitle, si ) ;
            hflat_ldp[si][bbi]->SetFillColor(11) ;

            sprintf( hname, "hflat_0lepldp_ratio_%db_%s", bbi+1, samplename[si] ) ;
            sprintf( htitle, "QCD 0lep/LDP ratio, nb=%d, %s", bbi+1, samplename[si] ) ;
            hflat_0lepldp_ratio[si][bbi] = bookHist( hname, htitle, si ) ;

            for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
               for ( int hbi=0; hbi<nBinsHT; hbi++ ) {

                  int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

                  float n0lep     = h0lep[si][bbi] -> GetBinContent( mbi+1, hbi+1 ) ;
                  float n0lep_err = h0lep[si][bbi] -> GetBinError( mbi+1, hbi+1 ) ;

                  float nldp      = hldp [si][bbi] -> GetBinContent( mbi+1, hbi+1 ) ;
                  float nldp_err  = hldp [si][bbi] -> GetBinError( mbi+1, hbi+1 ) ;

                  hflat_0lep[si][bbi] -> SetBinContent( histbin, n0lep ) ;
                  hflat_0lep[si][bbi] -> SetBinError( histbin, n0lep_err ) ;

                  hflat_ldp[si][bbi]  -> SetBinContent( histbin, nldp ) ;
                  hflat_ldp[si][bbi]  -> SetBinError( histbin, nldp_err ) ;

                  double ratio(0.) ;
                  if ( nldp > 0 ) { ratio = n0lep / nldp ; }
                  double ratio_err(0.) ;
                  if ( n0lep_err > 0. && nldp_err > 0. ) {
                     ratio_err = ratio * sqrt( pow(n0lep_err/n0lep,2) + pow(nldp_err/nldp,2) ) ;
                  }

                  hflat_0lepldp_ratio[si][bbi] -> SetBinContent( histbin, ratio ) ;
                  hflat_0lepldp_ratio[si][bbi] -> SetBinError( histbin, ratio_err ) ;


               } // hbi.
            } // mbi.

            hflat_0lep_all -> Add( hflat_0lep[si][bbi] ) ;
            hflat_ldp_all  -> Add( hflat_ldp [si][bbi] ) ;

            cqcd->Clear() ;
            cqcd->Divide(1,3) ;

            cqcd->cd(1) ;
            hflat_0lep[si][bbi] -> DrawCopy("histe") ;
            hflat_0lep[si][bbi] -> DrawCopy("samee") ;

            cqcd->cd(2) ;
            hflat_ldp[si][bbi] -> DrawCopy("histe") ;
            hflat_ldp[si][bbi] -> DrawCopy("samee") ;

            cqcd->cd(3) ;
            hflat_0lepldp_ratio[si][bbi]->SetMinimum(-0.1) ;
            hflat_0lepldp_ratio[si][bbi]->SetMaximum(0.6) ;
            hflat_0lepldp_ratio[si][bbi]->SetMarkerStyle(20) ;
            hflat_0lepldp_ratio[si][bbi]->SetLineWidth(2) ;
            hflat_0lepldp_ratio[si][bbi]->DrawCopy() ;
            gPad->SetGridx(1) ;
            gPad->SetGridy(1) ;
            line->DrawLine(0.5,0,nbins+0.5,0) ;


            cqcd->Update() ; cqcd->Draw() ;

         } // si.
      } // bbi.

      TH1F* hflat_0lep_all_sqrtNerrs = (TH1F*) hflat_0lep_all->Clone("hflat_0lep_all_sqrtNerrs") ;
      TH1F* hflat_ldp_all_sqrtNerrs  = (TH1F*) hflat_ldp_all ->Clone("hflat_ldp_all_sqrtNerrs") ;

      for ( int bi=1; bi<=nbins; bi++ ) {
         double val ;
         val = hflat_0lep_all_sqrtNerrs->GetBinContent( bi ) ;
         hflat_0lep_all_sqrtNerrs -> SetBinError( bi, sqrt(val) ) ;
         val = hflat_ldp_all_sqrtNerrs->GetBinContent( bi ) ;
         hflat_ldp_all_sqrtNerrs -> SetBinError( bi, sqrt(val) ) ;
      }












     //--- Compute sample average of 0lep/LDP ratio


      TH1F* hflat_0lepldp_ratio_ave[nBinsBjets] ;

      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

         sprintf( hname, "hflat_0lepldp_ratio_ave_%db", bbi+1 ) ;
         sprintf( htitle, "QCD 0lep/LDP average ratio, nb=%d", bbi+1 ) ;
         hflat_0lepldp_ratio_ave[bbi] = bookHist( hname, htitle, 5 ) ;

         hflat_0lepldp_ratio_ave[bbi] -> SetMarkerStyle(24) ;
         hflat_0lepldp_ratio_ave[bbi] -> SetMarkerSize(2.0) ;
         hflat_0lepldp_ratio_ave[bbi] -> SetLineWidth(3) ;

         for ( int binind=1; binind<=nbins; binind++ ) {

            double wvsum(0.) ;
            double wsum(0.) ;

            for ( int si=0; si<nQcdSamples; si++ ) {

               double val = hflat_0lepldp_ratio[si][bbi] -> GetBinContent( binind ) ;
               double err = hflat_0lepldp_ratio[si][bbi] -> GetBinError( binind ) ;
               if ( err > 0 ) {
                  double weight = 1./(err*err) ;
                  wvsum += val*weight ;
                  wsum  += weight ;
               }
            } // si.

            if ( wsum > 0. ) {
               hflat_0lepldp_ratio_ave[bbi] -> SetBinContent( binind, wvsum/wsum ) ;
               hflat_0lepldp_ratio_ave[bbi] -> SetBinError( binind, 1./sqrt(wsum) ) ;
            }

         } // binind

      } // bbi.








     //--- Compute RMS spread of samples and total uncertainty


      TH1F* hflat_0lepldp_ratio_ave_withRMSerror[nBinsBjets] ;

      printf("\n\n") ;
      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

         sprintf( hname, "hflat_0lepldp_ratio_ave_%db_withRMSerror", bbi+1 ) ;
         sprintf( htitle, "QCD 0lep/LDP average ratio, RMS included in error, nb=%d", bbi+1 ) ;
         hflat_0lepldp_ratio_ave_withRMSerror[bbi] = bookHist( hname, htitle, 5 ) ;

         hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetMarkerStyle(24) ;
         hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetMarkerSize(2.0) ;
         hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetLineWidth(3) ;

         for ( int binind=1; binind<=nbins; binind++ ) {

            double sumsquare(0.) ;
            int    nsum(0) ;

            double ave     = hflat_0lepldp_ratio_ave[bbi] -> GetBinContent( binind ) ;
            double ave_err = hflat_0lepldp_ratio_ave[bbi] -> GetBinError( binind ) ;

            for ( int si=0; si<nQcdSamples; si++ ) {

               double val = hflat_0lepldp_ratio[si][bbi] -> GetBinContent( binind ) ;
               double err = hflat_0lepldp_ratio[si][bbi] -> GetBinError( binind ) ;

               if ( err <= 0 ) { continue ; }
               if ( val <= 0 ) { continue ; }
               if ( err/val > 0.5 ) { continue ; } //-- do not include points with very large errors.


               double diff = val - ave ;

               sumsquare += diff*diff ;

               nsum ++ ;

            } // si.

            if ( nsum > 1 ) {
               double rms = sqrt(sumsquare/nsum) ;
               double totalerr = sqrt( ave_err*ave_err + rms*rms ) ;
               const char* bl = hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> GetXaxis() -> GetBinLabel( binind ) ;
               printf(" %s : %s : ratio = %5.3f, stat err = %5.3f, RMS = %5.3f,  total err = %5.3f\n",
                   hname, bl, ave, ave_err, rms, totalerr ) ;
               hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetBinContent( binind, ave ) ;
               hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetBinError( binind, totalerr ) ;
            } else {
               hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetBinContent( binind, ave ) ;
               hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> SetBinError( binind, ave_err ) ;
            }

         } // binind
         printf("\n") ;

      } // bbi.
      printf("\n\n") ;






     //--- compute average 0lep/LDP ratio globally and for each HT bin.

      double all0lep_ht[nBinsHT] ;
      double allldp_ht[nBinsHT] ;
      double all0lep(0.) ;
      double allldp(0.) ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         all0lep_ht[hbi] = 0. ;
         allldp_ht[hbi] = 0. ;
      } // hbi.

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

               all0lep_ht[hbi] += hflat_0lep_all -> GetBinContent( histbin ) ;
               allldp_ht[hbi]  += hflat_ldp_all  -> GetBinContent( histbin ) ;

               all0lep += hflat_0lep_all -> GetBinContent( histbin ) ;
               allldp  += hflat_ldp_all  -> GetBinContent( histbin ) ;

            } // bbi.
         } // hbi.
      } // mbi.


      printf("\n\n") ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         htbin_0lepldp_ratio[hbi] = 0. ;
         if ( allldp_ht[hbi] > 0. ) {
            htbin_0lepldp_ratio[hbi] = all0lep_ht[hbi] / allldp_ht[hbi] ;
         }
         printf(" HT bin %d : QCD :  0lep = %7.1f,   LDP = %7.1f,   0lep/LDP = %5.3f\n", hbi+1, all0lep_ht[hbi], allldp_ht[hbi], htbin_0lepldp_ratio[hbi] ) ;
      } // hbi.
      printf("\n\n") ;
      htbin_0lepldp_ratio_global = all0lep / allldp ;
      printf(" global    : QCD :  0lep = %7.1f,   LDP = %7.1f,   0lep/LDP = %5.3f\n\n\n", all0lep, allldp, htbin_0lepldp_ratio_global ) ;







     //--- Calculations for QCD Model 4 --------------

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

               double val, err ;

               val = hflat_0lepldp_ratio_ave[bbi] -> GetBinContent( histbin ) ;
               err = hflat_0lepldp_ratio_ave[bbi] -> GetBinError( histbin ) ;

               data_Rqcd[mbi][hbi][bbi] = val ;
               data_Rqcd_err[mbi][hbi][bbi] = err ;

            } // bbi.
         } // hbi.
      } // mbi.


      TMinuit *myMinuit = new TMinuit(nBinsHT+nBinsMET-1+nBinsBjets-1) ; // arg is # of parameters

      myMinuit->SetFCN( minuit_fcn ) ;

      Double_t arglist[10] ;
      Int_t ierflg = 0 ;

      arglist[0] = 1 ;
      myMinuit->mnexcm("SET ERR", arglist ,1,ierflg); //--- do this for chi2 fit.

      int parind(0) ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         char pname[1000] ;
         sprintf( pname, "Rqcd_HT%d", hbi+1 ) ;
         myMinuit->mnparm( parind, pname, htbin_0lepldp_ratio[hbi], 0.03, 0., 2., ierflg ) ;
         parind++ ;
      } // hbi.
      for ( int mbi=1; mbi<nBinsHT; mbi++ ) {
         char pname[1000] ;
         sprintf( pname, "SFqcd_MET%d", mbi+1 ) ;
         myMinuit->mnparm( parind, pname, 1.0, 0.10, 0., 2., ierflg ) ;
         parind++ ;
      } // mbi.
      for ( int bbi=1; bbi<nBinsHT; bbi++ ) {
         char pname[1000] ;
         sprintf( pname, "SFqcd_nb%d", bbi+1 ) ;
         myMinuit->mnparm( parind, pname, 1.0, 0.10, 0., 2., ierflg ) ;
         parind++ ;
      } // mbi.

      myMinuit->Migrad() ;
      myMinuit->mncomd("hesse",ierflg) ;

      printf("\n\n") ;
      parind = 0 ;
      for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
         char pname[1000] ;
         double val, err ;
         sprintf( pname, "Rqcd_HT%d", hbi+1 ) ;
         myMinuit->GetParameter( parind, val, err ) ;
         printf(" %11s : %6.3f +/- %5.3f\n", pname, val, err ) ;
         fit_Rqcd_HT[hbi] = val ;
         parind++ ;
      } // hbi.
      for ( int mbi=1; mbi<nBinsHT; mbi++ ) {
         char pname[1000] ;
         double val, err ;
         sprintf( pname, "SFqcd_MET%d", mbi+1 ) ;
         myMinuit->GetParameter( parind, val, err ) ;
         printf(" %11s : %6.3f +/- %5.3f\n", pname, val, err ) ;
         fit_SFqcd_MET[mbi] = val ;
         parind++ ;
      } // mbi.
      for ( int bbi=1; bbi<nBinsHT; bbi++ ) {
         char pname[1000] ;
         double val, err ;
         sprintf( pname, "SFqcd_nb%d", bbi+1 ) ;
         myMinuit->GetParameter( parind, val, err ) ;
         printf(" %11s : %6.3f +/- %5.3f\n", pname, val, err ) ;
         fit_SFqcd_nb[bbi] = val ;
         parind++ ;
      } // mbi.
      printf("\n\n") ;




      TH1F* hflat_0lepldp_ratio_model[nBinsBjets] ;

      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
         sprintf( hname, "hflat_0lepldp_ratio_model_%db", bbi+1 ) ;
         hflat_0lepldp_ratio_model[bbi] = (TH1F*) hflat_0lepldp_ratio_ave[bbi]->Clone( hname ) ;
      } // bbi.


      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

               hflat_0lepldp_ratio_model[bbi] -> SetBinContent( histbin, fit_Rqcd_HT[hbi] * fit_SFqcd_MET[mbi] * fit_SFqcd_nb[bbi] ) ;

            } // bbi.
         } // hbi.
      } // mbi.










     //--- Compute QCD scale factors

      TH1F* hflat_scale_factor[nBinsBjets] ;
      TH1F* hflat_scale_factor_withRMSerror[nBinsBjets] ;

      for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

         sprintf( hname, "hflat_scale_factor_%db", bbi+1 ) ;
         sprintf( htitle, "QCD scale factor, nb=%d", bbi+1 ) ;
         hflat_scale_factor[bbi] = (TH1F*) hflat_0lepldp_ratio_ave[bbi]->Clone( hname ) ;
         hflat_scale_factor[bbi] -> SetTitle( htitle ) ;

         sprintf( hname, "hflat_scale_factor_%db_withRMSerror", bbi+1 ) ;
         sprintf( htitle, "QCD scale factor, nb=%d, RMS error included ", bbi+1 ) ;
         hflat_scale_factor_withRMSerror[bbi] = (TH1F*) hflat_0lepldp_ratio_ave_withRMSerror[bbi]->Clone( hname ) ;
         hflat_scale_factor_withRMSerror[bbi] -> SetTitle( htitle ) ;

      } // bbi.

      for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

               int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;

               double val, err ;
               double valwrmse, errwrmse ;

               val = hflat_0lepldp_ratio_ave[bbi] -> GetBinContent( histbin ) ;
               err = hflat_0lepldp_ratio_ave[bbi] -> GetBinError( histbin ) ;

               valwrmse = hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> GetBinContent( histbin ) ;
               errwrmse = hflat_0lepldp_ratio_ave_withRMSerror[bbi] -> GetBinError( histbin ) ;


               if ( qcdModelIndex == 4 ) {

                  //-- QCD Model 4.  SF_ijk = (0lep/LDP)_ijk / (Rqcd_HTj * SF_METi * SF_nbk).

                  hflat_scale_factor[bbi] -> SetBinContent( histbin, val / ( fit_Rqcd_HT[hbi] * fit_SFqcd_MET[mbi] * fit_SFqcd_nb[bbi] ) ) ;
                  hflat_scale_factor[bbi] -> SetBinError(   histbin, err / ( fit_Rqcd_HT[hbi] * fit_SFqcd_MET[mbi] * fit_SFqcd_nb[bbi] ) ) ;

                  hflat_scale_factor_withRMSerror[bbi] -> SetBinContent( histbin, valwrmse / ( fit_Rqcd_HT[hbi] * fit_SFqcd_MET[mbi] * fit_SFqcd_nb[bbi] ) ) ;
                  hflat_scale_factor_withRMSerror[bbi] -> SetBinError(   histbin, errwrmse / ( fit_Rqcd_HT[hbi] * fit_SFqcd_MET[mbi] * fit_SFqcd_nb[bbi] ) ) ;

               } else if ( qcdModelIndex == 3 ) {

                  //-- QCD Model 3.  SF_i = (0lep/LDP)_i / global_(0lep/LDP)

                  hflat_scale_factor[bbi] -> SetBinContent( histbin, val / htbin_0lepldp_ratio_global ) ;
                  hflat_scale_factor[bbi] -> SetBinError(   histbin, err / htbin_0lepldp_ratio_global ) ;

                  hflat_scale_factor_withRMSerror[bbi] -> SetBinContent( histbin, valwrmse / htbin_0lepldp_ratio_global ) ;
                  hflat_scale_factor_withRMSerror[bbi] -> SetBinError(   histbin, errwrmse / htbin_0lepldp_ratio_global ) ;

               } else if ( qcdModelIndex == 2 ) {

                  //-- QCD Model 2.  SF_i = (0lep/LDP)_i / HTbin_(0lep/LDP)

                  hflat_scale_factor[bbi] -> SetBinContent( histbin, val / htbin_0lepldp_ratio[hbi] ) ;
                  hflat_scale_factor[bbi] -> SetBinError(   histbin, err / htbin_0lepldp_ratio[hbi] ) ;

                  hflat_scale_factor_withRMSerror[bbi] -> SetBinContent( histbin, valwrmse / htbin_0lepldp_ratio[hbi] ) ;
                  hflat_scale_factor_withRMSerror[bbi] -> SetBinError(   histbin, errwrmse / htbin_0lepldp_ratio[hbi] ) ;

               }

            } // bbi.
         } // hbi.
      } // mbi.










     //--- set Scale Factor values in dat file if one is provided.

      if ( strcmp( datfile, "null" ) != 0 ) {

         printf("\n\n\n Setting QCD scale factors in %s\n\n", datfile ) ;

         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
               for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {
                  char parname[1000] ;
                  int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
                  double val = hflat_scale_factor_withRMSerror[bbi]->GetBinContent( histbin ) ;
                  double err = hflat_scale_factor_withRMSerror[bbi]->GetBinError( histbin ) ;
                  if ( err <= 0 ) {
                     val = 1.0 ;
                     err = 3.0 ;
                  }
                  sprintf( parname, "sf_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  updateFileValue( datfile, parname, val ) ;
                  sprintf( parname, "sf_qcd_M%d_H%d_%db_err", mbi+1, hbi+1, bbi+1 ) ;
                  updateFileValue( datfile, parname, err ) ;
               } // bbi.
            } // hbi.
         } // mbi.

      }


     //===========  End QCD calculations.  QCD plots and other output below here.  ===================================


     //--- Plot all samples together

      TCanvas* cqcd2 = (TCanvas*) gDirectory->FindObject("cqcd2") ;
      if ( cqcd2 == 0x0 ) {
         cqcd2 = new TCanvas("cqcd2", "qcd study", 1700, 600 ) ;
      }

      TLegend* l2 = new TLegend( 0.79, 0.70,  0.99, 0.99 ) ;
      l2->AddEntry( hflat_0lepldp_ratio[0][0], "120 to 170" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[1][0], "170 to 300" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[2][0], "300 to 470" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[3][0], "470 to 600" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[4][0], "600 to 800" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[5][0], "800 to 1000" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[6][0], "1000 to 1400" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[7][0], "1400 to 1800" ) ;
      l2->AddEntry( hflat_0lepldp_ratio[8][0], "> 1800" ) ;
      l2->AddEntry( hflat_0lepldp_ratio_ave[0], "Average ratio" ) ;


      cqcd2 -> Clear() ;
      cqcd2 -> Divide(3,1) ;

      cqcd2 -> cd(1) ;
      hflat_0lepldp_ratio[0][0] -> SetTitle("QCD 0lep/LDP ratio, nb=1") ;
      for ( int si=0; si<nQcdSamples; si++ ) {
         hflat_0lepldp_ratio[si][0] -> SetMarkerStyle(20+si) ;
         hflat_0lepldp_ratio[si][0] -> SetLineColor(samplecolor[si]) ;
         hflat_0lepldp_ratio[si][0] -> SetMarkerColor(samplecolor[si]) ;
         if ( si == 0 ) { hflat_0lepldp_ratio[si][0] -> DrawCopy() ; } else { hflat_0lepldp_ratio[si][0] -> DrawCopy("same") ; }
      }
      hflat_0lepldp_ratio_ave[0]->DrawCopy("same") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      l2->Draw() ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd2 -> cd(2) ;
      hflat_0lepldp_ratio[0][1] -> SetTitle("QCD 0lep/LDP ratio, nb=2") ;
      for ( int si=0; si<nQcdSamples; si++ ) {
         hflat_0lepldp_ratio[si][1] -> SetMarkerStyle(20+si) ;
         hflat_0lepldp_ratio[si][1] -> SetLineColor(samplecolor[si]) ;
         hflat_0lepldp_ratio[si][1] -> SetMarkerColor(samplecolor[si]) ;
         if ( si == 0 ) { hflat_0lepldp_ratio[si][1] -> DrawCopy() ; } else { hflat_0lepldp_ratio[si][1] -> DrawCopy("same") ; }
      }
      hflat_0lepldp_ratio_ave[1]->DrawCopy("same") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      l2->Draw() ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd2 -> cd(3) ;
      hflat_0lepldp_ratio[0][2] -> SetTitle("QCD 0lep/LDP ratio, nb>=3") ;
      for ( int si=0; si<nQcdSamples; si++ ) {
         hflat_0lepldp_ratio[si][2] -> SetMarkerStyle(20+si) ;
         hflat_0lepldp_ratio[si][2] -> SetLineColor(samplecolor[si]) ;
         hflat_0lepldp_ratio[si][2] -> SetMarkerColor(samplecolor[si]) ;
         if ( si == 0 ) { hflat_0lepldp_ratio[si][2] -> DrawCopy() ; } else { hflat_0lepldp_ratio[si][2] -> DrawCopy("same") ; }
      }
      hflat_0lepldp_ratio_ave[2]->DrawCopy("same") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      l2->Draw() ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;


      cqcd2->Update() ; cqcd2->Draw() ;

      if ( savePlots ) {
         cqcd2 -> SaveAs("outputfiles/mcclosure4-qcd-allonone.png") ;
         cqcd2 -> SaveAs("outputfiles/mcclosure4-qcd-allonone.pdf") ;
      }








     //--- Go through the samples, 1 by 1, and draw numerator, denominator, and ratio

      TCanvas* cqcd3 = (TCanvas*) gDirectory->FindObject("cqcd3") ;
      if ( cqcd3 == 0x0 ) {
         cqcd3 = new TCanvas("cqcd3", "qcd study", 1200, 900 ) ;
      }

      char oldtitle[1000] ;
      sprintf( oldtitle, "QCD 0lep/LDP ratio, nb=1, %s", samplename[0] ) ;
      hflat_0lepldp_ratio[0][0] -> SetTitle(oldtitle) ;
      sprintf( oldtitle, "QCD 0lep/LDP ratio, nb=2, %s", samplename[0] ) ;
      hflat_0lepldp_ratio[0][1] -> SetTitle(oldtitle) ;
      sprintf( oldtitle, "QCD 0lep/LDP ratio, nb>=3, %s", samplename[0] ) ;
      hflat_0lepldp_ratio[0][2] -> SetTitle(oldtitle) ;

      cqcd3->Clear() ;
      cqcd3->Divide(3,3) ;

      for ( int si=0; si<nQcdSamples; si++ ) {
         TLegend* l3 = new TLegend( 0.65, 0.78, 0.99, 0.92) ;
         l3->AddEntry( hflat_0lepldp_ratio[si][0], samplename[si] ) ;
         l3->AddEntry( hflat_0lepldp_ratio_ave[0], "Average") ;
         for ( int bbi=0; bbi<nBinsBjets; bbi++ ) {

            cqcd3->cd(1+bbi) ;
            hflat_0lep[si][bbi] -> DrawCopy("histe") ;
            hflat_0lep[si][bbi] -> DrawCopy("samee") ;

            cqcd3->cd(4+bbi) ;
            hflat_ldp[si][bbi] -> DrawCopy("histe") ;
            hflat_ldp[si][bbi] -> DrawCopy("samee") ;


            cqcd3->cd(7+bbi) ;
            gPad->SetGridx(1) ;
            gPad->SetGridy(1) ;
            hflat_0lepldp_ratio[si][bbi]->DrawCopy() ;
            hflat_0lepldp_ratio_ave[bbi]->DrawCopy("same") ;
            line->DrawLine(0.5,0,nbins+0.5,0) ;
            l3->Draw() ;

         } // bbi.

         cqcd3->Update() ; cqcd3->Draw() ;

         if ( savePlots ) {
            char filename[1000] ;
            sprintf( filename, "outputfiles/mcclosure4-qcd-%s.pdf", samplename[si] ) ;
            cqcd3 -> SaveAs( filename ) ;
            sprintf( filename, "outputfiles/mcclosure4-qcd-%s.png", samplename[si] ) ;
            cqcd3 -> SaveAs( filename ) ;
         }

      } // si.







     //--- Draw the average 0lep/LDP ratio with MC stat error and total error.

      TCanvas* cqcd4 = (TCanvas*) gDirectory->FindObject("cqcd4") ;
      if ( cqcd4 == 0x0 ) {
         cqcd4 = new TCanvas("cqcd4", "qcd study", 1700, 600 ) ;
      }

      gStyle->SetEndErrorSize(5) ;

      cqcd4->Clear() ;
      cqcd4->Divide(3,1) ;

      cqcd4 -> cd(1) ;
      hflat_0lepldp_ratio_ave[0] -> SetTitle("Average QCD 0lep/LDP ratio, nb=1") ;
      hflat_0lepldp_ratio_ave[0] -> SetMarkerSize(1.2) ;
      hflat_0lepldp_ratio_ave[0] -> SetMarkerStyle(20) ;
      hflat_0lepldp_ratio_ave[0] -> SetMaximum(0.6) ;
      hflat_0lepldp_ratio_ave[0] -> SetMinimum(-0.1) ;
      hflat_0lepldp_ratio_ave[0]->DrawCopy("e1") ;
      hflat_0lepldp_ratio_ave_withRMSerror[0]->SetMarkerStyle() ;
      hflat_0lepldp_ratio_ave_withRMSerror[0]->SetLineColor(4) ;
      hflat_0lepldp_ratio_ave_withRMSerror[0]->DrawCopy("samee1") ;
      hflat_0lepldp_ratio_ave[0]->DrawCopy("samee1") ;
      hflat_0lepldp_ratio_model[0] -> SetLineColor(2) ;
      hflat_0lepldp_ratio_model[0] -> SetLineWidth(2) ;
      hflat_0lepldp_ratio_model[0] -> Draw("histsame") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd4 -> cd(2) ;
      hflat_0lepldp_ratio_ave[1] -> SetTitle("Average QCD 0lep/LDP ratio, nb=2") ;
      hflat_0lepldp_ratio_ave[1] -> SetMarkerSize(1.2) ;
      hflat_0lepldp_ratio_ave[1] -> SetMarkerStyle(20) ;
      hflat_0lepldp_ratio_ave[1] -> SetMaximum(0.6) ;
      hflat_0lepldp_ratio_ave[1] -> SetMinimum(-0.1) ;
      hflat_0lepldp_ratio_ave[1]->DrawCopy("e1") ;
      hflat_0lepldp_ratio_ave_withRMSerror[1]->SetMarkerStyle() ;
      hflat_0lepldp_ratio_ave_withRMSerror[1]->SetLineColor(4) ;
      hflat_0lepldp_ratio_ave_withRMSerror[1]->DrawCopy("samee1") ;
      hflat_0lepldp_ratio_ave[1]->DrawCopy("samee1") ;
      hflat_0lepldp_ratio_model[1] -> SetLineColor(2) ;
      hflat_0lepldp_ratio_model[1] -> SetLineWidth(2) ;
      hflat_0lepldp_ratio_model[1] -> Draw("histsame") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd4 -> cd(3) ;
      hflat_0lepldp_ratio_ave[2] -> SetTitle("Average QCD 0lep/LDP ratio, nb>=3") ;
      hflat_0lepldp_ratio_ave[2] -> SetMarkerSize(1.2) ;
      hflat_0lepldp_ratio_ave[2] -> SetMarkerStyle(20) ;
      hflat_0lepldp_ratio_ave[2] -> SetMaximum(0.6) ;
      hflat_0lepldp_ratio_ave[2] -> SetMinimum(-0.1) ;
      hflat_0lepldp_ratio_ave[2]->DrawCopy("e1") ;
      hflat_0lepldp_ratio_ave_withRMSerror[2]->SetMarkerStyle() ;
      hflat_0lepldp_ratio_ave_withRMSerror[2]->SetLineColor(4) ;
      hflat_0lepldp_ratio_ave_withRMSerror[2]->DrawCopy("samee1") ;
      hflat_0lepldp_ratio_ave[2]->DrawCopy("samee1") ;
      hflat_0lepldp_ratio_model[2] -> SetLineColor(2) ;
      hflat_0lepldp_ratio_model[2] -> SetLineWidth(2) ;
      hflat_0lepldp_ratio_model[2] -> Draw("histsame") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;


      cqcd4->Update() ; cqcd4->Draw() ;

      if ( savePlots ) {
         cqcd4 -> SaveAs("outputfiles/mcclosure4-qcd-averatio.png") ;
         cqcd4 -> SaveAs("outputfiles/mcclosure4-qcd-averatio.pdf") ;
      }












     //--- Draw the QCD scale factor with MC stat error and total error.

      TCanvas* cqcd5 = (TCanvas*) gDirectory->FindObject("cqcd5") ;
      if ( cqcd5 == 0x0 ) {
         cqcd5 = new TCanvas("cqcd5", "qcd study", 1700, 600 ) ;
      }

      gStyle->SetEndErrorSize(5) ;

      cqcd5->Clear() ;
      cqcd5->Divide(3,1) ;

      cqcd5 -> cd(1) ;
      hflat_scale_factor[0] -> SetTitle("QCD Scale Factor, nb=1") ;
      hflat_scale_factor[0] -> SetMarkerSize(1.2) ;
      hflat_scale_factor[0] -> SetMarkerStyle(20) ;
      hflat_scale_factor[0] -> SetMaximum(2.6) ;
      hflat_scale_factor[0] -> SetMinimum(-0.1) ;
      hflat_scale_factor[0]->DrawCopy("e1") ;
      hflat_scale_factor_withRMSerror[0]->SetMarkerStyle() ;
      hflat_scale_factor_withRMSerror[0]->SetLineColor(4) ;
      hflat_scale_factor_withRMSerror[0]->DrawCopy("samee1") ;
      hflat_scale_factor[0]->DrawCopy("samee1") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd5 -> cd(2) ;
      hflat_scale_factor[1] -> SetTitle("QCD Scale Factor, nb=2") ;
      hflat_scale_factor[1] -> SetMarkerSize(1.2) ;
      hflat_scale_factor[1] -> SetMarkerStyle(20) ;
      hflat_scale_factor[1] -> SetMaximum(2.6) ;
      hflat_scale_factor[1] -> SetMinimum(-0.1) ;
      hflat_scale_factor[1]->DrawCopy("e1") ;
      hflat_scale_factor_withRMSerror[1]->SetMarkerStyle() ;
      hflat_scale_factor_withRMSerror[1]->SetLineColor(4) ;
      hflat_scale_factor_withRMSerror[1]->DrawCopy("samee1") ;
      hflat_scale_factor[1]->DrawCopy("samee1") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;

      cqcd5 -> cd(3) ;
      hflat_scale_factor[2] -> SetTitle("QCD Scale Factor, nb>=3") ;
      hflat_scale_factor[2] -> SetMarkerSize(1.2) ;
      hflat_scale_factor[2] -> SetMarkerStyle(20) ;
      hflat_scale_factor[2] -> SetMaximum(2.6) ;
      hflat_scale_factor[2] -> SetMinimum(-0.1) ;
      hflat_scale_factor[2]->DrawCopy("e1") ;
      hflat_scale_factor_withRMSerror[2]->SetMarkerStyle() ;
      hflat_scale_factor_withRMSerror[2]->SetLineColor(4) ;
      hflat_scale_factor_withRMSerror[2]->DrawCopy("samee1") ;
      hflat_scale_factor[2]->DrawCopy("samee1") ;
      line->DrawLine(0.5,0,nbins+0.5,0) ;
      gPad->SetGridx(1) ;
      gPad->SetGridy(1) ;


      cqcd5->Update() ; cqcd5->Draw() ;

      if ( savePlots ) {
         char filename[1000] ;
         sprintf( filename, "outputfiles/mcclosure4-qcd-scalefactor-model%d.png", qcdModelIndex ) ;
         cqcd5 -> SaveAs( filename ) ;
         sprintf( filename, "outputfiles/mcclosure4-qcd-scalefactor-model%d.pdf", qcdModelIndex ) ;
         cqcd5 -> SaveAs( filename ) ;
      }



     } //-- scoping bracket for QCD part.

























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













     //--- Write MC values for background fit parameters to a file.

      TString mcvalfile( datfile ) ;
      char newend[100] ;
      sprintf( newend, ".mc-fitpar-vals-qcdmodel%d", qcdModelIndex ) ;
      mcvalfile.ReplaceAll( ".dat", newend ) ;
      printf("\n\n Creating %s\n\n", mcvalfile.Data() ) ;
      FILE* mcval_file = fopen( mcvalfile.Data(),"w") ;


      fprintf( mcval_file, "ttwj_0lep1lep_ratio  %5.3f\n", simpleAveR_0over1 ) ;

      if ( qcdModelIndex == 2 ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            fprintf( mcval_file, "qcd_0lepLDP_ratio_H%d  %5.3f\n", hbi+1, htbin_0lepldp_ratio[hbi] ) ;
         } // hbi.
      } else if ( qcdModelIndex == 3 ) {
         fprintf( mcval_file, "qcd_0lepLDP_ratio  %5.3f\n", htbin_0lepldp_ratio_global ) ;
      } else if ( qcdModelIndex == 4 ) {
         for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
            fprintf( mcval_file, "qcd_0lepLDP_ratio_H%d  %5.3f\n", hbi+1, htbin_0lepldp_ratio[hbi] ) ;
         } // hbi.
         for ( int mbi=1; mbi<nBinsHT; mbi++ ) {
            fprintf( mcval_file, "SFqcd_met%d  %5.3f\n", mbi+1, fit_SFqcd_MET[mbi] ) ;
         } // mbi.
         for ( int bbi=1; bbi<nBinsHT; bbi++ ) {
            fprintf( mcval_file, "SFqcd_nb%d  %5.3f\n", bbi+1, fit_SFqcd_nb[bbi] ) ;
         } // bbi.
      }

      fclose( mcval_file ) ;













     //--- Reset the 0lep observables to the expected values from the model, if requested.

      if ( setExpected0lepObs ) {

         printf("\n\n\n ============== Resetting 0lep observables to expectations from model in %s ================\n\n", datfile ) ;

         char pname[1000] ;

         sprintf( pname, "Z_ee_pur" ) ;
         float Z_ee_pur ;
         getFileValue( datfile, pname, Z_ee_pur ) ;

         sprintf( pname, "Z_mm_pur" ) ;
         float Z_mm_pur ;
         getFileValue( datfile, pname, Z_mm_pur ) ;

         sprintf( pname, "Z_ee_eff" ) ;
         float Z_ee_eff ;
         getFileValue( datfile, pname, Z_ee_eff ) ;

         sprintf( pname, "Z_mm_eff" ) ;
         float Z_mm_eff ;
         getFileValue( datfile, pname, Z_mm_eff ) ;

         sprintf( pname, "acc_Zee_M1" ) ;
         float acc_Zee ;
         getFileValue( datfile, pname, acc_Zee ) ;

         sprintf( pname, "acc_Zmm_M1" ) ;
         float acc_Zmm ;
         getFileValue( datfile, pname, acc_Zmm ) ;


         float Zee_factor[10] ;
         float Zmm_factor[10] ;
         for ( int bbi=0; bbi<nBinsMET; bbi++ ) {

            sprintf( pname, "knn_%db", bbi+1 ) ;
            float knn ;
            getFileValue( datfile, pname, knn ) ;

            Zee_factor[bbi] = (5.94 * Z_ee_pur * knn ) / ( acc_Zee * Z_ee_eff ) ;
            Zmm_factor[bbi] = (5.94 * Z_mm_pur * knn ) / ( acc_Zmm * Z_mm_eff ) ;

         } // bbi.

         for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
            for ( int hbi=0; hbi<nBinsMET; hbi++ ) {

               sprintf( pname, "N_Zee_M%d_H%d", mbi+1, hbi+1 ) ;
               float NZee ;
               getFileValue( datfile, pname, NZee ) ;

               sprintf( pname, "N_Zmm_M%d_H%d", mbi+1, hbi+1 ) ;
               float NZmm ;
               getFileValue( datfile, pname, NZmm ) ;

               for ( int bbi=0; bbi<nBinsMET; bbi++ ) {


                  sprintf( pname, "N_1lep_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  float N1lep ;
                  getFileValue( datfile, pname, N1lep ) ;

                  sprintf( pname, "N_ldp_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  float Nldp ;
                  getFileValue( datfile, pname, Nldp ) ;

                  sprintf( pname, "sf_qcd_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  float sf_qcd ;
                  getFileValue( datfile, pname, sf_qcd ) ;

                  sprintf( pname, "sf_ttwj_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  float sf_ttwj ;
                  getFileValue( datfile, pname, sf_ttwj ) ;

                  sprintf( pname, "ttwj_mc_ldpover0lep_ratio_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  float ttwj_mc_ldpover0lep_ratio ;
                  getFileValue( datfile, pname, ttwj_mc_ldpover0lep_ratio ) ;

                  sprintf( pname, "znn_mc_ldpover0lep_ratio_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  float znn_mc_ldpover0lep_ratio ;
                  getFileValue( datfile, pname, znn_mc_ldpover0lep_ratio ) ;


                  double exp_0lep_ttwj, exp_0lep_znn, exp_0lep_qcd ;

                  exp_0lep_ttwj = sf_ttwj * N1lep * simpleAveR_0over1 ;

                  exp_0lep_znn = 0.5 * ( NZee * Zee_factor[bbi] + NZmm * Zmm_factor[bbi] ) ;

                  double qcdratio(0.) ;

                  if ( qcdModelIndex == 2 ) {
                     qcdratio = htbin_0lepldp_ratio[hbi] ;
                  } else if ( qcdModelIndex == 3 ) {
                     qcdratio = htbin_0lepldp_ratio_global ;
                  } else if ( qcdModelIndex == 4 ) {
                     qcdratio = fit_Rqcd_HT[hbi] * fit_SFqcd_MET[mbi] * fit_SFqcd_nb[bbi] ;
                  }
                  exp_0lep_qcd = sf_qcd * ( Nldp - exp_0lep_ttwj * ttwj_mc_ldpover0lep_ratio - exp_0lep_znn * znn_mc_ldpover0lep_ratio ) * qcdratio ;

                  double newN0lep = exp_0lep_ttwj + exp_0lep_znn + exp_0lep_qcd ;


                  sprintf( pname, "N_0lep_M%d_H%d_%db", mbi+1, hbi+1, bbi+1 ) ;
                  float oldN0lep ;
                  getFileValue( datfile, pname, oldN0lep ) ;

                  updateFileValue( datfile, pname, newN0lep ) ;

                  printf("\n M%d,H%d,b%d : ttwj = %7.1f, qcd = %7.1f, Znn = %7.1f, expected model total = %7.1f  (MC val %7.1f)\n\n",
                     mbi+1, hbi+1, bbi+1, exp_0lep_ttwj, exp_0lep_qcd, exp_0lep_znn, newN0lep, oldN0lep ) ;


               } // bbi.
            } // hbi.
         } // mbi.


      } // setExpected0lepObs?









   } // mcclosure4


  //==========================================================================================
  
  
    TH1F* bookHist(const char* hname, const char* htitle, int sampleIndex ) {
  
       int nbins = nBinsMET*(nBinsHT+1) + 1 ;
  
       TH1F* retVal = new TH1F( hname, htitle, nbins, 0.5+0.05*(sampleIndex-5), nbins+0.5+0.05*(sampleIndex-5) ) ;
       TAxis* xaxis = retVal->GetXaxis() ;
  
       for ( int mbi=0; mbi<nBinsMET; mbi++ ) {
          for ( int hbi=0; hbi<nBinsHT; hbi++ ) {
             int histbin = 1 + (nBinsHT+1)*mbi + hbi + 1 ;
             char binlabel[1000] ;
             sprintf( binlabel, "M%d_H%d", mbi+1, hbi+1 ) ;
             xaxis->SetBinLabel( histbin, binlabel ) ;
          } // hbi.
       } // mbi.
  
       retVal->SetLabelSize(0.055,"x") ;
       xaxis->LabelsOption("v") ;
  
       return retVal ;
  
    }
  
  //==========================================================================================


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

