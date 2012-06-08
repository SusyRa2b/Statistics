
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TAxis.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TText.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TDirectory.h"

#include <string.h>
#include <sstream>

#include <iostream>

 using std::cout ;


   void loadHist1(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;

//----------------

   void draw_mctruth_plots3D( const char* histfile = "gi-plots.root",
                              bool logy=false,
                              bool doNorm=false,
                              double normmax=2.0,
                              int metgroupzoom=1 ) {


     if ( doNorm ) { logy = false ; }

     gStyle->SetOptStat(0) ;
     gStyle->SetPadTopMargin(0.03) ;
     gStyle->SetPadBottomMargin(0.30) ;
     gStyle->SetPadRightMargin(0.10) ;
     gStyle->SetTitleX(0.95) ;
     gStyle->SetTitleAlign(33) ;

     gDirectory->Delete("hmctruth*") ;

     loadHist1( histfile ) ;
     gDirectory->ls() ;

     if ( logy ) {
        gStyle->SetOptLogy(1) ;
     } else {
        gStyle->SetOptLogy(0) ;
     }


     TH1F* hmctruth_allsm_0lep_1b = (TH1F*) gDirectory->FindObject("hmctruth_allsm_0lep_1b") ;
     TH1F* hmctruth_ttwj_0lep_1b = (TH1F*) gDirectory->FindObject("hmctruth_ttwj_0lep_1b") ;
     TH1F* hmctruth_qcd_0lep_1b  = (TH1F*) gDirectory->FindObject("hmctruth_qcd_0lep_1b") ;
     TH1F* hmctruth_znn_0lep_1b  = (TH1F*) gDirectory->FindObject("hmctruth_znn_0lep_1b") ;
     if ( hmctruth_allsm_0lep_1b == 0 ) { printf("\n\n*** hmctruth_allsm_0lep_1b missing.\n\n") ; return ; }
     if ( hmctruth_ttwj_0lep_1b == 0 ) { printf("\n\n*** hmctruth_ttwj_0lep_1b missing.\n\n") ; return ; }
     if ( hmctruth_qcd_0lep_1b == 0 ) { printf("\n\n*** hmctruth_qcd_0lep_1b missing.\n\n") ; return ; }
     if ( hmctruth_znn_0lep_1b == 0 ) { printf("\n\n*** hmctruth_znn_0lep_1b missing.\n\n") ; return ; }

     TH1F* hmctruth_allsm_0lep_2b = (TH1F*) gDirectory->FindObject("hmctruth_allsm_0lep_2b") ;
     TH1F* hmctruth_ttwj_0lep_2b = (TH1F*) gDirectory->FindObject("hmctruth_ttwj_0lep_2b") ;
     TH1F* hmctruth_qcd_0lep_2b  = (TH1F*) gDirectory->FindObject("hmctruth_qcd_0lep_2b") ;
     TH1F* hmctruth_znn_0lep_2b  = (TH1F*) gDirectory->FindObject("hmctruth_znn_0lep_2b") ;
     if ( hmctruth_allsm_0lep_2b == 0 ) { printf("\n\n*** hmctruth_allsm_0lep_2b missing.\n\n") ; return ; }
     if ( hmctruth_ttwj_0lep_2b == 0 ) { printf("\n\n*** hmctruth_ttwj_0lep_2b missing.\n\n") ; return ; }
     if ( hmctruth_qcd_0lep_2b == 0 ) { printf("\n\n*** hmctruth_qcd_0lep_2b missing.\n\n") ; return ; }
     if ( hmctruth_znn_0lep_2b == 0 ) { printf("\n\n*** hmctruth_znn_0lep_2b missing.\n\n") ; return ; }

     TH1F* hmctruth_allsm_0lep_3b =  (TH1F*) gDirectory->FindObject("hmctruth_allsm_0lep_3b") ;
     TH1F* hmctruth_ttwj_0lep_3b =  (TH1F*) gDirectory->FindObject("hmctruth_ttwj_0lep_3b") ;
     TH1F* hmctruth_qcd_0lep_3b  =  (TH1F*) gDirectory->FindObject("hmctruth_qcd_0lep_3b") ;
     TH1F* hmctruth_znn_0lep_3b  =  (TH1F*) gDirectory->FindObject("hmctruth_znn_0lep_3b") ;
     if ( hmctruth_allsm_0lep_3b == 0 ) { printf("\n\n*** hmctruth_allsm_0lep_3b missing.\n\n") ; return ; }
     if ( hmctruth_ttwj_0lep_3b == 0 ) { printf("\n\n*** hmctruth_ttwj_0lep_3b missing.\n\n") ; return ; }
     if ( hmctruth_qcd_0lep_3b == 0 ) { printf("\n\n*** hmctruth_qcd_0lep_3b missing.\n\n") ; return ; }
     if ( hmctruth_znn_0lep_3b == 0 ) { printf("\n\n*** hmctruth_znn_0lep_3b missing.\n\n") ; return ; }



     TH1F* hmctruth_allsm_1lep_1b =  (TH1F*) gDirectory->FindObject("hmctruth_allsm_1lep_1b") ;
     TH1F* hmctruth_ttwj_1lep_1b =  (TH1F*) gDirectory->FindObject("hmctruth_ttwj_1lep_1b") ;
     TH1F* hmctruth_qcd_1lep_1b  =  (TH1F*) gDirectory->FindObject("hmctruth_qcd_1lep_1b") ;
     TH1F* hmctruth_znn_1lep_1b  =  (TH1F*) gDirectory->FindObject("hmctruth_znn_1lep_1b") ;

     TH1F* hmctruth_allsm_1lep_2b =  (TH1F*) gDirectory->FindObject("hmctruth_allsm_1lep_2b") ;
     TH1F* hmctruth_ttwj_1lep_2b =  (TH1F*) gDirectory->FindObject("hmctruth_ttwj_1lep_2b") ;
     TH1F* hmctruth_qcd_1lep_2b  =  (TH1F*) gDirectory->FindObject("hmctruth_qcd_1lep_2b") ;
     TH1F* hmctruth_znn_1lep_2b  =  (TH1F*) gDirectory->FindObject("hmctruth_znn_1lep_2b") ;

     TH1F* hmctruth_allsm_1lep_3b =  (TH1F*) gDirectory->FindObject("hmctruth_allsm_1lep_3b") ;
     TH1F* hmctruth_ttwj_1lep_3b =  (TH1F*) gDirectory->FindObject("hmctruth_ttwj_1lep_3b") ;
     TH1F* hmctruth_qcd_1lep_3b  =  (TH1F*) gDirectory->FindObject("hmctruth_qcd_1lep_3b") ;
     TH1F* hmctruth_znn_1lep_3b  =  (TH1F*) gDirectory->FindObject("hmctruth_znn_1lep_3b") ;



     TH1F* hmctruth_allsm_ldp_1b =  (TH1F*) gDirectory->FindObject("hmctruth_allsm_ldp_1b") ;
     TH1F* hmctruth_ttwj_ldp_1b =  (TH1F*) gDirectory->FindObject("hmctruth_ttwj_ldp_1b") ;
     TH1F* hmctruth_qcd_ldp_1b  =  (TH1F*) gDirectory->FindObject("hmctruth_qcd_ldp_1b") ;
     TH1F* hmctruth_znn_ldp_1b  =  (TH1F*) gDirectory->FindObject("hmctruth_znn_ldp_1b") ;

     TH1F* hmctruth_allsm_ldp_2b =  (TH1F*) gDirectory->FindObject("hmctruth_allsm_ldp_2b") ;
     TH1F* hmctruth_ttwj_ldp_2b =  (TH1F*) gDirectory->FindObject("hmctruth_ttwj_ldp_2b") ;
     TH1F* hmctruth_qcd_ldp_2b  =  (TH1F*) gDirectory->FindObject("hmctruth_qcd_ldp_2b") ;
     TH1F* hmctruth_znn_ldp_2b  =  (TH1F*) gDirectory->FindObject("hmctruth_znn_ldp_2b") ;

     TH1F* hmctruth_allsm_ldp_3b =  (TH1F*) gDirectory->FindObject("hmctruth_allsm_ldp_3b") ;
     TH1F* hmctruth_ttwj_ldp_3b =  (TH1F*) gDirectory->FindObject("hmctruth_ttwj_ldp_3b") ;
     TH1F* hmctruth_qcd_ldp_3b  =  (TH1F*) gDirectory->FindObject("hmctruth_qcd_ldp_3b") ;
     TH1F* hmctruth_znn_ldp_3b  =  (TH1F*) gDirectory->FindObject("hmctruth_znn_ldp_3b") ;


     TH1F* hmctruth_fit_zee_1b  =  (TH1F*) gDirectory->FindObject("hmctruth_fit_zee_1b") ;
     TH1F* hmctruth_fit_zmm_1b  =  (TH1F*) gDirectory->FindObject("hmctruth_fit_zmm_1b") ;

     TH1F* hmctruth_np   =  (TH1F*) gDirectory->FindObject("hmctruth_np") ;




     TH1F* hmctruth_susy_0lep_1b = (TH1F*) gDirectory->FindObject("hmctruth_susy_0lep_1b") ;
     TH1F* hmctruth_susy_0lep_2b = (TH1F*) gDirectory->FindObject("hmctruth_susy_0lep_2b") ;
     TH1F* hmctruth_susy_0lep_3b =  (TH1F*) gDirectory->FindObject("hmctruth_susy_0lep_3b") ;
     TH1F* hmctruth_susy_1lep_1b =  (TH1F*) gDirectory->FindObject("hmctruth_susy_1lep_1b") ;
     TH1F* hmctruth_susy_1lep_2b =  (TH1F*) gDirectory->FindObject("hmctruth_susy_1lep_2b") ;
     TH1F* hmctruth_susy_1lep_3b =  (TH1F*) gDirectory->FindObject("hmctruth_susy_1lep_3b") ;
     TH1F* hmctruth_susy_ldp_1b =  (TH1F*) gDirectory->FindObject("hmctruth_susy_ldp_1b") ;
     TH1F* hmctruth_susy_ldp_2b =  (TH1F*) gDirectory->FindObject("hmctruth_susy_ldp_2b") ;
     TH1F* hmctruth_susy_ldp_3b =  (TH1F*) gDirectory->FindObject("hmctruth_susy_ldp_3b") ;
     if ( hmctruth_susy_0lep_1b == 0 ) { printf("\n\n*** hmctruth_susy_0lep_1b missing.\n\n") ; return ; }
     if ( hmctruth_susy_0lep_2b == 0 ) { printf("\n\n*** hmctruth_susy_0lep_2b missing.\n\n") ; return ; }
     if ( hmctruth_susy_0lep_3b == 0 ) { printf("\n\n*** hmctruth_susy_0lep_3b missing.\n\n") ; return ; }




     int ncomp(5) ;
     char compname[5][100] = { "allsm", "susy", "ttwj", "qcd", "znn" } ;

     int nsel(3) ;
     char selname[3][100] = { "0lep", "1lep", "ldp" } ;

     int nbtagmult(3) ;
     char btagmultname[3][100] = { "1b", "2b", "3b" } ;




     if ( doNorm ) {

       //--- Divide all bins by the overall model prediction.
       //    This is slightly different than the way its done
       //    in the other one, where the normalization is the allsm entries.
       //    Doing it this way because we will have several observables
       //    with no events in the allsm.  Model should always be non-zero(?)
       //

        printf("\n\n Renormalizing histograms.\n\n") ;

        int nbins = hmctruth_allsm_0lep_1b -> GetNbinsX() ;

        for ( int bini=1; bini<=nbins; bini++ ) {
           for ( int seli=0; seli<nsel; seli++ ) {
              for ( int nbji=0; nbji<nbtagmult; nbji++ ) {

                 double modelsum(0.) ;

                 for ( int ci=1; ci<ncomp; ci++ ) {

                    char hname[100] ;
                    sprintf( hname, "hmctruth_%s_%s_%s", compname[ci], selname[seli], btagmultname[nbji] ) ;

                    TH1F* hist = (TH1F*) gDirectory->FindObject( hname ) ;
                    if ( hist == 0 ) { printf("\n\n\n *** Can't find histogram %s\n\n", hname ) ; return ; }

                    modelsum += hist->GetBinContent( bini ) ;


                 } // ci


                 if ( modelsum > 0. ) {
                    for ( int ci=0; ci<ncomp; ci++ ) {

                       char hname[100] ;
                       sprintf( hname, "hmctruth_%s_%s_%s", compname[ci], selname[seli], btagmultname[nbji] ) ;


                       TH1F* hist = (TH1F*) gDirectory->FindObject( hname ) ;
                       if ( hist == 0 ) { printf("\n\n\n *** Can't find histogram %s\n\n", hname ) ; return ; }

                       if ( ci==0 ) { // this is the allsm
                          hist -> SetBinError( bini, (hist->GetBinError(bini))/modelsum ) ;
                       }
                       hist -> SetBinContent( bini, (hist->GetBinContent(bini))/modelsum ) ;

                       ///// printf(" %d (%d,%d,%d) : %4.2f  %s\n", bini, seli, nbji, ci, hist->GetBinContent( bini ), hname ) ;

                       hist -> SetMaximum(normmax) ;

                    } // ci
                    /////// printf("\n") ;
                 }


              } // nbji
           } // seli

           printf("\n\n debug 1 : after loops.\n\n") ;

           if ( hmctruth_fit_zee_1b->GetBinContent( bini ) > 0. ) {
              hmctruth_fit_zee_1b ->SetBinContent( bini, 1. ) ;
           }

           if ( hmctruth_fit_zmm_1b->GetBinContent( bini ) > 0. ) {
              hmctruth_fit_zmm_1b ->SetBinContent( bini, 1. ) ;
           }

           hmctruth_fit_zee_1b -> SetMaximum(normmax) ;
           hmctruth_fit_zmm_1b -> SetMaximum(normmax) ;

        } // bini


     } // doNorm?




     hmctruth_fit_zee_1b->SetFillColor(kGreen-3) ;
     hmctruth_fit_zmm_1b->SetFillColor(kGreen-3) ;



     hmctruth_ttwj_0lep_1b  -> SetFillColor(kBlue-9) ;
     hmctruth_qcd_0lep_1b   -> SetFillColor(2) ;
     hmctruth_znn_0lep_1b   -> SetFillColor(kGreen-3) ;
     hmctruth_susy_0lep_1b  -> SetFillColor(6) ;
     
     hmctruth_ttwj_0lep_2b  -> SetFillColor(kBlue-9) ;
     hmctruth_qcd_0lep_2b   -> SetFillColor(2) ;
     hmctruth_znn_0lep_2b   -> SetFillColor(kGreen-3) ;
     hmctruth_susy_0lep_2b  -> SetFillColor(6) ;
     
     hmctruth_ttwj_0lep_3b  -> SetFillColor(kBlue-9) ;
     hmctruth_qcd_0lep_3b   -> SetFillColor(2) ;
     hmctruth_znn_0lep_3b   -> SetFillColor(kGreen-3) ;
     hmctruth_susy_0lep_3b  -> SetFillColor(6) ;



     
     hmctruth_ttwj_1lep_1b  -> SetFillColor(kBlue-9) ;
     hmctruth_qcd_1lep_1b   -> SetFillColor(2) ;
     hmctruth_znn_1lep_1b   -> SetFillColor(kGreen-3) ;
     hmctruth_susy_1lep_1b  -> SetFillColor(6) ;

     hmctruth_ttwj_1lep_2b  -> SetFillColor(kBlue-9) ;
     hmctruth_qcd_1lep_2b   -> SetFillColor(2) ;
     hmctruth_znn_1lep_2b   -> SetFillColor(kGreen-3) ;
     hmctruth_susy_1lep_2b  -> SetFillColor(6) ;

     hmctruth_ttwj_1lep_3b  -> SetFillColor(kBlue-9) ;
     hmctruth_qcd_1lep_3b   -> SetFillColor(2) ;
     hmctruth_znn_1lep_3b   -> SetFillColor(kGreen-3) ;
     hmctruth_susy_1lep_3b  -> SetFillColor(6) ;




     
     hmctruth_ttwj_ldp_1b  -> SetFillColor(kBlue-9) ;
     hmctruth_qcd_ldp_1b   -> SetFillColor(2) ;
     hmctruth_znn_ldp_1b   -> SetFillColor(kGreen-3) ;
     hmctruth_susy_ldp_1b  -> SetFillColor(6) ;

     hmctruth_ttwj_ldp_2b  -> SetFillColor(kBlue-9) ;
     hmctruth_qcd_ldp_2b   -> SetFillColor(2) ;
     hmctruth_znn_ldp_2b   -> SetFillColor(kGreen-3) ;
     hmctruth_susy_ldp_2b  -> SetFillColor(6) ;

     hmctruth_ttwj_ldp_3b  -> SetFillColor(kBlue-9) ;
     hmctruth_qcd_ldp_3b   -> SetFillColor(2) ;
     hmctruth_znn_ldp_3b   -> SetFillColor(kGreen-3) ;
     hmctruth_susy_ldp_3b  -> SetFillColor(6) ;






     printf("\n\n Making stacks...\n") ; cout << flush ;

     THStack* hmctruth_fit_0lep_1b = new THStack( "hmctruth_fit_0lep_1b", "RA2b likelihood fit results, fit" ) ;
     THStack* hmctruth_fit_0lep_2b = new THStack( "hmctruth_fit_0lep_2b", "RA2b likelihood fit results, fit" ) ;
     THStack* hmctruth_fit_0lep_3b = new THStack( "hmctruth_fit_0lep_3b", "RA2b likelihood fit results, fit" ) ;

     THStack* hmctruth_fit_1lep_1b = new THStack( "hmctruth_fit_1lep_1b", "RA2b likelihood fit results, fit" ) ;
     THStack* hmctruth_fit_1lep_2b = new THStack( "hmctruth_fit_1lep_2b", "RA2b likelihood fit results, fit" ) ;
     THStack* hmctruth_fit_1lep_3b = new THStack( "hmctruth_fit_1lep_3b", "RA2b likelihood fit results, fit" ) ;

     THStack* hmctruth_fit_ldp_1b  = new THStack( "hmctruth_fit_ldp_1b",  "RA2b likelihood fit results, fit" ) ;
     THStack* hmctruth_fit_ldp_2b  = new THStack( "hmctruth_fit_ldp_2b",  "RA2b likelihood fit results, fit" ) ;
     THStack* hmctruth_fit_ldp_3b  = new THStack( "hmctruth_fit_ldp_3b",  "RA2b likelihood fit results, fit" ) ;

     hmctruth_fit_0lep_1b->Add( hmctruth_znn_0lep_1b ) ;
     hmctruth_fit_0lep_1b->Add( hmctruth_qcd_0lep_1b ) ;
     hmctruth_fit_0lep_1b->Add( hmctruth_ttwj_0lep_1b ) ;
     hmctruth_fit_0lep_1b->Add( hmctruth_susy_0lep_1b ) ;

     hmctruth_fit_0lep_2b->Add( hmctruth_znn_0lep_2b ) ;
     hmctruth_fit_0lep_2b->Add( hmctruth_qcd_0lep_2b ) ;
     hmctruth_fit_0lep_2b->Add( hmctruth_ttwj_0lep_2b ) ;
     hmctruth_fit_0lep_2b->Add( hmctruth_susy_0lep_2b ) ;

     hmctruth_fit_0lep_3b->Add( hmctruth_znn_0lep_3b ) ;
     hmctruth_fit_0lep_3b->Add( hmctruth_qcd_0lep_3b ) ;
     hmctruth_fit_0lep_3b->Add( hmctruth_ttwj_0lep_3b ) ;
     hmctruth_fit_0lep_3b->Add( hmctruth_susy_0lep_3b ) ;



     hmctruth_fit_1lep_1b->Add( hmctruth_znn_1lep_1b ) ;
     hmctruth_fit_1lep_1b->Add( hmctruth_qcd_1lep_1b ) ;
     hmctruth_fit_1lep_1b->Add( hmctruth_ttwj_1lep_1b ) ;
     hmctruth_fit_1lep_1b->Add( hmctruth_susy_1lep_1b ) ;

     hmctruth_fit_1lep_2b->Add( hmctruth_znn_1lep_2b ) ;
     hmctruth_fit_1lep_2b->Add( hmctruth_qcd_1lep_2b ) ;
     hmctruth_fit_1lep_2b->Add( hmctruth_ttwj_1lep_2b ) ;
     hmctruth_fit_1lep_2b->Add( hmctruth_susy_1lep_2b ) ;

     hmctruth_fit_1lep_3b->Add( hmctruth_znn_1lep_3b ) ;
     hmctruth_fit_1lep_3b->Add( hmctruth_qcd_1lep_3b ) ;
     hmctruth_fit_1lep_3b->Add( hmctruth_ttwj_1lep_3b ) ;
     hmctruth_fit_1lep_3b->Add( hmctruth_susy_1lep_3b ) ;




     hmctruth_fit_ldp_1b->Add( hmctruth_znn_ldp_1b ) ;
     hmctruth_fit_ldp_1b->Add( hmctruth_qcd_ldp_1b ) ;
     hmctruth_fit_ldp_1b->Add( hmctruth_ttwj_ldp_1b ) ;
     hmctruth_fit_ldp_1b->Add( hmctruth_susy_ldp_1b ) ;

     hmctruth_fit_ldp_2b->Add( hmctruth_znn_ldp_2b ) ;
     hmctruth_fit_ldp_2b->Add( hmctruth_qcd_ldp_2b ) ;
     hmctruth_fit_ldp_2b->Add( hmctruth_ttwj_ldp_2b ) ;
     hmctruth_fit_ldp_2b->Add( hmctruth_susy_ldp_2b ) ;

     hmctruth_fit_ldp_3b->Add( hmctruth_znn_ldp_3b ) ;
     hmctruth_fit_ldp_3b->Add( hmctruth_qcd_ldp_3b ) ;
     hmctruth_fit_ldp_3b->Add( hmctruth_ttwj_ldp_3b ) ;
     hmctruth_fit_ldp_3b->Add( hmctruth_susy_ldp_3b ) ;



     printf("\n\n Done making stacks.\n\n") ; cout << flush ;

     printf("\n Making legend...\n") ; cout << flush ;
     TLegend* legend = new TLegend(0.4,0.35,0.7,0.85) ;

     legend->AddEntry( hmctruth_allsm_0lep_1b, "allsm" ) ;
     legend->AddEntry( hmctruth_susy_0lep_1b, "SUSY" ) ;
     legend->AddEntry( hmctruth_ttwj_0lep_1b, "ttwj" ) ;
     legend->AddEntry( hmctruth_qcd_0lep_1b,  "QCD" ) ;
     legend->AddEntry( hmctruth_znn_0lep_1b,  "Znunu" ) ;
     legend->AddEntry( hmctruth_np,           "Eff PG" ) ;

     printf("\n\n Done making legend.\n\n\n") ; cout << flush ;

     TCanvas* cmctruth = (TCanvas*) gDirectory->FindObject("cmctruth") ;
     if ( cmctruth==0 ) {
        printf("\n\n Creating cmctruth canvas.\n\n") ;
        cmctruth = new TCanvas("cmctruth","RA2b fit quality", 850, 1000 ) ;
     } else {
        printf("\n\n Found existing cmctruth canvas.\n\n") ;
        cmctruth->Clear() ;
     }





     if ( metgroupzoom>1 && !logy ) {

        int nhistbins = hmctruth_ttwj_0lep_1b->GetNbinsX() ;

        int nBinsHT(0) ;
        for ( int bi=2; bi<nhistbins; bi++ ) {
           if ( hmctruth_ttwj_0lep_1b->GetBinContent( bi ) <= 0. ) break ;
           nBinsHT++ ;
        }

        int zoomRefBin = 1 + (nBinsHT+1)*(metgroupzoom-1) + 1 ;

        printf("\n\n Number of HT bins=%d, met group ref bin=%d\n\n", nBinsHT, zoomRefBin ) ;

        double maxSF(1.3) ;

        hmctruth_allsm_0lep_1b->SetMaximum( maxSF*(hmctruth_allsm_0lep_1b->GetBinContent( zoomRefBin )) ) ;
        hmctruth_allsm_0lep_2b->SetMaximum( maxSF*(hmctruth_allsm_0lep_2b->GetBinContent( zoomRefBin )) ) ;
        hmctruth_allsm_0lep_3b->SetMaximum( maxSF*(hmctruth_allsm_0lep_3b->GetBinContent( zoomRefBin )) ) ;

        hmctruth_allsm_1lep_1b->SetMaximum( maxSF*(hmctruth_allsm_1lep_1b->GetBinContent( zoomRefBin )) ) ;
        hmctruth_allsm_1lep_2b->SetMaximum( maxSF*(hmctruth_allsm_1lep_2b->GetBinContent( zoomRefBin )) ) ;
        hmctruth_allsm_1lep_3b->SetMaximum( maxSF*(hmctruth_allsm_1lep_3b->GetBinContent( zoomRefBin )) ) ;

        hmctruth_allsm_ldp_1b->SetMaximum( maxSF*(hmctruth_allsm_ldp_1b->GetBinContent( zoomRefBin+2 )) ) ;
        hmctruth_allsm_ldp_2b->SetMaximum( maxSF*(hmctruth_allsm_ldp_2b->GetBinContent( zoomRefBin+2 )) ) ;
        hmctruth_allsm_ldp_3b->SetMaximum( maxSF*(hmctruth_allsm_ldp_3b->GetBinContent( zoomRefBin+2 )) ) ;

        hmctruth_fit_zee_1b->SetMaximum( maxSF*(hmctruth_fit_zee_1b->GetBinContent( zoomRefBin )) ) ;
        hmctruth_fit_zmm_1b->SetMaximum( maxSF*(hmctruth_fit_zmm_1b->GetBinContent( zoomRefBin )) ) ;

     }

     hmctruth_allsm_0lep_1b->SetLineWidth(2) ;
     hmctruth_allsm_0lep_2b->SetLineWidth(2) ;
     hmctruth_allsm_0lep_3b->SetLineWidth(2) ;

     hmctruth_allsm_1lep_1b->SetLineWidth(2) ;
     hmctruth_allsm_1lep_2b->SetLineWidth(2) ;
     hmctruth_allsm_1lep_3b->SetLineWidth(2) ;

     hmctruth_allsm_ldp_1b->SetLineWidth(2) ;
     hmctruth_allsm_ldp_2b->SetLineWidth(2) ;
     hmctruth_allsm_ldp_3b->SetLineWidth(2) ;


     cmctruth->Divide(3,4);

     gPad->SetTicks(1,0) ;


     printf(" pad 1\n") ; cout << flush ;
     cmctruth->cd(1);
     hmctruth_allsm_0lep_1b->Draw("histpe") ;
     hmctruth_fit_0lep_1b->Draw("histsame") ;
     hmctruth_allsm_0lep_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     printf(" pad 2\n") ; cout << flush ;
     cmctruth->cd(2);
     hmctruth_allsm_0lep_2b->Draw("histpe") ;
     hmctruth_fit_0lep_2b->Draw("histsame") ;
     hmctruth_allsm_0lep_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     printf(" pad 3\n") ; cout << flush ;
     cmctruth->cd(3);
     hmctruth_allsm_0lep_3b->Draw("histpe") ;
     hmctruth_fit_0lep_3b->Draw("histsame") ;
     hmctruth_allsm_0lep_3b->Draw("same") ;
     gPad->SetGridy(1) ;
     


     printf(" pad 4\n") ; cout << flush ;
     cmctruth->cd(4);
     hmctruth_allsm_1lep_1b->Draw("histpe") ;
     hmctruth_fit_1lep_1b->Draw("histsame") ;
     hmctruth_allsm_1lep_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     printf(" pad 5\n") ; cout << flush ;
     cmctruth->cd(5);
     hmctruth_allsm_1lep_2b->Draw("histpe") ;
     hmctruth_fit_1lep_2b->Draw("histsame") ;
     hmctruth_allsm_1lep_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     printf(" pad 6\n") ; cout << flush ;
     cmctruth->cd(6);
     hmctruth_allsm_1lep_3b->Draw("histpe") ;
     hmctruth_fit_1lep_3b->Draw("histsame") ;
     hmctruth_allsm_1lep_3b->Draw("same") ;
     gPad->SetGridy(1) ;




     
     printf(" pad 7\n") ; cout << flush ;
     cmctruth->cd(7);
     hmctruth_allsm_ldp_1b->Draw("histpe") ;
     hmctruth_fit_ldp_1b->Draw("histsame") ;
     hmctruth_allsm_ldp_1b->Draw("same") ;
     gPad->SetGridy(1) ;

     printf(" pad 8\n") ; cout << flush ;
     cmctruth->cd(8);
     hmctruth_allsm_ldp_2b->Draw("histpe") ;
     hmctruth_fit_ldp_2b->Draw("histsame") ;
     hmctruth_allsm_ldp_2b->Draw("same") ;
     gPad->SetGridy(1) ;

     printf(" pad 9\n") ; cout << flush ;
     cmctruth->cd(9);
     hmctruth_allsm_ldp_3b->Draw("histpe") ;
     hmctruth_fit_ldp_3b->Draw("histsame") ;
     hmctruth_allsm_ldp_3b->Draw("same") ;
     gPad->SetGridy(1) ;




     printf(" pad 10\n") ; cout << flush ;
     cmctruth->cd(10) ;
     hmctruth_fit_zee_1b->Draw("hist") ;
     hmctruth_fit_zee_1b->Draw("esame") ;
     gPad->SetGridy(1) ;

     printf(" pad 11\n") ; cout << flush ;
     cmctruth->cd(11) ;
     hmctruth_fit_zmm_1b->Draw("histe") ;
     hmctruth_fit_zmm_1b->Draw("esame") ;
     gPad->SetGridy(1) ;




     cmctruth->cd(12);
     legend->Draw() ;

     cmctruth->Update() ;




     cmctruth->SaveAs("mctruth.gif") ;




   }



//==========================================================================================


void loadHist1(const char* filename, const char* pfx, const char* pat, Bool_t doAdd, Double_t scaleFactor)
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

