
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


   void loadHist(const char* filename="in.root", const char* pfx=0, const char* pat="*", Bool_t doAdd=kFALSE, Double_t scaleFactor=-1.0) ;

//----------------

   void ws_redraw_fitqual_plots3D( const char* histfile = "fitqual-hists-ws-met3-ht3-v1.root",
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

     gDirectory->Delete("hfitqual*") ;

     loadHist( histfile ) ;
     gDirectory->ls() ;

     if ( logy ) {
        gStyle->SetOptLogy(1) ;
     } else {
        gStyle->SetOptLogy(0) ;
     }


     TH1F* hfitqual_data_0lep_1b = (TH1F*) gDirectory->FindObject("hfitqual_data_0lep_1b") ;
     TH1F* hfitqual_susy_0lep_1b = (TH1F*) gDirectory->FindObject("hfitqual_susy_0lep_1b") ;
     TH1F* hfitqual_ttwj_0lep_1b = (TH1F*) gDirectory->FindObject("hfitqual_ttwj_0lep_1b") ;
     TH1F* hfitqual_qcd_0lep_1b  = (TH1F*) gDirectory->FindObject("hfitqual_qcd_0lep_1b") ;
     TH1F* hfitqual_znn_0lep_1b  = (TH1F*) gDirectory->FindObject("hfitqual_znn_0lep_1b") ;
     if ( hfitqual_data_0lep_1b == 0 ) { printf("\n\n*** hfitqual_data_0lep_1b missing.\n\n") ; return ; }
     if ( hfitqual_susy_0lep_1b == 0 ) { printf("\n\n*** hfitqual_susy_0lep_1b missing.\n\n") ; return ; }
     if ( hfitqual_ttwj_0lep_1b == 0 ) { printf("\n\n*** hfitqual_ttwj_0lep_1b missing.\n\n") ; return ; }
     if ( hfitqual_qcd_0lep_1b == 0 ) { printf("\n\n*** hfitqual_qcd_0lep_1b missing.\n\n") ; return ; }
     if ( hfitqual_znn_0lep_1b == 0 ) { printf("\n\n*** hfitqual_znn_0lep_1b missing.\n\n") ; return ; }

     TH1F* hfitqual_data_0lep_2b = (TH1F*) gDirectory->FindObject("hfitqual_data_0lep_2b") ;
     TH1F* hfitqual_susy_0lep_2b = (TH1F*) gDirectory->FindObject("hfitqual_susy_0lep_2b") ;
     TH1F* hfitqual_ttwj_0lep_2b = (TH1F*) gDirectory->FindObject("hfitqual_ttwj_0lep_2b") ;
     TH1F* hfitqual_qcd_0lep_2b  = (TH1F*) gDirectory->FindObject("hfitqual_qcd_0lep_2b") ;
     TH1F* hfitqual_znn_0lep_2b  = (TH1F*) gDirectory->FindObject("hfitqual_znn_0lep_2b") ;
     if ( hfitqual_data_0lep_2b == 0 ) { printf("\n\n*** hfitqual_data_0lep_2b missing.\n\n") ; return ; }
     if ( hfitqual_susy_0lep_2b == 0 ) { printf("\n\n*** hfitqual_susy_0lep_2b missing.\n\n") ; return ; }
     if ( hfitqual_ttwj_0lep_2b == 0 ) { printf("\n\n*** hfitqual_ttwj_0lep_2b missing.\n\n") ; return ; }
     if ( hfitqual_qcd_0lep_2b == 0 ) { printf("\n\n*** hfitqual_qcd_0lep_2b missing.\n\n") ; return ; }
     if ( hfitqual_znn_0lep_2b == 0 ) { printf("\n\n*** hfitqual_znn_0lep_2b missing.\n\n") ; return ; }

     TH1F* hfitqual_data_0lep_3b =  (TH1F*) gDirectory->FindObject("hfitqual_data_0lep_3b") ;
     TH1F* hfitqual_susy_0lep_3b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_0lep_3b") ;
     TH1F* hfitqual_ttwj_0lep_3b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_0lep_3b") ;
     TH1F* hfitqual_qcd_0lep_3b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_0lep_3b") ;
     TH1F* hfitqual_znn_0lep_3b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_0lep_3b") ;
     if ( hfitqual_data_0lep_3b == 0 ) { printf("\n\n*** hfitqual_data_0lep_3b missing.\n\n") ; return ; }
     if ( hfitqual_susy_0lep_3b == 0 ) { printf("\n\n*** hfitqual_susy_0lep_3b missing.\n\n") ; return ; }
     if ( hfitqual_ttwj_0lep_3b == 0 ) { printf("\n\n*** hfitqual_ttwj_0lep_3b missing.\n\n") ; return ; }
     if ( hfitqual_qcd_0lep_3b == 0 ) { printf("\n\n*** hfitqual_qcd_0lep_3b missing.\n\n") ; return ; }
     if ( hfitqual_znn_0lep_3b == 0 ) { printf("\n\n*** hfitqual_znn_0lep_3b missing.\n\n") ; return ; }



     TH1F* hfitqual_data_1lep_1b =  (TH1F*) gDirectory->FindObject("hfitqual_data_1lep_1b") ;
     TH1F* hfitqual_susy_1lep_1b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_1lep_1b") ;
     TH1F* hfitqual_ttwj_1lep_1b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_1lep_1b") ;
     TH1F* hfitqual_qcd_1lep_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_1lep_1b") ;
     TH1F* hfitqual_znn_1lep_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_1lep_1b") ;

     TH1F* hfitqual_data_1lep_2b =  (TH1F*) gDirectory->FindObject("hfitqual_data_1lep_2b") ;
     TH1F* hfitqual_susy_1lep_2b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_1lep_2b") ;
     TH1F* hfitqual_ttwj_1lep_2b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_1lep_2b") ;
     TH1F* hfitqual_qcd_1lep_2b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_1lep_2b") ;
     TH1F* hfitqual_znn_1lep_2b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_1lep_2b") ;

     TH1F* hfitqual_data_1lep_3b =  (TH1F*) gDirectory->FindObject("hfitqual_data_1lep_3b") ;
     TH1F* hfitqual_susy_1lep_3b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_1lep_3b") ;
     TH1F* hfitqual_ttwj_1lep_3b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_1lep_3b") ;
     TH1F* hfitqual_qcd_1lep_3b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_1lep_3b") ;
     TH1F* hfitqual_znn_1lep_3b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_1lep_3b") ;



     TH1F* hfitqual_data_ldp_1b =  (TH1F*) gDirectory->FindObject("hfitqual_data_ldp_1b") ;
     TH1F* hfitqual_susy_ldp_1b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_ldp_1b") ;
     TH1F* hfitqual_ttwj_ldp_1b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_ldp_1b") ;
     TH1F* hfitqual_qcd_ldp_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_ldp_1b") ;
     TH1F* hfitqual_znn_ldp_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_ldp_1b") ;

     TH1F* hfitqual_data_ldp_2b =  (TH1F*) gDirectory->FindObject("hfitqual_data_ldp_2b") ;
     TH1F* hfitqual_susy_ldp_2b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_ldp_2b") ;
     TH1F* hfitqual_ttwj_ldp_2b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_ldp_2b") ;
     TH1F* hfitqual_qcd_ldp_2b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_ldp_2b") ;
     TH1F* hfitqual_znn_ldp_2b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_ldp_2b") ;

     TH1F* hfitqual_data_ldp_3b =  (TH1F*) gDirectory->FindObject("hfitqual_data_ldp_3b") ;
     TH1F* hfitqual_susy_ldp_3b =  (TH1F*) gDirectory->FindObject("hfitqual_susy_ldp_3b") ;
     TH1F* hfitqual_ttwj_ldp_3b =  (TH1F*) gDirectory->FindObject("hfitqual_ttwj_ldp_3b") ;
     TH1F* hfitqual_qcd_ldp_3b  =  (TH1F*) gDirectory->FindObject("hfitqual_qcd_ldp_3b") ;
     TH1F* hfitqual_znn_ldp_3b  =  (TH1F*) gDirectory->FindObject("hfitqual_znn_ldp_3b") ;


     TH1F* hfitqual_data_zee_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_data_zee_1b") ;
     TH1F* hfitqual_data_zmm_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_data_zmm_1b") ;
     TH1F* hfitqual_fit_zee_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_fit_zee_1b") ;
     TH1F* hfitqual_fit_zmm_1b  =  (TH1F*) gDirectory->FindObject("hfitqual_fit_zmm_1b") ;

     TH1F* hfitqual_np   =  (TH1F*) gDirectory->FindObject("hfitqual_np") ;


     int ncomp(5) ;
     char compname[5][100] = { "data", "susy", "ttwj", "qcd", "znn" } ;

     int nsel(3) ;
     char selname[3][100] = { "0lep", "1lep", "ldp" } ;

     int nbtagmult(3) ;
     char btagmultname[3][100] = { "1b", "2b", "3b" } ;




     if ( doNorm ) {

       //--- Divide all bins by the overall model prediction.
       //    This is slightly different than the way its done
       //    in the other one, where the normalization is the data entries.
       //    Doing it this way because we will have several observables
       //    with no events in the data.  Model should always be non-zero(?)
       //

        printf("\n\n Renormalizing histograms.\n\n") ;

        int nbins = hfitqual_data_0lep_1b -> GetNbinsX() ;

        for ( int bini=1; bini<=nbins; bini++ ) {

           for ( int seli=0; seli<nsel; seli++ ) {
              for ( int nbji=0; nbji<nbtagmult; nbji++ ) {

                 double modelsum(0.) ;

                 for ( int ci=1; ci<ncomp; ci++ ) {

                    char hname[100] ;
                    sprintf( hname, "hfitqual_%s_%s_%s", compname[ci], selname[seli], btagmultname[nbji] ) ;

                    TH1F* hist = (TH1F*) gDirectory->FindObject( hname ) ;
                    if ( hist == 0 ) { printf("\n\n\n *** Can't find histogram %s\n\n", hname ) ; return ; }

                    modelsum += hist->GetBinContent( bini ) ;


                 } // ci

                 if ( modelsum > 0. ) {
                    for ( int ci=0; ci<ncomp; ci++ ) {

                       char hname[100] ;
                       sprintf( hname, "hfitqual_%s_%s_%s", compname[ci], selname[seli], btagmultname[nbji] ) ;


                       TH1F* hist = (TH1F*) gDirectory->FindObject( hname ) ;
                       if ( hist == 0 ) { printf("\n\n\n *** Can't find histogram %s\n\n", hname ) ; return ; }

                       if ( ci==0 ) { // this is the data
                          hist -> SetBinError( bini, (hist->GetBinError(bini))/modelsum ) ;
                       }
                       hist -> SetBinContent( bini, (hist->GetBinContent(bini))/modelsum ) ;

                       //// printf(" %d (%d,%d,%d) : %4.2f  %s\n", bini, seli, nbji, ci, hist->GetBinContent( bini ), hname ) ;

                       hist -> SetMaximum(normmax) ;

                    } // ci
                    //// printf("\n") ;
                 }

              } // nbji
           } // seli

           if ( hfitqual_fit_zee_1b->GetBinContent( bini ) > 0. ) {
              hfitqual_data_zee_1b->SetBinError( bini, (hfitqual_data_zee_1b->GetBinError(bini))/(hfitqual_fit_zee_1b->GetBinContent( bini )) ) ;
              hfitqual_data_zee_1b->SetBinContent( bini, (hfitqual_data_zee_1b->GetBinContent(bini))/(hfitqual_fit_zee_1b->GetBinContent( bini )) ) ;
              hfitqual_fit_zee_1b ->SetBinContent( bini, 1. ) ;
           }

           if ( hfitqual_fit_zmm_1b->GetBinContent( bini ) > 0. ) {
              hfitqual_data_zmm_1b->SetBinError( bini, (hfitqual_data_zmm_1b->GetBinError(bini))/(hfitqual_fit_zmm_1b->GetBinContent( bini )) ) ;
              hfitqual_data_zmm_1b->SetBinContent( bini, (hfitqual_data_zmm_1b->GetBinContent(bini))/(hfitqual_fit_zmm_1b->GetBinContent( bini )) ) ;
              hfitqual_fit_zmm_1b ->SetBinContent( bini, 1. ) ;
           }

           hfitqual_data_zee_1b -> SetMaximum(normmax) ;
           hfitqual_data_zmm_1b -> SetMaximum(normmax) ;
           hfitqual_fit_zee_1b -> SetMaximum(normmax) ;
           hfitqual_fit_zmm_1b -> SetMaximum(normmax) ;

        } // bini


     } // doNorm?








     printf("\n\n Making stacks...\n") ; cout << flush ;

     THStack* hfitqual_fit_0lep_1b = new THStack( "hfitqual_fit_0lep_1b", "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_0lep_2b = new THStack( "hfitqual_fit_0lep_2b", "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_0lep_3b = new THStack( "hfitqual_fit_0lep_3b", "RA2b likelihood fit results, fit" ) ;

     THStack* hfitqual_fit_1lep_1b = new THStack( "hfitqual_fit_1lep_1b", "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_1lep_2b = new THStack( "hfitqual_fit_1lep_2b", "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_1lep_3b = new THStack( "hfitqual_fit_1lep_3b", "RA2b likelihood fit results, fit" ) ;

     THStack* hfitqual_fit_ldp_1b  = new THStack( "hfitqual_fit_ldp_1b",  "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_ldp_2b  = new THStack( "hfitqual_fit_ldp_2b",  "RA2b likelihood fit results, fit" ) ;
     THStack* hfitqual_fit_ldp_3b  = new THStack( "hfitqual_fit_ldp_3b",  "RA2b likelihood fit results, fit" ) ;

     hfitqual_fit_0lep_1b->Add( hfitqual_znn_0lep_1b ) ;
     hfitqual_fit_0lep_1b->Add( hfitqual_qcd_0lep_1b ) ;
     hfitqual_fit_0lep_1b->Add( hfitqual_ttwj_0lep_1b ) ;
     hfitqual_fit_0lep_1b->Add( hfitqual_susy_0lep_1b ) ;

     hfitqual_fit_0lep_2b->Add( hfitqual_znn_0lep_2b ) ;
     hfitqual_fit_0lep_2b->Add( hfitqual_qcd_0lep_2b ) ;
     hfitqual_fit_0lep_2b->Add( hfitqual_ttwj_0lep_2b ) ;
     hfitqual_fit_0lep_2b->Add( hfitqual_susy_0lep_2b ) ;

     hfitqual_fit_0lep_3b->Add( hfitqual_znn_0lep_3b ) ;
     hfitqual_fit_0lep_3b->Add( hfitqual_qcd_0lep_3b ) ;
     hfitqual_fit_0lep_3b->Add( hfitqual_ttwj_0lep_3b ) ;
     hfitqual_fit_0lep_3b->Add( hfitqual_susy_0lep_3b ) ;



     hfitqual_fit_1lep_1b->Add( hfitqual_znn_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_qcd_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_ttwj_1lep_1b ) ;
     hfitqual_fit_1lep_1b->Add( hfitqual_susy_1lep_1b ) ;

     hfitqual_fit_1lep_2b->Add( hfitqual_znn_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_qcd_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_ttwj_1lep_2b ) ;
     hfitqual_fit_1lep_2b->Add( hfitqual_susy_1lep_2b ) ;

     hfitqual_fit_1lep_3b->Add( hfitqual_znn_1lep_3b ) ;
     hfitqual_fit_1lep_3b->Add( hfitqual_qcd_1lep_3b ) ;
     hfitqual_fit_1lep_3b->Add( hfitqual_ttwj_1lep_3b ) ;
     hfitqual_fit_1lep_3b->Add( hfitqual_susy_1lep_3b ) ;




     hfitqual_fit_ldp_1b->Add( hfitqual_znn_ldp_1b ) ;
     hfitqual_fit_ldp_1b->Add( hfitqual_qcd_ldp_1b ) ;
     hfitqual_fit_ldp_1b->Add( hfitqual_ttwj_ldp_1b ) ;
     hfitqual_fit_ldp_1b->Add( hfitqual_susy_ldp_1b ) ;

     hfitqual_fit_ldp_2b->Add( hfitqual_znn_ldp_2b ) ;
     hfitqual_fit_ldp_2b->Add( hfitqual_qcd_ldp_2b ) ;
     hfitqual_fit_ldp_2b->Add( hfitqual_ttwj_ldp_2b ) ;
     hfitqual_fit_ldp_2b->Add( hfitqual_susy_ldp_2b ) ;

     hfitqual_fit_ldp_3b->Add( hfitqual_znn_ldp_3b ) ;
     hfitqual_fit_ldp_3b->Add( hfitqual_qcd_ldp_3b ) ;
     hfitqual_fit_ldp_3b->Add( hfitqual_ttwj_ldp_3b ) ;
     hfitqual_fit_ldp_3b->Add( hfitqual_susy_ldp_3b ) ;



     printf("\n\n Done making stacks.\n\n") ; cout << flush ;

     printf("\n Making legend...\n") ; cout << flush ;
     TLegend* legend = new TLegend(0.4,0.35,0.7,0.85) ;

     legend->AddEntry( hfitqual_data_0lep_1b, "data" ) ;
     legend->AddEntry( hfitqual_susy_0lep_1b, "SUSY" ) ;
     legend->AddEntry( hfitqual_ttwj_0lep_1b, "ttwj" ) ;
     legend->AddEntry( hfitqual_qcd_0lep_1b,  "QCD" ) ;
     legend->AddEntry( hfitqual_znn_0lep_1b,  "Znunu" ) ;
     legend->AddEntry( hfitqual_np,           "Eff PG" ) ;

     printf("\n\n Done making legend.\n\n\n") ; cout << flush ;

     TCanvas* cfitqual = (TCanvas*) gDirectory->FindObject("cfitqual") ;
     if ( cfitqual==0 ) {
        printf("\n\n Creating cfitqual canvas.\n\n") ;
        cfitqual = new TCanvas("cfitqual","RA2b fit quality", 850, 1000 ) ;
     } else {
        printf("\n\n Found existing cfitqual canvas.\n\n") ;
        cfitqual->Clear() ;
     }



     if ( metgroupzoom>1 && !logy ) {

        int nhistbins = hfitqual_ttwj_0lep_1b->GetNbinsX() ;

        int nBinsHT(0) ;
        for ( int bi=2; bi<nhistbins; bi++ ) {
           if ( hfitqual_ttwj_0lep_1b->GetBinContent( bi ) <= 0. ) break ;
           nBinsHT++ ;
        }

        int zoomRefBin = 1 + (nBinsHT+1)*(metgroupzoom-1) + 1 ;

        printf("\n\n Number of HT bins=%d, met group ref bin=%d\n\n", nBinsHT, zoomRefBin ) ;

        double maxSF(1.3) ;

        hfitqual_data_0lep_1b->SetMaximum( maxSF*(hfitqual_data_0lep_1b->GetBinContent( zoomRefBin )) ) ;
        hfitqual_data_0lep_2b->SetMaximum( maxSF*(hfitqual_data_0lep_2b->GetBinContent( zoomRefBin )) ) ;
        hfitqual_data_0lep_3b->SetMaximum( maxSF*(hfitqual_data_0lep_3b->GetBinContent( zoomRefBin )) ) ;

        hfitqual_data_1lep_1b->SetMaximum( maxSF*(hfitqual_data_1lep_1b->GetBinContent( zoomRefBin )) ) ;
        hfitqual_data_1lep_2b->SetMaximum( maxSF*(hfitqual_data_1lep_2b->GetBinContent( zoomRefBin )) ) ;
        hfitqual_data_1lep_3b->SetMaximum( maxSF*(hfitqual_data_1lep_3b->GetBinContent( zoomRefBin )) ) ;

        hfitqual_data_ldp_1b->SetMaximum( maxSF*(hfitqual_data_ldp_1b->GetBinContent( zoomRefBin+2 )) ) ;
        hfitqual_data_ldp_2b->SetMaximum( maxSF*(hfitqual_data_ldp_2b->GetBinContent( zoomRefBin+2 )) ) ;
        hfitqual_data_ldp_3b->SetMaximum( maxSF*(hfitqual_data_ldp_3b->GetBinContent( zoomRefBin+2 )) ) ;

        hfitqual_data_zee_1b->SetMaximum( maxSF*(hfitqual_data_zee_1b->GetBinContent( zoomRefBin )) ) ;
        hfitqual_data_zmm_1b->SetMaximum( maxSF*(hfitqual_data_zmm_1b->GetBinContent( zoomRefBin )) ) ;

     }


     cfitqual->Divide(3,4);

     gPad->SetTicks(1,0) ;


     printf(" pad 1\n") ; cout << flush ;
     cfitqual->cd(1);
     hfitqual_data_0lep_1b->Draw("histpe") ;
     hfitqual_fit_0lep_1b->Draw("same") ;
     hfitqual_data_0lep_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     printf(" pad 2\n") ; cout << flush ;
     cfitqual->cd(2);
     hfitqual_data_0lep_2b->Draw("histpe") ;
     hfitqual_fit_0lep_2b->Draw("same") ;
     hfitqual_data_0lep_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     printf(" pad 3\n") ; cout << flush ;
     cfitqual->cd(3);
     hfitqual_data_0lep_3b->Draw("histpe") ;
     hfitqual_fit_0lep_3b->Draw("same") ;
     hfitqual_data_0lep_3b->Draw("same") ;
     gPad->SetGridy(1) ;
     


     printf(" pad 4\n") ; cout << flush ;
     cfitqual->cd(4);
     hfitqual_data_1lep_1b->Draw("histpe") ;
     hfitqual_fit_1lep_1b->Draw("same") ;
     hfitqual_data_1lep_1b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     printf(" pad 5\n") ; cout << flush ;
     cfitqual->cd(5);
     hfitqual_data_1lep_2b->Draw("histpe") ;
     hfitqual_fit_1lep_2b->Draw("same") ;
     hfitqual_data_1lep_2b->Draw("same") ;
     gPad->SetGridy(1) ;
     
     printf(" pad 6\n") ; cout << flush ;
     cfitqual->cd(6);
     hfitqual_data_1lep_3b->Draw("histpe") ;
     hfitqual_fit_1lep_3b->Draw("same") ;
     hfitqual_data_1lep_3b->Draw("same") ;
     gPad->SetGridy(1) ;




     
     printf(" pad 7\n") ; cout << flush ;
     cfitqual->cd(7);
     hfitqual_data_ldp_1b->Draw("histpe") ;
     hfitqual_fit_ldp_1b->Draw("same") ;
     hfitqual_data_ldp_1b->Draw("same") ;
     gPad->SetGridy(1) ;

     printf(" pad 8\n") ; cout << flush ;
     cfitqual->cd(8);
     hfitqual_data_ldp_2b->Draw("histpe") ;
     hfitqual_fit_ldp_2b->Draw("same") ;
     hfitqual_data_ldp_2b->Draw("same") ;
     gPad->SetGridy(1) ;

     printf(" pad 9\n") ; cout << flush ;
     cfitqual->cd(9);
     hfitqual_data_ldp_3b->Draw("histpe") ;
     hfitqual_fit_ldp_3b->Draw("same") ;
     hfitqual_data_ldp_3b->Draw("same") ;
     gPad->SetGridy(1) ;




     cfitqual->cd(10) ;
     hfitqual_data_zee_1b->Draw("histpe") ;
     hfitqual_fit_zee_1b->Draw("same") ;
     hfitqual_data_zee_1b->Draw("same") ;
     gPad->SetGridy(1) ;

     cfitqual->cd(11) ;
     hfitqual_data_zmm_1b->Draw("histpe") ;
     hfitqual_fit_zmm_1b->Draw("same") ;
     hfitqual_data_zmm_1b->Draw("same") ;
     gPad->SetGridy(1) ;




     cfitqual->cd(12);
     legend->Draw() ;

     cfitqual->Update() ;



     //--- Efficiency scale factor primary Gaussian value.

     hfitqual_np->SetMinimum(-5.) ;
     hfitqual_np->SetMaximum( 5.) ;
     
     hfitqual_np->SetNdivisions(101,"x") ;
     hfitqual_np->SetNdivisions(101,"y") ;
     hfitqual_np->SetLabelOffset(99,"y") ;
     

     TPad* tp = new TPad("tp","tp",0.09,0.,0.18,1.0) ;

     tp->SetLogy(0);
     tp->SetRightMargin(0.4) ;

     tp->Draw() ;
     tp->cd() ;
     hfitqual_np->SetLabelSize(0.5,"x") ;
     TAxis *xaxis ;
     xaxis = hfitqual_np->GetXaxis() ;
     xaxis->SetBinLabel(1,"Eff PG") ;
     hfitqual_np->GetXaxis()->LabelsOption("v") ;
     hfitqual_np->Draw() ;

     cfitqual->Update() ;


     TGaxis* axis = new TGaxis() ;
     axis->SetLabelOffset(0.1) ;
     axis->SetLabelSize(0.30) ;
     axis->SetTickSize(0.2) ;
     axis->DrawAxis( 1.0, -5., 1.0, 5., -5., 5., 510, "+LS") ;
     
     cfitqual->Update() ;

     cfitqual->SaveAs("fitqual.gif") ;






     //--- some numerical output below here.

     int nbins = hfitqual_ttwj_0lep_1b -> GetNbinsX() ;




     double ttwj0lep1btotal(0.) ;
     double ttwj0lep2btotal(0.) ;
     double ttwj0lep3btotal(0.) ;
     double ttwj0leptotal(0.) ;

     double ttwj1lep1btotal(0.) ;
     double ttwj1lep2btotal(0.) ;
     double ttwj1lep3btotal(0.) ;
     double ttwj1leptotal(0.) ;

     double susy0leptotal(0.) ;

     for ( int bi=1; bi<=nbins; bi++ ) {
        ttwj0lep1btotal += hfitqual_ttwj_0lep_1b->GetBinContent(bi) ;
        ttwj0lep2btotal += hfitqual_ttwj_0lep_2b->GetBinContent(bi) ;
        ttwj0lep3btotal += hfitqual_ttwj_0lep_3b->GetBinContent(bi) ;
        ttwj0leptotal += hfitqual_ttwj_0lep_1b->GetBinContent(bi) + hfitqual_ttwj_0lep_2b->GetBinContent(bi) + hfitqual_ttwj_0lep_3b->GetBinContent(bi) ;
        ttwj1lep1btotal += hfitqual_ttwj_1lep_1b->GetBinContent(bi) ;
        ttwj1lep2btotal += hfitqual_ttwj_1lep_2b->GetBinContent(bi) ;
        ttwj1lep3btotal += hfitqual_ttwj_1lep_3b->GetBinContent(bi) ;
        ttwj1leptotal += hfitqual_ttwj_1lep_1b->GetBinContent(bi) + hfitqual_ttwj_1lep_2b->GetBinContent(bi) + hfitqual_ttwj_1lep_3b->GetBinContent(bi) ;
        susy0leptotal += hfitqual_susy_0lep_1b->GetBinContent(bi) + hfitqual_susy_0lep_2b->GetBinContent(bi) + hfitqual_susy_0lep_3b->GetBinContent(bi) ;
     }

     printf("\n\n 1b, ttwj 0lep/1lep ratio: %6.1f / %6.1f = %6.3f\n", ttwj0lep1btotal, ttwj1lep1btotal, ttwj0lep1btotal/ttwj1lep1btotal ) ;
     printf(" 2b, ttwj 0lep/1lep ratio: %6.1f / %6.1f = %6.3f\n", ttwj0lep2btotal, ttwj1lep2btotal, ttwj0lep2btotal/ttwj1lep2btotal ) ;
     printf(" 3b, ttwj 0lep/1lep ratio: %6.1f / %6.1f = %6.3f\n", ttwj0lep3btotal, ttwj1lep3btotal, ttwj0lep3btotal/ttwj1lep3btotal ) ;
     printf("     ttwj 0lep/1lep ratio: %6.1f / %6.1f = %6.3f\n\n\n", ttwj0leptotal, ttwj1leptotal, ttwj0leptotal/ttwj1leptotal ) ;

     printf("\n\n SUSY 0lep total : %7.2f\n", susy0leptotal ) ;






     double qcd0lep1btotal(0.) ;
     double qcd0lep2btotal(0.) ;
     double qcd0lep3btotal(0.) ;
     double qcd0leptotal(0.) ;

     double qcdldp1btotal(0.) ;
     double qcdldp2btotal(0.) ;
     double qcdldp3btotal(0.) ;
     double qcdldptotal(0.) ;

     for ( int bi=1; bi<=nbins; bi++ ) {
        qcd0lep1btotal += hfitqual_qcd_0lep_1b->GetBinContent(bi) ;
        qcd0lep2btotal += hfitqual_qcd_0lep_2b->GetBinContent(bi) ;
        qcd0lep3btotal += hfitqual_qcd_0lep_3b->GetBinContent(bi) ;
        qcd0leptotal += hfitqual_qcd_0lep_1b->GetBinContent(bi) + hfitqual_qcd_0lep_2b->GetBinContent(bi) + hfitqual_qcd_0lep_3b->GetBinContent(bi) ;
        qcdldp1btotal += hfitqual_qcd_ldp_1b->GetBinContent(bi) ;
        qcdldp2btotal += hfitqual_qcd_ldp_2b->GetBinContent(bi) ;
        qcdldp3btotal += hfitqual_qcd_ldp_3b->GetBinContent(bi) ;
        qcdldptotal += hfitqual_qcd_ldp_1b->GetBinContent(bi) + hfitqual_qcd_ldp_2b->GetBinContent(bi) + hfitqual_qcd_ldp_3b->GetBinContent(bi) ;
     }

     printf("\n\n 1b, qcd 0lep/ldp ratio: %6.1f / %6.1f = %6.3f\n", qcd0lep1btotal, qcdldp1btotal, qcd0lep1btotal/qcdldp1btotal ) ;
     printf(" 2b, qcd 0lep/ldp ratio: %6.1f / %6.1f = %6.3f\n", qcd0lep2btotal, qcdldp2btotal, qcd0lep2btotal/qcdldp2btotal ) ;
     printf(" 3b, qcd 0lep/ldp ratio: %6.1f / %6.1f = %6.3f\n", qcd0lep3btotal, qcdldp3btotal, qcd0lep3btotal/qcdldp3btotal ) ;
     printf("     qcd 0lep/ldp ratio: %6.1f / %6.1f = %6.3f\n\n\n", qcd0leptotal, qcdldptotal, qcd0leptotal/qcdldptotal ) ;









   //--- nuisance parameters
     TH1F* hnp_qcd_1b_val = (TH1F*) gDirectory->FindObject("hnp_qcd_1b_val") ;
     if ( hnp_qcd_1b_val == 0x0 ) { printf("\n\n *** No nuisance parameter plots.\n\n") ; return ; }
     TH1F* hnp_qcd_2b_val = (TH1F*) gDirectory->FindObject("hnp_qcd_2b_val") ;
     TH1F* hnp_qcd_3b_val = (TH1F*) gDirectory->FindObject("hnp_qcd_3b_val") ;
     TH1F* hnp_qcd_1b_pull = (TH1F*) gDirectory->FindObject("hnp_qcd_1b_pull") ;
     TH1F* hnp_qcd_2b_pull = (TH1F*) gDirectory->FindObject("hnp_qcd_2b_pull") ;
     TH1F* hnp_qcd_3b_pull = (TH1F*) gDirectory->FindObject("hnp_qcd_3b_pull") ;

     TH1F* hnp_ttwj_1b_val = (TH1F*) gDirectory->FindObject("hnp_ttwj_1b_val") ;
     TH1F* hnp_ttwj_2b_val = (TH1F*) gDirectory->FindObject("hnp_ttwj_2b_val") ;
     TH1F* hnp_ttwj_3b_val = (TH1F*) gDirectory->FindObject("hnp_ttwj_3b_val") ;
     TH1F* hnp_ttwj_1b_pull = (TH1F*) gDirectory->FindObject("hnp_ttwj_1b_pull") ;
     TH1F* hnp_ttwj_2b_pull = (TH1F*) gDirectory->FindObject("hnp_ttwj_2b_pull") ;
     TH1F* hnp_ttwj_3b_pull = (TH1F*) gDirectory->FindObject("hnp_ttwj_3b_pull") ;

     gStyle->SetPadGridY(1) ;
     TCanvas* cnp = (TCanvas*) gDirectory->FindObject("cnp") ;
     if ( cnp == 0x0 ) {
        cnp = new TCanvas("cnp","RA2b nuisance pars", 850, 1000 ) ;
     }
     cnp->Divide(3,4) ;

   //---
     cnp->cd(1) ;
     hnp_ttwj_1b_val->Draw() ;

     cnp->cd(2) ;
     hnp_ttwj_2b_val->Draw() ;

     cnp->cd(3) ;
     hnp_ttwj_3b_val->Draw() ;


   //---
     cnp->cd(4) ;
     hnp_ttwj_1b_pull->Draw() ;

     cnp->cd(5) ;
     hnp_ttwj_2b_pull->Draw() ;

     cnp->cd(6) ;
     hnp_ttwj_3b_pull->Draw() ;


   //---
     cnp->cd(7) ;
     hnp_qcd_1b_val->Draw() ;

     cnp->cd(8) ;
     hnp_qcd_2b_val->Draw() ;

     cnp->cd(9) ;
     hnp_qcd_3b_val->Draw() ;


   //---
     cnp->cd(10) ;
     hnp_qcd_1b_pull->Draw() ;

     cnp->cd(11) ;
     hnp_qcd_2b_pull->Draw() ;

     cnp->cd(12) ;
     hnp_qcd_3b_pull->Draw() ;







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

