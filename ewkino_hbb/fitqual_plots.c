
#include "TMath.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2.h"
#include "TAxis.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TText.h"


#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooConstVar.h"
#include "RooStats/ModelConfig.h"



  using namespace RooFit;
  using namespace RooStats;



   void fitqual_plots( const char* wsfile = "outputfiles/ws.root" ) {

      gStyle -> SetOptStat(0) ;

      TFile* wstf = new TFile( wsfile ) ;

      RooWorkspace* ws = dynamic_cast<RooWorkspace*>( wstf->Get("ws") );
      ws->Print() ;

      int bins_of_met = TMath::Nint( ws->var("bins_of_met")->getVal()  ) ;
      printf("\n\n Bins of MET : %d\n\n", bins_of_met ) ;

      int bins_of_nb = TMath::Nint( ws->var("bins_of_nb")->getVal()  ) ;
      printf("\n\n Bins of nb : %d\n\n", bins_of_nb ) ;

      ModelConfig* modelConfig = (ModelConfig*) ws->obj( "SbModel" ) ;

      RooDataSet* rds = (RooDataSet*) ws->obj( "hbb_observed_rds" ) ;

      rds->Print() ;
      rds->printMultiline(cout, 1, kTRUE, "") ;

      RooAbsPdf* likelihood = modelConfig->GetPdf() ;

      RooFitResult* fitResult = likelihood->fitTo( *rds, Save(true), PrintLevel(0) ) ;
      fitResult->Print() ;


      char hname[1000] ;
      char htitle[1000] ;
      char pname[1000] ;




     //-- unpack observables.

      int obs_N_msig[10][50] ; // first index is n btags, second is met bin.
      int obs_N_msb[10][50]  ; // first index is n btags, second is met bin.

      const RooArgSet* dsras = rds->get() ;
      TIterator* obsIter = dsras->createIterator() ;
      while ( RooRealVar* obs = (RooRealVar*) obsIter->Next() ) {
         for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
            for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
               sprintf( pname, "N_%db_msig_met%d", nbi+2, mbi+1 ) ;
               if ( strcmp( obs->GetName(), pname ) == 0 ) { obs_N_msig[nbi][mbi] = TMath::Nint( obs -> getVal() ) ; }
               sprintf( pname, "N_%db_msb_met%d", nbi+2, mbi+1 ) ;
               if ( strcmp( obs->GetName(), pname ) == 0 ) { obs_N_msb[nbi][mbi] = TMath::Nint( obs -> getVal() ) ; }
            } // mbi.
         } // nbi.
      } // obs iterator.


      printf("\n\n") ;
      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {
         printf(" nb=%d :  ", nbi+2 ) ;
         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
            printf("  sig=%3d, sb=%3d  |", obs_N_msig[nbi][mbi], obs_N_msb[nbi][mbi] ) ;
         } // mbi.
         printf("\n") ;
      } // nbi.
      printf("\n\n") ;




      TCanvas* can1 = (TCanvas*) gDirectory->FindObject("can1") ;
      if ( can1 == 0x0 ) {
         can1 = new TCanvas("can1","hbb fit", 700, 1000 ) ;
      }
      can1->Clear() ;
      can1->Divide( 2, bins_of_nb+1 ) ;

      int pad(1) ;

      for ( int nbi=0; nbi<bins_of_nb; nbi++ ) {


         sprintf( hname, "h_bg_%db_msig_met", nbi+2 ) ;
         sprintf( htitle, "mass sig, %db, MET", nbi+2 ) ;
         TH1F* hist_bg_msig = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;
         hist_bg_msig -> SetFillColor( kBlue-9 ) ;

         sprintf( hname, "h_bg_%db_msb_met", nbi+2 ) ;
         sprintf( htitle, "mass sig, %db, MET", nbi+2 ) ;
         TH1F* hist_bg_msb = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;
         hist_bg_msb -> SetFillColor( kBlue-9 ) ;

         sprintf( hname, "h_sig_%db_msig_met", nbi+2 ) ;
         sprintf( htitle, "mass sig, %db, MET", nbi+2 ) ;
         TH1F* hist_sig_msig = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;
         hist_sig_msig -> SetFillColor( kMagenta+2 ) ;

         sprintf( hname, "h_sig_%db_msb_met", nbi+2 ) ;
         sprintf( htitle, "mass sig, %db, MET", nbi+2 ) ;
         TH1F* hist_sig_msb = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;
         hist_sig_msb -> SetFillColor( kMagenta+2 ) ;

         sprintf( hname, "h_all_%db_msig_met", nbi+2 ) ;
         sprintf( htitle, "mass sig, %db, MET", nbi+2 ) ;
         TH1F* hist_all_msig = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;

         sprintf( hname, "h_all_%db_msb_met", nbi+2 ) ;
         sprintf( htitle, "mass sig, %db, MET", nbi+2 ) ;
         TH1F* hist_all_msb = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;

         sprintf( hname, "h_data_%db_msig_met", nbi+2 ) ;
         sprintf( htitle, "mass sig, %db, MET", nbi+2 ) ;
         TH1F* hist_data_msig = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;
         hist_data_msig -> SetLineWidth(2) ;
         hist_data_msig -> SetMarkerStyle(20) ;

         sprintf( hname, "h_data_%db_msb_met", nbi+2 ) ;
         sprintf( htitle, "mass sig, %db, MET", nbi+2 ) ;
         TH1F* hist_data_msb = new TH1F( hname, htitle, bins_of_met, 0.5, bins_of_met+0.5 ) ;
         hist_data_msb -> SetLineWidth(2) ;
         hist_data_msb -> SetMarkerStyle(20) ;

         for ( int mbi=0; mbi<bins_of_met; mbi++ ) {



            sprintf( pname, "mu_bg_%db_msig_met%d", nbi+2, mbi+1 ) ;
            RooAbsReal* mu_bg_msig = ws->function( pname ) ;
            if ( mu_bg_msig == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
            hist_bg_msig -> SetBinContent( mbi+1, mu_bg_msig->getVal() ) ;

            sprintf( pname, "mu_sig_%db_msig_met%d", nbi+2, mbi+1 ) ;
            RooAbsReal* mu_sig_msig = ws->function( pname ) ;
            if ( mu_sig_msig == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
            hist_sig_msig -> SetBinContent( mbi+1, mu_sig_msig->getVal() ) ;

            hist_all_msig -> SetBinContent( mbi+1, mu_bg_msig->getVal() + mu_sig_msig->getVal() ) ;

            hist_data_msig -> SetBinContent( mbi+1, obs_N_msig[nbi][mbi] ) ;



            sprintf( pname, "mu_bg_%db_msb_met%d", nbi+2, mbi+1 ) ;
            RooAbsReal* mu_bg_msb = ws->function( pname ) ;
            if ( mu_bg_msb == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
            hist_bg_msb -> SetBinContent( mbi+1, mu_bg_msb->getVal() ) ;

            sprintf( pname, "mu_sig_%db_msb_met%d", nbi+2, mbi+1 ) ;
            RooAbsReal* mu_sig_msb = ws->function( pname ) ;
            if ( mu_sig_msb == 0x0 ) { printf("\n\n *** ws missing %s\n\n", pname ) ; return ; }
            hist_sig_msb -> SetBinContent( mbi+1, mu_sig_msb->getVal() ) ;

            hist_all_msb -> SetBinContent( mbi+1, mu_bg_msb->getVal() + mu_sig_msb->getVal() ) ;

            hist_data_msb -> SetBinContent( mbi+1, obs_N_msb[nbi][mbi] ) ;



         } // mbi.

         can1->cd( pad ) ;

         sprintf( hname, "h_stack_%db_msig_met", nbi+2 ) ;
         sprintf( htitle, "mass sig, %db, MET", nbi+2 ) ;
         THStack* hstack_msig = new THStack( hname, htitle ) ;
         hstack_msig -> Add( hist_bg_msig ) ;
         hstack_msig -> Add( hist_sig_msig ) ;

         hist_data_msig -> Draw("e") ;
         hstack_msig -> Draw("same") ;
         hist_data_msig -> Draw("same e") ;

         pad++ ;



         can1->cd( pad ) ;

         sprintf( hname, "h_stack_%db_msb_met", nbi+2 ) ;
         sprintf( htitle, "mass sig, %db, MET", nbi+2 ) ;
         THStack* hstack_msb = new THStack( hname, htitle ) ;
         hstack_msb -> Add( hist_bg_msb ) ;
         hstack_msb -> Add( hist_sig_msb ) ;

         hist_data_msb -> Draw("e") ;
         hstack_msb -> Draw("same") ;
         hist_data_msb -> Draw("same e") ;

         pad++ ;



      } // nbi.




      TH1F* hist_R_msigmsb = new TH1F( "h_R_msigmsb", "R msig/msb vs met bin", bins_of_met, 0.5, 0.5+bins_of_met ) ;
      hist_R_msigmsb -> SetLineWidth(2) ;
      hist_R_msigmsb -> SetMarkerStyle(20) ;
      hist_R_msigmsb -> SetXTitle("MET bin") ;
      hist_R_msigmsb -> SetYTitle("R msig/msb") ;


      for ( int mbi=0; mbi<bins_of_met; mbi++ ) {
         sprintf( pname, "R_msigmsb_met%d", mbi+1 ) ;
         RooRealVar* rrv_R = ws->var( pname ) ;
         if ( rrv_R == 0x0 ) { printf("\n\n *** Can't find %s in ws.\n\n", pname ) ; return ; }
         hist_R_msigmsb -> SetBinContent( mbi+1, rrv_R -> getVal() ) ;
         hist_R_msigmsb -> SetBinError( mbi+1, rrv_R -> getError() ) ;
      } // mbi.

      can1->cd( pad ) ;

      hist_R_msigmsb -> Draw("e") ;







   } // fitqual_plots


















