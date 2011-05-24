
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TH1.h"

#include <stdio.h>


   void plottoysb( const char* infile = "toy-fits-2011-ge1-sbvars-model.root" ) {

      TChain tf("toyfittree") ;
      int nt = tf.Add( infile ) ;
      if ( nt < 1 ) {
         printf("\n\n *** Didn't find toyfittree TTree in %s\n\n", infile ) ;
         return ;
      }

      double toy_mu0_ttbar_sb ;
      double toy_mu0_qcd_sb ;

      TBranch *b_toy_mu0_ttbar_sb ;
      TBranch *b_toy_mu0_qcd_sb ;

      tf.SetBranchAddress("mu0_ttbar_sb" , &toy_mu0_ttbar_sb , &b_toy_mu0_ttbar_sb);
      tf.SetBranchAddress("mu0_qcd_sb"   , &toy_mu0_qcd_sb   , &b_toy_mu0_qcd_sb);

      tf.GetEntry(0) ;
      printf("True ttbar, sb %5.1f\n", toy_mu0_ttbar_sb ) ;
      printf("True qcd,   sb %5.1f\n", toy_mu0_qcd_sb ) ;

      TCanvas* ctoy_tt = new TCanvas("ctoy_sb","Toy study, sb region", 1150, 650) ;
      ctoy_tt->Clear() ;
      ctoy_tt->Divide(3,2) ;

      gStyle->SetOptStat("erm") ;
      gStyle->SetStatW(0.35) ;
      gStyle->SetStatH(0.25) ;

      int nbins(60) ;

      TH1F* hsb_ttbar_val = new TH1F("hsb_ttbar_val","Nttbar,sb value", nbins, 0., 260. ) ;
      TH1F* hsb_ttbar_err = new TH1F("hsb_ttbar_err","Nttbar,sb error", nbins, 0., 80. ) ;
      TH1F* hsb_ttbar_pull = new TH1F("hsb_ttbar_pull","Nttbar,sb pull", nbins, -8., 8. ) ;

      TH1F* hsb_qcd_val = new TH1F("hsb_qcd_val","Nqcd,sb value", nbins, 0., 260. ) ;
      TH1F* hsb_qcd_err = new TH1F("hsb_qcd_err","Nqcd,sb error", nbins, 0., 80. ) ;
      TH1F* hsb_qcd_pull = new TH1F("hsb_qcd_pull","Nqcd,sb pull", nbins, -8., 8. ) ;

      hsb_ttbar_val->SetFillColor(11) ;
      hsb_ttbar_err->SetFillColor(11) ;
      hsb_ttbar_pull->SetFillColor(11) ;

      hsb_qcd_val->SetFillColor(11) ;
      hsb_qcd_err->SetFillColor(11) ;
      hsb_qcd_pull->SetFillColor(11) ;

      TArrow* arrow = new TArrow() ;
      arrow->SetLineWidth(2) ;
      arrow->SetLineColor(4) ;


   //-- ttbar
      ctoy_tt->cd(1) ;
      tf.Draw( "mu_ttbar_sb>>hsb_ttbar_val", "" ) ;
      arrow->DrawArrow( toy_mu0_ttbar_sb, 0.4*(hsb_ttbar_val->GetMaximum()), toy_mu0_ttbar_sb, 0., 0.02, ">" ) ;
      ctoy_tt->Update() ;

      ctoy_tt->cd(2) ;
      tf.Draw( "mu_ttbar_sb_err>>hsb_ttbar_err", "" ) ;
      ctoy_tt->Update() ;

      ctoy_tt->cd(3) ;
      tf.Draw( "(mu_ttbar_sb-mu0_ttbar_sb)/mu_ttbar_sb_err>>hsb_ttbar_pull", "" ) ;
      ctoy_tt->Update() ;

   //-- qcd
      ctoy_tt->cd(4) ;
      tf.Draw( "mu_qcd_sb>>hsb_qcd_val", "" ) ;
      arrow->DrawArrow( toy_mu0_qcd_sb, 0.4*(hsb_qcd_val->GetMaximum()), toy_mu0_qcd_sb, 0., 0.02, ">" ) ;
      ctoy_tt->Update() ;

      ctoy_tt->cd(5) ;
      tf.Draw( "mu_qcd_sb_err>>hsb_qcd_err", "" ) ;
      ctoy_tt->Update() ;

      ctoy_tt->cd(6) ;
      tf.Draw( "(mu_qcd_sb-mu0_qcd_sb)/mu_qcd_sb_err>>hsb_qcd_pull", "" ) ;
      ctoy_tt->Update() ;


   }

