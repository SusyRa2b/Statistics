
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TH1.h"

#include <stdio.h>


   void plottoysig( const char* infile = "toy-fits-2011-ge1-sigvars-mc.root", double xmax=60. ) {

      TChain tf("toyfittree") ;
      int nt = tf.Add( infile ) ;
      if ( nt < 1 ) {
         printf("\n\n *** Didn't find toyfittree TTree in %s\n\n", infile ) ;
         return ;
      }

      double toy_mu0_ttbar_sig ;
      double toy_mu0_qcd_sig ;
      double toy_mu0_ttbar_sb ;
      double toy_mu0_qcd_sb ;
      double toy_mu0_susy_sig ;
      double toy_mu0_allbg_sig ;

      TBranch *b_toy_mu0_ttbar_sig ;
      TBranch *b_toy_mu0_qcd_sig ;
      TBranch *b_toy_mu0_ttbar_sb ;
      TBranch *b_toy_mu0_qcd_sb ;
      TBranch *b_toy_mu0_susy_sig ;
      TBranch *b_toy_mu0_allbg_sig ;

      tf.SetBranchAddress("mu0_ttbar_sig" , &toy_mu0_ttbar_sig , &b_toy_mu0_ttbar_sig);
      tf.SetBranchAddress("mu0_qcd_sig"   , &toy_mu0_qcd_sig   , &b_toy_mu0_qcd_sig);
      tf.SetBranchAddress("mu0_ttbar_sb"  , &toy_mu0_ttbar_sb  , &b_toy_mu0_ttbar_sb);
      tf.SetBranchAddress("mu0_qcd_sb"    , &toy_mu0_qcd_sb    , &b_toy_mu0_qcd_sb);
      tf.SetBranchAddress("mu0_susy_sig"  , &toy_mu0_susy_sig  , &b_toy_mu0_susy_sig);
      tf.SetBranchAddress("mu0_allbg_sig" , &toy_mu0_allbg_sig , &b_toy_mu0_allbg_sig);

      tf.GetEntry(0) ;
      printf("True ttbar, sig %5.1f\n", toy_mu0_ttbar_sig ) ;
      printf("True qcd,   sig %5.1f\n", toy_mu0_qcd_sig ) ;
      printf("True allbg, sig %5.1f\n", toy_mu0_allbg_sig ) ;

      TCanvas* ctoy_tt = new TCanvas("ctoy_sig","Toy study, SIG region", 1150, 950) ;
      ctoy_tt->Clear() ;
      ctoy_tt->Divide(3,3) ;

      gStyle->SetOptStat("erm") ;
      gStyle->SetStatW(0.35) ;
      gStyle->SetStatH(0.25) ;

      int nbins(60) ;

      TH1F* hsig_ttbar_val = new TH1F("hsig_ttbar_val","Nttbar,SIG value", nbins, 0., xmax ) ;
      TH1F* hsig_ttbar_err = new TH1F("hsig_ttbar_err","Nttbar,SIG error", nbins, 0., 50. ) ;
      TH1F* hsig_ttbar_pull = new TH1F("hsig_ttbar_pull","Nttbar,SIG pull", nbins, -8., 8. ) ;

      TH1F* hsig_qcd_val = new TH1F("hsig_qcd_val","Nqcd,SIG value", nbins, 0., xmax ) ;
      TH1F* hsig_qcd_err = new TH1F("hsig_qcd_err","Nqcd,SIG error", nbins, 0., 50. ) ;
      TH1F* hsig_qcd_pull = new TH1F("hsig_qcd_pull","Nqcd,SIG pull", nbins, -8., 8. ) ;

      TH1F* hsig_allbg_val = new TH1F("hsig_allbg_val","Nallbg,SIG value", nbins, 0., xmax ) ;
      TH1F* hsig_allbg_err = new TH1F("hsig_allbg_err","Nallbg,SIG error", nbins, 0., 50. ) ;
      TH1F* hsig_allbg_pull = new TH1F("hsig_allbg_pull","Nallbg,SIG pull", nbins, -8., 8. ) ;

      hsig_ttbar_val->SetFillColor(11) ;
      hsig_ttbar_err->SetFillColor(11) ;
      hsig_ttbar_pull->SetFillColor(11) ;

      hsig_qcd_val->SetFillColor(11) ;
      hsig_qcd_err->SetFillColor(11) ;
      hsig_qcd_pull->SetFillColor(11) ;

      hsig_allbg_val->SetFillColor(11) ;
      hsig_allbg_err->SetFillColor(11) ;
      hsig_allbg_pull->SetFillColor(11) ;


      TArrow* arrow = new TArrow() ;
      arrow->SetLineWidth(2) ;
      arrow->SetLineColor(4) ;


   //-- ttbar
      ctoy_tt->cd(1) ;
      tf.Draw( "mu_ttbar_sig>>hsig_ttbar_val", "" ) ;
      arrow->DrawArrow( toy_mu0_ttbar_sig, 0.4*(hsig_ttbar_val->GetMaximum()), toy_mu0_ttbar_sig, 0., 0.02, ">" ) ;
      ctoy_tt->Update() ;

      ctoy_tt->cd(2) ;
      tf.Draw( "mu_ttbar_sig_err>>hsig_ttbar_err", "" ) ;
      ctoy_tt->Update() ;

      ctoy_tt->cd(3) ;
      tf.Draw( "(mu_ttbar_sig-mu0_ttbar_sig)/mu_ttbar_sig_err>>hsig_ttbar_pull", "" ) ;
      ctoy_tt->Update() ;

   //-- qcd
      ctoy_tt->cd(4) ;
      tf.Draw( "mu_qcd_sig>>hsig_qcd_val", "" ) ;
      arrow->DrawArrow( toy_mu0_qcd_sig, 0.4*(hsig_qcd_val->GetMaximum()), toy_mu0_qcd_sig, 0., 0.02, ">" ) ;
      ctoy_tt->Update() ;

      ctoy_tt->cd(5) ;
      tf.Draw( "mu_qcd_sig_err>>hsig_qcd_err", "" ) ;
      ctoy_tt->Update() ;

      ctoy_tt->cd(6) ;
      tf.Draw( "(mu_qcd_sig-mu0_qcd_sig)/mu_qcd_sig_err>>hsig_qcd_pull", "" ) ;
      ctoy_tt->Update() ;


   //-- allbg
      ctoy_tt->cd(7) ;
      tf.Draw( "mu_allbg_sig>>hsig_allbg_val", "" ) ;
      arrow->DrawArrow( toy_mu0_allbg_sig, 0.4*(hsig_allbg_val->GetMaximum()), toy_mu0_allbg_sig, 0., 0.02, ">" ) ;
      ctoy_tt->Update() ;

      ctoy_tt->cd(8) ;
      tf.Draw( "mu_allbg_sig_err>>hsig_allbg_err", "" ) ;
      ctoy_tt->Update() ;

      ctoy_tt->cd(9) ;
      tf.Draw( "(mu_allbg_sig-mu0_allbg_sig)/mu_allbg_sig_err>>hsig_allbg_pull", "" ) ;
      ctoy_tt->Update() ;


   }

