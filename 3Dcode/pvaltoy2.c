
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"


   void pvaltoy2( double Nobs_sl = 0, double Nobs_ldp = 1, double Nobs_zl = 4 ) {

      gStyle->SetOptStat(0) ;

      double mu_znn_zl = 0.20 ;
      double mu_znn_ldp = 0.06 ;
      double Rttwj_zloversl = 0.93 ;
      double Rttwj_ldpoverzl = 0.66 ;
      double Sttwj = 1.59 ;
      double Sqcd = 1.00 ;
      double K3qcd = 0.14 ;
      double etrig_had = 1.00 ;
      double etrig_sl = 1.00 ;




      double min_mu_qcd_ldp(0.) ;
      double max_mu_qcd_ldp(15.0) ;
      //double max_mu_qcd_ldp(6.0) ;

      double min_mu_ttwj_sl(0.) ;
      double max_mu_ttwj_sl(6.0) ;
      //double max_mu_ttwj_sl(3.0) ;

      const int nbins_mu_qcd_ldp(120) ;
      const int nbins_mu_ttwj_sl(60) ;
      double mu_qcd_ldp_val[nbins_mu_qcd_ldp] ;
      double mu_ttwj_sl_val[nbins_mu_ttwj_sl] ;
      double lh_2bin_val[nbins_mu_qcd_ldp][nbins_mu_ttwj_sl] ;
      double lh_2bin_norm[nbins_mu_qcd_ldp][nbins_mu_ttwj_sl] ;
      double lh_2bin_cdf[nbins_mu_qcd_ldp][nbins_mu_ttwj_sl] ;

      TH2F* hlh_2bin_data_obs = new TH2F("hlh_2bin_data_obs","SL*LDP 2-bin likelihood, data Nobs, mu_ttwj_sl vs mu_qcd_ldp",
                nbins_mu_qcd_ldp, min_mu_qcd_ldp, max_mu_qcd_ldp,
                nbins_mu_ttwj_sl, min_mu_ttwj_sl, max_mu_ttwj_sl ) ;

      TH2F* hlh_3bin_data_obs = new TH2F("hlh_3bin_data_obs","SL*LDP 3-bin likelihood, data Nobs, mu_ttwj_sl vs mu_qcd_ldp",
                nbins_mu_qcd_ldp, min_mu_qcd_ldp, max_mu_qcd_ldp,
                nbins_mu_ttwj_sl, min_mu_ttwj_sl, max_mu_ttwj_sl ) ;

      TH2F* hlh_2bin_gen_data_obs = new TH2F("hlh_2bin_gen_data_obs","2-bin likelihood Generated mu_ttwj_sl vs mu_qcd_ldp",
                nbins_mu_qcd_ldp, min_mu_qcd_ldp, max_mu_qcd_ldp,
                nbins_mu_ttwj_sl, min_mu_ttwj_sl, max_mu_ttwj_sl ) ;

      TH1F* h_Nzl_gen = new TH1F("h_Nzl_gen","N_zl, generated", 16, -0.5, 15.5 ) ;
      TH1F* h_Nsl_gen = new TH1F("h_Nsl_gen","N_sl, generated", 16, -0.5, 15.5 ) ;
      TH1F* h_Nldp_gen = new TH1F("h_Nldp_gen","N_ldp, generated", 16, -0.5, 15.5 ) ;

      TH1F* h_toy_3bin_max_lh = new TH1F("h_toy_3bin_max_lh", "Toy Min 3-bin -lnL", 60, 0., 10. ) ;

     //--- Compute the 2-bin and 3-bin likelihood for the observed data.

      double data_2bin_lh_max(0.) ;
      double data_2bin_lh_max_mu_qcd_ldp(0.) ;
      double data_2bin_lh_max_mu_ttwj_sl(0.) ;
      double data_2bin_lh_max_n_zl(0.) ;

      double data_3bin_lh_max(0.) ;
      double data_3bin_lh_max_mu_qcd_ldp(0.) ;
      double data_3bin_lh_max_mu_ttwj_sl(0.) ;
      double data_3bin_lh_max_n_zl(0.) ;

      double lh_2bin_sum(0.) ;
      for ( int qbi=0; qbi<nbins_mu_qcd_ldp; qbi++ ) {
         double mu_qcd_ldp = min_mu_qcd_ldp + (qbi+0.5)*(max_mu_qcd_ldp-min_mu_qcd_ldp)/nbins_mu_qcd_ldp ;
         mu_qcd_ldp_val[qbi] = mu_qcd_ldp ;
         for ( int tbi=0; tbi<nbins_mu_ttwj_sl; tbi++ ) {
            double mu_ttwj_sl = min_mu_ttwj_sl + (tbi+0.5)*(max_mu_ttwj_sl-min_mu_ttwj_sl)/nbins_mu_ttwj_sl ;
            mu_ttwj_sl_val[tbi] = mu_ttwj_sl ;
            double mu_ttwj_zl = Sttwj * Rttwj_zloversl * mu_ttwj_sl ;
            double mu_ttwj_ldp = Rttwj_ldpoverzl * mu_ttwj_zl ;
            double mu_qcd_zl = Sqcd * K3qcd * mu_qcd_ldp ;
            double n_ldp = etrig_had * mu_qcd_ldp + etrig_sl * ( mu_ttwj_ldp + mu_znn_ldp ) ;
            double n_sl = etrig_sl * mu_ttwj_sl ;
            double n_zl = etrig_sl * mu_ttwj_zl + etrig_had * mu_qcd_zl + etrig_sl * mu_znn_zl ;
            double lh_2bin = TMath::Poisson( Nobs_sl, n_sl) * TMath::Poisson( Nobs_ldp, n_ldp ) ;
            lh_2bin_val[qbi][tbi] = lh_2bin ;
            lh_2bin_sum += lh_2bin ;
            //// printf(" qbi,tbi : %2d,%2d : n_sl = %6.2f, n_ldp = %6.2f, lh = %g\n", qbi, tbi, n_sl, n_ldp, lh_2bin ) ;
            hlh_2bin_data_obs -> SetBinContent( qbi+1, tbi+1, lh_2bin ) ;
            if ( lh_2bin > data_2bin_lh_max ) {
               data_2bin_lh_max = lh_2bin ;
               data_2bin_lh_max_mu_qcd_ldp = mu_qcd_ldp ;
               data_2bin_lh_max_mu_ttwj_sl = mu_ttwj_sl ;
               data_2bin_lh_max_n_zl = n_zl ;
            }
            double lh_3bin = TMath::Poisson( Nobs_zl, n_zl) * TMath::Poisson( Nobs_sl, n_sl) * TMath::Poisson( Nobs_ldp, n_ldp ) ;
            hlh_3bin_data_obs -> SetBinContent( qbi+1, tbi+1, lh_3bin ) ;
            if ( lh_3bin > data_3bin_lh_max ) {
               data_3bin_lh_max = lh_3bin ;
               data_3bin_lh_max_mu_qcd_ldp = mu_qcd_ldp ;
               data_3bin_lh_max_mu_ttwj_sl = mu_ttwj_sl ;
               data_3bin_lh_max_n_zl = n_zl ;
            }
         } // tbi
      } // qbi

      printf("\n\n Max 2-bin likelihood for Nsl_data = %.0f, Nldp_data = %.0f: mu_ttwj_sl = %6.2f, mu_qcd_ldp = %6.2f, n_zl = %6.2f, -lnL = %g\n\n",
      Nobs_sl, Nobs_ldp, data_2bin_lh_max_mu_ttwj_sl, data_2bin_lh_max_mu_qcd_ldp, data_2bin_lh_max_n_zl, -log(data_2bin_lh_max) ) ;

      printf("\n\n Max 3-bin likelihood for Nsl_data = %.0f, Nldp_data = %.0f, Nzl_data = %.0f: mu_ttwj_sl = %6.2f, mu_qcd_ldp = %6.2f, n_zl = %6.2f, -lnL = %g\n\n",
      Nobs_sl, Nobs_ldp, Nobs_zl, data_3bin_lh_max_mu_ttwj_sl, data_3bin_lh_max_mu_qcd_ldp, data_3bin_lh_max_n_zl, -log(data_3bin_lh_max) ) ;

      double normsum(0.) ;
      for ( int qbi=0; qbi<nbins_mu_qcd_ldp; qbi++ ) {
         for ( int tbi=0; tbi<nbins_mu_ttwj_sl; tbi++ ) {
            lh_2bin_norm[qbi][tbi] = lh_2bin_val[qbi][tbi] / lh_2bin_sum ;
            normsum += lh_2bin_norm[qbi][tbi] ;
            lh_2bin_cdf[qbi][tbi] = normsum ;
         } // tbi
      } // qbi

      TCanvas* can1 = new TCanvas("can1","can1", 1000, 500 ) ;
      can1->Divide(2,1) ;

      can1->cd(1) ;
      hlh_2bin_data_obs -> Draw("colz") ;
      can1->Draw() ;

      can1->cd(2) ;
      hlh_3bin_data_obs -> Draw("colz") ;
      can1->Draw() ;

      can1->Update() ;



      TFile outfile("rootfiles/pvaltoy2.root","recreate") ;

      TTree tt("pvaltoytt", "p-val toy study") ;

      //tt.Branch("


      TRandom tran(12345) ;

      int ngen(100000) ;

      double gen_mu_qcd_ldp ;
      double gen_mu_ttwj_sl ;
      int    Nzl_gen ;
      int    Nsl_gen ;
      int    Nldp_gen ;
      double toy_3bin_lh_max(0.) ;
      double toy_3bin_nll_min(0.) ;
      double toy_3bin_lh_max_mu_qcd_ldp(0.) ;
      double toy_3bin_lh_max_mu_ttwj_sl(0.) ;
      double toy_3bin_lh_max_n_zl(0.) ;

      tt.Branch( "gen_mu_qcd_ldp", &gen_mu_qcd_ldp, "gen_mu_qcd_ldp/D" ) ;
      tt.Branch( "gen_mu_ttwj_sl", &gen_mu_ttwj_sl, "gen_mu_ttwj_sl/D" ) ;
      tt.Branch( "Nzl_gen"       , &Nzl_gen       , "Nzl_gen/I"        ) ;
      tt.Branch( "Nsl_gen"       , &Nsl_gen       , "Nsl_gen/I"        ) ;
      tt.Branch( "Nldp_gen"      , &Nldp_gen      , "Nldp_gen/I"       ) ;
      tt.Branch( "toy_3bin_lh_max"           , &toy_3bin_lh_max           , "toy_3bin_lh_max/D"            ) ;
      tt.Branch( "toy_3bin_nll_min"          , &toy_3bin_nll_min          , "toy_3bin_nll_min/D"           ) ;
      tt.Branch( "toy_3bin_lh_max_mu_qcd_ldp", &toy_3bin_lh_max_mu_qcd_ldp, "toy_3bin_lh_max_mu_qcd_ldp/D" ) ;
      tt.Branch( "toy_3bin_lh_max_mu_ttwj_sl", &toy_3bin_lh_max_mu_ttwj_sl, "toy_3bin_lh_max_mu_ttwj_sl/D" ) ;
      tt.Branch( "toy_3bin_lh_max_n_zl"      , &toy_3bin_lh_max_n_zl      , "toy_3bin_lh_max_n_zl/D"       ) ;

      for ( int ti=0; ti<ngen; ti++ ) {

         if ( ti%(ngen/10) == 0 ) { printf(" %7d, (%4.0f%%)\n", ti, 100*(1.*ti)/(1.*ngen) ) ; }

         double rn = tran.Uniform() ;
         bool done(false) ;
         double mu_qcd_ldp(0.) ;
         double mu_ttwj_sl(0.) ;
         double prev_sum(0.) ;

         for ( int qbi=0; qbi<nbins_mu_qcd_ldp; qbi++ ) {
            for ( int tbi=0; tbi<nbins_mu_ttwj_sl; tbi++ ) {
               if ( rn >= prev_sum && rn < lh_2bin_cdf[qbi][tbi] ) {
                  mu_qcd_ldp = mu_qcd_ldp_val[qbi] ;
                  mu_ttwj_sl = mu_ttwj_sl_val[tbi] ;
                  gen_mu_qcd_ldp = mu_qcd_ldp ;
                  gen_mu_ttwj_sl = mu_ttwj_sl ;
                  done = true ;
                  break ;
               }
               prev_sum = lh_2bin_cdf[qbi][tbi] ;
            } // tbi
            if ( done ) break ;
         } // qbi

      // printf(" ti=%6d : rn = %7.5f, mu_qcd_ldp = %6.2f, mu_ttwj_sl = %6.2f\n", ti, rn, mu_qcd_ldp, mu_ttwj_sl ) ;
         if ( done ) { hlh_2bin_gen_data_obs->Fill( mu_qcd_ldp, mu_ttwj_sl ) ; }

         double mu_qcd_zl = Sqcd * K3qcd * mu_qcd_ldp ;
         double mu_ttwj_zl = Sttwj * Rttwj_zloversl * mu_ttwj_sl ;
         double mu_ttwj_ldp = Rttwj_ldpoverzl * mu_ttwj_zl ;

         double n_zl = etrig_sl * mu_ttwj_zl + etrig_had * mu_qcd_zl + etrig_sl * mu_znn_zl ;
         double n_sl = etrig_sl * mu_ttwj_sl ;
         double n_ldp = etrig_had * mu_qcd_ldp + etrig_sl * ( mu_ttwj_ldp + mu_znn_ldp ) ;
         Nzl_gen = tran.Poisson( n_zl ) ;
         Nsl_gen = tran.Poisson( n_sl ) ;
         Nldp_gen = tran.Poisson( n_ldp ) ;

         h_Nsl_gen  -> Fill( Nsl_gen ) ;
         h_Nldp_gen -> Fill( Nldp_gen ) ;
         h_Nzl_gen  -> Fill( Nzl_gen ) ;

      //-- find the maximum 3-bin likelihood for this generated set of observables.

         toy_3bin_lh_max = 0. ;
         toy_3bin_lh_max_mu_qcd_ldp = 0. ;
         toy_3bin_lh_max_mu_ttwj_sl = 0. ;
         toy_3bin_lh_max_n_zl = 0. ;

         for ( int qbi=0; qbi<nbins_mu_qcd_ldp; qbi++ ) {
            double lhs_mu_qcd_ldp = min_mu_qcd_ldp + (qbi+0.5)*(max_mu_qcd_ldp-min_mu_qcd_ldp)/nbins_mu_qcd_ldp ;
            for ( int tbi=0; tbi<nbins_mu_ttwj_sl; tbi++ ) {
               double lhs_mu_ttwj_sl = min_mu_ttwj_sl + (tbi+0.5)*(max_mu_ttwj_sl-min_mu_ttwj_sl)/nbins_mu_ttwj_sl ;
               double lhs_mu_ttwj_zl = Sttwj * Rttwj_zloversl * lhs_mu_ttwj_sl ;
               double lhs_mu_ttwj_ldp = Rttwj_ldpoverzl * lhs_mu_ttwj_zl ;
               double lhs_mu_qcd_zl = Sqcd * K3qcd * lhs_mu_qcd_ldp ;
               double lhs_n_ldp = etrig_had * lhs_mu_qcd_ldp + etrig_sl * ( lhs_mu_ttwj_ldp + mu_znn_ldp ) ;
               double lhs_n_sl = etrig_sl * lhs_mu_ttwj_sl ;
               double lhs_n_zl = etrig_sl * lhs_mu_ttwj_zl + etrig_had * lhs_mu_qcd_zl + etrig_sl * mu_znn_zl ;
               double lhs_lh_3bin = TMath::Poisson( Nzl_gen, lhs_n_zl) * TMath::Poisson( Nsl_gen, lhs_n_sl) * TMath::Poisson( Nldp_gen, lhs_n_ldp ) ;
               if ( lhs_lh_3bin > toy_3bin_lh_max ) {
                  toy_3bin_lh_max = lhs_lh_3bin ;
                  toy_3bin_lh_max_mu_qcd_ldp = lhs_mu_qcd_ldp ;
                  toy_3bin_lh_max_mu_ttwj_sl = lhs_mu_ttwj_sl ;
                  toy_3bin_lh_max_n_zl = lhs_n_zl ;
               }
            } // tbi
         } // qbi


         toy_3bin_nll_min = -log(toy_3bin_lh_max) ;
         h_toy_3bin_max_lh->Fill( toy_3bin_nll_min ) ;

         tt.Fill() ;

      } // ti.

      TCanvas* can2 = new TCanvas("can2","can2", 500, 500 ) ;
      hlh_2bin_gen_data_obs -> Draw("colz") ;
      can2->Draw() ;

      gStyle->SetOptStat("emrou") ;
      h_Nsl_gen->UseCurrentStyle() ;
      h_Nzl_gen->UseCurrentStyle() ;
      h_Nldp_gen->UseCurrentStyle() ;
      TCanvas* can3 = new TCanvas("can3","can3", 1000, 400 ) ;
      can3->Divide(3,1) ;

      can3->cd(1) ;
      h_Nsl_gen->SetFillColor(11) ;
      h_Nsl_gen->Draw() ;

      can3->cd(2) ;
      h_Nldp_gen->SetFillColor(11) ;
      h_Nldp_gen->Draw() ;

      can3->cd(3) ;
      h_Nzl_gen->SetFillColor(11) ;
      h_Nzl_gen->Draw() ;




      TCanvas* can4 = new TCanvas("can4","can4", 500, 400 ) ;
      h_toy_3bin_max_lh->UseCurrentStyle() ;
      h_toy_3bin_max_lh->SetFillColor(11) ;
      h_toy_3bin_max_lh->Draw() ;
      can4->Draw() ;


      outfile.cd() ;
      tt.Write() ;
      hlh_2bin_data_obs -> Write() ;
      hlh_3bin_data_obs -> Write() ;
      hlh_2bin_gen_data_obs -> Write() ;
      h_Nsl_gen -> Write() ;
      h_Nzl_gen -> Write() ;
      h_Nldp_gen -> Write() ;
      h_toy_3bin_max_lh->Write() ;
      outfile.Close() ;





   } // pvaltoy2



