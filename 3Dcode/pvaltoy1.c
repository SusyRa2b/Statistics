
#include "TRandom.h"
#include "TH1F.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TArrow.h"

   void pvaltoy( int nToy = 1000, double mu_sl_true = 0.99, double mu_smear = 0.1 ) {

    //--- M4_H4_3b

      int N_sl = 0 ;
      int N_ldp = 0 ;
      int N_zl = 4 ;

      ////double R_ldp = 1.14 ;
      ////double R_zl = 1.68 ;
      ////double N_Znn_zl = 0.4 ;

      double R_ldp = 1.06 ;
      double R_zl = 1.26 ;
      double N_Znn_zl = 0.3 ;



      TRandom rangen(12345) ;

      double mu_sl_data = 0.99 ;

      gDirectory->Delete("h*") ;

      int nbins(60) ;
      TH1F* h_best_mu = new TH1F("h_best_mu", "best mu", nbins, 0., 10. ) ;
      TH1F* h_logL = new TH1F("h_logL", "-lnL", nbins, 0., 10. ) ;
      TH1F* h_Ngen_sl = new TH1F("h_Ngen_sl", "Ngen, SL", 11, -0.5, 10.5 ) ;
      TH1F* h_Ngen_zl = new TH1F("h_Ngen_zl", "Ngen, ZL", 11, -0.5, 10.5 ) ;
      TH1F* h_Ngen_ldp = new TH1F("h_Ngen_ldp", "Ngen, LDP", 11, -0.5, 10.5 ) ;
      TH1F* h_gen_mu = new TH1F("h_gen_mu", "generated mu", 40, 0., 5. ) ;

      double scan_mu_max(10.) ;
      double scan_deltamu(0.01) ;
      int nScan = TMath::Nint( scan_mu_max / scan_deltamu ) ;

      for ( int ti=0; ti<nToy; ti++ ) {

     //------------------
         double mu_sl_gen = rangen.Gaus(mu_sl_true, mu_smear) ;
         if ( mu_sl_gen < 0. ) continue ;
     //------------------
      // /////// double mu_sl_gen = rangen.Exp( 0.467 ) ;  // 1./(1+1.14) = 0.467
      // double mu_sl_gen = rangen.Exp( 0.546 ) ;  // 1./(1+0.83) = 0.546
     //------------------

         h_gen_mu -> Fill( mu_sl_gen ) ;


         int gen_N_sl  = rangen.Poisson( mu_sl_gen ) ;
         int gen_N_ldp = rangen.Poisson( R_ldp * mu_sl_gen ) ;
         int gen_N_zl  = rangen.Poisson( N_Znn_zl + R_zl * mu_sl_gen ) ;

         h_Ngen_sl -> Fill( gen_N_sl ) ;
         h_Ngen_zl -> Fill( gen_N_zl ) ;
         h_Ngen_ldp -> Fill( gen_N_ldp ) ;

         double minNll(1.e9) ;
         double best_mu_sl(0.) ;
         for ( int si=0; si<nScan; si++ ) {
            double mu_sl = si * scan_deltamu ;
            double nll = -1. * log( TMath::Poisson( gen_N_sl, mu_sl )
                                   * TMath::Poisson( gen_N_ldp, R_ldp * mu_sl )
                                   * TMath::Poisson( gen_N_zl, N_Znn_zl + R_zl * mu_sl ) ) ;
            if ( nll < minNll ) {
               minNll = nll ;
               best_mu_sl = mu_sl ;
            }
         } // si.

         h_best_mu -> Fill( best_mu_sl ) ;
         h_logL -> Fill( minNll ) ;

      } // ti.

      double dataNll = -1. * log( TMath::Poisson( N_sl, mu_sl_data )
                                * TMath::Poisson( N_ldp, R_ldp * mu_sl_data )
                                * TMath::Poisson( N_zl, N_Znn_zl + R_zl * mu_sl_data ) ) ;

      printf("\n\n Data -lnL = %g\n\n", dataNll ) ;

      int thebin = h_logL->FindBin( dataNll ) ;
      int nworse = h_logL -> Integral( thebin, nbins ) ;
      printf(" pval = %d / %d = %g\n\n", nworse, nToy, (1.*nworse)/(1.*nToy) ) ;

      TH1F* hshaded = (TH1F*) h_logL->Clone("hshaded") ;
      for ( int bi=1; bi<thebin; bi++ ) { hshaded->SetBinContent(bi,0.) ; }
      hshaded->SetFillColor(4) ;




      TArrow* arrow = new TArrow() ;
      arrow->SetLineWidth(2) ;
      arrow->SetLineColor(4) ;


      TCanvas* cpval = (TCanvas*) gDirectory->FindObject("cpval") ;
      if ( cpval == 0x0 ) {
          cpval = new TCanvas("cpval", "pval", 800, 1000 ) ;
      }

      cpval->Clear() ;
      cpval->Divide(2,3) ;

      cpval->cd(1) ;
      h_best_mu->Draw() ;

      cpval->cd(2) ;
      h_logL->Draw() ;
      hshaded->Draw("same") ;

      cpval->cd(3) ;
      h_Ngen_sl->Draw() ;
      arrow->DrawArrow( N_sl, 0.4*h_Ngen_sl->GetMaximum(), N_sl, 0., 0.02 ) ;

      cpval->cd(4) ;
      h_Ngen_zl->Draw() ;
      arrow->DrawArrow( N_zl, 0.4*h_Ngen_zl->GetMaximum(), N_zl, 0., 0.02 ) ;

      cpval->cd(5) ;
      h_Ngen_ldp->Draw() ;
      arrow->DrawArrow( N_ldp, 0.4*h_Ngen_ldp->GetMaximum(), N_ldp, 0., 0.02 ) ;

      cpval->cd(6) ;
      h_gen_mu->Draw() ;

      cpval->Update() ;
      cpval->Draw() ;

   }

