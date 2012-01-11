#ifndef ra2bRoostatsClass9_h
#define ra2bRoostatsClass9_h

#include "TRandom2.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

//
//   Owen Long, UCR
//   Harrison Prosper, FSU
//   Sezen Sekmen, FSU
//


   class ra2bRoostatsClass9 {

     public :

       ra2bRoostatsClass9( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

       virtual ~ra2bRoostatsClass9();

       bool initialize( const char* infile = "an-11-257-v2-files/input-files/byhand-data-ge2b-loose.txt",
                        const char* inputScanFile = "an-11-257-v3-files/input-files/signalSyst.T1bbbb.ge1bTight.dat",
                        double m0 = 875., double m12 = 525., bool isT1bbbb = false, double t1bbbbXsec=0.,
                        const char* inputSusy_deff_dbtageff_file = "paper-2011-files/input-files/likelihood-newfit-syst-deff_dbtageff-LM9-HT400-SIGMET250.txt"
                        ) ;
       bool setSusyScanPoint( const char* inputScanFile,
                              double m0, double m12, bool isT1bbbb = false, double t1bbbbXsec=0.,
                              const char* inputSusy_deff_dbtageff_file = "paper-2011-files/input-files/likelihood-newfit-syst-deff_dbtageff-LM9-HT400-SIGMET250.txt"
                            ) ;

     private :

       char initializeFile[10000] ;

       bool varsAtFitVals ;
       bool useSigTtwjVar ;
       bool useLdpVars ;
       bool initialized ;

       double toy_mu0_ttbar_sig ;
       double toy_mu0_qcd_sig ;
       double toy_mu0_ttbar_sb ;
       double toy_mu0_qcd_sb ;
       double toy_mu0_susy_sig ;
       double toy_mu0_allbg_sig ;


      //======= Observables ==============================================================

       //-- data counts in signal region, A, and D, SB, SIG-SL, SB-SL.

       RooRealVar* rv_Nsig_1b        ;
       RooRealVar* rv_Nsb_1b         ;
       RooRealVar* rv_Nsig_sl_1b     ;
       RooRealVar* rv_Nsb_sl_1b      ;
       RooRealVar* rv_Nsig_ldp_1b    ;
       RooRealVar* rv_Nsb_ldp_1b     ;

       RooRealVar* rv_Nsig_2b        ;
       RooRealVar* rv_Nsb_2b         ;
       RooRealVar* rv_Nsig_sl_2b     ;
       RooRealVar* rv_Nsb_sl_2b      ;
       RooRealVar* rv_Nsig_ldp_2b    ;
       RooRealVar* rv_Nsb_ldp_2b     ;

       RooRealVar* rv_Nsig_3b        ;
       RooRealVar* rv_Nsb_3b         ;
       RooRealVar* rv_Nsig_sl_3b     ;
       RooRealVar* rv_Nsb_sl_3b      ;
       RooRealVar* rv_Nsig_ldp_3b    ;
       RooRealVar* rv_Nsb_ldp_3b     ;

       RooRealVar* rv_Nsb_ee      ;
       RooRealVar* rv_Nsig_ee     ;
       RooRealVar* rv_Nsb_mm      ;
       RooRealVar* rv_Nsig_mm     ;





       //========= Parameters ==============================================================


      //-- ttwj

       RooRealVar*    rrv_mu_ttwj_sig_1b ;
       RooRealVar*    rrv_mu_ttwj_sb_1b  ;

       RooFormulaVar* rfv_mu_ttwj_sig_1b ;
       RooFormulaVar* rfv_mu_ttwj_sb_1b  ;

       RooAbsArg*     rv_mu_ttwj_sig_1b ;
       RooAbsArg*     rv_mu_ttwj_sb_1b  ;

       RooRealVar*    rv_mu_ttwj_sig_sl_1b ;
       RooRealVar*    rv_mu_ttwj_sb_sl_1b  ;
       RooFormulaVar* rv_mu_ttwj_sig_ldp_1b ;
       RooFormulaVar* rv_mu_ttwj_sb_ldp_1b  ;




       RooRealVar*    rrv_mu_ttwj_sig_2b ;
       RooRealVar*    rrv_mu_ttwj_sb_2b  ;

       RooFormulaVar* rfv_mu_ttwj_sig_2b ;
       RooFormulaVar* rfv_mu_ttwj_sb_2b  ;

       RooAbsArg*     rv_mu_ttwj_sig_2b ;
       RooAbsArg*     rv_mu_ttwj_sb_2b  ;

       RooRealVar*    rv_mu_ttwj_sig_sl_2b ;
       RooRealVar*    rv_mu_ttwj_sb_sl_2b  ;
       RooFormulaVar* rv_mu_ttwj_sig_ldp_2b ;
       RooFormulaVar* rv_mu_ttwj_sb_ldp_2b  ;




       RooRealVar*    rrv_mu_ttwj_sig_3b ;
       RooRealVar*    rrv_mu_ttwj_sb_3b  ;

       RooFormulaVar* rfv_mu_ttwj_sig_3b ;
       RooFormulaVar* rfv_mu_ttwj_sb_3b  ;

       RooAbsArg*     rv_mu_ttwj_sig_3b ;
       RooAbsArg*     rv_mu_ttwj_sb_3b  ;

       RooRealVar*    rv_mu_ttwj_sig_sl_3b ;
       RooRealVar*    rv_mu_ttwj_sb_sl_3b  ;
       RooFormulaVar* rv_mu_ttwj_sig_ldp_3b ;
       RooFormulaVar* rv_mu_ttwj_sb_ldp_3b  ;




      //-- susy

       RooRealVar*    rv_mu_susy_sig_1b     ;
       RooFormulaVar* rv_mu_susy_sb_1b      ;
       RooFormulaVar* rv_mu_susy_sig_sl_1b  ;
       RooFormulaVar* rv_mu_susy_sb_sl_1b   ;
       RooFormulaVar* rv_mu_susy_sig_ldp_1b ;
       RooFormulaVar* rv_mu_susy_sb_ldp_1b  ;



       RooFormulaVar* rv_mu_susy_sig_2b     ;
       RooFormulaVar* rv_mu_susy_sb_2b      ;
       RooFormulaVar* rv_mu_susy_sig_sl_2b  ;
       RooFormulaVar* rv_mu_susy_sb_sl_2b   ;
       RooFormulaVar* rv_mu_susy_sig_ldp_2b ;
       RooFormulaVar* rv_mu_susy_sb_ldp_2b  ;



       RooFormulaVar* rv_mu_susy_sig_3b     ;
       RooFormulaVar* rv_mu_susy_sb_3b      ;
       RooFormulaVar* rv_mu_susy_sig_sl_3b  ;
       RooFormulaVar* rv_mu_susy_sb_sl_3b   ;
       RooFormulaVar* rv_mu_susy_sig_ldp_3b ;
       RooFormulaVar* rv_mu_susy_sb_ldp_3b  ;



      //-- QCD

       RooRealVar*    rrv_mu_qcd_sig_1b       ;
       RooRealVar*    rrv_mu_qcd_sb_1b        ;
       RooRealVar*    rrv_mu_qcd_sig_ldp_1b   ;
       RooRealVar*    rrv_mu_qcd_sb_ldp_1b    ;

       RooFormulaVar* rfv_mu_qcd_sig_1b       ;
       RooFormulaVar* rfv_mu_qcd_sb_1b        ;
       RooFormulaVar* rfv_mu_qcd_sig_ldp_1b   ;
       RooFormulaVar* rfv_mu_qcd_sb_ldp_1b    ;

       RooAbsArg*     rv_mu_qcd_sig_1b        ;
       RooAbsArg*     rv_mu_qcd_sb_1b         ;
       RooAbsArg*     rv_mu_qcd_sig_ldp_1b    ;
       RooAbsArg*     rv_mu_qcd_sb_ldp_1b     ;




       RooRealVar*    rrv_mu_qcd_sig_2b       ;
       RooRealVar*    rrv_mu_qcd_sb_2b        ;
       RooRealVar*    rrv_mu_qcd_sig_ldp_2b   ;
       RooRealVar*    rrv_mu_qcd_sb_ldp_2b    ;

       RooFormulaVar* rfv_mu_qcd_sig_2b       ;
       RooFormulaVar* rfv_mu_qcd_sb_2b        ;
       RooFormulaVar* rfv_mu_qcd_sig_ldp_2b   ;
       RooFormulaVar* rfv_mu_qcd_sb_ldp_2b    ;

       RooAbsArg*     rv_mu_qcd_sig_2b        ;
       RooAbsArg*     rv_mu_qcd_sb_2b         ;
       RooAbsArg*     rv_mu_qcd_sig_ldp_2b    ;
       RooAbsArg*     rv_mu_qcd_sb_ldp_2b     ;




       RooRealVar*    rrv_mu_qcd_sig_3b       ;
       RooRealVar*    rrv_mu_qcd_sb_3b        ;
       RooRealVar*    rrv_mu_qcd_sig_ldp_3b   ;
       RooRealVar*    rrv_mu_qcd_sb_ldp_3b    ;

       RooFormulaVar* rfv_mu_qcd_sig_3b       ;
       RooFormulaVar* rfv_mu_qcd_sb_3b        ;
       RooFormulaVar* rfv_mu_qcd_sig_ldp_3b   ;
       RooFormulaVar* rfv_mu_qcd_sb_ldp_3b    ;

       RooAbsArg*     rv_mu_qcd_sig_3b        ;
       RooAbsArg*     rv_mu_qcd_sb_3b         ;
       RooAbsArg*     rv_mu_qcd_sig_ldp_3b    ;
       RooAbsArg*     rv_mu_qcd_sb_ldp_3b     ;








      //-- Z to nunu, general
       RooRealVar*    rv_mu_znn_sig_1b     ;
       RooAbsArg*     rv_mu_znn_sb_1b ;
       RooFormulaVar* rv_mu_znn_sig_ldp_1b ;
       RooFormulaVar* rv_mu_znn_sb_ldp_1b  ;
       RooRealVar*    rrv_mu_znn_sb_1b      ;
       RooFormulaVar* rfv_mu_znn_sb_1b ;



       RooFormulaVar* rv_mu_znn_sig_2b     ;

       RooAbsArg*     rv_mu_znn_sb_2b ;
       RooFormulaVar* rv_mu_znn_sig_ldp_2b ;
       RooFormulaVar* rv_mu_znn_sb_ldp_2b  ;
       RooRealVar*    rrv_mu_znn_sb_2b      ;
       RooFormulaVar* rfv_mu_znn_sb_2b ;



       RooFormulaVar* rv_mu_znn_sig_3b     ;

       RooAbsArg*     rv_mu_znn_sb_3b ;
       RooFormulaVar* rv_mu_znn_sig_ldp_3b ;
       RooFormulaVar* rv_mu_znn_sb_ldp_3b  ;
       RooRealVar*    rrv_mu_znn_sb_3b      ;
       RooFormulaVar* rfv_mu_znn_sb_3b ;




       RooFormulaVar* rv_mu_zee_sig ;
       RooFormulaVar* rv_mu_zee_sb  ;
       RooFormulaVar* rv_mu_zmm_sig ;
       RooFormulaVar* rv_mu_zmm_sb  ;










      //-- Gaussian constraints and constants.


       RooFormulaVar* rv_eff_sf_sig_1b ;
       RooFormulaVar* rv_eff_sf_sb_1b ;
       RooFormulaVar* rv_eff_sf_sig_sl_1b ;
       RooFormulaVar* rv_eff_sf_sb_sl_1b ;
       RooFormulaVar* rv_eff_sf_sig_ldp_1b ;
       RooFormulaVar* rv_eff_sf_sb_ldp_1b ;

       RooFormulaVar* rv_eff_sf_sig_2b ;
       RooFormulaVar* rv_eff_sf_sb_2b ;
       RooFormulaVar* rv_eff_sf_sig_sl_2b ;
       RooFormulaVar* rv_eff_sf_sb_sl_2b ;
       RooFormulaVar* rv_eff_sf_sig_ldp_2b ;
       RooFormulaVar* rv_eff_sf_sb_ldp_2b ;

       RooFormulaVar* rv_eff_sf_sig_3b ;
       RooFormulaVar* rv_eff_sf_sb_3b ;
       RooFormulaVar* rv_eff_sf_sig_sl_3b ;
       RooFormulaVar* rv_eff_sf_sb_sl_3b ;
       RooFormulaVar* rv_eff_sf_sig_ldp_3b ;
       RooFormulaVar* rv_eff_sf_sb_ldp_3b ;



       RooRealVar* rv_width_eff_sf_sig_1b ;
       RooRealVar* rv_width_eff_sf_sb_1b ;
       RooRealVar* rv_width_eff_sf_sig_sl_1b ;
       RooRealVar* rv_width_eff_sf_sb_sl_1b ;
       RooRealVar* rv_width_eff_sf_sig_ldp_1b ;
       RooRealVar* rv_width_eff_sf_sb_ldp_1b ;

       RooRealVar* rv_width_eff_sf_sig_2b ;
       RooRealVar* rv_width_eff_sf_sb_2b ;
       RooRealVar* rv_width_eff_sf_sig_sl_2b ;
       RooRealVar* rv_width_eff_sf_sb_sl_2b ;
       RooRealVar* rv_width_eff_sf_sig_ldp_2b ;
       RooRealVar* rv_width_eff_sf_sb_ldp_2b ;

       RooRealVar* rv_width_eff_sf_sig_3b ;
       RooRealVar* rv_width_eff_sf_sb_3b ;
       RooRealVar* rv_width_eff_sf_sig_sl_3b ;
       RooRealVar* rv_width_eff_sf_sb_sl_3b ;
       RooRealVar* rv_width_eff_sf_sig_ldp_3b ;
       RooRealVar* rv_width_eff_sf_sb_ldp_3b ;




       RooRealVar* rv_mean_eff_sf_sig_1b ;
       RooRealVar* rv_mean_eff_sf_sb_1b ;
       RooRealVar* rv_mean_eff_sf_sig_sl_1b ;
       RooRealVar* rv_mean_eff_sf_sb_sl_1b ;
       RooRealVar* rv_mean_eff_sf_sig_ldp_1b ;
       RooRealVar* rv_mean_eff_sf_sb_ldp_1b ;

       RooRealVar* rv_mean_eff_sf_sig_2b ;
       RooRealVar* rv_mean_eff_sf_sb_2b ;
       RooRealVar* rv_mean_eff_sf_sig_sl_2b ;
       RooRealVar* rv_mean_eff_sf_sb_sl_2b ;
       RooRealVar* rv_mean_eff_sf_sig_ldp_2b ;
       RooRealVar* rv_mean_eff_sf_sb_ldp_2b ;

       RooRealVar* rv_mean_eff_sf_sig_3b ;
       RooRealVar* rv_mean_eff_sf_sb_3b ;
       RooRealVar* rv_mean_eff_sf_sig_sl_3b ;
       RooRealVar* rv_mean_eff_sf_sb_sl_3b ;
       RooRealVar* rv_mean_eff_sf_sig_ldp_3b ;
       RooRealVar* rv_mean_eff_sf_sb_ldp_3b ;





       RooRealVar* rv_deff_dbtageff_sig_1b ;
       RooRealVar* rv_deff_dbtageff_sb_1b ;
       RooRealVar* rv_deff_dbtageff_sig_sl_1b ;
       RooRealVar* rv_deff_dbtageff_sb_sl_1b ;
       RooRealVar* rv_deff_dbtageff_sig_ldp_1b ;
       RooRealVar* rv_deff_dbtageff_sb_ldp_1b ;

       RooRealVar* rv_deff_dbtageff_sig_2b ;
       RooRealVar* rv_deff_dbtageff_sb_2b ;
       RooRealVar* rv_deff_dbtageff_sig_sl_2b ;
       RooRealVar* rv_deff_dbtageff_sb_sl_2b ;
       RooRealVar* rv_deff_dbtageff_sig_ldp_2b ;
       RooRealVar* rv_deff_dbtageff_sb_ldp_2b ;

       RooRealVar* rv_deff_dbtageff_sig_3b ;
       RooRealVar* rv_deff_dbtageff_sb_3b ;
       RooRealVar* rv_deff_dbtageff_sig_sl_3b ;
       RooRealVar* rv_deff_dbtageff_sb_sl_3b ;
       RooRealVar* rv_deff_dbtageff_sig_ldp_3b ;
       RooRealVar* rv_deff_dbtageff_sb_ldp_3b ;




       RooFormulaVar* rv_btageff_sf_sig_1b ;
       RooFormulaVar* rv_btageff_sf_sb_1b ;
       RooFormulaVar* rv_btageff_sf_sig_sl_1b ;
       RooFormulaVar* rv_btageff_sf_sb_sl_1b ;
       RooFormulaVar* rv_btageff_sf_sig_ldp_1b ;
       RooFormulaVar* rv_btageff_sf_sb_ldp_1b ;

       RooFormulaVar* rv_btageff_sf_sig_2b ;
       RooFormulaVar* rv_btageff_sf_sb_2b ;
       RooFormulaVar* rv_btageff_sf_sig_sl_2b ;
       RooFormulaVar* rv_btageff_sf_sb_sl_2b ;
       RooFormulaVar* rv_btageff_sf_sig_ldp_2b ;
       RooFormulaVar* rv_btageff_sf_sb_ldp_2b ;

       RooFormulaVar* rv_btageff_sf_sig_3b ;
       RooFormulaVar* rv_btageff_sf_sb_3b ;
       RooFormulaVar* rv_btageff_sf_sig_sl_3b ;
       RooFormulaVar* rv_btageff_sf_sb_sl_3b ;
       RooFormulaVar* rv_btageff_sf_sig_ldp_3b ;
       RooFormulaVar* rv_btageff_sf_sb_ldp_3b ;


       RooRealVar* rv_znnoverll_bfratio    ;
       RooRealVar* rv_dataoverll_lumiratio ;


    //++++ MC inputs +++++++++++++++++++++++

       RooRealVar*    rv_mu_susymc_sig_1b      ;
       RooRealVar*    rv_mu_susymc_sb_1b       ;
       RooRealVar*    rv_mu_susymc_sig_sl_1b   ;
       RooRealVar*    rv_mu_susymc_sb_sl_1b    ;
       RooRealVar*    rv_mu_susymc_sig_ldp_1b  ;
       RooRealVar*    rv_mu_susymc_sb_ldp_1b   ;

       RooRealVar*    rv_mu_susymc_sig_2b      ;
       RooRealVar*    rv_mu_susymc_sb_2b       ;
       RooRealVar*    rv_mu_susymc_sig_sl_2b   ;
       RooRealVar*    rv_mu_susymc_sb_sl_2b    ;
       RooRealVar*    rv_mu_susymc_sig_ldp_2b  ;
       RooRealVar*    rv_mu_susymc_sb_ldp_2b   ;

       RooRealVar*    rv_mu_susymc_sig_3b      ;
       RooRealVar*    rv_mu_susymc_sb_3b       ;
       RooRealVar*    rv_mu_susymc_sig_sl_3b   ;
       RooRealVar*    rv_mu_susymc_sb_sl_3b    ;
       RooRealVar*    rv_mu_susymc_sig_ldp_3b  ;
       RooRealVar*    rv_mu_susymc_sb_ldp_3b   ;





     //-- SIG, LDP

       RooRealVar*    rv_mu_ttbarsingletopzjetsmc_sig_ldp_1b ;
       RooRealVar*    rv_mu_WJmc_sig_ldp_1b    ;
       RooRealVar*    rv_mu_Znnmc_sig_ldp_1b   ;

       RooRealVar*    rv_mu_ttbarsingletopzjetsmc_sig_ldp_2b ;
       RooRealVar*    rv_mu_WJmc_sig_ldp_2b    ;
       RooRealVar*    rv_mu_Znnmc_sig_ldp_2b   ;

       RooRealVar*    rv_mu_ttbarsingletopzjetsmc_sig_ldp_3b ;
       RooRealVar*    rv_mu_WJmc_sig_ldp_3b    ;
       RooRealVar*    rv_mu_Znnmc_sig_ldp_3b   ;




     //-- SB, LDP

       RooRealVar*    rv_mu_ttbarsingletopzjetsmc_sb_ldp_1b  ;
       RooRealVar*    rv_mu_WJmc_sb_ldp_1b     ;
       RooRealVar*    rv_mu_Znnmc_sb_ldp_1b    ;

       RooRealVar*    rv_mu_ttbarsingletopzjetsmc_sb_ldp_2b  ;
       RooRealVar*    rv_mu_WJmc_sb_ldp_2b     ;
       RooRealVar*    rv_mu_Znnmc_sb_ldp_2b    ;

       RooRealVar*    rv_mu_ttbarsingletopzjetsmc_sb_ldp_3b  ;
       RooRealVar*    rv_mu_WJmc_sb_ldp_3b     ;
       RooRealVar*    rv_mu_Znnmc_sb_ldp_3b    ;




       //========= Expected counts for observables in terms of parameters ============================

       RooFormulaVar* rv_n_sig_1b        ;
       RooFormulaVar* rv_n_sb_1b         ;
       RooFormulaVar* rv_n_sig_ldp_1b    ;
       RooFormulaVar* rv_n_sb_ldp_1b     ;
       RooFormulaVar* rv_n_sig_sl_1b     ;
       RooFormulaVar* rv_n_sb_sl_1b      ;

       RooFormulaVar* rv_n_sig_2b        ;
       RooFormulaVar* rv_n_sb_2b         ;
       RooFormulaVar* rv_n_sig_ldp_2b    ;
       RooFormulaVar* rv_n_sb_ldp_2b     ;
       RooFormulaVar* rv_n_sig_sl_2b     ;
       RooFormulaVar* rv_n_sb_sl_2b      ;

       RooFormulaVar* rv_n_sig_3b        ;
       RooFormulaVar* rv_n_sb_3b         ;
       RooFormulaVar* rv_n_sig_ldp_3b    ;
       RooFormulaVar* rv_n_sb_ldp_3b     ;
       RooFormulaVar* rv_n_sig_sl_3b     ;
       RooFormulaVar* rv_n_sb_sl_3b      ;

       RooFormulaVar* rv_n_sig_ee     ;
       RooFormulaVar* rv_n_sb_ee      ;
       RooFormulaVar* rv_n_sig_mm     ;
       RooFormulaVar* rv_n_sb_mm      ;







       //=========== PDFs for the likelihood ============================================================

       RooPoisson*  pdf_Nsig_1b        ;
       RooPoisson*  pdf_Nsb_1b         ;
       RooPoisson*  pdf_Nsig_ldp_1b    ;
       RooPoisson*  pdf_Nsb_ldp_1b     ;
       RooPoisson*  pdf_Nsig_sl_1b     ;
       RooPoisson*  pdf_Nsb_sl_1b      ;

       RooPoisson*  pdf_Nsig_2b        ;
       RooPoisson*  pdf_Nsb_2b         ;
       RooPoisson*  pdf_Nsig_ldp_2b    ;
       RooPoisson*  pdf_Nsb_ldp_2b     ;
       RooPoisson*  pdf_Nsig_sl_2b     ;
       RooPoisson*  pdf_Nsb_sl_2b      ;

       RooPoisson*  pdf_Nsig_3b        ;
       RooPoisson*  pdf_Nsb_3b         ;
       RooPoisson*  pdf_Nsig_ldp_3b    ;
       RooPoisson*  pdf_Nsb_ldp_3b     ;
       RooPoisson*  pdf_Nsig_sl_3b     ;
       RooPoisson*  pdf_Nsb_sl_3b      ;


       RooPoisson*  pdf_Nsig_ee     ;
       RooPoisson*  pdf_Nsb_ee      ;
       RooPoisson*  pdf_Nsig_mm     ;
       RooPoisson*  pdf_Nsb_mm      ;

       RooProdPdf*  likelihood ;



       //============ Other things =================================================================


       float EffScaleFactor, EffScaleFactorErr ;

       float lsf_WJmc, lsf_WJmc_err ;
       float lsf_Znnmc, lsf_Znnmc_err ;
       float sf_ttbarmc, sf_ttbarmc_err ;

       float  DataLumi                ;
       float  acc_ee_sig_mean             ;
       float  acc_ee_sig_err              ;
       float  acc_ee_sb_mean             ;
       float  acc_ee_sb_err              ;
       float  acc_mm_sig_mean             ;
       float  acc_mm_sig_err              ;
       float  acc_mm_sb_mean             ;
       float  acc_mm_sb_err              ;
       float  eff_ee_mean             ;
       float  eff_ee_err              ;
       float  eff_mm_mean             ;
       float  eff_mm_err              ;
       float  knn_sig_1b_mean            ;
       float  knn_sig_1b_err             ;
       float  knn_sig_2b_mean            ;
       float  knn_sig_2b_err             ;
       float  knn_sig_3b_mean            ;
       float  knn_sig_3b_err             ;
       float  knn_sb_1b_mean             ;
       float  knn_sb_1b_err              ;
       float  knn_sb_2b_mean             ;
       float  knn_sb_2b_err              ;
       float  knn_sb_3b_mean             ;
       float  knn_sb_3b_err              ;
       float  Ztoll_lumi              ;
       float  fsig_ee_mean            ;
       float  fsig_ee_err             ;
       float  fsig_mm_mean            ;
       float  fsig_mm_err             ;


       RooArgSet observedParametersList ;
       RooDataSet* dsObserved ;


       RooFitResult* fitResult ;




   } ;




#endif
