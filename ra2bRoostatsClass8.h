#ifndef ra2bRoostatsClass8_h
#define ra2bRoostatsClass8_h

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


   class ra2bRoostatsClass8 {

     public :

       ra2bRoostatsClass8( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

       virtual ~ra2bRoostatsClass8();

       bool initialize( const char* infile = "an-11-257-v2-files/input-files/byhand-data-ge2b-loose.txt",
                        const char* inputScanFile = "an-11-257-v3-files/input-files/signalSyst.T1bbbb.ge1bTight.dat", double m0 = 875., double m12 = 525., bool isT1bbbb = false, double t1bbbbXsec=0. ) ;
       bool setSusyScanPoint( const char* inputScanFile, double m0, double m12, bool isT1bbbb = false, double t1bbbbXsec=0. ) ;

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

       RooRealVar* rv_Nsig        ;
       RooRealVar* rv_Nsb         ;
       RooRealVar* rv_Nsig_sl     ;
       RooRealVar* rv_Nsb_sl      ;
       RooRealVar* rv_Nsig_ldp    ;
       RooRealVar* rv_Nsb_ldp     ;
       RooRealVar* rv_Nlsb_0b     ;
       RooRealVar* rv_Nlsb_0b_ldp ;

       //-- Znn model 1
       RooRealVar* rv_Nsb_ee      ;
       RooRealVar* rv_Nsig_ee     ;
       RooRealVar* rv_Nsb_mm      ;
       RooRealVar* rv_Nsig_mm     ;

       //-- Znn model 2
       RooRealVar* rv_Nsigsb_ee   ;
       RooRealVar* rv_Nsigsb_mm   ;





       //========= Parameters ==============================================================


      //-- ttwj

       RooRealVar*    rrv_mu_ttwj_sig ;
       RooRealVar*    rrv_mu_ttwj_sb  ;

       RooFormulaVar* rfv_mu_ttwj_sig ;
       RooFormulaVar* rfv_mu_ttwj_sb  ;

       RooAbsArg*     rv_mu_ttwj_sig ;
       RooAbsArg*     rv_mu_ttwj_sb  ;

       RooRealVar*    rv_mu_ttwj_sig_sl ;
       RooRealVar*    rv_mu_ttwj_sb_sl  ;
       RooFormulaVar* rv_mu_ttwj_sig_ldp ;
       RooFormulaVar* rv_mu_ttwj_sb_ldp  ;




      //-- susy

       RooRealVar*    rv_mu_susy_sig ;
       RooFormulaVar* rv_mu_susy_sb      ;
       RooFormulaVar* rv_mu_susy_sig_sl  ;
       RooFormulaVar* rv_mu_susy_sb_sl   ;
       RooFormulaVar* rv_mu_susy_sig_ldp ;
       RooFormulaVar* rv_mu_susy_sb_ldp  ;



      //-- QCD

       RooRealVar*    rrv_mu_qcd_sig       ;
       RooRealVar*    rrv_mu_qcd_sb        ;
       RooRealVar*    rrv_mu_qcd_sig_ldp   ;
       RooRealVar*    rrv_mu_qcd_sb_ldp    ;

       RooFormulaVar* rfv_mu_qcd_sig       ;
       RooFormulaVar* rfv_mu_qcd_sb        ;
       RooFormulaVar* rfv_mu_qcd_sig_ldp   ;
       RooFormulaVar* rfv_mu_qcd_sb_ldp    ;

       RooAbsArg*     rv_mu_qcd_sig        ;
       RooAbsArg*     rv_mu_qcd_sb         ;
       RooAbsArg*     rv_mu_qcd_sig_ldp    ;
       RooAbsArg*     rv_mu_qcd_sb_ldp     ;


   //  RooRealVar*    rv_mu_qcd_lsb_0b     ;
   //  RooRealVar*    rv_mu_qcd_lsb_0b_ldp ;






      //-- Z to nunu, general
       RooRealVar*    rv_mu_znn_sig     ;
       RooAbsArg*     rv_mu_znn_sb ;
       RooFormulaVar* rv_mu_znn_sig_ldp ;
       RooFormulaVar* rv_mu_znn_sb_ldp  ;


      //-- Z to nunu, model 1
       RooRealVar*    rrv_mu_znn_sb      ;
       RooFormulaVar* rv_mu_zee_sig_ee ;
       RooFormulaVar* rv_mu_zee_sb_ee  ;
       RooFormulaVar* rv_mu_zmm_sig_mm ;
       RooFormulaVar* rv_mu_zmm_sb_mm  ;

      //-- Z to nunu, model 2
       RooFormulaVar* rfv_mu_znn_sb ;
       RooFormulaVar* rv_mu_zee_sigsb_ee ;
       RooFormulaVar* rv_mu_zmm_sigsb_mm ;






      //-- EW0

       RooFormulaVar* rv_mu_ewo_sig ;
       RooFormulaVar* rv_mu_ewo_sb      ;
       RooFormulaVar* rv_mu_ewo_sig_ldp ;
       RooFormulaVar* rv_mu_ewo_sb_ldp  ;


      //-- Gaussian constraints and constants.

       RooRealVar*    rv_eff_sf_m ;
       RooFormulaVar* rv_eff_sf_sig ;
       RooFormulaVar* rv_eff_sf_sb ;
       RooFormulaVar* rv_eff_sf_sig_sl ;
       RooFormulaVar* rv_eff_sf_sb_sl ;
       RooFormulaVar* rv_eff_sf_sig_ldp ;
       RooFormulaVar* rv_eff_sf_sb_ldp ;

       RooRealVar* rv_width_eff_sf_sig ;
       RooRealVar* rv_width_eff_sf_sb ;
       RooRealVar* rv_width_eff_sf_sig_sl ;
       RooRealVar* rv_width_eff_sf_sb_sl ;
       RooRealVar* rv_width_eff_sf_sig_ldp ;
       RooRealVar* rv_width_eff_sf_sb_ldp ;

       RooRealVar* rv_mean_eff_sf_sig ;
       RooRealVar* rv_mean_eff_sf_sb ;
       RooRealVar* rv_mean_eff_sf_sig_sl ;
       RooRealVar* rv_mean_eff_sf_sb_sl ;
       RooRealVar* rv_mean_eff_sf_sig_ldp ;
       RooRealVar* rv_mean_eff_sf_sb_ldp ;



       RooRealVar* rv_znnoverll_bfratio    ;
       RooRealVar* rv_dataoverll_lumiratio ;


    //++++ MC inputs +++++++++++++++++++++++

       RooRealVar*    rv_mu_susymc_sig      ;
       RooRealVar*    rv_mu_susymc_sb       ;
       RooRealVar*    rv_mu_susymc_sig_sl   ;
       RooRealVar*    rv_mu_susymc_sb_sl    ;
       RooRealVar*    rv_mu_susymc_sig_ldp  ;
       RooRealVar*    rv_mu_susymc_sb_ldp   ;



     //-- SIG

       RooRealVar*    rv_mu_Ewomc_sig       ;

     //-- SB

       RooRealVar*    rv_mu_Ewomc_sb        ;


     //-- SIG, LDP

       RooRealVar*    rv_mu_ttbarmc_sig_ldp ;
       RooRealVar*    rv_mu_WJmc_sig_ldp    ;
       RooRealVar*    rv_mu_Znnmc_sig_ldp   ;
       RooRealVar*    rv_mu_Ewomc_sig_ldp   ;


     //-- SB, LDP

       RooRealVar*    rv_mu_ttbarmc_sb_ldp  ;
       RooRealVar*    rv_mu_WJmc_sb_ldp     ;
       RooRealVar*    rv_mu_Znnmc_sb_ldp    ;
       RooRealVar*    rv_mu_Ewomc_sb_ldp    ;




       //========= Expected counts for observables in terms of parameters ============================

       RooFormulaVar* rv_n_sig        ;
       RooFormulaVar* rv_n_sb         ;
       RooFormulaVar* rv_n_sig_ldp    ;
       RooFormulaVar* rv_n_sb_ldp     ;
       RooFormulaVar* rv_n_sig_sl     ;
       RooFormulaVar* rv_n_sb_sl      ;
       RooFormulaVar* rv_n_lsb_0b     ;
       RooFormulaVar* rv_n_lsb_0b_ldp ;

       //-- Znn model 1
       RooFormulaVar* rv_n_sig_ee     ;
       RooFormulaVar* rv_n_sb_ee      ;
       RooFormulaVar* rv_n_sig_mm     ;
       RooFormulaVar* rv_n_sb_mm      ;

       //-- Znn model 2
       RooFormulaVar* rv_n_sigsb_ee   ;
       RooFormulaVar* rv_n_sigsb_mm   ;




       //=========== Mean values and sigmas for nuisance parameters (things with Gaussian PDFs in likelihood).

       double np_m_Eff_sf_m ;
       double np_m_sf_mc ;
       double np_m_sf_qcd_sb ;
       double np_m_sf_qcd_sig ;
       double np_m_sf_ttwj_sig ;
       double np_m_sf_ee ;
       double np_m_sf_mm ;
       double np_m_acc_ee ;
       double np_m_acc_mm ;
       double np_m_eff_ee ;
       double np_m_eff_mm ;
       double np_m_fsig_ee ;
       double np_m_fsig_mm ;
       double np_m_knn_sig ;
       double np_m_knn_sb ;
       double np_m_Rlsb_passfail ;


       double np_s_Eff_sf_m ;
       double np_s_sf_mc ;
       double np_s_sf_qcd_sb ;
       double np_s_sf_qcd_sig ;
       double np_s_sf_ttwj_sig ;
       double np_s_sf_ee ;
       double np_s_sf_mm ;
       double np_s_acc_ee ;
       double np_s_acc_mm ;
       double np_s_eff_ee ;
       double np_s_eff_mm ;
       double np_s_fsig_ee ;
       double np_s_fsig_mm ;
       double np_s_knn_sig ;
       double np_s_knn_sb ;
       double np_s_Rlsb_passfail ;


       //=========== PDFs for the likelihood ============================================================

       RooPoisson*  pdf_Nsig        ;
       RooPoisson*  pdf_Nsb         ;
       RooPoisson*  pdf_Nsig_ldp    ;
       RooPoisson*  pdf_Nsb_ldp     ;
       RooPoisson*  pdf_Nsig_sl     ;
       RooPoisson*  pdf_Nsb_sl      ;
       RooPoisson*  pdf_Nlsb_0b     ;
       RooPoisson*  pdf_Nlsb_0b_ldp ;


       //-- Znn model 1
       RooPoisson*  pdf_Nsig_ee     ;
       RooPoisson*  pdf_Nsb_ee      ;
       RooPoisson*  pdf_Nsig_mm     ;
       RooPoisson*  pdf_Nsb_mm      ;

       //-- Znn model 2
       RooPoisson*  pdf_Nsigsb_ee   ;
       RooPoisson*  pdf_Nsigsb_mm   ;

       RooProdPdf*  likelihood ;



       //============ Other things =================================================================


       float EffScaleFactor, EffScaleFactorErr ;

       float lsf_WJmc, lsf_WJmc_err ;
       float lsf_Znnmc, lsf_Znnmc_err ;
       float lsf_Ewomc, lsf_Ewomc_err ;
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
       float  knn_ee_sig_mean            ;
       float  knn_ee_sig_err             ;
       float  knn_ee_sb_mean             ;
       float  knn_ee_sb_err              ;
       float  knn_mm_sig_mean            ;
       float  knn_mm_sig_err             ;
       float  knn_mm_sb_mean             ;
       float  knn_mm_sb_err              ;
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
