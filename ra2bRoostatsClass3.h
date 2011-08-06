#ifndef ra2bRoostatsClass3_h
#define ra2bRoostatsClass3_h

#include "Trandom2.h"
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


   class ra2bRoostatsClass3 {

     public :

       ra2bRoostatsClass3( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true, int ArgznnModel=1 ) ;

       virtual ~ra2bRoostatsClass3();

       bool initialize( const char* infile ) ;

       bool reinitialize( ) ;

       bool doFit() ;


       bool susyScanNoContam( const char* inputScanFile, const char* outputFilebase="output-files/scanplot-nocontam" ) ;
       bool susyScanWithContam( const char* inputScanFile, const char* outputFilebase="output-files/scanplot-withcontam"  ) ;
       bool discoveryScanWithContam( const char* inputScanFile ) ;
       bool nosusysignifScanWithContam( const char* inputScanFile ) ;


       bool setSusyScanPoint( const char* inputScanFile, double m0, double m12 ) ;


       bool profileSusySig( float& susySigLow, float& susySigHigh, bool makePlot=true, const char* plotname="output-files/prof_susy_sig.png", double scanMax=-1. ) ;

       bool profileZnnSig( float& znnSigLow, float& znnSigHigh, bool makePlot=true, const char* plotname="output-files/prof_znn_sig.png", double scanMax=-1. ) ;
       bool profileZnnSb( float& znnSbLow, float& znnSbHigh, bool makePlot=true, const char* plotname="output-files/prof_znn_sb.png", double scanMax=-1. ) ;

       bool profilettwjSig( float& ttwjSigLow, float& ttwjSigHigh, bool makePlot=true, const char* plotname="output-files/prof_ttwj_sig.png", double scanMax=-1. ) ;
       bool profilettwjSb( float& ttwjSbLow, float& ttwjSbHigh, bool makePlot=true, const char* plotname="output-files/prof_ttwj_sb.png", double scanMax=-1. ) ;

       bool profileqcdSig( float& qcdSigLow, float& qcdSigHigh, bool makePlot=true, const char* plotname="output-files/prof_qcd_sig.png", double scanMax=-1. ) ;
       bool profileqcdSb( float& qcdSbLow, float& qcdSbHigh, bool makePlot=true, const char* plotname="output-files/prof_qcd_sb.png", double scanMax=-1. ) ;



       void parameterSnapshot() ;

       bool fitQualityPlot( bool doNorm=false, const char* plotname="output-files/fit_qual.png", double hmax=1.5 ) ;

       void setAndFixSusySig( double setVal = 0. ) ;
       void freeSusySig() ;


      //--- things for CLs

       void saveToymeanSnapshot() ;

       void genToyExperiment() ;

       double doToyStudy( int nToys=1000, bool isBgonlyStudy=true, double data_q = 0. ) ;

       void setAndFixSusySigToPredictedValue() ;

       double getLogLikelihoodValue() ;

     private :

       char initializeFile[10000] ;

       bool varsAtFitVals ;
       bool useSigTtwjVar ;
       bool useLdpVars ;
       int  znnModel ;
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


       RooRealVar*    rv_mu_qcd_lsb_0b     ;
       RooRealVar*    rv_mu_qcd_lsb_0b_ldp ;






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

       RooRealVar* rv_eff_sf ;


       RooRealVar* rv_sf_mc ;
       RooRealVar* rv_sf_qcd_sb ;
       RooRealVar* rv_sf_qcd_sig ;
       RooRealVar* rv_sf_ttwj_sig ;
       RooRealVar* rv_sf_ee ;
       RooRealVar* rv_sf_mm ;

       RooRealVar* rv_acc_ee ;
       RooRealVar* rv_acc_mm ;
       RooRealVar* rv_eff_ee ;
       RooRealVar* rv_eff_mm ;

       RooRealVar* rv_fsig_ee ;
       RooRealVar* rv_fsig_mm ;

       RooRealVar* rv_knn_sig ;
       RooRealVar* rv_knn_sb ;

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






       //=========== PDFs for the likelihood ============================================================

       RooPoisson*  pdf_Nsig        ;
       RooPoisson*  pdf_Nsb         ;
       RooPoisson*  pdf_Nsig_ldp    ;
       RooPoisson*  pdf_Nsb_ldp     ;
       RooPoisson*  pdf_Nsig_sl     ;
       RooPoisson*  pdf_Nsb_sl      ;
       RooPoisson*  pdf_Nlsb_0b     ;
       RooPoisson*  pdf_Nlsb_0b_ldp ;

       RooGaussian* pdf_sf_mc ;
       RooGaussian* pdf_sf_qcd_sb ;
       RooGaussian* pdf_sf_qcd_sig ;
       RooGaussian* pdf_sf_ttwj_sig ;
       RooGaussian* pdf_sf_ee ;
       RooGaussian* pdf_sf_mm ;


       RooGaussian* pdf_acc_ee      ;
       RooGaussian* pdf_acc_mm      ;
       RooGaussian* pdf_eff_ee      ;
       RooGaussian* pdf_eff_mm      ;
       RooGaussian* pdf_fsig_ee      ;
       RooGaussian* pdf_fsig_mm      ;
       RooGaussian* pdf_Eff_sf      ;

       //-- Znn model 1
       RooPoisson*  pdf_Nsig_ee     ;
       RooPoisson*  pdf_Nsb_ee      ;
       RooPoisson*  pdf_Nsig_mm     ;
       RooPoisson*  pdf_Nsb_mm      ;

       //-- Znn model 2
       RooPoisson*  pdf_Nsigsb_ee   ;
       RooPoisson*  pdf_Nsigsb_mm   ;
       RooGaussian* pdf_knn_sig     ;
       RooGaussian* pdf_knn_sb      ;



       RooProdPdf*  likelihood ;



       //============ Other things =================================================================


       float EffScaleFactor, EffScaleFactorErr ;

       float lsf_WJmc, lsf_WJmc_err ;
       float lsf_Znnmc, lsf_Znnmc_err ;
       float lsf_Ewomc, lsf_Ewomc_err ;
       float sf_ttbarmc, sf_ttbarmc_err ;

       float  DataLumi                ;
       float  acc_ee_mean             ;
       float  acc_ee_err              ;
       float  acc_mm_mean             ;
       float  acc_mm_err              ;
       float  eff_ee_mean             ;
       float  eff_ee_err              ;
       float  eff_mm_mean             ;
       float  eff_mm_err              ;
       float  knn_sig_mean            ;
       float  knn_sig_err             ;
       float  knn_sb_mean             ;
       float  knn_sb_err              ;
       float  Ztoll_lumi              ;
       float  fsig_ee_mean            ;
       float  fsig_ee_err             ;
       float  fsig_mm_mean            ;
       float  fsig_mm_err             ;


       RooArgSet observedParametersList ;
       RooDataSet* dsObserved ;

       RooWorkspace* workspace ;

       RooFitResult* fitResult ;


       //========= Stuff needed for CLs below here ==================================================

     //++++ Model predictions from fitting data to BG-only hypothesis (susy fixed to zero in fit)

       double toymean_n_sig        ;
       double toymean_n_sb         ;
       double toymean_n_sig_ldp    ;
       double toymean_n_sb_ldp     ;
       double toymean_n_sig_sl     ;
       double toymean_n_sb_sl      ;
       double toymean_n_lsb_0b     ;
       double toymean_n_lsb_0b_ldp ;

       //-- Znn model 1
       double toymean_n_sig_ee     ;
       double toymean_n_sb_ee      ;
       double toymean_n_sig_mm     ;
       double toymean_n_sb_mm      ;

       //-- Znn model 2
       double toymean_n_sigsb_ee   ;
       double toymean_n_sigsb_mm   ;

       TTree* tt_bgonly ;
       TTree* tt_splusb ;


       TRandom2* trandom_cls ;

   } ;




#endif
