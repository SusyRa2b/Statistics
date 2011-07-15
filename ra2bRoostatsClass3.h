#ifndef ra2bRoostatsClass3_h
#define ra2bRoostatsClass3_h

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

       ra2bRoostatsClass3( bool ArgUseSigTtwjVar=false, bool ArgUseLdpVars=true ) ;

       virtual ~ra2bRoostatsClass3();

       bool initialize( const char* infile ) ;

//     bool reinitialize( ) ;

//     bool readTextDataset( const char* inputTextFile ) ;

//     bool generateToyDatasetsFromLikelihood( const char* outputRootFile, int nToys ) ;

//     bool generateToyDatasetsFromNVals( const char* outputRootFile, int nToys ) ;

//     bool readToyDataset( const char* inputRootFile, int dsIndex ) ;

       bool doFit() ;

//     bool doToyStudy( const char* inputRootFile, const char* outputRootFile, int dsFirst, int nToys ) ;

//     bool susyScanNoContam( const char* inputScanFile, double dataLumi ) ;
//     bool susyScanWithContam( const char* inputScanFile, double dataLumi ) ;
//     bool setSusyScanPoint( const char* inputScanFile, double dataLumi, double m0, double m12 ) ;


//     bool sbPlotsUniformBins( const char* plotBaseName ) ;
//     bool sbPlotsVariableBins( const char* plotBaseName ) ;

       bool profileSusySig( float& susySigLow, float& susySigHigh, bool makePlot=true ) ;

       bool profileZnnSig( float& znnSigLow, float& znnSigHigh, bool makePlot=true ) ;
       bool profileZnnSb( float& znnSbLow, float& znnSbHigh, bool makePlot=true ) ;

       bool profilettwjSig( float& ttwjSigLow, float& ttwjSigHigh, bool makePlot=true ) ;
       bool profilettwjSb( float& ttwjSbLow, float& ttwjSbHigh, bool makePlot=true ) ;

       bool profileqcdSig( float& qcdSigLow, float& qcdSigHigh, bool makePlot=true ) ;
       bool profileqcdSb( float& qcdSbLow, float& qcdSbHigh, bool makePlot=true ) ;



       void parameterSnapshot() ;

//     bool fitQualityPlot( bool doNorm=false, double hmax=1.5 ) ;

       void setAndFixSusySig( double setVal = 0. ) ;
       void freeSusySig() ;


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
       RooRealVar* rv_Nsb_ee      ;
       RooRealVar* rv_Nsig_ee     ;
       RooRealVar* rv_Nsb_mm      ;
       RooRealVar* rv_Nsig_mm     ;






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


      //-- Z to nunu

       RooRealVar*    rv_mu_znn_sig     ;
       RooRealVar*    rv_mu_znn_sb      ;
       RooFormulaVar* rv_mu_znn_sig_ldp ;
       RooFormulaVar* rv_mu_znn_sb_ldp  ;
       RooFormulaVar* rv_mu_znn_sig_ee  ;
       RooFormulaVar* rv_mu_znn_sb_ee   ;
       RooFormulaVar* rv_mu_znn_sig_mm  ;
       RooFormulaVar* rv_mu_znn_sb_mm   ;


       RooFormulaVar* rv_mu_zee_sig_ee ;
       RooFormulaVar* rv_mu_zee_sb_ee  ;
       RooFormulaVar* rv_mu_zmm_sig_mm ;
       RooFormulaVar* rv_mu_zmm_sb_mm  ;



      //-- EW0

       RooFormulaVar* rv_mu_ewo_sig ;
       RooFormulaVar* rv_mu_ewo_sb      ;
       RooFormulaVar* rv_mu_ewo_sig_ldp ;
       RooFormulaVar* rv_mu_ewo_sb_ldp  ;


      //-- Gaussian constraints and constants.

       RooRealVar* rv_eff_sf ;

       RooRealVar* rv_lsf_WJmc   ;
       RooRealVar* rv_lsf_Znnmc  ;
       RooRealVar* rv_lsf_Ewomc  ;
       RooRealVar* rv_sf_ttbarmc ;

       RooRealVar* rv_acc_ee ;
       RooRealVar* rv_acc_mm ;
       RooRealVar* rv_eff_ee ;
       RooRealVar* rv_eff_mm ;

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
       RooFormulaVar* rv_n_sig_ee     ;
       RooFormulaVar* rv_n_sb_ee      ;
       RooFormulaVar* rv_n_sig_mm     ;
       RooFormulaVar* rv_n_sb_mm      ;
       RooFormulaVar* rv_n_lsb_0b     ;
       RooFormulaVar* rv_n_lsb_0b_ldp ;


       //=========== PDFs for the likelihood ============================================================

       RooPoisson*  pdf_Nsig        ;
       RooPoisson*  pdf_Nsb         ;
       RooPoisson*  pdf_Nsig_ldp    ;
       RooPoisson*  pdf_Nsb_ldp     ;
       RooPoisson*  pdf_Nsig_sl     ;
       RooPoisson*  pdf_Nsb_sl      ;
       RooPoisson*  pdf_Nsig_ee     ;
       RooPoisson*  pdf_Nsb_ee      ;
       RooPoisson*  pdf_Nsig_mm     ;
       RooPoisson*  pdf_Nsb_mm      ;
       RooPoisson*  pdf_Nlsb_0b     ;
       RooPoisson*  pdf_Nlsb_0b_ldp ;
       RooGaussian* pdf_lsf_WJmc    ;
       RooGaussian* pdf_lsf_Znnmc   ;
       RooGaussian* pdf_lsf_Ewomc   ;
       RooGaussian* pdf_sf_ttbarmc  ;
       RooGaussian* pdf_acc_ee      ;
       RooGaussian* pdf_acc_mm      ;
       RooGaussian* pdf_eff_ee      ;
       RooGaussian* pdf_eff_mm      ;
       RooGaussian* pdf_Eff_sf      ;

       RooProdPdf*  likelihood ;



       //============ Other things =================================================================


       float EffScaleFactor, EffScaleFactorErr ;

       float lsf_WJmc, lsf_WJmc_err ;
       float lsf_Znnmc, lsf_Znnmc_err ;
       float lsf_Ewomc, lsf_Ewomc_err ;
       float sf_ttbarmc, sf_ttbarmc_err ;

       float  acc_ee_mean             ;
       float  acc_ee_err              ;
       float  acc_mm_mean             ;
       float  acc_mm_err              ;
       float  eff_ee_mean             ;
       float  eff_ee_err              ;
       float  eff_mm_mean             ;
       float  eff_mm_err              ;
       float  Ztoll_lumi              ;
       float  Ztoll_tight_sf          ;
       float  Ztoll_tight_sf_err      ;
       float  DataLumi                ;

       RooArgSet observedParametersList ;
       RooDataSet* dsObserved ;

       RooWorkspace* workspace ;

       RooFitResult* fitResult ;


   } ;




#endif
