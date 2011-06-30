#ifndef ra2bRoostatsClass2_h
#define ra2bRoostatsClass2_h

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


   class ra2bRoostatsClass2 {

     public :

       ra2bRoostatsClass2( bool ArgUseSigBgVars=false ) ;

       virtual ~ra2bRoostatsClass2();

       bool initialize( const char* infile ) ;
       bool reinitialize( ) ;

//     bool readTextDataset( const char* inputTextFile ) ;

//     bool generateToyDatasetsFromLikelihood( const char* outputRootFile, int nToys ) ;

//     bool generateToyDatasetsFromNVals( const char* outputRootFile, int nToys ) ;

//     bool readToyDataset( const char* inputRootFile, int dsIndex ) ;

       bool doFit() ;

//     bool doToyStudy( const char* inputRootFile, const char* outputRootFile, int dsFirst, int nToys ) ;

       bool susyScanNoContam( const char* inputScanFile, double dataLumi ) ;
       bool susyScanWithContam( const char* inputScanFile, double dataLumi ) ;
       bool setSusyScanPoint( const char* inputScanFile, double dataLumi, double m0, double m12 ) ;


//     bool sbPlotsUniformBins( const char* plotBaseName ) ;
//     bool sbPlotsVariableBins( const char* plotBaseName ) ;

       bool profileSusySig( float& susySigLow, float& susySigHigh, bool makePlot=true ) ;
       bool profileTtbarSig( float& ttbarSigLow, float& ttbarSigHigh ) ;
       bool profileQcdSig( float& qcdSigLow, float& qcdSigHigh ) ;

       void parameterSnapshot() ;

       bool fitQualityPlot( bool doNorm=false, double hmax=1.5 ) ;

       void setAndFixSusySig( double setVal = 0. ) ;
       void freeSusySig() ;


     private :

       char initializeFile[10000] ;

       bool varsAtFitVals ;
       bool useSigBgVars ;
       bool initialized ;

       double toy_mu0_ttbar_sig ;
       double toy_mu0_qcd_sig ;
       double toy_mu0_ttbar_sb ;
       double toy_mu0_qcd_sb ;
       double toy_mu0_susy_sig ;
       double toy_mu0_allbg_sig ;


      //======= Observables ==============================================================

       //-- data counts in signal region, A, and D, SB, SIG-SL, SB-SL.

       RooRealVar* rv_Nsig ;
       RooRealVar* rv_Na ;
       RooRealVar* rv_Nd ;
       RooRealVar* rv_Nsb ;
       RooRealVar* rv_Nslsig ;
       RooRealVar* rv_Nslsb ;


       //-- QCD MC counts

       RooRealVar* rv_Nqcdmca ;
       RooRealVar* rv_Nqcdmcd ;
       RooRealVar* rv_Nqcdmcsig ;
       RooRealVar* rv_Nqcdmcsb ;


       //-- EW MC counts

       RooRealVar* rv_NWJmcsig ;
       RooRealVar* rv_NWJmca ;
       RooRealVar* rv_NWJmcd ;
       RooRealVar* rv_NWJmcsb ;
       RooRealVar* rv_NWJmcslsig ;
       RooRealVar* rv_NWJmcslsb ;

       RooRealVar* rv_NZnnmcsig ;
       RooRealVar* rv_NZnnmca ;
       RooRealVar* rv_NZnnmcd ;
       RooRealVar* rv_NZnnmcsb ;
       RooRealVar* rv_NZnnmcslsig ;
       RooRealVar* rv_NZnnmcslsb ;

       RooRealVar* rv_NEwomcsig ;
       RooRealVar* rv_NEwomca ;
       RooRealVar* rv_NEwomcd ;
       RooRealVar* rv_NEwomcsb ;
       RooRealVar* rv_NEwomcslsig ;
       RooRealVar* rv_NEwomcslsb ;






       //========= Parameters ==============================================================




       //-- Counts in SIG, signal selection.

       RooAbsArg*     rv_mu_ttbar_sig ;
       RooRealVar*    rrv_mu_ttbar_sig ;
       RooFormulaVar* rfv_mu_ttbar_sig ;

       RooRealVar*    rv_mu_qcd_sig ;
       RooFormulaVar* rv_mu_ew_sig ;
       RooRealVar*    rv_mu_susymc_sig ;
       RooRealVar*    rv_mu_susy_sig ;


       //-- Counts in SB, signal selection.

       RooAbsArg*     rv_mu_ttbar_sb ;
       RooRealVar*    rrv_mu_ttbar_sb ;
       RooFormulaVar* rfv_mu_ttbar_sb ;

       RooRealVar*    rv_mu_qcd_sb ;
       RooFormulaVar* rv_mu_ew_sb ;
       RooRealVar*    rv_mu_susymc_sb ;
       RooFormulaVar* rv_mu_susy_sb ;


       //-- Counts in A, signal selection

       RooRealVar*    rv_mu_ttbar_a ;
       RooFormulaVar* rv_mu_qcd_a ;
       RooFormulaVar* rv_mu_ew_a ;
       RooRealVar*    rv_mu_susymc_a ;
       RooFormulaVar* rv_mu_susy_a ;


       //-- Counts in D, signal selection

       RooRealVar*    rv_mu_ttbar_d ;
       RooFormulaVar* rv_mu_qcd_d ;
       RooFormulaVar* rv_mu_ew_d ;
       RooRealVar*    rv_mu_susymc_d ;
       RooFormulaVar* rv_mu_susy_d ;


       //-- Single Lepton (SL) counts, SIG region

       RooRealVar*    rv_mu_sl_ttbar_sig ;
       //-- ignoring qcd in SL, SIG
       RooFormulaVar* rv_mu_sl_ew_sig ;
       RooRealVar*    rv_mu_sl_susymc_sig ;
       RooFormulaVar* rv_mu_sl_susy_sig ;


       //-- Single Lepton (SL) counts, SIG region

       RooRealVar*    rv_mu_sl_ttbar_sb ;
       //-- ignoring qcd in SL, SIG
       RooFormulaVar* rv_mu_sl_ew_sb ;
       RooRealVar*    rv_mu_sl_susymc_sb ;
       RooFormulaVar* rv_mu_sl_susy_sb ;






       //-- QCD MC, SIG,SB,A,D counts, signal selection

       RooRealVar* rv_mu_qcdmc_sig ;
       RooRealVar* rv_mu_qcdmc_sb ;
       RooRealVar* rv_mu_qcdmc_a ;
       RooRealVar* rv_mu_qcdmc_d ;



       RooRealVar* rv_mu_wjmc_sig ;
       RooRealVar* rv_mu_wjmc_a ;
       RooRealVar* rv_mu_wjmc_d ;
       RooRealVar* rv_mu_wjmc_sb ;
       RooRealVar* rv_mu_wjmc_slsig ;
       RooRealVar* rv_mu_wjmc_slsb ;

       RooRealVar* rv_mu_znnmc_sig ;
       RooRealVar* rv_mu_znnmc_a ;
       RooRealVar* rv_mu_znnmc_d ;
       RooRealVar* rv_mu_znnmc_sb ;
       RooRealVar* rv_mu_znnmc_slsig ;
       RooRealVar* rv_mu_znnmc_slsb ;

       RooRealVar* rv_mu_ewomc_sig ;
       RooRealVar* rv_mu_ewomc_a ;
       RooRealVar* rv_mu_ewomc_d ;
       RooRealVar* rv_mu_ewomc_sb ;
       RooRealVar* rv_mu_ewomc_slsig ;
       RooRealVar* rv_mu_ewomc_slsb ;







       //-- Efficiency scale factor

       RooRealVar* rv_eff_sf ;

       RooRealVar* rv_lsf_wjmc ;
       RooRealVar* rv_lsf_znnmc ;
       RooRealVar* rv_lsf_ewomc ;







       //========= Expected counts for observables in terms of parameters ============================


       RooFormulaVar* rv_n_sig ;
       RooFormulaVar* rv_n_sb ;
       RooFormulaVar* rv_n_a ;
       RooFormulaVar* rv_n_d ;
       RooFormulaVar* rv_n_sl_sig ;
       RooFormulaVar* rv_n_sl_sb ;

       RooFormulaVar* rv_n_wjmc_sig ;
       RooFormulaVar* rv_n_wjmc_a ;
       RooFormulaVar* rv_n_wjmc_d ;
       RooFormulaVar* rv_n_wjmc_sb ;
       RooFormulaVar* rv_n_wjmc_sl_sig ;
       RooFormulaVar* rv_n_wjmc_sl_sb ;

       RooFormulaVar* rv_n_znnmc_sig ;
       RooFormulaVar* rv_n_znnmc_a ;
       RooFormulaVar* rv_n_znnmc_d ;
       RooFormulaVar* rv_n_znnmc_sb ;
       RooFormulaVar* rv_n_znnmc_sl_sig ;
       RooFormulaVar* rv_n_znnmc_sl_sb ;

       RooFormulaVar* rv_n_ewomc_sig ;
       RooFormulaVar* rv_n_ewomc_a ;
       RooFormulaVar* rv_n_ewomc_d ;
       RooFormulaVar* rv_n_ewomc_sb ;
       RooFormulaVar* rv_n_ewomc_sl_sig ;
       RooFormulaVar* rv_n_ewomc_sl_sb ;



       //=========== PDFs for the likelihood ============================================================


       RooPoisson* pdf_Nsig ;
       RooPoisson* pdf_Na ;
       RooPoisson* pdf_Nd ;
       RooPoisson* pdf_Nsb ;
       RooPoisson* pdf_Nsl_sig ;
       RooPoisson* pdf_Nsl_sb ;

       RooPoisson* pdf_NWJmcsig ;
       RooPoisson* pdf_NWJmca ;
       RooPoisson* pdf_NWJmcd ;
       RooPoisson* pdf_NWJmcsb ;
       RooPoisson* pdf_NWJmcsl_sig ;
       RooPoisson* pdf_NWJmcsl_sb ;

       RooPoisson* pdf_NZnnmcsig ;
       RooPoisson* pdf_NZnnmca ;
       RooPoisson* pdf_NZnnmcd ;
       RooPoisson* pdf_NZnnmcsb ;
       RooPoisson* pdf_NZnnmcsl_sig ;
       RooPoisson* pdf_NZnnmcsl_sb ;

       RooPoisson* pdf_NEwomcsig ;
       RooPoisson* pdf_NEwomca ;
       RooPoisson* pdf_NEwomcd ;
       RooPoisson* pdf_NEwomcsb ;
       RooPoisson* pdf_NEwomcsl_sig ;
       RooPoisson* pdf_NEwomcsl_sb ;

       RooGaussian* pdf_Nqcdmc_sig ;
       RooGaussian* pdf_Nqcdmc_sb ;
       RooGaussian* pdf_Nqcdmc_a ;
       RooGaussian* pdf_Nqcdmc_d ;

       RooGaussian* pdf_Eff_sf ;
       RooGaussian* pdf_lsf_WJmc ;
       RooGaussian* pdf_lsf_Znnmc ;
       RooGaussian* pdf_lsf_Ewomc ;

       RooProdPdf* likelihood ;


       //============ Other things =================================================================


       float Nqcdmcsigerr, Nqcdmcsberr, Nqcdmcaerr, Nqcdmcderr ; //-- QCD MC uncertainties in SIG, SB, A, and D.
       float EffScaleFactor, EffScaleFactorErr ;

       float lsf_WJmc, lsf_WJmc_err ;
       float lsf_Znnmc, lsf_Znnmc_err ;
       float lsf_Ewomc, lsf_Ewomc_err ;

       RooArgSet observedParametersList ;
       RooDataSet* dsObserved ;

       RooWorkspace* workspace ;

       RooFitResult* fitResult ;

       float qcdCorrection ;
       float qcdCorrectionErr ;

   } ;




#endif
