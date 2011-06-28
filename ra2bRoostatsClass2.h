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

       ra2bRoostatsClass2( bool ArgUseSigBgVars=true ) ;

       virtual ~ra2bRoostatsClass2();

       bool initialize( const char* infile ) ;
//     bool reinitialize( ) ;

//     bool readTextDataset( const char* inputTextFile ) ;

//     bool generateToyDatasetsFromLikelihood( const char* outputRootFile, int nToys ) ;

//     bool generateToyDatasetsFromNVals( const char* outputRootFile, int nToys ) ;

//     bool readToyDataset( const char* inputRootFile, int dsIndex ) ;

//     bool doFit() ;

//     bool doToyStudy( const char* inputRootFile, const char* outputRootFile, int dsFirst, int nToys ) ;

//     bool susyScanNoContam( const char* inputScanFile, double dataLumi ) ;
//     bool susyScanWithContam( const char* inputScanFile, double dataLumi ) ;
//     bool setSusyScanPoint( const char* inputScanFile, double dataLumi, double m0, double m12 ) ;


//     bool sbPlotsUniformBins( const char* plotBaseName ) ;
//     bool sbPlotsVariableBins( const char* plotBaseName ) ;

//     bool profileSusySig( float& susySigLow, float& susySigHigh, bool makePlot=true ) ;
//     bool profileTtbarSig( float& ttbarSigLow, float& ttbarSigHigh ) ;
//     bool profileQcdSig( float& qcdSigLow, float& qcdSigHigh ) ;
//     bool profileTtbarSb( float& ttbarSbLow, float& ttbarSbHigh ) ;
//     bool profileQcdSb( float& qcdSbLow, float& qcdSbHigh ) ;

//     void parameterSnapshot() ;

//     bool fitQualityPlot( bool doNorm=false, double hmax=1.5 ) ;

//     void setAndFixSusySig( double setVal = 0. ) ;
//     void freeSusySig() ;


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




       //========= Parameters ==============================================================




       //-- Counts in SIG, signal selection.

       RooRealVar* rv_mu_ttbar_sig ;
       RooRealVar* rv_mu_qcd_sig ;
       RooRealVar* rv_mu_ew_sig ;
       RooRealVar* rv_mu_susy_sig ;
       RooRealVar* rv_mu_susymc_sig ;

       RooRealVar* rv_eff_sf ;


       //-- Counts in SB, EW, signal selection.

       RooFormulaVar* rv_mu_ttbar_sb ;
       RooFormulaVar* rv_mu_qcd_sb ;
       RooRealVar* rv_mu_ew_sb ;


       //-- Counts in SB, SUSY, signal selection.

       RooRealVar* rv_mu_susymc_sb ;


       //-- Counts in A, signal selection

       RooRealVar* rv_mu_ttbar_a ;
       RooRealVar* rv_mu_qcd_a ;
       RooRealVar* rv_mu_ew_a ;
       RooRealVar* rv_mu_susymc_a ;


       //-- Counts in D, signal selection

       RooRealVar* rv_mu_ttbar_d ;
       RooRealVar* rv_mu_qcd_d ;
       RooRealVar* rv_mu_ew_d ;
       RooRealVar* rv_mu_susymc_d ;


       //-- QCD MC, SIG,SB,A,D counts, signal selection

       RooRealVar* rv_mu_qcdmc_sig ;
       RooRealVar* rv_mu_qcdmc_sb ;
       RooRealVar* rv_mu_qcdmc_a ;
       RooRealVar* rv_mu_qcdmc_d ;


       //-- Single Lepton (SL) counts, ttbar, SIG region

       RooRealVar* rv_mu_sl_ttbar_sig ;


       //-- Single Lepton (SL) counts, ttbar, SB region

       RooRealVar* rv_mu_sl_ttbar_sb ;




       //-- Single Lepton (SL) counts, ew, SIG region

       RooRealVar* rv_mu_sl_ew_sig ;


       //-- Single Lepton (SL) counts, ew, SB region

       RooRealVar* rv_mu_sl_ew_sb ;


       //-- Single Lepton (SL) counts, susy, SIG region

       RooRealVar* rv_mu_sl_susymc_sig ;


       //-- Single Lepton (SL) counts, susy, SB region

       RooRealVar* rv_mu_sl_susymc_sb ;





      //========= Relationships between parameters ==================================================


       RooFormulaVar* rv_mu_sl_ttbar_sb ;
       RooFormulaVar* rv_mu_sl_ttbar_sig ;

       RooFormulaVar* rv_mu_qcd_sb ;

       RooFormulaVar* rv_mu_sl_ttbar ;

       RooFormulaVar* rv_mu_sl_ttbar ;

       RooFormulaVar* rv_f_sl_ttbar ;

       RooFormulaVar* rv_mu_ttbar_sb ;


       RooFormulaVar* rv_mu_susy_a ;
       RooFormulaVar* rv_mu_susy_d ;

       RooFormulaVar* rv_mu_susy_sb ;

       RooFormulaVar* rv_mu_sl_susy_sig ;

       RooFormulaVar* rv_mu_sl_susy_sb ;



       //========= Expected counts for observables in terms of parameters ============================


       RooFormulaVar* rv_n_sig ;
       RooFormulaVar* rv_n_sb ;
       RooFormulaVar* rv_n_a ;
       RooFormulaVar* rv_n_d ;
       RooFormulaVar* rv_n_sl_sig ;
       RooFormulaVar* rv_n_sl_sb ;




       //=========== PDFs for the likelihood ============================================================


       RooPoisson* pdf_Nsig ;
       RooPoisson* pdf_Na ;
       RooPoisson* pdf_Nd ;
       RooPoisson* pdf_Nsb ;
       RooPoisson* pdf_Nsl_sig ;
       RooPoisson* pdf_Nsl_sb ;

       RooGaussian* pdf_Nqcdmc_sig ;
       RooGaussian* pdf_Nqcdmc_sb ;
       RooGaussian* pdf_Nqcdmc_a ;
       RooGaussian* pdf_Nqcdmc_d ;

       RooGaussian* pdf_Eff_sf ;

       RooProdPdf* likelihood ;


       //============ Other things =================================================================


       float Nqcdmcsigerr, Nqcdmcsberr, Nqcdmcaerr, Nqcdmcderr ; //-- QCD MC uncertainties in SIG, SB, A, and D.
       float EffScaleFactor, EffScaleFactorErr ;

       RooArgSet observedParametersList ;
       RooDataSet* dsObserved ;

       RooWorkspace* workspace ;

       RooFitResult* fitResult ;

       float qcdCorrection ;
       float qcdCorrectionErr ;

   } ;




#endif
