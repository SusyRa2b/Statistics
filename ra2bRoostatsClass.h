#ifndef ra2bRoostatsClass_h
#define ra2bRoostatsClass_h

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


   class ra2bRoostatsClass {

     public :

       ra2bRoostatsClass( bool ArgUseSigBgVars=true ) ;

       virtual ~ra2bRoostatsClass();

       bool initialize( const char* infile ) ;
       bool reinitialize( ) ;

       bool readTextDataset( const char* inputTextFile ) ;

       bool generateToyDatasetsFromLikelihood( const char* outputRootFile, int nToys ) ;

       bool generateToyDatasetsFromNVals( const char* outputRootFile, int nToys ) ;

       bool readToyDataset( const char* inputRootFile, int dsIndex ) ;

       bool doFit() ;

       bool doToyStudy( const char* inputRootFile, const char* outputRootFile, int dsFirst, int nToys ) ;

       bool susyScanNoContam( const char* inputScanFile, double dataLumi ) ;


       bool sbPlotsUniformBins( const char* plotBaseName ) ;
       bool sbPlotsVariableBins( const char* plotBaseName ) ;

       bool profileSusySig( float& susySigLow, float& susySigHigh, bool makePlot=true ) ;
       bool profileTtbarSig( float& ttbarSigLow, float& ttbarSigHigh ) ;
       bool profileQcdSig( float& qcdSigLow, float& qcdSigHigh ) ;
       bool profileTtbarSb( float& ttbarSbLow, float& ttbarSbHigh ) ;
       bool profileQcdSb( float& qcdSbLow, float& qcdSbHigh ) ;





     private :

       char initializeFile[10000] ;

       bool varsAtFitVals ;
       bool useSigBgVars ;
       bool initialized ;


      //======= Observables ==============================================================

       //-- data counts in signal region, A, and D.

       RooRealVar* rv_Nsig ;
       RooRealVar* rv_Na ;
       RooRealVar* rv_Nd ;

       double toy_mu0_ttbar_sig ;
       double toy_mu0_qcd_sig ;
       double toy_mu0_ttbar_sb ;
       double toy_mu0_qcd_sb ;
       double toy_mu0_susy_sig ;
       double toy_mu0_allbg_sig ;


       //-- data counts in 5 3-jet mass bins of SB.

       RooRealVar* rv_Nsb1 ;
       RooRealVar* rv_Nsb2 ;
       RooRealVar* rv_Nsb3 ;
       RooRealVar* rv_Nsb4 ;
       RooRealVar* rv_Nsb5 ;


       //-- data counts in 5 3-jet mass bins of LSB.

       RooRealVar* rv_Nlsb1 ;
       RooRealVar* rv_Nlsb2 ;
       RooRealVar* rv_Nlsb3 ;
       RooRealVar* rv_Nlsb4 ;
       RooRealVar* rv_Nlsb5 ;


       //-- data counts in 5 3-jet mass bins of SL, SIG.

       RooRealVar* rv_Nslsig1 ;
       RooRealVar* rv_Nslsig2 ;
       RooRealVar* rv_Nslsig3 ;
       RooRealVar* rv_Nslsig4 ;
       RooRealVar* rv_Nslsig5 ;


       //-- data counts in 5 3-jet mass bins of SL, SB.

       RooRealVar* rv_Nslsb1 ;
       RooRealVar* rv_Nslsb2 ;
       RooRealVar* rv_Nslsb3 ;
       RooRealVar* rv_Nslsb4 ;
       RooRealVar* rv_Nslsb5 ;


       //-- data counts in 5 3-jet mass bins of SL, MSB.

       RooRealVar* rv_Nslmsb1 ;
       RooRealVar* rv_Nslmsb2 ;
       RooRealVar* rv_Nslmsb3 ;
       RooRealVar* rv_Nslmsb4 ;
       RooRealVar* rv_Nslmsb5 ;



       //-- QCD MC counts

       RooRealVar* rv_Nqcdmca ;
       RooRealVar* rv_Nqcdmcd ;
       RooRealVar* rv_Nqcdmcsig ;
       RooRealVar* rv_Nqcdmcsb ;




       //========= Parameters ==============================================================

       RooAbsArg* rv_mu_ttbar_sig ;
       RooAbsArg* rv_mu_qcd_sig ;

       RooRealVar* rrv_mu_ttbar_sig ;    // copies of above for conveinence (to avoid casting)
       RooRealVar* rrv_mu_qcd_sig ;      // copies of above for conveinence (to avoid casting)
       RooFormulaVar* rfv_mu_ttbar_sig ; // copies of above for conveinence (to avoid casting)
       RooFormulaVar* rfv_mu_qcd_sig ;   // copies of above for conveinence (to avoid casting)

       RooAbsArg* rv_mu_ttbar_sb ;
       RooAbsArg* rv_mu_qcd_sb ;

       RooRealVar* rrv_mu_ttbar_sb ;    // copies of above for conveinence (to avoid casting)
       RooRealVar* rrv_mu_qcd_sb ;      // copies of above for conveinence (to avoid casting)
       RooFormulaVar* rfv_mu_ttbar_sb ; // copies of above for conveinence (to avoid casting)
       RooFormulaVar* rfv_mu_qcd_sb ;   // copies of above for conveinence (to avoid casting)


       //-- Counts in SIG, signal selection.

       RooRealVar* rv_mu_ew_sig ;
       RooRealVar* rv_mu_susy_sig ;

       RooRealVar* rv_eff_sf ;


       //-- Counts in SB, bins of 3-jet mass, EW, signal selection.

       RooRealVar* rv_mu_ew_sb1 ;
       RooRealVar* rv_mu_ew_sb2 ;
       RooRealVar* rv_mu_ew_sb3 ;
       RooRealVar* rv_mu_ew_sb4 ;
       RooRealVar* rv_mu_ew_sb5 ;


       //-- Counts in SB, bins of 3-jet mass, SUSY, signal selection.

       RooRealVar* rv_mu_susy_sb1 ;
       RooRealVar* rv_mu_susy_sb2 ;
       RooRealVar* rv_mu_susy_sb3 ;
       RooRealVar* rv_mu_susy_sb4 ;
       RooRealVar* rv_mu_susy_sb5 ;


       //-- Counts in A, signal selection

       RooRealVar* rv_mu_ttbar_a ;
       RooRealVar* rv_mu_qcd_a ;
       RooRealVar* rv_mu_ew_a ;
       RooRealVar* rv_mu_susy_a ;


       //-- Counts in D, signal selection

       RooRealVar* rv_mu_ttbar_d ;
       RooRealVar* rv_mu_qcd_d ;
       RooRealVar* rv_mu_ew_d ;
       RooRealVar* rv_mu_susy_d ;


       //-- Counts in LSB, bins of 3-jet mass, qcd, signal selection

       RooRealVar* rv_mu_qcd_lsb1 ;
       RooRealVar* rv_mu_qcd_lsb2 ;
       RooRealVar* rv_mu_qcd_lsb3 ;
       RooRealVar* rv_mu_qcd_lsb4 ;
       RooRealVar* rv_mu_qcd_lsb5 ;


       //-- QCD MC, SIG,SB,A,D counts, signal selection

       RooRealVar* rv_mu_qcdmc_sig ;
       RooRealVar* rv_mu_qcdmc_sb ;
       RooRealVar* rv_mu_qcdmc_a ;
       RooRealVar* rv_mu_qcdmc_d ;


       //-- Single Lepton (SL) counts, ttbar, SIG region, bins of 3-jet mass

       RooRealVar* rv_mu_sl_ttbar_sig1 ;
       RooRealVar* rv_mu_sl_ttbar_sig2 ;
       RooRealVar* rv_mu_sl_ttbar_sig3 ;
       RooRealVar* rv_mu_sl_ttbar_sig4 ;
       RooRealVar* rv_mu_sl_ttbar_sig5 ;


       //-- Single Lepton (SL) counts, ttbar, SB region, bins of 3-jet mass

       RooRealVar* rv_mu_sl_ttbar_sb1 ;
       RooRealVar* rv_mu_sl_ttbar_sb2 ;
       RooRealVar* rv_mu_sl_ttbar_sb3 ;
       RooRealVar* rv_mu_sl_ttbar_sb4 ;
       RooRealVar* rv_mu_sl_ttbar_sb5 ;


       //-- Single Lepton (SL) counts, ttbar, MSB region, bins of 3-jet mass

       RooRealVar* rv_mu_sl_ttbar_msb1 ;
       RooRealVar* rv_mu_sl_ttbar_msb2 ;
       RooRealVar* rv_mu_sl_ttbar_msb3 ;
       RooRealVar* rv_mu_sl_ttbar_msb4 ;
       RooRealVar* rv_mu_sl_ttbar_msb5 ;



       //-- Single Lepton (SL) counts, ew, SIG region, bins of 3-jet mass

       RooRealVar* rv_mu_sl_ew_sig1 ;
       RooRealVar* rv_mu_sl_ew_sig2 ;
       RooRealVar* rv_mu_sl_ew_sig3 ;
       RooRealVar* rv_mu_sl_ew_sig4 ;
       RooRealVar* rv_mu_sl_ew_sig5 ;


       //-- Single Lepton (SL) counts, ew, SB region, bins of 3-jet mass

       RooRealVar* rv_mu_sl_ew_sb1 ;
       RooRealVar* rv_mu_sl_ew_sb2 ;
       RooRealVar* rv_mu_sl_ew_sb3 ;
       RooRealVar* rv_mu_sl_ew_sb4 ;
       RooRealVar* rv_mu_sl_ew_sb5 ;


       //-- Single Lepton (SL) counts, ew, MSB region, bins of 3-jet mass

       RooRealVar* rv_mu_sl_ew_msb1 ;
       RooRealVar* rv_mu_sl_ew_msb2 ;
       RooRealVar* rv_mu_sl_ew_msb3 ;
       RooRealVar* rv_mu_sl_ew_msb4 ;
       RooRealVar* rv_mu_sl_ew_msb5 ;


       //-- Single Lepton (SL) counts, susy, SIG region, bins of 3-jet mass

       RooRealVar* rv_mu_sl_susy_sig1 ;
       RooRealVar* rv_mu_sl_susy_sig2 ;
       RooRealVar* rv_mu_sl_susy_sig3 ;
       RooRealVar* rv_mu_sl_susy_sig4 ;
       RooRealVar* rv_mu_sl_susy_sig5 ;


       //-- Single Lepton (SL) counts, susy, SB region, bins of 3-jet mass

       RooRealVar* rv_mu_sl_susy_sb1 ;
       RooRealVar* rv_mu_sl_susy_sb2 ;
       RooRealVar* rv_mu_sl_susy_sb3 ;
       RooRealVar* rv_mu_sl_susy_sb4 ;
       RooRealVar* rv_mu_sl_susy_sb5 ;


       //-- Single Lepton (SL) counts, susy, MSB region, bins of 3-jet mass

       RooRealVar* rv_mu_sl_susy_msb1 ;
       RooRealVar* rv_mu_sl_susy_msb2 ;
       RooRealVar* rv_mu_sl_susy_msb3 ;
       RooRealVar* rv_mu_sl_susy_msb4 ;
       RooRealVar* rv_mu_sl_susy_msb5 ;



      //========= Relationships between parameters ==================================================


       RooFormulaVar* rv_mu_sl_ttbar_sb ;
       RooFormulaVar* rv_mu_sl_ttbar_sig ;

       RooFormulaVar* rv_mu_qcd_lsb ;

       RooFormulaVar* rv_f_qcd_lsb1 ;
       RooFormulaVar* rv_f_qcd_lsb2 ;
       RooFormulaVar* rv_f_qcd_lsb3 ;
       RooFormulaVar* rv_f_qcd_lsb4 ;
       RooFormulaVar* rv_f_qcd_lsb5 ;

       RooFormulaVar* rv_mu_qcd_sb1 ;
       RooFormulaVar* rv_mu_qcd_sb2 ;
       RooFormulaVar* rv_mu_qcd_sb3 ;
       RooFormulaVar* rv_mu_qcd_sb4 ;
       RooFormulaVar* rv_mu_qcd_sb5 ;

       RooFormulaVar* rv_mu_sl_ttbar1 ;
       RooFormulaVar* rv_mu_sl_ttbar2 ;
       RooFormulaVar* rv_mu_sl_ttbar3 ;
       RooFormulaVar* rv_mu_sl_ttbar4 ;
       RooFormulaVar* rv_mu_sl_ttbar5 ;

       RooFormulaVar* rv_mu_sl_ttbar ;

       RooFormulaVar* rv_f_sl_ttbar1 ;
       RooFormulaVar* rv_f_sl_ttbar2 ;
       RooFormulaVar* rv_f_sl_ttbar3 ;
       RooFormulaVar* rv_f_sl_ttbar4 ;
       RooFormulaVar* rv_f_sl_ttbar5 ;

       RooFormulaVar* rv_mu_ttbar_sb1 ;
       RooFormulaVar* rv_mu_ttbar_sb2 ;
       RooFormulaVar* rv_mu_ttbar_sb3 ;
       RooFormulaVar* rv_mu_ttbar_sb4 ;
       RooFormulaVar* rv_mu_ttbar_sb5 ;


       //========= Expected counts for observables in terms of parameters ============================


       RooFormulaVar* rv_n_sig ;

       RooFormulaVar* rv_n_sb1 ;
       RooFormulaVar* rv_n_sb2 ;
       RooFormulaVar* rv_n_sb3 ;
       RooFormulaVar* rv_n_sb4 ;
       RooFormulaVar* rv_n_sb5 ;

       RooFormulaVar* rv_n_a ;
       RooFormulaVar* rv_n_d ;

       RooFormulaVar* rv_n_lsb1 ;
       RooFormulaVar* rv_n_lsb2 ;
       RooFormulaVar* rv_n_lsb3 ;
       RooFormulaVar* rv_n_lsb4 ;
       RooFormulaVar* rv_n_lsb5 ;

       RooFormulaVar* rv_n_sl_sig1 ;
       RooFormulaVar* rv_n_sl_sig2 ;
       RooFormulaVar* rv_n_sl_sig3 ;
       RooFormulaVar* rv_n_sl_sig4 ;
       RooFormulaVar* rv_n_sl_sig5 ;

       RooFormulaVar* rv_n_sl_sb1 ;
       RooFormulaVar* rv_n_sl_sb2 ;
       RooFormulaVar* rv_n_sl_sb3 ;
       RooFormulaVar* rv_n_sl_sb4 ;
       RooFormulaVar* rv_n_sl_sb5 ;

       RooFormulaVar* rv_n_sl_msb1 ;
       RooFormulaVar* rv_n_sl_msb2 ;
       RooFormulaVar* rv_n_sl_msb3 ;
       RooFormulaVar* rv_n_sl_msb4 ;
       RooFormulaVar* rv_n_sl_msb5 ;



       //=========== PDFs for the likelihood ============================================================


       RooPoisson* pdf_Nsig ;

       RooPoisson* pdf_Na ;
       RooPoisson* pdf_Nd ;

       RooPoisson* pdf_Nsb1 ;
       RooPoisson* pdf_Nsb2 ;
       RooPoisson* pdf_Nsb3 ;
       RooPoisson* pdf_Nsb4 ;
       RooPoisson* pdf_Nsb5 ;

       RooPoisson* pdf_Nlsb1 ;
       RooPoisson* pdf_Nlsb2 ;
       RooPoisson* pdf_Nlsb3 ;
       RooPoisson* pdf_Nlsb4 ;
       RooPoisson* pdf_Nlsb5 ;

       RooPoisson* pdf_Nsl_sig1 ;
       RooPoisson* pdf_Nsl_sig2 ;
       RooPoisson* pdf_Nsl_sig3 ;
       RooPoisson* pdf_Nsl_sig4 ;
       RooPoisson* pdf_Nsl_sig5 ;

       RooPoisson* pdf_Nsl_sb1 ;
       RooPoisson* pdf_Nsl_sb2 ;
       RooPoisson* pdf_Nsl_sb3 ;
       RooPoisson* pdf_Nsl_sb4 ;
       RooPoisson* pdf_Nsl_sb5 ;

       RooPoisson* pdf_Nsl_msb1 ;
       RooPoisson* pdf_Nsl_msb2 ;
       RooPoisson* pdf_Nsl_msb3 ;
       RooPoisson* pdf_Nsl_msb4 ;
       RooPoisson* pdf_Nsl_msb5 ;

       RooGaussian* pdf_Nqcdmc_sig ;
       RooGaussian* pdf_Nqcdmc_sb ;
       RooGaussian* pdf_Nqcdmc_a ;
       RooGaussian* pdf_Nqcdmc_d ;

       RooGaussian* pdf_Eff_sf ;

       RooProdPdf* likelihood ;


       //============ Other things =================================================================


       RooArgSet observedParametersList ;
       RooDataSet* dsObserved ;

       RooWorkspace* workspace ;

       RooFitResult* fitResult ;

       float qcdCorrection ;
       float qcdCorrectionErr ;

   } ;




#endif
