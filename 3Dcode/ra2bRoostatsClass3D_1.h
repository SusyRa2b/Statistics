#ifndef ra2bRoostatsClass3D_1_h
#define ra2bRoostatsClass3D_1_h

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


   class ra2bRoostatsClass3D_1 {

     public :

       ra2bRoostatsClass3D_1() ;

       virtual ~ra2bRoostatsClass3D_1();

       bool initialize( const char* infile = "inFile.txt",
                        const char* inputScanFile = "scanFile.txt",
                        double m0 = 875., double m12 = 525., bool isT1bbbb = false, double t1bbbbXsec=0.,
                        const char* inputSusy_deff_dbtageff_file = "deff_dbtag_file.txt"
                        ) ;
       bool setSusyScanPoint( const char* inputScanFile,
                              double m0, double m12, bool isT1bbbb = false, double t1bbbbXsec=0.,
                              const char* inputSusy_deff_dbtageff_file = "deff_dbtag_file.txt"
                            ) ;

       void mismatchErr(char* label, TString inPar);

     private :

       char initializeFile[10000] ;

       bool varsAtFitVals ;
       bool initialized ;

       RooArgSet observedParametersList ;
       RooDataSet* dsObserved ;


       RooFitResult* fitResult ;


       // number of bins of the analysis

       static const int nBinsMET  = 3 ;
       static const int nBinsHT   = 3 ;
       static const int nBinsBtag = 3 ;    // this must always be 3


       // luminosity

       static const float DataLumi = 9999. ;     // integrated luminosity (in pb-1)


       // observables

       RooRealVar* rv_0lep[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_1lep[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_Zee[nBinsMET][nBinsHT] ;
       RooRealVar* rv_Zmm[nBinsMET][nBinsHT] ;


       // likelihood parameters

       RooRealVar* rrv_mu_ttwj[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*  rv_mu_ttwj[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_ttwj_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*  rv_mu_ttwj_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rfv_mu_ttwj[nBinsMET][nBinsHT][nBinsBtag] ;

       RooRealVar* rrv_mu_qcd[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*  rv_mu_qcd[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rrv_mu_qcd_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_qcd_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rfv_mu_qcd[nBinsMET][nBinsHT][nBinsBtag] ;

       RooRealVar* rrv_mu_znn[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*  rv_mu_znn[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rfv_mu_znn[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_mu_znn_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_mu_zee[nBinsMET][nBinsHT] ;
       RooFormulaVar* rv_mu_zmm[nBinsMET][nBinsHT] ;

       RooRealVar*  rv_mu_susy_M1_H1_1b ;
       RooAbsArg*   rv_mu_susy[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*   rv_mu_susy_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*   rv_mu_susy_ldp[nBinsMET][nBinsHT][nBinsBtag] ;


       // MC inputs

       RooRealVar* rv_mu_susymc[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_susymc_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_susymc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       RooRealVar* rv_mu_ttbarsingletopzjetsmc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_WJmc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_Znnmc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       
       // gaussian constraints

       RooRealVar* rv_mean_eff_sf[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_width_eff_sf[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mean_eff_sf_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_width_eff_sf_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mean_eff_sf_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_width_eff_sf_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       RooFormulaVar* rv_eff_sf[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_eff_sf_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_eff_sf_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       
       // btag efficiency derivatives

       RooRealVar* rv_deff_dbtageff[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_deff_dbtageff_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_deff_dbtageff_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       RooFormulaVar* rv_btageff_sf[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_btageff_sf_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_btageff_sf_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       
       // expected event counts

       RooFormulaVar* rv_n[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_n_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_n_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       RooFormulaVar* rv_n_ee[nBinsMET][nBinsHT] ;
       RooFormulaVar* rv_n_mm[nBinsMET][nBinsHT] ;


       // pdf's
       
       RooPoisson* pdf_N_0lep[nBinsMET][nBinsHT][nBinsBtag] ;
       RooPoisson* pdf_N_1lep[nBinsMET][nBinsHT][nBinsBtag] ;
       RooPoisson* pdf_N_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       RooPoisson* pdf_N_Zee[nBinsMET][nBinsHT] ;
       RooPoisson* pdf_N_Zmm[nBinsMET][nBinsHT] ;

       RooProdPdf* likelihood;


   } ;




#endif
