#ifndef ra2bRoostatsClass3D_3b_h
#define ra2bRoostatsClass3D_3b_h

#include "TRandom2.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooStats/ModelConfig.h"


   class ra2bRoostatsClass3D_3b {

     public :

       ra2bRoostatsClass3D_3b() ;

       virtual ~ra2bRoostatsClass3D_3b();

       bool initialize( const char* infile = "inFile.txt",
                        const char* inputScanFile = "scanFile.txt",
                        double m0 = 1225., double m12 = 225., bool isT1bbbb = false, double t1bbbbXsec=0.,
                        const char* inputSusy_deff_dbtageff_file = "deff_dbtag_file.txt",
			const char* inputSusy_deff_dbtageff_lightflavor_file = "deff_dbtag_file.txt",
                        int   qcdModelIndex = 4,
                        const char* wsrootfilename = "ws.root",
                        const char* blindBinsList = "null",
			bool constrainBjetShape = false,
			bool floatSLSigRatios = false,
			const char* systFile1 = "systFile1.txt",
			const char* pdf_syst_file = "foo",
			const char* isr_syst_file = "foo2",
			const char* wjets_xsec_shapesyst_file = "blah.txt",
			const char* singletop_xsec_shapesyst_file = "blah.txt"
                        ) ;
       bool setSusyScanPoint( const char* inputScanFile,
                              double m0, double m12
                            ) ;

       void mismatchErr(char* label, TString inPar);
       void SetConstants(RooWorkspace pWs, RooStats::ModelConfig pMc);
       void SetConstant(const RooArgSet * vars, Bool_t value );

     private :

       bool setupShapeSyst( const char* infile, const char* systName,
                            int constraintType, // 1=gaussian, 2=...
                            double target_mgl, double target_mlsp,
                            RooWorkspace& workspace
                             ) ;

       RooAbsReal* makeBetaPrimeConstraint( const char* NP_name, double NP_val, double NP_err ) ;
       RooAbsReal* makeBetaConstraint( const char* NP_name, double NP_val, double NP_err, RooWorkspace& workspace ) ;
       RooAbsReal* makeCorrelatedBetaPrimeConstraint( const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign=false ) ;
       RooAbsReal* makeCorrelatedBetaConstraint( const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign=false ) ;

       RooAbsReal* makeGaussianConstraint( const char* NP_name, double NP_val, double NP_err, bool allowNegative = false ) ;
       RooAbsReal* makeLognormalConstraint( const char* NP_name, double NP_val, double NP_err ) ;
       RooAbsReal* makeCorrelatedGaussianConstraint( const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign=false, bool allowNegative = false ) ;
       RooAbsReal* makeCorrelatedLognormalConstraint( const char* NP_name, double NP_val, double NP_err, const char* NP_base_name, bool changeSign=false ) ;

       char initializeFile[10000] ;

       bool varsAtFitVals ;
       bool initialized ;


       RooFitResult* fitResult ;


       RooArgSet* globalObservables ;
       RooArgSet* allNuisances ;
       RooArgSet* allNuisancePdfs ;
       RooArgSet* observedParametersList ;



       // number of bins of the analysis

       static const int nBinsMET  = 3 ;
       static const int nBinsHT   = 3 ;

       static const int nBinsBtag = 2 ;
       
     
       // luminosity

       /// static const float DataLumi(9999.) ;     // integrated luminosity (in pb-1)

       bool ignoreBin[nBinsMET][nBinsHT] ;

       // observables

       RooRealVar* rv_0lep[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_1lepSig[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_1lep[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_Zee[nBinsMET][nBinsHT] ;
       RooRealVar* rv_Zmm[nBinsMET][nBinsHT] ;


       // likelihood parameters

       RooRealVar* rrv_mu_ttwj[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rrv_mu_ttwj_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rrv_mu_ttwj_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*  rv_mu_ttwj[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*  rv_mu_ttwj_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*  rv_mu_ttwj_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*  rv_mu_ttwj_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rfv_mu_ttwj[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rfv_mu_ttwj_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rfv_mu_ttwj_sl[nBinsMET][nBinsHT][nBinsBtag] ;

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

       ///// RooRealVar*  rv_mu_susy_M1_H1_1b ;
       RooRealVar*  rv_mu_susy_all0lep ;
       RooAbsArg*   rv_mu_susy[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*   rv_mu_susy_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*   rv_mu_susy_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsArg*   rv_mu_susy_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       RooFormulaVar* rv_mu_vv[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_mu_vv_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_mu_vv_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_mu_vv_ldp[nBinsMET][nBinsHT][nBinsBtag] ;


       // MC inputs

       RooRealVar* rv_mu_susymc[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_susymc_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_susymc_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_susymc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       RooRealVar* rv_mu_ttbarsingletopzjetsmc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_WJmc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_Znnmc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       RooRealVar* rv_mu_vvmc[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_vvmc_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_vvmc_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mu_vvmc_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       
       // gaussian constraints

       RooRealVar* rv_mean_eff_sf[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_width_eff_sf[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mean_eff_sf_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_width_eff_sf_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mean_eff_sf_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_width_eff_sf_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_mean_eff_sf_ldp[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_width_eff_sf_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

  //   RooFormulaVar* rv_eff_sf[nBinsMET][nBinsHT][nBinsBtag] ;
  //   RooFormulaVar* rv_eff_sf_sl[nBinsMET][nBinsHT][nBinsBtag] ;
  //   RooFormulaVar* rv_eff_sf_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       RooRealVar* rv_mean_vv_sf;
       RooRealVar* rv_width_vv_sf;
       
       // btag efficiency derivatives

       RooRealVar* rv_deff_dbtageff[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_deff_dbtageff_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooRealVar* rv_deff_dbtageff_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       RooFormulaVar* rv_btageff_sf[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_btageff_sf_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_btageff_sf_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       
       // expected event counts

       RooFormulaVar* rv_n[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_n_slSig[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_n_sl[nBinsMET][nBinsHT][nBinsBtag] ;
       RooFormulaVar* rv_n_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       RooFormulaVar* rv_n_ee[nBinsMET][nBinsHT] ;
       RooFormulaVar* rv_n_mm[nBinsMET][nBinsHT] ;


       // pdf's
       
       ///// RooPoisson* pdf_N_0lep[nBinsMET][nBinsHT][nBinsBtag] ;
       RooAbsReal* pdf_N_0lep[nBinsMET][nBinsHT][nBinsBtag] ;
       RooPoisson* pdf_N_1lepSig[nBinsMET][nBinsHT][nBinsBtag] ;
       RooPoisson* pdf_N_1lep[nBinsMET][nBinsHT][nBinsBtag] ;
       RooPoisson* pdf_N_ldp[nBinsMET][nBinsHT][nBinsBtag] ;

       RooPoisson* pdf_N_Zee[nBinsMET][nBinsHT] ;
       RooPoisson* pdf_N_Zmm[nBinsMET][nBinsHT] ;

       RooProdPdf* likelihood;


   } ;




#endif
