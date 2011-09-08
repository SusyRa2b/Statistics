//
//   Owen Long, UCR
//   Harrison Prosper, FSU
//   Sezen Sekmen, FSU
//
//


#include "ra2bRoostatsClass7.h"

#include <iostream>
#include <string.h>


#include "TCanvas.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TText.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TH2F.h"
#include "TGaxis.h"
#include "TLine.h"
#include "TStringLong.h"

#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooTrace.h"
#include "RooUniform.h"

#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"

#include "LikelihoodIntervalPlot.cxx"

  using namespace RooFit ;
  using namespace RooStats ;


  //=====================================================================================================


   ra2bRoostatsClass7::ra2bRoostatsClass7( bool ArgUseSigTtwjVar, bool ArgUseLdpVars, int ArgznnModel ) {

      gStyle->SetOptStat(0) ;

      useSigTtwjVar = ArgUseSigTtwjVar ;
      useLdpVars = ArgUseLdpVars ;
      znnModel = ArgznnModel ;

      if ( znnModel > 2 ) {
         printf("\n\n\n *** Unrecognized znn model number %d.  Should be 1 or 2.\n\n\n", znnModel ) ;
      }

     //--- Tell RooFit to shut up about anything less important than an ERROR.
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
  /// RooMsgService::instance().addStream(DEBUG,Topic(Tracing),OutputFile("debug.log")) ;
      printf("\n\n ==== RooFit output configuration ===============\n") ;
      RooMsgService::instance().Print("v") ;
      printf("\n\n ================================================\n") ;

      varsAtFitVals = false ;
      initialized = false ;




      //--- CLs stuff below here
      tt_bgonly = 0x0 ;
      tt_splusb = 0x0 ;

      toymean_n_sig         = 0. ;
      toymean_n_sb          = 0. ;
      toymean_n_sig_ldp     = 0. ;
      toymean_n_sb_ldp      = 0. ;
      toymean_n_sig_sl      = 0. ;
      toymean_n_sb_sl       = 0. ;
      toymean_n_lsb_0b      = 0. ;
      toymean_n_lsb_0b_ldp  = 0. ;

       //-- Znn model 1
      toymean_n_sig_ee      = 0. ;
      toymean_n_sb_ee       = 0. ;
      toymean_n_sig_mm      = 0. ;
      toymean_n_sb_mm       = 0. ;

       //-- Znn model 2
      toymean_n_sigsb_ee    = 0. ;
      toymean_n_sigsb_mm    = 0. ;





      trandom_cls = new TRandom2(12345) ;

   }





  //===================================================================

   ra2bRoostatsClass7::~ra2bRoostatsClass7() {

      if ( !initialized ) return ;

//    delete dsObserved ;
//    //delete plinterval_susy_sig ;

//    delete rv_Nsig ;
//    delete rv_Na ;
//    delete rv_Nd ;
//    delete rv_Nsb1 ;
//    delete rv_Nsb2 ;
//    delete rv_Nsb3 ;
//    delete rv_Nsb4 ;
//    delete rv_Nsb5 ;
//    delete rv_Nlsb1 ;
//    delete rv_Nlsb2 ;
//    delete rv_Nlsb3 ;
//    delete rv_Nlsb4 ;
//    delete rv_Nlsb5 ;
//    delete rv_Nslsig1 ;
//    delete rv_Nslsig2 ;
//    delete rv_Nslsig3 ;
//    delete rv_Nslsig4 ;
//    delete rv_Nslsig5 ;
//    delete rv_Nslsb1 ;
//    delete rv_Nslsb2 ;
//    delete rv_Nslsb3 ;
//    delete rv_Nslsb4 ;
//    delete rv_Nslsb5 ;
//    delete rv_Nslmsb1 ;
//    delete rv_Nslmsb2 ;
//    delete rv_Nslmsb3 ;
//    delete rv_Nslmsb4 ;
//    delete rv_Nslmsb5 ;
//    delete rv_mu_ttbar_sig ;
//    delete rv_mu_qcd_sig ;
//    delete rv_mu_ttbar_sb ;
//    delete rv_mu_qcd_sb ;
//    delete rv_mu_ew_sig    ;
//    delete rv_mu_susy_sig  ;
//    delete rv_mu_ew_sb1 ;
//    delete rv_mu_ew_sb2 ;
//    delete rv_mu_ew_sb3 ;
//    delete rv_mu_ew_sb4 ;
//    delete rv_mu_ew_sb5 ;

//    delete rv_mu_susy_sb1 ;
//    delete rv_mu_susy_sb2 ;
//    delete rv_mu_susy_sb3 ;
//    delete rv_mu_susy_sb4 ;
//    delete rv_mu_susy_sb5 ;
//    delete rv_mu_ttbar_a ;
//    delete rv_mu_qcd_a   ;
//    delete rv_mu_ew_a    ;
//    delete rv_mu_susy_a  ;
//    delete rv_mu_ttbar_d ;
//    delete rv_mu_qcd_d   ;
//    delete rv_mu_ew_d    ;
//    delete rv_mu_susy_d  ;
//    delete rv_mu_qcd_lsb1 ;
//    delete rv_mu_qcd_lsb2 ;
//    delete rv_mu_qcd_lsb3 ;
//    delete rv_mu_qcd_lsb4 ;
//    delete rv_mu_qcd_lsb5 ;
//    delete rv_mu_qcdmc_sig ;
//    delete rv_mu_qcdmc_sb  ;
//    delete rv_mu_qcdmc_a   ;
//    delete rv_mu_qcdmc_d   ;

//    delete rv_mu_sl_ttbar_sig1 ;
//    delete rv_mu_sl_ttbar_sig2 ;
//    delete rv_mu_sl_ttbar_sig3 ;
//    delete rv_mu_sl_ttbar_sig4 ;
//    delete rv_mu_sl_ttbar_sig5 ;
//    delete rv_mu_sl_ttbar_sb1 ;
//    delete rv_mu_sl_ttbar_sb2 ;
//    delete rv_mu_sl_ttbar_sb3 ;
//    delete rv_mu_sl_ttbar_sb4 ;
//    delete rv_mu_sl_ttbar_sb5 ;
//    delete rv_mu_sl_ttbar_msb1 ;
//    delete rv_mu_sl_ttbar_msb2 ;
//    delete rv_mu_sl_ttbar_msb3 ;
//    delete rv_mu_sl_ttbar_msb4 ;
//    delete rv_mu_sl_ttbar_msb5 ;
//    delete rv_mu_sl_ew_sig1 ;
//    delete rv_mu_sl_ew_sig2 ;
//    delete rv_mu_sl_ew_sig3 ;
//    delete rv_mu_sl_ew_sig4 ;
//    delete rv_mu_sl_ew_sig5 ;
//    delete rv_mu_sl_ew_sb1 ;
//    delete rv_mu_sl_ew_sb2 ;
//    delete rv_mu_sl_ew_sb3 ;
//    delete rv_mu_sl_ew_sb4 ;
//    delete rv_mu_sl_ew_sb5 ;
//    delete rv_mu_sl_ew_msb1 ;
//    delete rv_mu_sl_ew_msb2 ;
//    delete rv_mu_sl_ew_msb3 ;
//    delete rv_mu_sl_ew_msb4 ;
//    delete rv_mu_sl_ew_msb5 ;

//    delete rv_mu_sl_susy_sig1 ;
//    delete rv_mu_sl_susy_sig2 ;
//    delete rv_mu_sl_susy_sig3;
//    delete rv_mu_sl_susy_sig4 ;
//    delete rv_mu_sl_susy_sig5 ;
//    delete rv_mu_sl_susy_sb1 ;
//    delete rv_mu_sl_susy_sb2 ;
//    delete rv_mu_sl_susy_sb3 ;
//    delete rv_mu_sl_susy_sb4 ;
//    delete rv_mu_sl_susy_sb5 ;
//    delete rv_mu_sl_susy_msb1 ;
//    delete rv_mu_sl_susy_msb2 ;
//    delete rv_mu_sl_susy_msb3 ;
//    delete rv_mu_sl_susy_msb4 ;
//    delete rv_mu_sl_susy_msb5 ;

//    delete rv_mu_sl_ttbar_sb ;
//    delete rv_mu_sl_ttbar_sig ;
//    delete rv_mu_qcd_lsb ;

//    delete rv_f_qcd_lsb1 ;
//    delete rv_f_qcd_lsb2 ;
//    delete rv_f_qcd_lsb3 ;
//    delete rv_f_qcd_lsb4 ;
//    delete rv_f_qcd_lsb5 ;
//    delete rv_mu_qcd_sb1 ;
//    delete rv_mu_qcd_sb2 ;
//    delete rv_mu_qcd_sb3 ;
//    delete rv_mu_qcd_sb4 ;
//    delete rv_mu_qcd_sb5 ;
//    delete rv_mu_sl_ttbar1 ;
//    delete rv_mu_sl_ttbar2 ;
//    delete rv_mu_sl_ttbar3 ;
//    delete rv_mu_sl_ttbar4 ;
//    delete rv_mu_sl_ttbar5 ;
//    delete rv_mu_sl_ttbar ;
//    delete rv_f_sl_ttbar1 ;
//    delete rv_f_sl_ttbar2 ;
//    delete rv_f_sl_ttbar3 ;
//    delete rv_f_sl_ttbar4 ;
//    delete rv_f_sl_ttbar5 ;
//    delete rv_mu_ttbar_sb1 ;
//    delete rv_mu_ttbar_sb2 ;
//    delete rv_mu_ttbar_sb3 ;
//    delete rv_mu_ttbar_sb4 ;
//    delete rv_mu_ttbar_sb5 ;
//    delete rv_n_sig ;
//    delete rv_n_sb1 ;
//    delete rv_n_sb2 ;
//    delete rv_n_sb3 ;
//    delete rv_n_sb4 ;
//    delete rv_n_sb5 ;
//    delete rv_n_a ;
//    delete rv_n_d ;
//    delete rv_n_lsb1 ;
//    delete rv_n_lsb2 ;
//    delete rv_n_lsb3 ;
//    delete rv_n_lsb4 ;
//    delete rv_n_lsb5 ;
//    delete rv_n_sl_sig1 ;
//    delete rv_n_sl_sig2 ;
//    delete rv_n_sl_sig3 ;
//    delete rv_n_sl_sig4 ;
//    delete rv_n_sl_sig5 ;
//    delete rv_n_sl_sb1 ;
//    delete rv_n_sl_sb2 ;
//    delete rv_n_sl_sb3 ;
//    delete rv_n_sl_sb4 ;
//    delete rv_n_sl_sb5 ;
//    delete rv_n_sl_msb1 ;
//    delete rv_n_sl_msb2 ;
//    delete rv_n_sl_msb3 ;
//    delete rv_n_sl_msb4 ;
//    delete rv_n_sl_msb5 ;

//    delete pdf_Nsig ;
//    delete pdf_Na ;
//    delete pdf_Nd ;
//    delete pdf_Nsb1 ;
//    delete pdf_Nsb2 ;
//    delete pdf_Nsb3 ;
//    delete pdf_Nsb4 ;
//    delete pdf_Nsb5 ;
//    delete pdf_Nlsb1 ;
//    delete pdf_Nlsb2 ;
//    delete pdf_Nlsb3 ;
//    delete pdf_Nlsb4 ;
//    delete pdf_Nlsb5 ;
//    delete pdf_Nsl_sig1 ;
//    delete pdf_Nsl_sig2 ;
//    delete pdf_Nsl_sig3 ;
//    delete pdf_Nsl_sig4 ;
//    delete pdf_Nsl_sig5 ;
//    delete pdf_Nsl_sb1 ;
//    delete pdf_Nsl_sb2 ;
//    delete pdf_Nsl_sb3 ;
//    delete pdf_Nsl_sb4 ;
//    delete pdf_Nsl_sb5 ;
//    delete pdf_Nsl_msb1 ;
//    delete pdf_Nsl_msb2 ;
//    delete pdf_Nsl_msb3 ;
//    delete pdf_Nsl_msb4 ;
//    delete pdf_Nsl_msb5 ;
//    delete pdf_Nqcdmc_sig ;
//    delete pdf_Nqcdmc_sb  ;
//    delete pdf_Nqcdmc_a   ;
//    delete pdf_Nqcdmc_d   ;

//    delete likelihood ;

//    delete workspace ;



   }



  //===================================================================

    bool ra2bRoostatsClass7::initialize( const char* infile ,
                                         const char* inputScanFile, double m0, double m12, bool isT1bbbb, double t1bbbbXsec ) {


       printf( "\n\n Opening input file : %s\n\n", infile ) ;

       FILE* infp ;
       if ( (infp=fopen( infile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
          return false ;
       }

       sprintf( initializeFile, "%s", infile ) ;


       int    Nsig                  ;
       int    Nsb                   ;
       int    Nsig_sl               ;
       int    Nsb_sl                ;
       int    Nsig_ldp              ;
       int    Nsb_ldp               ;
       int    Nlsb                  ;
       int    Nlsb_ldp              ;
       int    Nlsb_0b               ;
       int    Nlsb_0b_ldp           ;
       float  Nqcdmc_sig            ;
       float  Nqcdmc_sig_err        ;
       float  Nqcdmc_sb             ;
       float  Nqcdmc_sb_err         ;
       float  Nqcdmc_sig_sl         ;
       float  Nqcdmc_sig_sl_err     ;
       float  Nqcdmc_sb_sl          ;
       float  Nqcdmc_sb_sl_err      ;
       float  Nqcdmc_sig_ldp        ;
       float  Nqcdmc_sig_ldp_err    ;
       float  Nqcdmc_sb_ldp         ;
       float  Nqcdmc_sb_ldp_err     ;
       float  Nqcdmc_lsb            ;
       float  Nqcdmc_lsb_err        ;
       float  Nqcdmc_lsb_ldp        ;
       float  Nqcdmc_lsb_ldp_err    ;
       float  Nqcdmc_lsb_0b         ;
       float  Nqcdmc_lsb_0b_err     ;
       float  Nqcdmc_lsb_0b_ldp     ;
       float  Nqcdmc_lsb_0b_ldp_err ;
       float  Nttbarmc_sig          ;
       float  Nttbarmc_sb           ;
       float  Nttbarmc_sig_sl       ;
       float  Nttbarmc_sb_sl        ;
       float  Nttbarmc_sig_ldp      ;
       float  Nttbarmc_sb_ldp       ;
       float  Nttbarmc_lsb          ;
       float  Nttbarmc_lsb_ldp      ;
       float  NWJmc_sig             ;
       float  NWJmc_sb              ;
       float  NWJmc_sig_sl          ;
       float  NWJmc_sb_sl           ;
       float  NWJmc_sig_ldp         ;
       float  NWJmc_sb_ldp          ;
       float  NWJmc_lsb             ;
       float  NWJmc_lsb_ldp         ;
       float  NZnnmc_sig            ;
       float  NZnnmc_sb             ;
       float  NZnnmc_sig_sl         ;
       float  NZnnmc_sb_sl          ;
       float  NZnnmc_sig_ldp        ;
       float  NZnnmc_sb_ldp         ;
       float  NZnnmc_lsb            ;
       float  NZnnmc_lsb_ldp        ;
       float  NEwomc_sig            ;
       float  NEwomc_sb             ;
       float  NEwomc_sig_sl         ;
       float  NEwomc_sb_sl          ;
       float  NEwomc_sig_ldp        ;
       float  NEwomc_sb_ldp         ;
       float  NEwomc_lsb            ;
       float  NEwomc_lsb_ldp        ;
       float  Nsusymc_sig           ;
       float  Nsusymc_sb            ;
       float  Nsusymc_sig_sl        ;
       float  Nsusymc_sb_sl         ;
       float  Nsusymc_sig_ldp       ;
       float  Nsusymc_sb_ldp        ;
       float  Nsusymc_lsb           ;
       float  Nsusymc_lsb_ldp       ;
       float  Nsusymc_lsb_0b        ;
       float  Nsusymc_lsb_0b_ldp    ;
       int    Nhtonlytrig_lsb_0b      ;
       int    Nhtonlytrig_lsb_0b_ldp  ;
       int    Nsb_ee                  ;
       int    Nsig_ee                 ;
       int    Nsb_mm                  ;
       int    Nsig_mm                 ;

       float  sf_mc            ;
       float  sf_mc_err        ;
       float  sf_qcd_sb        ;
       float  sf_qcd_sb_err    ;
       float  sf_qcd_sig       ;
       float  sf_qcd_sig_err   ;
       float  sf_ttwj_sig      ;
       float  sf_ttwj_sig_err  ;
       float  sf_ee            ;
       float  sf_ee_err        ;
       float  sf_mm            ;
       float  sf_mm_err        ;











       //--- read in description line.
       printf("\n\n") ;
       char c(0) ;
       while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
       printf("\n\n") ;


       char label[1000] ;

      //--- Inputs generated with gen_roostats_input4.c
      //    The order here must be consistent with the order there!

       fscanf( infp, "%s %g", label, &EffScaleFactor        ) ;   printf( "%s %g\n", label, EffScaleFactor        ) ;
       fscanf( infp, "%s %g", label, &EffScaleFactorErr     ) ;   printf( "%s %g\n", label, EffScaleFactorErr     ) ;
       fscanf( infp, "%s %d", label, &Nsig                  ) ;   printf( "%s %d\n", label, Nsig                  ) ;
       fscanf( infp, "%s %d", label, &Nsb                   ) ;   printf( "%s %d\n", label, Nsb                   ) ;
       fscanf( infp, "%s %d", label, &Nsig_sl               ) ;   printf( "%s %d\n", label, Nsig_sl               ) ;
       fscanf( infp, "%s %d", label, &Nsb_sl                ) ;   printf( "%s %d\n", label, Nsb_sl                ) ;
       fscanf( infp, "%s %d", label, &Nsig_ldp              ) ;   printf( "%s %d\n", label, Nsig_ldp              ) ;
       fscanf( infp, "%s %d", label, &Nsb_ldp               ) ;   printf( "%s %d\n", label, Nsb_ldp               ) ;
       fscanf( infp, "%s %d", label, &Nlsb                  ) ;   printf( "%s %d\n", label, Nlsb                  ) ;
       fscanf( infp, "%s %d", label, &Nlsb_ldp              ) ;   printf( "%s %d\n", label, Nlsb_ldp              ) ;
       fscanf( infp, "%s %d", label, &Nlsb_0b               ) ;   printf( "%s %d\n", label, Nlsb_0b               ) ;
       fscanf( infp, "%s %d", label, &Nlsb_0b_ldp           ) ;   printf( "%s %d\n", label, Nlsb_0b_ldp           ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig            ) ;   printf( "%s %g\n", label, Nqcdmc_sig            ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_err        ) ;   printf( "%s %g\n", label, Nqcdmc_sig_err        ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb             ) ;   printf( "%s %g\n", label, Nqcdmc_sb             ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_err         ) ;   printf( "%s %g\n", label, Nqcdmc_sb_err         ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_sl         ) ;   printf( "%s %g\n", label, Nqcdmc_sig_sl         ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_sl_err     ) ;   printf( "%s %g\n", label, Nqcdmc_sig_sl_err     ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_sl          ) ;   printf( "%s %g\n", label, Nqcdmc_sb_sl          ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_sl_err      ) ;   printf( "%s %g\n", label, Nqcdmc_sb_sl_err      ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_ldp        ) ;   printf( "%s %g\n", label, Nqcdmc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sig_ldp_err    ) ;   printf( "%s %g\n", label, Nqcdmc_sig_ldp_err    ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_ldp         ) ;   printf( "%s %g\n", label, Nqcdmc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_sb_ldp_err     ) ;   printf( "%s %g\n", label, Nqcdmc_sb_ldp_err     ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb            ) ;   printf( "%s %g\n", label, Nqcdmc_lsb            ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_err        ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_err        ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_ldp        ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_ldp        ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_ldp_err    ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_ldp_err    ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_0b         ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_0b         ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_0b_err     ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_0b_err     ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_0b_ldp     ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_0b_ldp     ) ;
       fscanf( infp, "%s %g", label, &Nqcdmc_lsb_0b_ldp_err ) ;   printf( "%s %g\n", label, Nqcdmc_lsb_0b_ldp_err ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sig          ) ;   printf( "%s %g\n", label, Nttbarmc_sig          ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sb           ) ;   printf( "%s %g\n", label, Nttbarmc_sb           ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sig_sl       ) ;   printf( "%s %g\n", label, Nttbarmc_sig_sl       ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sb_sl        ) ;   printf( "%s %g\n", label, Nttbarmc_sb_sl        ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sig_ldp      ) ;   printf( "%s %g\n", label, Nttbarmc_sig_ldp      ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_sb_ldp       ) ;   printf( "%s %g\n", label, Nttbarmc_sb_ldp       ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_lsb          ) ;   printf( "%s %g\n", label, Nttbarmc_lsb          ) ;
       fscanf( infp, "%s %g", label, &Nttbarmc_lsb_ldp      ) ;   printf( "%s %g\n", label, Nttbarmc_lsb_ldp      ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig             ) ;   printf( "%s %g\n", label, NWJmc_sig             ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb              ) ;   printf( "%s %g\n", label, NWJmc_sb              ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig_sl          ) ;   printf( "%s %g\n", label, NWJmc_sig_sl          ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb_sl           ) ;   printf( "%s %g\n", label, NWJmc_sb_sl           ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig_ldp         ) ;   printf( "%s %g\n", label, NWJmc_sig_ldp         ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb_ldp          ) ;   printf( "%s %g\n", label, NWJmc_sb_ldp          ) ;
       fscanf( infp, "%s %g", label, &NWJmc_lsb             ) ;   printf( "%s %g\n", label, NWJmc_lsb             ) ;
       fscanf( infp, "%s %g", label, &NWJmc_lsb_ldp         ) ;   printf( "%s %g\n", label, NWJmc_lsb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig            ) ;   printf( "%s %g\n", label, NZnnmc_sig            ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb             ) ;   printf( "%s %g\n", label, NZnnmc_sb             ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig_sl         ) ;   printf( "%s %g\n", label, NZnnmc_sig_sl         ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb_sl          ) ;   printf( "%s %g\n", label, NZnnmc_sb_sl          ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig_ldp        ) ;   printf( "%s %g\n", label, NZnnmc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb_ldp         ) ;   printf( "%s %g\n", label, NZnnmc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_lsb            ) ;   printf( "%s %g\n", label, NZnnmc_lsb            ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_lsb_ldp        ) ;   printf( "%s %g\n", label, NZnnmc_lsb_ldp        ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig            ) ;   printf( "%s %g\n", label, NEwomc_sig            ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb             ) ;   printf( "%s %g\n", label, NEwomc_sb             ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig_sl         ) ;   printf( "%s %g\n", label, NEwomc_sig_sl         ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb_sl          ) ;   printf( "%s %g\n", label, NEwomc_sb_sl          ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sig_ldp        ) ;   printf( "%s %g\n", label, NEwomc_sig_ldp        ) ;
       fscanf( infp, "%s %g", label, &NEwomc_sb_ldp         ) ;   printf( "%s %g\n", label, NEwomc_sb_ldp         ) ;
       fscanf( infp, "%s %g", label, &NEwomc_lsb            ) ;   printf( "%s %g\n", label, NEwomc_lsb            ) ;
       fscanf( infp, "%s %g", label, &NEwomc_lsb_ldp        ) ;   printf( "%s %g\n", label, NEwomc_lsb_ldp        ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sig           ) ;   printf( "%s %g\n", label, Nsusymc_sig           ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sb            ) ;   printf( "%s %g\n", label, Nsusymc_sb            ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sig_sl        ) ;   printf( "%s %g\n", label, Nsusymc_sig_sl        ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sb_sl         ) ;   printf( "%s %g\n", label, Nsusymc_sb_sl         ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sig_ldp       ) ;   printf( "%s %g\n", label, Nsusymc_sig_ldp       ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_sb_ldp        ) ;   printf( "%s %g\n", label, Nsusymc_sb_ldp        ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_lsb           ) ;   printf( "%s %g\n", label, Nsusymc_lsb           ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_lsb_ldp       ) ;   printf( "%s %g\n", label, Nsusymc_lsb_ldp       ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_lsb_0b        ) ;   printf( "%s %g\n", label, Nsusymc_lsb_0b        ) ;
       fscanf( infp, "%s %g", label, &Nsusymc_lsb_0b_ldp    ) ;   printf( "%s %g\n", label, Nsusymc_lsb_0b_ldp    ) ;
       fscanf( infp, "%s %d", label, &Nhtonlytrig_lsb_0b        ) ;   printf( "%s %d\n", label, Nhtonlytrig_lsb_0b        ) ;
       fscanf( infp, "%s %d", label, &Nhtonlytrig_lsb_0b_ldp    ) ;   printf( "%s %d\n", label, Nhtonlytrig_lsb_0b_ldp    ) ;
       fscanf( infp, "%s %g", label, &DataLumi                  ) ;   printf( "%s %g\n", label, DataLumi                  ) ;
       fscanf( infp, "%s %d", label, &Nsb_ee                    ) ;   printf( "%s %d\n", label, Nsb_ee                    ) ;
       fscanf( infp, "%s %d", label, &Nsig_ee                   ) ;   printf( "%s %d\n", label, Nsig_ee                   ) ;
       fscanf( infp, "%s %d", label, &Nsb_mm                    ) ;   printf( "%s %d\n", label, Nsb_mm                    ) ;
       fscanf( infp, "%s %d", label, &Nsig_mm                   ) ;   printf( "%s %d\n", label, Nsig_mm                   ) ;
       fscanf( infp, "%s %g", label, &acc_ee_mean               ) ;   printf( "%s %g\n", label, acc_ee_mean               ) ;
       fscanf( infp, "%s %g", label, &acc_ee_err                ) ;   printf( "%s %g\n", label, acc_ee_err                ) ;
       fscanf( infp, "%s %g", label, &acc_mm_mean               ) ;   printf( "%s %g\n", label, acc_mm_mean               ) ;
       fscanf( infp, "%s %g", label, &acc_mm_err                ) ;   printf( "%s %g\n", label, acc_mm_err                ) ;
       fscanf( infp, "%s %g", label, &eff_ee_mean               ) ;   printf( "%s %g\n", label, eff_ee_mean               ) ;
       fscanf( infp, "%s %g", label, &eff_ee_err                ) ;   printf( "%s %g\n", label, eff_ee_err                ) ;
       fscanf( infp, "%s %g", label, &eff_mm_mean               ) ;   printf( "%s %g\n", label, eff_mm_mean               ) ;
       fscanf( infp, "%s %g", label, &eff_mm_err                ) ;   printf( "%s %g\n", label, eff_mm_err                ) ;
       fscanf( infp, "%s %g", label, &Ztoll_lumi                ) ;   printf( "%s %g\n", label, Ztoll_lumi                ) ;
       fscanf( infp, "%s %g", label, &knn_sig_mean              ) ;   printf( "%s %g\n", label, knn_sig_mean              ) ;
       fscanf( infp, "%s %g", label, &knn_sig_err               ) ;   printf( "%s %g\n", label, knn_sig_err               ) ;
       fscanf( infp, "%s %g", label, &knn_sb_mean               ) ;   printf( "%s %g\n", label, knn_sb_mean               ) ;
       fscanf( infp, "%s %g", label, &knn_sb_err                ) ;   printf( "%s %g\n", label, knn_sb_err                ) ;

       fscanf( infp, "%s %g", label, &fsig_ee_mean              ) ;   printf( "fsig_ee_mean %s %g\n", label, fsig_ee_mean              ) ;
       fscanf( infp, "%s %g", label, &fsig_ee_err               ) ;   printf( "fsig_ee_err %s %g\n", label, fsig_ee_err               ) ;
       fscanf( infp, "%s %g", label, &fsig_mm_mean              ) ;   printf( "fsig_mm_mean %s %g\n", label, fsig_mm_mean              ) ;
       fscanf( infp, "%s %g", label, &fsig_mm_err               ) ;   printf( "fsig_mm_err %s %g\n", label, fsig_mm_err               ) ;

       fscanf( infp, "%s %g", label, &sf_mc                     ) ;   printf( "sf_mc %s %g\n", label, sf_mc                     ) ;
       fscanf( infp, "%s %g", label, &sf_mc_err                 ) ;   printf( "sf_mc_err %s %g\n", label, sf_mc_err                 ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb                 ) ;   printf( "sf_qcd_sb %s %g\n", label, sf_qcd_sb                 ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb_err             ) ;   printf( "sf_qcd_sb_err %s %g\n", label, sf_qcd_sb_err             ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig                ) ;   printf( "sf_qcd_sig%s %g\n", label, sf_qcd_sig                ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig_err            ) ;   printf( "sf_qcd_sig_err%s %g\n", label, sf_qcd_sig_err            ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig               ) ;   printf( "sf_ttwj_sig%s %g\n", label, sf_ttwj_sig               ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig_err           ) ;   printf( "sf_ttwj_sig_err%s %g\n", label, sf_ttwj_sig_err           ) ;
       fscanf( infp, "%s %g", label, &sf_ee                     ) ;   printf( "sf_ee%s %g\n", label, sf_ee                     ) ;
       fscanf( infp, "%s %g", label, &sf_ee_err                 ) ;   printf( "sf_ee_err%s %g\n", label, sf_ee_err                 ) ;
       fscanf( infp, "%s %g", label, &sf_mm                     ) ;   printf( "sf_mm%s %g\n", label, sf_mm                     ) ;
       fscanf( infp, "%s %g", label, &sf_mm_err                 ) ;   printf( "sf_mm_err%s %g\n", label, sf_mm_err                 ) ;



       printf("\n Done reading in %s\n\n", infile ) ;
       fclose( infp ) ;


       //+++++ Owen, Aug 6.  set all lsf_ factors to 1.  Input MC counts are now all weighted.
       //                    Calculations of stat errors below are now wrong.  Fix later.
       lsf_WJmc = 1.0 ;
       lsf_Znnmc = 1.0 ;
       lsf_Ewomc = 0.0 ; // don't use ewo.
       sf_ttbarmc = 1.0 ;

       //--- Print out a nice summary of the inputs.

       float Nsm_sig         = Nttbarmc_sig          +  lsf_WJmc*NWJmc_sig          +  Nqcdmc_sig          +  lsf_Znnmc*NZnnmc_sig          +  lsf_Ewomc*NEwomc_sig         ;
       float Nsm_sb          = Nttbarmc_sb           +  lsf_WJmc*NWJmc_sb           +  Nqcdmc_sb           +  lsf_Znnmc*NZnnmc_sb           +  lsf_Ewomc*NEwomc_sb          ;
       float Nsm_sig_sl      = Nttbarmc_sig_sl       +  lsf_WJmc*NWJmc_sig_sl       +  Nqcdmc_sig_sl       +  lsf_Znnmc*NZnnmc_sig_sl       +  lsf_Ewomc*NEwomc_sig_sl      ;
       float Nsm_sb_sl       = Nttbarmc_sb_sl        +  lsf_WJmc*NWJmc_sb_sl        +  Nqcdmc_sb_sl        +  lsf_Znnmc*NZnnmc_sb_sl        +  lsf_Ewomc*NEwomc_sb_sl       ;
       float Nsm_sig_ldp     = Nttbarmc_sig_ldp      +  lsf_WJmc*NWJmc_sig_ldp      +  Nqcdmc_sig_ldp      +  lsf_Znnmc*NZnnmc_sig_ldp      +  lsf_Ewomc*NEwomc_sig_ldp     ;
       float Nsm_sb_ldp      = Nttbarmc_sb_ldp       +  lsf_WJmc*NWJmc_sb_ldp       +  Nqcdmc_sb_ldp       +  lsf_Znnmc*NZnnmc_sb_ldp       +  lsf_Ewomc*NEwomc_sb_ldp      ;
       float Nsm_lsb         = Nttbarmc_lsb          +  lsf_WJmc*NWJmc_lsb          +  Nqcdmc_lsb          +  lsf_Znnmc*NZnnmc_lsb          +  lsf_Ewomc*NEwomc_lsb         ;
       float Nsm_lsb_ldp     = Nttbarmc_lsb_ldp      +  lsf_WJmc*NWJmc_lsb_ldp      +  Nqcdmc_lsb_ldp      +  lsf_Znnmc*NZnnmc_lsb_ldp      +  lsf_Ewomc*NEwomc_lsb_ldp     ;
       float Nsm_lsb_0b      =                                                         Nqcdmc_lsb_0b                                                                        ;
       float Nsm_lsb_0b_ldp  =                                                         Nqcdmc_lsb_0b_ldp                                                                    ;



       //--- QCD min Delta phi N  pass / fail ratios.

       float Rqcd_sig        =  0. ;
       if ( Nqcdmc_sig_ldp > 0. ) { Rqcd_sig        =  Nqcdmc_sig        / Nqcdmc_sig_ldp ;  }
       float Rqcd_sig_err    =  0. ;
       if ( Nqcdmc_sig > 0. && Nqcdmc_sig_ldp > 0. ) { Rqcd_sig_err  =  Rqcd_sig * sqrt( pow(Nqcdmc_sig_err/Nqcdmc_sig,2) + pow(Nqcdmc_sig_ldp_err/Nqcdmc_sig_ldp,2) ) ; }


       float Rqcd_sb        =  0. ;
       if ( Nqcdmc_sb_ldp > 0. ) { Rqcd_sb        =  Nqcdmc_sb        / Nqcdmc_sb_ldp ;  }
       float Rqcd_sb_err    =  0. ;
       if ( Nqcdmc_sb > 0. && Nqcdmc_sb_ldp > 0. ) { Rqcd_sb_err  =  Rqcd_sb * sqrt( pow(Nqcdmc_sb_err/Nqcdmc_sb,2) + pow(Nqcdmc_sb_ldp_err/Nqcdmc_sb_ldp,2) ) ; }


       float Rqcd_lsb        =  0. ;
       if ( Nqcdmc_lsb_ldp > 0. ) { Rqcd_lsb        =  Nqcdmc_lsb        / Nqcdmc_lsb_ldp ;  }
       float Rqcd_lsb_err    =  0. ;
       if ( Nqcdmc_lsb > 0. && Nqcdmc_lsb_ldp > 0. ) { Rqcd_lsb_err  =  Rqcd_lsb * sqrt( pow(Nqcdmc_lsb_err/Nqcdmc_lsb,2) + pow(Nqcdmc_lsb_ldp_err/Nqcdmc_lsb_ldp,2) ) ; }


       float Nlsb_corrected     = Nlsb     - Nttbarmc_lsb     - lsf_WJmc*NWJmc_lsb     - lsf_Znnmc*NZnnmc_lsb     - lsf_Ewomc*NEwomc_lsb ;
       float Nlsb_ldp_corrected = Nlsb_ldp - Nttbarmc_lsb_ldp - lsf_WJmc*NWJmc_lsb_ldp - lsf_Znnmc*NZnnmc_lsb_ldp - lsf_Ewomc*NEwomc_lsb_ldp ;
       float Rdata_lsb = 0. ;
       if ( Nlsb_ldp_corrected > 0 ) { Rdata_lsb = (Nlsb_corrected) / (Nlsb_ldp_corrected) ; }
       float Rdata_lsb_err = 0. ;
       if ( Nlsb_corrected > 0 && Nlsb_ldp_corrected > 0 ) { Rdata_lsb_err = Rdata_lsb * sqrt( 1.0/(Nlsb_corrected) + 1.0/(Nlsb_ldp_corrected) ) ; }

       float RdataHTOT_lsb_0b = 0. ;
       if ( Nhtonlytrig_lsb_0b_ldp > 0 ) { RdataHTOT_lsb_0b = (1.0*Nhtonlytrig_lsb_0b) / (1.0*Nhtonlytrig_lsb_0b_ldp) ; }
       float RdataHTOT_lsb_0b_err = 0. ;
       if ( Nhtonlytrig_lsb_0b_ldp > 0 && Nhtonlytrig_lsb_0b > 0 ) { RdataHTOT_lsb_0b_err = RdataHTOT_lsb_0b * sqrt( 1.0/(1.0*Nhtonlytrig_lsb_0b) + 1.0/(1.0*Nhtonlytrig_lsb_0b_ldp) ) ; }


       float Rqcd_lsb_0b        =  0. ;
       if ( Nqcdmc_lsb_0b_ldp > 0. ) { Rqcd_lsb_0b        =  Nqcdmc_lsb_0b        / Nqcdmc_lsb_0b_ldp ;  }
       float Rqcd_lsb_0b_err    =  0. ;
       if ( Nqcdmc_lsb_0b > 0. && Nqcdmc_lsb_0b_ldp > 0. ) { Rqcd_lsb_0b_err  =  Rqcd_lsb_0b * sqrt( pow(Nqcdmc_lsb_0b_err/Nqcdmc_lsb_0b,2) + pow(Nqcdmc_lsb_0b_ldp_err/Nqcdmc_lsb_0b_ldp,2) ) ; }


       float Rdata_lsb_0b = 0. ;
       if ( Nlsb_0b_ldp > 0 ) { Rdata_lsb_0b = (1.0*Nlsb_0b) / (1.0*Nlsb_0b_ldp) ; }
       float Rdata_lsb_0b_err = 0. ;
       if ( Nlsb_0b > 0 && Nlsb_0b_ldp > 0 ) { Rdata_lsb_0b_err = Rdata_lsb_0b * sqrt( 1.0/(1.0*Nlsb_0b) + 1.0/(1.0*Nlsb_0b_ldp) ) ; }


       //--- ttbar + Wjets  MET  SIG / SB ratios

       float Rttwj = 0. ;
       if ( (Nttbarmc_sb+lsf_WJmc*NWJmc_sb) > 0. ) { Rttwj    =   (Nttbarmc_sig+lsf_WJmc*NWJmc_sig) / (Nttbarmc_sb+lsf_WJmc*NWJmc_sb) ; }
       float Rttwj_err2 = pow( lsf_WJmc/ (Nttbarmc_sb+lsf_WJmc*NWJmc_sb),2)*NWJmc_sig
                        + pow( (Nttbarmc_sig+lsf_WJmc*NWJmc_sig) * lsf_WJmc / pow((Nttbarmc_sb+lsf_WJmc*NWJmc_sb),2), 2)*NWJmc_sb ;
       float Rttwj_err = sqrt(Rttwj_err2) ;

       float Rttwj_sl = 0. ;
       if ( (Nttbarmc_sb_sl+lsf_WJmc*NWJmc_sb_sl) > 0. ) { Rttwj_sl    =   (Nttbarmc_sig_sl+lsf_WJmc*NWJmc_sig_sl) / (Nttbarmc_sb_sl+lsf_WJmc*NWJmc_sb_sl) ; }
       float Rttwj_sl_err2 = pow( lsf_WJmc / (Nttbarmc_sb_sl+lsf_WJmc*NWJmc_sb_sl), 2) * NWJmc_sig_sl
                           + pow( lsf_WJmc * (Nttbarmc_sig_sl+lsf_WJmc*NWJmc_sig_sl) / pow((Nttbarmc_sb_sl+lsf_WJmc*NWJmc_sb_sl),2), 2)*NWJmc_sb_sl ;
       float Rttwj_sl_err = sqrt(Rttwj_sl_err2) ;

       float Rdata_sl = 0. ;
       if ( Nsb_sl > 0 ) { Rdata_sl = (1.0*Nsig_sl)/(1.0*Nsb_sl) ; }
       float Rdata_sl_err = 0. ;
       if ( Nsig_sl>0 && Nsb_sl>0 ) { Rdata_sl_err = Rdata_sl * sqrt( 1.0/(1.0*Nsig_sl) + 1.0/(1.0*Nsb_sl) ) ; }


       //--- Simple MC closure tests.

       float comp_mc_qcd_sb = Nqcdmc_sb_ldp * ( Nqcdmc_lsb_0b / Nqcdmc_lsb_0b_ldp ) ;
       float comp_mc_qcd_sb_err = Nqcdmc_sb_ldp_err * ( Nqcdmc_lsb_0b / Nqcdmc_lsb_0b_ldp ) ;

       float comp_mc_qcd_sig = Nqcdmc_sig_ldp * ( Nqcdmc_lsb_0b / Nqcdmc_lsb_0b_ldp ) ;
       float comp_mc_qcd_sig_err = Nqcdmc_sig_ldp_err * ( Nqcdmc_lsb_0b / Nqcdmc_lsb_0b_ldp ) ;

       float comp_mc_ttwj_sig = (Nttbarmc_sb + lsf_WJmc*NWJmc_sb) * ( (Nttbarmc_sig_sl + lsf_WJmc*NWJmc_sig_sl) / (Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl) ) ;

       //--- below ignores ttbar errors, uses sqrt(N) on raw WJ counts.
       float comp_mc_ttwj_sig_err2 = pow( lsf_WJmc*( (Nttbarmc_sig_sl + lsf_WJmc*NWJmc_sig_sl) / (Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl) ),2)*NWJmc_sb
                                   + pow( (Nttbarmc_sb + lsf_WJmc*NWJmc_sb) *(lsf_WJmc/ (Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl) ), 2)*NWJmc_sig_sl
                                   + pow( (Nttbarmc_sb + lsf_WJmc*NWJmc_sb) *  (Nttbarmc_sig_sl + lsf_WJmc*NWJmc_sig_sl) / pow((Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl) ,2)*lsf_WJmc, 2 )*NWJmc_sb_sl ;
       float comp_mc_ttwj_sig_err = sqrt(comp_mc_ttwj_sig_err2) ;

       //--- Simple data calculations.


       float Nsb_ldp_corrected  = Nsb_ldp  - (Nttbarmc_sb_ldp  + lsf_WJmc*NWJmc_sb_ldp  + lsf_Znnmc*NZnnmc_sb_ldp  + lsf_Ewomc*NEwomc_sb_ldp  ) ;
       float comp_data_qcd_sb = Nsb_ldp_corrected * ( (1.0*Nlsb_0b)/(1.0*Nlsb_0b_ldp) ) ;
       float comp_data_qcd_sb_err = sqrt( Nsb_ldp + Nttbarmc_sb_ldp  + lsf_WJmc*NWJmc_sb_ldp  + lsf_Znnmc*NZnnmc_sb_ldp  + lsf_Ewomc*NEwomc_sb_ldp ) * ( (1.0*Nlsb_0b)/(1.0*Nlsb_0b_ldp) ) ;

       float Nsig_ldp_corrected = Nsig_ldp - (Nttbarmc_sig_ldp + lsf_WJmc*NWJmc_sig_ldp + lsf_Znnmc*NZnnmc_sig_ldp + lsf_Ewomc*NEwomc_sig_ldp ) ;
       float comp_data_qcd_sig = Nsig_ldp_corrected * ( (1.0*Nlsb_0b)/(1.0*Nlsb_0b_ldp) ) ;
       float comp_data_qcd_sig_err = sqrt( Nsig_ldp + Nttbarmc_sig_ldp  + lsf_WJmc*NWJmc_sig_ldp  + lsf_Znnmc*NZnnmc_sig_ldp  + lsf_Ewomc*NEwomc_sig_ldp ) * ( (1.0*Nlsb_0b)/(1.0*Nlsb_0b_ldp) ) ;

       float comp_data_ttwj_sig = (Nsb - (comp_data_qcd_sb + lsf_Znnmc*NZnnmc_sb + lsf_Ewomc*NEwomc_sb)) * ( (1.0*Nsig_sl) / (1.0*Nsb_sl) ) ;
       float comp_data_ttwj_sig_err2 = pow( ( (1.0*Nsig_sl) / (1.0*Nsb_sl) ), 2)*Nsb
                                     + pow( ( (1.0*Nsig_sl) / (1.0*Nsb_sl) ), 2)*pow(comp_data_qcd_sb_err,2)
                                     + pow( (Nsb - (comp_data_qcd_sb + lsf_Znnmc*NZnnmc_sb + lsf_Ewomc*NEwomc_sb)) / (1.0*Nsb_sl), 2 )*Nsig_sl
                                     + pow( (Nsb - (comp_data_qcd_sb + lsf_Znnmc*NZnnmc_sb + lsf_Ewomc*NEwomc_sb)) * (1.0*Nsig_sl) / pow(1.0*Nsb_sl,2), 2 ) * Nsb_sl ;
       float comp_data_ttwj_sig_err = sqrt(comp_data_ttwj_sig_err2) ;


       float comp_znn_sig_ee(2.) ;
       float comp_znn_sig_mm(2.) ;
       float comp_znn_sig(2.) ;
       float comp_znn_sb_ee(2.) ;
       float comp_znn_sb_mm(2.) ;
       float comp_znn_sb(2.) ;

       float comp_znn_sig_ee_err(0.) ;
       float comp_znn_sig_mm_err(0.) ;
       float comp_znn_sig_err(0.) ;
       float comp_znn_sb_ee_err(0.) ;
       float comp_znn_sb_mm_err(0.) ;
       float comp_znn_sb_err(0.) ;


       if ( znnModel == 1 ) {

         comp_znn_sig_ee = Nsig_ee * ( 5.95 * DataLumi * fsig_ee_mean ) / ( acc_ee_mean * eff_ee_mean * Ztoll_lumi ) ;
         comp_znn_sb_ee  = Nsb_ee  * ( 5.95 * DataLumi * fsig_ee_mean ) / ( acc_ee_mean * eff_ee_mean * Ztoll_lumi ) ;

         comp_znn_sig_mm = Nsig_mm * ( 5.95 * DataLumi * fsig_mm_mean ) / ( acc_mm_mean * eff_mm_mean * Ztoll_lumi ) ;
         comp_znn_sb_mm  = Nsb_mm  * ( 5.95 * DataLumi * fsig_mm_mean ) / ( acc_mm_mean * eff_mm_mean * Ztoll_lumi ) ;

         if ( Nsig_ee > 0 ) { comp_znn_sig_ee_err = comp_znn_sig_ee * sqrt( 1.0/(1.0*Nsig_ee) + pow(fsig_ee_err/fsig_ee_mean,2) + pow(acc_ee_err/acc_ee_mean,2) + pow(eff_ee_err/eff_ee_mean,2) ) ; }
         if ( Nsb_ee  > 0 ) { comp_znn_sb_ee_err  = comp_znn_sb_ee  * sqrt( 1.0/(1.0*Nsb_ee)  + pow(fsig_ee_err/fsig_ee_mean,2) + pow(acc_ee_err/acc_ee_mean,2) + pow(eff_ee_err/eff_ee_mean,2) ) ; }

         if ( Nsig_mm > 0 ) { comp_znn_sig_mm_err = comp_znn_sig_mm * sqrt( 1.0/(1.0*Nsig_mm) + pow(fsig_mm_err/fsig_mm_mean,2) + pow(acc_mm_err/acc_mm_mean,2) + pow(eff_mm_err/eff_mm_mean,2) ) ; }
         if ( Nsb_mm  > 0 ) { comp_znn_sb_mm_err  = comp_znn_sb_mm  * sqrt( 1.0/(1.0*Nsb_mm)  + pow(fsig_mm_err/fsig_mm_mean,2) + pow(acc_mm_err/acc_mm_mean,2) + pow(eff_mm_err/eff_mm_mean,2) ) ; }



       } else if ( znnModel == 2 ) {

         comp_znn_sig_ee = (Nsig_ee + Nsb_ee) * ( 5.95 * DataLumi * fsig_ee_mean * knn_sig_mean ) / ( acc_ee_mean * eff_ee_mean * Ztoll_lumi ) ;
         comp_znn_sb_ee  = (Nsig_ee + Nsb_ee) * ( 5.95 * DataLumi * fsig_ee_mean * knn_sb_mean  ) / ( acc_ee_mean * eff_ee_mean * Ztoll_lumi ) ;

         comp_znn_sig_mm = (Nsig_mm + Nsb_mm) * ( 5.95 * DataLumi * fsig_mm_mean * knn_sig_mean ) / ( acc_mm_mean * eff_mm_mean * Ztoll_lumi ) ;
         comp_znn_sb_mm  = (Nsig_mm + Nsb_mm) * ( 5.95 * DataLumi * fsig_mm_mean * knn_sb_mean  ) / ( acc_mm_mean * eff_mm_mean * Ztoll_lumi ) ;

         if ( (Nsig_ee+Nsb_ee) > 0 ) { comp_znn_sig_ee_err = comp_znn_sig_ee * sqrt( 1.0/(1.0*(Nsig_ee+Nsb_ee)) + pow(fsig_ee_err/fsig_ee_mean,2) + pow(acc_ee_err/acc_ee_mean,2) + pow(eff_ee_err/eff_ee_mean,2) + pow(knn_sig_err/knn_sig_mean,2) ) ; }
         if ( (Nsig_ee+Nsb_ee) > 0 ) { comp_znn_sb_ee_err  = comp_znn_sb_ee  * sqrt( 1.0/(1.0*(Nsig_ee+Nsb_ee)) + pow(fsig_ee_err/fsig_ee_mean,2) + pow(acc_ee_err/acc_ee_mean,2) + pow(eff_ee_err/eff_ee_mean,2) + pow(knn_sig_err/knn_sig_mean,2) ) ; }

         if ( (Nsig_mm+Nsb_mm) > 0 ) { comp_znn_sig_mm_err = comp_znn_sig_mm * sqrt( 1.0/(1.0*(Nsig_mm+Nsb_mm)) + pow(fsig_mm_err/fsig_mm_mean,2) + pow(acc_mm_err/acc_mm_mean,2) + pow(eff_mm_err/eff_mm_mean,2) + pow(knn_sig_err/knn_sig_mean,2) ) ; }
         if ( (Nsig_mm+Nsb_mm) > 0 ) { comp_znn_sb_mm_err  = comp_znn_sb_mm  * sqrt( 1.0/(1.0*(Nsig_mm+Nsb_mm)) + pow(fsig_mm_err/fsig_mm_mean,2) + pow(acc_mm_err/acc_mm_mean,2) + pow(eff_mm_err/eff_mm_mean,2) + pow(knn_sig_err/knn_sig_mean,2) ) ; }


       }

       //-- really dumb ave.
       comp_znn_sig = 0.5 * ( comp_znn_sig_ee + comp_znn_sig_mm ) ;
       comp_znn_sig_err = comp_znn_sig_ee_err ;
       if ( comp_znn_sig_mm_err > comp_znn_sig_ee_err ) comp_znn_sig_err = comp_znn_sig_mm_err ;

       comp_znn_sb = 0.5 * ( comp_znn_sb_ee + comp_znn_sb_mm ) ;
       comp_znn_sb_err = comp_znn_sb_ee_err ;
       if ( comp_znn_sb_mm_err > comp_znn_sb_ee_err ) comp_znn_sb_err = comp_znn_sb_mm_err ;


       printf("\n\n\n") ;

       printf("------------+----------+----------+------------------------+----------+----------+------------+------------+---------------\n") ;
       printf("   Sample   |  ttbar   |  W+jets  |          QCD           |  Z to nn |   other  |   All SM   |    Data    |     SUSY      \n") ;
       printf("------------+----------+----------+------------------------+----------+----------+------------+------------+---------------\n") ;
       printf("            |          |          |                        |          |          |            |            |               \n") ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "sig", Nttbarmc_sig, (lsf_WJmc*NWJmc_sig), Nqcdmc_sig, Nqcdmc_sig_err, lsf_Znnmc*NZnnmc_sig, lsf_Ewomc*NEwomc_sig, Nsm_sig, Nsig, Nsusymc_sig ) ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "sb", Nttbarmc_sb, (lsf_WJmc*NWJmc_sb), Nqcdmc_sb, Nqcdmc_sb_err, lsf_Znnmc*NZnnmc_sb, lsf_Ewomc*NEwomc_sb, Nsm_sb, Nsb, Nsusymc_sb ) ;
       printf("            |          |          |                        |          |          |            |            |               \n") ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "sig_sl", Nttbarmc_sig_sl, (lsf_WJmc*NWJmc_sig_sl), Nqcdmc_sig_sl, Nqcdmc_sig_sl_err, lsf_Znnmc*NZnnmc_sig_sl, lsf_Ewomc*NEwomc_sig_sl, Nsm_sig_sl, Nsig_sl, Nsusymc_sig_sl ) ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "sb_sl", Nttbarmc_sb_sl, (lsf_WJmc*NWJmc_sb_sl), Nqcdmc_sb_sl, Nqcdmc_sb_sl_err, lsf_Znnmc*NZnnmc_sb_sl, lsf_Ewomc*NEwomc_sb_sl, Nsm_sb_sl, Nsb_sl, Nsusymc_sb_sl ) ;
       printf("            |          |          |                        |          |          |            |            |               \n") ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "sig_ldp", Nttbarmc_sig_ldp, (lsf_WJmc*NWJmc_sig_ldp), Nqcdmc_sig_ldp, Nqcdmc_sig_ldp_err, lsf_Znnmc*NZnnmc_sig_ldp, lsf_Ewomc*NEwomc_sig_ldp, Nsm_sig_ldp, Nsig_ldp, Nsusymc_sig_ldp ) ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "sb_ldp", Nttbarmc_sb_ldp, (lsf_WJmc*NWJmc_sb_ldp), Nqcdmc_sb_ldp, Nqcdmc_sb_ldp_err, lsf_Znnmc*NZnnmc_sb_ldp, lsf_Ewomc*NEwomc_sb_ldp, Nsm_sb_ldp, Nsb_ldp, Nsusymc_sb_ldp ) ;
       printf("            |          |          |                        |          |          |            |            |               \n") ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "lsb", Nttbarmc_lsb, (lsf_WJmc*NWJmc_lsb), Nqcdmc_lsb, Nqcdmc_lsb_err, lsf_Znnmc*NZnnmc_lsb, lsf_Ewomc*NEwomc_lsb, Nsm_lsb, Nlsb, Nsusymc_lsb ) ;
       printf(" %10s | %8.1f | %8.1f | %10.1f +/- %7.1f | %8.1f | %8.1f | %10.1f | %10d | %8.1f\n",
            "lsb_ldp", Nttbarmc_lsb_ldp, (lsf_WJmc*NWJmc_lsb_ldp), Nqcdmc_lsb_ldp, Nqcdmc_lsb_ldp_err, lsf_Znnmc*NZnnmc_lsb_ldp, lsf_Ewomc*NEwomc_lsb_ldp, Nsm_lsb_ldp, Nlsb_ldp, Nsusymc_lsb_ldp ) ;
       printf("            |          |          |                        |          |          |            |            |               \n") ;
       printf(" %10s | -------- | -------- | %10.1f +/- %7.1f | -------- | -------- | %10.1f | %10d | %8.1f\n",
            "lsb_0b",  Nqcdmc_lsb_0b, Nqcdmc_lsb_0b_err, Nsm_lsb_0b, Nlsb_0b, Nsusymc_lsb_0b ) ;
       printf(" %10s | -------- | -------- | %10.1f +/- %7.1f | -------- | -------- | %10.1f | %10d | %8.1f\n",
            "lsb_0b_ldp",  Nqcdmc_lsb_0b_ldp, Nqcdmc_lsb_0b_ldp_err, Nsm_lsb_0b_ldp, Nlsb_0b_ldp, Nsusymc_lsb_0b_ldp ) ;
       printf("            |          |          |                        |          |          |            |            |               \n") ;
       printf("------------+----------+----------+------------------------+----------+----------+------------+------------+---------------\n") ;


       printf("\n\n\n") ;

       printf(" R QCD  :  sig    /    sig_ldp  : ( %10.1f +/- %6.1f ) / ( %10.1f +/- %7.1f )    =    %5.3f +/- %5.3f\n",
            Nqcdmc_sig, Nqcdmc_sig_err, Nqcdmc_sig_ldp, Nqcdmc_sig_ldp_err, Rqcd_sig, Rqcd_sig_err ) ;

       printf(" R QCD  :  sb     /     sb_ldp  : ( %10.1f +/- %6.1f ) / ( %10.1f +/- %7.1f )    =    %5.3f +/- %5.3f\n",
            Nqcdmc_sb, Nqcdmc_sb_err, Nqcdmc_sb_ldp, Nqcdmc_sb_ldp_err, Rqcd_sb, Rqcd_sb_err ) ;

       printf(" R QCD  :  lsb    /    lsb_ldp  : ( %10.1f +/- %6.1f ) / ( %10.1f +/- %7.1f )    =    %5.3f +/- %5.3f\n",
            Nqcdmc_lsb, Nqcdmc_lsb_err, Nqcdmc_lsb_ldp, Nqcdmc_lsb_ldp_err, Rqcd_lsb, Rqcd_lsb_err ) ;

       printf(" R data :  lsb    /    lsb_ldp  : ( %10.1f            ) / ( %10.1f             )    =    %5.3f +/- %5.3f\n",
            Nlsb_corrected, Nlsb_ldp_corrected, Rdata_lsb, Rdata_lsb_err ) ;

       printf(" R QCD  :  lsb_0b / lsb_0b_ldp  : ( %10.1f +/- %6.1f ) / ( %10.1f +/- %7.1f )    =    %5.3f +/- %5.3f\n",
            Nqcdmc_lsb_0b, Nqcdmc_lsb_0b_err, Nqcdmc_lsb_0b_ldp, Nqcdmc_lsb_0b_ldp_err, Rqcd_lsb_0b, Rqcd_lsb_0b_err ) ;

       printf(" R data :  lsb_0b / lsb_0b_ldp  : ( %8d              ) / ( %8d               )    =    %5.3f +/- %5.3f\n",
            Nlsb_0b, Nlsb_0b_ldp, Rdata_lsb_0b, Rdata_lsb_0b_err ) ;

       printf(" R data :  lsb_0b / lsb_0b_ldp  : ( %8d              ) / ( %8d               )    =    %5.3f +/- %5.3f (HT-only trigger)\n",
            Nhtonlytrig_lsb_0b, Nhtonlytrig_lsb_0b_ldp, RdataHTOT_lsb_0b, RdataHTOT_lsb_0b_err ) ;



       printf("\n\n\n") ;


       printf(" R ttwj :  sig    /    sb    :  ( %5.1f +/- %4.1f ) / ( %5.1f +/- %4.1f )   =   %5.3f +/- %5.3f\n",
                 (Nttbarmc_sig+lsf_WJmc*NWJmc_sig), lsf_WJmc*sqrt(NWJmc_sig),
                 (Nttbarmc_sb+lsf_WJmc*NWJmc_sb), lsf_WJmc*sqrt(NWJmc_sb),
                 Rttwj, Rttwj_err ) ;

       printf(" R ttwj :  sig_sl /    sb_sl :  ( %5.1f +/- %4.1f ) / ( %5.1f +/- %4.1f )   =   %5.3f +/- %5.3f\n",
                 (Nttbarmc_sig_sl+lsf_WJmc*NWJmc_sig_sl), lsf_WJmc*sqrt(NWJmc_sig_sl),
                 (Nttbarmc_sb_sl+lsf_WJmc*NWJmc_sb_sl), lsf_WJmc*sqrt(NWJmc_sb_sl),
                 Rttwj_sl, Rttwj_sl_err ) ;

       printf(" R data :  sig_sl /    sb_sl :    %3d              /   %3d                =   %5.3f +/- %5.3f\n",
                 Nsig_sl, Nsb_sl, Rdata_sl, Rdata_sl_err ) ;



       printf("\n\n\n") ;


       printf(" ----- Simple MC closure tests\n\n") ;

       printf("  QCD :  sb =  sb_ldp * (lsb_0b/lsb_0b_ldp) : ( %6.1f +/- %4.1f ) * ( %8.1f / %10.1f ) = %8.1f +/- %6.1f\n",
             Nqcdmc_sb_ldp, Nqcdmc_sb_ldp_err, Nqcdmc_lsb_0b, Nqcdmc_lsb_0b_ldp, comp_mc_qcd_sb, comp_mc_qcd_sb_err ) ;
       printf("                                                                               MC truth value : %8.1f +/- %6.1f\n",
                   Nqcdmc_sb, Nqcdmc_sb_err ) ;

       printf("\n") ;
       printf("  QCD : sig = sig_ldp * (lsb_0b/lsb_0b_ldp) : ( %6.1f +/- %4.1f ) * ( %8.1f / %10.1f ) = %8.1f +/- %6.1f\n",
             Nqcdmc_sig_ldp, Nqcdmc_sig_ldp_err, Nqcdmc_lsb_0b, Nqcdmc_lsb_0b_ldp, comp_mc_qcd_sig, comp_mc_qcd_sig_err ) ;
       printf("                                                                               MC truth value : %8.1f +/- %6.1f\n",
                   Nqcdmc_sig, Nqcdmc_sig_err ) ;


       printf("\n") ;
       printf("  ttwj : sig = ttwj_sb * ( ttwj_sig_sl / ttwj_sb_sl ) :  %6.1f * ( %6.1f / %6.1f ) = %6.1f +/- %4.1f\n",
              (Nttbarmc_sb + lsf_WJmc*NWJmc_sb), (Nttbarmc_sig_sl + lsf_WJmc*NWJmc_sig_sl), (Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl),
              comp_mc_ttwj_sig, comp_mc_ttwj_sig_err ) ;
       printf("                                                                       MC truth value : %6.1f +/- %4.1f\n",
                (Nttbarmc_sig + lsf_WJmc*NWJmc_sig), lsf_WJmc*sqrt(NWJmc_sig) ) ;


       printf("\n\n\n") ;

       printf(" ----- Simple Data calculations\n\n") ;


       printf(" QCD :  sb  = (  sb_ldp - (ttwj+znn+other)_ldp_mc ) * (lsb_0b/lsb_0b_ldp) : ( %4d - %5.1f ) * ( %5d / %7d ) = %6.1f +/- %4.1f\n",
           Nsb_ldp, (Nttbarmc_sb_ldp  + lsf_WJmc*NWJmc_sb_ldp  + lsf_Znnmc*NZnnmc_sb_ldp  + lsf_Ewomc*NEwomc_sb_ldp  ),
           Nlsb_0b, Nlsb_0b_ldp, comp_data_qcd_sb, comp_data_qcd_sb_err ) ;
       printf("                                                                                                  MC truth value : %8.1f +/- %4.1f\n",
                   Nqcdmc_sb, Nqcdmc_sb_err ) ;

       printf("\n") ;
       printf(" QCD :  sig = ( sig_ldp - (ttwj+znn+other)_ldp_mc ) * (lsb_0b/lsb_0b_ldp) : ( %4d - %5.1f ) * ( %5d / %7d ) = %6.1f +/- %4.1f\n",
           Nsig_ldp, (Nttbarmc_sig_ldp  + lsf_WJmc*NWJmc_sig_ldp  + lsf_Znnmc*NZnnmc_sig_ldp  + lsf_Ewomc*NEwomc_sig_ldp  ),
           Nlsb_0b, Nlsb_0b_ldp, comp_data_qcd_sig, comp_data_qcd_sig_err ) ;
       printf("                                                                                                  MC truth value : %8.1f +/- %4.1f\n",
                   Nqcdmc_sig, Nqcdmc_sig_err ) ;


       printf("\n") ;
       printf(" ttwj : sig = ( sb - (qcd_sb + (znn+other)_sb_mc) ) * (sig_sl/sb_sl) : ( %5d - ( %5.1f + %5.1f))*(%3d/%3d) = %5.1f +/- %4.1f\n",
              Nsb, comp_data_qcd_sb, (lsf_Znnmc*NZnnmc_sb + lsf_Ewomc*NEwomc_sb), Nsig_sl, Nsb_sl, comp_data_ttwj_sig, comp_data_ttwj_sig_err ) ;
       printf("                                                                                              MC truth value : %5.1f +/- %4.1f\n",
                 (Nttbarmc_sig + lsf_WJmc*NWJmc_sig), lsf_WJmc*sqrt(NWJmc_sig) ) ;


       printf("\n\n\n") ;

       printf(" ---- Z to nunu calculations\n\n" ) ;

       printf("  Znn SIG : %6.1f +/- %5.1f\n", comp_znn_sig, comp_znn_sig_err ) ;
       printf("  Znn SB  : %6.1f +/- %5.1f\n", comp_znn_sb , comp_znn_sb_err ) ;


       printf("\n\n\n") ;



     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


       printf(" --- Defining observables.\n" ) ;


      rv_Nsig        = new RooRealVar( "Nsig"        , "Nsig"        , 0.0, 1000000. ) ;
      rv_Nsb         = new RooRealVar( "Nsb"         , "Nsb"         , 0.0, 1000000. ) ;

      rv_Nsig_sl     = new RooRealVar( "Nsig_sl"     , "Nsig_sl"     , 0.0, 1000000. ) ;
      rv_Nsb_sl      = new RooRealVar( "Nsb_sl"      , "Nsb_sl"      , 0.0, 1000000. ) ;

      rv_Nsig_ldp    = new RooRealVar( "Nsig_ldp"    , "Nsig_ldp"    , 0.0, 1000000. ) ;
      rv_Nsb_ldp     = new RooRealVar( "Nsb_ldp"     , "Nsb_ldp"     , 0.0, 1000000. ) ;

      rv_Nlsb_0b     = new RooRealVar( "Nlsb_0b"     , "Nlsb_0b"     , 0.0, 1000000. ) ;
      rv_Nlsb_0b_ldp = new RooRealVar( "Nlsb_0b_ldp" , "Nlsb_0b_ldp" , 0.0, 1000000. ) ;

      if ( znnModel == 1 ) {

         rv_Nsb_ee      = new RooRealVar( "Nsb_ee"      ,"Nsb_ee"       , 0., 100. ) ;
         rv_Nsig_ee     = new RooRealVar( "Nsig_ee"     ,"Nsig_ee"      , 0., 100. ) ;

         rv_Nsb_mm      = new RooRealVar( "Nsb_mm"      ,"Nsb_mm"       , 0., 100. ) ;
         rv_Nsig_mm     = new RooRealVar( "Nsig_mm"     ,"Nsig_mm"      , 0., 100. ) ;

      } else if ( znnModel == 2 ) {

         rv_Nsigsb_ee     = new RooRealVar( "Nsigsb_ee"     ,"Nsigsb_ee"      , 0., 100. ) ;
         rv_Nsigsb_mm     = new RooRealVar( "Nsigsb_mm"     ,"Nsigsb_mm"      , 0., 100. ) ;

      }





      if ( Nsig < 0 ) {
         printf("\n\n *** Negative value for Nsig in input file.  Will set Nsig to MC expectation, which is %d.\n\n",
             TMath::Nint( Nsm_sig ) ) ;
         Nsig = TMath::Nint( Nsm_sig ) ;
      }

      rv_Nsig        -> setVal( Nsig ) ;
      rv_Nsb         -> setVal( Nsb ) ;

      rv_Nsig_sl     -> setVal( Nsig_sl ) ;
      rv_Nsb_sl      -> setVal( Nsb_sl ) ;

      rv_Nsig_ldp    -> setVal( Nsig_ldp ) ;
      rv_Nsb_ldp     -> setVal( Nsb_ldp ) ;

      rv_Nlsb_0b     -> setVal( Nhtonlytrig_lsb_0b ) ;
      rv_Nlsb_0b_ldp -> setVal( Nhtonlytrig_lsb_0b_ldp ) ;

      if ( znnModel == 1 ) {

         rv_Nsb_ee      -> setVal( Nsb_ee ) ;
         rv_Nsig_ee     -> setVal( Nsig_ee ) ;

         rv_Nsb_mm      -> setVal( Nsb_mm ) ;
         rv_Nsig_mm     -> setVal( Nsig_mm ) ;

      } else if ( znnModel == 2 ) {

         rv_Nsigsb_ee     -> setVal( Nsig_ee + Nsb_ee ) ;
         rv_Nsigsb_mm     -> setVal( Nsig_mm + Nsb_mm ) ;

      }







    //++++++++ Parameters of the likelihood +++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining parameters.\n" ) ;



    //____ Counts in SIG ______________________

      if ( useSigTtwjVar ) {
         rrv_mu_ttwj_sig = new RooRealVar( "mu_ttwj_sig"   , "mu_ttwj_sig"   , 0.0, 10000. ) ;
         rv_mu_ttwj_sig = rrv_mu_ttwj_sig ;
         rrv_mu_ttwj_sig   -> setVal( Nttbarmc_sig + lsf_WJmc*NWJmc_sig ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sig  = new RooRealVar( "mu_qcd_sig"    , "mu_qcd_sig"    , 0.0, 50. ) ;
         rv_mu_qcd_sig = rrv_mu_qcd_sig ;
         rrv_mu_qcd_sig  -> setVal( Nqcdmc_sig ) ; //-- this is a starting value only.
      }


      //-- Note: Ewo is rfv

      rv_mu_znn_sig      = new RooRealVar( "mu_znn_sig"    , "mu_znn_sig"    , 0.0, 80. ) ;

      float maxSusySig = 4.0*Nsig ;
      rv_mu_susy_sig     = new RooRealVar( "mu_susy_sig"   , "mu_susy_sig"   , 0.0, maxSusySig ) ;


      rv_mu_znn_sig   -> setVal( comp_znn_sig ) ;  //-- this is a starting value only.
      rv_mu_susy_sig    -> setVal( 0. ) ;  //-- this is a starting value only.





    //____ Counts in SB  ______________________


      if ( !useSigTtwjVar ) {
         rrv_mu_ttwj_sb  = new RooRealVar( "mu_ttwj_sb"    , "mu_ttwj_sb"    , 0.0, 10000. ) ;
         rv_mu_ttwj_sb = rrv_mu_ttwj_sb ;
         rrv_mu_ttwj_sb   -> setVal( Nttbarmc_sb + lsf_WJmc*NWJmc_sb ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sb  = new RooRealVar( "mu_qcd_sb"    , "mu_qcd_sb"    , 0.0, 100. ) ;
         rv_mu_qcd_sb = rrv_mu_qcd_sb ;
         rrv_mu_qcd_sb  -> setVal( Nqcdmc_sb ) ; //-- this is a starting value only.
      }

      //-- Note: QCD is rfv
      //-- Note: Ewo is rfv
      //-- Note: SUSY is rfv

      if ( znnModel == 1 ) {

         rrv_mu_znn_sb       = new RooRealVar( "mu_znn_sb"     , "mu_znn_sb"     , 0.0, 150. ) ;

         rrv_mu_znn_sb   -> setVal( comp_znn_sb ) ;  //-- this is a starting value only.

         rv_mu_znn_sb = rrv_mu_znn_sb ;
      }
      //-- Note: Znn is rfv in Znn model 2.






    //____ Counts in SIG, SL  ______________________

      rv_mu_ttwj_sig_sl  = new RooRealVar( "mu_ttwj_sig_sl"    , "mu_ttwj_sig_sl"    , 0.0, 500. ) ;

      //-- Note: QCD, Ewo, and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      rv_mu_ttwj_sig_sl  -> setVal( Nttbarmc_sig_sl + lsf_WJmc*NWJmc_sig_sl ) ;  //-- this is a starting value only.







    //____ Counts in SB, SL  ______________________

      rv_mu_ttwj_sb_sl  = new RooRealVar( "mu_ttwj_sb_sl"    , "mu_ttwj_sb_sl"    , 0.0, 3000. ) ;

      //-- Note: QCD, Ewo, and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      rv_mu_ttwj_sb_sl  -> setVal( Nttbarmc_sb_sl + lsf_WJmc*NWJmc_sb_sl ) ;  //-- this is a starting value only.







    //____ Counts in SIG, LDP  ______________________

      if ( useLdpVars ) {
         rrv_mu_qcd_sig_ldp  = new RooRealVar( "mu_qcd_sig_ldp"    , "mu_qcd_sig_ldp"    , 0.0, 500. ) ;
         rv_mu_qcd_sig_ldp = rrv_mu_qcd_sig_ldp ;
         rrv_mu_qcd_sig_ldp  -> setVal( Nqcdmc_sig_ldp ) ; //-- this is a starting value only.
      }

      //-- Note: Ewo is assumed to be negligible and is not explicitly included.
      //-- Note: Znn is rfv (MC)
      //-- Note: ttwj is rfv (MC)
      //-- Note: SUSY is rfv








    //____ Counts in SB, LDP  ______________________

      if ( useLdpVars ) {
         rrv_mu_qcd_sb_ldp  = new RooRealVar( "mu_qcd_sb_ldp"    , "mu_qcd_sb_ldp"    , 0.0, 1000. ) ;
         rv_mu_qcd_sb_ldp = rrv_mu_qcd_sb_ldp ;
         rrv_mu_qcd_sb_ldp  -> setVal( Nqcdmc_sb_ldp ) ; //-- this is a starting value only.
      }


      //-- Note: Ewo is assumed to be negligible and is not explicitly included.
      //-- Note: Znn is rfv (MC)
      //-- Note: ttwj is rfv (MC)
      //-- Note: SUSY is rfv








    //____ Counts in LSB, 0b  ______________________

      rv_mu_qcd_lsb_0b  = new RooRealVar( "mu_qcd_lsb_0b"    ,  "mu_qcd_lsb_0b" ,  0.    ,  10000. ) ;

      //-- Note: The 0btag LSB is assumed to be 100% QCD.

      rv_mu_qcd_lsb_0b  -> setVal( Nhtonlytrig_lsb_0b ) ;  //-- this is a starting value only.






    //____ Counts in LSB, 0b, LDP  ______________________

      rv_mu_qcd_lsb_0b_ldp  = new RooRealVar( "mu_qcd_lsb_0b_ldp"    ,  "mu_qcd_lsb_0b_ldp",   0. ,  10000. ) ;

      //-- Note: The 0btag LSB, LDP is assumed to be 100% QCD.

      rv_mu_qcd_lsb_0b_ldp  -> setVal( Nhtonlytrig_lsb_0b_ldp ) ;  //-- this is a starting value only.






    //____ Counts in SB, ee  ______________________


      //-- Note: The Z to ee sample is assumed to be 100% Z to ee.
      //-- Note: zee is rfv






    //____ Counts in SIG, ee  ______________________


      //-- Note: The Z to ee sample is assumed to be 100% Z to ee.
      //-- Note: zee is rfv








    //____ Counts in SB, mm  ______________________


      //-- Note: The Z to mm sample is assumed to be 100% Z to mm.
      //-- Note: zmm is rfv






    //____ Counts in SIG, mm  ______________________


      //-- Note: The Z to mm sample is assumed to be 100% Z to mm.
      //-- Note: zmm is rfv






    //____ MC inputs _______________________________

      printf(" --- Defining MC parameters.\n" ) ;

     //-- SUSY

      rv_mu_susymc_sig      = new RooRealVar( "mu_susymc_sig"     , "mu_susymc_sig"     , 0.0, 100000. ) ;
      rv_mu_susymc_sb       = new RooRealVar( "mu_susymc_sb"      , "mu_susymc_sb"      , 0.0, 100000. ) ;
      rv_mu_susymc_sig_sl   = new RooRealVar( "mu_susymc_sig_sl"  , "mu_susymc_sig_sl"  , 0.0, 100000. ) ;
      rv_mu_susymc_sb_sl    = new RooRealVar( "mu_susymc_sb_sl"   , "mu_susymc_sb_sl"   , 0.0, 100000. ) ;
      rv_mu_susymc_sig_ldp  = new RooRealVar( "mu_susymc_sig_ldp" , "mu_susymc_sig_ldp" , 0.0, 100000. ) ;
      rv_mu_susymc_sb_ldp   = new RooRealVar( "mu_susymc_sb_ldp"  , "mu_susymc_sb_ldp"  , 0.0, 100000. ) ;

      rv_mu_susymc_sig     -> setVal( 0.1 ) ;
      rv_mu_susymc_sb      -> setVal( 0. ) ;
      rv_mu_susymc_sig_sl  -> setVal( 0. ) ;
      rv_mu_susymc_sb_sl   -> setVal( 0. ) ;
      rv_mu_susymc_sig_ldp -> setVal( 0. ) ;
      rv_mu_susymc_sb_ldp  -> setVal( 0. ) ;

      rv_mu_susymc_sig     -> setConstant(kTRUE) ;
      rv_mu_susymc_sb      -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_sl  -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_sl   -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_ldp -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_ldp  -> setConstant(kTRUE) ;



     //-- SIG, LDP

      rv_mu_ttbarmc_sig_ldp   = new RooRealVar( "mu_ttbarmc_sig_ldp" ,"mu_ttbarmc_sig_ldp" , 0., 1000. ) ;
      rv_mu_WJmc_sig_ldp      = new RooRealVar( "mu_WJmc_sig_ldp"    ,"mu_WJmc_sig_ldp"    , 0., 1000. ) ;
      rv_mu_Znnmc_sig_ldp     = new RooRealVar( "mu_Znnmc_sig_ldp"   ,"mu_Znnmc_sig_ldp"   , 0., 1000. ) ;
      rv_mu_Ewomc_sig_ldp     = new RooRealVar( "mu_Ewomc_sig_ldp"   ,"mu_Ewomc_sig_ldp"   , 0., 1000. ) ;

      rv_mu_ttbarmc_sig_ldp  -> setVal( Nttbarmc_sig_ldp ) ;
      rv_mu_WJmc_sig_ldp     -> setVal( NWJmc_sig_ldp ) ;
      rv_mu_Znnmc_sig_ldp    -> setVal( NZnnmc_sig_ldp ) ;
      rv_mu_Ewomc_sig_ldp    -> setVal( NEwomc_sig_ldp ) ;

      rv_mu_ttbarmc_sig_ldp  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sig_ldp     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sig_ldp    -> setConstant( kTRUE ) ;
      rv_mu_Ewomc_sig_ldp    -> setConstant( kTRUE ) ;


     //-- SB, LDP

      rv_mu_ttbarmc_sb_ldp   = new RooRealVar( "mu_ttbarmc_sb_ldp" ,"mu_ttbarmc_sb_ldp" , 0., 1000. ) ;
      rv_mu_WJmc_sb_ldp      = new RooRealVar( "mu_WJmc_sb_ldp"    ,"mu_WJmc_sb_ldp"    , 0., 1000. ) ;
      rv_mu_Znnmc_sb_ldp     = new RooRealVar( "mu_Znnmc_sb_ldp"   ,"mu_Znnmc_sb_ldp"   , 0., 1000. ) ;
      rv_mu_Ewomc_sb_ldp     = new RooRealVar( "mu_Ewomc_sb_ldp"   ,"mu_Ewomc_sb_ldp"   , 0., 1000. ) ;

      rv_mu_ttbarmc_sb_ldp  -> setVal( Nttbarmc_sb_ldp ) ;
      rv_mu_WJmc_sb_ldp     -> setVal( NWJmc_sb_ldp ) ;
      rv_mu_Znnmc_sb_ldp    -> setVal( NZnnmc_sb_ldp ) ;
      rv_mu_Ewomc_sb_ldp    -> setVal( NEwomc_sb_ldp ) ;

      rv_mu_ttbarmc_sb_ldp  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sb_ldp     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sb_ldp    -> setConstant( kTRUE ) ;
      rv_mu_Ewomc_sb_ldp    -> setConstant( kTRUE ) ;


     //-- SIG

      rv_mu_Ewomc_sig     = new RooRealVar( "mu_Ewomc_sig"   ,"mu_Ewomc_sig"   , 0., 1000. ) ;
      rv_mu_Ewomc_sig    -> setVal( NEwomc_sig ) ;
      rv_mu_Ewomc_sig    -> setConstant( kTRUE ) ;

     //-- SB

      rv_mu_Ewomc_sb     = new RooRealVar( "mu_Ewomc_sb"   ,"mu_Ewomc_sb"   , 0., 1000. ) ;
      rv_mu_Ewomc_sb    -> setVal( NEwomc_sb ) ;
      rv_mu_Ewomc_sb    -> setConstant( kTRUE ) ;





    //+++++++ Gaussian constraints ++++++++++++++++++++++++++++++++

      printf(" --- Defining Gaussian constraint and constant parameters.\n" ) ;

    //_______ Efficiency scale factor.  Applied to SUSY and all MC inputs _______________

    //   August 23, 2011:  switching to correlated log-normal PDFs
    //
    //
    //  double pmin, pmax ;


    //--- mean parameters.
      rv_mean_eff_sf_sig     = new RooRealVar( "mean_eff_sf_sig"    , "mean_eff_sf_sig", 0., 10. ) ;
      rv_mean_eff_sf_sb      = new RooRealVar( "mean_eff_sf_sb"     , "mean_eff_sf_sb", 0., 10. ) ;
      rv_mean_eff_sf_sig_sl  = new RooRealVar( "mean_eff_sf_sig_sl" , "mean_eff_sf_sig_sl", 0., 10. ) ;
      rv_mean_eff_sf_sb_sl   = new RooRealVar( "mean_eff_sf_sb_sl"  , "mean_eff_sf_sb_sl", 0., 10. ) ;
      rv_mean_eff_sf_sig_ldp = new RooRealVar( "mean_eff_sf_sig_ldp", "mean_eff_sf_sig_ldp", 0., 10. ) ;
      rv_mean_eff_sf_sb_ldp  = new RooRealVar( "mean_eff_sf_sb_ldp" , "mean_eff_sf_sb_ldp", 0., 10. ) ;

    //--- width parameters.
      rv_width_eff_sf_sig     = new RooRealVar( "width_eff_sf_sig"    , "width_eff_sf_sig", 0., 10. ) ;
      rv_width_eff_sf_sb      = new RooRealVar( "width_eff_sf_sb"     , "width_eff_sf_sb", 0., 10. ) ;
      rv_width_eff_sf_sig_sl  = new RooRealVar( "width_eff_sf_sig_sl" , "width_eff_sf_sig_sl", 0., 10. ) ;
      rv_width_eff_sf_sb_sl   = new RooRealVar( "width_eff_sf_sb_sl"  , "width_eff_sf_sb_sl", 0., 10. ) ;
      rv_width_eff_sf_sig_ldp = new RooRealVar( "width_eff_sf_sig_ldp", "width_eff_sf_sig_ldp", 0., 10. ) ;
      rv_width_eff_sf_sb_ldp  = new RooRealVar( "width_eff_sf_sb_ldp" , "width_eff_sf_sb_ldp", 0., 10. ) ;


      rv_mean_eff_sf_sig     -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb      -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_sl  -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_sl   -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_ldp -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_ldp  -> setVal( 1.00 ) ;

      rv_mean_eff_sf_sig     -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb      -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_sl  -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_sl   -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_ldp -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_ldp  -> setConstant( kTRUE ) ;


      //--- Initialize all width parameters to 15%.  They will be reset in the susy scan.
      rv_width_eff_sf_sig     -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb      -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_sl  -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_sl   -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_ldp -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_ldp  -> setVal( 0.15 ) ;

      rv_width_eff_sf_sig     -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb      -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_sl  -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_sl   -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_ldp -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_ldp  -> setConstant( kTRUE ) ;







   //
   // Owen : Sept 7, 2011: need to set these here so that it gets into the workspace.
   //
   //++++++++++
      setSusyScanPoint( inputScanFile,  m0,  m12,  isT1bbbb,  t1bbbbXsec ) ;
   //++++++++++









    //--- Underlying Gaussian variable for log-normal.

   //--- Owen, Aug 30: increase range from +-5 to +-15
      rv_eff_sf_prim = new RooRealVar("eff_sf_prim", "eff_sf_prim", 0., -15., 15. ) ;
      rv_eff_sf_nom  = new RooRealVar("eff_sf_nom" , "eff_sf_nom" , 0., -15., 15. ) ;
      rv_eff_sf_nom->setConstant() ;



    //--- Systematics
      char formula[1024];
      RooArgSet globalObservables ("globalObservables");
      RooArgSet allNuisances ("allNuisances");
      RooArgSet allNuisancePdfs ("allNuisancePdfs");

      RooRealVar acc_ee_prim ( "acc_ee_prim", "acc_ee_prim", 0, -5, 5);
      RooRealVar acc_ee_nom ( "acc_ee_nom", "acc_ee_nom", 0, -5, 5);
      RooGaussian pdf_acc_ee ("pdf_acc_ee" , "pdf_acc_ee", acc_ee_prim, acc_ee_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", acc_ee_mean, exp(acc_ee_err/acc_ee_mean));
      RooFormulaVar fv_acc_ee ("acc_ee", formula, RooArgList(acc_ee_prim));
      acc_ee_nom.setConstant();
      globalObservables.add (acc_ee_nom);
      allNuisances.add (acc_ee_prim);
      allNuisancePdfs.add (pdf_acc_ee);

      RooRealVar acc_mm_prim ( "acc_mm_prim", "acc_mm_prim", 0, -5, 5);
      RooRealVar acc_mm_nom ( "acc_mm_nom", "acc_mm_nom", 0, -5, 5);
      RooGaussian pdf_acc_mm ("pdf_acc_mm" , "pdf_acc_mm", acc_mm_prim, acc_mm_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", acc_mm_mean, exp(acc_mm_err/acc_mm_mean));
      RooFormulaVar fv_acc_mm ("acc_mm", formula, RooArgList(acc_mm_prim));
      acc_mm_nom.setConstant();
      globalObservables.add (acc_mm_nom);
      allNuisances.add (acc_mm_prim);
      allNuisancePdfs.add (pdf_acc_mm);

      RooRealVar eff_ee_prim ( "eff_ee_prim", "eff_ee_prim", 0, -5, 5);
      RooRealVar eff_ee_nom ( "eff_ee_nom", "eff_ee_nom", 0, -5, 5);
      RooGaussian pdf_eff_ee ("pdf_eff_ee" , "pdf_eff_ee", eff_ee_prim, eff_ee_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", eff_ee_mean, exp(eff_ee_err/eff_ee_mean));
      RooFormulaVar fv_eff_ee ("eff_ee", formula, RooArgList(eff_ee_prim));
      eff_ee_nom.setConstant();
      globalObservables.add (eff_ee_nom);
      allNuisances.add (eff_ee_prim);
      allNuisancePdfs.add (pdf_eff_ee);

      RooRealVar eff_mm_prim ( "eff_mm_prim", "eff_mm_prim", 0, -5, 5);
      RooRealVar eff_mm_nom ( "eff_mm_nom", "eff_mm_nom", 0, -5, 5);
      RooGaussian pdf_eff_mm ("pdf_eff_mm" , "pdf_eff_mm", eff_mm_prim, eff_mm_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", eff_mm_mean, exp(eff_mm_err/eff_mm_mean));
      RooFormulaVar fv_eff_mm ("eff_mm", formula, RooArgList(eff_mm_prim));
      eff_mm_nom.setConstant();
      globalObservables.add (eff_mm_nom);
      allNuisances.add (eff_mm_prim);
      allNuisancePdfs.add (pdf_eff_mm);

      RooRealVar fsig_ee_prim ( "fsig_ee_prim", "fsig_ee_prim", 0, -5, 5);
      RooRealVar fsig_ee_nom ( "fsig_ee_nom", "fsig_ee_nom", 0, -5, 5);
      RooGaussian pdf_fsig_ee ("pdf_fsig_ee" , "pdf_fsig_ee", fsig_ee_prim, fsig_ee_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", fsig_ee_mean, exp(fsig_ee_err/fsig_ee_mean));
      RooFormulaVar fv_fsig_ee ("fsig_ee", formula, RooArgList(fsig_ee_prim));
      fsig_ee_nom.setConstant();
      globalObservables.add (fsig_ee_nom);
      allNuisances.add (fsig_ee_prim);
      allNuisancePdfs.add (pdf_fsig_ee);

      RooRealVar fsig_mm_prim ( "fsig_mm_prim", "fsig_mm_prim", 0, -5, 5);
      RooRealVar fsig_mm_nom ( "fsig_mm_nom", "fsig_mm_nom", 0, -5, 5);
      RooGaussian pdf_fsig_mm ("pdf_fsig_mm" , "pdf_fsig_mm", fsig_mm_prim, fsig_mm_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", fsig_mm_mean, exp(fsig_mm_err/fsig_mm_mean));
      RooFormulaVar fv_fsig_mm ("fsig_mm", formula, RooArgList(fsig_mm_prim));
      fsig_mm_nom.setConstant();
      globalObservables.add (fsig_mm_nom);
      allNuisances.add (fsig_mm_prim);
      allNuisancePdfs.add (pdf_fsig_mm);

      RooRealVar sf_ee_prim ( "sf_ee_prim", "sf_ee_prim", 0, -5, 5);
      RooRealVar sf_ee_nom ( "sf_ee_nom", "sf_ee_nom", 0, -5, 5);
      RooGaussian pdf_sf_ee ("pdf_sf_ee" , "pdf_sf_ee", sf_ee_prim, sf_ee_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_ee, exp(sf_ee_err/sf_ee));
      RooFormulaVar fv_sf_ee ("sf_ee", formula, RooArgList(sf_ee_prim));
      sf_ee_nom.setConstant();
      globalObservables.add (sf_ee_nom);
      allNuisances.add (sf_ee_prim);
      allNuisancePdfs.add (pdf_sf_ee);

      RooRealVar sf_mc_prim ( "sf_mc_prim", "sf_mc_prim", 0, -5, 5);
      RooRealVar sf_mc_nom ( "sf_mc_nom", "sf_mc_nom", 0, -5, 5);
      RooGaussian pdf_sf_mc ("pdf_sf_mc" , "pdf_sf_mc", sf_mc_prim, sf_mc_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_mc, exp(sf_mc_err/sf_mc));
      RooFormulaVar fv_sf_mc ("sf_mc", formula, RooArgList(sf_mc_prim));
      sf_mc_nom.setConstant();
      globalObservables.add (sf_mc_nom);
      allNuisances.add (sf_mc_prim);
      allNuisancePdfs.add (pdf_sf_mc);

      RooRealVar sf_mm_prim ( "sf_mm_prim", "sf_mm_prim", 0, -5, 5);
      RooRealVar sf_mm_nom ( "sf_mm_nom", "sf_mm_nom", 0, -5, 5);
      RooGaussian pdf_sf_mm ("pdf_sf_mm" , "pdf_sf_mm", sf_mm_prim, sf_mm_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_mm, exp(sf_mm_err/sf_mm));
      RooFormulaVar fv_sf_mm ("sf_mm", formula, RooArgList(sf_mm_prim));
      sf_mm_nom.setConstant();
      globalObservables.add (sf_mm_nom);
      allNuisances.add (sf_mm_prim);
      allNuisancePdfs.add (pdf_sf_mm);

      RooRealVar sf_qcd_sb_prim ( "sf_qcd_sb_prim", "sf_qcd_sb_prim", 0, -5, 5);
      RooRealVar sf_qcd_sb_nom ( "sf_qcd_sb_nom", "sf_qcd_sb_nom", 0, -5, 5);
      RooGaussian pdf_sf_qcd_sb ("pdf_sf_qcd_sb" , "pdf_sf_qcd_sb", sf_qcd_sb_prim, sf_qcd_sb_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_qcd_sb, exp(sf_qcd_sb_err/sf_qcd_sb));
      RooFormulaVar fv_sf_qcd_sb ("sf_qcd_sb", formula, RooArgList(sf_qcd_sb_prim));
      sf_qcd_sb_nom.setConstant();
      globalObservables.add (sf_qcd_sb_nom);
      allNuisances.add (sf_qcd_sb_prim);
      allNuisancePdfs.add (pdf_sf_qcd_sb);

      RooRealVar sf_qcd_sig_prim ( "sf_qcd_sig_prim", "sf_qcd_sig_prim", 0, -5, 5);
      RooRealVar sf_qcd_sig_nom ( "sf_qcd_sig_nom", "sf_qcd_sig_nom", 0, -5, 5);
      RooGaussian pdf_sf_qcd_sig ("pdf_sf_qcd_sig" , "pdf_sf_qcd_sig", sf_qcd_sig_prim, sf_qcd_sig_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_qcd_sig, exp(sf_qcd_sig_err/sf_qcd_sig));
      RooFormulaVar fv_sf_qcd_sig ("sf_qcd_sig", formula, RooArgList(sf_qcd_sig_prim));
      sf_qcd_sig_nom.setConstant();
      globalObservables.add (sf_qcd_sig_nom);
      allNuisances.add (sf_qcd_sig_prim);
      allNuisancePdfs.add (pdf_sf_qcd_sig);

      RooRealVar sf_ttwj_sig_prim ( "sf_ttwj_sig_prim", "sf_ttwj_sig_prim", 0, -5, 5);
      RooRealVar sf_ttwj_sig_nom ( "sf_ttwj_sig_nom", "sf_ttwj_sig_nom", 0, -5, 5);
      RooGaussian pdf_sf_ttwj_sig ("pdf_sf_ttwj_sig" , "pdf_sf_ttwj_sig", sf_ttwj_sig_prim, sf_ttwj_sig_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_ttwj_sig, exp(sf_ttwj_sig_err/sf_ttwj_sig));
      RooFormulaVar fv_sf_ttwj_sig ("sf_ttwj_sig", formula, RooArgList(sf_ttwj_sig_prim));
      sf_ttwj_sig_nom.setConstant();
      globalObservables.add (sf_ttwj_sig_nom);
      allNuisances.add (sf_ttwj_sig_prim);
      allNuisancePdfs.add (pdf_sf_ttwj_sig);

      RooRealVar knn_sig_prim ( "knn_sig_prim", "knn_sig_prim", 0, -5, 5);
      RooRealVar knn_sig_nom ( "knn_sig_nom", "knn_sig_nom", 0, -5, 5);
      RooGaussian pdf_knn_sig ("pdf_knn_sig" , "pdf_knn_sig", knn_sig_prim, knn_sig_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", knn_sig_mean, exp(knn_sig_err/knn_sig_mean));
      RooFormulaVar fv_knn_sig ("knn_sig", formula, RooArgList(knn_sig_prim));
      knn_sig_nom.setConstant();
      if (znnModel == 2 ) {
	globalObservables.add (knn_sig_nom);
	allNuisances.add (knn_sig_prim);
	allNuisancePdfs.add (pdf_knn_sig);
      }

      RooRealVar knn_sb_prim ( "knn_sb_prim", "knn_sb_prim", 0, -5, 5);
      RooRealVar knn_sb_nom ( "knn_sb_nom", "knn_sb_nom", 0, -5, 5);
      RooGaussian pdf_knn_sb ("pdf_knn_sb" , "pdf_knn_sb", knn_sb_prim, knn_sb_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", knn_sb_mean, exp(knn_sb_err/knn_sb_mean));
      RooFormulaVar fv_knn_sb ("knn_sb", formula, RooArgList(knn_sb_prim));
      knn_sb_nom.setConstant();
      if (znnModel == 2 ) {
	globalObservables.add (knn_sb_nom);
	allNuisances.add (knn_sb_prim);
	allNuisancePdfs.add (pdf_knn_sb);
      }


      RooRealVar eff_sf_prim ( "eff_sf_prim", "eff_sf_prim", 0, -5, 5);
      RooRealVar eff_sf_nom ( "eff_sf_nom", "eff_sf_nom", 0, -5, 5);
      RooGaussian pdf_eff_sf_prim ("pdf_eff_sf_prim" , "master pdf_eff_sf_prim", *rv_eff_sf_prim, *rv_eff_sf_nom, RooConst(1));
      sprintf (formula, "pow(%f,@0)", exp(0.10));
      RooFormulaVar fv_Eff_sf_m ("Eff_sf_m", formula, RooArgList(eff_sf_prim));
      rv_eff_sf_nom->setConstant();
      globalObservables.add (*rv_eff_sf_nom);
      allNuisances.add (eff_sf_prim);
      allNuisancePdfs.add (pdf_eff_sf_prim);

/*       RooRealVar xxx_prim ( "xxx_prim", "xxx_prim", 0, -5, 5); */
/*       RooRealVar xxx_nom ( "xxx_nom", "xxx_nom", 0, -5, 5); */
/*       RooGaussian pdf_xxx ("pdf_xxx" , "pdf_xxx", *xxx_prim, *xxx_nom, RooConst(1)); */
/*       sprintf (formula, "%f*pow(%f,@0)", xxx, exp(xxx_err/xxx)); */
/*       RooFormulaVar fv_xxx ("xxx", formula, RooArgList(xxx_prim)); */
/*       globalObservables.add (xxx_nom); */
/*       allNuisances.add (xxx_prim); */
/*       allNuisancePdfs.add (pdf_xxx); */



    //-- Z to nunu stuff



      rv_znnoverll_bfratio = new RooRealVar( "znnoverll_bfratio", "znnoverll_bfratio", 0.1, 10. ) ;
      rv_znnoverll_bfratio -> setVal( 5.95 ) ;
      rv_znnoverll_bfratio -> setConstant( kTRUE ) ;

      rv_dataoverll_lumiratio = new RooRealVar( "dataoverll_lumiratio", "dataoverll_lumiratio", 0.1, 10.0 ) ;
      rv_dataoverll_lumiratio  -> setVal( DataLumi / Ztoll_lumi ) ;
      rv_dataoverll_lumiratio  -> setConstant( kTRUE ) ;




     //+++++++++++++++++ Relationships between parameters ++++++++++++++++++++++++++++++++++++++++++++

       printf(" --- Defining relationships between parameters.\n" ) ;




    //-- ttwj

      if ( useSigTtwjVar ) {
         rfv_mu_ttwj_sb = new RooFormulaVar( "mu_ttwj_sb",
                                                       "mu_ttwj_sig * (1.0/sf_ttwj_sig) * (mu_ttwj_sb_sl/mu_ttwj_sig_sl)",
                                                       RooArgSet( *rv_mu_ttwj_sig, fv_sf_ttwj_sig, *rv_mu_ttwj_sb_sl, *rv_mu_ttwj_sig_sl ) ) ;
         rv_mu_ttwj_sb = rfv_mu_ttwj_sb ;
      } else {
         rfv_mu_ttwj_sig = new RooFormulaVar( "mu_ttwj_sig",
                                                       "mu_ttwj_sb * sf_ttwj_sig * (mu_ttwj_sig_sl/mu_ttwj_sb_sl)",
                                                       RooArgSet( *rv_mu_ttwj_sb, fv_sf_ttwj_sig, *rv_mu_ttwj_sig_sl, *rv_mu_ttwj_sb_sl ) ) ;
         rv_mu_ttwj_sig = rfv_mu_ttwj_sig ;
      }




      rv_mu_ttwj_sig_ldp = new RooFormulaVar( "mu_ttwj_sig_ldp",
                              "mu_ttbarmc_sig_ldp + mu_WJmc_sig_ldp",
                              RooArgSet( *rv_mu_ttbarmc_sig_ldp,
                                         *rv_mu_WJmc_sig_ldp ) ) ;


      rv_mu_ttwj_sb_ldp = new RooFormulaVar( "mu_ttwj_sb_ldp",
                              "mu_ttbarmc_sb_ldp + mu_WJmc_sb_ldp",
                              RooArgSet( *rv_mu_ttbarmc_sb_ldp,
                                         *rv_mu_WJmc_sb_ldp ) ) ;






    //-- QCD

      if ( useLdpVars ) {

         rfv_mu_qcd_sig = new RooFormulaVar( "mu_qcd_sig",
                                     "mu_qcd_sig_ldp * sf_qcd_sig * ( mu_qcd_lsb_0b / mu_qcd_lsb_0b_ldp )",
                                     RooArgSet( *rv_mu_qcd_sig_ldp, fv_sf_qcd_sig, *rv_mu_qcd_lsb_0b, *rv_mu_qcd_lsb_0b_ldp ) ) ;
         rv_mu_qcd_sig = rfv_mu_qcd_sig ;

         rfv_mu_qcd_sb = new RooFormulaVar( "mu_qcd_sb",
                                     "mu_qcd_sb_ldp * sf_qcd_sb * ( mu_qcd_lsb_0b / mu_qcd_lsb_0b_ldp )",
                                     RooArgSet( *rv_mu_qcd_sb_ldp, fv_sf_qcd_sb, *rv_mu_qcd_lsb_0b, *rv_mu_qcd_lsb_0b_ldp ) ) ;
         rv_mu_qcd_sb = rfv_mu_qcd_sb ;

      } else {

         rfv_mu_qcd_sig_ldp = new RooFormulaVar( "mu_qcd_sig_ldp",
                                     "mu_qcd_sig * (1.0/sf_qcd_sig) * ( mu_qcd_lsb_0b_ldp / mu_qcd_lsb_0b )",
                                     RooArgSet( *rv_mu_qcd_sig, fv_sf_qcd_sig, *rv_mu_qcd_lsb_0b_ldp, *rv_mu_qcd_lsb_0b ) ) ;
         rv_mu_qcd_sig_ldp = rfv_mu_qcd_sig_ldp ;

         rfv_mu_qcd_sb_ldp = new RooFormulaVar( "mu_qcd_sb_ldp",
                                     "mu_qcd_sb * (1.0/sf_qcd_sb) * ( mu_qcd_lsb_0b_ldp / mu_qcd_lsb_0b )",
                                     RooArgSet( *rv_mu_qcd_sb, fv_sf_qcd_sb, *rv_mu_qcd_lsb_0b_ldp, *rv_mu_qcd_lsb_0b ) ) ;
         rv_mu_qcd_sb_ldp = rfv_mu_qcd_sb_ldp ;

      }




    //-- SUSY

      rv_mu_susy_sb = new RooFormulaVar( "mu_susy_sb",
                                        "mu_susymc_sb * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sb, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sig_sl = new RooFormulaVar( "mu_susy_sig_sl",
                                        "mu_susymc_sig_sl * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sig_sl, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sb_sl = new RooFormulaVar( "mu_susy_sb_sl",
                                        "mu_susymc_sb_sl * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sb_sl, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sig_ldp = new RooFormulaVar( "mu_susy_sig_ldp",
                                        "mu_susymc_sig_ldp * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sig_ldp, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;

      rv_mu_susy_sb_ldp = new RooFormulaVar( "mu_susy_sb_ldp",
                                        "mu_susymc_sb_ldp * (mu_susy_sig/mu_susymc_sig)",
                                        RooArgSet( *rv_mu_susymc_sb_ldp, *rv_mu_susy_sig, *rv_mu_susymc_sig ) ) ;





    //-- Z to nu nu

      rv_mu_znn_sig_ldp = new RooFormulaVar( "mu_znn_sig_ldp",
                              "mu_Znnmc_sig_ldp",
                              RooArgSet( *rv_mu_Znnmc_sig_ldp ) ) ;

      rv_mu_znn_sb_ldp = new RooFormulaVar( "mu_znn_sb_ldp",
                              "mu_Znnmc_sb_ldp",
                              RooArgSet( *rv_mu_Znnmc_sb_ldp ) ) ;


      if ( znnModel == 1 ) {

         rv_mu_zee_sb_ee = new RooFormulaVar( "mu_zee_sb_ee",
                                         "mu_znn_sb * sf_ee * ( acc_ee * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
                                         RooArgSet( *rv_mu_znn_sb, fv_sf_ee, fv_acc_ee, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

         rv_mu_zee_sig_ee = new RooFormulaVar( "mu_zee_sig_ee",
                                         "mu_znn_sig * sf_ee * ( acc_ee * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
                                         RooArgSet( *rv_mu_znn_sig, fv_sf_ee, fv_acc_ee, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

         rv_mu_zmm_sb_mm = new RooFormulaVar( "mu_zmm_sb_mm",
                                         "mu_znn_sb * sf_mm * ( acc_mm * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
                                         RooArgSet( *rv_mu_znn_sb, fv_sf_mm, fv_acc_mm, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

         rv_mu_zmm_sig_mm = new RooFormulaVar( "mu_zmm_sig_mm",
                                         "mu_znn_sig * sf_mm * ( acc_mm * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio )",
                                         RooArgSet( *rv_mu_znn_sig, fv_sf_mm, fv_acc_mm, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

      } else if ( znnModel == 2 ) {

         rfv_mu_znn_sb = new RooFormulaVar( "mu_znn_sb",
                                          "mu_znn_sig * ( knn_sb / knn_sig )",
                                          RooArgSet( *rv_mu_znn_sig, fv_knn_sb, fv_knn_sig ) ) ;

         rv_mu_znn_sb = rfv_mu_znn_sb ;

         rv_mu_zee_sigsb_ee = new RooFormulaVar( "mu_zee_sigsb_ee",
                                       "( mu_znn_sig / knn_sig ) * sf_ee * ( (acc_ee * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                         RooArgSet( *rv_mu_znn_sig, fv_knn_sig, fv_sf_ee, fv_acc_ee, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

         rv_mu_zmm_sigsb_mm = new RooFormulaVar( "mu_zmm_sigsb_mm",
                                       "( mu_znn_sig / knn_sig ) * sf_mm * ( (acc_mm * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                         RooArgSet( *rv_mu_znn_sig, fv_knn_sig, fv_sf_mm, fv_acc_mm, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

      }

    //-- EWO

      rv_mu_ewo_sig     = new RooFormulaVar( "mu_ewo_sig"     , "mu_Ewomc_sig"     , RooArgSet( *rv_mu_Ewomc_sig     ) ) ;
      rv_mu_ewo_sb      = new RooFormulaVar( "mu_ewo_sb"      , "mu_Ewomc_sb"      , RooArgSet( *rv_mu_Ewomc_sb      ) ) ;
      rv_mu_ewo_sig_ldp = new RooFormulaVar( "mu_ewo_sig_ldp" , "mu_Ewomc_sig_ldp" , RooArgSet( *rv_mu_Ewomc_sig_ldp ) ) ;
      rv_mu_ewo_sb_ldp  = new RooFormulaVar( "mu_ewo_sb_ldp"  , "mu_Ewomc_sb_ldp"  , RooArgSet( *rv_mu_Ewomc_sb_ldp  ) ) ;




    //-- Parametric relations between correlated efficiency scale factors.



      rv_eff_sf_sig = new RooFormulaVar( "eff_sf_sig",
                                         "mean_eff_sf_sig * pow( exp( width_eff_sf_sig/mean_eff_sf_sig ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig, *rv_width_eff_sf_sig, *rv_mean_eff_sf_sig, *rv_eff_sf_prim ) ) ;

      rv_eff_sf_sb = new RooFormulaVar( "eff_sf_sb",
                                         "mean_eff_sf_sb * pow( exp( width_eff_sf_sb/mean_eff_sf_sb ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb, *rv_width_eff_sf_sb, *rv_mean_eff_sf_sb, *rv_eff_sf_prim ) ) ;

      rv_eff_sf_sig_sl = new RooFormulaVar( "eff_sf_sig_sl",
                                         "mean_eff_sf_sig_sl * pow( exp( width_eff_sf_sig_sl/mean_eff_sf_sig_sl ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_sl, *rv_width_eff_sf_sig_sl, *rv_mean_eff_sf_sig_sl, *rv_eff_sf_prim ) ) ;

      rv_eff_sf_sb_sl = new RooFormulaVar( "eff_sf_sb_sl",
                                         "mean_eff_sf_sb_sl * pow( exp( width_eff_sf_sb_sl/mean_eff_sf_sb_sl ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_sl, *rv_width_eff_sf_sb_sl, *rv_mean_eff_sf_sb_sl, *rv_eff_sf_prim ) ) ;

      rv_eff_sf_sig_ldp = new RooFormulaVar( "eff_sf_sig_ldp",
                                         "mean_eff_sf_sig_ldp * pow( exp( width_eff_sf_sig_ldp/mean_eff_sf_sig_ldp ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_ldp, *rv_width_eff_sf_sig_ldp, *rv_mean_eff_sf_sig_ldp, *rv_eff_sf_prim ) ) ;

      rv_eff_sf_sb_ldp = new RooFormulaVar( "eff_sf_sb_ldp",
                                         "mean_eff_sf_sb_ldp * pow( exp( width_eff_sf_sb_ldp/mean_eff_sf_sb_ldp ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_ldp, *rv_width_eff_sf_sb_ldp, *rv_mean_eff_sf_sb_ldp, *rv_eff_sf_prim ) ) ;


    //+++++++++++++ Expected counts for observables in terms of parameters ++++++++++++++++++

       printf(" --- Defining expected counts in terms of parameters.\n" ) ;


      rv_n_sig         = new RooFormulaVar( "n_sig",
                                     "mu_ttwj_sig + mu_qcd_sig + mu_znn_sig + eff_sf_sig*( mu_ewo_sig + mu_susy_sig)",
                                     RooArgSet( *rv_mu_ttwj_sig, *rv_mu_qcd_sig, *rv_mu_znn_sig, *rv_eff_sf_sig, *rv_mu_ewo_sig, *rv_mu_susy_sig ) ) ;

      rv_n_sb          = new RooFormulaVar( "n_sb",
                                     "mu_ttwj_sb  + mu_qcd_sb  + mu_znn_sb  + eff_sf_sb*( mu_ewo_sb  + mu_susy_sb )",
                                     RooArgSet( *rv_mu_ttwj_sb , *rv_mu_qcd_sb , *rv_mu_znn_sb , *rv_eff_sf_sb, *rv_mu_ewo_sb , *rv_mu_susy_sb  ) ) ;

      rv_n_sig_ldp     = new RooFormulaVar( "n_sig_ldp",
                                     "mu_qcd_sig_ldp + eff_sf_sig_ldp*( sf_mc * (mu_ttwj_sig_ldp + mu_znn_sig_ldp + mu_ewo_sig_ldp) + mu_susy_sig_ldp)",
                                     RooArgSet( *rv_mu_qcd_sig_ldp, *rv_eff_sf_sig_ldp, fv_sf_mc, *rv_mu_ttwj_sig_ldp, *rv_mu_znn_sig_ldp, *rv_mu_ewo_sig_ldp, *rv_mu_susy_sig_ldp ) ) ;

      rv_n_sb_ldp      = new RooFormulaVar( "n_sb_ldp",
                                     "mu_qcd_sb_ldp + eff_sf_sb_ldp*( sf_mc * (mu_ttwj_sb_ldp + mu_znn_sb_ldp + mu_ewo_sb_ldp) + mu_susy_sb_ldp)",
                                     RooArgSet( *rv_mu_qcd_sb_ldp, *rv_eff_sf_sb_ldp, fv_sf_mc, *rv_mu_ttwj_sb_ldp, *rv_mu_znn_sb_ldp, *rv_mu_ewo_sb_ldp, *rv_mu_susy_sb_ldp ) ) ;

      rv_n_sig_sl      = new RooFormulaVar( "n_sig_sl",
                                     "mu_ttwj_sig_sl + eff_sf_sig_sl*mu_susy_sig_sl",
                                     RooArgSet( *rv_mu_ttwj_sig_sl, *rv_eff_sf_sig_sl, *rv_mu_susy_sig_sl ) ) ;

      rv_n_sb_sl       = new RooFormulaVar( "n_sb_sl",
                                     "mu_ttwj_sb_sl + eff_sf_sb_sl*mu_susy_sb_sl",
                                     RooArgSet( *rv_mu_ttwj_sb_sl, *rv_eff_sf_sb_sl, *rv_mu_susy_sb_sl ) ) ;

      rv_n_lsb_0b      = new RooFormulaVar( "n_lsb_0b",
                                     "mu_qcd_lsb_0b",
                                     RooArgSet( *rv_mu_qcd_lsb_0b ) ) ;

      rv_n_lsb_0b_ldp  = new RooFormulaVar( "n_lsb_0b_ldp",
                                     "mu_qcd_lsb_0b_ldp",
                                     RooArgSet( *rv_mu_qcd_lsb_0b_ldp ) ) ;

      if ( znnModel == 1 ) {

         rv_n_sig_ee      = new RooFormulaVar( "n_sig_ee",
                                        "mu_zee_sig_ee / fsig_ee",
                                        RooArgSet( *rv_mu_zee_sig_ee, fv_fsig_ee ) ) ;

         rv_n_sb_ee       = new RooFormulaVar( "n_sb_ee",
                                        "mu_zee_sb_ee / fsig_ee",
                                        RooArgSet( *rv_mu_zee_sb_ee, fv_fsig_ee ) ) ;

         rv_n_sig_mm      = new RooFormulaVar( "n_sig_mm",
                                        "mu_zmm_sig_mm / fsig_mm",
                                        RooArgSet( *rv_mu_zmm_sig_mm, fv_fsig_mm ) ) ;

         rv_n_sb_mm       = new RooFormulaVar( "n_sb_mm",
                                        "mu_zmm_sb_mm / fsig_mm",
                                        RooArgSet( *rv_mu_zmm_sb_mm, fv_fsig_mm ) ) ;

      } else if ( znnModel == 2 ) {

         rv_n_sigsb_ee      = new RooFormulaVar( "n_sigsb_ee",
                                        "mu_zee_sigsb_ee / fsig_ee",
                                        RooArgSet( *rv_mu_zee_sigsb_ee, fv_fsig_ee ) ) ;

         rv_n_sigsb_mm      = new RooFormulaVar( "n_sigsb_mm",
                                        "mu_zmm_sigsb_mm / fsig_mm",
                                        RooArgSet( *rv_mu_zmm_sigsb_mm, fv_fsig_mm ) ) ;

      }

   //++++++++++++ PDFs for the likelihood +++++++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining PDFs of the likelihood.\n" ) ;

      pdf_Nsig        = new RooPoisson( "pdf_Nsig"        , "Nsig Poisson PDF"        , *rv_Nsig        , *rv_n_sig ) ;
      pdf_Nsb         = new RooPoisson( "pdf_Nsb"         , "Nsb Poisson PDF"         , *rv_Nsb         , *rv_n_sb ) ;
      pdf_Nsig_ldp    = new RooPoisson( "pdf_Nsig_ldp"    , "Nsig_ldp Poisson PDF"    , *rv_Nsig_ldp    , *rv_n_sig_ldp ) ;
      pdf_Nsb_ldp     = new RooPoisson( "pdf_Nsb_ldp"     , "Nsb_ldp Poisson PDF"     , *rv_Nsb_ldp     , *rv_n_sb_ldp ) ;
      pdf_Nsig_sl     = new RooPoisson( "pdf_Nsig_sl"     , "Nsig_sl Poisson PDF"     , *rv_Nsig_sl     , *rv_n_sig_sl ) ;
      pdf_Nsb_sl      = new RooPoisson( "pdf_Nsb_sl"      , "Nsb_sl Poisson PDF"      , *rv_Nsb_sl      , *rv_n_sb_sl ) ;
      pdf_Nlsb_0b     = new RooPoisson( "pdf_Nlsb_0b"     , "Nlsb_0b Poisson PDF"     , *rv_Nlsb_0b     , *rv_n_lsb_0b ) ;
      pdf_Nlsb_0b_ldp = new RooPoisson( "pdf_Nlsb_0b_ldp" , "Nlsb_0b_ldp Poisson PDF" , *rv_Nlsb_0b_ldp , *rv_n_lsb_0b_ldp ) ;

      if ( znnModel == 1 ) {
         pdf_Nsig_ee     = new RooPoisson( "pdf_Nsig_ee"     , "Nsig_ee Poisson PDF"     , *rv_Nsig_ee     , *rv_n_sig_ee ) ;
         pdf_Nsb_ee      = new RooPoisson( "pdf_Nsb_ee"      , "Nsb_ee Poisson PDF"      , *rv_Nsb_ee      , *rv_n_sb_ee ) ;
         pdf_Nsig_mm     = new RooPoisson( "pdf_Nsig_mm"     , "Nsig_mm Poisson PDF"     , *rv_Nsig_mm     , *rv_n_sig_mm ) ;
         pdf_Nsb_mm      = new RooPoisson( "pdf_Nsb_mm"      , "Nsb_mm Poisson PDF"      , *rv_Nsb_mm      , *rv_n_sb_mm ) ;
      } else if ( znnModel == 2 ) {
         pdf_Nsigsb_ee     = new RooPoisson( "pdf_Nsigsb_ee"     , "Nsigsb_ee Poisson PDF"     , *rv_Nsigsb_ee     , *rv_n_sigsb_ee ) ;
         pdf_Nsigsb_mm     = new RooPoisson( "pdf_Nsigsb_mm"     , "Nsigsb_mm Poisson PDF"     , *rv_Nsigsb_mm     , *rv_n_sigsb_mm ) ;
      }


      {
         RooArgSet pdflist ;
         pdflist.add( *pdf_Nsig        ) ;
         pdflist.add( *pdf_Nsb         ) ;
         pdflist.add( *pdf_Nsig_ldp    ) ;
         pdflist.add( *pdf_Nsb_ldp     ) ;
         pdflist.add( *pdf_Nsig_sl     ) ;
         pdflist.add( *pdf_Nsb_sl      ) ;
         pdflist.add( *pdf_Nlsb_0b     ) ;
         pdflist.add( *pdf_Nlsb_0b_ldp ) ;
         if ( znnModel == 1 ) {
            pdflist.add( *pdf_Nsig_ee     ) ;
            pdflist.add( *pdf_Nsb_ee      ) ;
            pdflist.add( *pdf_Nsig_mm     ) ;
            pdflist.add( *pdf_Nsb_mm      ) ;
         } else if ( znnModel == 2 ) {
            pdflist.add( *pdf_Nsigsb_ee   ) ;
            pdflist.add( *pdf_Nsigsb_mm   ) ;
         }

	 pdflist.add(allNuisancePdfs);

         likelihood = new RooProdPdf( "likelihood", "ra2b likelihood", pdflist ) ;
      }


     //---- Define the list of observables.

       observedParametersList.add( *rv_Nsig        ) ;
       observedParametersList.add( *rv_Nsb         ) ;
       observedParametersList.add( *rv_Nsig_sl     ) ;
       observedParametersList.add( *rv_Nsb_sl      ) ;
       observedParametersList.add( *rv_Nsig_ldp    ) ;
       observedParametersList.add( *rv_Nsb_ldp     ) ;
       observedParametersList.add( *rv_Nlsb_0b     ) ;
       observedParametersList.add( *rv_Nlsb_0b_ldp ) ;
       if ( znnModel == 1 ) {
          observedParametersList.add( *rv_Nsb_ee      ) ;
          observedParametersList.add( *rv_Nsig_ee     ) ;
          observedParametersList.add( *rv_Nsb_mm      ) ;
          observedParametersList.add( *rv_Nsig_mm     ) ;
       } else if ( znnModel == 2 ) {
          observedParametersList.add( *rv_Nsigsb_ee      ) ;
          observedParametersList.add( *rv_Nsigsb_mm      ) ;
       }



       dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
                                      observedParametersList ) ;
       dsObserved->add( observedParametersList ) ;


       RooWorkspace workspace ("ws") ;

       workspace.import(*dsObserved);

       // parameters of interest
       RooArgSet poi(*rv_mu_susy_sig, "poi");
       // flat prior for POI
       RooUniform signal_prior ("signal_prior","signal_prior",*rv_mu_susy_sig);


       // signal+background model
       ModelConfig sbModel ("SbModel");
       sbModel.SetWorkspace(workspace);
       sbModel.SetPdf(*likelihood);
       sbModel.SetParametersOfInterest(poi);
       sbModel.SetPriorPdf(signal_prior);
       sbModel.SetNuisanceParameters(allNuisances);
       sbModel.SetObservables(observedParametersList);
       sbModel.SetGlobalObservables(globalObservables);

       // find global maximum with the signal+background model
       // with conditional MLEs for nuisance parameters
       // and save the parameter point snapshot in the Workspace
       //  - safer to keep a default name because some RooStats calculators
       //    will anticipate it
       RooAbsReal * pNll = sbModel.GetPdf()->createNLL(*dsObserved);
       RooAbsReal * pProfile = pNll->createProfile(RooArgSet());
       pProfile->getVal(); // this will do fit and set POI and nuisance parameters to fitted values
       RooArgSet * pPoiAndNuisance = new RooArgSet();
       pPoiAndNuisance->add(*sbModel.GetParametersOfInterest());
       if(sbModel.GetNuisanceParameters()) pPoiAndNuisance->add(*sbModel.GetNuisanceParameters());
       cout << "\nWill save these parameter points that correspond to the fit to data" << endl;
       pPoiAndNuisance->Print("v");
       sbModel.SetSnapshot(*pPoiAndNuisance);
       workspace.import (sbModel);

       delete pProfile;
       delete pNll;
       delete pPoiAndNuisance;


       // background-only model
       // use the same PDF as s+b, with xsec=0
       // POI value under the background hypothesis
       ModelConfig bModel (*(RooStats::ModelConfig *)workspace.obj("SbModel"));
       bModel.SetName("BModel");
       bModel.SetWorkspace(workspace);

       // Find a parameter point for generating pseudo-data
       // with the background-only data.
       // Save the parameter point snapshot in the Workspace
       pNll = bModel.GetPdf()->createNLL(*dsObserved);
       // bug discovered by Fedor on Sep 6th 2011:
       //pProfile = pNll->createProfile(poi);
       //((RooRealVar *)poi.first())->setVal(0.); // set signal = 0
       pProfile = pNll->createProfile(*bModel.GetParametersOfInterest());
       ((RooRealVar *)(bModel.GetParametersOfInterest()->first()))->setVal(0.); // set signal = 0
       pProfile->getVal(); // this will do fit and set nuisance parameters to profiled values
       pPoiAndNuisance = new RooArgSet();
       pPoiAndNuisance->add(*bModel.GetParametersOfInterest());
       if(bModel.GetNuisanceParameters()) pPoiAndNuisance->add(*bModel.GetNuisanceParameters());
       cout << "\nShould use these parameter points to generate pseudo data for bkg only" << endl;
       pPoiAndNuisance->Print("v");
       bModel.SetSnapshot(*pPoiAndNuisance);
       workspace.import (bModel);

       delete pProfile;
       delete pNll;
       delete pPoiAndNuisance;


       workspace.Print() ;
       workspace.writeToFile("ws.root");

       return true ;


    } // initialize.


   //===================================================================================================================================


    bool ra2bRoostatsClass7::setSusyScanPoint( const char* inputScanFile, double m0, double m12, bool isT1bbbb, double t1bbbbXsec ) {


       //--- Aug 15, 2011: updated to new format for AN, v3.


       printf("\n\n Opening SUSY scan input file : %s\n", inputScanFile ) ;

       FILE* infp ;
       if ( (infp=fopen( inputScanFile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", inputScanFile ) ;
          return false ;
       }

       double deltaM0(0.) ;
       double deltaM12(0.) ;

       if ( !isT1bbbb ) {
          deltaM0 = 20 ;
          deltaM12 = 20 ;
       } else {
          deltaM0 = 25 ;
          deltaM12 = 25 ;
       }

       bool found(false) ;

       //--- Loop over the scan points.
       while ( !feof( infp ) ) {

          float pointM0 ;
          float pointM12 ;



   //+++ ORL: Aug 14, 2011 ++++++++++++++++++++++++++
          float n_sig_raw ;
          float n_sb_raw ;
          float n_sig_sl_raw ;
          float n_sb_sl_raw ;
          float n_sig_ldp_raw ;
          float n_sb_ldp_raw ;

          float n_sig_correction ;
          float n_sb_correction ;
          float n_sig_sl_correction ;
          float n_sb_sl_correction ;
          float n_sig_ldp_correction ;
          float n_sb_ldp_correction ;

          float n_sig_error ;
          float n_sb_error ;
          float n_sig_sl_error ;
          float n_sb_sl_error ;
          float n_sig_ldp_error ;
          float n_sb_ldp_error ;

          int nGen ;

          fscanf( infp, "%f %f %d  %f %f %f %f %f %f  %f %f %f %f %f %f %f %f %f %f %f %f",
            &pointM0, &pointM12, &nGen,
            &n_sig_raw, &n_sb_raw, &n_sig_sl_raw, &n_sb_sl_raw, &n_sig_ldp_raw, &n_sb_ldp_raw,
            &n_sig_correction, &n_sb_correction, &n_sig_sl_correction, &n_sb_sl_correction, &n_sig_ldp_correction, &n_sb_ldp_correction,
            &n_sig_error, &n_sb_error, &n_sig_sl_error, &n_sb_sl_error, &n_sig_ldp_error, &n_sb_ldp_error ) ;

          if ( feof(infp) ) break ;
          if ( n_sig_raw < 0.00001 ) continue ;
   //--If you are asking for it, I'll assume it's good.  Josh is using 0 for ngen dummy in LM9.
   /////  if ( nGen != 10000 ) continue ; // get rid of bad scan points.

          if (    fabs( pointM0 - m0 ) <= deltaM0/2.
               && fabs( pointM12 - m12 ) <= deltaM12/2. ) {

              double nGenPerPoint = 10000 ; // for t1bbbb

             printf("\n\n") ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_error     = %6.1f %%\n", n_sig_error     ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_error      = %6.1f %%\n", n_sb_error      ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_sl_error  = %6.1f %%\n", n_sig_sl_error  ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_sl_error   = %6.1f %%\n", n_sb_sl_error   ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_ldp_error = %6.1f %%\n", n_sig_ldp_error ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_ldp_error  = %6.1f %%\n", n_sb_ldp_error  ) ;
             printf("\n\n") ;

      //--Include the stat error on the efficiency for t1bbbb.
             if ( isT1bbbb ) {

              //-- absolute raw eff
                 float n_sig_raw_eff     = n_sig_raw     / nGenPerPoint ;
                 float n_sb_raw_eff      = n_sb_raw      / nGenPerPoint ;
                 float n_sig_sl_raw_eff  = n_sig_sl_raw  / nGenPerPoint ;
                 float n_sb_sl_raw_eff   = n_sb_sl_raw   / nGenPerPoint ;
                 float n_sig_ldp_raw_eff = n_sig_ldp_raw / nGenPerPoint ;
                 float n_sb_ldp_raw_eff  = n_sb_ldp_raw  / nGenPerPoint ;


               //-- absolute stat err.
                 float n_sig_stat_error     =  sqrt(  n_sig_raw_eff    * ( 1.0 - n_sig_raw_eff     ) / nGenPerPoint ) ;
                 float n_sb_stat_error      =  sqrt(  n_sb_raw_eff     * ( 1.0 - n_sb_raw_eff      ) / nGenPerPoint ) ;
                 float n_sig_sl_stat_error  =  sqrt(  n_sig_sl_raw_eff * ( 1.0 - n_sig_sl_raw_eff  ) / nGenPerPoint ) ;
                 float n_sb_sl_stat_error   =  sqrt(  n_sb_sl_raw_eff  * ( 1.0 - n_sb_sl_raw_eff   ) / nGenPerPoint ) ;
                 float n_sig_ldp_stat_error =  sqrt(  n_sig_ldp_raw_eff* ( 1.0 - n_sig_ldp_raw_eff ) / nGenPerPoint ) ;
                 float n_sb_ldp_stat_error  =  sqrt(  n_sb_ldp_raw_eff * ( 1.0 - n_sb_ldp_raw_eff  ) / nGenPerPoint ) ;

               //-- relative stat err in percent.
                 if ( n_sig_raw_eff     > 0 ) { n_sig_stat_error     = 100.* n_sig_stat_error     / n_sig_raw_eff     ; } else { n_sig_stat_error     = 0. ; }
                 if ( n_sb_raw_eff      > 0 ) { n_sb_stat_error      = 100.* n_sb_stat_error      / n_sb_raw_eff      ; } else { n_sb_stat_error      = 0. ; }
                 if ( n_sig_sl_raw_eff  > 0 ) { n_sig_sl_stat_error  = 100.* n_sig_sl_stat_error  / n_sig_sl_raw_eff  ; } else { n_sig_sl_stat_error  = 0. ; }
                 if ( n_sb_sl_raw_eff   > 0 ) { n_sb_sl_stat_error   = 100.* n_sb_sl_stat_error   / n_sb_sl_raw_eff   ; } else { n_sb_sl_stat_error   = 0. ; }
                 if ( n_sig_ldp_raw_eff > 0 ) { n_sig_ldp_stat_error = 100.* n_sig_ldp_stat_error / n_sig_ldp_raw_eff ; } else { n_sig_ldp_stat_error = 0. ; }
                 if ( n_sb_ldp_raw_eff  > 0 ) { n_sb_ldp_stat_error  = 100.* n_sb_ldp_stat_error  / n_sb_ldp_raw_eff  ; } else { n_sb_ldp_stat_error  = 0. ; }

                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_stat_error     = %6.1f %%\n", n_sig_stat_error     ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_stat_error      = %6.1f %%\n", n_sb_stat_error      ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_sl_stat_error  = %6.1f %%\n", n_sig_sl_stat_error  ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_sl_stat_error   = %6.1f %%\n", n_sb_sl_stat_error   ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_ldp_stat_error = %6.1f %%\n", n_sig_ldp_stat_error ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_ldp_stat_error  = %6.1f %%\n", n_sb_ldp_stat_error  ) ;

               //-- total err in percent.
                 n_sig_error     = sqrt( pow( n_sig_error    , 2) + pow( n_sig_stat_error    , 2) ) ;
                 n_sb_error      = sqrt( pow( n_sb_error     , 2) + pow( n_sb_stat_error     , 2) ) ;
                 n_sig_sl_error  = sqrt( pow( n_sig_sl_error , 2) + pow( n_sig_sl_stat_error , 2) ) ;
                 n_sb_sl_error   = sqrt( pow( n_sb_sl_error  , 2) + pow( n_sb_sl_stat_error  , 2) ) ;
                 n_sig_ldp_error = sqrt( pow( n_sig_ldp_error, 2) + pow( n_sig_ldp_stat_error, 2) ) ;
                 n_sb_ldp_error  = sqrt( pow( n_sb_ldp_error , 2) + pow( n_sb_ldp_stat_error , 2) ) ;

                 printf("\n\n") ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_error     = %6.1f %%\n", n_sig_error     ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_error      = %6.1f %%\n", n_sb_error      ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_sl_error  = %6.1f %%\n", n_sig_sl_error  ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_sl_error   = %6.1f %%\n", n_sb_sl_error   ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_ldp_error = %6.1f %%\n", n_sig_ldp_error ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_ldp_error  = %6.1f %%\n", n_sb_ldp_error  ) ;
                 printf("\n\n") ;

             }

       //--- Not needed with log-normal
        ///  //-- enforce a maximum efficiency uncertainty (to avoid negative scale factors).
        ///  if ( n_sig_error     > 35. ) { n_sig_error     = 35. ; }
        ///  if ( n_sb_error      > 35. ) { n_sb_error      = 35. ; }
        ///  if ( n_sig_sl_error  > 35. ) { n_sig_sl_error  = 35. ; }
        ///  if ( n_sb_sl_error   > 35. ) { n_sb_sl_error   = 35. ; }
        ///  if ( n_sig_ldp_error > 35. ) { n_sig_ldp_error = 35. ; }
        ///  if ( n_sb_ldp_error  > 35. ) { n_sb_ldp_error  = 35. ; }

      //++++++++++++++++++++++++++++++++++++++++++++++++



             double setVal_n_sig(0.) ;
             double setVal_n_sb(0.) ;
             double setVal_n_sig_sl(0.) ;
             double setVal_n_sb_sl(0.) ;
             double setVal_n_sig_ldp(0.) ;
             double setVal_n_sb_ldp(0.) ;


             if ( !isT1bbbb ) {
                //-- tanb40
                setVal_n_sig     = n_sig_raw     * n_sig_correction     ;
                setVal_n_sb      = n_sb_raw      * n_sb_correction      ;
                setVal_n_sig_sl  = n_sig_sl_raw  * n_sig_sl_correction  ;
                setVal_n_sb_sl   = n_sb_sl_raw   * n_sb_sl_correction   ;
                setVal_n_sig_ldp = n_sig_ldp_raw * n_sig_ldp_correction ;
                setVal_n_sb_ldp  = n_sb_ldp_raw  * n_sb_ldp_correction  ;
             } else {
                //-- t1bbbb
                setVal_n_sig     = DataLumi * t1bbbbXsec * (( n_sig_raw     * n_sig_correction     )/ nGenPerPoint ) ;
                setVal_n_sb      = DataLumi * t1bbbbXsec * (( n_sb_raw      * n_sb_correction      )/ nGenPerPoint ) ;
                setVal_n_sig_sl  = DataLumi * t1bbbbXsec * (( n_sig_sl_raw  * n_sig_sl_correction  )/ nGenPerPoint ) ;
                setVal_n_sb_sl   = DataLumi * t1bbbbXsec * (( n_sb_sl_raw   * n_sb_sl_correction   )/ nGenPerPoint ) ;
                setVal_n_sig_ldp = DataLumi * t1bbbbXsec * (( n_sig_ldp_raw * n_sig_ldp_correction )/ nGenPerPoint ) ;
                setVal_n_sb_ldp  = DataLumi * t1bbbbXsec * (( n_sb_ldp_raw  * n_sb_ldp_correction  )/ nGenPerPoint ) ;
             }

             rv_mu_susymc_sig       -> setVal( setVal_n_sig      ) ;
             rv_mu_susymc_sb        -> setVal( setVal_n_sb       ) ;
             rv_mu_susymc_sig_sl    -> setVal( setVal_n_sig_sl   ) ;
             rv_mu_susymc_sb_sl     -> setVal( setVal_n_sb_sl    ) ;
             rv_mu_susymc_sig_ldp   -> setVal( setVal_n_sig_ldp  ) ;
             rv_mu_susymc_sb_ldp    -> setVal( setVal_n_sb_ldp   ) ;

             rv_width_eff_sf_sig     -> setVal( n_sig_error     / 100. ) ;
             rv_width_eff_sf_sb      -> setVal( n_sb_error      / 100. ) ;
             rv_width_eff_sf_sig_sl  -> setVal( n_sig_sl_error  / 100. ) ;
             rv_width_eff_sf_sb_sl   -> setVal( n_sb_sl_error   / 100. ) ;
             rv_width_eff_sf_sig_ldp -> setVal( n_sig_ldp_error / 100. ) ;
             rv_width_eff_sf_sb_ldp  -> setVal( n_sb_ldp_error  / 100. ) ;



             if ( !isT1bbbb ) {
                printf("\n\n Found point m0 = %4.0f,  m1/2 = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, setVal_n_sig ) ;
             } else {
                printf("\n\n Found point mGluino = %4.0f,  mLSP = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, setVal_n_sig ) ;
             }


             printf("\n\n") ;
             printf(" Setting susy N_sig     to  %7.1f\n", setVal_n_sig       ) ;
             printf(" Setting susy N_sb      to  %7.1f\n", setVal_n_sb        ) ;
             printf(" Setting susy N_sig_sl  to  %7.1f\n", setVal_n_sig_sl    ) ;
             printf(" Setting susy N_sb_sl   to  %7.1f\n", setVal_n_sb_sl     ) ;
             printf(" Setting susy N_sig_ldp to  %7.1f\n", setVal_n_sig_ldp   ) ;
             printf(" Setting susy N_sb_ldp  to  %7.1f\n", setVal_n_sb_ldp    ) ;
             printf("\n\n") ;

             found = true ;

             break ;

          } // point match?

       } // not eof ?

       fclose( infp ) ;

       if ( found ) {
          return true ;
       } else {
          printf("\n\n *** Point not found in scan.\n\n" ) ;
          return false ;
       }

    } // setSusyScanPoint.
