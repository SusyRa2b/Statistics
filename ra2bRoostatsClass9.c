//
//   Owen Long, UCR
//   Harrison Prosper, FSU
//   Sezen Sekmen, FSU
//
//


#include "ra2bRoostatsClass9.h"

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

//////////#include "LikelihoodIntervalPlot.cxx"

  using namespace RooFit ;
  using namespace RooStats ;


  //=====================================================================================================


   ra2bRoostatsClass9::ra2bRoostatsClass9( bool ArgUseSigTtwjVar, bool ArgUseLdpVars ) {

      gStyle->SetOptStat(0) ;

      useSigTtwjVar = ArgUseSigTtwjVar ;
      useLdpVars = ArgUseLdpVars ;


     //--- Tell RooFit to shut up about anything less important than an ERROR.
      RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR) ;
  /// RooMsgService::instance().addStream(DEBUG,Topic(Tracing),OutputFile("debug.log")) ;
      printf("\n\n ==== RooFit output configuration ===============\n") ;
      RooMsgService::instance().Print("v") ;
      printf("\n\n ================================================\n") ;

      varsAtFitVals = false ;
      initialized = false ;




   }





  //===================================================================

   ra2bRoostatsClass9::~ra2bRoostatsClass9() { }



  //===================================================================

    bool ra2bRoostatsClass9::initialize( const char* infile ,
                                         const char* inputScanFile,
                                         double m0, double m12, bool isT1bbbb, double t1bbbbXsec,
                                         const char* inputSusy_deff_dbtageff_file
                                       ) {


       printf( "\n\n Opening input file : %s\n\n", infile ) ;

       FILE* infp ;
       if ( (infp=fopen( infile,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", infile ) ;
          return false ;
       }

       sprintf( initializeFile, "%s", infile ) ;


       int    Nsig_1b                  ;
       int    Nsig_2b                  ;
       int    Nsig_3b                  ;
       int    Nsb_1b                   ;
       int    Nsb_2b                   ;
       int    Nsb_3b                   ;
       int    Nsig_sl_1b               ;
       int    Nsig_sl_2b               ;
       int    Nsig_sl_3b               ;
       int    Nsb_sl_1b                ;
       int    Nsb_sl_2b                ;
       int    Nsb_sl_3b                ;
       int    Nsig_ldp_1b              ;
       int    Nsig_ldp_2b              ;
       int    Nsig_ldp_3b              ;
       int    Nsb_ldp_1b               ;
       int    Nsb_ldp_2b               ;
       int    Nsb_ldp_3b               ;
       int    Nsb_ee                  ;
       int    Nsig_ee                 ;
       int    Nsb_mm                  ;
       int    Nsig_mm                 ;
       float  Nttbarsingletopzjetsmc_sig_ldp_1b      ;
       float  Nttbarsingletopzjetsmc_sig_ldp_2b      ;
       float  Nttbarsingletopzjetsmc_sig_ldp_3b      ;
       float  Nttbarsingletopzjetsmc_sb_ldp_1b       ;
       float  Nttbarsingletopzjetsmc_sb_ldp_2b       ;
       float  Nttbarsingletopzjetsmc_sb_ldp_3b       ;
       float  NWJmc_sig_ldp_1b         ;
       float  NWJmc_sig_ldp_2b         ;
       float  NWJmc_sig_ldp_3b         ;
       float  NWJmc_sb_ldp_1b          ;
       float  NWJmc_sb_ldp_2b          ;
       float  NWJmc_sb_ldp_3b          ;
       float  NZnnmc_sig_ldp_1b        ;
       float  NZnnmc_sig_ldp_2b        ;
       float  NZnnmc_sig_ldp_3b        ;
       float  NZnnmc_sb_ldp_1b         ;
       float  NZnnmc_sb_ldp_2b         ;
       float  NZnnmc_sb_ldp_3b         ;


       float  sf_mc            ;
       float  sf_mc_err        ;
       float  sf_qcd_sb_1b        ;
       float  sf_qcd_sb_1b_err    ;
       float  sf_qcd_sb_2b        ;
       float  sf_qcd_sb_2b_err    ;
       float  sf_qcd_sb_3b        ;
       float  sf_qcd_sb_3b_err    ;
       float  sf_qcd_sig_1b       ;
       float  sf_qcd_sig_1b_err   ;
       float  sf_qcd_sig_2b       ;
       float  sf_qcd_sig_2b_err   ;
       float  sf_qcd_sig_3b       ;
       float  sf_qcd_sig_3b_err   ;
       float  sf_ttwj_sig_1b      ;
       float  sf_ttwj_sig_1b_err  ;
       float  sf_ttwj_sig_2b      ;
       float  sf_ttwj_sig_2b_err  ;
       float  sf_ttwj_sig_3b      ;
       float  sf_ttwj_sig_3b_err  ;
       float  sf_ttwj_sb_2b      ;
       float  sf_ttwj_sb_2b_err  ;
       float  sf_ttwj_sb_3b      ;
       float  sf_ttwj_sb_3b_err  ;
       float  sf_ee            ;
       float  sf_ee_err        ;
       float  sf_mm            ;
       float  sf_mm_err        ;

       float  Rlsb_passfail     ;
       float  Rlsb_passfail_err ;


       float  btageff_err ;






       //--- read in description line.
       printf("\n\n") ;
       char c(0) ;
       while ( c!=10  ) { c = fgetc( infp ) ; printf("%c", c ) ; }
       printf("\n\n") ;


       char label[1000] ;




       fscanf( infp, "%s %d", label, &Nsig_1b                                ) ;   printf( "Nsig_1b                                %s %d\n", label, Nsig_1b                               ) ;
       fscanf( infp, "%s %d", label, &Nsig_2b                                ) ;   printf( "Nsig_2b                                %s %d\n", label, Nsig_2b                               ) ;
       fscanf( infp, "%s %d", label, &Nsig_3b                                ) ;   printf( "Nsig_3b                                %s %d\n", label, Nsig_3b                               ) ;
       fscanf( infp, "%s %d", label, &Nsb_1b                                 ) ;   printf( "Nsb_1b                                 %s %d\n", label, Nsb_1b                                ) ;
       fscanf( infp, "%s %d", label, &Nsb_2b                                 ) ;   printf( "Nsb_2b                                 %s %d\n", label, Nsb_2b                                ) ;
       fscanf( infp, "%s %d", label, &Nsb_3b                                 ) ;   printf( "Nsb_3b                                 %s %d\n", label, Nsb_3b                                ) ;
       fscanf( infp, "%s %d", label, &Nsig_sl_1b                             ) ;   printf( "Nsig_sl_1b                             %s %d\n", label, Nsig_sl_1b                            ) ;
       fscanf( infp, "%s %d", label, &Nsig_sl_2b                             ) ;   printf( "Nsig_sl_2b                             %s %d\n", label, Nsig_sl_2b                            ) ;
       fscanf( infp, "%s %d", label, &Nsig_sl_3b                             ) ;   printf( "Nsig_sl_3b                             %s %d\n", label, Nsig_sl_3b                            ) ;
       fscanf( infp, "%s %d", label, &Nsb_sl_1b                              ) ;   printf( "Nsb_sl_1b                              %s %d\n", label, Nsb_sl_1b                             ) ;
       fscanf( infp, "%s %d", label, &Nsb_sl_2b                              ) ;   printf( "Nsb_sl_2b                              %s %d\n", label, Nsb_sl_2b                             ) ;
       fscanf( infp, "%s %d", label, &Nsb_sl_3b                              ) ;   printf( "Nsb_sl_3b                              %s %d\n", label, Nsb_sl_3b                             ) ;
       fscanf( infp, "%s %d", label, &Nsig_ldp_1b                            ) ;   printf( "Nsig_ldp_1b                            %s %d\n", label, Nsig_ldp_1b                           ) ;
       fscanf( infp, "%s %d", label, &Nsig_ldp_2b                            ) ;   printf( "Nsig_ldp_2b                            %s %d\n", label, Nsig_ldp_2b                           ) ;
       fscanf( infp, "%s %d", label, &Nsig_ldp_3b                            ) ;   printf( "Nsig_ldp_3b                            %s %d\n", label, Nsig_ldp_3b                           ) ;
       fscanf( infp, "%s %d", label, &Nsb_ldp_1b                             ) ;   printf( "Nsb_ldp_1b                             %s %d\n", label, Nsb_ldp_1b                            ) ;
       fscanf( infp, "%s %d", label, &Nsb_ldp_2b                             ) ;   printf( "Nsb_ldp_2b                             %s %d\n", label, Nsb_ldp_2b                            ) ;
       fscanf( infp, "%s %d", label, &Nsb_ldp_3b                             ) ;   printf( "Nsb_ldp_3b                             %s %d\n", label, Nsb_ldp_3b                            ) ;
       fscanf( infp, "%s %d", label, &Nsb_ee                                 ) ;   printf( "Nsb_ee                                 %s %d\n", label, Nsb_ee                                ) ;
       fscanf( infp, "%s %d", label, &Nsig_ee                                ) ;   printf( "Nsig_ee                                %s %d\n", label, Nsig_ee                               ) ;
       fscanf( infp, "%s %d", label, &Nsb_mm                                 ) ;   printf( "Nsb_mm                                 %s %d\n", label, Nsb_mm                                ) ;
       fscanf( infp, "%s %d", label, &Nsig_mm                                ) ;   printf( "Nsig_mm                                %s %d\n", label, Nsig_mm                               ) ;
       fscanf( infp, "%s %g", label, &Rlsb_passfail                          ) ;   printf( "Rlsb_passfail                          %s %g\n", label, Rlsb_passfail                         ) ;
       fscanf( infp, "%s %g", label, &Rlsb_passfail_err                      ) ;   printf( "Rlsb_passfail_err                      %s %g\n", label, Rlsb_passfail_err                     ) ;
       fscanf( infp, "%s %g", label, &Nttbarsingletopzjetsmc_sig_ldp_1b      ) ;   printf( "Nttbarsingletopzjetsmc_sig_ldp_1b      %s %g\n", label, Nttbarsingletopzjetsmc_sig_ldp_1b     ) ;
       fscanf( infp, "%s %g", label, &Nttbarsingletopzjetsmc_sig_ldp_2b      ) ;   printf( "Nttbarsingletopzjetsmc_sig_ldp_2b      %s %g\n", label, Nttbarsingletopzjetsmc_sig_ldp_2b     ) ;
       fscanf( infp, "%s %g", label, &Nttbarsingletopzjetsmc_sig_ldp_3b      ) ;   printf( "Nttbarsingletopzjetsmc_sig_ldp_3b      %s %g\n", label, Nttbarsingletopzjetsmc_sig_ldp_3b     ) ;
       fscanf( infp, "%s %g", label, &Nttbarsingletopzjetsmc_sb_ldp_1b       ) ;   printf( "Nttbarsingletopzjetsmc_sb_ldp_1b       %s %g\n", label, Nttbarsingletopzjetsmc_sb_ldp_1b      ) ;
       fscanf( infp, "%s %g", label, &Nttbarsingletopzjetsmc_sb_ldp_2b       ) ;   printf( "Nttbarsingletopzjetsmc_sb_ldp_2b       %s %g\n", label, Nttbarsingletopzjetsmc_sb_ldp_2b      ) ;
       fscanf( infp, "%s %g", label, &Nttbarsingletopzjetsmc_sb_ldp_3b       ) ;   printf( "Nttbarsingletopzjetsmc_sb_ldp_3b       %s %g\n", label, Nttbarsingletopzjetsmc_sb_ldp_3b      ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig_ldp_1b                       ) ;   printf( "NWJmc_sig_ldp_1b                       %s %g\n", label, NWJmc_sig_ldp_1b                      ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig_ldp_2b                       ) ;   printf( "NWJmc_sig_ldp_2b                       %s %g\n", label, NWJmc_sig_ldp_2b                      ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sig_ldp_3b                       ) ;   printf( "NWJmc_sig_ldp_3b                       %s %g\n", label, NWJmc_sig_ldp_3b                      ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb_ldp_1b                        ) ;   printf( "NWJmc_sb_ldp_1b                        %s %g\n", label, NWJmc_sb_ldp_1b                       ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb_ldp_2b                        ) ;   printf( "NWJmc_sb_ldp_2b                        %s %g\n", label, NWJmc_sb_ldp_2b                       ) ;
       fscanf( infp, "%s %g", label, &NWJmc_sb_ldp_3b                        ) ;   printf( "NWJmc_sb_ldp_3b                        %s %g\n", label, NWJmc_sb_ldp_3b                       ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig_ldp_1b                      ) ;   printf( "NZnnmc_sig_ldp_1b                      %s %g\n", label, NZnnmc_sig_ldp_1b                     ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig_ldp_2b                      ) ;   printf( "NZnnmc_sig_ldp_2b                      %s %g\n", label, NZnnmc_sig_ldp_2b                     ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sig_ldp_3b                      ) ;   printf( "NZnnmc_sig_ldp_3b                      %s %g\n", label, NZnnmc_sig_ldp_3b                     ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb_ldp_1b                       ) ;   printf( "NZnnmc_sb_ldp_1b                       %s %g\n", label, NZnnmc_sb_ldp_1b                      ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb_ldp_2b                       ) ;   printf( "NZnnmc_sb_ldp_2b                       %s %g\n", label, NZnnmc_sb_ldp_2b                      ) ;
       fscanf( infp, "%s %g", label, &NZnnmc_sb_ldp_3b                       ) ;   printf( "NZnnmc_sb_ldp_3b                       %s %g\n", label, NZnnmc_sb_ldp_3b                      ) ;
       fscanf( infp, "%s %g", label, &DataLumi                               ) ;   printf( "DataLumi                               %s %g\n", label, DataLumi                              ) ;
       fscanf( infp, "%s %g", label, &acc_ee_sig_mean                        ) ;   printf( "acc_ee_sig_mean                        %s %g\n", label, acc_ee_sig_mean                       ) ;
       fscanf( infp, "%s %g", label, &acc_ee_sig_err                         ) ;   printf( "acc_ee_sig_err                         %s %g\n", label, acc_ee_sig_err                        ) ;
       fscanf( infp, "%s %g", label, &acc_ee_sb_mean                         ) ;   printf( "acc_ee_sb_mean                         %s %g\n", label, acc_ee_sb_mean                        ) ;
       fscanf( infp, "%s %g", label, &acc_ee_sb_err                          ) ;   printf( "acc_ee_sb_err                          %s %g\n", label, acc_ee_sb_err                         ) ;
       fscanf( infp, "%s %g", label, &acc_mm_sig_mean                        ) ;   printf( "acc_mm_sig_mean                        %s %g\n", label, acc_mm_sig_mean                       ) ;
       fscanf( infp, "%s %g", label, &acc_mm_sig_err                         ) ;   printf( "acc_mm_sig_err                         %s %g\n", label, acc_mm_sig_err                        ) ;
       fscanf( infp, "%s %g", label, &acc_mm_sb_mean                         ) ;   printf( "acc_mm_sb_mean                         %s %g\n", label, acc_mm_sb_mean                        ) ;
       fscanf( infp, "%s %g", label, &acc_mm_sb_err                          ) ;   printf( "acc_mm_sb_err                          %s %g\n", label, acc_mm_sb_err                         ) ;
       fscanf( infp, "%s %g", label, &eff_ee_mean                            ) ;   printf( "eff_ee_mean                            %s %g\n", label, eff_ee_mean                           ) ;
       fscanf( infp, "%s %g", label, &eff_ee_err                             ) ;   printf( "eff_ee_err                             %s %g\n", label, eff_ee_err                            ) ;
       fscanf( infp, "%s %g", label, &eff_mm_mean                            ) ;   printf( "eff_mm_mean                            %s %g\n", label, eff_mm_mean                           ) ;
       fscanf( infp, "%s %g", label, &eff_mm_err                             ) ;   printf( "eff_mm_err                             %s %g\n", label, eff_mm_err                            ) ;
       fscanf( infp, "%s %g", label, &Ztoll_lumi                             ) ;   printf( "Ztoll_lumi                             %s %g\n", label, Ztoll_lumi                            ) ;
       fscanf( infp, "%s %g", label, &knn_sig_1b_mean                        ) ;   printf( "knn_sig_1b_mean                        %s %g\n", label, knn_sig_1b_mean                       ) ;
       fscanf( infp, "%s %g", label, &knn_sig_1b_err                         ) ;   printf( "knn_sig_1b_err                         %s %g\n", label, knn_sig_1b_err                        ) ;
       fscanf( infp, "%s %g", label, &knn_sig_2b_mean                        ) ;   printf( "knn_sig_2b_mean                        %s %g\n", label, knn_sig_2b_mean                       ) ;
       fscanf( infp, "%s %g", label, &knn_sig_2b_err                         ) ;   printf( "knn_sig_2b_err                         %s %g\n", label, knn_sig_2b_err                        ) ;
       fscanf( infp, "%s %g", label, &knn_sig_3b_mean                        ) ;   printf( "knn_sig_3b_mean                        %s %g\n", label, knn_sig_3b_mean                       ) ;
       fscanf( infp, "%s %g", label, &knn_sig_3b_err                         ) ;   printf( "knn_sig_3b_err                         %s %g\n", label, knn_sig_3b_err                        ) ;
       fscanf( infp, "%s %g", label, &knn_sb_1b_mean                         ) ;   printf( "knn_sb_1b_mean                         %s %g\n", label, knn_sb_1b_mean                        ) ;
       fscanf( infp, "%s %g", label, &knn_sb_1b_err                          ) ;   printf( "knn_sb_1b_err                          %s %g\n", label, knn_sb_1b_err                         ) ;
       fscanf( infp, "%s %g", label, &knn_sb_2b_mean                         ) ;   printf( "knn_sb_2b_mean                         %s %g\n", label, knn_sb_2b_mean                        ) ;
       fscanf( infp, "%s %g", label, &knn_sb_2b_err                          ) ;   printf( "knn_sb_2b_err                          %s %g\n", label, knn_sb_2b_err                         ) ;
       fscanf( infp, "%s %g", label, &knn_sb_3b_mean                         ) ;   printf( "knn_sb_3b_mean                         %s %g\n", label, knn_sb_3b_mean                        ) ;
       fscanf( infp, "%s %g", label, &knn_sb_3b_err                          ) ;   printf( "knn_sb_3b_err                          %s %g\n", label, knn_sb_3b_err                         ) ;
       fscanf( infp, "%s %g", label, &fsig_ee_mean                           ) ;   printf( "fsig_ee_mean                           %s %g\n", label, fsig_ee_mean                          ) ;
       fscanf( infp, "%s %g", label, &fsig_ee_err                            ) ;   printf( "fsig_ee_err                            %s %g\n", label, fsig_ee_err                           ) ;
       fscanf( infp, "%s %g", label, &fsig_mm_mean                           ) ;   printf( "fsig_mm_mean                           %s %g\n", label, fsig_mm_mean                          ) ;
       fscanf( infp, "%s %g", label, &fsig_mm_err                            ) ;   printf( "fsig_mm_err                            %s %g\n", label, fsig_mm_err                           ) ;
       fscanf( infp, "%s %g", label, &sf_mc                                  ) ;   printf( "sf_mc                                  %s %g\n", label, sf_mc                                 ) ;
       fscanf( infp, "%s %g", label, &sf_mc_err                              ) ;   printf( "sf_mc_err                              %s %g\n", label, sf_mc_err                             ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb_1b                           ) ;   printf( "sf_qcd_sb_1b                           %s %g\n", label, sf_qcd_sb_1b                          ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb_1b_err                       ) ;   printf( "sf_qcd_sb_1b_err                       %s %g\n", label, sf_qcd_sb_1b_err                      ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb_2b                           ) ;   printf( "sf_qcd_sb_2b                           %s %g\n", label, sf_qcd_sb_2b                          ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb_2b_err                       ) ;   printf( "sf_qcd_sb_2b_err                       %s %g\n", label, sf_qcd_sb_2b_err                      ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb_3b                           ) ;   printf( "sf_qcd_sb_3b                           %s %g\n", label, sf_qcd_sb_3b                          ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sb_3b_err                       ) ;   printf( "sf_qcd_sb_3b_err                       %s %g\n", label, sf_qcd_sb_3b_err                      ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig_1b                          ) ;   printf( "sf_qcd_sig_1b                          %s %g\n", label, sf_qcd_sig_1b                         ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig_1b_err                      ) ;   printf( "sf_qcd_sig_1b_err                      %s %g\n", label, sf_qcd_sig_1b_err                     ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig_2b                          ) ;   printf( "sf_qcd_sig_2b                          %s %g\n", label, sf_qcd_sig_2b                         ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig_2b_err                      ) ;   printf( "sf_qcd_sig_2b_err                      %s %g\n", label, sf_qcd_sig_2b_err                     ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig_3b                          ) ;   printf( "sf_qcd_sig_3b                          %s %g\n", label, sf_qcd_sig_3b                         ) ;
       fscanf( infp, "%s %g", label, &sf_qcd_sig_3b_err                      ) ;   printf( "sf_qcd_sig_3b_err                      %s %g\n", label, sf_qcd_sig_3b_err                     ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig_1b                         ) ;   printf( "sf_ttwj_sig_1b                         %s %g\n", label, sf_ttwj_sig_1b                        ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig_1b_err                     ) ;   printf( "sf_ttwj_sig_1b_err                     %s %g\n", label, sf_ttwj_sig_1b_err                    ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig_2b                         ) ;   printf( "sf_ttwj_sig_2b                         %s %g\n", label, sf_ttwj_sig_2b                        ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig_2b_err                     ) ;   printf( "sf_ttwj_sig_2b_err                     %s %g\n", label, sf_ttwj_sig_2b_err                    ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig_3b                         ) ;   printf( "sf_ttwj_sig_3b                         %s %g\n", label, sf_ttwj_sig_3b                        ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sig_3b_err                     ) ;   printf( "sf_ttwj_sig_3b_err                     %s %g\n", label, sf_ttwj_sig_3b_err                    ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sb_2b                          ) ;   printf( "sf_ttwj_sb_2b                          %s %g\n", label, sf_ttwj_sb_2b                         ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sb_2b_err                      ) ;   printf( "sf_ttwj_sb_2b_err                      %s %g\n", label, sf_ttwj_sb_2b_err                     ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sb_3b                          ) ;   printf( "sf_ttwj_sb_3b                          %s %g\n", label, sf_ttwj_sb_3b                         ) ;
       fscanf( infp, "%s %g", label, &sf_ttwj_sb_3b_err                      ) ;   printf( "sf_ttwj_sb_3b_err                      %s %g\n", label, sf_ttwj_sb_3b_err                     ) ;
       fscanf( infp, "%s %g", label, &sf_ee                               ) ;   printf( "sf_ee                               %s %g\n", label, sf_ee                              ) ;
       fscanf( infp, "%s %g", label, &sf_ee_err                           ) ;   printf( "sf_ee_err                           %s %g\n", label, sf_ee_err                          ) ;
       fscanf( infp, "%s %g", label, &sf_mm                               ) ;   printf( "sf_mm                               %s %g\n", label, sf_mm                              ) ;
       fscanf( infp, "%s %g", label, &sf_mm_err                           ) ;   printf( "sf_mm_err                           %s %g\n", label, sf_mm_err                          ) ;
       fscanf( infp, "%s %g", label, &btageff_err                           ) ;   printf( "btageff_err                           %s %g\n", label, btageff_err                          ) ;





       printf("\n Done reading in %s\n\n", infile ) ;
       fclose( infp ) ;

    //---- calculations for determining initial values for floating parameters.

     //-- Znunu stuff

       float initialval_znn_sig_ee_1b(2.) ;
       float initialval_znn_sig_ee_2b(2.) ;
       float initialval_znn_sig_ee_3b(2.) ;
       float initialval_znn_sig_mm_1b(2.) ;
       float initialval_znn_sig_mm_2b(2.) ;
       float initialval_znn_sig_mm_3b(2.) ;
       float initialval_znn_sb_ee_1b(2.) ;
       float initialval_znn_sb_ee_2b(2.) ;
       float initialval_znn_sb_ee_3b(2.) ;
       float initialval_znn_sb_mm_1b(2.) ;
       float initialval_znn_sb_mm_2b(2.) ;
       float initialval_znn_sb_mm_3b(2.) ;
       float initialval_znn_sig_1b(2.) ;
       float initialval_znn_sig_2b(2.) ;
       float initialval_znn_sig_3b(2.) ;
       float initialval_znn_sb_1b(2.) ;
       float initialval_znn_sb_2b(2.) ;
       float initialval_znn_sb_3b(2.) ;


     //--- 1b

       initialval_znn_sig_ee_1b = (Nsig_ee) * ( 5.95 * DataLumi * fsig_ee_mean * knn_sig_1b_mean ) / ( acc_ee_sig_mean * eff_ee_mean * Ztoll_lumi ) ;
       initialval_znn_sb_ee_1b  = (Nsb_ee ) * ( 5.95 * DataLumi * fsig_ee_mean * knn_sb_1b_mean  ) / ( acc_ee_sb_mean  * eff_ee_mean * Ztoll_lumi ) ;

       initialval_znn_sig_mm_1b = (Nsig_mm) * ( 5.95 * DataLumi * fsig_mm_mean * knn_sig_1b_mean ) / ( acc_mm_sig_mean * eff_mm_mean * Ztoll_lumi ) ;
       initialval_znn_sb_mm_1b  = (Nsb_mm ) * ( 5.95 * DataLumi * fsig_mm_mean * knn_sb_1b_mean  ) / ( acc_mm_sb_mean  * eff_mm_mean * Ztoll_lumi ) ;


       //-- really dumb ave.
       initialval_znn_sig_1b = 0.5 * ( initialval_znn_sig_ee_1b + initialval_znn_sig_mm_1b ) ;
       initialval_znn_sb_1b  = 0.5 * ( initialval_znn_sb_ee_1b  + initialval_znn_sb_mm_1b ) ;




     //--- 2b

       initialval_znn_sig_ee_2b = (Nsig_ee) * ( 5.95 * DataLumi * fsig_ee_mean * knn_sig_2b_mean ) / ( acc_ee_sig_mean * eff_ee_mean * Ztoll_lumi ) ;
       initialval_znn_sb_ee_2b  = (Nsb_ee ) * ( 5.95 * DataLumi * fsig_ee_mean * knn_sb_2b_mean  ) / ( acc_ee_sb_mean  * eff_ee_mean * Ztoll_lumi ) ;

       initialval_znn_sig_mm_2b = (Nsig_mm) * ( 5.95 * DataLumi * fsig_mm_mean * knn_sig_2b_mean ) / ( acc_mm_sig_mean * eff_mm_mean * Ztoll_lumi ) ;
       initialval_znn_sb_mm_2b  = (Nsb_mm ) * ( 5.95 * DataLumi * fsig_mm_mean * knn_sb_2b_mean  ) / ( acc_mm_sb_mean  * eff_mm_mean * Ztoll_lumi ) ;


       //-- really dumb ave.
       initialval_znn_sig_2b = 0.5 * ( initialval_znn_sig_ee_2b + initialval_znn_sig_mm_2b ) ;
       initialval_znn_sb_2b  = 0.5 * ( initialval_znn_sb_ee_2b  + initialval_znn_sb_mm_2b ) ;




     //--- 3b

       initialval_znn_sig_ee_3b = (Nsig_ee) * ( 5.95 * DataLumi * fsig_ee_mean * knn_sig_3b_mean ) / ( acc_ee_sig_mean * eff_ee_mean * Ztoll_lumi ) ;
       initialval_znn_sb_ee_3b  = (Nsb_ee ) * ( 5.95 * DataLumi * fsig_ee_mean * knn_sb_3b_mean  ) / ( acc_ee_sb_mean  * eff_ee_mean * Ztoll_lumi ) ;

       initialval_znn_sig_mm_3b = (Nsig_mm) * ( 5.95 * DataLumi * fsig_mm_mean * knn_sig_3b_mean ) / ( acc_mm_sig_mean * eff_mm_mean * Ztoll_lumi ) ;
       initialval_znn_sb_mm_3b  = (Nsb_mm ) * ( 5.95 * DataLumi * fsig_mm_mean * knn_sb_3b_mean  ) / ( acc_mm_sb_mean  * eff_mm_mean * Ztoll_lumi ) ;


       //-- really dumb ave.
       initialval_znn_sig_3b = 0.5 * ( initialval_znn_sig_ee_3b + initialval_znn_sig_mm_3b ) ;
       initialval_znn_sb_3b  = 0.5 * ( initialval_znn_sb_ee_3b  + initialval_znn_sb_mm_3b ) ;









     //--- QCD ldp stuff

       double initialval_qcd_sig_ldp_1b = Nsig_ldp_1b - ( Nttbarsingletopzjetsmc_sig_ldp_1b + NWJmc_sig_ldp_1b + NZnnmc_sig_ldp_1b ) ;
       double initialval_qcd_sb_ldp_1b  = Nsb_ldp_1b  - ( Nttbarsingletopzjetsmc_sb_ldp_1b  + NWJmc_sb_ldp_1b  + NZnnmc_sb_ldp_1b  ) ;

       double initialval_qcd_sig_ldp_2b = Nsig_ldp_2b - ( Nttbarsingletopzjetsmc_sig_ldp_2b + NWJmc_sig_ldp_2b + NZnnmc_sig_ldp_2b ) ;
       double initialval_qcd_sb_ldp_2b  = Nsb_ldp_2b  - ( Nttbarsingletopzjetsmc_sb_ldp_2b  + NWJmc_sb_ldp_2b  + NZnnmc_sb_ldp_2b  ) ;

       double initialval_qcd_sig_ldp_3b = Nsig_ldp_3b - ( Nttbarsingletopzjetsmc_sig_ldp_3b + NWJmc_sig_ldp_3b + NZnnmc_sig_ldp_3b ) ;
       double initialval_qcd_sb_ldp_3b  = Nsb_ldp_3b  - ( Nttbarsingletopzjetsmc_sb_ldp_3b  + NWJmc_sb_ldp_3b  + NZnnmc_sb_ldp_3b  ) ;


     //--- QCD sig and sb stuff

       double initialval_qcd_sig_1b = Rlsb_passfail * initialval_qcd_sig_ldp_1b ;
       double initialval_qcd_sb_1b  = Rlsb_passfail * initialval_qcd_sb_ldp_1b  ;

       double initialval_qcd_sig_2b = Rlsb_passfail * initialval_qcd_sig_ldp_2b ;
       double initialval_qcd_sb_2b  = Rlsb_passfail * initialval_qcd_sb_ldp_2b  ;

       double initialval_qcd_sig_3b = Rlsb_passfail * initialval_qcd_sig_ldp_3b ;
       double initialval_qcd_sb_3b  = Rlsb_passfail * initialval_qcd_sb_ldp_3b  ;


     //--- ttwj SL

       double initialval_ttwj_sig_sl_1b = Nsig_sl_1b ;
       double initialval_ttwj_sb_sl_1b  = Nsb_sl_1b ;

       double initialval_ttwj_sig_sl_2b = Nsig_sl_2b ;
       double initialval_ttwj_sb_sl_2b  = Nsb_sl_2b ;

       double initialval_ttwj_sig_sl_3b = Nsig_sl_3b ;
       double initialval_ttwj_sb_sl_3b  = Nsb_sl_3b ;


     //--- ttwj sig and sb

       //-- 1b
       if ( initialval_ttwj_sb_sl_1b <= 0 ) {
          printf("\n\n\n *** initialval_ttwj_sb_sl_1b has a crazy value (%.1f).  I quit.\n\n\n", initialval_ttwj_sb_sl_1b ) ;
          return false ;
       }
       double initialval_ttwj_sb_1b     = Nsb_1b - ( initialval_qcd_sb_1b + initialval_znn_sb_1b ) ;
       double initialval_ttwj_sig_1b    = initialval_ttwj_sb_1b * ( initialval_ttwj_sig_sl_1b / initialval_ttwj_sb_sl_1b ) ;

       //-- 2b
       if ( initialval_ttwj_sb_sl_2b <= 0 ) {
          printf("\n\n\n *** initialval_ttwj_sb_sl_2b has a crazy value (%.1f).  I quit.\n\n\n", initialval_ttwj_sb_sl_2b ) ;
          return false ;
       }
       double initialval_ttwj_sb_2b     = Nsb_2b - ( initialval_qcd_sb_2b + initialval_znn_sb_2b ) ;
       double initialval_ttwj_sig_2b    = initialval_ttwj_sb_2b * ( initialval_ttwj_sig_sl_2b / initialval_ttwj_sb_sl_2b ) ;

       //-- 3b
       if ( initialval_ttwj_sb_sl_3b <= 0 ) {
          printf("\n\n\n *** initialval_ttwj_sb_sl_3b has a crazy value (%.1f).  I quit.\n\n\n", initialval_ttwj_sb_sl_3b ) ;
          return false ;
       }
       double initialval_ttwj_sb_3b     = Nsb_3b - ( initialval_qcd_sb_3b + initialval_znn_sb_3b ) ;
       double initialval_ttwj_sig_3b    = initialval_ttwj_sb_3b * ( initialval_ttwj_sig_sl_3b / initialval_ttwj_sb_sl_3b ) ;






       printf("\n\n\n --------- Observables and floating parameter initial values. ------------\n\n") ;

       printf("             |  Nobs   ||  ttwj  |  QCD  |  Znn  |\n") ;
       printf("-------------------------------------------------------------\n") ;
       printf(" 1b SIG      | %5d   || %5.1f | %5.1f | %5.1f |\n", Nsig_1b, initialval_ttwj_sig_1b, initialval_qcd_sig_1b, initialval_znn_sig_1b ) ;
       printf(" 1b SB       | %5d   || %5.1f | %5.1f | %5.1f |\n", Nsb_1b , initialval_ttwj_sb_1b , initialval_qcd_sb_1b , initialval_znn_sb_1b  ) ;
       printf(" 1b SIG,LDP  | %5d   ||*%5.1f | %5.1f |*%5.1f |\n", Nsig_ldp_1b, (Nttbarsingletopzjetsmc_sig_ldp_1b+NWJmc_sig_ldp_1b), initialval_qcd_sig_ldp_1b, NZnnmc_sig_ldp_1b ) ;
       printf(" 1b SB ,LDP  | %5d   ||*%5.1f | %5.1f |*%5.1f |\n", Nsb_ldp_1b , (Nttbarsingletopzjetsmc_sb_ldp_1b +NWJmc_sb_ldp_1b ) , initialval_qcd_sb_ldp_1b , NZnnmc_sb_ldp_1b  ) ;
       printf("-------------------------------------------------------------\n") ;
       printf(" 2b SIG      | %5d   || %5.1f | %5.1f | %5.1f |\n", Nsig_2b, initialval_ttwj_sig_2b, initialval_qcd_sig_2b, initialval_znn_sig_2b ) ;
       printf(" 2b SB       | %5d   || %5.1f | %5.1f | %5.1f |\n", Nsb_2b , initialval_ttwj_sb_2b , initialval_qcd_sb_2b , initialval_znn_sb_2b  ) ;
       printf(" 2b SIG,LDP  | %5d   ||*%5.1f | %5.1f |*%5.1f |\n", Nsig_ldp_2b, (Nttbarsingletopzjetsmc_sig_ldp_2b+NWJmc_sig_ldp_2b), initialval_qcd_sig_ldp_2b, NZnnmc_sig_ldp_2b ) ;
       printf(" 2b SB ,LDP  | %5d   ||*%5.1f | %5.1f |*%5.1f |\n", Nsb_ldp_2b , (Nttbarsingletopzjetsmc_sb_ldp_2b +NWJmc_sb_ldp_2b ) , initialval_qcd_sb_ldp_2b , NZnnmc_sb_ldp_2b  ) ;
       printf("-------------------------------------------------------------\n") ;
       printf(" 3b SIG      | %5d   || %5.1f | %5.1f | %5.1f |\n", Nsig_3b, initialval_ttwj_sig_3b, initialval_qcd_sig_3b, initialval_znn_sig_3b ) ;
       printf(" 3b SB       | %5d   || %5.1f | %5.1f | %5.1f |\n", Nsb_3b , initialval_ttwj_sb_3b , initialval_qcd_sb_3b , initialval_znn_sb_3b  ) ;
       printf(" 3b SIG,LDP  | %5d   ||*%5.1f | %5.1f |*%5.1f |\n", Nsig_ldp_3b, (Nttbarsingletopzjetsmc_sig_ldp_3b+NWJmc_sig_ldp_3b), initialval_qcd_sig_ldp_3b, NZnnmc_sig_ldp_3b ) ;
       printf(" 3b SB ,LDP  | %5d   ||*%5.1f | %5.1f |*%5.1f |\n", Nsb_ldp_3b , (Nttbarsingletopzjetsmc_sb_ldp_3b +NWJmc_sb_ldp_3b ) , initialval_qcd_sb_ldp_3b , NZnnmc_sb_ldp_3b  ) ;
       printf("-------------------------------------------------------------\n") ;
       printf("\n * means fixed MC value.\n\n\n") ;




     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


       printf(" --- Defining observables.\n" ) ;


      rv_Nsig_1b        = new RooRealVar( "Nsig_1b"        , "Nsig_1b"        , 0.0, 1000000. ) ;
      rv_Nsb_1b         = new RooRealVar( "Nsb_1b"         , "Nsb_1b"         , 0.0, 1000000. ) ;

      rv_Nsig_sl_1b     = new RooRealVar( "Nsig_sl_1b"     , "Nsig_sl_1b"     , 0.0, 1000000. ) ;
      rv_Nsb_sl_1b      = new RooRealVar( "Nsb_sl_1b"      , "Nsb_sl_1b"      , 0.0, 1000000. ) ;

      rv_Nsig_ldp_1b    = new RooRealVar( "Nsig_ldp_1b"    , "Nsig_ldp_1b"    , 0.0, 1000000. ) ;
      rv_Nsb_ldp_1b     = new RooRealVar( "Nsb_ldp_1b"     , "Nsb_ldp_1b"     , 0.0, 1000000. ) ;


      rv_Nsig_2b        = new RooRealVar( "Nsig_2b"        , "Nsig_2b"        , 0.0, 1000000. ) ;
      rv_Nsb_2b         = new RooRealVar( "Nsb_2b"         , "Nsb_2b"         , 0.0, 1000000. ) ;

      rv_Nsig_sl_2b     = new RooRealVar( "Nsig_sl_2b"     , "Nsig_sl_2b"     , 0.0, 1000000. ) ;
      rv_Nsb_sl_2b      = new RooRealVar( "Nsb_sl_2b"      , "Nsb_sl_2b"      , 0.0, 1000000. ) ;

      rv_Nsig_ldp_2b    = new RooRealVar( "Nsig_ldp_2b"    , "Nsig_ldp_2b"    , 0.0, 1000000. ) ;
      rv_Nsb_ldp_2b     = new RooRealVar( "Nsb_ldp_2b"     , "Nsb_ldp_2b"     , 0.0, 1000000. ) ;


      rv_Nsig_3b        = new RooRealVar( "Nsig_3b"        , "Nsig_3b"        , 0.0, 1000000. ) ;
      rv_Nsb_3b         = new RooRealVar( "Nsb_3b"         , "Nsb_3b"         , 0.0, 1000000. ) ;

      rv_Nsig_sl_3b     = new RooRealVar( "Nsig_sl_3b"     , "Nsig_sl_3b"     , 0.0, 1000000. ) ;
      rv_Nsb_sl_3b      = new RooRealVar( "Nsb_sl_3b"      , "Nsb_sl_3b"      , 0.0, 1000000. ) ;

      rv_Nsig_ldp_3b    = new RooRealVar( "Nsig_ldp_3b"    , "Nsig_ldp_3b"    , 0.0, 1000000. ) ;
      rv_Nsb_ldp_3b     = new RooRealVar( "Nsb_ldp_3b"     , "Nsb_ldp_3b"     , 0.0, 1000000. ) ;



      rv_Nsb_ee      = new RooRealVar( "Nsb_ee"      ,"Nsb_ee"       , 0., 10000. ) ;
      rv_Nsig_ee     = new RooRealVar( "Nsig_ee"     ,"Nsig_ee"      , 0., 10000. ) ;

      rv_Nsb_mm      = new RooRealVar( "Nsb_mm"      ,"Nsb_mm"       , 0., 10000. ) ;
      rv_Nsig_mm     = new RooRealVar( "Nsig_mm"     ,"Nsig_mm"      , 0., 10000. ) ;






      rv_Nsig_1b        -> setVal( Nsig_1b ) ;
      rv_Nsb_1b         -> setVal( Nsb_1b ) ;

      rv_Nsig_sl_1b     -> setVal( Nsig_sl_1b ) ;
      rv_Nsb_sl_1b      -> setVal( Nsb_sl_1b ) ;

      rv_Nsig_ldp_1b    -> setVal( Nsig_ldp_1b ) ;
      rv_Nsb_ldp_1b     -> setVal( Nsb_ldp_1b ) ;



      rv_Nsig_2b        -> setVal( Nsig_2b ) ;
      rv_Nsb_2b         -> setVal( Nsb_2b ) ;

      rv_Nsig_sl_2b     -> setVal( Nsig_sl_2b ) ;
      rv_Nsb_sl_2b      -> setVal( Nsb_sl_2b ) ;

      rv_Nsig_ldp_2b    -> setVal( Nsig_ldp_2b ) ;
      rv_Nsb_ldp_2b     -> setVal( Nsb_ldp_2b ) ;



      rv_Nsig_3b        -> setVal( Nsig_3b ) ;
      rv_Nsb_3b         -> setVal( Nsb_3b ) ;

      rv_Nsig_sl_3b     -> setVal( Nsig_sl_3b ) ;
      rv_Nsb_sl_3b      -> setVal( Nsb_sl_3b ) ;

      rv_Nsig_ldp_3b    -> setVal( Nsig_ldp_3b ) ;
      rv_Nsb_ldp_3b     -> setVal( Nsb_ldp_3b ) ;



      rv_Nsb_ee      -> setVal( Nsb_ee ) ;
      rv_Nsig_ee     -> setVal( Nsig_ee ) ;

      rv_Nsb_mm      -> setVal( Nsb_mm ) ;
      rv_Nsig_mm     -> setVal( Nsig_mm ) ;








    //++++++++ Parameters of the likelihood +++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining parameters.\n" ) ;



    //____ Counts in SIG ______________________

      if ( useSigTtwjVar ) {
         rrv_mu_ttwj_sig_1b = new RooRealVar( "mu_ttwj_sig_1b"   , "mu_ttwj_sig_1b"   , 0.0, 10000. ) ;
         rv_mu_ttwj_sig_1b = rrv_mu_ttwj_sig_1b ;
         rrv_mu_ttwj_sig_1b   -> setVal( initialval_ttwj_sig_1b ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sig_1b  = new RooRealVar( "mu_qcd_sig_1b"    , "mu_qcd_sig_1b"    , 0.0, 20000. ) ;
         rv_mu_qcd_sig_1b = rrv_mu_qcd_sig_1b ;
         rrv_mu_qcd_sig_1b  -> setVal( initialval_qcd_sig_1b ) ; //-- this is a starting value only.
      }


      if ( useSigTtwjVar ) {
         rrv_mu_ttwj_sig_2b = new RooRealVar( "mu_ttwj_sig_2b"   , "mu_ttwj_sig_2b"   , 0.0, 10000. ) ;
         rv_mu_ttwj_sig_2b = rrv_mu_ttwj_sig_2b ;
         rrv_mu_ttwj_sig_2b   -> setVal( initialval_ttwj_sig_2b ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sig_2b  = new RooRealVar( "mu_qcd_sig_2b"    , "mu_qcd_sig_2b"    , 0.0, 20000. ) ;
         rv_mu_qcd_sig_2b = rrv_mu_qcd_sig_2b ;
         rrv_mu_qcd_sig_2b  -> setVal( initialval_qcd_sig_2b ) ; //-- this is a starting value only.
      }


      if ( useSigTtwjVar ) {
         rrv_mu_ttwj_sig_3b = new RooRealVar( "mu_ttwj_sig_3b"   , "mu_ttwj_sig_3b"   , 0.0, 10000. ) ;
         rv_mu_ttwj_sig_3b = rrv_mu_ttwj_sig_3b ;
         rrv_mu_ttwj_sig_3b   -> setVal( initialval_ttwj_sig_3b ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sig_3b  = new RooRealVar( "mu_qcd_sig_3b"    , "mu_qcd_sig_3b"    , 0.0, 20000. ) ;
         rv_mu_qcd_sig_3b = rrv_mu_qcd_sig_3b ;
         rrv_mu_qcd_sig_3b  -> setVal( initialval_qcd_sig_3b ) ; //-- this is a starting value only.
      }






      rv_mu_znn_sig_1b      = new RooRealVar( "mu_znn_sig_1b"    , "mu_znn_sig_1b"    , 0.0, 3000. ) ;

      float maxSusySig_1b = 4.0*Nsig_1b ;
      rv_mu_susy_sig_1b     = new RooRealVar( "mu_susy_sig_1b"   , "mu_susy_sig_1b"   , 0.0, maxSusySig_1b ) ;


      rv_mu_znn_sig_1b   -> setVal( initialval_znn_sig_1b ) ;  //-- this is a starting value only.
      rv_mu_susy_sig_1b    -> setVal( 0. ) ;  //-- this is a starting value only.







    //____ Counts in SB  ______________________


      if ( !useSigTtwjVar ) {
         rrv_mu_ttwj_sb_1b  = new RooRealVar( "mu_ttwj_sb_1b"    , "mu_ttwj_sb_1b"    , 0.0, 10000. ) ;
         rv_mu_ttwj_sb_1b = rrv_mu_ttwj_sb_1b ;
         rrv_mu_ttwj_sb_1b   -> setVal( initialval_ttwj_sb_1b ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sb_1b  = new RooRealVar( "mu_qcd_sb_1b"    , "mu_qcd_sb_1b"    , 0.0, 500. ) ;
         rv_mu_qcd_sb_1b = rrv_mu_qcd_sb_1b ;
         rrv_mu_qcd_sb_1b  -> setVal( initialval_qcd_sb_1b ) ; //-- this is a starting value only.
      }



      if ( !useSigTtwjVar ) {
         rrv_mu_ttwj_sb_2b  = new RooRealVar( "mu_ttwj_sb_2b"    , "mu_ttwj_sb_2b"    , 0.0, 10000. ) ;
         rv_mu_ttwj_sb_2b = rrv_mu_ttwj_sb_2b ;
         rrv_mu_ttwj_sb_2b   -> setVal( initialval_ttwj_sb_2b ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sb_2b  = new RooRealVar( "mu_qcd_sb_2b"    , "mu_qcd_sb_2b"    , 0.0, 500. ) ;
         rv_mu_qcd_sb_2b = rrv_mu_qcd_sb_2b ;
         rrv_mu_qcd_sb_2b  -> setVal( initialval_qcd_sb_2b ) ; //-- this is a starting value only.
      }



      if ( !useSigTtwjVar ) {
         rrv_mu_ttwj_sb_3b  = new RooRealVar( "mu_ttwj_sb_3b"    , "mu_ttwj_sb_3b"    , 0.0, 10000. ) ;
         rv_mu_ttwj_sb_3b = rrv_mu_ttwj_sb_3b ;
         rrv_mu_ttwj_sb_3b   -> setVal( initialval_ttwj_sb_3b ) ;  //-- this is a starting value only.
      }
      if ( !useLdpVars ) {
         rrv_mu_qcd_sb_3b  = new RooRealVar( "mu_qcd_sb_3b"    , "mu_qcd_sb_3b"    , 0.0, 500. ) ;
         rv_mu_qcd_sb_3b = rrv_mu_qcd_sb_3b ;
         rrv_mu_qcd_sb_3b  -> setVal( initialval_qcd_sb_3b ) ; //-- this is a starting value only.
      }




      //-- Note: QCD is rfv
      //-- Note: SUSY is rfv


      rrv_mu_znn_sb_1b       = new RooRealVar( "mu_znn_sb_1b"     , "mu_znn_sb_1b"     , 0.0, 3500. ) ;

      rrv_mu_znn_sb_1b   -> setVal( initialval_znn_sb_1b ) ;  //-- this is a starting value only.

      rv_mu_znn_sb_1b = rrv_mu_znn_sb_1b ;




      rrv_mu_znn_sb_2b       = new RooRealVar( "mu_znn_sb_2b"     , "mu_znn_sb_2b"     , 0.0, 3500. ) ;

      rrv_mu_znn_sb_2b   -> setVal( initialval_znn_sb_2b ) ;  //-- this is a starting value only.

      rv_mu_znn_sb_2b = rrv_mu_znn_sb_2b ;




      rrv_mu_znn_sb_3b       = new RooRealVar( "mu_znn_sb_3b"     , "mu_znn_sb_3b"     , 0.0, 3500. ) ;

      rrv_mu_znn_sb_3b   -> setVal( initialval_znn_sb_3b ) ;  //-- this is a starting value only.

      rv_mu_znn_sb_3b = rrv_mu_znn_sb_3b ;









      //-- Note: Znn is rfv in Znn model 2.






    //____ Counts in SIG, SL  ______________________

      rv_mu_ttwj_sig_sl_1b  = new RooRealVar( "mu_ttwj_sig_sl_1b"    , "mu_ttwj_sig_sl_1b"    , 0.0, 2500. ) ;

      //-- Note: QCD and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      rv_mu_ttwj_sig_sl_1b  -> setVal( initialval_ttwj_sig_sl_1b ) ;  //-- this is a starting value only.





      rv_mu_ttwj_sig_sl_2b  = new RooRealVar( "mu_ttwj_sig_sl_2b"    , "mu_ttwj_sig_sl_2b"    , 0.0, 2500. ) ;

      //-- Note: QCD and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      rv_mu_ttwj_sig_sl_2b  -> setVal( initialval_ttwj_sig_sl_2b ) ;  //-- this is a starting value only.





      rv_mu_ttwj_sig_sl_3b  = new RooRealVar( "mu_ttwj_sig_sl_3b"    , "mu_ttwj_sig_sl_3b"    , 0.0, 2500. ) ;

      //-- Note: QCD and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      rv_mu_ttwj_sig_sl_3b  -> setVal( initialval_ttwj_sig_sl_3b ) ;  //-- this is a starting value only.







    //____ Counts in SB, SL  ______________________

      rv_mu_ttwj_sb_sl_1b  = new RooRealVar( "mu_ttwj_sb_sl_1b"    , "mu_ttwj_sb_sl_1b"    , 0.0, 30000. ) ;

      //-- Note: QCD and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      rv_mu_ttwj_sb_sl_1b  -> setVal( initialval_ttwj_sb_sl_1b ) ;  //-- this is a starting value only.





      rv_mu_ttwj_sb_sl_2b  = new RooRealVar( "mu_ttwj_sb_sl_2b"    , "mu_ttwj_sb_sl_2b"    , 0.0, 30000. ) ;

      //-- Note: QCD and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      rv_mu_ttwj_sb_sl_2b  -> setVal( initialval_ttwj_sb_sl_2b ) ;  //-- this is a starting value only.





      rv_mu_ttwj_sb_sl_3b  = new RooRealVar( "mu_ttwj_sb_sl_3b"    , "mu_ttwj_sb_sl_3b"    , 0.0, 30000. ) ;

      //-- Note: QCD and Znn are assumed to be negligible and are not explicitly included.
      //-- Note: SUSY is rfv

      rv_mu_ttwj_sb_sl_3b  -> setVal( initialval_ttwj_sb_sl_3b ) ;  //-- this is a starting value only.







    //____ Counts in SIG, LDP  ______________________

      if ( useLdpVars ) {
         rrv_mu_qcd_sig_ldp_1b  = new RooRealVar( "mu_qcd_sig_ldp_1b"    , "mu_qcd_sig_ldp_1b"    , 0.0, 3500. ) ;
         rv_mu_qcd_sig_ldp_1b = rrv_mu_qcd_sig_ldp_1b ;
         rrv_mu_qcd_sig_ldp_1b  -> setVal( initialval_qcd_sig_ldp_1b ) ; //-- this is a starting value only.
      }



      if ( useLdpVars ) {
         rrv_mu_qcd_sig_ldp_2b  = new RooRealVar( "mu_qcd_sig_ldp_2b"    , "mu_qcd_sig_ldp_2b"    , 0.0, 3500. ) ;
         rv_mu_qcd_sig_ldp_2b = rrv_mu_qcd_sig_ldp_2b ;
         rrv_mu_qcd_sig_ldp_2b  -> setVal( initialval_qcd_sig_ldp_2b ) ; //-- this is a starting value only.
      }



      if ( useLdpVars ) {
         rrv_mu_qcd_sig_ldp_3b  = new RooRealVar( "mu_qcd_sig_ldp_3b"    , "mu_qcd_sig_ldp_3b"    , 0.0, 3500. ) ;
         rv_mu_qcd_sig_ldp_3b = rrv_mu_qcd_sig_ldp_3b ;
         rrv_mu_qcd_sig_ldp_3b  -> setVal( initialval_qcd_sig_ldp_3b ) ; //-- this is a starting value only.
      }





      //-- Note: Znn is rfv (MC)
      //-- Note: ttwj is rfv (MC)
      //-- Note: SUSY is rfv








    //____ Counts in SB, LDP  ______________________

      if ( useLdpVars ) {
         rrv_mu_qcd_sb_ldp_1b  = new RooRealVar( "mu_qcd_sb_ldp_1b"    , "mu_qcd_sb_ldp_1b"    , 0.0, 3000. ) ;
         rv_mu_qcd_sb_ldp_1b = rrv_mu_qcd_sb_ldp_1b ;
         rrv_mu_qcd_sb_ldp_1b  -> setVal( initialval_qcd_sb_ldp_1b ) ; //-- this is a starting value only.
      }



      if ( useLdpVars ) {
         rrv_mu_qcd_sb_ldp_2b  = new RooRealVar( "mu_qcd_sb_ldp_2b"    , "mu_qcd_sb_ldp_2b"    , 0.0, 3000. ) ;
         rv_mu_qcd_sb_ldp_2b = rrv_mu_qcd_sb_ldp_2b ;
         rrv_mu_qcd_sb_ldp_2b  -> setVal( initialval_qcd_sb_ldp_2b ) ; //-- this is a starting value only.
      }



      if ( useLdpVars ) {
         rrv_mu_qcd_sb_ldp_3b  = new RooRealVar( "mu_qcd_sb_ldp_3b"    , "mu_qcd_sb_ldp_3b"    , 0.0, 3000. ) ;
         rv_mu_qcd_sb_ldp_3b = rrv_mu_qcd_sb_ldp_3b ;
         rrv_mu_qcd_sb_ldp_3b  -> setVal( initialval_qcd_sb_ldp_3b ) ; //-- this is a starting value only.
      }





      //-- Note: Znn is rfv (MC)
      //-- Note: ttwj is rfv (MC)
      //-- Note: SUSY is rfv














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

      rv_mu_susymc_sig_1b      = new RooRealVar( "mu_susymc_sig_1b"     , "mu_susymc_sig_1b"     , 0.0, 100000. ) ;
      rv_mu_susymc_sb_1b       = new RooRealVar( "mu_susymc_sb_1b"      , "mu_susymc_sb_1b"      , 0.0, 100000. ) ;
      rv_mu_susymc_sig_sl_1b   = new RooRealVar( "mu_susymc_sig_sl_1b"  , "mu_susymc_sig_sl_1b"  , 0.0, 100000. ) ;
      rv_mu_susymc_sb_sl_1b    = new RooRealVar( "mu_susymc_sb_sl_1b"   , "mu_susymc_sb_sl_1b"   , 0.0, 100000. ) ;
      rv_mu_susymc_sig_ldp_1b  = new RooRealVar( "mu_susymc_sig_ldp_1b" , "mu_susymc_sig_ldp_1b" , 0.0, 100000. ) ;
      rv_mu_susymc_sb_ldp_1b   = new RooRealVar( "mu_susymc_sb_ldp_1b"  , "mu_susymc_sb_ldp_1b"  , 0.0, 100000. ) ;

      rv_mu_susymc_sig_2b      = new RooRealVar( "mu_susymc_sig_2b"     , "mu_susymc_sig_2b"     , 0.0, 100000. ) ;
      rv_mu_susymc_sb_2b       = new RooRealVar( "mu_susymc_sb_2b"      , "mu_susymc_sb_2b"      , 0.0, 100000. ) ;
      rv_mu_susymc_sig_sl_2b   = new RooRealVar( "mu_susymc_sig_sl_2b"  , "mu_susymc_sig_sl_2b"  , 0.0, 100000. ) ;
      rv_mu_susymc_sb_sl_2b    = new RooRealVar( "mu_susymc_sb_sl_2b"   , "mu_susymc_sb_sl_2b"   , 0.0, 100000. ) ;
      rv_mu_susymc_sig_ldp_2b  = new RooRealVar( "mu_susymc_sig_ldp_2b" , "mu_susymc_sig_ldp_2b" , 0.0, 100000. ) ;
      rv_mu_susymc_sb_ldp_2b   = new RooRealVar( "mu_susymc_sb_ldp_2b"  , "mu_susymc_sb_ldp_2b"  , 0.0, 100000. ) ;

      rv_mu_susymc_sig_3b      = new RooRealVar( "mu_susymc_sig_3b"     , "mu_susymc_sig_3b"     , 0.0, 100000. ) ;
      rv_mu_susymc_sb_3b       = new RooRealVar( "mu_susymc_sb_3b"      , "mu_susymc_sb_3b"      , 0.0, 100000. ) ;
      rv_mu_susymc_sig_sl_3b   = new RooRealVar( "mu_susymc_sig_sl_3b"  , "mu_susymc_sig_sl_3b"  , 0.0, 100000. ) ;
      rv_mu_susymc_sb_sl_3b    = new RooRealVar( "mu_susymc_sb_sl_3b"   , "mu_susymc_sb_sl_3b"   , 0.0, 100000. ) ;
      rv_mu_susymc_sig_ldp_3b  = new RooRealVar( "mu_susymc_sig_ldp_3b" , "mu_susymc_sig_ldp_3b" , 0.0, 100000. ) ;
      rv_mu_susymc_sb_ldp_3b   = new RooRealVar( "mu_susymc_sb_ldp_3b"  , "mu_susymc_sb_ldp_3b"  , 0.0, 100000. ) ;



      rv_mu_susymc_sig_1b     -> setVal( 0.1 ) ;
      rv_mu_susymc_sb_1b      -> setVal( 0. ) ;
      rv_mu_susymc_sig_sl_1b  -> setVal( 0. ) ;
      rv_mu_susymc_sb_sl_1b   -> setVal( 0. ) ;
      rv_mu_susymc_sig_ldp_1b -> setVal( 0. ) ;
      rv_mu_susymc_sb_ldp_1b  -> setVal( 0. ) ;

      rv_mu_susymc_sig_2b     -> setVal( 0.1 ) ;
      rv_mu_susymc_sb_2b      -> setVal( 0. ) ;
      rv_mu_susymc_sig_sl_2b  -> setVal( 0. ) ;
      rv_mu_susymc_sb_sl_2b   -> setVal( 0. ) ;
      rv_mu_susymc_sig_ldp_2b -> setVal( 0. ) ;
      rv_mu_susymc_sb_ldp_2b  -> setVal( 0. ) ;

      rv_mu_susymc_sig_3b     -> setVal( 0.1 ) ;
      rv_mu_susymc_sb_3b      -> setVal( 0. ) ;
      rv_mu_susymc_sig_sl_3b  -> setVal( 0. ) ;
      rv_mu_susymc_sb_sl_3b   -> setVal( 0. ) ;
      rv_mu_susymc_sig_ldp_3b -> setVal( 0. ) ;
      rv_mu_susymc_sb_ldp_3b  -> setVal( 0. ) ;



      rv_mu_susymc_sig_1b     -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_1b      -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_sl_1b  -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_sl_1b   -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_ldp_1b -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_ldp_1b  -> setConstant(kTRUE) ;

      rv_mu_susymc_sig_2b     -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_2b      -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_sl_2b  -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_sl_2b   -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_ldp_2b -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_ldp_2b  -> setConstant(kTRUE) ;

      rv_mu_susymc_sig_3b     -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_3b      -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_sl_3b  -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_sl_3b   -> setConstant(kTRUE) ;
      rv_mu_susymc_sig_ldp_3b -> setConstant(kTRUE) ;
      rv_mu_susymc_sb_ldp_3b  -> setConstant(kTRUE) ;







     //-- SIG, LDP

      rv_mu_ttbarsingletopzjetsmc_sig_ldp_1b   = new RooRealVar( "mu_ttbarsingletopzjetsmc_sig_ldp_1b" ,"mu_ttbarsingletopzjetsmc_sig_ldp_1b" , 0., 100000. ) ;
      rv_mu_WJmc_sig_ldp_1b      = new RooRealVar( "mu_WJmc_sig_ldp_1b"    ,"mu_WJmc_sig_ldp_1b"    , 0., 100000. ) ;
      rv_mu_Znnmc_sig_ldp_1b     = new RooRealVar( "mu_Znnmc_sig_ldp_1b"   ,"mu_Znnmc_sig_ldp_1b"   , 0., 100000. ) ;

      rv_mu_ttbarsingletopzjetsmc_sig_ldp_2b   = new RooRealVar( "mu_ttbarsingletopzjetsmc_sig_ldp_2b" ,"mu_ttbarsingletopzjetsmc_sig_ldp_2b" , 0., 100000. ) ;
      rv_mu_WJmc_sig_ldp_2b      = new RooRealVar( "mu_WJmc_sig_ldp_2b"    ,"mu_WJmc_sig_ldp_2b"    , 0., 100000. ) ;
      rv_mu_Znnmc_sig_ldp_2b     = new RooRealVar( "mu_Znnmc_sig_ldp_2b"   ,"mu_Znnmc_sig_ldp_2b"   , 0., 100000. ) ;

      rv_mu_ttbarsingletopzjetsmc_sig_ldp_3b   = new RooRealVar( "mu_ttbarsingletopzjetsmc_sig_ldp_3b" ,"mu_ttbarsingletopzjetsmc_sig_ldp_3b" , 0., 100000. ) ;
      rv_mu_WJmc_sig_ldp_3b      = new RooRealVar( "mu_WJmc_sig_ldp_3b"    ,"mu_WJmc_sig_ldp_3b"    , 0., 100000. ) ;
      rv_mu_Znnmc_sig_ldp_3b     = new RooRealVar( "mu_Znnmc_sig_ldp_3b"   ,"mu_Znnmc_sig_ldp_3b"   , 0., 100000. ) ;



      rv_mu_ttbarsingletopzjetsmc_sig_ldp_1b  -> setVal( Nttbarsingletopzjetsmc_sig_ldp_1b ) ;
      rv_mu_WJmc_sig_ldp_1b     -> setVal( NWJmc_sig_ldp_1b ) ;
      rv_mu_Znnmc_sig_ldp_1b    -> setVal( NZnnmc_sig_ldp_1b ) ;

      rv_mu_ttbarsingletopzjetsmc_sig_ldp_2b  -> setVal( Nttbarsingletopzjetsmc_sig_ldp_2b ) ;
      rv_mu_WJmc_sig_ldp_2b     -> setVal( NWJmc_sig_ldp_2b ) ;
      rv_mu_Znnmc_sig_ldp_2b    -> setVal( NZnnmc_sig_ldp_2b ) ;

      rv_mu_ttbarsingletopzjetsmc_sig_ldp_3b  -> setVal( Nttbarsingletopzjetsmc_sig_ldp_3b ) ;
      rv_mu_WJmc_sig_ldp_3b     -> setVal( NWJmc_sig_ldp_3b ) ;
      rv_mu_Znnmc_sig_ldp_3b    -> setVal( NZnnmc_sig_ldp_3b ) ;



      rv_mu_ttbarsingletopzjetsmc_sig_ldp_1b  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sig_ldp_1b     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sig_ldp_1b    -> setConstant( kTRUE ) ;

      rv_mu_ttbarsingletopzjetsmc_sig_ldp_2b  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sig_ldp_2b     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sig_ldp_2b    -> setConstant( kTRUE ) ;

      rv_mu_ttbarsingletopzjetsmc_sig_ldp_3b  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sig_ldp_3b     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sig_ldp_3b    -> setConstant( kTRUE ) ;






     //-- SB, LDP

      rv_mu_ttbarsingletopzjetsmc_sb_ldp_1b   = new RooRealVar( "mu_ttbarsingletopzjetsmc_sb_ldp_1b" ,"mu_ttbarsingletopzjetsmc_sb_ldp_1b" , 0., 100000. ) ;
      rv_mu_WJmc_sb_ldp_1b      = new RooRealVar( "mu_WJmc_sb_ldp_1b"    ,"mu_WJmc_sb_ldp_1b"    , 0., 100000. ) ;
      rv_mu_Znnmc_sb_ldp_1b     = new RooRealVar( "mu_Znnmc_sb_ldp_1b"   ,"mu_Znnmc_sb_ldp_1b"   , 0., 100000. ) ;

      rv_mu_ttbarsingletopzjetsmc_sb_ldp_2b   = new RooRealVar( "mu_ttbarsingletopzjetsmc_sb_ldp_2b" ,"mu_ttbarsingletopzjetsmc_sb_ldp_2b" , 0., 100000. ) ;
      rv_mu_WJmc_sb_ldp_2b      = new RooRealVar( "mu_WJmc_sb_ldp_2b"    ,"mu_WJmc_sb_ldp_2b"    , 0., 100000. ) ;
      rv_mu_Znnmc_sb_ldp_2b     = new RooRealVar( "mu_Znnmc_sb_ldp_2b"   ,"mu_Znnmc_sb_ldp_2b"   , 0., 100000. ) ;

      rv_mu_ttbarsingletopzjetsmc_sb_ldp_3b   = new RooRealVar( "mu_ttbarsingletopzjetsmc_sb_ldp_3b" ,"mu_ttbarsingletopzjetsmc_sb_ldp_3b" , 0., 100000. ) ;
      rv_mu_WJmc_sb_ldp_3b      = new RooRealVar( "mu_WJmc_sb_ldp_3b"    ,"mu_WJmc_sb_ldp_3b"    , 0., 100000. ) ;
      rv_mu_Znnmc_sb_ldp_3b     = new RooRealVar( "mu_Znnmc_sb_ldp_3b"   ,"mu_Znnmc_sb_ldp_3b"   , 0., 100000. ) ;



      rv_mu_ttbarsingletopzjetsmc_sb_ldp_1b  -> setVal( Nttbarsingletopzjetsmc_sb_ldp_1b ) ;
      rv_mu_WJmc_sb_ldp_1b     -> setVal( NWJmc_sb_ldp_1b ) ;
      rv_mu_Znnmc_sb_ldp_1b    -> setVal( NZnnmc_sb_ldp_1b ) ;

      rv_mu_ttbarsingletopzjetsmc_sb_ldp_2b  -> setVal( Nttbarsingletopzjetsmc_sb_ldp_2b ) ;
      rv_mu_WJmc_sb_ldp_2b     -> setVal( NWJmc_sb_ldp_2b ) ;
      rv_mu_Znnmc_sb_ldp_2b    -> setVal( NZnnmc_sb_ldp_2b ) ;

      rv_mu_ttbarsingletopzjetsmc_sb_ldp_3b  -> setVal( Nttbarsingletopzjetsmc_sb_ldp_3b ) ;
      rv_mu_WJmc_sb_ldp_3b     -> setVal( NWJmc_sb_ldp_3b ) ;
      rv_mu_Znnmc_sb_ldp_3b    -> setVal( NZnnmc_sb_ldp_3b ) ;



      rv_mu_ttbarsingletopzjetsmc_sb_ldp_1b  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sb_ldp_1b     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sb_ldp_1b    -> setConstant( kTRUE ) ;

      rv_mu_ttbarsingletopzjetsmc_sb_ldp_2b  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sb_ldp_2b     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sb_ldp_2b    -> setConstant( kTRUE ) ;

      rv_mu_ttbarsingletopzjetsmc_sb_ldp_3b  -> setConstant( kTRUE ) ;
      rv_mu_WJmc_sb_ldp_3b     -> setConstant( kTRUE ) ;
      rv_mu_Znnmc_sb_ldp_3b    -> setConstant( kTRUE ) ;






    //+++++++ Gaussian constraints ++++++++++++++++++++++++++++++++

      printf(" --- Defining Gaussian constraint and constant parameters.\n" ) ;

    //_______ Efficiency scale factor.  Applied to SUSY and all MC inputs _______________

    //   August 23, 2011:  switching to correlated log-normal PDFs
    //
    //
    //  double pmin, pmax ;


    //--- mean parameters.
      rv_mean_eff_sf_sig_1b     = new RooRealVar( "mean_eff_sf_sig_1b"    , "mean_eff_sf_sig_1b", 0., 10. ) ;
      rv_mean_eff_sf_sb_1b      = new RooRealVar( "mean_eff_sf_sb_1b"     , "mean_eff_sf_sb_1b", 0., 10. ) ;
      rv_mean_eff_sf_sig_sl_1b  = new RooRealVar( "mean_eff_sf_sig_sl_1b" , "mean_eff_sf_sig_sl_1b", 0., 10. ) ;
      rv_mean_eff_sf_sb_sl_1b   = new RooRealVar( "mean_eff_sf_sb_sl_1b"  , "mean_eff_sf_sb_sl_1b", 0., 10. ) ;
      rv_mean_eff_sf_sig_ldp_1b = new RooRealVar( "mean_eff_sf_sig_ldp_1b", "mean_eff_sf_sig_ldp_1b", 0., 10. ) ;
      rv_mean_eff_sf_sb_ldp_1b  = new RooRealVar( "mean_eff_sf_sb_ldp_1b" , "mean_eff_sf_sb_ldp_1b", 0., 10. ) ;

      rv_mean_eff_sf_sig_2b     = new RooRealVar( "mean_eff_sf_sig_2b"    , "mean_eff_sf_sig_2b", 0., 10. ) ;
      rv_mean_eff_sf_sb_2b      = new RooRealVar( "mean_eff_sf_sb_2b"     , "mean_eff_sf_sb_2b", 0., 10. ) ;
      rv_mean_eff_sf_sig_sl_2b  = new RooRealVar( "mean_eff_sf_sig_sl_2b" , "mean_eff_sf_sig_sl_2b", 0., 10. ) ;
      rv_mean_eff_sf_sb_sl_2b   = new RooRealVar( "mean_eff_sf_sb_sl_2b"  , "mean_eff_sf_sb_sl_2b", 0., 10. ) ;
      rv_mean_eff_sf_sig_ldp_2b = new RooRealVar( "mean_eff_sf_sig_ldp_2b", "mean_eff_sf_sig_ldp_2b", 0., 10. ) ;
      rv_mean_eff_sf_sb_ldp_2b  = new RooRealVar( "mean_eff_sf_sb_ldp_2b" , "mean_eff_sf_sb_ldp_2b", 0., 10. ) ;

      rv_mean_eff_sf_sig_3b     = new RooRealVar( "mean_eff_sf_sig_3b"    , "mean_eff_sf_sig_3b", 0., 10. ) ;
      rv_mean_eff_sf_sb_3b      = new RooRealVar( "mean_eff_sf_sb_3b"     , "mean_eff_sf_sb_3b", 0., 10. ) ;
      rv_mean_eff_sf_sig_sl_3b  = new RooRealVar( "mean_eff_sf_sig_sl_3b" , "mean_eff_sf_sig_sl_3b", 0., 10. ) ;
      rv_mean_eff_sf_sb_sl_3b   = new RooRealVar( "mean_eff_sf_sb_sl_3b"  , "mean_eff_sf_sb_sl_3b", 0., 10. ) ;
      rv_mean_eff_sf_sig_ldp_3b = new RooRealVar( "mean_eff_sf_sig_ldp_3b", "mean_eff_sf_sig_ldp_3b", 0., 10. ) ;
      rv_mean_eff_sf_sb_ldp_3b  = new RooRealVar( "mean_eff_sf_sb_ldp_3b" , "mean_eff_sf_sb_ldp_3b", 0., 10. ) ;



    //--- width parameters.
      rv_width_eff_sf_sig_1b     = new RooRealVar( "width_eff_sf_sig_1b"    , "width_eff_sf_sig_1b", 0., 10. ) ;
      rv_width_eff_sf_sb_1b      = new RooRealVar( "width_eff_sf_sb_1b"     , "width_eff_sf_sb_1b", 0., 10. ) ;
      rv_width_eff_sf_sig_sl_1b  = new RooRealVar( "width_eff_sf_sig_sl_1b" , "width_eff_sf_sig_sl_1b", 0., 10. ) ;
      rv_width_eff_sf_sb_sl_1b   = new RooRealVar( "width_eff_sf_sb_sl_1b"  , "width_eff_sf_sb_sl_1b", 0., 10. ) ;
      rv_width_eff_sf_sig_ldp_1b = new RooRealVar( "width_eff_sf_sig_ldp_1b", "width_eff_sf_sig_ldp_1b", 0., 10. ) ;
      rv_width_eff_sf_sb_ldp_1b  = new RooRealVar( "width_eff_sf_sb_ldp_1b" , "width_eff_sf_sb_ldp_1b", 0., 10. ) ;

      rv_width_eff_sf_sig_2b     = new RooRealVar( "width_eff_sf_sig_2b"    , "width_eff_sf_sig_2b", 0., 10. ) ;
      rv_width_eff_sf_sb_2b      = new RooRealVar( "width_eff_sf_sb_2b"     , "width_eff_sf_sb_2b", 0., 10. ) ;
      rv_width_eff_sf_sig_sl_2b  = new RooRealVar( "width_eff_sf_sig_sl_2b" , "width_eff_sf_sig_sl_2b", 0., 10. ) ;
      rv_width_eff_sf_sb_sl_2b   = new RooRealVar( "width_eff_sf_sb_sl_2b"  , "width_eff_sf_sb_sl_2b", 0., 10. ) ;
      rv_width_eff_sf_sig_ldp_2b = new RooRealVar( "width_eff_sf_sig_ldp_2b", "width_eff_sf_sig_ldp_2b", 0., 10. ) ;
      rv_width_eff_sf_sb_ldp_2b  = new RooRealVar( "width_eff_sf_sb_ldp_2b" , "width_eff_sf_sb_ldp_2b", 0., 10. ) ;

      rv_width_eff_sf_sig_3b     = new RooRealVar( "width_eff_sf_sig_3b"    , "width_eff_sf_sig_3b", 0., 10. ) ;
      rv_width_eff_sf_sb_3b      = new RooRealVar( "width_eff_sf_sb_3b"     , "width_eff_sf_sb_3b", 0., 10. ) ;
      rv_width_eff_sf_sig_sl_3b  = new RooRealVar( "width_eff_sf_sig_sl_3b" , "width_eff_sf_sig_sl_3b", 0., 10. ) ;
      rv_width_eff_sf_sb_sl_3b   = new RooRealVar( "width_eff_sf_sb_sl_3b"  , "width_eff_sf_sb_sl_3b", 0., 10. ) ;
      rv_width_eff_sf_sig_ldp_3b = new RooRealVar( "width_eff_sf_sig_ldp_3b", "width_eff_sf_sig_ldp_3b", 0., 10. ) ;
      rv_width_eff_sf_sb_ldp_3b  = new RooRealVar( "width_eff_sf_sb_ldp_3b" , "width_eff_sf_sb_ldp_3b", 0., 10. ) ;


      rv_mean_eff_sf_sig_1b     -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_1b      -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_sl_1b  -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_sl_1b   -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_ldp_1b -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_ldp_1b  -> setVal( 1.00 ) ;

      rv_mean_eff_sf_sig_2b     -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_2b      -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_sl_2b  -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_sl_2b   -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_ldp_2b -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_ldp_2b  -> setVal( 1.00 ) ;

      rv_mean_eff_sf_sig_3b     -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_3b      -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_sl_3b  -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_sl_3b   -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sig_ldp_3b -> setVal( 1.00 ) ;
      rv_mean_eff_sf_sb_ldp_3b  -> setVal( 1.00 ) ;



      rv_mean_eff_sf_sig_1b     -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_1b      -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_sl_1b  -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_sl_1b   -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_ldp_1b -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_ldp_1b  -> setConstant( kTRUE ) ;

      rv_mean_eff_sf_sig_2b     -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_2b      -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_sl_2b  -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_sl_2b   -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_ldp_2b -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_ldp_2b  -> setConstant( kTRUE ) ;

      rv_mean_eff_sf_sig_3b     -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_3b      -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_sl_3b  -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_sl_3b   -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sig_ldp_3b -> setConstant( kTRUE ) ;
      rv_mean_eff_sf_sb_ldp_3b  -> setConstant( kTRUE ) ;





      //--- Initialize all width parameters to 15%.  They will be reset in the susy scan.
      rv_width_eff_sf_sig_1b     -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_1b      -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_sl_1b  -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_sl_1b   -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_ldp_1b -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_ldp_1b  -> setVal( 0.15 ) ;

      rv_width_eff_sf_sig_2b     -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_2b      -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_sl_2b  -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_sl_2b   -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_ldp_2b -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_ldp_2b  -> setVal( 0.15 ) ;

      rv_width_eff_sf_sig_3b     -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_3b      -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_sl_3b  -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_sl_3b   -> setVal( 0.15 ) ;
      rv_width_eff_sf_sig_ldp_3b -> setVal( 0.15 ) ;
      rv_width_eff_sf_sb_ldp_3b  -> setVal( 0.15 ) ;



      rv_width_eff_sf_sig_1b     -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_1b      -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_sl_1b  -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_sl_1b   -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_ldp_1b -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_ldp_1b  -> setConstant( kTRUE ) ;

      rv_width_eff_sf_sig_2b     -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_2b      -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_sl_2b  -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_sl_2b   -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_ldp_2b -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_ldp_2b  -> setConstant( kTRUE ) ;

      rv_width_eff_sf_sig_3b     -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_3b      -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_sl_3b  -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_sl_3b   -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sig_ldp_3b -> setConstant( kTRUE ) ;
      rv_width_eff_sf_sb_ldp_3b  -> setConstant( kTRUE ) ;



    //--- Jan 10, 2012: btag efficiency derivatives.

      rv_deff_dbtageff_sig_1b     = new RooRealVar( "deff_dbtageff_sig_1b"    , "deff_dbtageff_sig_1b", -10., 10. ) ;
      rv_deff_dbtageff_sb_1b      = new RooRealVar( "deff_dbtageff_sb_1b"     , "deff_dbtageff_sb_1b", -10., 10. ) ;
      rv_deff_dbtageff_sig_sl_1b  = new RooRealVar( "deff_dbtageff_sig_sl_1b" , "deff_dbtageff_sig_sl_1b", -10., 10. ) ;
      rv_deff_dbtageff_sb_sl_1b   = new RooRealVar( "deff_dbtageff_sb_sl_1b"  , "deff_dbtageff_sb_sl_1b", -10., 10. ) ;
      rv_deff_dbtageff_sig_ldp_1b = new RooRealVar( "deff_dbtageff_sig_ldp_1b", "deff_dbtageff_sig_ldp_1b", -10., 10. ) ;
      rv_deff_dbtageff_sb_ldp_1b  = new RooRealVar( "deff_dbtageff_sb_ldp_1b" , "deff_dbtageff_sb_ldp_1b", -10., 10. ) ;

      rv_deff_dbtageff_sig_2b     = new RooRealVar( "deff_dbtageff_sig_2b"    , "deff_dbtageff_sig_2b", -10., 10. ) ;
      rv_deff_dbtageff_sb_2b      = new RooRealVar( "deff_dbtageff_sb_2b"     , "deff_dbtageff_sb_2b", -10., 10. ) ;
      rv_deff_dbtageff_sig_sl_2b  = new RooRealVar( "deff_dbtageff_sig_sl_2b" , "deff_dbtageff_sig_sl_2b", -10., 10. ) ;
      rv_deff_dbtageff_sb_sl_2b   = new RooRealVar( "deff_dbtageff_sb_sl_2b"  , "deff_dbtageff_sb_sl_2b", -10., 10. ) ;
      rv_deff_dbtageff_sig_ldp_2b = new RooRealVar( "deff_dbtageff_sig_ldp_2b", "deff_dbtageff_sig_ldp_2b", -10., 10. ) ;
      rv_deff_dbtageff_sb_ldp_2b  = new RooRealVar( "deff_dbtageff_sb_ldp_2b" , "deff_dbtageff_sb_ldp_2b", -10., 10. ) ;

      rv_deff_dbtageff_sig_3b     = new RooRealVar( "deff_dbtageff_sig_3b"    , "deff_dbtageff_sig_3b", -10., 10. ) ;
      rv_deff_dbtageff_sb_3b      = new RooRealVar( "deff_dbtageff_sb_3b"     , "deff_dbtageff_sb_3b", -10., 10. ) ;
      rv_deff_dbtageff_sig_sl_3b  = new RooRealVar( "deff_dbtageff_sig_sl_3b" , "deff_dbtageff_sig_sl_3b", -10., 10. ) ;
      rv_deff_dbtageff_sb_sl_3b   = new RooRealVar( "deff_dbtageff_sb_sl_3b"  , "deff_dbtageff_sb_sl_3b", -10., 10. ) ;
      rv_deff_dbtageff_sig_ldp_3b = new RooRealVar( "deff_dbtageff_sig_ldp_3b", "deff_dbtageff_sig_ldp_3b", -10., 10. ) ;
      rv_deff_dbtageff_sb_ldp_3b  = new RooRealVar( "deff_dbtageff_sb_ldp_3b" , "deff_dbtageff_sb_ldp_3b", -10., 10. ) ;


      rv_deff_dbtageff_sig_1b     -> setVal( -1.0 ) ;
      rv_deff_dbtageff_sb_1b      -> setVal( -1.0 ) ;
      rv_deff_dbtageff_sig_sl_1b  -> setVal( -1.0 ) ;
      rv_deff_dbtageff_sb_sl_1b   -> setVal( -1.0 ) ;
      rv_deff_dbtageff_sig_ldp_1b -> setVal( -1.0 ) ;
      rv_deff_dbtageff_sb_ldp_1b  -> setVal( -1.0 ) ;

      rv_deff_dbtageff_sig_2b     -> setVal( 0.5 ) ;
      rv_deff_dbtageff_sb_2b      -> setVal( 0.5 ) ;
      rv_deff_dbtageff_sig_sl_2b  -> setVal( 0.5 ) ;
      rv_deff_dbtageff_sb_sl_2b   -> setVal( 0.5 ) ;
      rv_deff_dbtageff_sig_ldp_2b -> setVal( 0.5 ) ;
      rv_deff_dbtageff_sb_ldp_2b  -> setVal( 0.5 ) ;

      rv_deff_dbtageff_sig_3b     -> setVal( 2.0 ) ;
      rv_deff_dbtageff_sb_3b      -> setVal( 2.0 ) ;
      rv_deff_dbtageff_sig_sl_3b  -> setVal( 2.0 ) ;
      rv_deff_dbtageff_sb_sl_3b   -> setVal( 2.0 ) ;
      rv_deff_dbtageff_sig_ldp_3b -> setVal( 2.0 ) ;
      rv_deff_dbtageff_sb_ldp_3b  -> setVal( 2.0 ) ;



      rv_deff_dbtageff_sig_1b     -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sb_1b      -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sig_sl_1b  -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sb_sl_1b   -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sig_ldp_1b -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sb_ldp_1b  -> setConstant( kTRUE ) ;

      rv_deff_dbtageff_sig_2b     -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sb_2b      -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sig_sl_2b  -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sb_sl_2b   -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sig_ldp_2b -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sb_ldp_2b  -> setConstant( kTRUE ) ;

      rv_deff_dbtageff_sig_3b     -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sb_3b      -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sig_sl_3b  -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sb_sl_3b   -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sig_ldp_3b -> setConstant( kTRUE ) ;
      rv_deff_dbtageff_sb_ldp_3b  -> setConstant( kTRUE ) ;




   //
   // Owen : Sept 7, 2011: need to set these here so that it gets into the workspace.
   //
   //++++++++++
      setSusyScanPoint( inputScanFile,  m0,  m12,  isT1bbbb,  t1bbbbXsec, inputSusy_deff_dbtageff_file ) ;
   //++++++++++





      RooRealVar rv_btageff_err ( "btageff_err", "btageff_err", 0., 10. ) ;
      rv_btageff_err. setVal( btageff_err ) ;
      rv_btageff_err. setConstant( kTRUE ) ;





    //--- Systematics and other nuisance parameters
      char formula[1024];
      RooArgSet globalObservables ("globalObservables");
      RooArgSet allNuisances ("allNuisances");
      RooArgSet allNuisancePdfs ("allNuisancePdfs");

      RooRealVar Rlsb_passfail_prim ( "Rlsb_passfail_prim", "Rlsb_passfail_prim", 0, -5, 5);
      RooRealVar Rlsb_passfail_nom ( "Rlsb_passfail_nom", "Rlsb_passfail_nom", 0, -5, 5);
      RooGaussian pdf_Rlsb_passfail ("pdf_Rlsb_passfail" , "pdf_Rlsb_passfail", Rlsb_passfail_prim, Rlsb_passfail_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", Rlsb_passfail, exp(Rlsb_passfail_err/Rlsb_passfail));
      RooFormulaVar fv_Rlsb_passfail ("Rlsb_passfail", formula, RooArgList(Rlsb_passfail_prim));
      Rlsb_passfail_nom.setConstant();
      globalObservables.add (Rlsb_passfail_nom);
      allNuisances.add (Rlsb_passfail_prim);
      allNuisancePdfs.add (pdf_Rlsb_passfail);

      RooRealVar acc_ee_sig_prim ( "acc_ee_sig_prim", "acc_ee_sig_prim", 0, -5, 5);
      RooRealVar acc_ee_sig_nom ( "acc_ee_sig_nom", "acc_ee_sig_nom", 0, -5, 5);
      RooGaussian pdf_acc_ee_sig ("pdf_acc_ee_sig" , "pdf_acc_ee_sig", acc_ee_sig_prim, acc_ee_sig_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", acc_ee_sig_mean, exp(acc_ee_sig_err/acc_ee_sig_mean));
      RooFormulaVar fv_acc_ee_sig ("acc_ee_sig", formula, RooArgList(acc_ee_sig_prim));
      acc_ee_sig_nom.setConstant();
      globalObservables.add (acc_ee_sig_nom);
      allNuisances.add (acc_ee_sig_prim);
      allNuisancePdfs.add (pdf_acc_ee_sig);

      RooRealVar acc_ee_sb_prim ( "acc_ee_sb_prim", "acc_ee_sb_prim", 0, -5, 5);
      RooRealVar acc_ee_sb_nom ( "acc_ee_sb_nom", "acc_ee_sb_nom", 0, -5, 5);
      RooGaussian pdf_acc_ee_sb ("pdf_acc_ee_sb" , "pdf_acc_ee_sb", acc_ee_sb_prim, acc_ee_sb_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", acc_ee_sb_mean, exp(acc_ee_sb_err/acc_ee_sb_mean));
      RooFormulaVar fv_acc_ee_sb ("acc_ee_sb", formula, RooArgList(acc_ee_sb_prim));
      acc_ee_sb_nom.setConstant();
      globalObservables.add (acc_ee_sb_nom);
      allNuisances.add (acc_ee_sb_prim);
      allNuisancePdfs.add (pdf_acc_ee_sb);

      RooRealVar acc_mm_sig_prim ( "acc_mm_sig_prim", "acc_mm_sig_prim", 0, -5, 5);
      RooRealVar acc_mm_sig_nom ( "acc_mm_sig_nom", "acc_mm_sig_nom", 0, -5, 5);
      RooGaussian pdf_acc_mm_sig ("pdf_acc_mm_sig" , "pdf_acc_mm_sig", acc_mm_sig_prim, acc_mm_sig_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", acc_mm_sig_mean, exp(acc_mm_sig_err/acc_mm_sig_mean));
      RooFormulaVar fv_acc_mm_sig ("acc_mm_sig", formula, RooArgList(acc_mm_sig_prim));
      acc_mm_sig_nom.setConstant();
      globalObservables.add (acc_mm_sig_nom);
      allNuisances.add (acc_mm_sig_prim);
      allNuisancePdfs.add (pdf_acc_mm_sig);

      RooRealVar acc_mm_sb_prim ( "acc_mm_sb_prim", "acc_mm_sb_prim", 0, -5, 5);
      RooRealVar acc_mm_sb_nom ( "acc_mm_sb_nom", "acc_mm_sb_nom", 0, -5, 5);
      RooGaussian pdf_acc_mm_sb ("pdf_acc_mm_sb" , "pdf_acc_mm_sb", acc_mm_sb_prim, acc_mm_sb_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", acc_mm_sb_mean, exp(acc_mm_sb_err/acc_mm_sb_mean));
      RooFormulaVar fv_acc_mm_sb ("acc_mm_sb", formula, RooArgList(acc_mm_sb_prim));
      acc_mm_sb_nom.setConstant();
      globalObservables.add (acc_mm_sb_nom);
      allNuisances.add (acc_mm_sb_prim);
      allNuisancePdfs.add (pdf_acc_mm_sb);


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


    //--- Nov 24, 2011:  sf_ee and sf_mm are now derived from a common underlying Gaussian.

      RooRealVar sf_ll_prim ( "sf_ll_prim", "sf_ll_prim", 0, -5, 5);
      RooRealVar sf_ll_nom ( "sf_ll_nom", "sf_ll_nom", 0, -5, 5);
      RooGaussian pdf_sf_ll ("pdf_sf_ll" , "pdf_sf_ll", sf_ll_prim, sf_ll_nom, RooConst(1));
      sf_ll_nom.setConstant();
      globalObservables.add (sf_ll_nom);
      allNuisances.add (sf_ll_prim);
      allNuisancePdfs.add (pdf_sf_ll);

      sprintf (formula, "%f*pow(%f,@0)", sf_ee, exp(sf_ee_err/sf_ee));
      RooFormulaVar fv_sf_ee ("sf_ee", formula, RooArgList(sf_ll_prim));

      sprintf (formula, "%f*pow(%f,@0)", sf_mm, exp(sf_mm_err/sf_mm));
      RooFormulaVar fv_sf_mm ("sf_mm", formula, RooArgList(sf_ll_prim));










      RooRealVar sf_mc_prim ( "sf_mc_prim", "sf_mc_prim", 0, -5, 5);
      RooRealVar sf_mc_nom ( "sf_mc_nom", "sf_mc_nom", 0, -5, 5);
      RooGaussian pdf_sf_mc ("pdf_sf_mc" , "pdf_sf_mc", sf_mc_prim, sf_mc_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_mc, exp(sf_mc_err/sf_mc));
      RooFormulaVar fv_sf_mc ("sf_mc", formula, RooArgList(sf_mc_prim));
      sf_mc_nom.setConstant();
      globalObservables.add (sf_mc_nom);
      allNuisances.add (sf_mc_prim);
      allNuisancePdfs.add (pdf_sf_mc);








      RooRealVar sf_qcd_sb_1b_prim ( "sf_qcd_sb_1b_prim", "sf_qcd_sb_1b_prim", 0, -5, 5);
      RooRealVar sf_qcd_sb_1b_nom ( "sf_qcd_sb_1b_nom", "sf_qcd_sb_1b_nom", 0, -5, 5);
      RooGaussian pdf_sf_qcd_sb_1b ("pdf_sf_qcd_sb_1b" , "pdf_sf_qcd_sb_1b", sf_qcd_sb_1b_prim, sf_qcd_sb_1b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_qcd_sb_1b, exp(sf_qcd_sb_1b_err/sf_qcd_sb_1b));
      RooFormulaVar fv_sf_qcd_sb_1b ("sf_qcd_sb_1b", formula, RooArgList(sf_qcd_sb_1b_prim));
      sf_qcd_sb_1b_nom.setConstant();
      globalObservables.add (sf_qcd_sb_1b_nom);
      allNuisances.add (sf_qcd_sb_1b_prim);
      allNuisancePdfs.add (pdf_sf_qcd_sb_1b);


      RooRealVar sf_qcd_sb_2b_prim ( "sf_qcd_sb_2b_prim", "sf_qcd_sb_2b_prim", 0, -5, 5);
      RooRealVar sf_qcd_sb_2b_nom ( "sf_qcd_sb_2b_nom", "sf_qcd_sb_2b_nom", 0, -5, 5);
      RooGaussian pdf_sf_qcd_sb_2b ("pdf_sf_qcd_sb_2b" , "pdf_sf_qcd_sb_2b", sf_qcd_sb_2b_prim, sf_qcd_sb_2b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_qcd_sb_2b, exp(sf_qcd_sb_2b_err/sf_qcd_sb_2b));
      RooFormulaVar fv_sf_qcd_sb_2b ("sf_qcd_sb_2b", formula, RooArgList(sf_qcd_sb_2b_prim));
      sf_qcd_sb_2b_nom.setConstant();
      globalObservables.add (sf_qcd_sb_2b_nom);
      allNuisances.add (sf_qcd_sb_2b_prim);
      allNuisancePdfs.add (pdf_sf_qcd_sb_2b);


      RooRealVar sf_qcd_sb_3b_prim ( "sf_qcd_sb_3b_prim", "sf_qcd_sb_3b_prim", 0, -5, 5);
      RooRealVar sf_qcd_sb_3b_nom ( "sf_qcd_sb_3b_nom", "sf_qcd_sb_3b_nom", 0, -5, 5);
      RooGaussian pdf_sf_qcd_sb_3b ("pdf_sf_qcd_sb_3b" , "pdf_sf_qcd_sb_3b", sf_qcd_sb_3b_prim, sf_qcd_sb_3b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_qcd_sb_3b, exp(sf_qcd_sb_3b_err/sf_qcd_sb_3b));
      RooFormulaVar fv_sf_qcd_sb_3b ("sf_qcd_sb_3b", formula, RooArgList(sf_qcd_sb_3b_prim));
      sf_qcd_sb_3b_nom.setConstant();
      globalObservables.add (sf_qcd_sb_3b_nom);
      allNuisances.add (sf_qcd_sb_3b_prim);
      allNuisancePdfs.add (pdf_sf_qcd_sb_3b);








      RooRealVar sf_qcd_sig_1b_prim ( "sf_qcd_sig_1b_prim", "sf_qcd_sig_1b_prim", 0, -5, 5);
      RooRealVar sf_qcd_sig_1b_nom ( "sf_qcd_sig_1b_nom", "sf_qcd_sig_1b_nom", 0, -5, 5);
      RooGaussian pdf_sf_qcd_sig_1b ("pdf_sf_qcd_sig_1b" , "pdf_sf_qcd_sig_1b", sf_qcd_sig_1b_prim, sf_qcd_sig_1b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_qcd_sig_1b, exp(sf_qcd_sig_1b_err/sf_qcd_sig_1b));
      RooFormulaVar fv_sf_qcd_sig_1b ("sf_qcd_sig_1b", formula, RooArgList(sf_qcd_sig_1b_prim));
      sf_qcd_sig_1b_nom.setConstant();
      globalObservables.add (sf_qcd_sig_1b_nom);
      allNuisances.add (sf_qcd_sig_1b_prim);
      allNuisancePdfs.add (pdf_sf_qcd_sig_1b);


      RooRealVar sf_qcd_sig_2b_prim ( "sf_qcd_sig_2b_prim", "sf_qcd_sig_2b_prim", 0, -5, 5);
      RooRealVar sf_qcd_sig_2b_nom ( "sf_qcd_sig_2b_nom", "sf_qcd_sig_2b_nom", 0, -5, 5);
      RooGaussian pdf_sf_qcd_sig_2b ("pdf_sf_qcd_sig_2b" , "pdf_sf_qcd_sig_2b", sf_qcd_sig_2b_prim, sf_qcd_sig_2b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_qcd_sig_2b, exp(sf_qcd_sig_2b_err/sf_qcd_sig_2b));
      RooFormulaVar fv_sf_qcd_sig_2b ("sf_qcd_sig_2b", formula, RooArgList(sf_qcd_sig_2b_prim));
      sf_qcd_sig_2b_nom.setConstant();
      globalObservables.add (sf_qcd_sig_2b_nom);
      allNuisances.add (sf_qcd_sig_2b_prim);
      allNuisancePdfs.add (pdf_sf_qcd_sig_2b);


      RooRealVar sf_qcd_sig_3b_prim ( "sf_qcd_sig_3b_prim", "sf_qcd_sig_3b_prim", 0, -5, 5);
      RooRealVar sf_qcd_sig_3b_nom ( "sf_qcd_sig_3b_nom", "sf_qcd_sig_3b_nom", 0, -5, 5);
      RooGaussian pdf_sf_qcd_sig_3b ("pdf_sf_qcd_sig_3b" , "pdf_sf_qcd_sig_3b", sf_qcd_sig_3b_prim, sf_qcd_sig_3b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_qcd_sig_3b, exp(sf_qcd_sig_3b_err/sf_qcd_sig_3b));
      RooFormulaVar fv_sf_qcd_sig_3b ("sf_qcd_sig_3b", formula, RooArgList(sf_qcd_sig_3b_prim));
      sf_qcd_sig_3b_nom.setConstant();
      globalObservables.add (sf_qcd_sig_3b_nom);
      allNuisances.add (sf_qcd_sig_3b_prim);
      allNuisancePdfs.add (pdf_sf_qcd_sig_3b);














      RooRealVar sf_ttwj_sig_1b_prim ( "sf_ttwj_sig_1b_prim", "sf_ttwj_sig_1b_prim", 0, -5, 5);
      RooRealVar sf_ttwj_sig_1b_nom ( "sf_ttwj_sig_1b_nom", "sf_ttwj_sig_1b_nom", 0, -5, 5);
      RooGaussian pdf_sf_ttwj_sig_1b ("pdf_sf_ttwj_sig_1b" , "pdf_sf_ttwj_sig_1b", sf_ttwj_sig_1b_prim, sf_ttwj_sig_1b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_ttwj_sig_1b, exp(sf_ttwj_sig_1b_err/sf_ttwj_sig_1b));
      RooFormulaVar fv_sf_ttwj_sig_1b ("sf_ttwj_sig_1b", formula, RooArgList(sf_ttwj_sig_1b_prim));
      sf_ttwj_sig_1b_nom.setConstant();
      globalObservables.add (sf_ttwj_sig_1b_nom);
      allNuisances.add (sf_ttwj_sig_1b_prim);
      allNuisancePdfs.add (pdf_sf_ttwj_sig_1b);


      RooRealVar sf_ttwj_sig_2b_prim ( "sf_ttwj_sig_2b_prim", "sf_ttwj_sig_2b_prim", 0, -5, 5);
      RooRealVar sf_ttwj_sig_2b_nom ( "sf_ttwj_sig_2b_nom", "sf_ttwj_sig_2b_nom", 0, -5, 5);
      RooGaussian pdf_sf_ttwj_sig_2b ("pdf_sf_ttwj_sig_2b" , "pdf_sf_ttwj_sig_2b", sf_ttwj_sig_2b_prim, sf_ttwj_sig_2b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_ttwj_sig_2b, exp(sf_ttwj_sig_2b_err/sf_ttwj_sig_2b));
      RooFormulaVar fv_sf_ttwj_sig_2b ("sf_ttwj_sig_2b", formula, RooArgList(sf_ttwj_sig_2b_prim));
      sf_ttwj_sig_2b_nom.setConstant();
      globalObservables.add (sf_ttwj_sig_2b_nom);
      allNuisances.add (sf_ttwj_sig_2b_prim);
      allNuisancePdfs.add (pdf_sf_ttwj_sig_2b);


      RooRealVar sf_ttwj_sig_3b_prim ( "sf_ttwj_sig_3b_prim", "sf_ttwj_sig_3b_prim", 0, -5, 5);
      RooRealVar sf_ttwj_sig_3b_nom ( "sf_ttwj_sig_3b_nom", "sf_ttwj_sig_3b_nom", 0, -5, 5);
      RooGaussian pdf_sf_ttwj_sig_3b ("pdf_sf_ttwj_sig_3b" , "pdf_sf_ttwj_sig_3b", sf_ttwj_sig_3b_prim, sf_ttwj_sig_3b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_ttwj_sig_3b, exp(sf_ttwj_sig_3b_err/sf_ttwj_sig_3b));
      RooFormulaVar fv_sf_ttwj_sig_3b ("sf_ttwj_sig_3b", formula, RooArgList(sf_ttwj_sig_3b_prim));
      sf_ttwj_sig_3b_nom.setConstant();
      globalObservables.add (sf_ttwj_sig_3b_nom);
      allNuisances.add (sf_ttwj_sig_3b_prim);
      allNuisancePdfs.add (pdf_sf_ttwj_sig_3b);


      RooRealVar sf_ttwj_sb_2b_prim ( "sf_ttwj_sb_2b_prim", "sf_ttwj_sb_2b_prim", 0, -5, 5);
      RooRealVar sf_ttwj_sb_2b_nom ( "sf_ttwj_sb_2b_nom", "sf_ttwj_sb_2b_nom", 0, -5, 5);
      RooGaussian pdf_sf_ttwj_sb_2b ("pdf_sf_ttwj_sb_2b" , "pdf_sf_ttwj_sb_2b", sf_ttwj_sb_2b_prim, sf_ttwj_sb_2b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_ttwj_sb_2b, exp(sf_ttwj_sb_2b_err/sf_ttwj_sb_2b));
      RooFormulaVar fv_sf_ttwj_sb_2b ("sf_ttwj_sb_2b", formula, RooArgList(sf_ttwj_sb_2b_prim));
      sf_ttwj_sb_2b_nom.setConstant();
      globalObservables.add (sf_ttwj_sb_2b_nom);
      allNuisances.add (sf_ttwj_sb_2b_prim);
      allNuisancePdfs.add (pdf_sf_ttwj_sb_2b);


      RooRealVar sf_ttwj_sb_3b_prim ( "sf_ttwj_sb_3b_prim", "sf_ttwj_sb_3b_prim", 0, -5, 5);
      RooRealVar sf_ttwj_sb_3b_nom ( "sf_ttwj_sb_3b_nom", "sf_ttwj_sb_3b_nom", 0, -5, 5);
      RooGaussian pdf_sf_ttwj_sb_3b ("pdf_sf_ttwj_sb_3b" , "pdf_sf_ttwj_sb_3b", sf_ttwj_sb_3b_prim, sf_ttwj_sb_3b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", sf_ttwj_sb_3b, exp(sf_ttwj_sb_3b_err/sf_ttwj_sb_3b));
      RooFormulaVar fv_sf_ttwj_sb_3b ("sf_ttwj_sb_3b", formula, RooArgList(sf_ttwj_sb_3b_prim));
      sf_ttwj_sb_3b_nom.setConstant();
      globalObservables.add (sf_ttwj_sb_3b_nom);
      allNuisances.add (sf_ttwj_sb_3b_prim);
      allNuisancePdfs.add (pdf_sf_ttwj_sb_3b);




   //--------------------------------------

      RooRealVar knn_sig_1b_prim ( "knn_sig_1b_prim", "knn_sig_1b_prim", 0, -5, 5);
      RooRealVar knn_sig_1b_nom ( "knn_sig_1b_nom", "knn_sig_1b_nom", 0, -5, 5);
      RooGaussian pdf_knn_sig_1b ("pdf_knn_sig_1b" , "pdf_knn_sig_1b", knn_sig_1b_prim, knn_sig_1b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", knn_sig_1b_mean, exp(knn_sig_1b_err/knn_sig_1b_mean));
      RooFormulaVar fv_knn_sig_1b ("knn_sig_1b", formula, RooArgList(knn_sig_1b_prim));
      knn_sig_1b_nom.setConstant();
      globalObservables.add (knn_sig_1b_nom);
      allNuisances.add (knn_sig_1b_prim);
      allNuisancePdfs.add (pdf_knn_sig_1b);


      RooRealVar knn_sig_2b_prim ( "knn_sig_2b_prim", "knn_sig_2b_prim", 0, -5, 5);
      RooRealVar knn_sig_2b_nom ( "knn_sig_2b_nom", "knn_sig_2b_nom", 0, -5, 5);
      RooGaussian pdf_knn_sig_2b ("pdf_knn_sig_2b" , "pdf_knn_sig_2b", knn_sig_2b_prim, knn_sig_2b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", knn_sig_2b_mean, exp(knn_sig_2b_err/knn_sig_2b_mean));
      RooFormulaVar fv_knn_sig_2b ("knn_sig_2b", formula, RooArgList(knn_sig_2b_prim));
      knn_sig_2b_nom.setConstant();
      globalObservables.add (knn_sig_2b_nom);
      allNuisances.add (knn_sig_2b_prim);
      allNuisancePdfs.add (pdf_knn_sig_2b);


      RooRealVar knn_sig_3b_prim ( "knn_sig_3b_prim", "knn_sig_3b_prim", 0, -5, 5);
      RooRealVar knn_sig_3b_nom ( "knn_sig_3b_nom", "knn_sig_3b_nom", 0, -5, 5);
      RooGaussian pdf_knn_sig_3b ("pdf_knn_sig_3b" , "pdf_knn_sig_3b", knn_sig_3b_prim, knn_sig_3b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", knn_sig_3b_mean, exp(knn_sig_3b_err/knn_sig_3b_mean));
      RooFormulaVar fv_knn_sig_3b ("knn_sig_3b", formula, RooArgList(knn_sig_3b_prim));
      knn_sig_3b_nom.setConstant();
      globalObservables.add (knn_sig_3b_nom);
      allNuisances.add (knn_sig_3b_prim);
      allNuisancePdfs.add (pdf_knn_sig_3b);







      RooRealVar knn_sb_1b_prim ( "knn_sb_1b_prim", "knn_sb_1b_prim", 0, -5, 5);
      RooRealVar knn_sb_1b_nom ( "knn_sb_1b_nom", "knn_sb_1b_nom", 0, -5, 5);
      RooGaussian pdf_knn_sb_1b ("pdf_knn_sb_1b" , "pdf_knn_sb_1b", knn_sb_1b_prim, knn_sb_1b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", knn_sb_1b_mean, exp(knn_sb_1b_err/knn_sb_1b_mean));
      RooFormulaVar fv_knn_sb_1b ("knn_sb_1b", formula, RooArgList(knn_sb_1b_prim));
      knn_sb_1b_nom.setConstant();
      globalObservables.add (knn_sb_1b_nom);
      allNuisances.add (knn_sb_1b_prim);
      allNuisancePdfs.add (pdf_knn_sb_1b);


      RooRealVar knn_sb_2b_prim ( "knn_sb_2b_prim", "knn_sb_2b_prim", 0, -5, 5);
      RooRealVar knn_sb_2b_nom ( "knn_sb_2b_nom", "knn_sb_2b_nom", 0, -5, 5);
      RooGaussian pdf_knn_sb_2b ("pdf_knn_sb_2b" , "pdf_knn_sb_2b", knn_sb_2b_prim, knn_sb_2b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", knn_sb_2b_mean, exp(knn_sb_2b_err/knn_sb_2b_mean));
      RooFormulaVar fv_knn_sb_2b ("knn_sb_2b", formula, RooArgList(knn_sb_2b_prim));
      knn_sb_2b_nom.setConstant();
      globalObservables.add (knn_sb_2b_nom);
      allNuisances.add (knn_sb_2b_prim);
      allNuisancePdfs.add (pdf_knn_sb_2b);


      RooRealVar knn_sb_3b_prim ( "knn_sb_3b_prim", "knn_sb_3b_prim", 0, -5, 5);
      RooRealVar knn_sb_3b_nom ( "knn_sb_3b_nom", "knn_sb_3b_nom", 0, -5, 5);
      RooGaussian pdf_knn_sb_3b ("pdf_knn_sb_3b" , "pdf_knn_sb_3b", knn_sb_3b_prim, knn_sb_3b_nom, RooConst(1));
      sprintf (formula, "%f*pow(%f,@0)", knn_sb_3b_mean, exp(knn_sb_3b_err/knn_sb_3b_mean));
      RooFormulaVar fv_knn_sb_3b ("knn_sb_3b", formula, RooArgList(knn_sb_3b_prim));
      knn_sb_3b_nom.setConstant();
      globalObservables.add (knn_sb_3b_nom);
      allNuisances.add (knn_sb_3b_prim);
      allNuisancePdfs.add (pdf_knn_sb_3b);










    //--------------------------------------


      RooRealVar eff_sf_prim ( "eff_sf_prim", "eff_sf_prim", 0, -5, 5);
      RooRealVar eff_sf_nom ( "eff_sf_nom", "eff_sf_nom", 0, -5, 5);
      RooGaussian pdf_eff_sf ("pdf_eff_sf" , "pdf_eff_sf", eff_sf_prim, eff_sf_nom, RooConst(1));
      eff_sf_nom.setConstant();
      globalObservables.add (eff_sf_nom);
      allNuisances.add (eff_sf_prim);
      allNuisancePdfs.add (pdf_eff_sf);



    //--------------------------------------


      RooRealVar btageff_sf_prim ( "btageff_sf_prim", "btageff_sf_prim", 0, -5, 5);
      RooRealVar btageff_sf_nom ( "btageff_sf_nom", "btageff_sf_nom", 0, -5, 5);
      RooGaussian pdf_btageff_sf ("pdf_btageff_sf" , "pdf_btageff_sf", btageff_sf_prim, btageff_sf_nom, RooConst(1));
      btageff_sf_nom.setConstant();
      globalObservables.add (btageff_sf_nom);
      allNuisances.add (btageff_sf_prim);
      allNuisancePdfs.add (pdf_btageff_sf);




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

         //-- Here, float 1b SIG bin of ttwj and use SL ratios for other 5 ttwj SIG and SB bins.
         //
         //  It is STRONGLY recommended NOT to use this option, except for making profile plots of the ttwj SIG variable.
         //  The reason is that the SIG-SL component can be small or even zero, so the ratio can blow up.
         //
         rfv_mu_ttwj_sb_1b = new RooFormulaVar( "mu_ttwj_sb_1b",
                                                       "mu_ttwj_sig_1b * (1.0/sf_ttwj_sig_1b) * (mu_ttwj_sb_sl_1b/mu_ttwj_sig_sl_1b)",
                                                       RooArgSet( *rv_mu_ttwj_sig_1b, fv_sf_ttwj_sig_1b, *rv_mu_ttwj_sb_sl_1b, *rv_mu_ttwj_sig_sl_1b ) ) ;

         rfv_mu_ttwj_sb_2b = new RooFormulaVar( "mu_ttwj_sb_2b",
                                                       "mu_ttwj_sig_1b * (1.0/sf_ttwj_sig_2b) * (mu_ttwj_sb_sl_2b/mu_ttwj_sig_sl_1b)",
                                                       RooArgSet( *rv_mu_ttwj_sig_1b, fv_sf_ttwj_sig_2b, *rv_mu_ttwj_sb_sl_2b, *rv_mu_ttwj_sig_sl_1b ) ) ;

         rfv_mu_ttwj_sb_3b = new RooFormulaVar( "mu_ttwj_sb_3b",
                                                       "mu_ttwj_sig_1b * (1.0/sf_ttwj_sig_3b) * (mu_ttwj_sb_sl_3b/mu_ttwj_sig_sl_1b)",
                                                       RooArgSet( *rv_mu_ttwj_sig_1b, fv_sf_ttwj_sig_3b, *rv_mu_ttwj_sb_sl_3b, *rv_mu_ttwj_sig_sl_1b ) ) ;

         rfv_mu_ttwj_sig_2b = new RooFormulaVar( "mu_ttwj_sig_2b",
                                                       "mu_ttwj_sig_1b * (1.0/sf_ttwj_sb_2b) * (mu_ttwj_sig_sl_2b/mu_ttwj_sig_sl_1b)",
                                                       RooArgSet( *rv_mu_ttwj_sig_1b, fv_sf_ttwj_sb_2b, *rv_mu_ttwj_sig_sl_2b, *rv_mu_ttwj_sig_sl_1b ) ) ;

         rfv_mu_ttwj_sig_3b = new RooFormulaVar( "mu_ttwj_sig_3b",
                                                       "mu_ttwj_sig_1b * (1.0/sf_ttwj_sb_3b) * (mu_ttwj_sig_sl_3b/mu_ttwj_sig_sl_1b)",
                                                       RooArgSet( *rv_mu_ttwj_sig_1b, fv_sf_ttwj_sb_3b, *rv_mu_ttwj_sig_sl_3b, *rv_mu_ttwj_sig_sl_1b ) ) ;



         rv_mu_ttwj_sb_1b = rfv_mu_ttwj_sb_1b ;
         rv_mu_ttwj_sb_2b = rfv_mu_ttwj_sb_2b ;
         rv_mu_ttwj_sb_3b = rfv_mu_ttwj_sb_3b ;
         rv_mu_ttwj_sig_2b = rfv_mu_ttwj_sig_2b ;
         rv_mu_ttwj_sig_3b = rfv_mu_ttwj_sig_3b ;


      } else {

         //-- Here, float 1b SB bin of ttwj and use SL ratios for other 5 ttwj SIG and SB bins.
         //
         //  It is STRONGLY recommended to use this option.  See comment above.
         //
         rfv_mu_ttwj_sig_1b = new RooFormulaVar( "mu_ttwj_sig_1b",
                                                       "mu_ttwj_sb_1b * sf_ttwj_sig_1b * (mu_ttwj_sig_sl_1b/mu_ttwj_sb_sl_1b)",
                                                       RooArgSet( *rv_mu_ttwj_sb_1b, fv_sf_ttwj_sig_1b, *rv_mu_ttwj_sig_sl_1b, *rv_mu_ttwj_sb_sl_1b ) ) ;

         rfv_mu_ttwj_sig_2b = new RooFormulaVar( "mu_ttwj_sig_2b",
                                                       "mu_ttwj_sb_1b * sf_ttwj_sig_2b * (mu_ttwj_sig_sl_2b/mu_ttwj_sb_sl_1b)",
                                                       RooArgSet( *rv_mu_ttwj_sb_1b, fv_sf_ttwj_sig_2b, *rv_mu_ttwj_sig_sl_2b, *rv_mu_ttwj_sb_sl_1b ) ) ;

         rfv_mu_ttwj_sig_3b = new RooFormulaVar( "mu_ttwj_sig_3b",
                                                       "mu_ttwj_sb_1b * sf_ttwj_sig_3b * (mu_ttwj_sig_sl_3b/mu_ttwj_sb_sl_1b)",
                                                       RooArgSet( *rv_mu_ttwj_sb_1b, fv_sf_ttwj_sig_3b, *rv_mu_ttwj_sig_sl_3b, *rv_mu_ttwj_sb_sl_1b ) ) ;

         rfv_mu_ttwj_sb_2b  = new RooFormulaVar( "mu_ttwj_sb_2b",
                                                       "mu_ttwj_sb_1b * sf_ttwj_sb_2b * (mu_ttwj_sb_sl_2b/mu_ttwj_sb_sl_1b)",
                                                       RooArgSet( *rv_mu_ttwj_sb_1b, fv_sf_ttwj_sb_2b, *rv_mu_ttwj_sb_sl_2b, *rv_mu_ttwj_sb_sl_1b ) ) ;

         rfv_mu_ttwj_sb_3b  = new RooFormulaVar( "mu_ttwj_sb_3b",
                                                       "mu_ttwj_sb_1b * sf_ttwj_sb_3b * (mu_ttwj_sb_sl_3b/mu_ttwj_sb_sl_1b)",
                                                       RooArgSet( *rv_mu_ttwj_sb_1b, fv_sf_ttwj_sb_3b, *rv_mu_ttwj_sb_sl_3b, *rv_mu_ttwj_sb_sl_1b ) ) ;

         rv_mu_ttwj_sig_1b = rfv_mu_ttwj_sig_1b ;
         rv_mu_ttwj_sig_2b = rfv_mu_ttwj_sig_2b ;
         rv_mu_ttwj_sig_3b = rfv_mu_ttwj_sig_3b ;
         rv_mu_ttwj_sb_2b  = rfv_mu_ttwj_sb_2b  ;
         rv_mu_ttwj_sb_3b  = rfv_mu_ttwj_sb_3b  ;

      }




      rv_mu_ttwj_sig_ldp_1b = new RooFormulaVar( "mu_ttwj_sig_ldp_1b",
                              "mu_ttbarsingletopzjetsmc_sig_ldp_1b + mu_WJmc_sig_ldp_1b",
                              RooArgSet( *rv_mu_ttbarsingletopzjetsmc_sig_ldp_1b,
                                         *rv_mu_WJmc_sig_ldp_1b ) ) ;

      rv_mu_ttwj_sig_ldp_2b = new RooFormulaVar( "mu_ttwj_sig_ldp_2b",
                              "mu_ttbarsingletopzjetsmc_sig_ldp_2b + mu_WJmc_sig_ldp_2b",
                              RooArgSet( *rv_mu_ttbarsingletopzjetsmc_sig_ldp_2b,
                                         *rv_mu_WJmc_sig_ldp_2b ) ) ;

      rv_mu_ttwj_sig_ldp_3b = new RooFormulaVar( "mu_ttwj_sig_ldp_3b",
                              "mu_ttbarsingletopzjetsmc_sig_ldp_3b + mu_WJmc_sig_ldp_3b",
                              RooArgSet( *rv_mu_ttbarsingletopzjetsmc_sig_ldp_3b,
                                         *rv_mu_WJmc_sig_ldp_3b ) ) ;





      rv_mu_ttwj_sb_ldp_1b = new RooFormulaVar( "mu_ttwj_sb_ldp_1b",
                              "mu_ttbarsingletopzjetsmc_sb_ldp_1b + mu_WJmc_sb_ldp_1b",
                              RooArgSet( *rv_mu_ttbarsingletopzjetsmc_sb_ldp_1b,
                                         *rv_mu_WJmc_sb_ldp_1b ) ) ;

      rv_mu_ttwj_sb_ldp_2b = new RooFormulaVar( "mu_ttwj_sb_ldp_2b",
                              "mu_ttbarsingletopzjetsmc_sb_ldp_2b + mu_WJmc_sb_ldp_2b",
                              RooArgSet( *rv_mu_ttbarsingletopzjetsmc_sb_ldp_2b,
                                         *rv_mu_WJmc_sb_ldp_2b ) ) ;

      rv_mu_ttwj_sb_ldp_3b = new RooFormulaVar( "mu_ttwj_sb_ldp_3b",
                              "mu_ttbarsingletopzjetsmc_sb_ldp_3b + mu_WJmc_sb_ldp_3b",
                              RooArgSet( *rv_mu_ttbarsingletopzjetsmc_sb_ldp_3b,
                                         *rv_mu_WJmc_sb_ldp_3b ) ) ;






    //-- QCD


      if ( useLdpVars ) {

         //--- Here, float the LSB variable and derive the SIG or SB variable using the LSB ratio.
         //
         //    This option is recommended because Minuit will have an easier time with this floating variable.
         //

         rfv_mu_qcd_sig_1b = new RooFormulaVar( "mu_qcd_sig_1b",
                                     "mu_qcd_sig_ldp_1b * sf_qcd_sig_1b * Rlsb_passfail",
                                     RooArgSet( *rv_mu_qcd_sig_ldp_1b, fv_sf_qcd_sig_1b, fv_Rlsb_passfail ) ) ;

         rfv_mu_qcd_sig_2b = new RooFormulaVar( "mu_qcd_sig_2b",
                                     "mu_qcd_sig_ldp_2b * sf_qcd_sig_2b * Rlsb_passfail",
                                     RooArgSet( *rv_mu_qcd_sig_ldp_2b, fv_sf_qcd_sig_2b, fv_Rlsb_passfail ) ) ;

         rfv_mu_qcd_sig_3b = new RooFormulaVar( "mu_qcd_sig_3b",
                                     "mu_qcd_sig_ldp_3b * sf_qcd_sig_3b * Rlsb_passfail",
                                     RooArgSet( *rv_mu_qcd_sig_ldp_3b, fv_sf_qcd_sig_3b, fv_Rlsb_passfail ) ) ;


         rfv_mu_qcd_sb_1b = new RooFormulaVar( "mu_qcd_sb_1b",
                                     "mu_qcd_sb_ldp_1b * sf_qcd_sb_1b * Rlsb_passfail",
                                     RooArgSet( *rv_mu_qcd_sb_ldp_1b, fv_sf_qcd_sb_1b, fv_Rlsb_passfail ) ) ;

         rfv_mu_qcd_sb_2b = new RooFormulaVar( "mu_qcd_sb_2b",
                                     "mu_qcd_sb_ldp_2b * sf_qcd_sb_2b * Rlsb_passfail",
                                     RooArgSet( *rv_mu_qcd_sb_ldp_2b, fv_sf_qcd_sb_2b, fv_Rlsb_passfail ) ) ;

         rfv_mu_qcd_sb_3b = new RooFormulaVar( "mu_qcd_sb_3b",
                                     "mu_qcd_sb_ldp_3b * sf_qcd_sb_3b * Rlsb_passfail",
                                     RooArgSet( *rv_mu_qcd_sb_ldp_3b, fv_sf_qcd_sb_3b, fv_Rlsb_passfail ) ) ;



         rv_mu_qcd_sig_1b = rfv_mu_qcd_sig_1b ;
         rv_mu_qcd_sig_2b = rfv_mu_qcd_sig_2b ;
         rv_mu_qcd_sig_3b = rfv_mu_qcd_sig_3b ;

         rv_mu_qcd_sb_1b = rfv_mu_qcd_sb_1b ;
         rv_mu_qcd_sb_2b = rfv_mu_qcd_sb_2b ;
         rv_mu_qcd_sb_3b = rfv_mu_qcd_sb_3b ;


      } else {

         //--- Here, float the SIG or SB variable and derive the LSB variable with the LSB ratio.
         //
         //    This option is not recommended unless you are making profile plots of the SIG or SB variables.
         //

         rfv_mu_qcd_sig_ldp_1b = new RooFormulaVar( "mu_qcd_sig_ldp_1b",
                                     "mu_qcd_sig_1b * (1.0/sf_qcd_sig_1b) * ( 1.0 / Rlsb_passfail )",
                                     RooArgSet( *rv_mu_qcd_sig_1b, fv_sf_qcd_sig_1b, fv_Rlsb_passfail ) ) ;

         rfv_mu_qcd_sig_ldp_2b = new RooFormulaVar( "mu_qcd_sig_ldp_2b",
                                     "mu_qcd_sig_2b * (1.0/sf_qcd_sig_2b) * ( 1.0 / Rlsb_passfail )",
                                     RooArgSet( *rv_mu_qcd_sig_2b, fv_sf_qcd_sig_2b, fv_Rlsb_passfail ) ) ;

         rfv_mu_qcd_sig_ldp_3b = new RooFormulaVar( "mu_qcd_sig_ldp_3b",
                                     "mu_qcd_sig_3b * (1.0/sf_qcd_sig_3b) * ( 1.0 / Rlsb_passfail )",
                                     RooArgSet( *rv_mu_qcd_sig_3b, fv_sf_qcd_sig_3b, fv_Rlsb_passfail ) ) ;



         rfv_mu_qcd_sb_ldp_1b = new RooFormulaVar( "mu_qcd_sb_ldp_1b",
                                     "mu_qcd_sb_1b * (1.0/sf_qcd_sb_1b) * ( 1.0 / Rlsb_passfail )",
                                     RooArgSet( *rv_mu_qcd_sb_1b, fv_sf_qcd_sb_1b, fv_Rlsb_passfail ) ) ;

         rfv_mu_qcd_sb_ldp_2b = new RooFormulaVar( "mu_qcd_sb_ldp_2b",
                                     "mu_qcd_sb_2b * (1.0/sf_qcd_sb_2b) * ( 1.0 / Rlsb_passfail )",
                                     RooArgSet( *rv_mu_qcd_sb_2b, fv_sf_qcd_sb_2b, fv_Rlsb_passfail ) ) ;

         rfv_mu_qcd_sb_ldp_3b = new RooFormulaVar( "mu_qcd_sb_ldp_3b",
                                     "mu_qcd_sb_3b * (1.0/sf_qcd_sb_3b) * ( 1.0 / Rlsb_passfail )",
                                     RooArgSet( *rv_mu_qcd_sb_3b, fv_sf_qcd_sb_3b, fv_Rlsb_passfail ) ) ;


         rv_mu_qcd_sig_ldp_1b = rfv_mu_qcd_sig_ldp_1b ;
         rv_mu_qcd_sig_ldp_2b = rfv_mu_qcd_sig_ldp_2b ;
         rv_mu_qcd_sig_ldp_3b = rfv_mu_qcd_sig_ldp_3b ;

         rv_mu_qcd_sb_ldp_1b = rfv_mu_qcd_sb_ldp_1b ;
         rv_mu_qcd_sb_ldp_2b = rfv_mu_qcd_sb_ldp_2b ;
         rv_mu_qcd_sb_ldp_3b = rfv_mu_qcd_sb_ldp_3b ;

      }

   //-------------------------------



    //-- SUSY : float mu_susy_sig_1b and derive other 17 parameters from it and MC ratios.

      rv_mu_susy_sb_1b = new RooFormulaVar( "mu_susy_sb_1b",
                                        "mu_susymc_sb_1b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sb_1b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sig_sl_1b = new RooFormulaVar( "mu_susy_sig_sl_1b",
                                        "mu_susymc_sig_sl_1b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sig_sl_1b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sb_sl_1b = new RooFormulaVar( "mu_susy_sb_sl_1b",
                                        "mu_susymc_sb_sl_1b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sb_sl_1b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sig_ldp_1b = new RooFormulaVar( "mu_susy_sig_ldp_1b",
                                        "mu_susymc_sig_ldp_1b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sig_ldp_1b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sb_ldp_1b = new RooFormulaVar( "mu_susy_sb_ldp_1b",
                                        "mu_susymc_sb_ldp_1b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sb_ldp_1b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;



      rv_mu_susy_sig_2b = new RooFormulaVar( "mu_susy_sig_2b",
                                        "mu_susymc_sig_2b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sig_2b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sb_2b = new RooFormulaVar( "mu_susy_sb_2b",
                                        "mu_susymc_sb_2b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sb_2b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sig_sl_2b = new RooFormulaVar( "mu_susy_sig_sl_2b",
                                        "mu_susymc_sig_sl_2b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sig_sl_2b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sb_sl_2b = new RooFormulaVar( "mu_susy_sb_sl_2b",
                                        "mu_susymc_sb_sl_2b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sb_sl_2b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sig_ldp_2b = new RooFormulaVar( "mu_susy_sig_ldp_2b",
                                        "mu_susymc_sig_ldp_2b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sig_ldp_2b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sb_ldp_2b = new RooFormulaVar( "mu_susy_sb_ldp_2b",
                                        "mu_susymc_sb_ldp_2b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sb_ldp_2b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;




      rv_mu_susy_sig_3b = new RooFormulaVar( "mu_susy_sig_3b",
                                        "mu_susymc_sig_3b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sig_3b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sb_3b = new RooFormulaVar( "mu_susy_sb_3b",
                                        "mu_susymc_sb_3b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sb_3b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sig_sl_3b = new RooFormulaVar( "mu_susy_sig_sl_3b",
                                        "mu_susymc_sig_sl_3b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sig_sl_3b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sb_sl_3b = new RooFormulaVar( "mu_susy_sb_sl_3b",
                                        "mu_susymc_sb_sl_3b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sb_sl_3b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sig_ldp_3b = new RooFormulaVar( "mu_susy_sig_ldp_3b",
                                        "mu_susymc_sig_ldp_3b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sig_ldp_3b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;

      rv_mu_susy_sb_ldp_3b = new RooFormulaVar( "mu_susy_sb_ldp_3b",
                                        "mu_susymc_sb_ldp_3b * (mu_susy_sig_1b/mu_susymc_sig_1b)",
                                        RooArgSet( *rv_mu_susymc_sb_ldp_3b, *rv_mu_susy_sig_1b, *rv_mu_susymc_sig_1b ) ) ;





    //-- Z to nu nu

      rv_mu_znn_sig_ldp_1b = new RooFormulaVar( "mu_znn_sig_ldp_1b",
                              "mu_Znnmc_sig_ldp_1b",
                              RooArgSet( *rv_mu_Znnmc_sig_ldp_1b ) ) ;

      rv_mu_znn_sig_ldp_2b = new RooFormulaVar( "mu_znn_sig_ldp_2b",
                              "mu_Znnmc_sig_ldp_2b",
                              RooArgSet( *rv_mu_Znnmc_sig_ldp_2b ) ) ;

      rv_mu_znn_sig_ldp_3b = new RooFormulaVar( "mu_znn_sig_ldp_3b",
                              "mu_Znnmc_sig_ldp_3b",
                              RooArgSet( *rv_mu_Znnmc_sig_ldp_3b ) ) ;




      rv_mu_znn_sb_ldp_1b = new RooFormulaVar( "mu_znn_sb_ldp_1b",
                              "mu_Znnmc_sb_ldp_1b",
                              RooArgSet( *rv_mu_Znnmc_sb_ldp_1b ) ) ;

      rv_mu_znn_sb_ldp_2b = new RooFormulaVar( "mu_znn_sb_ldp_2b",
                              "mu_Znnmc_sb_ldp_2b",
                              RooArgSet( *rv_mu_Znnmc_sb_ldp_2b ) ) ;

      rv_mu_znn_sb_ldp_3b = new RooFormulaVar( "mu_znn_sb_ldp_3b",
                              "mu_Znnmc_sb_ldp_3b",
                              RooArgSet( *rv_mu_Znnmc_sb_ldp_3b ) ) ;



     //-- Float the Znn 1b vars and derive 2b and 3b vars using knn ratios.



      rv_mu_zee_sig = new RooFormulaVar( "mu_zee_sig",
                                    "( mu_znn_sig_1b / knn_sig_1b ) * sf_ee * ( (acc_ee_sig * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                      RooArgSet( *rv_mu_znn_sig_1b, fv_knn_sig_1b, fv_sf_ee, fv_acc_ee_sig, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

      rv_mu_zee_sb = new RooFormulaVar( "mu_zee_sb",
                                    "( mu_znn_sb_1b / knn_sb_1b ) * sf_ee * ( (acc_ee_sb * eff_ee ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                      RooArgSet( *rv_mu_znn_sb_1b, fv_knn_sb_1b, fv_sf_ee, fv_acc_ee_sb, fv_eff_ee, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

      rv_mu_zmm_sig = new RooFormulaVar( "mu_zmm_sig",
                                    "( mu_znn_sig_1b / knn_sig_1b ) * sf_mm * ( (acc_mm_sig * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                      RooArgSet( *rv_mu_znn_sig_1b, fv_knn_sig_1b, fv_sf_mm, fv_acc_mm_sig, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;

      rv_mu_zmm_sb = new RooFormulaVar( "mu_zmm_sb",
                                    "( mu_znn_sb_1b / knn_sb_1b ) * sf_mm * ( (acc_mm_sb * eff_mm ) / ( znnoverll_bfratio * dataoverll_lumiratio ) )",
                                      RooArgSet( *rv_mu_znn_sb_1b, fv_knn_sb_1b, fv_sf_mm, fv_acc_mm_sb, fv_eff_mm, *rv_znnoverll_bfratio, *rv_dataoverll_lumiratio ) ) ;







      rv_mu_znn_sig_2b = new RooFormulaVar( "mu_znn_sig_2b",
                                               "mu_znn_sig_1b * ( knn_sig_2b / knn_sig_1b )",
                                         RooArgSet( *rv_mu_znn_sig_1b, fv_knn_sig_2b, fv_knn_sig_1b ) ) ;

      rv_mu_znn_sig_3b = new RooFormulaVar( "mu_znn_sig_3b",
                                               "mu_znn_sig_1b * ( knn_sig_3b / knn_sig_1b )",
                                         RooArgSet( *rv_mu_znn_sig_1b, fv_knn_sig_3b, fv_knn_sig_1b ) ) ;



      rv_mu_znn_sb_2b = new RooFormulaVar( "mu_znn_sb_2b",
                                               "mu_znn_sb_1b * ( knn_sb_2b / knn_sb_1b )",
                                         RooArgSet( *rv_mu_znn_sb_1b, fv_knn_sb_2b, fv_knn_sb_1b ) ) ;

      rv_mu_znn_sb_3b = new RooFormulaVar( "mu_znn_sb_3b",
                                               "mu_znn_sb_1b * ( knn_sb_3b / knn_sb_1b )",
                                         RooArgSet( *rv_mu_znn_sb_1b, fv_knn_sb_3b, fv_knn_sb_1b ) ) ;








    //-- Parametric relations between correlated efficiency scale factors.



      rv_eff_sf_sig_1b = new RooFormulaVar( "eff_sf_sig_1b",
                                         "mean_eff_sf_sig_1b * pow( exp( width_eff_sf_sig_1b/mean_eff_sf_sig_1b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_1b, *rv_width_eff_sf_sig_1b, *rv_mean_eff_sf_sig_1b, eff_sf_prim ) ) ;

      rv_eff_sf_sb_1b = new RooFormulaVar( "eff_sf_sb_1b",
                                         "mean_eff_sf_sb_1b * pow( exp( width_eff_sf_sb_1b/mean_eff_sf_sb_1b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_1b, *rv_width_eff_sf_sb_1b, *rv_mean_eff_sf_sb_1b, eff_sf_prim ) ) ;

      rv_eff_sf_sig_sl_1b = new RooFormulaVar( "eff_sf_sig_sl_1b",
                                         "mean_eff_sf_sig_sl_1b * pow( exp( width_eff_sf_sig_sl_1b/mean_eff_sf_sig_sl_1b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_sl_1b, *rv_width_eff_sf_sig_sl_1b, *rv_mean_eff_sf_sig_sl_1b, eff_sf_prim ) ) ;

      rv_eff_sf_sb_sl_1b = new RooFormulaVar( "eff_sf_sb_sl_1b",
                                         "mean_eff_sf_sb_sl_1b * pow( exp( width_eff_sf_sb_sl_1b/mean_eff_sf_sb_sl_1b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_sl_1b, *rv_width_eff_sf_sb_sl_1b, *rv_mean_eff_sf_sb_sl_1b, eff_sf_prim ) ) ;

      rv_eff_sf_sig_ldp_1b = new RooFormulaVar( "eff_sf_sig_ldp_1b",
                                         "mean_eff_sf_sig_ldp_1b * pow( exp( width_eff_sf_sig_ldp_1b/mean_eff_sf_sig_ldp_1b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_ldp_1b, *rv_width_eff_sf_sig_ldp_1b, *rv_mean_eff_sf_sig_ldp_1b, eff_sf_prim ) ) ;

      rv_eff_sf_sb_ldp_1b = new RooFormulaVar( "eff_sf_sb_ldp_1b",
                                         "mean_eff_sf_sb_ldp_1b * pow( exp( width_eff_sf_sb_ldp_1b/mean_eff_sf_sb_ldp_1b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_ldp_1b, *rv_width_eff_sf_sb_ldp_1b, *rv_mean_eff_sf_sb_ldp_1b, eff_sf_prim ) ) ;






      rv_eff_sf_sig_2b = new RooFormulaVar( "eff_sf_sig_2b",
                                         "mean_eff_sf_sig_2b * pow( exp( width_eff_sf_sig_2b/mean_eff_sf_sig_2b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_2b, *rv_width_eff_sf_sig_2b, *rv_mean_eff_sf_sig_2b, eff_sf_prim ) ) ;

      rv_eff_sf_sb_2b = new RooFormulaVar( "eff_sf_sb_2b",
                                         "mean_eff_sf_sb_2b * pow( exp( width_eff_sf_sb_2b/mean_eff_sf_sb_2b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_2b, *rv_width_eff_sf_sb_2b, *rv_mean_eff_sf_sb_2b, eff_sf_prim ) ) ;

      rv_eff_sf_sig_sl_2b = new RooFormulaVar( "eff_sf_sig_sl_2b",
                                         "mean_eff_sf_sig_sl_2b * pow( exp( width_eff_sf_sig_sl_2b/mean_eff_sf_sig_sl_2b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_sl_2b, *rv_width_eff_sf_sig_sl_2b, *rv_mean_eff_sf_sig_sl_2b, eff_sf_prim ) ) ;

      rv_eff_sf_sb_sl_2b = new RooFormulaVar( "eff_sf_sb_sl_2b",
                                         "mean_eff_sf_sb_sl_2b * pow( exp( width_eff_sf_sb_sl_2b/mean_eff_sf_sb_sl_2b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_sl_2b, *rv_width_eff_sf_sb_sl_2b, *rv_mean_eff_sf_sb_sl_2b, eff_sf_prim ) ) ;

      rv_eff_sf_sig_ldp_2b = new RooFormulaVar( "eff_sf_sig_ldp_2b",
                                         "mean_eff_sf_sig_ldp_2b * pow( exp( width_eff_sf_sig_ldp_2b/mean_eff_sf_sig_ldp_2b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_ldp_2b, *rv_width_eff_sf_sig_ldp_2b, *rv_mean_eff_sf_sig_ldp_2b, eff_sf_prim ) ) ;

      rv_eff_sf_sb_ldp_2b = new RooFormulaVar( "eff_sf_sb_ldp_2b",
                                         "mean_eff_sf_sb_ldp_2b * pow( exp( width_eff_sf_sb_ldp_2b/mean_eff_sf_sb_ldp_2b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_ldp_2b, *rv_width_eff_sf_sb_ldp_2b, *rv_mean_eff_sf_sb_ldp_2b, eff_sf_prim ) ) ;






      rv_eff_sf_sig_3b = new RooFormulaVar( "eff_sf_sig_3b",
                                         "mean_eff_sf_sig_3b * pow( exp( width_eff_sf_sig_3b/mean_eff_sf_sig_3b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_3b, *rv_width_eff_sf_sig_3b, *rv_mean_eff_sf_sig_3b, eff_sf_prim ) ) ;

      rv_eff_sf_sb_3b = new RooFormulaVar( "eff_sf_sb_3b",
                                         "mean_eff_sf_sb_3b * pow( exp( width_eff_sf_sb_3b/mean_eff_sf_sb_3b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_3b, *rv_width_eff_sf_sb_3b, *rv_mean_eff_sf_sb_3b, eff_sf_prim ) ) ;

      rv_eff_sf_sig_sl_3b = new RooFormulaVar( "eff_sf_sig_sl_3b",
                                         "mean_eff_sf_sig_sl_3b * pow( exp( width_eff_sf_sig_sl_3b/mean_eff_sf_sig_sl_3b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_sl_3b, *rv_width_eff_sf_sig_sl_3b, *rv_mean_eff_sf_sig_sl_3b, eff_sf_prim ) ) ;

      rv_eff_sf_sb_sl_3b = new RooFormulaVar( "eff_sf_sb_sl_3b",
                                         "mean_eff_sf_sb_sl_3b * pow( exp( width_eff_sf_sb_sl_3b/mean_eff_sf_sb_sl_3b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_sl_3b, *rv_width_eff_sf_sb_sl_3b, *rv_mean_eff_sf_sb_sl_3b, eff_sf_prim ) ) ;

      rv_eff_sf_sig_ldp_3b = new RooFormulaVar( "eff_sf_sig_ldp_3b",
                                         "mean_eff_sf_sig_ldp_3b * pow( exp( width_eff_sf_sig_ldp_3b/mean_eff_sf_sig_ldp_3b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sig_ldp_3b, *rv_width_eff_sf_sig_ldp_3b, *rv_mean_eff_sf_sig_ldp_3b, eff_sf_prim ) ) ;

      rv_eff_sf_sb_ldp_3b = new RooFormulaVar( "eff_sf_sb_ldp_3b",
                                         "mean_eff_sf_sb_ldp_3b * pow( exp( width_eff_sf_sb_ldp_3b/mean_eff_sf_sb_ldp_3b ), eff_sf_prim )",
                                         RooArgSet( *rv_mean_eff_sf_sb_ldp_3b, *rv_width_eff_sf_sb_ldp_3b, *rv_mean_eff_sf_sb_ldp_3b, eff_sf_prim ) ) ;





     //--- btag efficiency scale factors.

      if ( rv_deff_dbtageff_sig_1b->getVal() > 0. ) {
         rv_btageff_sf_sig_1b = new RooFormulaVar( "btageff_sf_sig_1b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_1b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_1b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sig_1b->setVal( -1.0*(rv_deff_dbtageff_sig_1b->getVal()) ) ;
         rv_btageff_sf_sig_1b = new RooFormulaVar( "btageff_sf_sig_1b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_1b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_1b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sb_1b->getVal() > 0. ) {
         rv_btageff_sf_sb_1b = new RooFormulaVar( "btageff_sf_sb_1b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_1b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_1b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sb_1b->setVal( -1.0*(rv_deff_dbtageff_sb_1b->getVal()) ) ;
         rv_btageff_sf_sb_1b = new RooFormulaVar( "btageff_sf_sb_1b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_1b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_1b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sig_sl_1b->getVal() > 0. ) {
         rv_btageff_sf_sig_sl_1b = new RooFormulaVar( "btageff_sf_sig_sl_1b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_sl_1b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_sl_1b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sig_sl_1b->setVal( -1.0*(rv_deff_dbtageff_sig_sl_1b->getVal()) ) ;
         rv_btageff_sf_sig_sl_1b = new RooFormulaVar( "btageff_sf_sig_sl_1b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_sl_1b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_sl_1b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sb_sl_1b->getVal() > 0. ) {
         rv_btageff_sf_sb_sl_1b = new RooFormulaVar( "btageff_sf_sb_sl_1b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_sl_1b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_sl_1b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sb_sl_1b->setVal( -1.0*(rv_deff_dbtageff_sb_sl_1b->getVal()) ) ;
         rv_btageff_sf_sb_sl_1b = new RooFormulaVar( "btageff_sf_sb_sl_1b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_sl_1b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_sl_1b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sig_ldp_1b->getVal() > 0. ) {
         rv_btageff_sf_sig_ldp_1b = new RooFormulaVar( "btageff_sf_sig_ldp_1b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_ldp_1b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_ldp_1b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sig_ldp_1b->setVal( -1.0*(rv_deff_dbtageff_sig_ldp_1b->getVal()) ) ;
         rv_btageff_sf_sig_ldp_1b = new RooFormulaVar( "btageff_sf_sig_ldp_1b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_ldp_1b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_ldp_1b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sb_ldp_1b->getVal() > 0. ) {
         rv_btageff_sf_sb_ldp_1b = new RooFormulaVar( "btageff_sf_sb_ldp_1b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_ldp_1b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_ldp_1b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sb_ldp_1b->setVal( -1.0*(rv_deff_dbtageff_sb_ldp_1b->getVal()) ) ;
         rv_btageff_sf_sb_ldp_1b = new RooFormulaVar( "btageff_sf_sb_ldp_1b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_ldp_1b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_ldp_1b, btageff_sf_prim ) ) ;
      }









      if ( rv_deff_dbtageff_sig_2b->getVal() > 0. ) {
         rv_btageff_sf_sig_2b = new RooFormulaVar( "btageff_sf_sig_2b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_2b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_2b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sig_2b->setVal( -1.0*(rv_deff_dbtageff_sig_2b->getVal()) ) ;
         rv_btageff_sf_sig_2b = new RooFormulaVar( "btageff_sf_sig_2b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_2b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_2b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sb_2b->getVal() > 0. ) {
         rv_btageff_sf_sb_2b = new RooFormulaVar( "btageff_sf_sb_2b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_2b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_2b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sb_2b->setVal( -1.0*(rv_deff_dbtageff_sb_2b->getVal()) ) ;
         rv_btageff_sf_sb_2b = new RooFormulaVar( "btageff_sf_sb_2b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_2b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_2b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sig_sl_2b->getVal() > 0. ) {
         rv_btageff_sf_sig_sl_2b = new RooFormulaVar( "btageff_sf_sig_sl_2b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_sl_2b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_sl_2b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sig_sl_2b->setVal( -1.0*(rv_deff_dbtageff_sig_sl_2b->getVal()) ) ;
         rv_btageff_sf_sig_sl_2b = new RooFormulaVar( "btageff_sf_sig_sl_2b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_sl_2b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_sl_2b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sb_sl_2b->getVal() > 0. ) {
         rv_btageff_sf_sb_sl_2b = new RooFormulaVar( "btageff_sf_sb_sl_2b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_sl_2b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_sl_2b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sb_sl_2b->setVal( -1.0*(rv_deff_dbtageff_sb_sl_2b->getVal()) ) ;
         rv_btageff_sf_sb_sl_2b = new RooFormulaVar( "btageff_sf_sb_sl_2b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_sl_2b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_sl_2b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sig_ldp_2b->getVal() > 0. ) {
         rv_btageff_sf_sig_ldp_2b = new RooFormulaVar( "btageff_sf_sig_ldp_2b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_ldp_2b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_ldp_2b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sig_ldp_2b->setVal( -1.0*(rv_deff_dbtageff_sig_ldp_2b->getVal()) ) ;
         rv_btageff_sf_sig_ldp_2b = new RooFormulaVar( "btageff_sf_sig_ldp_2b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_ldp_2b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_ldp_2b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sb_ldp_2b->getVal() > 0. ) {
         rv_btageff_sf_sb_ldp_2b = new RooFormulaVar( "btageff_sf_sb_ldp_2b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_ldp_2b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_ldp_2b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sb_ldp_2b->setVal( -1.0*(rv_deff_dbtageff_sb_ldp_2b->getVal()) ) ;
         rv_btageff_sf_sb_ldp_2b = new RooFormulaVar( "btageff_sf_sb_ldp_2b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_ldp_2b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_ldp_2b, btageff_sf_prim ) ) ;
      }









      if ( rv_deff_dbtageff_sig_3b->getVal() > 0. ) {
         rv_btageff_sf_sig_3b = new RooFormulaVar( "btageff_sf_sig_3b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_3b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_3b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sig_3b->setVal( -1.0*(rv_deff_dbtageff_sig_3b->getVal()) ) ;
         rv_btageff_sf_sig_3b = new RooFormulaVar( "btageff_sf_sig_3b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_3b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_3b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sb_3b->getVal() > 0. ) {
         rv_btageff_sf_sb_3b = new RooFormulaVar( "btageff_sf_sb_3b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_3b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_3b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sb_3b->setVal( -1.0*(rv_deff_dbtageff_sb_3b->getVal()) ) ;
         rv_btageff_sf_sb_3b = new RooFormulaVar( "btageff_sf_sb_3b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_3b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_3b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sig_sl_3b->getVal() > 0. ) {
         rv_btageff_sf_sig_sl_3b = new RooFormulaVar( "btageff_sf_sig_sl_3b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_sl_3b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_sl_3b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sig_sl_3b->setVal( -1.0*(rv_deff_dbtageff_sig_sl_3b->getVal()) ) ;
         rv_btageff_sf_sig_sl_3b = new RooFormulaVar( "btageff_sf_sig_sl_3b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_sl_3b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_sl_3b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sb_sl_3b->getVal() > 0. ) {
         rv_btageff_sf_sb_sl_3b = new RooFormulaVar( "btageff_sf_sb_sl_3b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_sl_3b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_sl_3b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sb_sl_3b->setVal( -1.0*(rv_deff_dbtageff_sb_sl_3b->getVal()) ) ;
         rv_btageff_sf_sb_sl_3b = new RooFormulaVar( "btageff_sf_sb_sl_3b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_sl_3b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_sl_3b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sig_ldp_3b->getVal() > 0. ) {
         rv_btageff_sf_sig_ldp_3b = new RooFormulaVar( "btageff_sf_sig_ldp_3b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_ldp_3b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_ldp_3b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sig_ldp_3b->setVal( -1.0*(rv_deff_dbtageff_sig_ldp_3b->getVal()) ) ;
         rv_btageff_sf_sig_ldp_3b = new RooFormulaVar( "btageff_sf_sig_ldp_3b",
                                            "pow( exp( btageff_err * deff_dbtageff_sig_ldp_3b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_ldp_3b, btageff_sf_prim ) ) ;
      }

      if ( rv_deff_dbtageff_sb_ldp_3b->getVal() > 0. ) {
         rv_btageff_sf_sb_ldp_3b = new RooFormulaVar( "btageff_sf_sb_ldp_3b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_ldp_3b ), btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_ldp_3b, btageff_sf_prim ) ) ;
      } else {
         rv_deff_dbtageff_sb_ldp_3b->setVal( -1.0*(rv_deff_dbtageff_sb_ldp_3b->getVal()) ) ;
         rv_btageff_sf_sb_ldp_3b = new RooFormulaVar( "btageff_sf_sb_ldp_3b",
                                            "pow( exp( btageff_err * deff_dbtageff_sb_ldp_3b ), -1.0*btageff_sf_prim )",
                                            RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_ldp_3b, btageff_sf_prim ) ) ;
      }











///   rv_btageff_sf_sb_1b = new RooFormulaVar( "btageff_sf_sb_1b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sb_1b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_1b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sig_sl_1b = new RooFormulaVar( "btageff_sf_sig_sl_1b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sig_sl_1b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_sl_1b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sb_sl_1b = new RooFormulaVar( "btageff_sf_sb_sl_1b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sb_sl_1b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_sl_1b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sig_ldp_1b = new RooFormulaVar( "btageff_sf_sig_ldp_1b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sig_ldp_1b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_ldp_1b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sb_ldp_1b = new RooFormulaVar( "btageff_sf_sb_ldp_1b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sb_ldp_1b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_ldp_1b, btageff_sf_prim ) ) ;



///   rv_btageff_sf_sig_2b = new RooFormulaVar( "btageff_sf_sig_2b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sig_2b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_2b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sb_2b = new RooFormulaVar( "btageff_sf_sb_2b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sb_2b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_2b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sig_sl_2b = new RooFormulaVar( "btageff_sf_sig_sl_2b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sig_sl_2b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_sl_2b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sb_sl_2b = new RooFormulaVar( "btageff_sf_sb_sl_2b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sb_sl_2b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_sl_2b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sig_ldp_2b = new RooFormulaVar( "btageff_sf_sig_ldp_2b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sig_ldp_2b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_ldp_2b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sb_ldp_2b = new RooFormulaVar( "btageff_sf_sb_ldp_2b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sb_ldp_2b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_ldp_2b, btageff_sf_prim ) ) ;



///   rv_btageff_sf_sig_3b = new RooFormulaVar( "btageff_sf_sig_3b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sig_3b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_3b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sb_3b = new RooFormulaVar( "btageff_sf_sb_3b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sb_3b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_3b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sig_sl_3b = new RooFormulaVar( "btageff_sf_sig_sl_3b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sig_sl_3b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_sl_3b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sb_sl_3b = new RooFormulaVar( "btageff_sf_sb_sl_3b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sb_sl_3b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_sl_3b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sig_ldp_3b = new RooFormulaVar( "btageff_sf_sig_ldp_3b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sig_ldp_3b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sig_ldp_3b, btageff_sf_prim ) ) ;

///   rv_btageff_sf_sb_ldp_3b = new RooFormulaVar( "btageff_sf_sb_ldp_3b",
///                                      "pow( exp( btageff_err * deff_dbtageff_sb_ldp_3b ), btageff_sf_prim )",
///                                      RooArgSet( rv_btageff_err, *rv_deff_dbtageff_sb_ldp_3b, btageff_sf_prim ) ) ;







    //+++++++++++++ Expected counts for observables in terms of parameters ++++++++++++++++++

       printf(" --- Defining expected counts in terms of parameters.\n" ) ;


      rv_n_sig_1b         = new RooFormulaVar( "n_sig_1b",
                                     "mu_ttwj_sig_1b + mu_qcd_sig_1b + mu_znn_sig_1b + btageff_sf_sig_1b*eff_sf_sig_1b*( mu_susy_sig_1b)",
                                     RooArgSet( *rv_mu_ttwj_sig_1b, *rv_mu_qcd_sig_1b, *rv_mu_znn_sig_1b, *rv_btageff_sf_sig_1b, *rv_eff_sf_sig_1b, *rv_mu_susy_sig_1b ) ) ;

      rv_n_sb_1b          = new RooFormulaVar( "n_sb_1b",
                                     "mu_ttwj_sb_1b  + mu_qcd_sb_1b  + mu_znn_sb_1b  + btageff_sf_sb_1b*eff_sf_sb_1b*( mu_susy_sb_1b )",
                                     RooArgSet( *rv_mu_ttwj_sb_1b , *rv_mu_qcd_sb_1b , *rv_mu_znn_sb_1b , *rv_btageff_sf_sb_1b, *rv_eff_sf_sb_1b, *rv_mu_susy_sb_1b  ) ) ;

      rv_n_sig_ldp_1b     = new RooFormulaVar( "n_sig_ldp_1b",
                                     "mu_qcd_sig_ldp_1b + btageff_sf_sig_ldp_1b*eff_sf_sig_ldp_1b*( sf_mc * (mu_ttwj_sig_ldp_1b + mu_znn_sig_ldp_1b ) + mu_susy_sig_ldp_1b)",
                                     RooArgSet( *rv_mu_qcd_sig_ldp_1b, *rv_btageff_sf_sig_ldp_1b, *rv_eff_sf_sig_ldp_1b, fv_sf_mc, *rv_mu_ttwj_sig_ldp_1b, *rv_mu_znn_sig_ldp_1b, *rv_mu_susy_sig_ldp_1b ) ) ;

      rv_n_sb_ldp_1b      = new RooFormulaVar( "n_sb_ldp_1b",
                                     "mu_qcd_sb_ldp_1b + btageff_sf_sb_ldp_1b*eff_sf_sb_ldp_1b*( sf_mc * (mu_ttwj_sb_ldp_1b + mu_znn_sb_ldp_1b ) + mu_susy_sb_ldp_1b)",
                                     RooArgSet( *rv_mu_qcd_sb_ldp_1b, *rv_btageff_sf_sb_ldp_1b, *rv_eff_sf_sb_ldp_1b, fv_sf_mc, *rv_mu_ttwj_sb_ldp_1b, *rv_mu_znn_sb_ldp_1b, *rv_mu_susy_sb_ldp_1b ) ) ;

      rv_n_sig_sl_1b      = new RooFormulaVar( "n_sig_sl_1b",
                                     "mu_ttwj_sig_sl_1b + btageff_sf_sig_sl_1b*eff_sf_sig_sl_1b*mu_susy_sig_sl_1b",
                                     RooArgSet( *rv_mu_ttwj_sig_sl_1b, *rv_btageff_sf_sig_sl_1b, *rv_eff_sf_sig_sl_1b, *rv_mu_susy_sig_sl_1b ) ) ;

      rv_n_sb_sl_1b       = new RooFormulaVar( "n_sb_sl_1b",
                                     "mu_ttwj_sb_sl_1b + btageff_sf_sb_sl_1b*eff_sf_sb_sl_1b*mu_susy_sb_sl_1b",
                                     RooArgSet( *rv_mu_ttwj_sb_sl_1b, *rv_btageff_sf_sb_sl_1b, *rv_eff_sf_sb_sl_1b, *rv_mu_susy_sb_sl_1b ) ) ;




      rv_n_sig_2b         = new RooFormulaVar( "n_sig_2b",
                                     "mu_ttwj_sig_2b + mu_qcd_sig_2b + mu_znn_sig_2b + btageff_sf_sig_2b*eff_sf_sig_2b*( mu_susy_sig_2b)",
                                     RooArgSet( *rv_mu_ttwj_sig_2b, *rv_mu_qcd_sig_2b, *rv_mu_znn_sig_2b, *rv_btageff_sf_sig_2b, *rv_eff_sf_sig_2b, *rv_mu_susy_sig_2b ) ) ;

      rv_n_sb_2b          = new RooFormulaVar( "n_sb_2b",
                                     "mu_ttwj_sb_2b  + mu_qcd_sb_2b  + mu_znn_sb_2b  + btageff_sf_sb_2b*eff_sf_sb_2b*( mu_susy_sb_2b )",
                                     RooArgSet( *rv_mu_ttwj_sb_2b , *rv_mu_qcd_sb_2b , *rv_mu_znn_sb_2b , *rv_btageff_sf_sb_2b, *rv_eff_sf_sb_2b, *rv_mu_susy_sb_2b  ) ) ;

      rv_n_sig_ldp_2b     = new RooFormulaVar( "n_sig_ldp_2b",
                                     "mu_qcd_sig_ldp_2b + btageff_sf_sig_ldp_2b*eff_sf_sig_ldp_2b*( sf_mc * (mu_ttwj_sig_ldp_2b + mu_znn_sig_ldp_2b ) + mu_susy_sig_ldp_2b)",
                                     RooArgSet( *rv_mu_qcd_sig_ldp_2b, *rv_btageff_sf_sig_ldp_2b, *rv_eff_sf_sig_ldp_2b, fv_sf_mc, *rv_mu_ttwj_sig_ldp_2b, *rv_mu_znn_sig_ldp_2b, *rv_mu_susy_sig_ldp_2b ) ) ;

      rv_n_sb_ldp_2b      = new RooFormulaVar( "n_sb_ldp_2b",
                                     "mu_qcd_sb_ldp_2b + btageff_sf_sb_ldp_2b*eff_sf_sb_ldp_2b*( sf_mc * (mu_ttwj_sb_ldp_2b + mu_znn_sb_ldp_2b ) + mu_susy_sb_ldp_2b)",
                                     RooArgSet( *rv_mu_qcd_sb_ldp_2b, *rv_btageff_sf_sb_ldp_2b, *rv_eff_sf_sb_ldp_2b, fv_sf_mc, *rv_mu_ttwj_sb_ldp_2b, *rv_mu_znn_sb_ldp_2b, *rv_mu_susy_sb_ldp_2b ) ) ;

      rv_n_sig_sl_2b      = new RooFormulaVar( "n_sig_sl_2b",
                                     "mu_ttwj_sig_sl_2b + btageff_sf_sig_sl_2b*eff_sf_sig_sl_2b*mu_susy_sig_sl_2b",
                                     RooArgSet( *rv_mu_ttwj_sig_sl_2b, *rv_btageff_sf_sig_sl_2b, *rv_eff_sf_sig_sl_2b, *rv_mu_susy_sig_sl_2b ) ) ;

      rv_n_sb_sl_2b       = new RooFormulaVar( "n_sb_sl_2b",
                                     "mu_ttwj_sb_sl_2b + btageff_sf_sb_sl_2b*eff_sf_sb_sl_2b*mu_susy_sb_sl_2b",
                                     RooArgSet( *rv_mu_ttwj_sb_sl_2b, *rv_btageff_sf_sb_sl_2b, *rv_eff_sf_sb_sl_2b, *rv_mu_susy_sb_sl_2b ) ) ;




      rv_n_sig_3b         = new RooFormulaVar( "n_sig_3b",
                                     "mu_ttwj_sig_3b + mu_qcd_sig_3b + mu_znn_sig_3b + btageff_sf_sig_3b*eff_sf_sig_3b*( mu_susy_sig_3b)",
                                     RooArgSet( *rv_mu_ttwj_sig_3b, *rv_mu_qcd_sig_3b, *rv_mu_znn_sig_3b, *rv_btageff_sf_sig_3b, *rv_eff_sf_sig_3b, *rv_mu_susy_sig_3b ) ) ;

      rv_n_sb_3b          = new RooFormulaVar( "n_sb_3b",
                                     "mu_ttwj_sb_3b  + mu_qcd_sb_3b  + mu_znn_sb_3b  + btageff_sf_sb_3b*eff_sf_sb_3b*( mu_susy_sb_3b )",
                                     RooArgSet( *rv_mu_ttwj_sb_3b , *rv_mu_qcd_sb_3b , *rv_mu_znn_sb_3b , *rv_btageff_sf_sb_3b, *rv_eff_sf_sb_3b, *rv_mu_susy_sb_3b  ) ) ;

      rv_n_sig_ldp_3b     = new RooFormulaVar( "n_sig_ldp_3b",
                                     "mu_qcd_sig_ldp_3b + btageff_sf_sig_ldp_3b*eff_sf_sig_ldp_3b*( sf_mc * (mu_ttwj_sig_ldp_3b + mu_znn_sig_ldp_3b ) + mu_susy_sig_ldp_3b)",
                                     RooArgSet( *rv_mu_qcd_sig_ldp_3b, *rv_btageff_sf_sig_ldp_3b, *rv_eff_sf_sig_ldp_3b, fv_sf_mc, *rv_mu_ttwj_sig_ldp_3b, *rv_mu_znn_sig_ldp_3b, *rv_mu_susy_sig_ldp_3b ) ) ;

      rv_n_sb_ldp_3b      = new RooFormulaVar( "n_sb_ldp_3b",
                                     "mu_qcd_sb_ldp_3b + btageff_sf_sb_ldp_3b*eff_sf_sb_ldp_3b*( sf_mc * (mu_ttwj_sb_ldp_3b + mu_znn_sb_ldp_3b ) + mu_susy_sb_ldp_3b)",
                                     RooArgSet( *rv_mu_qcd_sb_ldp_3b, *rv_btageff_sf_sb_ldp_3b, *rv_eff_sf_sb_ldp_3b, fv_sf_mc, *rv_mu_ttwj_sb_ldp_3b, *rv_mu_znn_sb_ldp_3b, *rv_mu_susy_sb_ldp_3b ) ) ;

      rv_n_sig_sl_3b      = new RooFormulaVar( "n_sig_sl_3b",
                                     "mu_ttwj_sig_sl_3b + btageff_sf_sig_sl_3b*eff_sf_sig_sl_3b*mu_susy_sig_sl_3b",
                                     RooArgSet( *rv_mu_ttwj_sig_sl_3b, *rv_btageff_sf_sig_sl_3b, *rv_eff_sf_sig_sl_3b, *rv_mu_susy_sig_sl_3b ) ) ;

      rv_n_sb_sl_3b       = new RooFormulaVar( "n_sb_sl_3b",
                                     "mu_ttwj_sb_sl_3b + btageff_sf_sb_sl_3b*eff_sf_sb_sl_3b*mu_susy_sb_sl_3b",
                                     RooArgSet( *rv_mu_ttwj_sb_sl_3b, *rv_btageff_sf_sb_sl_3b, *rv_eff_sf_sb_sl_3b, *rv_mu_susy_sb_sl_3b ) ) ;





      rv_n_sig_ee      = new RooFormulaVar( "n_sig_ee",
                                     "mu_zee_sig / fsig_ee",
                                     RooArgSet( *rv_mu_zee_sig, fv_fsig_ee ) ) ;

      rv_n_sb_ee       = new RooFormulaVar( "n_sb_ee",
                                     "mu_zee_sb / fsig_ee",
                                     RooArgSet( *rv_mu_zee_sb, fv_fsig_ee ) ) ;

      rv_n_sig_mm      = new RooFormulaVar( "n_sig_mm",
                                     "mu_zmm_sig / fsig_mm",
                                     RooArgSet( *rv_mu_zmm_sig, fv_fsig_mm ) ) ;

      rv_n_sb_mm       = new RooFormulaVar( "n_sb_mm",
                                     "mu_zmm_sb / fsig_mm",
                                     RooArgSet( *rv_mu_zmm_sb, fv_fsig_mm ) ) ;


   //++++++++++++ PDFs for the likelihood +++++++++++++++++++++++++++++++++++++++++++++

      printf(" --- Defining PDFs of the likelihood.\n" ) ;

      pdf_Nsig_1b        = new RooPoisson( "pdf_Nsig_1b"        , "Nsig Poisson PDF, 1b"        , *rv_Nsig_1b        , *rv_n_sig_1b ) ;
      pdf_Nsb_1b         = new RooPoisson( "pdf_Nsb_1b"         , "Nsb Poisson PDF, 1b"         , *rv_Nsb_1b         , *rv_n_sb_1b ) ;
      pdf_Nsig_ldp_1b    = new RooPoisson( "pdf_Nsig_ldp_1b"    , "Nsig_ldp Poisson PDF, 1b"    , *rv_Nsig_ldp_1b    , *rv_n_sig_ldp_1b ) ;
      pdf_Nsb_ldp_1b     = new RooPoisson( "pdf_Nsb_ldp_1b"     , "Nsb_ldp Poisson PDF, 1b"     , *rv_Nsb_ldp_1b     , *rv_n_sb_ldp_1b ) ;
      pdf_Nsig_sl_1b     = new RooPoisson( "pdf_Nsig_sl_1b"     , "Nsig_sl Poisson PDF, 1b"     , *rv_Nsig_sl_1b     , *rv_n_sig_sl_1b ) ;
      pdf_Nsb_sl_1b      = new RooPoisson( "pdf_Nsb_sl_1b"      , "Nsb_sl Poisson PDF, 1b"      , *rv_Nsb_sl_1b      , *rv_n_sb_sl_1b ) ;

      pdf_Nsig_2b        = new RooPoisson( "pdf_Nsig_2b"        , "Nsig Poisson PDF, 2b"        , *rv_Nsig_2b        , *rv_n_sig_2b ) ;
      pdf_Nsb_2b         = new RooPoisson( "pdf_Nsb_2b"         , "Nsb Poisson PDF, 2b"         , *rv_Nsb_2b         , *rv_n_sb_2b ) ;
      pdf_Nsig_ldp_2b    = new RooPoisson( "pdf_Nsig_ldp_2b"    , "Nsig_ldp Poisson PDF, 2b"    , *rv_Nsig_ldp_2b    , *rv_n_sig_ldp_2b ) ;
      pdf_Nsb_ldp_2b     = new RooPoisson( "pdf_Nsb_ldp_2b"     , "Nsb_ldp Poisson PDF, 2b"     , *rv_Nsb_ldp_2b     , *rv_n_sb_ldp_2b ) ;
      pdf_Nsig_sl_2b     = new RooPoisson( "pdf_Nsig_sl_2b"     , "Nsig_sl Poisson PDF, 2b"     , *rv_Nsig_sl_2b     , *rv_n_sig_sl_2b ) ;
      pdf_Nsb_sl_2b      = new RooPoisson( "pdf_Nsb_sl_2b"      , "Nsb_sl Poisson PDF, 2b"      , *rv_Nsb_sl_2b      , *rv_n_sb_sl_2b ) ;

      pdf_Nsig_3b        = new RooPoisson( "pdf_Nsig_3b"        , "Nsig Poisson PDF, 3b"        , *rv_Nsig_3b        , *rv_n_sig_3b ) ;
      pdf_Nsb_3b         = new RooPoisson( "pdf_Nsb_3b"         , "Nsb Poisson PDF, 3b"         , *rv_Nsb_3b         , *rv_n_sb_3b ) ;
      pdf_Nsig_ldp_3b    = new RooPoisson( "pdf_Nsig_ldp_3b"    , "Nsig_ldp Poisson PDF, 3b"    , *rv_Nsig_ldp_3b    , *rv_n_sig_ldp_3b ) ;
      pdf_Nsb_ldp_3b     = new RooPoisson( "pdf_Nsb_ldp_3b"     , "Nsb_ldp Poisson PDF, 3b"     , *rv_Nsb_ldp_3b     , *rv_n_sb_ldp_3b ) ;
      pdf_Nsig_sl_3b     = new RooPoisson( "pdf_Nsig_sl_3b"     , "Nsig_sl Poisson PDF, 3b"     , *rv_Nsig_sl_3b     , *rv_n_sig_sl_3b ) ;
      pdf_Nsb_sl_3b      = new RooPoisson( "pdf_Nsb_sl_3b"      , "Nsb_sl Poisson PDF, 3b"      , *rv_Nsb_sl_3b      , *rv_n_sb_sl_3b ) ;


      pdf_Nsig_ee     = new RooPoisson( "pdf_Nsig_ee"     , "Nsig_ee Poisson PDF"     , *rv_Nsig_ee     , *rv_n_sig_ee ) ;
      pdf_Nsb_ee      = new RooPoisson( "pdf_Nsb_ee"      , "Nsb_ee Poisson PDF"      , *rv_Nsb_ee      , *rv_n_sb_ee ) ;
      pdf_Nsig_mm     = new RooPoisson( "pdf_Nsig_mm"     , "Nsig_mm Poisson PDF"     , *rv_Nsig_mm     , *rv_n_sig_mm ) ;
      pdf_Nsb_mm      = new RooPoisson( "pdf_Nsb_mm"      , "Nsb_mm Poisson PDF"      , *rv_Nsb_mm      , *rv_n_sb_mm ) ;


      {
         RooArgSet pdflist ;

         pdflist.add( *pdf_Nsig_1b        ) ;
         pdflist.add( *pdf_Nsb_1b         ) ;
         pdflist.add( *pdf_Nsig_ldp_1b    ) ;
         pdflist.add( *pdf_Nsb_ldp_1b     ) ;
         pdflist.add( *pdf_Nsig_sl_1b     ) ;
         pdflist.add( *pdf_Nsb_sl_1b      ) ;

         pdflist.add( *pdf_Nsig_2b        ) ;
         pdflist.add( *pdf_Nsb_2b         ) ;
         pdflist.add( *pdf_Nsig_ldp_2b    ) ;
         pdflist.add( *pdf_Nsb_ldp_2b     ) ;
         pdflist.add( *pdf_Nsig_sl_2b     ) ;
         pdflist.add( *pdf_Nsb_sl_2b      ) ;

         pdflist.add( *pdf_Nsig_3b        ) ;
         pdflist.add( *pdf_Nsb_3b         ) ;
         pdflist.add( *pdf_Nsig_ldp_3b    ) ;
         pdflist.add( *pdf_Nsb_ldp_3b     ) ;
         pdflist.add( *pdf_Nsig_sl_3b     ) ;
         pdflist.add( *pdf_Nsb_sl_3b      ) ;

         pdflist.add( *pdf_Nsig_ee     ) ;
         pdflist.add( *pdf_Nsb_ee      ) ;
         pdflist.add( *pdf_Nsig_mm     ) ;
         pdflist.add( *pdf_Nsb_mm      ) ;

         pdflist.add(allNuisancePdfs);

         likelihood = new RooProdPdf( "likelihood", "ra2b likelihood", pdflist ) ;
      }


     //---- Define the list of observables.

       observedParametersList.add( *rv_Nsig_1b        ) ;
       observedParametersList.add( *rv_Nsb_1b         ) ;
       observedParametersList.add( *rv_Nsig_sl_1b     ) ;
       observedParametersList.add( *rv_Nsb_sl_1b      ) ;
       observedParametersList.add( *rv_Nsig_ldp_1b    ) ;
       observedParametersList.add( *rv_Nsb_ldp_1b     ) ;

       observedParametersList.add( *rv_Nsig_2b        ) ;
       observedParametersList.add( *rv_Nsb_2b         ) ;
       observedParametersList.add( *rv_Nsig_sl_2b     ) ;
       observedParametersList.add( *rv_Nsb_sl_2b      ) ;
       observedParametersList.add( *rv_Nsig_ldp_2b    ) ;
       observedParametersList.add( *rv_Nsb_ldp_2b     ) ;

       observedParametersList.add( *rv_Nsig_3b        ) ;
       observedParametersList.add( *rv_Nsb_3b         ) ;
       observedParametersList.add( *rv_Nsig_sl_3b     ) ;
       observedParametersList.add( *rv_Nsb_sl_3b      ) ;
       observedParametersList.add( *rv_Nsig_ldp_3b    ) ;
       observedParametersList.add( *rv_Nsb_ldp_3b     ) ;

       observedParametersList.add( *rv_Nsb_ee      ) ;
       observedParametersList.add( *rv_Nsig_ee     ) ;
       observedParametersList.add( *rv_Nsb_mm      ) ;
       observedParametersList.add( *rv_Nsig_mm     ) ;



       dsObserved = new RooDataSet("ra2b_observed_rds", "RA2b observed data values",
                                      observedParametersList ) ;
       dsObserved->add( observedParametersList ) ;


       RooWorkspace workspace ("ws") ;

       workspace.import(*dsObserved);

       // parameters of interest
       RooArgSet poi(*rv_mu_susy_sig_1b, "poi");
       // flat prior for POI
       RooUniform signal_prior ("signal_prior","signal_prior",*rv_mu_susy_sig_1b);


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


    bool ra2bRoostatsClass9::setSusyScanPoint( const char* inputScanFile,
                                               double m0, double m12, bool isT1bbbb, double t1bbbbXsec,
                                               const char* inputSusy_deff_dbtageff_file
                                             ) {


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
          float n_sig_1b_raw ;
          float n_sb_1b_raw ;
          float n_sig_sl_1b_raw ;
          float n_sb_sl_1b_raw ;
          float n_sig_ldp_1b_raw ;
          float n_sb_ldp_1b_raw ;

          float n_sig_2b_raw ;
          float n_sb_2b_raw ;
          float n_sig_sl_2b_raw ;
          float n_sb_sl_2b_raw ;
          float n_sig_ldp_2b_raw ;
          float n_sb_ldp_2b_raw ;

          float n_sig_3b_raw ;
          float n_sb_3b_raw ;
          float n_sig_sl_3b_raw ;
          float n_sb_sl_3b_raw ;
          float n_sig_ldp_3b_raw ;
          float n_sb_ldp_3b_raw ;



          float n_sig_1b_correction ;
          float n_sb_1b_correction ;
          float n_sig_sl_1b_correction ;
          float n_sb_sl_1b_correction ;
          float n_sig_ldp_1b_correction ;
          float n_sb_ldp_1b_correction ;

          float n_sig_2b_correction ;
          float n_sb_2b_correction ;
          float n_sig_sl_2b_correction ;
          float n_sb_sl_2b_correction ;
          float n_sig_ldp_2b_correction ;
          float n_sb_ldp_2b_correction ;

          float n_sig_3b_correction ;
          float n_sb_3b_correction ;
          float n_sig_sl_3b_correction ;
          float n_sb_sl_3b_correction ;
          float n_sig_ldp_3b_correction ;
          float n_sb_ldp_3b_correction ;



          float n_sig_1b_error ;
          float n_sb_1b_error ;
          float n_sig_sl_1b_error ;
          float n_sb_sl_1b_error ;
          float n_sig_ldp_1b_error ;
          float n_sb_ldp_1b_error ;

          float n_sig_2b_error ;
          float n_sb_2b_error ;
          float n_sig_sl_2b_error ;
          float n_sb_sl_2b_error ;
          float n_sig_ldp_2b_error ;
          float n_sb_ldp_2b_error ;

          float n_sig_3b_error ;
          float n_sb_3b_error ;
          float n_sig_sl_3b_error ;
          float n_sb_sl_3b_error ;
          float n_sig_ldp_3b_error ;
          float n_sb_ldp_3b_error ;




          int nGen ;

          fscanf( infp, "%f %f %d  %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
            &pointM0, &pointM12, &nGen,
            &n_sig_1b_raw, &n_sb_1b_raw, &n_sig_sl_1b_raw, &n_sb_sl_1b_raw, &n_sig_ldp_1b_raw, &n_sb_ldp_1b_raw,
            &n_sig_2b_raw, &n_sb_2b_raw, &n_sig_sl_2b_raw, &n_sb_sl_2b_raw, &n_sig_ldp_2b_raw, &n_sb_ldp_2b_raw,
            &n_sig_3b_raw, &n_sb_3b_raw, &n_sig_sl_3b_raw, &n_sb_sl_3b_raw, &n_sig_ldp_3b_raw, &n_sb_ldp_3b_raw,
            &n_sig_1b_correction, &n_sb_1b_correction, &n_sig_sl_1b_correction, &n_sb_sl_1b_correction, &n_sig_ldp_1b_correction, &n_sb_ldp_1b_correction,
            &n_sig_2b_correction, &n_sb_2b_correction, &n_sig_sl_2b_correction, &n_sb_sl_2b_correction, &n_sig_ldp_2b_correction, &n_sb_ldp_2b_correction,
            &n_sig_3b_correction, &n_sb_3b_correction, &n_sig_sl_3b_correction, &n_sb_sl_3b_correction, &n_sig_ldp_3b_correction, &n_sb_ldp_3b_correction,
            &n_sig_1b_error, &n_sb_1b_error, &n_sig_sl_1b_error, &n_sb_sl_1b_error, &n_sig_ldp_1b_error, &n_sb_ldp_1b_error,
            &n_sig_2b_error, &n_sb_2b_error, &n_sig_sl_2b_error, &n_sb_sl_2b_error, &n_sig_ldp_2b_error, &n_sb_ldp_2b_error,
            &n_sig_3b_error, &n_sb_3b_error, &n_sig_sl_3b_error, &n_sb_sl_3b_error, &n_sig_ldp_3b_error, &n_sb_ldp_3b_error ) ;

          if ( feof(infp) ) break ;
          if ( n_sig_1b_raw < 0.00001 ) continue ;
          if ( n_sig_2b_raw < 0.00001 ) continue ;
          if ( n_sig_3b_raw < 0.00001 ) continue ;
   //--If you are asking for it, I'll assume it's good.  Josh is using 0 for ngen dummy in LM9.
   /////  if ( nGen != 10000 ) continue ; // get rid of bad scan points.

          printf("  pointM0 = %g , pointM12 = %g\n", pointM0, pointM12 ) ;
          if (    fabs( pointM0 - m0 ) <= deltaM0/2.
               && fabs( pointM12 - m12 ) <= deltaM12/2. ) {

              double nGenPerPoint = 10000 ; // for t1bbbb

             printf("\n\n") ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_1b_error     = %6.1f %%\n", n_sig_1b_error     ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_1b_error      = %6.1f %%\n", n_sb_1b_error      ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_sl_1b_error  = %6.1f %%\n", n_sig_sl_1b_error  ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_sl_1b_error   = %6.1f %%\n", n_sb_sl_1b_error   ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_ldp_1b_error = %6.1f %%\n", n_sig_ldp_1b_error ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_ldp_1b_error  = %6.1f %%\n", n_sb_ldp_1b_error  ) ;
             printf("\n\n") ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_2b_error     = %6.1f %%\n", n_sig_2b_error     ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_2b_error      = %6.1f %%\n", n_sb_2b_error      ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_sl_2b_error  = %6.1f %%\n", n_sig_sl_2b_error  ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_sl_2b_error   = %6.1f %%\n", n_sb_sl_2b_error   ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_ldp_2b_error = %6.1f %%\n", n_sig_ldp_2b_error ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_ldp_2b_error  = %6.1f %%\n", n_sb_ldp_2b_error  ) ;
             printf("\n\n") ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_3b_error     = %6.1f %%\n", n_sig_3b_error     ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_3b_error      = %6.1f %%\n", n_sb_3b_error      ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_sl_3b_error  = %6.1f %%\n", n_sig_sl_3b_error  ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_sl_3b_error   = %6.1f %%\n", n_sb_sl_3b_error   ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sig_ldp_3b_error = %6.1f %%\n", n_sig_ldp_3b_error ) ;
             printf("  SUSY efficiency  systematic uncertainty,   n_sb_ldp_3b_error  = %6.1f %%\n", n_sb_ldp_3b_error  ) ;
             printf("\n\n") ;

      //--Include the stat error on the efficiency for t1bbbb.
             if ( isT1bbbb ) {

              //-- absolute raw eff
                 float n_sig_1b_raw_eff     = n_sig_1b_raw     / nGenPerPoint ;
                 float n_sb_1b_raw_eff      = n_sb_1b_raw      / nGenPerPoint ;
                 float n_sig_sl_1b_raw_eff  = n_sig_sl_1b_raw  / nGenPerPoint ;
                 float n_sb_sl_1b_raw_eff   = n_sb_sl_1b_raw   / nGenPerPoint ;
                 float n_sig_ldp_1b_raw_eff = n_sig_ldp_1b_raw / nGenPerPoint ;
                 float n_sb_ldp_1b_raw_eff  = n_sb_ldp_1b_raw  / nGenPerPoint ;

                 float n_sig_2b_raw_eff     = n_sig_2b_raw     / nGenPerPoint ;
                 float n_sb_2b_raw_eff      = n_sb_2b_raw      / nGenPerPoint ;
                 float n_sig_sl_2b_raw_eff  = n_sig_sl_2b_raw  / nGenPerPoint ;
                 float n_sb_sl_2b_raw_eff   = n_sb_sl_2b_raw   / nGenPerPoint ;
                 float n_sig_ldp_2b_raw_eff = n_sig_ldp_2b_raw / nGenPerPoint ;
                 float n_sb_ldp_2b_raw_eff  = n_sb_ldp_2b_raw  / nGenPerPoint ;

                 float n_sig_3b_raw_eff     = n_sig_3b_raw     / nGenPerPoint ;
                 float n_sb_3b_raw_eff      = n_sb_3b_raw      / nGenPerPoint ;
                 float n_sig_sl_3b_raw_eff  = n_sig_sl_3b_raw  / nGenPerPoint ;
                 float n_sb_sl_3b_raw_eff   = n_sb_sl_3b_raw   / nGenPerPoint ;
                 float n_sig_ldp_3b_raw_eff = n_sig_ldp_3b_raw / nGenPerPoint ;
                 float n_sb_ldp_3b_raw_eff  = n_sb_ldp_3b_raw  / nGenPerPoint ;


               //-- absolute stat err.
                 float n_sig_1b_stat_error     =  sqrt(  n_sig_1b_raw_eff    * ( 1.0 - n_sig_1b_raw_eff     ) / nGenPerPoint ) ;
                 float n_sb_1b_stat_error      =  sqrt(  n_sb_1b_raw_eff     * ( 1.0 - n_sb_1b_raw_eff      ) / nGenPerPoint ) ;
                 float n_sig_sl_1b_stat_error  =  sqrt(  n_sig_sl_1b_raw_eff * ( 1.0 - n_sig_sl_1b_raw_eff  ) / nGenPerPoint ) ;
                 float n_sb_sl_1b_stat_error   =  sqrt(  n_sb_sl_1b_raw_eff  * ( 1.0 - n_sb_sl_1b_raw_eff   ) / nGenPerPoint ) ;
                 float n_sig_ldp_1b_stat_error =  sqrt(  n_sig_ldp_1b_raw_eff* ( 1.0 - n_sig_ldp_1b_raw_eff ) / nGenPerPoint ) ;
                 float n_sb_ldp_1b_stat_error  =  sqrt(  n_sb_ldp_1b_raw_eff * ( 1.0 - n_sb_ldp_1b_raw_eff  ) / nGenPerPoint ) ;

                 float n_sig_2b_stat_error     =  sqrt(  n_sig_2b_raw_eff    * ( 1.0 - n_sig_2b_raw_eff     ) / nGenPerPoint ) ;
                 float n_sb_2b_stat_error      =  sqrt(  n_sb_2b_raw_eff     * ( 1.0 - n_sb_2b_raw_eff      ) / nGenPerPoint ) ;
                 float n_sig_sl_2b_stat_error  =  sqrt(  n_sig_sl_2b_raw_eff * ( 1.0 - n_sig_sl_2b_raw_eff  ) / nGenPerPoint ) ;
                 float n_sb_sl_2b_stat_error   =  sqrt(  n_sb_sl_2b_raw_eff  * ( 1.0 - n_sb_sl_2b_raw_eff   ) / nGenPerPoint ) ;
                 float n_sig_ldp_2b_stat_error =  sqrt(  n_sig_ldp_2b_raw_eff* ( 1.0 - n_sig_ldp_2b_raw_eff ) / nGenPerPoint ) ;
                 float n_sb_ldp_2b_stat_error  =  sqrt(  n_sb_ldp_2b_raw_eff * ( 1.0 - n_sb_ldp_2b_raw_eff  ) / nGenPerPoint ) ;

                 float n_sig_3b_stat_error     =  sqrt(  n_sig_3b_raw_eff    * ( 1.0 - n_sig_3b_raw_eff     ) / nGenPerPoint ) ;
                 float n_sb_3b_stat_error      =  sqrt(  n_sb_3b_raw_eff     * ( 1.0 - n_sb_3b_raw_eff      ) / nGenPerPoint ) ;
                 float n_sig_sl_3b_stat_error  =  sqrt(  n_sig_sl_3b_raw_eff * ( 1.0 - n_sig_sl_3b_raw_eff  ) / nGenPerPoint ) ;
                 float n_sb_sl_3b_stat_error   =  sqrt(  n_sb_sl_3b_raw_eff  * ( 1.0 - n_sb_sl_3b_raw_eff   ) / nGenPerPoint ) ;
                 float n_sig_ldp_3b_stat_error =  sqrt(  n_sig_ldp_3b_raw_eff* ( 1.0 - n_sig_ldp_3b_raw_eff ) / nGenPerPoint ) ;
                 float n_sb_ldp_3b_stat_error  =  sqrt(  n_sb_ldp_3b_raw_eff * ( 1.0 - n_sb_ldp_3b_raw_eff  ) / nGenPerPoint ) ;

               //-- relative stat err in percent.
                 if ( n_sig_1b_raw_eff     > 0 ) { n_sig_1b_stat_error     = 100.* n_sig_1b_stat_error     / n_sig_1b_raw_eff     ; } else { n_sig_1b_stat_error     = 0. ; }
                 if ( n_sb_1b_raw_eff      > 0 ) { n_sb_1b_stat_error      = 100.* n_sb_1b_stat_error      / n_sb_1b_raw_eff      ; } else { n_sb_1b_stat_error      = 0. ; }
                 if ( n_sig_sl_1b_raw_eff  > 0 ) { n_sig_sl_1b_stat_error  = 100.* n_sig_sl_1b_stat_error  / n_sig_sl_1b_raw_eff  ; } else { n_sig_sl_1b_stat_error  = 0. ; }
                 if ( n_sb_sl_1b_raw_eff   > 0 ) { n_sb_sl_1b_stat_error   = 100.* n_sb_sl_1b_stat_error   / n_sb_sl_1b_raw_eff   ; } else { n_sb_sl_1b_stat_error   = 0. ; }
                 if ( n_sig_ldp_1b_raw_eff > 0 ) { n_sig_ldp_1b_stat_error = 100.* n_sig_ldp_1b_stat_error / n_sig_ldp_1b_raw_eff ; } else { n_sig_ldp_1b_stat_error = 0. ; }
                 if ( n_sb_ldp_1b_raw_eff  > 0 ) { n_sb_ldp_1b_stat_error  = 100.* n_sb_ldp_1b_stat_error  / n_sb_ldp_1b_raw_eff  ; } else { n_sb_ldp_1b_stat_error  = 0. ; }

                 if ( n_sig_2b_raw_eff     > 0 ) { n_sig_2b_stat_error     = 100.* n_sig_2b_stat_error     / n_sig_2b_raw_eff     ; } else { n_sig_2b_stat_error     = 0. ; }
                 if ( n_sb_2b_raw_eff      > 0 ) { n_sb_2b_stat_error      = 100.* n_sb_2b_stat_error      / n_sb_2b_raw_eff      ; } else { n_sb_2b_stat_error      = 0. ; }
                 if ( n_sig_sl_2b_raw_eff  > 0 ) { n_sig_sl_2b_stat_error  = 100.* n_sig_sl_2b_stat_error  / n_sig_sl_2b_raw_eff  ; } else { n_sig_sl_2b_stat_error  = 0. ; }
                 if ( n_sb_sl_2b_raw_eff   > 0 ) { n_sb_sl_2b_stat_error   = 100.* n_sb_sl_2b_stat_error   / n_sb_sl_2b_raw_eff   ; } else { n_sb_sl_2b_stat_error   = 0. ; }
                 if ( n_sig_ldp_2b_raw_eff > 0 ) { n_sig_ldp_2b_stat_error = 100.* n_sig_ldp_2b_stat_error / n_sig_ldp_2b_raw_eff ; } else { n_sig_ldp_2b_stat_error = 0. ; }
                 if ( n_sb_ldp_2b_raw_eff  > 0 ) { n_sb_ldp_2b_stat_error  = 100.* n_sb_ldp_2b_stat_error  / n_sb_ldp_2b_raw_eff  ; } else { n_sb_ldp_2b_stat_error  = 0. ; }

                 if ( n_sig_3b_raw_eff     > 0 ) { n_sig_3b_stat_error     = 100.* n_sig_3b_stat_error     / n_sig_3b_raw_eff     ; } else { n_sig_3b_stat_error     = 0. ; }
                 if ( n_sb_3b_raw_eff      > 0 ) { n_sb_3b_stat_error      = 100.* n_sb_3b_stat_error      / n_sb_3b_raw_eff      ; } else { n_sb_3b_stat_error      = 0. ; }
                 if ( n_sig_sl_3b_raw_eff  > 0 ) { n_sig_sl_3b_stat_error  = 100.* n_sig_sl_3b_stat_error  / n_sig_sl_3b_raw_eff  ; } else { n_sig_sl_3b_stat_error  = 0. ; }
                 if ( n_sb_sl_3b_raw_eff   > 0 ) { n_sb_sl_3b_stat_error   = 100.* n_sb_sl_3b_stat_error   / n_sb_sl_3b_raw_eff   ; } else { n_sb_sl_3b_stat_error   = 0. ; }
                 if ( n_sig_ldp_3b_raw_eff > 0 ) { n_sig_ldp_3b_stat_error = 100.* n_sig_ldp_3b_stat_error / n_sig_ldp_3b_raw_eff ; } else { n_sig_ldp_3b_stat_error = 0. ; }
                 if ( n_sb_ldp_3b_raw_eff  > 0 ) { n_sb_ldp_3b_stat_error  = 100.* n_sb_ldp_3b_stat_error  / n_sb_ldp_3b_raw_eff  ; } else { n_sb_ldp_3b_stat_error  = 0. ; }

                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_1b_stat_error     = %6.1f %%\n", n_sig_1b_stat_error     ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_1b_stat_error      = %6.1f %%\n", n_sb_1b_stat_error      ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_sl_1b_stat_error  = %6.1f %%\n", n_sig_sl_1b_stat_error  ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_sl_1b_stat_error   = %6.1f %%\n", n_sb_sl_1b_stat_error   ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_ldp_1b_stat_error = %6.1f %%\n", n_sig_ldp_1b_stat_error ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_ldp_1b_stat_error  = %6.1f %%\n", n_sb_ldp_1b_stat_error  ) ;

                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_2b_stat_error     = %6.1f %%\n", n_sig_2b_stat_error     ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_2b_stat_error      = %6.1f %%\n", n_sb_2b_stat_error      ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_sl_2b_stat_error  = %6.1f %%\n", n_sig_sl_2b_stat_error  ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_sl_2b_stat_error   = %6.1f %%\n", n_sb_sl_2b_stat_error   ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_ldp_2b_stat_error = %6.1f %%\n", n_sig_ldp_2b_stat_error ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_ldp_2b_stat_error  = %6.1f %%\n", n_sb_ldp_2b_stat_error  ) ;

                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_3b_stat_error     = %6.1f %%\n", n_sig_3b_stat_error     ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_3b_stat_error      = %6.1f %%\n", n_sb_3b_stat_error      ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_sl_3b_stat_error  = %6.1f %%\n", n_sig_sl_3b_stat_error  ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_sl_3b_stat_error   = %6.1f %%\n", n_sb_sl_3b_stat_error   ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sig_ldp_3b_stat_error = %6.1f %%\n", n_sig_ldp_3b_stat_error ) ;
                 printf("  SUSY efficiency  statistical uncertainty,   n_sb_ldp_3b_stat_error  = %6.1f %%\n", n_sb_ldp_3b_stat_error  ) ;

               //-- total err in percent.
                 n_sig_1b_error     = sqrt( pow( n_sig_1b_error    , 2) + pow( n_sig_1b_stat_error    , 2) ) ;
                 n_sb_1b_error      = sqrt( pow( n_sb_1b_error     , 2) + pow( n_sb_1b_stat_error     , 2) ) ;
                 n_sig_sl_1b_error  = sqrt( pow( n_sig_sl_1b_error , 2) + pow( n_sig_sl_1b_stat_error , 2) ) ;
                 n_sb_sl_1b_error   = sqrt( pow( n_sb_sl_1b_error  , 2) + pow( n_sb_sl_1b_stat_error  , 2) ) ;
                 n_sig_ldp_1b_error = sqrt( pow( n_sig_ldp_1b_error, 2) + pow( n_sig_ldp_1b_stat_error, 2) ) ;
                 n_sb_ldp_1b_error  = sqrt( pow( n_sb_ldp_1b_error , 2) + pow( n_sb_ldp_1b_stat_error , 2) ) ;

                 n_sig_2b_error     = sqrt( pow( n_sig_2b_error    , 2) + pow( n_sig_2b_stat_error    , 2) ) ;
                 n_sb_2b_error      = sqrt( pow( n_sb_2b_error     , 2) + pow( n_sb_2b_stat_error     , 2) ) ;
                 n_sig_sl_2b_error  = sqrt( pow( n_sig_sl_2b_error , 2) + pow( n_sig_sl_2b_stat_error , 2) ) ;
                 n_sb_sl_2b_error   = sqrt( pow( n_sb_sl_2b_error  , 2) + pow( n_sb_sl_2b_stat_error  , 2) ) ;
                 n_sig_ldp_2b_error = sqrt( pow( n_sig_ldp_2b_error, 2) + pow( n_sig_ldp_2b_stat_error, 2) ) ;
                 n_sb_ldp_2b_error  = sqrt( pow( n_sb_ldp_2b_error , 2) + pow( n_sb_ldp_2b_stat_error , 2) ) ;

                 n_sig_3b_error     = sqrt( pow( n_sig_3b_error    , 2) + pow( n_sig_3b_stat_error    , 2) ) ;
                 n_sb_3b_error      = sqrt( pow( n_sb_3b_error     , 2) + pow( n_sb_3b_stat_error     , 2) ) ;
                 n_sig_sl_3b_error  = sqrt( pow( n_sig_sl_3b_error , 2) + pow( n_sig_sl_3b_stat_error , 2) ) ;
                 n_sb_sl_3b_error   = sqrt( pow( n_sb_sl_3b_error  , 2) + pow( n_sb_sl_3b_stat_error  , 2) ) ;
                 n_sig_ldp_3b_error = sqrt( pow( n_sig_ldp_3b_error, 2) + pow( n_sig_ldp_3b_stat_error, 2) ) ;
                 n_sb_ldp_3b_error  = sqrt( pow( n_sb_ldp_3b_error , 2) + pow( n_sb_ldp_3b_stat_error , 2) ) ;

                 printf("\n\n") ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_1b_error     = %6.1f %%\n", n_sig_1b_error     ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_1b_error      = %6.1f %%\n", n_sb_1b_error      ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_sl_1b_error  = %6.1f %%\n", n_sig_sl_1b_error  ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_sl_1b_error   = %6.1f %%\n", n_sb_sl_1b_error   ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_ldp_1b_error = %6.1f %%\n", n_sig_ldp_1b_error ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_ldp_1b_error  = %6.1f %%\n", n_sb_ldp_1b_error  ) ;
                 printf("\n\n") ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_2b_error     = %6.1f %%\n", n_sig_2b_error     ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_2b_error      = %6.1f %%\n", n_sb_2b_error      ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_sl_2b_error  = %6.1f %%\n", n_sig_sl_2b_error  ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_sl_2b_error   = %6.1f %%\n", n_sb_sl_2b_error   ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_ldp_2b_error = %6.1f %%\n", n_sig_ldp_2b_error ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_ldp_2b_error  = %6.1f %%\n", n_sb_ldp_2b_error  ) ;
                 printf("\n\n") ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_3b_error     = %6.1f %%\n", n_sig_3b_error     ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_3b_error      = %6.1f %%\n", n_sb_3b_error      ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_sl_3b_error  = %6.1f %%\n", n_sig_sl_3b_error  ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_sl_3b_error   = %6.1f %%\n", n_sb_sl_3b_error   ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sig_ldp_3b_error = %6.1f %%\n", n_sig_ldp_3b_error ) ;
                 printf("  SUSY efficiency  total uncertainty,   n_sb_ldp_3b_error  = %6.1f %%\n", n_sb_ldp_3b_error  ) ;
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



             double setVal_n_sig_1b(0.) ;
             double setVal_n_sb_1b(0.) ;
             double setVal_n_sig_sl_1b(0.) ;
             double setVal_n_sb_sl_1b(0.) ;
             double setVal_n_sig_ldp_1b(0.) ;
             double setVal_n_sb_ldp_1b(0.) ;

             double setVal_n_sig_2b(0.) ;
             double setVal_n_sb_2b(0.) ;
             double setVal_n_sig_sl_2b(0.) ;
             double setVal_n_sb_sl_2b(0.) ;
             double setVal_n_sig_ldp_2b(0.) ;
             double setVal_n_sb_ldp_2b(0.) ;

             double setVal_n_sig_3b(0.) ;
             double setVal_n_sb_3b(0.) ;
             double setVal_n_sig_sl_3b(0.) ;
             double setVal_n_sb_sl_3b(0.) ;
             double setVal_n_sig_ldp_3b(0.) ;
             double setVal_n_sb_ldp_3b(0.) ;


             if ( !isT1bbbb ) {

                //-- tanb40
                setVal_n_sig_1b     = n_sig_1b_raw     * n_sig_1b_correction     ;
                setVal_n_sb_1b      = n_sb_1b_raw      * n_sb_1b_correction      ;
                setVal_n_sig_sl_1b  = n_sig_sl_1b_raw  * n_sig_sl_1b_correction  ;
                setVal_n_sb_sl_1b   = n_sb_sl_1b_raw   * n_sb_sl_1b_correction   ;
                setVal_n_sig_ldp_1b = n_sig_ldp_1b_raw * n_sig_ldp_1b_correction ;
                setVal_n_sb_ldp_1b  = n_sb_ldp_1b_raw  * n_sb_ldp_1b_correction  ;

                setVal_n_sig_2b     = n_sig_2b_raw     * n_sig_2b_correction     ;
                setVal_n_sb_2b      = n_sb_2b_raw      * n_sb_2b_correction      ;
                setVal_n_sig_sl_2b  = n_sig_sl_2b_raw  * n_sig_sl_2b_correction  ;
                setVal_n_sb_sl_2b   = n_sb_sl_2b_raw   * n_sb_sl_2b_correction   ;
                setVal_n_sig_ldp_2b = n_sig_ldp_2b_raw * n_sig_ldp_2b_correction ;
                setVal_n_sb_ldp_2b  = n_sb_ldp_2b_raw  * n_sb_ldp_2b_correction  ;

                setVal_n_sig_3b     = n_sig_3b_raw     * n_sig_3b_correction     ;
                setVal_n_sb_3b      = n_sb_3b_raw      * n_sb_3b_correction      ;
                setVal_n_sig_sl_3b  = n_sig_sl_3b_raw  * n_sig_sl_3b_correction  ;
                setVal_n_sb_sl_3b   = n_sb_sl_3b_raw   * n_sb_sl_3b_correction   ;
                setVal_n_sig_ldp_3b = n_sig_ldp_3b_raw * n_sig_ldp_3b_correction ;
                setVal_n_sb_ldp_3b  = n_sb_ldp_3b_raw  * n_sb_ldp_3b_correction  ;

             } else {

                //-- t1bbbb
                setVal_n_sig_1b     = DataLumi * t1bbbbXsec * (( n_sig_1b_raw     * n_sig_1b_correction     )/ nGenPerPoint ) ;
                setVal_n_sb_1b      = DataLumi * t1bbbbXsec * (( n_sb_1b_raw      * n_sb_1b_correction      )/ nGenPerPoint ) ;
                setVal_n_sig_sl_1b  = DataLumi * t1bbbbXsec * (( n_sig_sl_1b_raw  * n_sig_sl_1b_correction  )/ nGenPerPoint ) ;
                setVal_n_sb_sl_1b   = DataLumi * t1bbbbXsec * (( n_sb_sl_1b_raw   * n_sb_sl_1b_correction   )/ nGenPerPoint ) ;
                setVal_n_sig_ldp_1b = DataLumi * t1bbbbXsec * (( n_sig_ldp_1b_raw * n_sig_ldp_1b_correction )/ nGenPerPoint ) ;
                setVal_n_sb_ldp_1b  = DataLumi * t1bbbbXsec * (( n_sb_ldp_1b_raw  * n_sb_ldp_1b_correction  )/ nGenPerPoint ) ;

                setVal_n_sig_2b     = DataLumi * t1bbbbXsec * (( n_sig_2b_raw     * n_sig_2b_correction     )/ nGenPerPoint ) ;
                setVal_n_sb_2b      = DataLumi * t1bbbbXsec * (( n_sb_2b_raw      * n_sb_2b_correction      )/ nGenPerPoint ) ;
                setVal_n_sig_sl_2b  = DataLumi * t1bbbbXsec * (( n_sig_sl_2b_raw  * n_sig_sl_2b_correction  )/ nGenPerPoint ) ;
                setVal_n_sb_sl_2b   = DataLumi * t1bbbbXsec * (( n_sb_sl_2b_raw   * n_sb_sl_2b_correction   )/ nGenPerPoint ) ;
                setVal_n_sig_ldp_2b = DataLumi * t1bbbbXsec * (( n_sig_ldp_2b_raw * n_sig_ldp_2b_correction )/ nGenPerPoint ) ;
                setVal_n_sb_ldp_2b  = DataLumi * t1bbbbXsec * (( n_sb_ldp_2b_raw  * n_sb_ldp_2b_correction  )/ nGenPerPoint ) ;

                setVal_n_sig_3b     = DataLumi * t1bbbbXsec * (( n_sig_3b_raw     * n_sig_3b_correction     )/ nGenPerPoint ) ;
                setVal_n_sb_3b      = DataLumi * t1bbbbXsec * (( n_sb_3b_raw      * n_sb_3b_correction      )/ nGenPerPoint ) ;
                setVal_n_sig_sl_3b  = DataLumi * t1bbbbXsec * (( n_sig_sl_3b_raw  * n_sig_sl_3b_correction  )/ nGenPerPoint ) ;
                setVal_n_sb_sl_3b   = DataLumi * t1bbbbXsec * (( n_sb_sl_3b_raw   * n_sb_sl_3b_correction   )/ nGenPerPoint ) ;
                setVal_n_sig_ldp_3b = DataLumi * t1bbbbXsec * (( n_sig_ldp_3b_raw * n_sig_ldp_3b_correction )/ nGenPerPoint ) ;
                setVal_n_sb_ldp_3b  = DataLumi * t1bbbbXsec * (( n_sb_ldp_3b_raw  * n_sb_ldp_3b_correction  )/ nGenPerPoint ) ;

             }




             rv_mu_susymc_sig_1b       -> setVal( setVal_n_sig_1b      ) ;
             rv_mu_susymc_sb_1b        -> setVal( setVal_n_sb_1b       ) ;
             rv_mu_susymc_sig_sl_1b    -> setVal( setVal_n_sig_sl_1b   ) ;
             rv_mu_susymc_sb_sl_1b     -> setVal( setVal_n_sb_sl_1b    ) ;
             rv_mu_susymc_sig_ldp_1b   -> setVal( setVal_n_sig_ldp_1b  ) ;
             rv_mu_susymc_sb_ldp_1b    -> setVal( setVal_n_sb_ldp_1b   ) ;

             rv_mu_susymc_sig_2b       -> setVal( setVal_n_sig_2b      ) ;
             rv_mu_susymc_sb_2b        -> setVal( setVal_n_sb_2b       ) ;
             rv_mu_susymc_sig_sl_2b    -> setVal( setVal_n_sig_sl_2b   ) ;
             rv_mu_susymc_sb_sl_2b     -> setVal( setVal_n_sb_sl_2b    ) ;
             rv_mu_susymc_sig_ldp_2b   -> setVal( setVal_n_sig_ldp_2b  ) ;
             rv_mu_susymc_sb_ldp_2b    -> setVal( setVal_n_sb_ldp_2b   ) ;

             rv_mu_susymc_sig_3b       -> setVal( setVal_n_sig_3b      ) ;
             rv_mu_susymc_sb_3b        -> setVal( setVal_n_sb_3b       ) ;
             rv_mu_susymc_sig_sl_3b    -> setVal( setVal_n_sig_sl_3b   ) ;
             rv_mu_susymc_sb_sl_3b     -> setVal( setVal_n_sb_sl_3b    ) ;
             rv_mu_susymc_sig_ldp_3b   -> setVal( setVal_n_sig_ldp_3b  ) ;
             rv_mu_susymc_sb_ldp_3b    -> setVal( setVal_n_sb_ldp_3b   ) ;




             rv_width_eff_sf_sig_1b     -> setVal( n_sig_1b_error     / 100. ) ;
             rv_width_eff_sf_sb_1b      -> setVal( n_sb_1b_error      / 100. ) ;
             rv_width_eff_sf_sig_sl_1b  -> setVal( n_sig_sl_1b_error  / 100. ) ;
             rv_width_eff_sf_sb_sl_1b   -> setVal( n_sb_sl_1b_error   / 100. ) ;
             rv_width_eff_sf_sig_ldp_1b -> setVal( n_sig_ldp_1b_error / 100. ) ;
             rv_width_eff_sf_sb_ldp_1b  -> setVal( n_sb_ldp_1b_error  / 100. ) ;

             rv_width_eff_sf_sig_2b     -> setVal( n_sig_2b_error     / 100. ) ;
             rv_width_eff_sf_sb_2b      -> setVal( n_sb_2b_error      / 100. ) ;
             rv_width_eff_sf_sig_sl_2b  -> setVal( n_sig_sl_2b_error  / 100. ) ;
             rv_width_eff_sf_sb_sl_2b   -> setVal( n_sb_sl_2b_error   / 100. ) ;
             rv_width_eff_sf_sig_ldp_2b -> setVal( n_sig_ldp_2b_error / 100. ) ;
             rv_width_eff_sf_sb_ldp_2b  -> setVal( n_sb_ldp_2b_error  / 100. ) ;

             rv_width_eff_sf_sig_3b     -> setVal( n_sig_3b_error     / 100. ) ;
             rv_width_eff_sf_sb_3b      -> setVal( n_sb_3b_error      / 100. ) ;
             rv_width_eff_sf_sig_sl_3b  -> setVal( n_sig_sl_3b_error  / 100. ) ;
             rv_width_eff_sf_sb_sl_3b   -> setVal( n_sb_sl_3b_error   / 100. ) ;
             rv_width_eff_sf_sig_ldp_3b -> setVal( n_sig_ldp_3b_error / 100. ) ;
             rv_width_eff_sf_sb_ldp_3b  -> setVal( n_sb_ldp_3b_error  / 100. ) ;





             if ( !isT1bbbb ) {
                printf("\n\n Found point m0 = %4.0f,  m1/2 = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, setVal_n_sig_1b ) ;
             } else {
                printf("\n\n Found point mGluino = %4.0f,  mLSP = %4.0f,  Npred = %7.1f\n\n\n", pointM0, pointM12, setVal_n_sig_1b ) ;
             }


             printf("\n\n") ;
             printf(" Setting susy N_sig_1b     to  %7.1f\n", setVal_n_sig_1b       ) ;
             printf(" Setting susy N_sb_1b      to  %7.1f\n", setVal_n_sb_1b        ) ;
             printf(" Setting susy N_sig_sl_1b  to  %7.1f\n", setVal_n_sig_sl_1b    ) ;
             printf(" Setting susy N_sb_sl_1b   to  %7.1f\n", setVal_n_sb_sl_1b     ) ;
             printf(" Setting susy N_sig_ldp_1b to  %7.1f\n", setVal_n_sig_ldp_1b   ) ;
             printf(" Setting susy N_sb_ldp_1b  to  %7.1f\n", setVal_n_sb_ldp_1b    ) ;
             printf("\n\n") ;
             printf(" Setting susy N_sig_2b     to  %7.1f\n", setVal_n_sig_2b       ) ;
             printf(" Setting susy N_sb_2b      to  %7.1f\n", setVal_n_sb_2b        ) ;
             printf(" Setting susy N_sig_sl_2b  to  %7.1f\n", setVal_n_sig_sl_2b    ) ;
             printf(" Setting susy N_sb_sl_2b   to  %7.1f\n", setVal_n_sb_sl_2b     ) ;
             printf(" Setting susy N_sig_ldp_2b to  %7.1f\n", setVal_n_sig_ldp_2b   ) ;
             printf(" Setting susy N_sb_ldp_2b  to  %7.1f\n", setVal_n_sb_ldp_2b    ) ;
             printf("\n\n") ;
             printf(" Setting susy N_sig_3b     to  %7.1f\n", setVal_n_sig_3b       ) ;
             printf(" Setting susy N_sb_3b      to  %7.1f\n", setVal_n_sb_3b        ) ;
             printf(" Setting susy N_sig_sl_3b  to  %7.1f\n", setVal_n_sig_sl_3b    ) ;
             printf(" Setting susy N_sb_sl_3b   to  %7.1f\n", setVal_n_sb_sl_3b     ) ;
             printf(" Setting susy N_sig_ldp_3b to  %7.1f\n", setVal_n_sig_ldp_3b   ) ;
             printf(" Setting susy N_sb_ldp_3b  to  %7.1f\n", setVal_n_sb_ldp_3b    ) ;
             printf("\n\n") ;

             found = true ;

             break ;

          } // point match?

       } // not eof ?

       fclose( infp ) ;

       if ( !found ) {
          printf("\n\n *** Point not found in scan.\n\n" ) ;
          return false ;
       }






       //----- Now, read in the deff_dbtageff numbers.


       printf("\n\n Opening SUSY deff_dbtageff input file : %s\n", inputSusy_deff_dbtageff_file ) ;

       if ( (infp=fopen( inputSusy_deff_dbtageff_file,"r"))==NULL ) {
          printf("\n\n *** Problem opening input file: %s.\n\n", inputSusy_deff_dbtageff_file ) ;
          return false ;
       }

       while ( !feof( infp ) ) {

          float pointM0, pointM12 ;

          float  deff_dbtageff_sig_1b, deff_dbtageff_sb_1b, deff_dbtageff_sig_sl_1b, deff_dbtageff_sb_sl_1b, deff_dbtageff_sig_ldp_1b, deff_dbtageff_sb_ldp_1b ;
          float  deff_dbtageff_sig_2b, deff_dbtageff_sb_2b, deff_dbtageff_sig_sl_2b, deff_dbtageff_sb_sl_2b, deff_dbtageff_sig_ldp_2b, deff_dbtageff_sb_ldp_2b ;
          float  deff_dbtageff_sig_3b, deff_dbtageff_sb_3b, deff_dbtageff_sig_sl_3b, deff_dbtageff_sb_sl_3b, deff_dbtageff_sig_ldp_3b, deff_dbtageff_sb_ldp_3b ;

          fscanf( infp, "%f %f   %f %f %f %f %f %f   %f %f %f %f %f %f   %f %f %f %f %f %f",
            &pointM0, &pointM12,
            &deff_dbtageff_sig_1b, &deff_dbtageff_sb_1b, &deff_dbtageff_sig_sl_1b, &deff_dbtageff_sb_sl_1b, &deff_dbtageff_sig_ldp_1b, &deff_dbtageff_sb_ldp_1b,
            &deff_dbtageff_sig_2b, &deff_dbtageff_sb_2b, &deff_dbtageff_sig_sl_2b, &deff_dbtageff_sb_sl_2b, &deff_dbtageff_sig_ldp_2b, &deff_dbtageff_sb_ldp_2b,
            &deff_dbtageff_sig_3b, &deff_dbtageff_sb_3b, &deff_dbtageff_sig_sl_3b, &deff_dbtageff_sb_sl_3b, &deff_dbtageff_sig_ldp_3b, &deff_dbtageff_sb_ldp_3b
          ) ;

          printf("  pointM0 = %g , pointM12 = %g\n", pointM0, pointM12 ) ;
          if (    fabs( pointM0 - m0 ) <= deltaM0/2.
               && fabs( pointM12 - m12 ) <= deltaM12/2. ) {

             printf("\n\n") ;
             printf(" Setting susy deff_dbtageff_sig_1b     to  %7.1f\n", deff_dbtageff_sig_1b       ) ;
             printf(" Setting susy deff_dbtageff_sb_1b      to  %7.1f\n", deff_dbtageff_sb_1b        ) ;
             printf(" Setting susy deff_dbtageff_sig_sl_1b  to  %7.1f\n", deff_dbtageff_sig_sl_1b    ) ;
             printf(" Setting susy deff_dbtageff_sb_sl_1b   to  %7.1f\n", deff_dbtageff_sb_sl_1b     ) ;
             printf(" Setting susy deff_dbtageff_sig_ldp_1b to  %7.1f\n", deff_dbtageff_sig_ldp_1b   ) ;
             printf(" Setting susy deff_dbtageff_sb_ldp_1b  to  %7.1f\n", deff_dbtageff_sb_ldp_1b    ) ;
             printf("\n\n") ;
             printf(" Setting susy deff_dbtageff_sig_2b     to  %7.1f\n", deff_dbtageff_sig_2b       ) ;
             printf(" Setting susy deff_dbtageff_sb_2b      to  %7.1f\n", deff_dbtageff_sb_2b        ) ;
             printf(" Setting susy deff_dbtageff_sig_sl_2b  to  %7.1f\n", deff_dbtageff_sig_sl_2b    ) ;
             printf(" Setting susy deff_dbtageff_sb_sl_2b   to  %7.1f\n", deff_dbtageff_sb_sl_2b     ) ;
             printf(" Setting susy deff_dbtageff_sig_ldp_2b to  %7.1f\n", deff_dbtageff_sig_ldp_2b   ) ;
             printf(" Setting susy deff_dbtageff_sb_ldp_2b  to  %7.1f\n", deff_dbtageff_sb_ldp_2b    ) ;
             printf("\n\n") ;
             printf(" Setting susy deff_dbtageff_sig_3b     to  %7.1f\n", deff_dbtageff_sig_3b       ) ;
             printf(" Setting susy deff_dbtageff_sb_3b      to  %7.1f\n", deff_dbtageff_sb_3b        ) ;
             printf(" Setting susy deff_dbtageff_sig_sl_3b  to  %7.1f\n", deff_dbtageff_sig_sl_3b    ) ;
             printf(" Setting susy deff_dbtageff_sb_sl_3b   to  %7.1f\n", deff_dbtageff_sb_sl_3b     ) ;
             printf(" Setting susy deff_dbtageff_sig_ldp_3b to  %7.1f\n", deff_dbtageff_sig_ldp_3b   ) ;
             printf(" Setting susy deff_dbtageff_sb_ldp_3b  to  %7.1f\n", deff_dbtageff_sb_ldp_3b    ) ;
             printf("\n\n") ;

             rv_deff_dbtageff_sig_1b      -> setVal( deff_dbtageff_sig_1b       ) ;
             rv_deff_dbtageff_sb_1b       -> setVal( deff_dbtageff_sb_1b        ) ;
             rv_deff_dbtageff_sig_sl_1b   -> setVal( deff_dbtageff_sig_sl_1b    ) ;
             rv_deff_dbtageff_sb_sl_1b    -> setVal( deff_dbtageff_sb_sl_1b     ) ;
             rv_deff_dbtageff_sig_ldp_1b  -> setVal( deff_dbtageff_sig_ldp_1b   ) ;
             rv_deff_dbtageff_sb_ldp_1b   -> setVal( deff_dbtageff_sb_ldp_1b    ) ;

             rv_deff_dbtageff_sig_2b      -> setVal( deff_dbtageff_sig_2b       ) ;
             rv_deff_dbtageff_sb_2b       -> setVal( deff_dbtageff_sb_2b        ) ;
             rv_deff_dbtageff_sig_sl_2b   -> setVal( deff_dbtageff_sig_sl_2b    ) ;
             rv_deff_dbtageff_sb_sl_2b    -> setVal( deff_dbtageff_sb_sl_2b     ) ;
             rv_deff_dbtageff_sig_ldp_2b  -> setVal( deff_dbtageff_sig_ldp_2b   ) ;
             rv_deff_dbtageff_sb_ldp_2b   -> setVal( deff_dbtageff_sb_ldp_2b    ) ;

             rv_deff_dbtageff_sig_3b      -> setVal( deff_dbtageff_sig_3b       ) ;
             rv_deff_dbtageff_sb_3b       -> setVal( deff_dbtageff_sb_3b        ) ;
             rv_deff_dbtageff_sig_sl_3b   -> setVal( deff_dbtageff_sig_sl_3b    ) ;
             rv_deff_dbtageff_sb_sl_3b    -> setVal( deff_dbtageff_sb_sl_3b     ) ;
             rv_deff_dbtageff_sig_ldp_3b  -> setVal( deff_dbtageff_sig_ldp_3b   ) ;
             rv_deff_dbtageff_sb_ldp_3b   -> setVal( deff_dbtageff_sb_ldp_3b    ) ;



             found = true ;

             break ;

          }

          if ( feof(infp) ) break ;

       } // not eof



       if ( found ) {
          return true ;
       } else {
          printf("\n\n *** Point not found in scan.\n\n" ) ;
          return false ;
       }




    } // setSusyScanPoint.
